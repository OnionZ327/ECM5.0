/* The copyright in this software is being made available under the BSD
* License, included below. This software may be subject to other third party
* and contributor rights, including patent rights, and no such rights are
* granted under this license.
*
* Copyright (c) 2010-2022, ITU/ISO/IEC
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
*  * Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*/

/** \file     Buffer.cpp
 *  \brief    Low-overhead class describing 2D memory layout
 */

#define DONT_UNDEF_SIZE_AWARE_PER_EL_OP

// unit needs to come first due to a forward declaration
#include "Unit.h"
#include "Buffer.h"
#include "InterpolationFilter.h"

void applyPROFCore(Pel* dst, int dstStride, const Pel* src, int srcStride, int width, int height, const Pel* gradX, const Pel* gradY, int gradStride, const int* dMvX, const int* dMvY, int dMvStride, const bool& bi, int shiftNum, Pel offset, const ClpRng& clpRng)
{
  int idx = 0;

  const int dILimit = 1 << std::max<int>(clpRng.bd + 1, 13);
  for (int h = 0; h < height; h++)
  {
    for (int w = 0; w < width; w++)
    {
      int32_t dI = dMvX[idx] * gradX[w] + dMvY[idx] * gradY[w];
      dI = Clip3(-dILimit, dILimit - 1, dI);
      dst[w] = src[w] + dI;
      if (!bi)
      {
        dst[w] = (dst[w] + offset) >> shiftNum;
        dst[w] = ClipPel(dst[w], clpRng);
      }

      idx++;
    }
    gradX += gradStride;
    gradY += gradStride;
    dst += dstStride;
    src += srcStride;
  }
}

#if TM_AMVP || TM_MRG
int64_t getSumOfDifferenceCore(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, int width, int height, int rowSubShift, int bitDepth)
{
  height     >>= rowSubShift;
  src0Stride <<= rowSubShift;
  src1Stride <<= rowSubShift;

  int64_t sum = 0;
#define GET_SUM_DIFF_CORE_OP( ADDR ) sum += ( src0[ADDR] - src1[ADDR] )
#define GET_SUM_DIFF_CORE_INC    \
  src0 += src0Stride;            \
  src1 += src1Stride;            \

  SIZE_AWARE_PER_EL_OP(GET_SUM_DIFF_CORE_OP, GET_SUM_DIFF_CORE_INC);

#undef GET_SUM_DIFF_CORE_OP
#undef GET_SUM_DIFF_CORE_INC

  return sum;
}
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
void getAbsoluteDifferencePerSampleCore(Pel* dst, int dstStride, const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, int width, int height)
{
#define GET_ABS_DIFF_PER_SAMPLE_CORE_OP( ADDR ) dst[ADDR] = std::abs( src0[ADDR] - src1[ADDR] )
#define GET_ABS_DIFF_PER_SAMPLE_CORE_INC    \
  src0 += src0Stride;                       \
  src1 += src1Stride;                       \
  dst  += dstStride;                        \

  SIZE_AWARE_PER_EL_OP(GET_ABS_DIFF_PER_SAMPLE_CORE_OP, GET_ABS_DIFF_PER_SAMPLE_CORE_INC);

#undef GET_ABS_DIFF_PER_SAMPLE_CORE_OP
#undef GET_ABS_DIFF_PER_SAMPLE_CORE_INC
}

template <uint8_t maskType>
int64_t getMaskedSampleSumCore(Pel* src, int srcStride, int width, int height, int bitDepth, short* weightMask, int maskStepX, int maskStride, int maskStride2)
{
  const Pel* mask      = weightMask;
  const int  cols      = width;
        int  rows      = height;

  int64_t sum = 0;
  if (maskType == 1) // 1: Use mask
  {
    for (; rows != 0; rows--)
    {
      for (int n = 0; n < cols; n++)
      {
        sum  += (src[n]) * (*mask);
        mask += maskStepX;
      }
      src  += srcStride;
      mask += (maskStride + maskStride2);
    }
  }
  else if (maskType == 2 || maskType == 3) // 2: Use binary mask that contains only 0's and 1's, 3: Inverse the input binary mask before use
  {
    for (; rows != 0; rows--)
    {
      for (int n = 0; n < cols; n++)
      {
        sum += (src[n]) & (maskType == 3 ? ((*mask) - 1) : (-(*mask)));
        mask += maskStepX;
      }
      src  += srcStride;
      mask += (maskStride + maskStride2);
    }
  }
  else // No mask
  {
    for (; rows != 0; rows--)
    {
      for (int n = 0; n < cols; n++)
      {
        sum += src[n];
      }
      src  += srcStride;
    }
  }

  return sum;
}
#endif

#if JVET_W0097_GPM_MMVD_TM
void roundBDCore(const Pel* srcp, const int srcStride, Pel* dest, const int destStride, int width, int height, const ClpRng& clpRng)
{
  const int32_t clipbd = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int32_t shiftDefault = IF_INTERNAL_FRAC_BITS(clipbd);
#else
  const int32_t shiftDefault = std::max<int>(2, (IF_INTERNAL_PREC - clipbd));
#endif
  const int32_t offsetDefault = (1 << (shiftDefault - 1)) + IF_INTERNAL_OFFS;

  if (width == 1)
  {
    THROW("Blocks of width = 1 not supported");
  }
  else
  {
#define RND_OP( ADDR ) dest[ADDR] = ClipPel( rightShift( srcp[ADDR] + offsetDefault, shiftDefault), clpRng )
#define RND_INC        \
    srcp += srcStride;  \
    dest += destStride; \

    SIZE_AWARE_PER_EL_OP(RND_OP, RND_INC);

#undef RND_OP
#undef RND_INC
  }
}

void weightedAvgCore(const Pel* src0, const unsigned src0Stride, const Pel* src1, const unsigned src1Stride, Pel* dest, const unsigned destStride, const int8_t w0, const int8_t w1, int width, int height, const ClpRng& clpRng)
{
  const int8_t log2WeightBase = g_BcwLog2WeightBase;
  const int clipbd = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int shiftNum = IF_INTERNAL_FRAC_BITS(clipbd) + log2WeightBase;
#else
  const int shiftNum = std::max<int>(2, (IF_INTERNAL_PREC - clipbd)) + log2WeightBase;
#endif
  const int offset = (1 << (shiftNum - 1)) + (IF_INTERNAL_OFFS << log2WeightBase);

#define ADD_AVG_OP( ADDR ) dest[ADDR] = ClipPel( rightShift( ( src0[ADDR]*w0 + src1[ADDR]*w1 + offset ), shiftNum ), clpRng )
#define ADD_AVG_INC     \
    src0 += src0Stride; \
    src1 += src1Stride; \
    dest += destStride; \

  SIZE_AWARE_PER_EL_OP(ADD_AVG_OP, ADD_AVG_INC);

#undef ADD_AVG_OP
#undef ADD_AVG_INC
}

void copyClipCore(const Pel* srcp, const unsigned srcStride, Pel* dest, const unsigned destStride, int width, int height, const ClpRng& clpRng)
{
#define RECO_OP( ADDR ) dest[ADDR] = ClipPel( srcp[ADDR], clpRng )
#define RECO_INC        \
  srcp += srcStride;  \
  dest += destStride; \

  SIZE_AWARE_PER_EL_OP(RECO_OP, RECO_INC);

#undef RECO_OP
#undef RECO_INC
}
#endif
template< typename T >
#if JVET_Z0136_OOB
void addAvgCore( const T* src1, int src1Stride, const T* src2, int src2Stride, T* dest, int dstStride, int width, int height, int rshift, int offset, const ClpRng& clpRng, bool *mcMask[2], int mcStride, bool *isOOB)
#else
void addAvgCore( const T* src1, int src1Stride, const T* src2, int src2Stride, T* dest, int dstStride, int width, int height, int rshift, int offset, const ClpRng& clpRng )
#endif
{
#define ADD_AVG_CORE_OP( ADDR ) dest[ADDR] = ClipPel( rightShift( ( src1[ADDR] + src2[ADDR] + offset ), rshift ), clpRng )
#define ADD_AVG_CORE_INC    \
  src1 += src1Stride;       \
  src2 += src2Stride;       \
  dest +=  dstStride;       \

  SIZE_AWARE_PER_EL_OP( ADD_AVG_CORE_OP, ADD_AVG_CORE_INC );

#undef ADD_AVG_CORE_OP
#undef ADD_AVG_CORE_INC
}

void addBIOAvgCore(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, Pel *dst, int dstStride, const Pel *gradX0, const Pel *gradX1, const Pel *gradY0, const Pel*gradY1, int gradStride, int width, int height, int tmpx, int tmpy, int shift, int offset, const ClpRng& clpRng)
{
  int b = 0;

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x += 4)
    {
      b = tmpx * (gradX0[x] - gradX1[x]) + tmpy * (gradY0[x] - gradY1[x]);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      dst[x] = ClipPel(rightShift((src0[x] + src1[x] + b + offset), shift), clpRng);
#else
      dst[x] = ClipPel((int16_t)rightShift((src0[x] + src1[x] + b + offset), shift), clpRng);
#endif

      b = tmpx * (gradX0[x + 1] - gradX1[x + 1]) + tmpy * (gradY0[x + 1] - gradY1[x + 1]);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      dst[x + 1] = ClipPel(rightShift((src0[x + 1] + src1[x + 1] + b + offset), shift), clpRng);
#else
      dst[x + 1] = ClipPel((int16_t)rightShift((src0[x + 1] + src1[x + 1] + b + offset), shift), clpRng);
#endif

      b = tmpx * (gradX0[x + 2] - gradX1[x + 2]) + tmpy * (gradY0[x + 2] - gradY1[x + 2]);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      dst[x + 2] = ClipPel(rightShift((src0[x + 2] + src1[x + 2] + b + offset), shift), clpRng);
#else
      dst[x + 2] = ClipPel((int16_t)rightShift((src0[x + 2] + src1[x + 2] + b + offset), shift), clpRng);
#endif

      b = tmpx * (gradX0[x + 3] - gradX1[x + 3]) + tmpy * (gradY0[x + 3] - gradY1[x + 3]);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      dst[x + 3] = ClipPel(rightShift((src0[x + 3] + src1[x + 3] + b + offset), shift), clpRng);
#else
      dst[x + 3] = ClipPel((int16_t)rightShift((src0[x + 3] + src1[x + 3] + b + offset), shift), clpRng);
#endif
    }
    dst += dstStride;       src0 += src0Stride;     src1 += src1Stride;
    gradX0 += gradStride; gradX1 += gradStride; gradY0 += gradStride; gradY1 += gradStride;
  }
}

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
void calcBIOParameterCore(const Pel* srcY0Tmp, const Pel* srcY1Tmp, Pel* gradX0, Pel* gradX1, Pel* gradY0, Pel* gradY1, int width, int height, const int src0Stride, const int src1Stride, const int widthG, const int bitDepth, Pel* absGX, Pel* absGY, Pel* dIX, Pel* dIY, Pel* signGY_GX, Pel* dI)
{
  width -= 2;
  height -= 2;
  const int bioParamOffset = widthG + 1;
  srcY0Tmp += src0Stride + 1;
  srcY1Tmp += src1Stride + 1;
  gradX0 += bioParamOffset;  gradX1 += bioParamOffset;
  gradY0 += bioParamOffset;  gradY1 += bioParamOffset;
  absGX  += bioParamOffset;  absGY  += bioParamOffset;
  dIX    += bioParamOffset;  dIY    += bioParamOffset;
  signGY_GX += bioParamOffset;
  int shift4 = 4;
  int shift5 = 1;
  if (dI)
  {
    dI += bioParamOffset;
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        int tmpGX = (gradX0[x] + gradX1[x]) >> shift5;
        int tmpGY = (gradY0[x] + gradY1[x]) >> shift5;
        int tmpDI = (int)((srcY1Tmp[x] >> shift4) - (srcY0Tmp[x] >> shift4));
        dI[x] = tmpDI;
        absGX[x] = (tmpGX < 0 ? -tmpGX : tmpGX);
        absGY[x] = (tmpGY < 0 ? -tmpGY : tmpGY);
        dIX[x] = (tmpGX < 0 ? -tmpDI : (tmpGX == 0 ? 0 : tmpDI));
        dIY[x] = (tmpGY < 0 ? -tmpDI : (tmpGY == 0 ? 0 : tmpDI));
        signGY_GX[x] = (tmpGY < 0 ? -tmpGX : (tmpGY == 0 ? 0 : tmpGX));
      }
      srcY0Tmp += src0Stride;
      srcY1Tmp += src1Stride;
      gradX0 += widthG;
      gradX1 += widthG;
      gradY0 += widthG;
      gradY1 += widthG;
      absGX += widthG;
      absGY += widthG;
      dI += widthG;
      dIX += widthG;
      dIY += widthG;
      signGY_GX += widthG;
    }

    return;
  }

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      int tmpGX = (gradX0[x] + gradX1[x]) >> shift5;
      int tmpGY = (gradY0[x] + gradY1[x]) >> shift5;
      int tmpDI = (int)((srcY1Tmp[x] >> shift4) - (srcY0Tmp[x] >> shift4));
      absGX[x] = (tmpGX < 0 ? -tmpGX : tmpGX);
      absGY[x] = (tmpGY < 0 ? -tmpGY : tmpGY);
      dIX[x] = (tmpGX < 0 ? -tmpDI : (tmpGX == 0 ? 0 : tmpDI));
      dIY[x] = (tmpGY < 0 ? -tmpDI : (tmpGY == 0 ? 0 : tmpDI));
      signGY_GX[x] = (tmpGY < 0 ? -tmpGX : (tmpGY == 0 ? 0 : tmpGX));
    }
    srcY0Tmp += src0Stride;
    srcY1Tmp += src1Stride;
    gradX0 += widthG;
    gradX1 += widthG;
    gradY0 += widthG;
    gradY1 += widthG;
    absGX += widthG;
    absGY += widthG;
    dIX += widthG;
    dIY += widthG;
    signGY_GX += widthG;
  }
}

void calcBIOParamSum5Core(Pel* absGX, Pel* absGY, Pel* dIX, Pel* dIY, Pel* signGY_GX, const int widthG, const int width, const int height, int* sumAbsGX, int* sumAbsGY, int* sumDIX, int* sumDIY, int* sumSignGY_GX)
{
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      const int pixel_idx = y * width + x;
      sumAbsGX[pixel_idx] = 0;
      sumAbsGY[pixel_idx] = 0;
      sumDIX[pixel_idx] = 0;
      sumDIY[pixel_idx] = 0;
      sumSignGY_GX[pixel_idx] = 0;
      for (int yy = 0; yy < 5; yy++)
      {
        for (int xx = 0; xx < 5; xx++)
        {
          sumAbsGX[pixel_idx] += absGX[xx];
          sumAbsGY[pixel_idx] += absGY[xx];
          sumDIX[pixel_idx] += dIX[xx];
          sumDIY[pixel_idx] += dIY[xx];
          sumSignGY_GX[pixel_idx] += signGY_GX[xx];
        }
        absGX += widthG;
        absGY += widthG;
        dIX += widthG;
        dIY += widthG;
        signGY_GX += widthG;
      }
      sumDIX[pixel_idx] <<= 2;
      sumDIY[pixel_idx] <<= 2;
      absGX += (1 - 5 * widthG);
      absGY += (1 - 5 * widthG);
      dIX += (1 - 5 * widthG);
      dIY += (1 - 5 * widthG);
      signGY_GX += (1 - 5 * widthG);
    }
    absGX += (widthG - width);
    absGY += (widthG - width);
    dIX += (widthG - width);
    dIY += (widthG - width);
    signGY_GX += (widthG - width);
  }
}

void calcBIOParamSum4Core(Pel* absGX, Pel* absGY, Pel* dIX, Pel* dIY, Pel* signGY_GX, int width, int height, const int widthG, int* sumAbsGX, int* sumAbsGY, int* sumDIX, int* sumDIY, int* sumSignGY_GX)
{
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      *sumAbsGX += absGX[x];
      *sumAbsGY += absGY[x];
      *sumDIX += dIX[x];
      *sumDIY += dIY[x];
      *sumSignGY_GX += signGY_GX[x];
    }
    absGX += widthG;
    absGY += widthG;
    dIX += widthG;
    dIY += widthG;
    signGY_GX += widthG;
  }
}

void calcBIOClippedVxVyCore(int* sumDIX_pixel_32bit, int* sumAbsGX_pixel_32bit, int* sumDIY_pixel_32bit, int* sumAbsGY_pixel_32bit, int* sumSignGY_GX_pixel_32bit, const int limit, const int bioSubblockSize, int* tmpx_pixel_32bit, int* tmpy_pixel_32bit)
{
  for (int idx = 0; idx < bioSubblockSize; idx++)
  {
    *tmpx_pixel_32bit = Clip3(-limit, limit, (*sumDIX_pixel_32bit) >> (*sumAbsGX_pixel_32bit));
    int tmpData = ((*sumSignGY_GX_pixel_32bit) * (*tmpx_pixel_32bit)) >> 1;
    *tmpy_pixel_32bit = Clip3(-limit, limit, (((*sumDIY_pixel_32bit) - tmpData) >> (*sumAbsGY_pixel_32bit)));
    sumDIX_pixel_32bit++;
    sumAbsGX_pixel_32bit++;
    sumDIY_pixel_32bit++;
    sumAbsGY_pixel_32bit++;
    sumSignGY_GX_pixel_32bit++;
    tmpx_pixel_32bit++;
    tmpy_pixel_32bit++;
  }
}
#if JVET_Z0136_OOB
void addBIOAvgNCore(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, Pel *dst, int dstStride, const Pel *gradX0, const Pel *gradX1, const Pel *gradY0, const Pel *gradY1, int gradStride, int width, int height, int* tmpx, int* tmpy, int shift, int offset, const ClpRng& clpRng, bool *mcMask[2], int mcStride, bool *isOOB)
#else
void addBIOAvgNCore(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, Pel *dst, int dstStride, const Pel *gradX0, const Pel *gradX1, const Pel *gradY0, const Pel *gradY1, int gradStride, int width, int height, int* tmpx, int* tmpy, int shift, int offset, const ClpRng& clpRng)
#endif
{
  int b = 0;
#if JVET_Z0136_OOB
  int offset2 = offset >> 1;
  int shift2 = shift - 1;
  bool *pMcMask0 = mcMask[0];
  bool *pMcMask1 = mcMask[1];
  if (isOOB[0] || isOOB[1])
  {
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        b = (int)tmpx[x] * (gradX0[x] - gradX1[x]) + (int)tmpy[x] * (gradY0[x] - gradY1[x]);
        bool oob0 = pMcMask0[x];
        bool oob1 = pMcMask1[x];
        if (oob0 && !oob1)
        {
          dst[x] = ClipPel(rightShift(src1[x] + offset2, shift2), clpRng);
        }
        else if (!oob0 && oob1)
        {
          dst[x] = ClipPel(rightShift(src0[x] + offset2, shift2), clpRng);
        }
        else
        {
          dst[x] = ClipPel(rightShift((src0[x] + src1[x] + b + offset), shift), clpRng);
        }
      }
      pMcMask0 += mcStride;
      pMcMask1 += mcStride;
      tmpx += width;
      tmpy += width;
      dst += dstStride;
      src0 += src0Stride;
      src1 += src1Stride;
      gradX0 += gradStride;
      gradX1 += gradStride;
      gradY0 += gradStride;
      gradY1 += gradStride;
    }
  }
  else
  {
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        b = (int)tmpx[x] * (gradX0[x] - gradX1[x]) + (int)tmpy[x] * (gradY0[x] - gradY1[x]);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
        dst[x] = ClipPel(rightShift((src0[x] + src1[x] + b + offset), shift), clpRng);
#else
        dst[x] = ClipPel((int16_t)rightShift((src0[x] + src1[x] + b + offset), shift), clpRng);
#endif
      }
      tmpx += width;
      tmpy += width;
      dst += dstStride;
      src0 += src0Stride;
      src1 += src1Stride;
      gradX0 += gradStride;
      gradX1 += gradStride;
      gradY0 += gradStride;
      gradY1 += gradStride;
    }
  }
#else
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      b = (int)tmpx[x] * (gradX0[x] - gradX1[x]) + (int)tmpy[x] * (gradY0[x] - gradY1[x]);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      dst[x] = ClipPel(rightShift((src0[x] + src1[x] + b + offset), shift), clpRng);
#else
      dst[x] = ClipPel((int16_t)rightShift((src0[x] + src1[x] + b + offset), shift), clpRng);
#endif
    }
    tmpx += width;    tmpy += width;
    dst += dstStride;       src0 += src0Stride;     src1 += src1Stride;
    gradX0 += gradStride; gradX1 += gradStride; gradY0 += gradStride; gradY1 += gradStride;
  }
#endif
  return;
}

void calAbsSumCore(const Pel* diff, int stride, int width, int height, int* absSum)
{
  *absSum = 0;
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      *absSum += ::abs(diff[x]);
    }
    diff += stride;
  }
}
#endif

template<bool PAD = true>
void gradFilterCore(Pel* pSrc, int srcStride, int width, int height, int gradStride, Pel* gradX, Pel* gradY, const int bitDepth)
{
  Pel* srcTmp = pSrc + srcStride + 1;
  Pel* gradXTmp = gradX + gradStride + 1;
  Pel* gradYTmp = gradY + gradStride + 1;
  int  shift1 = 6;

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  for (int y = 0; y < (height - 2); y++)
  {
    for (int x = 0; x < (width - 2); x++)
#else
  for (int y = 0; y < (height - 2 * BIO_EXTEND_SIZE); y++)
  {
    for (int x = 0; x < (width - 2 * BIO_EXTEND_SIZE); x++)
#endif
    {
      gradYTmp[x] = ( srcTmp[x + srcStride] >> shift1 ) - ( srcTmp[x - srcStride] >> shift1 );
      gradXTmp[x] = ( srcTmp[x + 1] >> shift1 ) - ( srcTmp[x - 1] >> shift1 );
    }
    gradXTmp += gradStride;
    gradYTmp += gradStride;
    srcTmp += srcStride;
  }

#if !MULTI_PASS_DMVR && !SAMPLE_BASED_BDOF
  if (PAD)
  {
  gradXTmp = gradX + gradStride + 1;
  gradYTmp = gradY + gradStride + 1;
  for (int y = 0; y < (height - 2 * BIO_EXTEND_SIZE); y++)
  {
    gradXTmp[-1] = gradXTmp[0];
    gradXTmp[width - 2 * BIO_EXTEND_SIZE] = gradXTmp[width - 2 * BIO_EXTEND_SIZE - 1];
    gradXTmp += gradStride;

    gradYTmp[-1] = gradYTmp[0];
    gradYTmp[width - 2 * BIO_EXTEND_SIZE] = gradYTmp[width - 2 * BIO_EXTEND_SIZE - 1];
    gradYTmp += gradStride;
  }

  gradXTmp = gradX + gradStride;
  gradYTmp = gradY + gradStride;
  ::memcpy(gradXTmp - gradStride, gradXTmp, sizeof(Pel)*(width));
  ::memcpy(gradXTmp + (height - 2 * BIO_EXTEND_SIZE)*gradStride, gradXTmp + (height - 2 * BIO_EXTEND_SIZE - 1)*gradStride, sizeof(Pel)*(width));
  ::memcpy(gradYTmp - gradStride, gradYTmp, sizeof(Pel)*(width));
  ::memcpy(gradYTmp + (height - 2 * BIO_EXTEND_SIZE)*gradStride, gradYTmp + (height - 2 * BIO_EXTEND_SIZE - 1)*gradStride, sizeof(Pel)*(width));
  }
#endif
}

void calcBIOSumsCore(const Pel* srcY0Tmp, const Pel* srcY1Tmp, Pel* gradX0, Pel* gradX1, Pel* gradY0, Pel* gradY1, int xu, int yu, const int src0Stride, const int src1Stride, const int widthG, const int bitDepth, int* sumAbsGX, int* sumAbsGY, int* sumDIX, int* sumDIY, int* sumSignGY_GX)
{
  int shift4 = 4;
  int shift5 = 1;

  for (int y = 0; y < 6; y++)
  {
    for (int x = 0; x < 6; x++)
    {
      int tmpGX = (gradX0[x] + gradX1[x]) >> shift5;
      int tmpGY = (gradY0[x] + gradY1[x]) >> shift5;
      int tmpDI = (int)((srcY1Tmp[x] >> shift4) - (srcY0Tmp[x] >> shift4));
      *sumAbsGX += (tmpGX < 0 ? -tmpGX : tmpGX);
      *sumAbsGY += (tmpGY < 0 ? -tmpGY : tmpGY);
      *sumDIX += (tmpGX < 0 ? -tmpDI : (tmpGX == 0 ? 0 : tmpDI));
      *sumDIY += (tmpGY < 0 ? -tmpDI : (tmpGY == 0 ? 0 : tmpDI));
      *sumSignGY_GX += (tmpGY < 0 ? -tmpGX : (tmpGY == 0 ? 0 : tmpGX));

    }
    srcY1Tmp += src1Stride;
    srcY0Tmp += src0Stride;
    gradX0 += widthG;
    gradX1 += widthG;
    gradY0 += widthG;
    gradY1 += widthG;
  }
}


void calcBlkGradientCore(int sx, int sy, int     *arraysGx2, int     *arraysGxGy, int     *arraysGxdI, int     *arraysGy2, int     *arraysGydI, int     &sGx2, int     &sGy2, int     &sGxGy, int     &sGxdI, int     &sGydI, int width, int height, int unitSize)
{
  int     *Gx2 = arraysGx2;
  int     *Gy2 = arraysGy2;
  int     *GxGy = arraysGxGy;
  int     *GxdI = arraysGxdI;
  int     *GydI = arraysGydI;

  // set to the above row due to JVET_K0485_BIO_EXTEND_SIZE
  Gx2 -= (BIO_EXTEND_SIZE*width);
  Gy2 -= (BIO_EXTEND_SIZE*width);
  GxGy -= (BIO_EXTEND_SIZE*width);
  GxdI -= (BIO_EXTEND_SIZE*width);
  GydI -= (BIO_EXTEND_SIZE*width);

  for (int y = -BIO_EXTEND_SIZE; y < unitSize + BIO_EXTEND_SIZE; y++)
  {
    for (int x = -BIO_EXTEND_SIZE; x < unitSize + BIO_EXTEND_SIZE; x++)
    {
      sGx2 += Gx2[x];
      sGy2 += Gy2[x];
      sGxGy += GxGy[x];
      sGxdI += GxdI[x];
      sGydI += GydI[x];
    }
    Gx2 += width;
    Gy2 += width;
    GxGy += width;
    GxdI += width;
    GydI += width;
  }
}

#if ENABLE_SIMD_OPT_BCW && defined(TARGET_SIMD_X86)
void removeWeightHighFreq(int16_t* dst, int dstStride, const int16_t* src, int srcStride, int width, int height, int shift, int bcwWeight)
{
  int normalizer = ((1 << 16) + (bcwWeight > 0 ? (bcwWeight >> 1) : -(bcwWeight >> 1))) / bcwWeight;
  int weight0 = normalizer << g_BcwLog2WeightBase;
  int weight1 = (g_BcwWeightBase - bcwWeight)*normalizer;
#define REM_HF_INC  \
  src += srcStride; \
  dst += dstStride; \

#define REM_HF_OP( ADDR )      dst[ADDR] =             (dst[ADDR]*weight0 - src[ADDR]*weight1 + (1<<15))>>16

  SIZE_AWARE_PER_EL_OP(REM_HF_OP, REM_HF_INC);

#undef REM_HF_INC
#undef REM_HF_OP
#undef REM_HF_OP_CLIP
}

void removeHighFreq(int16_t* dst, int dstStride, const int16_t* src, int srcStride, int width, int height)
{
#define REM_HF_INC  \
  src += srcStride; \
  dst += dstStride; \

#define REM_HF_OP( ADDR )      dst[ADDR] =             2 * dst[ADDR] - src[ADDR]

  SIZE_AWARE_PER_EL_OP(REM_HF_OP, REM_HF_INC);

#undef REM_HF_INC
#undef REM_HF_OP
#undef REM_HF_OP_CLIP
}
#endif

template<typename T>
void reconstructCore( const T* src1, int src1Stride, const T* src2, int src2Stride, T* dest, int dstStride, int width, int height, const ClpRng& clpRng )
{
#define RECO_CORE_OP( ADDR ) dest[ADDR] = ClipPel( src1[ADDR] + src2[ADDR], clpRng )
#define RECO_CORE_INC     \
  src1 += src1Stride;     \
  src2 += src2Stride;     \
  dest +=  dstStride;     \

  SIZE_AWARE_PER_EL_OP( RECO_CORE_OP, RECO_CORE_INC );

#undef RECO_CORE_OP
#undef RECO_CORE_INC
}


template<typename T>
void linTfCore( const T* src, int srcStride, Pel *dst, int dstStride, int width, int height, int scale, int shift, int offset, const ClpRng& clpRng, bool bClip )
{
#define LINTF_CORE_OP( ADDR ) dst[ADDR] = ( Pel ) bClip ? ClipPel( rightShift( scale * src[ADDR], shift ) + offset, clpRng ) : ( rightShift( scale * src[ADDR], shift ) + offset )
#define LINTF_CORE_INC  \
  src += srcStride;     \
  dst += dstStride;     \

  SIZE_AWARE_PER_EL_OP( LINTF_CORE_OP, LINTF_CORE_INC );

#undef LINTF_CORE_OP
#undef LINTF_CORE_INC
}

#if JVET_Z0136_OOB
bool isMvOOBCore(const Mv& rcMv, const struct Position pos, const struct Size size, const SPS* sps, const PPS* pps, bool *mcMask, bool *mcMaskChroma, bool lumaOnly, ChromaFormat componentID)
{
  int chromaScale = getComponentScaleX(COMPONENT_Cb, componentID);
  const int mvstep = 1 << MV_FRACTIONAL_BITS_INTERNAL;
  const int mvstepHalf = mvstep >> 1;

  int horMax = (((int)pps->getPicWidthInLumaSamples() - 1) << MV_FRACTIONAL_BITS_INTERNAL) + mvstepHalf;
  int horMin = -mvstepHalf;
  int verMax = (((int)pps->getPicHeightInLumaSamples() - 1) << MV_FRACTIONAL_BITS_INTERNAL) + mvstepHalf;
  int verMin = -mvstepHalf;

  int offsetX = (pos.x << MV_FRACTIONAL_BITS_INTERNAL) + rcMv.getHor();
  int offsetY = (pos.y << MV_FRACTIONAL_BITS_INTERNAL) + rcMv.getVer();
  bool isOOB = false;
  if ((offsetX <= horMin)
    || ((offsetX + (size.width << MV_FRACTIONAL_BITS_INTERNAL) - 1) >= horMax)
    || (offsetY <= verMin)
    || ((offsetY + (size.height << MV_FRACTIONAL_BITS_INTERNAL) - 1) >= verMax))
  {
    isOOB = true;
  }
  if (isOOB)
  {
    int baseOffsetX = offsetX;
    bool *pMcMask = mcMask;

    for (int y = 0; y < size.height; y++, offsetY += mvstep)
    {
      offsetX = baseOffsetX;
      bool checkY = (offsetY <= verMin) || (offsetY >= verMax);
      for (int x = 0; x < size.width; x++, offsetX += mvstep)
      {
        pMcMask[x] = (offsetX <= horMin) || (offsetX >= horMax) || checkY;
      }
      pMcMask += size.width;
    }

    if (!lumaOnly)
    {
      bool *pMcMaskChroma = mcMaskChroma;
      pMcMask = mcMask;
      int widthChroma = (size.width) >> chromaScale;
      int heightChroma = (size.height) >> chromaScale;
      int widthLuma2 = size.width << chromaScale;
      for (int y = 0; y < heightChroma; y++)
      {
        for (int x = 0; x < widthChroma; x++)
        {
          pMcMaskChroma[x] = pMcMask[x << chromaScale];
        }
        pMcMaskChroma += widthChroma;
        pMcMask += widthLuma2;
      }
    }
  }
  else
  {
    bool *pMcMask = mcMask;
    memset(pMcMask, false, size.width * size.height);

    bool *pMcMaskChroma = mcMaskChroma;
    int widthChroma = (size.width) >> chromaScale;
    int heightChroma = (size.height) >> chromaScale;
    memset(pMcMaskChroma, false, widthChroma * heightChroma);
  }
  return isOOB;
}

bool isMvOOBSubBlkCore(const Mv& rcMv, const struct Position pos, const struct Size size, const SPS* sps, const PPS* pps, bool *mcMask, int mcStride, bool *mcMaskChroma, int mcCStride, bool lumaOnly, ChromaFormat componentID)
{
  int chromaScale = getComponentScaleX(COMPONENT_Cb, componentID);
  const int mvstep = 1 << MV_FRACTIONAL_BITS_INTERNAL;
  const int mvstepHalf = mvstep >> 1;

  int horMax = (((int)pps->getPicWidthInLumaSamples() - 1) << MV_FRACTIONAL_BITS_INTERNAL) + mvstepHalf;
  int horMin = -mvstepHalf;
  int verMax = (((int)pps->getPicHeightInLumaSamples() - 1) << MV_FRACTIONAL_BITS_INTERNAL) + mvstepHalf;
  int verMin = -mvstepHalf;

  int offsetX = (pos.x << MV_FRACTIONAL_BITS_INTERNAL) + rcMv.getHor();
  int offsetY = (pos.y << MV_FRACTIONAL_BITS_INTERNAL) + rcMv.getVer();
  bool isOOB = false;
  if ((offsetX <= horMin)
    || ((offsetX + (size.width << MV_FRACTIONAL_BITS_INTERNAL) - 1) >= horMax)
    || (offsetY <= verMin)
    || ((offsetY + (size.height << MV_FRACTIONAL_BITS_INTERNAL) - 1) >= verMax))
  {
    isOOB = true;
  }
  if (isOOB)
  {
    int baseOffsetX = offsetX;
    bool *pMcMask = mcMask;
    for (int y = 0; y < size.height; y++, offsetY += mvstep)
    {
      offsetX = baseOffsetX;
      bool checkY = (offsetY <= verMin) || (offsetY >= verMax);;
      for (int x = 0; x < size.width; x++, offsetX += mvstep)
      {
        pMcMask[x] = (offsetX <= horMin) || (offsetX >= horMax) || checkY;
      }
      pMcMask += mcStride;
    }

    if (!lumaOnly)
    {
      bool *pMcMaskChroma = mcMaskChroma;
      pMcMask = mcMask;
      int widthChroma = (size.width) >> chromaScale;
      int heightChroma = (size.height) >> chromaScale;
      int strideLuma2 = mcStride << chromaScale;
      for (int y = 0; y < heightChroma; y++)
      {
        for (int x = 0; x < widthChroma; x++)
        {
          pMcMaskChroma[x] = pMcMask[x << chromaScale];
        }
        pMcMaskChroma += mcCStride;
        pMcMask += strideLuma2;
      }
    }
  }
  else
  {
    bool *pMcMask = mcMask;
    for (int y = 0; y < size.height; y++)
    {
      memset(pMcMask, false, size.width);
      pMcMask += mcStride;
    }

    bool *pMcMaskChroma = mcMaskChroma;
    int widthChroma = (size.width) >> chromaScale;
    int heightChroma = (size.height) >> chromaScale;
    for (int y = 0; y < heightChroma; y++)
    {
      memset(pMcMaskChroma, false, widthChroma);
      pMcMaskChroma += mcCStride;
    }
  }
  return isOOB;
}
#endif

PelBufferOps::PelBufferOps()
{
#if JVET_W0097_GPM_MMVD_TM
  roundBD = roundBDCore;
  weightedAvg = weightedAvgCore;
  copyClip = copyClipCore;
#endif
  addAvg4 = addAvgCore<Pel>;
  addAvg8 = addAvgCore<Pel>;

  reco4 = reconstructCore<Pel>;
  reco8 = reconstructCore<Pel>;

  linTf4 = linTfCore<Pel>;
  linTf8 = linTfCore<Pel>;

  addBIOAvg4      = addBIOAvgCore;
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  calcBIOParameter   = calcBIOParameterCore;
  calcBIOParamSum5   = calcBIOParamSum5Core;
  calcBIOParamSum4   = calcBIOParamSum4Core;
  calcBIOClippedVxVy = calcBIOClippedVxVyCore;
  addBIOAvgN         = addBIOAvgNCore;
  calAbsSum          = calAbsSumCore;
  bioGradFilter      = gradFilterCore <false>;
#else
  bioGradFilter   = gradFilterCore;
#endif
  calcBIOSums = calcBIOSumsCore;

  copyBuffer = copyBufferCore;
  padding = paddingCore;
#if ENABLE_SIMD_OPT_BCW && defined(TARGET_SIMD_X86)
  removeWeightHighFreq8 = removeWeightHighFreq;
  removeWeightHighFreq4 = removeWeightHighFreq;
  removeHighFreq8 = removeHighFreq;
  removeHighFreq4 = removeHighFreq;
#endif

  profGradFilter = gradFilterCore <false>;
  applyPROF      = applyPROFCore;
  roundIntVector = nullptr;
#if TM_AMVP || TM_MRG
  getSumOfDifference = getSumOfDifferenceCore;
#endif
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  getAbsoluteDifferencePerSample = getAbsoluteDifferencePerSampleCore;
  getSampleSumFunc[0] = getMaskedSampleSumCore<0>;
  getSampleSumFunc[1] = getMaskedSampleSumCore<1>;
  getSampleSumFunc[2] = getMaskedSampleSumCore<2>;
  getSampleSumFunc[3] = getMaskedSampleSumCore<3>;
#endif
#if JVET_Z0136_OOB
  isMvOOB = isMvOOBCore;
  isMvOOBSubBlk = isMvOOBSubBlkCore;
#endif
}

PelBufferOps g_pelBufOP = PelBufferOps();

void copyBufferCore(Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height)
{
  int numBytes = width * sizeof(Pel);
  for (int i = 0; i < height; i++)
  {
    memcpy(dst + i * dstStride, src + i * srcStride, numBytes);
  }
}

void paddingCore(Pel *ptr, int stride, int width, int height, int padSize)
{
  /*left and right padding*/
  Pel *ptrTemp1 = ptr;
  Pel *ptrTemp2 = ptr + (width - 1);
  int offset = 0;
  for (int i = 0; i < height; i++)
  {
    offset = stride * i;
    for (int j = 1; j <= padSize; j++)
    {
      *(ptrTemp1 - j + offset) = *(ptrTemp1 + offset);
      *(ptrTemp2 + j + offset) = *(ptrTemp2 + offset);
    }
  }
  /*Top and Bottom padding*/
  int numBytes = (width + padSize + padSize) * sizeof(Pel);
  ptrTemp1 = (ptr - padSize);
  ptrTemp2 = (ptr + (stride * (height - 1)) - padSize);
  for (int i = 1; i <= padSize; i++)
  {
    memcpy(ptrTemp1 - (i * stride), (ptrTemp1), numBytes);
    memcpy(ptrTemp2 + (i * stride), (ptrTemp2), numBytes);
  }
}

#if MULTI_HYP_PRED
template<>
void AreaBuf<Pel>::addHypothesisAndClip(const AreaBuf<const Pel> &other, const int weight, const ClpRng& clpRng)
{
  CHECK(width != other.width, "Incompatible size");
  CHECK(height != other.height, "Incompatible size");

  Pel* dest = buf;
  const Pel* src = other.buf;
  const int counterweight = (1 << MULTI_HYP_PRED_WEIGHT_BITS) - weight;
  const int add = 1 << (MULTI_HYP_PRED_WEIGHT_BITS - 1);

#define ADD_HYP_OP( ADDR ) dest[ADDR] = ClipPel( ( counterweight*dest[ADDR] + weight*src[ADDR] + add ) >> MULTI_HYP_PRED_WEIGHT_BITS, clpRng )
#define ADD_HYP_INC     \
    dest += stride; \
    src += other.stride;

  SIZE_AWARE_PER_EL_OP(ADD_HYP_OP, ADD_HYP_INC);

#undef ADD_HYP_OP
#undef ADD_HYP_INC
}
#endif

template<>
#if JVET_Z0136_OOB
void AreaBuf<Pel>::addWeightedAvg(const AreaBuf<const Pel> &other1, const AreaBuf<const Pel> &other2, const ClpRng& clpRng, const int8_t bcwIdx, bool *mcMask[2], int mcStride, bool* isOOB)
#else
void AreaBuf<Pel>::addWeightedAvg(const AreaBuf<const Pel> &other1, const AreaBuf<const Pel> &other2, const ClpRng& clpRng, const int8_t bcwIdx)
#endif
{
#if JVET_W0097_GPM_MMVD_TM
#if JVET_Z0136_OOB
  int8_t w0 = getBcwWeight(bcwIdx, REF_PIC_LIST_0);
  int8_t w1 = getBcwWeight(bcwIdx, REF_PIC_LIST_1);

  const int8_t log2WeightBase = g_BcwLog2WeightBase;
  const Pel* src1 = other1.buf;
  const Pel* src2 = other2.buf;
  Pel* dest = buf;

  const unsigned src1Stride = other1.stride;
  const unsigned src2Stride = other2.stride;
  const unsigned destStride = stride;
  const int clipbd = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int shiftNum = IF_INTERNAL_FRAC_BITS(clipbd) + log2WeightBase;
#else
  const int shiftNum = std::max<int>(2, (IF_INTERNAL_PREC - clipbd)) + log2WeightBase;
#endif
  const int offset = (1 << (shiftNum - 1)) + (IF_INTERNAL_OFFS << log2WeightBase);
  if (!isOOB[0] && !isOOB[1])
  {
    g_pelBufOP.weightedAvg(src1, src1Stride, src2, src2Stride, dest, destStride, w0, w1, width, height, clpRng);
  }
  else
  {
    int shiftNum2 = IF_INTERNAL_FRAC_BITS(clipbd);
    const int offset2 = (1 << (shiftNum2 - 1)) + IF_INTERNAL_OFFS;
    bool *pMcMask0 = mcMask[0];
    bool *pMcMask1 = mcMask[1];

    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        bool oob0 = pMcMask0[x];
        bool oob1 = pMcMask1[x];
        if (oob0 && !oob1)
        {
          dest[x] = ClipPel(rightShift(src2[x] + offset2, shiftNum2), clpRng);
        }
        else if (!oob0 && oob1)
        {
          dest[x] = ClipPel(rightShift(src1[x] + offset2, shiftNum2), clpRng);
        }
        else
        {
          dest[x + 0] = ClipPel(rightShift((src1[x] * w0 + src2[x + 0] * w1 + offset), shiftNum), clpRng);
        }
      }
      pMcMask0 += mcStride;
      pMcMask1 += mcStride;
      src1 += src1Stride;
      src2 += src2Stride;
      dest += destStride;
    }
  }
#else
  const int8_t w0 = getBcwWeight(bcwIdx, REF_PIC_LIST_0);
  const int8_t w1 = getBcwWeight(bcwIdx, REF_PIC_LIST_1);

  const Pel*            src0 = other1.buf;
  const Pel*            src1 = other2.buf;
  Pel*                 dest = buf;
  const unsigned src0Stride = other1.stride;
  const unsigned src1Stride = other2.stride;
  const unsigned destStride = stride;

  g_pelBufOP.weightedAvg(src0, src0Stride, src1, src1Stride, dest, destStride, w0, w1, width, height, clpRng);
#endif
#else
  const int8_t w0 = getBcwWeight(bcwIdx, REF_PIC_LIST_0);
  const int8_t w1 = getBcwWeight(bcwIdx, REF_PIC_LIST_1);
  const int8_t log2WeightBase = g_BcwLog2WeightBase;

  const Pel* src0 = other1.buf;
  const Pel* src2 = other2.buf;
  Pel* dest = buf;

  const unsigned src1Stride = other1.stride;
  const unsigned src2Stride = other2.stride;
  const unsigned destStride = stride;
  const int clipbd = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int shiftNum = IF_INTERNAL_FRAC_BITS(clipbd) + log2WeightBase;
#else
  const int shiftNum = std::max<int>(2, (IF_INTERNAL_PREC - clipbd)) + log2WeightBase;
#endif
  const int offset = (1 << (shiftNum - 1)) + (IF_INTERNAL_OFFS << log2WeightBase);

#define ADD_AVG_OP( ADDR ) dest[ADDR] = ClipPel( rightShift( ( src0[ADDR]*w0 + src2[ADDR]*w1 + offset ), shiftNum ), clpRng )
#define ADD_AVG_INC     \
    src0 += src1Stride; \
    src2 += src2Stride; \
    dest += destStride; \

  SIZE_AWARE_PER_EL_OP(ADD_AVG_OP, ADD_AVG_INC);

#undef ADD_AVG_OP
#undef ADD_AVG_INC
#endif
}

template<>
void AreaBuf<Pel>::rspSignal(std::vector<Pel>& pLUT)
{
  Pel* dst = buf;
  Pel* src = buf;
    for (unsigned y = 0; y < height; y++)
    {
      for (unsigned x = 0; x < width; x++)
      {
        dst[x] = pLUT[src[x]];
      }
      dst += stride;
      src += stride;
    }
}

template<>
void AreaBuf<Pel>::rspSignal(const AreaBuf<const Pel>& other, std::vector<Pel>& pLUT)
{
  CHECK( width != other.width, "Incompatible size" );
  CHECK( height != other.height, "Incompatible size" );

  Pel* dst = buf;
  const Pel* src = other.buf;
  for (unsigned y = 0; y < height; y++)
  {
    for (unsigned x = 0; x < width; x++)
    {
      dst[x] = pLUT[src[x]];
    }
    dst += stride;
    src += other.stride;
  }
}

template<>
void AreaBuf<Pel>::rspSignal( const AreaBuf<Pel> &toReshape, std::vector<Pel>& pLUT )
{
  CHECK( width != toReshape.width, "Incompatible size" );
  CHECK( height != toReshape.height, "Incompatible size" );

  Pel* dst = buf;
  Pel* src = toReshape.buf;
  const int srcStride = toReshape.stride;

  for( unsigned y = 0; y < height; y++ )
  {
    for( unsigned x = 0; x < width; x++ )
    {
      dst[x] = pLUT[src[x]];
    }
    dst += stride;
    src += srcStride;
  }
}

template<>
void AreaBuf<Pel>::rspSignalAllAndSubtract( const AreaBuf<Pel> &buffer1, const AreaBuf<Pel> &buffer2, std::vector<Pel>& pLUT )
{
  CHECK( width != buffer1.width, "Incompatible size in buffer1" );
  CHECK( height != buffer1.height, "Incompatible size in buffer1" );
  CHECK( width != buffer2.width, "Incompatible size in buffer2" );
  CHECK( height != buffer2.height, "Incompatible size in buffer2" );

  Pel* dest = buf;
  const Pel* buf1 = buffer1.buf;
  const Pel* buf2 = buffer2.buf;

#define SUBS_INC           \
  dest +=          stride; \
  buf1 +=  buffer1.stride; \
  buf2 +=  buffer2.stride; \

#define SUBS_OP( ADDR ) dest[ADDR] = pLUT[buf1[ADDR]] - pLUT[buf2[ADDR]]

  SIZE_AWARE_PER_EL_OP( SUBS_OP, SUBS_INC );

#undef SUBS_OP
#undef SUBS_INC
}

template<>
void AreaBuf<Pel>::rspSignalAndSubtract( const AreaBuf<Pel> &buffer1, const AreaBuf<Pel> &buffer2, std::vector<Pel>& pLUT )
{
  CHECK( width != buffer1.width, "Incompatible size in buffer1" );
  CHECK( height != buffer1.height, "Incompatible size in buffer1" );
  CHECK( width != buffer2.width, "Incompatible size in buffer2" );
  CHECK( height != buffer2.height, "Incompatible size in buffer2" );

  Pel* dest = buf;
  const Pel* buf1 = buffer1.buf;
  const Pel* buf2 = buffer2.buf;

#define SUBS_INC           \
  dest +=          stride; \
  buf1 +=  buffer1.stride; \
  buf2 +=  buffer2.stride; \

#define SUBS_OP( ADDR ) dest[ADDR] = pLUT[buf1[ADDR]] - buf2[ADDR]

  SIZE_AWARE_PER_EL_OP( SUBS_OP, SUBS_INC );

#undef SUBS_OP
#undef SUBS_INC
}

template<>
void AreaBuf<Pel>::scaleSignal(const int scale, const bool dir, const ClpRng& clpRng)
{
  Pel* dst = buf;
  Pel* src = buf;
  int sign, absval;
  int maxAbsclipBD = (1<<clpRng.bd) - 1;

  if (dir) // forward
  {
    if (width == 1)
    {
      THROW("Blocks of width = 1 not supported");
    }
    else
    {
      for (unsigned y = 0; y < height; y++)
      {
        for (unsigned x = 0; x < width; x++)
        {
          sign = src[x] >= 0 ? 1 : -1;
          absval = sign * src[x];
          dst[x] = (Pel)Clip3(-maxAbsclipBD, maxAbsclipBD, sign * (((absval << CSCALE_FP_PREC) + (scale >> 1)) / scale));
        }
        dst += stride;
        src += stride;
      }
    }
  }
  else // inverse
  {
    for (unsigned y = 0; y < height; y++)
    {
      for (unsigned x = 0; x < width; x++)
      {
        src[x] = (Pel)Clip3((Pel)(-maxAbsclipBD - 1), (Pel)maxAbsclipBD, src[x]);
        sign = src[x] >= 0 ? 1 : -1;
        absval = sign * src[x];
        int val = sign * ((absval * scale + (1 << (CSCALE_FP_PREC - 1))) >> CSCALE_FP_PREC);
        if (sizeof(Pel) == 2) // avoid overflow when storing data
        {
           val = Clip3<int>(-32768, 32767, val);
        }
        dst[x] = (Pel)val;
      }
      dst += stride;
      src += stride;
    }
  }
}

template<>
#if JVET_Z0136_OOB
void AreaBuf<Pel>::addAvg(const AreaBuf<const Pel> &other1, const AreaBuf<const Pel> &other2, const ClpRng& clpRng, bool *mcMask[2], int mcStride, bool* isOOB)
#else
void AreaBuf<Pel>::addAvg( const AreaBuf<const Pel> &other1, const AreaBuf<const Pel> &other2, const ClpRng& clpRng)
#endif
{
  const Pel* src0 = other1.buf;
  const Pel* src2 = other2.buf;
        Pel* dest =        buf;

  const unsigned src1Stride = other1.stride;
  const unsigned src2Stride = other2.stride;
  const unsigned destStride =        stride;
  const int     clipbd      = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int shiftNum = IF_INTERNAL_FRAC_BITS(clipbd) + 1;
#else
  const int     shiftNum    = std::max<int>(2, (IF_INTERNAL_PREC - clipbd)) + 1;
#endif
  const int     offset      = (1 << (shiftNum - 1)) + 2 * IF_INTERNAL_OFFS;

#if JVET_Z0136_OOB
  if (mcMask == NULL || (!isOOB[0] && !isOOB[1]))
  {
#if ENABLE_SIMD_OPT_BUFFER && defined(TARGET_SIMD_X86)
    if ((width & 7) == 0)
    {
      g_pelBufOP.addAvg8(src0, src1Stride, src2, src2Stride, dest, destStride, width, height, shiftNum, offset, clpRng, mcMask, mcStride, isOOB);
    }
    else if ((width & 3) == 0)
    {
      g_pelBufOP.addAvg4(src0, src1Stride, src2, src2Stride, dest, destStride, width, height, shiftNum, offset, clpRng, mcMask, mcStride, isOOB);
    }
    else
#endif
    {
#define ADD_AVG_OP( ADDR ) dest[ADDR] = ClipPel( rightShift( ( src0[ADDR] + src2[ADDR] + offset ), shiftNum ), clpRng )
#define ADD_AVG_INC     \
    src0 += src1Stride; \
    src2 += src2Stride; \
    dest += destStride; \

    SIZE_AWARE_PER_EL_OP(ADD_AVG_OP, ADD_AVG_INC);

#undef ADD_AVG_OP
#undef ADD_AVG_INC
    }
  }
  else
  {
    int shiftNum2 = IF_INTERNAL_FRAC_BITS(clipbd);
    const int offset2 = (1 << (shiftNum2 - 1)) + IF_INTERNAL_OFFS;
    bool *pMcMask0 = mcMask[0];
    bool *pMcMask1 = mcMask[1];
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        bool oob0 = pMcMask0[x];
        bool oob1 = pMcMask1[x];
        if (oob0 && !oob1)
        {
          dest[x] = ClipPel(rightShift(src2[x] + offset2, shiftNum2), clpRng);
        }
        else if (!oob0 && oob1)
        {
          dest[x] = ClipPel(rightShift(src0[x] + offset2, shiftNum2), clpRng);
        }
        else
        {
          dest[x] = ClipPel(rightShift((src0[x] + src2[x] + offset), shiftNum), clpRng);
        }
      }
      pMcMask0 += mcStride;
      pMcMask1 += mcStride;
      src0 += src1Stride;
      src2 += src2Stride;
      dest += destStride;
    }
  }
#else
#if ENABLE_SIMD_OPT_BUFFER && defined(TARGET_SIMD_X86)
  if( ( width & 7 ) == 0 )
  {
    g_pelBufOP.addAvg8( src0, src1Stride, src2, src2Stride, dest, destStride, width, height, shiftNum, offset, clpRng );
  }
  else if( ( width & 3 ) == 0 )
  {
    g_pelBufOP.addAvg4( src0, src1Stride, src2, src2Stride, dest, destStride, width, height, shiftNum, offset, clpRng );
  }
  else
#endif
  {
#define ADD_AVG_OP( ADDR ) dest[ADDR] = ClipPel( rightShift( ( src0[ADDR] + src2[ADDR] + offset ), shiftNum ), clpRng )
#define ADD_AVG_INC     \
    src0 += src1Stride; \
    src2 += src2Stride; \
    dest += destStride; \

    SIZE_AWARE_PER_EL_OP( ADD_AVG_OP, ADD_AVG_INC );

#undef ADD_AVG_OP
#undef ADD_AVG_INC
  }
#endif
}

template<>
void AreaBuf<Pel>::toLast( const ClpRng& clpRng )
{
        Pel* src       = buf;
  const uint32_t srcStride = stride;

  const int  clipbd    = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int shiftNum = IF_INTERNAL_FRAC_BITS(clipbd);
#else
  const int  shiftNum  = std::max<int>(2, (IF_INTERNAL_PREC - clipbd));
#endif
  const int  offset    = ( 1 << ( shiftNum - 1 ) ) + IF_INTERNAL_OFFS;

  if (width == 1)
  {
    THROW( "Blocks of width = 1 not supported" );
  }
  else if (width&2)
  {
    for ( int y = 0; y < height; y++ )
    {
      for (int x=0 ; x < width; x+=2 )
      {
        src[x + 0] = ClipPel( rightShift( ( src[x + 0] + offset ), shiftNum ), clpRng );
        src[x + 1] = ClipPel( rightShift( ( src[x + 1] + offset ), shiftNum ), clpRng );
      }
      src += srcStride;
    }
  }
  else
  {
    for ( int y = 0; y < height; y++ )
    {
      for (int x=0 ; x < width; x+=4 )
      {
        src[x + 0] = ClipPel( rightShift( ( src[x + 0] + offset ), shiftNum ), clpRng );
        src[x + 1] = ClipPel( rightShift( ( src[x + 1] + offset ), shiftNum ), clpRng );
        src[x + 2] = ClipPel( rightShift( ( src[x + 2] + offset ), shiftNum ), clpRng );
        src[x + 3] = ClipPel( rightShift( ( src[x + 3] + offset ), shiftNum ), clpRng );

      }
      src += srcStride;
    }
  }
}


template<>
void AreaBuf<Pel>::copyClip( const AreaBuf<const Pel> &src, const ClpRng& clpRng )
{
  const Pel* srcp = src.buf;
        Pel* dest =     buf;

  const unsigned srcStride  = src.stride;
  const unsigned destStride = stride;

#if !JVET_W0090_ARMC_TM && !JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  if( width == 1 )
  {
    THROW( "Blocks of width = 1 not supported" );
  }
  else
#endif
  {
#if JVET_W0097_GPM_MMVD_TM
    g_pelBufOP.copyClip(srcp, srcStride, dest, destStride, width, height, clpRng);
#else
#define RECO_OP( ADDR ) dest[ADDR] = ClipPel( srcp[ADDR], clpRng )
#define RECO_INC        \
    srcp += srcStride;  \
    dest += destStride; \

    SIZE_AWARE_PER_EL_OP( RECO_OP, RECO_INC );

#undef RECO_OP
#undef RECO_INC
#endif
  }
}

template<>
void AreaBuf<Pel>::roundToOutputBitdepth( const AreaBuf<const Pel> &src, const ClpRng& clpRng )
{
  const Pel* srcp = src.buf;
        Pel* dest =     buf;
  const unsigned srcStride  = src.stride;
  const unsigned destStride = stride;
#if !JVET_W0097_GPM_MMVD_TM
  const int32_t clipbd            = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int32_t shiftDefault      = IF_INTERNAL_FRAC_BITS(clipbd);
#else
  const int32_t shiftDefault      = std::max<int>(2, (IF_INTERNAL_PREC - clipbd));
#endif
  const int32_t offsetDefault     = (1<<(shiftDefault-1)) + IF_INTERNAL_OFFS;
#endif
  if( width == 1 )
  {
    THROW( "Blocks of width = 1 not supported" );
  }
  else
  {
#if JVET_W0097_GPM_MMVD_TM
    g_pelBufOP.roundBD(srcp, srcStride, dest, destStride, width, height, clpRng);
#else
#define RND_OP( ADDR ) dest[ADDR] = ClipPel( rightShift( srcp[ADDR] + offsetDefault, shiftDefault), clpRng )
#define RND_INC        \
    srcp += srcStride;  \
    dest += destStride; \

    SIZE_AWARE_PER_EL_OP( RND_OP, RND_INC );

#undef RND_OP
#undef RND_INC
#endif
  }
}


template<>
void AreaBuf<Pel>::reconstruct( const AreaBuf<const Pel> &pred, const AreaBuf<const Pel> &resi, const ClpRng& clpRng )
{
  const Pel* src1 = pred.buf;
  const Pel* src2 = resi.buf;
        Pel* dest =      buf;

  const unsigned src1Stride = pred.stride;
  const unsigned src2Stride = resi.stride;
  const unsigned destStride =      stride;

#if ENABLE_SIMD_OPT_BUFFER && defined(TARGET_SIMD_X86)
  if( ( width & 7 ) == 0 )
  {
    g_pelBufOP.reco8( src1, src1Stride, src2, src2Stride, dest, destStride, width, height, clpRng );
  }
  else if( ( width & 3 ) == 0 )
  {
    g_pelBufOP.reco4( src1, src1Stride, src2, src2Stride, dest, destStride, width, height, clpRng );
  }
  else
#endif
  {
#define RECO_OP( ADDR ) dest[ADDR] = ClipPel( src1[ADDR] + src2[ADDR], clpRng )
#define RECO_INC        \
    src1 += src1Stride; \
    src2 += src2Stride; \
    dest += destStride; \

    SIZE_AWARE_PER_EL_OP( RECO_OP, RECO_INC );

#undef RECO_OP
#undef RECO_INC
  }
}

template<>
void AreaBuf<Pel>::linearTransform( const int scale, const int shift, const int offset, bool bClip, const ClpRng& clpRng )
{
  const Pel* src = buf;
        Pel* dst = buf;

#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  if (width == 0)
  {
    THROW("Blocks of width = 0 not supported");
  }
#else
  if( width == 1 )
  {
    THROW( "Blocks of width = 1 not supported" );
  }
#endif
#if ENABLE_SIMD_OPT_BUFFER && defined(TARGET_SIMD_X86)
  else if( ( width & 7 ) == 0 )
  {
    g_pelBufOP.linTf8( src, stride, dst, stride, width, height, scale, shift, offset, clpRng, bClip );
  }
  else if( ( width & 3 ) == 0 )
  {
    g_pelBufOP.linTf4( src, stride, dst, stride, width, height, scale, shift, offset, clpRng, bClip );
  }
#endif
  else
  {
#define LINTF_OP( ADDR ) dst[ADDR] = ( Pel ) bClip ? ClipPel( rightShift( scale * src[ADDR], shift ) + offset, clpRng ) : ( rightShift( scale * src[ADDR], shift ) + offset )
#define LINTF_INC        \
    src += stride;       \
    dst += stride;       \

    SIZE_AWARE_PER_EL_OP( LINTF_OP, LINTF_INC );

#undef RECO_OP
#undef RECO_INC
  }
}

#if ENABLE_SIMD_OPT_BUFFER && defined(TARGET_SIMD_X86)
template<>
void AreaBuf<Pel>::subtract( const Pel val )
{
  ClpRng clpRngDummy;

  clpRngDummy.min = 0;
  clpRngDummy.max = 0;
  clpRngDummy.bd = 0;
  clpRngDummy.n = 0;

  linearTransform( 1, 0, -val, false, clpRngDummy );
}
#endif


PelStorage::PelStorage()
{
  for( uint32_t i = 0; i < MAX_NUM_COMPONENT; i++ )
  {
    m_origin[i] = nullptr;
  }
}

PelStorage::~PelStorage()
{
  destroy();
}

void PelStorage::create( const UnitArea &_UnitArea )
{
  create( _UnitArea.chromaFormat, _UnitArea.blocks[0] );
}

void PelStorage::create( const ChromaFormat &_chromaFormat, const Area& _area, const unsigned _maxCUSize, const unsigned _margin, const unsigned _alignment, const bool _scaleChromaMargin )
{
  CHECK( !bufs.empty(), "Trying to re-create an already initialized buffer" );

  chromaFormat = _chromaFormat;

  const uint32_t numCh = getNumberValidComponents( _chromaFormat );

  unsigned extHeight = _area.height;
  unsigned extWidth  = _area.width;

  if( _maxCUSize )
  {
    extHeight = ( ( _area.height + _maxCUSize - 1 ) / _maxCUSize ) * _maxCUSize;
    extWidth  = ( ( _area.width  + _maxCUSize - 1 ) / _maxCUSize ) * _maxCUSize;
  }

  for( uint32_t i = 0; i < numCh; i++ )
  {
    const ComponentID compID = ComponentID( i );
    const unsigned scaleX = ::getComponentScaleX( compID, _chromaFormat );
    const unsigned scaleY = ::getComponentScaleY( compID, _chromaFormat );

    unsigned scaledHeight = extHeight >> scaleY;
    unsigned scaledWidth  = extWidth  >> scaleX;
    unsigned ymargin      = _margin >> (_scaleChromaMargin?scaleY:0);
    unsigned xmargin      = _margin >> (_scaleChromaMargin?scaleX:0);
    unsigned totalWidth   = scaledWidth + 2*xmargin;
    unsigned totalHeight  = scaledHeight +2*ymargin;

    if( _alignment )
    {
      // make sure buffer lines are align
      CHECK( _alignment != MEMORY_ALIGN_DEF_SIZE, "Unsupported alignment" );
      totalWidth = ( ( totalWidth + _alignment - 1 ) / _alignment ) * _alignment;
    }
    uint32_t area = totalWidth * totalHeight;
    CHECK( !area, "Trying to create a buffer with zero area" );

    m_origin[i] = ( Pel* ) xMalloc( Pel, area );
    Pel* topLeft = m_origin[i] + totalWidth * ymargin + xmargin;
    bufs.push_back( PelBuf( topLeft, totalWidth, _area.width >> scaleX, _area.height >> scaleY ) );
  }
}

void PelStorage::createFromBuf( PelUnitBuf buf )
{
  chromaFormat = buf.chromaFormat;

  const uint32_t numCh = ::getNumberValidComponents( chromaFormat );

  bufs.resize(numCh);

  for( uint32_t i = 0; i < numCh; i++ )
  {
    PelBuf cPelBuf = buf.get( ComponentID( i ) );
    bufs[i] = PelBuf( cPelBuf.bufAt( 0, 0 ), cPelBuf.stride, cPelBuf.width, cPelBuf.height );
  }
}

void PelStorage::swap( PelStorage& other )
{
  const uint32_t numCh = ::getNumberValidComponents( chromaFormat );

  for( uint32_t i = 0; i < numCh; i++ )
  {
    // check this otherwise it would turn out to get very weird
    CHECK( chromaFormat                   != other.chromaFormat                  , "Incompatible formats" );
    CHECK( get( ComponentID( i ) )        != other.get( ComponentID( i ) )       , "Incompatible formats" );
    CHECK( get( ComponentID( i ) ).stride != other.get( ComponentID( i ) ).stride, "Incompatible formats" );

    std::swap( bufs[i].buf,    other.bufs[i].buf );
    std::swap( bufs[i].stride, other.bufs[i].stride );
    std::swap( m_origin[i],    other.m_origin[i] );
  }
}

void PelStorage::destroy()
{
  chromaFormat = NUM_CHROMA_FORMAT;
  for( uint32_t i = 0; i < MAX_NUM_COMPONENT; i++ )
  {
    if( m_origin[i] )
    {
      xFree( m_origin[i] );
      m_origin[i] = nullptr;
    }
  }
  bufs.clear();
}

PelBuf PelStorage::getBuf( const ComponentID CompID )
{
  return bufs[CompID];
}

const CPelBuf PelStorage::getBuf( const ComponentID CompID ) const
{
  return bufs[CompID];
}

PelBuf PelStorage::getBuf( const CompArea &blk )
{
  const PelBuf& r = bufs[blk.compID];

  CHECKD( rsAddr( blk.bottomRight(), r.stride ) >= ( ( r.height - 1 ) * r.stride + r.width ), "Trying to access a buf outside of bound!" );

  return PelBuf( r.buf + rsAddr( blk, r.stride ), r.stride, blk );
}

const CPelBuf PelStorage::getBuf( const CompArea &blk ) const
{
  const PelBuf& r = bufs[blk.compID];
  return CPelBuf( r.buf + rsAddr( blk, r.stride ), r.stride, blk );
}

PelUnitBuf PelStorage::getBuf( const UnitArea &unit )
{
  return ( chromaFormat == CHROMA_400 ) ? PelUnitBuf( chromaFormat, getBuf( unit.Y() ) ) : PelUnitBuf( chromaFormat, getBuf( unit.Y() ), getBuf( unit.Cb() ), getBuf( unit.Cr() ) );
}

const CPelUnitBuf PelStorage::getBuf( const UnitArea &unit ) const
{
  return ( chromaFormat == CHROMA_400 ) ? CPelUnitBuf( chromaFormat, getBuf( unit.Y() ) ) : CPelUnitBuf( chromaFormat, getBuf( unit.Y() ), getBuf( unit.Cb() ), getBuf( unit.Cr() ) );
}

template<>
void UnitBuf<Pel>::colorSpaceConvert(const UnitBuf<Pel> &other, const bool forward, const ClpRng& clpRng)
{
  const Pel* pOrg0 = bufs[COMPONENT_Y].buf;
  const Pel* pOrg1 = bufs[COMPONENT_Cb].buf;
  const Pel* pOrg2 = bufs[COMPONENT_Cr].buf;
  const int  strideOrg = bufs[COMPONENT_Y].stride;

  Pel* pDst0 = other.bufs[COMPONENT_Y].buf;
  Pel* pDst1 = other.bufs[COMPONENT_Cb].buf;
  Pel* pDst2 = other.bufs[COMPONENT_Cr].buf;
  const int strideDst = other.bufs[COMPONENT_Y].stride;

  int width = bufs[COMPONENT_Y].width;
  int height = bufs[COMPONENT_Y].height;
  int maxAbsclipBD = (1 << (clpRng.bd + 1)) - 1;
  int r, g, b;
  int y0, cg, co;

  CHECK(bufs[COMPONENT_Y].stride != bufs[COMPONENT_Cb].stride || bufs[COMPONENT_Y].stride != bufs[COMPONENT_Cr].stride, "unequal stride for 444 content");
  CHECK(other.bufs[COMPONENT_Y].stride != other.bufs[COMPONENT_Cb].stride || other.bufs[COMPONENT_Y].stride != other.bufs[COMPONENT_Cr].stride, "unequal stride for 444 content");
  CHECK(bufs[COMPONENT_Y].width != other.bufs[COMPONENT_Y].width || bufs[COMPONENT_Y].height != other.bufs[COMPONENT_Y].height, "unequal block size")

    if (forward)
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          r = pOrg2[x];
          g = pOrg0[x];
          b = pOrg1[x];

          co = r - b;
          int t = b + (co >> 1);
          cg = g - t;
          pDst0[x] = t + (cg >> 1);
          pDst1[x] = cg;
          pDst2[x] = co;
        }
        pOrg0 += strideOrg;
        pOrg1 += strideOrg;
        pOrg2 += strideOrg;
        pDst0 += strideDst;
        pDst1 += strideDst;
        pDst2 += strideDst;
      }
    }
    else
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          y0 = pOrg0[x];
          cg = pOrg1[x];
          co = pOrg2[x];

          y0 = Clip3((-maxAbsclipBD - 1), maxAbsclipBD, y0);
          cg = Clip3((-maxAbsclipBD - 1), maxAbsclipBD, cg);
          co = Clip3((-maxAbsclipBD - 1), maxAbsclipBD, co);

          int t = y0 - (cg >> 1);
          pDst0[x] = cg + t;
          pDst1[x] = t - (co >> 1);
          pDst2[x] = co + pDst1[x];
        }

        pOrg0 += strideOrg;
        pOrg1 += strideOrg;
        pOrg2 += strideOrg;
        pDst0 += strideDst;
        pDst1 += strideDst;
        pDst2 += strideDst;
      }
    }
}
