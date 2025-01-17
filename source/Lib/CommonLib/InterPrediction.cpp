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

/** \file     Prediction.cpp
    \brief    prediction class
*/

#include "InterPrediction.h"

#include "Buffer.h"
#include "UnitTools.h"
#include "MCTS.h"

#include <memory.h>
#include <algorithm>

#if INTER_LIC || (TM_AMVP || TM_MRG) || JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
#include "Reshape.h"
#endif

#if ENABLE_SIMD_TMP
#include "CommonDefX86.h"
#endif

//! \ingroup CommonLib
//! \{

// ====================================================================================================================
// Constructor / destructor / initialize
// ====================================================================================================================
#if JVET_Z0136_OOB
bool InterPrediction::isMvOOB(const Mv& rcMv, const struct Position pos, const struct Size size, const SPS* sps, const PPS* pps, bool *mcMask, bool *mcMaskChroma, bool lumaOnly)
{
  return g_pelBufOP.isMvOOB(rcMv, pos, size, sps, pps, mcMask, mcMaskChroma, lumaOnly, m_currChromaFormat);
}
bool InterPrediction::isMvOOBSubBlk(const Mv& rcMv, const struct Position pos, const struct Size size, const SPS* sps, const PPS* pps, bool *mcMask, int mcStride, bool *mcMaskChroma, int mcCStride, bool lumaOnly)
{
  return g_pelBufOP.isMvOOBSubBlk(rcMv, pos, size, sps, pps, mcMask, mcStride, mcMaskChroma, mcCStride, lumaOnly, m_currChromaFormat);
}
#endif

InterPrediction::InterPrediction()
:
#if INTER_LIC
  m_storeBeforeLIC  (false),
#endif
#if INTER_LIC || (TM_AMVP || TM_MRG) // note: already refactor
  m_pcReshape            ( nullptr ),
#endif
#if INTER_LIC || JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  m_pcLICRefLeftTemplate ( nullptr ),
  m_pcLICRefAboveTemplate( nullptr ),
  m_pcLICRecLeftTemplate ( nullptr ),
  m_pcLICRecAboveTemplate( nullptr ),
#endif
#if TM_AMVP || TM_MRG
  m_pcCurTplLeft ( nullptr ),
  m_pcCurTplAbove( nullptr ),
  m_pcRefTplLeft ( nullptr ),
  m_pcRefTplAbove( nullptr ),
#endif
  m_currChromaFormat( NUM_CHROMA_FORMAT )
, m_maxCompIDToPred ( MAX_NUM_COMPONENT )
, m_pcRdCost        ( nullptr )
, m_storedMv        ( nullptr )
, m_skipPROF (false)
, m_encOnly  (false)
, m_isBi     (false)
, m_gradX0(nullptr)
, m_gradY0(nullptr)
, m_gradX1(nullptr)
, m_gradY1(nullptr)
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
, m_absGx(nullptr)
, m_absGy(nullptr)
, m_dIx(nullptr)
, m_dIy(nullptr)
, m_dI(nullptr)
, m_signGxGy(nullptr)
, m_tmpx_pixel_32bit(nullptr)
, m_tmpy_pixel_32bit(nullptr)
, m_sumAbsGX_pixel_32bit(nullptr)
, m_sumAbsGY_pixel_32bit(nullptr)
, m_sumDIX_pixel_32bit(nullptr)
, m_sumDIY_pixel_32bit(nullptr)
, m_sumSignGY_GX_pixel_32bit(nullptr)
#endif
, m_subPuMC(false)
{
  for( uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++ )
  {
    for( uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++ )
    {
      m_acYuvPred[refList][ch] = nullptr;
    }
  }

  for( uint32_t c = 0; c < MAX_NUM_COMPONENT; c++ )
  {
    for( uint32_t i = 0; i < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL; i++ )
    {
      for( uint32_t j = 0; j < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL; j++ )
      {
        m_filteredBlock[i][j][c] = nullptr;
      }

      m_filteredBlockTmp[i][c] = nullptr;
    }
  }
  m_cYuvPredTempDMVRL1 = nullptr;
  m_cYuvPredTempDMVRL0 = nullptr;
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    m_cRefSamplesDMVRL0[ch] = nullptr;
    m_cRefSamplesDMVRL1[ch] = nullptr;
  }

#if INTER_LIC
  m_LICMultApprox[0] = 0;
  for (int k = 1; k < 64; k++)
  {
    m_LICMultApprox[k] = ((1 << 15) + (k >> 1)) / k;
  }
#endif

#if MULTI_PASS_DMVR
  int mvSearchIdx_bilMrg = 0;
#if JVET_X0049_BDMVR_SW_OPT
  uint16_t currtPrio = 0, currIdx = 0;
  ::memset(m_searchEnlargeOffsetNum, 0, sizeof(m_searchEnlargeOffsetNum));
#endif
  for (int y = -BDMVR_INTME_RANGE; y <= BDMVR_INTME_RANGE; y++)
  {
    for (int x = -BDMVR_INTME_RANGE; x <= BDMVR_INTME_RANGE; x++)
    {
#if JVET_X0049_BDMVR_SW_OPT
#else
      m_searchEnlargeOffsetBilMrg[mvSearchIdx_bilMrg] = Mv(x, y);
#endif
      if ( (abs(x) + abs(y)) == 0 )
      {
#if JVET_X0049_BDMVR_SW_OPT
        currtPrio = 0;
        currIdx = m_searchEnlargeOffsetNum[currtPrio];
        m_searchEnlargeOffsetToIdx[currtPrio][currIdx] = mvSearchIdx_bilMrg;
#else
        m_searchPriorityBilMrg[mvSearchIdx_bilMrg] = 0;
#endif
        m_costShiftBilMrg1[mvSearchIdx_bilMrg] = 63;
        m_costShiftBilMrg2[mvSearchIdx_bilMrg++] = 63;
      }
      else if ( (abs(x) + abs(y)) < 4 )
      {
#if JVET_X0049_BDMVR_SW_OPT
        currtPrio = 1;
        currIdx = m_searchEnlargeOffsetNum[currtPrio];
        m_searchEnlargeOffsetToIdx[currtPrio][currIdx] = mvSearchIdx_bilMrg;
#else
        m_searchPriorityBilMrg[mvSearchIdx_bilMrg] = 1;
#endif
        m_costShiftBilMrg1[mvSearchIdx_bilMrg] = 63;
        m_costShiftBilMrg2[mvSearchIdx_bilMrg++] = 63;
      }
      else if ((abs(x) + abs(y)) < 7)
      {
#if JVET_X0049_BDMVR_SW_OPT
        currtPrio = 2;
        currIdx = m_searchEnlargeOffsetNum[currtPrio];
        m_searchEnlargeOffsetToIdx[currtPrio][currIdx] = mvSearchIdx_bilMrg;
#else
        m_searchPriorityBilMrg[mvSearchIdx_bilMrg] = 2;
#endif
        m_costShiftBilMrg1[mvSearchIdx_bilMrg] = 2;
        m_costShiftBilMrg2[mvSearchIdx_bilMrg++] = 63;
      }
      else if ((abs(x) + abs(y)) < 11)
      {
#if JVET_X0049_BDMVR_SW_OPT
        currtPrio = 3;
        currIdx = m_searchEnlargeOffsetNum[currtPrio];
        m_searchEnlargeOffsetToIdx[currtPrio][currIdx] = mvSearchIdx_bilMrg;
#else
        m_searchPriorityBilMrg[mvSearchIdx_bilMrg] = 3;
#endif
        m_costShiftBilMrg1[mvSearchIdx_bilMrg] = 1;
        m_costShiftBilMrg2[mvSearchIdx_bilMrg++] = 63;
      }
      else
      {
#if JVET_X0049_BDMVR_SW_OPT
        currtPrio = 4;
        currIdx = m_searchEnlargeOffsetNum[currtPrio];
        m_searchEnlargeOffsetToIdx[currtPrio][currIdx] = mvSearchIdx_bilMrg;
#else
        m_searchPriorityBilMrg[mvSearchIdx_bilMrg] = 4;
#endif
        m_costShiftBilMrg1[mvSearchIdx_bilMrg] = 1;
        m_costShiftBilMrg2[mvSearchIdx_bilMrg++] = 2;
      }
#if JVET_X0049_BDMVR_SW_OPT
      m_searchEnlargeOffsetBilMrg[currtPrio][currIdx] = Mv(x, y);
      m_searchEnlargeOffsetNum[currtPrio]++;
#endif
    }
  }
  CHECK(mvSearchIdx_bilMrg != (2 * BDMVR_INTME_RANGE + 1) * (2 * BDMVR_INTME_RANGE + 1),
      "this is wrong, mvSearchIdx_bilMrg != (2 * BDMVR_INTME_RANGE + 1) * (2 * BDMVR_INTME_RANGE + 1)");
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    for (uint32_t tmplt = 0; tmplt < 2; tmplt++)
    {
      m_acYuvCurAMLTemplate[tmplt][ch] = nullptr;
      m_acYuvRefAboveTemplate[tmplt][ch] = nullptr;
      m_acYuvRefLeftTemplate[tmplt][ch] = nullptr;
      m_acYuvRefAMLTemplate[tmplt][ch] = nullptr;
    }
  }
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  for (uint32_t tmplt = 0; tmplt < 2 + 2; tmplt++)
  {
    m_acYuvRefAMLTemplatePart0[tmplt] = nullptr;
    m_acYuvRefAMLTemplatePart1[tmplt] = nullptr;
  }
  m_tplWeightTblInitialized = false;
#endif
#endif
#if JVET_Z0061_TM_OBMC
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    for (uint32_t tmplt = 0; tmplt < 2; tmplt++)
    {
      m_acYuvRefAboveTemplateOBMC[tmplt][ch] = nullptr;
      m_acYuvRefLeftTemplateOBMC[tmplt][ch]  = nullptr;
      m_acYuvBlendTemplateOBMC[tmplt][ch]    = nullptr;
    }
  }
#endif
}

InterPrediction::~InterPrediction()
{
  destroy();
}

void InterPrediction::destroy()
{
  for( uint32_t i = 0; i < NUM_REF_PIC_LIST_01; i++ )
  {
    for( uint32_t c = 0; c < MAX_NUM_COMPONENT; c++ )
    {
      xFree( m_acYuvPred[i][c] );
      m_acYuvPred[i][c] = nullptr;
    }
  }

  for( uint32_t c = 0; c < MAX_NUM_COMPONENT; c++ )
  {
    for( uint32_t i = 0; i < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL; i++ )
    {
      for( uint32_t j = 0; j < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL; j++ )
      {
        xFree( m_filteredBlock[i][j][c] );
        m_filteredBlock[i][j][c] = nullptr;
      }

      xFree( m_filteredBlockTmp[i][c] );
      m_filteredBlockTmp[i][c] = nullptr;
    }
  }

  m_geoPartBuf[0].destroy();
  m_geoPartBuf[1].destroy();
  m_colorTransResiBuf[0].destroy();
  m_colorTransResiBuf[1].destroy();
  m_colorTransResiBuf[2].destroy();

  if (m_storedMv != nullptr)
  {
    delete[]m_storedMv;
    m_storedMv = nullptr;
  }

  xFree(m_gradX0);   m_gradX0 = nullptr;
  xFree(m_gradY0);   m_gradY0 = nullptr;
  xFree(m_gradX1);   m_gradX1 = nullptr;
  xFree(m_gradY1);   m_gradY1 = nullptr;
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  xFree(m_absGx);    m_absGx = nullptr;
  xFree(m_absGy);    m_absGy = nullptr;
  xFree(m_dIx);      m_dIx = nullptr;
  xFree(m_dIy);      m_dIy = nullptr;
  xFree(m_dI);    m_dI = nullptr;
  xFree(m_signGxGy); m_signGxGy = nullptr;
  xFree(m_tmpx_pixel_32bit); m_tmpx_pixel_32bit = nullptr;
  xFree(m_tmpy_pixel_32bit); m_tmpy_pixel_32bit = nullptr;
  xFree(m_sumAbsGX_pixel_32bit);     m_sumAbsGX_pixel_32bit = nullptr;
  xFree(m_sumAbsGY_pixel_32bit);     m_sumAbsGY_pixel_32bit = nullptr;
  xFree(m_sumDIX_pixel_32bit);       m_sumDIX_pixel_32bit = nullptr;
  xFree(m_sumDIY_pixel_32bit);       m_sumDIY_pixel_32bit = nullptr;
  xFree(m_sumSignGY_GX_pixel_32bit); m_sumSignGY_GX_pixel_32bit = nullptr;
#endif
#if ENABLE_OBMC
  m_tmpObmcBufL0.destroy();
  m_tmpObmcBufT0.destroy();
  m_tmpSubObmcBuf.destroy();
#endif
  xFree(m_cYuvPredTempDMVRL0);
  m_cYuvPredTempDMVRL0 = nullptr;
  xFree(m_cYuvPredTempDMVRL1);
  m_cYuvPredTempDMVRL1 = nullptr;
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    xFree(m_cRefSamplesDMVRL0[ch]);
    m_cRefSamplesDMVRL0[ch] = nullptr;
    xFree(m_cRefSamplesDMVRL1[ch]);
    m_cRefSamplesDMVRL1[ch] = nullptr;
  }
#if JVET_Z0118_GDR
  m_IBCBuffer0.destroy();
  m_IBCBuffer1.destroy();
#else
  m_IBCBuffer.destroy();
#endif

#if TM_AMVP || TM_MRG
  xFree(m_pcCurTplLeft ); m_pcCurTplLeft  = nullptr;
  xFree(m_pcCurTplAbove); m_pcCurTplAbove = nullptr;
  xFree(m_pcRefTplLeft ); m_pcRefTplLeft  = nullptr;
  xFree(m_pcRefTplAbove); m_pcRefTplAbove = nullptr;
#endif
#if INTER_LIC || JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  xFree(m_pcLICRefLeftTemplate);  m_pcLICRefLeftTemplate  = nullptr;
  xFree(m_pcLICRefAboveTemplate); m_pcLICRefAboveTemplate = nullptr;
  xFree(m_pcLICRecLeftTemplate);  m_pcLICRecLeftTemplate  = nullptr;
  xFree(m_pcLICRecAboveTemplate); m_pcLICRecAboveTemplate = nullptr;
#endif
#if MULTI_HYP_PRED
  m_additionalHypothesisStorage.destroy();
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    for (uint32_t tmplt = 0; tmplt < 2; tmplt++)
    {
      xFree(m_acYuvCurAMLTemplate[tmplt][ch]);
      xFree(m_acYuvRefAboveTemplate[tmplt][ch]);
      xFree(m_acYuvRefLeftTemplate[tmplt][ch]);
      xFree(m_acYuvRefAMLTemplate[tmplt][ch]);
      m_acYuvCurAMLTemplate[tmplt][ch] = nullptr;
      m_acYuvRefAboveTemplate[tmplt][ch] = nullptr;
      m_acYuvRefLeftTemplate[tmplt][ch] = nullptr;
      m_acYuvRefAMLTemplate[tmplt][ch] = nullptr;
  }
}
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  for (uint32_t tmplt = 0; tmplt < 2 + 2; tmplt++)
  {
    xFree(m_acYuvRefAMLTemplatePart0[tmplt]);
    xFree(m_acYuvRefAMLTemplatePart1[tmplt]);
    m_acYuvRefAMLTemplatePart0[tmplt] = nullptr;
    m_acYuvRefAMLTemplatePart1[tmplt] = nullptr;
  }
#endif
#endif
#if JVET_Z0061_TM_OBMC
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    for (uint32_t tmplt = 0; tmplt < 2; tmplt++)
    {
      xFree(m_acYuvRefAboveTemplateOBMC[tmplt][ch]);
      xFree(m_acYuvRefLeftTemplateOBMC[tmplt][ch]);
      xFree(m_acYuvBlendTemplateOBMC[tmplt][ch]);
      m_acYuvRefAboveTemplateOBMC[tmplt][ch] = nullptr;
      m_acYuvRefLeftTemplateOBMC[tmplt][ch]  = nullptr;
      m_acYuvBlendTemplateOBMC[tmplt][ch]    = nullptr;
    }
  }
#endif
}

#if INTER_LIC || (TM_AMVP || TM_MRG) || JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
#if JVET_Z0153_IBC_EXT_REF
void InterPrediction::init( RdCost* pcRdCost, ChromaFormat chromaFormatIDC, const int ctuSize, Reshape* reshape, const int picWidth )
#else
void InterPrediction::init( RdCost* pcRdCost, ChromaFormat chromaFormatIDC, const int ctuSize, Reshape* reshape )
#endif
#else
void InterPrediction::init( RdCost* pcRdCost, ChromaFormat chromaFormatIDC, const int ctuSize )
#endif
{
  m_pcRdCost = pcRdCost;
#if INTER_LIC || (TM_AMVP || TM_MRG) || JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
  m_pcReshape = reshape;
#endif

  // if it has been initialised before, but the chroma format has changed, release the memory and start again.
  if( m_acYuvPred[REF_PIC_LIST_0][COMPONENT_Y] != nullptr && m_currChromaFormat != chromaFormatIDC )
  {
    destroy();
  }

  m_currChromaFormat = chromaFormatIDC;
  if( m_acYuvPred[REF_PIC_LIST_0][COMPONENT_Y] == nullptr ) // check if first is null (in which case, nothing initialised yet)
  {
    for( uint32_t c = 0; c < MAX_NUM_COMPONENT; c++ )
    {
#if IF_12TAP || MULTI_PASS_DMVR
#if MULTI_PASS_DMVR
      int extendSize = std::max(2 * BIO_EXTEND_SIZE + 2, 2 * BDMVR_INTME_RANGE);
#else
      int extendSize = std::max(2 * BIO_EXTEND_SIZE + 2, 2 * DMVR_NUM_ITERATION);
#endif
#if IF_12TAP
      int extWidth   = MAX_CU_SIZE + extendSize + 32;
#else
      int extWidth   = MAX_CU_SIZE + extendSize + 16;
#endif
      int extHeight  = MAX_CU_SIZE + extendSize + 1;
#else
      int extWidth = MAX_CU_SIZE + (2 * BIO_EXTEND_SIZE + 2) + 16;
      int extHeight = MAX_CU_SIZE + (2 * BIO_EXTEND_SIZE + 2) + 1;
      extWidth = extWidth > (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + 16) ? extWidth : MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + 16;
      extHeight = extHeight > (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + 1) ? extHeight : MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + 1;
#endif
      for( uint32_t i = 0; i < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL; i++ )
      {
#if IF_12TAP
        m_filteredBlockTmp[i][c] = ( Pel* ) xMalloc( Pel, ( extWidth + 4 ) * ( extHeight + 15 + 4 ) );
#else
        m_filteredBlockTmp[i][c] = ( Pel* ) xMalloc( Pel, ( extWidth + 4 ) * ( extHeight + 7 + 4 ) );
#endif

        for( uint32_t j = 0; j < LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL; j++ )
        {
          m_filteredBlock[i][j][c] = ( Pel* ) xMalloc( Pel, extWidth * extHeight );
        }
      }

      // new structure
      for( uint32_t i = 0; i < NUM_REF_PIC_LIST_01; i++ )
      {
        m_acYuvPred[i][c] = ( Pel* ) xMalloc( Pel, MAX_CU_SIZE * MAX_CU_SIZE );
      }
    }

    m_geoPartBuf[0].create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
    m_geoPartBuf[1].create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
    m_colorTransResiBuf[0].create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
    m_colorTransResiBuf[1].create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
    m_colorTransResiBuf[2].create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
#if MULTI_HYP_PRED
    m_additionalHypothesisStorage.create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
#endif

    m_iRefListIdx = -1;

    m_gradX0 = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_gradY0 = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_gradX1 = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_gradY1 = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
    m_absGx = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_absGy = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_dIx = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_dIy = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_dI = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_signGxGy = (Pel*)xMalloc(Pel, BIO_TEMP_BUFFER_SIZE);
    m_tmpx_pixel_32bit = (int*)xMalloc(int, BDOF_SUBPU_SIZE);
    m_tmpy_pixel_32bit = (int*)xMalloc(int, BDOF_SUBPU_SIZE);
    m_sumAbsGX_pixel_32bit = (int*)xMalloc(int, BDOF_SUBPU_SIZE);
    m_sumAbsGY_pixel_32bit = (int*)xMalloc(int, BDOF_SUBPU_SIZE);
    m_sumDIX_pixel_32bit = (int*)xMalloc(int, BDOF_SUBPU_SIZE);
    m_sumDIY_pixel_32bit = (int*)xMalloc(int, BDOF_SUBPU_SIZE);
    m_sumSignGY_GX_pixel_32bit = (int*)xMalloc(int, BDOF_SUBPU_SIZE);
#endif
#if ENABLE_OBMC
    m_tmpObmcBufL0.create(UnitArea(chromaFormatIDC, Area(0, 0, 4, MAX_CU_SIZE)));
    m_tmpObmcBufT0.create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, 4)));
    m_tmpSubObmcBuf.create(UnitArea(chromaFormatIDC, Area(0, 0, 20, 4)));
    m_tmpSubObmcBuf.bufs[0].memset(0);
    m_tmpSubObmcBuf.bufs[1].memset(0);
    m_tmpSubObmcBuf.bufs[2].memset(0);
#endif
  }

  if (m_cYuvPredTempDMVRL0 == nullptr && m_cYuvPredTempDMVRL1 == nullptr)
  {
    m_cYuvPredTempDMVRL0 = (Pel*)xMalloc(Pel, (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION)) * (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION)));
    m_cYuvPredTempDMVRL1 = (Pel*)xMalloc(Pel, (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION)) * (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION)));
    for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
    {
#if IF_12TAP
      m_cRefSamplesDMVRL0[ch] = (Pel*)xMalloc(Pel, (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA(0)) * (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA(0)));
      m_cRefSamplesDMVRL1[ch] = (Pel*)xMalloc(Pel, (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA(0)) * (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA(0)));
#else
      m_cRefSamplesDMVRL0[ch] = (Pel*)xMalloc(Pel, (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA) * (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA));
      m_cRefSamplesDMVRL1[ch] = (Pel*)xMalloc(Pel, (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA) * (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA));
#endif
    }
  }
#if !JVET_J0090_MEMORY_BANDWITH_MEASURE
  m_if.initInterpolationFilter( true );
#endif
#if TM_AMVP || TM_MRG
  if (m_pcCurTplLeft == nullptr)
  {
    m_pcCurTplLeft  = (Pel*)xMalloc(Pel, TM_TPL_SIZE * MAX_CU_SIZE);
    m_pcCurTplAbove = (Pel*)xMalloc(Pel, TM_TPL_SIZE * MAX_CU_SIZE);
    m_pcRefTplLeft  = (Pel*)xMalloc(Pel, TM_TPL_SIZE * MAX_CU_SIZE);
    m_pcRefTplAbove = (Pel*)xMalloc(Pel, TM_TPL_SIZE * MAX_CU_SIZE);
  }
#endif
#if INTER_LIC || JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  if (m_pcLICRefLeftTemplate == nullptr)
  {
    m_pcLICRefLeftTemplate  = (Pel*)xMalloc(Pel, MAX_CU_SIZE);
    m_pcLICRefAboveTemplate = (Pel*)xMalloc(Pel, MAX_CU_SIZE);
    m_pcLICRecLeftTemplate  = (Pel*)xMalloc(Pel, MAX_CU_SIZE);
    m_pcLICRecAboveTemplate = (Pel*)xMalloc(Pel, MAX_CU_SIZE);
  }
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    for (uint32_t tmplt = 0; tmplt < 2; tmplt++)
    {
      m_acYuvCurAMLTemplate[tmplt][ch] = (Pel*)xMalloc(Pel, MAX_CU_SIZE * MAX_CU_SIZE);
      m_acYuvRefAboveTemplate[tmplt][ch] = (Pel*)xMalloc(Pel, MAX_CU_SIZE * MAX_CU_SIZE);
      m_acYuvRefLeftTemplate[tmplt][ch] = (Pel*)xMalloc(Pel, MAX_CU_SIZE * MAX_CU_SIZE);
      m_acYuvRefAMLTemplate[tmplt][ch] = (Pel*)xMalloc(Pel, MAX_CU_SIZE * MAX_CU_SIZE);
    }
  }
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  for (uint32_t tmplt = 0; tmplt < 2 + 2; tmplt++)
  {
    m_acYuvRefAMLTemplatePart0[tmplt] = (Pel*)xMalloc(Pel, GEO_MAX_CU_SIZE * GEO_MODE_SEL_TM_SIZE);
    m_acYuvRefAMLTemplatePart1[tmplt] = (Pel*)xMalloc(Pel, GEO_MAX_CU_SIZE * GEO_MODE_SEL_TM_SIZE);
  }
#endif
#endif
#if JVET_Z0061_TM_OBMC
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++)
  {
    for (uint32_t tmplt = 0; tmplt < 2; tmplt++)
    {
      m_acYuvRefAboveTemplateOBMC[tmplt][ch] = (Pel *) xMalloc(Pel, MAX_CU_SIZE * MAX_CU_SIZE);
      m_acYuvRefLeftTemplateOBMC[tmplt][ch]  = (Pel *) xMalloc(Pel, MAX_CU_SIZE * MAX_CU_SIZE);
      m_acYuvBlendTemplateOBMC[tmplt][ch]    = (Pel *) xMalloc(Pel, MAX_CU_SIZE * MAX_CU_SIZE);
    }
  }
#endif

  if (m_storedMv == nullptr)
  {
    const int MVBUFFER_SIZE = MAX_CU_SIZE / MIN_PU_SIZE;
    m_storedMv = new Mv[MVBUFFER_SIZE*MVBUFFER_SIZE];
  }

#if JVET_Z0118_GDR
#if JVET_Z0153_IBC_EXT_REF
  m_IBCBufferWidth  = (picWidth + ctuSize - 1) / ctuSize * ctuSize;
  m_IBCBufferHeight = 3 * ctuSize;

  if (m_IBCBuffer0.bufs.empty())
  {    
    m_IBCBuffer0.create(UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, m_IBCBufferHeight)));
  }
  if (m_IBCBuffer1.bufs.empty())
  {   
    m_IBCBuffer1.create(UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, m_IBCBufferHeight)));
  }  
#else
  m_IBCBufferWidth = g_IBCBufferSize / ctuSize;

  if (m_IBCBuffer0.bufs.empty())
  {    
    m_IBCBuffer0.create(UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, ctuSize)));
  }
  if (m_IBCBuffer1.bufs.empty())
  {   
    m_IBCBuffer1.create(UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, ctuSize)));
  }
#endif
#else
  if (m_IBCBuffer.bufs.empty())
  {
#if JVET_Z0153_IBC_EXT_REF
    m_IBCBufferWidth = (picWidth + ctuSize - 1) / ctuSize * ctuSize;
    m_IBCBufferHeight = 3 * ctuSize;
    m_IBCBuffer.create(UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, m_IBCBufferHeight)));
#else
    m_IBCBufferWidth = g_IBCBufferSize / ctuSize;
    m_IBCBuffer.create(UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, ctuSize)));
#endif
  }
#endif
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

bool InterPrediction::xCheckIdenticalMotion( const PredictionUnit &pu )
{
  const Slice &slice = *pu.cs->slice;

  if( slice.isInterB() && !pu.cs->pps->getWPBiPred() )
  {
    if( pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0 )
    {
      int RefPOCL0 = slice.getRefPic( REF_PIC_LIST_0, pu.refIdx[0] )->getPOC();
      int RefPOCL1 = slice.getRefPic( REF_PIC_LIST_1, pu.refIdx[1] )->getPOC();

      if( RefPOCL0 == RefPOCL1 )
      {
        if( !pu.cu->affine )
        {
          if( pu.mv[0] == pu.mv[1] )
          {
            return true;
          }
        }
        else
        {
          if ( (pu.cu->affineType == AFFINEMODEL_4PARAM && (pu.mvAffi[0][0] == pu.mvAffi[1][0]) && (pu.mvAffi[0][1] == pu.mvAffi[1][1]))
            || (pu.cu->affineType == AFFINEMODEL_6PARAM && (pu.mvAffi[0][0] == pu.mvAffi[1][0]) && (pu.mvAffi[0][1] == pu.mvAffi[1][1]) && (pu.mvAffi[0][2] == pu.mvAffi[1][2])) )
          {
            return true;
          }
        }
      }
    }
  }

  return false;
}

void InterPrediction::xSubPuMC( PredictionUnit& pu, PelUnitBuf& predBuf, const RefPicList &eRefPicList /*= REF_PIC_LIST_X*/, const bool luma /*= true*/, const bool chroma /*= true*/)
{
#if MULTI_HYP_PRED
  CHECK(!pu.addHypData.empty(), "Multi Hyp: !pu.addHypData.empty()");
#endif
  // compute the location of the current PU
  Position puPos    = pu.lumaPos();
  Size puSize       = pu.lumaSize();

  int numPartLine, numPartCol, puHeight, puWidth;
  {
    numPartLine = std::max(puSize.width >> ATMVP_SUB_BLOCK_SIZE, 1u);
    numPartCol = std::max(puSize.height >> ATMVP_SUB_BLOCK_SIZE, 1u);
    puHeight = numPartCol == 1 ? puSize.height : 1 << ATMVP_SUB_BLOCK_SIZE;
    puWidth = numPartLine == 1 ? puSize.width : 1 << ATMVP_SUB_BLOCK_SIZE;
  }

  PredictionUnit subPu;

  subPu.cs        = pu.cs;
  subPu.cu        = pu.cu;
  subPu.mergeType = MRG_TYPE_DEFAULT_N;

  bool isAffine = pu.cu->affine;
  subPu.cu->affine = false;

  // join sub-pus containing the same motion
  bool verMC = puSize.height > puSize.width;
  int  fstStart = (!verMC ? puPos.y : puPos.x);
  int  secStart = (!verMC ? puPos.x : puPos.y);
  int  fstEnd = (!verMC ? puPos.y + puSize.height : puPos.x + puSize.width);
  int  secEnd = (!verMC ? puPos.x + puSize.width : puPos.y + puSize.height);
  int  fstStep = (!verMC ? puHeight : puWidth);
  int  secStep = (!verMC ? puWidth : puHeight);

  const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();

  const bool scaled = isResamplingPossible && ( pu.cu->slice->getRefPic( REF_PIC_LIST_0, 0 )->isRefScaled( pu.cs->pps ) || ( pu.cs->slice->getSliceType() == B_SLICE ? pu.cu->slice->getRefPic( REF_PIC_LIST_1, 0 )->isRefScaled( pu.cs->pps ) : false ) );
  m_subPuMC = true;

  for (int fstDim = fstStart; fstDim < fstEnd; fstDim += fstStep)
  {
    for (int secDim = secStart; secDim < secEnd; secDim += secStep)
    {
      int x = !verMC ? secDim : fstDim;
      int y = !verMC ? fstDim : secDim;
      const MotionInfo &curMi = pu.getMotionInfo(Position{ x, y });

      int length = secStep;
      int later  = secDim + secStep;

      while (later < secEnd)
      {
        const MotionInfo &laterMi = !verMC ? pu.getMotionInfo(Position{ later, fstDim }) : pu.getMotionInfo(Position{ fstDim, later });
        if (!scaled && laterMi == curMi
#if INTER_LIC
          && laterMi.usesLIC == curMi.usesLIC
#endif
          )
        {
          length += secStep;
        }
        else
        {
          break;
        }
        later += secStep;
      }
      int dx = !verMC ? length : puWidth;
      int dy = !verMC ? puHeight : length;

      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, dy)));
      subPu = curMi;
      PelUnitBuf subPredBuf = predBuf.subBuf(UnitAreaRelative(pu, subPu));
      subPu.mmvdEncOptMode = 0;
      subPu.mvRefine = false;
      motionCompensation(subPu, subPredBuf, eRefPicList, luma, chroma);
      secDim = later - secStep;
    }
  }
  m_subPuMC = false;

  pu.cu->affine = isAffine;
}
#if !BDOF_RM_CONSTRAINTS
void InterPrediction::xSubPuBio(PredictionUnit& pu, PelUnitBuf& predBuf, const RefPicList &eRefPicList /*= REF_PIC_LIST_X*/, PelUnitBuf* yuvDstTmp /*= NULL*/)
{
  // compute the location of the current PU
  Position puPos = pu.lumaPos();
  Size puSize = pu.lumaSize();

#if JVET_J0090_MEMORY_BANDWITH_MEASURE
  JVET_J0090_SET_CACHE_ENABLE(true);
  int mvShift = (MV_FRACTIONAL_BITS_INTERNAL);
  for (int k = 0; k < NUM_REF_PIC_LIST_01; k++)
  {
    RefPicList refId = (RefPicList)k;
    const Picture* refPic = pu.cu->slice->getRefPic(refId, pu.refIdx[refId]);
    for (int compID = 0; compID < MAX_NUM_COMPONENT; compID++)
    {
      Mv cMv = pu.mv[refId];
      int mvshiftTemp = mvShift + getComponentScaleX((ComponentID)compID, pu.chromaFormat);
      int filtersize = (compID == (COMPONENT_Y)) ? NTAPS_LUMA : NTAPS_CHROMA;
      cMv += Mv(-(((filtersize >> 1) - 1) << mvshiftTemp), -(((filtersize >> 1) - 1) << mvshiftTemp));
      bool wrapRef = false;
      if ( pu.cu->slice->getRefPic(refId, pu.refIdx[refId])->isWrapAroundEnabled( pu.cs->pps ) )
      {
        wrapRef = wrapClipMv(cMv, pu.blocks[0].pos(), pu.blocks[0].size(), pu.cs->sps, pu.cs->pps);
      }
      else
      {
        clipMv(cMv, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);
      }

      int width = predBuf.bufs[compID].width + (filtersize - 1);
      int height = predBuf.bufs[compID].height + (filtersize - 1);

      CPelBuf refBuf;
      Position recOffset = pu.blocks[compID].pos().offset(cMv.getHor() >> mvshiftTemp, cMv.getVer() >> mvshiftTemp);
      refBuf = refPic->getRecoBuf(CompArea((ComponentID)compID, pu.chromaFormat, recOffset, pu.blocks[compID].size()), wrapRef);

      JVET_J0090_SET_REF_PICTURE(refPic, (ComponentID)compID);
      for (int row = 0; row < height; row++)
      {
        for (int col = 0; col < width; col++)
        {
          JVET_J0090_CACHE_ACCESS(((Pel *)refBuf.buf) + row * refBuf.stride + col, __FILE__, __LINE__);
        }
      }
    }
  }
  JVET_J0090_SET_CACHE_ENABLE(false);
#endif
  PredictionUnit subPu;

  subPu.cs = pu.cs;
  subPu.cu = pu.cu;
  subPu.mergeType = pu.mergeType;
  subPu.mmvdMergeFlag = pu.mmvdMergeFlag;
  subPu.mmvdEncOptMode = pu.mmvdEncOptMode;
  subPu.mergeFlag = pu.mergeFlag;
  subPu.ciipFlag = pu.ciipFlag;
  subPu.mvRefine = pu.mvRefine;
#if TM_MRG
  subPu.tmMergeFlag = pu.tmMergeFlag;
#endif
  subPu.refIdx[0] = pu.refIdx[0];
  subPu.refIdx[1] = pu.refIdx[1];
  int  fstStart = puPos.y;
  int  secStart = puPos.x;
  int  fstEnd = puPos.y + puSize.height;
  int  secEnd = puPos.x + puSize.width;
  int  fstStep = std::min((int)MAX_BDOF_APPLICATION_REGION, (int)puSize.height);
  int  secStep = std::min((int)MAX_BDOF_APPLICATION_REGION, (int)puSize.width);
  for (int fstDim = fstStart; fstDim < fstEnd; fstDim += fstStep)
  {
    for (int secDim = secStart; secDim < secEnd; secDim += secStep)
    {
      int x = secDim;
      int y = fstDim;
      int dx = secStep;
      int dy = fstStep;

#if !JVET_W0097_GPM_MMVD_TM
      const MotionInfo &curMi = pu.getMotionInfo(Position{ x, y });
#endif

      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, dy)));
#if JVET_W0097_GPM_MMVD_TM
      subPu.interDir = pu.interDir;
      for (uint32_t i = 0; i < NUM_REF_PIC_LIST_01; i++)
      {
        subPu.refIdx[i] = pu.refIdx[i];
        subPu.mv[i] = pu.mv[i];
      }
#else
      subPu = curMi;
#endif
      PelUnitBuf subPredBuf = predBuf.subBuf(UnitAreaRelative(pu, subPu));

      if (yuvDstTmp)
      {
        PelUnitBuf subPredBufTmp = yuvDstTmp->subBuf(UnitAreaRelative(pu, subPu));
        motionCompensation(subPu, subPredBuf, eRefPicList, true, true, &subPredBufTmp);
      }
      else
      motionCompensation(subPu, subPredBuf, eRefPicList);
    }
  }
  JVET_J0090_SET_CACHE_ENABLE(true);
}
#endif

#if MULTI_PASS_DMVR
void InterPrediction::xPredInterUni(const PredictionUnit &pu, const RefPicList &eRefPicList, PelUnitBuf &pcYuvPred, const bool &bi,
                                    const bool &bioApplied, const bool luma, const bool chroma, const bool isBdofMvRefine)
#else
void InterPrediction::xPredInterUni(const PredictionUnit &pu, const RefPicList &eRefPicList, PelUnitBuf &pcYuvPred,
                                    const bool &bi, const bool &bioApplied, const bool luma, const bool chroma)
#endif
{
  const SPS &sps = *pu.cs->sps;

  int iRefIdx = pu.refIdx[eRefPicList];
  Mv mv[3];
  bool isIBC = false;
#if !INTER_RM_SIZE_CONSTRAINTS
#if ENABLE_OBMC
  if (pu.cu->isobmcMC == false)
#endif
  CHECK( !CU::isIBC( *pu.cu ) && pu.lwidth() == 4 && pu.lheight() == 4, "invalid 4x4 inter blocks" );
#endif
  if (CU::isIBC(*pu.cu))
  {
    isIBC = true;
  }
  if( pu.cu->affine )
  {
    CHECK( iRefIdx < 0, "iRefIdx incorrect." );

    mv[0] = pu.mvAffi[eRefPicList][0];
    mv[1] = pu.mvAffi[eRefPicList][1];
    mv[2] = pu.mvAffi[eRefPicList][2];
  }
  else
  {
    mv[0] = pu.mv[eRefPicList];
  }

  if( !pu.cu->affine )
  {
    const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();

    if( !isIBC && ( !isResamplingPossible || !pu.cu->slice->getRefPic( eRefPicList, iRefIdx )->isRefScaled( pu.cs->pps ) ) )
    {
      if( !pu.cs->pps->getWrapAroundEnabledFlag() )
      {
        clipMv( mv[0], pu.cu->lumaPos(), pu.cu->lumaSize(), sps, *pu.cs->pps );
      }
    }
  }

  for( uint32_t comp = COMPONENT_Y; comp < pcYuvPred.bufs.size() && comp <= m_maxCompIDToPred; comp++ )
  {
    const ComponentID compID = ComponentID( comp );
    if (compID == COMPONENT_Y && !luma)
    {
      continue;
    }
    if (compID != COMPONENT_Y && !chroma)
    {
      continue;
    }
#if MULTI_PASS_DMVR
    if (compID != COMPONENT_Y && bioApplied && isBdofMvRefine)
    {
      continue;
    }
#endif
    if ( pu.cu->affine )
    {
      CHECK( bioApplied, "BIO is not allowed with affine" );
      m_iRefListIdx = eRefPicList;
      bool genChromaMv = (!luma && chroma && compID == COMPONENT_Cb);
#if JVET_Z0136_OOB
      xPredAffineBlk(compID, pu, pu.cu->slice->getRefPic(eRefPicList, iRefIdx)->unscaledPic, mv, pcYuvPred, bi, pu.cu->slice->clpRng(compID), eRefPicList, genChromaMv, pu.cu->slice->getScalingRatio(eRefPicList, iRefIdx));
#else
      xPredAffineBlk(compID, pu, pu.cu->slice->getRefPic(eRefPicList, iRefIdx)->unscaledPic, mv, pcYuvPred, bi, pu.cu->slice->clpRng(compID), genChromaMv, pu.cu->slice->getScalingRatio(eRefPicList, iRefIdx));
#endif
    }
    else
    {
      if (isIBC)
      {
        xPredInterBlk(compID, pu, pu.cu->slice->getPic(), mv[0], pcYuvPred, bi, pu.cu->slice->clpRng(compID),
                      bioApplied, isIBC);
      }
      else
      {
        xPredInterBlk( compID, pu, pu.cu->slice->getRefPic( eRefPicList, iRefIdx )->unscaledPic, mv[0], pcYuvPred, bi, pu.cu->slice->clpRng( compID ), bioApplied, isIBC, pu.cu->slice->getScalingRatio( eRefPicList, iRefIdx ) );
      }
    }
  }
}

#if MULTI_PASS_DMVR
void InterPrediction::xPredInterBiSubPuBDOF(PredictionUnit &pu, PelUnitBuf &pcYuvPred, const bool luma, const bool chroma)
{
  const Slice &slice = *pu.cs->slice;
  bool bioApplied = true;
  // common variable for all subPu
  const bool lumaOnly = (luma && !chroma), chromaOnly = (!luma && chroma);
  const int bioDy = std::min<int>(pu.lumaSize().height, BDOF_SUBPU_DIM);
  const int bioDx = std::min<int>(pu.lumaSize().width,  BDOF_SUBPU_DIM);
  const int scaleX = getComponentScaleX(COMPONENT_Cb, pu.chromaFormat);
  const int scaleY = getComponentScaleY(COMPONENT_Cb, pu.chromaFormat);
  CPelUnitBuf srcPred0 = ( pu.chromaFormat == CHROMA_400 ?
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y())) :
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[0][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[0][2], pcYuvPred.Cr())) );
  CPelUnitBuf srcPred1 = ( pu.chromaFormat == CHROMA_400 ?
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y())) :
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[1][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[1][2], pcYuvPred.Cr())) );
  Position puPos = pu.lumaPos();
  PredictionUnit subPu = pu;
  CHECK(subPu.refIdx[0] < 0, "this is not possible for BDOF");
  CHECK(subPu.refIdx[1] < 0, "this is not possible for BDOF");
  int bioSubPuIdx = 0;
  const int bioSubPuStrideIncr = BDOF_SUBPU_STRIDE - std::max(1, (int)(pu.lumaSize().width >> BDOF_SUBPU_DIM_LOG2));
  for (int y = puPos.y, yStart = 0; y < (puPos.y + pu.lumaSize().height); y = y + bioDy, yStart = yStart + bioDy)
  {
    for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + bioDx, xStart = xStart + bioDx)
    {
      Mv bioMv = m_bdofSubPuMvOffset[bioSubPuIdx];

      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, bioDx, bioDy)));
      if (pu.bdmvrRefine)
      {
        const int bdmvrSubPuIdx = (yStart >> DMVR_SUBCU_HEIGHT_LOG2) * DMVR_SUBPU_STRIDE + (xStart >> DMVR_SUBCU_WIDTH_LOG2);
        subPu.mv[0] = m_bdmvrSubPuMvBuf[0][bdmvrSubPuIdx] + bioMv;
        subPu.mv[1] = m_bdmvrSubPuMvBuf[1][bdmvrSubPuIdx] - bioMv;
      }
      else
      {
        subPu.mv[0] = pu.mv[0] + bioMv;
        subPu.mv[1] = pu.mv[1] - bioMv;
      }

      // inter pred to generate buf data
      for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
      {
        if( subPu.refIdx[refList] < 0)
        {
          continue;
        }
        RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

        CHECK(CU::isIBC(*subPu.cu) && eRefPicList != REF_PIC_LIST_0, "Invalid interdir for ibc mode");
        CHECK(CU::isIBC(*subPu.cu) && subPu.refIdx[refList] != MAX_NUM_REF, "Invalid reference index for ibc mode");
        CHECK((CU::isInter(*subPu.cu) && subPu.refIdx[refList] >= slice.getNumRefIdx(eRefPicList)), "Invalid reference index");
        m_iRefListIdx = refList;

        PelUnitBuf pcMbBuf = ( subPu.chromaFormat == CHROMA_400 ?
                               PelUnitBuf(subPu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y())) :
                               PelUnitBuf(subPu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[refList][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[refList][2], pcYuvPred.Cr())) );
        pcMbBuf = pcMbBuf.subBuf(UnitAreaRelative(pu, subPu));
        if (bioMv.hor == 0 && bioMv.ver == 0)
        {
          // only chroma MC
          if (!lumaOnly)
            xPredInterUni ( subPu, eRefPicList, pcMbBuf, true, bioApplied, false, chroma, false );
        }
        else
        {
          xPredInterUni ( subPu, eRefPicList, pcMbBuf, true, bioApplied, luma, chroma, false );
        }
      }
      // prepare dst sub buf
      PelUnitBuf subYuvPredBuf = pcYuvPred.subBuf(UnitAreaRelative(pu, subPu));
      int dstStride[MAX_NUM_COMPONENT] = { pcYuvPred.bufs[COMPONENT_Y].stride,
                                           isChromaEnabled(pu.chromaFormat) ? pcYuvPred.bufs[COMPONENT_Cb].stride : 0,
                                           isChromaEnabled(pu.chromaFormat) ? pcYuvPred.bufs[COMPONENT_Cr].stride : 0};
      subYuvPredBuf.bufs[COMPONENT_Y].buf = pcYuvPred.bufs[COMPONENT_Y].buf + xStart + yStart * dstStride[COMPONENT_Y];
      if (isChromaEnabled(pu.chromaFormat))
      {
        subYuvPredBuf.bufs[COMPONENT_Cb].buf = pcYuvPred.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cb]);
        subYuvPredBuf.bufs[COMPONENT_Cr].buf = pcYuvPred.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cr]);
      }
      // prepare src sub buf
      int srcStride[MAX_NUM_COMPONENT] = { srcPred0.bufs[COMPONENT_Y].stride,
                                           isChromaEnabled(pu.chromaFormat) ? srcPred0.bufs[COMPONENT_Cb].stride : 0,
                                           isChromaEnabled(pu.chromaFormat) ? srcPred0.bufs[COMPONENT_Cr].stride : 0};
      CPelUnitBuf srcSubPred0 = srcPred0.subBuf(UnitAreaRelative(pu, subPu));
      CPelUnitBuf srcSubPred1 = srcPred1.subBuf(UnitAreaRelative(pu, subPu));
      srcSubPred0.bufs[COMPONENT_Y].buf = srcPred0.bufs[COMPONENT_Y].buf + xStart + yStart * srcStride[COMPONENT_Y];
      if (isChromaEnabled(pu.chromaFormat))
      {
        srcSubPred0.bufs[COMPONENT_Cb].buf = srcPred0.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cb]);
        srcSubPred0.bufs[COMPONENT_Cr].buf = srcPred0.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cr]);
      }
      srcSubPred1.bufs[COMPONENT_Y].buf = srcPred1.bufs[COMPONENT_Y].buf + xStart + yStart * srcStride[COMPONENT_Y];
      if (isChromaEnabled(pu.chromaFormat))
      {
        srcSubPred1.bufs[COMPONENT_Cb].buf = srcPred1.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cb]);
        srcSubPred1.bufs[COMPONENT_Cr].buf = srcPred1.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cr]);
      }
      // generate the dst buf
      {
        if (bioMv.hor == 0 && bioMv.ver == 0)
        {
          // only derive chroma prediction
#if JVET_Z0136_OOB
          if (!lumaOnly)
          {
            bool isOOB[2] = { false,false };
            if (pu.interDir == 3)
            {
              isOOB[0] = isMvOOB(subPu.mv[0], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
              isOOB[1] = isMvOOB(subPu.mv[1], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
              xWeightedAverage(false/*isBdofMvRefine*/, 0/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, true, NULL, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
            }
            else
            {
              xWeightedAverage(false/*isBdofMvRefine*/, 0/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, true, NULL, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
            }
          }
#else
          if (!lumaOnly)
            xWeightedAverage( false/*isBdofMvRefine*/, 0/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), false/*bioApplied*/, lumaOnly, true/*chromaOnly*/, NULL/*yuvPredTmp*/ );
#endif
        }
        else
        {
#if JVET_Z0136_OOB
          bool isOOB[2] = { false,false };
          if (pu.interDir == 3)
          {
            isOOB[0] = isMvOOB(subPu.mv[0], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
            isOOB[1] = isMvOOB(subPu.mv[1], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
            xWeightedAverage(false/*isBdofMvRefine*/, 0/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, NULL, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
          }
          else
          {
            xWeightedAverage(false/*isBdofMvRefine*/, 0/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, NULL, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
          }
#else
          xWeightedAverage( false/*isBdofMvRefine*/, 0/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, NULL/*yuvPredTmp*/ );
#endif
        }
      }
      bioSubPuIdx += 1;
    }
    bioSubPuIdx += bioSubPuStrideIncr;
  }
}
#endif

#if MULTI_PASS_DMVR
void InterPrediction::xPredInterBiBDMVR(PredictionUnit &pu, PelUnitBuf &pcYuvPred, const bool luma, const bool chroma, PelUnitBuf *yuvPredTmp /*= NULL*/)
{
  const PPS   &pps   = *pu.cs->pps;
  const Slice &slice = *pu.cs->slice;
#if !INTER_RM_SIZE_CONSTRAINTS
#if ENABLE_OBMC
  if (pu.cu->isobmcMC == false)
#endif
  CHECK( !pu.cu->affine && pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0 && ( pu.lwidth() + pu.lheight() == 12 ), "invalid 4x8/8x4 bi-predicted blocks" );
#endif

  int refIdx0 = pu.refIdx[REF_PIC_LIST_0];
  int refIdx1 = pu.refIdx[REF_PIC_LIST_1];

  const WPScalingParam *wp0 = pu.cs->slice->getWpScaling(REF_PIC_LIST_0, refIdx0);
  const WPScalingParam *wp1 = pu.cs->slice->getWpScaling(REF_PIC_LIST_1, refIdx1);

  bool bioApplied = false;
  if (pu.cs->sps->getBDOFEnabledFlag() && (!pu.cs->picHeader->getDisBdofFlag()))
  {
#if INTER_LIC
    if (pu.cu->affine || m_subPuMC || pu.cu->LICFlag
#if ENABLE_OBMC
      || pu.cu->isobmcMC
#endif
      )
#else
    if (pu.cu->affine || m_subPuMC)
#endif
    {
      bioApplied = false;
    }
    else
    {
      const bool biocheck0 =
        !((WPScalingParam::isWeighted(wp0) || WPScalingParam::isWeighted(wp1)) && slice.getSliceType() == B_SLICE);
      const bool biocheck1 = !(pps.getUseWP() && slice.getSliceType() == P_SLICE);
      if (biocheck0
        && biocheck1
        && PU::isBiPredFromDifferentDirEqDistPoc(pu)
#if !BDOF_RM_CONSTRAINTS
        && (pu.Y().height >= 8)
        && (pu.Y().width >= 8)
        && ((pu.Y().height * pu.Y().width) >= 128)
#endif
       )
      {
        bioApplied = true;
      }
    }

    if (bioApplied && pu.ciipFlag)
    {
      bioApplied = false;
    }

    if (bioApplied && pu.cu->smvdMode)
    {
      bioApplied = false;
    }

    if (pu.cu->cs->sps->getUseBcw() && bioApplied && pu.cu->BcwIdx != BCW_DEFAULT)
    {
      bioApplied = false;
    }
  }
  if (pu.mmvdEncOptMode == 2 && pu.mmvdMergeFlag)
  {
    bioApplied = false;
  }
#if ENABLE_OBMC
  if (pu.cu->isobmcMC)
  {
    bioApplied = false;
  }
#endif
  const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();
  bool dmvrApplied = false;
  dmvrApplied = (pu.mvRefine) && PU::checkDMVRCondition(pu);
  const bool refIsScaled = isResamplingPossible && ( ( refIdx0 < 0 ? false : pu.cu->slice->getRefPic( REF_PIC_LIST_0, refIdx0 )->isRefScaled( pu.cs->pps ) ) || ( refIdx1 < 0 ? false : pu.cu->slice->getRefPic( REF_PIC_LIST_1, refIdx1 )->isRefScaled( pu.cs->pps ) ) );
  dmvrApplied = dmvrApplied && !refIsScaled;
  bioApplied = bioApplied && !refIsScaled;
  // common variable for all subPu
  const bool lumaOnly = (luma && !chroma), chromaOnly = (!luma && chroma);
  const int dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
  const int dx = std::min<int>(pu.lumaSize().width,  DMVR_SUBCU_WIDTH);
  const int scaleX = getComponentScaleX(COMPONENT_Cb, pu.chromaFormat);
  const int scaleY = getComponentScaleY(COMPONENT_Cb, pu.chromaFormat);
  CPelUnitBuf srcPred0 = ( pu.chromaFormat == CHROMA_400 ?
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y())) :
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[0][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[0][2], pcYuvPred.Cr())) );
  CPelUnitBuf srcPred1 = ( pu.chromaFormat == CHROMA_400 ?
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y())) :
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[1][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[1][2], pcYuvPred.Cr())) );
  Position puPos = pu.lumaPos();
  PredictionUnit subPu = pu;

  int subPuIdx = 0;
  const int dmvrSubPuStrideIncr = DMVR_SUBPU_STRIDE - std::max(1, (int)(pu.lumaSize().width >> DMVR_SUBCU_WIDTH_LOG2));
  int length = 0, later = 0;
  int width = pu.lwidth(), height = pu.lheight();
  int subPuIdxColumn = 0;
  if (height > width)
  {
    for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + dx, xStart = xStart + dx)
    {
      subPuIdx = subPuIdxColumn;
      for (int y = puPos.y, yStart = 0; y < (puPos.y + pu.lumaSize().height); y = y + dy, yStart = yStart + dy)
      {
        subPu.mv[0] = m_bdmvrSubPuMvBuf[0][subPuIdx];
        subPu.mv[1] = m_bdmvrSubPuMvBuf[1][subPuIdx];
        length = dy;
        later = yStart + dy;
        subPuIdx += DMVR_SUBPU_STRIDE;
        while (later < width)
        {
          Mv nextMv[2] = { m_bdmvrSubPuMvBuf[0][subPuIdx] , m_bdmvrSubPuMvBuf[1][subPuIdx] };
          if (nextMv[0] == subPu.mv[0] && nextMv[1] == subPu.mv[1])
          {
            length += dy;
          }
          else
          {
            break;
          }

          later += dy;
          subPuIdx += DMVR_SUBPU_STRIDE;
        }
        subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, length)));

        // inter pred to generate buf data
        for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
        {
          if (subPu.refIdx[refList] < 0)
          {
            continue;
          }

          RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

          CHECK(CU::isIBC(*subPu.cu) && eRefPicList != REF_PIC_LIST_0, "Invalid interdir for ibc mode");
          CHECK(CU::isIBC(*subPu.cu) && subPu.refIdx[refList] != MAX_NUM_REF, "Invalid reference index for ibc mode");
          CHECK((CU::isInter(*subPu.cu) && subPu.refIdx[refList] >= slice.getNumRefIdx(eRefPicList)), "Invalid reference index");
          m_iRefListIdx = refList;

          PelUnitBuf pcMbBuf = (subPu.chromaFormat == CHROMA_400 ?
            PelUnitBuf(subPu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y())) :
            PelUnitBuf(subPu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[refList][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[refList][2], pcYuvPred.Cr())));
          pcMbBuf = pcMbBuf.subBuf(UnitAreaRelative(pu, subPu));
          if (subPu.refIdx[0] >= 0 && subPu.refIdx[1] >= 0)
          {
            bool isBdofMvRefineSkipChromaMC = (yuvPredTmp == NULL);
            xPredInterUni(subPu, eRefPicList, pcMbBuf, true
              , bioApplied, luma, chroma, isBdofMvRefineSkipChromaMC);
          }
          else
          {
            if (((pps.getUseWP() && slice.getSliceType() == P_SLICE) || (pps.getWPBiPred() && slice.getSliceType() == B_SLICE))
#if INTER_LIC
              && !subPu.cu->LICFlag
#endif
              )
            {
              xPredInterUni(subPu, eRefPicList, pcMbBuf, true
                , bioApplied
                , luma, chroma
              );
            }
            else
            {
              xPredInterUni(subPu, eRefPicList, pcMbBuf, subPu.cu->geoFlag
                , bioApplied
                , luma, chroma
              );
            }
          }
        }
        // prepare dst sub buf
        PelUnitBuf subYuvPredBuf = pcYuvPred.subBuf(UnitAreaRelative(pu, subPu));
        int dstStride[MAX_NUM_COMPONENT] = { pcYuvPred.bufs[COMPONENT_Y].stride,
                                             isChromaEnabled(pu.chromaFormat) ? pcYuvPred.bufs[COMPONENT_Cb].stride : 0,
                                             isChromaEnabled(pu.chromaFormat) ? pcYuvPred.bufs[COMPONENT_Cr].stride : 0 };
        subYuvPredBuf.bufs[COMPONENT_Y].buf = pcYuvPred.bufs[COMPONENT_Y].buf + xStart + yStart * dstStride[COMPONENT_Y];
        if (isChromaEnabled(pu.chromaFormat))
        {
          subYuvPredBuf.bufs[COMPONENT_Cb].buf = pcYuvPred.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cb]);
          subYuvPredBuf.bufs[COMPONENT_Cr].buf = pcYuvPred.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cr]);
        }
        // prepare src sub buf
        int srcStride[MAX_NUM_COMPONENT] = { srcPred0.bufs[COMPONENT_Y].stride,
                                             isChromaEnabled(pu.chromaFormat) ? srcPred0.bufs[COMPONENT_Cb].stride : 0,
                                             isChromaEnabled(pu.chromaFormat) ? srcPred0.bufs[COMPONENT_Cr].stride : 0 };
        CPelUnitBuf srcSubPred0 = srcPred0.subBuf(UnitAreaRelative(pu, subPu));
        CPelUnitBuf srcSubPred1 = srcPred1.subBuf(UnitAreaRelative(pu, subPu));
        srcSubPred0.bufs[COMPONENT_Y].buf = srcPred0.bufs[COMPONENT_Y].buf + xStart + yStart * srcStride[COMPONENT_Y];
        if (isChromaEnabled(pu.chromaFormat))
        {
          srcSubPred0.bufs[COMPONENT_Cb].buf = srcPred0.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cb]);
          srcSubPred0.bufs[COMPONENT_Cr].buf = srcPred0.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cr]);
        }
        srcSubPred1.bufs[COMPONENT_Y].buf = srcPred1.bufs[COMPONENT_Y].buf + xStart + yStart * srcStride[COMPONENT_Y];
        if (isChromaEnabled(pu.chromaFormat))
        {
          srcSubPred1.bufs[COMPONENT_Cb].buf = srcPred1.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cb]);
          srcSubPred1.bufs[COMPONENT_Cr].buf = srcPred1.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cr]);
        }
        // generate the dst buf
        {
          const int bioSubPuOffset = (xStart >> BDOF_SUBPU_DIM_LOG2) + (yStart >> BDOF_SUBPU_DIM_LOG2) * BDOF_SUBPU_STRIDE;
#if JVET_Z0136_OOB
          bool isOOB[2] = { false,false };
          if (pu.interDir == 3)
          {
            isOOB[0] = isMvOOB(subPu.mv[0], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
            isOOB[1] = isMvOOB(subPu.mv[1], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
            xWeightedAverage(true/*isBdofMvRefine*/, bioSubPuOffset/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
          }
          else
          {
            xWeightedAverage(true/*isBdofMvRefine*/, bioSubPuOffset/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
          }
#else
          xWeightedAverage(true/*isBdofMvRefine*/, bioSubPuOffset/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf,
              slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp);
#endif
        }

        yStart = later - dy;
        y = puPos.y + yStart;
      }
      subPuIdxColumn++;
    }
  }
  else
  for (int y = puPos.y, yStart = 0; y < (puPos.y + pu.lumaSize().height); y = y + dy, yStart = yStart + dy)
  {
    for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + dx, xStart = xStart + dx)
    {
      subPu.mv[0] = m_bdmvrSubPuMvBuf[0][subPuIdx];
      subPu.mv[1] = m_bdmvrSubPuMvBuf[1][subPuIdx];
      length = dx;
      later = xStart + dx;
      subPuIdx++;
      while (later < width)
      {
        Mv nextMv[2] = { m_bdmvrSubPuMvBuf[0][subPuIdx] , m_bdmvrSubPuMvBuf[1][subPuIdx] };
        if (nextMv[0] == subPu.mv[0] && nextMv[1] == subPu.mv[1])
        {
          length += dx;
        }
        else
        {
          break;
        }

        later += dx;
        subPuIdx++;
      }
      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, length, dy)));
      // inter pred to generate buf data
      for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
      {
        if( subPu.refIdx[refList] < 0)
        {
          continue;
        }

        RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

        CHECK(CU::isIBC(*subPu.cu) && eRefPicList != REF_PIC_LIST_0, "Invalid interdir for ibc mode");
        CHECK(CU::isIBC(*subPu.cu) && subPu.refIdx[refList] != MAX_NUM_REF, "Invalid reference index for ibc mode");
        CHECK((CU::isInter(*subPu.cu) && subPu.refIdx[refList] >= slice.getNumRefIdx(eRefPicList)), "Invalid reference index");
        m_iRefListIdx = refList;

        PelUnitBuf pcMbBuf = ( subPu.chromaFormat == CHROMA_400 ?
                               PelUnitBuf(subPu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y())) :
                               PelUnitBuf(subPu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[refList][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[refList][2], pcYuvPred.Cr())) );
        pcMbBuf = pcMbBuf.subBuf(UnitAreaRelative(pu, subPu));
        if (subPu.refIdx[0] >= 0 && subPu.refIdx[1] >= 0)
        {
          bool isBdofMvRefineSkipChromaMC = (yuvPredTmp == NULL);
          xPredInterUni ( subPu, eRefPicList, pcMbBuf, true
            , bioApplied, luma, chroma, isBdofMvRefineSkipChromaMC);
        }
        else
        {
          if( ( (pps.getUseWP() && slice.getSliceType() == P_SLICE) || (pps.getWPBiPred() && slice.getSliceType() == B_SLICE) )
    #if INTER_LIC
            && !subPu.cu->LICFlag
    #endif
            )
          {
            xPredInterUni ( subPu, eRefPicList, pcMbBuf, true
              , bioApplied
              , luma, chroma
            );
          }
          else
          {
            xPredInterUni(subPu, eRefPicList, pcMbBuf, subPu.cu->geoFlag
              , bioApplied
              , luma, chroma
            );
          }
        }
      }
      // prepare dst sub buf
      PelUnitBuf subYuvPredBuf = pcYuvPred.subBuf(UnitAreaRelative(pu, subPu));
      int dstStride[MAX_NUM_COMPONENT] = { pcYuvPred.bufs[COMPONENT_Y].stride,
                                           isChromaEnabled(pu.chromaFormat) ? pcYuvPred.bufs[COMPONENT_Cb].stride : 0,
                                           isChromaEnabled(pu.chromaFormat) ? pcYuvPred.bufs[COMPONENT_Cr].stride : 0};
      subYuvPredBuf.bufs[COMPONENT_Y].buf = pcYuvPred.bufs[COMPONENT_Y].buf + xStart + yStart * dstStride[COMPONENT_Y];
      if (isChromaEnabled(pu.chromaFormat))
      {
        subYuvPredBuf.bufs[COMPONENT_Cb].buf = pcYuvPred.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cb]);
        subYuvPredBuf.bufs[COMPONENT_Cr].buf = pcYuvPred.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cr]);
      }
      // prepare src sub buf
      int srcStride[MAX_NUM_COMPONENT] = { srcPred0.bufs[COMPONENT_Y].stride,
                                           isChromaEnabled(pu.chromaFormat) ? srcPred0.bufs[COMPONENT_Cb].stride : 0,
                                           isChromaEnabled(pu.chromaFormat) ? srcPred0.bufs[COMPONENT_Cr].stride : 0};
      CPelUnitBuf srcSubPred0 = srcPred0.subBuf(UnitAreaRelative(pu, subPu));
      CPelUnitBuf srcSubPred1 = srcPred1.subBuf(UnitAreaRelative(pu, subPu));
      srcSubPred0.bufs[COMPONENT_Y].buf = srcPred0.bufs[COMPONENT_Y].buf + xStart + yStart * srcStride[COMPONENT_Y];
      if (isChromaEnabled(pu.chromaFormat))
      {
        srcSubPred0.bufs[COMPONENT_Cb].buf = srcPred0.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cb]);
        srcSubPred0.bufs[COMPONENT_Cr].buf = srcPred0.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cr]);
      }
      srcSubPred1.bufs[COMPONENT_Y].buf = srcPred1.bufs[COMPONENT_Y].buf + xStart + yStart * srcStride[COMPONENT_Y];
      if (isChromaEnabled(pu.chromaFormat))
      {
        srcSubPred1.bufs[COMPONENT_Cb].buf = srcPred1.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cb]);
        srcSubPred1.bufs[COMPONENT_Cr].buf = srcPred1.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * srcStride[COMPONENT_Cr]);
      }
      // generate the dst buf
      {
        const int bioSubPuOffset = (xStart >> BDOF_SUBPU_DIM_LOG2) + (yStart >> BDOF_SUBPU_DIM_LOG2) * BDOF_SUBPU_STRIDE;
#if JVET_Z0136_OOB
        bool isOOB[2] = { false,false };
        if (pu.interDir == 3)
        {
          isOOB[0] = isMvOOB(subPu.mv[0], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
          isOOB[1] = isMvOOB(subPu.mv[1], subPu.Y().topLeft(), subPu.lumaSize(), subPu.cu->slice->getSPS(), subPu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
          xWeightedAverage(true/*isBdofMvRefine*/, bioSubPuOffset/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
        }
        else
        {
          xWeightedAverage(true/*isBdofMvRefine*/, bioSubPuOffset/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, subYuvPredBuf.Y().width, pu.cs->mcMaskChroma, subYuvPredBuf.Cb().width, isOOB);
        }
#else
        xWeightedAverage( true/*isBdofMvRefine*/, bioSubPuOffset/*bdofBlockOffset*/, subPu, srcSubPred0, srcSubPred1, subYuvPredBuf,
            slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp );
#endif
      }

      xStart = later - dx;
      x = puPos.x + xStart;
    }
    subPuIdx += dmvrSubPuStrideIncr;
  }
}
#endif

void InterPrediction::xPredInterBi(PredictionUnit &pu, PelUnitBuf &pcYuvPred, const bool luma, const bool chroma, PelUnitBuf *yuvPredTmp /*= NULL*/)
{
  const PPS   &pps = *pu.cs->pps;
  const Slice &slice = *pu.cs->slice;
#if MULTI_PASS_DMVR
  if ( pu.bdmvrRefine )
  {
    if (yuvPredTmp && (pu.lwidth() > DMVR_SUBCU_WIDTH || pu.lheight() > DMVR_SUBCU_HEIGHT)) // pre-do MC for yuvPredTmp to avoid MC for yuvPredTmp within the subblock loop
    {
      for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
      {
        CHECK(pu.refIdx[refList] == NOT_VALID, "pu.refIdx[refList] shouldn't be NOT_VALID.")
          RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
        m_iRefListIdx = refList;

        PelUnitBuf pcMbBuf = (pu.chromaFormat == CHROMA_400 ?
          PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y())) :
          PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[refList][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[refList][2], pcYuvPred.Cr())));

        xPredInterUni(pu, eRefPicList, pcMbBuf, true, false, luma, chroma);
      }
      CPelUnitBuf srcPred0 = (pu.chromaFormat == CHROMA_400 ?
        CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y())) :
        CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[0][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[0][2], pcYuvPred.Cr())));
      CPelUnitBuf srcPred1 = (pu.chromaFormat == CHROMA_400 ?
        CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y())) :
        CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[1][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[1][2], pcYuvPred.Cr())));
      const bool lumaOnly = luma && !chroma;
      const bool chromaOnly = !luma && chroma;
      if (pps.getWPBiPred() && slice.getSliceType() == B_SLICE && pu.cu->BcwIdx == BCW_DEFAULT)
      {
        xWeightedPredictionBi(pu, srcPred0, srcPred1, *yuvPredTmp, m_maxCompIDToPred, lumaOnly, chromaOnly);
      }
      else if (pps.getUseWP() && slice.getSliceType() == P_SLICE)
      {
        xWeightedPredictionUni(pu, srcPred0, REF_PIC_LIST_0, *yuvPredTmp, -1, m_maxCompIDToPred, lumaOnly, chromaOnly);
      }
      else
      {
#if JVET_Z0136_OOB
        bool isOOB[2] = { false,false };
        if (pu.interDir == 3)
        {
          isOOB[0] = isMvOOB(pu.mv[0], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
          isOOB[1] = isMvOOB(pu.mv[1], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
          xWeightedAverage(false/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, *yuvPredTmp, slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, chromaOnly, NULL, pu.cs->mcMask, yuvPredTmp->Y().width, pu.cs->mcMaskChroma, yuvPredTmp->Cb().width, isOOB);
        }
        else
        {
          xWeightedAverage(false/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, *yuvPredTmp, slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, chromaOnly, NULL, pu.cs->mcMask, yuvPredTmp->Y().width, pu.cs->mcMaskChroma, yuvPredTmp->Cb().width, isOOB);
        }
#else
        xWeightedAverage(false/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, *yuvPredTmp,
            slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, chromaOnly);
#endif
      }

      yuvPredTmp = nullptr;
    }
    xPredInterBiBDMVR(pu, pcYuvPred, luma, chroma, yuvPredTmp);
    return;
  }
#endif
#if !INTER_RM_SIZE_CONSTRAINTS
#if ENABLE_OBMC
  if (pu.cu->isobmcMC == false)
#endif
  CHECK( !pu.cu->affine && pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0 && ( pu.lwidth() + pu.lheight() == 12 ), "invalid 4x8/8x4 bi-predicted blocks" );
#endif
  int refIdx0 = pu.refIdx[REF_PIC_LIST_0];
  int refIdx1 = pu.refIdx[REF_PIC_LIST_1];

  const WPScalingParam *wp0 = pu.cs->slice->getWpScaling(REF_PIC_LIST_0, refIdx0);
  const WPScalingParam *wp1 = pu.cs->slice->getWpScaling(REF_PIC_LIST_1, refIdx1);

  bool bioApplied = false;
  if (pu.cs->sps->getBDOFEnabledFlag() && (!pu.cs->picHeader->getDisBdofFlag()))
  {
#if INTER_LIC
    if (pu.cu->affine || m_subPuMC || pu.cu->LICFlag)
#else
    if (pu.cu->affine || m_subPuMC)
#endif
    {
      bioApplied = false;
    }
    else
    {
      const bool biocheck0 =
        !((WPScalingParam::isWeighted(wp0) || WPScalingParam::isWeighted(wp1)) && slice.getSliceType() == B_SLICE);
      const bool biocheck1 = !(pps.getUseWP() && slice.getSliceType() == P_SLICE);      if (biocheck0
        && biocheck1
        && PU::isBiPredFromDifferentDirEqDistPoc(pu)
#if !BDOF_RM_CONSTRAINTS
        && (pu.Y().height >= 8)
        && (pu.Y().width >= 8)
        && ((pu.Y().height * pu.Y().width) >= 128)
#endif
       )
      {
        bioApplied = true;
      }
    }

    if (bioApplied && pu.ciipFlag)
    {
      bioApplied = false;
    }

    if (bioApplied && pu.cu->smvdMode)
    {
      bioApplied = false;
    }

    if (pu.cu->cs->sps->getUseBcw() && bioApplied && pu.cu->BcwIdx != BCW_DEFAULT)
    {
      bioApplied = false;
    }
  }
  if (pu.mmvdEncOptMode == 2 && pu.mmvdMergeFlag)
  {
    bioApplied = false;
  }
#if ENABLE_OBMC
  if (pu.cu->isobmcMC)
  {
    bioApplied = false;
  }
#endif
  bool dmvrApplied = false;
  dmvrApplied = (pu.mvRefine) && PU::checkDMVRCondition(pu);

  const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();
  const bool refIsScaled = isResamplingPossible && ( ( refIdx0 < 0 ? false : pu.cu->slice->getRefPic( REF_PIC_LIST_0, refIdx0 )->isRefScaled( pu.cs->pps ) ) || ( refIdx1 < 0 ? false : pu.cu->slice->getRefPic( REF_PIC_LIST_1, refIdx1 )->isRefScaled( pu.cs->pps ) ) );
  dmvrApplied = dmvrApplied && !refIsScaled;
  bioApplied = bioApplied && !refIsScaled;

#if MULTI_PASS_DMVR
  if (yuvPredTmp && bioApplied && (pu.lwidth() > BDOF_SUBPU_DIM || pu.lheight() > BDOF_SUBPU_DIM)) // pre-do MC for yuvPredTmp to avoid MC for yuvPredTmp within the subblock loop
  {
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      CHECK(pu.refIdx[refList] == NOT_VALID, "pu.refIdx[refList] shouldn't be NOT_VALID.")
        RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
      m_iRefListIdx = refList;

      PelUnitBuf pcMbBuf = (pu.chromaFormat == CHROMA_400 ?
        PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y())) :
        PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[refList][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[refList][2], pcYuvPred.Cr())));

      xPredInterUni(pu, eRefPicList, pcMbBuf, true, false, luma, chroma);
    }
    CPelUnitBuf srcPred0 = (pu.chromaFormat == CHROMA_400 ?
      CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y())) :
      CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[0][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[0][2], pcYuvPred.Cr())));
    CPelUnitBuf srcPred1 = (pu.chromaFormat == CHROMA_400 ?
      CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y())) :
      CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[1][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[1][2], pcYuvPred.Cr())));
    const bool lumaOnly = luma && !chroma;
    const bool chromaOnly = !luma && chroma;
    if (pps.getWPBiPred() && slice.getSliceType() == B_SLICE && pu.cu->BcwIdx == BCW_DEFAULT)
    {
      xWeightedPredictionBi(pu, srcPred0, srcPred1, *yuvPredTmp, m_maxCompIDToPred, lumaOnly, chromaOnly);
    }
    else if (pps.getUseWP() && slice.getSliceType() == P_SLICE)
    {
      xWeightedPredictionUni(pu, srcPred0, REF_PIC_LIST_0, *yuvPredTmp, -1, m_maxCompIDToPred, lumaOnly, chromaOnly);
    }
    else
    {
#if JVET_Z0136_OOB
      bool isOOB[2] = { false,false };
      if (pu.interDir == 3)
      {
        isOOB[0] = isMvOOB(pu.mv[0], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
        isOOB[1] = isMvOOB(pu.mv[1], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
        xWeightedAverage(false/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, *yuvPredTmp, slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, chromaOnly, NULL, pu.cs->mcMask, yuvPredTmp->Y().width, pu.cs->mcMaskChroma, yuvPredTmp->Cb().width, isOOB);
      }
      else
      {
        xWeightedAverage(false/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, *yuvPredTmp, slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, chromaOnly, NULL, pu.cs->mcMask, yuvPredTmp->Y().width, pu.cs->mcMaskChroma, yuvPredTmp->Cb().width, isOOB);
      }
#else
      xWeightedAverage(false, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, *yuvPredTmp, slice.getSPS()->getBitDepths(), slice.clpRngs(), false, lumaOnly, chromaOnly);
#endif
    }

    yuvPredTmp = nullptr;
  }
#endif

  for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
  {
    if( pu.refIdx[refList] < 0)
    {
      continue;
    }

    RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

    CHECK(CU::isIBC(*pu.cu) && eRefPicList != REF_PIC_LIST_0, "Invalid interdir for ibc mode");
    CHECK(CU::isIBC(*pu.cu) && pu.refIdx[refList] != MAX_NUM_REF, "Invalid reference index for ibc mode");
    CHECK((CU::isInter(*pu.cu) && pu.refIdx[refList] >= slice.getNumRefIdx(eRefPicList)), "Invalid reference index");
    m_iRefListIdx = refList;

    PelUnitBuf pcMbBuf = ( pu.chromaFormat == CHROMA_400 ?
                           PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y())) :
                           PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[refList][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[refList][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[refList][2], pcYuvPred.Cr())) );

    if (pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0)
    {
      if (dmvrApplied)
      {
        if (yuvPredTmp)
        {
          xPredInterUni(pu, eRefPicList, pcMbBuf, true, false, luma, chroma);
        }
        continue;
      }
#if MULTI_PASS_DMVR
      bool isBdofMvRefineSkipChromaMC = (yuvPredTmp == NULL);
      xPredInterUni(pu, eRefPicList, pcMbBuf, true, bioApplied, luma, chroma, isBdofMvRefineSkipChromaMC);
#else
      xPredInterUni(pu, eRefPicList, pcMbBuf, true, bioApplied, luma, chroma);
#endif
    }
    else
    {
      if( ( (pps.getUseWP() && slice.getSliceType() == P_SLICE) || (pps.getWPBiPred() && slice.getSliceType() == B_SLICE) )
#if INTER_LIC
        && !pu.cu->LICFlag
#endif
        )
      {
        xPredInterUni(pu, eRefPicList, pcMbBuf, true, bioApplied, luma, chroma);
      }
      else
      {
        xPredInterUni(pu, eRefPicList, pcMbBuf, pu.cu->geoFlag, bioApplied, luma, chroma);
      }
    }
  }
  CPelUnitBuf srcPred0 = ( pu.chromaFormat == CHROMA_400 ?
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y())) :
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[0][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[0][2], pcYuvPred.Cr())) );
  CPelUnitBuf srcPred1 = ( pu.chromaFormat == CHROMA_400 ?
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y())) :
                           CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvPred.Y()), PelBuf(m_acYuvPred[1][1], pcYuvPred.Cb()), PelBuf(m_acYuvPred[1][2], pcYuvPred.Cr())) );
  const bool lumaOnly   = luma && !chroma;
  const bool chromaOnly = !luma && chroma;
  if( !pu.cu->geoFlag && (!dmvrApplied) && (!bioApplied) && pps.getWPBiPred() && slice.getSliceType() == B_SLICE && pu.cu->BcwIdx == BCW_DEFAULT)
  {
    xWeightedPredictionBi( pu, srcPred0, srcPred1, pcYuvPred, m_maxCompIDToPred, lumaOnly, chromaOnly );
    if (yuvPredTmp)
    {
      yuvPredTmp->copyFrom(pcYuvPred);
    }
  }
  else if( !pu.cu->geoFlag && pps.getUseWP() && slice.getSliceType() == P_SLICE )
  {
    xWeightedPredictionUni( pu, srcPred0, REF_PIC_LIST_0, pcYuvPred, -1, m_maxCompIDToPred, lumaOnly, chromaOnly );
    if (yuvPredTmp)
    {
      yuvPredTmp->copyFrom(pcYuvPred);
    }
  }
  else
  {
    if (dmvrApplied)
    {
      if (yuvPredTmp)
      {
        yuvPredTmp->addAvg(srcPred0, srcPred1, slice.clpRngs(), false);
      }
      xProcessDMVR(pu, pcYuvPred, slice.clpRngs(), bioApplied);
    }
    else
    {
#if MULTI_PASS_DMVR
#if JVET_Z0136_OOB
      bool isOOB[2] = { false,false };
      if (pu.interDir == 3)
      {
        if (pu.cu->affine && pu.mergeType != MRG_TYPE_SUBPU_ATMVP)  // affine
        {
          bool *pMcMask0 = pu.cs->mcMask[0];
          bool *pMcMask1 = pu.cs->mcMask[1];
          for (int h = 0; h < (int)pu.lumaSize().height && (!isOOB[0] || !isOOB[1]); h++)
          {
            for (int w = 0; w < (int)pu.lumaSize().width && (!isOOB[0] || !isOOB[1]); w++)
            {
              isOOB[0] |= pMcMask0[w];
              isOOB[1] |= pMcMask1[w];
            }
            pMcMask0 += (int)pu.lumaSize().width; pMcMask1 += (int)pu.lumaSize().width;
          }
        }
        else
        {
          isOOB[0] = isMvOOB(pu.mv[0], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
          isOOB[1] = isMvOOB(pu.mv[1], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
        }
        xWeightedAverage(true/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, pcYuvPred, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, pcYuvPred.Y().width, pu.cs->mcMaskChroma, pcYuvPred.Cb().width, isOOB);
      }
      else
      {
        xWeightedAverage(true/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, pcYuvPred, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, pcYuvPred.Y().width, pu.cs->mcMaskChroma, pcYuvPred.Cb().width, isOOB);
      }
#else
      xWeightedAverage( true/*isBdofMvRefine*/, 0/*bioSubPuOffset*/, pu, srcPred0, srcPred1, pcYuvPred,
          slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp );
#endif
#else
#if JVET_Z0136_OOB
      bool isOOB[2] = { false,false };
      if (pu.interDir == 3)
      {
        if (pu.cu->affine && pu.mergeType != MRG_TYPE_SUBPU_ATMVP)  // affine
        {
          bool *pMcMask0 = pu.cs->mcMask[0];
          bool *pMcMask1 = pu.cs->mcMask[1];
          for (int h = 0; h < (int)pu.lumaSize().height && (!isOOB[0] || !isOOB[1]); h++)
          {
            for (int w = 0; w < (int)pu.lumaSize().width && (!isOOB[0] || !isOOB[1]); w++)
            {
              isOOB[0] |= pMcMask0[w];
              isOOB[1] |= pMcMask1[w];
            }
            pMcMask0 += (int)pu.lumaSize().width; pMcMask1 += (int)pu.lumaSize().width;
          }
        }
        else
        {
          isOOB[0] = isMvOOB(pu.mv[0], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[0], pu.cs->mcMaskChroma[0]);
          isOOB[1] = isMvOOB(pu.mv[1], pu.Y().topLeft(), pu.lumaSize(), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pu.cs->mcMask[1], pu.cs->mcMaskChroma[1]);
        }
        xWeightedAverage(pu, srcPred0, srcPred1, pcYuvPred, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, pcYuvPred.Y().width, pu.cs->mcMaskChroma, pcYuvPred.Cb().width, isOOB);
      }
      else
      {
        xWeightedAverage(pu, srcPred0, srcPred1, pcYuvPred, slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp, pu.cs->mcMask, pcYuvPred.Y().width, pu.cs->mcMaskChroma, pcYuvPred.Cb().width, isOOB);
      }
#else
      xWeightedAverage( pu, srcPred0, srcPred1, pcYuvPred,
          slice.getSPS()->getBitDepths(), slice.clpRngs(), bioApplied, lumaOnly, chromaOnly, yuvPredTmp );
#endif
#endif
    }
  }
}

void InterPrediction::xPredInterBlk ( const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv& _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng
                                     , const bool& bioApplied
                                     , bool isIBC
                                     , const std::pair<int, int> scalingRatio
                                     , SizeType dmvrWidth
                                     , SizeType dmvrHeight
                                     , bool bilinearMC
                                     , Pel *srcPadBuf
                                     , int32_t srcPadStride
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                                     , bool isAML
#if INTER_LIC
                                     , bool doLic
                                     , Mv   mvCurr
#endif
#endif
#if JVET_Z0061_TM_OBMC
                                    , bool fastOBMC
#endif
                                    )
{
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
#if  JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  int filterIdx = isAML && pu.mmvdMergeFlag ? 1 : 0;
#else
  int filterIdx = 0;
#endif
  if (bilinearMC)
  {
    filterIdx = 1;
  }
#if JVET_Z0061_TM_OBMC
  if (fastOBMC)
  {
    filterIdx = 1;
  }
#endif
#endif

  JVET_J0090_SET_REF_PICTURE( refPic, compID );
  const ChromaFormat  chFmt = pu.chromaFormat;
  const bool          rndRes = !bi;

  int shiftHor = MV_FRACTIONAL_BITS_INTERNAL + ::getComponentScaleX(compID, chFmt);
  int shiftVer = MV_FRACTIONAL_BITS_INTERNAL + ::getComponentScaleY(compID, chFmt);

  bool  wrapRef = false;
  Mv    mv(_mv);
  if( !isIBC && refPic->isWrapAroundEnabled( pu.cs->pps ) )
  {
    wrapRef = wrapClipMv( mv, pu.blocks[0].pos(), pu.blocks[0].size(), pu.cs->sps, pu.cs->pps );
  }

  bool useAltHpelIf = pu.cu->imv == IMV_HPEL;

  const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();

  if( isResamplingPossible && !isIBC && xPredInterBlkRPR( scalingRatio, *pu.cs->pps, CompArea( compID, chFmt, pu.blocks[compID], Size( dstPic.bufs[compID].width, dstPic.bufs[compID].height ) ), refPic, mv, dstPic.bufs[compID].buf, dstPic.bufs[compID].stride, bi, wrapRef, clpRng, 0, useAltHpelIf ) )
  {
    CHECK( bilinearMC, "DMVR should be disabled with RPR" );
    CHECK( bioApplied, "BDOF should be disabled with RPR" );
  }
  else
  {
    int xFrac = mv.hor & ((1 << shiftHor) - 1);
    int yFrac = mv.ver & ((1 << shiftVer) - 1);
    if (isIBC)
    {
      xFrac = yFrac = 0;
      JVET_J0090_SET_CACHE_ENABLE(false);
    }

    PelBuf & dstBuf = dstPic.bufs[compID];
    unsigned width  = dstBuf.width;
    unsigned height = dstBuf.height;

    CPelBuf refBuf;
    {
      Position offset = pu.blocks[compID].pos().offset(mv.getHor() >> shiftHor, mv.getVer() >> shiftVer);
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
      int refBufExtendSize = 0;
      if (bioApplied && compID == COMPONENT_Y)
      {
        refBufExtendSize = ((BIO_EXTEND_SIZE + 1) << 1);  // trick to use SIMD filter
        offset.x -= (BIO_EXTEND_SIZE + 1);
        offset.y -= (BIO_EXTEND_SIZE + 1);
      }
      if (dmvrWidth)
      {
        refBuf = refPic->getRecoBuf(CompArea(compID, chFmt, offset, Size(dmvrWidth + refBufExtendSize, dmvrHeight + refBufExtendSize)), wrapRef);
      }
      else
      {
        refBuf = refPic->getRecoBuf( CompArea( compID, chFmt, offset, Size(pu.blocks[compID].width + refBufExtendSize, pu.blocks[compID].height + refBufExtendSize) ), wrapRef);
      }
#else
      if (dmvrWidth)
      {
        refBuf = refPic->getRecoBuf(CompArea(compID, chFmt, offset, Size(dmvrWidth, dmvrHeight)), wrapRef);
      }
      else
      {
        refBuf = refPic->getRecoBuf(CompArea(compID, chFmt, offset, pu.blocks[compID].size()), wrapRef);
      }
#endif
    }

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
    if (NULL != srcPadBuf && bioApplied == false)
#else
    if (NULL != srcPadBuf)
#endif
    {
      refBuf.buf    = srcPadBuf;
      refBuf.stride = srcPadStride;
    }
    if (dmvrWidth)
    {
      width  = dmvrWidth;
      height = dmvrHeight;
    }
    // backup data
    int  backupWidth        = width;
    int  backupHeight       = height;
    Pel *backupDstBufPtr    = dstBuf.buf;
    int  backupDstBufStride = dstBuf.stride;

    if (bioApplied && compID == COMPONENT_Y)
    {
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
      backupWidth += ((BIO_EXTEND_SIZE + 1) << 1);
      backupHeight += ((BIO_EXTEND_SIZE + 1) << 1);
      dstBuf.stride = backupWidth;
      dstBuf.buf = m_filteredBlockTmp[2 + m_iRefListIdx][compID];
#else
      width  = width + 2 * BIO_EXTEND_SIZE + 2;
      height = height + 2 * BIO_EXTEND_SIZE + 2;

      // change MC output
      dstBuf.stride = width;
      dstBuf.buf = m_filteredBlockTmp[2 + m_iRefListIdx][compID] + 2 * dstBuf.stride + 2;
#endif
    }

    if( yFrac == 0 )
    {
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
      m_if.filterHor( compID, (Pel*)refBuf.buf, refBuf.stride, dstBuf.buf, dstBuf.stride, backupWidth, backupHeight, xFrac, rndRes, chFmt, clpRng, filterIdx, bilinearMC, useAltHpelIf );
#else
      m_if.filterHor( compID, ( Pel* ) refBuf.buf, refBuf.stride, dstBuf.buf, dstBuf.stride, backupWidth, backupHeight, xFrac, rndRes, chFmt, clpRng, bilinearMC, bilinearMC, useAltHpelIf);
#endif
    }
    else if( xFrac == 0 )
    {
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
      m_if.filterVer( compID, (Pel*)refBuf.buf, refBuf.stride, dstBuf.buf, dstBuf.stride, backupWidth, backupHeight, yFrac, true, rndRes, chFmt, clpRng, filterIdx, bilinearMC, useAltHpelIf );
#else
      m_if.filterVer( compID, ( Pel* ) refBuf.buf, refBuf.stride, dstBuf.buf, dstBuf.stride, backupWidth, backupHeight, yFrac, true, rndRes, chFmt, clpRng, bilinearMC, bilinearMC, useAltHpelIf);
#endif
    }
    else
    {
#if SIMD_4x4_12 && defined(TARGET_SIMD_X86)
      //use 4x4 if possible
      if( compID == COMPONENT_Y
          && backupWidth == 4
          && backupHeight == 4
          && !( (xFrac == 8 || yFrac == 8) && useAltHpelIf ) //to avoid (8,12 or 12,8 passes)
          && dmvrWidth == 0                                  //seems to conflict with DMVR, not sure //kolya
        )
        m_if.filter4x4(clpRng,  (Pel*)refBuf.buf, refBuf.stride,  dstBuf.buf, dstBuf.stride, xFrac, yFrac, rndRes);
      else
      {
#endif
        PelBuf tmpBuf = dmvrWidth ? PelBuf( m_filteredBlockTmp[0][compID], Size( dmvrWidth, dmvrHeight ) ) : PelBuf( m_filteredBlockTmp[0][compID], pu.blocks[compID] );
        if( dmvrWidth == 0 )
        {
          tmpBuf.stride = dstBuf.stride;
        }
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
        if (bioApplied && compID == COMPONENT_Y)
        {
          tmpBuf = PelBuf(m_filteredBlockTmp[0][compID], Size(backupWidth, backupWidth));
          tmpBuf.stride = dstBuf.stride;
        }
#endif
#if IF_12TAP
        int vFilterSize = isLuma( compID ) ? NTAPS_LUMA( 0 ) : NTAPS_CHROMA;
#else
        int vFilterSize = isLuma( compID ) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
        if (isLuma(compID) && filterIdx == 1)
#else
        if (bilinearMC)
#endif
        {
          vFilterSize = NTAPS_BILINEAR;
        }
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
        m_if.filterHor(compID, (Pel*)refBuf.buf - ((vFilterSize >> 1) - 1) * refBuf.stride, refBuf.stride, tmpBuf.buf, tmpBuf.stride, backupWidth, backupHeight + vFilterSize - 1, xFrac, false, chFmt, clpRng, filterIdx, bilinearMC, useAltHpelIf);
        JVET_J0090_SET_CACHE_ENABLE(false);
        m_if.filterVer(compID, (Pel*)tmpBuf.buf + ((vFilterSize >> 1) - 1) * tmpBuf.stride, tmpBuf.stride, dstBuf.buf, dstBuf.stride, backupWidth, backupHeight, yFrac, false, rndRes, chFmt, clpRng, filterIdx, bilinearMC, useAltHpelIf);
#else
        m_if.filterHor( compID, ( Pel* ) refBuf.buf - ( ( vFilterSize >> 1 ) - 1 ) * refBuf.stride, refBuf.stride, tmpBuf.buf, tmpBuf.stride, backupWidth, backupHeight + vFilterSize - 1, xFrac, false, chFmt, clpRng, bilinearMC, bilinearMC, useAltHpelIf);
        JVET_J0090_SET_CACHE_ENABLE( false );
        m_if.filterVer( compID, ( Pel* ) tmpBuf.buf + ( ( vFilterSize >> 1 ) - 1 ) * tmpBuf.stride, tmpBuf.stride, dstBuf.buf, dstBuf.stride, backupWidth, backupHeight, yFrac, false, rndRes, chFmt, clpRng, bilinearMC, bilinearMC, useAltHpelIf);
#endif
#if SIMD_4x4_12 && defined(TARGET_SIMD_X86)
      }
#endif
    }
    JVET_J0090_SET_CACHE_ENABLE(
      ( srcPadStride == 0 )
      && ( bioApplied
           == false ) );   // Enabled only in non-DMVR-non-BDOF process, In DMVR process, srcPadStride is always non-zero

    if( bioApplied && compID == COMPONENT_Y )
    {
#if !MULTI_PASS_DMVR && !SAMPLE_BASED_BDOF
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      const int shift = IF_INTERNAL_FRAC_BITS( clpRng.bd );
#else
      const int shift = std::max<int>( 2, ( IF_INTERNAL_PREC - clpRng.bd ) );
#endif
      int xOffset = ( xFrac < 8 ) ? 1 : 0;
      int yOffset = ( yFrac < 8 ) ? 1 : 0;
      const Pel* refPel = refBuf.buf - yOffset * refBuf.stride - xOffset;
      Pel* dstPel = m_filteredBlockTmp[2 + m_iRefListIdx][compID] + dstBuf.stride + 1;
      for( int w = 0; w < ( width - 2 * BIO_EXTEND_SIZE ); w++ )
      {
        Pel val = leftShift_round( refPel[w], shift );
        dstPel[w] = val - ( Pel ) IF_INTERNAL_OFFS;
      }

      refPel = refBuf.buf + ( 1 - yOffset )*refBuf.stride - xOffset;
      dstPel = m_filteredBlockTmp[2 + m_iRefListIdx][compID] + 2 * dstBuf.stride + 1;
      for( int h = 0; h < ( height - 2 * BIO_EXTEND_SIZE - 2 ); h++ )
      {
        Pel val = leftShift_round( refPel[0], shift );
        dstPel[0] = val - ( Pel ) IF_INTERNAL_OFFS;

        val = leftShift_round( refPel[width - 3], shift );
        dstPel[width - 3] = val - ( Pel ) IF_INTERNAL_OFFS;

        refPel += refBuf.stride;
        dstPel += dstBuf.stride;
      }

      refPel = refBuf.buf + ( height - 2 * BIO_EXTEND_SIZE - 2 + 1 - yOffset )*refBuf.stride - xOffset;
      dstPel = m_filteredBlockTmp[2 + m_iRefListIdx][compID] + ( height - 2 * BIO_EXTEND_SIZE )*dstBuf.stride + 1;
      for( int w = 0; w < ( width - 2 * BIO_EXTEND_SIZE ); w++ )
      {
        Pel val = leftShift_round( refPel[w], shift );
        dstPel[w] = val - ( Pel ) IF_INTERNAL_OFFS;
      }

      // restore data
      width = backupWidth;
      height = backupHeight;
#endif
      dstBuf.buf = backupDstBufPtr;
      dstBuf.stride = backupDstBufStride;
    }
#if RPR_ENABLE
  }
#endif
#if INTER_LIC
#if RPR_ENABLE
  PelBuf& dstBuf = dstPic.bufs[compID];
#endif

    if( m_storeBeforeLIC )
    {
      m_predictionBeforeLIC.bufs[compID].copyFrom( dstBuf );
    }

#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
    if (pu.cu->LICFlag && (!pu.ciipFlag || doLic))
#else
    if( pu.cu->LICFlag && !pu.ciipFlag )
#endif
    {
      CHECK( pu.cu->geoFlag, "Geometric mode is not used with LIC" );
      CHECK( CU::isIBC( *pu.cu ), "IBC mode is not used with LIC" );
      CHECK( pu.interDir == 3, "Bi-prediction is not used with LIC" );
#if !JVET_W0090_ARMC_TM && !JVET_Z0056_GPM_SPLIT_MODE_REORDERING
      CHECK( pu.ciipFlag, "CIIP mode is not used with LIC" );
#endif
#if RPR_ENABLE
      if (PU::checkRprLicCondition(pu))
      {
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
      if( isAML )
      {
        xLocalIlluComp(pu, compID, *refPic, mvCurr, bi, dstBuf);
      }
      else
#endif
      xLocalIlluComp( pu, compID, *refPic, _mv, bi, dstBuf );
#if RPR_ENABLE
      }
#endif
    }
#endif
#if !RPR_ENABLE
  }
#endif
}

#if !AFFINE_RM_CONSTRAINTS_AND_OPT
bool InterPrediction::isSubblockVectorSpreadOverLimit( int a, int b, int c, int d, int predType )
{
  int s4 = ( 4 << 11 );
  int filterTap = 6;

  if ( predType == 3 )
  {
    int refBlkWidth  = std::max( std::max( 0, 4 * a + s4 ), std::max( 4 * c, 4 * a + 4 * c + s4 ) ) - std::min( std::min( 0, 4 * a + s4 ), std::min( 4 * c, 4 * a + 4 * c + s4 ) );
    int refBlkHeight = std::max( std::max( 0, 4 * b ), std::max( 4 * d + s4, 4 * b + 4 * d + s4 ) ) - std::min( std::min( 0, 4 * b ), std::min( 4 * d + s4, 4 * b + 4 * d + s4 ) );
    refBlkWidth  = ( refBlkWidth >> 11 ) + filterTap + 3;
    refBlkHeight = ( refBlkHeight >> 11 ) + filterTap + 3;

    if ( refBlkWidth * refBlkHeight > ( filterTap + 9 ) * ( filterTap + 9 ) )
    {
      return true;
    }
  }
  else
  {
    int refBlkWidth  = std::max( 0, 4 * a + s4 ) - std::min( 0, 4 * a + s4 );
    int refBlkHeight = std::max( 0, 4 * b ) - std::min( 0, 4 * b );
    refBlkWidth  = ( refBlkWidth >> 11 ) + filterTap + 3;
    refBlkHeight = ( refBlkHeight >> 11 ) + filterTap + 3;
    if ( refBlkWidth * refBlkHeight > ( filterTap + 9 ) * ( filterTap + 5 ) )
    {
      return true;
    }

    refBlkWidth  = std::max( 0, 4 * c ) - std::min( 0, 4 * c );
    refBlkHeight = std::max( 0, 4 * d + s4 ) - std::min( 0, 4 * d + s4 );
    refBlkWidth  = ( refBlkWidth >> 11 ) + filterTap + 3;
    refBlkHeight = ( refBlkHeight >> 11 ) + filterTap + 3;
    if ( refBlkWidth * refBlkHeight > ( filterTap + 5 ) * ( filterTap + 9 ) )
    {
      return true;
    }
  }
  return false;
}
#endif

#if AFFINE_ENC_OPT
#if JVET_Z0136_OOB
void InterPrediction::xPredAffineBlk(const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv* _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng, RefPicList eRefPicList, const bool genChromaMv, const std::pair<int, int> scalingRatio, const bool calGradient)
#else
void InterPrediction::xPredAffineBlk(const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv* _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng, const bool genChromaMv, const std::pair<int, int> scalingRatio, const bool calGradient)
#endif
#else
#if JVET_Z0136_OOB
void InterPrediction::xPredAffineBlk(const ComponentID &compID, const PredictionUnit &pu, const Picture *refPic, const Mv *_mv, PelUnitBuf &dstPic, const bool &bi, const ClpRng &clpRng, RefPicList eRefPicList, bool genChromaMv, const std::pair<int, int> scalingRatio)
#else
void InterPrediction::xPredAffineBlk(const ComponentID &compID, const PredictionUnit &pu, const Picture *refPic, const Mv *_mv, PelUnitBuf &dstPic, const bool &bi, const ClpRng &clpRng, bool genChromaMv, const std::pair<int, int> scalingRatio)
#endif
#endif
{

  JVET_J0090_SET_REF_PICTURE( refPic, compID );
  const ChromaFormat chFmt = pu.chromaFormat;
  int iScaleX = ::getComponentScaleX( compID, chFmt );
  int iScaleY = ::getComponentScaleY( compID, chFmt );

  Mv mvLT =_mv[0];
  Mv mvRT =_mv[1];
  Mv mvLB =_mv[2];
#if INTER_LIC
  Pel* refLeftTemplate = m_pcLICRefLeftTemplate;
  Pel* refAboveTemplate = m_pcLICRefAboveTemplate;
  Pel* recLeftTemplate = m_pcLICRecLeftTemplate;
  Pel* recAboveTemplate = m_pcLICRecAboveTemplate;
  int numTemplate[2] = { 0 , 0 }; // 0:Above, 1:Left
#endif

  // get affine sub-block width and height
  const int width  = pu.Y().width;
  const int height = pu.Y().height;
  int blockWidth = AFFINE_MIN_BLOCK_SIZE;
  int blockHeight = AFFINE_MIN_BLOCK_SIZE;

  CHECK(blockWidth  > (width >> iScaleX ), "Sub Block width  > Block width");
  CHECK(blockHeight > (height >> iScaleY), "Sub Block height > Block height");
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  const int MVBUFFER_SIZE = MAX_CU_SIZE / MIN_PU_SIZE;
#endif

  const int cxWidth  = width  >> iScaleX;
  const int cxHeight = height >> iScaleY;
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  const int iHalfBW  = blockWidth  >> 1;
  const int iHalfBH  = blockHeight >> 1;
#endif

  const int iBit = MAX_CU_DEPTH;
  int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;
#if AFFINE_RM_CONSTRAINTS_AND_OPT
  iDMvHorX = (mvRT - mvLT).getHor() << (iBit - floorLog2(width));
  iDMvHorY = (mvRT - mvLT).getVer() << (iBit - floorLog2(width));
  if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
  {
    iDMvVerX = (mvLB - mvLT).getHor() << (iBit - floorLog2(height));
    iDMvVerY = (mvLB - mvLT).getVer() << (iBit - floorLog2(height));
  }
#else
  iDMvHorX = (mvRT - mvLT).getHor() << (iBit - floorLog2(cxWidth));
  iDMvHorY = (mvRT - mvLT).getVer() << (iBit - floorLog2(cxWidth));
  if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
  {
    iDMvVerX = (mvLB - mvLT).getHor() << (iBit - floorLog2(cxHeight));
    iDMvVerY = (mvLB - mvLT).getVer() << (iBit - floorLog2(cxHeight));
  }
#endif
  else
  {
    iDMvVerX = -iDMvHorY;
    iDMvVerY = iDMvHorX;
  }

  int iMvScaleHor = mvLT.getHor() << iBit;
  int iMvScaleVer = mvLT.getVer() << iBit;
  const SPS &sps    = *pu.cs->sps;

#if IF_12TAP
  const int vFilterSize = isLuma(compID) ? NTAPS_LUMA( 0 ) : NTAPS_CHROMA;
#else
  const int vFilterSize = isLuma(compID) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif

  const int shift = iBit - 4 + MV_FRACTIONAL_BITS_INTERNAL;
  bool      wrapRef = false;
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  const bool subblkMVSpreadOverLimit = isSubblockVectorSpreadOverLimit( iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY, pu.interDir );
#endif

  bool enablePROF = (sps.getUsePROF()) && (!m_skipPROF) && (compID == COMPONENT_Y);
  enablePROF &= (!pu.cs->picHeader->getDisProfFlag());
  enablePROF &= !((pu.cu->affineType == AFFINEMODEL_6PARAM && _mv[0] == _mv[1] && _mv[0] == _mv[2]) || (pu.cu->affineType == AFFINEMODEL_4PARAM && _mv[0] == _mv[1]));
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  enablePROF &= !subblkMVSpreadOverLimit;
#endif
  const int profThres = 1 << (iBit + (m_isBi ? 1 : 0));
  enablePROF &= !m_encOnly || pu.cu->slice->getCheckLDC() || iDMvHorX > profThres || iDMvHorY > profThres || iDMvVerX > profThres || iDMvVerY > profThres || iDMvHorX < -profThres || iDMvHorY < -profThres || iDMvVerX < -profThres || iDMvVerY < -profThres;
  enablePROF &= (refPic->isRefScaled( pu.cs->pps ) == false);
#if AFFINE_MMVD
  enablePROF &= ((pu.mmvdEncOptMode & 3) != 3); // encoder-only
#endif

#if AFFINE_ENC_OPT
  bool isLast = (enablePROF || calGradient) ? false : !bi;
#else
  bool isLast = enablePROF ? false : !bi;
#endif

#if AFFINE_RM_CONSTRAINTS_AND_OPT
  const int cuExtW = width + PROF_BORDER_EXT_W * 2;
  const int cuExtH = height + PROF_BORDER_EXT_H * 2;
#else
  const int cuExtW = AFFINE_MIN_BLOCK_SIZE + PROF_BORDER_EXT_W * 2;
  const int cuExtH = AFFINE_MIN_BLOCK_SIZE + PROF_BORDER_EXT_H * 2;
#endif

  PelBuf gradXExt(m_gradBuf[0], cuExtW, cuExtH);
  PelBuf gradYExt(m_gradBuf[1], cuExtW, cuExtH);
#if IF_12TAP
  const int MAX_FILTER_SIZE = std::max<int>(NTAPS_LUMA(0), NTAPS_CHROMA);
#else
  const int MAX_FILTER_SIZE = std::max<int>(NTAPS_LUMA, NTAPS_CHROMA);
#endif
#if AFFINE_RM_CONSTRAINTS_AND_OPT
  const int dstExtW = ((width + PROF_BORDER_EXT_W * 2 + 7) >> 3) << 3;
  const int dstExtH = cuExtH;
  PelBuf dstExtBuf(m_filteredBlockTmp[1][compID], cuExtW, cuExtH);
#else
  const int dstExtW = ((blockWidth + PROF_BORDER_EXT_W * 2 + 7) >> 3) << 3;
  const int dstExtH = blockHeight + PROF_BORDER_EXT_H * 2;
  PelBuf dstExtBuf(m_filteredBlockTmp[1][compID], dstExtW, dstExtH);
#endif

  const int refExtH = dstExtH + MAX_FILTER_SIZE - 1;
  PelBuf tmpBuf = PelBuf(m_filteredBlockTmp[0][compID], dstExtW, refExtH);

  PelBuf &dstBuf = dstPic.bufs[compID];

  int *dMvScaleHor = m_dMvBuf[m_iRefListIdx];
  int *dMvScaleVer = m_dMvBuf[m_iRefListIdx] + 16;

  if (enablePROF)
  {
    int* dMvH = dMvScaleHor;
    int* dMvV = dMvScaleVer;
    int quadHorX = iDMvHorX << 2;
    int quadHorY = iDMvHorY << 2;
    int quadVerX = iDMvVerX << 2;
    int quadVerY = iDMvVerY << 2;

    dMvH[0] = ((iDMvHorX + iDMvVerX) << 1) - ((quadHorX + quadVerX) << 1);
    dMvV[0] = ((iDMvHorY + iDMvVerY) << 1) - ((quadHorY + quadVerY) << 1);

    for (int w = 1; w < blockWidth; w++)
    {
      dMvH[w] = dMvH[w - 1] + quadHorX;
      dMvV[w] = dMvV[w - 1] + quadHorY;
    }

    dMvH += blockWidth;
    dMvV += blockWidth;
    for (int h = 1; h < blockHeight; h++)
    {
      for (int w = 0; w < blockWidth; w++)
      {
        dMvH[w] = dMvH[w - blockWidth] + quadVerX;
        dMvV[w] = dMvV[w - blockWidth] + quadVerY;
      }
      dMvH += blockWidth;
      dMvV += blockWidth;
    }

#if CTU_256
    const int mvShift = MAX_CU_DEPTH + 1;
    const int dmvLimit = (1 << 5) - 1; // this means the maximum magnitude of dmv is half pel. The target MV precision is 1/64, thus the bit shift is 5
#else
    const int mvShift  = 8;
    const int dmvLimit = ( 1 << 5 ) - 1;
#endif

    if (!g_pelBufOP.roundIntVector)
    {
      for (int idx = 0; idx < blockWidth * blockHeight; idx++)
      {
        roundAffineMv(dMvScaleHor[idx], dMvScaleVer[idx], mvShift);
        dMvScaleHor[idx] = Clip3( -dmvLimit, dmvLimit, dMvScaleHor[idx] );
        dMvScaleVer[idx] = Clip3( -dmvLimit, dmvLimit, dMvScaleVer[idx] );
      }
    }
    else
    {
      int sz = blockWidth * blockHeight;
      g_pelBufOP.roundIntVector(dMvScaleHor, sz, mvShift, dmvLimit);
      g_pelBufOP.roundIntVector(dMvScaleVer, sz, mvShift, dmvLimit);
    }
  }
#if AFFINE_ENC_OPT
  else if (calGradient)
  {
    ::memset(m_dMvBuf, 0, sizeof(m_dMvBuf));
  }
#endif
#if AFFINE_RM_CONSTRAINTS_AND_OPT
  if (compID == COMPONENT_Y)
  {
    if (iDMvHorX == 0 && iDMvHorY == 0)
      blockWidth = width;
    else
    {
      int maxDmv = std::max(abs(iDMvHorX), abs(iDMvHorY)) * blockWidth;
      int TH = 1 << (iBit - 1); // Half pel
      while (maxDmv < TH && blockWidth < width)
      {
        blockWidth <<= 1;
        maxDmv <<= 1;
      }
    }
    if (iDMvVerX == 0 && iDMvVerY == 0)
      blockHeight = height;
    else
    {
      int maxDmv = std::max(abs(iDMvVerX), abs(iDMvVerY)) * blockHeight;
      int TH = 1 << (iBit - 1); // Half pel
      while (maxDmv < TH && blockHeight < height)
      {
        blockHeight <<= 1;
        maxDmv <<= 1;
      }
    }
  }
  CMotionBuf mb = pu.getMotionBuf();
  const MotionInfo  *miLine = mb.buf;
  const MotionInfo  *miLine2 = mb.buf + iScaleX + iScaleY * mb.stride;
  int stride = ((blockHeight << iScaleY) >> 2) * mb.stride;

  int iMvScaleTmpHor0 = iMvScaleHor + ((iDMvHorX * blockWidth + iDMvVerX * blockHeight) >> 1);
  int iMvScaleTmpVer0 = iMvScaleVer + ((iDMvHorY * blockWidth + iDMvVerY * blockHeight) >> 1);
#endif
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  int scaleXLuma = ::getComponentScaleX(COMPONENT_Y, chFmt);
  int scaleYLuma = ::getComponentScaleY(COMPONENT_Y, chFmt);

  if (genChromaMv && pu.chromaFormat != CHROMA_444)
  {
    CHECK(compID == COMPONENT_Y, "Chroma only subblock MV calculation should not apply to Luma");
    int lumaBlockWidth  = AFFINE_MIN_BLOCK_SIZE;
    int lumaBlockHeight = AFFINE_MIN_BLOCK_SIZE;

    CHECK(lumaBlockWidth > (width >> scaleXLuma), "Sub Block width  > Block width");
    CHECK(lumaBlockHeight > (height >> scaleYLuma), "Sub Block height > Block height");

    const int cxWidthLuma  = width >> scaleXLuma;
    const int cxHeightLuma = height >> scaleYLuma;
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
    const int halfBWLuma  = lumaBlockWidth >> 1;
    const int halfBHLuma  = lumaBlockHeight >> 1;

    int dMvHorXLuma, dMvHorYLuma, dMvVerXLuma, dMvVerYLuma;
    dMvHorXLuma = (mvRT - mvLT).getHor() << (iBit - floorLog2(cxWidthLuma));
    dMvHorYLuma = (mvRT - mvLT).getVer() << (iBit - floorLog2(cxWidthLuma));
    if (pu.cu->affineType == AFFINEMODEL_6PARAM)
    {
      dMvVerXLuma = (mvLB - mvLT).getHor() << (iBit - floorLog2(cxHeightLuma));
      dMvVerYLuma = (mvLB - mvLT).getVer() << (iBit - floorLog2(cxHeightLuma));
    }
    else
    {
      dMvVerXLuma = -dMvHorYLuma;
      dMvVerYLuma = dMvHorXLuma;
    }
#endif

#if !AFFINE_RM_CONSTRAINTS_AND_OPT
    const bool subblkMVSpreadOverLimitLuma = isSubblockVectorSpreadOverLimit(dMvHorXLuma, dMvHorYLuma, dMvVerXLuma, dMvVerYLuma, pu.interDir);
#endif

    // get luma MV block by block
    for (int h = 0; h < cxHeightLuma; h += lumaBlockHeight)
    {
      for (int w = 0; w < cxWidthLuma; w += lumaBlockWidth)
      {
        int mvScaleTmpHor, mvScaleTmpVer;
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
        if (!subblkMVSpreadOverLimitLuma)
#endif
        {
#if AFFINE_RM_CONSTRAINTS_AND_OPT
          mvScaleTmpHor = iMvScaleTmpHor0 + iDMvHorX * w + iDMvVerX * h;
          mvScaleTmpVer = iMvScaleTmpVer0 + iDMvHorY * w + iDMvVerY * h;
#else
          mvScaleTmpHor = iMvScaleHor + dMvHorXLuma * (halfBWLuma + w) + dMvVerXLuma * (halfBHLuma + h);
          mvScaleTmpVer = iMvScaleVer + dMvHorYLuma * (halfBWLuma + w) + dMvVerYLuma * (halfBHLuma + h);
#endif
        }
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
        else
        {
          mvScaleTmpHor = iMvScaleHor + dMvHorXLuma * (cxWidthLuma >> 1) + dMvVerXLuma * (cxHeightLuma >> 1);
          mvScaleTmpVer = iMvScaleVer + dMvHorYLuma * (cxWidthLuma >> 1) + dMvVerYLuma * (cxHeightLuma >> 1);
        }
#endif

        roundAffineMv(mvScaleTmpHor, mvScaleTmpVer, shift);
        Mv tmpMv(mvScaleTmpHor, mvScaleTmpVer);
        tmpMv.clipToStorageBitDepth();
        mvScaleTmpHor = tmpMv.getHor();
        mvScaleTmpVer = tmpMv.getVer();

        m_storedMv[h / AFFINE_MIN_BLOCK_SIZE * MVBUFFER_SIZE + w / AFFINE_MIN_BLOCK_SIZE].set(mvScaleTmpHor, mvScaleTmpVer);
      }
    }
  }
#endif
#if AFFINE_ENC_OPT
  int gradLineOffset = 0, gradOffset = 0;
  int gradSubBlkStride = blockHeight * width;
#elif AFFINE_RM_CONSTRAINTS_AND_OPT
  int gradLineOffset = 0, gradOffset = 0;
  int gradSubBlkStride = blockHeight * cuExtW;
#endif
#if JVET_Z0136_OOB
  if (compID == COMPONENT_Y && pu.interDir == 3)
  {
    bool *pMcMask = pu.cs->mcMask[int(eRefPicList)];
    memset(pMcMask, false, cxWidth * cxHeight);
    bool *pMcMaskChroma = pu.cs->mcMaskChroma[int(eRefPicList)];
    int chromaScale = getComponentScaleX(COMPONENT_Cb, m_currChromaFormat);
    int cxWidthChroma = cxWidth >> chromaScale;
    int cxHeightChroma = cxHeight >> chromaScale;
    memset(pMcMaskChroma, false, cxWidthChroma * cxHeightChroma);
  }
#endif
  // get prediction block by block
  for ( int h = 0; h < cxHeight; h += blockHeight )
  {
    for ( int w = 0; w < cxWidth; w += blockWidth )
    {

      int iMvScaleTmpHor, iMvScaleTmpVer;
      if (compID == COMPONENT_Y || pu.chromaFormat == CHROMA_444)
      {
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
        if ( !subblkMVSpreadOverLimit )
#endif
        {
#if AFFINE_RM_CONSTRAINTS_AND_OPT
          iMvScaleTmpHor = iMvScaleTmpHor0 + iDMvHorX * w + iDMvVerX * h;
          iMvScaleTmpVer = iMvScaleTmpVer0 + iDMvHorY * w + iDMvVerY * h;
#else
          iMvScaleTmpHor = iMvScaleHor + iDMvHorX * (iHalfBW + w) + iDMvVerX * (iHalfBH + h);
          iMvScaleTmpVer = iMvScaleVer + iDMvHorY * (iHalfBW + w) + iDMvVerY * (iHalfBH + h);
#endif
        }
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
        else
        {
          iMvScaleTmpHor = iMvScaleHor + iDMvHorX * ( cxWidth >> 1 ) + iDMvVerX * ( cxHeight >> 1 );
          iMvScaleTmpVer = iMvScaleVer + iDMvHorY * ( cxWidth >> 1 ) + iDMvVerY * ( cxHeight >> 1 );
        }
#endif
        roundAffineMv(iMvScaleTmpHor, iMvScaleTmpVer, shift);
        Mv tmpMv(iMvScaleTmpHor, iMvScaleTmpVer);
        tmpMv.clipToStorageBitDepth();
        iMvScaleTmpHor = tmpMv.getHor();
        iMvScaleTmpVer = tmpMv.getVer();

        // clip and scale
        if ( refPic->isWrapAroundEnabled( pu.cs->pps ) )
        {
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
          m_storedMv[h / AFFINE_MIN_BLOCK_SIZE * MVBUFFER_SIZE + w / AFFINE_MIN_BLOCK_SIZE].set(iMvScaleTmpHor, iMvScaleTmpVer);
#endif
          Mv tmpMv(iMvScaleTmpHor, iMvScaleTmpVer);
          wrapRef = wrapClipMv( tmpMv, Position( pu.Y().x + w, pu.Y().y + h ), Size( blockWidth, blockHeight ), &sps, pu.cs->pps );
          iMvScaleTmpHor = tmpMv.getHor();
          iMvScaleTmpVer = tmpMv.getVer();
        }
        else
        {
          wrapRef = false;
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
          m_storedMv[h / AFFINE_MIN_BLOCK_SIZE * MVBUFFER_SIZE + w / AFFINE_MIN_BLOCK_SIZE].set(iMvScaleTmpHor, iMvScaleTmpVer);
#endif
          if (refPic->isRefScaled(pu.cs->pps) == false)
          {
            clipMv(tmpMv, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);
            iMvScaleTmpHor = tmpMv.getHor();
            iMvScaleTmpVer = tmpMv.getVer();
          }
        }
#if JVET_Z0136_OOB
        if (compID == COMPONENT_Y && pu.interDir == 3)
        {
          int chromaScale = getComponentScaleX(COMPONENT_Cb, m_currChromaFormat);
          bool *pMcMask = pu.cs->mcMask[int(eRefPicList)] + w + h * cxWidth;
          bool *pMcMaskChroma = pu.cs->mcMaskChroma[int(eRefPicList)] + (w >> chromaScale) + (h >> chromaScale) * (cxWidth >> chromaScale);
          int cxWidthChroma = cxWidth >> chromaScale;

          isMvOOBSubBlk(tmpMv, Position(pu.Y().x + w, pu.Y().y + h), Size(blockWidth, blockHeight), pu.cu->slice->getSPS(), pu.cu->slice->getPPS(), pMcMask, cxWidth, pMcMaskChroma, cxWidthChroma);
        }
#endif
      }
      else
      {
#if AFFINE_RM_CONSTRAINTS_AND_OPT
        Mv curMv = miLine[(w << iScaleX) >> 2].mv[m_iRefListIdx] + miLine2[(w << iScaleX) >> 2].mv[m_iRefListIdx];
#else
        Mv curMv = m_storedMv[((h << iScaleY) / AFFINE_MIN_BLOCK_SIZE) * MVBUFFER_SIZE + ((w << iScaleX) / AFFINE_MIN_BLOCK_SIZE)] +
          m_storedMv[((h << iScaleY) / AFFINE_MIN_BLOCK_SIZE + iScaleY)* MVBUFFER_SIZE + ((w << iScaleX) / AFFINE_MIN_BLOCK_SIZE + iScaleX)];
#endif
        roundAffineMv(curMv.hor, curMv.ver, 1);
        if ( refPic->isWrapAroundEnabled( pu.cs->pps ) )
        {
          wrapRef = wrapClipMv( curMv, Position( pu.Y().x + ( w << iScaleX ), pu.Y().y + ( h << iScaleY ) ), Size( blockWidth << iScaleX, blockHeight << iScaleY ), &sps, pu.cs->pps );
        }
        else
        {
          wrapRef = false;
          if (refPic->isRefScaled(pu.cs->pps) == false)
          {
            clipMv(curMv, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);
          }
        }
        iMvScaleTmpHor = curMv.hor;
        iMvScaleTmpVer = curMv.ver;
      }

      if( xPredInterBlkRPR( scalingRatio, *pu.cs->pps, CompArea( compID, chFmt, pu.blocks[compID].offset( w, h ), Size( blockWidth, blockHeight ) ), refPic, Mv( iMvScaleTmpHor, iMvScaleTmpVer ), dstBuf.buf + w + h * dstBuf.stride, dstBuf.stride, bi, wrapRef, clpRng, 2 ) )
      {
        CHECK( enablePROF, "PROF should be disabled with RPR" );
      }
      else
      {
#if INTER_LIC
        if (pu.cu->LICFlag && (w == 0 || h == 0))
        {
          xGetSublkTemplate(*pu.cu, compID, *refPic, Mv(iMvScaleTmpHor, iMvScaleTmpVer), blockWidth, blockHeight, w, h, numTemplate, refLeftTemplate, refAboveTemplate, recLeftTemplate, recAboveTemplate);
        }
#endif
        // get the MV in high precision
        int xFrac, yFrac, xInt, yInt;

        if (!iScaleX)
        {
          xInt  = iMvScaleTmpHor >> 4;
          xFrac = iMvScaleTmpHor & 15;
        }
        else
        {
          xInt  = iMvScaleTmpHor >> 5;
          xFrac = iMvScaleTmpHor & 31;
        }
        if (!iScaleY)
        {
          yInt  = iMvScaleTmpVer >> 4;
          yFrac = iMvScaleTmpVer & 15;
        }
        else
        {
          yInt  = iMvScaleTmpVer >> 5;
          yFrac = iMvScaleTmpVer & 31;
        }

        const CPelBuf refBuf = refPic->getRecoBuf(
          CompArea(compID, chFmt, pu.blocks[compID].offset(xInt + w, yInt + h), pu.blocks[compID]), wrapRef);

        Pel *ref = (Pel *) refBuf.buf;
        Pel *dst = dstBuf.buf + w + h * dstBuf.stride;

        int refStride = refBuf.stride;
        int dstStride = dstBuf.stride;

        int bw = blockWidth;
        int bh = blockHeight;

#if AFFINE_ENC_OPT
        if (enablePROF || calGradient)
#else
        if (enablePROF)
#endif
        {
          dst       = dstExtBuf.bufAt(PROF_BORDER_EXT_W, PROF_BORDER_EXT_H);
          dstStride = dstExtBuf.stride;
        }

#if IF_12TAP
        if (yFrac == 0)
        {
          m_if.filterHor(compID, (Pel*)ref, refStride, dst, dstStride, bw, bh, xFrac, isLast, chFmt, clpRng
            , 0, false, false);
        }
        else if (xFrac == 0)
        {
          m_if.filterVer(compID, (Pel*)ref, refStride, dst, dstStride, bw, bh, yFrac, true, isLast, chFmt, clpRng
            , 0, false, false);
        }
        else
        {
#if SIMD_4x4_12 && defined(TARGET_SIMD_X86)
#if AFFINE_RM_CONSTRAINTS_AND_OPT
          if (compID == COMPONENT_Y && bw == 4 && bh == 4)
#else
          if (compID == COMPONENT_Y)
#endif
            m_if.filter4x4(clpRng, (Pel*)ref, refStride, dst, dstStride, xFrac, yFrac, isLast);
          else {
#endif
            m_if.filterHor(compID, (Pel*)ref - ((vFilterSize >> 1) - 1)*refStride, refStride, tmpBuf.buf, tmpBuf.stride, bw, bh + vFilterSize - 1, xFrac, false, chFmt, clpRng
              , 0, false, false);
            JVET_J0090_SET_CACHE_ENABLE(false);
            m_if.filterVer(compID, tmpBuf.buf + ((vFilterSize >> 1) - 1)*tmpBuf.stride, tmpBuf.stride, dst, dstStride, bw, bh, yFrac, false, isLast, chFmt, clpRng
              , 0, false, false);
            JVET_J0090_SET_CACHE_ENABLE(true);
#if SIMD_4x4_12 && defined(TARGET_SIMD_X86)
          }
#endif
        }
#else
        if (yFrac == 0)
        {
          m_if.filterHor(compID, (Pel *) ref, refStride, dst, dstStride, bw, bh, xFrac, isLast, chFmt, clpRng);
        }
        else if (xFrac == 0)
        {
          m_if.filterVer(compID, (Pel *) ref, refStride, dst, dstStride, bw, bh, yFrac, true, isLast, chFmt, clpRng);
        }
        else
        {
          m_if.filterHor(compID, (Pel *) ref - ((vFilterSize >> 1) - 1) * refStride, refStride, tmpBuf.buf,
                         tmpBuf.stride, bw, bh + vFilterSize - 1, xFrac, false, chFmt, clpRng);
          JVET_J0090_SET_CACHE_ENABLE(false);
          m_if.filterVer(compID, tmpBuf.buf + ((vFilterSize >> 1) - 1) * tmpBuf.stride, tmpBuf.stride, dst, dstStride,
                         bw, bh, yFrac, false, isLast, chFmt, clpRng);
          JVET_J0090_SET_CACHE_ENABLE(true);
        }
#endif
#if AFFINE_ENC_OPT
        if (enablePROF || calGradient)
#else
        if (enablePROF)
#endif
        {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
          const int shift = IF_INTERNAL_FRAC_BITS(clpRng.bd);
#else
          const int shift    = std::max<int>(2, (IF_INTERNAL_PREC - clpRng.bd));
#endif
          const int xOffset = xFrac >> 3;
          const int yOffset = yFrac >> 3;

          const int refOffset = (blockHeight + 1) * refStride;
          const int dstOffset = (blockHeight + 1) * dstStride;

          const Pel *refPel = ref - (1 - yOffset) * refStride + xOffset - 1;
          Pel *      dstPel = dst - dstStride - 1;
          for (int pw = 0; pw < blockWidth + 2; pw++)
          {
            dstPel[pw]             = leftShift_round(refPel[pw], shift) - (Pel) IF_INTERNAL_OFFS;
            dstPel[pw + dstOffset] = leftShift_round(refPel[pw + refOffset], shift) - (Pel) IF_INTERNAL_OFFS;
          }

          refPel = ref + yOffset * refBuf.stride + xOffset;
          dstPel = dst;
          for (int ph = 0; ph < blockHeight; ph++, refPel += refStride, dstPel += dstStride)
          {
            dstPel[-1]         = leftShift_round(refPel[-1], shift) - (Pel) IF_INTERNAL_OFFS;
            dstPel[blockWidth] = leftShift_round(refPel[blockWidth], shift) - (Pel) IF_INTERNAL_OFFS;
          }

#if AFFINE_RM_CONSTRAINTS_AND_OPT
          gradOffset = gradLineOffset + w;
#if AFFINE_ENC_OPT
          g_pelBufOP.profGradFilter(dstExtBuf.buf, dstExtBuf.stride, blockWidth + 2, blockHeight + 2, width, m_gradX0 + gradOffset, m_gradY0 + gradOffset, clpRng.bd);
#else
          g_pelBufOP.profGradFilter(dstExtBuf.buf, dstExtBuf.stride, blockWidth + 2, blockHeight + 2, cuExtW, m_gradX0 + gradOffset, m_gradY0 + gradOffset, clpRng.bd);
#endif
#else
#if AFFINE_ENC_OPT
          gradOffset = gradLineOffset + w;
          g_pelBufOP.profGradFilter(dstExtBuf.buf, dstExtBuf.stride, blockWidth + 2, blockHeight + 2, width, m_gradX0 + gradOffset, m_gradY0 + gradOffset, clpRng.bd);
#else
          PelBuf gradXBuf = gradXExt.subBuf(0, 0, blockWidth + 2, blockHeight + 2);
          PelBuf gradYBuf = gradYExt.subBuf(0, 0, blockWidth + 2, blockHeight + 2);
          g_pelBufOP.profGradFilter(dstExtBuf.buf, dstExtBuf.stride, blockWidth + 2, blockHeight + 2, gradXBuf.stride,
                                    gradXBuf.buf, gradYBuf.buf, clpRng.bd);
#endif
#endif

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
          const Pel offset = (1 << (shift - 1)) + IF_INTERNAL_OFFS;
#else
          const int shiftNum = std::max<int>(2, (IF_INTERNAL_PREC - clpRng.bd));
          const Pel offset   = (1 << (shiftNum - 1)) + IF_INTERNAL_OFFS;
#endif
#if AFFINE_RM_CONSTRAINTS_AND_OPT
          Pel* src = dst;
#if AFFINE_ENC_OPT
          Pel* gX = m_gradX0 + gradOffset + width + 1;
          Pel* gY = m_gradY0 + gradOffset + width + 1;
#else
          Pel* gX = m_gradX0 + gradOffset + cuExtW + 1;
          Pel* gY = m_gradY0 + gradOffset + cuExtW + 1;
#endif
          Pel * dstY = dstBuf.buf + w + h * dstBuf.stride;
          for (int sh = 0; sh < blockHeight; sh += AFFINE_MIN_BLOCK_SIZE)
          {
            for (int sw = 0; sw < blockWidth; sw += AFFINE_MIN_BLOCK_SIZE)
            {
#if AFFINE_ENC_OPT
              g_pelBufOP.applyPROF(dstY + sw, dstBuf.stride, src + sw, dstExtBuf.stride, AFFINE_MIN_BLOCK_SIZE, AFFINE_MIN_BLOCK_SIZE, gX + sw, gY + sw, width, dMvScaleHor, dMvScaleVer, AFFINE_MIN_BLOCK_SIZE, bi, shift, offset, clpRng);
#else
              g_pelBufOP.applyPROF(dstY + sw, dstBuf.stride, src + sw, dstExtBuf.stride, AFFINE_MIN_BLOCK_SIZE, AFFINE_MIN_BLOCK_SIZE, gX + sw, gY + sw, cuExtW, dMvScaleHor, dMvScaleVer, AFFINE_MIN_BLOCK_SIZE, bi, shift, offset, clpRng);
#endif
            }
            src += (dstStride << 2);
#if AFFINE_ENC_OPT
            gX += (width << 2);
            gY += (width << 2);
#else
            gX += (cuExtW << 2);
            gY += (cuExtW << 2);
#endif
            dstY += (dstBuf.stride << 2);
          }
#else
          Pel *src = dstExtBuf.bufAt(PROF_BORDER_EXT_W, PROF_BORDER_EXT_H);
#if AFFINE_ENC_OPT
          Pel* gX = m_gradX0 + gradOffset + width + 1;
          Pel* gY = m_gradY0 + gradOffset + width + 1;
          Pel * dstY = dstBuf.bufAt(w, h);
          g_pelBufOP.applyPROF(dstY, dstBuf.stride, src, dstExtBuf.stride, blockWidth, blockHeight, gX, gY, width, dMvScaleHor, dMvScaleVer, blockWidth, bi, shift, offset, clpRng);
#else
          Pel *gX  = gradXBuf.bufAt(PROF_BORDER_EXT_W, PROF_BORDER_EXT_H);
          Pel *gY  = gradYBuf.bufAt(PROF_BORDER_EXT_W, PROF_BORDER_EXT_H);

          Pel *dstY = dstBuf.bufAt(w, h);

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
          g_pelBufOP.applyPROF(dstY, dstBuf.stride, src, dstExtBuf.stride, blockWidth, blockHeight, gX, gY,
                               gradXBuf.stride, dMvScaleHor, dMvScaleVer, blockWidth, bi, shift, offset, clpRng);
#else
          g_pelBufOP.applyPROF(dstY, dstBuf.stride, src, dstExtBuf.stride, blockWidth, blockHeight, gX, gY,
                               gradXBuf.stride, dMvScaleHor, dMvScaleVer, blockWidth, bi, shiftNum, offset, clpRng);
#endif
#endif
#endif
        }
      }
    }
#if AFFINE_RM_CONSTRAINTS_AND_OPT || AFFINE_ENC_OPT
    gradLineOffset += gradSubBlkStride;
#endif
#if AFFINE_RM_CONSTRAINTS_AND_OPT
    miLine += stride;
    miLine2 += stride;
#endif
  }
#if INTER_LIC
  if (m_storeBeforeLIC)
  {
    m_predictionBeforeLIC.bufs[compID].copyFrom(dstBuf);
  }

#if RPR_ENABLE
  if( pu.cu->LICFlag && PU::checkRprLicCondition( pu ) )
#else
  if (pu.cu->LICFlag)
#endif
  {
    PelBuf &dstBuf = dstPic.bufs[compID];
    int LICshift = 0, scale = 0, offset = 0;
    xGetLICParamGeneral(*pu.cu, compID, numTemplate, refLeftTemplate, refAboveTemplate, recLeftTemplate, recAboveTemplate, LICshift, scale, offset);

    const ClpRng& clpRng = pu.cu->cs->slice->clpRng(compID);
    dstBuf.linearTransform(scale, LICshift, offset, true, clpRng);
  }
#endif
}

#if MULTI_PASS_DMVR
#if JVET_Z0136_OOB
void InterPrediction::applyBiOptFlow(const bool isBdofMvRefine, const int bdofBlockOffset, const PredictionUnit &pu, const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1, PelUnitBuf &yuvDst, const BitDepths &clipBitDepths, bool *mcMask[2], bool *mcMaskChroma[2], bool *isOOB)
#else
void InterPrediction::applyBiOptFlow(const bool isBdofMvRefine, const int bdofBlockOffset, const PredictionUnit &pu, const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1, PelUnitBuf &yuvDst, const BitDepths &clipBitDepths)
#endif
#else
#if JVET_Z0136_OOB
void InterPrediction::applyBiOptFlow(const PredictionUnit &pu, const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1, PelUnitBuf &yuvDst, const BitDepths &clipBitDepths, bool *mcMask[2], bool *mcMaskChroma[2], bool *isOOB)
#else
void InterPrediction::applyBiOptFlow(const PredictionUnit &pu, const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1, PelUnitBuf &yuvDst, const BitDepths &clipBitDepths)
#endif
#endif
{
  const int     height = yuvDst.Y().height;
  const int     width = yuvDst.Y().width;
  int           heightG = height + 2 * BIO_EXTEND_SIZE;
  int           widthG = width + 2 * BIO_EXTEND_SIZE;
  int           offsetPos = widthG*BIO_EXTEND_SIZE + BIO_EXTEND_SIZE;

  Pel*          gradX0 = m_gradX0;
  Pel*          gradX1 = m_gradX1;
  Pel*          gradY0 = m_gradY0;
  Pel*          gradY1 = m_gradY1;

  int           stridePredMC = widthG + 2;
  const Pel*    srcY0 = m_filteredBlockTmp[2][COMPONENT_Y] + stridePredMC + 1;
  const Pel*    srcY1 = m_filteredBlockTmp[3][COMPONENT_Y] + stridePredMC + 1;
  const int     src0Stride = stridePredMC;
  const int     src1Stride = stridePredMC;

  Pel*          dstY = yuvDst.Y().buf;
  const int     dstStride = yuvDst.Y().stride;
  const Pel*    srcY0Temp = srcY0;
  const Pel*    srcY1Temp = srcY1;

  for (int refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
  {
    Pel* dstTempPtr = m_filteredBlockTmp[2 + refList][COMPONENT_Y] + stridePredMC + 1;
    Pel* gradY = (refList == 0) ? m_gradY0 : m_gradY1;
    Pel* gradX = (refList == 0) ? m_gradX0 : m_gradX1;

    xBioGradFilter(dstTempPtr, stridePredMC, widthG, heightG, widthG, gradX, gradY, clipBitDepths.recon[toChannelType(COMPONENT_Y)]);
#if !MULTI_PASS_DMVR && !SAMPLE_BASED_BDOF
    Pel* padStr = m_filteredBlockTmp[2 + refList][COMPONENT_Y] + 2 * stridePredMC + 2;
    for (int y = 0; y< height; y++)
    {
      padStr[-1] = padStr[0];
      padStr[width] = padStr[width - 1];
      padStr += stridePredMC;
    }

    padStr = m_filteredBlockTmp[2 + refList][COMPONENT_Y] + 2 * stridePredMC + 1;
    ::memcpy(padStr - stridePredMC, padStr, sizeof(Pel)*(widthG));
    ::memcpy(padStr + height*stridePredMC, padStr + (height - 1)*stridePredMC, sizeof(Pel)*(widthG));
#endif
  }

  const ClpRng& clpRng = pu.cu->cs->slice->clpRng(COMPONENT_Y);
  const int   bitDepth = clipBitDepths.recon[toChannelType(COMPONENT_Y)];
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int   shiftNum = IF_INTERNAL_FRAC_BITS(bitDepth) + 1;
#else
  const int   shiftNum = IF_INTERNAL_PREC + 1 - bitDepth;
#endif
  const int   offset = (1 << (shiftNum - 1)) + 2 * IF_INTERNAL_OFFS;
  const int   limit = ( 1 << 4 ) - 1;

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  int srcBlockOffset = (stridePredMC + 1) * BIO_EXTEND_SIZE;
  int bioBlockParamOffset = (widthG + 1);
  int dstBlockOffset = 0;
  const int bioDx = (width < BDOF_SUBPU_DIM) ? width : BDOF_SUBPU_DIM;
  const int bioDy = (height < BDOF_SUBPU_DIM) ? height : BDOF_SUBPU_DIM;
  const int srcBlockOffsetIncrementY = (stridePredMC << BDOF_SUBPU_DIM_LOG2) - width;
  const int dstBlockOffsetIncrementY = (dstStride << BDOF_SUBPU_DIM_LOG2) - width;
  const int bioBlockParamOffsetIncrementY = (widthG << BDOF_SUBPU_DIM_LOG2) - width;
#endif
#if MULTI_PASS_DMVR
  if (isBdofMvRefine)
  {
    g_pelBufOP.calcBIOParameter(srcY0, srcY1, gradX0, gradX1, gradY0, gradY1, widthG, heightG, src0Stride, src1Stride, widthG,
                                bitDepth, m_absGx, m_absGy, m_dIx, m_dIy, m_signGxGy, m_dI);
    m_bdofMvRefined = true;
    int bioSubPuMvIndex = 0;
    const int bioSubPuMvIndexIncrementY = BDOF_SUBPU_STRIDE - std::max(1, (width >> BDOF_SUBPU_DIM_LOG2));
    const int   bioBlockDistTh = (bioDx * bioDy) << (5 - 4); //4 is to compensate the shift4 of dI in calcBIOParameter
    Pel* dI = m_dI + 2 + 2 * widthG;
    for (int yBlock = 0; yBlock < height; yBlock += bioDy)
    {
      for (int xBlock = 0; xBlock < width; xBlock += bioDx)
      {
        srcY0Temp = srcY0 + srcBlockOffset;
        srcY1Temp = srcY1 + srcBlockOffset;

        int costSubblockSAD = 0;
        Pel* tmp = dI + bioBlockParamOffset;
        g_pelBufOP.calAbsSum(tmp, widthG, bioDx, bioDy, &costSubblockSAD);

        if (costSubblockSAD < bioBlockDistTh)
        {
#if JVET_Z0136_OOB
          int maskOffset = yBlock * width + xBlock;
          bool *pSubMcMask[2] = { pu.cs->mcMask[0] + maskOffset, pu.cs->mcMask[1] + maskOffset };
          bool isOOBTmp[2] = { false, false };

          if (isOOB[0] || isOOB[1])
          {
            for (int dir = 0; dir < 2; dir++)
            {
              bool *pMcMask = (dir == 0) ? pSubMcMask[0] : pSubMcMask[1];
              for (int y = 0; y < bioDy && !isOOBTmp[dir]; y++)
              {
                for (int x = 0; x < bioDx && !isOOBTmp[dir]; x++)
                {
                  isOOBTmp[dir] |= pMcMask[x];
                }
                pMcMask += width;
              }
            }
          }
          m_bdofSubPuMvOffset[bdofBlockOffset + bioSubPuMvIndex].setZero();
          if (bioDx == 4)
          {
            g_pelBufOP.addAvg4(srcY0Temp, src0Stride, srcY1Temp, src1Stride, dstY + dstBlockOffset,
              dstStride, bioDx, bioDy, shiftNum, offset, clpRng, pSubMcMask, width, isOOBTmp);
          }
          else
          {
            g_pelBufOP.addAvg8(srcY0Temp, src0Stride, srcY1Temp, src1Stride, dstY + dstBlockOffset,
              dstStride, bioDx, bioDy, shiftNum, offset, clpRng, pSubMcMask, width, isOOBTmp);
          }
#else
          m_bdofSubPuMvOffset[bdofBlockOffset + bioSubPuMvIndex].setZero();
          if (bioDx == 4)
          {
            g_pelBufOP.addAvg4(srcY0Temp, src0Stride, srcY1Temp, src1Stride, dstY + dstBlockOffset,
                               dstStride, bioDx, bioDy, shiftNum, offset, clpRng);
          }
          else
          {
            g_pelBufOP.addAvg8(srcY0Temp, src0Stride, srcY1Temp, src1Stride, dstY + dstBlockOffset,
                               dstStride, bioDx, bioDy, shiftNum, offset, clpRng);
          }
#endif
          srcBlockOffset += bioDx;
          dstBlockOffset += bioDx;
          bioBlockParamOffset += bioDx;
          bioSubPuMvIndex += 1;
          continue;
        }
        if (!pu.bdmvrRefine)
        {
          m_bdofSubPuMvOffset[bdofBlockOffset + bioSubPuMvIndex].setZero();
#if JVET_Z0136_OOB
          int maskOffset = yBlock * width + xBlock;
          bool *pSubMcMask[2] = { pu.cs->mcMask[0] + maskOffset, pu.cs->mcMask[1] + maskOffset };
          bool isOOBTmp[2] = { false, false };
          if (isOOB[0] || isOOB[1])
          {
            for (int dir = 0; dir < 2; dir++)
            {
              bool *pMcMask = (dir == 0) ? pSubMcMask[0] : pSubMcMask[1];
              for (int y = 0; y < bioDy && !isOOBTmp[dir]; y++)
              {
                for (int x = 0; x < bioDx && !isOOBTmp[dir]; x++)
                {
                  isOOBTmp[dir] |= pMcMask[x];
                }
                pMcMask += width;
              }
            }
          }
          subBlockBiOptFlow(dstY + dstBlockOffset, dstStride, srcY0Temp, src0Stride, srcY1Temp, src1Stride,
            bioBlockParamOffset, widthG, bioDx, bioDy, clpRng, shiftNum, offset, limit, pSubMcMask, width, isOOBTmp);
#else
          subBlockBiOptFlow(dstY + dstBlockOffset, dstStride, srcY0Temp, src0Stride, srcY1Temp, src1Stride,
                            bioBlockParamOffset, widthG, bioDx, bioDy, clpRng, shiftNum, offset, limit);
#endif
          srcBlockOffset += bioDx;
          dstBlockOffset += bioDx;
          bioBlockParamOffset += bioDx;
          bioSubPuMvIndex += 1;
          continue;
        }

        int sumAbsGX_block = 0, sumAbsGY_block = 0, sumDIX_block = 0, sumDIY_block = 0, sumSignGY_GX_block = 0;
        g_pelBufOP.calcBIOParamSum4(m_absGx + bioBlockParamOffset, m_absGy + bioBlockParamOffset, m_dIx + bioBlockParamOffset,
                                    m_dIy + bioBlockParamOffset, m_signGxGy + bioBlockParamOffset, bioDx + 4, bioDy + 4, widthG,
                                    &sumAbsGX_block, &sumAbsGY_block, &sumDIX_block, &sumDIY_block, &sumSignGY_GX_block);

        int tmpx_block = (sumAbsGX_block == 0 ? 0 : rightShiftMSB(sumDIX_block << 3, sumAbsGX_block));
        int tmpData_block = ((tmpx_block * sumSignGY_GX_block) >> 1);
        int tmpy_block = (sumAbsGY_block == 0 ? 0 : rightShiftMSB(((sumDIY_block << 3) - tmpData_block), sumAbsGY_block));
        tmpx_block = Clip3(-256, 256, tmpx_block);
        tmpy_block = Clip3(-256, 256, tmpy_block);

        Mv bioMv;
        if (tmpx_block >= 0)
          bioMv.hor = ((tmpx_block + 4) >> 3);
        else
          bioMv.hor = (-1) * ((((-1) * tmpx_block) + 4) >> 3);
        if (tmpy_block >= 0)
          bioMv.ver = ((tmpy_block + 4) >> 3);
        else
          bioMv.ver = (-1) * ((((-1) * tmpy_block) + 4) >> 3);

        m_bdofSubPuMvOffset[bdofBlockOffset + bioSubPuMvIndex] = bioMv;
        if (bioMv.hor == 0 && bioMv.ver == 0)
        {
          // by doing this, we do not need to do second LUMA MC
#if JVET_Z0136_OOB
          int maskOffset = yBlock * width + xBlock;
          bool *pSubMcMask[2] = { pu.cs->mcMask[0] + maskOffset, pu.cs->mcMask[1] + maskOffset };
          bool isOOBTmp[2] = { false, false };
          if (isOOB[0] || isOOB[1])
          {
            for (int dir = 0; dir < 2; dir++)
            {
              bool *pMcMask = (dir == 0) ? pSubMcMask[0] : pSubMcMask[1];
              for (int y = 0; y < bioDy && !isOOBTmp[dir]; y++)
              {
                for (int x = 0; x < bioDx && !isOOBTmp[dir]; x++)
                {
                  isOOBTmp[dir] |= pMcMask[x];
                }
                pMcMask += width;
              }
            }
          }
          subBlockBiOptFlow(dstY + dstBlockOffset, dstStride, srcY0Temp, src0Stride, srcY1Temp, src1Stride,
                            bioBlockParamOffset, widthG, bioDx, bioDy, clpRng, shiftNum, offset, limit, pSubMcMask, width, isOOBTmp);
#else
          subBlockBiOptFlow(dstY + dstBlockOffset, dstStride, srcY0Temp, src0Stride, srcY1Temp, src1Stride,
                            bioBlockParamOffset, widthG, bioDx, bioDy, clpRng, shiftNum, offset, limit);
#endif
        }
        srcBlockOffset += bioDx;
        dstBlockOffset += bioDx;
        bioBlockParamOffset += bioDx;
        bioSubPuMvIndex += 1;
      }
      srcBlockOffset += srcBlockOffsetIncrementY;
      dstBlockOffset += dstBlockOffsetIncrementY;
      bioBlockParamOffset += bioBlockParamOffsetIncrementY;
      bioSubPuMvIndex += bioSubPuMvIndexIncrementY;
    }
    return;
  }
#endif
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  g_pelBufOP.calcBIOParameter(srcY0, srcY1, gradX0, gradX1, gradY0, gradY1, widthG, heightG, src0Stride, src1Stride, widthG,
                              bitDepth, m_absGx, m_absGy, m_dIx, m_dIy, m_signGxGy, nullptr);
  for (int yBlock = 0; yBlock < height; yBlock += bioDy)
  {
    for (int xBlock = 0; xBlock < width; xBlock += bioDx)
    {
      srcY0Temp = srcY0 + srcBlockOffset;
      srcY1Temp = srcY1 + srcBlockOffset;
#if JVET_Z0136_OOB
      int maskOffset = yBlock * width + xBlock;
      bool *pSubMcMask[2] = { pu.cs->mcMask[0] + maskOffset, pu.cs->mcMask[1] + maskOffset };
      bool isOOBTmp[2] = { false, false };
      if (isOOB[0] || isOOB[1])
      {
        for (int dir = 0; dir < 2; dir++)
        {
          bool *pMcMask = (dir == 0) ? pSubMcMask[0] : pSubMcMask[1];
          for (int y = 0; y < bioDy && !isOOBTmp[dir]; y++)
          {
            for (int x = 0; x < bioDx && !isOOBTmp[dir]; x++)
            {
              isOOBTmp[dir] |= pMcMask[x];
            }
            pMcMask += width;
          }
        }
      }
      subBlockBiOptFlow(dstY + dstBlockOffset, dstStride, srcY0Temp, src0Stride, srcY1Temp, src1Stride,
                        bioBlockParamOffset, widthG, bioDx, bioDy, clpRng, shiftNum, offset, limit, pSubMcMask, width, isOOBTmp);
#else
      subBlockBiOptFlow(dstY + dstBlockOffset, dstStride, srcY0Temp, src0Stride, srcY1Temp, src1Stride,
                        bioBlockParamOffset, widthG, bioDx, bioDy, clpRng, shiftNum, offset, limit);
#endif
      srcBlockOffset += bioDx;
      dstBlockOffset += bioDx;
      bioBlockParamOffset += bioDx;
    }
    srcBlockOffset += srcBlockOffsetIncrementY;
    dstBlockOffset += dstBlockOffsetIncrementY;
    bioBlockParamOffset += bioBlockParamOffsetIncrementY;
  }
  return;
#endif

  int xUnit = (width >> 2);
  int yUnit = (height >> 2);

  Pel *dstY0 = dstY;
  gradX0 = m_gradX0; gradX1 = m_gradX1;
  gradY0 = m_gradY0; gradY1 = m_gradY1;

  Pel *pGradX0Tmp, *pGradX1Tmp, *pGradY0Tmp, *pGradY1Tmp;
  const Pel *SrcY0Tmp, *SrcY1Tmp;
  int tmpx = 0, tmpy = 0;
  int sumAbsGX = 0, sumAbsGY = 0, sumDIX = 0, sumDIY = 0, sumSignGY_GX = 0;
  int gradOfst, srcOfst, dstOfst, gradLineOfst = 0, srcLineOfst = 0, dstLineOfst = 0;

  for (int yu = 0; yu < yUnit; yu++)
  {
    gradOfst = gradLineOfst;
    srcOfst = srcLineOfst;
    dstOfst = dstLineOfst;

    for (int xu = 0; xu < xUnit; xu++)
    {
      sumAbsGX = 0; sumAbsGY = 0; sumDIX = 0; sumDIY = 0, sumSignGY_GX = 0;
      pGradX0Tmp = m_gradX0 + gradOfst;
      pGradX1Tmp = m_gradX1 + gradOfst;
      pGradY0Tmp = m_gradY0 + gradOfst;
      pGradY1Tmp = m_gradY1 + gradOfst;
      SrcY1Tmp = srcY1 + srcOfst;
      SrcY0Tmp = srcY0 + srcOfst;

      g_pelBufOP.calcBIOSums(SrcY0Tmp, SrcY1Tmp, pGradX0Tmp, pGradX1Tmp, pGradY0Tmp, pGradY1Tmp, xu, yu, src0Stride, src1Stride, widthG, bitDepth, &sumAbsGX, &sumAbsGY, &sumDIX, &sumDIY, &sumSignGY_GX);
      tmpx = (sumAbsGX == 0 ? 0 : rightShiftMSB(sumDIX << 2, sumAbsGX));
      tmpx = Clip3(-limit, limit, tmpx);

      int     mainsGxGy = sumSignGY_GX >> 12;
      int     secsGxGy = sumSignGY_GX & ((1 << 12) - 1);
      int     tmpData = tmpx * mainsGxGy;
      tmpData = ((tmpData << 12) + tmpx*secsGxGy) >> 1;
      tmpy = (sumAbsGY == 0 ? 0 : rightShiftMSB(((sumDIY << 2) - tmpData), sumAbsGY));
      tmpy = Clip3(-limit, limit, tmpy);

      srcY0Temp = SrcY0Tmp + ( stridePredMC + 1 );
      srcY1Temp = SrcY1Tmp + ( stridePredMC + 1 );
      gradX0 = pGradX0Tmp + offsetPos;
      gradX1 = pGradX1Tmp + offsetPos;
      gradY0 = pGradY0Tmp + offsetPos;
      gradY1 = pGradY1Tmp + offsetPos;
      dstY0 = dstY + dstOfst;

      gradOfst += 4;
      srcOfst += 4;
      dstOfst += 4;

      xAddBIOAvg4(srcY0Temp, src0Stride, srcY1Temp, src1Stride, dstY0, dstStride, gradX0, gradX1, gradY0, gradY1, widthG, (1 << 2), (1 << 2), (int)tmpx, (int)tmpy, shiftNum, offset, clpRng);
    }  // xu

    gradLineOfst += ( widthG << 2 );
    srcLineOfst += ( src0Stride << 2 );
    dstLineOfst += ( dstStride << 2 );
  }  // yu
}

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
#if JVET_Z0136_OOB
void InterPrediction::subBlockBiOptFlow(Pel* dstY, const int dstStride, const Pel* src0, const int src0Stride, const Pel* src1, const int src1Stride, int bioParamOffset, const int bioParamStride, int width, int height, const ClpRng& clpRng, const int shiftNum, const int offset, const int limit, bool *mcMask[2], int mcStride, bool *isOOB)
#else
void InterPrediction::subBlockBiOptFlow(Pel* dstY, const int dstStride, const Pel* src0, const int src0Stride, const Pel* src1, const int src1Stride, int bioParamOffset, const int bioParamStride, int width, int height, const ClpRng& clpRng, const int shiftNum, const int offset, const int limit)
#endif
{
#if SAMPLE_BASED_BDOF
  g_pelBufOP.calcBIOParamSum5(m_absGx + bioParamOffset, m_absGy + bioParamOffset, m_dIx + bioParamOffset,
                              m_dIy + bioParamOffset, m_signGxGy + bioParamOffset, bioParamStride, width, height,
                              m_sumAbsGX_pixel_32bit, m_sumAbsGY_pixel_32bit, m_sumDIX_pixel_32bit, m_sumDIY_pixel_32bit, m_sumSignGY_GX_pixel_32bit);
  // sumDIX and sumDIY left shift by 2 is calculated in previous step
  const int bioSubblockSize = width * height;
  for (int pixel_index = 0; pixel_index < bioSubblockSize; pixel_index++)
  {
    if (m_sumAbsGX_pixel_32bit[pixel_index] == 0)
    {
      m_sumDIX_pixel_32bit[pixel_index] = 0;
      m_sumAbsGX_pixel_32bit[pixel_index] = 32;
    }
    else
    {
      m_sumAbsGX_pixel_32bit[pixel_index] = floorLog2(m_sumAbsGX_pixel_32bit[pixel_index]);
    }
    if (m_sumAbsGY_pixel_32bit[pixel_index] == 0)
    {
      m_sumDIY_pixel_32bit[pixel_index] = 0;
      m_sumSignGY_GX_pixel_32bit[pixel_index] = 0;
      m_sumAbsGY_pixel_32bit[pixel_index] = 32;
    }
    else
    {
      m_sumAbsGY_pixel_32bit[pixel_index] = floorLog2(m_sumAbsGY_pixel_32bit[pixel_index]);
    }
  }
  g_pelBufOP.calcBIOClippedVxVy(m_sumDIX_pixel_32bit, m_sumAbsGX_pixel_32bit, m_sumDIY_pixel_32bit, m_sumAbsGY_pixel_32bit, m_sumSignGY_GX_pixel_32bit, limit, bioSubblockSize, m_tmpx_pixel_32bit, m_tmpy_pixel_32bit);
  bioParamOffset += ((bioParamStride + 1) << 1);
#else
  bioParamOffset += ((bioParamStride + 1) << 1);
  int unitSize = 4, extendSize = 1; // unitSize = 1, extendSize = 2 gives same results as per-pixel BDOF
  for (int yUnit = 0; yUnit < height; yUnit += unitSize)
  {
    for (int xUnit = 0; xUnit < width; xUnit += unitSize)
    {
      int subTmpx = 0, subTmpy = 0;
      int subSumGx = 0, subSumGy = 0, subSumDIX = 0, subSumDIY = 0, subSumSignGY_GX = 0;
      int subBioParamOffset = bioParamOffset + (yUnit - extendSize) * bioParamStride + xUnit;
      for (int ySub = -extendSize; ySub < (extendSize + unitSize); ySub++)
      {
        for (int xSub = -extendSize; xSub < (extendSize + unitSize); xSub++)
        {
          subSumGx += m_absGx[subBioParamOffset + xSub];
          subSumGy += m_absGy[subBioParamOffset + xSub];
          subSumDIX += m_dIx[subBioParamOffset + xSub];
          subSumDIY += m_dIy[subBioParamOffset + xSub];
          subSumSignGY_GX += m_signGxGy[subBioParamOffset + xSub];
        }
        subBioParamOffset += bioParamStride;
      }
      subTmpx = (subSumGx == 0 ? 0 : rightShiftMSB(subSumDIX << 2, subSumGx));
      subTmpx = Clip3(-limit, limit, subTmpx);

      int     mainsGxGy = subSumSignGY_GX >> 12;
      int     secsGxGy = subSumSignGY_GX & ((1 << 12) - 1);
      int     tmpData = subTmpx * mainsGxGy;
      tmpData = ((tmpData << 12) + subTmpx*secsGxGy) >> 1;
      subTmpy = (subSumGy == 0 ? 0 : rightShiftMSB(((subSumDIY << 2) - tmpData), subSumGy));
      subTmpy = Clip3(-limit, limit, subTmpy);
      int curSubIdx = yUnit * width + xUnit;
      for (int ySub = 0; ySub < unitSize; ySub++)
      {
        for (int xSub = 0; xSub < unitSize; xSub++)
        {
          m_tmpx_pixel_32bit[curSubIdx + xSub] = subTmpx;
          m_tmpy_pixel_32bit[curSubIdx + xSub] = subTmpy;
        }
        curSubIdx += width;
      }
    }
  }
#endif
#if JVET_Z0136_OOB
  g_pelBufOP.addBIOAvgN(src0, src0Stride, src1, src1Stride, dstY, dstStride, m_gradX0 + bioParamOffset, m_gradX1 + bioParamOffset, m_gradY0 + bioParamOffset, m_gradY1 + bioParamOffset, bioParamStride, width, height, m_tmpx_pixel_32bit, m_tmpy_pixel_32bit, shiftNum, offset, clpRng, mcMask, mcStride, isOOB);
#else
  g_pelBufOP.addBIOAvgN(src0, src0Stride, src1, src1Stride, dstY, dstStride, m_gradX0 + bioParamOffset, m_gradX1 + bioParamOffset, m_gradY0 + bioParamOffset, m_gradY1 + bioParamOffset, bioParamStride, width, height, m_tmpx_pixel_32bit, m_tmpy_pixel_32bit, shiftNum, offset, clpRng);
#endif
}
#endif

void InterPrediction::xAddBIOAvg4(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, Pel *dst, int dstStride, const Pel *gradX0, const Pel *gradX1, const Pel *gradY0, const Pel*gradY1, int gradStride, int width, int height, int tmpx, int tmpy, int shift, int offset, const ClpRng& clpRng)
{
  g_pelBufOP.addBIOAvg4(src0, src0Stride, src1, src1Stride, dst, dstStride, gradX0, gradX1, gradY0, gradY1, gradStride, width, height, tmpx, tmpy, shift, offset, clpRng);
}

void InterPrediction::xBioGradFilter(Pel* pSrc, int srcStride, int width, int height, int gradStride, Pel* gradX, Pel* gradY, int bitDepth)
{
  g_pelBufOP.bioGradFilter(pSrc, srcStride, width, height, gradStride, gradX, gradY, bitDepth);
}

void InterPrediction::xCalcBIOPar(const Pel* srcY0Temp, const Pel* srcY1Temp, const Pel* gradX0, const Pel* gradX1, const Pel* gradY0, const Pel* gradY1, int* dotProductTemp1, int* dotProductTemp2, int* dotProductTemp3, int* dotProductTemp5, int* dotProductTemp6, const int src0Stride, const int src1Stride, const int gradStride, const int widthG, const int heightG, int bitDepth)
{
  g_pelBufOP.calcBIOPar(srcY0Temp, srcY1Temp, gradX0, gradX1, gradY0, gradY1, dotProductTemp1, dotProductTemp2, dotProductTemp3, dotProductTemp5, dotProductTemp6, src0Stride, src1Stride, gradStride, widthG, heightG, bitDepth);
}

void InterPrediction::xCalcBlkGradient(int sx, int sy, int    *arraysGx2, int     *arraysGxGy, int     *arraysGxdI, int     *arraysGy2, int     *arraysGydI, int     &sGx2, int     &sGy2, int     &sGxGy, int     &sGxdI, int     &sGydI, int width, int height, int unitSize)
{
  g_pelBufOP.calcBlkGradient(sx, sy, arraysGx2, arraysGxGy, arraysGxdI, arraysGy2, arraysGydI, sGx2, sGy2, sGxGy, sGxdI, sGydI, width, height, unitSize);
}

#if MULTI_PASS_DMVR
void InterPrediction::xWeightedAverage(const bool isBdofMvRefine, const int bdofBlockOffset,
#else
void InterPrediction::xWeightedAverage(
#endif
    const PredictionUnit& pu, const CPelUnitBuf& pcYuvSrc0, const CPelUnitBuf& pcYuvSrc1, PelUnitBuf& pcYuvDst, const BitDepths& clipBitDepths,
#if JVET_Z0136_OOB
    const ClpRngs& clpRngs, const bool& bioApplied, bool lumaOnly, bool chromaOnly, PelUnitBuf* yuvDstTmp /*= NULL*/, bool *mcMask[2], int mcStride, bool *mcMaskChroma[2], int mcCStride, bool *isOOB)
#else
    const ClpRngs& clpRngs, const bool& bioApplied, bool lumaOnly, bool chromaOnly, PelUnitBuf* yuvDstTmp /*= NULL*/)
#endif
{
  CHECK( (chromaOnly && lumaOnly), "should not happen" );

  const int iRefIdx0 = pu.refIdx[0];
  const int iRefIdx1 = pu.refIdx[1];

  if( iRefIdx0 >= 0 && iRefIdx1 >= 0 )
  {
    if( pu.cu->BcwIdx != BCW_DEFAULT && (yuvDstTmp || !pu.ciipFlag) )
    {
      CHECK(bioApplied, "Bcw is disallowed with BIO");
#if JVET_Z0136_OOB
      pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, pu.cu->BcwIdx, chromaOnly, lumaOnly, mcMask, mcStride, mcMaskChroma, mcCStride, isOOB);
#else
      pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, pu.cu->BcwIdx, chromaOnly, lumaOnly);
#endif
#if JVET_Z0136_OOB
      if (yuvDstTmp)
      {
        yuvDstTmp->addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, BCW_DEFAULT, chromaOnly, lumaOnly, mcMask, mcStride, mcMaskChroma, mcCStride, isOOB);
      }
#else
      if (yuvDstTmp)
        yuvDstTmp->addAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, chromaOnly, lumaOnly);
#endif
      return;
    }
    if (bioApplied)
    {
#if !JVET_Z0136_OOB
      const int  src0Stride = pu.lwidth() + 2 * BIO_EXTEND_SIZE + 2;
      const int  src1Stride = pu.lwidth() + 2 * BIO_EXTEND_SIZE + 2;
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
      const Pel* pSrcY0 = m_filteredBlockTmp[2][COMPONENT_Y] + (1 + BIO_EXTEND_SIZE) * (src0Stride + 1);
      const Pel* pSrcY1 = m_filteredBlockTmp[3][COMPONENT_Y] + (1 + BIO_EXTEND_SIZE) * (src1Stride + 1);
#else
      const Pel* pSrcY0 = m_filteredBlockTmp[2][COMPONENT_Y] + 2 * src0Stride + 2;
      const Pel* pSrcY1 = m_filteredBlockTmp[3][COMPONENT_Y] + 2 * src1Stride + 2;
#endif
#endif

      bool bioEnabled = true;
      if (bioEnabled)
      {
#if MULTI_PASS_DMVR
#if JVET_Z0136_OOB
        applyBiOptFlow(isBdofMvRefine, bdofBlockOffset, pu, pcYuvSrc0, pcYuvSrc1, iRefIdx0, iRefIdx1, pcYuvDst, clipBitDepths, mcMask, mcMaskChroma, isOOB);
#else
        applyBiOptFlow(isBdofMvRefine, bdofBlockOffset, pu, pcYuvSrc0, pcYuvSrc1, iRefIdx0, iRefIdx1, pcYuvDst, clipBitDepths);
#endif
#else
#if JVET_Z0136_OOB
        applyBiOptFlow(pu, pcYuvSrc0, pcYuvSrc1, iRefIdx0, iRefIdx1, pcYuvDst, clipBitDepths, mcMask, mcMaskChroma, isOOB);
#else
        applyBiOptFlow(pu, pcYuvSrc0, pcYuvSrc1, iRefIdx0, iRefIdx1, pcYuvDst, clipBitDepths);
#endif
#endif
#if JVET_Z0136_OOB
        if (yuvDstTmp)
        {
          yuvDstTmp->addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, pu.cu->BcwIdx, false, true, mcMask, mcStride, mcMaskChroma, mcCStride, isOOB);
        }
#else
        if (yuvDstTmp)
          yuvDstTmp->bufs[0].addAvg(CPelBuf(pSrcY0, src0Stride, pu.lumaSize()), CPelBuf(pSrcY1, src1Stride, pu.lumaSize()), clpRngs.comp[0]);
#endif
      }
      else
      {
#if JVET_Z0136_OOB
        pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, BCW_DEFAULT, chromaOnly, lumaOnly, mcMask, mcStride, mcMaskChroma, mcCStride, isOOB);
        if (yuvDstTmp)
        {
          yuvDstTmp->bufs[0].copyFrom(pcYuvDst.bufs[0]);
        }
#else
        pcYuvDst.bufs[0].addAvg(CPelBuf(pSrcY0, src0Stride, pu.lumaSize()), CPelBuf(pSrcY1, src1Stride, pu.lumaSize()), clpRngs.comp[0]);
        if (yuvDstTmp)
          yuvDstTmp->bufs[0].copyFrom(pcYuvDst.bufs[0]);
#endif
      }
    }
    if (!bioApplied && (lumaOnly || chromaOnly))
    {
#if JVET_Z0136_OOB
      pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, BCW_DEFAULT, chromaOnly, lumaOnly, mcMask, mcStride, mcMaskChroma, mcCStride, isOOB);
#else
      pcYuvDst.addAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, chromaOnly, lumaOnly);
#endif
    }
#if MULTI_PASS_DMVR
      // this part is to derive the chroma dst pred
    else if (!isBdofMvRefine || !bioApplied || yuvDstTmp != NULL)
#else
    else
#endif
    {
#if JVET_Z0136_OOB
      if (bioApplied)
      {
        pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, BCW_DEFAULT, true, false, mcMask, mcStride, mcMaskChroma, mcCStride, isOOB);
      }
      else
      {
        pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, BCW_DEFAULT, chromaOnly, lumaOnly, mcMask, mcStride, mcMaskChroma, mcCStride, isOOB);
      }
#else
      pcYuvDst.addAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, bioApplied);
#endif
    }
    if (yuvDstTmp)
    {
      if (bioApplied)
      {
        if (isChromaEnabled(yuvDstTmp->chromaFormat))
        {
          yuvDstTmp->bufs[1].copyFrom(pcYuvDst.bufs[1]);
          yuvDstTmp->bufs[2].copyFrom(pcYuvDst.bufs[2]);
        }
      }
      else
      {
        yuvDstTmp->copyFrom(pcYuvDst, lumaOnly, chromaOnly);
      }
    }
  }
  else if( iRefIdx0 >= 0 && iRefIdx1 < 0 )
  {
    if( pu.cu->geoFlag )
    {
#if JVET_W0097_GPM_MMVD_TM
      pcYuvDst.copyFrom(pcYuvSrc0, lumaOnly, chromaOnly);
#else
      pcYuvDst.copyFrom( pcYuvSrc0 );
#endif
    }
    else
    {
      pcYuvDst.copyClip( pcYuvSrc0, clpRngs, lumaOnly, chromaOnly );
    }
    if (yuvDstTmp)
    {
      yuvDstTmp->copyFrom( pcYuvDst, lumaOnly, chromaOnly );
    }
  }
  else if( iRefIdx0 < 0 && iRefIdx1 >= 0 )
  {
    if( pu.cu->geoFlag )
    {
#if JVET_W0097_GPM_MMVD_TM
      pcYuvDst.copyFrom(pcYuvSrc1, lumaOnly, chromaOnly);
#else
      pcYuvDst.copyFrom( pcYuvSrc1 );
#endif
    }
    else
    {
      pcYuvDst.copyClip( pcYuvSrc1, clpRngs, lumaOnly, chromaOnly );
    }
    if (yuvDstTmp)
    {
      yuvDstTmp->copyFrom(pcYuvDst, lumaOnly, chromaOnly);
    }
  }
}

#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
#if !INTER_LIC
template <bool trueAfalseL>
void InterPrediction::xGetPredBlkTpl(const CodingUnit& cu, const ComponentID compID, const CPelBuf& refBuf, const Mv& mv, const int posW, const int posH, const int tplSize, Pel* predBlkTpl
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
                                     , bool AML
#endif
                                     )
{
  const int lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
  const int horShift = (lumaShift + ::getComponentScaleX(compID, cu.chromaFormat));
  const int verShift = (lumaShift + ::getComponentScaleY(compID, cu.chromaFormat));

  const int xInt = mv.getHor() >> horShift;
  const int yInt = mv.getVer() >> verShift;
  const int xFrac = mv.getHor() & ((1 << horShift) - 1);
  const int yFrac = mv.getVer() & ((1 << verShift) - 1);

  const Pel* ref;
  Pel* dst;
  int refStride, dstStride, bw, bh;
  if (trueAfalseL)
  {
    ref = refBuf.bufAt(cu.blocks[compID].pos().offset(xInt + posW, yInt + posH - 1));
    dst = predBlkTpl + posW;
    refStride = refBuf.stride;
    dstStride = tplSize;
    bw = tplSize;
    bh = 1;
  }
  else
  {
    ref = refBuf.bufAt(cu.blocks[compID].pos().offset(xInt + posW - 1, yInt + posH));
    dst = predBlkTpl + posH;
    refStride = refBuf.stride;
    dstStride = 1;
    bw = 1;
    bh = tplSize;
  }

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  int   nFilterIdx =  AML  ? 1 : 0;
#else
  const int  nFilterIdx   = 0;
#endif
  const bool useAltHpelIf = false;

  if (yFrac == 0)
  {
    m_if.filterHor(compID, (Pel*)ref, refStride, dst, dstStride, bw, bh, xFrac, true, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
  }
  else if (xFrac == 0)
  {
    m_if.filterVer(compID, (Pel*)ref, refStride, dst, dstStride, bw, bh, yFrac, true, true, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
  }
  else
  {
      
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if IF_12TAP
    int vFilterSize = isLuma(compID) ? NTAPS_LUMA(0) : NTAPS_CHROMA;
#else
    int vFilterSize = isLuma(compID) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif
    if (isLuma(compID) && nFilterIdx == 1)
    {
      vFilterSize = NTAPS_BILINEAR;
    }
#else
#if IF_12TAP
    const int vFilterSize = isLuma(compID) ? NTAPS_LUMA(0) : NTAPS_CHROMA;
#else
    const int vFilterSize = isLuma(compID) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif
#endif
    PelBuf tmpBuf = PelBuf(m_filteredBlockTmp[0][compID], Size(bw, bh + vFilterSize - 1));

    m_if.filterHor(compID, (Pel*)ref - ((vFilterSize >> 1) - 1)*refStride, refStride, tmpBuf.buf, tmpBuf.stride, bw, bh + vFilterSize - 1, xFrac, false, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
    JVET_J0090_SET_CACHE_ENABLE(false);
    m_if.filterVer(compID, tmpBuf.buf + ((vFilterSize >> 1) - 1)*tmpBuf.stride, tmpBuf.stride, dst, dstStride, bw, bh, yFrac, false, true, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
    JVET_J0090_SET_CACHE_ENABLE(true);
  }
}
#endif
void InterPrediction::xWeightedAverageY(const PredictionUnit& pu, const CPelUnitBuf& pcYuvSrc0, const CPelUnitBuf& pcYuvSrc1, PelUnitBuf& pcYuvDst, const BitDepths& clipBitDepths, const ClpRngs& clpRngs)
{
  const int iRefIdx0 = pu.refIdx[0];
  const int iRefIdx1 = pu.refIdx[1];

  if (iRefIdx0 >= 0 && iRefIdx1 >= 0)
  {
    if (pu.cu->BcwIdx != BCW_DEFAULT)
    {
#if JVET_Z0136_OOB
      bool isOOB[2] = { false,false };
      pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, pu.cu->BcwIdx, false, true, pu.cs->mcMask, -1, pu.cs->mcMaskChroma, -1, isOOB);
#else
      pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, pu.cu->BcwIdx, false, true);
#endif
    }
    else
    {
      pcYuvDst.addAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, false, true);
    }
  }
  else if (iRefIdx0 >= 0 && iRefIdx1 < 0)
  {
    pcYuvDst.copyClip(pcYuvSrc0, clpRngs, true);
  }
  else if (iRefIdx0 < 0 && iRefIdx1 >= 0)
  {
    pcYuvDst.copyClip(pcYuvSrc1, clpRngs, true);
  }
}
#endif

#if JVET_W0090_ARMC_TM
void InterPrediction::xPredAffineTpl(const PredictionUnit &pu, const RefPicList &eRefPicList, int* numTemplate, Pel* refLeftTemplate, Pel* refAboveTemplate)
{
  int iRefIdx = pu.refIdx[eRefPicList];
  CHECK(iRefIdx < 0, "iRefIdx incorrect.");
  const Picture* refPic = pu.cu->slice->getRefPic(eRefPicList, iRefIdx)->unscaledPic;
  Mv mvLT = pu.mvAffi[eRefPicList][0];
  Mv mvRT = pu.mvAffi[eRefPicList][1];
  Mv mvLB = pu.mvAffi[eRefPicList][2];
  // get affine sub-block width and height
  const int width = pu.Y().width;
  const int height = pu.Y().height;
  int blockWidth = AFFINE_MIN_BLOCK_SIZE;
  int blockHeight = AFFINE_MIN_BLOCK_SIZE;

  CHECK(blockWidth > width, "Sub Block width  > Block width");
  CHECK(blockHeight > height, "Sub Block height > Block height");

  const int cxWidth = width;
  const int cxHeight = height;
  const int iBit = MAX_CU_DEPTH;
  int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;
  iDMvHorX = (mvRT - mvLT).getHor() << (iBit - floorLog2(width));
  iDMvHorY = (mvRT - mvLT).getVer() << (iBit - floorLog2(width));
  if (pu.cu->affineType == AFFINEMODEL_6PARAM)
  {
    iDMvVerX = (mvLB - mvLT).getHor() << (iBit - floorLog2(height));
    iDMvVerY = (mvLB - mvLT).getVer() << (iBit - floorLog2(height));
  }
  else
  {
    iDMvVerX = -iDMvHorY;
    iDMvVerY = iDMvHorX;
  }
  int iMvScaleHor = mvLT.getHor() << iBit;
  int iMvScaleVer = mvLT.getVer() << iBit;

  const int shift = iBit - 4 + MV_FRACTIONAL_BITS_INTERNAL;
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  const bool subblkMVSpreadOverLimit = isSubblockVectorSpreadOverLimit(iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY, pu.interDir);
#endif
#if AFFINE_RM_CONSTRAINTS_AND_OPT
  if (iDMvHorX == 0 && iDMvHorY == 0)
    blockWidth = width;
  else
  {
    int maxDmv = std::max(abs(iDMvHorX), abs(iDMvHorY)) * blockWidth;
    int TH = 1 << (iBit - 1); // Half pel
    while (maxDmv < TH && blockWidth < width)
    {
      blockWidth <<= 1;
      maxDmv <<= 1;
    }
  }
  if (iDMvVerX == 0 && iDMvVerY == 0)
    blockHeight = height;
  else
  {
    int maxDmv = std::max(abs(iDMvVerX), abs(iDMvVerY)) * blockHeight;
    int TH = 1 << (iBit - 1); // Half pel
    while (maxDmv < TH && blockHeight < height)
    {
      blockHeight <<= 1;
      maxDmv <<= 1;
    }
  }
#endif
  int iMvScaleTmpHor0 = iMvScaleHor + ((iDMvHorX * blockWidth + iDMvVerX * blockHeight) >> 1);
  int iMvScaleTmpVer0 = iMvScaleVer + ((iDMvHorY * blockWidth + iDMvVerY * blockHeight) >> 1);

#if JVET_Z0139_NA_AFF
  const CodingUnit *const cuAbove = pu.cu->cs->getCU(pu.cu->blocks[COMPONENT_Y].pos().offset(0, -1), toChannelType(COMPONENT_Y));
  const CodingUnit *const cuLeft  = pu.cu->cs->getCU(pu.cu->blocks[COMPONENT_Y].pos().offset(-1, 0), toChannelType(COMPONENT_Y));

  // get prediction block by block
  for (int h = 0; (cuLeft && h < cxHeight) || h < 1; h += blockHeight)
  {
    for (int w = 0; (cuAbove && w < cxWidth) || w < 1; w += blockWidth)
#else
  for (int h = 0; h < cxHeight; h += blockHeight)
  {
    for (int w = 0; w < cxWidth; w += blockWidth)
#endif
    {
      if (w == 0 || h == 0)
      {
        int iMvScaleTmpHor, iMvScaleTmpVer;

#if !AFFINE_RM_CONSTRAINTS_AND_OPT
        if (!subblkMVSpreadOverLimit)
#endif
        {
          iMvScaleTmpHor = iMvScaleTmpHor0 + iDMvHorX * w + iDMvVerX * h;
          iMvScaleTmpVer = iMvScaleTmpVer0 + iDMvHorY * w + iDMvVerY * h;
        }
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
        else
        {
          iMvScaleTmpHor = iMvScaleHor + iDMvHorX * (cxWidth >> 1) + iDMvVerX * (cxHeight >> 1);
          iMvScaleTmpVer = iMvScaleVer + iDMvHorY * (cxWidth >> 1) + iDMvVerY * (cxHeight >> 1);
        }
#endif
        roundAffineMv(iMvScaleTmpHor, iMvScaleTmpVer, shift);
        Mv tmpMv(iMvScaleTmpHor, iMvScaleTmpVer);
        tmpMv.clipToStorageBitDepth();
        iMvScaleTmpHor = tmpMv.getHor();
        iMvScaleTmpVer = tmpMv.getVer();

        // clip and scale
        if (refPic->isRefScaled(pu.cs->pps) == false)
        {
          clipMv(tmpMv, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);
          iMvScaleTmpHor = tmpMv.getHor();
          iMvScaleTmpVer = tmpMv.getVer();
        }
        xGetSublkAMLTemplate(*pu.cu, COMPONENT_Y, *refPic, Mv(iMvScaleTmpHor, iMvScaleTmpVer), blockWidth, blockHeight, w, h, numTemplate, refLeftTemplate, refAboveTemplate
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
                             , pu.afMmvdFlag
#endif
                             );
      }
    }
  }
}
#endif

void InterPrediction::motionCompensation( PredictionUnit &pu, PelUnitBuf &predBuf, const RefPicList &eRefPicList
  , const bool luma, const bool chroma
  , PelUnitBuf* predBufWOBIO /*= NULL*/
)
{
  // Note: there appears to be an interaction with weighted prediction that
  // makes the code follow different paths if chroma is on or off (in the encoder).
  // Therefore for 4:0:0, "chroma" is not changed to false.
#if MULTI_HYP_PRED
  if (!pu.addHypData.empty())
  {
    CHECK(eRefPicList != REF_PIC_LIST_X, "Multi Hyp: eRefPicList != REF_PIC_LIST_X");
    CHECK(!luma, "Multi Hyp: !luma");
    xAddHypMC(pu, predBuf, predBufWOBIO, !chroma);
    return;
  }
#endif
  CHECK(predBufWOBIO && pu.ciipFlag, "the case should not happen!");

  if (!pu.cs->pcv->isEncoder)
  {
    if (CU::isIBC(*pu.cu))
    {
      CHECK(!luma, "IBC only for Chroma is not allowed.");
      xIntraBlockCopy(pu, predBuf, COMPONENT_Y);
      if (chroma && isChromaEnabled(pu.chromaFormat))
      {
        xIntraBlockCopy(pu, predBuf, COMPONENT_Cb);
        xIntraBlockCopy(pu, predBuf, COMPONENT_Cr);
      }
      return;
    }
  }
  // dual tree handling for IBC as the only ref
  if ((!luma || !chroma) && eRefPicList == REF_PIC_LIST_0)
  {
    xPredInterUni(pu, eRefPicList, predBuf, false, false, luma, chroma);
    return;
  }
  // else, go with regular MC below
        CodingStructure &cs = *pu.cs;
  const PPS &pps            = *cs.pps;
  const SliceType sliceType =  cs.slice->getSliceType();

  if( eRefPicList != REF_PIC_LIST_X )
  {
    CHECK(predBufWOBIO != NULL, "the case should not happen!");
    if ((CU::isIBC(*pu.cu) == false) && ((sliceType == P_SLICE && pps.getUseWP()) || (sliceType == B_SLICE && pps.getWPBiPred()))
#if INTER_LIC
      && !pu.cu->LICFlag
#endif
      )
    {
      xPredInterUni(pu, eRefPicList, predBuf, true, false, luma, chroma);
      xWeightedPredictionUni(pu, predBuf, eRefPicList, predBuf, -1, m_maxCompIDToPred, (luma && !chroma),
                             (!luma && chroma));
    }
    else
    {
      xPredInterUni(pu, eRefPicList, predBuf, false, false, luma, chroma);
    }
  }
  else
  {
#if !INTER_RM_SIZE_CONSTRAINTS
#if ENABLE_OBMC
  if (pu.cu->isobmcMC == false)
#endif
    CHECK( !pu.cu->affine && pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0 && ( pu.lwidth() + pu.lheight() == 12 ), "invalid 4x8/8x4 bi-predicted blocks" );
#endif
#if !BDOF_RM_CONSTRAINTS
    int refIdx0 = pu.refIdx[REF_PIC_LIST_0];
    int refIdx1 = pu.refIdx[REF_PIC_LIST_1];

    const WPScalingParam *wp0 = pu.cs->slice->getWpScaling(REF_PIC_LIST_0, refIdx0);
    const WPScalingParam *wp1 = pu.cs->slice->getWpScaling(REF_PIC_LIST_1, refIdx1);

    bool bioApplied = false;
    const Slice &slice = *pu.cs->slice;
    if (pu.cs->sps->getBDOFEnabledFlag() && (!pu.cs->picHeader->getDisBdofFlag()))
    {
      if (pu.cu->affine || m_subPuMC)
      {
        bioApplied = false;
      }
      else
      {
        const bool biocheck0 =
          !((WPScalingParam::isWeighted(wp0) || WPScalingParam::isWeighted(wp1)) && slice.getSliceType() == B_SLICE);
        const bool biocheck1 = !(pps.getUseWP() && slice.getSliceType() == P_SLICE);
        if (biocheck0
          && biocheck1
          && PU::isBiPredFromDifferentDirEqDistPoc(pu)
          && (pu.Y().height >= 8)
          && (pu.Y().width >= 8)
          && ((pu.Y().height * pu.Y().width) >= 128)
          )
        {
          bioApplied = true;
        }
      }

      if (bioApplied && pu.ciipFlag)
      {
        bioApplied = false;
      }

      if (bioApplied && pu.cu->smvdMode)
      {
        bioApplied = false;
      }
      if (pu.cu->cs->sps->getUseBcw() && bioApplied && pu.cu->BcwIdx != BCW_DEFAULT)
      {
        bioApplied = false;
      }
      if (pu.mmvdEncOptMode == 2 && pu.mmvdMergeFlag)
      {
        bioApplied = false;
      }
    }
#if ENABLE_OBMC
      if (pu.cu->isobmcMC)
      {
        bioApplied = false;
      }
#endif
    bool refIsScaled = ( refIdx0 < 0 ? false : pu.cu->slice->getRefPic( REF_PIC_LIST_0, refIdx0 )->isRefScaled( pu.cs->pps ) ) ||
                       ( refIdx1 < 0 ? false : pu.cu->slice->getRefPic( REF_PIC_LIST_1, refIdx1 )->isRefScaled( pu.cs->pps ) );
    bioApplied = refIsScaled ? false : bioApplied;

    bool dmvrApplied = false;
    dmvrApplied = (pu.mvRefine) && PU::checkDMVRCondition(pu);
#if MULTI_PASS_DMVR
    if ((pu.lumaSize().width > MAX_BDOF_APPLICATION_REGION || pu.lumaSize().height > MAX_BDOF_APPLICATION_REGION) && pu.mergeType != MRG_TYPE_SUBPU_ATMVP && (bioApplied && !dmvrApplied && !pu.bdmvrRefine))
#else
    if ((pu.lumaSize().width > MAX_BDOF_APPLICATION_REGION || pu.lumaSize().height > MAX_BDOF_APPLICATION_REGION) && pu.mergeType != MRG_TYPE_SUBPU_ATMVP && (bioApplied && !dmvrApplied))
#endif
    {
      xSubPuBio(pu, predBuf, eRefPicList, predBufWOBIO);
    }
    else
#endif
      if( pu.mergeType != MRG_TYPE_DEFAULT_N && pu.mergeType != MRG_TYPE_IBC )
      {
        CHECK( predBufWOBIO != NULL, "the case should not happen!" );
        xSubPuMC( pu, predBuf, eRefPicList, luma, chroma );
      }
      else if( xCheckIdenticalMotion( pu ) )
      {
        xPredInterUni( pu, REF_PIC_LIST_0, predBuf, false, false, luma, chroma );
        if( predBufWOBIO )
        {
          predBufWOBIO->copyFrom( predBuf, ( luma && !chroma ), ( chroma && !luma ) );
        }
      }
      else
      {
#if MULTI_PASS_DMVR
        m_bdofMvRefined = false;
#if !BDOF_RM_CONSTRAINTS
        if (pu.bdmvrRefine && !bioApplied)
        {
          for (int bdofSubPuIdx = 0; bdofSubPuIdx < BDOF_SUBPU_MAX_NUM; bdofSubPuIdx++)
          {
            m_bdofSubPuMvOffset[bdofSubPuIdx].setZero();
          }
        }
#endif
        xPredInterBi(pu, predBuf, luma, chroma, predBufWOBIO);
        if (m_bdofMvRefined)
        {
          xPredInterBiSubPuBDOF(pu, predBuf, luma, chroma);  // do not change the predBufWOBIO
          m_bdofMvRefined = false;
        }
#else
        xPredInterBi( pu, predBuf, luma, chroma, predBufWOBIO );
#endif
      }
  }
  return;
}

void InterPrediction::motionCompensation( CodingUnit &cu, const RefPicList &eRefPicList
  , const bool luma, const bool chroma
)
{
  for( auto &pu : CU::traversePUs( cu ) )
  {
    PelUnitBuf predBuf = cu.cs->getPredBuf( pu );
    pu.mvRefine = true;
    motionCompensation(pu, predBuf, eRefPicList, luma, chroma);
    pu.mvRefine = false;
  }
}

void InterPrediction::motionCompensation( PredictionUnit &pu, const RefPicList &eRefPicList /*= REF_PIC_LIST_X*/
  , const bool luma, const bool chroma
)
{
  PelUnitBuf predBuf = pu.cs->getPredBuf( pu );
  motionCompensation(pu, predBuf, eRefPicList, luma, chroma);
}

#if ENABLE_OBMC
/** Function for sub-block based Overlapped Block Motion Compensation (OBMC).
*
* This function can:
* 1. Perform sub-block OBMC for a CU.
* 2. Before motion estimation, subtract (scaled) predictors generated by applying neighboring motions to current CU/PU from the original signal of current CU/PU,
*    to make the motion estimation biased to OBMC.
*/
void InterPrediction::subBlockOBMC(PredictionUnit  &pu, PelUnitBuf* pDst)
{
  if (
    pu.cs->sps->getUseOBMC() == false
    || pu.cu->obmcFlag == false
#if INTER_LIC
    || pu.cu->LICFlag
#endif
    || pu.lwidth() * pu.lheight() < 32
    )
  {
    return;
  }

  const UnitArea   orgPuArea = pu;
  PredictionUnit subPu = pu;

  const uint32_t uiWidth = pu.lwidth();
  const uint32_t uiHeight = pu.lheight();

  const uint32_t uiMinCUW = pu.cs->pcv->minCUWidth;

  const uint32_t uiHeightInBlock = uiHeight / uiMinCUW;
  const uint32_t uiWidthInBlock = uiWidth / uiMinCUW;

#if MULTI_PASS_DMVR
  const bool bSubMotion = pu.cu->affine || pu.bdmvrRefine;
#else
  const bool bSubMotion = pu.cu->affine || PU::checkDMVRCondition(pu);
#endif

  MotionInfo NeighMi = MotionInfo();

  int BcwIdx = pu.cu->BcwIdx;
  bool affine = pu.cu->affine;
  bool geo = pu.cu->geoFlag;
  subPu.cu->affine = false;
  subPu.cu->BcwIdx = BCW_DEFAULT;
  subPu.cu->geoFlag = false;
#if INTER_LIC
  subPu.cu->LICFlag = false;
#endif
  subPu.ciipFlag = false;
#if TM_MRG
  subPu.tmMergeFlag = false;
#endif
#if MULTI_PASS_DMVR
  subPu.bdmvrRefine = false;
#endif
  subPu.mvRefine = false;
  subPu.mmvdMergeFlag = false;
  PelUnitBuf pcYuvPred = pDst == nullptr ? pu.cs->getPredBuf(pu) : *pDst;

  PelUnitBuf pcYuvTmpPredL0 = m_tmpObmcBufL0.subBuf(UnitAreaRelative(*pu.cu, pu));
  PelUnitBuf pcYuvTmpPredT0 = m_tmpObmcBufT0.subBuf(UnitAreaRelative(*pu.cu, pu));

  for (int iBlkBoundary = 0; iBlkBoundary < 2; iBlkBoundary++)  // 0 - top; 1 - left
  {
    unsigned int uiLengthInBlock = ((iBlkBoundary == 0) ? uiWidthInBlock : uiHeightInBlock);

    int iSub = 0, iState = 0;

    while (iSub < uiLengthInBlock)
    {
      int      iLength = 0;
      Position curOffset = (iBlkBoundary == 0) ? Position(iSub * uiMinCUW, 0) : Position(0, iSub * uiMinCUW);

      iState = PU::getSameNeigMotion(pu, NeighMi, curOffset, iBlkBoundary, iLength, uiLengthInBlock - iSub);

      if (iState == 2)  // do OBMC
      {
        subPu = NeighMi;
        if (iBlkBoundary == 0)
        {
          subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(orgPuArea.lumaPos().offset(iSub * uiMinCUW, 0), Size{ iLength*uiMinCUW, uiMinCUW })));
        }
        else
        {
          subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(orgPuArea.lumaPos().offset(0, iSub * uiMinCUW), Size{ uiMinCUW, iLength*uiMinCUW })));
        }

        const UnitArea predArea = UnitAreaRelative(orgPuArea, subPu);
        PelUnitBuf     cPred = pcYuvPred.subBuf(predArea);
        PelUnitBuf cTmp1;

        if (iBlkBoundary == 0)//above
        {
          cTmp1 = pcYuvTmpPredT0.subBuf(predArea);
        }
        else//left
        {
          cTmp1 = pcYuvTmpPredL0.subBuf(predArea);
        }
#if JVET_Z0061_TM_OBMC
        bool isAbove   = (iBlkBoundary == 0) ? 1 : 0;
        int  iOBMCmode = selectOBMCmode(pu, subPu, isAbove, iLength, uiMinCUW, curOffset);

        if (iOBMCmode == 1)   // 1: current;
        {
          iSub += iLength;
        }
        else if (iOBMCmode == 2)   // 2: neighbour;
        {
          xSubBlockMotionCompensation(subPu, cTmp1);

          for (int compID = 0; compID < MAX_NUM_COMPONENT; compID++)
          {
            xSubblockTMOBMC(ComponentID(compID), subPu, cPred, cTmp1, iBlkBoundary, iOBMCmode);
          }
          iSub += iLength;
        }
        else   // 3: blend (OBMC) or default 0: best has not been found;
        {
          xSubBlockMotionCompensation(subPu, cTmp1);

          for (int compID = 0; compID < MAX_NUM_COMPONENT; compID++)
          {
            xSubblockTMOBMC(ComponentID(compID), subPu, cPred, cTmp1, iBlkBoundary, iOBMCmode);
          }
          iSub += iLength;
        }
#else
        xSubBlockMotionCompensation(subPu, cTmp1);

        for (int compID = 0; compID < MAX_NUM_COMPONENT; compID++)
        {
          xSubblockOBMC(ComponentID(compID), subPu, cPred, cTmp1, iBlkBoundary);
        }

        iSub += iLength;
#endif
      }
      else if (iState == 1 || iState == 3)   // consecutive intra neighbors or   skip OBMC based on MV similarity
      {
        iSub += iLength;
      }
      else                    // unavailable neighbors
      {
        iSub += uiLengthInBlock;
        break;
      }
    }
    CHECK(iSub != uiLengthInBlock, "not all sub-blocks are merged");
  }


  if (!bSubMotion)
  {
    pu.cu->BcwIdx = BcwIdx;
    pu.cu->affine = affine;
    pu.cu->geoFlag = geo;
    return;
  }

  PelUnitBuf pcYuvTmpPred = m_tmpSubObmcBuf;

  PelUnitBuf cTmp1 = pcYuvTmpPred.subBuf(UnitArea(pu.chromaFormat, Area(0, 0, uiMinCUW, uiMinCUW)));
  PelUnitBuf cTmp2 = pcYuvTmpPred.subBuf(UnitArea(pu.chromaFormat, Area(4, 0, uiMinCUW, uiMinCUW)));
  PelUnitBuf cTmp3 = pcYuvTmpPred.subBuf(UnitArea(pu.chromaFormat, Area(8, 0, uiMinCUW, uiMinCUW)));
  PelUnitBuf cTmp4 = pcYuvTmpPred.subBuf(UnitArea(pu.chromaFormat, Area(12, 0, uiMinCUW, uiMinCUW)));
  PelUnitBuf zero  = pcYuvTmpPred.subBuf(UnitArea(pu.chromaFormat, Area(16, 0, uiMinCUW, uiMinCUW)));

  for (int iSubX = 0; iSubX < uiWidthInBlock; iSubX += 1)
  {
    for (int iSubY = 0; iSubY < uiHeightInBlock; iSubY += 1)
    {
      bool bCURBoundary = (iSubX == uiWidthInBlock - 1);
      bool bCUBBoundary = (iSubY == uiHeightInBlock - 1);

      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(orgPuArea.lumaPos().offset(iSubX * uiMinCUW, iSubY * uiMinCUW), Size{ uiMinCUW, uiMinCUW })));
      const UnitArea predArea = UnitAreaRelative(orgPuArea, subPu);
      PelUnitBuf     cPred = pcYuvPred.subBuf(predArea);

      bool isAboveAvail = false, isLeftAvail = false, isBelowAvail = false, isRightAvail = false;

      // above
      if (iSubY)
      {
        isAboveAvail = PU::getNeighborMotion(pu, NeighMi, Position(iSubX * uiMinCUW, iSubY * uiMinCUW), Size(uiMinCUW, uiMinCUW), 0);
        if (isAboveAvail)
        {
          subPu = NeighMi;
          xSubBlockMotionCompensation(subPu, cTmp1);
        }
      }

      // left
      if (iSubX)
      {
        isLeftAvail = PU::getNeighborMotion(pu, NeighMi, Position(iSubX * uiMinCUW, iSubY * uiMinCUW), Size(uiMinCUW, uiMinCUW), 1);
        if (isLeftAvail)
        {
          subPu = NeighMi;
          xSubBlockMotionCompensation(subPu, cTmp2);
        }
      }

      // below
      if (!bCUBBoundary)
      {
        isBelowAvail = PU::getNeighborMotion(pu, NeighMi, Position(iSubX * uiMinCUW, iSubY * uiMinCUW), Size(uiMinCUW, uiMinCUW), 2);
        if (isBelowAvail)
        {
          subPu = NeighMi;
          xSubBlockMotionCompensation(subPu, cTmp3);
        }
      }

      // right
      if (!bCURBoundary)
      {
        isRightAvail = PU::getNeighborMotion(pu, NeighMi, Position(iSubX * uiMinCUW, iSubY * uiMinCUW), Size(uiMinCUW, uiMinCUW), 3);
        if (isRightAvail)
        {
          subPu = NeighMi;
          xSubBlockMotionCompensation(subPu, cTmp4);
        }
      }

      if( isAboveAvail || isLeftAvail || isBelowAvail || isRightAvail )
      {
        for( int compID = 0; compID < MAX_NUM_COMPONENT; compID++ )
        {
          xSubblockOBMCBlending( ComponentID( compID ), subPu, cPred, isAboveAvail ? cTmp1: zero, isLeftAvail ? cTmp2: zero, isBelowAvail ? cTmp3: zero, isRightAvail ? cTmp4: zero, isAboveAvail, isLeftAvail, isBelowAvail, isRightAvail, true );
        }
      }
    }
  }
  pu.cu->BcwIdx = BcwIdx;
  pu.cu->affine = affine;
  pu.cu->geoFlag = geo;

  return;
}

// Function for (weighted) averaging predictors of current block and predictors generated by applying neighboring motions to current block.
void InterPrediction::xSubblockOBMC(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst, PelUnitBuf &pcYuvPredSrc, int iDir, bool bSubMotion)
{
  int iWidth = pu.blocks[eComp].width;
  int iHeight = pu.blocks[eComp].height;

  if (iWidth == 0 || iHeight == 0)
  {
    return;
  }

  Pel* pOrgDst = pcYuvPredDst.bufs[eComp].buf;
  Pel* pOrgSrc = pcYuvPredSrc.bufs[eComp].buf;
  const int strideDst = pcYuvPredDst.bufs[eComp].stride;
  const int strideSrc = pcYuvPredSrc.bufs[eComp].stride;

  if (iDir == 0) //above
  {
    for (int i = 0; i < iWidth; i++)
    {
      Pel* pDst = pOrgDst;
      Pel* pSrc = pOrgSrc;

      pDst[i] = bSubMotion ? (3 * pDst[i] + pSrc[i] + 2) >> 2 : (26 * pDst[i] + 6 * pSrc[i] + 16) >> 5;

      if (eComp == COMPONENT_Y)
      {
        pDst += strideDst;
        pSrc += strideSrc;
        pDst[i] = (7 * pDst[i] + pSrc[i] + 4) >> 3;
        pDst += strideDst;
        pSrc += strideSrc;
        pDst[i] = (15 * pDst[i] + pSrc[i] + 8) >> 4;
        if (!bSubMotion)
        {
          pDst += strideDst;
          pSrc += strideSrc;
          pDst[i] = (31 * pDst[i] + pSrc[i] + 16) >> 5;
        }
      }
    }
  }

  if (iDir == 1) //left
  {
    Pel* pDst = pOrgDst;
    Pel* pSrc = pOrgSrc;
    for (int i = 0; i < iHeight; i++)
    {
      pDst[0] = bSubMotion ? (3 * pDst[0] + pSrc[0] + 2) >> 2 : (26 * pDst[0] + 6 * pSrc[0] + 16) >> 5;

      if (eComp == COMPONENT_Y)
      {
        pDst[1] = (7 * pDst[1] + pSrc[1] + 4) >> 3;
        pDst[2] = (15 * pDst[2] + pSrc[2] + 8) >> 4;
        if (!bSubMotion)
        {
          pDst[3] = (31 * pDst[3] + pSrc[3] + 16) >> 5;
        }
      }

      pDst += strideDst;
      pSrc += strideSrc;
    }
  }

  if (iDir == 2) //below
  {
    for (int i = 0; i < iWidth; i++)
    {
      Pel* pDst = pOrgDst + (iHeight - 1) * strideDst;
      Pel* pSrc = pOrgSrc + (iHeight - 1) * strideSrc;

      pDst[i] = (3 * pDst[i] + pSrc[i] + 2) >> 2;

      if (eComp == COMPONENT_Y)
      {
        pDst -= strideDst;
        pSrc -= strideSrc;
        pDst[i] = (7 * pDst[i] + pSrc[i] + 4) >> 3;
        pDst -= strideDst;
        pSrc -= strideSrc;
        pDst[i] = (15 * pDst[i] + pSrc[i] + 8) >> 4;
      }
    }
  }

  if (iDir == 3) //right
  {
    Pel* pDst = pOrgDst + (iWidth - 4);
    Pel* pSrc = pOrgSrc + (iWidth - 4);
    for (int i = 0; i < iHeight; i++)
    {
      pDst[3] = (3 * pDst[3] + pSrc[3] + 2) >> 2;
      if (eComp == COMPONENT_Y)
      {
        pDst[2] = (7 * pDst[2] + pSrc[2] + 4) >> 3;
        pDst[1] = (15 * pDst[1] + pSrc[1] + 8) >> 4;
      }
      pDst += strideDst;
      pSrc += strideSrc;
    }
  }
}

#if ENABLE_OBMC
void InterPrediction::xSubblockOBMCBlending(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst, PelUnitBuf &pcYuvPredSrc1, PelUnitBuf &pcYuvPredSrc2, PelUnitBuf &pcYuvPredSrc3, PelUnitBuf &pcYuvPredSrc4, bool isAboveAvail, bool isLeftAvail, bool isBelowAvail, bool isRightAvail, bool bSubMotion)
{
  int iWidth = pu.blocks[eComp].width;
  int iHeight = pu.blocks[eComp].height;

  if (iWidth == 0 || iHeight == 0)
  {
    return;
  }

  Pel* pOrgDst = pcYuvPredDst.bufs[eComp].buf;
  Pel* pOrgSrc1 = pcYuvPredSrc1.bufs[eComp].buf;
  Pel* pOrgSrc2 = pcYuvPredSrc2.bufs[eComp].buf;
  Pel* pOrgSrc3 = pcYuvPredSrc3.bufs[eComp].buf;
  Pel* pOrgSrc4 = pcYuvPredSrc4.bufs[eComp].buf;

  const int strideDst = pcYuvPredDst.bufs[eComp].stride;
  const int strideSrc = pcYuvPredSrc1.bufs[eComp].stride;

  unsigned int isChroma = !isLuma( eComp );
  unsigned int aboveWeight[4], leftWeight[4], belowWeight[4], rightWeight[4];

  if( isAboveAvail )
  {
    memcpy( aboveWeight, defaultWeight[isChroma], sizeof( aboveWeight ) );
  }
  else
  {
    memset( aboveWeight, 0, sizeof( aboveWeight ) );
  }

  if( isLeftAvail )
  {
    memcpy( leftWeight, defaultWeight[isChroma], sizeof( leftWeight ) );
  }
  else
  {
    memset( leftWeight, 0, sizeof( leftWeight ) );
  }

  if( isBelowAvail )
  {
    memcpy( belowWeight, defaultWeight[isChroma], sizeof( belowWeight ) );
  }
  else
  {
    memset( belowWeight, 0, sizeof( belowWeight ) );
  }

  if( isRightAvail )
  {
    memcpy( rightWeight, defaultWeight[isChroma], sizeof( rightWeight ) );
  }
  else
  {
    memset( rightWeight, 0, sizeof( rightWeight ) );
  }

  unsigned int shift = 7;
  unsigned int sumWeight = 1 << shift;
  unsigned int add = 1 << (shift - 1);

  Pel* pDst = pOrgDst;
  Pel* pSrc1 = pOrgSrc1;
  Pel* pSrc2 = pOrgSrc2;
  Pel* pSrc3 = pOrgSrc3;
  Pel* pSrc4 = pOrgSrc4;

  if( isLuma( eComp ) )
  {
    for( int j = 0; j < iHeight; j++ )
    {
      unsigned int idx_h = iHeight - 1 - j;

      for( int i = 0; i < iWidth; i++ )
      {
        unsigned int idx_w = iWidth - 1 - i;

        unsigned int sumOBMCWeight = aboveWeight[j] + leftWeight[i] + belowWeight[idx_h] + rightWeight[idx_w];

        if( sumOBMCWeight == 0 )
        {
          continue;
        }

        unsigned int currentWeight = sumWeight - sumOBMCWeight;

        pDst[i] = (currentWeight * pDst[i] + aboveWeight[j] * pSrc1[i] + leftWeight[i] * pSrc2[i] + belowWeight[idx_h] * pSrc3[i] + rightWeight[idx_w] * pSrc4[i] + add) >> shift;
      }
      pDst += strideDst;
      pSrc1 += strideSrc;
      pSrc2 += strideSrc;
      pSrc3 += strideSrc;
      pSrc4 += strideSrc;
    }
  }
  else
  {
    pDst[0] = ((sumWeight - aboveWeight[0] - leftWeight[0]) * pDst[0] + aboveWeight[0] * pSrc1[0] + leftWeight[0] * pSrc2[0] + add) >> shift;
    pDst[1] = ((sumWeight - aboveWeight[0] - rightWeight[0]) * pDst[1] + aboveWeight[0] * pSrc1[1] + rightWeight[0] * pSrc4[1] + add) >> shift;

    pDst += strideDst;
    pSrc2 += strideSrc;
    pSrc3 += strideSrc;
    pSrc4 += strideSrc;

    pDst[0] = ((sumWeight - leftWeight[0] - belowWeight[0]) * pDst[0] + leftWeight[0] * pSrc2[0] + belowWeight[0] * pSrc3[0] + add) >> shift;
    pDst[1] = ((sumWeight - belowWeight[0] - rightWeight[0]) * pDst[1] + belowWeight[0] * pSrc3[1] + rightWeight[0] * pSrc4[1] + add) >> shift;
  }
}
#endif

void InterPrediction::xSubBlockMotionCompensation(PredictionUnit &pu, PelUnitBuf &pcYuvPred)
{
  if (xCheckIdenticalMotion(pu))
  {
    xPredInterUni(pu, REF_PIC_LIST_0, pcYuvPred, false, false, true, true);
  }
  else
  {
    xPredInterBi(pu, pcYuvPred, true);
  }
}
#endif

int InterPrediction::rightShiftMSB(int numer, int denom)
{
  return numer >> floorLog2(denom);
}

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
void InterPrediction::initTplWeightTable()
{
  if (m_tplWeightTblInitialized)
  {
    return;
  }
  m_tplWeightTblInitialized = true;

  for (int hIdx = 0; hIdx < GEO_NUM_CU_SIZE; hIdx++)
  {
    int height = 1 << ( hIdx + GEO_MIN_CU_LOG2);

    for (int wIdx = 0; wIdx < GEO_NUM_CU_SIZE; wIdx++)
    {
      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t(&offset)[2] = g_weightOffset[splitDir][hIdx][wIdx];

        int16_t angle = g_GeoParams[splitDir][0];
        Pel* weight = &g_globalGeoWeightsTpl[g_angle2mask[angle]][GEO_TM_ADDED_WEIGHT_MASK_SIZE * GEO_WEIGHT_MASK_SIZE_EXT + GEO_TM_ADDED_WEIGHT_MASK_SIZE];
        if (g_angle2mirror[angle] == 2)
        {
          weight += ((GEO_WEIGHT_MASK_SIZE - 1 - offset[1]) * GEO_WEIGHT_MASK_SIZE_EXT + offset[0]);
        }
        else if (g_angle2mirror[angle] == 1)
        {
          weight += (offset[1] * GEO_WEIGHT_MASK_SIZE_EXT + (GEO_WEIGHT_MASK_SIZE - 1 - offset[0]));
        }
        else
        {
          weight += (offset[1] * GEO_WEIGHT_MASK_SIZE_EXT + offset[0]);
        }

        m_tplWeightTblDict[hIdx][wIdx][splitDir] = weight;
        m_tplWeightTbl = m_tplWeightTblDict[hIdx][wIdx];

        // Transpose weights for left template
        weight = getTplWeightTableCU<false, 1>(splitDir);
        int verticalOffset = g_angle2mirror[angle] == 2 ? -GEO_WEIGHT_MASK_SIZE_EXT : GEO_WEIGHT_MASK_SIZE_EXT;
        for (int h = 0; h < height; ++h)
        {
          m_tplColWeightTblDict[hIdx][wIdx][splitDir][h] = weight[0];
          weight += verticalOffset;
        }
      }
    }
  }

  m_tplWeightTbl    = nullptr;
  m_tplColWeightTbl = nullptr;
}
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
void InterPrediction::deriveGpmSplitMode(PredictionUnit& pu, MergeCtx &geoMrgCtx
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
                                       , MergeCtx(&geoTmMrgCtx)[GEO_NUM_TM_MV_CAND]
#endif
#if JVET_Y0065_GPM_INTRA
                                       , IntraPrediction* pcIntraPred
#endif
)
{
  if ( pu.cu->cs->pcv->isEncoder || !pu.cs->slice->getSPS()->getUseAltGPMSplitModeCode())
  {
    return;
  }

  uint8_t numValidInList = 0;
  uint8_t modeList[GEO_NUM_SIG_PARTMODE];
  bool refinedSplitMode = !PU::checkRprRefExistingInGpm(pu, geoMrgCtx, pu.geoMergeIdx0, geoMrgCtx, pu.geoMergeIdx1)
                       && xAMLGetCurBlkTemplate(pu, pu.lwidth(), pu.lheight());

  if (refinedSplitMode)
  {
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
    if (pu.tmMergeFlag)
    {
      Pel* pRefTopPart0 [GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart0[0], m_acYuvRefAMLTemplatePart0[2], nullptr                      }; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_A]
      Pel* pRefLeftPart0[GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart0[1], m_acYuvRefAMLTemplatePart0[3], nullptr                      }; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_A]
      Pel* pRefTopPart1 [GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart1[0], nullptr,                       m_acYuvRefAMLTemplatePart1[2]}; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_L]
      Pel* pRefLeftPart1[GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart1[1], nullptr,                       m_acYuvRefAMLTemplatePart1[3]}; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_L]
      fillPartGPMRefTemplate<0, false>(pu, geoTmMrgCtx[GEO_TM_SHAPE_AL], pu.geoMergeIdx0, -1, pRefTopPart0[GEO_TM_SHAPE_AL], pRefLeftPart0[GEO_TM_SHAPE_AL]);
      fillPartGPMRefTemplate<0, false>(pu, geoTmMrgCtx[GEO_TM_SHAPE_A ], pu.geoMergeIdx0, -1, pRefTopPart0[GEO_TM_SHAPE_A ], pRefLeftPart0[GEO_TM_SHAPE_A ]);
      fillPartGPMRefTemplate<1, false>(pu, geoTmMrgCtx[GEO_TM_SHAPE_AL], pu.geoMergeIdx1, -1, pRefTopPart1[GEO_TM_SHAPE_AL], pRefLeftPart1[GEO_TM_SHAPE_AL]);
      fillPartGPMRefTemplate<1, false>(pu, geoTmMrgCtx[GEO_TM_SHAPE_L ], pu.geoMergeIdx1, -1, pRefTopPart1[GEO_TM_SHAPE_L ], pRefLeftPart1[GEO_TM_SHAPE_L ]);

      // Update split mode
      getBestGeoTMModeList(pu, numValidInList, modeList, pRefTopPart0, pRefLeftPart0, pRefTopPart1, pRefLeftPart1);
      CHECK(pu.geoSyntaxMode < 0 || pu.geoSyntaxMode >= GEO_NUM_SIG_PARTMODE || pu.geoSyntaxMode >= numValidInList, "Invalid GEO split direction!");
      CHECK(numValidInList <= 0 || numValidInList > GEO_NUM_SIG_PARTMODE, "Error occurs");
      pu.geoSplitDir = modeList[pu.geoSyntaxMode];
      return;
    }
    else
#endif
    {
#if JVET_W0097_GPM_MMVD_TM
      int geoMmvdIdx0 = (pu.geoMMVDFlag0 ? pu.geoMMVDIdx0 : -1);
      int geoMmvdIdx1 = (pu.geoMMVDFlag1 ? pu.geoMMVDIdx1 : -1);
#else
      int geoMmvdIdx0 = -1;
      int geoMmvdIdx1 = -1;
#endif
      fillPartGPMRefTemplate<0>(pu, geoMrgCtx, pu.geoMergeIdx0, geoMmvdIdx0);
      fillPartGPMRefTemplate<1>(pu, geoMrgCtx, pu.geoMergeIdx1, geoMmvdIdx1);

#if JVET_Y0065_GPM_INTRA
      if (pu.gpmIntraFlag)
      {
        std::vector<Pel>* LUT = m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag() ? &m_pcReshape->getInvLUT() : nullptr;
        pcIntraPred->clearPrefilledIntraGPMRefTemplate();
        pcIntraPred->fillIntraGPMRefTemplateAll(pu, m_bAMLTemplateAvailabe[0], m_bAMLTemplateAvailabe[1], false, true, true, LUT, pu.geoMergeIdx0, pu.geoMergeIdx1);
      }
#endif
    }
  }
  else
  {
    m_bAMLTemplateAvailabe[0] = false;
    m_bAMLTemplateAvailabe[1] = false;
  }

  // Update split mode
#if JVET_Y0065_GPM_INTRA
  bool isIintra[2] = {false, false};
  Pel* pIntraRefTop [2][GEO_NUM_PARTITION_MODE];
  Pel* pIntraRefLeft[2][GEO_NUM_PARTITION_MODE];
  if (refinedSplitMode && pu.gpmIntraFlag)
  {
    isIintra[0] = pu.geoMergeIdx0 >= GEO_MAX_NUM_UNI_CANDS;
    isIintra[1] = pu.geoMergeIdx1 >= GEO_MAX_NUM_UNI_CANDS;
    for (uint8_t partIdx = 0; partIdx < 2; ++partIdx)
    {
      if (isIintra[partIdx])
      {
        uint8_t realCandIdx = (partIdx == 0 ? pu.geoMergeIdx0 : pu.geoMergeIdx1) - GEO_MAX_NUM_UNI_CANDS;
        for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
        {
          pIntraRefTop [partIdx][splitDir] = pcIntraPred->getPrefilledIntraGPMRefTemplate(partIdx, splitDir, realCandIdx, 0);
          pIntraRefLeft[partIdx][splitDir] = pcIntraPred->getPrefilledIntraGPMRefTemplate(partIdx, splitDir, realCandIdx, 1);
        }
      }
    }
  }
#endif

  getBestGeoModeList(pu, numValidInList, modeList, m_acYuvRefAMLTemplatePart0[0], m_acYuvRefAMLTemplatePart0[1], m_acYuvRefAMLTemplatePart1[0], m_acYuvRefAMLTemplatePart1[1]
#if JVET_Y0065_GPM_INTRA
                   , isIintra[0] ? pIntraRefTop [0] : nullptr
                   , isIintra[0] ? pIntraRefLeft[0] : nullptr
                   , isIintra[1] ? pIntraRefTop [1] : nullptr
                   , isIintra[1] ? pIntraRefLeft[1] : nullptr
#endif
  );
  CHECK(pu.geoSyntaxMode < 0 || pu.geoSyntaxMode >= GEO_NUM_SIG_PARTMODE || pu.geoSyntaxMode >= numValidInList, "Invalid GEO split direction!");
  CHECK(numValidInList <= 0 || numValidInList > GEO_NUM_SIG_PARTMODE, "Error occurs");
  pu.geoSplitDir = modeList[pu.geoSyntaxMode];
}
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
void InterPrediction::motionCompensationGeo( CodingUnit &cu, MergeCtx &geoMrgCtx
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
                                           , MergeCtx(&geoTmMrgCtx)[GEO_NUM_TM_MV_CAND]
#endif
#if JVET_Y0065_GPM_INTRA
                                           , IntraPrediction* pcIntraPred, std::vector<Pel>* reshapeLUT
#endif
)
#else
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
#if JVET_Y0065_GPM_INTRA
void InterPrediction::motionCompensationGeo( CodingUnit &cu, MergeCtx &geoMrgCtx, MergeCtx &geoTmMrgCtx0, MergeCtx &geoTmMrgCtx1, IntraPrediction* pcIntraPred, std::vector<Pel>* reshapeLUT )
#else
void InterPrediction::motionCompensationGeo(CodingUnit &cu, MergeCtx &geoMrgCtx, MergeCtx &geoTmMrgCtx0, MergeCtx &geoTmMrgCtx1)
#endif
#else
#if JVET_Y0065_GPM_INTRA
void InterPrediction::motionCompensationGeo( CodingUnit &cu, MergeCtx &geoMrgCtx, IntraPrediction* pcIntraPred, std::vector<Pel>* reshapeLUT )
#else
void InterPrediction::motionCompensationGeo( CodingUnit &cu, MergeCtx &geoMrgCtx )
#endif
#endif
#endif
{
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  deriveGpmSplitMode(*cu.firstPU, geoMrgCtx
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
                   , geoTmMrgCtx
#endif
#if JVET_Y0065_GPM_INTRA
                   , pcIntraPred
#endif
  );
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
  MergeCtx& geoTmMrgCtx0 = geoTmMrgCtx[g_geoTmShape[0][g_GeoParams[cu.firstPU->geoSplitDir][0]]];
  MergeCtx& geoTmMrgCtx1 = geoTmMrgCtx[g_geoTmShape[1][g_GeoParams[cu.firstPU->geoSplitDir][0]]];
#endif
#endif

  const uint8_t splitDir = cu.firstPU->geoSplitDir;
  const uint8_t candIdx0 = cu.firstPU->geoMergeIdx0;
  const uint8_t candIdx1 = cu.firstPU->geoMergeIdx1;
#if JVET_W0097_GPM_MMVD_TM
  const bool    geoMMVDFlag0 = cu.firstPU->geoMMVDFlag0;
  const uint8_t geoMMVDIdx0 = cu.firstPU->geoMMVDIdx0;
  const bool    geoMMVDFlag1 = cu.firstPU->geoMMVDFlag1;
  const uint8_t geoMMVDIdx1 = cu.firstPU->geoMMVDIdx1;
#if TM_MRG
  const bool    geoTmFlag0 = cu.firstPU->geoTmFlag0;
  const bool    geoTmFlag1 = cu.firstPU->geoTmFlag1;
#endif
#endif

  for( auto &pu : CU::traversePUs( cu ) )
  {
    const UnitArea localUnitArea( cu.cs->area.chromaFormat, Area( 0, 0, pu.lwidth(), pu.lheight() ) );
    PelUnitBuf tmpGeoBuf0 = m_geoPartBuf[0].getBuf( localUnitArea );
    PelUnitBuf tmpGeoBuf1 = m_geoPartBuf[1].getBuf( localUnitArea );
    PelUnitBuf predBuf    = cu.cs->getPredBuf( pu );
#if JVET_Y0065_GPM_INTRA
    bool isIntra0 = candIdx0 >= GEO_MAX_NUM_UNI_CANDS;
    bool isIntra1 = candIdx1 >= GEO_MAX_NUM_UNI_CANDS;
    if (isIntra0)
    {
      PU::getGeoIntraMPMs(pu, pu.intraMPM, splitDir, g_geoTmShape[0][g_GeoParams[pu.geoSplitDir][0]]);
      pu.intraDir[0] = pu.intraMPM[candIdx0 - GEO_MAX_NUM_UNI_CANDS];
      pcIntraPred->initIntraPatternChType(cu, pu.Y());
      pcIntraPred->predIntraAng(COMPONENT_Y, tmpGeoBuf0.Y(), pu);
      if (isChromaEnabled(pu.chromaFormat))
      {
        pu.intraDir[1] = pu.intraDir[0];
        pcIntraPred->initIntraPatternChType(cu, pu.Cb());
        pcIntraPred->predIntraAng(COMPONENT_Cb, tmpGeoBuf0.Cb(), pu);
        pcIntraPred->initIntraPatternChType(cu, pu.Cr());
        pcIntraPred->predIntraAng(COMPONENT_Cr, tmpGeoBuf0.Cr(), pu);
      }
    }
    else
    {
#endif
#if JVET_W0097_GPM_MMVD_TM
#if TM_MRG
    if (geoTmFlag0)
    {
      geoTmMrgCtx0.setMergeInfo(pu, candIdx0);
    }
    else
#endif
    if (geoMMVDFlag0)
    {
      geoMrgCtx.setGeoMmvdMergeInfo(pu, candIdx0, geoMMVDIdx0);
    }
    else
#endif
    geoMrgCtx.setMergeInfo( pu, candIdx0 );

    motionCompensation(pu, tmpGeoBuf0, REF_PIC_LIST_X, true, isChromaEnabled(pu.chromaFormat)); // TODO: check 4:0:0 interaction with weighted prediction.
    if( g_mctsDecCheckEnabled && !MCTSHelper::checkMvBufferForMCTSConstraint( pu, true ) )
    {
      printf( "DECODER_GEO_PU: pu motion vector across tile boundaries (%d,%d,%d,%d)\n", pu.lx(), pu.ly(), pu.lwidth(), pu.lheight() );
    }
#if JVET_Y0065_GPM_INTRA
      if (isIntra1)
      {
        tmpGeoBuf0.roundToOutputBitdepth(tmpGeoBuf0, cu.slice->clpRngs());
#if ENABLE_OBMC
#if JVET_W0123_TIMD_FUSION
        PU::spanMotionInfo2(pu);
#else
        PU::spanMotionInfo(pu);
#endif
        cu.isobmcMC = true;
        subBlockOBMC(pu, &tmpGeoBuf0);
        cu.isobmcMC = false;
#endif
      }
    }

    if (isIntra1)
    {
      PU::getGeoIntraMPMs(pu, pu.intraMPM+GEO_MAX_NUM_INTRA_CANDS, splitDir, g_geoTmShape[1][g_GeoParams[pu.geoSplitDir][0]]);
      pu.intraDir[0] = pu.intraMPM[candIdx1 - GEO_MAX_NUM_UNI_CANDS + GEO_MAX_NUM_INTRA_CANDS];
      pcIntraPred->initIntraPatternChType(cu, pu.Y());
      pcIntraPred->predIntraAng(COMPONENT_Y, tmpGeoBuf1.Y(), pu);
      if (isChromaEnabled(pu.chromaFormat))
      {
        pu.intraDir[1] = pu.intraDir[0];
        pcIntraPred->initIntraPatternChType(cu, pu.Cb());
        pcIntraPred->predIntraAng(COMPONENT_Cb, tmpGeoBuf1.Cb(), pu);
        pcIntraPred->initIntraPatternChType(cu, pu.Cr());
        pcIntraPred->predIntraAng(COMPONENT_Cr, tmpGeoBuf1.Cr(), pu);
      }
    }
    else
    {
#endif
#if JVET_W0097_GPM_MMVD_TM
#if TM_MRG
    if (geoTmFlag1)
    {
      geoTmMrgCtx1.setMergeInfo(pu, candIdx1);
    }
    else
#endif
    if (geoMMVDFlag1)
    {
      geoMrgCtx.setGeoMmvdMergeInfo(pu, candIdx1, geoMMVDIdx1);
    }
    else
#endif
    geoMrgCtx.setMergeInfo( pu, candIdx1 );

    motionCompensation(pu, tmpGeoBuf1, REF_PIC_LIST_X, true, isChromaEnabled(pu.chromaFormat)); // TODO: check 4:0:0 interaction with weighted prediction.
    if( g_mctsDecCheckEnabled && !MCTSHelper::checkMvBufferForMCTSConstraint( pu, true ) )
    {
      printf( "DECODER_GEO_PU: pu motion vector across tile boundaries (%d,%d,%d,%d)\n", pu.lx(), pu.ly(), pu.lwidth(), pu.lheight() );
    }
#if JVET_Y0065_GPM_INTRA
      if (isIntra0)
      {
        tmpGeoBuf1.roundToOutputBitdepth(tmpGeoBuf1, cu.slice->clpRngs());
#if ENABLE_OBMC
#if JVET_W0123_TIMD_FUSION
        PU::spanMotionInfo2(pu);
#else
        PU::spanMotionInfo(pu);
#endif
        cu.isobmcMC = true;
        subBlockOBMC(pu, &tmpGeoBuf1);
        cu.isobmcMC = false;
#endif
      }
    }
    if (pu.gpmIntraFlag)
    {
      if (reshapeLUT)
      {
        if (!isIntra1)
        {
          tmpGeoBuf1.Y().rspSignal(*reshapeLUT);
        }
        else if (!isIntra0)
        {
          tmpGeoBuf0.Y().rspSignal(*reshapeLUT);
        }
      }
      weightedGeoBlkRounded(pu, splitDir, isChromaEnabled(pu.chromaFormat)? MAX_NUM_CHANNEL_TYPE : CHANNEL_TYPE_LUMA, predBuf, tmpGeoBuf0, tmpGeoBuf1);
    }
    else
#endif
    weightedGeoBlk(pu, splitDir, isChromaEnabled(pu.chromaFormat)? MAX_NUM_CHANNEL_TYPE : CHANNEL_TYPE_LUMA, predBuf, tmpGeoBuf0, tmpGeoBuf1);
  }
}

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
void InterPrediction::getBestGeoTMModeList(PredictionUnit &pu, uint8_t& numValidInList, uint8_t(&modeList)[GEO_NUM_SIG_PARTMODE], Pel* (&pRefTopPart0)[GEO_NUM_TM_MV_CAND], Pel* (&pRefLeftPart0)[GEO_NUM_TM_MV_CAND], Pel* (&pRefTopPart1)[GEO_NUM_TM_MV_CAND], Pel* (&pRefLeftPart1)[GEO_NUM_TM_MV_CAND])
{
  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    for (int i = 0; i < GEO_NUM_SIG_PARTMODE; ++i)
    {
      modeList[i] = i;
    }
    numValidInList = GEO_NUM_SIG_PARTMODE;
    return;
  }

  // Check split mode cost
  uint32_t uiCost[GEO_NUM_PARTITION_MODE] = { 0, };

  if (m_bAMLTemplateAvailabe[0])
  {
    SizeType   szPerLine       = pu.lwidth();
    PelUnitBuf pcBufPredCurTop = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTop = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTopPart0[GEO_NUM_TM_MV_CAND] = {PelUnitBuf(), 
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart0[GEO_TM_SHAPE_AL], szPerLine, GEO_MODE_SEL_TM_SIZE)),
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart0[GEO_TM_SHAPE_A ], szPerLine, GEO_MODE_SEL_TM_SIZE)),
                                                           PelUnitBuf()};
    PelUnitBuf pcBufPredRefTopPart1[GEO_NUM_TM_MV_CAND] = {PelUnitBuf(),
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart1[GEO_TM_SHAPE_AL], szPerLine, GEO_MODE_SEL_TM_SIZE)),
                                                           PelUnitBuf(),
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart1[GEO_TM_SHAPE_L ], szPerLine, GEO_MODE_SEL_TM_SIZE))};

    DistParam cDistParam;
    cDistParam.applyWeight = false;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
    {
      uint8_t shapeIdx0 = g_geoTmShape[0][g_GeoParams[splitDir][0]];
      uint8_t shapeIdx1 = g_geoTmShape[1][g_GeoParams[splitDir][0]];
      weightedGeoTpl<true>(pu, splitDir, pcBufPredRefTop, pcBufPredRefTopPart0[shapeIdx0], pcBufPredRefTopPart1[shapeIdx1]);
      uint32_t tempDist = (uint32_t)cDistParam.distFunc(cDistParam);
      uiCost[splitDir] += tempDist;
    }
  }

  if (m_bAMLTemplateAvailabe[1])
  {
    SizeType   szPerLine          = pu.lheight();
    PelUnitBuf pcBufPredCurLeftTr = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufPredRefLeftTr = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufPredRefLeft   = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], GEO_MODE_SEL_TM_SIZE, szPerLine));
    PelUnitBuf pcBufPredRefLeftPart0[GEO_NUM_TM_MV_CAND] = {PelUnitBuf(),
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart0[GEO_TM_SHAPE_AL], GEO_MODE_SEL_TM_SIZE, szPerLine)),
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart0[GEO_TM_SHAPE_A ], GEO_MODE_SEL_TM_SIZE, szPerLine)),
                                                            PelUnitBuf()};
    PelUnitBuf pcBufPredRefLeftPart1[GEO_NUM_TM_MV_CAND] = {PelUnitBuf(),
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart1[GEO_TM_SHAPE_AL], GEO_MODE_SEL_TM_SIZE, szPerLine)),
                                                            PelUnitBuf(),
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart1[GEO_TM_SHAPE_L ], GEO_MODE_SEL_TM_SIZE, szPerLine))};

    DistParam cDistParam;
    cDistParam.applyWeight = false;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeftTr.Y(), pcBufPredRefLeftTr.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false); // To enable SIMD for cost computation
    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
    {
      uint8_t shapeIdx0  = g_geoTmShape[0][g_GeoParams[splitDir][0]];
      uint8_t shapeIdx1  = g_geoTmShape[1][g_GeoParams[splitDir][0]];
      weightedGeoTpl<false>(pu, splitDir, pcBufPredRefLeft, pcBufPredRefLeftPart0[shapeIdx0], pcBufPredRefLeftPart1[shapeIdx1]);
      uint32_t tempDist = (uint32_t)cDistParam.distFunc(cDistParam);
      uiCost[splitDir] += tempDist;
    }
  }

  // Find best N candidates
  numValidInList = (uint8_t)getIndexMappingTableToSortedArray1D<uint32_t, GEO_NUM_PARTITION_MODE, uint8_t, GEO_NUM_SIG_PARTMODE>(uiCost, modeList);

}
#endif

void InterPrediction::getBestGeoModeList(PredictionUnit &pu, uint8_t& numValidInList, uint8_t(&modeList)[GEO_NUM_SIG_PARTMODE], Pel* pRefTopPart0, Pel* pRefLeftPart0, Pel* pRefTopPart1, Pel* pRefLeftPart1
#if JVET_Y0065_GPM_INTRA
                                       , Pel** pIntraRefTopPart0, Pel** pIntraRefLeftPart0, Pel** pIntraRefTopPart1, Pel** pIntraRefLeftPart1
#endif
)
{
  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    for (int i = 0; i < GEO_NUM_SIG_PARTMODE; ++i)
    {
      modeList[i] = i;
    }
    numValidInList = GEO_NUM_SIG_PARTMODE;
    return;
  }

  // Check split mode cost
  uint32_t uiCost[GEO_NUM_PARTITION_MODE] = { 0, };

  if (m_bAMLTemplateAvailabe[0])
  {
    SizeType   szPerLine            = pu.lwidth();
    PelUnitBuf pcBufPredCurTop      = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTop      = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTopPart0 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart0,                szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTopPart1 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart1,                szPerLine, GEO_MODE_SEL_TM_SIZE));

    DistParam cDistParam;
    cDistParam.applyWeight = false;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
    {
#if JVET_Y0065_GPM_INTRA
      pcBufPredRefTopPart0.Y().buf = pIntraRefTopPart0 == nullptr ? pcBufPredRefTopPart0.Y().buf : pIntraRefTopPart0[splitDir];
      pcBufPredRefTopPart1.Y().buf = pIntraRefTopPart1 == nullptr ? pcBufPredRefTopPart1.Y().buf : pIntraRefTopPart1[splitDir];
#endif
      weightedGeoTpl<true>(pu, splitDir, pcBufPredRefTop, pcBufPredRefTopPart0, pcBufPredRefTopPart1);
      uint32_t tempDist = (uint32_t)cDistParam.distFunc(cDistParam);
      uiCost[splitDir] += tempDist;
    }
  }

  if (m_bAMLTemplateAvailabe[1])
  {
    SizeType   szPerLine             = pu.lheight();
    PelUnitBuf pcBufPredCurLeftTr    = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufPredRefLeftTr    = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufPredRefLeft      = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], GEO_MODE_SEL_TM_SIZE, szPerLine));
    PelUnitBuf pcBufPredRefLeftPart0 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart0,               GEO_MODE_SEL_TM_SIZE, szPerLine));
    PelUnitBuf pcBufPredRefLeftPart1 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart1,               GEO_MODE_SEL_TM_SIZE, szPerLine));

    DistParam cDistParam;
    cDistParam.applyWeight = false;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeftTr.Y(), pcBufPredRefLeftTr.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false); // To enable SIMD for cost computation
    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
    {
#if JVET_Y0065_GPM_INTRA
      pcBufPredRefLeftPart0.Y().buf = pIntraRefLeftPart0 == nullptr ? pcBufPredRefLeftPart0.Y().buf : pIntraRefLeftPart0[splitDir];
      pcBufPredRefLeftPart1.Y().buf = pIntraRefLeftPart1 == nullptr ? pcBufPredRefLeftPart1.Y().buf : pIntraRefLeftPart1[splitDir];
#endif
      weightedGeoTpl<false>(pu, splitDir, pcBufPredRefLeft, pcBufPredRefLeftPart0, pcBufPredRefLeftPart1);
      uint32_t tempDist = (uint32_t)cDistParam.distFunc(cDistParam);
      uiCost[splitDir] += tempDist;
    }
  }

  // Find best N candidates
  numValidInList = (uint8_t)getIndexMappingTableToSortedArray1D<uint32_t, GEO_NUM_PARTITION_MODE, uint8_t, GEO_NUM_SIG_PARTMODE>(uiCost, modeList);

}

template <bool trueTFalseL>
void InterPrediction::weightedGeoTpl( PredictionUnit &pu, const uint8_t splitDir, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1)
{
  m_if.weightedGeoTpl<trueTFalseL>( pu, splitDir, predDst, predSrc0, predSrc1 );
}
#endif

void InterPrediction::weightedGeoBlk( PredictionUnit &pu, const uint8_t splitDir, int32_t channel, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1)
{
  if( channel == CHANNEL_TYPE_LUMA )
  {
    m_if.weightedGeoBlk( pu, pu.lumaSize().width, pu.lumaSize().height, COMPONENT_Y, splitDir, predDst, predSrc0, predSrc1 );
  }
  else if( channel == CHANNEL_TYPE_CHROMA )
  {
    m_if.weightedGeoBlk( pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cb, splitDir, predDst, predSrc0, predSrc1 );
    m_if.weightedGeoBlk( pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cr, splitDir, predDst, predSrc0, predSrc1 );
  }
  else
  {
    m_if.weightedGeoBlk( pu, pu.lumaSize().width,   pu.lumaSize().height,   COMPONENT_Y,  splitDir, predDst, predSrc0, predSrc1 );
    if (isChromaEnabled(pu.chromaFormat))
    {
      m_if.weightedGeoBlk(pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cb, splitDir, predDst, predSrc0,
                          predSrc1);
      m_if.weightedGeoBlk(pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cr, splitDir, predDst, predSrc0,
                          predSrc1);
    }
  }
}

#if JVET_Y0065_GPM_INTRA
void InterPrediction::weightedGeoBlkRounded( PredictionUnit &pu, const uint8_t splitDir, int32_t channel, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1)
{
  if( channel == CHANNEL_TYPE_LUMA )
  {
    m_if.weightedGeoBlkRounded( pu, pu.lumaSize().width, pu.lumaSize().height, COMPONENT_Y, splitDir, predDst, predSrc0, predSrc1 );
  }
  else if( channel == CHANNEL_TYPE_CHROMA )
  {
    m_if.weightedGeoBlkRounded( pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cb, splitDir, predDst, predSrc0, predSrc1 );
    m_if.weightedGeoBlkRounded( pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cr, splitDir, predDst, predSrc0, predSrc1 );
  }
  else
  {
    m_if.weightedGeoBlkRounded( pu, pu.lumaSize().width,   pu.lumaSize().height,   COMPONENT_Y,  splitDir, predDst, predSrc0, predSrc1 );
    if (isChromaEnabled(pu.chromaFormat))
    {
      m_if.weightedGeoBlkRounded(pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cb, splitDir, predDst, predSrc0,
                          predSrc1);
      m_if.weightedGeoBlkRounded(pu, pu.chromaSize().width, pu.chromaSize().height, COMPONENT_Cr, splitDir, predDst, predSrc0,
                          predSrc1);
    }
  }
}
#endif

void InterPrediction::xPrefetch(PredictionUnit& pu, PelUnitBuf &pcPad, RefPicList refId, bool forLuma)
{
  int offset, width, height;
  Mv cMv;
  const Picture* refPic = pu.cu->slice->getRefPic( refId, pu.refIdx[refId] )->unscaledPic;
  int mvShift = (MV_FRACTIONAL_BITS_INTERNAL);

  int start = 0;
  int end = MAX_NUM_COMPONENT;

  start = forLuma ? 0 : 1;
  end = forLuma ? 1 : MAX_NUM_COMPONENT;

  for (int compID = start; compID < end; compID++)
  {
    cMv = Mv(pu.mv[refId].getHor(), pu.mv[refId].getVer());
#if IF_12TAP
    pcPad.bufs[compID].stride = (pcPad.bufs[compID].width + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA(0));
    int filtersize = (compID == (COMPONENT_Y)) ? NTAPS_LUMA(0) : NTAPS_CHROMA;
#else
    pcPad.bufs[compID].stride = (pcPad.bufs[compID].width + (2 * DMVR_NUM_ITERATION) + NTAPS_LUMA);
    int filtersize = (compID == (COMPONENT_Y)) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif
    width = pcPad.bufs[compID].width;
    height = pcPad.bufs[compID].height;
    offset = (DMVR_NUM_ITERATION) * (pcPad.bufs[compID].stride + 1);

    int mvshiftTempHor = mvShift + getComponentScaleX((ComponentID)compID, pu.chromaFormat);
    int mvshiftTempVer = mvShift + getComponentScaleY((ComponentID)compID, pu.chromaFormat);
    width += (filtersize - 1);
    height += (filtersize - 1);
    cMv += Mv(-(((filtersize >> 1) - 1) << mvshiftTempHor),
      -(((filtersize >> 1) - 1) << mvshiftTempVer));
    bool wrapRef = false;
    if( refPic->isWrapAroundEnabled( pu.cs->pps ) )
    {
      wrapRef = wrapClipMv( cMv, pu.blocks[0].pos(), pu.blocks[0].size(), pu.cs->sps, pu.cs->pps );
    }
    else
    {
      clipMv( cMv, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps );
    }
    /* Pre-fetch similar to HEVC*/
    {
      CPelBuf refBuf;
      Position Rec_offset = pu.blocks[compID].pos().offset(cMv.getHor() >> mvshiftTempHor, cMv.getVer() >> mvshiftTempVer);
      refBuf = refPic->getRecoBuf(CompArea((ComponentID)compID, pu.chromaFormat, Rec_offset, pu.blocks[compID].size()), wrapRef);
      PelBuf &dstBuf = pcPad.bufs[compID];
      g_pelBufOP.copyBuffer((Pel *)refBuf.buf, refBuf.stride, ((Pel *)dstBuf.buf) + offset, dstBuf.stride, width, height);
    }
  }
}

void InterPrediction::xPad(PredictionUnit& pu, PelUnitBuf &pcPad, RefPicList refId)
{
  int offset = 0, width, height;
  int padsize;
  Mv cMv;
  for (int compID = 0; compID < getNumberValidComponents(pu.chromaFormat); compID++)
  {
#if IF_12TAP
    int filtersize = (compID == (COMPONENT_Y)) ? NTAPS_LUMA(0) : NTAPS_CHROMA;
#else
    int filtersize = (compID == (COMPONENT_Y)) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif
    width = pcPad.bufs[compID].width;
    height = pcPad.bufs[compID].height;
    offset = (DMVR_NUM_ITERATION) * (pcPad.bufs[compID].stride + 1);
    /*using the larger padsize for 422*/
    padsize = (DMVR_NUM_ITERATION) >> getComponentScaleY((ComponentID)compID, pu.chromaFormat);
    width += (filtersize - 1);
    height += (filtersize - 1);
    /*padding on all side of size DMVR_PAD_LENGTH*/
    g_pelBufOP.padding(pcPad.bufs[compID].buf + offset, pcPad.bufs[compID].stride, width, height, padsize);
  }
}

inline int32_t div_for_maxq7(int64_t N, int64_t D)
{
  int32_t sign, q;
  sign = 0;
  if (N < 0)
  {
    sign = 1;
    N = -N;
  }

  q = 0;
  D = (D << 3);
  if (N >= D)
  {
    N -= D;
    q++;
  }
  q = (q << 1);

  D = (D >> 1);
  if (N >= D)
  {
    N -= D;
    q++;
  }
  q = (q << 1);

  if (N >= (D >> 1))
  {
    q++;
  }
  if (sign)
  {
    return (-q);
  }
  return(q);
}

void xSubPelErrorSrfc(uint64_t *sadBuffer, int32_t *deltaMv)
{
  int64_t numerator, denominator;
  int32_t mvDeltaSubPel;
  int32_t mvSubPelLvl = 4;/*1: half pel, 2: Qpel, 3:1/8, 4: 1/16*/
                                                        /*horizontal*/
  numerator   = (int64_t)((sadBuffer[1] - sadBuffer[3]) << mvSubPelLvl);
  denominator = (int64_t)((sadBuffer[1] + sadBuffer[3] - (sadBuffer[0] << 1)));

#if MULTI_PASS_DMVR
    if (denominator > 0)
    {
      if ((sadBuffer[1] != sadBuffer[0]) && (sadBuffer[3] != sadBuffer[0]))
      {
        mvDeltaSubPel = div_for_maxq7(numerator, denominator);
        deltaMv[0] = (mvDeltaSubPel);
      }
      else
      {
        if (sadBuffer[1] == sadBuffer[0])
        {
          deltaMv[0] = -8;// half pel
        }
        else
        {
          deltaMv[0] = 8;// half pel
        }
      }
    }
    else
    {
      if (sadBuffer[1] < sadBuffer[3])
      {
        deltaMv[0] = -8;
      }
      else if (sadBuffer[1] == sadBuffer[3])
      {
        deltaMv[0] = 0;
      }
      else
      {
        deltaMv[0] = 8;
      }
    }
#else
    if (0 != denominator)
    {
      if ((sadBuffer[1] != sadBuffer[0]) && (sadBuffer[3] != sadBuffer[0]))
      {
        mvDeltaSubPel = div_for_maxq7(numerator, denominator);
        deltaMv[0] = (mvDeltaSubPel);
      }
      else
      {
        if (sadBuffer[1] == sadBuffer[0])
        {
          deltaMv[0] = -8;// half pel
        }
        else
        {
          deltaMv[0] = 8;// half pel
        }
      }
    }
#endif

    /*vertical*/
    numerator = (int64_t)((sadBuffer[2] - sadBuffer[4]) << mvSubPelLvl);
    denominator = (int64_t)((sadBuffer[2] + sadBuffer[4] - (sadBuffer[0] << 1)));
#if MULTI_PASS_DMVR
    if (denominator > 0)
    {
      if ((sadBuffer[2] != sadBuffer[0]) && (sadBuffer[4] != sadBuffer[0]))
      {
        mvDeltaSubPel = div_for_maxq7(numerator, denominator);
        deltaMv[1] = (mvDeltaSubPel);
      }
      else
      {
        if (sadBuffer[2] == sadBuffer[0])
        {
          deltaMv[1] = -8;// half pel
        }
        else
        {
          deltaMv[1] = 8;// half pel
        }
      }
    }
    else
    {
      if (sadBuffer[2] < sadBuffer[4])
      {
        deltaMv[1] = -8;
      }
      else if (sadBuffer[2] == sadBuffer[4])
      {
        deltaMv[1] = 0;
      }
      else
      {
        deltaMv[1] = 8;
      }
    }
#else
    if (0 != denominator)
    {
      if ((sadBuffer[2] != sadBuffer[0]) && (sadBuffer[4] != sadBuffer[0]))
      {
        mvDeltaSubPel = div_for_maxq7(numerator, denominator);
        deltaMv[1] = (mvDeltaSubPel);
      }
      else
      {
        if (sadBuffer[2] == sadBuffer[0])
        {
          deltaMv[1] = -8;// half pel
        }
        else
        {
          deltaMv[1] = 8;// half pel
        }
      }
    }
#endif
  return;
}

void InterPrediction::xBIPMVRefine(int bd, Pel *pRefL0, Pel *pRefL1, uint64_t& minCost, int16_t *deltaMV, uint64_t *pSADsArray, int width, int height)
{
  const int32_t refStrideL0 = m_biLinearBufStride;
  const int32_t refStrideL1 = m_biLinearBufStride;
  Pel *pRefL0Orig = pRefL0;
  Pel *pRefL1Orig = pRefL1;
  for (int nIdx = 0; (nIdx < 25); ++nIdx)
  {
    int32_t sadOffset = ((m_pSearchOffset[nIdx].getVer() * ((2 * DMVR_NUM_ITERATION) + 1)) + m_pSearchOffset[nIdx].getHor());
    pRefL0 = pRefL0Orig + m_pSearchOffset[nIdx].hor + (m_pSearchOffset[nIdx].ver * refStrideL0);
    pRefL1 = pRefL1Orig - m_pSearchOffset[nIdx].hor - (m_pSearchOffset[nIdx].ver * refStrideL1);
    if (*(pSADsArray + sadOffset) == MAX_UINT64)
    {
      const uint64_t cost = xDMVRCost(bd, pRefL0, refStrideL0, pRefL1, refStrideL1, width, height);
      *(pSADsArray + sadOffset) = cost;
    }
    if (*(pSADsArray + sadOffset) < minCost)
    {
      minCost = *(pSADsArray + sadOffset);
      deltaMV[0] = m_pSearchOffset[nIdx].getHor();
      deltaMV[1] = m_pSearchOffset[nIdx].getVer();
    }
  }
}

void InterPrediction::xFinalPaddedMCForDMVR(PredictionUnit &pu, PelUnitBuf &pcYuvSrc0, PelUnitBuf &pcYuvSrc1,
                                            PelUnitBuf &pcPad0, PelUnitBuf &pcPad1, const bool bioApplied,
                                            const Mv mergeMV[NUM_REF_PIC_LIST_01], bool blockMoved)
{
  int offset, deltaIntMvX, deltaIntMvY;

  PelUnitBuf pcYUVTemp = pcYuvSrc0;
  PelUnitBuf pcPadTemp = pcPad0;
  /*always high precision MVs are used*/
  int mvShift = MV_FRACTIONAL_BITS_INTERNAL;

  for (int k = 0; k < NUM_REF_PIC_LIST_01; k++)
  {
    RefPicList refId = (RefPicList)k;
    Mv cMv = pu.mv[refId];
    m_iRefListIdx = refId;
    const Picture* refPic = pu.cu->slice->getRefPic( refId, pu.refIdx[refId] )->unscaledPic;
    Mv cMvClipped = cMv;
    if( !pu.cs->pps->getWrapAroundEnabledFlag() )
    {
      clipMv( cMvClipped, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps );
    }

    Mv startMv = mergeMV[refId];

    if( g_mctsDecCheckEnabled && !MCTSHelper::checkMvForMCTSConstraint( pu, startMv, MV_PRECISION_INTERNAL ) )
    {
      const Area& tileArea = pu.cs->picture->mctsInfo.getTileArea();
      printf( "Attempt an access over tile boundary at block %d,%d %d,%d with MV %d,%d (in Tile TL: %d,%d BR: %d,%d)\n",
        pu.lx(), pu.ly(), pu.lwidth(), pu.lheight(), startMv.getHor(), startMv.getVer(), tileArea.topLeft().x, tileArea.topLeft().y, tileArea.bottomRight().x, tileArea.bottomRight().y );
      THROW( "MCTS constraint failed!" );
    }
    for (int compID = 0; compID < getNumberValidComponents(pu.chromaFormat); compID++)
    {
      Pel *srcBufPelPtr = NULL;
      int pcPadstride = 0;
      if (blockMoved || (compID == 0))
      {
        pcPadstride = pcPadTemp.bufs[compID].stride;
        int mvshiftTempHor = mvShift + getComponentScaleX((ComponentID)compID, pu.chromaFormat);
        int mvshiftTempVer = mvShift + getComponentScaleY((ComponentID)compID, pu.chromaFormat);
        int leftPixelExtra;
        if (compID == COMPONENT_Y)
        {
#if IF_12TAP
          leftPixelExtra = (NTAPS_LUMA(0) >> 1) - 1;
#else
          leftPixelExtra = (NTAPS_LUMA >> 1) - 1;
#endif
        }
        else
        {
          leftPixelExtra = (NTAPS_CHROMA >> 1) - 1;
        }
        PelBuf &srcBuf = pcPadTemp.bufs[compID];
        deltaIntMvX = (cMv.getHor() >> mvshiftTempHor) -
          (startMv.getHor() >> mvshiftTempHor);
        deltaIntMvY = (cMv.getVer() >> mvshiftTempVer) -
          (startMv.getVer() >> mvshiftTempVer);

        CHECK((abs(deltaIntMvX) > DMVR_NUM_ITERATION) || (abs(deltaIntMvY) > DMVR_NUM_ITERATION), "not expected DMVR movement");
        offset = (DMVR_NUM_ITERATION + leftPixelExtra) * (pcPadTemp.bufs[compID].stride + 1);
        offset += (deltaIntMvY)* pcPadTemp.bufs[compID].stride;
        offset += (deltaIntMvX);
        srcBufPelPtr = (srcBuf.buf + offset);
      }
      JVET_J0090_SET_CACHE_ENABLE(false);
      xPredInterBlk((ComponentID) compID, pu, refPic, cMvClipped, pcYUVTemp, true,
                    pu.cs->slice->getClpRngs().comp[compID], bioApplied, false,
                    pu.cu->slice->getScalingRatio(refId, pu.refIdx[refId]), 0, 0, 0, srcBufPelPtr, pcPadstride);
      JVET_J0090_SET_CACHE_ENABLE(false);
    }
    pcYUVTemp = pcYuvSrc1;
    pcPadTemp = pcPad1;
  }
}

uint64_t InterPrediction::xDMVRCost(int bitDepth, Pel* pOrg, uint32_t refStride, const Pel* pRef, uint32_t orgStride, int width, int height)
{
  DistParam cDistParam;
  cDistParam.applyWeight = false;
  cDistParam.useMR = false;
  m_pcRdCost->setDistParam(cDistParam, pOrg, pRef, orgStride, refStride, bitDepth, COMPONENT_Y, width, height, 1);
  uint64_t uiCost = cDistParam.distFunc(cDistParam);
  return uiCost>>1;
}

void xDMVRSubPixelErrorSurface(bool notZeroCost, int16_t *totalDeltaMV, int16_t *deltaMV, uint64_t *pSADsArray)
{
  int sadStride = (((2 * DMVR_NUM_ITERATION) + 1));
  uint64_t sadbuffer[5];
  if (notZeroCost && (abs(totalDeltaMV[0]) != (2 << MV_FRACTIONAL_BITS_INTERNAL))
    && (abs(totalDeltaMV[1]) != (2 << MV_FRACTIONAL_BITS_INTERNAL)))
  {
    int32_t tempDeltaMv[2] = { 0,0 };
    sadbuffer[0] = pSADsArray[0];
    sadbuffer[1] = pSADsArray[-1];
    sadbuffer[2] = pSADsArray[-sadStride];
    sadbuffer[3] = pSADsArray[1];
    sadbuffer[4] = pSADsArray[sadStride];
    xSubPelErrorSrfc(sadbuffer, tempDeltaMv);
    totalDeltaMV[0] += tempDeltaMv[0];
    totalDeltaMV[1] += tempDeltaMv[1];
  }
}

void InterPrediction::xinitMC(PredictionUnit& pu, const ClpRngs &clpRngs)
{
  const int refIdx0 = pu.refIdx[0];
  const int refIdx1 = pu.refIdx[1];
  /*use merge MV as starting MV*/
  Mv mergeMVL0(pu.mv[REF_PIC_LIST_0]);
  Mv mergeMVL1(pu.mv[REF_PIC_LIST_1]);

  /*Clip the starting MVs*/
  if( !pu.cs->pps->getWrapAroundEnabledFlag() )
  {
    clipMv( mergeMVL0, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps );
    clipMv( mergeMVL1, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps );
  }

  /*L0 MC for refinement*/
  {
    int offset;
#if IF_12TAP
    int leftPixelExtra = (NTAPS_LUMA(0) >> 1) - 1;
#else
    int leftPixelExtra = (NTAPS_LUMA >> 1) - 1;
#endif
    offset = (DMVR_NUM_ITERATION + leftPixelExtra) * (m_cYuvRefBuffDMVRL0.bufs[COMPONENT_Y].stride + 1);
    offset += (-(int)DMVR_NUM_ITERATION)* (int)m_cYuvRefBuffDMVRL0.bufs[COMPONENT_Y].stride;
    offset += (-(int)DMVR_NUM_ITERATION);
    PelBuf srcBuf = m_cYuvRefBuffDMVRL0.bufs[COMPONENT_Y];
    PelUnitBuf yuvPredTempL0 = PelUnitBuf(pu.chromaFormat, PelBuf(m_cYuvPredTempDMVRL0,
      m_biLinearBufStride
      , pu.lwidth() + (2 * DMVR_NUM_ITERATION), pu.lheight() + (2 * DMVR_NUM_ITERATION)));

    xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( REF_PIC_LIST_0, refIdx0 )->unscaledPic, mergeMVL0, yuvPredTempL0, true, clpRngs.comp[COMPONENT_Y],
      false, false, pu.cu->slice->getScalingRatio( REF_PIC_LIST_0, refIdx0 ), pu.lwidth() + ( 2 * DMVR_NUM_ITERATION ), pu.lheight() + ( 2 * DMVR_NUM_ITERATION ), true, ( (Pel *)srcBuf.buf ) + offset, srcBuf.stride );
  }

  /*L1 MC for refinement*/
  {
    int offset;
#if IF_12TAP
    int leftPixelExtra = (NTAPS_LUMA(0) >> 1) - 1;
#else
    int leftPixelExtra = (NTAPS_LUMA >> 1) - 1;
#endif
    offset = (DMVR_NUM_ITERATION + leftPixelExtra) * (m_cYuvRefBuffDMVRL1.bufs[COMPONENT_Y].stride + 1);
    offset += (-(int)DMVR_NUM_ITERATION)* (int)m_cYuvRefBuffDMVRL1.bufs[COMPONENT_Y].stride;
    offset += (-(int)DMVR_NUM_ITERATION);
    PelBuf srcBuf = m_cYuvRefBuffDMVRL1.bufs[COMPONENT_Y];
    PelUnitBuf yuvPredTempL1 = PelUnitBuf(pu.chromaFormat, PelBuf(m_cYuvPredTempDMVRL1,
      m_biLinearBufStride
      , pu.lwidth() + (2 * DMVR_NUM_ITERATION), pu.lheight() + (2 * DMVR_NUM_ITERATION)));

    xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( REF_PIC_LIST_1, refIdx1 )->unscaledPic, mergeMVL1, yuvPredTempL1, true, clpRngs.comp[COMPONENT_Y],
      false, false, pu.cu->slice->getScalingRatio( REF_PIC_LIST_1, refIdx1 ), pu.lwidth() + ( 2 * DMVR_NUM_ITERATION ), pu.lheight() + ( 2 * DMVR_NUM_ITERATION ), true, ( (Pel *)srcBuf.buf ) + offset, srcBuf.stride );
  }
}

void InterPrediction::xProcessDMVR(PredictionUnit& pu, PelUnitBuf &pcYuvDst, const ClpRngs &clpRngs, const bool bioApplied)
{
#if MULTI_PASS_DMVR
  CHECK( true, "DMVR is removed when MULTI_PASS_DMVR is turned on." );
#else
  int iterationCount = 1;
  /*Always High Precision*/
  int mvShift = MV_FRACTIONAL_BITS_INTERNAL;

  /*use merge MV as starting MV*/
  Mv mergeMv[] = { pu.mv[REF_PIC_LIST_0] , pu.mv[REF_PIC_LIST_1] };

  m_biLinearBufStride = (MAX_CU_SIZE + (2 * DMVR_NUM_ITERATION));
  int dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
  int dx = std::min<int>(pu.lumaSize().width,  DMVR_SUBCU_WIDTH);
  Position puPos = pu.lumaPos();

  int bd = pu.cs->slice->getClpRngs().comp[COMPONENT_Y].bd;

  int            bioEnabledThres = 2 * dy * dx;
  bool           bioAppliedType[MAX_NUM_SUBCU_DMVR];

#if JVET_J0090_MEMORY_BANDWITH_MEASURE
  JVET_J0090_SET_CACHE_ENABLE(true);
  for (int k = 0; k < NUM_REF_PIC_LIST_01; k++)
  {
    RefPicList refId = (RefPicList)k;
    const Picture* refPic = pu.cu->slice->getRefPic(refId, pu.refIdx[refId]);
    for (int compID = 0; compID < MAX_NUM_COMPONENT; compID++)
    {
      Mv cMv = pu.mv[refId];
      int mvshiftTemp = mvShift + getComponentScaleX((ComponentID)compID, pu.chromaFormat);
      int filtersize = (compID == (COMPONENT_Y)) ? NTAPS_LUMA : NTAPS_CHROMA;
      cMv += Mv(-(((filtersize >> 1) - 1) << mvshiftTemp), -(((filtersize >> 1) - 1) << mvshiftTemp));
      bool wrapRef = false;
      if ( pu.cs->pps->getWrapAroundEnabledFlag() )
      {
        wrapRef = wrapClipMv(cMv, pu.blocks[0].pos(), pu.blocks[0].size(), pu.cs->sps, pu.cs->pps);
      }
      else
      {
        clipMv(cMv, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);
      }

      int width = pcYuvDst.bufs[compID].width + (filtersize - 1);
      int height = pcYuvDst.bufs[compID].height + (filtersize - 1);

      CPelBuf refBuf;
      Position recOffset = pu.blocks[compID].pos().offset(cMv.getHor() >> mvshiftTemp, cMv.getVer() >> mvshiftTemp);
      refBuf = refPic->getRecoBuf(CompArea((ComponentID)compID, pu.chromaFormat, recOffset, pu.blocks[compID].size()), wrapRef);

      JVET_J0090_SET_REF_PICTURE(refPic, (ComponentID)compID);
      for (int row = 0; row < height; row++)
      {
        for (int col = 0; col < width; col++)
        {
          JVET_J0090_CACHE_ACCESS(((Pel *)refBuf.buf) + row * refBuf.stride + col, __FILE__, __LINE__);
        }
      }
    }
  }
  JVET_J0090_SET_CACHE_ENABLE(false);
#endif

  {
    int num = 0;

    int scaleX = getComponentScaleX(COMPONENT_Cb, pu.chromaFormat);
    int scaleY = getComponentScaleY(COMPONENT_Cb, pu.chromaFormat);
    m_biLinearBufStride = (dx + (2 * DMVR_NUM_ITERATION));
    // point mc buffer to cetre point to avoid multiplication to reach each iteration to the begining
    Pel *biLinearPredL0 = m_cYuvPredTempDMVRL0 + (DMVR_NUM_ITERATION * m_biLinearBufStride) + DMVR_NUM_ITERATION;
    Pel *biLinearPredL1 = m_cYuvPredTempDMVRL1 + (DMVR_NUM_ITERATION * m_biLinearBufStride) + DMVR_NUM_ITERATION;

    PredictionUnit subPu = pu;
    subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(puPos.x, puPos.y, dx, dy)));
    m_cYuvRefBuffDMVRL0 = (pu.chromaFormat == CHROMA_400 ?
      PelUnitBuf(pu.chromaFormat, PelBuf(m_cRefSamplesDMVRL0[0], pcYuvDst.Y())) :
      PelUnitBuf(pu.chromaFormat, PelBuf(m_cRefSamplesDMVRL0[0], pcYuvDst.Y()),
        PelBuf(m_cRefSamplesDMVRL0[1], pcYuvDst.Cb()), PelBuf(m_cRefSamplesDMVRL0[2], pcYuvDst.Cr())));
    m_cYuvRefBuffDMVRL0 = m_cYuvRefBuffDMVRL0.subBuf(UnitAreaRelative(pu, subPu));

    m_cYuvRefBuffDMVRL1 = (pu.chromaFormat == CHROMA_400 ?
      PelUnitBuf(pu.chromaFormat, PelBuf(m_cRefSamplesDMVRL1[0], pcYuvDst.Y())) :
      PelUnitBuf(pu.chromaFormat, PelBuf(m_cRefSamplesDMVRL1[0], pcYuvDst.Y()), PelBuf(m_cRefSamplesDMVRL1[1], pcYuvDst.Cb()),
        PelBuf(m_cRefSamplesDMVRL1[2], pcYuvDst.Cr())));
    m_cYuvRefBuffDMVRL1 = m_cYuvRefBuffDMVRL1.subBuf(UnitAreaRelative(pu, subPu));

    PelUnitBuf srcPred0 = (pu.chromaFormat == CHROMA_400 ?
      PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvDst.Y())) :
      PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[0][0], pcYuvDst.Y()), PelBuf(m_acYuvPred[0][1], pcYuvDst.Cb()), PelBuf(m_acYuvPred[0][2], pcYuvDst.Cr())));
    PelUnitBuf srcPred1 = (pu.chromaFormat == CHROMA_400 ?
      PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvDst.Y())) :
      PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvPred[1][0], pcYuvDst.Y()), PelBuf(m_acYuvPred[1][1], pcYuvDst.Cb()), PelBuf(m_acYuvPred[1][2], pcYuvDst.Cr())));

    srcPred0 = srcPred0.subBuf(UnitAreaRelative(pu, subPu));
    srcPred1 = srcPred1.subBuf(UnitAreaRelative(pu, subPu));

    int yStart = 0;

    for (int y = puPos.y; y < (puPos.y + pu.lumaSize().height); y = y + dy, yStart = yStart + dy)
    {
      for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + dx, xStart = xStart + dx)
      {
        subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, dy)));
        subPu.mv[0] = mergeMv[0];
        subPu.mv[1] = mergeMv[1];
        xPrefetch(subPu, m_cYuvRefBuffDMVRL0, REF_PIC_LIST_0, 1);
        xPrefetch(subPu, m_cYuvRefBuffDMVRL1, REF_PIC_LIST_1, 1);

        xinitMC(subPu, clpRngs);

        uint64_t minCost = MAX_UINT64;
        bool notZeroCost = true;
        int16_t totalDeltaMV[2] = { 0,0 };
        int16_t deltaMV[2] = { 0, 0 };
        uint64_t  *pSADsArray;

        for (int i = 0; i < (((2 * DMVR_NUM_ITERATION) + 1) * ((2 * DMVR_NUM_ITERATION) + 1)); i++)
        {
          m_SADsArray[i] = MAX_UINT64;
        }
        pSADsArray = &m_SADsArray[(((2 * DMVR_NUM_ITERATION) + 1) * ((2 * DMVR_NUM_ITERATION) + 1)) >> 1];

        for (int i = 0; i < iterationCount; i++)
        {
          deltaMV[0] = 0;
          deltaMV[1] = 0;
          Pel *addrL0 = biLinearPredL0 + totalDeltaMV[0] + (totalDeltaMV[1] * m_biLinearBufStride);
          Pel *addrL1 = biLinearPredL1 - totalDeltaMV[0] - (totalDeltaMV[1] * m_biLinearBufStride);
          if (i == 0)
          {
            minCost = xDMVRCost(clpRngs.comp[COMPONENT_Y].bd, addrL0, m_biLinearBufStride, addrL1, m_biLinearBufStride, dx, dy);
            minCost -= (minCost >>2);
            pSADsArray[0] = minCost;
            if (minCost < (dx * dy))
            {
              notZeroCost = false;
              break;
            }
          }
          if (!minCost)
          {
            notZeroCost = false;
            break;
          }

          xBIPMVRefine(bd, addrL0, addrL1, minCost, deltaMV, pSADsArray, dx, dy);

          if (deltaMV[0] == 0 && deltaMV[1] == 0)
          {
            break;
          }
          totalDeltaMV[0] += deltaMV[0];
          totalDeltaMV[1] += deltaMV[1];
          pSADsArray += ((deltaMV[1] * (((2 * DMVR_NUM_ITERATION) + 1))) + deltaMV[0]);
        }

        bioAppliedType[num] = (minCost < bioEnabledThres) ? false : bioApplied;
        totalDeltaMV[0] = (totalDeltaMV[0] << mvShift);
        totalDeltaMV[1] = (totalDeltaMV[1] << mvShift);
        xDMVRSubPixelErrorSurface(notZeroCost, totalDeltaMV, deltaMV, pSADsArray);

        pu.mvdL0SubPu[num] = Mv(totalDeltaMV[0], totalDeltaMV[1]);
        PelUnitBuf subPredBuf = pcYuvDst.subBuf(UnitAreaRelative(pu, subPu));

        bool blockMoved = false;

        if (pu.mvdL0SubPu[num] != Mv(0, 0))
        {
          blockMoved = true;
          if (isChromaEnabled(pu.chromaFormat))
          {
          xPrefetch(subPu, m_cYuvRefBuffDMVRL0, REF_PIC_LIST_0, 0);
          xPrefetch(subPu, m_cYuvRefBuffDMVRL1, REF_PIC_LIST_1, 0);
          }
          xPad(subPu, m_cYuvRefBuffDMVRL0, REF_PIC_LIST_0);
          xPad(subPu, m_cYuvRefBuffDMVRL1, REF_PIC_LIST_1);
        }

        int dstStride[MAX_NUM_COMPONENT] = { pcYuvDst.bufs[COMPONENT_Y].stride,
                                             isChromaEnabled(pu.chromaFormat) ? pcYuvDst.bufs[COMPONENT_Cb].stride : 0,
                                             isChromaEnabled(pu.chromaFormat) ? pcYuvDst.bufs[COMPONENT_Cr].stride : 0};
        subPu.mv[0] = mergeMv[REF_PIC_LIST_0] + pu.mvdL0SubPu[num];
        subPu.mv[1] = mergeMv[REF_PIC_LIST_1] - pu.mvdL0SubPu[num];

        subPu.mv[0].clipToStorageBitDepth();
        subPu.mv[1].clipToStorageBitDepth();

        xFinalPaddedMCForDMVR(subPu, srcPred0, srcPred1, m_cYuvRefBuffDMVRL0, m_cYuvRefBuffDMVRL1, bioAppliedType[num],
                              mergeMv, blockMoved);

        subPredBuf.bufs[COMPONENT_Y].buf = pcYuvDst.bufs[COMPONENT_Y].buf + xStart + yStart * dstStride[COMPONENT_Y];

        if (isChromaEnabled(pu.chromaFormat))
        {
        subPredBuf.bufs[COMPONENT_Cb].buf = pcYuvDst.bufs[COMPONENT_Cb].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cb]);

        subPredBuf.bufs[COMPONENT_Cr].buf = pcYuvDst.bufs[COMPONENT_Cr].buf + (xStart >> scaleX) + ((yStart >> scaleY) * dstStride[COMPONENT_Cr]);
        }
        xWeightedAverage(subPu, srcPred0, srcPred1, subPredBuf, subPu.cu->slice->getSPS()->getBitDepths(), subPu.cu->slice->clpRngs(), bioAppliedType[num]);
        num++;
      }
    }
  }
  JVET_J0090_SET_CACHE_ENABLE(true);
#endif
}

#if JVET_J0090_MEMORY_BANDWITH_MEASURE
void InterPrediction::cacheAssign( CacheModel *cache )
{
  m_cacheModel = cache;
  m_if.cacheAssign( cache );
  m_if.initInterpolationFilter( !cache->isCacheEnable() );
}
#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
void  InterPrediction::sortInterMergeMMVDCandidates(PredictionUnit &pu, MergeCtx& mrgCtx, uint32_t * mmvdLUT, uint32_t MMVDIdx)
{
  
  const int tempNum = (const int) (std::min<int>(MMVD_BASE_MV_NUM, mrgCtx.numValidMergeCand) * MMVD_MAX_REFINE_NUM);
  const int groupSize = std::min<int>(tempNum, ADAPTIVE_SUB_GROUP_SIZE_MMVD);
#if _WINDOWS
  Distortion candCostList[MMVD_BASE_MV_NUM* MMVD_MAX_REFINE_NUM];
#else
  Distortion candCostList[tempNum] ;
#endif
  
  for (uint32_t i = 0; i < tempNum; i++)
  {
    mmvdLUT[i] = i;
    candCostList[i] = MAX_UINT;
  }
  Distortion uiCost;
  DistParam cDistParam;
  cDistParam.applyWeight = false;
  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;
  if (!xAMLGetCurBlkTemplate(pu, nWidth, nHeight))
  {
    return;
  }
  
  int startMMVDIdx = 0;
  int endMMVDIdx = tempNum;
  if(MMVDIdx != -1)
  {
    uint32_t gpId = MMVDIdx/groupSize;
    startMMVDIdx = gpId * groupSize;
    endMMVDIdx = (gpId+1) * groupSize;
  }
  
  int shiftEnc = MMVD_SIZE_SHIFT;
  int encGrpSize = groupSize >> shiftEnc;
  for (int mmvdMergeCand = startMMVDIdx; mmvdMergeCand < endMMVDIdx; mmvdMergeCand++)
  {
    mrgCtx.setMmvdMergeCandiInfo(pu, mmvdMergeCand, mmvdMergeCand);
    
    for (int refList = 0; refList < 2; refList++)
    {
      if (pu.refIdx[refList] >= 0)
      {
        pu.mv[refList].roundToPrecision(MV_PRECISION_QUARTER, MV_PRECISION_INT);
      }
    }
    
    uiCost = 0;
    
    PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
    PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
    
    getBlkAMLRefTemplate(pu, pcBufPredRefTop, pcBufPredRefLeft);
    
    if (m_bAMLTemplateAvailabe[0])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
      
      uiCost += cDistParam.distFunc(cDistParam);
    }
    
    if (m_bAMLTemplateAvailabe[1])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
      
      uiCost += cDistParam.distFunc(cDistParam);
    }
    // update part
    uint32_t i;
    uint32_t shift = 0;
    uint32_t gpIdx = mmvdMergeCand/groupSize;
    uint32_t endIdx = gpIdx * groupSize + encGrpSize;
    while (shift < encGrpSize && uiCost < candCostList[endIdx - 1 - shift])
    {
      shift++;
    }
    if (shift != 0)
    {
      for (i = 1; i < shift; i++)
      {
        mmvdLUT[endIdx - i] = mmvdLUT[endIdx - 1 - i];
        candCostList[endIdx - i] = candCostList[endIdx - 1 - i];
      }
      mmvdLUT[endIdx - shift] = mmvdMergeCand;
      candCostList[endIdx - shift] = uiCost;
    }
  }
  
}
#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
void  InterPrediction::sortAffineMergeCandidates(PredictionUnit pu, AffineMergeCtx& affMrgCtx, uint32_t * affMmvdLUT, uint32_t afMMVDIdx)
{
  const int tempNum = AF_MMVD_NUM;
  int baseIdxToMergeIdxOffset = (int)PU::getMergeIdxFromAfMmvdBaseIdx(affMrgCtx, 0);
  int baseCount               = std::min<int>((int)AF_MMVD_BASE_NUM, affMrgCtx.numValidMergeCand - baseIdxToMergeIdxOffset);
  const int groupSize = std::min<int>(tempNum, ADAPTIVE_SUB_GROUP_SIZE_MMVD_AFF);
  Distortion candCostList[tempNum];
  for (uint32_t i = 0; i < tempNum; i++)
  {
    affMmvdLUT[i] = i;
    candCostList[i] = MAX_UINT;
  }
  
  if (baseCount < 1)
  {
    return;
  }
  Distortion uiCost;
  
  DistParam cDistParam;
  cDistParam.applyWeight = false;
  
  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;
  
  if (!xAMLGetCurBlkTemplate(pu, nWidth, nHeight))
  {
    return;
  }
  
  int startMMVDIdx = 0;
  int endMMVDIdx = tempNum;
  if(afMMVDIdx != -1)
  {
    uint32_t gpId = afMMVDIdx/groupSize;
    startMMVDIdx = gpId * groupSize;
    endMMVDIdx = (gpId+1) * groupSize;
  }
  int shiftEnc = AFFINE_MMVD_SIZE_SHIFT;
  int encGrpSize = groupSize >> shiftEnc;
  for (int mmvdMergeCand = startMMVDIdx; mmvdMergeCand < endMMVDIdx; mmvdMergeCand++)
  {
    pu.afMmvdMergeIdx = (uint8_t)mmvdMergeCand;
    
    int baseIdx = (int)mmvdMergeCand / AF_MMVD_MAX_REFINE_NUM;
    int stepIdx = (int)mmvdMergeCand - baseIdx * AF_MMVD_MAX_REFINE_NUM;
    int dirIdx  = stepIdx % AF_MMVD_OFFSET_DIR;
    stepIdx = stepIdx / AF_MMVD_OFFSET_DIR;
    
    pu.cu->affine = true;
    pu.cu->imv    = IMV_OFF;
    pu.cu->mmvdSkip         = false;
    pu.regularMergeFlag = false;
    pu.mmvdMergeFlag    = false;
    pu.mergeFlag      = true;
    pu.afMmvdFlag     = true;
    pu.afMmvdBaseIdx  = (uint8_t)baseIdx;
    pu.afMmvdDir      = (uint8_t)dirIdx;
    pu.afMmvdStep     = (uint8_t)stepIdx;
    pu.mergeIdx       = (uint8_t)(baseIdxToMergeIdxOffset + baseIdx);
    pu.mergeType = affMrgCtx.mergeType[pu.mergeIdx];
#if INTER_LIC
    pu.cu->LICFlag = affMrgCtx.LICFlags[pu.mergeIdx];
    pu.cu->LICFlag = false;
#endif
    pu.interDir = affMrgCtx.interDirNeighbours[pu.mergeIdx];
    pu.cu->affineType = affMrgCtx.affineType[pu.mergeIdx];
    pu.cu->BcwIdx = affMrgCtx.BcwIdx[pu.mergeIdx];
    pu.ciipFlag = false;
    MvField mvfMmvd[2][3];
    PU::getAfMmvdMvf(pu, affMrgCtx, mvfMmvd, pu.mergeIdx, pu.afMmvdStep, pu.afMmvdDir);
#if JVET_Z0067_RPR_ENABLE
    bool  bIsRefScaled = false;
#endif
    for (int i = 0; i < 2; i++)
    {
      if( pu.cs->slice->getNumRefIdx( RefPicList( i ) ) > 0 )
      {
        pu.mvpIdx[i] = 0;
        pu.mvpNum[i] = 0;
        pu.mvd[i]    = Mv();
        pu.refIdx[i] = mvfMmvd[i][0].refIdx;
        pu.mvAffi[i][0] = mvfMmvd[i][0].mv;
        pu.mvAffi[i][1] = mvfMmvd[i][1].mv;
        pu.mvAffi[i][2] = mvfMmvd[i][2].mv;
      }
#if JVET_Z0067_RPR_ENABLE
      if ( !bIsRefScaled && pu.refIdx[i]>=0 && pu.cu->slice->getRefPic(i ? REF_PIC_LIST_1 : REF_PIC_LIST_0, pu.refIdx[i])->isRefScaled(pu.cs->pps) )
      {
        bIsRefScaled = true;
      }
#endif
    }
    uiCost = 0;
#if JVET_Z0067_RPR_ENABLE
    if ( bIsRefScaled )
    {
      uiCost = std::numeric_limits<Distortion>::max();
    }
    else
    {
#endif
    PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
    PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
    getAffAMLRefTemplate(pu, pcBufPredRefTop, pcBufPredRefLeft);
    
    if (m_bAMLTemplateAvailabe[0])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
      
      uiCost += cDistParam.distFunc(cDistParam);
    }
    
    if (m_bAMLTemplateAvailabe[1])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
      
      uiCost += cDistParam.distFunc(cDistParam);
    }
#if JVET_Z0067_RPR_ENABLE
    }
#endif

    // update part
    uint32_t i;
    uint32_t shift = 0;
    uint32_t gpIdx = mmvdMergeCand/groupSize;
    uint32_t endIdx = gpIdx * groupSize + encGrpSize;
    while (shift < encGrpSize && uiCost < candCostList[endIdx - 1 - shift])
    {
      shift++;
    }
    
    if (shift != 0)
    {
      for (i = 1; i < shift; i++)
      {
        affMmvdLUT[endIdx - i] = affMmvdLUT[endIdx - 1 - i];
        candCostList[endIdx - i] = candCostList[endIdx - 1 - i];
      }
      affMmvdLUT[endIdx - shift] = mmvdMergeCand;
      candCostList[endIdx - shift] = uiCost;
    }
    
  }
}
#endif

#if JVET_W0090_ARMC_TM
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
void InterPrediction::adjustMergeCandidatesInOneCandidateGroup(PredictionUnit &pu, MergeCtx& mvpMergeCandCtx, int numRetrievedMergeCand, int mrgCandIdx)
{
  if (mvpMergeCandCtx.numValidMergeCand <= 1)
  {
    return;
  }

  const int numCandInCategory = std::min(numRetrievedMergeCand, mvpMergeCandCtx.numValidMergeCand);

  uint32_t RdCandList[MRG_MAX_NUM_CANDS];
  Distortion candCostList[MRG_MAX_NUM_CANDS];

  for (uint32_t j = 0; j < MRG_MAX_NUM_CANDS; j++)
  {
    RdCandList[j] = j;
    candCostList[j] = MAX_UINT;
  }

  Distortion uiCost;

  DistParam cDistParam;
  cDistParam.applyWeight = false;

  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;

  auto origMergeIdx = pu.mergeIdx;
  for (uint32_t uiMergeCand = 0; uiMergeCand < mvpMergeCandCtx.numValidMergeCand; uiMergeCand++)
  {
    if (mvpMergeCandCtx.candCost[uiMergeCand] == MAX_UINT64)
    {
      uiCost = 0;

      mvpMergeCandCtx.setMergeInfo(pu, uiMergeCand);

      PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
      PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
      PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
      PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));

#if JVET_Z0067_RPR_ENABLE
    bool bRefIsRescaled = false;
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      const RefPicList eRefPicList = refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
      bRefIsRescaled |= (pu.refIdx[refList] >= 0) ? pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->isRefScaled(pu.cs->pps) : false;
    }
    if (bRefIsRescaled)
    {
      uiCost = std::numeric_limits<Distortion>::max();
    }
    else
    {
#endif
      getBlkAMLRefTemplate(pu, pcBufPredRefTop, pcBufPredRefLeft);

      if (m_bAMLTemplateAvailabe[0])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }

      if (m_bAMLTemplateAvailabe[1])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }
#if JVET_Z0067_RPR_ENABLE
    }
#endif
    }
    else
    {
      uiCost = mvpMergeCandCtx.candCost[uiMergeCand];
    }

    updateCandList(uiMergeCand, uiCost, numCandInCategory, RdCandList, candCostList);
  }
  pu.mergeIdx = origMergeIdx;

  updateCandInOneCandidateGroup(mvpMergeCandCtx, RdCandList, numCandInCategory);

  mvpMergeCandCtx.numValidMergeCand = numCandInCategory;

  for (int idx = 0; idx < numCandInCategory; idx++)
  {
    mvpMergeCandCtx.candCost[idx] = candCostList[idx];
  }
}

#if JVET_Z0102_NO_ARMC_FOR_ZERO_CAND
void InterPrediction::adjustMergeCandidates(PredictionUnit& pu, MergeCtx& mvpMergeCandCtx, int numRetrievedMergeCand)
{
  if (mvpMergeCandCtx.numValidMergeCand <= 1)
  {
    return;
  }

  const int numCandInCategory = std::min(numRetrievedMergeCand, mvpMergeCandCtx.numValidMergeCand);

  uint32_t RdCandList[MRG_MAX_NUM_CANDS];
  Distortion candCostList[MRG_MAX_NUM_CANDS];

  for (uint32_t j = 0; j < MRG_MAX_NUM_CANDS; j++)
  {
    RdCandList[j] = j;
    candCostList[j] = MAX_UINT64;
  }

  Distortion uiCost;

  DistParam cDistParam;
  cDistParam.applyWeight = false;

  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;

  auto origMergeIdx = pu.mergeIdx;

  PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));

  for (uint32_t uiMergeCand = 0; uiMergeCand < mvpMergeCandCtx.numValidMergeCand; uiMergeCand++)
  {
    if (mvpMergeCandCtx.numCandToTestEnc != mvpMergeCandCtx.numValidMergeCand)
    {
      if (uiMergeCand >= mvpMergeCandCtx.numCandToTestEnc)
      {
        mvpMergeCandCtx.candCost[uiMergeCand] = MAX_UINT64 - 1;
      }
    }

    if (mvpMergeCandCtx.candCost[uiMergeCand] == MAX_UINT64)
    {
      uiCost = 0;

      mvpMergeCandCtx.setMergeInfo(pu, uiMergeCand);

      PelUnitBuf pcBufPredRefTop =
        (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
      PelUnitBuf pcBufPredRefLeft =
        (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));

#if JVET_Z0067_RPR_ENABLE
    bool bRefIsRescaled = false;
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      const RefPicList eRefPicList = refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
      bRefIsRescaled |= (pu.refIdx[refList] >= 0) ? pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->isRefScaled(pu.cs->pps) : false;
    }
    if (bRefIsRescaled)
    {
      uiCost = std::numeric_limits<Distortion>::max();
    }
    else
    {
#endif
      getBlkAMLRefTemplate(pu, pcBufPredRefTop, pcBufPredRefLeft);

      if (m_bAMLTemplateAvailabe[0])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(),
          pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }

      if (m_bAMLTemplateAvailabe[1])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(),
          pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }
#if JVET_Z0067_RPR_ENABLE
    }
#endif
    }
    else
    {
      uiCost = mvpMergeCandCtx.candCost[uiMergeCand];
    }

    updateCandList(uiMergeCand, uiCost, mvpMergeCandCtx.numValidMergeCand, RdCandList, candCostList);

  }
  pu.mergeIdx = origMergeIdx;

  updateCandInOneCandidateGroup(mvpMergeCandCtx, RdCandList, mvpMergeCandCtx.numValidMergeCand);

  for (int idx = 0; idx < mvpMergeCandCtx.numValidMergeCand; idx++)
  {
    mvpMergeCandCtx.candCost[idx] = candCostList[idx];
  }

  mvpMergeCandCtx.numValidMergeCand = numCandInCategory;

  for (int idx = 0; idx < numCandInCategory; idx++)
  {
    mvpMergeCandCtx.candCost[idx] = candCostList[idx];
  }
}
#endif

void  InterPrediction::updateCandInOneCandidateGroup(MergeCtx& mrgCtx, uint32_t* RdCandList, int numCandInCategory)
{
  MergeCtx mrgCtxTmp;
  for (uint32_t uiMergeCand = 0; uiMergeCand < mrgCtx.numValidMergeCand; uiMergeCand++)
  {
    mrgCtxTmp.BcwIdx[uiMergeCand] = mrgCtx.BcwIdx[uiMergeCand];
    mrgCtxTmp.interDirNeighbours[uiMergeCand] = mrgCtx.interDirNeighbours[uiMergeCand];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1];
    mrgCtxTmp.useAltHpelIf[uiMergeCand] = mrgCtx.useAltHpelIf[uiMergeCand];
#if INTER_LIC 
    mrgCtxTmp.LICFlags[uiMergeCand] = mrgCtx.LICFlags[uiMergeCand];
#endif
#if MULTI_HYP_PRED
    mrgCtxTmp.addHypNeighbours[uiMergeCand] = mrgCtx.addHypNeighbours[uiMergeCand];
#endif
  }
  //update
  for (uint32_t uiMergeCand = 0; uiMergeCand < numCandInCategory; uiMergeCand++)
  {
    mrgCtx.BcwIdx[uiMergeCand] = mrgCtxTmp.BcwIdx[RdCandList[uiMergeCand]];
    mrgCtx.interDirNeighbours[uiMergeCand] = mrgCtxTmp.interDirNeighbours[RdCandList[uiMergeCand]];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand] << 1)];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand] << 1) + 1];
    mrgCtx.useAltHpelIf[uiMergeCand] = mrgCtxTmp.useAltHpelIf[RdCandList[uiMergeCand]];
#if INTER_LIC
    mrgCtx.LICFlags[uiMergeCand] = mrgCtxTmp.LICFlags[RdCandList[uiMergeCand]];
#endif
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[uiMergeCand] = mrgCtxTmp.addHypNeighbours[RdCandList[uiMergeCand]];
#endif
  }
}
#endif

void  InterPrediction::adjustInterMergeCandidates(PredictionUnit &pu, MergeCtx& mrgCtx, int mrgCandIdx)
{
  uint32_t RdCandList[MRG_MAX_NUM_CANDS][MRG_MAX_NUM_CANDS];
  Distortion candCostList[MRG_MAX_NUM_CANDS][MRG_MAX_NUM_CANDS];

  for (uint32_t i = 0; i < MRG_MAX_NUM_CANDS; i++)
  {
    for (uint32_t j = 0; j < MRG_MAX_NUM_CANDS; j++)
    {
      RdCandList[i][j] = j;
      candCostList[i][j] = MAX_UINT;
    }
  }

  Distortion uiCost;

  DistParam cDistParam;
  cDistParam.applyWeight = false;

  /*const SPS &sps = *pu.cs->sps;
  Position puPos = pu.lumaPos();*/
  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;

  if (!xAMLGetCurBlkTemplate(pu, nWidth, nHeight))
  {
    return;
  }

#if JVET_X0049_ADAPT_DMVR
  uint8_t origMergeIdx = pu.mergeIdx;
#endif
  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE)*ADAPTIVE_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE > mrgCtx.numValidMergeCand)) ? mrgCtx.numValidMergeCand : ((mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE >= mrgCtx.numValidMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    uiCost = 0;

    mrgCtx.setMergeInfo(pu, uiMergeCand);
    PU::spanMotionInfo(pu, mrgCtx);

    PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
    PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));

#if JVET_Y0128_NON_CTC
    bool bRefIsRescaled = false;
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      const RefPicList eRefPicList = refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
      bRefIsRescaled |= (pu.refIdx[refList] >= 0) ? pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->isRefScaled(pu.cs->pps) : false;
    }
    if ( !bRefIsRescaled )
    {
#endif
    getBlkAMLRefTemplate(pu, pcBufPredRefTop, pcBufPredRefLeft);

    if (m_bAMLTemplateAvailabe[0])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

      uiCost += cDistParam.distFunc(cDistParam);
    }

    if (m_bAMLTemplateAvailabe[1])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

      uiCost += cDistParam.distFunc(cDistParam);
    }
#if JVET_Y0128_NON_CTC
    }
#endif

    updateCandList(uiMergeCand, uiCost, ADAPTIVE_SUB_GROUP_SIZE, RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE], candCostList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE]);
  }
#if JVET_X0049_ADAPT_DMVR
  pu.mergeIdx = origMergeIdx;
#else
  pu.mergeIdx = mrgCandIdx;    //restore the merge index
#endif
  updateCandInfo(mrgCtx, RdCandList
    , mrgCandIdx
  );

}
#endif

#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
bool InterPrediction::xAMLGetCurBlkTemplate(PredictionUnit& pu, int nCurBlkWidth, int nCurBlkHeight)
{
  m_bAMLTemplateAvailabe[0] = xAMLIsTopTempAvailable(pu);
  m_bAMLTemplateAvailabe[1] = xAMLIsLeftTempAvailable(pu);

  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    return false;
  }

  /* const int       lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
   const int       horShift = (lumaShift + ::getComponentScaleX(COMPONENT_Y, pu.chromaFormat));
   const int       verShift = (lumaShift + ::getComponentScaleY(COMPONENT_Y, pu.chromaFormat));*/
  const Picture&  currPic = *pu.cs->picture;
  const CPelBuf recBuf = currPic.getRecoBuf(pu.cs->picture->blocks[COMPONENT_Y]);
  std::vector<Pel>& invLUT = m_pcReshape->getInvLUT();

#if JVET_Z0054_BLK_REF_PIC_REORDER
  if(!m_fillCurTplAboveARMC)
#endif
  if (m_bAMLTemplateAvailabe[0])
  {
    const Pel*    rec = recBuf.bufAt(pu.blocks[COMPONENT_Y].pos().offset(0, -AML_MERGE_TEMPLATE_SIZE));
    PelBuf pcYBuf = PelBuf(m_acYuvCurAMLTemplate[0][0], nCurBlkWidth, AML_MERGE_TEMPLATE_SIZE);
    Pel*   pcY = pcYBuf.bufAt(0, 0);
    for (int k = 0; k < nCurBlkWidth; k++)
    {
      for (int l = 0; l < AML_MERGE_TEMPLATE_SIZE; l++)
      {
        int recVal = rec[k + l * recBuf.stride];

        if (m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag())
        {
          recVal = invLUT[recVal];
        }

        pcY[k + l * nCurBlkWidth] = recVal;
      }
    }
#if JVET_Z0054_BLK_REF_PIC_REORDER
    m_fillCurTplAboveARMC = true;
#endif
  }

#if JVET_Z0054_BLK_REF_PIC_REORDER
  if(!m_fillCurTplLeftARMC)
#endif
  if (m_bAMLTemplateAvailabe[1])
  {
    PelBuf pcYBuf = PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nCurBlkHeight);
    Pel*   pcY = pcYBuf.bufAt(0, 0);
    const Pel*    rec = recBuf.bufAt(pu.blocks[COMPONENT_Y].pos().offset(-AML_MERGE_TEMPLATE_SIZE, 0));
    for (int k = 0; k < nCurBlkHeight; k++)
    {
      for (int l = 0; l < AML_MERGE_TEMPLATE_SIZE; l++)
      {
        int recVal = rec[recBuf.stride * k + l];

        if (m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag())
        {
          recVal = invLUT[recVal];
        }

        pcY[AML_MERGE_TEMPLATE_SIZE * k + l] = recVal;
      }
    }
#if JVET_Z0054_BLK_REF_PIC_REORDER
    m_fillCurTplLeftARMC = true;
#endif
  }

  return true;
}

bool InterPrediction::xAMLIsTopTempAvailable(PredictionUnit& pu)
{
  const CodingStructure &cs = *pu.cs;
  Position posRT = pu.Y().topRight();
  const PredictionUnit *puAbove = cs.getPURestricted(posRT.offset(0, -1), pu, pu.chType);

  return (puAbove && pu.cu != puAbove->cu);
}

bool InterPrediction::xAMLIsLeftTempAvailable(PredictionUnit& pu)
{
  const CodingStructure &cs = *pu.cs;
  Position posLB = pu.Y().bottomLeft();
  const PredictionUnit *puLeft = cs.getPURestricted(posLB.offset(-1, 0), pu, pu.chType);

  return (puLeft && pu.cu != puLeft->cu);
}
#endif

#if JVET_W0090_ARMC_TM
void InterPrediction::updateCandList(uint32_t uiCand, Distortion uiCost, uint32_t uiMrgCandNum, uint32_t* RdCandList, Distortion* CandCostList)
{
  uint32_t i;
  uint32_t shift = 0;

  while (shift < uiMrgCandNum && uiCost < CandCostList[uiMrgCandNum - 1 - shift])
  {
    shift++;
  }

  if (shift != 0)
  {
    for (i = 1; i < shift; i++)
    {
      RdCandList[uiMrgCandNum - i] = RdCandList[uiMrgCandNum - 1 - i];
      CandCostList[uiMrgCandNum - i] = CandCostList[uiMrgCandNum - 1 - i];
    }
    RdCandList[uiMrgCandNum - shift] = uiCand;
    CandCostList[uiMrgCandNum - shift] = uiCost;
  }
}
void  InterPrediction::updateCandInfo(MergeCtx& mrgCtx, uint32_t(*RdCandList)[MRG_MAX_NUM_CANDS], int mrgCandIdx)
{
  MergeCtx mrgCtxTmp;
  for (uint32_t ui = 0; ui < MRG_MAX_NUM_CANDS; ++ui)
  {
    mrgCtxTmp.BcwIdx[ui] = BCW_DEFAULT;
    mrgCtxTmp.interDirNeighbours[ui] = 0;
    mrgCtxTmp.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtxTmp.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtxTmp.useAltHpelIf[ui] = false;
#if INTER_LIC
    mrgCtxTmp.LICFlags[ui] = false;
#endif
#if MULTI_HYP_PRED
    mrgCtxTmp.addHypNeighbours[ui].clear();
#endif
  }
  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE)*ADAPTIVE_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE > mrgCtx.numValidMergeCand)) ? mrgCtx.numValidMergeCand : ((mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE >= mrgCtx.numValidMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    mrgCtxTmp.BcwIdx[uiMergeCand] = mrgCtx.BcwIdx[uiMergeCand];
    mrgCtxTmp.interDirNeighbours[uiMergeCand] = mrgCtx.interDirNeighbours[uiMergeCand];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1];
    mrgCtxTmp.useAltHpelIf[uiMergeCand] = mrgCtx.useAltHpelIf[uiMergeCand];
#if INTER_LIC 
    mrgCtxTmp.LICFlags[uiMergeCand] = mrgCtx.LICFlags[uiMergeCand];
#endif
#if MULTI_HYP_PRED
    mrgCtxTmp.addHypNeighbours[uiMergeCand] = mrgCtx.addHypNeighbours[uiMergeCand];
#endif
  }
  //update
  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE)*ADAPTIVE_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE > mrgCtx.numValidMergeCand)) ? mrgCtx.numValidMergeCand : ((mrgCandIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE >= mrgCtx.numValidMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }

    mrgCtx.BcwIdx[uiMergeCand] = mrgCtxTmp.BcwIdx[RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_SUB_GROUP_SIZE]];
    mrgCtx.interDirNeighbours[uiMergeCand] = mrgCtxTmp.interDirNeighbours[RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_SUB_GROUP_SIZE]];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_SUB_GROUP_SIZE] << 1)];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_SUB_GROUP_SIZE] << 1) + 1];
    mrgCtx.useAltHpelIf[uiMergeCand] = mrgCtxTmp.useAltHpelIf[RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_SUB_GROUP_SIZE]];
#if INTER_LIC
    mrgCtx.LICFlags[uiMergeCand] = mrgCtxTmp.LICFlags[RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_SUB_GROUP_SIZE]];
#endif
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[uiMergeCand] = mrgCtxTmp.addHypNeighbours[RdCandList[uiMergeCand / ADAPTIVE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_SUB_GROUP_SIZE]];
#endif
  }
}
#endif

#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
void InterPrediction::getBlkAMLRefTemplate(PredictionUnit &pu, PelUnitBuf &pcBufPredRefTop, PelUnitBuf &pcBufPredRefLeft)
{
  Mv mvCurr;
  const int lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
  const int horShift  = (lumaShift + ::getComponentScaleX(COMPONENT_Y, pu.chromaFormat));
  const int verShift  = (lumaShift + ::getComponentScaleY(COMPONENT_Y, pu.chromaFormat));

  if (xCheckIdenticalMotion(pu))
  {
    mvCurr = pu.mv[0];
    /*const int horIntMv = (mvCurr.getHor() + ((1 << horShift) >> 1)) >> horShift;
    const int verIntMv = (mvCurr.getVer() + ((1 << verShift) >> 1)) >> verShift;
    Mv    subPelMv(horIntMv << horShift, verIntMv << verShift);*/
    Mv subPelMv = mvCurr;
    clipMv(mvCurr, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);
    CHECK(pu.refIdx[0] < 0, "invalid ref idx");

    if (m_bAMLTemplateAvailabe[0])
    {
      Mv mvTop(0, -(AML_MERGE_TEMPLATE_SIZE << verShift));
      mvTop += subPelMv;

      clipMv(mvTop, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);

#if RPR_ENABLE
      const Picture* picRef = pu.cu->slice->getRefPic(REF_PIC_LIST_0, pu.refIdx[0])->unscaledPic;
      const std::pair<int, int>& scalingRatio = pu.cu->slice->getScalingRatio(REF_PIC_LIST_0, pu.refIdx[0]);
#if INTER_LIC
      xPredInterBlk( COMPONENT_Y, pu, picRef, mvTop, pcBufPredRefTop,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
      xPredInterBlk( COMPONENT_Y, pu, picRef, mvTop, pcBufPredRefTop,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true );
#endif
#else
#if INTER_LIC
      xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( REF_PIC_LIST_0, pu.refIdx[0] ), mvTop, pcBufPredRefTop,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true, true,
                     mvCurr );
#else
      xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( REF_PIC_LIST_0, pu.refIdx[0] ), mvTop, pcBufPredRefTop,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true );
#endif
#endif
    }
    if (m_bAMLTemplateAvailabe[1])
    {
      Mv mvLeft(-(AML_MERGE_TEMPLATE_SIZE << horShift), 0);
      mvLeft += subPelMv;

      clipMv(mvLeft, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);

#if RPR_ENABLE
      const Picture* picRef = pu.cu->slice->getRefPic(REF_PIC_LIST_0, pu.refIdx[0])->unscaledPic;
      const std::pair<int, int>& scalingRatio = pu.cu->slice->getScalingRatio(REF_PIC_LIST_0, pu.refIdx[0]);
#if INTER_LIC
      xPredInterBlk( COMPONENT_Y, pu, picRef, mvLeft, pcBufPredRefLeft,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
      xPredInterBlk( COMPONENT_Y, pu, picRef, mvLeft, pcBufPredRefLeft,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true );
#endif
#else
#if INTER_LIC
      xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( REF_PIC_LIST_0, pu.refIdx[0] ), mvLeft, pcBufPredRefLeft,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true, true,
                     mvCurr );
#else
      xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( REF_PIC_LIST_0, pu.refIdx[0] ), mvLeft, pcBufPredRefLeft,
                     false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true );
#endif
#endif
    }
  }
  else
  {
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      if (pu.refIdx[refList] < 0)
      {
        continue;
      }
      RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
      CHECK(pu.refIdx[refList] >= pu.cu->slice->getNumRefIdx(eRefPicList), "Invalid reference index");

      m_iRefListIdx = refList;
      mvCurr        = pu.mv[refList];
      /*const int horIntMv = (mvCurr.getHor() + ((1 << horShift) >> 1)) >> horShift;
      const int verIntMv = (mvCurr.getVer() + ((1 << verShift) >> 1)) >> verShift;
      Mv    subPelMv(horIntMv << horShift, verIntMv << verShift);*/
      Mv subPelMv = mvCurr;
      clipMv(mvCurr, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);

      if (m_bAMLTemplateAvailabe[0])
      {
        Mv mvTop(0, -(AML_MERGE_TEMPLATE_SIZE << verShift));
        mvTop += subPelMv;

        clipMv(mvTop, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);

        PelUnitBuf pcMbBuf =
          PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplate[refList][0], pcBufPredRefTop.Y()));

        if (pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0)
        {
#if RPR_ENABLE
          const Picture* picRef = pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->unscaledPic;
          const std::pair<int, int>& scalingRatio = pu.cu->slice->getScalingRatio(eRefPicList, pu.refIdx[refList]);
#if INTER_LIC
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvTop, pcMbBuf, true,
                         pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvTop, pcMbBuf, true,
                         pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true );
#endif
#else
#if INTER_LIC
          xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( eRefPicList, pu.refIdx[refList] ), mvTop, pcMbBuf, true,
                         pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true, true,
                         mvCurr );
#else
          xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( eRefPicList, pu.refIdx[refList] ), mvTop, pcMbBuf, true,
                         pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true );
#endif
#endif
        }
        else
        {
#if RPR_ENABLE
          const Picture* picRef = pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->unscaledPic;
          const std::pair<int, int>& scalingRatio = pu.cu->slice->getScalingRatio(eRefPicList, pu.refIdx[refList]);
#if INTER_LIC
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvTop, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvTop, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true );
#endif
#else
#if INTER_LIC
          xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( eRefPicList, pu.refIdx[refList] ), mvTop, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
          xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( eRefPicList, pu.refIdx[refList] ), mvTop, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true );
#endif
#endif
        }
      }
      if (m_bAMLTemplateAvailabe[1])
      {
        Mv mvLeft(-(AML_MERGE_TEMPLATE_SIZE << horShift), 0);
        mvLeft += subPelMv;

        clipMv(mvLeft, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);

        PelUnitBuf pcMbBuf =
          PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplate[refList][0], pcBufPredRefLeft.Y()));

        if (pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0)
        {
#if RPR_ENABLE
          const Picture* picRef = pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->unscaledPic;
          const std::pair<int, int>& scalingRatio = pu.cu->slice->getScalingRatio( eRefPicList, pu.refIdx[refList] );
#if INTER_LIC
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvLeft, pcMbBuf,
                         true, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvLeft, pcMbBuf,
                         true, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true );
#endif
#else
  #if INTER_LIC
          xPredInterBlk(COMPONENT_Y, pu, pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList]), mvLeft, pcMbBuf,
                        true, pu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X, 0, 0, false, NULL, 0, true, true, mvCurr);
  #else
          xPredInterBlk(COMPONENT_Y, pu, pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList]), mvLeft, pcMbBuf,
                        true, pu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X, 0, 0, false, NULL, 0, true);
  #endif
#endif
        }
        else
        {
#if RPR_ENABLE
          const Picture* picRef = pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->unscaledPic;
          const std::pair<int, int>& scalingRatio = pu.cu->slice->getScalingRatio(eRefPicList, pu.refIdx[refList]);
#if INTER_LIC
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvLeft, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
          xPredInterBlk( COMPONENT_Y, pu, picRef, mvLeft, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, scalingRatio, 0, 0, false, NULL, 0, true );
#endif
#else
#if INTER_LIC
          xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( eRefPicList, pu.refIdx[refList] ), mvLeft, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true, true, mvCurr );
#else
          xPredInterBlk( COMPONENT_Y, pu, pu.cu->slice->getRefPic( eRefPicList, pu.refIdx[refList] ), mvLeft, pcMbBuf,
                         false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false, SCALE_1X, 0, 0, false, NULL, 0, true );
#endif
#endif
        }
      }
    }
    if (m_bAMLTemplateAvailabe[0])
    {
      CPelUnitBuf srcPred0 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplate[0][0], pcBufPredRefTop.Y()));
      CPelUnitBuf srcPred1 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplate[1][0], pcBufPredRefTop.Y()));
      xWeightedAverageY(pu, srcPred0, srcPred1, pcBufPredRefTop, pu.cu->slice->getSPS()->getBitDepths(),
                        pu.cu->slice->clpRngs());
    }
    if (m_bAMLTemplateAvailabe[1])
    {
      CPelUnitBuf srcPred0 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplate[0][0], pcBufPredRefLeft.Y()));
      CPelUnitBuf srcPred1 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplate[1][0], pcBufPredRefLeft.Y()));
      xWeightedAverageY(pu, srcPred0, srcPred1, pcBufPredRefLeft, pu.cu->slice->getSPS()->getBitDepths(),
                        pu.cu->slice->clpRngs());
    }
  }
}
#endif

#if JVET_W0090_ARMC_TM
void  InterPrediction::adjustAffineMergeCandidates(PredictionUnit &pu, AffineMergeCtx& affMrgCtx, int mrgCandIdx
#if JVET_Z0139_NA_AFF
    , int sortedCandNum
#endif
)
{
#if JVET_Z0139_NA_AFF
  const uint32_t maxNumAffineMergeCand = (sortedCandNum > 0) ? sortedCandNum: pu.cs->slice->getPicHeader()->getMaxNumAffineMergeCand();
#endif
  uint32_t RdCandList[AFFINE_MRG_MAX_NUM_CANDS][AFFINE_MRG_MAX_NUM_CANDS];
  Distortion candCostList[AFFINE_MRG_MAX_NUM_CANDS][AFFINE_MRG_MAX_NUM_CANDS];

  for (uint32_t i = 0; i < AFFINE_MRG_MAX_NUM_CANDS; i++)
  {
    for (uint32_t j = 0; j < AFFINE_MRG_MAX_NUM_CANDS; j++)
    {
      RdCandList[i][j] = j;
      candCostList[i][j] = MAX_UINT;
    }
  }

  Distortion uiCost;

  DistParam cDistParam;
  cDistParam.applyWeight = false;

  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;

  if (!xAMLGetCurBlkTemplate(pu, nWidth, nHeight))
  {
    return;
  }

  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE > affMrgCtx.maxNumMergeCand)) ? affMrgCtx.maxNumMergeCand : ((mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE >= affMrgCtx.maxNumMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    uiCost = 0;

    // set merge information
    pu.interDir = affMrgCtx.interDirNeighbours[uiMergeCand];
    pu.mergeFlag = true;
    pu.regularMergeFlag = false;
    pu.mergeIdx = uiMergeCand;
    pu.cu->affine = true;
    pu.cu->affineType = affMrgCtx.affineType[uiMergeCand];
#if AFFINE_MMVD
    pu.afMmvdFlag = false;
#endif
    pu.cu->BcwIdx = affMrgCtx.BcwIdx[uiMergeCand];
#if INTER_LIC
    pu.cu->LICFlag = affMrgCtx.LICFlags[uiMergeCand];
#endif

    pu.mergeType = affMrgCtx.mergeType[uiMergeCand];
#if JVET_Z0139_NA_AFF
    if (pu.mergeType == MRG_TYPE_DEFAULT_N)
#else
    if (pu.mergeType == MRG_TYPE_SUBPU_ATMVP)
    {
      pu.refIdx[0] = affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0][0].refIdx;
      pu.refIdx[1] = affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1][0].refIdx;
      PU::spanMotionInfo(pu, *affMrgCtx.mrgCtx);
    }
    else
#endif
    {
      for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
      {
        for (int i = 0; i < 3; i++)
        {
          pu.mvAffi[refList][i] = affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + refList][i].mv;
        }
        pu.refIdx[refList] = affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + refList][0].refIdx;
      }

      PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
      PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
      PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
      PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));

#if RPR_ENABLE
      bool bRefIsRescaled = false;
      for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
      {
        const RefPicList eRefPicList = refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
#if JVET_Z0118_GDR             
        Picture *refPic = pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList]);
        if (refPic)
        {
          bRefIsRescaled |= (pu.refIdx[refList] >= 0) ? pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->isRefScaled(pu.cs->pps) : false;
        }
#else
        bRefIsRescaled |= (pu.refIdx[refList] >= 0) ? pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[refList])->isRefScaled(pu.cs->pps) : false;
#endif        
      }

      if ( !bRefIsRescaled )
      {
#endif

      getAffAMLRefTemplate(pu, pcBufPredRefTop, pcBufPredRefLeft);

      if (m_bAMLTemplateAvailabe[0])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }

      if (m_bAMLTemplateAvailabe[1])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }

#if RPR_ENABLE
      }
#endif
    }
#if JVET_Z0139_NA_AFF
    updateCandList(uiMergeCand, uiCost, maxNumAffineMergeCand, RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE], candCostList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE]);
#else
    updateCandList(uiMergeCand, uiCost, ADAPTIVE_AFFINE_SUB_GROUP_SIZE, RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE], candCostList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE]);
#endif
  }
  pu.mergeIdx = mrgCandIdx;    //restore the merge index
  updateAffineCandInfo(pu, affMrgCtx, RdCandList
    , mrgCandIdx
  );

}
void  InterPrediction::updateAffineCandInfo(PredictionUnit &pu, AffineMergeCtx& affMrgCtx, uint32_t(*RdCandList)[AFFINE_MRG_MAX_NUM_CANDS], int mrgCandIdx)
{
  AffineMergeCtx affMrgCtxTmp;
  const uint32_t maxNumAffineMergeCand = pu.cs->slice->getPicHeader()->getMaxNumAffineMergeCand();
  for (int i = 0; i < maxNumAffineMergeCand; i++)
  {
    for (int mvNum = 0; mvNum < 3; mvNum++)
    {
      affMrgCtxTmp.mvFieldNeighbours[(i << 1) + 0][mvNum].setMvField(Mv(), -1);
      affMrgCtxTmp.mvFieldNeighbours[(i << 1) + 1][mvNum].setMvField(Mv(), -1);
    }
    affMrgCtxTmp.interDirNeighbours[i] = 0;
    affMrgCtxTmp.affineType[i] = AFFINEMODEL_4PARAM;
    affMrgCtxTmp.mergeType[i] = MRG_TYPE_DEFAULT_N;
    affMrgCtxTmp.BcwIdx[i] = BCW_DEFAULT;
#if INTER_LIC
    affMrgCtxTmp.LICFlags[i] = false;
#endif
  }
  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE > affMrgCtx.maxNumMergeCand)) ? affMrgCtx.maxNumMergeCand : ((mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE >= affMrgCtx.maxNumMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    for (int mvNum = 0; mvNum < 3; mvNum++)
    {
      affMrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1) + 0][mvNum] = affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0][mvNum];
      affMrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1) + 1][mvNum] = affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1][mvNum];
    }
    affMrgCtxTmp.interDirNeighbours[uiMergeCand] = affMrgCtx.interDirNeighbours[uiMergeCand];
    affMrgCtxTmp.affineType[uiMergeCand] = affMrgCtx.affineType[uiMergeCand];
    affMrgCtxTmp.mergeType[uiMergeCand] = affMrgCtx.mergeType[uiMergeCand];
    affMrgCtxTmp.BcwIdx[uiMergeCand] = affMrgCtx.BcwIdx[uiMergeCand];
#if INTER_LIC                                                   
    affMrgCtxTmp.LICFlags[uiMergeCand] = affMrgCtx.LICFlags[uiMergeCand];
#endif
  }
  //update
  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE > affMrgCtx.maxNumMergeCand)) ? affMrgCtx.maxNumMergeCand : ((mrgCandIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE >= affMrgCtx.maxNumMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    for (int mvNum = 0; mvNum < 3; mvNum++)
    {
      affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0][mvNum] = affMrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_AFFINE_SUB_GROUP_SIZE] << 1) + 0][mvNum];
      affMrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1][mvNum] = affMrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_AFFINE_SUB_GROUP_SIZE] << 1) + 1][mvNum];
    }
    affMrgCtx.interDirNeighbours[uiMergeCand] = affMrgCtxTmp.interDirNeighbours[RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_AFFINE_SUB_GROUP_SIZE]];
    affMrgCtx.affineType[uiMergeCand] = affMrgCtxTmp.affineType[RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_AFFINE_SUB_GROUP_SIZE]];
    affMrgCtx.mergeType[uiMergeCand] = affMrgCtxTmp.mergeType[RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_AFFINE_SUB_GROUP_SIZE]];
    affMrgCtx.BcwIdx[uiMergeCand] = affMrgCtxTmp.BcwIdx[RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_AFFINE_SUB_GROUP_SIZE]];
#if INTER_LIC 
    affMrgCtx.LICFlags[uiMergeCand] = affMrgCtxTmp.LICFlags[RdCandList[uiMergeCand / ADAPTIVE_AFFINE_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_AFFINE_SUB_GROUP_SIZE]];
#endif
  }
}
void InterPrediction::xGetSublkAMLTemplate(const CodingUnit& cu,
  const ComponentID compID,
  const Picture&    refPic,
  const Mv&         mv,
  const int         sublkWidth,
  const int         sublkHeight,
  const int         posW,
  const int         posH,
  int*              numTemplate,
  Pel*              refLeftTemplate,
  Pel*              refAboveTemplate
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
     , bool afMMVD
#endif
                                           )
{
  const int       bitDepth = cu.cs->sps->getBitDepth(toChannelType(compID));
  const int       precShift = std::max(0, bitDepth - 12);

  const CodingUnit* const cuAbove = cu.cs->getCU(cu.blocks[compID].pos().offset(0, -1), toChannelType(compID));
  const CodingUnit* const cuLeft = cu.cs->getCU(cu.blocks[compID].pos().offset(-1, 0), toChannelType(compID));
  const CPelBuf refBuf = cuAbove || cuLeft ? refPic.getRecoBuf(refPic.blocks[compID]) : CPelBuf();

  // above
  if (cuAbove && posH == 0)
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    xGetPredBlkTpl<true>(cu, compID, refBuf, mv, posW, posH, sublkWidth, refAboveTemplate, afMMVD);
#else
    xGetPredBlkTpl<true>(cu, compID, refBuf, mv, posW, posH, sublkWidth, refAboveTemplate);
#endif

    for (int k = posW; k < posW + sublkWidth; k++)
    {
      int refVal = refAboveTemplate[k];
      refVal >>= precShift;
      refAboveTemplate[k] = refVal;
      numTemplate[0]++;
    }
  }

  // left
  if (cuLeft && posW == 0)
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    xGetPredBlkTpl<false>(cu, compID, refBuf, mv, posW, posH, sublkHeight, refLeftTemplate, afMMVD);
#else
    xGetPredBlkTpl<false>(cu, compID, refBuf, mv, posW, posH, sublkHeight, refLeftTemplate);
#endif

    for (int k = posH; k < posH + sublkHeight; k++)
    {
      int refVal = refLeftTemplate[k];
      refVal >>= precShift;
      refLeftTemplate[k] = refVal;
      numTemplate[1]++;
    }
  }
}
void InterPrediction::getAffAMLRefTemplate(PredictionUnit &pu, PelUnitBuf &pcBufPredRefTop, PelUnitBuf &pcBufPredRefLeft)
{
#if INTER_LIC
  int LICshift[2] = { 0 };
  int scale[2]    = { 0 };
  int offset[2]   = { 0 };
#endif
  const int bitDepth = pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA);
  if (xCheckIdenticalMotion(pu))
  {
    Pel *            refLeftTemplate  = m_acYuvRefAMLTemplate[1][0];
    Pel *            refAboveTemplate = m_acYuvRefAMLTemplate[0][0];
    int              numTemplate[2]   = { 0, 0 };   // 0:Above, 1:Left
    const RefPicList eRefPicList      = REF_PIC_LIST_0;
#if JVET_Z0067_RPR_ENABLE
    CHECK(pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[eRefPicList])->isRefScaled(pu.cs->pps), "getAffAMLRefTemplate not supported with ref scaled.");
#endif
    xPredAffineTpl(pu, eRefPicList, numTemplate, refLeftTemplate, refAboveTemplate);
#if INTER_LIC
    if (pu.cu->LICFlag)
    {
      Pel *recLeftTemplate  = m_acYuvCurAMLTemplate[1][0];
      Pel *recAboveTemplate = m_acYuvCurAMLTemplate[0][0];
      xGetLICParamGeneral(*pu.cu, COMPONENT_Y, numTemplate, refLeftTemplate, refAboveTemplate, recLeftTemplate,
                          recAboveTemplate, LICshift[0], scale[0], offset[0]);
      if (m_bAMLTemplateAvailabe[0])
      {
        PelBuf &      dstBuf = pcBufPredRefTop.bufs[0];
        const ClpRng &clpRng = pu.cu->cs->slice->clpRng(COMPONENT_Y);
        dstBuf.linearTransform(scale[0], LICshift[0], offset[0], true, clpRng);
      }
      if (m_bAMLTemplateAvailabe[1])
      {
        PelBuf &      dstBuf = pcBufPredRefLeft.bufs[0];
        const ClpRng &clpRng = pu.cu->cs->slice->clpRng(COMPONENT_Y);
        dstBuf.linearTransform(scale[0], LICshift[0], offset[0], true, clpRng);
      }
    }
#endif
  }
  else
  {
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      if (pu.refIdx[refList] < 0)
      {
        continue;
      }
      RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
      CHECK(pu.refIdx[refList] >= pu.cu->slice->getNumRefIdx(eRefPicList), "Invalid reference index");
#if JVET_Z0067_RPR_ENABLE
      CHECK(pu.cu->slice->getRefPic(eRefPicList, pu.refIdx[eRefPicList])->isRefScaled(pu.cs->pps), "getAffAMLRefTemplate not supported with ref scaled.");
#endif
      Pel *refLeftTemplate  = m_acYuvRefLeftTemplate[refList][0];
      Pel *refAboveTemplate = m_acYuvRefAboveTemplate[refList][0];
      int  numTemplate[2]   = { 0, 0 };   // 0:Above, 1:Left
      xPredAffineTpl(pu, eRefPicList, numTemplate, refLeftTemplate, refAboveTemplate);
#if INTER_LIC
      if (pu.cu->LICFlag)
      {
        Pel *recLeftTemplate  = m_acYuvCurAMLTemplate[1][0];
        Pel *recAboveTemplate = m_acYuvCurAMLTemplate[0][0];
        xGetLICParamGeneral(*pu.cu, COMPONENT_Y, numTemplate, refLeftTemplate, refAboveTemplate, recLeftTemplate,
                            recAboveTemplate, LICshift[refList], scale[refList], offset[refList]);
      }
#endif
    }
    if (m_bAMLTemplateAvailabe[0])
    {
      PelUnitBuf srcPred[2];
      srcPred[0] = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplate[0][0], pcBufPredRefTop.Y()));
      srcPred[1] = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplate[1][0], pcBufPredRefTop.Y()));
#if INTER_LIC
      if (pu.cu->LICFlag)
      {
        for (int i = 0; i < 2; i++)
        {
          if (pu.refIdx[i] >= 0)
          {
            PelBuf &      dstBuf = srcPred[i].bufs[0];
            const ClpRng &clpRng = pu.cu->cs->slice->clpRng(COMPONENT_Y);
            dstBuf.linearTransform(scale[i], LICshift[i], offset[i], true, clpRng);
          }
        }
      }
#endif
      const int iRefIdx0 = pu.refIdx[0];
      const int iRefIdx1 = pu.refIdx[1];
      if (iRefIdx0 >= 0 && iRefIdx1 >= 0)
      {
        for (int i = 0; i < 2; i++)
        {
          PelBuf &  dstBuf   = srcPred[i].bufs[0];
          const int biShift  = IF_INTERNAL_PREC - bitDepth;
          const Pel biOffset = -IF_INTERNAL_OFFS;
          ClpRng    clpRngDummy;
          dstBuf.linearTransform(1, -biShift, biOffset, false, clpRngDummy);
        }
      }
      xWeightedAverageY(pu, srcPred[0], srcPred[1], pcBufPredRefTop, pu.cu->slice->getSPS()->getBitDepths(),
                        pu.cu->slice->clpRngs());
    }
    if (m_bAMLTemplateAvailabe[1])
    {
      PelUnitBuf srcPred[2];
      srcPred[0] = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplate[0][0], pcBufPredRefLeft.Y()));
      srcPred[1] = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplate[1][0], pcBufPredRefLeft.Y()));
#if INTER_LIC
      if (pu.cu->LICFlag)
      {
        for (int i = 0; i < 2; i++)
        {
          if (pu.refIdx[i] >= 0)
          {
            PelBuf &      dstBuf = srcPred[i].bufs[0];
            const ClpRng &clpRng = pu.cu->cs->slice->clpRng(COMPONENT_Y);
            dstBuf.linearTransform(scale[i], LICshift[i], offset[i], true, clpRng);
          }
        }
      }
#endif
      const int iRefIdx0 = pu.refIdx[0];
      const int iRefIdx1 = pu.refIdx[1];
      if (iRefIdx0 >= 0 && iRefIdx1 >= 0)
      {
        for (int i = 0; i < 2; i++)
        {
          PelBuf &  dstBuf   = srcPred[i].bufs[0];
          const int biShift  = IF_INTERNAL_PREC - bitDepth;
          const Pel biOffset = -IF_INTERNAL_OFFS;
          ClpRng    clpRngDummy;
          dstBuf.linearTransform(1, -biShift, biOffset, false, clpRngDummy);
        }
      }
      xWeightedAverageY(pu, srcPred[0], srcPred[1], pcBufPredRefLeft, pu.cu->slice->getSPS()->getBitDepths(),
                        pu.cu->slice->clpRngs());
    }
  }
}
#if JVET_Y0058_IBC_LIST_MODIFY
void  InterPrediction::adjustIBCMergeCandidates(PredictionUnit &pu, MergeCtx& mrgCtx, int mrgCandIdx)
{
#if JVET_Z0084_IBC_TM
  if (mrgCtx.numValidMergeCand <= 1)
  {
    return;
  }
#endif

  uint32_t RdCandList[IBC_MRG_MAX_NUM_CANDS][IBC_MRG_MAX_NUM_CANDS];
  Distortion candCostList[IBC_MRG_MAX_NUM_CANDS][IBC_MRG_MAX_NUM_CANDS];

  for (uint32_t i = 0; i < IBC_MRG_MAX_NUM_CANDS; i++)
  {
    for (uint32_t j = 0; j < IBC_MRG_MAX_NUM_CANDS; j++)
    {
      RdCandList[i][j] = j;
      candCostList[i][j] = MAX_UINT;
    }
  }

  Distortion uiCost;

  DistParam cDistParam;
  cDistParam.applyWeight = false;

  /*const SPS &sps = *pu.cs->sps;
  Position puPos = pu.lumaPos();*/
  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;

  if (!xAMLIBCGetCurBlkTemplate(pu, nWidth, nHeight))
  {
    return;
  }

  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE)*ADAPTIVE_IBC_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE > mrgCtx.numValidMergeCand)) ? mrgCtx.numValidMergeCand : ((mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE >= mrgCtx.numValidMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    uiCost = 0;

    mrgCtx.setMergeInfo(pu, uiMergeCand);

    if (pu.bv == Mv(0, 0))
    {
      break;
    }

    PelBuf pcBufPredRefTop = PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE);
    PelBuf pcBufPredCurTop = PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE);
    PelBuf pcBufPredRefLeft = PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight);
    PelBuf pcBufPredCurLeft = PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight);

    getIBCAMLRefTemplate(pu, nWidth, nHeight);

    if (m_bAMLTemplateAvailabe[0])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop, pcBufPredRefTop, pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

      uiCost += cDistParam.distFunc(cDistParam);
    }

    if (m_bAMLTemplateAvailabe[1])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft, pcBufPredRefLeft, pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

      uiCost += cDistParam.distFunc(cDistParam);
    }
    updateCandList(uiMergeCand, uiCost, ADAPTIVE_IBC_SUB_GROUP_SIZE, RdCandList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE], candCostList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE]);
  }

  updateIBCCandInfo(pu, mrgCtx, RdCandList
    , mrgCandIdx
  );
  pu.mergeIdx = mrgCandIdx;    //restore the merge index
}
void  InterPrediction::updateIBCCandInfo(PredictionUnit &pu, MergeCtx& mrgCtx, uint32_t(*RdCandList)[IBC_MRG_MAX_NUM_CANDS], int mrgCandIdx)
{
  MergeCtx mrgCtxTmp;
  for (uint32_t ui = 0; ui < IBC_MRG_MAX_NUM_CANDS; ++ui)
  {
    mrgCtxTmp.BcwIdx[ui] = BCW_DEFAULT;
    mrgCtxTmp.interDirNeighbours[ui] = 0;
    mrgCtxTmp.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtxTmp.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtxTmp.useAltHpelIf[ui] = false;
#if INTER_LIC
    mrgCtxTmp.LICFlags[ui] = false;
#endif
  }
  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE)*ADAPTIVE_IBC_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE > mrgCtx.numValidMergeCand)) ? mrgCtx.numValidMergeCand : ((mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE >= mrgCtx.numValidMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    mrgCtx.setMergeInfo(pu, uiMergeCand);
    if (pu.bv == Mv(0, 0))
    {
      break;
    }
    mrgCtxTmp.BcwIdx[uiMergeCand] = mrgCtx.BcwIdx[uiMergeCand];
    mrgCtxTmp.interDirNeighbours[uiMergeCand] = mrgCtx.interDirNeighbours[uiMergeCand];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1];
    mrgCtxTmp.useAltHpelIf[uiMergeCand] = mrgCtx.useAltHpelIf[uiMergeCand];
#if INTER_LIC 
    mrgCtxTmp.LICFlags[uiMergeCand] = mrgCtx.LICFlags[uiMergeCand];
#endif
  }
  //update
  for (uint32_t uiMergeCand = ((mrgCandIdx < 0) ? 0 : (mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE)*ADAPTIVE_IBC_SUB_GROUP_SIZE); uiMergeCand < (((mrgCandIdx < 0) || ((mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE > mrgCtx.numValidMergeCand)) ? mrgCtx.numValidMergeCand : ((mrgCandIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE)); ++uiMergeCand)
  {
    bool firstGroup = (uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE) == 0 ? true : false;
    bool lastGroup = ((uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE >= mrgCtx.numValidMergeCand) ? true : false;
    if (lastGroup && !firstGroup)
    {
      break;
    }
    mrgCtx.setMergeInfo(pu, uiMergeCand);
    if (pu.bv == Mv(0, 0))
    {
      break;
    }
    mrgCtx.BcwIdx[uiMergeCand] = mrgCtxTmp.BcwIdx[RdCandList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_IBC_SUB_GROUP_SIZE]];
    mrgCtx.interDirNeighbours[uiMergeCand] = mrgCtxTmp.interDirNeighbours[RdCandList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_IBC_SUB_GROUP_SIZE]];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_IBC_SUB_GROUP_SIZE] << 1)];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_IBC_SUB_GROUP_SIZE] << 1) + 1];
    mrgCtx.useAltHpelIf[uiMergeCand] = mrgCtxTmp.useAltHpelIf[RdCandList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_IBC_SUB_GROUP_SIZE]];
#if INTER_LIC
    mrgCtx.LICFlags[uiMergeCand] = mrgCtxTmp.LICFlags[RdCandList[uiMergeCand / ADAPTIVE_IBC_SUB_GROUP_SIZE][uiMergeCand%ADAPTIVE_IBC_SUB_GROUP_SIZE]];
#endif
  }
}
bool InterPrediction::xAMLIBCGetCurBlkTemplate(PredictionUnit& pu, int nCurBlkWidth, int nCurBlkHeight)
{
  m_bAMLTemplateAvailabe[0] = xAMLIsTopTempAvailable(pu);
  m_bAMLTemplateAvailabe[1] = xAMLIsLeftTempAvailable(pu);

  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    return false;
  }

  const Picture&  currPic = *pu.cs->picture;
  const CPelBuf recBuf = currPic.getRecoBuf(pu.cs->picture->blocks[COMPONENT_Y]);
  /* std::vector<Pel>& invLUT = m_pcReshape->getInvLUT();*/

  if (m_bAMLTemplateAvailabe[0])
  {
    const Pel*    rec = recBuf.bufAt(pu.blocks[COMPONENT_Y].pos().offset(0, -AML_MERGE_TEMPLATE_SIZE));
    PelBuf pcYBuf = PelBuf(m_acYuvCurAMLTemplate[0][0], nCurBlkWidth, AML_MERGE_TEMPLATE_SIZE);
    Pel*   pcY = pcYBuf.bufAt(0, 0);
    for (int k = 0; k < nCurBlkWidth; k++)
    {
      for (int l = 0; l < AML_MERGE_TEMPLATE_SIZE; l++)
      {
        int recVal = rec[k + l * recBuf.stride];

        //if (m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag())
        //{
        //  recVal = invLUT[recVal];
        //}

        pcY[k + l * nCurBlkWidth] = recVal;
      }
    }
  }

  if (m_bAMLTemplateAvailabe[1])
  {
    PelBuf pcYBuf = PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nCurBlkHeight);
    Pel*   pcY = pcYBuf.bufAt(0, 0);
    const Pel*    rec = recBuf.bufAt(pu.blocks[COMPONENT_Y].pos().offset(-AML_MERGE_TEMPLATE_SIZE, 0));
    for (int k = 0; k < nCurBlkHeight; k++)
    {
      for (int l = 0; l < AML_MERGE_TEMPLATE_SIZE; l++)
      {
        int recVal = rec[recBuf.stride * k + l];

        //if (m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag())
        //{
        //  recVal = invLUT[recVal];
        //}

        pcY[AML_MERGE_TEMPLATE_SIZE * k + l] = recVal;
      }
    }
  }

  return true;
}
void InterPrediction::getIBCAMLRefTemplate(PredictionUnit &pu, int nCurBlkWidth, int nCurBlkHeight)
{
  Mv mvCurr;
  mvCurr = pu.bv;
  const int lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
  const int horShift  = (lumaShift + ::getComponentScaleX(COMPONENT_Y, pu.chromaFormat));
  const int verShift  = (lumaShift + ::getComponentScaleY(COMPONENT_Y, pu.chromaFormat));
  const Picture&  currPic = *pu.cs->picture;
  const CPelBuf recBuf = currPic.getRecoBuf(pu.cs->picture->blocks[COMPONENT_Y]);
  /* std::vector<Pel>& invLUT = m_pcReshape->getInvLUT();*/
  if (m_bAMLTemplateAvailabe[0])
  {
    Mv mvTop(0, -AML_MERGE_TEMPLATE_SIZE);
    mvTop += mvCurr;

    MotionInfo miTop;
    miTop.mv[0] = Mv(mvTop.hor <<horShift , mvTop.ver<< verShift);
    miTop.refIdx[0] = MAX_NUM_REF;
    if (!PU::checkIsIBCCandidateValid(pu, miTop))
    {
      mvTop = mvCurr;
    }
    const Pel*    rec = recBuf.bufAt(pu.blocks[COMPONENT_Y].pos().offset(mvTop.hor, mvTop.ver));
    PelBuf pcYBuf = PelBuf(m_acYuvRefAMLTemplate[0][0], nCurBlkWidth, AML_MERGE_TEMPLATE_SIZE);
    Pel*   pcY = pcYBuf.bufAt(0, 0);
    for (int k = 0; k < nCurBlkWidth; k++)
    {
      for (int l = 0; l < AML_MERGE_TEMPLATE_SIZE; l++)
      {
        int recVal = rec[k + l * recBuf.stride];

        //if (m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag())
        //{
        //  recVal = invLUT[recVal];
        //}

        pcY[k + l * nCurBlkWidth] = recVal;
      }
    }
  }

  if (m_bAMLTemplateAvailabe[1])
  {
    Mv mvLeft(-AML_MERGE_TEMPLATE_SIZE, 0);
    mvLeft += mvCurr;

    MotionInfo miLeft;
    miLeft.mv[0] = Mv(mvLeft.hor <<horShift , mvLeft.ver<< verShift);
    miLeft.refIdx[0] = MAX_NUM_REF;
    if (!PU::checkIsIBCCandidateValid(pu, miLeft))
    {
      mvLeft = mvCurr;
    }
    PelBuf pcYBuf = PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nCurBlkHeight);
    Pel*   pcY = pcYBuf.bufAt(0, 0);
    const Pel*    rec = recBuf.bufAt(pu.blocks[COMPONENT_Y].pos().offset( mvLeft.hor,  mvLeft.ver));
    for (int k = 0; k < nCurBlkHeight; k++)
    {
      for (int l = 0; l < AML_MERGE_TEMPLATE_SIZE; l++)
      {
        int recVal = rec[recBuf.stride * k + l];

        //if (m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag())
        //{
        //  recVal = invLUT[recVal];
        //}

        pcY[AML_MERGE_TEMPLATE_SIZE * k + l] = recVal;
      }
    }
  }
}
#endif
#if JVET_Z0075_IBC_HMVP_ENLARGE
void  InterPrediction::adjustIBCMergeCandidates(PredictionUnit &pu, MergeCtx& mrgCtx,uint32_t startPos,uint32_t endPos)
{
#if JVET_Z0084_IBC_TM
  if (mrgCtx.numValidMergeCand <= 1)
  {
    return;
  }
#endif

  uint32_t RdCandList[IBC_MRG_MAX_NUM_CANDS_MEM];
  Distortion candCostList[IBC_MRG_MAX_NUM_CANDS_MEM];

  for (uint32_t i = 0; i < IBC_MRG_MAX_NUM_CANDS_MEM; i++)
  {
    RdCandList[i] = i;
    candCostList[i] = MAX_UINT;
  }

  Distortion uiCost;

  DistParam cDistParam;
  cDistParam.applyWeight = false;

  /*const SPS &sps = *pu.cs->sps;
  Position puPos = pu.lumaPos();*/
  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;

  if (!xAMLIBCGetCurBlkTemplate(pu, nWidth, nHeight))
  {
    return;
  }

  for (uint32_t uiMergeCand = startPos; uiMergeCand < endPos; ++uiMergeCand)
  {
    uiCost = 0;

    mrgCtx.setMergeInfo(pu, uiMergeCand);

    if (pu.bv == Mv(0, 0))
    {
      break;
    }

    PelBuf pcBufPredRefTop = PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE);
    PelBuf pcBufPredCurTop = PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE);
    PelBuf pcBufPredRefLeft = PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight);
    PelBuf pcBufPredCurLeft = PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight);

    getIBCAMLRefTemplate(pu, nWidth, nHeight);

    if (m_bAMLTemplateAvailabe[0])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop, pcBufPredRefTop, pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

      uiCost += cDistParam.distFunc(cDistParam);
    }

    if (m_bAMLTemplateAvailabe[1])
    {
      m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft, pcBufPredRefLeft, pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

      uiCost += cDistParam.distFunc(cDistParam);
    }

    updateCandList(uiMergeCand, uiCost, IBC_MRG_MAX_NUM_CANDS_MEM, RdCandList, candCostList);
  }

  updateIBCCandInfo(pu, mrgCtx, RdCandList, startPos, endPos);
}
void  InterPrediction::updateIBCCandInfo(PredictionUnit &pu, MergeCtx& mrgCtx, uint32_t* RdCandList,uint32_t startPos,uint32_t endPos)
{
  MergeCtx mrgCtxTmp;
  for (uint32_t ui = 0; ui < IBC_MRG_MAX_NUM_CANDS_MEM; ++ui)
  {
    mrgCtxTmp.BcwIdx[ui] = BCW_DEFAULT;
    mrgCtxTmp.interDirNeighbours[ui] = 0;
    mrgCtxTmp.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtxTmp.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtxTmp.useAltHpelIf[ui] = false;
#if INTER_LIC
    mrgCtxTmp.LICFlags[ui] = false;
#endif
  }
  for (uint32_t uiMergeCand = startPos; uiMergeCand < endPos; ++uiMergeCand)
  {
    mrgCtx.setMergeInfo(pu, uiMergeCand);
    if (pu.bv == Mv(0, 0))
    {
      break;
    }
    mrgCtxTmp.BcwIdx[uiMergeCand] = mrgCtx.BcwIdx[uiMergeCand];
    mrgCtxTmp.interDirNeighbours[uiMergeCand] = mrgCtx.interDirNeighbours[uiMergeCand];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)];
    mrgCtxTmp.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1];
    mrgCtxTmp.useAltHpelIf[uiMergeCand] = mrgCtx.useAltHpelIf[uiMergeCand];
#if INTER_LIC 
    mrgCtxTmp.LICFlags[uiMergeCand] = mrgCtx.LICFlags[uiMergeCand];
#endif
  }
  //update
  for (uint32_t uiMergeCand = startPos; uiMergeCand < endPos; ++uiMergeCand)
  {
    mrgCtx.setMergeInfo(pu, uiMergeCand);
    if (pu.bv == Mv(0, 0))
    {
      break;
    }
    mrgCtx.BcwIdx[uiMergeCand] = mrgCtxTmp.BcwIdx[RdCandList[uiMergeCand -startPos]];
    mrgCtx.interDirNeighbours[uiMergeCand] = mrgCtxTmp.interDirNeighbours[RdCandList[uiMergeCand -startPos]];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1)] = mrgCtxTmp.mvFieldNeighbours[RdCandList[uiMergeCand -startPos] << 1];
    mrgCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1] = mrgCtxTmp.mvFieldNeighbours[(RdCandList[uiMergeCand -startPos] << 1) + 1];
    mrgCtx.useAltHpelIf[uiMergeCand] = mrgCtxTmp.useAltHpelIf[RdCandList[uiMergeCand -startPos]];
#if INTER_LIC
    mrgCtx.LICFlags[uiMergeCand] = mrgCtxTmp.LICFlags[RdCandList[uiMergeCand -startPos]];
#endif
  }
}
#endif
#endif
#if JVET_Z0061_TM_OBMC
void InterPrediction::xOBMCWeightedAverageY(const PredictionUnit &pu, const CPelUnitBuf &pcYuvSrc0,
                                            const CPelUnitBuf &pcYuvSrc1, PelUnitBuf &pcYuvDst,
                                            const BitDepths &clipBitDepths, const ClpRngs &clpRngs, MotionInfo currMi)
{
  const int iRefIdx0 = currMi.refIdx[0];
  const int iRefIdx1 = currMi.refIdx[1];

  if (iRefIdx0 >= 0 && iRefIdx1 >= 0)
  {
    if (pu.cu->BcwIdx != BCW_DEFAULT)
    {
#if JVET_Z0136_OOB
      bool isOOB[2] = { false, false };
      pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, pu.cu->BcwIdx, false, true, pu.cs->mcMask, -1,
                              pu.cs->mcMaskChroma, -1, isOOB);
#else
      pcYuvDst.addWeightedAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, pu.cu->BcwIdx, false, true);
#endif
    }
    else
    {
      pcYuvDst.addAvg(pcYuvSrc0, pcYuvSrc1, clpRngs, false, true);
    }
  }
  else if (iRefIdx0 >= 0 && iRefIdx1 < 0)
  {
    pcYuvDst.copyClip(pcYuvSrc0, clpRngs, true);
  }
  else if (iRefIdx0 < 0 && iRefIdx1 >= 0)
  {
    pcYuvDst.copyClip(pcYuvSrc1, clpRngs, true);
  }
}
int InterPrediction::selectOBMCmode(PredictionUnit &pu, PredictionUnit &subblockPu, bool isAbove, int iLength,
                                    uint32_t uiMinCUW, Position off)
{
  const Position posSubBlock(pu.lumaPos().offset(off));
  Position       posNeighborMotion = Position(0, 0);
  if (isAbove)
  {
    posNeighborMotion = posSubBlock.offset(0, -1);
  }
  else if (!isAbove)
  {
    posNeighborMotion = posSubBlock.offset(-1, 0);
  }
  PredictionUnit *tmpPu = nullptr;
  tmpPu                 = pu.cs->getPU(posNeighborMotion, pu.chType);

  if (!tmpPu)
  {
    return 0;
  }

  MotionInfo neigMi = tmpPu->getMotionInfo(posNeighborMotion);
  MotionInfo currMi = pu.getMotionInfo(posSubBlock);
  Distortion candCostList[3];
  for (uint32_t i = 0; i < 3; i++)
  {
    candCostList[i] = MAX_UINT;
  }

  Distortion uiCost;
  DistParam  cDistParam;
  cDistParam.applyWeight = false;
  int nWidth             = isAbove ? (iLength * uiMinCUW) : uiMinCUW;
  int nHeight            = isAbove ? uiMinCUW : (iLength * uiMinCUW);

  if (!xAMLGetCurBlkTemplate(pu, pu.lumaSize().width, pu.lumaSize().height))
  {
    return 0;
  }

  // Process above boundary
  PelUnitBuf pcBufPredCurTopTmp =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], pu.lumaSize().width, TM_OBMC_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredCurTop =
    pcBufPredCurTopTmp.subBuf(UnitArea(pu.chromaFormat, Area(off.x, 0, nWidth, TM_OBMC_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredRefTop0 =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplateOBMC[0][0], nWidth, TM_OBMC_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredRefTop1 =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplateOBMC[1][0], nWidth, TM_OBMC_TEMPLATE_SIZE)));
  PelUnitBuf pcBufBlendDstAbove =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvBlendTemplateOBMC[0][0], nWidth, TM_OBMC_TEMPLATE_SIZE)));

  if ((isAbove) && (m_bAMLTemplateAvailabe[0]))
  {
    // 0: use current mv to do MC for template;
    getBlkOBMCRefTemplate(subblockPu, pcBufPredRefTop0, isAbove, currMi);
    uiCost = 0;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop0.Y(),
                             subblockPu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    uiCost          = cDistParam.distFunc(cDistParam);
    candCostList[0] = uiCost;

    // 1: use neighbour mv to do MC for template;
    getBlkOBMCRefTemplate(subblockPu, pcBufPredRefTop1, isAbove, neigMi);
    uiCost = 0;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop1.Y(),
                             subblockPu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    uiCost          = cDistParam.distFunc(cDistParam);
    candCostList[1] = uiCost;

    // 2: Now calculate the blending template
    uiCost               = 0;
    CPelUnitBuf srcPred0 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplateOBMC[0][0], pcBufPredCurTop.Y()));
    CPelUnitBuf srcPred1 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAboveTemplateOBMC[1][0], pcBufPredCurTop.Y()));

    for (int i = 0; i < nWidth; i++)
    {
      pcBufBlendDstAbove.Y().buf[i] = (26 * srcPred0.Y().buf[i] + 6 * srcPred1.Y().buf[i] + 16) >> 5;
    }

    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufBlendDstAbove.Y(),
                             subblockPu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    uiCost          = cDistParam.distFunc(cDistParam);
    candCostList[2] = uiCost;

    int bestOBMCmode = 0;
    if ((candCostList[0] < candCostList[1]) && (candCostList[0] < candCostList[2]))
    {
      bestOBMCmode = 1;
    }
    else
    {
      candCostList[0] = candCostList[0] << 3;
      candCostList[1] = candCostList[1] << 3;
      if ((candCostList[1] + (candCostList[1] >> 2) + (candCostList[1] >> 3)) <= candCostList[0])
      {
        bestOBMCmode = 2;
      }
      else if (candCostList[0] <= candCostList[1])
      {
        bestOBMCmode = 3;
      }
      else if (candCostList[1] <= candCostList[0])
      {
        bestOBMCmode = 4;
      }
    }
    return bestOBMCmode;
  }

  // Process left boundary
  PelUnitBuf pcBufPredCurLeftTmp =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], TM_OBMC_TEMPLATE_SIZE, pu.lumaSize().height)));
  PelUnitBuf pcBufPredCurLeft =
    pcBufPredCurLeftTmp.subBuf(UnitArea(pu.chromaFormat, Area(0, off.y, TM_OBMC_TEMPLATE_SIZE, nHeight)));
  PelUnitBuf pcBufPredRefLeft0 =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplateOBMC[0][0], TM_OBMC_TEMPLATE_SIZE, nHeight)));
  PelUnitBuf pcBufPredRefLeft1 =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplateOBMC[1][0], TM_OBMC_TEMPLATE_SIZE, nHeight)));
  PelUnitBuf pcBufBlendDstLeft =
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvBlendTemplateOBMC[1][0], TM_OBMC_TEMPLATE_SIZE, nHeight)));

  if ((!isAbove) && (m_bAMLTemplateAvailabe[1]))
  {
    // 0: use current mv to do MC for template;
    getBlkOBMCRefTemplate(subblockPu, pcBufPredRefLeft0, isAbove, currMi);
    uiCost = 0;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft0.Y(),
                             subblockPu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    uiCost          = cDistParam.distFunc(cDistParam);
    candCostList[0] = uiCost;

    // 1: use neighbour mv to do MC for template;
    getBlkOBMCRefTemplate(subblockPu, pcBufPredRefLeft1, isAbove, neigMi);
    uiCost = 0;
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft1.Y(),
                             subblockPu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    uiCost          = cDistParam.distFunc(cDistParam);
    candCostList[1] = uiCost;

    // 2: Now calculate the blending template
    uiCost               = 0;
    CPelUnitBuf srcPred0 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplateOBMC[0][0], pcBufPredCurLeft.Y()));
    CPelUnitBuf srcPred1 = CPelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefLeftTemplateOBMC[1][0], pcBufPredCurLeft.Y()));
    int         idx      = 0;
    for (int i = 0; i < nHeight; i++)
    {
      pcBufBlendDstLeft.Y().buf[idx] = (26 * srcPred0.Y().buf[idx] + 6 * srcPred1.Y().buf[idx] + 16) >> 5;
      idx += pcBufBlendDstLeft.bufs[COMPONENT_Y].stride;
    }
    m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufBlendDstLeft.Y(),
                             subblockPu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);
    uiCost          = cDistParam.distFunc(cDistParam);
    candCostList[2] = uiCost;

    int bestOBMCmode = 0;
    if ((candCostList[0] < candCostList[1]) && (candCostList[0] < candCostList[2]))
    {
      bestOBMCmode = 1;
    }
    else
    {
      candCostList[0] = candCostList[0] << 3;
      candCostList[1] = candCostList[1] << 3;
      if ((candCostList[1] + (candCostList[1] >> 2) + (candCostList[1] >> 3)) <= candCostList[0])
      {
        bestOBMCmode = 2;
      }
      else if (candCostList[0] <= candCostList[1])
      {
        bestOBMCmode = 3;
      }
      else if (candCostList[1] <= candCostList[0])
      {
        bestOBMCmode = 4;
      }
    }
    return bestOBMCmode;
  }
  else
  {
    return 0;
  }
}
bool InterPrediction::xCheckIdenticalMotionOBMC(PredictionUnit &pu, MotionInfo tryMi)
{
  const Slice &slice = *pu.cs->slice;

  if (slice.isInterB() && !pu.cs->pps->getWPBiPred())
  {
    if (tryMi.refIdx[0] >= 0 && tryMi.refIdx[1] >= 0)
    {
      int RefPOCL0 = slice.getRefPic(REF_PIC_LIST_0, tryMi.refIdx[0])->getPOC();
      int RefPOCL1 = slice.getRefPic(REF_PIC_LIST_1, tryMi.refIdx[1])->getPOC();

      if (RefPOCL0 == RefPOCL1)
      {
        if (!pu.cu->affine)
        {
          if (tryMi.mv[0] == tryMi.mv[1])
          {
            return true;
          }
        }
        else
        {
          if ((pu.cu->affineType == AFFINEMODEL_4PARAM && (pu.mvAffi[0][0] == pu.mvAffi[1][0])
               && (pu.mvAffi[0][1] == pu.mvAffi[1][1]))
              || (pu.cu->affineType == AFFINEMODEL_6PARAM && (pu.mvAffi[0][0] == pu.mvAffi[1][0])
                  && (pu.mvAffi[0][1] == pu.mvAffi[1][1]) && (pu.mvAffi[0][2] == pu.mvAffi[1][2])))
          {
            return true;
          }
        }
      }
    }
  }
  return false;
}
void InterPrediction::getBlkOBMCRefTemplate(PredictionUnit &subblockPu, PelUnitBuf &pcBufPredRef, bool isAbove,
                                            MotionInfo tryMi)
{
  Mv        mvCurr;
  const int lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
  const int horShift  = (lumaShift + ::getComponentScaleX(COMPONENT_Y, subblockPu.chromaFormat));
  const int verShift  = (lumaShift + ::getComponentScaleY(COMPONENT_Y, subblockPu.chromaFormat));

  if (xCheckIdenticalMotionOBMC(subblockPu, tryMi))
  {
    mvCurr      = tryMi.mv[0];
    Mv subPelMv = mvCurr;
    clipMv(mvCurr, subblockPu.lumaPos(), subblockPu.lumaSize(), *subblockPu.cs->sps, *subblockPu.cs->pps);
    CHECK(tryMi.refIdx[0] < 0, "invalid ref idx");

    if ((isAbove) && (m_bAMLTemplateAvailabe[0]))
    {
      Mv mvTop(0, -(TM_OBMC_TEMPLATE_SIZE << verShift));
      mvTop += subPelMv;
      clipMv(mvTop, subblockPu.lumaPos(), subblockPu.lumaSize(), *subblockPu.cs->sps, *subblockPu.cs->pps);

#if RPR_ENABLE
      const Picture *            picRef = subblockPu.cu->slice->getRefPic(REF_PIC_LIST_0, tryMi.refIdx[0])->unscaledPic;
      const std::pair<int, int> &scalingRatio = subblockPu.cu->slice->getScalingRatio(REF_PIC_LIST_0, tryMi.refIdx[0]);
      xPredInterBlk(COMPONENT_Y, subblockPu, picRef, mvTop,
                    pcBufPredRef, false, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, scalingRatio, 0, 0,
                    false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                    , false
#if INTER_LIC
                    , false, Mv(0, 0)
#endif
#endif
                    , true);
#else
      xPredInterBlk(COMPONENT_Y, subblockPu, subblockPu.cu->slice->getRefPic(REF_PIC_LIST_0, tryMi.refIdx[0]), mvTop,
                    pcBufPredRef, false, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X, 0, 0, false,
                    NULL, 0 
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING       
        , false
#if INTER_LIC        
        , false, Mv(0, 0)
#endif
#endif
                    , true);
#endif
    }

    if ((!isAbove) && (m_bAMLTemplateAvailabe[1]))
    {
      Mv mvLeft(-(TM_OBMC_TEMPLATE_SIZE << horShift), 0);
      mvLeft += subPelMv;
      clipMv(mvLeft, subblockPu.lumaPos(), subblockPu.lumaSize(), *subblockPu.cs->sps, *subblockPu.cs->pps);
#if RPR_ENABLE
      const Picture *            picRef = subblockPu.cu->slice->getRefPic(REF_PIC_LIST_0, tryMi.refIdx[0])->unscaledPic;
      const std::pair<int, int> &scalingRatio = subblockPu.cu->slice->getScalingRatio(REF_PIC_LIST_0, tryMi.refIdx[0]);
      xPredInterBlk(COMPONENT_Y, subblockPu, picRef, mvLeft,
                    pcBufPredRef, false, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, scalingRatio, 0, 0,
                    false,
                    NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                    , false
#if INTER_LIC
                    , false, Mv(0, 0)
#endif
#endif
                    , true);
#else
      xPredInterBlk(COMPONENT_Y, subblockPu, subblockPu.cu->slice->getRefPic(REF_PIC_LIST_0, tryMi.refIdx[0]), mvLeft,
                    pcBufPredRef, false, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X, 0, 0, false,
                    NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING         
        , false
#if INTER_LIC        
        , false, Mv(0, 0)
#endif
#endif
                    , true);
#endif
    }
  }
  else
  {
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      if (tryMi.refIdx[refList] < 0)
      {
        continue;
      }
      RefPicList eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
      CHECK(tryMi.refIdx[refList] >= subblockPu.cu->slice->getNumRefIdx(eRefPicList), "Invalid reference index");

      m_iRefListIdx = refList;
      mvCurr        = tryMi.mv[refList];
      Mv subPelMv   = mvCurr;
      clipMv(mvCurr, subblockPu.lumaPos(), subblockPu.lumaSize(), *subblockPu.cs->sps, *subblockPu.cs->pps);

      if ((isAbove) && (m_bAMLTemplateAvailabe[0]))
      {
        Mv mvTop(0, -(TM_OBMC_TEMPLATE_SIZE << verShift));
        mvTop += subPelMv;

        clipMv(mvTop, subblockPu.lumaPos(), subblockPu.lumaSize(), *subblockPu.cs->sps, *subblockPu.cs->pps);

        PelUnitBuf pcMbBuf =
          PelUnitBuf(subblockPu.chromaFormat,
                     PelBuf(m_acYuvRefAboveTemplate[refList][0], pcBufPredRef.Y()));
        if (tryMi.refIdx[0] >= 0 && tryMi.refIdx[1] >= 0)
        {
#if RPR_ENABLE
          const Picture *picRef = subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList])->unscaledPic;
          const std::pair<int, int> &scalingRatio =
            subblockPu.cu->slice->getScalingRatio(eRefPicList, tryMi.refIdx[refList]);
          xPredInterBlk(COMPONENT_Y, subblockPu, picRef, mvTop, pcMbBuf, true,
                        subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, scalingRatio, 0, 0,
                        false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                        , false
#if INTER_LIC
                        , false, Mv(0, 0)
#endif
#endif
                        , true);
#else
          xPredInterBlk(COMPONENT_Y, subblockPu, subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList]),
                        mvTop, pcMbBuf, true, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X, 0, 0,
                        false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING            
            , false
#if INTER_LIC             
            , false, Mv(0, 0)
#endif
#endif
                        , true);
#endif
          if (refList == 1)
          {
            CPelUnitBuf srcPred0 =
              CPelUnitBuf(subblockPu.chromaFormat,
                          PelBuf(m_acYuvRefAboveTemplate[0][0], pcBufPredRef.Y()));
            CPelUnitBuf srcPred1 =
              CPelUnitBuf(subblockPu.chromaFormat,
                          PelBuf(m_acYuvRefAboveTemplate[1][0], pcBufPredRef.Y()));
            xOBMCWeightedAverageY(subblockPu, srcPred0, srcPred1, pcBufPredRef,
                                  subblockPu.cu->slice->getSPS()->getBitDepths(), subblockPu.cu->slice->clpRngs(),
                                  tryMi);
          }
        }
        else
        {
#if RPR_ENABLE
          const Picture *picRef = subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList])->unscaledPic;
          const std::pair<int, int> &scalingRatio =
            subblockPu.cu->slice->getScalingRatio(eRefPicList, tryMi.refIdx[refList]);
          xPredInterBlk(COMPONENT_Y, subblockPu, picRef, mvTop, pcBufPredRef, false,
                        subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, scalingRatio,
                        0, 0, false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                        , false
#if INTER_LIC
                        , false, Mv(0, 0)
#endif
#endif
                        , true);
#else
          xPredInterBlk(COMPONENT_Y, subblockPu, subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList]),
                        mvTop, pcBufPredRef, false, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X,
                        0, 0, false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING            
            , false
#if INTER_LIC             
            , false, Mv(0, 0)
#endif
#endif
                        , true);
#endif
        }
      }
      if ((!isAbove) && (m_bAMLTemplateAvailabe[1]))
      {
        Mv mvLeft(-(TM_OBMC_TEMPLATE_SIZE << horShift), 0);
        mvLeft += subPelMv;

        clipMv(mvLeft, subblockPu.lumaPos(), subblockPu.lumaSize(), *subblockPu.cs->sps, *subblockPu.cs->pps);

        PelUnitBuf pcMbBuf =
          PelUnitBuf(subblockPu.chromaFormat,
                     PelBuf(m_acYuvRefLeftTemplate[refList][0], pcBufPredRef.Y()));
        if (tryMi.refIdx[0] >= 0 && tryMi.refIdx[1] >= 0)
        {
#if RPR_ENABLE
          const Picture *picRef = subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList])->unscaledPic;
          const std::pair<int, int> &scalingRatio =
            subblockPu.cu->slice->getScalingRatio(eRefPicList, tryMi.refIdx[refList]);
          xPredInterBlk(COMPONENT_Y, subblockPu, picRef, mvLeft, pcMbBuf, true,
                        subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, scalingRatio, 0, 0,
                        false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                        , false
#if INTER_LIC
                        , false, Mv(0, 0)
#endif
#endif
                        , true);
#else
          xPredInterBlk(COMPONENT_Y, subblockPu, subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList]),
                        mvLeft, pcMbBuf, true, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X, 0, 0,
                        false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING              
            , false
#if INTER_LIC             
            , false, Mv(0, 0)
#endif
#endif
                        , true);
#endif
          if (refList == 1)
          {
            CPelUnitBuf srcPred0 =
              CPelUnitBuf(subblockPu.chromaFormat,
                          PelBuf(m_acYuvRefLeftTemplate[0][0], pcBufPredRef.Y()));
            CPelUnitBuf srcPred1 =
              CPelUnitBuf(subblockPu.chromaFormat,
                          PelBuf(m_acYuvRefLeftTemplate[1][0], pcBufPredRef.Y()));
            xOBMCWeightedAverageY(subblockPu, srcPred0, srcPred1, pcBufPredRef,
                                  subblockPu.cu->slice->getSPS()->getBitDepths(), subblockPu.cu->slice->clpRngs(),
                                  tryMi);
          }
        }
        else
        {
#if RPR_ENABLE
          const Picture *picRef = subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList])->unscaledPic;
          const std::pair<int, int> &scalingRatio =
            subblockPu.cu->slice->getScalingRatio(eRefPicList, tryMi.refIdx[refList]);
          xPredInterBlk(COMPONENT_Y, subblockPu, picRef, mvLeft, pcBufPredRef, false,
                        subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, scalingRatio,
                        0, 0, false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                        , false
#if INTER_LIC
                        , false, Mv(0, 0)
#endif
#endif
                        , true);
#else
          xPredInterBlk(COMPONENT_Y, subblockPu, subblockPu.cu->slice->getRefPic(eRefPicList, tryMi.refIdx[refList]),
                        mvLeft, pcBufPredRef, false, subblockPu.cu->slice->clpRng(COMPONENT_Y), false, false, SCALE_1X,
                        0, 0, false, NULL, 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING            
            , false
#if INTER_LIC             
            , false, Mv(0, 0)
#endif
#endif
                        , true);
#endif
        }
      }
    }
  }
}

void InterPrediction::xSubblockOBMCCopy(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst,
                                        PelUnitBuf &pcYuvPredSrc, int iDir)
{
  int iWidth  = pu.blocks[eComp].width;
  int iHeight = pu.blocks[eComp].height;

  if (iWidth == 0 || iHeight == 0)
  {
    return;
  }

  Pel *     pOrgDst   = pcYuvPredDst.bufs[eComp].buf;
  Pel *     pOrgSrc   = pcYuvPredSrc.bufs[eComp].buf;
  const int strideDst = pcYuvPredDst.bufs[eComp].stride;
  const int strideSrc = pcYuvPredSrc.bufs[eComp].stride;

  if (iDir == 0)   // above
  {
    for (int i = 0; i < iWidth; i++)
    {
      Pel *pDst = pOrgDst;
      Pel *pSrc = pOrgSrc;
      pDst[i]   = pSrc[i];
    }
  }

  if (iDir == 1)   // left
  {
    Pel *pDst = pOrgDst;
    Pel *pSrc = pOrgSrc;
    for (int i = 0; i < iHeight; i++)
    {
      pDst[0] = pSrc[0];
      pDst += strideDst;
      pSrc += strideSrc;
    }
  }
}

void InterPrediction::xSubblockTMOBMC(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst,
                                      PelUnitBuf &pcYuvPredSrc, int iDir, int iOBMCmode)
{
  int iWidth  = pu.blocks[eComp].width;
  int iHeight = pu.blocks[eComp].height;

  if (iWidth == 0 || iHeight == 0)
  {
    return;
  }

  Pel *     pOrgDst   = pcYuvPredDst.bufs[eComp].buf;
  Pel *     pOrgSrc   = pcYuvPredSrc.bufs[eComp].buf;
  const int strideDst = pcYuvPredDst.bufs[eComp].stride;
  const int strideSrc = pcYuvPredSrc.bufs[eComp].stride;

  if (iDir == 0)   // above
  {
    for (int i = 0; i < iWidth; i++)
    {
      Pel *pDst = pOrgDst;
      Pel *pSrc = pOrgSrc;
      if (iOBMCmode == 2)   // neighbor is best
      {
        pDst[i] = (26 * pDst[i] + 6 * pSrc[i] + 16) >> 5;

        if (eComp == COMPONENT_Y)
        {
          pDst += strideDst;
          pSrc += strideSrc;
          pDst[i] = (7 * pDst[i] + pSrc[i] + 4) >> 3;

          pDst += strideDst;
          pSrc += strideSrc;
          pDst[i] = (15 * pDst[i] + pSrc[i] + 8) >> 4;

          pDst += strideDst;
          pSrc += strideSrc;
          pDst[i] = (31 * pDst[i] + pSrc[i] + 16) >> 5;
        }
      }
      else if (iOBMCmode == 4)
      {
        pDst[i] = (7 * pDst[i] + pSrc[i] + 4) >> 3;

        if (eComp == COMPONENT_Y)
        {
          pDst += strideDst;
          pSrc += strideSrc;
          pDst[i] = (15 * pDst[i] + pSrc[i] + 8) >> 4;

          pDst += strideDst;
          pSrc += strideSrc;
          pDst[i] = (31 * pDst[i] + pSrc[i] + 16) >> 5;
        }
      }
      else   // blending is best
      {
        pDst[i] = (15 * pDst[i] + pSrc[i] + 8) >> 4;

        // luma blend 3 lines
        if (eComp == COMPONENT_Y)
        {
          pDst += strideDst;
          pSrc += strideSrc;
          pDst[i] = (31 * pDst[i] + pSrc[i] + 16) >> 5;
        }
      }
    }
  }

  if (iDir == 1)   // left
  {
    Pel *pDst = pOrgDst;
    Pel *pSrc = pOrgSrc;
    for (int i = 0; i < iHeight; i++)
    {
      if (iOBMCmode == 2)   // neighbor is best
      {
        pDst[0] = (26 * pDst[0] + 6 * pSrc[0] + 16) >> 5;

        if (eComp == COMPONENT_Y)
        {
          pDst[1] = (7 * pDst[1] + pSrc[1] + 4) >> 3;
          pDst[2] = (15 * pDst[2] + pSrc[2] + 8) >> 4;
          pDst[3] = (31 * pDst[3] + pSrc[3] + 16) >> 5;
        }
      }
      else if (iOBMCmode == 4)   // neighbor is best
      {
        pDst[0] = (7 * pDst[0] + pSrc[0] + 4) >> 3;

        if (eComp == COMPONENT_Y)
        {
          pDst[1] = (15 * pDst[1] + pSrc[1] + 8) >> 4;
          pDst[2] = (31 * pDst[2] + pSrc[2] + 16) >> 5;
        }
      }
      else   // blending is best
      {
        pDst[0] = (15 * pDst[0] + pSrc[0] + 8) >> 4;

        // luma blend 3 lines
        if (eComp == COMPONENT_Y)
        {
          pDst[1] = (31 * pDst[1] + pSrc[1] + 16) >> 5;
        }
      }
      pDst += strideDst;
      pSrc += strideSrc;
    }
  }
}
#endif
void InterPrediction::xFillIBCBuffer(CodingUnit &cu)
{
#if JVET_Z0118_GDR
  bool isCleanCu              = cu.cs->isClean(cu);
  bool useCleanIBCBuffer      = cu.cs->isInGdrIntervalOrRecoveryPoc() && isCleanCu;
#endif

  for (auto &currPU : CU::traverseTUs(cu))
  {
    for (const CompArea &area : currPU.blocks)
    {
      if (!area.valid())
      {
        continue;
      }

#if JVET_Z0153_IBC_EXT_REF
      const int shiftSampleHor = ::getComponentScaleX(area.compID, cu.chromaFormat);
      const int shiftSampleVer = ::getComponentScaleY(area.compID, cu.chromaFormat);
      const int pux = area.x % (m_IBCBufferWidth  >> shiftSampleHor);
      const int puy = area.y % (m_IBCBufferHeight >> shiftSampleVer);
#else
      const unsigned int lcuWidth = cu.cs->slice->getSPS()->getMaxCUWidth();
      const int shiftSampleHor = ::getComponentScaleX(area.compID, cu.chromaFormat);
      const int shiftSampleVer = ::getComponentScaleY(area.compID, cu.chromaFormat);
      const int ctuSizeLog2Ver = floorLog2(lcuWidth) - shiftSampleVer;
      const int pux = area.x & ((m_IBCBufferWidth >> shiftSampleHor) - 1);
      const int puy = area.y & (( 1 << ctuSizeLog2Ver ) - 1);
#endif

      const CompArea dstArea = CompArea(area.compID, cu.chromaFormat, Position(pux, puy), Size(area.width, area.height));
      CPelBuf srcBuf = cu.cs->getRecoBuf(area);

#if JVET_Z0118_GDR
      PelBuf dstBuf;

      // 1. copy to Dirty IBC Buffer
      dstBuf = m_IBCBuffer0.getBuf(dstArea);
      dstBuf.copyFrom(srcBuf);      

      // 2. copy to Clean IBC Buffer 
      if (useCleanIBCBuffer)
      {
        dstBuf = m_IBCBuffer1.getBuf(dstArea);
        dstBuf.copyFrom(srcBuf);
      }      
#else
      PelBuf dstBuf = m_IBCBuffer.getBuf(dstArea);

      dstBuf.copyFrom(srcBuf);
#endif      
    }
  }
}

void InterPrediction::xIntraBlockCopy(PredictionUnit &pu, PelUnitBuf &predBuf, const ComponentID compID)
{
#if JVET_Z0118_GDR  
  bool isCleanCu           = pu.cs->isClean(pu);
  bool useCleanIBCBuffer   = pu.cs->isInGdrIntervalOrRecoveryPoc() && isCleanCu;
#endif

#if JVET_Z0153_IBC_EXT_REF
  const int shiftSampleHor = ::getComponentScaleX(compID, pu.chromaFormat);
  const int shiftSampleVer = ::getComponentScaleY(compID, pu.chromaFormat);
#else
  const unsigned int lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
  const int shiftSampleHor = ::getComponentScaleX(compID, pu.chromaFormat);
  const int shiftSampleVer = ::getComponentScaleY(compID, pu.chromaFormat);
  const int ctuSizeLog2Ver = floorLog2(lcuWidth) - shiftSampleVer;
#endif

  pu.bv = pu.mv[REF_PIC_LIST_0];
  pu.bv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
  int refx, refy;

  if (compID == COMPONENT_Y)
  {
    refx = pu.Y().x + pu.bv.hor;
    refy = pu.Y().y + pu.bv.ver;
  }
  else // compID == COMPONENT_Cb, COMPONENT_Cr
  { 
    refx = pu.Cb().x + (pu.bv.hor >> shiftSampleHor);
    refy = pu.Cb().y + (pu.bv.ver >> shiftSampleVer);
  }

#if JVET_Z0153_IBC_EXT_REF
  refx = refx % (m_IBCBufferWidth  >> shiftSampleHor);
  refy = refy % (m_IBCBufferHeight >> shiftSampleVer);
#else
  refx &= ((m_IBCBufferWidth >> shiftSampleHor) - 1);
  refy &= ((1 << ctuSizeLog2Ver) - 1);
#endif

#if JVET_Z0153_IBC_EXT_REF
  if (refy + predBuf.bufs[compID].height <= (m_IBCBufferHeight >> shiftSampleVer))
#else
  if (refx + predBuf.bufs[compID].width <= (m_IBCBufferWidth >> shiftSampleHor))
#endif
  {
    const CompArea srcArea = CompArea(compID, pu.chromaFormat, Position(refx, refy), Size(predBuf.bufs[compID].width, predBuf.bufs[compID].height));

#if JVET_Z0118_GDR
    if (useCleanIBCBuffer)
    {
      const CPelBuf refBuf = m_IBCBuffer1.getBuf(srcArea);
      predBuf.bufs[compID].copyFrom(refBuf);
    }
    else
    {
      const CPelBuf refBuf = m_IBCBuffer0.getBuf(srcArea);
      predBuf.bufs[compID].copyFrom(refBuf);
    }
#else
    const CPelBuf refBuf = m_IBCBuffer.getBuf(srcArea);
    predBuf.bufs[compID].copyFrom(refBuf);
#endif
  }
  else
  { //wrap around
#if JVET_Z0153_IBC_EXT_REF
#if JVET_Z0118_GDR
    if (useCleanIBCBuffer)
    {
      int height = (m_IBCBufferHeight >> shiftSampleVer) - refy;
      CompArea srcArea = CompArea(compID, pu.chromaFormat, Position(refx, refy), Size(predBuf.bufs[compID].width, height));
      CPelBuf srcBuf = m_IBCBuffer1.getBuf(srcArea);
      PelBuf  dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, 0)), predBuf.bufs[compID].stride, Size(predBuf.bufs[compID].width, height));
      dstBuf.copyFrom(srcBuf);

      height = refy + predBuf.bufs[compID].height - (m_IBCBufferHeight >> shiftSampleVer);
      srcArea = CompArea(compID, pu.chromaFormat, Position(refx, 0), Size(predBuf.bufs[compID].width, height));
      srcBuf = m_IBCBuffer1.getBuf(srcArea);
      dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, (m_IBCBufferHeight >> shiftSampleVer) - refy)), predBuf.bufs[compID].stride, Size(predBuf.bufs[compID].width, height));
      dstBuf.copyFrom(srcBuf);     
    }
    else
    {
      int height = (m_IBCBufferHeight >> shiftSampleVer) - refy;
      CompArea srcArea = CompArea(compID, pu.chromaFormat, Position(refx, refy), Size(predBuf.bufs[compID].width, height));
      CPelBuf srcBuf = m_IBCBuffer0.getBuf(srcArea);
      PelBuf dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, 0)), predBuf.bufs[compID].stride, Size(predBuf.bufs[compID].width, height));
      dstBuf.copyFrom(srcBuf);

      height = refy + predBuf.bufs[compID].height - (m_IBCBufferHeight >> shiftSampleVer);
      srcArea = CompArea(compID, pu.chromaFormat, Position(refx, 0), Size(predBuf.bufs[compID].width, height));
      srcBuf = m_IBCBuffer0.getBuf(srcArea);
      dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, (m_IBCBufferHeight >> shiftSampleVer) - refy)), predBuf.bufs[compID].stride, Size(predBuf.bufs[compID].width, height));
      dstBuf.copyFrom(srcBuf);
    }
#else
    int height = (m_IBCBufferHeight >> shiftSampleVer) - refy;
    CompArea srcArea = CompArea(compID, pu.chromaFormat, Position(refx, refy), Size(predBuf.bufs[compID].width, height));
    CPelBuf srcBuf = m_IBCBuffer.getBuf(srcArea);
    PelBuf dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, 0)), predBuf.bufs[compID].stride, Size(predBuf.bufs[compID].width, height));
    dstBuf.copyFrom(srcBuf);

    height = refy + predBuf.bufs[compID].height - (m_IBCBufferHeight >> shiftSampleVer);
    srcArea = CompArea(compID, pu.chromaFormat, Position(refx, 0), Size(predBuf.bufs[compID].width, height));
    srcBuf = m_IBCBuffer.getBuf(srcArea);
    dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, (m_IBCBufferHeight >> shiftSampleVer) - refy)), predBuf.bufs[compID].stride, Size(predBuf.bufs[compID].width, height));
    dstBuf.copyFrom(srcBuf);
#endif
#else
    int width = (m_IBCBufferWidth >> shiftSampleHor) - refx;
    CompArea srcArea = CompArea(compID, pu.chromaFormat, Position(refx, refy), Size(width, predBuf.bufs[compID].height));

#if JVET_Z0118_GDR
    if (useCleanIBCBuffer)
    {
      CPelBuf srcBuf = m_IBCBuffer1.getBuf(srcArea);
      PelBuf dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, 0)), predBuf.bufs[compID].stride, Size(width, predBuf.bufs[compID].height));
      dstBuf.copyFrom(srcBuf);

      width = refx + predBuf.bufs[compID].width - (m_IBCBufferWidth >> shiftSampleHor);
      srcArea = CompArea(compID, pu.chromaFormat, Position(0, refy), Size(width, predBuf.bufs[compID].height));
      srcBuf = m_IBCBuffer1.getBuf(srcArea);
      dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position((m_IBCBufferWidth >> shiftSampleHor) - refx, 0)), predBuf.bufs[compID].stride, Size(width, predBuf.bufs[compID].height));
      dstBuf.copyFrom(srcBuf);
    }
    else
    {
      CPelBuf srcBuf = m_IBCBuffer0.getBuf(srcArea);
      PelBuf dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, 0)), predBuf.bufs[compID].stride, Size(width, predBuf.bufs[compID].height));
      dstBuf.copyFrom(srcBuf);

      width = refx + predBuf.bufs[compID].width - (m_IBCBufferWidth >> shiftSampleHor);
      srcArea = CompArea(compID, pu.chromaFormat, Position(0, refy), Size(width, predBuf.bufs[compID].height));
      srcBuf = m_IBCBuffer0.getBuf(srcArea);
      dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position((m_IBCBufferWidth >> shiftSampleHor) - refx, 0)), predBuf.bufs[compID].stride, Size(width, predBuf.bufs[compID].height));
      dstBuf.copyFrom(srcBuf);
    }
#else
    CPelBuf srcBuf = m_IBCBuffer.getBuf(srcArea);
    PelBuf dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position(0, 0)), predBuf.bufs[compID].stride, Size(width, predBuf.bufs[compID].height));
    dstBuf.copyFrom(srcBuf);

    width = refx + predBuf.bufs[compID].width - (m_IBCBufferWidth >> shiftSampleHor);
    srcArea = CompArea(compID, pu.chromaFormat, Position(0, refy), Size(width, predBuf.bufs[compID].height));
    srcBuf = m_IBCBuffer.getBuf(srcArea);
    dstBuf = PelBuf(predBuf.bufs[compID].bufAt(Position((m_IBCBufferWidth >> shiftSampleHor) - refx, 0)), predBuf.bufs[compID].stride, Size(width, predBuf.bufs[compID].height));
    dstBuf.copyFrom(srcBuf);
#endif    
#endif
  }
}

void InterPrediction::resetIBCBuffer(const ChromaFormat chromaFormatIDC, const int ctuSize)
{
#if JVET_Z0153_IBC_EXT_REF
  const UnitArea area = UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, m_IBCBufferHeight));
#else
  const UnitArea area = UnitArea(chromaFormatIDC, Area(0, 0, m_IBCBufferWidth, ctuSize));
#endif
#if JVET_Z0118_GDR  
  m_IBCBuffer0.getBuf(area).fill(-1);  
#else
  m_IBCBuffer.getBuf(area).fill(-1);
#endif
}

#if JVET_Z0118_GDR
void InterPrediction::resetCurIBCBuffer(const ChromaFormat chromaFormatIDC, const Area ctuArea, const int ctuSize, const Pel dirtyPel)
{
#if JVET_Z0153_IBC_EXT_REF
  const int shiftSampleHor = ::getComponentScaleX(COMPONENT_Y, chromaFormatIDC);
  const int shiftSampleVer = ::getComponentScaleY(COMPONENT_Y, chromaFormatIDC);
  const int pux = ctuArea.x % (m_IBCBufferWidth  >> shiftSampleHor);
  const int puy = ctuArea.y % (m_IBCBufferHeight >> shiftSampleVer);
#else
  const int shiftSampleHor = ::getComponentScaleX(COMPONENT_Y, chromaFormatIDC);
  const int shiftSampleVer = ::getComponentScaleY(COMPONENT_Y, chromaFormatIDC);
  const int ctuSizeLog2Ver = floorLog2(ctuSize) - shiftSampleVer;
  const int pux = ctuArea.x & ((m_IBCBufferWidth >> shiftSampleHor) - 1);
  const int puy = ctuArea.y & ((1 << ctuSizeLog2Ver) - 1);
#endif

  const UnitArea area = UnitArea(chromaFormatIDC, Area(pux, puy, ctuSize, ctuSize));

  m_IBCBuffer1.getBuf(area).fill(dirtyPel);
}
#endif

void InterPrediction::resetVPDUforIBC(const ChromaFormat chromaFormatIDC, const int ctuSize, const int vSize, const int xPos, const int yPos)
{
#if JVET_Z0153_IBC_EXT_REF
  if(xPos == 0)
  {
    const UnitArea area = UnitArea(chromaFormatIDC, Area(0, yPos % m_IBCBufferHeight, m_IBCBufferWidth, ctuSize));
#if JVET_Z0118_GDR
    m_IBCBuffer0.getBuf(area).fill(-1);
#else
    m_IBCBuffer.getBuf(area).fill(-1);
#endif
  }

  if(xPos - 3 * ctuSize >= 0)
  {
    const UnitArea area = UnitArea(chromaFormatIDC, Area((xPos - 3 * ctuSize) % m_IBCBufferWidth, (yPos + ctuSize) % m_IBCBufferHeight, ctuSize, ctuSize));
#if JVET_Z0118_GDR
    m_IBCBuffer0.getBuf(area).fill(-1);
#else
    m_IBCBuffer.getBuf(area).fill(-1);
#endif
  }
#else
  const UnitArea area = UnitArea(chromaFormatIDC, Area(xPos & (m_IBCBufferWidth - 1), yPos & (ctuSize - 1), vSize, vSize));

#if JVET_Z0118_GDR
  m_IBCBuffer0.getBuf(area).fill(-1);
#else
  m_IBCBuffer.getBuf(area).fill(-1);
#endif
#endif
}

bool InterPrediction::isLumaBvValid(const int ctuSize, const int xCb, const int yCb, const int width, const int height, const int xBv, const int yBv)
{
#if JVET_Z0153_IBC_EXT_REF
  int refTLx = xCb + xBv;
  int refTLy = yCb + yBv;
#else
  if(((yCb + yBv) & (ctuSize - 1)) + height > ctuSize)
  {
    return false;
  }
  int refTLx = xCb + xBv;
  int refTLy = (yCb + yBv) & (ctuSize - 1);
#endif

#if JVET_Z0118_GDR
  PelBuf buf = m_IBCBuffer0.Y();
#else
  PelBuf buf = m_IBCBuffer.Y();
#endif

  for(int x = 0; x < width; x += 4)
  {
    for(int y = 0; y < height; y += 4)
    {
#if JVET_Z0153_IBC_EXT_REF
      if(buf.at((x + refTLx) % m_IBCBufferWidth, (y + refTLy) % m_IBCBufferHeight) == -1)
      {
        return false;
      }
      if(buf.at((x + 3 + refTLx) % m_IBCBufferWidth, (y + refTLy) % m_IBCBufferHeight) == -1)
      {
        return false;
      }
      if(buf.at((x + refTLx) % m_IBCBufferWidth, (y + 3 + refTLy) % m_IBCBufferHeight) == -1)
      {
        return false;
      }
      if(buf.at((x + 3 + refTLx) % m_IBCBufferWidth, (y + 3 + refTLy) % m_IBCBufferHeight) == -1)
      {
        return false;
      }
#else
      if(buf.at((x + refTLx) & (m_IBCBufferWidth - 1), y + refTLy) == -1) return false;
      if(buf.at((x + 3 + refTLx) & (m_IBCBufferWidth - 1), y + refTLy) == -1) return false;
      if(buf.at((x + refTLx) & (m_IBCBufferWidth - 1), y + 3 + refTLy) == -1) return false;
      if(buf.at((x + 3 + refTLx) & (m_IBCBufferWidth - 1), y + 3 + refTLy) == -1) return false;
#endif
    }
  }

  return true;
}

bool InterPrediction::xPredInterBlkRPR( const std::pair<int, int>& scalingRatio, const PPS& pps, const CompArea &blk, const Picture* refPic, const Mv& mv, Pel* dst, const int dstStride, const bool bi, const bool wrapRef, const ClpRng& clpRng, const int filterIndex, const bool useAltHpelIf )
{
  const ChromaFormat  chFmt = blk.chromaFormat;
  const ComponentID compID = blk.compID;
  const bool          rndRes = !bi;

  int shiftHor = MV_FRACTIONAL_BITS_INTERNAL + ::getComponentScaleX( compID, chFmt );
  int shiftVer = MV_FRACTIONAL_BITS_INTERNAL + ::getComponentScaleY( compID, chFmt );

  int width = blk.width;
  int height = blk.height;
  CPelBuf refBuf;

  const bool scaled = refPic->isRefScaled( &pps );

  if( scaled )
  {
    int row, col;
    int refPicWidth = refPic->getPicWidthInLumaSamples();
    int refPicHeight = refPic->getPicHeightInLumaSamples();

    int xFilter = filterIndex;
    int yFilter = filterIndex;
    const int rprThreshold1 = ( 1 << SCALE_RATIO_BITS ) * 5 / 4;
    const int rprThreshold2 = ( 1 << SCALE_RATIO_BITS ) * 7 / 4;
    if( filterIndex == 0 )
    {
      if( scalingRatio.first > rprThreshold2 )
      {
        xFilter = 4;
      }
      else if( scalingRatio.first > rprThreshold1 )
      {
        xFilter = 3;
      }

      if( scalingRatio.second > rprThreshold2 )
      {
        yFilter = 4;
      }
      else if( scalingRatio.second > rprThreshold1 )
      {
        yFilter = 3;
      }
    }
    if (filterIndex == 2)
    {
      if (isLuma(compID))
      {
        if (scalingRatio.first > rprThreshold2)
        {
          xFilter = 6;
        }
        else if (scalingRatio.first > rprThreshold1)
        {
          xFilter = 5;
        }

        if (scalingRatio.second > rprThreshold2)
        {
          yFilter = 6;
        }
        else if (scalingRatio.second > rprThreshold1)
        {
          yFilter = 5;
        }
      }
      else
      {
        if (scalingRatio.first > rprThreshold2)
        {
          xFilter = 4;
        }
        else if (scalingRatio.first > rprThreshold1)
        {
          xFilter = 3;
        }

        if (scalingRatio.second > rprThreshold2)
        {
          yFilter = 4;
        }
        else if (scalingRatio.second > rprThreshold1)
        {
          yFilter = 3;
        }
      }
    }

    const int posShift = SCALE_RATIO_BITS - 4;
    int stepX = ( scalingRatio.first + 8 ) >> 4;
    int stepY = ( scalingRatio.second + 8 ) >> 4;
    int64_t x0Int;
    int64_t y0Int;
    int offX = 1 << ( posShift - shiftHor - 1 );
    int offY = 1 << ( posShift - shiftVer - 1 );

    const int64_t posX = ( ( blk.pos().x << ::getComponentScaleX( compID, chFmt ) ) - ( pps.getScalingWindow().getWindowLeftOffset() * SPS::getWinUnitX( chFmt ) ) ) >> ::getComponentScaleX( compID, chFmt );
    const int64_t posY = ( ( blk.pos().y << ::getComponentScaleY( compID, chFmt ) ) - ( pps.getScalingWindow().getWindowTopOffset()  * SPS::getWinUnitY( chFmt ) ) ) >> ::getComponentScaleY( compID, chFmt );

    int addX = isLuma( compID ) ? 0 : int( 1 - refPic->cs->sps->getHorCollocatedChromaFlag() ) * 8 * ( scalingRatio.first - SCALE_1X.first );
    int addY = isLuma( compID ) ? 0 : int( 1 - refPic->cs->sps->getVerCollocatedChromaFlag() ) * 8 * ( scalingRatio.second - SCALE_1X.second );

    x0Int = ( ( posX << ( 4 + ::getComponentScaleX( compID, chFmt ) ) ) + mv.getHor() ) * (int64_t)scalingRatio.first + addX;
    x0Int = SIGN( x0Int ) * ( ( llabs( x0Int ) + ( (long long)1 << ( 7 + ::getComponentScaleX( compID, chFmt ) ) ) ) >> ( 8 + ::getComponentScaleX( compID, chFmt ) ) ) + ( ( refPic->getScalingWindow().getWindowLeftOffset() * SPS::getWinUnitX( chFmt ) ) << ( ( posShift - ::getComponentScaleX( compID, chFmt ) ) ) );

    y0Int = ( ( posY << ( 4 + ::getComponentScaleY( compID, chFmt ) ) ) + mv.getVer() ) * (int64_t)scalingRatio.second + addY;
    y0Int = SIGN( y0Int ) * ( ( llabs( y0Int ) + ( (long long)1 << ( 7 + ::getComponentScaleY( compID, chFmt ) ) ) ) >> ( 8 + ::getComponentScaleY( compID, chFmt ) ) ) + ( ( refPic->getScalingWindow().getWindowTopOffset() * SPS::getWinUnitY( chFmt ) ) << ( ( posShift - ::getComponentScaleY( compID, chFmt ) ) ) );

    const int extSize = isLuma( compID ) ? 1 : 2;
#if IF_12TAP
#if RPR_ENABLE
    const int iTap = 0;
#else
    const int iTap = 1;
#endif
    int vFilterSize = isLuma(compID) ? NTAPS_LUMA(iTap) : NTAPS_CHROMA;
#else
    int vFilterSize = isLuma( compID ) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif

    int yInt0 = ( (int32_t)y0Int + offY ) >> posShift;
#if IF_12TAP
    yInt0 = std::min(std::max(-(NTAPS_LUMA(iTap) / 2), yInt0), (refPicHeight >> ::getComponentScaleY(compID, chFmt)) + (NTAPS_LUMA(iTap) / 2));
#else
    yInt0 = std::min( std::max( -(NTAPS_LUMA / 2), yInt0 ), ( refPicHeight >> ::getComponentScaleY( compID, chFmt ) ) + (NTAPS_LUMA / 2) );
#endif

    int xInt0 = ( (int32_t)x0Int + offX ) >> posShift;
#if IF_12TAP
    xInt0 = std::min(std::max(-(NTAPS_LUMA(iTap) / 2), xInt0), (refPicWidth >> ::getComponentScaleX(compID, chFmt)) + (NTAPS_LUMA(iTap) / 2));
#else
    xInt0 = std::min( std::max( -(NTAPS_LUMA / 2), xInt0 ), ( refPicWidth >> ::getComponentScaleX( compID, chFmt ) ) + (NTAPS_LUMA / 2) );
#endif

    int refHeight = ((((int32_t)y0Int + (height-1) * stepY) + offY ) >> posShift) - ((((int32_t)y0Int + 0 * stepY) + offY ) >> posShift) + 1;
    refHeight = std::max<int>( 1, refHeight );

    CHECK( MAX_CU_SIZE * MAX_SCALING_RATIO + 16 < refHeight + vFilterSize - 1 + extSize, "Buffer is not large enough, increase MAX_SCALING_RATIO" );

    Pel buffer[( MAX_CU_SIZE + 16 ) * ( MAX_CU_SIZE * MAX_SCALING_RATIO + 16 )];
    int tmpStride = width;
    int xInt = 0, yInt = 0;

    for( col = 0; col < width; col++ )
    {
      int posX = (int32_t)x0Int + col * stepX;
      xInt = ( posX + offX ) >> posShift;
#if IF_12TAP
      xInt = std::min(std::max(-(NTAPS_LUMA(iTap) / 2), xInt), (refPicWidth >> ::getComponentScaleX(compID, chFmt)) + (NTAPS_LUMA(iTap) / 2));
#else
      xInt = std::min( std::max( -(NTAPS_LUMA / 2), xInt ), ( refPicWidth >> ::getComponentScaleX( compID, chFmt ) ) + (NTAPS_LUMA / 2) );
#endif
      int xFrac = ( ( posX + offX ) >> ( posShift - shiftHor ) ) & ( ( 1 << shiftHor ) - 1 );

      CHECK( xInt0 > xInt, "Wrong horizontal starting point" );

      Position offset = Position( xInt, yInt0 );
      refBuf = refPic->getRecoBuf( CompArea( compID, chFmt, offset, Size( 1, refHeight ) ), wrapRef );
      Pel* tempBuf = buffer + col;

      m_if.filterHor( compID, (Pel*)refBuf.buf - ( ( vFilterSize >> 1 ) - 1 ) * refBuf.stride, refBuf.stride, tempBuf, tmpStride, 1, refHeight + vFilterSize - 1 + extSize, xFrac, false, chFmt, clpRng, xFilter, false, useAltHpelIf && scalingRatio.first == 1 << SCALE_RATIO_BITS );
    }

    for( row = 0; row < height; row++ )
    {
      int posY = (int32_t)y0Int + row * stepY;
      yInt = ( posY + offY ) >> posShift;
#if IF_12TAP
      yInt = std::min(std::max(-(NTAPS_LUMA(iTap) / 2), yInt), (refPicHeight >> ::getComponentScaleY(compID, chFmt)) + (NTAPS_LUMA(iTap) / 2));
#else
      yInt = std::min( std::max( -(NTAPS_LUMA / 2), yInt ), ( refPicHeight >> ::getComponentScaleY( compID, chFmt ) ) + (NTAPS_LUMA / 2) );
#endif
      int yFrac = ( ( posY + offY ) >> ( posShift - shiftVer ) ) & ( ( 1 << shiftVer ) - 1 );

      CHECK( yInt0 > yInt, "Wrong vertical starting point" );

      Pel* tempBuf = buffer + ( yInt - yInt0 ) * tmpStride;

      JVET_J0090_SET_CACHE_ENABLE( false );
      m_if.filterVer( compID, tempBuf + ( ( vFilterSize >> 1 ) - 1 ) * tmpStride, tmpStride, dst + row * dstStride, dstStride, width, 1, yFrac, false, rndRes, chFmt, clpRng, yFilter, false, useAltHpelIf && scalingRatio.second == 1 << SCALE_RATIO_BITS );
      JVET_J0090_SET_CACHE_ENABLE( true );
    }
  }

  return scaled;
}


#if INTER_LIC
void InterPrediction::xLocalIlluComp(const PredictionUnit& pu,
                                     const ComponentID     compID,
                                     const Picture&        refPic,
                                     const Mv&             mv,
                                     const bool            biPred,
                                     PelBuf&               dstBuf
)
{
  Pel* refLeftTemplate  = m_pcLICRefLeftTemplate;
  Pel* refAboveTemplate = m_pcLICRefAboveTemplate;
  Pel* recLeftTemplate  = m_pcLICRecLeftTemplate;
  Pel* recAboveTemplate = m_pcLICRecAboveTemplate;
  int numTemplate[2] = { 0 , 0 }; // 0:Above, 1:Left
  xGetSublkTemplate(*pu.cu, compID, refPic, mv, pu.blocks[compID].width, pu.blocks[compID].height, 0, 0, numTemplate, refLeftTemplate, refAboveTemplate, recLeftTemplate, recAboveTemplate);

  int shift = 0, scale = 0, offset = 0;
  xGetLICParamGeneral(*pu.cu, compID, numTemplate, refLeftTemplate, refAboveTemplate, recLeftTemplate, recAboveTemplate, shift, scale, offset);

  const ClpRng& clpRng = pu.cu->cs->slice->clpRng(compID);
  dstBuf.linearTransform(scale, shift, offset, true, clpRng);
}

void InterPrediction::xGetSublkTemplate(const CodingUnit& cu,
                                        const ComponentID compID,
                                        const Picture&    refPic,
                                        const Mv&         mv,
                                        const int         sublkWidth,
                                        const int         sublkHeight,
                                        const int         posW,
                                        const int         posH,
                                        int*              numTemplate,
                                        Pel*              refLeftTemplate,
                                        Pel*              refAboveTemplate,
                                        Pel*              recLeftTemplate,
                                        Pel*              recAboveTemplate)
{
  const int       bitDepth = cu.cs->sps->getBitDepth(toChannelType(compID));
  const int       precShift = std::max(0, bitDepth - 12);

  const Picture&  currPic = *cu.cs->picture;
  const CodingUnit* const cuAbove = cu.cs->getCU(cu.blocks[compID].pos().offset(0, -1), toChannelType(compID));
  const CodingUnit* const cuLeft = cu.cs->getCU(cu.blocks[compID].pos().offset(-1, 0), toChannelType(compID));
  const CPelBuf recBuf = cuAbove || cuLeft ? currPic.getRecoBuf(cu.cs->picture->blocks[compID]) : CPelBuf();
  const CPelBuf refBuf = cuAbove || cuLeft ? refPic.getRecoBuf(refPic.blocks[compID]) : CPelBuf();

  std::vector<Pel>& invLUT = m_pcReshape->getInvLUT();

  // above
  if (cuAbove && posH == 0)
  {
    xGetPredBlkTpl<true>(cu, compID, refBuf, mv, posW, posH, sublkWidth, refAboveTemplate);
    const Pel*    rec = recBuf.bufAt(cu.blocks[compID].pos().offset(0, -1));

    for (int k = posW; k < posW + sublkWidth; k++)
    {
      int refVal = refAboveTemplate[k];
      int recVal = rec[k];

      if (isLuma(compID) && cu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
      {
        recVal = invLUT[recVal];
      }

      recVal >>= precShift;
      refVal >>= precShift;

      refAboveTemplate[k] = refVal;
      recAboveTemplate[k] = recVal;
      numTemplate[0]++;
    }
  }

  // left
  if (cuLeft && posW == 0)
  {
    xGetPredBlkTpl<false>(cu, compID, refBuf, mv, posW, posH, sublkHeight, refLeftTemplate);
    const Pel*    rec = recBuf.bufAt(cu.blocks[compID].pos().offset(-1, 0));

    for (int k = posH; k < posH + sublkHeight; k++)
    {
      int refVal = refLeftTemplate[k];
      int recVal = rec[recBuf.stride * k];

      if (isLuma(compID) && cu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
      {
        recVal = invLUT[recVal];
      }

      recVal >>= precShift;
      refVal >>= precShift;

      refLeftTemplate[k] = refVal;
      recLeftTemplate[k] = recVal;
      numTemplate[1]++;
    }
  }
}

void InterPrediction::xGetLICParamGeneral(const CodingUnit& cu,
                                          const ComponentID compID,
                                          int*              numTemplate,
                                          Pel*              refLeftTemplate,
                                          Pel*              refAboveTemplate,
                                          Pel*              recLeftTemplate,
                                          Pel*              recAboveTemplate,
                                          int&              shift,
                                          int&              scale,
                                          int&              offset
)
{
  const int       cuWidth = cu.blocks[compID].width;
  const int       cuHeight = cu.blocks[compID].height;

  const int       bitDepth = cu.cs->sps->getBitDepth(toChannelType(compID));
  const int       precShift = std::max(0, bitDepth - 12);
  const int       maxNumMinus1 = 30 - 2 * std::min(bitDepth, 12) - 1;
  const int       minDimBit = floorLog2(std::min(cuHeight, cuWidth));
  const int       minDim = 1 << minDimBit;
  int       minStepBit = minDim > 8 ? 1 : 0;
  while (minDimBit > minStepBit + maxNumMinus1) { minStepBit++; } //make sure log2(2*minDim/tmpStep) + 2*min(bitDepth,12) <= 30
  const int       numSteps = minDim >> minStepBit;
  const int       dimShift = minDimBit - minStepBit;

  //----- get correlation data -----
  int x = 0, y = 0, xx = 0, xy = 0, cntShift = 0;

  // above
  if (numTemplate[0] != 0)
  {
    for (int k = 0; k < numSteps; k++)
    {
      CHECK(((k * cuWidth) >> dimShift) >= cuWidth, "Out of range");

      int refVal = refAboveTemplate[((k * cuWidth) >> dimShift)];
      int recVal = recAboveTemplate[((k * cuWidth) >> dimShift)];

      x += refVal;
      y += recVal;
      xx += refVal * refVal;
      xy += refVal * recVal;
    }

    cntShift = dimShift;
  }

  // left
  if (numTemplate[1] != 0)
  {
    for (int k = 0; k < numSteps; k++)
    {
      CHECK(((k * cuHeight) >> dimShift) >= cuHeight, "Out of range");

      int refVal = refLeftTemplate[((k * cuHeight) >> dimShift)];
      int recVal = recLeftTemplate[((k * cuHeight) >> dimShift)];

      x += refVal;
      y += recVal;
      xx += refVal * refVal;
      xy += refVal * recVal;
    }

    cntShift += (cntShift ? 1 : dimShift);
  }

  //----- determine scale and offset -----
  shift = m_LICShift;
  if (cntShift == 0)
  {
    scale = (1 << shift);
    offset = 0;
    return;
  }

  const int cropShift = std::max(0, bitDepth - precShift + cntShift - 15);
  const int xzOffset = (xx >> m_LICRegShift);
  const int sumX = x << precShift;
  const int sumY = y << precShift;
  const int sumXX = ((xx + xzOffset) >> (cropShift << 1)) << cntShift;
  const int sumXY = ((xy + xzOffset) >> (cropShift << 1)) << cntShift;
  const int sumXsumX = (x >> cropShift) * (x >> cropShift);
  const int sumXsumY = (x >> cropShift) * (y >> cropShift);
  int a1 = sumXY - sumXsumY;
  int a2 = sumXX - sumXsumX;
  int scaleShiftA2 = getMSB(abs(a2)) - 6;
  int scaleShiftA1 = scaleShiftA2 - m_LICShiftDiff;
  scaleShiftA2 = std::max(0, scaleShiftA2);
  scaleShiftA1 = std::max(0, scaleShiftA1);
  const int scaleShiftA = scaleShiftA2 + 15 - shift - scaleShiftA1;
  a1 = a1 >> scaleShiftA1;
  a2 = Clip3(0, 63, a2 >> scaleShiftA2);
  scale = int((int64_t(a1) * int64_t(m_LICMultApprox[a2])) >> scaleShiftA);
  scale = Clip3(0, 1 << (shift + 2), scale);
  const int maxOffset = (1 << (bitDepth - 1)) - 1;
  const int minOffset = -1 - maxOffset;
  offset = (sumY - ((scale * sumX) >> shift) + ((1 << (cntShift)) >> 1)) >> cntShift;
  offset = Clip3(minOffset, maxOffset, offset);
}

template <bool trueAfalseL>
void InterPrediction::xGetPredBlkTpl(const CodingUnit& cu, const ComponentID compID, const CPelBuf& refBuf, const Mv& mv, const int posW, const int posH, const int tplSize, Pel* predBlkTpl
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
                      , bool AML
#endif
                                     )
{
  const int lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
  const int horShift  = (lumaShift + ::getComponentScaleX(compID, cu.chromaFormat));
  const int verShift  = (lumaShift + ::getComponentScaleY(compID, cu.chromaFormat));

  const int xInt      = mv.getHor() >> horShift;
  const int yInt      = mv.getVer() >> verShift;
  const int xFrac     = mv.getHor() & ((1 << horShift) - 1);
  const int yFrac     = mv.getVer() & ((1 << verShift) - 1);

  const Pel* ref;
        Pel* dst;
        int refStride, dstStride, bw, bh;
  if( trueAfalseL )
  {
    ref       = refBuf.bufAt(cu.blocks[compID].pos().offset(xInt + posW, yInt + posH - 1));
    dst       = predBlkTpl + posW;
    refStride = refBuf.stride;
    dstStride = tplSize;
    bw        = tplSize;
    bh        = 1;
  }
  else
  {
    ref       = refBuf.bufAt(cu.blocks[compID].pos().offset(xInt + posW - 1, yInt + posH));
    dst       = predBlkTpl + posH;
    refStride = refBuf.stride;
    dstStride = 1;
    bw        = 1;
    bh        = tplSize;
  }

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  int   nFilterIdx =  AML  ? 1 : 0;
#else
  const int  nFilterIdx   = 0;
#endif
  const bool useAltHpelIf = false;

  if ( yFrac == 0 )
  {
    m_if.filterHor( compID, (Pel*) ref, refStride, dst, dstStride, bw, bh, xFrac, true, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
  }
  else if ( xFrac == 0 )
  {
    m_if.filterVer( compID, (Pel*) ref, refStride, dst, dstStride, bw, bh, yFrac, true, true, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
  }
  else
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if IF_12TAP
    int vFilterSize = isLuma(compID) ? NTAPS_LUMA(0) : NTAPS_CHROMA;
#else
    int vFilterSize = isLuma(compID) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif
    if (isLuma(compID) && nFilterIdx == 1)
    {
      vFilterSize = NTAPS_BILINEAR;
    }
#else
#if IF_12TAP
    const int vFilterSize = isLuma(compID) ? NTAPS_LUMA(0) : NTAPS_CHROMA;
#else
    const int vFilterSize = isLuma(compID) ? NTAPS_LUMA : NTAPS_CHROMA;
#endif
#endif
    PelBuf tmpBuf = PelBuf(m_filteredBlockTmp[0][compID], Size(bw, bh+vFilterSize-1));

    m_if.filterHor( compID, (Pel*)ref - ((vFilterSize>>1) -1)*refStride, refStride, tmpBuf.buf, tmpBuf.stride, bw, bh+vFilterSize-1, xFrac, false, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
    JVET_J0090_SET_CACHE_ENABLE( false );
    m_if.filterVer( compID, tmpBuf.buf + ((vFilterSize>>1) -1)*tmpBuf.stride, tmpBuf.stride, dst, dstStride, bw, bh, yFrac, false, true, cu.chromaFormat, cu.slice->clpRng(compID), nFilterIdx, false, useAltHpelIf);
    JVET_J0090_SET_CACHE_ENABLE( true );
  }
}
#endif // INTER_LIC

#if TM_AMVP || TM_MRG
Distortion InterPrediction::deriveTMMv(const PredictionUnit& pu, bool fillCurTpl, Distortion curBestCost, RefPicList eRefList, int refIdx, int maxSearchRounds, Mv& mv, const MvField* otherMvf)
{
  CHECK(refIdx < 0, "Invalid reference index for TM");
  const CodingUnit& cu   = *pu.cu;
#if JVET_Z0084_IBC_TM
#if JVET_Y0128_NON_CTC
  if ( !CU::isIBC(cu) && cu.slice->getRefPic(eRefList, refIdx)->isRefScaled( pu.cs->pps ) )
  {
    return std::numeric_limits<Distortion>::max();
  }
#endif
  CHECK(CU::isIBC(cu) && otherMvf != nullptr, "IBC TM for bidir is not allowed.");
  const Picture& refPic  = CU::isIBC(cu) ? *cu.slice->getPic() : *cu.slice->getRefPic(eRefList, refIdx)->unscaledPic;
#else
#if JVET_Y0128_NON_CTC
  if ( cu.slice->getRefPic(eRefList, refIdx)->isRefScaled( pu.cs->pps ) )
  {
    return std::numeric_limits<Distortion>::max();
  }
#endif
  const Picture& refPic  = *cu.slice->getRefPic(eRefList, refIdx)->unscaledPic;
#endif
  bool doSimilarityCheck = otherMvf == nullptr ? false : cu.slice->getRefPOC((RefPicList)eRefList, refIdx) == cu.slice->getRefPOC((RefPicList)(1 - eRefList), otherMvf->refIdx);

  InterPredResources interRes(m_pcReshape, m_pcRdCost, m_if, m_filteredBlockTmp[0][COMPONENT_Y]
                           ,  m_filteredBlock[3][1][0], m_filteredBlock[3][0][0]
  );
  TplMatchingCtrl tplCtrl(pu, interRes, refPic, fillCurTpl, COMPONENT_Y, true, maxSearchRounds, m_pcCurTplAbove, m_pcCurTplLeft, m_pcRefTplAbove, m_pcRefTplLeft, mv, (doSimilarityCheck ? &(otherMvf->mv) : nullptr), curBestCost);
  if (!tplCtrl.getTemplatePresentFlag())
  {
    return std::numeric_limits<Distortion>::max();
  }

  if (otherMvf == nullptr) // uni prediction
  {
    tplCtrl.deriveMvUni<TM_TPL_SIZE>();
    mv = tplCtrl.getFinalMv();
    return tplCtrl.getMinCost();
  }
  else // bi prediction
  {
#if JVET_Y0128_NON_CTC
    if ( cu.slice->getRefPic((RefPicList)(1 - eRefList), otherMvf->refIdx)->isRefScaled(pu.cs->pps) )
    {
      return std::numeric_limits<Distortion>::max();
    }
#endif
    const Picture& otherRefPic = *cu.slice->getRefPic((RefPicList)(1-eRefList), otherMvf->refIdx)->unscaledPic;
    tplCtrl.removeHighFreq<TM_TPL_SIZE>(otherRefPic, otherMvf->mv, getBcwWeight(cu.BcwIdx, eRefList));
    tplCtrl.deriveMvUni<TM_TPL_SIZE>();
    mv = tplCtrl.getFinalMv();

    int8_t intWeight = getBcwWeight(cu.BcwIdx, eRefList);
    return (tplCtrl.getMinCost() * intWeight + (g_BcwWeightBase >> 1)) >> g_BcwLog2WeightBase;
  }
}

#if TM_MRG
void InterPrediction::deriveTMMv(PredictionUnit& pu)
{
  if( !pu.tmMergeFlag )
  {
    return;
  }

  Distortion minCostUni[NUM_REF_PIC_LIST_01] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max() };

  for (int iRefList = 0; iRefList < ( pu.cu->slice->isInterB() ? NUM_REF_PIC_LIST_01 : 1 ) ; ++iRefList)
  {
    if (pu.interDir & (iRefList + 1))
    {
      minCostUni[iRefList] = deriveTMMv(pu, true, std::numeric_limits<Distortion>::max(), (RefPicList)iRefList, pu.refIdx[iRefList], TM_MAX_NUM_OF_ITERATIONS, pu.mv[iRefList]);
    }
  }

  if (pu.cu->slice->isInterB() && pu.interDir == 3
#if MULTI_PASS_DMVR
    && !PU::checkBDMVRCondition(pu)
#endif
    )
  {
    if (minCostUni[0] == std::numeric_limits<Distortion>::max() || minCostUni[1] == std::numeric_limits<Distortion>::max())
    {
      return;
    }

    RefPicList eTargetPicList = (minCostUni[0] <= minCostUni[1]) ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
    MvField    mvfBetterUni(pu.mv[1 - eTargetPicList], pu.refIdx[1 - eTargetPicList]);
    Distortion minCostBi = deriveTMMv(pu, true, std::numeric_limits<Distortion>::max(), eTargetPicList, pu.refIdx[eTargetPicList], TM_MAX_NUM_OF_ITERATIONS, pu.mv[eTargetPicList], &mvfBetterUni);

    if (minCostBi > (minCostUni[1 - eTargetPicList] + (minCostUni[1 - eTargetPicList] >> 3)))
    {
      pu.interDir = 1 + (1 - eTargetPicList);
      pu.mv    [eTargetPicList] = Mv();
      pu.refIdx[eTargetPicList] = NOT_VALID;
    }
  }
}
#endif // TM_MRG
#endif // TM_AMVP || TM_MRG

#if TM_AMVP || TM_MRG
TplMatchingCtrl::TplMatchingCtrl( const PredictionUnit&     pu,
                                        InterPredResources& interRes,
                                  const Picture&            refPic,
                                  const bool                fillCurTpl,
                                  const ComponentID         compID,
                                  const bool                useWeight,
                                  const int                 maxSearchRounds,
                                        Pel*                curTplAbove,
                                        Pel*                curTplLeft,
                                        Pel*                refTplAbove,
                                        Pel*                refTplLeft,
                                  const Mv&                 mvStart,
                                  const Mv*                 otherRefListMv,
                                  const Distortion          curBestCost
)
: m_cu              (*pu.cu)
, m_pu              (pu)
, m_interRes        (interRes)
, m_refPic          (refPic)
, m_mvStart         (mvStart)
, m_mvFinal         (mvStart)
, m_otherRefListMv  (otherRefListMv)
, m_minCost         (curBestCost)
, m_useWeight       (useWeight)
, m_maxSearchRounds (maxSearchRounds)
, m_compID          (compID)
{
#if JVET_Z0067_RPR_ENABLE
  if ( refPic.isRefScaled(pu.cs->pps) )
  {
    return;
  }
#endif
  // Initialization
  const bool tplAvalableAbove = xFillCurTemplate<TM_TPL_SIZE, true >((fillCurTpl ? curTplAbove : nullptr));
  const bool tplAvalableLeft  = xFillCurTemplate<TM_TPL_SIZE, false>((fillCurTpl ? curTplLeft  : nullptr));
  m_curTplAbove = tplAvalableAbove ? PelBuf(curTplAbove, pu.lwidth(),   TM_TPL_SIZE ) : PelBuf();
  m_curTplLeft  = tplAvalableLeft  ? PelBuf(curTplLeft , TM_TPL_SIZE,   pu.lheight()) : PelBuf();
  m_refTplAbove = tplAvalableAbove ? PelBuf(refTplAbove, m_curTplAbove              ) : PelBuf();
  m_refTplLeft  = tplAvalableLeft  ? PelBuf(refTplLeft , m_curTplLeft               ) : PelBuf();
#if JVET_X0056_DMVD_EARLY_TERMINATION
  m_earlyTerminateTh = TM_TPL_SIZE * ((tplAvalableAbove ? m_pu.lwidth() : 0) + (tplAvalableLeft ? m_pu.lheight() : 0));
#endif

  // Pre-interpolate samples on search area
#if JVET_Z0084_IBC_TM
  m_refSrAbove = tplAvalableAbove && maxSearchRounds > 0 && !CU::isIBC(m_cu) ? PelBuf(interRes.m_preFillBufA, m_curTplAbove.width + 2 * TM_SEARCH_RANGE, m_curTplAbove.height + 2 * TM_SEARCH_RANGE) : PelBuf();
#else
  m_refSrAbove = tplAvalableAbove && maxSearchRounds > 0 ? PelBuf(interRes.m_preFillBufA, m_curTplAbove.width + 2 * TM_SEARCH_RANGE, m_curTplAbove.height + 2 * TM_SEARCH_RANGE) : PelBuf();
#endif
  if (m_refSrAbove.buf != nullptr)
  {
    m_refSrAbove = xGetRefTemplate<TM_TPL_SIZE, true, TM_SEARCH_RANGE>(m_pu, m_refPic, mvStart, m_refSrAbove);
    m_refSrAbove = m_refSrAbove.subBuf(Position(TM_SEARCH_RANGE, TM_SEARCH_RANGE), m_curTplAbove);
  }

#if JVET_Z0084_IBC_TM
  m_refSrLeft  = tplAvalableLeft  && maxSearchRounds > 0 && !CU::isIBC(m_cu) ? PelBuf(interRes.m_preFillBufL, m_curTplLeft .width + 2 * TM_SEARCH_RANGE, m_curTplLeft .height + 2 * TM_SEARCH_RANGE) : PelBuf();
#else
  m_refSrLeft  = tplAvalableLeft  && maxSearchRounds > 0 ? PelBuf(interRes.m_preFillBufL, m_curTplLeft .width + 2 * TM_SEARCH_RANGE, m_curTplLeft .height + 2 * TM_SEARCH_RANGE) : PelBuf();
#endif
  if (m_refSrLeft.buf != nullptr)
  {
    m_refSrLeft = xGetRefTemplate<TM_TPL_SIZE, false, TM_SEARCH_RANGE>(m_pu, m_refPic, mvStart, m_refSrLeft);
    m_refSrLeft = m_refSrLeft.subBuf(Position(TM_SEARCH_RANGE, TM_SEARCH_RANGE), m_curTplLeft);
  }
}

int TplMatchingCtrl::getDeltaMean(const PelBuf& bufCur, const PelBuf& bufRef, const int rowSubShift, const int bd)
{
  int64_t deltaSum = g_pelBufOP.getSumOfDifference(bufCur.buf, bufCur.stride, bufRef.buf, bufRef.stride, bufCur.width, bufCur.height, rowSubShift, bd);
  return int(deltaSum / (int64_t)bufCur.area());
}

template <int tplSize>
void TplMatchingCtrl::deriveMvUni()
{
  if (m_minCost == std::numeric_limits<Distortion>::max())
  {
    m_minCost = xGetTempMatchError<tplSize>(m_mvStart);
  }

  if (m_maxSearchRounds <= 0)
  {
    return;
  }

  int searchStepShift = (m_cu.imv == IMV_4PEL ? MV_FRACTIONAL_BITS_INTERNAL + 2 : MV_FRACTIONAL_BITS_INTERNAL);
  xRefineMvSearch<tplSize, TplMatchingCtrl::TMSEARCH_DIAMOND>(m_maxSearchRounds, searchStepShift);
  xRefineMvSearch<tplSize, TplMatchingCtrl::TMSEARCH_CROSS  >(                1, searchStepShift);
  xRefineMvSearch<tplSize, TplMatchingCtrl::TMSEARCH_CROSS  >(                1, searchStepShift - 1);
#if MULTI_PASS_DMVR
  if (!m_pu.bdmvrRefine)
  {
#endif
  xRefineMvSearch<tplSize, TplMatchingCtrl::TMSEARCH_CROSS  >(                1, searchStepShift - 2);
  xRefineMvSearch<tplSize, TplMatchingCtrl::TMSEARCH_CROSS  >(                1, searchStepShift - 3);
#if MULTI_PASS_DMVR
  }
  else
  {
    xDeriveCostBasedMv<TplMatchingCtrl::TMSEARCH_CROSS>();
  }
#endif
}

template <int tplSize>
void TplMatchingCtrl::removeHighFreq(const Picture& otherRefPic, const Mv& otherRefMv, const uint8_t curRefBcwWeight)
{
  xRemoveHighFreq<tplSize,  true>(otherRefPic, otherRefMv, curRefBcwWeight);
  xRemoveHighFreq<tplSize, false>(otherRefPic, otherRefMv, curRefBcwWeight);
}

template <int tplSize, bool trueAfalseL>
bool TplMatchingCtrl::xFillCurTemplate(Pel* tpl)
{
  const Position          posOffset = trueAfalseL ? Position(0, -tplSize) : Position(-tplSize, 0);
  const CodingUnit* const cuNeigh   = m_cu.cs->getCU(m_pu.blocks[m_compID].pos().offset(posOffset), toChannelType(m_compID));

  if (cuNeigh == nullptr)
  {
    return false;
  }

  if (tpl == nullptr)
  {
    return true;
  }

#if JVET_Z0084_IBC_TM
  // Stay in reference region for IBC
  if( CU::isIBC(m_cu) )
  {
    const int cuPelX       = m_pu.lx();
    const int cuPelY       = m_pu.ly();
    const int roiWidth     = trueAfalseL ? m_pu.lwidth() : tplSize;
    const int roiHeight    = trueAfalseL ? tplSize       : m_pu.lheight();
    const int picWidth     = m_pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
    const int picHeight    = m_pu.cs->slice->getPPS()->getPicHeightInLumaSamples();
    const uint32_t ctuSize = m_pu.cs->slice->getSPS()->getMaxCUWidth();
    const Mv  tempBv       = trueAfalseL ? Mv(0, -tplSize) : Mv(-tplSize, 0);

    if (!PU::searchBv(m_pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, tempBv.getHor(), tempBv.getVer(), ctuSize))
    {
      return false;
    }
  }
#endif

  const Picture&          currPic = *m_cu.cs->picture;
  const CPelBuf           recBuf  = currPic.getRecoBuf(m_cu.cs->picture->blocks[m_compID]);
        std::vector<Pel>& invLUT  = m_interRes.m_pcReshape->getInvLUT();
  const bool              useLUT  = isLuma(m_compID) && m_cu.cs->picHeader->getLmcsEnabledFlag() && m_interRes.m_pcReshape->getCTUFlag();
#if JVET_W0097_GPM_MMVD_TM & TM_MRG
  if (m_cu.geoFlag)
  {
    CHECK(m_pu.geoTmType == GEO_TM_OFF, "invalid geo template type value");
    if (m_pu.geoTmType == GEO_TM_SHAPE_A)
    {
      if (trueAfalseL == 0)
      {
        return false;
      }
    }
    if (m_pu.geoTmType == GEO_TM_SHAPE_L)
    {
      if (trueAfalseL == 1)
      {
        return false;
      }
    }
  }
#endif
  const Size dstSize = (trueAfalseL ? Size(m_pu.lwidth(), tplSize) : Size(tplSize, m_pu.lheight()));
  for (int h = 0; h < (int)dstSize.height; h++)
  {
    const Position recPos = trueAfalseL ? Position(0, -tplSize + h) : Position(-tplSize, h);
    const Pel*     rec    = recBuf.bufAt(m_pu.blocks[m_compID].pos().offset(recPos));
          Pel*     dst    = tpl + h * dstSize.width;

    for (int w = 0; w < (int)dstSize.width; w++)
    {
      int recVal = rec[w];
      dst[w] = useLUT ? invLUT[recVal] : recVal;
    }
  }

  return true;
}

template <int tplSize, bool trueAfalseL, int sr>
PelBuf TplMatchingCtrl::xGetRefTemplate(const PredictionUnit& curPu, const Picture& refPic, const Mv& _mv, PelBuf& dstBuf)
{
#if JVET_Z0084_IBC_TM
  // Stay in reference region for IBC
  if( CU::isIBC(m_cu) )
  {
    const int cuPelX       = m_pu.lx();
    const int cuPelY       = m_pu.ly();
    const int roiWidth     = trueAfalseL ? m_pu.lwidth() : tplSize;
    const int roiHeight    = trueAfalseL ? tplSize       : m_pu.lheight();
    const int picWidth     = m_pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
    const int picHeight    = m_pu.cs->slice->getPPS()->getPicHeightInLumaSamples();
    const uint32_t ctuSize = m_pu.cs->slice->getSPS()->getMaxCUWidth();

    Mv tempBv = _mv;
    tempBv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
    tempBv += trueAfalseL ? Mv(0, -tplSize) : Mv(-tplSize, 0);

    if (!PU::searchBv(m_pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, tempBv.getHor(), tempBv.getVer(), ctuSize))
    {
      return PelBuf();
    }
  }
#endif

  // read from pre-interpolated buffer
  PelBuf& refSrBuf = trueAfalseL ? m_refSrAbove : m_refSrLeft;
#if JVET_Z0084_IBC_TM
  if (!CU::isIBC(m_cu) && sr == 0 && refPic.getPOC() == m_refPic.getPOC() && refSrBuf.buf != nullptr)
#else
  if (sr == 0 && refPic.getPOC() == m_refPic.getPOC() && refSrBuf.buf != nullptr)
#endif
  {
    Mv mvDiff = _mv - m_mvStart;
    if ((mvDiff.getAbsHor() & ((1 << MV_FRACTIONAL_BITS_INTERNAL) - 1)) == 0 && (mvDiff.getAbsVer() & ((1 << MV_FRACTIONAL_BITS_INTERNAL) - 1)) == 0)
    {
      mvDiff >>= MV_FRACTIONAL_BITS_INTERNAL;
      if (mvDiff.getAbsHor() <= TM_SEARCH_RANGE && mvDiff.getAbsVer() <= TM_SEARCH_RANGE)
      {
        return refSrBuf.subBuf(Position(mvDiff.getHor(), mvDiff.getVer()), dstBuf);
      }
    }
  }

  // Do interpolation on the fly
  Position blkPos  = ( trueAfalseL ? Position(curPu.lx(), curPu.ly() - tplSize) : Position(curPu.lx() - tplSize, curPu.ly()) );
  Size     blkSize = Size(dstBuf.width, dstBuf.height);
  Mv       mv      = _mv - Mv(sr << MV_FRACTIONAL_BITS_INTERNAL, sr << MV_FRACTIONAL_BITS_INTERNAL);
#if JVET_Z0084_IBC_TM
  if( !CU::isIBC(m_cu) )
#endif
  clipMv( mv, blkPos, blkSize, *m_cu.cs->sps, *m_cu.cs->pps );

  const int lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
  const int horShift  = (lumaShift + ::getComponentScaleX(m_compID, m_cu.chromaFormat));
  const int verShift  = (lumaShift + ::getComponentScaleY(m_compID, m_cu.chromaFormat));

  const int xInt  = mv.getHor() >> horShift;
  const int yInt  = mv.getVer() >> verShift;
  const int xFrac = mv.getHor() & ((1 << horShift) - 1);
  const int yFrac = mv.getVer() & ((1 << verShift) - 1);

  const CPelBuf refBuf = refPic.getRecoBuf(refPic.blocks[m_compID]);
  const Pel* ref       = refBuf.bufAt(blkPos.offset(xInt, yInt));
        Pel* dst       = dstBuf.buf;
        int  refStride = refBuf.stride;
        int  dstStride = dstBuf.stride;
        int  bw        = (int)blkSize.width;
        int  bh        = (int)blkSize.height;

  const int  nFilterIdx   = 1;
  const bool useAltHpelIf = false;
  const bool biMCForDMVR  = false;

  if ( yFrac == 0 )
  {
    m_interRes.m_if.filterHor( m_compID, (Pel*) ref, refStride, dst, dstStride, bw, bh, xFrac, true, m_cu.chromaFormat, m_cu.slice->clpRng(m_compID), nFilterIdx, biMCForDMVR, useAltHpelIf );
  }
  else if ( xFrac == 0 )
  {
    m_interRes.m_if.filterVer( m_compID, (Pel*) ref, refStride, dst, dstStride, bw, bh, yFrac, true, true, m_cu.chromaFormat, m_cu.slice->clpRng(m_compID), nFilterIdx, biMCForDMVR, useAltHpelIf );
  }
  else
  {
    const int vFilterSize = isLuma(m_compID) ? NTAPS_BILINEAR : NTAPS_CHROMA;
    PelBuf tmpBuf = PelBuf(m_interRes.m_ifBuf, Size(bw, bh+vFilterSize-1));

    m_interRes.m_if.filterHor( m_compID, (Pel*)ref - ((vFilterSize>>1) -1)*refStride, refStride, tmpBuf.buf, tmpBuf.stride, bw, bh+vFilterSize-1, xFrac, false, m_cu.chromaFormat, m_cu.slice->clpRng(m_compID), nFilterIdx, biMCForDMVR, useAltHpelIf );
    JVET_J0090_SET_CACHE_ENABLE( false );
    m_interRes.m_if.filterVer( m_compID, tmpBuf.buf + ((vFilterSize>>1) -1)*tmpBuf.stride, tmpBuf.stride, dst, dstStride, bw, bh, yFrac, false, true, m_cu.chromaFormat, m_cu.slice->clpRng(m_compID), nFilterIdx, biMCForDMVR, useAltHpelIf );
    JVET_J0090_SET_CACHE_ENABLE( true );
  }

  return dstBuf;
}

template <int tplSize, bool trueAfalseL>
void TplMatchingCtrl::xRemoveHighFreq(const Picture& otherRefPic, const Mv& otherRefMv, const uint8_t curRefBcwWeight)
{
  PelBuf& curTplBuf = trueAfalseL ? m_curTplAbove : m_curTplLeft;
  PelBuf  refTplBuf = trueAfalseL ? m_refTplAbove : m_refTplLeft;

  if (curTplBuf.buf != nullptr)
  {
    refTplBuf = xGetRefTemplate<tplSize, trueAfalseL, 0>(m_pu, otherRefPic, otherRefMv, refTplBuf);
    if (curRefBcwWeight != g_BcwWeights[BCW_DEFAULT])
    {
      curTplBuf.removeWeightHighFreq(refTplBuf, false, m_cu.slice->clpRng(m_compID), curRefBcwWeight);
    }
    else
    {
      curTplBuf.removeHighFreq(refTplBuf, false, m_cu.slice->clpRng(m_compID));
    }
  }
}

template <int tplSize, int searchPattern>
void TplMatchingCtrl::xRefineMvSearch(int maxSearchRounds, int searchStepShift)
{
  static const int finestMvdPrec[NUM_IMV_MODES] = { MV_FRACTIONAL_BITS_INTERNAL - 2, MV_FRACTIONAL_BITS_INTERNAL, MV_FRACTIONAL_BITS_INTERNAL + 2, MV_FRACTIONAL_BITS_INTERNAL - 1 };
  if (searchStepShift < finestMvdPrec[m_cu.imv] && (!m_pu.mergeFlag || m_cu.imv == IMV_HPEL))
  {
    return;
  }

#if JVET_Z0084_IBC_TM
  // Limit to integer pel search for IBC
  if( CU::isIBC(m_cu) && (searchStepShift < MV_FRACTIONAL_BITS_INTERNAL) )
  {
    return;
  }
#endif
  // Search pattern configuration
  static const Mv patternCross  [4] = { Mv(0, 1), Mv(1, 0), Mv(0, -1), Mv(-1, 0) };
  static const Mv patternDiamond[8] = { Mv(0, 2), Mv(1, 1), Mv(2, 0), Mv(1, -1), Mv(0, -2), Mv(-1, -1), Mv(-2, 0), Mv(-1, 1) };

  int       directStart = 0, directEnd = 0, directRounding = 0, directMask = 0;
  const Mv *pSearchOffset = nullptr;
#if MULTI_PASS_DMVR
  Distortion *costArray   = nullptr; 
#endif
  if (searchPattern == TMSEARCH_CROSS)
  {
    directEnd      = 3;
    directRounding = 4;
    directMask     = 0x03;
    pSearchOffset  = patternCross;
#if MULTI_PASS_DMVR
    memset(m_tmCostArrayCross, -1, sizeof(m_tmCostArrayCross));
    costArray      = m_tmCostArrayCross;
    costArray[4]   = m_minCost;
#endif
  }
  else if (searchPattern == TMSEARCH_DIAMOND)
  {
    directEnd      = 7;
    directRounding = 8;
    directMask     = 0x07;
    pSearchOffset  = patternDiamond;
#if MULTI_PASS_DMVR
    memset(m_tmCostArrayDiamond, -1, sizeof(m_tmCostArrayDiamond));
    costArray      = m_tmCostArrayDiamond;
    costArray[8]   = m_minCost;
#endif
  }
  else
  {
    CHECK(true, "Unknown search method for TM");
  }

#if JVET_Z0084_IBC_TM
  const int cuPelX       = m_pu.lx();
  const int cuPelY       = m_pu.ly();
  const int roiWidth     = m_pu.lwidth();
  const int roiHeight    = m_pu.lheight();
  const int picWidth     = m_pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
  const int picHeight    = m_pu.cs->slice->getPPS()->getPicHeightInLumaSamples();
  const uint32_t ctuSize = m_pu.cs->slice->getSPS()->getMaxCUWidth();
#endif

  // Iterative search
  for (int uiRound = 0; uiRound < maxSearchRounds; uiRound++)
  {
    int directBest = -1;
    Mv  mvCurCenter(m_mvFinal);
#if JVET_X0056_DMVD_EARLY_TERMINATION
    Distortion prevMinCost = m_minCost;
#endif
    for (int nIdx = directStart; nIdx <= directEnd; nIdx++)
    {
      int nDirect = (nIdx + directRounding) & directMask;

      Mv mvOffset = pSearchOffset[nDirect];
      mvOffset  <<= searchStepShift;
      Mv mvCand   = mvCurCenter + mvOffset;
#if JVET_Z0084_IBC_TM
      // Stay in reference region for IBC
      if( CU::isIBC(m_cu) )
      {
        Mv tempBv = mvCand;
        tempBv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);

        if (!PU::searchBv(m_pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, tempBv.getHor(), tempBv.getVer(), ctuSize)
          || (m_curTplAbove.buf != nullptr && !PU::searchBv(m_pu, cuPelX, cuPelY, roiWidth, tplSize,   picWidth, picHeight, tempBv.getHor(),         tempBv.getVer()-tplSize, ctuSize))
          || (m_curTplLeft.buf  != nullptr && !PU::searchBv(m_pu, cuPelX, cuPelY, tplSize,  roiHeight, picWidth, picHeight, tempBv.getHor()-tplSize, tempBv.getVer(),         ctuSize)))
        {
          continue;
        }
      }
#endif
      Distortion cost = InterPrediction::getDecoderSideDerivedMvCost(m_mvStart, mvCand, TM_SEARCH_RANGE, DECODER_SIDE_MV_WEIGHT); // MV cost is used just for skipping search
      if (cost >= m_minCost || (m_otherRefListMv != nullptr && *m_otherRefListMv == mvCand))
      {
        continue;
      }

      cost = xGetTempMatchError<tplSize>(mvCand);
#if MULTI_PASS_DMVR
      costArray[nDirect] = cost;
#endif

      if (cost < m_minCost)
      {
        m_minCost  = cost;
        m_mvFinal  = mvCand;
        directBest = nDirect;
      }
    }

    if (directBest == -1)
    {
      break;
    }
#if JVET_X0056_DMVD_EARLY_TERMINATION
    if (uiRound > 0 && prevMinCost < m_minCost + m_earlyTerminateTh)
    {
      break;
    }
#endif

    int nStep   = searchPattern == TMSEARCH_DIAMOND ? (2 - (directBest & 0x01)) : 1;
    directStart = directBest - nStep;
    directEnd   = directBest + nStep;

#if MULTI_PASS_DMVR
    if ((uiRound + 1) < maxSearchRounds)
    {
      xNextTmCostAarray<searchPattern>(directBest);
    }
#endif
  }
}

#if MULTI_PASS_DMVR
template <int searchPattern>
void TplMatchingCtrl::xNextTmCostAarray(int bestDirect)
{
  Distortion *costLog = searchPattern == TMSEARCH_CROSS ? m_tmCostArrayCross : (searchPattern == TMSEARCH_DIAMOND ? m_tmCostArrayDiamond : nullptr);

  if (searchPattern == TMSEARCH_CROSS)
  {
    CHECK(bestDirect < 0 || bestDirect > 3, "Error: Unknown bestDirect");

    int prevCenter = (bestDirect + 2) & 0x3;
    costLog[prevCenter] = costLog[4];
    costLog[4]          = costLog[bestDirect];

    for (int offset = 1; offset < 4; ++offset)
    {
      costLog[(prevCenter + offset + 4) & 0x3] = std::numeric_limits<Distortion>::max();
    }
  }
  else if (searchPattern == TMSEARCH_DIAMOND)
  {
  }
  else
  {
    CHECK(true, "Unknown search method for TM");
  }
}

template <int searchPattern>
void TplMatchingCtrl::xDeriveCostBasedMv()
{
  if (m_minCost == 0)
  {
    return;
  }

  if (searchPattern == TMSEARCH_CROSS)
  {
    xDeriveCostBasedOffset<true >(m_tmCostArrayCross[3], m_tmCostArrayCross[4], m_tmCostArrayCross[1], 0);
    xDeriveCostBasedOffset<false>(m_tmCostArrayCross[2], m_tmCostArrayCross[4], m_tmCostArrayCross[0], 0);
  }

  else
  {
    CHECK(true, "Unknown search method for TM");
  }
}

template <bool TrueX_FalseY>
void TplMatchingCtrl::xDeriveCostBasedOffset(Distortion costLorA, Distortion costCenter, Distortion costRorB, int log2StepSize)
{
  if (!m_pu.mergeFlag || m_cu.imv != IMV_OFF)
  {
    return;
  }
  if (costLorA == std::numeric_limits<Distortion>::max() || costRorB == std::numeric_limits<Distortion>::max() || (costCenter > costLorA || costCenter > costRorB))
  {
    return;
  }

  const int extraMvFracBit = MV_FRACTIONAL_BITS_INTERNAL - 1;

  int&    mvComp      = TrueX_FalseY ? m_mvFinal.hor : m_mvFinal.ver;
  int64_t numerator   = (int64_t)(costLorA - costRorB);
  int64_t denominator = (int64_t)((costLorA + costRorB - (costCenter << 1)) << 1);
  
  if (denominator != 0)
  {
    if (costCenter != costLorA && costCenter != costRorB)
    {
      if (extraMvFracBit > 1 || log2StepSize > 1)
      {
        mvComp += xBinaryDivision(numerator, denominator, extraMvFracBit + log2StepSize);
      }
    }
    else
    {
      const int off = 1 << (extraMvFracBit - 1);
      mvComp += ((costCenter == costLorA ? -off : off) << log2StepSize);
    }
  }
}

int TplMatchingCtrl::xBinaryDivision(int64_t numerator, int64_t denominator, int fracBits)
{
  if (fracBits < 2) // Because the result of division is assumed to be less than 0.5
  {
    return 0;
  }

  int sign = 0;
  if (numerator < 0)
  {
    sign = 1;
    numerator = -numerator;
  }

  numerator   <<=  fracBits;
  denominator <<= (fracBits - 2); // This "-2" is by the assumption that the result of division is always less than 0.5

  int quotient = 0;
  for (int binIdx = 0; binIdx < fracBits - 2; ++binIdx)
  {
    if (numerator >= denominator)
    {
      numerator -= denominator;
      ++quotient;
    }
    quotient    <<= 1;
    denominator >>= 1;
  }

  if (numerator >= denominator)
  {
    ++quotient;
  }

  return sign ? -quotient : quotient;
}
#endif

template <int tplSize>
Distortion TplMatchingCtrl::xGetTempMatchError(const Mv& mv)
{
  if (!getTemplatePresentFlag())
  {
    return std::numeric_limits<Distortion>::max();
  }

  Distortion sum = 0;
  sum += xGetTempMatchError<tplSize, true >(mv);
  sum += xGetTempMatchError<tplSize, false>(mv);
  return sum;
}

template <int tplSize, bool trueAfalseL>
Distortion TplMatchingCtrl::xGetTempMatchError(const Mv& mv)
{
  PelBuf& curTplBuf = trueAfalseL ? m_curTplAbove : m_curTplLeft;
  PelBuf  refTplBuf = trueAfalseL ? m_refTplAbove : m_refTplLeft;

  if (curTplBuf.buf == nullptr)
  {
    return 0;
  }

  const int rowSubShift   = 0;
  const int bitDepth      = m_cu.slice->clpRng(m_compID).bd;

  // fetch reference template block
  refTplBuf = xGetRefTemplate<tplSize, trueAfalseL, 0>(m_pu, m_refPic, mv, refTplBuf);
#if JVET_Z0084_IBC_TM
  if (refTplBuf.buf == nullptr)
  {
    return std::numeric_limits<Distortion>::max();
  }
#endif

  // compute matching cost
  Distortion partSum = 0;
  if (m_useWeight)
  {
    DistParam cDistParam;
    cDistParam.applyWeight = false;
#if INTER_LIC
    cDistParam.useMR = m_cu.LICFlag;
#endif
    int tmWeightIdx  = (m_pu.lwidth() >= TM_MIN_CU_SIZE_FOR_ALT_WEIGHTED_COST && m_pu.lheight() >= TM_MIN_CU_SIZE_FOR_ALT_WEIGHTED_COST ? 1 : 0);
    m_interRes.m_pcRdCost->setDistParam( cDistParam, curTplBuf, refTplBuf, bitDepth, trueAfalseL, tmWeightIdx, rowSubShift, m_compID );
    CHECK(TM_TPL_SIZE != 4, "The distortion function of template matching is implemetned currently only for size=4.");
    partSum = cDistParam.distFunc( cDistParam );
  }
  else
  {
    DistParam cDistParam;
    cDistParam.applyWeight = false;
#if INTER_LIC
    cDistParam.useMR = m_cu.LICFlag;
#endif
    m_interRes.m_pcRdCost->setDistParam( cDistParam, curTplBuf, refTplBuf, bitDepth, m_compID, false );
    cDistParam.subShift = rowSubShift;

    partSum = cDistParam.distFunc( cDistParam );
#if FULL_NBIT
    partSum >>= (bitDepth > 8 ? bitDepth - 8 : 0);
#endif
  }

  return partSum;
}
#endif // TM_AMVP || TM_MRG

#if TM_AMVP || TM_MRG || MULTI_PASS_DMVR
Distortion InterPrediction::getDecoderSideDerivedMvCost(const Mv& mvStart, const Mv& mvCur, int searchRangeInFullPel, int weight)
{
  int searchRange = searchRangeInFullPel << MV_FRACTIONAL_BITS_INTERNAL;
  Mv mvDist = mvStart - mvCur;
  Distortion cost = std::numeric_limits<Distortion>::max();
  if (mvDist.getAbsHor() <= searchRange && mvDist.getAbsVer() <= searchRange)
  {
    cost = (mvDist.getAbsHor() + mvDist.getAbsVer()) * weight;
    cost >>= MV_FRACTIONAL_BITS_DIFF;
  }

  return cost;
}

#if MULTI_PASS_DMVR
void InterPrediction::xBDMVRUpdateSquareSearchCostLog( Distortion* costLog, int bestDirect )
{
  CHECK(bestDirect < 0 || bestDirect > 7, "Error: Unknown bestDirect");

  int prevCenter      = ( bestDirect + 4 ) & 0x7;
  costLog[prevCenter] = costLog[8];
  costLog[8]          = costLog[bestDirect];

  if( prevCenter & 0x1 )
  {
    costLog[( prevCenter - 1 + 8 ) & 0x7] = costLog[( prevCenter - 2 + 8 ) & 0x7];
    costLog[( prevCenter + 1 + 8 ) & 0x7] = costLog[( prevCenter + 2 + 8 ) & 0x7];
    costLog[( prevCenter - 2 + 8 ) & 0x7] = costLog[( prevCenter - 3 + 8 ) & 0x7];
    costLog[( prevCenter + 2 + 8 ) & 0x7] = costLog[( prevCenter + 3 + 8 ) & 0x7];
    for( int offset = 3 ; offset < 6 ; ++ offset )
    {
      costLog[( prevCenter + offset + 8 ) & 0x7] = std::numeric_limits<Distortion>::max();
    }
  }
  else
  {
    costLog[( prevCenter - 1 + 8 ) & 0x7] = costLog[( prevCenter - 3 + 8 ) & 0x7];
    costLog[( prevCenter + 1 + 8 ) & 0x7] = costLog[( prevCenter + 3 + 8 ) & 0x7];
    for( int offset = 2 ; offset < 7 ; ++ offset )
    {
      costLog[( prevCenter + offset + 8 ) & 0x7] = std::numeric_limits<Distortion>::max();
    }
  }
}
#endif
#endif

#if TM_AMVP
void InterPrediction::clearTplAmvpBuffer()
{
  for (int imv = 0; imv < NUM_IMV_MODES; ++imv)
  {
    for (int refIdx = 0; refIdx < MAX_NUM_REF; ++refIdx)
    {
      m_tplAmvpInfo   [imv][0][refIdx] = AMVPInfo();
      m_tplAmvpInfo   [imv][1][refIdx] = AMVPInfo();
#if INTER_LIC
      m_tplAmvpInfoLIC[imv][0][refIdx] = AMVPInfo();
      m_tplAmvpInfoLIC[imv][1][refIdx] = AMVPInfo();
#endif
    }
  }
}

void InterPrediction::writeTplAmvpBuffer(const AMVPInfo& src, const CodingUnit& cu, RefPicList eRefList, int refIdx)
{
#if INTER_LIC
  AMVPInfo& dst = (cu.LICFlag ? m_tplAmvpInfoLIC : m_tplAmvpInfo)[cu.imv][eRefList][refIdx];
#else
  AMVPInfo& dst = m_tplAmvpInfo[cu.imv][eRefList][refIdx];
#endif
  dst = src;
}

bool InterPrediction::readTplAmvpBuffer(AMVPInfo& dst, const CodingUnit& cu, RefPicList eRefList, int refIdx)
{
#if INTER_LIC
  AMVPInfo& src = (cu.LICFlag ? m_tplAmvpInfoLIC : m_tplAmvpInfo)[cu.imv][eRefList][refIdx];
#else
  AMVPInfo& src = m_tplAmvpInfo[cu.imv][eRefList][refIdx];
#endif
  if (src.numCand > 0)
  {
    dst = src;
    return true;
  }
  return false;
}
#endif

#if MULTI_PASS_DMVR
#if JVET_X0049_ADAPT_DMVR
bool InterPrediction::processBDMVRPU2Dir(PredictionUnit& pu, bool subPURefine[2], Mv(&finalMvDir)[2])
{
  const int lumaArea = pu.lumaSize().area();
  bool       bUseMR = lumaArea > 64;
#if JVET_Y0089_DMVR_BCW
  bUseMR |= (pu.cu->BcwIdx != BCW_DEFAULT);
#endif
  subPURefine[0] = subPURefine[1] = true;

  Distortion minCost = std::numeric_limits<Distortion>::max();
  Mv         mvInitial_PU[2] = { pu.mv[0], pu.mv[1] };
  Mv         mvFinal[2] = { pu.mv[0], pu.mv[1] };

  Distortion initCost = xBDMVRGetMatchingError(pu, mvInitial_PU, bUseMR, false);
  if (initCost < lumaArea)
  {
    subPURefine[0] = false;
    subPURefine[1] = false;
    finalMvDir[0] = mvFinal[0];
    finalMvDir[1] = mvFinal[1];
    return false;
  }

  minCost = xBDMVRMvOneTemplateHPelSquareSearch<1>(mvFinal, initCost, pu, mvInitial_PU, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, bUseMR, false);
  subPURefine[0] = minCost >= lumaArea;
  finalMvDir[0] = mvFinal[0];

  mvFinal[0] = mvInitial_PU[0];
  mvFinal[1] = mvInitial_PU[1];
  minCost = xBDMVRMvOneTemplateHPelSquareSearch<2>(mvFinal, initCost, pu, mvInitial_PU, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, bUseMR, false);
  subPURefine[1] = minCost >= lumaArea;
  finalMvDir[1] = mvFinal[1];

  return true;
}

void InterPrediction::processBDMVRSubPU(PredictionUnit& pu, bool subPURefine)
{
  if (!subPURefine)
  {
    // span motion to subPU
    const int dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
    const int dx = std::min<int>(pu.lumaSize().width, DMVR_SUBCU_WIDTH);
    Position puPos = pu.lumaPos();

    int subPuIdx = 0;
    const int dmvrSubPuStrideIncr = DMVR_SUBPU_STRIDE - std::max(1, (int)(pu.lumaSize().width >> DMVR_SUBCU_WIDTH_LOG2));
    for (int y = puPos.y, yStart = 0; y < (puPos.y + pu.lumaSize().height); y = y + dy, yStart = yStart + dy)
    {
      for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + dx, xStart = xStart + dx)
      {
        m_bdmvrSubPuMvBuf[REF_PIC_LIST_0][subPuIdx] = pu.mv[0];
        m_bdmvrSubPuMvBuf[REF_PIC_LIST_1][subPuIdx] = pu.mv[1];
        subPuIdx++;
      }
      subPuIdx += dmvrSubPuStrideIncr;
    }
    return;
  }

  const int dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
  const int dx = std::min<int>(pu.lumaSize().width, DMVR_SUBCU_WIDTH);
  Position puPos = pu.lumaPos();
  PredictionUnit subPu = pu;

  int subPuIdx = 0;
  const int dmvrSubPuStrideIncr = DMVR_SUBPU_STRIDE - std::max(1, (int)(pu.lumaSize().width >> DMVR_SUBCU_WIDTH_LOG2));
  Distortion minCost = std::numeric_limits<Distortion>::max();
  const Mv   mvInitial[2] = { pu.mv[0], pu.mv[1] };
  Mv         mvFinal[2] = { pu.mv[0], pu.mv[1] };
  Mv         mvOffset;

  const Distortion earlyTerminateTh = dx * dy;
  const int adaptiveSearchRangeHor = (dx >> 1) < BDMVR_INTME_RANGE ? (dx >> 1) : BDMVR_INTME_RANGE;
  const int adaptiveSearchRangeVer = (dy >> 1) < BDMVR_INTME_RANGE ? (dy >> 1) : BDMVR_INTME_RANGE;
  const bool adaptRange = (adaptiveSearchRangeHor != BDMVR_INTME_RANGE || adaptiveSearchRangeVer != BDMVR_INTME_RANGE);
  const int maxSearchRound = pu.bmMergeFlag ? BM_MRG_SUB_PU_INT_MAX_SRCH_ROUND : BDMVR_INTME_FULL_SEARCH_MAX_NUM_ITERATIONS;

  // prepare cDistParam for cost calculation
  DistParam cDistParam;
  cDistParam.applyWeight = false;
  cDistParam.useMR = false;
#if JVET_Y0089_DMVR_BCW
  cDistParam.useMR |= (pu.cu->BcwIdx != BCW_DEFAULT);
#endif

  Pel* pelBuffer[2] = { nullptr, nullptr };
  pelBuffer[0] = m_filteredBlock[3][REF_PIC_LIST_0][0] + BDMVR_CENTER_POSITION;
  pelBuffer[1] = m_filteredBlock[3][REF_PIC_LIST_1][0] + BDMVR_CENTER_POSITION;

  PelUnitBuf predBuf[2] = { PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_0], BDMVR_BUF_STRIDE, dx, dy)),
    PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_1], BDMVR_BUF_STRIDE, dx, dy)) };

  bool useHadamard = true;          // STAD cost function
  m_pcRdCost->setDistParam(cDistParam, predBuf[0].Y(), predBuf[1].Y(), pu.cu->slice->clpRng(COMPONENT_Y).bd, COMPONENT_Y, useHadamard);

  // prepare buffer for pre-interpolaction 
  const Picture& refPic0 = *pu.cu->slice->getRefPic(REF_PIC_LIST_0, pu.refIdx[REF_PIC_LIST_0])->unscaledPic;
  const Picture& refPic1 = *pu.cu->slice->getRefPic(REF_PIC_LIST_1, pu.refIdx[REF_PIC_LIST_1])->unscaledPic;

  int iWidthExt = dx + (BDMVR_INTME_RANGE << 1);
  int iHeightExt = dy + (BDMVR_INTME_RANGE << 1);
  int iWidthOffset = BDMVR_SIMD_IF_FACTOR - (iWidthExt & (BDMVR_SIMD_IF_FACTOR - 1));
  iWidthOffset &= (BDMVR_SIMD_IF_FACTOR - 1);
  iWidthExt += iWidthOffset; // This ensures that iWidthExt is a factor-of-n number, assuming BDMVR_SIMD_IF_FACTOR is equal to n

  PelUnitBuf predBufExt[2] = { (PelUnitBuf(pu.chromaFormat, PelBuf(m_filteredBlock[3][REF_PIC_LIST_0][0], BDMVR_BUF_STRIDE, iWidthExt, iHeightExt))),
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_filteredBlock[3][REF_PIC_LIST_1][0], BDMVR_BUF_STRIDE, iWidthExt, iHeightExt))) };

  Mv mvTopLeft[2] = { mvInitial[0] - Mv((BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL), (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL)),
    mvInitial[1] - Mv((BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL), (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL)) };
  for (int y = puPos.y, yStart = 0; y < (puPos.y + pu.lumaSize().height); y = y + dy, yStart = yStart + dy)
  {
    for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + dx, xStart = xStart + dx)
    {
      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, dy)));

      minCost = std::numeric_limits<Distortion>::max();

      // Pre-interpolation
      xBDMVRFillBlkPredPelBuffer(subPu, refPic0, mvTopLeft[0], predBufExt[0], pu.cs->slice->clpRng(COMPONENT_Y));
      xBDMVRFillBlkPredPelBuffer(subPu, refPic1, mvTopLeft[1], predBufExt[1], pu.cs->slice->clpRng(COMPONENT_Y));

      if (adaptRange)
      {
        minCost = xBDMVRMvIntPelFullSearch<true, true>(mvOffset, minCost, mvInitial,
          maxSearchRound,
          adaptiveSearchRangeHor, adaptiveSearchRangeVer,
          pu.bmMergeFlag,
          earlyTerminateTh, cDistParam,
          pelBuffer, BDMVR_BUF_STRIDE);
      }
      else
      {
        minCost = xBDMVRMvIntPelFullSearch<false, true>(mvOffset, minCost, mvInitial,
          maxSearchRound,
          adaptiveSearchRangeHor, adaptiveSearchRangeVer,
          pu.bmMergeFlag,
          earlyTerminateTh, cDistParam,
          pelBuffer, BDMVR_BUF_STRIDE);
      }
      if (minCost >= earlyTerminateTh)
      {
        int bestOffsetIdx = (mvOffset.getVer() + BDMVR_INTME_RANGE) * BDMVR_INTME_STRIDE + (mvOffset.getHor() + BDMVR_INTME_RANGE);
        mvOffset <<= MV_FRACTIONAL_BITS_INTERNAL;

        mvFinal[0] = mvInitial[0] + mvOffset;
        mvFinal[1] = mvInitial[1] - mvOffset;
        minCost = m_sadEnlargeArrayBilMrg[bestOffsetIdx];
        Distortion tmpCost = getDecoderSideDerivedMvCost(mvInitial[0], mvFinal[0], BDMVR_INTME_RANGE + 1, DECODER_SIDE_MV_WEIGHT);
        if (minCost >= tmpCost)
        {
          minCost += tmpCost;
          minCost = xBDMVRMvSquareSearch<true>(mvFinal, minCost/*std::numeric_limits<Distortion>::max()*/, subPu, mvInitial, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, false, true);
        }
      }
      else
      {
        mvOffset <<= MV_FRACTIONAL_BITS_INTERNAL;

        mvFinal[0] = mvInitial[0] + mvOffset;
        mvFinal[1] = mvInitial[1] - mvOffset;
      }

      m_bdmvrSubPuMvBuf[REF_PIC_LIST_0][subPuIdx] = mvFinal[0];
      m_bdmvrSubPuMvBuf[REF_PIC_LIST_1][subPuIdx] = mvFinal[1];
      subPuIdx++;
    }
    subPuIdx += dmvrSubPuStrideIncr;
  }
}

#endif

bool InterPrediction::processBDMVR(PredictionUnit& pu)
{
  if( !pu.cs->slice->getSPS()->getUseDMVDMode() || !pu.cs->slice->isInterB() )
  {
    return false;
  }
  CHECK( !pu.mergeFlag, "Merge mode must be used here" );
  CHECK( pu.refIdx[0] < 0 || pu.refIdx[1] < 0, "Bilateral DMVR is performed for bi-prediction" );

  const int lumaArea = pu.lumaSize().area();
  bool subPURefine = true;
  Mv puOrgMv[2] = { pu.mv[0], pu.mv[1] };
  {
    Distortion minCost = std::numeric_limits<Distortion>::max();
    bool       bUseMR = lumaArea > 64;
#if JVET_Y0089_DMVR_BCW
    bUseMR    |= (pu.cu->BcwIdx != BCW_DEFAULT);
#endif
    Mv         mvFinal_PU[2] = { pu.mv[0], pu.mv[1] };
    Mv         mvInitial_PU[2] = { pu.mv[0], pu.mv[1] };

#if JVET_X0049_BDMVR_SW_OPT
#if JVET_X0049_ADAPT_DMVR
    if (pu.bmDir == 1)
    {
      minCost = xBDMVRGetMatchingError(pu, mvInitial_PU, bUseMR, false);
      if (minCost >= lumaArea)
      {
        minCost = xBDMVRMvOneTemplateHPelSquareSearch<1>(mvFinal_PU, minCost, pu, mvInitial_PU, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, bUseMR, false);
      }
    }
    else if (pu.bmDir == 2)
    {
      minCost = xBDMVRGetMatchingError(pu, mvInitial_PU, bUseMR, false);
      if (minCost >= lumaArea)
      {
        minCost = xBDMVRMvOneTemplateHPelSquareSearch<2>(mvFinal_PU, minCost, pu, mvInitial_PU, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, bUseMR, false);
      }
    }
    else
#endif
    {
      minCost = xBDMVRMvSquareSearch<false>( mvFinal_PU, minCost, pu, mvInitial_PU, BDMVR_INTME_SQUARE_SEARCH_MAX_NUM_ITERATIONS, MV_FRACTIONAL_BITS_INTERNAL,     bUseMR, false );
      if (minCost > 0)
      {
        minCost = xBDMVRMvSquareSearch<true>(mvFinal_PU, minCost, pu, mvInitial_PU, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, bUseMR, false);
      }
    }
#else
    minCost = xBDMVRMvSquareSearch( mvFinal_PU, minCost, pu, mvInitial_PU, BDMVR_INTME_SQUARE_SEARCH_MAX_NUM_ITERATIONS, MV_FRACTIONAL_BITS_INTERNAL,     bUseMR, false );
    minCost = xBDMVRMvSquareSearch( mvFinal_PU, minCost, pu, mvInitial_PU, 2, MV_FRACTIONAL_BITS_INTERNAL - 1,     bUseMR, false );
#endif

    subPURefine = minCost >= lumaArea;
    pu.mv[REF_PIC_LIST_0] = mvFinal_PU[0];
    pu.mv[REF_PIC_LIST_1] = mvFinal_PU[1];
  }

#if TM_MRG
  if (pu.tmMergeFlag)
  {
    deriveTMMv(pu);
    if (pu.interDir != 3)
    {
      return false;
    }
  }
#endif

  if (!subPURefine)
  {
    // span motion to subPU
    const int dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
    const int dx = std::min<int>(pu.lumaSize().width, DMVR_SUBCU_WIDTH);
    Position puPos = pu.lumaPos();

    int subPuIdx = 0;
    const int dmvrSubPuStrideIncr = DMVR_SUBPU_STRIDE - std::max(1, (int)(pu.lumaSize().width >> DMVR_SUBCU_WIDTH_LOG2));
    for (int y = puPos.y, yStart = 0; y < (puPos.y + pu.lumaSize().height); y = y + dy, yStart = yStart + dy)
    {
      for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + dx, xStart = xStart + dx)
      {
        m_bdmvrSubPuMvBuf[REF_PIC_LIST_0][subPuIdx] = pu.mv[0];
        m_bdmvrSubPuMvBuf[REF_PIC_LIST_1][subPuIdx] = pu.mv[1];
        subPuIdx++;
      }
      subPuIdx += dmvrSubPuStrideIncr;
    }
    pu.mv[0] = puOrgMv[0];
    pu.mv[1] = puOrgMv[1];
    return true;
  }

  const int dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
  const int dx = std::min<int>(pu.lumaSize().width,  DMVR_SUBCU_WIDTH);
  Position puPos = pu.lumaPos();
  PredictionUnit subPu = pu;

  int subPuIdx = 0;
  const int dmvrSubPuStrideIncr = DMVR_SUBPU_STRIDE - std::max(1, (int)(pu.lumaSize().width >> DMVR_SUBCU_WIDTH_LOG2));
#if JVET_X0049_BDMVR_SW_OPT
  Distortion minCost = std::numeric_limits<Distortion>::max();
  const Mv   mvInitial[2] = { pu.mv[0], pu.mv[1] };
  Mv         mvFinal[2] = { pu.mv[0], pu.mv[1] };
  Mv         mvOffset;

  const Distortion earlyTerminateTh = dx * dy;
  const int adaptiveSearchRangeHor = (dx >> 1) < BDMVR_INTME_RANGE ? (dx >> 1) : BDMVR_INTME_RANGE;
  const int adaptiveSearchRangeVer = (dy >> 1) < BDMVR_INTME_RANGE ? (dy >> 1) : BDMVR_INTME_RANGE;
  const bool adaptRange = (adaptiveSearchRangeHor != BDMVR_INTME_RANGE || adaptiveSearchRangeVer != BDMVR_INTME_RANGE);
#if JVET_X0049_ADAPT_DMVR
  const int maxSearchRound = pu.bmMergeFlag ? BM_MRG_SUB_PU_INT_MAX_SRCH_ROUND : BDMVR_INTME_FULL_SEARCH_MAX_NUM_ITERATIONS;
#else
  const int maxSearchRound = BDMVR_INTME_FULL_SEARCH_MAX_NUM_ITERATIONS;
#endif

  // prepare cDistParam for cost calculation
  DistParam cDistParam;
  cDistParam.applyWeight = false;
  cDistParam.useMR = false;
#if JVET_Y0089_DMVR_BCW
  cDistParam.useMR |= (pu.cu->BcwIdx != BCW_DEFAULT);
#endif

  Pel* pelBuffer[2] = { nullptr, nullptr };
  pelBuffer[0] = m_filteredBlock[3][REF_PIC_LIST_0][0] + BDMVR_CENTER_POSITION;
  pelBuffer[1] = m_filteredBlock[3][REF_PIC_LIST_1][0] + BDMVR_CENTER_POSITION;

  PelUnitBuf predBuf[2] = { PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_0], BDMVR_BUF_STRIDE, dx, dy)),
    PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_1], BDMVR_BUF_STRIDE, dx, dy)) };

  bool useHadamard = true;          // STAD cost function
  m_pcRdCost->setDistParam(cDistParam, predBuf[0].Y(), predBuf[1].Y(), pu.cu->slice->clpRng(COMPONENT_Y).bd, COMPONENT_Y, useHadamard);

  // prepare buffer for pre-interpolaction 
  const Picture& refPic0 = *pu.cu->slice->getRefPic(REF_PIC_LIST_0, pu.refIdx[REF_PIC_LIST_0])->unscaledPic;
  const Picture& refPic1 = *pu.cu->slice->getRefPic(REF_PIC_LIST_1, pu.refIdx[REF_PIC_LIST_1])->unscaledPic;

  int iWidthExt = dx + (BDMVR_INTME_RANGE << 1);
  int iHeightExt = dy + (BDMVR_INTME_RANGE << 1);
  int iWidthOffset = BDMVR_SIMD_IF_FACTOR - (iWidthExt & (BDMVR_SIMD_IF_FACTOR - 1));
  iWidthOffset &= (BDMVR_SIMD_IF_FACTOR - 1);
  iWidthExt += iWidthOffset; // This ensures that iWidthExt is a factor-of-n number, assuming BDMVR_SIMD_IF_FACTOR is equal to n

  PelUnitBuf predBufExt[2] = { (PelUnitBuf(pu.chromaFormat, PelBuf(m_filteredBlock[3][REF_PIC_LIST_0][0], BDMVR_BUF_STRIDE, iWidthExt, iHeightExt))),
    (PelUnitBuf(pu.chromaFormat, PelBuf(m_filteredBlock[3][REF_PIC_LIST_1][0], BDMVR_BUF_STRIDE, iWidthExt, iHeightExt))) };

  Mv mvTopLeft[2] = { mvInitial[0] - Mv((BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL), (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL)),
                      mvInitial[1] - Mv((BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL), (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL)) };
#endif
  for (int y = puPos.y, yStart = 0; y < (puPos.y + pu.lumaSize().height); y = y + dy, yStart = yStart + dy)
  {
    for (int x = puPos.x, xStart = 0; x < (puPos.x + pu.lumaSize().width); x = x + dx, xStart = xStart + dx)
    {
#if JVET_X0049_BDMVR_SW_OPT
      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, dy)));

      minCost = std::numeric_limits<Distortion>::max();

      // Pre-interpolation
      xBDMVRFillBlkPredPelBuffer(subPu, refPic0, mvTopLeft[0], predBufExt[0], pu.cs->slice->clpRng(COMPONENT_Y));
      xBDMVRFillBlkPredPelBuffer(subPu, refPic1, mvTopLeft[1], predBufExt[1], pu.cs->slice->clpRng(COMPONENT_Y));
      
      if (adaptRange)
      {
        minCost = xBDMVRMvIntPelFullSearch<true, true>(mvOffset, minCost, mvInitial,
          maxSearchRound,
          adaptiveSearchRangeHor, adaptiveSearchRangeVer,
#if JVET_X0056_DMVD_EARLY_TERMINATION
          true,
#elif JVET_X0049_ADAPT_DMVR
          pu.bmMergeFlag,
#else
          false,
#endif
          earlyTerminateTh, cDistParam,
          pelBuffer, BDMVR_BUF_STRIDE);
      }
      else
      {
        minCost = xBDMVRMvIntPelFullSearch<false, true>(mvOffset, minCost, mvInitial,
          maxSearchRound,
          adaptiveSearchRangeHor, adaptiveSearchRangeVer,
#if JVET_X0056_DMVD_EARLY_TERMINATION
          true,
#elif JVET_X0049_ADAPT_DMVR
          pu.bmMergeFlag,
#else
          false,
#endif
          earlyTerminateTh, cDistParam,
          pelBuffer, BDMVR_BUF_STRIDE);
      }
      if (minCost >= earlyTerminateTh)
      {
        int bestOffsetIdx = (mvOffset.getVer() + BDMVR_INTME_RANGE) * BDMVR_INTME_STRIDE + (mvOffset.getHor() + BDMVR_INTME_RANGE);
        mvOffset <<= MV_FRACTIONAL_BITS_INTERNAL;

        mvFinal[0] = mvInitial[0] + mvOffset;
        mvFinal[1] = mvInitial[1] - mvOffset;
        minCost = m_sadEnlargeArrayBilMrg[bestOffsetIdx];
        Distortion tmpCost = getDecoderSideDerivedMvCost(mvInitial[0], mvFinal[0], BDMVR_INTME_RANGE + 1, DECODER_SIDE_MV_WEIGHT);
        if (minCost >= tmpCost)
        {
          minCost += tmpCost;
          minCost = xBDMVRMvSquareSearch<true>(mvFinal, minCost/*std::numeric_limits<Distortion>::max()*/, subPu, mvInitial, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, false, true);
        }
      }
      else
      {
        mvOffset <<= MV_FRACTIONAL_BITS_INTERNAL;

        mvFinal[0] = mvInitial[0] + mvOffset;
        mvFinal[1] = mvInitial[1] - mvOffset;
      }
#else
      subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, dy)));
      Distortion minCost = std::numeric_limits<Distortion>::max();
      bool       bUseMR = subPu.lumaSize().area() > 64;
      Mv         mvInitial[2] = { pu.mv[0], pu.mv[1] };
      Mv         mvFinal[2] = { pu.mv[0], pu.mv[1] };
      const int subPuBufOffset = 0;  // will do interpolation inside search

      minCost = xBDMVRMvIntPelFullSearch( mvFinal, minCost, subPu, mvInitial, BDMVR_INTME_SQUARE_SEARCH_MAX_NUM_ITERATIONS, MV_FRACTIONAL_BITS_INTERNAL, bUseMR, subPuBufOffset );
      minCost = (minCost < dx * dy) ? 0 : std::numeric_limits<Distortion>::max();
      minCost = xBDMVRMvSquareSearch( mvFinal, minCost, subPu, mvInitial, 2, MV_FRACTIONAL_BITS_INTERNAL - 1, false,  true);
#endif

      m_bdmvrSubPuMvBuf[REF_PIC_LIST_0][subPuIdx] = mvFinal[0];
      m_bdmvrSubPuMvBuf[REF_PIC_LIST_1][subPuIdx] = mvFinal[1];
      subPuIdx++;
    }
    subPuIdx += dmvrSubPuStrideIncr;
  }

  pu.mv[0] = puOrgMv[0];
  pu.mv[1] = puOrgMv[1];
  return true;
}

void InterPrediction::xBDMVRFillBlkPredPelBuffer(const PredictionUnit& pu, const Picture& refPic, const Mv &_mv, PelUnitBuf &dstBuf, const ClpRng& clpRng)
{
  const ComponentID compID = COMPONENT_Y;
  const CPelBuf     refBuf = refPic.getRecoBuf(refPic.blocks[compID]);

  const int lumaShift = 2 + MV_FRACTIONAL_BITS_DIFF;
  const int horShift = (lumaShift + ::getComponentScaleX(compID, pu.chromaFormat));
  const int verShift = (lumaShift + ::getComponentScaleY(compID, pu.chromaFormat));

  Mv mv(_mv);
  clipMv(mv, pu.lumaPos(), pu.lumaSize(), *pu.cu->cs->sps, *pu.cu->cs->pps);
  const int xInt  = mv.getHor() >> horShift;
  const int yInt  = mv.getVer() >> verShift;
  const int xFrac = mv.getHor() & ((1 << horShift) - 1);
  const int yFrac = mv.getVer() & ((1 << verShift) - 1);

  const Pel* ref       = refBuf.bufAt(pu.blocks[compID].pos().offset(xInt, yInt));
        Pel* dst       = dstBuf.bufs[compID].buf;
        int  refStride = refBuf.stride;
        int  dstStride = dstBuf.bufs[compID].stride;
        int  bw        = (int)dstBuf.bufs[compID].width;
        int  bh        = (int)dstBuf.bufs[compID].height;

  const int  nFilterIdx   = 0;
  const bool useAltHpelIf = pu.cu->imv == IMV_HPEL;
  const bool biMCForDMVR  = true;

  if ( yFrac == 0 )
  {
    m_if.filterHor( compID, (Pel*) ref, refStride, dst, dstStride, bw, bh, xFrac, false/*rndRes=!bi*/,
                    pu.chromaFormat, pu.cu->slice->clpRng(compID), biMCForDMVR, biMCForDMVR, useAltHpelIf );
  }
  else if ( xFrac == 0 )
  {
    m_if.filterVer( compID, (Pel*) ref, refStride, dst, dstStride, bw, bh, yFrac, true, false/*rndRes=!bi*/,
        pu.chromaFormat, pu.cu->slice->clpRng(compID), biMCForDMVR, biMCForDMVR, useAltHpelIf );
  }
  else
  {
#if IF_12TAP
    int vFilterSize = isLuma(compID) ? (nFilterIdx == 1 ? NTAPS_BILINEAR : NTAPS_LUMA(0)) : NTAPS_CHROMA;
#else
    int vFilterSize = isLuma(compID) ? (nFilterIdx == 1 ? NTAPS_BILINEAR : NTAPS_LUMA) : NTAPS_CHROMA;
#endif
    if (biMCForDMVR)
    {
      vFilterSize = NTAPS_BILINEAR;
    }
    PelBuf tmpBuf = PelBuf(m_filteredBlockTmp[0][compID], Size(bw + 2 * BDMVR_INTME_RANGE, bh + 2 * BDMVR_INTME_RANGE));

    m_if.filterHor( compID, (Pel*)ref - ((vFilterSize>>1) -1)*refStride, refStride, tmpBuf.buf, tmpBuf.stride, bw, bh+vFilterSize-1, xFrac,
                    false, pu.chromaFormat, pu.cu->slice->clpRng(compID), biMCForDMVR, biMCForDMVR, useAltHpelIf );
    JVET_J0090_SET_CACHE_ENABLE( false );
    m_if.filterVer( compID, tmpBuf.buf + ((vFilterSize>>1) -1)*tmpBuf.stride, tmpBuf.stride, dst, dstStride, bw, bh, yFrac,
                    false, false/*rndRes=!bi*/, pu.chromaFormat, pu.cu->slice->clpRng(compID), biMCForDMVR, biMCForDMVR, useAltHpelIf );
    JVET_J0090_SET_CACHE_ENABLE( true );
  }
}

#if JVET_X0049_ADAPT_DMVR
template <uint8_t dir>
#endif
void InterPrediction::xBDMVRPreInterpolation(const PredictionUnit& pu, const Mv (&mvCenter)[2], bool doPreInterpolationFP, bool doPreInterpolationHP)
{
  if (doPreInterpolationFP)
  {
    for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
#if JVET_X0049_ADAPT_DMVR
      if (!(dir & (1 << refList)))
      {
        continue;
      }
#endif
      const Picture& refPic  = *pu.cu->slice->getRefPic((RefPicList)refList, pu.refIdx[refList])->unscaledPic;

      int dstStride    = MAX_CU_SIZE + ( BDMVR_INTME_RANGE << 1 ) + ( BDMVR_SIMD_IF_FACTOR - 2 );
      int iWidthExt    = (int)pu.lwidth () + ( BDMVR_INTME_RANGE << 1 );
      int iHeightExt   = (int)pu.lheight() + ( BDMVR_INTME_RANGE << 1 );
      int iWidthOffset = BDMVR_SIMD_IF_FACTOR - ( iWidthExt & ( BDMVR_SIMD_IF_FACTOR - 1 ) );
      iWidthOffset    &= ( BDMVR_SIMD_IF_FACTOR - 1 );
      iWidthExt       += iWidthOffset; // This ensures that iWidthExt is a factor-of-n number, assuming BDMVR_SIMD_IF_FACTOR is equal to n

      Mv mv = mvCenter[refList] - Mv((BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL), (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL));
      PelUnitBuf predBuf = ( PelUnitBuf(pu.chromaFormat, PelBuf(m_filteredBlock[3][refList][0], dstStride, iWidthExt, iHeightExt ) ) );
      xBDMVRFillBlkPredPelBuffer(pu, refPic, mv, predBuf, pu.cs->slice->clpRng(COMPONENT_Y));
    }
  }

  if (doPreInterpolationHP)
  {
    const int offset = 0 - ( 1 << ( MV_FRACTIONAL_BITS_INTERNAL - 1 ) );
    const Mv  cPhaseOffset[3] = { Mv( offset , 0 ), Mv( offset, offset ), Mv( 0 , offset ) };

    for (int refList = 0; refList < NUM_REF_PIC_LIST_01 ; refList++)
    {
#if JVET_X0049_ADAPT_DMVR
      if (!(dir & (1 << refList)))
      {
        continue;
      }
#endif
      const Picture& refPic = *pu.cu->slice->getRefPic((RefPicList)refList, pu.refIdx[refList])->unscaledPic;

      for (int phaseIdx = 0 ; phaseIdx < 3 ; phaseIdx++)
      {
        int iRefStride   = MAX_CU_SIZE  + ( BDMVR_INTME_RANGE << 1 ) + ( BDMVR_SIMD_IF_FACTOR - 2 );
        int iWidthExt    = (int)pu.lwidth () + 1 - (   phaseIdx  >> 1);
        int iHeightExt   = (int)pu.lheight() + 1 - ((2-phaseIdx) >> 1);
        int iWidthOffset = BDMVR_SIMD_IF_FACTOR - ( iWidthExt & ( BDMVR_SIMD_IF_FACTOR - 1 ) );
        iWidthOffset    &= ( BDMVR_SIMD_IF_FACTOR - 1 );
        iWidthExt       += iWidthOffset; // This ensures that iWidthExt is a factor-of-n number, assuming BDMVR_SIMD_IF_FACTOR is equal to n

        Mv mv = mvCenter[refList] + cPhaseOffset[phaseIdx];
        PelUnitBuf predBuf = PelUnitBuf( pu.chromaFormat, PelBuf( m_filteredBlock[phaseIdx][refList][0], iRefStride, iWidthExt, iHeightExt ) );
        xBDMVRFillBlkPredPelBuffer( pu, refPic, mv, predBuf, pu.cs->slice->clpRng(COMPONENT_Y) );
      }
    }
  }
}

#if JVET_X0049_BDMVR_SW_OPT
template <bool adaptRange, bool useHadamard>
Distortion InterPrediction::xBDMVRMvIntPelFullSearch(Mv&mvOffset, Distortion curBestCost, const Mv(&initialMv)[2], const int32_t maxSearchRounds, const int maxHorOffset, const int maxVerOffset, const bool earlySkip, const Distortion earlyTerminateTh, DistParam &cDistParam, Pel* pelBuffer[2], const int stride)
{
  // check initial cost  
  mvOffset.setZero();
  cDistParam.org.buf = pelBuffer[0];
  cDistParam.cur.buf = pelBuffer[1];

#if FULL_NBIT
  if (useHadamard)
  {
    curBestCost = cDistParam.distFunc(cDistParam) >> 1;  // magic shift, benefit for early terminate
  }
  else
  {
    int32_t precisionAdj = cDistParam.bitDepth > 8 ? cDistParam.bitDepth - 8 : 0;
    curBestCost = cDistParam.distFunc(cDistParam) >> precisionAdj;
  }
#else
  curBestCost = cDistParam.distFunc(cDistParam);
#endif

  m_sadEnlargeArrayBilMrg[BDMVR_INTME_CENTER] = curBestCost;
  curBestCost = curBestCost - (curBestCost >> 2);  // cost tuning

  if (curBestCost < earlyTerminateTh)
  {
    return curBestCost;
  }

  Distortion tmCost = MAX_UINT64;
  Distortion prevMinCost = MAX_UINT64;

  for (int searchPrio = 1; searchPrio < maxSearchRounds; searchPrio++)
  {
    prevMinCost = curBestCost;
    for (int currIdx = 0; currIdx < m_searchEnlargeOffsetNum[searchPrio]; currIdx++)
    {
      tmCost = 0;
      int horOffset = m_searchEnlargeOffsetBilMrg[searchPrio][currIdx].getHor();
      int verOffset = m_searchEnlargeOffsetBilMrg[searchPrio][currIdx].getVer();
      int searchOffsetIdx = m_searchEnlargeOffsetToIdx[searchPrio][currIdx];

      if (adaptRange)
      {
        if (abs(horOffset) > maxHorOffset || abs(verOffset) > maxVerOffset)
        {
          continue;
        }
      }
        
      int bufOffset = verOffset * stride + horOffset;
      cDistParam.org.buf = pelBuffer[0] + bufOffset;
      cDistParam.cur.buf = pelBuffer[1] - bufOffset;

#if FULL_NBIT
      if (useHadamard)
      {
        m_sadEnlargeArrayBilMrg[searchOffsetIdx] = cDistParam.distFunc(cDistParam) >> 1;  // magic shift, benefit for early terminate
      }
      else
      {
        int32_t precisionAdj = cDistParam.bitDepth > 8 ? cDistParam.bitDepth - 8 : 0;
        m_sadEnlargeArrayBilMrg[searchOffsetIdx] = cDistParam.distFunc(cDistParam) >> precisionAdj;
      }
#else
      m_sadEnlargeArrayBilMrg[searchOffsetIdx] = cDistParam.distFunc(cDistParam);
#endif

      tmCost += m_sadEnlargeArrayBilMrg[searchOffsetIdx];
      tmCost += (m_sadEnlargeArrayBilMrg[searchOffsetIdx] >> m_costShiftBilMrg1[searchOffsetIdx]);
      tmCost += (m_sadEnlargeArrayBilMrg[searchOffsetIdx] >> m_costShiftBilMrg2[searchOffsetIdx]);

      if (tmCost < curBestCost)
      {
        mvOffset = Mv(horOffset, verOffset);
        curBestCost = tmCost;
      }
    }

    if (curBestCost < earlyTerminateTh)
    {
      break;
    }
    if (earlySkip && searchPrio > 1 && prevMinCost - curBestCost < earlyTerminateTh)
    {
      break;
    }
  }
  return curBestCost;
}
#else
Distortion InterPrediction::xBDMVRMvIntPelFullSearch(Mv(&curBestMv)[2], Distortion curBestCost, PredictionUnit& pu, const Mv(&initialMv)[2], int32_t maxSearchRounds, int32_t searchStepShift, bool useMR, const int subPuBufOffset)
{
  bool doPreInterpolation = true;
  bool useHadamard = true;          // STAD cost function
  useMR = false;                    // STAD cost function
  const int adaptiveSearchRangeHor = (pu.lwidth() >> 1) < BDMVR_INTME_RANGE ? (pu.lwidth() >> 1) : BDMVR_INTME_RANGE;
  const int adaptiveSearchRangeVer = (pu.lheight() >> 1) < BDMVR_INTME_RANGE ? (pu.lheight() >> 1) : BDMVR_INTME_RANGE;

  // Calculate TM cost of initial MVs, if it is not set
  if (curBestCost == std::numeric_limits<Distortion>::max())
  {
    curBestCost = xBDMVRGetMatchingError( pu, curBestMv, subPuBufOffset, useHadamard, useMR,
                                          doPreInterpolation, searchStepShift, curBestMv, initialMv, -1 );
  }

  for (int i = 0; i < BDMVR_INTME_AREA; i++)
  {
    m_sadEnlargeArrayBilMrg[i] = MAX_UINT64;
  }
  m_sadEnlargeArrayBilMrg[BDMVR_INTME_CENTER] = curBestCost;
  curBestCost = curBestCost - (curBestCost >> 2);  // cost tuning
  const Distortion earlyTerminateTh = pu.lumaSize().area();
  Distortion tmCost = MAX_UINT64;
#if JVET_X0056_DMVD_EARLY_TERMINATION
  Distortion prevMinCost = MAX_UINT64;
#endif

  for( int searchPrio = 0 ; searchPrio < BDMVR_INTME_FULL_SEARCH_MAX_NUM_ITERATIONS; searchPrio++ )
  {
    if( curBestCost < earlyTerminateTh )
    {
      break;
    }

#if JVET_X0056_DMVD_EARLY_TERMINATION
    prevMinCost = curBestCost;
#endif
    for (int searchOffsetIdx = 0; searchOffsetIdx < BDMVR_INTME_AREA; searchOffsetIdx++)
    {
      tmCost = 0;
      if( m_searchPriorityBilMrg[searchOffsetIdx] != searchPrio )
      {
        continue;
      }
      // adaptive search area base on block dimension
      if( m_searchEnlargeOffsetBilMrg[searchOffsetIdx].getAbsVer() > adaptiveSearchRangeVer )
      {
        continue;
      }
      if( m_searchEnlargeOffsetBilMrg[searchOffsetIdx].getAbsHor() > adaptiveSearchRangeHor )
      {
        continue;
      }

      Mv mvOffset(m_searchEnlargeOffsetBilMrg[searchOffsetIdx].getHor() << searchStepShift, m_searchEnlargeOffsetBilMrg[searchOffsetIdx].getVer() << searchStepShift);
      Mv mvCand[2] = {initialMv[0] + mvOffset, initialMv[1] - mvOffset};

      if ( m_sadEnlargeArrayBilMrg[searchOffsetIdx] == MAX_UINT64 )
      {
        m_sadEnlargeArrayBilMrg[searchOffsetIdx] = xBDMVRGetMatchingError( pu, mvCand, subPuBufOffset, useHadamard, useMR,
                                                                             doPreInterpolation, searchStepShift, initialMv, initialMv, -1 );
      }
      tmCost += m_sadEnlargeArrayBilMrg[searchOffsetIdx];
      tmCost += (m_sadEnlargeArrayBilMrg[searchOffsetIdx] >> m_costShiftBilMrg1[searchOffsetIdx]);
      tmCost += (m_sadEnlargeArrayBilMrg[searchOffsetIdx] >> m_costShiftBilMrg2[searchOffsetIdx]);

      if( tmCost < curBestCost )
      {
        curBestCost  = tmCost;
        curBestMv[0] = mvCand[0];
        curBestMv[1] = mvCand[1];
      }
    }
#if JVET_X0056_DMVD_EARLY_TERMINATION
    if (searchPrio > 1 && prevMinCost - curBestCost < earlyTerminateTh)
    {
      break;
    }
#endif
  }

  return curBestCost;
}
#endif

#if JVET_X0049_BDMVR_SW_OPT
template<bool hPel>
#endif
Distortion InterPrediction::xBDMVRMvSquareSearch(Mv (&curBestMv)[2], Distortion curBestCost, PredictionUnit& pu, const Mv (&initialMv)[2], int32_t maxSearchRounds, int32_t searchStepShift, bool useMR, bool useHadmard)
{
#if !JVET_X0049_BDMVR_SW_OPT
  if (curBestCost == 0)
  {
    return 0;
  }
#endif
  static const Mv   cSearchOffset[8] = { Mv( -1 , 1 ) , Mv( 0 , 1 ) , Mv(  1 ,  1 ) , Mv(  1 ,  0 ) , Mv(  1 , -1 ) , Mv(  0 , -1 ) , Mv( -1 , -1 ) , Mv( -1 , 0 )  };
               int  nDirectStart     = 0;
               int  nDirectEnd       = 7;
         const int  nDirectRounding  = 8;
         const int  nDirectMask      = 0x07;
               bool doPreInterpolation = searchStepShift == MV_FRACTIONAL_BITS_INTERNAL;

  // Calculate TM cost of initial MVs, if it is not set
  if (curBestCost == std::numeric_limits<Distortion>::max())
  {
    CHECK(searchStepShift < MV_FRACTIONAL_BITS_INTERNAL - 1, "this is not possible");
#if JVET_X0049_BDMVR_SW_OPT
    if (hPel)
    {
      doPreInterpolation = true;
      Distortion tmCost = getDecoderSideDerivedMvCost(initialMv[0], curBestMv[0], BDMVR_INTME_RANGE + (MV_FRACTIONAL_BITS_INTERNAL - searchStepShift), DECODER_SIDE_MV_WEIGHT);
      curBestCost = xBDMVRGetMatchingError(pu, curBestMv, useMR, useHadmard);
#else
    if (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1)
    {
      doPreInterpolation = true;
      Distortion tmCost = getDecoderSideDerivedMvCost(initialMv[0], curBestMv[0], BDMVR_INTME_RANGE + (MV_FRACTIONAL_BITS_INTERNAL - searchStepShift), DECODER_SIDE_MV_WEIGHT);
      curBestCost = xBDMVRGetMatchingError( pu, curBestMv, 0/*subPuOffset*/, useHadmard, useMR,
                                            doPreInterpolation, MV_FRACTIONAL_BITS_INTERNAL, curBestMv, curBestMv, -1 );
#endif
      if( curBestCost < tmCost )
      {
        return curBestCost;
      }

      curBestCost += tmCost;
    }
    else
    {
#if JVET_X0049_ADAPT_DMVR
      curBestCost = xBDMVRGetMatchingError<3>(pu, curBestMv, 0/*subPuOffset*/, useHadmard, useMR,
        doPreInterpolation, searchStepShift, curBestMv, initialMv, -1);
#else
      curBestCost = xBDMVRGetMatchingError( pu, curBestMv, 0/*subPuOffset*/, useHadmard, useMR,
                                            doPreInterpolation, searchStepShift, curBestMv, initialMv, -1 );
#endif
    }
  }

  Distortion localCostArray[9] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(),
                                   std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(),
                                   std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(), curBestCost };

  // Iterative search process
  for( uint32_t uiRound = 0 ; uiRound < maxSearchRounds ; uiRound++ )
  {
    int nBestDirect    = -1;
    Mv  mvCurCenter[2] = {curBestMv[0], curBestMv[1]};
    doPreInterpolation |= (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1);

    for( int nIdx = nDirectStart ; nIdx <= nDirectEnd ; nIdx++ )
    {
      int nDirect = ( nIdx + nDirectRounding ) & nDirectMask;

      Mv mvOffset(cSearchOffset[nDirect].getHor() << searchStepShift, cSearchOffset[nDirect].getVer() << searchStepShift);

#if JVET_X0049_BDMVR_SW_OPT
      if(hPel && uiRound > 0)
#else
      if (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1 && uiRound > 0)
#endif
      {
        if( ( nDirect % 2 ) == 0 )
        {
          continue;
        }
      }
      Mv mvCand[2] = {mvCurCenter[0] + mvOffset, mvCurCenter[1] - mvOffset};
#if JVET_X0049_BDMVR_SW_OPT
      if(!hPel)
#else
      if (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL)
#endif
      {
        int currentIdx = BDMVR_INTME_CENTER + ((mvCand[0] -initialMv[0]).hor >> searchStepShift) +
                                              ((mvCand[0] -initialMv[0]).ver >> searchStepShift) * BDMVR_INTME_STRIDE;
        if( currentIdx < 0 || currentIdx >= BDMVR_INTME_AREA )
        {
          continue;
        }
      }

      Distortion tmCost = getDecoderSideDerivedMvCost(initialMv[0], mvCand[0], BDMVR_INTME_RANGE + (MV_FRACTIONAL_BITS_INTERNAL - searchStepShift), DECODER_SIDE_MV_WEIGHT);

      if (tmCost > curBestCost)
      {
        localCostArray[nDirect] = 2 * tmCost;
        continue;
      }

#if JVET_X0049_ADAPT_DMVR
      tmCost += xBDMVRGetMatchingError<3>(pu, mvCand, 0/*subPuOffset*/, useHadmard, useMR,
        doPreInterpolation, searchStepShift, mvCurCenter, initialMv, nDirect);
#else
      tmCost += xBDMVRGetMatchingError( pu, mvCand, 0/*subPuOffset*/, useHadmard, useMR,
                                        doPreInterpolation, searchStepShift, mvCurCenter, initialMv, nDirect );
#endif
      localCostArray[nDirect] = tmCost;
#if JVET_X0049_BDMVR_SW_OPT
      if(hPel && uiRound > 0)
#else
      if (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1 && uiRound > 0)
#endif
      {
        continue;
      }

      if( tmCost < curBestCost )
      {
        nBestDirect  = nDirect;
        curBestCost  = tmCost;
        curBestMv[0] = mvCand[0];
        curBestMv[1] = mvCand[1];
      }
    }

    if( nBestDirect == -1 )
    {
      break;
    }

    int nStep    = 2 - ( nBestDirect & 0x01 );
    nDirectStart = nBestDirect - nStep;
    nDirectEnd   = nBestDirect + nStep;

    if ((uiRound + 1) < maxSearchRounds)
    {
      xBDMVRUpdateSquareSearchCostLog(localCostArray, nBestDirect);
    }
  }

#if JVET_X0049_BDMVR_SW_OPT
  if(!hPel)
#else
  if (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL)
#endif
  {
    return curBestCost;
  }

  // Model-based fractional MVD optimization
  Mv mvDiff = curBestMv[0] - initialMv[0];
  if (localCostArray[8] > 0 && localCostArray[8] == curBestCost && mvDiff.getAbsHor() != (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL) && mvDiff.getAbsVer() != (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL))
  {
    uint64_t sadbuffer[5];
    sadbuffer[0] = (uint64_t)localCostArray[8]; // center
    sadbuffer[1] = (uint64_t)localCostArray[7]; // left
    sadbuffer[2] = (uint64_t)localCostArray[5]; // above
    sadbuffer[3] = (uint64_t)localCostArray[3]; // right
    sadbuffer[4] = (uint64_t)localCostArray[1]; // bottom

    int32_t tempDeltaMv[2] = {0, 0};
    xSubPelErrorSrfc(sadbuffer, tempDeltaMv);

    curBestMv[0] += Mv(tempDeltaMv[0], tempDeltaMv[1]);
    curBestMv[1] -= Mv(tempDeltaMv[0], tempDeltaMv[1]);
  }

  return curBestCost;
}

#if JVET_X0049_ADAPT_DMVR
template <uint8_t dir>
Distortion InterPrediction::xBDMVRMvOneTemplateHPelSquareSearch(Mv(&curBestMv)[2], Distortion curBestCost, PredictionUnit& pu,
  const Mv(&initialMv)[2], int32_t maxSearchRounds, int32_t searchStepShift,
  bool useMR, bool useHadmard)
{
  if (curBestCost == 0)
  {
    return 0;
  }

  static const Mv   cSearchOffset[8] = { Mv(-1 , 1) , Mv(0 , 1) , Mv(1 ,  1) , Mv(1 ,  0) , Mv(1 , -1) , Mv(0 , -1) , Mv(-1 , -1) , Mv(-1 , 0) };
  int  nDirectStart = 0;
  int  nDirectEnd = 7;
  const int  nDirectRounding = 8;
  const int  nDirectMask = 0x07;
  bool doPreInterpolation = searchStepShift == MV_FRACTIONAL_BITS_INTERNAL;
  const int curRefList = (dir >> 1);
  const int templateRefList = 1 - curRefList;
  // Calculate TM cost of initial MVs, if it is not set
  if (curBestCost == std::numeric_limits<Distortion>::max())
  {
    CHECK(searchStepShift < MV_FRACTIONAL_BITS_INTERNAL - 1, "this is not possible");
    Distortion tmCost = getDecoderSideDerivedMvCost(initialMv[curRefList], curBestMv[curRefList], BDMVR_INTME_RANGE + (MV_FRACTIONAL_BITS_INTERNAL - searchStepShift), DECODER_SIDE_MV_WEIGHT);
    curBestCost = xBDMVRGetMatchingError(pu, curBestMv, useMR, useHadmard);
    if (curBestCost < tmCost)
    {
      return curBestCost;
    }

    curBestCost += tmCost;
  }

  Distortion localCostArray[9] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(),
    std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(),
    std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max(), curBestCost };

  // Iterative search process
  for (uint32_t uiRound = 0; uiRound < maxSearchRounds; uiRound++)
  {
    int nBestDirect = -1;
    Mv  mvCurCenter[2] = { curBestMv[0], curBestMv[1] };
    doPreInterpolation |= (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1);

    for (int nIdx = nDirectStart; nIdx <= nDirectEnd; nIdx++)
    {
      int nDirect = (nIdx + nDirectRounding) & nDirectMask;

      Mv mvOffset(cSearchOffset[nDirect].getHor() << searchStepShift, cSearchOffset[nDirect].getVer() << searchStepShift);

      if (uiRound > 0)
      {
        if ((nDirect % 2) == 0)
        {
          continue;
        }
      }
      Mv mvCand[2] = { mvCurCenter[0] + mvOffset, mvCurCenter[1] - mvOffset };
      mvCand[templateRefList] = initialMv[templateRefList];
      Distortion tmCost = getDecoderSideDerivedMvCost(initialMv[curRefList], mvCand[curRefList], BDMVR_INTME_RANGE + (MV_FRACTIONAL_BITS_INTERNAL - searchStepShift), DECODER_SIDE_MV_WEIGHT);

      if (tmCost > curBestCost)
      {
        localCostArray[nDirect] = 2 * tmCost;
        continue;
      }

      tmCost += xBDMVRGetMatchingError<dir>(pu, mvCand, 0/*subPuOffset*/, useHadmard, useMR,
        doPreInterpolation, searchStepShift, mvCurCenter, initialMv, nDirect);
      localCostArray[nDirect] = tmCost;
      if (uiRound > 0)
      {
        continue;
      }

      if (tmCost < curBestCost)
      {
        nBestDirect = nDirect;
        curBestCost = tmCost;
        curBestMv[0] = mvCand[0];
        curBestMv[1] = mvCand[1];
      }
    }

    if (nBestDirect == -1)
    {
      break;
    }

    int nStep = 2 - (nBestDirect & 0x01);
    nDirectStart = nBestDirect - nStep;
    nDirectEnd = nBestDirect + nStep;

    if ((uiRound + 1) < maxSearchRounds)
    {
      xBDMVRUpdateSquareSearchCostLog(localCostArray, nBestDirect);
    }
  }

  CHECK(curBestMv[templateRefList] != initialMv[templateRefList], "this is not possible");
  // Model-based fractional MVD optimization
  Mv mvDiff = curBestMv[curRefList] - initialMv[curRefList];
  if (localCostArray[8] > 0 && localCostArray[8] == curBestCost && mvDiff.getAbsHor() != (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL) && mvDiff.getAbsVer() != (BDMVR_INTME_RANGE << MV_FRACTIONAL_BITS_INTERNAL))
  {
    uint64_t sadbuffer[5];
    sadbuffer[0] = (uint64_t)localCostArray[8]; // center
    sadbuffer[1] = (uint64_t)localCostArray[7]; // left
    sadbuffer[2] = (uint64_t)localCostArray[5]; // above
    sadbuffer[3] = (uint64_t)localCostArray[3]; // right
    sadbuffer[4] = (uint64_t)localCostArray[1]; // bottom

    int32_t tempDeltaMv[2] = { 0, 0 };
    xSubPelErrorSrfc(sadbuffer, tempDeltaMv);

    if (dir == 1)
      curBestMv[0] += Mv(tempDeltaMv[0], tempDeltaMv[1]);
    else
      curBestMv[1] -= Mv(tempDeltaMv[0], tempDeltaMv[1]);
  }

  return curBestCost;
}
#endif

#if JVET_X0049_BDMVR_SW_OPT
Distortion InterPrediction::xBDMVRGetMatchingError(const PredictionUnit& pu, const Mv(&mv)[2], bool useMR, bool useHadmard)
#else
Distortion InterPrediction::xBDMVRGetMatchingError(const PredictionUnit& pu, const Mv(&mv)[2], bool useMR)
#endif
{
  // Fill L0'a and L1's prediction blocks
#if JVET_X0049_BDMVR_SW_OPT
        Pel*     pelBuffer[2] = { m_filteredBlock[3][REF_PIC_LIST_0][0] + BDMVR_CENTER_POSITION, m_filteredBlock[3][REF_PIC_LIST_1][0] + BDMVR_CENTER_POSITION };
        const SizeType stride = BDMVR_BUF_STRIDE;
#else
        Pel*     pelBuffer[2] = { m_filteredBlock[3][REF_PIC_LIST_0][0], m_filteredBlock[3][REF_PIC_LIST_1][0] };
  const SizeType stride       = pu.lwidth();
#endif
  PelUnitBuf     predBuf[2]   = { PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_0], stride, pu.lwidth(), pu.lheight())),
                                  PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_1], stride, pu.lwidth(), pu.lheight())) };

  for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
  {
#if JVET_X0083_BM_AMVP_MERGE_MODE
    if (pu.amvpMergeModeFlag[1 - refList])
    {
      continue;
    }
#endif
    const Picture&   refPic  = *pu.cu->slice->getRefPic((RefPicList)refList, pu.refIdx[refList])->unscaledPic;
    xBDMVRFillBlkPredPelBuffer( pu, refPic, mv[refList] , predBuf[refList], pu.cs->slice->clpRng(COMPONENT_Y) );
  }

  // Compute distortion between L0'a and L1's prediction blocks
  DistParam cDistParam;
  cDistParam.applyWeight = false;
  cDistParam.useMR       = useMR;

#if JVET_X0049_BDMVR_SW_OPT
  m_pcRdCost->setDistParam(cDistParam, predBuf[0].Y(), predBuf[1].Y(), pu.cu->slice->clpRng(COMPONENT_Y).bd, COMPONENT_Y, useHadmard);
#if FULL_NBIT
  if (useHadmard)
  {
    return cDistParam.distFunc(cDistParam) >> 1;  // magic shift, benefit for early terminate
  }
  else
  {
    int32_t precisionAdj = cDistParam.bitDepth > 8 ? cDistParam.bitDepth - 8 : 0;
    return cDistParam.distFunc(cDistParam) >> precisionAdj;
  }
#else
  return cDistParam.distFunc(cDistParam);
#endif
#else
  m_pcRdCost->setDistParam( cDistParam, predBuf[0].Y(), predBuf[1].Y(), pu.cu->slice->clpRng(COMPONENT_Y).bd, COMPONENT_Y, false );
#if FULL_NBIT
  int32_t precisionAdj = cDistParam.bitDepth > 8 ? cDistParam.bitDepth - 8 : 0;
  return cDistParam.distFunc( cDistParam ) >> precisionAdj;
#else
  return cDistParam.distFunc( cDistParam );
#endif
#endif
}

#if MULTI_PASS_DMVR
#if JVET_X0049_ADAPT_DMVR
template <uint8_t dir>
#endif
Distortion InterPrediction::xBDMVRGetMatchingError(const PredictionUnit& pu, const Mv (&mv)[2], const int subPuBufOffset, bool useHadmard, bool useMR
                                                 , bool& doPreInterpolation, int32_t searchStepShift, const Mv (&mvCenter)[2]
                                                 , const Mv (&mvInitial)[2]
                                                 , int nDirect
)
{
  // Pre-interpolation
  if (doPreInterpolation)
  {
#if JVET_X0049_ADAPT_DMVR
    xBDMVRPreInterpolation<dir>(pu, mvCenter, searchStepShift == MV_FRACTIONAL_BITS_INTERNAL, searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1);
#else
    xBDMVRPreInterpolation( pu, mvCenter, searchStepShift == MV_FRACTIONAL_BITS_INTERNAL, searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1 );
#endif
    doPreInterpolation = false;
  }

  // Locate L0'a and L1's prediction blocks in pre-interpolation buffer
#if JVET_X0049_BDMVR_SW_OPT
  const int32_t stride = BDMVR_BUF_STRIDE;
#else
  const int32_t stride = MAX_CU_SIZE + ( BDMVR_INTME_RANGE << 1 ) + ( BDMVR_SIMD_IF_FACTOR - 2 );
#endif
  Pel* pelBuffer[2] = { nullptr, nullptr };

  if (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL)
  {
    Mv mvDiff[2] = { mv[0] - mvInitial[0], mv[1] - mvInitial[1] };
    mvDiff[0]  >>= MV_FRACTIONAL_BITS_INTERNAL;
    mvDiff[1]  >>= MV_FRACTIONAL_BITS_INTERNAL;

#if JVET_X0049_ADAPT_DMVR
    if (dir == 1)
    {
      // fix template at refList1 
      CHECK(subPuBufOffset != 0, "this is not possible");
      pelBuffer[0] = m_filteredBlock[3][REF_PIC_LIST_0][0] + subPuBufOffset + BDMVR_CENTER_POSITION + mvDiff[0].getVer() * stride + mvDiff[0].getHor();
      pelBuffer[1] = m_filteredBlock[3][REF_PIC_LIST_1][0] + BDMVR_CENTER_POSITION;
    }
    else if (dir == 2)
    {
      CHECK(subPuBufOffset != 0, "this is not possible");
      pelBuffer[0] = m_filteredBlock[3][REF_PIC_LIST_0][0] + BDMVR_CENTER_POSITION;
      pelBuffer[1] = m_filteredBlock[3][REF_PIC_LIST_1][0] + subPuBufOffset + BDMVR_CENTER_POSITION + mvDiff[1].getVer() * stride + mvDiff[1].getHor();
    }
    else
    {
      pelBuffer[0] = m_filteredBlock[3][REF_PIC_LIST_0][0] + subPuBufOffset + BDMVR_CENTER_POSITION + mvDiff[0].getVer() * stride + mvDiff[0].getHor();
      pelBuffer[1] = m_filteredBlock[3][REF_PIC_LIST_1][0] + subPuBufOffset + BDMVR_CENTER_POSITION + mvDiff[1].getVer() * stride + mvDiff[1].getHor();
    }
#else
    pelBuffer[0] = m_filteredBlock[3][REF_PIC_LIST_0][0] + subPuBufOffset + ( BDMVR_INTME_RANGE + mvDiff[0].getVer() ) * stride + BDMVR_INTME_RANGE + mvDiff[0].getHor();
    pelBuffer[1] = m_filteredBlock[3][REF_PIC_LIST_1][0] + subPuBufOffset + ( BDMVR_INTME_RANGE + mvDiff[1].getVer() ) * stride + BDMVR_INTME_RANGE + mvDiff[1].getHor();
#endif
  }
  else if (searchStepShift == MV_FRACTIONAL_BITS_INTERNAL - 1)
  {
           const int32_t  cFracBufOffset[8] = { stride, stride, stride + 1, 1, 1, 0, 0, 0 };
    static const uint32_t phaseIdxList[4]   = { 1, 2, 1, 0 };

    int phaseIdx = phaseIdxList[ nDirect & 0x3 ];
#if JVET_X0049_ADAPT_DMVR 
    if (dir == 3)
    {
      pelBuffer[0] = m_filteredBlock[phaseIdx][REF_PIC_LIST_0][0] + cFracBufOffset[nDirect];
      pelBuffer[1] = m_filteredBlock[phaseIdx][REF_PIC_LIST_1][0] + cFracBufOffset[(nDirect + 4) & 0x7];
    }
    else if (dir == 1)
    {
      pelBuffer[0] = m_filteredBlock[phaseIdx][REF_PIC_LIST_0][0] + cFracBufOffset[nDirect];
      pelBuffer[1] = m_filteredBlock[3][REF_PIC_LIST_1][0] + BDMVR_CENTER_POSITION;
    }
    else
    {
      pelBuffer[0] = m_filteredBlock[3][REF_PIC_LIST_0][0] + BDMVR_CENTER_POSITION;
      pelBuffer[1] = m_filteredBlock[phaseIdx][REF_PIC_LIST_1][0] + cFracBufOffset[(nDirect + 4) & 0x7];
    }
#else
    pelBuffer[0] = m_filteredBlock[phaseIdx][REF_PIC_LIST_0][0] + cFracBufOffset[  nDirect            ];
    pelBuffer[1] = m_filteredBlock[phaseIdx][REF_PIC_LIST_1][0] + cFracBufOffset[( nDirect + 4 ) & 0x7];
#endif
  }
  else
  {
    return xBDMVRGetMatchingError(pu, mv, useMR);
  }

  PelUnitBuf predBuf[2] = { PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_0], stride, pu.lwidth(), pu.lheight())),
                            PelUnitBuf(pu.chromaFormat, PelBuf(pelBuffer[REF_PIC_LIST_1], stride, pu.lwidth(), pu.lheight())) };

  // Compute distortion between L0'a and L1's prediction blocks
  DistParam cDistParam;
  cDistParam.applyWeight = false;
  cDistParam.useMR       = useMR;

  m_pcRdCost->setDistParam( cDistParam, predBuf[0].Y(), predBuf[1].Y(), pu.cu->slice->clpRng(COMPONENT_Y).bd, COMPONENT_Y, useHadmard );
#if FULL_NBIT
  if (useHadmard)
  {
    return cDistParam.distFunc( cDistParam ) >> 1;  // magic shift, benefit for early terminate
  }
  else
  {
    int32_t precisionAdj = cDistParam.bitDepth > 8 ? cDistParam.bitDepth - 8 : 0;
    return cDistParam.distFunc( cDistParam ) >> precisionAdj;
  }
#else
  return cDistParam.distFunc( cDistParam );
#endif
}
#endif
#endif

#if MULTI_HYP_PRED
void InterPrediction::xAddHypMC(PredictionUnit& pu, PelUnitBuf& predBuf, PelUnitBuf* predBufWOBIO, const bool lumaOnly)
{
  CHECK(pu.Y().area() <= MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE, "Multi Hyp: Block too small!");
  CHECK(pu.cu->geoFlag, "multi-hyp does not work with geo");
  CHECK(pu.ciipFlag, "multi-hyp does not work with intra/inter");
  CHECK(!pu.mergeFlag && pu.interDir != 3, "multihyp selected for AMVP uni prediction");

  // get prediction for current additional hypothesis
  const UnitArea unitAreaFromPredBuf(predBuf.chromaFormat, Area(Position(0, 0), predBuf.Y()));
  PelUnitBuf tempBuf = m_additionalHypothesisStorage.getBuf(unitAreaFromPredBuf);
  const auto savedAffine = pu.cu->affine;
  const auto savedIMV = pu.cu->imv;
#if INTER_LIC
  auto savedLICFlag = pu.cu->LICFlag;
#endif
  MultiHypVec savedHypVec = pu.addHypData;
  pu.addHypData.clear();
  pu.mvRefine = true;
  motionCompensation(pu, predBuf, REF_PIC_LIST_X, true, !lumaOnly, predBufWOBIO);
  pu.mvRefine = false;
#if INTER_LIC
  m_storeBeforeLIC = false;
#endif
  PredictionUnit fakePredData = pu;
  fakePredData.cu->affine = false;
  fakePredData.mergeFlag = false;
  fakePredData.mergeType = MRG_TYPE_DEFAULT_N;
  fakePredData.mmvdMergeFlag = false;
  fakePredData.ciipFlag = false;
#if MULTI_PASS_DMVR
  fakePredData.bdmvrRefine = false;
#endif

  for (int i = 0; i < savedHypVec.size(); i++)
  {
    const MultiHypPredictionData mhData = savedHypVec[i];

    // get legacy ref list and ref idx
    const auto &MHRefPics = pu.cs->slice->getMultiHypRefPicList();
    CHECK(mhData.refIdx < 0, "Multi Hyp: mhData.refIdx < 0");
    const int iRefPicList = mhData.isMrg ? mhData.refList : MHRefPics[mhData.refIdx].refList;
    const int iRefIdx = mhData.isMrg ? mhData.refIdx : MHRefPics[mhData.refIdx].refIdx;

    // construct fake object using legacy indexing
    fakePredData.interDir = iRefPicList + 1;
    fakePredData.mv[iRefPicList] = mhData.mv;
    fakePredData.refIdx[iRefPicList] = iRefIdx;
    fakePredData.refIdx[1 - iRefPicList] = -1;
#if INTER_LIC
    fakePredData.cu->LICFlag = mhData.LICFlag;
#endif
    fakePredData.cu->imv = mhData.imv;
    fakePredData.mvRefine = true;
    motionCompensation(fakePredData, tempBuf, REF_PIC_LIST_X, true, !lumaOnly);
    fakePredData.mvRefine = false;

    CHECK(mhData.weightIdx < 0, "Multi Hyp: mhData.weightIdx < 0");
    CHECK(mhData.weightIdx >= MULTI_HYP_PRED_NUM_WEIGHTS, "Multi Hyp: mhData.weightIdx >= MULTI_HYP_PRED_NUM_WEIGHTS");
    predBuf.addHypothesisAndClip(tempBuf, g_addHypWeight[mhData.weightIdx], pu.cs->slice->clpRngs(), lumaOnly);
  }
#if INTER_LIC
  pu.cu->LICFlag = savedLICFlag;
#endif
  pu.cu->imv = savedIMV;
  pu.cu->affine = savedAffine;
  pu.addHypData = savedHypVec;
}
#endif

#if JVET_X0083_BM_AMVP_MERGE_MODE
void InterPrediction::getAmvpMergeModeMergeList(PredictionUnit& pu, MvField* mvFieldAmListCommon, const int decAmvpRefIdx)
{
  RefPicList refListMerge = pu.amvpMergeModeFlag[0] ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
  RefPicList refListAmvp = RefPicList(1 - refListMerge);
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
  for (int idx = 0; idx < pu.cu->slice->getNumRefIdx(refListAmvp) * AMVP_MAX_NUM_CANDS_MEM; idx++)
#else
  for (int idx = 0; idx < pu.cu->slice->getNumRefIdx(refListAmvp) * AMVP_MAX_NUM_CANDS; idx++)
#endif
  {
    mvFieldAmListCommon[idx] = MvField();
    mvFieldAmListCommon[MAX_NUM_AMVP_CANDS_MAX_REF + idx] = MvField();
  }
  int amvpRefIdxStart = 0;
  int amvpRefIdxEnd = pu.cu->slice->getNumRefIdx(refListAmvp);
  int decAmvpMvpIdx = -1;
  if (decAmvpRefIdx >= 0)
  {
    amvpRefIdxStart = decAmvpRefIdx;
    amvpRefIdxEnd = decAmvpRefIdx + 1;
    decAmvpMvpIdx = pu.mvpIdx[refListAmvp];
  }
#if !JVET_Y0128_NON_CTC
  const int curPoc = pu.cu->slice->getPOC();
#endif
  const bool useMR = pu.lumaSize().area() > 64;

  for (int refIdxAmvp = amvpRefIdxStart; refIdxAmvp < amvpRefIdxEnd; refIdxAmvp++)
  {
#if JVET_Y0128_NON_CTC
    if (pu.cu->slice->getAmvpMergeModeValidRefIdx(refListAmvp, refIdxAmvp) == false)
    {
      continue;
    }
    CHECK(pu.cu->slice->getRefPic(refListAmvp, refIdxAmvp)->isRefScaled(pu.cu->cs->pps), "this is not possible");
#else
    const int amvpRefPoc = pu.cu->slice->getRefPOC(refListAmvp, refIdxAmvp);
    bool findValidMergeRefPic = false;
    for (int refIdxCandMerge = 0; refIdxCandMerge < pu.cu->slice->getNumRefIdx(refListMerge); refIdxCandMerge++)
    {
      const int candMergePoc = pu.cu->slice->getRefPOC(refListMerge, refIdxCandMerge);
      if ((amvpRefPoc - curPoc) * (candMergePoc - curPoc) < 0)
      {
        findValidMergeRefPic = true;
        break;
      }
    }
    if (findValidMergeRefPic == false)
    {
      continue;
    }
#endif
    pu.refIdx[refListAmvp] = refIdxAmvp;
    AMVPInfo amvpInfo;
    PU::fillMvpCand( pu, refListAmvp, refIdxAmvp, amvpInfo
#if TM_AMVP
                   , this
#endif
    );
    MergeCtx bmMergeCtx;
    PU::getInterMergeCandidates(pu, bmMergeCtx, 0, AMVP_MERGE_MODE_MERGE_LIST_MAX_CANDS - 1);

    int bestMvpIdxLoopStart = 0;
    int bestMvpIdxLoopEnd = amvpInfo.numCand;
    if (decAmvpRefIdx >= 0)
    {
      bestMvpIdxLoopStart = decAmvpMvpIdx;
      bestMvpIdxLoopEnd = bestMvpIdxLoopStart + 1;
    }
    for (int bestMvpIdxToTest = bestMvpIdxLoopStart; bestMvpIdxToTest < bestMvpIdxLoopEnd; bestMvpIdxToTest++)
    {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      const int mvFieldMergeIdx = refIdxAmvp * AMVP_MAX_NUM_CANDS_MEM + bestMvpIdxToTest;
#else
      const int mvFieldMergeIdx = refIdxAmvp * AMVP_MAX_NUM_CANDS + bestMvpIdxToTest;
#endif
      const int mvFieldAmvpIdx = MAX_NUM_AMVP_CANDS_MAX_REF + mvFieldMergeIdx;
      pu.mv[refListAmvp] = amvpInfo.mvCand[bestMvpIdxToTest];

      // BM select merge candidate
      struct bmCostSort
      {
        int mergeIdx;
        Distortion bmCost;
      };
      bmCostSort temp;
      const auto CostIncSort = [](const bmCostSort &x, const bmCostSort &y) { return x.bmCost < y.bmCost; };
      std::vector<bmCostSort> input;
      // here to select the merge cand which has minimum BM cost, at each cand, the cost is derived by  minBMcost(mvpIdx0, mvpIdx1)
      if (bmMergeCtx.numValidMergeCand > 1)
      {
        // pre Fill AMVP prediction blocks
#if JVET_X0049_BDMVR_SW_OPT
        Pel* pelBufferAmvp = m_filteredBlock[3][refListAmvp][0] + BDMVR_CENTER_POSITION;
        const SizeType stride = BDMVR_BUF_STRIDE;
#else
        Pel* pelBufferAmvp        = m_filteredBlock[3][refListAmvp][0];
        const SizeType stride     = pu.lwidth();
#endif
        PelUnitBuf predBufAmvp    = PelUnitBuf(pu.chromaFormat, PelBuf(pelBufferAmvp, stride, pu.lwidth(), pu.lheight()));
        const Picture& refPicAmvp = *pu.cu->slice->getRefPic((RefPicList)refListAmvp, pu.refIdx[refListAmvp])->unscaledPic;
        xBDMVRFillBlkPredPelBuffer( pu, refPicAmvp, pu.mv[refListAmvp] , predBufAmvp, pu.cs->slice->clpRng(COMPONENT_Y) );
        Mv mvAmBdmvr[2/*refListId*/];
        for (int mergeIdx = 0; mergeIdx < bmMergeCtx.numValidMergeCand; mergeIdx++)
        {
          pu.refIdx[refListMerge] = bmMergeCtx.mvFieldNeighbours[(mergeIdx << 1) + refListMerge].refIdx;
          mvAmBdmvr[refListAmvp] = amvpInfo.mvCand[bestMvpIdxToTest];
          mvAmBdmvr[refListMerge] = bmMergeCtx.mvFieldNeighbours[(mergeIdx << 1) + refListMerge].mv;
#if JVET_Y0128_NON_CTC
          CHECK(pu.cu->slice->getRefPic((RefPicList)refListMerge, pu.refIdx[refListMerge])->isRefScaled(pu.cs->pps), "this is not possible");
#endif
#if JVET_Z0085_AMVPMERGE_DMVD_OFF
          if (pu.cu->cs->sps->getUseDMVDMode())
          {
#endif
            Distortion tmpBmCost = xBDMVRGetMatchingError(pu, mvAmBdmvr, useMR);
            temp.mergeIdx = mergeIdx;
            temp.bmCost = tmpBmCost;
#if JVET_Z0085_AMVPMERGE_DMVD_OFF
          }
          else
          {
            temp.mergeIdx = mergeIdx;
            temp.bmCost = std::numeric_limits<Distortion>::max();
          }
#endif
          input.push_back(temp);
        }
        stable_sort(input.begin(), input.end(), CostIncSort);
      }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      else
      {
        temp.mergeIdx = 0;
        temp.bmCost = 0;
        input.push_back(temp);
      }
#else
      if (bmMergeCtx.numValidMergeCand == 1)
      {
        pu.mv[refListMerge] = bmMergeCtx.mvFieldNeighbours[refListMerge].mv;
        pu.refIdx[refListMerge] = bmMergeCtx.mvFieldNeighbours[refListMerge].refIdx;
      }
      else
#endif
      {
        pu.mv[refListMerge] = bmMergeCtx.mvFieldNeighbours[(input[0].mergeIdx << 1) + refListMerge].mv;
        pu.refIdx[refListMerge] = bmMergeCtx.mvFieldNeighbours[(input[0].mergeIdx << 1) + refListMerge].refIdx;
      }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      if (bestMvpIdxToTest == 0 || bestMvpIdxToTest == 2)
      {
#endif
      amvpMergeModeMvRefinement(pu, mvFieldAmListCommon, mvFieldMergeIdx, mvFieldAmvpIdx);
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      }
      else if (bmMergeCtx.numValidMergeCand == 1)
      {
        mvFieldAmListCommon[mvFieldMergeIdx].refIdx = bmMergeCtx.mvFieldNeighbours[(input[0].mergeIdx << 1) + refListMerge].refIdx;
        mvFieldAmListCommon[mvFieldMergeIdx].mv = bmMergeCtx.mvFieldNeighbours[(input[0].mergeIdx << 1) + refListMerge].mv;
        mvFieldAmListCommon[mvFieldAmvpIdx].refIdx = refIdxAmvp;
        mvFieldAmListCommon[mvFieldAmvpIdx].mv = amvpInfo.mvCand[bestMvpIdxToTest];
      }
      else
      {
        pu.mv[refListAmvp] = amvpInfo.mvCand[bestMvpIdxToTest];
        pu.refIdx[refListAmvp] = refIdxAmvp;
        pu.mv[refListMerge] = bmMergeCtx.mvFieldNeighbours[(input[1].mergeIdx << 1) + refListMerge].mv;
        pu.refIdx[refListMerge] = bmMergeCtx.mvFieldNeighbours[(input[1].mergeIdx << 1) + refListMerge].refIdx;
        amvpMergeModeMvRefinement(pu, mvFieldAmListCommon, mvFieldMergeIdx, mvFieldAmvpIdx);
      }
      if (bestMvpIdxToTest == 2)
      {
        mvFieldAmListCommon[mvFieldAmvpIdx].mv.roundTransPrecInternal2Amvr(pu.cu->imv);
      }
#endif
    } // bestMvpIdxLoop
  }  // refIdxAmvp loop
}
void InterPrediction::amvpMergeModeMvRefinement(PredictionUnit& pu, MvField* mvFieldAmListCommon, const int mvFieldMergeIdx, const int mvFieldAmvpIdx)
{
  const RefPicList refListMerge = pu.amvpMergeModeFlag[0] ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
  const RefPicList refListAmvp = RefPicList(1 - refListMerge);
  const int curPoc = pu.cu->slice->getPOC();
  const int mergeRefPoc = pu.cu->slice->getRefPOC(refListMerge, pu.refIdx[refListMerge]);
  const bool useMR = pu.lumaSize().area() > 64;
  const int amvpRefPoc = pu.cu->slice->getRefPOC(refListAmvp, pu.refIdx[refListAmvp]);
#if JVET_Y0128_NON_CTC
  CHECK(pu.cu->slice->getRefPic(REF_PIC_LIST_0, pu.refIdx[0])->isRefScaled(pu.cs->pps), "this is not possible");
  CHECK(pu.cu->slice->getRefPic(REF_PIC_LIST_1, pu.refIdx[1])->isRefScaled(pu.cs->pps), "this is not possible");
#endif
#if JVET_Z0085_AMVPMERGE_DMVD_OFF
  if (pu.cu->cs->sps->getUseDMVDMode())
  {
#endif
  if ((mergeRefPoc - curPoc) == (curPoc - amvpRefPoc))
  {
    Mv         mvInitial[2];
    mvInitial[refListAmvp] = pu.mv[refListAmvp];;
    mvInitial[refListMerge] = pu.mv[refListMerge];
    Mv         mvFinal[2] = { mvInitial[0], mvInitial[1] };
    Distortion curBmCost = std::numeric_limits<Distortion>::max();
#if JVET_X0049_BDMVR_SW_OPT
    curBmCost = xBDMVRMvSquareSearch<false>(mvFinal, curBmCost, pu, mvInitial,
        AMVP_MERGE_MODE_REDUCED_MV_REFINE_SEARCH_ROUND, MV_FRACTIONAL_BITS_INTERNAL, useMR, false);
    curBmCost = xBDMVRMvSquareSearch<true>(mvFinal, curBmCost, pu, mvInitial,
      2, MV_FRACTIONAL_BITS_INTERNAL - 1, useMR, false);
#else
    curBmCost = xBDMVRMvSquareSearch( mvFinal, curBmCost, pu, mvInitial,
        AMVP_MERGE_MODE_REDUCED_MV_REFINE_SEARCH_ROUND, MV_FRACTIONAL_BITS_INTERNAL, useMR, false );
    curBmCost = xBDMVRMvSquareSearch( mvFinal, curBmCost, pu, mvInitial,
        2, MV_FRACTIONAL_BITS_INTERNAL - 1, useMR, false );
#endif
    pu.mv[refListMerge] = mvFinal[refListMerge];
    pu.mv[refListAmvp] = mvFinal[refListAmvp];
  }
#if TM_AMVP || TM_MRG
  else
  {
    Distortion tmCost[2];
    tmCost[refListMerge] = deriveTMMv(pu, true, std::numeric_limits<Distortion>::max(), refListMerge, pu.refIdx[refListMerge], 0, pu.mv[refListMerge]);
    tmCost[refListAmvp] = deriveTMMv(pu, true, std::numeric_limits<Distortion>::max(), refListAmvp, pu.refIdx[refListAmvp], 0, pu.mv[refListAmvp]);
    RefPicList refListToBeRefined = (tmCost[refListMerge] < tmCost[refListAmvp]) ? refListAmvp : refListMerge;
    MvField    mvfBetterUni(pu.mv[1 - refListToBeRefined], pu.refIdx[1 - refListToBeRefined]);
    deriveTMMv(pu, true, std::numeric_limits<Distortion>::max(), refListToBeRefined, pu.refIdx[refListToBeRefined],
        AMVP_MERGE_MODE_REDUCED_MV_REFINE_SEARCH_ROUND, pu.mv[refListToBeRefined], &mvfBetterUni);
    }
#endif
#if JVET_Z0085_AMVPMERGE_DMVD_OFF
  }
#endif
#if !JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
  pu.mv[refListAmvp].roundTransPrecInternal2Amvr(pu.cu->imv);
#endif
  mvFieldAmListCommon[mvFieldMergeIdx].refIdx = pu.refIdx[refListMerge];
  mvFieldAmListCommon[mvFieldMergeIdx].mv = pu.mv[refListMerge];
  mvFieldAmListCommon[mvFieldAmvpIdx].refIdx = pu.refIdx[refListAmvp];
  mvFieldAmListCommon[mvFieldAmvpIdx].mv = pu.mv[refListAmvp];
}
#endif


#if JVET_Z0054_BLK_REF_PIC_REORDER
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
void InterPrediction::deriveMVDCandVecFromMotionInforPredGeneral(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> &cMvdDerivedVec)
{
  cMvdDerivedVec.clear();

  deriveMVDcand(pu, eRefPicList, cMvdDerivedVec);
  if (!cMvdDerivedVec.empty() && !miPredList.empty())
  {
    std::stable_sort(miPredList.begin(), miPredList.end(), [](const MotionInfoPred & l, const MotionInfoPred & r) {return l.cost < r.cost; });

    std::vector<MotionInfoPred> miPredListSub;
    miPredListSub.clear();
    {
      for (std::vector<MotionInfoPred>::iterator it = miPredList.begin(); it != miPredList.end(); ++it)
      {
        if (it->interDir == pu.interDir && it->refIdx[eRefPicList] == pu.refIdx[eRefPicList])
        {
          if (it->interDir == 3 && eRefPicList == REF_PIC_LIST_1
            && (it->mvd[0] != pu.mvd[0]))
          {
            continue;
          }
          bool add = true;
          for (std::vector<MotionInfoPred>::iterator itSub = miPredListSub.begin(); itSub != miPredListSub.end(); ++itSub)
          {
            if (it->mvd[eRefPicList] == itSub->mvd[eRefPicList])
            {
              add = false;
              break;
            }
          }
          if (add)
          {
            miPredListSub.push_back(*it);
          }
        }
      }
      CHECK(!miPredListSub.empty() && cMvdDerivedVec.size() != miPredListSub.size(), "cMvdDerivedVec.size() != miPredListSub.size()");
      // when template doesn't exist, miPredListSub is empty but cMvdDerivedVec might not be.
      for (size_t i = 0; i < miPredListSub.size(); i++)
      {
        cMvdDerivedVec[i] = (miPredListSub[i].mvd[eRefPicList]);
      }
    }
  }
}

void InterPrediction::deriveAffineMVDCandVecFromMotionInforPredGeneral(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> cMvdDerivedVec[3])
{
  cMvdDerivedVec[0].clear();
  cMvdDerivedVec[1].clear();
  cMvdDerivedVec[2].clear();

  deriveMVDcandAffine(pu, eRefPicList, cMvdDerivedVec);
  if (!cMvdDerivedVec[0].empty() && !miPredList.empty())
  {
    std::stable_sort(miPredList.begin(), miPredList.end(), [](const MotionInfoPred & l, const MotionInfoPred & r) {return l.cost < r.cost; });

    std::vector<MotionInfoPred> miPredListSub;
    miPredListSub.clear();
    for (std::vector<MotionInfoPred>::iterator it = miPredList.begin(); it != miPredList.end(); ++it)
    {
      if (it->interDir == pu.interDir && it->refIdx[eRefPicList] == pu.refIdx[eRefPicList])
      {
        if (it->interDir == 3 && eRefPicList == REF_PIC_LIST_1
          && (it->mvdAffi[0][0] != pu.mvdAffi[0][0] || it->mvdAffi[0][1] != pu.mvdAffi[0][1] || (pu.cu->affineType == AFFINEMODEL_6PARAM && it->mvdAffi[0][2] != pu.mvdAffi[0][2])))
        {
          continue;
        }
        bool add = true;
        for (std::vector<MotionInfoPred>::iterator itSub = miPredListSub.begin(); itSub != miPredListSub.end(); ++itSub)
        {
          if (it->mvdAffi[eRefPicList][0] == itSub->mvdAffi[eRefPicList][0]
            && it->mvdAffi[eRefPicList][1] == itSub->mvdAffi[eRefPicList][1]
            && (pu.cu->affineType == AFFINEMODEL_4PARAM || it->mvdAffi[eRefPicList][2] == itSub->mvdAffi[eRefPicList][2]))
          {
            add = false;
            break;
          }
        }
        if (add)
        {
          miPredListSub.push_back(*it);
        }
      }
    }
    CHECK(cMvdDerivedVec[0].size() != miPredListSub.size(), "cMvdDerivedVec[0].size() != miPredListSub.size()");
    for (size_t i = 0; i < miPredListSub.size(); i++)
    {
      cMvdDerivedVec[0][i] = (miPredListSub[i].mvdAffi[eRefPicList][0]);
      cMvdDerivedVec[1][i] = (miPredListSub[i].mvdAffi[eRefPicList][1]);
      cMvdDerivedVec[2][i] = (miPredListSub[i].mvdAffi[eRefPicList][2]);
    }
  }
}

void InterPrediction::deriveMVDCandVecFromMotionInforPred(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> &cMvdDerivedVec)
{
  cMvdDerivedVec.clear();

  deriveMVDcand(pu, eRefPicList, cMvdDerivedVec);
  if (!cMvdDerivedVec.empty())
  {
    std::vector<MotionInfoPred> miPredListSub;
    miPredListSub.clear();
    if (!miPredList.empty())
    {
      for (std::vector<MotionInfoPred>::iterator it = miPredList.begin(); it != miPredList.end(); ++it)
      {
        if (it->interDir == pu.interDir && it->refIdx[eRefPicList] == pu.refIdx[eRefPicList])
        {
          miPredListSub.push_back(*it);
        }
      }
      std::stable_sort(miPredListSub.begin(), miPredListSub.end(), [](const MotionInfoPred & l, const MotionInfoPred & r) {return l.cost < r.cost; });
    }
    CHECK(!miPredListSub.empty() && cMvdDerivedVec.size() != miPredListSub.size(), "cMvdDerivedVec.size() != miPredListSub.size()");
    // when template doesn't exist, miPredListSub is empty but cMvdDerivedVec might not be.
    for (size_t i = 0; i < miPredListSub.size(); i++)
    {
      cMvdDerivedVec[i] = (miPredListSub[i].mvd[eRefPicList]);
    }
  }
}

void InterPrediction::deriveAffineMVDCandVecFromMotionInforPred(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> cMvdDerivedVec[3])
{
  cMvdDerivedVec[0].clear();
  cMvdDerivedVec[1].clear();
  cMvdDerivedVec[2].clear();

  deriveMVDcandAffine(pu, eRefPicList, cMvdDerivedVec);
  if (!cMvdDerivedVec[0].empty())
  {
    std::vector<MotionInfoPred> miPredListSub;
    miPredListSub.clear();
    if (!miPredList.empty())
    {
      for (std::vector<MotionInfoPred>::iterator it = miPredList.begin(); it != miPredList.end(); ++it)
      {
        if (it->interDir == pu.interDir && it->refIdx[eRefPicList] == pu.refIdx[eRefPicList])
        {
          miPredListSub.push_back(*it);
        }
      }
      std::stable_sort(miPredListSub.begin(), miPredListSub.end(), [](const MotionInfoPred & l, const MotionInfoPred & r) {return l.cost < r.cost; });
    }
    CHECK(!miPredListSub.empty() && cMvdDerivedVec[0].size() != miPredListSub.size(), "cMvdDerivedVec[0].size() != miPredListSub.size()");
    for (size_t i = 0; i < miPredListSub.size(); i++)
    {
      cMvdDerivedVec[0][i] = (miPredListSub[i].mvdAffi[eRefPicList][0]);
      cMvdDerivedVec[1][i] = (miPredListSub[i].mvdAffi[eRefPicList][1]);
      cMvdDerivedVec[2][i] = (miPredListSub[i].mvdAffi[eRefPicList][2]);
    }
  }
}
#endif

void InterPrediction::reorderRefCombList(PredictionUnit &pu, std::vector<RefListAndRefIdx> &refListComb
  , RefPicList currRefList
  , std::vector<MotionInfoPred> &miPredList
)
{
  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  if (!xAMLGetCurBlkTemplate(pu, nWidth, nHeight) && !pu.isMvsdApplicable())
#else
  if (refListComb.size() < 2 || !xAMLGetCurBlkTemplate(pu, nWidth, nHeight))
#endif
  {
    return;
  }

  PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
  PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
  PredictionUnit tmpPU = pu;

  DistParam cDistParam;
  cDistParam.applyWeight = false;
  Distortion uiCost;

  if (pu.cu->affine)
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    std::vector<Mv> cMvdCandList[3];
    cMvdCandList[0].push_back(tmpPU.mvdAffi[currRefList][0]);
    cMvdCandList[1].push_back(tmpPU.mvdAffi[currRefList][1]);
    cMvdCandList[2].push_back(tmpPU.mvdAffi[currRefList][2]);

    if (pu.isMvsdApplicable())
    {
      deriveMVDcandAffine(tmpPU, currRefList, cMvdCandList);
    }
#endif
    for (int idx = 0; idx < refListComb.size(); idx++)
    {
      RefPicList eRefList = refListComb[idx].refList;
      int refIdx = refListComb[idx].refIdx;
      tmpPU.interDir = 1 << eRefList;
      tmpPU.refIdx[1 - eRefList] = -1;
      tmpPU.refIdx[eRefList] = refIdx;
      tmpPU.mvpIdx[eRefList] = pu.mvpIdx[currRefList];
      tmpPU.mvdAffi[eRefList][0] = tmpPU.mvdAffi[currRefList][0];
      tmpPU.mvdAffi[eRefList][1] = tmpPU.mvdAffi[currRefList][1];
      tmpPU.mvdAffi[eRefList][2] = tmpPU.mvdAffi[currRefList][2];


      AffineAMVPInfo affineAMVPInfo;
      PU::fillAffineMvpCand(tmpPU, eRefList, tmpPU.refIdx[eRefList], affineAMVPInfo);

      const unsigned mvp_idx = tmpPU.mvpIdx[eRefList];

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      refListComb[idx].cost = std::numeric_limits<Distortion>::max();
      for (size_t i = 0; i < cMvdCandList[0].size(); i++)
      {
        Mv mvLT = affineAMVPInfo.mvCandLT[mvp_idx] + cMvdCandList[0][i];
        Mv mvRT = affineAMVPInfo.mvCandRT[mvp_idx] + cMvdCandList[1][i];
        mvRT += cMvdCandList[0][i];

        Mv mvLB;
        if (tmpPU.cu->affineType == AFFINEMODEL_6PARAM)
        {
          mvLB = affineAMVPInfo.mvCandLB[mvp_idx] + cMvdCandList[2][i];
          mvLB += cMvdCandList[0][i];
        }
#else
      Mv mvLT = affineAMVPInfo.mvCandLT[mvp_idx] + tmpPU.mvdAffi[eRefList][0];
      Mv mvRT = affineAMVPInfo.mvCandRT[mvp_idx] + tmpPU.mvdAffi[eRefList][1];
      mvRT += tmpPU.mvdAffi[eRefList][0];

      Mv mvLB;
      if (tmpPU.cu->affineType == AFFINEMODEL_6PARAM)
      {
        mvLB = affineAMVPInfo.mvCandLB[mvp_idx] + tmpPU.mvdAffi[eRefList][2];
        mvLB += tmpPU.mvdAffi[eRefList][0];
      }
#endif

      tmpPU.mvAffi[eRefList][0] = mvLT;
      tmpPU.mvAffi[eRefList][1] = mvRT;
      tmpPU.mvAffi[eRefList][2] = mvLB;

#if !JVET_Z0067_RPR_ENABLE
      getAffAMLRefTemplate(tmpPU, pcBufPredRefTop, pcBufPredRefLeft);
#endif

      uiCost = 0;
      bool bRefIsRescaled = (tmpPU.refIdx[eRefList] >= 0) ? tmpPU.cu->slice->getRefPic(eRefList, tmpPU.refIdx[eRefList])->isRefScaled(pu.cs->pps) : false;
      if (bRefIsRescaled)
      {
        uiCost = std::numeric_limits<Distortion>::max();
      }
      else
      {
#if JVET_Z0067_RPR_ENABLE
        getAffAMLRefTemplate(tmpPU, pcBufPredRefTop, pcBufPredRefLeft);
#endif
        if (m_bAMLTemplateAvailabe[0])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }

        if (m_bAMLTemplateAvailabe[1])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }
      }
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      MotionInfoPred miPred;
      miPred.interDir = 1 << eRefList;
      miPred.refIdx[eRefList] = refIdx;
      miPred.mvdAffi[eRefList][0] = cMvdCandList[0][i];
      miPred.mvdAffi[eRefList][1] = cMvdCandList[1][i];
      miPred.mvdAffi[eRefList][2] = cMvdCandList[2][i];
      miPred.mvAffi[eRefList][0] = tmpPU.mvAffi[eRefList][0];
      miPred.mvAffi[eRefList][1] = tmpPU.mvAffi[eRefList][1];
      miPred.mvAffi[eRefList][2] = tmpPU.mvAffi[eRefList][2];
      miPred.cost = uiCost;
      miPredList.push_back(miPred);

      if (uiCost < refListComb[idx].cost)
      {
        refListComb[idx].cost = uiCost;
      }
      }
#else
      refListComb[idx].cost = uiCost;
#endif
    }
  }
  else
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    std::vector<Mv> cMvdCandList;
    cMvdCandList.push_back(tmpPU.mvd[currRefList]);
    if (pu.isMvsdApplicable())
    {
      deriveMVDcand(tmpPU, currRefList, cMvdCandList);
    }
#endif
    for (int idx = 0; idx < refListComb.size(); idx++)
    {
      RefPicList eRefList = refListComb[idx].refList;
      int refIdx = refListComb[idx].refIdx;

      tmpPU.interDir = 1 << eRefList;
      tmpPU.refIdx[1 - eRefList] = -1;
      tmpPU.refIdx[eRefList] = refIdx;
      tmpPU.mvpIdx[eRefList] = tmpPU.mvpIdx[currRefList];

      AMVPInfo amvpInfo;
      PU::fillMvpCand(tmpPU, eRefList, tmpPU.refIdx[eRefList], amvpInfo
#if TM_AMVP
        , this
#endif
      );

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      refListComb[idx].cost = std::numeric_limits<Distortion>::max();
      for (std::vector<Mv>::iterator it = cMvdCandList.begin(); it != cMvdCandList.end(); ++it)
      {
        tmpPU.mvd[eRefList] = *it;
#else
      tmpPU.mvd[eRefList] = tmpPU.mvd[currRefList];
#endif
      tmpPU.mv[eRefList] = amvpInfo.mvCand[tmpPU.mvpIdx[eRefList]] + tmpPU.mvd[eRefList];
      tmpPU.mv[eRefList].mvCliptoStorageBitDepth();

      uiCost = 0;
      bool bRefIsRescaled = (tmpPU.refIdx[eRefList] >= 0) ? tmpPU.cu->slice->getRefPic(eRefList, tmpPU.refIdx[eRefList])->isRefScaled(pu.cs->pps) : false;
      if (bRefIsRescaled)
      {
        uiCost = std::numeric_limits<Distortion>::max();
      }
      else
      {
        getBlkAMLRefTemplate(tmpPU, pcBufPredRefTop, pcBufPredRefLeft);
        if (m_bAMLTemplateAvailabe[0])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }

        if (m_bAMLTemplateAvailabe[1])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }
      }

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      MotionInfoPred miPred;
      miPred.interDir = 1 << eRefList;
      miPred.refIdx[eRefList] = refIdx;
      miPred.mvd[eRefList] = tmpPU.mvd[eRefList];
      miPred.mv[eRefList] = tmpPU.mv[eRefList];
      miPred.cost = uiCost;
      miPredList.push_back(miPred);

      if (uiCost < refListComb[idx].cost)
      {
        refListComb[idx].cost = uiCost;
      }
      }
#else
      refListComb[idx].cost = uiCost;
#endif
    }
  }
  std::stable_sort(refListComb.begin(), refListComb.end(), [](const RefListAndRefIdx & l, const RefListAndRefIdx & r) {return l.cost < r.cost; });
}

void InterPrediction::setUniRefIdxLC(PredictionUnit &pu)
{
  RefPicList eRefPicList;
  std::vector<RefListAndRefIdx> refListComb;
#if JVET_X0083_BM_AMVP_MERGE_MODE
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    eRefPicList = pu.amvpMergeModeFlag[0] ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
    refListComb = pu.cs->slice->getRefPicCombinedListAmvpMerge();
  }
  else
#endif
  {
    eRefPicList = pu.interDir == 1 ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
    refListComb = pu.cs->slice->getRefPicCombinedList();
    std::vector<MotionInfoPred> miPredList;
    miPredList.clear();
    reorderRefCombList(pu, refListComb, RefPicList(pu.interDir >> 1), miPredList);

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    if (pu.isMvsdApplicable())
    {
      if (pu.cu->affine)
      {
        std::vector<Mv> cMvdDerivedVec[3];
        deriveAffineMVDCandVecFromMotionInforPred(pu, miPredList, (pu.interDir == 1 ? REF_PIC_LIST_0 : REF_PIC_LIST_1), cMvdDerivedVec);
        pu.mvsdIdx[eRefPicList] = deriveMVSDIdxFromMVDAffine(pu, eRefPicList, cMvdDerivedVec[0], cMvdDerivedVec[1], cMvdDerivedVec[2]);
      }
      else
      {
        std::vector<Mv> cMvdDerivedVec;
        deriveMVDCandVecFromMotionInforPred(pu, miPredList, eRefPicList, cMvdDerivedVec);
        pu.mvsdIdx[eRefPicList] = deriveMVSDIdxFromMVDTrans(pu.mvd[eRefPicList], cMvdDerivedVec);
      }
    }
#endif
  }

  int refIdx = pu.refIdx[eRefPicList];
  for (int8_t idx = 0; idx < refListComb.size(); idx++)
  {
    if (refListComb[idx].refList == eRefPicList && refListComb[idx].refIdx == refIdx)
    {
      pu.refIdxLC = idx;
      break;
    }
  }
}

void InterPrediction::setUniRefListAndIdx(PredictionUnit &pu)
{
  RefPicList eRefList;
  std::vector<RefListAndRefIdx> refListComb;
  std::vector<MotionInfoPred> miPredList;
  miPredList.clear();

#if JVET_X0083_BM_AMVP_MERGE_MODE
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    refListComb = pu.cs->slice->getRefPicCombinedListAmvpMerge();
    eRefList = refListComb[pu.refIdxLC].refList;
    pu.refIdx[eRefList] = refListComb[pu.refIdxLC].refIdx;
    pu.interDir = 3;
    pu.amvpMergeModeFlag[0] = (eRefList == REF_PIC_LIST_0 ? false : true);
    pu.amvpMergeModeFlag[1] = (eRefList == REF_PIC_LIST_0 ? true : false);
  }
  else
#endif
  {
    refListComb = pu.cs->slice->getRefPicCombinedList();
    reorderRefCombList(pu, refListComb, RefPicList(pu.interDir >> 1), miPredList);
    eRefList = refListComb[pu.refIdxLC].refList;
    pu.interDir = 1 << eRefList;
    pu.refIdx[eRefList] = refListComb[pu.refIdxLC].refIdx;
    pu.refIdx[1 - eRefList] = -1; // some code relies on whether refIdx equals to -1 instead of interDir to determine inter prediction direction
  }
  //move other motion informations from temporally locations to the correct ones
  if (pu.cu->affine)
  {
    pu.mvdAffi[eRefList][0] = pu.mvdAffi[0][0];
    pu.mvdAffi[eRefList][1] = pu.mvdAffi[0][1];
    pu.mvdAffi[eRefList][2] = pu.mvdAffi[0][2];
  }
  else
  {
    pu.mvd[eRefList] = pu.mvd[0];
  }
  pu.mvpIdx[eRefList] = pu.mvpIdx[0];
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  pu.mvsdIdx[eRefList] = pu.mvsdIdx[0];
#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if JVET_X0083_BM_AMVP_MERGE_MODE
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    return;
  }
#endif
  if (pu.isMvsdApplicable())
  {
    if (pu.cu->affine)
    {
      std::vector<Mv> cMvdDerivedVec[3];
      deriveAffineMVDCandVecFromMotionInforPred(pu, miPredList, (pu.interDir == 1 ? REF_PIC_LIST_0 : REF_PIC_LIST_1), cMvdDerivedVec);
      deriveMVDFromMVSDIdxAffine(pu, eRefList, cMvdDerivedVec[0], cMvdDerivedVec[1], cMvdDerivedVec[2]);
    }
    else
    {
      std::vector<Mv> cMvdDerivedVec;
      deriveMVDCandVecFromMotionInforPred(pu, miPredList, eRefList, cMvdDerivedVec);
      if (!cMvdDerivedVec.empty())
      {
        int mvsdIdx = pu.mvsdIdx[eRefList];
        pu.mvd[eRefList] = deriveMVDFromMVSDIdxTrans(mvsdIdx, cMvdDerivedVec);
      }
    }
  }
#endif
}

void InterPrediction::setBiRefPairIdx(PredictionUnit &pu)
{
  std::vector<RefPicPair>    refPicPairList = pu.cs->slice->getRefPicPairList();
  std::vector<MotionInfoPred> miPredList;
  miPredList.clear();
  reorderRefPairList(pu, refPicPairList, miPredList);
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  if (pu.isMvsdApplicable())
  {
    uint8_t numCandL0 = 0;
    for (int8_t refList = 0; refList < 2; refList++)
    {
      RefPicList eRefPicList = (RefPicList)refList;
      if (refList && numCandL0 >= 2)
      {
        continue;
      }
      if (pu.cu->affine)
      {
        std::vector<Mv> cMvdDerivedVec[3];
        deriveAffineMVDCandVecFromMotionInforPredGeneral(pu, miPredList, eRefPicList, cMvdDerivedVec);
        numCandL0 = (uint8_t)cMvdDerivedVec[0].size();
        pu.mvsdIdx[eRefPicList] = deriveMVSDIdxFromMVDAffine(pu, eRefPicList, cMvdDerivedVec[0], cMvdDerivedVec[1], cMvdDerivedVec[2]);
      }
      else
      {
        std::vector<Mv> cMvdDerivedVec;
        deriveMVDCandVecFromMotionInforPredGeneral(pu, miPredList, eRefPicList, cMvdDerivedVec);
        numCandL0 = (uint8_t)cMvdDerivedVec.size();
        pu.mvsdIdx[eRefPicList] = deriveMVSDIdxFromMVDTrans(pu.mvd[eRefPicList], cMvdDerivedVec);
      }
    }
  }
#endif
  for (int8_t idx = 0; idx < refPicPairList.size(); idx++)
  {
    if (refPicPairList[idx].refIdx[0] == pu.refIdx[0] && refPicPairList[idx].refIdx[1] == pu.refIdx[1])
    {
      pu.refPairIdx = idx;
      break;
    }
  }
  CHECK(pu.refPairIdx < 0, "");
}

void InterPrediction::setBiRefIdx(PredictionUnit &pu)
{
  std::vector<RefPicPair>    refPicPairList = pu.cs->slice->getRefPicPairList();
  std::vector<MotionInfoPred> miPredList;
  miPredList.clear();
  reorderRefPairList(pu, refPicPairList, miPredList);
  pu.refIdx[0] = refPicPairList[pu.refPairIdx].refIdx[0];
  pu.refIdx[1] = refPicPairList[pu.refPairIdx].refIdx[1];
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  if (pu.isMvsdApplicable())
  {
    uint8_t numCandL0 = 0;
    for (int8_t refList = 0; refList < 2; refList++)
    {
      RefPicList eRefList = (RefPicList)refList;
      if (refList && numCandL0 >= 2)
      {
        continue;
      }
      if (pu.cu->affine)
      {
        std::vector<Mv> cMvdDerivedVec[3];
        deriveAffineMVDCandVecFromMotionInforPredGeneral(pu, miPredList, eRefList, cMvdDerivedVec);
        numCandL0 = (uint8_t)cMvdDerivedVec[0].size();
        deriveMVDFromMVSDIdxAffine(pu, eRefList, cMvdDerivedVec[0], cMvdDerivedVec[1], cMvdDerivedVec[2]);
      }
      else
      {
        std::vector<Mv> cMvdDerivedVec;
        deriveMVDCandVecFromMotionInforPredGeneral(pu, miPredList, eRefList, cMvdDerivedVec);
        numCandL0 = (uint8_t)cMvdDerivedVec.size();
        if (!cMvdDerivedVec.empty())
        {
          int mvsdIdx = pu.mvsdIdx[eRefList];
          pu.mvd[eRefList] = deriveMVDFromMVSDIdxTrans(mvsdIdx, cMvdDerivedVec);
        }
      }
    }
  }
#endif
}

void InterPrediction::reorderRefPairList(PredictionUnit &pu, std::vector<RefPicPair> &refPairList
  , std::vector<MotionInfoPred> &miPredList
)
{
  int nWidth = pu.lumaSize().width;
  int nHeight = pu.lumaSize().height;
  if (refPairList.size() < 2 || !xAMLGetCurBlkTemplate(pu, nWidth, nHeight))
  {
    return;
  }

  PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
  PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
  PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
  PredictionUnit tmpPU = pu;

  DistParam cDistParam;
  cDistParam.applyWeight = false;
  Distortion uiCost;

  if (pu.cu->affine)
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    std::vector<Mv> cMvdCandList[2][3];
    for (int refList = 0; refList < 2; refList++)
    {
      cMvdCandList[refList][0].push_back(tmpPU.mvdAffi[refList][0]);
      cMvdCandList[refList][1].push_back(tmpPU.mvdAffi[refList][1]);
      cMvdCandList[refList][2].push_back(tmpPU.mvdAffi[refList][2]);
      if (pu.isMvsdApplicable() && (!refList || cMvdCandList[0][0].size() < 2))
      {
        deriveMVDcandAffine(tmpPU, RefPicList(refList), cMvdCandList[refList]);
      }
    }
#endif

    for (int idx = 0; idx < refPairList.size(); idx++)
    {
      tmpPU.refIdx[0] = refPairList[idx].refIdx[0];
      tmpPU.refIdx[1] = refPairList[idx].refIdx[1];

      AffineAMVPInfo affineAMVPInfo[2];
      for (int refList = 0; refList < 2; refList++)
      {
        RefPicList eRefList = (RefPicList)refList;
        PU::fillAffineMvpCand(tmpPU, eRefList, tmpPU.refIdx[eRefList], affineAMVPInfo[eRefList]);

#if !JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
        const unsigned mvp_idx = tmpPU.mvpIdx[eRefList];

        Mv mvLT = affineAMVPInfo[eRefList].mvCandLT[mvp_idx] + tmpPU.mvdAffi[eRefList][0];
        Mv mvRT = affineAMVPInfo[eRefList].mvCandRT[mvp_idx] + tmpPU.mvdAffi[eRefList][1];
        mvRT += tmpPU.mvdAffi[eRefList][0];

        Mv mvLB;
        if (tmpPU.cu->affineType == AFFINEMODEL_6PARAM)
        {
          mvLB = affineAMVPInfo[eRefList].mvCandLB[mvp_idx] + tmpPU.mvdAffi[eRefList][2];
          mvLB += tmpPU.mvdAffi[eRefList][0];
        }
 
        tmpPU.mvAffi[eRefList][0] = mvLT;
        tmpPU.mvAffi[eRefList][1] = mvRT;
        tmpPU.mvAffi[eRefList][2] = mvLB;
#endif
      }

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      refPairList[idx].cost = std::numeric_limits<Distortion>::max();
      for (size_t i = 0; i < cMvdCandList[0][0].size(); i++)
      {
        for (size_t j = 0; j < cMvdCandList[1][0].size(); j++)
        {
          MotionInfoPred miPred;
          miPred.interDir = 3;
          int refList = 0;
          unsigned mvp_idx = tmpPU.mvpIdx[refList];

          Mv mvLT = affineAMVPInfo[refList].mvCandLT[mvp_idx] + cMvdCandList[refList][0][i];
          Mv mvRT = affineAMVPInfo[refList].mvCandRT[mvp_idx] + cMvdCandList[refList][1][i];
          mvRT += cMvdCandList[refList][0][i];

          Mv mvLB;
          if (tmpPU.cu->affineType == AFFINEMODEL_6PARAM)
          {
            mvLB = affineAMVPInfo[refList].mvCandLB[mvp_idx] + cMvdCandList[refList][2][i];
            mvLB += cMvdCandList[refList][0][i];
          }

          tmpPU.mvAffi[refList][0] = mvLT;
          tmpPU.mvAffi[refList][1] = mvRT;
          tmpPU.mvAffi[refList][2] = mvLB;

          miPred.refIdx[refList] = tmpPU.refIdx[refList];
          miPred.mvdAffi[refList][0] = cMvdCandList[refList][0][i];
          miPred.mvdAffi[refList][1] = cMvdCandList[refList][1][i];
          miPred.mvdAffi[refList][2] = cMvdCandList[refList][2][i];
          miPred.mvAffi[refList][0] = tmpPU.mvAffi[refList][0];
          miPred.mvAffi[refList][1] = tmpPU.mvAffi[refList][1];
          miPred.mvAffi[refList][2] = tmpPU.mvAffi[refList][2];
          refList = 1;
          mvp_idx = tmpPU.mvpIdx[refList];

          mvLT = affineAMVPInfo[refList].mvCandLT[mvp_idx] + cMvdCandList[refList][0][j];
          mvRT = affineAMVPInfo[refList].mvCandRT[mvp_idx] + cMvdCandList[refList][1][j];
          mvRT += cMvdCandList[refList][0][j];

          if (tmpPU.cu->affineType == AFFINEMODEL_6PARAM)
          {
            mvLB = affineAMVPInfo[refList].mvCandLB[mvp_idx] + cMvdCandList[refList][2][j];
            mvLB += cMvdCandList[refList][0][j];
          }

          tmpPU.mvAffi[refList][0] = mvLT;
          tmpPU.mvAffi[refList][1] = mvRT;
          tmpPU.mvAffi[refList][2] = mvLB;

          miPred.refIdx[refList] = tmpPU.refIdx[refList];
          miPred.mvdAffi[refList][0] = cMvdCandList[refList][0][j];
          miPred.mvdAffi[refList][1] = cMvdCandList[refList][1][j];
          miPred.mvdAffi[refList][2] = cMvdCandList[refList][2][j];
          miPred.mvAffi[refList][0] = tmpPU.mvAffi[refList][0];
          miPred.mvAffi[refList][1] = tmpPU.mvAffi[refList][1];
          miPred.mvAffi[refList][2] = tmpPU.mvAffi[refList][2];

#endif
#if !JVET_Z0067_RPR_ENABLE
      getAffAMLRefTemplate(tmpPU, pcBufPredRefTop, pcBufPredRefLeft);
#endif

      uiCost = 0;
      bool bRefIsRescaled = false;
      for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
      {
        const RefPicList eRefPicList = refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
        bRefIsRescaled |= (tmpPU.refIdx[refList] >= 0) ? tmpPU.cu->slice->getRefPic(eRefPicList, tmpPU.refIdx[refList])->isRefScaled(tmpPU.cs->pps) : false;
      }
      if (bRefIsRescaled)
      {
        uiCost = std::numeric_limits<Distortion>::max();
      }
      else
      {
#if JVET_Z0067_RPR_ENABLE
        getAffAMLRefTemplate(tmpPU, pcBufPredRefTop, pcBufPredRefLeft);
#endif
        if (m_bAMLTemplateAvailabe[0])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }

        if (m_bAMLTemplateAvailabe[1])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }
      }
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
          miPred.cost = uiCost;
          miPredList.push_back(miPred);
          if (uiCost < refPairList[idx].cost)
          {
            refPairList[idx].cost = uiCost;
          }
        }
      }
#else
      refPairList[idx].cost = uiCost;
#endif
    }
  }
  else
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    std::vector<Mv> cMvdCandList[2];
    for (int refList = 0; refList < 2; refList++)
    {
      cMvdCandList[refList].push_back(tmpPU.mvd[refList]);
      if (pu.isMvsdApplicable() && (!refList || cMvdCandList[0].size() < 2))
      {
        deriveMVDcand(tmpPU, RefPicList(refList), cMvdCandList[refList]);
      }
    }
#endif
    for (int idx = 0; idx < refPairList.size(); idx++)
    {
      tmpPU.refIdx[0] = refPairList[idx].refIdx[0];
      tmpPU.refIdx[1] = refPairList[idx].refIdx[1];

      AMVPInfo amvpInfo[2];
      for (int refList = 0; refList < 2; refList++)
      {
        RefPicList eRefList = (RefPicList)refList;
        {
          PU::fillMvpCand(tmpPU, eRefList, tmpPU.refIdx[eRefList], amvpInfo[eRefList]
#if TM_AMVP
            , this
#endif
          );

#if !JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
          tmpPU.mv[eRefList] = amvpInfo[eRefList].mvCand[tmpPU.mvpIdx[eRefList]] + tmpPU.mvd[eRefList];
          tmpPU.mv[eRefList].mvCliptoStorageBitDepth();
#endif
        }
      }
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      refPairList[idx].cost = MAX_UINT;
      for (size_t i = 0; i < cMvdCandList[0].size(); i++)
      {
        for (size_t j = 0; j < cMvdCandList[1].size(); j++)
        {
          MotionInfoPred miPred;
          miPred.interDir = 3;
          tmpPU.mv[REF_PIC_LIST_0] = amvpInfo[REF_PIC_LIST_0].mvCand[tmpPU.mvpIdx[REF_PIC_LIST_0]] + cMvdCandList[REF_PIC_LIST_0][i];
          tmpPU.mv[REF_PIC_LIST_0].mvCliptoStorageBitDepth();

          miPred.refIdx[REF_PIC_LIST_0] = tmpPU.refIdx[REF_PIC_LIST_0];
          miPred.mvd[REF_PIC_LIST_0] = cMvdCandList[REF_PIC_LIST_0][i];
          miPred.mv[REF_PIC_LIST_0] = tmpPU.mv[REF_PIC_LIST_0];

          tmpPU.mv[REF_PIC_LIST_1] = amvpInfo[REF_PIC_LIST_1].mvCand[tmpPU.mvpIdx[REF_PIC_LIST_1]] + cMvdCandList[REF_PIC_LIST_1][j];
          tmpPU.mv[REF_PIC_LIST_1].mvCliptoStorageBitDepth();

          miPred.refIdx[REF_PIC_LIST_1] = tmpPU.refIdx[REF_PIC_LIST_1];
          miPred.mvd[REF_PIC_LIST_1] = cMvdCandList[REF_PIC_LIST_1][j];
          miPred.mv[REF_PIC_LIST_1] = tmpPU.mv[REF_PIC_LIST_1];
#endif

      uiCost = 0;
      bool bRefIsRescaled = false;
      for (uint32_t refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
      {
        const RefPicList eRefPicList = refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0;
        bRefIsRescaled |= (tmpPU.refIdx[refList] >= 0) ? tmpPU.cu->slice->getRefPic(eRefPicList, tmpPU.refIdx[refList])->isRefScaled(tmpPU.cs->pps) : false;
      }
      if (bRefIsRescaled)
      {
        uiCost = std::numeric_limits<Distortion>::max();
      }
      else
      {
        getBlkAMLRefTemplate(tmpPU, pcBufPredRefTop, pcBufPredRefLeft);
        if (m_bAMLTemplateAvailabe[0])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }

        if (m_bAMLTemplateAvailabe[1])
        {
          m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

          uiCost += cDistParam.distFunc(cDistParam);
        }
      }

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
          miPred.cost = uiCost;
          miPredList.push_back(miPred);
          if (uiCost < refPairList[idx].cost)
          {
            refPairList[idx].cost = uiCost;
          }
        }
      }
#else
      refPairList[idx].cost = uiCost;
#endif
    }
  }
  std::stable_sort(refPairList.begin(), refPairList.end(), [](const RefPicPair & l, const RefPicPair & r) {return l.cost < r.cost; });
}

#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if JVET_Z0054_BLK_REF_PIC_REORDER
void InterPrediction::deriveMVDcand(const PredictionUnit& pu, RefPicList eRefPicList, std::vector<Mv>& cMvdCandList)
{
  Mv cMvdKnownAtDecoder;
  cMvdKnownAtDecoder.set(pu.mvd[eRefPicList].getAbsHor(), pu.mvd[eRefPicList].getAbsVer());

  const static int patternsX[4][2] =
  {
    { +1, +1 },
    { +1, -1 },
  };

  const static int patternsY[4][2] =
  {
    { +1, +1 },
    { -1, +1 },
  };

  const static int patternsXY[4][2] =
  {
    { +1, +1 },
    { +1, -1 },
    { -1, +1 },
    { -1, -1 },
  };

  typedef int Int2[2];
  const Int2* patterns = 0;
  uint16_t patternsNum = 0;
  if (cMvdKnownAtDecoder.getHor() == 0 && cMvdKnownAtDecoder.getVer() == 0)
  {
    return;
  }
  if (cMvdKnownAtDecoder.getHor() == 0)
  {
    patterns = patternsX;
    patternsNum = 2;
  }
  else if (cMvdKnownAtDecoder.getVer() == 0)
  {
    patterns = patternsY;
    patternsNum = 2;
  }
  else
  {
    patterns = patternsXY;
    patternsNum = 4;
  }
  cMvdCandList.clear();
  for (int n = 0; n < patternsNum; ++n)
  {
    auto sign = patterns[n];
    auto cMv = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
    cMvdCandList.push_back(cMv);
  }
}

void InterPrediction::deriveMVDcandAffine(const PredictionUnit& pu, RefPicList eRefPicList, std::vector<Mv> cMvdCandList[3])
{
  Mv cMvdKnownAtDecoder[3];
  for (int i = 0; i < 3; i++)
  {
    cMvdKnownAtDecoder[i].set(pu.mvdAffi[eRefPicList][i].getAbsHor(), pu.mvdAffi[eRefPicList][i].getAbsVer());
  }

  const static int patterns2[2][1] =
  {
    { +1 },
    { -1 },
  };

  const static int patterns4[4][2] =
  {
    { +1, +1 },
    { +1, -1 },
    { -1, +1 },
    { -1, -1 },
  };

  const static int patterns8[8][3] =
  {
    { +1, +1, +1 },
    { +1, +1, -1 },
    { +1, -1, +1 },
    { +1, -1, -1 },
    { -1, +1, +1 },
    { -1, +1, -1 },
    { -1, -1, +1 },
    { -1, -1, -1 },
  };

  const static int patterns16[16][4] =
  {
    { +1, +1, +1, +1 },
    { +1, +1, +1, -1 },
    { +1, +1, -1, +1 },
    { +1, +1, -1, -1 },
    { +1, -1, +1, +1 },
    { +1, -1, +1, -1 },
    { +1, -1, -1, +1 },
    { +1, -1, -1, -1 },
    { -1, +1, +1, +1 },
    { -1, +1, +1, -1 },
    { -1, +1, -1, +1 },
    { -1, +1, -1, -1 },
    { -1, -1, +1, +1 },
    { -1, -1, +1, -1 },
    { -1, -1, -1, +1 },
    { -1, -1, -1, -1 },
  };

  const static int patterns32[32][5] =
  {
    { +1, +1, +1, +1, +1 },
    { +1, +1, +1, +1, -1 },
    { +1, +1, +1, -1, +1 },
    { +1, +1, +1, -1, -1 },
    { +1, +1, -1, +1, +1 },
    { +1, +1, -1, +1, -1 },
    { +1, +1, -1, -1, +1 },
    { +1, +1, -1, -1, -1 },
    { +1, -1, +1, +1, +1 },
    { +1, -1, +1, +1, -1 },
    { +1, -1, +1, -1, +1 },
    { +1, -1, +1, -1, -1 },
    { +1, -1, -1, +1, +1 },
    { +1, -1, -1, +1, -1 },
    { +1, -1, -1, -1, +1 },
    { +1, -1, -1, -1, -1 },
    { -1, +1, +1, +1, +1 },
    { -1, +1, +1, +1, -1 },
    { -1, +1, +1, -1, +1 },
    { -1, +1, +1, -1, -1 },
    { -1, +1, -1, +1, +1 },
    { -1, +1, -1, +1, -1 },
    { -1, +1, -1, -1, +1 },
    { -1, +1, -1, -1, -1 },
    { -1, -1, +1, +1, +1 },
    { -1, -1, +1, +1, -1 },
    { -1, -1, +1, -1, +1 },
    { -1, -1, +1, -1, -1 },
    { -1, -1, -1, +1, +1 },
    { -1, -1, -1, +1, -1 },
    { -1, -1, -1, -1, +1 },
    { -1, -1, -1, -1, -1 },
  };

  const static int patterns64[64][6] =
  {
    { +1, +1, +1, +1, +1, +1 },
    { +1, +1, +1, +1, +1, -1 },
    { +1, +1, +1, +1, -1, +1 },
    { +1, +1, +1, +1, -1, -1 },
    { +1, +1, +1, -1, +1, +1 },
    { +1, +1, +1, -1, +1, -1 },
    { +1, +1, +1, -1, -1, +1 },
    { +1, +1, +1, -1, -1, -1 },
    { +1, +1, -1, +1, +1, +1 },
    { +1, +1, -1, +1, +1, -1 },
    { +1, +1, -1, +1, -1, +1 },
    { +1, +1, -1, +1, -1, -1 },
    { +1, +1, -1, -1, +1, +1 },
    { +1, +1, -1, -1, +1, -1 },
    { +1, +1, -1, -1, -1, +1 },
    { +1, +1, -1, -1, -1, -1 },
    { +1, -1, +1, +1, +1, +1 },
    { +1, -1, +1, +1, +1, -1 },
    { +1, -1, +1, +1, -1, +1 },
    { +1, -1, +1, +1, -1, -1 },
    { +1, -1, +1, -1, +1, +1 },
    { +1, -1, +1, -1, +1, -1 },
    { +1, -1, +1, -1, -1, +1 },
    { +1, -1, +1, -1, -1, -1 },
    { +1, -1, -1, +1, +1, +1 },
    { +1, -1, -1, +1, +1, -1 },
    { +1, -1, -1, +1, -1, +1 },
    { +1, -1, -1, +1, -1, -1 },
    { +1, -1, -1, -1, +1, +1 },
    { +1, -1, -1, -1, +1, -1 },
    { +1, -1, -1, -1, -1, +1 },
    { +1, -1, -1, -1, -1, -1 },
    { -1, +1, +1, +1, +1, +1 },
    { -1, +1, +1, +1, +1, -1 },
    { -1, +1, +1, +1, -1, +1 },
    { -1, +1, +1, +1, -1, -1 },
    { -1, +1, +1, -1, +1, +1 },
    { -1, +1, +1, -1, +1, -1 },
    { -1, +1, +1, -1, -1, +1 },
    { -1, +1, +1, -1, -1, -1 },
    { -1, +1, -1, +1, +1, +1 },
    { -1, +1, -1, +1, +1, -1 },
    { -1, +1, -1, +1, -1, +1 },
    { -1, +1, -1, +1, -1, -1 },
    { -1, +1, -1, -1, +1, +1 },
    { -1, +1, -1, -1, +1, -1 },
    { -1, +1, -1, -1, -1, +1 },
    { -1, +1, -1, -1, -1, -1 },
    { -1, -1, +1, +1, +1, +1 },
    { -1, -1, +1, +1, +1, -1 },
    { -1, -1, +1, +1, -1, +1 },
    { -1, -1, +1, +1, -1, -1 },
    { -1, -1, +1, -1, +1, +1 },
    { -1, -1, +1, -1, +1, -1 },
    { -1, -1, +1, -1, -1, +1 },
    { -1, -1, +1, -1, -1, -1 },
    { -1, -1, -1, +1, +1, +1 },
    { -1, -1, -1, +1, +1, -1 },
    { -1, -1, -1, +1, -1, +1 },
    { -1, -1, -1, +1, -1, -1 },
    { -1, -1, -1, -1, +1, +1 },
    { -1, -1, -1, -1, +1, -1 },
    { -1, -1, -1, -1, -1, +1 },
    { -1, -1, -1, -1, -1, -1 },
  };
  std::vector<int> isZeroComp(6, 0);
  if (cMvdKnownAtDecoder[0].getHor() == 0)
  {
    isZeroComp[0] = 1;
  }
  if (cMvdKnownAtDecoder[0].getVer() == 0)
  {
    isZeroComp[1] = 1;
  }
  if (cMvdKnownAtDecoder[1].getHor() == 0)
  {
    isZeroComp[2] = 1;
  }
  if (cMvdKnownAtDecoder[1].getVer() == 0)
  {
    isZeroComp[3] = 1;
  }
  if (cMvdKnownAtDecoder[2].getHor() == 0)
  {
    isZeroComp[4] = 1;
  }
  if (cMvdKnownAtDecoder[2].getVer() == 0)
  {
    isZeroComp[5] = 1;
  }
  int nZeroComp = isZeroComp[0] + isZeroComp[1] + isZeroComp[2] + isZeroComp[3] + isZeroComp[4] + isZeroComp[5];
  if (nZeroComp == 6)
  {
    return;
  }

  uint16_t patternsNum = 0;
  if (nZeroComp == 0)
  {
    patternsNum = 64;
    cMvdCandList[0].resize(patternsNum);
    cMvdCandList[1].resize(patternsNum);
    cMvdCandList[2].resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto sign = patterns64[n];
      cMvdCandList[0][n] = Mv(sign[0] * cMvdKnownAtDecoder[0].getHor(), sign[1] * cMvdKnownAtDecoder[0].getVer());
      cMvdCandList[1][n] = Mv(sign[2] * cMvdKnownAtDecoder[1].getHor(), sign[3] * cMvdKnownAtDecoder[1].getVer());
      cMvdCandList[2][n] = Mv(sign[4] * cMvdKnownAtDecoder[2].getHor(), sign[5] * cMvdKnownAtDecoder[2].getVer());
    }
  }
  else if (nZeroComp == 1)
  {
    patternsNum = 32;
    cMvdCandList[0].resize(patternsNum);
    cMvdCandList[1].resize(patternsNum);
    cMvdCandList[2].resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto signR = patterns32[n];
      int sign[6];
      int k = 0;
      for (int i = 0; i < 6; i++)
      {
        if (isZeroComp[i])
        {
          sign[i] = +1;
        }
        else
        {
          sign[i] = signR[k];
          k++;
        }
      }
      cMvdCandList[0][n] = Mv(sign[0] * cMvdKnownAtDecoder[0].getHor(), sign[1] * cMvdKnownAtDecoder[0].getVer());
      cMvdCandList[1][n] = Mv(sign[2] * cMvdKnownAtDecoder[1].getHor(), sign[3] * cMvdKnownAtDecoder[1].getVer());
      cMvdCandList[2][n] = Mv(sign[4] * cMvdKnownAtDecoder[2].getHor(), sign[5] * cMvdKnownAtDecoder[2].getVer());
    }
  }
  else if (nZeroComp == 2)
  {
    patternsNum = 16;
    cMvdCandList[0].resize(patternsNum);
    cMvdCandList[1].resize(patternsNum);
    cMvdCandList[2].resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto signR = patterns16[n];
      int sign[6];
      int k = 0;
      for (int i = 0; i < 6; i++)
      {
        if (isZeroComp[i])
        {
          sign[i] = +1;
        }
        else
        {
          sign[i] = signR[k];
          k++;
        }
      }
      cMvdCandList[0][n] = Mv(sign[0] * cMvdKnownAtDecoder[0].getHor(), sign[1] * cMvdKnownAtDecoder[0].getVer());
      cMvdCandList[1][n] = Mv(sign[2] * cMvdKnownAtDecoder[1].getHor(), sign[3] * cMvdKnownAtDecoder[1].getVer());
      cMvdCandList[2][n] = Mv(sign[4] * cMvdKnownAtDecoder[2].getHor(), sign[5] * cMvdKnownAtDecoder[2].getVer());
    }
  }
  else if (nZeroComp == 3)
  {
    patternsNum = 8;
    cMvdCandList[0].resize(patternsNum);
    cMvdCandList[1].resize(patternsNum);
    cMvdCandList[2].resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto signR = patterns8[n];
      int sign[6];
      int k = 0;
      for (int i = 0; i < 6; i++)
      {
        if (isZeroComp[i])
        {
          sign[i] = +1;
        }
        else
        {
          sign[i] = signR[k];
          k++;
        }
      }
      cMvdCandList[0][n] = Mv(sign[0] * cMvdKnownAtDecoder[0].getHor(), sign[1] * cMvdKnownAtDecoder[0].getVer());
      cMvdCandList[1][n] = Mv(sign[2] * cMvdKnownAtDecoder[1].getHor(), sign[3] * cMvdKnownAtDecoder[1].getVer());
      cMvdCandList[2][n] = Mv(sign[4] * cMvdKnownAtDecoder[2].getHor(), sign[5] * cMvdKnownAtDecoder[2].getVer());
    }
  }
  else if (nZeroComp == 4)
  {
    patternsNum = 4;
    cMvdCandList[0].resize(patternsNum);
    cMvdCandList[1].resize(patternsNum);
    cMvdCandList[2].resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto signR = patterns4[n];
      int sign[6];
      int k = 0;
      for (int i = 0; i < 6; i++)
      {
        if (isZeroComp[i])
        {
          sign[i] = +1;
        }
        else
        {
          sign[i] = signR[k];
          k++;
        }
      }
      cMvdCandList[0][n] = Mv(sign[0] * cMvdKnownAtDecoder[0].getHor(), sign[1] * cMvdKnownAtDecoder[0].getVer());
      cMvdCandList[1][n] = Mv(sign[2] * cMvdKnownAtDecoder[1].getHor(), sign[3] * cMvdKnownAtDecoder[1].getVer());
      cMvdCandList[2][n] = Mv(sign[4] * cMvdKnownAtDecoder[2].getHor(), sign[5] * cMvdKnownAtDecoder[2].getVer());
    }
  }
  else
  {
    patternsNum = 2;
    cMvdCandList[0].resize(patternsNum);
    cMvdCandList[1].resize(patternsNum);
    cMvdCandList[2].resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto signR = patterns2[n];
      int sign[6];
      int k = 0;
      for (int i = 0; i < 6; i++)
      {
        if (isZeroComp[i])
        {
          sign[i] = +1;
        }
        else
        {
          sign[i] = signR[k];
          k++;
        }
      }
      cMvdCandList[0][n] = Mv(sign[0] * cMvdKnownAtDecoder[0].getHor(), sign[1] * cMvdKnownAtDecoder[0].getVer());
      cMvdCandList[1][n] = Mv(sign[2] * cMvdKnownAtDecoder[1].getHor(), sign[3] * cMvdKnownAtDecoder[1].getVer());
      cMvdCandList[2][n] = Mv(sign[4] * cMvdKnownAtDecoder[2].getHor(), sign[5] * cMvdKnownAtDecoder[2].getVer());
    }
  }
}
#endif

  void InterPrediction::deriveMvdSign(const Mv& cMvPred, const Mv& cMvdKnownAtDecoder, PredictionUnit& pu, RefPicList eRefList, int refIdx, std::vector<Mv>& cMvdDerived)
  {
#if JVET_Z0054_BLK_REF_PIC_REORDER
    deriveMVDcand(pu, eRefList, cMvdDerived);
#else
    const static int patternsX[4][2] =
    {
      { +1, +1 },
      { +1, -1 },
    };
    
    const static int patternsY[4][2] =
    {
      { +1, +1 },
      { -1, +1 },
    };
    
    const static int patternsXY[4][2] =
    {
      { +1, +1 },
      { +1, -1 },
      { -1, +1 },
      { -1, -1 },
    };
    
    typedef int Int2[2];
    const Int2* patterns = 0;
    uint16_t patternsNum = 0;
    if (cMvdKnownAtDecoder.getHor() == 0)
    {
      patterns = patternsX;
      patternsNum = 2;
    }
    else if (cMvdKnownAtDecoder.getVer() == 0)
    {
      patterns = patternsY;
      patternsNum = 2;
    }
    else
    {
      patterns = patternsXY;
      patternsNum = 4;
    }
    cMvdDerived.resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto sign = patterns[n];
      auto cMv = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
      cMvdDerived[n] = cMv;
    }
#endif
    if (!pu.lumaPos().x && !pu.lumaPos().y)
    {
      return;
    }
    CHECK(refIdx < 0, "Invalid reference index for FRUC");
    
    const Picture& refPic = *pu.cu->slice->getRefPic(eRefList, refIdx)->unscaledPic;
    InterPredResources interRes(m_pcReshape, m_pcRdCost, m_if, m_filteredBlockTmp[0][COMPONENT_Y]
                                , m_filteredBlock[3][1][0], m_filteredBlock[3][0][0]
                                );
    
    TplMatchingCtrl tplCtrl(pu, interRes, refPic, true, COMPONENT_Y, true, 0, m_pcCurTplAbove, m_pcCurTplLeft, m_pcRefTplAbove, m_pcRefTplLeft, Mv(0, 0), nullptr, 0);
    
#if JVET_Z0054_BLK_REF_PIC_REORDER
    size_t patternsNum = cMvdDerived.size();
#endif
    std::vector<std::pair<Mv, Distortion>> aMvCostVec(patternsNum);
    
#if JVET_Z0067_RPR_ENABLE
    bool  bIsRefScaled = pu.cu->slice->getRefPic(eRefList, refIdx)->isRefScaled( pu.cs->pps );
    for (int n = 0; n < patternsNum; ++n)
    {
      auto cMvdTest = cMvdDerived[n];
      Mv cMvTest = cMvPred + cMvdTest;
      Distortion uiCost = bIsRefScaled ? std::numeric_limits<Distortion>::max() : tplCtrl.xGetTempMatchError<TM_TPL_SIZE>(cMvTest);
      aMvCostVec[n] = { cMvdTest, uiCost };
    }
#else
    for (int n = 0; n < patternsNum; ++n)
    {
      auto cMvdTest = cMvdDerived[n];
      Mv cMvTest = cMvPred + cMvdTest;
      Distortion uiCost = tplCtrl.xGetTempMatchError<TM_TPL_SIZE>(cMvTest);
      aMvCostVec[n] = { cMvdTest, uiCost };
    }
#endif

    std::stable_sort(aMvCostVec.begin(), aMvCostVec.end(), [](const std::pair<Mv, Distortion> & l, const std::pair<Mv, Distortion> & r) {return l.second < r.second; });
    
    for (int n = 0; n < patternsNum; ++n)
    {
      cMvdDerived[n] = aMvCostVec[n].first;
    }
  }
  void InterPrediction::deriveMvdSignSMVD(const Mv& cMvPred, const Mv& cMvPred2, const Mv& cMvdKnownAtDecoder, PredictionUnit& pu, std::vector<Mv>& cMvdDerived)
  {
#if JVET_Z0054_BLK_REF_PIC_REORDER
    deriveMVDcand(pu, REF_PIC_LIST_0, cMvdDerived);
    int nWidth = pu.lumaSize().width;
    int nHeight = pu.lumaSize().height;
    if (cMvdDerived.size() < 2 || !xAMLGetCurBlkTemplate(pu, nWidth, nHeight))
    {
      return;
    }
    std::vector<std::pair<Mv, Distortion>> aMvCostVec;
    PelUnitBuf pcBufPredRefTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredCurTop = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], nWidth, AML_MERGE_TEMPLATE_SIZE)));
    PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvRefAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
    PelUnitBuf pcBufPredCurLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], AML_MERGE_TEMPLATE_SIZE, nHeight)));
    PredictionUnit tmpPU = pu;
    AMVPInfo amvpInfo;
    PU::fillMvpCand(tmpPU, REF_PIC_LIST_0, tmpPU.refIdx[REF_PIC_LIST_0], amvpInfo
#if TM_AMVP
      , this
#endif
    );
    AMVPInfo amvpInfo1;
    PU::fillMvpCand(tmpPU, REF_PIC_LIST_1, tmpPU.refIdx[REF_PIC_LIST_1], amvpInfo1
#if TM_AMVP
      , this
#endif
    );


    DistParam cDistParam;
    cDistParam.applyWeight = false;
    Distortion uiCost;
    for (std::vector<Mv>::iterator it = cMvdDerived.begin(); it != cMvdDerived.end(); ++it)
    {
      tmpPU.mvd[0] = *it;
      tmpPU.mv[0] = amvpInfo.mvCand[tmpPU.mvpIdx[0]] + tmpPU.mvd[0];
      tmpPU.mv[0].mvCliptoStorageBitDepth();
      tmpPU.mv[1] = amvpInfo.mvCand[tmpPU.mvpIdx[1]] - tmpPU.mvd[0];
      tmpPU.mv[1].mvCliptoStorageBitDepth();
      getBlkAMLRefTemplate(tmpPU, pcBufPredRefTop, pcBufPredRefLeft);

      uiCost = 0;
      if (m_bAMLTemplateAvailabe[0])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurTop.Y(), pcBufPredRefTop.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }

      if (m_bAMLTemplateAvailabe[1])
      {
        m_pcRdCost->setDistParam(cDistParam, pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, false);

        uiCost += cDistParam.distFunc(cDistParam);
      }

      aMvCostVec.push_back(std::pair<Mv, Distortion>(*it, uiCost));
    }
    std::stable_sort(aMvCostVec.begin(), aMvCostVec.end(), [](const std::pair<Mv, Distortion> & l, const std::pair<Mv, Distortion> & r) {return l.second < r.second; });
    for (size_t n = 0; n < cMvdDerived.size(); ++n)
    {
      cMvdDerived[n] = aMvCostVec[n].first;
    }
#else
    const static int patternsX[4][2] =
    {
      { +1, +1 },
      { +1, -1 },
    };
    
    const static int patternsY[4][2] =
    {
      { +1, +1 },
      { -1, +1 },
    };
    
    const static int patternsXY[4][2] =
    {
      { +1, +1 },
      { +1, -1 },
      { -1, +1 },
      { -1, -1 },
    };
    
    typedef int Int2[2];
    const Int2* patterns = 0;
    uint16_t patternsNum = 0;
    if (cMvdKnownAtDecoder.getHor() == 0)
    {
      patterns = patternsX;
      patternsNum = 2;
    }
    else if (cMvdKnownAtDecoder.getVer() == 0)
    {
      patterns = patternsY;
      patternsNum = 2;
    }
    else
    {
      patterns = patternsXY;
      patternsNum = 4;
    }
    
    cMvdDerived.resize(patternsNum);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto sign = patterns[n];
      auto cMv = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
      cMvdDerived[n] = cMv;
    }
    
    if (!pu.lumaPos().x && !pu.lumaPos().y)
    {
      return;
    }
    std::vector<std::pair<Mv, Distortion>> aMvCostVec(patternsNum);
    InterPredResources interRes(m_pcReshape, m_pcRdCost, m_if, m_filteredBlockTmp[0][COMPONENT_Y]
                                , m_filteredBlock[3][1][0], m_filteredBlock[3][0][0]
                                );
    
    // For L0
    int refIdx = pu.cs->slice->getSymRefIdx(REF_PIC_LIST_0);
    CHECK(refIdx < 0, "Invalid reference index for SMVD L0");
    const Picture& refPic = *pu.cu->slice->getRefPic(REF_PIC_LIST_0, refIdx)->unscaledPic;
    
    TplMatchingCtrl tplCtrl(pu, interRes, refPic, true, COMPONENT_Y, true, 0, m_pcCurTplAbove, m_pcCurTplLeft, m_pcRefTplAbove, m_pcRefTplLeft, Mv(0, 0), nullptr, 0);
    
#if JVET_Z0067_RPR_ENABLE
    bool  bIsRefScaled = pu.cu->slice->getRefPic(REF_PIC_LIST_0, refIdx)->isRefScaled(pu.cs->pps);
    for (int n = 0; n < patternsNum; ++n)
    {
      auto cMvdTest = cMvdDerived[n];
      Mv cMvTest = cMvPred + cMvdTest;
      Distortion uiCost = bIsRefScaled ? std::numeric_limits<Distortion>::max() : tplCtrl.xGetTempMatchError<TM_TPL_SIZE>(cMvTest);
      aMvCostVec[n] = { cMvdTest, uiCost };
    }
#else
    for (int n = 0; n < patternsNum; ++n)
    {
      auto cMvdTest = cMvdDerived[n];
      Mv cMvTest = cMvPred + cMvdTest;
      Distortion uiCost = tplCtrl.xGetTempMatchError<TM_TPL_SIZE>(cMvTest);
      aMvCostVec[n] = { cMvdTest, uiCost };
    }
#endif 

    std::stable_sort(aMvCostVec.begin(), aMvCostVec.end(), [](const std::pair<Mv, Distortion> & l, const std::pair<Mv, Distortion> & r) {return l.second < r.second; });
    
    for (int n = 0; n < patternsNum; ++n)
    {
      cMvdDerived[n] = aMvCostVec[n].first;
    }
#endif
  }

  void InterPrediction::deriveMvdSignAffine(const Mv& cMvPred, const Mv& cMvPred2, const Mv& cMvPred3,
#if JVET_Z0054_BLK_REF_PIC_REORDER
                                            const Mv cMvdKnownAtDecoder[3],
#else
                                            const Mv& cMvdKnownAtDecoder, const Mv& cMvdKnownAtDecoder2, const Mv& cMvdKnownAtDecoder3,
#endif
                                            PredictionUnit& pu, RefPicList eRefList, int refIdx, std::vector<Mv>& cMvdDerived, std::vector<Mv>& cMvdDerived2, std::vector<Mv>& cMvdDerived3)
  {
#if JVET_Z0054_BLK_REF_PIC_REORDER
    std::vector<Mv> MvdCand[3];
    deriveMVDcandAffine(pu, eRefList, MvdCand);
    size_t patternsNum = MvdCand[0].size();
#else
    int patterns2[2][1] =
    {
      { +1 },
      { -1 },
    };
    
    int patterns4[4][2] =
    {
      { +1, +1 },
      { +1, -1 },
      { -1, +1 },
      { -1, -1 },
    };
    
    int patterns8[8][3] =
    {
      {+1, +1, +1 },
      {+1, +1, -1 },
      {+1, -1, +1 },
      {+1, -1, -1 },
      {-1, +1, +1 },
      {-1, +1, -1 },
      {-1, -1, +1 },
      {-1, -1, -1 },
    };
    
    int patterns16[16][4] =
    {
      {+1, +1, +1, +1 },
      {+1, +1, +1, -1 },
      {+1, +1, -1, +1 },
      {+1, +1, -1, -1 },
      {+1, -1, +1, +1 },
      {+1, -1, +1, -1 },
      {+1, -1, -1, +1 },
      {+1, -1, -1, -1 },
      {-1, +1, +1, +1 },
      {-1, +1, +1, -1 },
      {-1, +1, -1, +1 },
      {-1, +1, -1, -1 },
      {-1, -1, +1, +1 },
      {-1, -1, +1, -1 },
      {-1, -1, -1, +1 },
      {-1, -1, -1, -1 },
    };
    
    int patterns32[32][5] =
    {
      {+1, +1, +1, +1, +1 },
      {+1, +1, +1, +1, -1 },
      {+1, +1, +1, -1, +1 },
      {+1, +1, +1, -1, -1 },
      {+1, +1, -1, +1, +1 },
      {+1, +1, -1, +1, -1 },
      {+1, +1, -1, -1, +1 },
      {+1, +1, -1, -1, -1 },
      {+1, -1, +1, +1, +1 },
      {+1, -1, +1, +1, -1 },
      {+1, -1, +1, -1, +1 },
      {+1, -1, +1, -1, -1 },
      {+1, -1, -1, +1, +1 },
      {+1, -1, -1, +1, -1 },
      {+1, -1, -1, -1, +1 },
      {+1, -1, -1, -1, -1 },
      {-1, +1, +1, +1, +1 },
      {-1, +1, +1, +1, -1 },
      {-1, +1, +1, -1, +1 },
      {-1, +1, +1, -1, -1 },
      {-1, +1, -1, +1, +1 },
      {-1, +1, -1, +1, -1 },
      {-1, +1, -1, -1, +1 },
      {-1, +1, -1, -1, -1 },
      {-1, -1, +1, +1, +1 },
      {-1, -1, +1, +1, -1 },
      {-1, -1, +1, -1, +1 },
      {-1, -1, +1, -1, -1 },
      {-1, -1, -1, +1, +1 },
      {-1, -1, -1, +1, -1 },
      {-1, -1, -1, -1, +1 },
      {-1, -1, -1, -1, -1 },
    };
    
    int patterns64[64][6] =
    {
      {+1, +1, +1, +1, +1, +1 },
      {+1, +1, +1, +1, +1, -1 },
      {+1, +1, +1, +1, -1, +1 },
      {+1, +1, +1, +1, -1, -1 },
      {+1, +1, +1, -1, +1, +1 },
      {+1, +1, +1, -1, +1, -1 },
      {+1, +1, +1, -1, -1, +1 },
      {+1, +1, +1, -1, -1, -1 },
      {+1, +1, -1, +1, +1, +1 },
      {+1, +1, -1, +1, +1, -1 },
      {+1, +1, -1, +1, -1, +1 },
      {+1, +1, -1, +1, -1, -1 },
      {+1, +1, -1, -1, +1, +1 },
      {+1, +1, -1, -1, +1, -1 },
      {+1, +1, -1, -1, -1, +1 },
      {+1, +1, -1, -1, -1, -1 },
      {+1, -1, +1, +1, +1, +1 },
      {+1, -1, +1, +1, +1, -1 },
      {+1, -1, +1, +1, -1, +1 },
      {+1, -1, +1, +1, -1, -1 },
      {+1, -1, +1, -1, +1, +1 },
      {+1, -1, +1, -1, +1, -1 },
      {+1, -1, +1, -1, -1, +1 },
      {+1, -1, +1, -1, -1, -1 },
      {+1, -1, -1, +1, +1, +1 },
      {+1, -1, -1, +1, +1, -1 },
      {+1, -1, -1, +1, -1, +1 },
      {+1, -1, -1, +1, -1, -1 },
      {+1, -1, -1, -1, +1, +1 },
      {+1, -1, -1, -1, +1, -1 },
      {+1, -1, -1, -1, -1, +1 },
      {+1, -1, -1, -1, -1, -1 },
      {-1, +1, +1, +1, +1, +1 },
      {-1, +1, +1, +1, +1, -1 },
      {-1, +1, +1, +1, -1, +1 },
      {-1, +1, +1, +1, -1, -1 },
      {-1, +1, +1, -1, +1, +1 },
      {-1, +1, +1, -1, +1, -1 },
      {-1, +1, +1, -1, -1, +1 },
      {-1, +1, +1, -1, -1, -1 },
      {-1, +1, -1, +1, +1, +1 },
      {-1, +1, -1, +1, +1, -1 },
      {-1, +1, -1, +1, -1, +1 },
      {-1, +1, -1, +1, -1, -1 },
      {-1, +1, -1, -1, +1, +1 },
      {-1, +1, -1, -1, +1, -1 },
      {-1, +1, -1, -1, -1, +1 },
      {-1, +1, -1, -1, -1, -1 },
      {-1, -1, +1, +1, +1, +1 },
      {-1, -1, +1, +1, +1, -1 },
      {-1, -1, +1, +1, -1, +1 },
      {-1, -1, +1, +1, -1, -1 },
      {-1, -1, +1, -1, +1, +1 },
      {-1, -1, +1, -1, +1, -1 },
      {-1, -1, +1, -1, -1, +1 },
      {-1, -1, +1, -1, -1, -1 },
      {-1, -1, -1, +1, +1, +1 },
      {-1, -1, -1, +1, +1, -1 },
      {-1, -1, -1, +1, -1, +1 },
      {-1, -1, -1, +1, -1, -1 },
      {-1, -1, -1, -1, +1, +1 },
      {-1, -1, -1, -1, +1, -1 },
      {-1, -1, -1, -1, -1, +1 },
      {-1, -1, -1, -1, -1, -1 },
    };
    std::vector<int> isZeroComp(6, 0);
    if (cMvdKnownAtDecoder.getHor() == 0)
    {
      isZeroComp[0] = 1;
    }
    if (cMvdKnownAtDecoder.getVer() == 0)
    {
      isZeroComp[1] = 1;
    }
    if (cMvdKnownAtDecoder2.getHor() == 0)
    {
      isZeroComp[2] = 1;
    }
    if (cMvdKnownAtDecoder2.getVer() == 0)
    {
      isZeroComp[3] = 1;
    }
    if (cMvdKnownAtDecoder3.getHor() == 0)
    {
      isZeroComp[4] = 1;
    }
    if (cMvdKnownAtDecoder3.getVer() == 0)
    {
      isZeroComp[5] = 1;
    }
    int nZeroComp = isZeroComp[0] + isZeroComp[1] + isZeroComp[2] + isZeroComp[3] + isZeroComp[4] + isZeroComp[5];
    CHECK(nZeroComp == 6, "nnZeroComp == 6");
    
    uint16_t patternsNum = 0;
    std::vector<Mv> MvdCand[3];
    if (nZeroComp == 0)
    {
      patternsNum = 64;
      MvdCand[0].resize(patternsNum);
      MvdCand[1].resize(patternsNum);
      MvdCand[2].resize(patternsNum);
      for (int n = 0; n < patternsNum; ++n)
      {
        auto sign = patterns64[n];
        MvdCand[0][n] = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
        MvdCand[1][n] = Mv(sign[2] * cMvdKnownAtDecoder2.getHor(), sign[3] * cMvdKnownAtDecoder2.getVer());
        MvdCand[2][n] = Mv(sign[4] * cMvdKnownAtDecoder3.getHor(), sign[5] * cMvdKnownAtDecoder3.getVer());
      }
    }
    else if (nZeroComp == 1)
    {
      patternsNum = 32;
      MvdCand[0].resize(patternsNum);
      MvdCand[1].resize(patternsNum);
      MvdCand[2].resize(patternsNum);
      for (int n = 0; n < patternsNum; ++n)
      {
        auto signR = patterns32[n];
        int sign[6];
        int k = 0;
        for (int i = 0; i < 6; i++)
        {
          if (isZeroComp[i])
          {
            sign[i] = +1;
          }
          else
          {
            sign[i] = signR[k];
            k++;
          }
        }
        MvdCand[0][n] = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
        MvdCand[1][n] = Mv(sign[2] * cMvdKnownAtDecoder2.getHor(), sign[3] * cMvdKnownAtDecoder2.getVer());
        MvdCand[2][n] = Mv(sign[4] * cMvdKnownAtDecoder3.getHor(), sign[5] * cMvdKnownAtDecoder3.getVer());
      }
    }
    else if (nZeroComp == 2)
    {
      patternsNum = 16;
      MvdCand[0].resize(patternsNum);
      MvdCand[1].resize(patternsNum);
      MvdCand[2].resize(patternsNum);
      for (int n = 0; n < patternsNum; ++n)
      {
        auto signR = patterns16[n];
        int sign[6];
        int k = 0;
        for (int i = 0; i < 6; i++)
        {
          if (isZeroComp[i])
          {
            sign[i] = +1;
          }
          else
          {
            sign[i] = signR[k];
            k++;
          }
        }
        MvdCand[0][n] = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
        MvdCand[1][n] = Mv(sign[2] * cMvdKnownAtDecoder2.getHor(), sign[3] * cMvdKnownAtDecoder2.getVer());
        MvdCand[2][n] = Mv(sign[4] * cMvdKnownAtDecoder3.getHor(), sign[5] * cMvdKnownAtDecoder3.getVer());
      }
    }
    else if (nZeroComp == 3)
    {
      patternsNum = 8;
      MvdCand[0].resize(patternsNum);
      MvdCand[1].resize(patternsNum);
      MvdCand[2].resize(patternsNum);
      for (int n = 0; n < patternsNum; ++n)
      {
        auto signR = patterns8[n];
        int sign[6];
        int k = 0;
        for (int i = 0; i < 6; i++)
        {
          if (isZeroComp[i])
          {
            sign[i] = +1;
          }
          else
          {
            sign[i] = signR[k];
            k++;
          }
        }
        MvdCand[0][n] = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
        MvdCand[1][n] = Mv(sign[2] * cMvdKnownAtDecoder2.getHor(), sign[3] * cMvdKnownAtDecoder2.getVer());
        MvdCand[2][n] = Mv(sign[4] * cMvdKnownAtDecoder3.getHor(), sign[5] * cMvdKnownAtDecoder3.getVer());
      }
    }
    else if (nZeroComp == 4)
    {
      patternsNum = 4;
      MvdCand[0].resize(patternsNum);
      MvdCand[1].resize(patternsNum);
      MvdCand[2].resize(patternsNum);
      for (int n = 0; n < patternsNum; ++n)
      {
        auto signR = patterns4[n];
        int sign[6];
        int k = 0;
        for (int i = 0; i < 6; i++)
        {
          if (isZeroComp[i])
          {
            sign[i] = +1;
          }
          else
          {
            sign[i] = signR[k];
            k++;
          }
        }
        MvdCand[0][n] = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
        MvdCand[1][n] = Mv(sign[2] * cMvdKnownAtDecoder2.getHor(), sign[3] * cMvdKnownAtDecoder2.getVer());
        MvdCand[2][n] = Mv(sign[4] * cMvdKnownAtDecoder3.getHor(), sign[5] * cMvdKnownAtDecoder3.getVer());
      }
    }
    else
    {
      patternsNum = 2;
      MvdCand[0].resize(patternsNum);
      MvdCand[1].resize(patternsNum);
      MvdCand[2].resize(patternsNum);
      for (int n = 0; n < patternsNum; ++n)
      {
        auto signR = patterns2[n];
        int sign[6];
        int k = 0;
        for (int i = 0; i < 6; i++)
        {
          if (isZeroComp[i])
          {
            sign[i] = +1;
          }
          else
          {
            sign[i] = signR[k];
            k++;
          }
        }
        MvdCand[0][n] = Mv(sign[0] * cMvdKnownAtDecoder.getHor(), sign[1] * cMvdKnownAtDecoder.getVer());
        MvdCand[1][n] = Mv(sign[2] * cMvdKnownAtDecoder2.getHor(), sign[3] * cMvdKnownAtDecoder2.getVer());
        MvdCand[2][n] = Mv(sign[4] * cMvdKnownAtDecoder3.getHor(), sign[5] * cMvdKnownAtDecoder3.getVer());
      }
    }
#endif
    
    cMvdDerived.resize(patternsNum);
    cMvdDerived2.resize(patternsNum);
    cMvdDerived3.resize(patternsNum);
    
    for (int n = 0; n < patternsNum; ++n)
    {
      cMvdDerived[n] = MvdCand[0][n];
      cMvdDerived2[n] = MvdCand[1][n];
      cMvdDerived3[n] = MvdCand[2][n];
    }
    
    if (!pu.lumaPos().x && !pu.lumaPos().y)
    {
      for (int n = 0; n < patternsNum; ++n)
      {
        cMvdDerived[n] = MvdCand[0][n];
        cMvdDerived2[n] = MvdCand[1][n];
        cMvdDerived3[n] = MvdCand[2][n];
      }
      return;
    }
    
    /////////////////////////////////////////////////////////////////
    Pel* refLeftTemplate = m_pcLICRefLeftTemplate;
    Pel* refAboveTemplate = m_pcLICRefAboveTemplate;
    Pel* recLeftTemplate = m_pcLICRecLeftTemplate;
    Pel* recAboveTemplate = m_pcLICRecAboveTemplate;
    int numTemplate[2] = { 0 , 0 }; // 0:Above, 1:Left
    
    const int width = pu.Y().width;
    const int height = pu.Y().height;
    int blockWidth = AFFINE_MIN_BLOCK_SIZE;
    int blockHeight = AFFINE_MIN_BLOCK_SIZE;
    
    const int iHalfBW = blockWidth >> 1;
    const int iHalfBH = blockHeight >> 1;
    
    const int iBit = MAX_CU_DEPTH;
    const int shift = iBit - 4 + MV_FRACTIONAL_BITS_INTERNAL;
    int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;
    
    CHECK(refIdx < 0, "Invalid reference index for FRUC");
    const Picture& refPic = *pu.cu->slice->getRefPic(eRefList, refIdx)->unscaledPic;
    std::vector<std::pair<int, Distortion>> aMvCostVec(patternsNum);
    Distortion uiCost = 0;
    
#if JVET_Z0067_RPR_ENABLE
    bool bIsRefScaled = pu.cu->slice->getRefPic(eRefList, refIdx)->isRefScaled( pu.cs->pps );
    if ( bIsRefScaled )
    {
      for (int n = 0; n < patternsNum; ++n)
      {
        aMvCostVec[n] = { n, uiCost };
      }
    }
    else
    {
#endif
    for (int n = 0; n < patternsNum; ++n)
    {
      uiCost = 0;
      //--------------------- (derive CPMVs)----------------------------------------------//
      Mv mvLT = cMvPred + MvdCand[0][n];
      Mv mvRT = cMvPred2 + MvdCand[1][n];
      mvRT += MvdCand[0][n];
      
      Mv mvLB;
      if (pu.cu->affineType == AFFINEMODEL_6PARAM)
      {
        mvLB = cMvPred3 + MvdCand[2][n];
        mvLB += MvdCand[0][n];
      }
      //--------------- Calculate dMVs ------------------------------------------//
      iDMvHorX = (mvRT - mvLT).getHor() << (iBit - floorLog2(width));
      iDMvHorY = (mvRT - mvLT).getVer() << (iBit - floorLog2(width));
      if (pu.cu->affineType == AFFINEMODEL_6PARAM)
      {
        iDMvVerX = (mvLB - mvLT).getHor() << (iBit - floorLog2(height));
        iDMvVerY = (mvLB - mvLT).getVer() << (iBit - floorLog2(height));
      }
      else
      {
        iDMvVerX = -iDMvHorY;
        iDMvVerY = iDMvHorX;
      }
      
      int iMvScaleHor = mvLT.getHor() << iBit;
      int iMvScaleVer = mvLT.getVer() << iBit;
      
      int mvScaleHorLine = iMvScaleHor + iDMvHorX * iHalfBW + iDMvVerX * iHalfBH;
      int mvScaleVerLine = iMvScaleVer + iDMvHorY * iHalfBW + iDMvVerY * iHalfBH;
      
      int deltaMvHorXBlk = iDMvHorX * blockWidth;
      int deltaMvHorYBlk = iDMvHorY * blockWidth;
      
      // get prediction block by block
      
      for (int h = 0; h < height; h += blockHeight)
      {
        int mvScaleHorBlk = mvScaleHorLine;
        int mvScaleVerBlk = mvScaleVerLine;
        
        for (int w = 0; w < width; w += blockWidth)
        {
          if (w != 0 && h != 0) continue; //applies only on boundary subblocks.
          int iMvScaleTmpHor, iMvScaleTmpVer;
          
          iMvScaleTmpHor = mvScaleHorBlk;
          iMvScaleTmpVer = mvScaleVerBlk;
          
          mvScaleHorBlk += deltaMvHorXBlk;
          mvScaleVerBlk += deltaMvHorYBlk;
          
          roundAffineMv(iMvScaleTmpHor, iMvScaleTmpVer, shift);
          Mv tmpMv(iMvScaleTmpHor, iMvScaleTmpVer);
          tmpMv.clipToStorageBitDepth();
          //clip
          clipMv(tmpMv, pu.lumaPos(), pu.lumaSize(), *pu.cs->sps, *pu.cs->pps);
          iMvScaleTmpHor = tmpMv.getHor();
          iMvScaleTmpVer = tmpMv.getVer();
          uiCost += xGetSublkTemplateCost(*pu.cu, COMPONENT_Y, refPic, Mv(iMvScaleTmpHor, iMvScaleTmpVer), blockWidth, blockHeight, w, h, numTemplate, refLeftTemplate, refAboveTemplate, recLeftTemplate, recAboveTemplate);
        }
      }
      aMvCostVec[n] = { n, uiCost };
    }
#if JVET_Z0067_RPR_ENABLE
    }
#endif
    //--------------------------------------------------------------------------------//
    /////////////////////////////////////////////////////////////////
    std::stable_sort(aMvCostVec.begin(), aMvCostVec.end(), [](const std::pair<int, Distortion> & l, const std::pair<int, Distortion> & r) {return l.second < r.second; });
    for (int n = 0; n < patternsNum; ++n)
    {
      int index = aMvCostVec[n].first;
      cMvdDerived[n] = MvdCand[0][index];
      cMvdDerived2[n] = MvdCand[1][index];
      cMvdDerived3[n] = MvdCand[2][index];
    }
  }
  
  Distortion InterPrediction::xGetSublkTemplateCost(const CodingUnit& cu,
                                                    const ComponentID compID,
                                                    const Picture&    refPic,
                                                    const Mv&         mv,
                                                    const int         sublkWidth,
                                                    const int         sublkHeight,
                                                    const int         posW,
                                                    const int         posH,
                                                    int*              numTemplate,
                                                    Pel*              refLeftTemplate,
                                                    Pel*              refAboveTemplate,
                                                    Pel*              recLeftTemplate,
                                                    Pel*              recAboveTemplate)
  {
#if JVET_Z0067_RPR_ENABLE
    CHECK(refPic.isRefScaled(cu.cs->pps), "xGetSublkTemplateCost ref Scaled not supported");
#endif
    const int       bitDepth = cu.cs->sps->getBitDepth(toChannelType(compID));
    const int       precShift = std::max(0, bitDepth - 12);
    Distortion cost = 0;
    
    const Picture&  currPic = *cu.cs->picture;
    const CodingUnit* const cuAbove = cu.cs->getCU(cu.blocks[compID].pos().offset(0, -1), toChannelType(compID));
    const CodingUnit* const cuLeft = cu.cs->getCU(cu.blocks[compID].pos().offset(-1, 0), toChannelType(compID));
    const CPelBuf recBuf = cuAbove || cuLeft ? currPic.getRecoBuf(cu.cs->picture->blocks[compID]) : CPelBuf();
    const CPelBuf refBuf = cuAbove || cuLeft ? refPic.getRecoBuf(refPic.blocks[compID]) : CPelBuf();
    
    std::vector<Pel>& invLUT = m_pcReshape->getInvLUT();
    
    // above
    if (cuAbove && posH == 0)
    {
      xGetPredBlkTpl<true>(cu, compID, refBuf, mv, posW, posH, sublkWidth, refAboveTemplate);
      const Pel*    rec = recBuf.bufAt(cu.blocks[compID].pos().offset(0, -1));
      for (int k = posW; k < posW + sublkWidth; k++)
      {
        int refVal = refAboveTemplate[k];
        int recVal = rec[k];
        
        if (isLuma(compID) && cu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
        {
          recVal = invLUT[recVal];
        }
        
        recVal >>= precShift;
        refVal >>= precShift;
        
        refAboveTemplate[k] = refVal;
        recAboveTemplate[k] = recVal;
        numTemplate[0]++;
        cost += (Distortion)(refVal - recVal) * (refVal - recVal);
      }
    }
    
    // left
    if (cuLeft && posW == 0)
    {
      xGetPredBlkTpl<false>(cu, compID, refBuf, mv, posW, posH, sublkHeight, refLeftTemplate);
      const Pel*    rec = recBuf.bufAt(cu.blocks[compID].pos().offset(-1, 0));
      for (int k = posH; k < posH + sublkHeight; k++)
      {
        int refVal = refLeftTemplate[k];
        int recVal = rec[recBuf.stride * k];
        
        if (isLuma(compID) && cu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
        {
          recVal = invLUT[recVal];
        }
        
        recVal >>= precShift;
        refVal >>= precShift;
        
        refLeftTemplate[k] = refVal;
        recLeftTemplate[k] = recVal;
        numTemplate[1]++;
        cost += (Distortion)(refVal - recVal) * (refVal - recVal);
      }
    }
    return cost;
  }
  int InterPrediction::deriveMVSDIdxFromMVDAffine(PredictionUnit& pu, RefPicList eRefList, std::vector<Mv>& cMvdDerived, std::vector<Mv>& cMvdDerived2, std::vector<Mv>& cMvdDerived3)
  {
    int mvsdIdx = 0;
    int shift = 0;
    int bin = 0;
    if (pu.mvdAffi[eRefList][0].getHor())
    {
      bin = (cMvdDerived[0].getHor() == pu.mvdAffi[eRefList][0].getHor()) ? 0 : 1;
      mvsdIdx += bin << shift;
      shift++;
    }
    if (pu.mvdAffi[eRefList][0].getVer())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i].getHor() == pu.mvdAffi[eRefList][0].getHor())
        {
          bin = (cMvdDerived[i].getVer() == pu.mvdAffi[eRefList][0].getVer()) ? 0 : 1;
          mvsdIdx += bin << shift;
          shift++;
          break;
        }
      }
    }
    if (pu.mvdAffi[eRefList][1].getHor())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i] == pu.mvdAffi[eRefList][0])
        {
          bin = (cMvdDerived2[i].getHor() == pu.mvdAffi[eRefList][1].getHor()) ? 0 : 1;
          mvsdIdx += bin << shift;
          shift++;
          break;
        }
      }
    }
    if (pu.mvdAffi[eRefList][1].getVer())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i] == pu.mvdAffi[eRefList][0] && cMvdDerived2[i].getHor() == pu.mvdAffi[eRefList][1].getHor())
        {
          bin = (cMvdDerived2[i].getVer() == pu.mvdAffi[eRefList][1].getVer()) ? 0 : 1;
          mvsdIdx += bin << shift;
          shift++;
          break;
        }
      }
    }
    if (pu.cu->affineType == AFFINEMODEL_6PARAM)
    {
      if (pu.mvdAffi[eRefList][2].getHor())
      {
        for (int i = 0; i < (int)cMvdDerived.size(); i++)
        {
          if (cMvdDerived[i] == pu.mvdAffi[eRefList][0] && cMvdDerived2[i] == pu.mvdAffi[eRefList][1])
          {
            bin = (cMvdDerived3[i].getHor() == pu.mvdAffi[eRefList][2].getHor()) ? 0 : 1;
            mvsdIdx += bin << shift;
            shift++;
            break;
          }
        }
      }
      if (pu.mvdAffi[eRefList][2].getVer())
      {
        for (int i = 0; i < (int)cMvdDerived.size(); i++)
        {
          if (cMvdDerived[i] == pu.mvdAffi[eRefList][0] && cMvdDerived2[i] == pu.mvdAffi[eRefList][1] && cMvdDerived3[i].getHor() == pu.mvdAffi[eRefList][2].getHor())
          {
            bin = (cMvdDerived3[i].getVer() == pu.mvdAffi[eRefList][2].getVer()) ? 0 : 1;
            mvsdIdx += bin << shift;
            shift++;
            break;
          }
        }
      }
    }
    return mvsdIdx;
  }
  
  void InterPrediction::deriveMVDFromMVSDIdxAffine(PredictionUnit& pu, RefPicList eRefList, std::vector<Mv>& cMvdDerived, std::vector<Mv>& cMvdDerived2, std::vector<Mv>& cMvdDerived3)
  {
    int mvsdIdx = pu.mvsdIdx[eRefList];
    int bin = 0;
    
    if (pu.mvdAffi[eRefList][0].getHor())
    {
      bin = mvsdIdx & 1;
      int val = bin ? -cMvdDerived[0].getHor() : cMvdDerived[0].getHor();
      pu.mvdAffi[eRefList][0].setHor(val);
      mvsdIdx >>= 1;
    }
    if (pu.mvdAffi[eRefList][0].getVer())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i].getHor() == pu.mvdAffi[eRefList][0].getHor())
        {
          bin = mvsdIdx & 1;
          int val = bin ? -cMvdDerived[i].getVer() : cMvdDerived[i].getVer();
          pu.mvdAffi[eRefList][0].setVer(val);
          mvsdIdx >>= 1;
          break;
        }
      }
    }
    if (pu.mvdAffi[eRefList][1].getHor())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i] == pu.mvdAffi[eRefList][0])
        {
          bin = mvsdIdx & 1;
          int val = bin ? -cMvdDerived2[i].getHor() : cMvdDerived2[i].getHor();
          pu.mvdAffi[eRefList][1].setHor(val);
          mvsdIdx >>= 1;
          break;
        }
      }
    }
    if (pu.mvdAffi[eRefList][1].getVer())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i] == pu.mvdAffi[eRefList][0] && cMvdDerived2[i].getHor() == pu.mvdAffi[eRefList][1].getHor())
        {
          bin = mvsdIdx & 1;
          int val = bin ? -cMvdDerived2[i].getVer() : cMvdDerived2[i].getVer();
          pu.mvdAffi[eRefList][1].setVer(val);
          mvsdIdx >>= 1;
          break;
        }
      }
    }
    if (pu.cu->affineType == AFFINEMODEL_6PARAM)
    {
      if (pu.mvdAffi[eRefList][2].getHor())
      {
        for (int i = 0; i < (int)cMvdDerived.size(); i++)
        {
          if (cMvdDerived[i] == pu.mvdAffi[eRefList][0] && cMvdDerived2[i] == pu.mvdAffi[eRefList][1])
          {
            bin = mvsdIdx & 1;
            int val = bin ? -cMvdDerived3[i].getHor() : cMvdDerived3[i].getHor();
            pu.mvdAffi[eRefList][2].setHor(val);
            mvsdIdx >>= 1;
            break;
          }
        }
      }
      if (pu.mvdAffi[eRefList][2].getVer())
      {
        for (int i = 0; i < (int)cMvdDerived.size(); i++)
        {
          if (cMvdDerived[i] == pu.mvdAffi[eRefList][0] && cMvdDerived2[i] == pu.mvdAffi[eRefList][1] && cMvdDerived3[i].getHor() == pu.mvdAffi[eRefList][2].getHor())
          {
            bin = mvsdIdx & 1;
            int val = bin ? -cMvdDerived3[i].getVer() : cMvdDerived3[i].getVer();
            pu.mvdAffi[eRefList][2].setVer(val);
            mvsdIdx >>= 1;
            break;
          }
        }
      }
    }
  }
  int InterPrediction::deriveMVSDIdxFromMVDTrans(Mv cMvd, std::vector<Mv>& cMvdDerived)
  {
    int mvsdIdx = 0;
    int shift = 0;
    int bin = 0;
    
    if (cMvd.getHor())
    {
      bin = (cMvdDerived[0].getHor() == cMvd.getHor()) ? 0 : 1;
      mvsdIdx += bin << shift;
      shift++;
    }
    if (cMvd.getVer())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i].getHor() == cMvd.getHor())
        {
          bin = (cMvdDerived[i].getVer() == cMvd.getVer()) ? 0 : 1;
          mvsdIdx += bin << shift;
          shift++;
          break;
        }
      }
    }
    return mvsdIdx;
  }
  Mv InterPrediction::deriveMVDFromMVSDIdxTrans(int mvsdIdx, std::vector<Mv>& cMvdDerived)
  {
    int bin = 0;
    Mv cMvd = Mv(0, 0);
    if (cMvdDerived[0].getHor())
    {
      bin = mvsdIdx & 1;
      int val = bin ? -cMvdDerived[0].getHor() : cMvdDerived[0].getHor();
      cMvd.setHor(val);
      mvsdIdx >>= 1;
    }
    if (cMvdDerived[0].getVer())
    {
      for (int i = 0; i < (int)cMvdDerived.size(); i++)
      {
        if (cMvdDerived[i].getHor() == cMvd.getHor())
        {
          bin = mvsdIdx & 1;
          int val = bin ? -cMvdDerived[i].getVer() : cMvdDerived[i].getVer();
          cMvd.setVer(val);
          mvsdIdx >>= 1;
          break;
        }
      }
    }
    return cMvd;
  }
#endif
