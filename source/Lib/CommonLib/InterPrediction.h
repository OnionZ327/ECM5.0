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

/** \file     InterPrediction.h
    \brief    inter prediction class (header)
*/

#ifndef __INTERPREDICTION__
#define __INTERPREDICTION__


// Include files
#include "InterpolationFilter.h"
#include "WeightPrediction.h"

#include "Buffer.h"
#include "Unit.h"
#include "Picture.h"

#include "RdCost.h"
#include "ContextModelling.h"
#if JVET_Y0065_GPM_INTRA
#include "IntraPrediction.h"
#endif
// forward declaration
class Mv;
#if INTER_LIC || (TM_AMVP || TM_MRG) || JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
class Reshape;
#endif

//! \ingroup CommonLib
//! \{


// ====================================================================================================================
// Class definition
// ====================================================================================================================

class InterPrediction : public WeightPrediction
{
#if INTER_LIC
public:
  PelUnitBuf        m_predictionBeforeLIC;
  bool              m_storeBeforeLIC;
#endif
#if JVET_Z0136_OOB
  bool isMvOOB(const Mv& rcMv, const struct Position pos, const struct Size size, const SPS* sps, const PPS* pps, bool *mcMask, bool *mcMaskChroma, bool lumaOnly = false);
  bool isMvOOBSubBlk(const Mv& rcMv, const struct Position pos, const struct Size size, const SPS* sps, const PPS* pps, bool *mcMask, int mcStride, bool *mcMaskChroma, int mcCStride, bool lumaOnly = false);
#endif
#if INTER_LIC || (TM_AMVP || TM_MRG) || JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC // note: already refactor
  Reshape*          m_pcReshape;
#endif

private:
#if INTER_LIC
  static const int  m_LICShift     = 5;
  static const int  m_LICRegShift  = 7;
  static const int  m_LICShiftDiff = 12;
  int               m_LICMultApprox[64];
#endif

#if INTER_LIC || JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  // buffer size for left/above current templates and left/above reference templates
  Pel* m_pcLICRefLeftTemplate;
  Pel* m_pcLICRefAboveTemplate;
  Pel* m_pcLICRecLeftTemplate;
  Pel* m_pcLICRecAboveTemplate;
#endif
#if TM_AMVP || TM_MRG
  // buffer size for left/above current templates
  Pel* m_pcCurTplLeft;
  Pel* m_pcCurTplAbove;
  Pel* m_pcRefTplLeft;
  Pel* m_pcRefTplAbove;
#endif
#if TM_AMVP
  AMVPInfo m_tplAmvpInfo[NUM_IMV_MODES][NUM_REF_PIC_LIST_01][MAX_NUM_REF];
#if INTER_LIC
  AMVPInfo m_tplAmvpInfoLIC[NUM_IMV_MODES][NUM_REF_PIC_LIST_01][MAX_NUM_REF];
#endif
#endif

protected:
  InterpolationFilter  m_if;

  Pel*                 m_acYuvPred            [NUM_REF_PIC_LIST_01][MAX_NUM_COMPONENT];
  Pel*                 m_filteredBlock        [LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL][LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL][MAX_NUM_COMPONENT];
  Pel*                 m_filteredBlockTmp     [LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL][MAX_NUM_COMPONENT];

#if MULTI_HYP_PRED
  PelStorage           m_additionalHypothesisStorage;
  int                  m_multiHypActive = 0;
#endif

  ChromaFormat         m_currChromaFormat;

  ComponentID          m_maxCompIDToPred;      ///< tells the predictor to only process the components up to (inklusive) this one - useful to skip chroma components during RD-search

  RdCost*              m_pcRdCost;

  int                  m_iRefListIdx;
  PelStorage           m_geoPartBuf[2];
  Mv*                  m_storedMv;

 /*buffers for bilinear Filter data for DMVR refinement*/
  Pel*                 m_cYuvPredTempDMVRL0;
  Pel*                 m_cYuvPredTempDMVRL1;
  int                  m_biLinearBufStride;
  /*buffers for padded data*/
  PelUnitBuf           m_cYuvRefBuffDMVRL0;
  PelUnitBuf           m_cYuvRefBuffDMVRL1;
  Pel*                 m_cRefSamplesDMVRL0[MAX_NUM_COMPONENT];
  Pel*                 m_cRefSamplesDMVRL1[MAX_NUM_COMPONENT];
  Mv m_pSearchOffset[25] = { Mv(-2,-2), Mv(-1,-2), Mv(0,-2), Mv(1,-2), Mv(2,-2),
                             Mv(-2,-1), Mv(-1,-1), Mv(0,-1), Mv(1,-1), Mv(2,-1),
                             Mv(-2, 0), Mv(-1, 0), Mv(0, 0), Mv(1, 0), Mv(2, 0),
                             Mv(-2, 1), Mv(-1, 1), Mv(0, 1), Mv(1, 1), Mv(2, 1),
                             Mv(-2, 2), Mv(-1, 2), Mv(0, 2), Mv(1, 2), Mv(2, 2) };
  uint64_t m_SADsArray[((2 * DMVR_NUM_ITERATION) + 1) * ((2 * DMVR_NUM_ITERATION) + 1)];
#if MULTI_PASS_DMVR
#if JVET_X0049_BDMVR_SW_OPT
  Mv                   m_searchEnlargeOffsetBilMrg[5][BDMVR_INTME_AREA];
  uint16_t             m_searchEnlargeOffsetToIdx[5][BDMVR_INTME_AREA];
  uint16_t             m_searchEnlargeOffsetNum[5];
  uint64_t             m_sadEnlargeArrayBilMrg[BDMVR_INTME_AREA];
#else
  Mv                   m_searchEnlargeOffsetBilMrg[BDMVR_INTME_AREA];
  uint64_t             m_sadEnlargeArrayBilMrg[BDMVR_INTME_AREA];
  int                  m_searchPriorityBilMrg[BDMVR_INTME_AREA];
#endif
  int                  m_costShiftBilMrg1[BDMVR_INTME_AREA];
  int                  m_costShiftBilMrg2[BDMVR_INTME_AREA];
#endif

  Pel                  m_gradBuf[2][(AFFINE_MIN_BLOCK_SIZE + 2) * (AFFINE_MIN_BLOCK_SIZE + 2)];
  int                  m_dMvBuf[2][16 * 2];
  bool                 m_skipPROF;
  bool                 m_encOnly;
  bool                 m_isBi;

  Pel*                 m_gradX0;
  Pel*                 m_gradY0;
  Pel*                 m_gradX1;
  Pel*                 m_gradY1;
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  Pel*                 m_absGx;
  Pel*                 m_absGy;
  Pel*                 m_dIx;
  Pel*                 m_dIy;
  Pel*                 m_dI;
  Pel*                 m_signGxGy;
  int*                 m_tmpx_pixel_32bit;
  int*                 m_tmpy_pixel_32bit;
  int*                 m_sumAbsGX_pixel_32bit;
  int*                 m_sumAbsGY_pixel_32bit;
  int*                 m_sumDIX_pixel_32bit;
  int*                 m_sumDIY_pixel_32bit;
  int*                 m_sumSignGY_GX_pixel_32bit;
  bool                 m_bdofMvRefined;
  Mv                   m_bdofSubPuMvOffset[BDOF_SUBPU_MAX_NUM];
#endif
  bool                 m_subPuMC;

  int                  m_IBCBufferWidth;
#if JVET_Z0153_IBC_EXT_REF
  int                  m_IBCBufferHeight;
#endif
#if JVET_Z0118_GDR
  PelStorage           m_IBCBuffer0; // for dirty
  PelStorage           m_IBCBuffer1; // for clean
#else
  PelStorage           m_IBCBuffer;
#endif
  void xIntraBlockCopy          (PredictionUnit &pu, PelUnitBuf &predBuf, const ComponentID compID);
  int             rightShiftMSB(int numer, int    denom);
#if MULTI_PASS_DMVR
  void            xPredInterUni            ( const PredictionUnit &pu, const RefPicList &eRefPicList, PelUnitBuf &pcYuvPred,
                                             const bool &bi, const bool &bioApplied, const bool luma, const bool chroma,
                                             const bool isBdofMvRefine = false);
  void            xPredInterBiSubPuBDOF    ( PredictionUnit &pu, PelUnitBuf &pcYuvPred, const bool luma, const bool chroma );
#if JVET_Z0136_OOB
  void            applyBiOptFlow           ( const bool isBdofMvRefine, const int bdofBlockOffset, const PredictionUnit &pu,
                                             const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1,
                                             PelUnitBuf &yuvDst, const BitDepths &clipBitDepths, bool *mcMask[2] = NULL, bool *mcMaskChroma[2] = NULL, bool *isOOB = NULL);
#else
  void            applyBiOptFlow           ( const bool isBdofMvRefine, const int bdofBlockOffset, const PredictionUnit &pu,
                                             const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1,
                                             PelUnitBuf &yuvDst, const BitDepths &clipBitDepths );
#endif
#else
#if JVET_Z0136_OOB
  void            applyBiOptFlow(const PredictionUnit &pu, const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1, PelUnitBuf &yuvDst, const BitDepths &clipBitDepths, bool *mcMask[2] = NULL, bool *mcMaskChroma[2] = NULL, bool *isOOB = NULL);
#else
  void            applyBiOptFlow(const PredictionUnit &pu, const CPelUnitBuf &yuvSrc0, const CPelUnitBuf &yuvSrc1, const int &refIdx0, const int &refIdx1, PelUnitBuf &yuvDst, const BitDepths &clipBitDepths);
#endif
  void            xPredInterUni ( const PredictionUnit& pu, const RefPicList& eRefPicList, PelUnitBuf& pcYuvPred, const bool& bi
                                  , const bool& bioApplied
                                  , const bool luma, const bool chroma
  );
#endif
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
#if JVET_Z0136_OOB
  void            subBlockBiOptFlow        ( Pel* dstY, const int dstStride, const Pel* src0, const int src0Stride, const Pel* src1,
                                             const int src1Stride, int bioParamOffset, const int bioParamStride, int width, int height,
                                             const ClpRng& clpRng, const int shiftNum, const int offset, const int limit, bool *mcMask[2], int mcStride, bool *isOOB = NULL);
#else
  void            subBlockBiOptFlow        ( Pel* dstY, const int dstStride, const Pel* src0, const int src0Stride, const Pel* src1,
                                             const int src1Stride, int bioParamOffset, const int bioParamStride, int width, int height,
                                             const ClpRng& clpRng, const int shiftNum, const int offset, const int limit );
#endif
#endif
#if ENABLE_OBMC
  PelStorage           m_tmpObmcBufL0;
  PelStorage           m_tmpObmcBufT0;
  PelStorage           m_tmpSubObmcBuf;
#endif
#if MULTI_PASS_DMVR
  void xPredInterBiBDMVR        ( PredictionUnit &pu, PelUnitBuf &pcYuvPred, const bool luma, const bool chroma, PelUnitBuf *yuvPredTmp = NULL );
#endif
  void xPredInterBi             ( PredictionUnit& pu, PelUnitBuf &pcYuvPred, const bool luma = true, const bool chroma = true, PelUnitBuf* yuvPredTmp = NULL );
  void xPredInterBlk            ( const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv& _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng
                                 , const bool& bioApplied
                                 , bool isIBC
                                 , const std::pair<int, int> scalingRatio = SCALE_1X
                                 , SizeType dmvrWidth = 0
                                 , SizeType dmvrHeight = 0
                                 , bool bilinearMC = false
                                 , Pel *srcPadBuf = NULL
                                 , int32_t srcPadStride = 0
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
                                 , bool isAML = false
#if INTER_LIC
                                 , bool doLic = false
                                 , Mv   mvCurr = Mv(0, 0)
#endif
#endif
#if JVET_Z0061_TM_OBMC
                                 , bool fastOBMC = false
#endif
                                 );

  void xAddBIOAvg4              (const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, Pel *dst, int dstStride, const Pel *gradX0, const Pel *gradX1, const Pel *gradY0, const Pel*gradY1, int gradStride, int width, int height, int tmpx, int tmpy, int shift, int offset, const ClpRng& clpRng);
  void xBioGradFilter           (Pel* pSrc, int srcStride, int width, int height, int gradStride, Pel* gradX, Pel* gradY, int bitDepth);
  void xCalcBIOPar              (const Pel* srcY0Temp, const Pel* srcY1Temp, const Pel* gradX0, const Pel* gradX1, const Pel* gradY0, const Pel* gradY1, int* dotProductTemp1, int* dotProductTemp2, int* dotProductTemp3, int* dotProductTemp5, int* dotProductTemp6, const int src0Stride, const int src1Stride, const int gradStride, const int widthG, const int heightG, int bitDepth);
  void xCalcBlkGradient         (int sx, int sy, int    *arraysGx2, int     *arraysGxGy, int     *arraysGxdI, int     *arraysGy2, int     *arraysGydI, int     &sGx2, int     &sGy2, int     &sGxGy, int     &sGxdI, int     &sGydI, int width, int height, int unitSize);
#if MULTI_PASS_DMVR
#if JVET_Z0136_OOB
  void xWeightedAverage         ( const bool isBdofMvRefine, const int bdofBlockOffset, const PredictionUnit& pu, const CPelUnitBuf& pcYuvSrc0, const CPelUnitBuf& pcYuvSrc1, PelUnitBuf& pcYuvDst, const BitDepths& clipBitDepths, const ClpRngs& clpRngs, const bool& bioApplied, const bool lumaOnly = false, const bool chromaOnly = false, PelUnitBuf* yuvDstTmp = NULL, bool *mcMask[2] = NULL, int mcStride = -1, bool *mcMaskChroma[2] = NULL, int mcCStride = -1, bool *isOOB = NULL );
#else
  void xWeightedAverage         ( const bool isBdofMvRefine, const int bdofBlockOffset, const PredictionUnit& pu, const CPelUnitBuf& pcYuvSrc0, const CPelUnitBuf& pcYuvSrc1, PelUnitBuf& pcYuvDst, const BitDepths& clipBitDepths, const ClpRngs& clpRngs, const bool& bioApplied, const bool lumaOnly = false, const bool chromaOnly = false, PelUnitBuf* yuvDstTmp = NULL );
#endif
#else
#if JVET_Z0136_OOB
  void xWeightedAverage         ( const PredictionUnit& pu, const CPelUnitBuf& pcYuvSrc0, const CPelUnitBuf& pcYuvSrc1, PelUnitBuf& pcYuvDst, const BitDepths& clipBitDepths, const ClpRngs& clpRngs, const bool& bioApplied, const bool lumaOnly = false, const bool chromaOnly = false, PelUnitBuf* yuvDstTmp = NULL, bool *mcMask[2] = NULL, int mcStride = -1, bool *mcMaskChroma[2] = NULL, int mcCStride = -1, bool *isOOB = NULL );
#else
  void xWeightedAverage         ( const PredictionUnit& pu, const CPelUnitBuf& pcYuvSrc0, const CPelUnitBuf& pcYuvSrc1, PelUnitBuf& pcYuvDst, const BitDepths& clipBitDepths, const ClpRngs& clpRngs, const bool& bioApplied, const bool lumaOnly = false, const bool chromaOnly = false, PelUnitBuf* yuvDstTmp = NULL );
#endif
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
#if !INTER_LIC
  template <bool trueAfalseL>
  void xGetPredBlkTpl(const CodingUnit& cu, const ComponentID compID, const CPelBuf& refBuf, const Mv& mv, const int posW, const int posH, const int tplSize, Pel* predBlkTpl
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
                      , bool AML = false
#endif
                      );
#endif
  void xWeightedAverageY(const PredictionUnit& pu, const CPelUnitBuf& pcYuvSrc0, const CPelUnitBuf& pcYuvSrc1, PelUnitBuf& pcYuvDst, const BitDepths& clipBitDepths, const ClpRngs& clpRngs);
#endif
#if JVET_W0090_ARMC_TM
  void xPredAffineTpl(const PredictionUnit &pu, const RefPicList &eRefPicList, int* numTemplate, Pel* refLeftTemplate, Pel* refAboveTemplate);
#endif
#if AFFINE_ENC_OPT
#if JVET_Z0136_OOB
  void xPredAffineBlk           ( const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv* _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng, RefPicList eRefPicList, const bool genChromaMv = false, const std::pair<int, int> scalingRatio = SCALE_1X, const bool calGradient = false);
#else
  void xPredAffineBlk           ( const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv* _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng, const bool genChromaMv = false, const std::pair<int, int> scalingRatio = SCALE_1X, const bool calGradient = false);
#endif
#else
#if JVET_Z0136_OOB
  void xPredAffineBlk           ( const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv* _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng, RefPicList eRefPicList, const bool genChromaMv = false, const std::pair<int, int> scalingRatio = SCALE_1X);
#else
  void xPredAffineBlk           ( const ComponentID& compID, const PredictionUnit& pu, const Picture* refPic, const Mv* _mv, PelUnitBuf& dstPic, const bool& bi, const ClpRng& clpRng, const bool genChromaMv = false, const std::pair<int, int> scalingRatio = SCALE_1X );
#endif
#endif

  static bool xCheckIdenticalMotion( const PredictionUnit& pu );

  void xSubPuMC(PredictionUnit& pu, PelUnitBuf& predBuf, const RefPicList &eRefPicList = REF_PIC_LIST_X, const bool luma = true, const bool chroma = true);
#if ENABLE_OBMC
  void xSubblockOBMC(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst, PelUnitBuf &pcYuvPredSrc, int iDir, bool bSubMotion = false);
  void xSubBlockMotionCompensation(PredictionUnit &pu, PelUnitBuf &pcYuvPred);
  void xSubblockOBMCBlending(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst, PelUnitBuf &pcYuvPredSrc1, PelUnitBuf &pcYuvPredSrc2, PelUnitBuf &pcYuvPredSrc3, PelUnitBuf &pcYuvPredSrc4, bool isAboveAvail = false, bool isLeftAvail = false, bool isBelowAvail = false, bool isRightAvail = false, bool bSubMotion = false);
#endif
#if !BDOF_RM_CONSTRAINTS
  void xSubPuBio(PredictionUnit& pu, PelUnitBuf& predBuf, const RefPicList &eRefPicList = REF_PIC_LIST_X, PelUnitBuf* yuvDstTmp = NULL);
#endif
  void destroy();


  MotionInfo      m_SubPuMiBuf[(MAX_CU_SIZE * MAX_CU_SIZE) >> (MIN_CU_LOG2 << 1)];
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
  Pel*   m_acYuvCurAMLTemplate[2][MAX_NUM_COMPONENT];   //0: top, 1: left
  bool   m_bAMLTemplateAvailabe[2];
  Pel*   m_acYuvRefAboveTemplate[2][MAX_NUM_COMPONENT];   //0: list0, 1: list1
  Pel*   m_acYuvRefLeftTemplate[2][MAX_NUM_COMPONENT];   //0: list0, 1: list1
  Pel*   m_acYuvRefAMLTemplate[2][MAX_NUM_COMPONENT];   //0: top, 1: left
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  Pel*   m_acYuvRefAMLTemplatePart0[2 + 2];   //0: top, 1: left, 2: top, 3: left
  Pel*   m_acYuvRefAMLTemplatePart1[2 + 2];   //0: top, 1: left, 2: top, 3: left
#endif
#endif
#if JVET_Z0061_TM_OBMC
  Pel *m_acYuvRefAboveTemplateOBMC[2][MAX_NUM_COMPONENT];   // 0: current motion, 1: neighbour motion
  Pel *m_acYuvRefLeftTemplateOBMC[2][MAX_NUM_COMPONENT];    // 0: current motion, 1: neighbour motion
  Pel *m_acYuvBlendTemplateOBMC[2][MAX_NUM_COMPONENT];      // 0: top, 1: left
#endif
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  bool   m_tplWeightTblInitialized;
  Pel*   m_tplWeightTblDict[GEO_NUM_CU_SIZE][GEO_NUM_CU_SIZE][GEO_NUM_PARTITION_MODE];
  Pel    m_tplColWeightTblDict[GEO_NUM_CU_SIZE][GEO_NUM_CU_SIZE][GEO_NUM_PARTITION_MODE][GEO_MAX_CU_SIZE * GEO_MODE_SEL_TM_SIZE];
  Pel**  m_tplWeightTbl;
  Pel    (*m_tplColWeightTbl)[GEO_MAX_CU_SIZE * GEO_MODE_SEL_TM_SIZE];
#endif

#if JVET_J0090_MEMORY_BANDWITH_MEASURE
  CacheModel      *m_cacheModel;
#endif
  PelStorage       m_colorTransResiBuf[3];  // 0-org; 1-act; 2-tmp
#if MULTI_HYP_PRED
  void xAddHypMC(PredictionUnit& pu, PelUnitBuf& predBuf, PelUnitBuf* predBufWOBIO, const bool lumaOnly = false);
#endif

public:
  InterPrediction();
  virtual ~InterPrediction();

#if INTER_LIC || (TM_AMVP || TM_MRG) || JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
#if JVET_Z0153_IBC_EXT_REF
  void    init                (RdCost* pcRdCost, ChromaFormat chromaFormatIDC, const int ctuSize, Reshape* reshape, const int picWidth);
#else
  void    init                (RdCost* pcRdCost, ChromaFormat chromaFormatIDC, const int ctuSize, Reshape* reshape);
#endif
#else
  void    init                (RdCost* pcRdCost, ChromaFormat chromaFormatIDC, const int ctuSize);
#endif

#if JVET_Z0054_BLK_REF_PIC_REORDER
  void     setUniRefIdxLC(PredictionUnit &pu);
  void     setUniRefListAndIdx(PredictionUnit &pu);
  void     reorderRefCombList(PredictionUnit &pu, std::vector<RefListAndRefIdx> &refListComb
    , RefPicList currRefList
    , std::vector<MotionInfoPred> &miPredList
  );
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  void      deriveMVDCandVecFromMotionInforPred(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> &cMvdDerivedVec);
  void      deriveAffineMVDCandVecFromMotionInforPred(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> cMvdDerivedVec[3]);
  void      deriveMVDCandVecFromMotionInforPredGeneral(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> &cMvdDerivedVec);
  void      deriveAffineMVDCandVecFromMotionInforPredGeneral(const PredictionUnit &pu, std::vector<MotionInfoPred> &miPredList, RefPicList eRefPicList, std::vector<Mv> cMvdDerivedVec[3]);
#endif
  void     setBiRefPairIdx(PredictionUnit &pu);
  void     setBiRefIdx(PredictionUnit &pu);
  void     reorderRefPairList(PredictionUnit &pu, std::vector<RefPicPair> &refPairList
    , std::vector<MotionInfoPred> &miPredList
  );
#endif

  // inter
  void    motionCompensation  (PredictionUnit &pu, PelUnitBuf& predBuf, const RefPicList &eRefPicList = REF_PIC_LIST_X
    , const bool luma = true, const bool chroma = true
    , PelUnitBuf* predBufWOBIO = NULL
  );
  void    motionCompensation  (PredictionUnit &pu, const RefPicList &eRefPicList = REF_PIC_LIST_X
    , const bool luma = true, const bool chroma = true
  );
  void    motionCompensation  (CodingUnit &cu,     const RefPicList &eRefPicList = REF_PIC_LIST_X
    , const bool luma = true, const bool chroma = true
  );
#if ENABLE_OBMC
  void    subBlockOBMC(PredictionUnit  &pu, PelUnitBuf *pDst = nullptr);
#endif
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  void    initTplWeightTable();

  template <bool regular, uint8_t top0Left1TrLeft2 = 0>
  Pel*    getTplWeightTableCU(int splitDir)
  {
    if (regular)
    {
      return m_tplWeightTbl[splitDir];
    }
    else if (top0Left1TrLeft2 == 2)
    {
      return m_tplColWeightTbl[splitDir];
    }
    else
    {
      int16_t angle = g_GeoParams[splitDir][0];
      if (g_angle2mirror[angle] == 2)
      {
        return m_tplWeightTbl[splitDir] + (top0Left1TrLeft2 == 0 ? GEO_WEIGHT_MASK_SIZE_EXT * GEO_MODE_SEL_TM_SIZE : -GEO_MODE_SEL_TM_SIZE); // Shift to template pos
      }
      else if (g_angle2mirror[angle] == 1)
      {
        return m_tplWeightTbl[splitDir] - (top0Left1TrLeft2 == 0 ? GEO_WEIGHT_MASK_SIZE_EXT * GEO_MODE_SEL_TM_SIZE : -GEO_MODE_SEL_TM_SIZE); // Shift to template pos;
      }
      else
      {
        return m_tplWeightTbl[splitDir] - (top0Left1TrLeft2 == 0 ? GEO_WEIGHT_MASK_SIZE_EXT * GEO_MODE_SEL_TM_SIZE : GEO_MODE_SEL_TM_SIZE); // Shift to template pos;
      }
    }
  }
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  void    deriveGpmSplitMode(PredictionUnit& pu, MergeCtx &geoMrgCtx
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
                           , MergeCtx(&geoTmMrgCtx)[GEO_NUM_TM_MV_CAND]
#endif
#if JVET_Y0065_GPM_INTRA
                           , IntraPrediction* pcIntraPred
#endif
  );
#endif
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  void    motionCompensationGeo(CodingUnit &cu, MergeCtx &geoMrgCtx
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
                              , MergeCtx (&geoTmMrgCtx)[GEO_NUM_TM_MV_CAND]
#endif
#if JVET_Y0065_GPM_INTRA
                              , IntraPrediction* pcIntraPred, std::vector<Pel>* reshapeLUT
#endif
  );
#else
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
#if JVET_Y0065_GPM_INTRA
  void    motionCompensationGeo(CodingUnit &cu, MergeCtx &geoMrgCtx, MergeCtx &geoTmMrgCtx0, MergeCtx &geoTmMrgCtx1, IntraPrediction* pcIntraPred, std::vector<Pel>* reshapeLUT);
#else
  void    motionCompensationGeo(CodingUnit &cu, MergeCtx &geoMrgCtx, MergeCtx &geoTmMrgCtx0, MergeCtx &geoTmMrgCtx1);
#endif
#else
#if JVET_Y0065_GPM_INTRA
  void    motionCompensationGeo( CodingUnit &cu, MergeCtx &geoMrgCtx, IntraPrediction* pcIntraPred, std::vector<Pel>* reshapeLUT );
#else
  void    motionCompensationGeo(CodingUnit &cu, MergeCtx &GeoMrgCtx);
#endif
#endif
#endif
  void    weightedGeoBlk(PredictionUnit &pu, const uint8_t splitDir, int32_t channel, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1);
#if JVET_Y0065_GPM_INTRA
  void    weightedGeoBlkRounded( PredictionUnit &pu, const uint8_t splitDir, int32_t channel, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1);
#endif
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
  void    getBestGeoTMModeList(PredictionUnit &pu, uint8_t& numValidInList, uint8_t(&modeList)[GEO_NUM_SIG_PARTMODE], Pel* (&pRefTopPart0)[GEO_NUM_TM_MV_CAND], Pel* (&pRefLeftPart0)[GEO_NUM_TM_MV_CAND], Pel* (&pRefTopPart1)[GEO_NUM_TM_MV_CAND], Pel* (&pRefLeftPart1)[GEO_NUM_TM_MV_CAND]);
#endif
  void    getBestGeoModeList  (PredictionUnit &pu, uint8_t& numValidInList, uint8_t(&modeList)[GEO_NUM_SIG_PARTMODE], Pel* pRefTopPart0, Pel* pRefLeftPart0, Pel* pRefTopPart1, Pel* pRefLeftPart1
#if JVET_Y0065_GPM_INTRA
                             , Pel** pIntraRefTopPart0 = nullptr, Pel** pIntraRefLeftPart0 = nullptr, Pel** pIntraRefTopPart1 = nullptr, Pel** pIntraRefLeftPart1 = nullptr
#endif
  );
  template <bool trueTFalseL>
  void    weightedGeoTpl(PredictionUnit &pu, const uint8_t splitDir, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1);
#endif
  void xPrefetch(PredictionUnit& pu, PelUnitBuf &pcPad, RefPicList refId, bool forLuma);
  void xPad(PredictionUnit& pu, PelUnitBuf &pcPad, RefPicList refId);
  void xFinalPaddedMCForDMVR(PredictionUnit& pu, PelUnitBuf &pcYuvSrc0, PelUnitBuf &pcYuvSrc1, PelUnitBuf &pcPad0, PelUnitBuf &pcPad1, const bool bioApplied
    , const Mv startMV[NUM_REF_PIC_LIST_01]
    , bool blockMoved
  );
  void xBIPMVRefine(int bd, Pel *pRefL0, Pel *pRefL1, uint64_t& minCost, int16_t *deltaMV, uint64_t *pSADsArray, int width, int height);
  uint64_t xDMVRCost(int bitDepth, Pel* pRef, uint32_t refStride, const Pel* pOrg, uint32_t orgStride, int width, int height);
  void xinitMC(PredictionUnit& pu, const ClpRngs &clpRngs);
  void xProcessDMVR(PredictionUnit& pu, PelUnitBuf &pcYuvDst, const ClpRngs &clpRngs, const bool bioApplied );
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if JVET_Z0054_BLK_REF_PIC_REORDER
  void deriveMVDcand(const PredictionUnit& pu, RefPicList eRefPicList, std::vector<Mv>& cMvdCandList);
  void deriveMVDcandAffine(const PredictionUnit& pu, RefPicList eRefPicList, std::vector<Mv> cMvdCandList[3]);
#endif
  void deriveMvdSign(const Mv& cMvPred, const Mv& cMvdKnownAtDecoder, PredictionUnit& pu, RefPicList eRefPicList, int iRefIdx, std::vector<Mv>& cMvdDerived);
  int deriveMVSDIdxFromMVDTrans(Mv cMvd, std::vector<Mv>& cMvdDerived);
  Mv deriveMVDFromMVSDIdxTrans(int mvsdIdx, std::vector<Mv>& cMvdDerived);
  void deriveMvdSignSMVD(const Mv& cMvPred, const Mv& cMvPred2, const Mv& cMvdKnownAtDecoder, PredictionUnit& pu, std::vector<Mv>& cMvdDerived);
  void deriveMvdSignAffine(const Mv& cMvPred, const Mv& cMvPred2, const Mv& cMvPred3,
#if JVET_Z0054_BLK_REF_PIC_REORDER
                           const Mv cMvdKnownAtDecoder[3],
#else
                           const Mv& cMvdKnownAtDecoder, const Mv& cMvdKnownAtDecoder2, const Mv& cMvdKnownAtDecoder3,
#endif
                           PredictionUnit& pu, RefPicList eRefList, int refIdx, std::vector<Mv>& cMvdDerived, std::vector<Mv>& cMvdDerived2, std::vector<Mv>& cMvdDerived3);
  Distortion xGetSublkTemplateCost(const CodingUnit& cu, const ComponentID compID, const Picture& refPic, const Mv& mv, const int sublkWidth, const int sublkHeight,
    const int posW, const int posH, int* numTemplate, Pel* refLeftTemplate, Pel* refAboveTemplate, Pel* recLeftTemplate, Pel* recAboveTemplate);
  int deriveMVSDIdxFromMVDAffine(PredictionUnit& pu, RefPicList eRefList, std::vector<Mv>& cMvdDerived, std::vector<Mv>& cMvdDerived2, std::vector<Mv>& cMvdDerived3);
  void deriveMVDFromMVSDIdxAffine(PredictionUnit& pu, RefPicList eRefList, std::vector<Mv>& cMvdDerived, std::vector<Mv>& cMvdDerived2, std::vector<Mv>& cMvdDerived3);
#endif
#if JVET_J0090_MEMORY_BANDWITH_MEASURE
  void    cacheAssign( CacheModel *cache );
#endif
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  static bool isSubblockVectorSpreadOverLimit( int a, int b, int c, int d, int predType );
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  void    sortInterMergeMMVDCandidates(PredictionUnit &pu, MergeCtx& mrgCtx, uint32_t * mmvdLUT, uint32_t MMVDIdx = -1);
#endif
    
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  void    sortAffineMergeCandidates(PredictionUnit pu, AffineMergeCtx& affMrgCtx, uint32_t * affMmvdLUT, uint32_t afMMVDIdx = -1);
#endif
    
#if JVET_W0090_ARMC_TM
  void    adjustInterMergeCandidates(PredictionUnit &pu, MergeCtx& mrgCtx, int mrgCandIdx = -1);
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING || JVET_Z0061_TM_OBMC
  bool    xAMLGetCurBlkTemplate(PredictionUnit& pu, int nCurBlkWidth, int nCurBlkHeight);
  bool    xAMLIsTopTempAvailable(PredictionUnit& pu);
  bool    xAMLIsLeftTempAvailable(PredictionUnit& pu);
#endif
#if JVET_Z0061_TM_OBMC
  void xOBMCWeightedAverageY(const PredictionUnit &pu, const CPelUnitBuf &pcYuvSrc0, const CPelUnitBuf &pcYuvSrc1,
                             PelUnitBuf &pcYuvDst, const BitDepths &clipBitDepths, const ClpRngs &clpRngs,
                             MotionInfo currMi);
  int  selectOBMCmode(PredictionUnit &curPu, PredictionUnit &neigPu, bool isAbove, int nLength, uint32_t minCUW,
                      Position curOffset);
  void getBlkOBMCRefTemplate(PredictionUnit &subblockPu, PelUnitBuf &pcBufPredRef, bool isAbove, MotionInfo currMi);
  bool xCheckIdenticalMotionOBMC(PredictionUnit &pu, MotionInfo tryMi);
  void xSubblockOBMCCopy(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst,
                         PelUnitBuf &pcYuvPredSrc, int iDir);
  void xSubblockTMOBMC(const ComponentID eComp, PredictionUnit &pu, PelUnitBuf &pcYuvPredDst, PelUnitBuf &pcYuvPredSrc,
                       int iDir, int iOBMCmode = 0);
#endif
#if JVET_W0090_ARMC_TM
  void    updateCandList(uint32_t uiCand, Distortion uiCost, uint32_t uiMrgCandNum, uint32_t* RdCandList, Distortion* CandCostList);
  void    updateCandInfo(MergeCtx& mrgCtx, uint32_t(*RdCandList)[MRG_MAX_NUM_CANDS], int mrgCandIdx = -1);
#endif
#if JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  void    getBlkAMLRefTemplate(PredictionUnit &pu, PelUnitBuf &pcBufPredRefTop, PelUnitBuf &pcBufPredRefLeft);
#endif
#if JVET_W0090_ARMC_TM
  void    adjustAffineMergeCandidates(PredictionUnit &pu, AffineMergeCtx& affMrgCtx, int mrgCandIdx = -1
#if JVET_Z0139_NA_AFF
    , int sortedCandNum = -1
#endif
  );
  void    updateAffineCandInfo(PredictionUnit &pu, AffineMergeCtx& affMrgCtx, uint32_t(*RdCandList)[AFFINE_MRG_MAX_NUM_CANDS], int mrgCandIdx = -1);
  void    xGetSublkAMLTemplate(const CodingUnit& cu, const ComponentID compID, const Picture& refPic, const Mv& mv, const int sublkWidth, const int sublkHeight, const int posW, const int posH, int* numTemplate, Pel* refLeftTemplate, Pel* refAboveTemplate
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
     , bool afMMVD = false
#endif
                               );
  void    getAffAMLRefTemplate(PredictionUnit &pu, PelUnitBuf &pcBufPredRefTop, PelUnitBuf &pcBufPredRefLeft);
#if JVET_Z0102_NO_ARMC_FOR_ZERO_CAND
  void adjustMergeCandidates(PredictionUnit& pu, MergeCtx& smvpMergeCandCtx, int numRetrievedMergeCand);
#endif
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
  void    adjustMergeCandidatesInOneCandidateGroup(PredictionUnit &pu, MergeCtx& smvpMergeCandCtx, int numRetrievedMergeCand, int mrgCandIdx = -1);
  void    updateCandInOneCandidateGroup(MergeCtx& mrgCtx, uint32_t* RdCandList, int numCandInCategory = -1);
#endif
#if JVET_Y0058_IBC_LIST_MODIFY
  void    adjustIBCMergeCandidates(PredictionUnit &pu, MergeCtx& mrgCtx, int mrgCandIdx = -1);
  void    updateIBCCandInfo(PredictionUnit &pu, MergeCtx& mrgCtx, uint32_t(*RdCandList)[IBC_MRG_MAX_NUM_CANDS], int mrgCandIdx = -1);
  bool    xAMLIBCGetCurBlkTemplate(PredictionUnit& pu, int nCurBlkWidth, int nCurBlkHeight);
  void    getIBCAMLRefTemplate(PredictionUnit &pu, int nCurBlkWidth, int nCurBlkHeight);
#endif
#if JVET_Z0075_IBC_HMVP_ENLARGE
  void    adjustIBCMergeCandidates(PredictionUnit &pu, MergeCtx& mrgCtx, uint32_t startPos,uint32_t endPos);
  void    updateIBCCandInfo(PredictionUnit &pu, MergeCtx& mrgCtx, uint32_t* RdCandList, uint32_t startPos,uint32_t endPos);
#endif
#endif
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  template <uint8_t partIdx, bool useDefaultPelBuffer = true>
  void    fillPartGPMRefTemplate(PredictionUnit &pu, Pel* bufTop = nullptr, Pel* bufLeft = nullptr)
  {
    if (useDefaultPelBuffer)
    {
      bufTop  = partIdx == 0 ? m_acYuvRefAMLTemplatePart0[0] : m_acYuvRefAMLTemplatePart1[0];
      bufLeft = partIdx == 0 ? m_acYuvRefAMLTemplatePart0[1] : m_acYuvRefAMLTemplatePart1[1];
    }
    PelUnitBuf pcBufPredRefTop  = (PelUnitBuf(pu.chromaFormat, PelBuf(bufTop,  pu.lwidth(),          GEO_MODE_SEL_TM_SIZE)));
    PelUnitBuf pcBufPredRefLeft = (PelUnitBuf(pu.chromaFormat, PelBuf(bufLeft, GEO_MODE_SEL_TM_SIZE, pu.lheight()        )));
    getBlkAMLRefTemplate(pu, pcBufPredRefTop, pcBufPredRefLeft);
  }

  template <uint8_t partIdx, bool useDefaultPelBuffer = true>
  void    fillPartGPMRefTemplate(PredictionUnit &pu, MergeCtx& geoMrgCtx, int candIdx, int geoMmvdIdx = -1, Pel* bufTop = nullptr, Pel* bufLeft = nullptr)
  {
#if JVET_Y0065_GPM_INTRA
    if (candIdx >= GEO_MAX_NUM_UNI_CANDS)
    {
      return;
    }
#endif
#if JVET_W0097_GPM_MMVD_TM
    if (geoMmvdIdx >= 0)
    {
      geoMrgCtx.setGeoMmvdMergeInfo(pu, candIdx, geoMmvdIdx);
    }
    else
#endif
    {
      geoMrgCtx.setMergeInfo(pu, candIdx);
    }
    
    fillPartGPMRefTemplate<partIdx, useDefaultPelBuffer>(pu, bufTop, bufLeft);
  }
#endif
#if INTER_LIC
  void xGetLICParamGeneral (const CodingUnit& cu, const ComponentID compID, int* numTemplate, Pel* refLeftTemplate, Pel* refAboveTemplate, Pel* recLeftTemplate, Pel* recAboveTemplate, int& shift, int& scale, int& offset);
  void xGetSublkTemplate   (const CodingUnit& cu, const ComponentID compID, const Picture& refPic, const Mv& mv, const int sublkWidth, const int sublkHeight, const int posW, const int posH, int* numTemplate, Pel* refLeftTemplate, Pel* refAboveTemplate, Pel* recLeftTemplate, Pel* recAboveTemplate);
  void xLocalIlluComp      (const PredictionUnit& pu, const ComponentID compID, const Picture& refPic, const Mv& mv, const bool biPred, PelBuf& dstBuf);

  template <bool trueAfalseL>
  void xGetPredBlkTpl(const CodingUnit& cu, const ComponentID compID, const CPelBuf& refBuf, const Mv& mv, const int posW, const int posH, const int tplSize, Pel* predBlkTpl
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
                      , bool AML = false
#endif
                      );
#endif

#if TM_AMVP || TM_MRG
#if TM_MRG
  void       deriveTMMv         (PredictionUnit& pu);
#endif
  Distortion deriveTMMv         (const PredictionUnit& pu, bool fillCurTpl, Distortion curBestCost, RefPicList eRefList, int refIdx, int maxSearchRounds, Mv& mv, const MvField* otherMvf = nullptr);
#if TM_AMVP
  void       clearTplAmvpBuffer ();
  void       writeTplAmvpBuffer (const AMVPInfo& src, const CodingUnit& cu, RefPicList eRefList, int refIdx);
  bool       readTplAmvpBuffer  (      AMVPInfo& dst, const CodingUnit& cu, RefPicList eRefList, int refIdx);
#endif
#endif // TM_AMVP || TM_MRG
#if TM_AMVP || TM_MRG || MULTI_PASS_DMVR
  static Distortion getDecoderSideDerivedMvCost (const Mv& mvStart, const Mv& mvCur, int searchRangeInFullPel, int weight);
  void       xBDMVRUpdateSquareSearchCostLog (Distortion* costLog, int bestDirect);
#endif
#if MULTI_PASS_DMVR
private:
  void       xBDMVRFillBlkPredPelBuffer(const PredictionUnit& pu, const Picture& refPic, const Mv &_mv, PelUnitBuf &dstBuf, const ClpRng& clpRng);
#if JVET_X0049_ADAPT_DMVR
  template <uint8_t dir>
#endif
  void       xBDMVRPreInterpolation    (const PredictionUnit& pu, const Mv (&mvCenter)[2], bool doPreInterpolationFP, bool doPreInterpolationHP);

#if JVET_X0049_BDMVR_SW_OPT
  Distortion xBDMVRGetMatchingError    (const PredictionUnit& pu, const Mv (&mv)[2], bool useMR, bool useHadmard = false );
#if JVET_X0049_ADAPT_DMVR
  template <uint8_t dir>
#endif
#else
  Distortion xBDMVRGetMatchingError    (const PredictionUnit& pu, const Mv (&mv)[2], bool useMR );
#endif
  Distortion xBDMVRGetMatchingError    (const PredictionUnit& pu, const Mv (&mv)[2], const int subPuOffset, bool useHadmard, bool useMR
                                      , bool& doPreInterpolation, int32_t searchStepShift, const Mv (&mvCenter)[2]
                                      , const Mv(&mvInitial)[2]  // only used for full-pel MVD
                                      , int nDirect              // only used for half-pel MVD
                                        );
#if JVET_X0049_BDMVR_SW_OPT
  template <bool adaptRange, bool useHadamard>
  Distortion xBDMVRMvIntPelFullSearch  (Mv&mvOffset, Distortion curBestCost, 
    const Mv(&initialMv)[2], 
    const int32_t maxSearchRounds, 
    const int maxHorOffset, const int maxVerOffset, 
    const bool earlySkip,
    const Distortion earlyTerminateTh, DistParam &cDistParam, Pel* pelBuffer[2], const int stride);
  template<bool hPel>
  Distortion xBDMVRMvSquareSearch(Mv(&curBestMv)[2], Distortion curBestCost, PredictionUnit& pu, const Mv(&initialMv)[2], int32_t maxSearchRounds, int32_t searchStepShift, bool useMR, bool useHadmard);
#if JVET_X0049_ADAPT_DMVR
  template <uint8_t dir>
  Distortion xBDMVRMvOneTemplateHPelSquareSearch(Mv(&curBestMv)[2], Distortion curBestCost, PredictionUnit& pu,
    const Mv(&initialMv)[2], int32_t maxSearchRounds, int32_t searchStepShift,
    bool useMR, bool useHadmard);
#endif
#else
  Distortion xBDMVRMvIntPelFullSearch  (Mv (&curBestMv)[2], Distortion curBestCost, PredictionUnit& pu, const Mv (&initialMv)[2], int32_t maxSearchRounds, int32_t searchStepShift, bool useMR, const int subPuBufOffset );
  Distortion xBDMVRMvSquareSearch      (Mv(&curBestMv)[2], Distortion curBestCost, PredictionUnit& pu, const Mv(&initialMv)[2], int32_t maxSearchRounds, int32_t searchStepShift, bool useMR, bool useHadmard);
#endif

  Mv*       m_bdmvrSubPuMvBuf[2];
#if JVET_X0083_BM_AMVP_MERGE_MODE
public:
  void      getAmvpMergeModeMergeList(PredictionUnit& pu, MvField* mvFieldAmListCommon, const int decAmvpRefIdx = -1);
  void      amvpMergeModeMvRefinement(PredictionUnit& pu, MvField* mvFieldAmListCommon, const int mvFieldMergeIdx, const int mvFieldAmvpIdx);
#endif

public:
  Mv*       getBdofSubPuMvOffset() {return m_bdofSubPuMvOffset;}
  void      setBdmvrSubPuMvBuf(Mv* mvBuf0, Mv* mvBuf1) { m_bdmvrSubPuMvBuf[0] = mvBuf0; m_bdmvrSubPuMvBuf[1] = mvBuf1; }
  bool      processBDMVR              (PredictionUnit& pu);
#if JVET_X0049_ADAPT_DMVR
  bool      processBDMVRPU2Dir        (PredictionUnit& pu, bool subPURefine[2], Mv(&finalMvDir)[2]);
  void      processBDMVRSubPU         (PredictionUnit& pu, bool subPURefine);
#endif
#endif
  void xFillIBCBuffer(CodingUnit &cu);
#if JVET_Z0118_GDR  
  void resetCurIBCBuffer(const ChromaFormat chromaFormatIDC, const Area ctuArea, const int ctuSize, const Pel dirtyPel);  
#endif
  void resetIBCBuffer(const ChromaFormat chromaFormatIDC, const int ctuSize);
  void resetVPDUforIBC(const ChromaFormat chromaFormatIDC, const int ctuSize, const int vSize, const int xPos, const int yPos);

  bool isLumaBvValid(const int ctuSize, const int xCb, const int yCb, const int width, const int height, const int xBv, const int yBv);

  bool xPredInterBlkRPR( const std::pair<int, int>& scalingRatio, const PPS& pps, const CompArea &blk, const Picture* refPic, const Mv& mv, Pel* dst, const int dstStride, const bool bi, const bool wrapRef, const ClpRng& clpRng, const int filterIndex, const bool useAltHpelIf = false);
#if JVET_Z0118_GDR
  void xPadIBCBuffer(const CodingStructure &cs, const UnitArea& ctuArea);
#endif

#if JVET_Z0054_BLK_REF_PIC_REORDER
private:
  bool m_fillCurTplLeftARMC;
  bool m_fillCurTplAboveARMC;
public:
  void setFillCurTplAboveARMC(bool b) { m_fillCurTplAboveARMC = b; }
  void setFillCurTplLeftARMC(bool b) { m_fillCurTplLeftARMC = b; }
#endif
};

#if TM_AMVP || TM_MRG
struct InterPredResources // Bridge required resource from InterPrediction
{
  Reshape*              m_pcReshape;
  RdCost*               m_pcRdCost;
  InterpolationFilter&  m_if;
  Pel*                  m_ifBuf;       // m_filteredBlockTmp[0][compID]: temp interpolation buffer used to buffer horizontal interpolation output before vertical one performs
  Pel*                  m_preFillBufA; // m_filteredBlock[3][1][0]: Pre-interpolation buffer used to store search area samples of above template
  Pel*                  m_preFillBufL; // m_filteredBlock[3][0][0]: Pre-interpolation buffer used to store search area samples of left  template

  InterPredResources( Reshape* pcReshape, RdCost* pcRdCost, InterpolationFilter& ifObj, Pel* ifBuf
                    , Pel* preFillBufA, Pel* preFillBufL
  )
  : m_pcReshape   (pcReshape)
  , m_pcRdCost    (pcRdCost)
  , m_if          (ifObj)
  , m_ifBuf       (ifBuf)
  , m_preFillBufA (preFillBufA)
  , m_preFillBufL (preFillBufL)
  {};
};

class TplMatchingCtrl
{
  enum TMSearchMethod
  {
    TMSEARCH_DIAMOND,
    TMSEARCH_CROSS,
    TMSEARCH_NUMBER_OF_METHODS
  };

  const CodingUnit&         m_cu;
  const PredictionUnit&     m_pu;
        InterPredResources& m_interRes;

  const Picture&    m_refPic;
  const Mv          m_mvStart;
        Mv          m_mvFinal;
  const Mv*         m_otherRefListMv;
        Distortion  m_minCost;
        bool        m_useWeight;
        int         m_maxSearchRounds;
        ComponentID m_compID;

  PelBuf m_curTplAbove;
  PelBuf m_curTplLeft;
  PelBuf m_refTplAbove;
  PelBuf m_refTplLeft;
  PelBuf m_refSrAbove; // pre-filled samples on search area
  PelBuf m_refSrLeft;  // pre-filled samples on search area

#if JVET_X0056_DMVD_EARLY_TERMINATION
  Distortion m_earlyTerminateTh;
#endif
#if MULTI_PASS_DMVR
  Distortion m_tmCostArrayDiamond[9];
  Distortion m_tmCostArrayCross[5];
#endif

public:
  TplMatchingCtrl(const PredictionUnit&     pu,
                        InterPredResources& interRes, // Bridge required resource from InterPrediction
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
  );

  bool       getTemplatePresentFlag() { return m_curTplAbove.buf != nullptr || m_curTplLeft.buf != nullptr; }
  Distortion getMinCost            () { return m_minCost; }
  Mv         getFinalMv            () { return m_mvFinal; }
  static int getDeltaMean          (const PelBuf& bufCur, const PelBuf& bufRef, const int rowSubShift, const int bd);

  template <int tplSize> void deriveMvUni    ();
  template <int tplSize> void removeHighFreq (const Picture& otherRefPic, const Mv& otherRefMv, const uint8_t curRefBcwWeight);

private:
  template <int tplSize, bool trueAfalseL>         bool       xFillCurTemplate   (Pel* tpl);
  template <int tplSize, bool trueAfalseL, int sr> PelBuf     xGetRefTemplate    (const PredictionUnit& curPu, const Picture& refPic, const Mv& _mv, PelBuf& dstBuf);
  template <int tplSize, bool trueAfalseL>         void       xRemoveHighFreq    (const Picture& otherRefPic, const Mv& otherRefMv, const uint8_t curRefBcwWeight);
  template <int tplSize, int searchPattern>         void       xRefineMvSearch    (int maxSearchRounds, int searchStepShift);
#if MULTI_PASS_DMVR
  template <int searchPattern>                      void       xNextTmCostAarray  (int bestDirect);
  template <int searchPattern>                      void       xDeriveCostBasedMv ();
  template <bool TrueX_FalseY>                      void       xDeriveCostBasedOffset (Distortion costLorA, Distortion costCenter, Distortion costRorB, int log2StepSize);
                                                    int        xBinaryDivision    (int64_t numerator, int64_t denominator, int fracBits);
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
public:
#endif
  template <int tplSize>                            Distortion xGetTempMatchError (const Mv& mv);
  template <int tplSize, bool trueAfalseL>         Distortion xGetTempMatchError (const Mv& mv);
};
#endif // TM_AMVP || TM_MRG

//! \}

#endif // __INTERPREDICTION__
