﻿/* The copyright in this software is being made available under the BSD
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

/** \file     EncSearch.cpp
 *  \brief    encoder inter search class
 */

#include "InterSearch.h"


#include "CommonLib/CommonDef.h"
#include "CommonLib/Rom.h"
#include "CommonLib/MotionInfo.h"
#include "CommonLib/Picture.h"
#include "CommonLib/UnitTools.h"
#include "CommonLib/dtrace_next.h"
#include "CommonLib/dtrace_buffer.h"
#if JVET_V0094_BILATERAL_FILTER || JVET_X0071_CHROMA_BILATERAL_FILTER
#include "CommonLib/BilateralFilter.h"
#endif
#include "CommonLib/MCTS.h"

#include "EncModeCtrl.h"
#include "EncLib.h"

#include <math.h>
#include <limits>


 //! \ingroup EncoderLib
 //! \{

static const Mv s_acMvRefineH[9] =
{
  Mv(  0,  0 ), // 0
  Mv(  0, -1 ), // 1
  Mv(  0,  1 ), // 2
  Mv( -1,  0 ), // 3
  Mv(  1,  0 ), // 4
  Mv( -1, -1 ), // 5
  Mv(  1, -1 ), // 6
  Mv( -1,  1 ), // 7
  Mv(  1,  1 )  // 8
};

static const Mv s_acMvRefineQ[9] =
{
  Mv(  0,  0 ), // 0
  Mv(  0, -1 ), // 1
  Mv(  0,  1 ), // 2
  Mv( -1, -1 ), // 5
  Mv(  1, -1 ), // 6
  Mv( -1,  0 ), // 3
  Mv(  1,  0 ), // 4
  Mv( -1,  1 ), // 7
  Mv(  1,  1 )  // 8
};

#if JVET_Z0131_IBC_BVD_BINARIZATION
void InterSearch::xEstBvdBitCosts(EstBvdBitsStruct *p)
{
  const FracBitsAccess& fracBits = m_CABACEstimator->getCtx().getFracBitsAcess();

  p->bitsGt0FlagH[0] = fracBits.getFracBitsArray(Ctx::Bvd(HOR_BVD_CTX_OFFSET)).intBits[0];
  p->bitsGt0FlagH[1] = fracBits.getFracBitsArray(Ctx::Bvd(HOR_BVD_CTX_OFFSET)).intBits[1];;

  p->bitsGt0FlagV[0] = fracBits.getFracBitsArray(Ctx::Bvd(VER_BVD_CTX_OFFSET)).intBits[0];
  p->bitsGt0FlagV[1] = fracBits.getFracBitsArray(Ctx::Bvd(VER_BVD_CTX_OFFSET)).intBits[1];

  const int epBitCost = 1 << SCALE_BITS;
  const int horCtxThre = NUM_HOR_BVD_CTX;
  const int verCtxThre = NUM_VER_BVD_CTX;

  const int horCtxOs = HOR_BVD_CTX_OFFSET;
  const int verCtxOs = VER_BVD_CTX_OFFSET;

  uint32_t singleBitH[2];
  uint32_t singleBitV[2];
  int bitsX = 0, bitsY = 0;

  for (int i = 0; i < BVD_IBC_MAX_PREFIX; i++)
  {
    if (i < horCtxThre)
    {
      const BinFracBits fracBitsPar = fracBits.getFracBitsArray(Ctx::Bvd(horCtxOs + i + 1));
      singleBitH[0] = fracBitsPar.intBits[0];
      singleBitH[1] = fracBitsPar.intBits[1];
    }
    else
    {
      singleBitH[0] = epBitCost;
      singleBitH[1] = epBitCost;
    }
    p->bitsH[i] = bitsX + singleBitH[0] + (i+BVD_CODING_GOLOMB_ORDER) * epBitCost;
    bitsX += singleBitH[1];
  }

  for (int i = 0; i < BVD_IBC_MAX_PREFIX; i++)
  {
    if (i < verCtxThre)
    {
      const BinFracBits fracBitsPar = fracBits.getFracBitsArray(Ctx::Bvd(verCtxOs + i + 1));
      singleBitV[0] = fracBitsPar.intBits[0];
      singleBitV[1] = fracBitsPar.intBits[1];
    }
    else
    {
      singleBitV[0] = epBitCost;
      singleBitV[1] = epBitCost;
    }
    p->bitsV[i] = bitsY + singleBitV[0] + (i+BVD_CODING_GOLOMB_ORDER) * epBitCost;
    bitsY += singleBitV[1];
  }

  p->bitsIdx[0] = fracBits.getFracBitsArray(Ctx::MVPIdx()).intBits[0];
  p->bitsIdx[1] = fracBits.getFracBitsArray(Ctx::MVPIdx()).intBits[1];
  p->bitsImv[0] = fracBits.getFracBitsArray(Ctx::ImvFlag(1)).intBits[0];
  p->bitsImv[1] = fracBits.getFracBitsArray(Ctx::ImvFlag(1)).intBits[1];
}
#endif

InterSearch::InterSearch()
  : m_modeCtrl                    (nullptr)
  , m_pSplitCS                    (nullptr)
  , m_pFullCS                     (nullptr)
  , m_pcEncCfg                    (nullptr)
#if JVET_V0094_BILATERAL_FILTER || JVET_X0071_CHROMA_BILATERAL_FILTER
, m_bilateralFilter             (nullptr)
#endif
  , m_pcTrQuant                   (nullptr)
  , m_pcReshape                   (nullptr)
  , m_iSearchRange                (0)
  , m_bipredSearchRange           (0)
  , m_motionEstimationSearchMethod(MESEARCH_FULL)
  , m_CABACEstimator              (nullptr)
  , m_CtxCache                    (nullptr)
  , m_pTempPel                    (nullptr)
  , m_isInitialized               (false)
{
  for (int i=0; i<MAX_NUM_REF_LIST_ADAPT_SR; i++)
  {
    memset (m_aaiAdaptSR[i], 0, MAX_IDX_ADAPT_SR * sizeof (int));
  }
  for (int i=0; i<AMVP_MAX_NUM_CANDS+1; i++)
  {
    memset (m_auiMVPIdxCost[i], 0, (AMVP_MAX_NUM_CANDS+1) * sizeof (uint32_t) );
  }

  setWpScalingDistParam( -1, REF_PIC_LIST_X, nullptr );
  m_affMVList = nullptr;
  m_affMVListSize = 0;
  m_affMVListIdx = 0;
  m_uniMvList = nullptr;
  m_uniMvListSize = 0;
  m_uniMvListIdx = 0;
#if INTER_LIC
  m_uniMvListLIC = nullptr;
  m_uniMvListSizeLIC = 0;
  m_uniMvListIdxLIC = 0;
#endif
  m_histBestSbt    = MAX_UCHAR;
  m_histBestMtsIdx = MAX_UCHAR;

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  m_tplWeightTblInitialized = false;
  initTplWeightTable();
#endif
}


void InterSearch::destroy()
{
  CHECK(!m_isInitialized, "Not initialized");
  if ( m_pTempPel )
  {
    delete [] m_pTempPel;
    m_pTempPel = NULL;
  }

  m_pSplitCS = m_pFullCS = nullptr;

  m_pSaveCS = nullptr;

  for(uint32_t i = 0; i < NUM_REF_PIC_LIST_01; i++)
  {
    m_tmpPredStorage[i].destroy();
  }
  m_tmpStorageLCU.destroy();
  m_tmpAffiStorage.destroy();

  if ( m_tmpAffiError != NULL )
  {
    delete[] m_tmpAffiError;
  }
  if ( m_tmpAffiDeri[0] != NULL )
  {
    delete[] m_tmpAffiDeri[0];
  }
  if ( m_tmpAffiDeri[1] != NULL )
  {
    delete[] m_tmpAffiDeri[1];
  }
  if (m_affMVList)
  {
    delete[] m_affMVList;
    m_affMVList = nullptr;
  }
  m_affMVListIdx = 0;
  m_affMVListSize = 0;
  if (m_uniMvList)
  {
    delete[] m_uniMvList;
    m_uniMvList = nullptr;
  }
  m_uniMvListIdx = 0;
  m_uniMvListSize = 0;
#if INTER_LIC
  if (m_uniMvListLIC)
  {
    delete[] m_uniMvListLIC;
    m_uniMvListLIC = nullptr;
  }
  m_uniMvListIdxLIC = 0;
  m_uniMvListSizeLIC = 0;
#endif
  m_isInitialized = false;
}

void InterSearch::setTempBuffers( CodingStructure ****pSplitCS, CodingStructure ****pFullCS, CodingStructure **pSaveCS )
{
  m_pSplitCS = pSplitCS;
  m_pFullCS  = pFullCS;
  m_pSaveCS  = pSaveCS;
}

#if ENABLE_SPLIT_PARALLELISM
void InterSearch::copyState( const InterSearch& other )
{
  memcpy( m_aaiAdaptSR, other.m_aaiAdaptSR, sizeof( m_aaiAdaptSR ) );
}
#endif

InterSearch::~InterSearch()
{
  if (m_isInitialized)
  {
    destroy();
  }
}

void InterSearch::init( EncCfg*        pcEncCfg,
#if JVET_V0094_BILATERAL_FILTER || JVET_X0071_CHROMA_BILATERAL_FILTER
                      BilateralFilter* bilateralFilter,
#endif
                        TrQuant*       pcTrQuant,
                        int            iSearchRange,
                        int            bipredSearchRange,
                        MESearchMethod motionEstimationSearchMethod,
                        bool           useCompositeRef,
                        const uint32_t     maxCUWidth,
                        const uint32_t     maxCUHeight,
                        const uint32_t     maxTotalCUDepth,
                        RdCost*        pcRdCost,
                        CABACWriter*   CABACEstimator,
                        CtxCache*      ctxCache
                      , EncReshape*    pcReshape
#if JVET_Z0153_IBC_EXT_REF
                      , const uint32_t curPicWidthY 
#endif
)
{
  CHECK(m_isInitialized, "Already initialized");
  m_numBVs = 0;
  for (int i = 0; i < IBC_NUM_CANDIDATES; i++)
  {
    m_defaultCachedBvs.m_bvCands[i].setZero();
  }
  m_defaultCachedBvs.currCnt = 0;
  m_pcEncCfg                     = pcEncCfg;
#if JVET_V0094_BILATERAL_FILTER || JVET_X0071_CHROMA_BILATERAL_FILTER
  m_bilateralFilter              = bilateralFilter;
#endif
  m_pcTrQuant                    = pcTrQuant;
  m_iSearchRange                 = iSearchRange;
  m_bipredSearchRange            = bipredSearchRange;
  m_motionEstimationSearchMethod = motionEstimationSearchMethod;
  m_CABACEstimator               = CABACEstimator;
  m_CtxCache                     = ctxCache;
  m_useCompositeRef              = useCompositeRef;
  m_pcReshape                    = pcReshape;

  for( uint32_t iDir = 0; iDir < MAX_NUM_REF_LIST_ADAPT_SR; iDir++ )
  {
    for( uint32_t iRefIdx = 0; iRefIdx < MAX_IDX_ADAPT_SR; iRefIdx++ )
    {
      m_aaiAdaptSR[iDir][iRefIdx] = iSearchRange;
    }
  }

  // initialize motion cost
  for( int iNum = 0; iNum < AMVP_MAX_NUM_CANDS + 1; iNum++ )
  {
    for( int iIdx = 0; iIdx < AMVP_MAX_NUM_CANDS; iIdx++ )
    {
      if( iIdx < iNum )
      {
        m_auiMVPIdxCost[iIdx][iNum] = xGetMvpIdxBits( iIdx, iNum );
      }
      else
      {
        m_auiMVPIdxCost[iIdx][iNum] = MAX_UINT;
      }
    }
  }

  const ChromaFormat cform = pcEncCfg->getChromaFormatIdc();
#if INTER_LIC || (TM_AMVP || TM_MRG) || JVET_W0090_ARMC_TM || JVET_Z0056_GPM_SPLIT_MODE_REORDERING
#if JVET_Z0153_IBC_EXT_REF
  InterPrediction::init( pcRdCost, cform, maxCUHeight, m_pcReshape, curPicWidthY );
#else
  InterPrediction::init( pcRdCost, cform, maxCUHeight, m_pcReshape );
#endif
#else 
  InterPrediction::init( pcRdCost, cform, maxCUHeight );
#endif

  for( uint32_t i = 0; i < NUM_REF_PIC_LIST_01; i++ )
  {
    m_tmpPredStorage[i].create( UnitArea( cform, Area( 0, 0, MAX_CU_SIZE, MAX_CU_SIZE ) ) );
  }
  m_tmpStorageLCU.create( UnitArea( cform, Area( 0, 0, MAX_CU_SIZE, MAX_CU_SIZE ) ) );
  m_tmpAffiStorage.create(UnitArea(cform, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
  m_tmpAffiError = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#if AFFINE_ENC_OPT
  m_tmpAffiDeri[0] = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
  m_tmpAffiDeri[1] = new Pel[MAX_CU_SIZE * MAX_CU_SIZE];
#else
  m_tmpAffiDeri[0] = new int[MAX_CU_SIZE * MAX_CU_SIZE];
  m_tmpAffiDeri[1] = new int[MAX_CU_SIZE * MAX_CU_SIZE];
#endif
  m_pTempPel = new Pel[maxCUWidth*maxCUHeight];
  m_affMVListMaxSize = (pcEncCfg->getIntraPeriod() == (uint32_t)-1) ? AFFINE_ME_LIST_SIZE_LD : AFFINE_ME_LIST_SIZE;
  if (!m_affMVList)
    m_affMVList = new AffineMVInfo[m_affMVListMaxSize];
  m_affMVListIdx = 0;
  m_affMVListSize = 0;
  m_uniMvListMaxSize = 15;
  if (!m_uniMvList)
  {
    m_uniMvList = new BlkUniMvInfo[m_uniMvListMaxSize];
  }
  m_uniMvListIdx = 0;
  m_uniMvListSize = 0;
#if INTER_LIC
  if (!m_uniMvListLIC)
  {
    m_uniMvListLIC = new BlkUniMvInfo[m_uniMvListMaxSize];
  }
  m_uniMvListIdxLIC = 0;
  m_uniMvListSizeLIC = 0;
#endif
  m_isInitialized = true;
}

void InterSearch::resetSavedAffineMotion()
{
  for ( int i = 0; i < 2; i++ )
  {
    for ( int j = 0; j < 2; j++ )
    {
      m_affineMotion.acMvAffine4Para[i][j] = Mv( 0, 0 );
      m_affineMotion.acMvAffine6Para[i][j] = Mv( 0, 0 );
    }
    m_affineMotion.acMvAffine6Para[i][2] = Mv( 0, 0 );

    m_affineMotion.affine4ParaRefIdx[i] = -1;
    m_affineMotion.affine6ParaRefIdx[i] = -1;
  }
  for ( int i = 0; i < 3; i++ )
  {
    m_affineMotion.hevcCost[i] = std::numeric_limits<Distortion>::max();
  }
  m_affineMotion.affine4ParaAvail = false;
  m_affineMotion.affine6ParaAvail = false;
}

void InterSearch::storeAffineMotion( Mv acAffineMv[2][3], int8_t affineRefIdx[2], EAffineModel affineType, int bcwIdx )
{
  if ( ( bcwIdx == BCW_DEFAULT || !m_affineMotion.affine6ParaAvail ) && affineType == AFFINEMODEL_6PARAM )
  {
    for ( int i = 0; i < 2; i++ )
    {
      for ( int j = 0; j < 3; j++ )
      {
        m_affineMotion.acMvAffine6Para[i][j] = acAffineMv[i][j];
      }
      m_affineMotion.affine6ParaRefIdx[i] = affineRefIdx[i];
    }
    m_affineMotion.affine6ParaAvail = true;
  }

  if ( ( bcwIdx == BCW_DEFAULT || !m_affineMotion.affine4ParaAvail ) && affineType == AFFINEMODEL_4PARAM )
  {
    for ( int i = 0; i < 2; i++ )
    {
      for ( int j = 0; j < 2; j++ )
      {
        m_affineMotion.acMvAffine4Para[i][j] = acAffineMv[i][j];
      }
      m_affineMotion.affine4ParaRefIdx[i] = affineRefIdx[i];
    }
    m_affineMotion.affine4ParaAvail = true;
  }
}

inline void InterSearch::xTZSearchHelp( IntTZSearchStruct& rcStruct, const int iSearchX, const int iSearchY, const uint8_t ucPointNr, const uint32_t uiDistance )
{
  Distortion  uiSad = 0;

//  CHECK(!( !( rcStruct.searchRange.left > iSearchX || rcStruct.searchRange.right < iSearchX || rcStruct.searchRange.top > iSearchY || rcStruct.searchRange.bottom < iSearchY )), "Unspecified error");

  const Pel* const  piRefSrch = rcStruct.piRefY + iSearchY * rcStruct.iRefStride + iSearchX;

  m_cDistParam.cur.buf = piRefSrch;

  if( 1 == rcStruct.subShiftMode )
  {
    // motion cost
    Distortion uiBitCost = m_pcRdCost->getCostOfVectorWithPredictor( iSearchX, iSearchY, rcStruct.imvShift );

    // Skip search if bit cost is already larger than best SAD
    if (uiBitCost < rcStruct.uiBestSad)
    {
      Distortion uiTempSad = m_cDistParam.distFunc( m_cDistParam );

      if((uiTempSad + uiBitCost) < rcStruct.uiBestSad)
      {
        // it's not supposed that any member of DistParams is manipulated beside cur.buf
        int subShift = m_cDistParam.subShift;
        const Pel* pOrgCpy = m_cDistParam.org.buf;
        uiSad += uiTempSad >> m_cDistParam.subShift;

        while( m_cDistParam.subShift > 0 )
        {
          int isubShift           = m_cDistParam.subShift -1;
          m_cDistParam.org.buf = rcStruct.pcPatternKey->buf + (rcStruct.pcPatternKey->stride << isubShift);
          m_cDistParam.cur.buf = piRefSrch + (rcStruct.iRefStride << isubShift);
          uiTempSad            = m_cDistParam.distFunc( m_cDistParam );
          uiSad               += uiTempSad >> m_cDistParam.subShift;

          if(((uiSad << isubShift) + uiBitCost) > rcStruct.uiBestSad)
          {
            break;
          }

          m_cDistParam.subShift--;
        }

        if(m_cDistParam.subShift == 0)
        {
          uiSad += uiBitCost;

          if( uiSad < rcStruct.uiBestSad )
          {
            rcStruct.uiBestSad      = uiSad;
            rcStruct.iBestX         = iSearchX;
            rcStruct.iBestY         = iSearchY;
            rcStruct.uiBestDistance = uiDistance;
            rcStruct.uiBestRound    = 0;
            rcStruct.ucPointNr      = ucPointNr;
            m_cDistParam.maximumDistortionForEarlyExit = uiSad;
          }
        }

        // restore org ptr
        m_cDistParam.org.buf  = pOrgCpy;
        m_cDistParam.subShift = subShift;
      }
    }
  }
  else
  {
    uiSad = m_cDistParam.distFunc( m_cDistParam );

    // only add motion cost if uiSad is smaller than best. Otherwise pointless
    // to add motion cost.
    if( uiSad < rcStruct.uiBestSad )
    {
      // motion cost
      uiSad += m_pcRdCost->getCostOfVectorWithPredictor( iSearchX, iSearchY, rcStruct.imvShift );

      if( uiSad < rcStruct.uiBestSad )
      {
        rcStruct.uiBestSad      = uiSad;
        rcStruct.iBestX         = iSearchX;
        rcStruct.iBestY         = iSearchY;
        rcStruct.uiBestDistance = uiDistance;
        rcStruct.uiBestRound    = 0;
        rcStruct.ucPointNr      = ucPointNr;
        m_cDistParam.maximumDistortionForEarlyExit = uiSad;
      }
    }
  }
}



inline void InterSearch::xTZ2PointSearch( IntTZSearchStruct& rcStruct )
{
  const SearchRange& sr = rcStruct.searchRange;

  static const int xOffset[2][9] = { {  0, -1, -1,  0, -1, +1, -1, -1, +1 }, {  0,  0, +1, +1, -1, +1,  0, +1,  0 } };
  static const int yOffset[2][9] = { {  0,  0, -1, -1, +1, -1,  0, +1,  0 }, {  0, -1, -1,  0, -1, +1, +1, +1, +1 } };

  // 2 point search,                   //   1 2 3
  // check only the 2 untested points  //   4 0 5
  // around the start point            //   6 7 8
  const int iX1 = rcStruct.iBestX + xOffset[0][rcStruct.ucPointNr];
  const int iX2 = rcStruct.iBestX + xOffset[1][rcStruct.ucPointNr];

  const int iY1 = rcStruct.iBestY + yOffset[0][rcStruct.ucPointNr];
  const int iY2 = rcStruct.iBestY + yOffset[1][rcStruct.ucPointNr];

  if( iX1 >= sr.left && iX1 <= sr.right && iY1 >= sr.top && iY1 <= sr.bottom )
  {
    xTZSearchHelp( rcStruct, iX1, iY1, 0, 2 );
  }

  if( iX2 >= sr.left && iX2 <= sr.right && iY2 >= sr.top && iY2 <= sr.bottom )
  {
    xTZSearchHelp( rcStruct, iX2, iY2, 0, 2 );
  }
}


inline void InterSearch::xTZ8PointSquareSearch( IntTZSearchStruct& rcStruct, const int iStartX, const int iStartY, const int iDist )
{
  const SearchRange& sr = rcStruct.searchRange;
  // 8 point search,                   //   1 2 3
  // search around the start point     //   4 0 5
  // with the required  distance       //   6 7 8
  CHECK( iDist == 0 , "Invalid distance");
  const int iTop        = iStartY - iDist;
  const int iBottom     = iStartY + iDist;
  const int iLeft       = iStartX - iDist;
  const int iRight      = iStartX + iDist;
  rcStruct.uiBestRound += 1;

  if ( iTop >= sr.top ) // check top
  {
    if ( iLeft >= sr.left ) // check top left
    {
      xTZSearchHelp( rcStruct, iLeft, iTop, 1, iDist );
    }
    // top middle
    xTZSearchHelp( rcStruct, iStartX, iTop, 2, iDist );

    if ( iRight <= sr.right ) // check top right
    {
      xTZSearchHelp( rcStruct, iRight, iTop, 3, iDist );
    }
  } // check top
  if ( iLeft >= sr.left ) // check middle left
  {
    xTZSearchHelp( rcStruct, iLeft, iStartY, 4, iDist );
  }
  if ( iRight <= sr.right ) // check middle right
  {
    xTZSearchHelp( rcStruct, iRight, iStartY, 5, iDist );
  }
  if ( iBottom <= sr.bottom ) // check bottom
  {
    if ( iLeft >= sr.left ) // check bottom left
    {
      xTZSearchHelp( rcStruct, iLeft, iBottom, 6, iDist );
    }
    // check bottom middle
    xTZSearchHelp( rcStruct, iStartX, iBottom, 7, iDist );

    if ( iRight <= sr.right ) // check bottom right
    {
      xTZSearchHelp( rcStruct, iRight, iBottom, 8, iDist );
    }
  } // check bottom
}




inline void InterSearch::xTZ8PointDiamondSearch( IntTZSearchStruct& rcStruct,
                                                 const int iStartX,
                                                 const int iStartY,
                                                 const int iDist,
                                                 const bool bCheckCornersAtDist1 )
{
  const SearchRange& sr = rcStruct.searchRange;
  // 8 point search,                   //   1 2 3
  // search around the start point     //   4 0 5
  // with the required  distance       //   6 7 8
  CHECK( iDist == 0, "Invalid distance" );
  const int iTop        = iStartY - iDist;
  const int iBottom     = iStartY + iDist;
  const int iLeft       = iStartX - iDist;
  const int iRight      = iStartX + iDist;
  rcStruct.uiBestRound += 1;

  if ( iDist == 1 )
  {
    if ( iTop >= sr.top ) // check top
    {
      if (bCheckCornersAtDist1)
      {
        if ( iLeft >= sr.left) // check top-left
        {
          xTZSearchHelp( rcStruct, iLeft, iTop, 1, iDist );
        }
        xTZSearchHelp( rcStruct, iStartX, iTop, 2, iDist );
        if ( iRight <= sr.right ) // check middle right
        {
          xTZSearchHelp( rcStruct, iRight, iTop, 3, iDist );
        }
      }
      else
      {
        xTZSearchHelp( rcStruct, iStartX, iTop, 2, iDist );
      }
    }
    if ( iLeft >= sr.left ) // check middle left
    {
      xTZSearchHelp( rcStruct, iLeft, iStartY, 4, iDist );
    }
    if ( iRight <= sr.right ) // check middle right
    {
      xTZSearchHelp( rcStruct, iRight, iStartY, 5, iDist );
    }
    if ( iBottom <= sr.bottom ) // check bottom
    {
      if (bCheckCornersAtDist1)
      {
        if ( iLeft >= sr.left) // check top-left
        {
          xTZSearchHelp( rcStruct, iLeft, iBottom, 6, iDist );
        }
        xTZSearchHelp( rcStruct, iStartX, iBottom, 7, iDist );
        if ( iRight <= sr.right ) // check middle right
        {
          xTZSearchHelp( rcStruct, iRight, iBottom, 8, iDist );
        }
      }
      else
      {
        xTZSearchHelp( rcStruct, iStartX, iBottom, 7, iDist );
      }
    }
  }
  else
  {
    if ( iDist <= 8 )
    {
      const int iTop_2      = iStartY - (iDist>>1);
      const int iBottom_2   = iStartY + (iDist>>1);
      const int iLeft_2     = iStartX - (iDist>>1);
      const int iRight_2    = iStartX + (iDist>>1);

      if (  iTop >= sr.top && iLeft >= sr.left &&
           iRight <= sr.right && iBottom <= sr.bottom ) // check border
      {
        xTZSearchHelp( rcStruct, iStartX,  iTop,      2, iDist    );
        xTZSearchHelp( rcStruct, iLeft_2,  iTop_2,    1, iDist>>1 );
        xTZSearchHelp( rcStruct, iRight_2, iTop_2,    3, iDist>>1 );
        xTZSearchHelp( rcStruct, iLeft,    iStartY,   4, iDist    );
        xTZSearchHelp( rcStruct, iRight,   iStartY,   5, iDist    );
        xTZSearchHelp( rcStruct, iLeft_2,  iBottom_2, 6, iDist>>1 );
        xTZSearchHelp( rcStruct, iRight_2, iBottom_2, 8, iDist>>1 );
        xTZSearchHelp( rcStruct, iStartX,  iBottom,   7, iDist    );
      }
      else // check border
      {
        if ( iTop >= sr.top ) // check top
        {
          xTZSearchHelp( rcStruct, iStartX, iTop, 2, iDist );
        }
        if ( iTop_2 >= sr.top ) // check half top
        {
          if ( iLeft_2 >= sr.left ) // check half left
          {
            xTZSearchHelp( rcStruct, iLeft_2, iTop_2, 1, (iDist>>1) );
          }
          if ( iRight_2 <= sr.right ) // check half right
          {
            xTZSearchHelp( rcStruct, iRight_2, iTop_2, 3, (iDist>>1) );
          }
        } // check half top
        if ( iLeft >= sr.left ) // check left
        {
          xTZSearchHelp( rcStruct, iLeft, iStartY, 4, iDist );
        }
        if ( iRight <= sr.right ) // check right
        {
          xTZSearchHelp( rcStruct, iRight, iStartY, 5, iDist );
        }
        if ( iBottom_2 <= sr.bottom ) // check half bottom
        {
          if ( iLeft_2 >= sr.left ) // check half left
          {
            xTZSearchHelp( rcStruct, iLeft_2, iBottom_2, 6, (iDist>>1) );
          }
          if ( iRight_2 <= sr.right ) // check half right
          {
            xTZSearchHelp( rcStruct, iRight_2, iBottom_2, 8, (iDist>>1) );
          }
        } // check half bottom
        if ( iBottom <= sr.bottom ) // check bottom
        {
          xTZSearchHelp( rcStruct, iStartX, iBottom, 7, iDist );
        }
      } // check border
    }
    else // iDist > 8
    {
      if ( iTop >= sr.top && iLeft >= sr.left &&
           iRight <= sr.right && iBottom <= sr.bottom ) // check border
      {
        xTZSearchHelp( rcStruct, iStartX, iTop,    0, iDist );
        xTZSearchHelp( rcStruct, iLeft,   iStartY, 0, iDist );
        xTZSearchHelp( rcStruct, iRight,  iStartY, 0, iDist );
        xTZSearchHelp( rcStruct, iStartX, iBottom, 0, iDist );
        for ( int index = 1; index < 4; index++ )
        {
          const int iPosYT = iTop    + ((iDist>>2) * index);
          const int iPosYB = iBottom - ((iDist>>2) * index);
          const int iPosXL = iStartX - ((iDist>>2) * index);
          const int iPosXR = iStartX + ((iDist>>2) * index);
          xTZSearchHelp( rcStruct, iPosXL, iPosYT, 0, iDist );
          xTZSearchHelp( rcStruct, iPosXR, iPosYT, 0, iDist );
          xTZSearchHelp( rcStruct, iPosXL, iPosYB, 0, iDist );
          xTZSearchHelp( rcStruct, iPosXR, iPosYB, 0, iDist );
        }
      }
      else // check border
      {
        if ( iTop >= sr.top ) // check top
        {
          xTZSearchHelp( rcStruct, iStartX, iTop, 0, iDist );
        }
        if ( iLeft >= sr.left ) // check left
        {
          xTZSearchHelp( rcStruct, iLeft, iStartY, 0, iDist );
        }
        if ( iRight <= sr.right ) // check right
        {
          xTZSearchHelp( rcStruct, iRight, iStartY, 0, iDist );
        }
        if ( iBottom <= sr.bottom ) // check bottom
        {
          xTZSearchHelp( rcStruct, iStartX, iBottom, 0, iDist );
        }
        for ( int index = 1; index < 4; index++ )
        {
          const int iPosYT = iTop    + ((iDist>>2) * index);
          const int iPosYB = iBottom - ((iDist>>2) * index);
          const int iPosXL = iStartX - ((iDist>>2) * index);
          const int iPosXR = iStartX + ((iDist>>2) * index);

          if ( iPosYT >= sr.top ) // check top
          {
            if ( iPosXL >= sr.left ) // check left
            {
              xTZSearchHelp( rcStruct, iPosXL, iPosYT, 0, iDist );
            }
            if ( iPosXR <= sr.right ) // check right
            {
              xTZSearchHelp( rcStruct, iPosXR, iPosYT, 0, iDist );
            }
          } // check top
          if ( iPosYB <= sr.bottom ) // check bottom
          {
            if ( iPosXL >= sr.left ) // check left
            {
              xTZSearchHelp( rcStruct, iPosXL, iPosYB, 0, iDist );
            }
            if ( iPosXR <= sr.right ) // check right
            {
              xTZSearchHelp( rcStruct, iPosXR, iPosYB, 0, iDist );
            }
          } // check bottom
        } // for ...
      } // check border
    } // iDist <= 8
  } // iDist == 1
}

Distortion InterSearch::xPatternRefinement( const CPelBuf* pcPatternKey,
                                            Mv baseRefMv,
                                            int iFrac, Mv& rcMvFrac,
                                            bool bAllowUseOfHadamard )
{
  Distortion  uiDist;
  Distortion  uiDistBest  = std::numeric_limits<Distortion>::max();
  uint32_t        uiDirecBest = 0;

  Pel*  piRefPos;
  int iRefStride = pcPatternKey->width + 1;
  m_pcRdCost->setDistParam( m_cDistParam, *pcPatternKey, m_filteredBlock[0][0][0], iRefStride, m_lumaClpRng.bd, COMPONENT_Y, 0, 1, m_pcEncCfg->getUseHADME() && bAllowUseOfHadamard );

  const Mv* pcMvRefine = (iFrac == 2 ? s_acMvRefineH : s_acMvRefineQ);
  for (uint32_t i = 0; i < 9; i++)
  {
    if (m_skipFracME && i > 0)
    {
      break;
    }
    Mv cMvTest = pcMvRefine[i];
    cMvTest += baseRefMv;

    int horVal = cMvTest.getHor() * iFrac;
    int verVal = cMvTest.getVer() * iFrac;
    piRefPos = m_filteredBlock[verVal & 3][horVal & 3][0];

    if (horVal == 2 && (verVal & 1) == 0)
    {
      piRefPos += 1;
    }
    if ((horVal & 1) == 0 && verVal == 2)
    {
      piRefPos += iRefStride;
    }
    cMvTest = pcMvRefine[i];
    cMvTest += rcMvFrac;


    m_cDistParam.cur.buf   = piRefPos;
    uiDist = m_cDistParam.distFunc( m_cDistParam );
    uiDist += m_pcRdCost->getCostOfVectorWithPredictor( cMvTest.getHor(), cMvTest.getVer(), 0 );

    if ( uiDist < uiDistBest )
    {
      uiDistBest  = uiDist;
      uiDirecBest = i;
      m_cDistParam.maximumDistortionForEarlyExit = uiDist;
    }
  }

  rcMvFrac = pcMvRefine[uiDirecBest];

  return uiDistBest;
}

Distortion InterSearch::xGetInterPredictionError( PredictionUnit& pu, PelUnitBuf& origBuf, const RefPicList &eRefPicList )
{
  PelUnitBuf predBuf = m_tmpStorageLCU.getBuf( UnitAreaRelative(*pu.cu, pu) );

  motionCompensation( pu, predBuf, eRefPicList );

  DistParam cDistParam;
  cDistParam.applyWeight = false;

  m_pcRdCost->setDistParam(cDistParam, origBuf.Y(), predBuf.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, m_pcEncCfg->getUseHADME() && !pu.cu->slice->getDisableSATDForRD());

  return (Distortion)cDistParam.distFunc( cDistParam );
}

/// add ibc search functions here

void InterSearch::xIBCSearchMVCandUpdate(Distortion  sad, int x, int y, Distortion* sadBestCand, Mv* cMVCand)
{
  int j = CHROMA_REFINEMENT_CANDIDATES - 1;

  if (sad < sadBestCand[CHROMA_REFINEMENT_CANDIDATES - 1])
  {
    for (int t = CHROMA_REFINEMENT_CANDIDATES - 1; t >= 0; t--)
    {
      if (sad < sadBestCand[t])
        j = t;
    }

    for (int k = CHROMA_REFINEMENT_CANDIDATES - 1; k > j; k--)
    {
      sadBestCand[k] = sadBestCand[k - 1];

      cMVCand[k].set(cMVCand[k - 1].getHor(), cMVCand[k - 1].getVer());
    }
    sadBestCand[j] = sad;
    cMVCand[j].set(x, y);
  }
}

int InterSearch::xIBCSearchMVChromaRefine(PredictionUnit& pu,
  int         roiWidth,
  int         roiHeight,
  int         cuPelX,
  int         cuPelY,
  Distortion* sadBestCand,
  Mv*     cMVCand

)
{
  if ( (!isChromaEnabled(pu.chromaFormat)) || (!pu.Cb().valid()) )
  {
    return 0;
  }

  int bestCandIdx = 0;
  Distortion  sadBest = std::numeric_limits<Distortion>::max();
  Distortion  tempSad;

  Pel* pRef;
  Pel* pOrg;
  int refStride, orgStride;
  int width, height;

  int picWidth = pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
  int picHeight = pu.cs->slice->getPPS()->getPicHeightInLumaSamples();

  UnitArea allCompBlocks(pu.chromaFormat, (Area)pu.block(COMPONENT_Y));
  for (int cand = 0; cand < CHROMA_REFINEMENT_CANDIDATES; cand++)
  {
    if (sadBestCand[cand] == std::numeric_limits<Distortion>::max())
    {
      continue;
    }

    if ((!cMVCand[cand].getHor()) && (!cMVCand[cand].getVer()))
      continue;

    if (((int)(cuPelY + cMVCand[cand].getVer() + roiHeight) >= picHeight) || ((cuPelY + cMVCand[cand].getVer()) < 0))
      continue;

    if (((int)(cuPelX + cMVCand[cand].getHor() + roiWidth) >= picWidth) || ((cuPelX + cMVCand[cand].getHor()) < 0))
      continue;

    tempSad = sadBestCand[cand];

    pu.mv[0] = cMVCand[cand];
    pu.mv[0].changePrecision(MV_PRECISION_INT, MV_PRECISION_INTERNAL);
    pu.interDir = 1;
    pu.refIdx[0] = pu.cs->slice->getNumRefIdx(REF_PIC_LIST_0); // last idx in the list

    PelUnitBuf predBufTmp = m_tmpPredStorage[REF_PIC_LIST_0].getBuf(UnitAreaRelative(*pu.cu, pu));
    motionCompensation(pu, predBufTmp, REF_PIC_LIST_0);

    for (unsigned int ch = COMPONENT_Cb; ch < ::getNumberValidComponents(pu.chromaFormat); ch++)
    {
      width = roiWidth >> ::getComponentScaleX(ComponentID(ch), pu.chromaFormat);
      height = roiHeight >> ::getComponentScaleY(ComponentID(ch), pu.chromaFormat);

      PelUnitBuf origBuf = pu.cs->getOrgBuf(allCompBlocks);
      PelUnitBuf* pBuf = &origBuf;
      CPelBuf  tmpPattern = pBuf->get(ComponentID(ch));
      pOrg = (Pel*)tmpPattern.buf;

      Picture* refPic = pu.cu->slice->getPic();
      const CPelBuf refBuf = refPic->getRecoBuf(allCompBlocks.blocks[ComponentID(ch)]);
      pRef = (Pel*)refBuf.buf;

      refStride = refBuf.stride;
      orgStride = tmpPattern.stride;

      //ComponentID compID = (ComponentID)ch;
      PelUnitBuf* pBufRef = &predBufTmp;
      CPelBuf  tmpPatternRef = pBufRef->get(ComponentID(ch));
      pRef = (Pel*)tmpPatternRef.buf;
      refStride = tmpPatternRef.stride;


      for (int row = 0; row < height; row++)
      {
        for (int col = 0; col < width; col++)
        {
          tempSad += ((abs(pRef[col] - pOrg[col])) >> (pu.cs->sps->getBitDepth(CHANNEL_TYPE_CHROMA) - 8));
        }
        pRef += refStride;
        pOrg += orgStride;
      }
    }

    if (tempSad < sadBest)
    {
      sadBest = tempSad;
      bestCandIdx = cand;
    }
  }

  return bestCandIdx;
}

static unsigned int xMergeCandLists(Mv *dst, unsigned int dn, unsigned int dstTotalLength, Mv *src, unsigned int sn)
{
  for (unsigned int cand = 0; cand < sn && dn < dstTotalLength; cand++)
  {
    if (src[cand] == Mv())
    {
      continue;
    }
    bool found = false;
    for (int j = 0; j<dn; j++)
    {
      if (src[cand] == dst[j])
      {
        found = true;
        break;
      }
    }

    if (!found)
    {
      dst[dn] = src[cand];
      dn++;
    }
  }

  return dn;
}

void InterSearch::xIntraPatternSearch(PredictionUnit& pu, IntTZSearchStruct&  cStruct, Mv& rcMv, Distortion&  ruiCost, Mv*  pcMvSrchRngLT, Mv*  pcMvSrchRngRB, Mv* pcMvPred)
{
  const int   srchRngHorLeft = pcMvSrchRngLT->getHor();
  const int   srchRngHorRight = pcMvSrchRngRB->getHor();
  const int   srchRngVerTop = pcMvSrchRngLT->getVer();
  const int   srchRngVerBottom = pcMvSrchRngRB->getVer();

  const unsigned int  lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
  const int   puPelOffsetX = 0;
  const int   puPelOffsetY = 0;
  const int   cuPelX = pu.Y().x;
  const int   cuPelY = pu.Y().y;

  int          roiWidth = pu.lwidth();
  int          roiHeight = pu.lheight();

  Distortion  sad;
  Distortion  sadBest = std::numeric_limits<Distortion>::max();
  int         bestX = 0;
  int         bestY = 0;

  const Pel*        piRefSrch = cStruct.piRefY;

  int         bestCandIdx = 0;

  Distortion  sadBestCand[CHROMA_REFINEMENT_CANDIDATES];
  Mv      cMVCand[CHROMA_REFINEMENT_CANDIDATES];


  for (int cand = 0; cand < CHROMA_REFINEMENT_CANDIDATES; cand++)
  {
    sadBestCand[cand] = std::numeric_limits<Distortion>::max();
    cMVCand[cand].set(0, 0);
  }

  m_cDistParam.useMR = false;
  m_pcRdCost->setDistParam(m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, cStruct.subShiftMode);

  const int picWidth = pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
  const int picHeight = pu.cs->slice->getPPS()->getPicHeightInLumaSamples();


  {
    m_cDistParam.subShift = 0;

    Distortion tempSadBest = 0;

    int srLeft = srchRngHorLeft, srRight = srchRngHorRight, srTop = srchRngVerTop, srBottom = srchRngVerBottom;
    m_numBVs = 0;
    m_numBVs = xMergeCandLists(m_acBVs, m_numBVs, (2 * IBC_NUM_CANDIDATES), m_defaultCachedBvs.m_bvCands, m_defaultCachedBvs.currCnt);

    Mv cMvPredEncOnly[IBC_NUM_CANDIDATES];
    int nbPreds = 0;
    PU::getIbcMVPsEncOnly(pu, cMvPredEncOnly, nbPreds);
    m_numBVs = xMergeCandLists(m_acBVs, m_numBVs, (2 * IBC_NUM_CANDIDATES), cMvPredEncOnly, nbPreds);

    for (unsigned int cand = 0; cand < m_numBVs; cand++)
    {
      int xPred = m_acBVs[cand].getHor();
      int yPred = m_acBVs[cand].getVer();

      if (!(xPred == 0 && yPred == 0)
        && !((yPred < srTop) || (yPred > srBottom))
        && !((xPred < srLeft) || (xPred > srRight)))
      {
#if JVET_Z0084_IBC_TM
        bool validCand = PU::searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, xPred, yPred, lcuWidth);
#else
        bool validCand = searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, xPred, yPred, lcuWidth);
#endif

        if (validCand)
        {
          sad = m_pcRdCost->getBvCostMultiplePreds(xPred, yPred, pu.cs->sps->getAMVREnabledFlag());
          m_cDistParam.cur.buf = piRefSrch + cStruct.iRefStride * yPred + xPred;
          sad += m_cDistParam.distFunc(m_cDistParam);

          xIBCSearchMVCandUpdate(sad, xPred, yPred, sadBestCand, cMVCand);
        }
      }
    }

    bestX = cMVCand[0].getHor();
    bestY = cMVCand[0].getVer();
    rcMv.set(bestX, bestY);
    sadBest = sadBestCand[0];

    const int boundY = (0 - roiHeight - puPelOffsetY);
    for (int y = std::max(srchRngVerTop, 0 - cuPelY); y <= boundY; ++y)
    {
#if JVET_Z0084_IBC_TM
      if (!PU::searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, 0, y, lcuWidth))
#else
      if (!searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, 0, y, lcuWidth))
#endif
      {
        continue;
      }

      sad = m_pcRdCost->getBvCostMultiplePreds(0, y, pu.cs->sps->getAMVREnabledFlag());
      m_cDistParam.cur.buf = piRefSrch + cStruct.iRefStride * y;
      sad += m_cDistParam.distFunc(m_cDistParam);

      xIBCSearchMVCandUpdate(sad, 0, y, sadBestCand, cMVCand);
      tempSadBest = sadBestCand[0];
      if (sadBestCand[0] <= 3)
      {
        bestX = cMVCand[0].getHor();
        bestY = cMVCand[0].getVer();
        sadBest = sadBestCand[0];
        rcMv.set(bestX, bestY);
        ruiCost = sadBest;
        goto end;
      }
    }

    const int boundX = std::max(srchRngHorLeft, -cuPelX);
    for (int x = 0 - roiWidth - puPelOffsetX; x >= boundX; --x)
    {
#if JVET_Z0084_IBC_TM
      if (!PU::searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, 0, lcuWidth))
#else
      if (!searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, 0, lcuWidth))
#endif
      {
        continue;
      }

      sad = m_pcRdCost->getBvCostMultiplePreds(x, 0, pu.cs->sps->getAMVREnabledFlag());
      m_cDistParam.cur.buf = piRefSrch + x;
      sad += m_cDistParam.distFunc(m_cDistParam);


      xIBCSearchMVCandUpdate(sad, x, 0, sadBestCand, cMVCand);
      tempSadBest = sadBestCand[0];
      if (sadBestCand[0] <= 3)
      {
        bestX = cMVCand[0].getHor();
        bestY = cMVCand[0].getVer();
        sadBest = sadBestCand[0];
        rcMv.set(bestX, bestY);
        ruiCost = sadBest;
        goto end;
      }
    }

    bestX = cMVCand[0].getHor();
    bestY = cMVCand[0].getVer();
    sadBest = sadBestCand[0];
    if ((!bestX && !bestY) || (sadBest - m_pcRdCost->getBvCostMultiplePreds(bestX, bestY, pu.cs->sps->getAMVREnabledFlag()) <= 32))
    {
      //chroma refine
      bestCandIdx = xIBCSearchMVChromaRefine(pu, roiWidth, roiHeight, cuPelX, cuPelY, sadBestCand, cMVCand);
      bestX = cMVCand[bestCandIdx].getHor();
      bestY = cMVCand[bestCandIdx].getVer();
      sadBest = sadBestCand[bestCandIdx];
      rcMv.set(bestX, bestY);
      ruiCost = sadBest;
      goto end;
    }


    if (pu.lwidth() < 16 && pu.lheight() < 16)
    {
#if JVET_Z0153_IBC_EXT_REF
      int verTop = - (int)lcuWidth;
      int verBottom = std::min((int)(lcuWidth>>2), (int)(lcuWidth - (cuPelY % lcuWidth) - roiHeight));
      int horLeft = - (int)lcuWidth*2;
      int horRight = lcuWidth>>2;

      for (int y = std::max(verTop, -cuPelY); y <= verBottom; y += 2)
      {
        if ((y == 0) || ((int)(cuPelY + y + roiHeight) >= picHeight))
        {
          continue;
        }        
        for (int x = std::max(horLeft, -cuPelX); x <= horRight; x++)
        {
#else
      for (int y = std::max(srchRngVerTop, -cuPelY); y <= srchRngVerBottom; y += 2)
      {
        if ((y == 0) || ((int)(cuPelY + y + roiHeight) >= picHeight))
          continue;

        for (int x = std::max(srchRngHorLeft, -cuPelX); x <= srchRngHorRight; x++)
        {
#endif
          if ((x == 0) || ((int)(cuPelX + x + roiWidth) >= picWidth))
            continue;

#if JVET_Z0084_IBC_TM
          if (!PU::searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, y, lcuWidth))
#else
          if (!searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, y, lcuWidth))
#endif
          {
            continue;
          }

          sad = m_pcRdCost->getBvCostMultiplePreds(x, y, pu.cs->sps->getAMVREnabledFlag());
          m_cDistParam.cur.buf = piRefSrch + cStruct.iRefStride * y + x;
          sad += m_cDistParam.distFunc(m_cDistParam);

          xIBCSearchMVCandUpdate(sad, x, y, sadBestCand, cMVCand);
        }
      }

      bestX = cMVCand[0].getHor();
      bestY = cMVCand[0].getVer();
      sadBest = sadBestCand[0];
      if (sadBest - m_pcRdCost->getBvCostMultiplePreds(bestX, bestY, pu.cs->sps->getAMVREnabledFlag()) <= 16)
      {
        //chroma refine
        bestCandIdx = xIBCSearchMVChromaRefine(pu, roiWidth, roiHeight, cuPelX, cuPelY, sadBestCand, cMVCand);

        bestX = cMVCand[bestCandIdx].getHor();
        bestY = cMVCand[bestCandIdx].getVer();
        sadBest = sadBestCand[bestCandIdx];
        rcMv.set(bestX, bestY);
        ruiCost = sadBest;
        goto end;
      }

#if JVET_Z0153_IBC_EXT_REF
      for (int y = (std::max(verTop, -cuPelY) + 1); y <= verBottom; y += 2)
      {
        if ((y == 0) || ((int)(cuPelY + y + roiHeight) >= picHeight))
        {
          continue;
        }

        for (int x = std::max(horLeft, -cuPelX); x <= horRight; x += 2)
#else
      for (int y = (std::max(srchRngVerTop, -cuPelY) + 1); y <= srchRngVerBottom; y += 2)
      {
        if ((y == 0) || ((int)(cuPelY + y + roiHeight) >= picHeight))
          continue;

        for (int x = std::max(srchRngHorLeft, -cuPelX); x <= srchRngHorRight; x += 2)
#endif
        {
          if ((x == 0) || ((int)(cuPelX + x + roiWidth) >= picWidth))
            continue;

#if JVET_Z0084_IBC_TM
          if (!PU::searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, y, lcuWidth))
#else
          if (!searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, y, lcuWidth))
#endif
          {
            continue;
          }

          sad = m_pcRdCost->getBvCostMultiplePreds(x, y, pu.cs->sps->getAMVREnabledFlag());
          m_cDistParam.cur.buf = piRefSrch + cStruct.iRefStride * y + x;
          sad += m_cDistParam.distFunc(m_cDistParam);


          xIBCSearchMVCandUpdate(sad, x, y, sadBestCand, cMVCand);
          if (sadBestCand[0] <= 5)
          {
            //chroma refine & return
            bestCandIdx = xIBCSearchMVChromaRefine(pu, roiWidth, roiHeight, cuPelX, cuPelY, sadBestCand, cMVCand);
            bestX = cMVCand[bestCandIdx].getHor();
            bestY = cMVCand[bestCandIdx].getVer();
            sadBest = sadBestCand[bestCandIdx];
            rcMv.set(bestX, bestY);
            ruiCost = sadBest;
            goto end;
          }
        }
      }

      bestX = cMVCand[0].getHor();
      bestY = cMVCand[0].getVer();
      sadBest = sadBestCand[0];

      if ((sadBest >= tempSadBest) || ((sadBest - m_pcRdCost->getBvCostMultiplePreds(bestX, bestY, pu.cs->sps->getAMVREnabledFlag())) <= 32))
      {
        //chroma refine
        bestCandIdx = xIBCSearchMVChromaRefine(pu, roiWidth, roiHeight, cuPelX, cuPelY, sadBestCand, cMVCand);
        bestX = cMVCand[bestCandIdx].getHor();
        bestY = cMVCand[bestCandIdx].getVer();
        sadBest = sadBestCand[bestCandIdx];
        rcMv.set(bestX, bestY);
        ruiCost = sadBest;
        goto end;
      }

      tempSadBest = sadBestCand[0];

#if JVET_Z0153_IBC_EXT_REF
      for (int y = (std::max(verTop, -cuPelY) + 1); y <= verBottom; y += 2)
      {
        if ((y == 0) || ((int)(cuPelY + y + roiHeight) >= picHeight))
        {
          continue;
        }
        for (int x = (std::max(horLeft, -cuPelX) + 1); x <= horRight; x += 2)
#else
      for (int y = (std::max(srchRngVerTop, -cuPelY) + 1); y <= srchRngVerBottom; y += 2)
      {
        if ((y == 0) || ((int)(cuPelY + y + roiHeight) >= picHeight))
          continue;

        for (int x = (std::max(srchRngHorLeft, -cuPelX) + 1); x <= srchRngHorRight; x += 2)
#endif
        {

          if ((x == 0) || ((int)(cuPelX + x + roiWidth) >= picWidth))
            continue;

#if JVET_Z0084_IBC_TM
          if (!PU::searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, y, lcuWidth))
#else
          if (!searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, x, y, lcuWidth))
#endif
          {
            continue;
          }

          sad = m_pcRdCost->getBvCostMultiplePreds(x, y, pu.cs->sps->getAMVREnabledFlag());
          m_cDistParam.cur.buf = piRefSrch + cStruct.iRefStride * y + x;
          sad += m_cDistParam.distFunc(m_cDistParam);


          xIBCSearchMVCandUpdate(sad, x, y, sadBestCand, cMVCand);
          if (sadBestCand[0] <= 5)
          {
            //chroma refine & return
            bestCandIdx = xIBCSearchMVChromaRefine(pu, roiWidth, roiHeight, cuPelX, cuPelY, sadBestCand, cMVCand);
            bestX = cMVCand[bestCandIdx].getHor();
            bestY = cMVCand[bestCandIdx].getVer();
            sadBest = sadBestCand[bestCandIdx];
            rcMv.set(bestX, bestY);
            ruiCost = sadBest;
            goto end;
          }
        }
      }
    }
  }

  bestCandIdx = xIBCSearchMVChromaRefine(pu, roiWidth, roiHeight, cuPelX, cuPelY, sadBestCand, cMVCand);

  bestX = cMVCand[bestCandIdx].getHor();
  bestY = cMVCand[bestCandIdx].getVer();
  sadBest = sadBestCand[bestCandIdx];
  rcMv.set(bestX, bestY);
  ruiCost = sadBest;

end:
  m_numBVs = 0;
  m_numBVs = xMergeCandLists(m_acBVs, m_numBVs, (2 * IBC_NUM_CANDIDATES), m_defaultCachedBvs.m_bvCands, m_defaultCachedBvs.currCnt);

  m_defaultCachedBvs.currCnt = 0;
  m_defaultCachedBvs.currCnt = xMergeCandLists(m_defaultCachedBvs.m_bvCands, m_defaultCachedBvs.currCnt, IBC_NUM_CANDIDATES, cMVCand, CHROMA_REFINEMENT_CANDIDATES);
  m_defaultCachedBvs.currCnt = xMergeCandLists(m_defaultCachedBvs.m_bvCands, m_defaultCachedBvs.currCnt, IBC_NUM_CANDIDATES, m_acBVs, m_numBVs);

  for (unsigned int cand = 0; cand < CHROMA_REFINEMENT_CANDIDATES; cand++)
  {
    if (cMVCand[cand].getHor() == 0 && cMVCand[cand].getVer() == 0)
    {
      continue;
    }
    m_ctuRecord[pu.lumaPos()][pu.lumaSize()].bvRecord[cMVCand[cand]] = sadBestCand[cand];
  }

  return;
}



// based on xMotionEstimation
void InterSearch::xIBCEstimation(PredictionUnit& pu, PelUnitBuf& origBuf,
  Mv     *pcMvPred,
  Mv     &rcMv,
  Distortion &ruiCost, const int localSearchRangeX, const int localSearchRangeY
)
{
  const int iPicWidth = pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
  const int iPicHeight = pu.cs->slice->getPPS()->getPicHeightInLumaSamples();
  const unsigned int  lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
  const int           cuPelX = pu.Y().x;
  const int           cuPelY = pu.Y().y;
  int                 iRoiWidth = pu.lwidth();
  int                 iRoiHeight = pu.lheight();

  PelUnitBuf* pBuf = &origBuf;

  //  Search key pattern initialization
  CPelBuf  tmpPattern = pBuf->Y();
  CPelBuf* pcPatternKey = &tmpPattern;
  PelBuf tmpOrgLuma;

  if ((pu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag()))
  {
    const CompArea &area = pu.blocks[COMPONENT_Y];
    CompArea    tmpArea(COMPONENT_Y, area.chromaFormat, Position(0, 0), area.size());
    tmpOrgLuma = m_tmpStorageLCU.getBuf(tmpArea);
    tmpOrgLuma.rspSignal( tmpPattern, m_pcReshape->getFwdLUT() );
    pcPatternKey = (CPelBuf*)&tmpOrgLuma;
  }

  m_lumaClpRng = pu.cs->slice->clpRng(COMPONENT_Y);
  Picture* refPic = pu.cu->slice->getPic();
  const CPelBuf refBuf = refPic->getRecoBuf(pu.blocks[COMPONENT_Y]);

  IntTZSearchStruct cStruct;
  cStruct.pcPatternKey = pcPatternKey;
  cStruct.iRefStride = refBuf.stride;
  cStruct.piRefY = refBuf.buf;
  CHECK(pu.cu->imv == IMV_HPEL, "IF_IBC");
  cStruct.imvShift = pu.cu->imv << 1;
  cStruct.subShiftMode = 0; // used by intra pattern search function

  // disable weighted prediction
  setWpScalingDistParam(-1, REF_PIC_LIST_X, pu.cs->slice);

  m_pcRdCost->getMotionCost(0);
  m_pcRdCost->setPredictors(pcMvPred);
  m_pcRdCost->setCostScale(0);

  m_cDistParam.useMR = false;
  m_pcRdCost->setDistParam(m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, cStruct.subShiftMode);
  bool buffered = false;
  if (m_pcEncCfg->getIBCFastMethod() & IBC_FAST_METHOD_BUFFERBV)
  {
    ruiCost = MAX_UINT;
    std::unordered_map<Mv, Distortion>& history = m_ctuRecord[pu.lumaPos()][pu.lumaSize()].bvRecord;
    for (std::unordered_map<Mv, Distortion>::iterator p = history.begin(); p != history.end(); p++)
    {
      const Mv& bv = p->first;

      int xBv = bv.hor;
      int yBv = bv.ver;
#if JVET_Z0084_IBC_TM
      if (PU::searchBv(pu, cuPelX, cuPelY, iRoiWidth, iRoiHeight, iPicWidth, iPicHeight, xBv, yBv, lcuWidth))
#else
      if (searchBv(pu, cuPelX, cuPelY, iRoiWidth, iRoiHeight, iPicWidth, iPicHeight, xBv, yBv, lcuWidth))
#endif
      {
        buffered = true;
        Distortion sad = m_pcRdCost->getBvCostMultiplePreds(xBv, yBv, pu.cs->sps->getAMVREnabledFlag());
        m_cDistParam.cur.buf = cStruct.piRefY + cStruct.iRefStride * yBv + xBv;
        sad += m_cDistParam.distFunc(m_cDistParam);
        if (sad < ruiCost)
        {
          rcMv = bv;
          ruiCost = sad;
        }
        else if (sad == ruiCost)
        {
          // stabilise the search through the unordered list
          if (bv.hor < rcMv.getHor()
            || (bv.hor == rcMv.getHor() && bv.ver < rcMv.getVer()))
          {
            // update the vector.
            rcMv = bv;
          }
        }
      }
    }

    if (buffered)
    {
      Mv cMvPredEncOnly[IBC_NUM_CANDIDATES];
      int nbPreds = 0;
      PU::getIbcMVPsEncOnly(pu, cMvPredEncOnly, nbPreds);

      for (unsigned int cand = 0; cand < nbPreds; cand++)
      {
        int xPred = cMvPredEncOnly[cand].getHor();
        int yPred = cMvPredEncOnly[cand].getVer();

#if JVET_Z0084_IBC_TM
        if (PU::searchBv(pu, cuPelX, cuPelY, iRoiWidth, iRoiHeight, iPicWidth, iPicHeight, xPred, yPred, lcuWidth))
#else
        if (searchBv(pu, cuPelX, cuPelY, iRoiWidth, iRoiHeight, iPicWidth, iPicHeight, xPred, yPred, lcuWidth))
#endif
        {
          Distortion sad = m_pcRdCost->getBvCostMultiplePreds(xPred, yPred, pu.cs->sps->getAMVREnabledFlag());
          m_cDistParam.cur.buf = cStruct.piRefY + cStruct.iRefStride * yPred + xPred;
          sad += m_cDistParam.distFunc(m_cDistParam);
          if (sad < ruiCost)
          {
            rcMv.set(xPred, yPred);
            ruiCost = sad;
          }
          else if (sad == ruiCost)
          {
            // stabilise the search through the unordered list
            if (xPred < rcMv.getHor()
              || (xPred == rcMv.getHor() && yPred < rcMv.getVer()))
            {
              // update the vector.
              rcMv.set(xPred, yPred);
            }
          }

          m_ctuRecord[pu.lumaPos()][pu.lumaSize()].bvRecord[Mv(xPred, yPred)] = sad;
        }
      }
    }
  }

  if (!buffered)
  {
    Mv        cMvSrchRngLT;
    Mv        cMvSrchRngRB;

    // assume that intra BV is integer-pel precision
    xSetIntraSearchRange(pu, pu.lwidth(), pu.lheight(), localSearchRangeX, localSearchRangeY, cMvSrchRngLT, cMvSrchRngRB);

    //  Do integer search
    xIntraPatternSearch(pu, cStruct, rcMv, ruiCost, &cMvSrchRngLT, &cMvSrchRngRB, pcMvPred);
  }
}

// based on xSetSearchRange
void InterSearch::xSetIntraSearchRange(PredictionUnit& pu, int iRoiWidth, int iRoiHeight, const int localSearchRangeX, const int localSearchRangeY, Mv& rcMvSrchRngLT, Mv& rcMvSrchRngRB)
{
  const SPS &sps = *pu.cs->sps;

  int srLeft, srRight, srTop, srBottom;

  const int cuPelX = pu.Y().x;
  const int cuPelY = pu.Y().y;

  const int lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
#if JVET_Z0153_IBC_EXT_REF
  const int picWidth = pu.cs->slice->getPPS()->getPicWidthInLumaSamples();

  srLeft = -cuPelX;
  srTop = - 2 * lcuWidth - (cuPelY % lcuWidth);
  srRight = picWidth - cuPelX - iRoiWidth;
  srBottom = lcuWidth - (cuPelY % lcuWidth) - iRoiHeight;  
#else
  const int ctuSizeLog2 = floorLog2(lcuWidth);
  int numLeftCTUs = (1 << ((7 - ctuSizeLog2) << 1)) - ((ctuSizeLog2 < 7) ? 1 : 0);

  srLeft = -(numLeftCTUs * lcuWidth + (cuPelX % lcuWidth));
  srTop = -(cuPelY % lcuWidth);

  srRight = lcuWidth - (cuPelX % lcuWidth) - iRoiWidth;
  srBottom = lcuWidth - (cuPelY % lcuWidth) - iRoiHeight;
#endif

  rcMvSrchRngLT.setHor(srLeft);
  rcMvSrchRngLT.setVer(srTop);
  rcMvSrchRngRB.setHor(srRight);
  rcMvSrchRngRB.setVer(srBottom);

  rcMvSrchRngLT <<= 2;
  rcMvSrchRngRB <<= 2;
  bool temp = m_clipMvInSubPic;
  m_clipMvInSubPic = true;
  xClipMv(rcMvSrchRngLT, pu.cu->lumaPos(),
         pu.cu->lumaSize(),
         sps
      , *pu.cs->pps
  );
  xClipMv(rcMvSrchRngRB, pu.cu->lumaPos(),
         pu.cu->lumaSize(),
         sps
      , *pu.cs->pps
  );
  m_clipMvInSubPic = temp;
  rcMvSrchRngLT >>= 2;
  rcMvSrchRngRB >>= 2;
}

bool InterSearch::predIBCSearch(CodingUnit& cu, Partitioner& partitioner, const int localSearchRangeX, const int localSearchRangeY, IbcHashMap& ibcHashMap)
{
  Mv           cMvSrchRngLT;
  Mv           cMvSrchRngRB;

  Mv           cMv;
  Mv           cMvPred;

#if JVET_Z0131_IBC_BVD_BINARIZATION
  xEstBvdBitCosts(m_pcRdCost->getBvdBitCosts());
#endif

  for (auto &pu : CU::traversePUs(cu))
  {
    m_maxCompIDToPred = MAX_NUM_COMPONENT;

    CHECK(pu.cu != &cu, "PU is contained in another CU");
    //////////////////////////////////////////////////////////
    /// ibc search
    pu.cu->imv = 2;
    AMVPInfo amvpInfo4Pel;
#if JVET_Z0084_IBC_TM && TM_AMVP
    PU::fillIBCMvpCand(pu, amvpInfo4Pel, this);
#else
    PU::fillIBCMvpCand(pu, amvpInfo4Pel);
#endif

    pu.cu->imv = 0;// (Int)cu.cs->sps->getUseIMV(); // set as IMV=0 initially
    Mv    cMv, cMvPred[2];
    AMVPInfo amvpInfo;
#if JVET_Z0084_IBC_TM && TM_AMVP
    PU::fillIBCMvpCand(pu, amvpInfo, this);
#else
    PU::fillIBCMvpCand(pu, amvpInfo);
#endif

    // store in full pel accuracy, shift before use in search
    cMvPred[0] = amvpInfo.mvCand[0];
    cMvPred[0].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
    cMvPred[1] = amvpInfo.mvCand[1];
    cMvPred[1].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);

    int iBvpNum = 2;
    int bvpIdxBest = 0;
    cMv.setZero();
    Distortion cost = 0;
    if (pu.cs->sps->getMaxNumIBCMergeCand() == 1)
    {
      iBvpNum = 1;
      cMvPred[1] = cMvPred[0];
    }

    if (m_pcEncCfg->getIBCHashSearch())
    {
      xxIBCHashSearch(pu, cMvPred, iBvpNum, cMv, bvpIdxBest, ibcHashMap);
    }

    if (cMv.getHor() == 0 && cMv.getVer() == 0)
    {
      // if hash search does not work or is not enabled
      PelUnitBuf origBuf = pu.cs->getOrgBuf(pu);
      xIBCEstimation(pu, origBuf, cMvPred, cMv, cost, localSearchRangeX, localSearchRangeY);
    }

    if (cMv.getHor() == 0 && cMv.getVer() == 0)
    {
      return false;
    }
    /// ibc search
    /////////////////////////////////////////////////////////
#if JVET_Z0131_IBC_BVD_BINARIZATION
    m_pcRdCost->setPredictors(cMvPred);
    m_pcRdCost->setCostScale(0);
#if JVET_Z0084_IBC_TM
    m_pcRdCost->getBvCostMultiplePreds(cMv.getHor(), cMv.getVer(), pu.cs->sps->getAMVREnabledFlag(), &pu.cu->imv, &bvpIdxBest, true, &amvpInfo4Pel);
#else
    m_pcRdCost->getBvCostMultiplePreds(cMv.getHor(), cMv.getVer(), pu.cs->sps->getAMVREnabledFlag(), &pu.cu->imv, &bvpIdxBest);
#endif
#else
    unsigned int bitsBVPBest, bitsBVPTemp;
    bitsBVPBest = MAX_INT;
    m_pcRdCost->setCostScale(0);

    for (int bvpIdxTemp = 0; bvpIdxTemp<iBvpNum; bvpIdxTemp++)
    {
      m_pcRdCost->setPredictor(cMvPred[bvpIdxTemp]);

      bitsBVPTemp = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer(), 0);

      if (bitsBVPTemp < bitsBVPBest)
      {
        bitsBVPBest = bitsBVPTemp;
        bvpIdxBest = bvpIdxTemp;

        if (cu.cs->sps->getAMVREnabledFlag() && cMv != cMvPred[bvpIdxTemp])
          pu.cu->imv = 1; // set as full-pel
        else
          pu.cu->imv = 0; // set as fractional-pel

      }

      unsigned int bitsBVPQP = MAX_UINT;


      Mv mvPredQuadPel;
      if ((cMv.getHor() % 4 == 0) && (cMv.getVer() % 4 == 0) && (pu.cs->sps->getAMVREnabledFlag()))
      {
        mvPredQuadPel = amvpInfo4Pel.mvCand[bvpIdxTemp];// cMvPred[bvpIdxTemp];

        mvPredQuadPel.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_4PEL);

        m_pcRdCost->setPredictor(mvPredQuadPel);

        bitsBVPQP = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor() >> 2, cMv.getVer() >> 2, 0);

      }
      mvPredQuadPel.changePrecision(MV_PRECISION_4PEL, MV_PRECISION_INT);
      if (bitsBVPQP < bitsBVPBest && cMv != mvPredQuadPel)
      {
        bitsBVPBest = bitsBVPQP;
        bvpIdxBest = bvpIdxTemp;

        if (cu.cs->sps->getAMVREnabledFlag())
          pu.cu->imv = 2; // set as quad-pel
      }

    }
#endif

    pu.bv = cMv; // bv is always at integer accuracy
    cMv.changePrecision(MV_PRECISION_INT, MV_PRECISION_INTERNAL);
    pu.mv[REF_PIC_LIST_0] = cMv; // store in fractional pel accuracy

    pu.mvpIdx[REF_PIC_LIST_0] = bvpIdxBest;

    if(pu.cu->imv == 2 && cMv != amvpInfo4Pel.mvCand[bvpIdxBest])
      pu.mvd[REF_PIC_LIST_0] = cMv - amvpInfo4Pel.mvCand[bvpIdxBest];
    else
      pu.mvd[REF_PIC_LIST_0] = cMv - amvpInfo.mvCand[bvpIdxBest];

    if (pu.mvd[REF_PIC_LIST_0] == Mv(0, 0))
      pu.cu->imv = 0;
    if (pu.cu->imv == 2)
      assert((cMv.getHor() % 16 == 0) && (cMv.getVer() % 16 == 0));
    if (cu.cs->sps->getAMVREnabledFlag())
      assert(pu.cu->imv>0 || pu.mvd[REF_PIC_LIST_0] == Mv());

    pu.refIdx[REF_PIC_LIST_0] = MAX_NUM_REF;

  }

  return true;
}

void InterSearch::xxIBCHashSearch(PredictionUnit& pu, Mv* mvPred, int numMvPred, Mv &mv, int& idxMvPred, IbcHashMap& ibcHashMap)
{
  mv.setZero();
  m_pcRdCost->setCostScale(0);

  std::vector<Position> candPos;
  if (ibcHashMap.ibcHashMatch(pu.Y(), candPos, *pu.cs, m_pcEncCfg->getIBCHashSearchMaxCand(), m_pcEncCfg->getIBCHashSearchRange4SmallBlk()))
  {
#if JVET_Z0131_IBC_BVD_BINARIZATION
    Distortion minCost = MAX_UINT64;
    m_pcRdCost->setPredictors(mvPred);
#else
    unsigned int minCost = MAX_UINT;
#endif

    const unsigned int  lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
    const int   cuPelX = pu.Y().x;
    const int   cuPelY = pu.Y().y;
    const int   picWidth = pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
    const int   picHeight = pu.cs->slice->getPPS()->getPicHeightInLumaSamples();
    int         roiWidth = pu.lwidth();
    int         roiHeight = pu.lheight();

    for (std::vector<Position>::iterator pos = candPos.begin(); pos != candPos.end(); pos++)
    {
      Position bottomRight = pos->offset(pu.Y().width - 1, pu.Y().height - 1);
      if (pu.cs->isDecomp(*pos, CHANNEL_TYPE_LUMA) && pu.cs->isDecomp(bottomRight, CHANNEL_TYPE_LUMA))
      {
        Position tmp = *pos - pu.Y().pos();
        Mv candMv;
        candMv.set(tmp.x, tmp.y);

#if JVET_Z0084_IBC_TM
        if (!PU::searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, candMv.getHor(), candMv.getVer(), lcuWidth))
#else
        if (!searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, candMv.getHor(), candMv.getVer(), lcuWidth))
#endif
        {
          continue;
        }

#if JVET_Z0131_IBC_BVD_BINARIZATION
        Distortion cost = m_pcRdCost->getBvCostMultiplePreds(candMv.getHor(), candMv.getVer(), pu.cs->sps->getAMVREnabledFlag());
        if (cost < minCost)
        {
          mv = candMv;
          minCost = cost;
        }
#else
        for (int n = 0; n < numMvPred; n++)
        {
          m_pcRdCost->setPredictor(mvPred[n]);

          unsigned int cost = m_pcRdCost->getBitsOfVectorWithPredictor(candMv.getHor(), candMv.getVer(), 0);

          if (cost < minCost)
          {
            mv = candMv;
            idxMvPred = n;
            minCost = cost;
          }

          int costQuadPel = MAX_UINT;
          if ((candMv.getHor() % 4 == 0) && (candMv.getVer() % 4 == 0) && (pu.cs->sps->getAMVREnabledFlag()))
          {
            Mv mvPredQuadPel;
            int imvShift = 2;
            int offset = 1 << (imvShift - 1);

            int x = (mvPred[n].hor + offset - (mvPred[n].hor >= 0)) >> 2;
            int y = (mvPred[n].ver + offset - (mvPred[n].ver >= 0)) >> 2;
            mvPredQuadPel.set(x, y);

            m_pcRdCost->setPredictor(mvPredQuadPel);

            costQuadPel = m_pcRdCost->getBitsOfVectorWithPredictor(candMv.getHor() >> 2, candMv.getVer() >> 2, 0);

          }
          if (costQuadPel < minCost)
          {
            mv = candMv;
            idxMvPred = n;
            minCost = costQuadPel;
          }

        }
#endif
      }
    }
  }

}


void InterSearch::addToSortList(std::list<BlockHash>& listBlockHash, std::list<int>& listCost, int cost, const BlockHash& blockHash)
{
  std::list<BlockHash>::iterator itBlockHash = listBlockHash.begin();
  std::list<int>::iterator itCost = listCost.begin();

  while (itCost != listCost.end())
  {
    if (cost < (*itCost))
    {
      listCost.insert(itCost, cost);
      listBlockHash.insert(itBlockHash, blockHash);
      return;
    }

    ++itCost;
    ++itBlockHash;
  }

  listCost.push_back(cost);
  listBlockHash.push_back(blockHash);
}

void InterSearch::selectMatchesInter(const MapIterator& itBegin, int count, std::list<BlockHash>& listBlockHash, const BlockHash& currBlockHash)
{
  const int maxReturnNumber = 5;

  listBlockHash.clear();
  std::list<int> listCost;
  listCost.clear();

  MapIterator it = itBegin;
  for (int i = 0; i < count; i++, it++)
  {
    if ((*it).hashValue2 != currBlockHash.hashValue2)
    {
      continue;
    }

    int currCost = RdCost::xGetExpGolombNumberOfBits((*it).x - currBlockHash.x) +
      RdCost::xGetExpGolombNumberOfBits((*it).y - currBlockHash.y);

    if (listBlockHash.size() < maxReturnNumber)
    {
      addToSortList(listBlockHash, listCost, currCost, (*it));
    }
    else if (!listCost.empty() && currCost < listCost.back())
    {
      listCost.pop_back();
      listBlockHash.pop_back();
      addToSortList(listBlockHash, listCost, currCost, (*it));
    }
  }
}
void InterSearch::selectRectangleMatchesInter(const MapIterator& itBegin, int count, std::list<BlockHash>& listBlockHash, const BlockHash& currBlockHash, int width, int height, int idxNonSimple, unsigned int* &hashValues, int baseNum, int picWidth, int picHeight, bool isHorizontal, uint16_t* curHashPic)
{
  const int maxReturnNumber = 5;
  int baseSize = min(width, height);
  unsigned int crcMask = 1 << 16;
  crcMask -= 1;

  listBlockHash.clear();
  std::list<int> listCost;
  listCost.clear();

  MapIterator it = itBegin;

  for (int i = 0; i < count; i++, it++)
  {
    if ((*it).hashValue2 != currBlockHash.hashValue2)
    {
      continue;
    }
    int xRef = (*it).x;
    int yRef = (*it).y;
    if (isHorizontal)
    {
      xRef -= idxNonSimple * baseSize;
    }
    else
    {
      yRef -= idxNonSimple * baseSize;
    }
    if (xRef < 0 || yRef < 0 || xRef + width >= picWidth || yRef + height >= picHeight)
    {
      continue;
    }
    //check Other baseSize hash values
    uint16_t* refHashValue = curHashPic + yRef * picWidth + xRef;
    bool isSame = true;

    for (int k = 0; k < baseNum; k++)
    {
      if ((*refHashValue) != (uint16_t)(hashValues[k] & crcMask))
      {
        isSame = false;
        break;
      }
      refHashValue += (isHorizontal ? baseSize : (baseSize*picWidth));
    }
    if (!isSame)
    {
      continue;
    }

    int currCost = RdCost::xGetExpGolombNumberOfBits(xRef - currBlockHash.x) +
      RdCost::xGetExpGolombNumberOfBits(yRef - currBlockHash.y);

    BlockHash refBlockHash;
    refBlockHash.hashValue2 = (*it).hashValue2;
    refBlockHash.x = xRef;
    refBlockHash.y = yRef;

    if (listBlockHash.size() < maxReturnNumber)
    {
      addToSortList(listBlockHash, listCost, currCost, refBlockHash);
    }
    else if (!listCost.empty() && currCost < listCost.back())
    {
      listCost.pop_back();
      listBlockHash.pop_back();
      addToSortList(listBlockHash, listCost, currCost, refBlockHash);
    }
  }
}

bool InterSearch::xRectHashInterEstimation(PredictionUnit& pu, RefPicList& bestRefPicList, int& bestRefIndex, Mv& bestMv, Mv& bestMvd, int& bestMVPIndex, bool& isPerfectMatch)
{
  int width = pu.cu->lumaSize().width;
  int height = pu.cu->lumaSize().height;

  int baseSize = min(width, height);
  bool isHorizontal = true;;
  int baseNum = 0;
  if (height < width)
  {
    isHorizontal = true;
    baseNum = 1 << (floorLog2(width) - floorLog2(height));
  }
  else
  {
    isHorizontal = false;
    baseNum = 1 << (floorLog2(height) - floorLog2(width));
  }

  int xPos = pu.cu->lumaPos().x;
  int yPos = pu.cu->lumaPos().y;
  const int currStride = pu.cs->picture->getOrigBuf().get(COMPONENT_Y).stride;
  const Pel* curPel = pu.cs->picture->getOrigBuf().get(COMPONENT_Y).buf + yPos * currStride + xPos;
  int picWidth = pu.cu->slice->getPPS()->getPicWidthInLumaSamples();
  int picHeight = pu.cu->slice->getPPS()->getPicHeightInLumaSamples();

  int xBase = xPos;
  int yBase = yPos;
  const Pel* basePel = curPel;
  int idxNonSimple = -1;
  unsigned int* hashValue1s = new unsigned int[baseNum];
  unsigned int* hashValue2s = new unsigned int[baseNum];

  for (int k = 0; k < baseNum; k++)
  {
    if (isHorizontal)
    {
      xBase = xPos + k * baseSize;
      basePel = curPel + k * baseSize;
    }
    else
    {
      yBase = yPos + k * baseSize;
      basePel = curPel + k * baseSize * currStride;
    }

    if (idxNonSimple == -1 && !TComHash::isHorizontalPerfectLuma(basePel, currStride, baseSize, baseSize) && !TComHash::isVerticalPerfectLuma(basePel, currStride, baseSize, baseSize))
    {
      idxNonSimple = k;
    }
    TComHash::getBlockHashValue((pu.cs->picture->getOrigBuf()), baseSize, baseSize, xBase, yBase, pu.cu->slice->getSPS()->getBitDepths(), hashValue1s[k], hashValue2s[k]);
  }
  if (idxNonSimple == -1)
  {
    idxNonSimple = 0;
  }

  Distortion bestCost = UINT64_MAX;

  BlockHash currBlockHash;
  currBlockHash.x = xPos;//still use the first base block location
  currBlockHash.y = yPos;

  currBlockHash.hashValue2 = hashValue2s[idxNonSimple];

  m_pcRdCost->setDistParam(m_cDistParam, pu.cs->getOrgBuf(pu).Y(), 0, 0, m_lumaClpRng.bd, COMPONENT_Y, 0, 1, false);

  int imvBest = 0;
  int numPredDir = pu.cu->slice->isInterP() ? 1 : 2;
  for (int refList = 0; refList < numPredDir; refList++)
  {
    RefPicList eRefPicList = (refList == 0) ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
    int refPicNumber = pu.cu->slice->getNumRefIdx(eRefPicList);

    for (int refIdx = 0; refIdx < refPicNumber; refIdx++)
    {
      int bitsOnRefIdx = 1;
      if (refPicNumber > 1)
      {
        bitsOnRefIdx += refIdx + 1;
        if (refIdx == refPicNumber - 1)
        {
          bitsOnRefIdx--;
        }
      }
      m_numHashMVStoreds[eRefPicList][refIdx] = 0;

      const std::pair<int, int>& scaleRatio = pu.cu->slice->getScalingRatio( eRefPicList, refIdx );
      if( scaleRatio != SCALE_1X )
      {
        continue;
      }

      CHECK( pu.cu->slice->getRefPic( eRefPicList, refIdx )->getHashMap() == nullptr, "Hash table is not initialized" );

      if (refList == 0 || pu.cu->slice->getList1IdxToList0Idx(refIdx) < 0)
      {
        int count = static_cast<int>(pu.cu->slice->getRefPic(eRefPicList, refIdx)->getHashMap()->count(hashValue1s[idxNonSimple]));
        if (count == 0)
        {
          continue;
        }

        list<BlockHash> listBlockHash;
        selectRectangleMatchesInter(pu.cu->slice->getRefPic(eRefPicList, refIdx)->getHashMap()->getFirstIterator(hashValue1s[idxNonSimple]), count, listBlockHash, currBlockHash, width, height, idxNonSimple, hashValue2s, baseNum, picWidth, picHeight, isHorizontal, pu.cu->slice->getRefPic(eRefPicList, refIdx)->getHashMap()->getHashPic(baseSize));

        m_numHashMVStoreds[eRefPicList][refIdx] = int(listBlockHash.size());
        if (listBlockHash.empty())
        {
          continue;
        }
        AMVPInfo currAMVPInfoPel;
        AMVPInfo currAMVPInfo4Pel;
        AMVPInfo currAMVPInfoQPel;
        pu.cu->imv = 2;
        PU::fillMvpCand(pu, eRefPicList, refIdx, currAMVPInfo4Pel
#if TM_AMVP
                      , this
#endif
        );
        pu.cu->imv = 1;
        PU::fillMvpCand(pu, eRefPicList, refIdx, currAMVPInfoPel
#if TM_AMVP
                      , this
#endif
        );
        pu.cu->imv = 0;
        PU::fillMvpCand(pu, eRefPicList, refIdx, currAMVPInfoQPel
#if TM_AMVP
                      , this
#endif
        );
#if TM_AMVP
        CHECK(currAMVPInfoPel.numCand != currAMVPInfoQPel.numCand, "The number of full-Pel AMVP candidates and that of Q-Pel should be identical");
        CHECK(currAMVPInfoPel.numCand != currAMVPInfo4Pel.numCand, "The number of full-Pel AMVP candidates and that of 4-Pel should be identical");
        const uint8_t amvpNumCand = currAMVPInfoPel.numCand;
        for (int mvpIdxTemp = 0; mvpIdxTemp < amvpNumCand; mvpIdxTemp++)
#else
        for (int mvpIdxTemp = 0; mvpIdxTemp < 2; mvpIdxTemp++)
#endif
        {
          currAMVPInfoQPel.mvCand[mvpIdxTemp].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
          currAMVPInfoPel.mvCand[mvpIdxTemp].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
          currAMVPInfo4Pel.mvCand[mvpIdxTemp].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
        }

        bool wrap = pu.cu->slice->getRefPic(eRefPicList, refIdx)->isWrapAroundEnabled( pu.cs->pps );
        const Pel* refBufStart = pu.cu->slice->getRefPic(eRefPicList, refIdx)->getRecoBuf(wrap).get(COMPONENT_Y).buf;
        const int refStride = pu.cu->slice->getRefPic(eRefPicList, refIdx)->getRecoBuf(wrap).get(COMPONENT_Y).stride;
        m_cDistParam.cur.stride = refStride;

        m_pcRdCost->selectMotionLambda( );
        m_pcRdCost->setCostScale(0);

        list<BlockHash>::iterator it;
        int countMV = 0;
        for (it = listBlockHash.begin(); it != listBlockHash.end(); ++it)
        {
          int curMVPIdx = 0;
          unsigned int curMVPbits = MAX_UINT;
          Mv cMv((*it).x - currBlockHash.x, (*it).y - currBlockHash.y);
          m_hashMVStoreds[eRefPicList][refIdx][countMV++] = cMv;
          cMv.changePrecision(MV_PRECISION_INT, MV_PRECISION_QUARTER);

#if TM_AMVP
          for (int mvpIdxTemp = 0; mvpIdxTemp < amvpNumCand; mvpIdxTemp++)
#else
          for (int mvpIdxTemp = 0; mvpIdxTemp < 2; mvpIdxTemp++)
#endif
          {
            Mv cMvPredPel = currAMVPInfoQPel.mvCand[mvpIdxTemp];
            m_pcRdCost->setPredictor(cMvPredPel);

            unsigned int tempMVPbits = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer(), 0);

            if (tempMVPbits < curMVPbits)
            {
              curMVPbits = tempMVPbits;
              curMVPIdx = mvpIdxTemp;
              pu.cu->imv = 0;
            }

            if (pu.cu->slice->getSPS()->getAMVREnabledFlag())
            {
              unsigned int bitsMVP1Pel = MAX_UINT;
              Mv mvPred1Pel = currAMVPInfoPel.mvCand[mvpIdxTemp];
              m_pcRdCost->setPredictor(mvPred1Pel);
              bitsMVP1Pel = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer(), 2);
              if (bitsMVP1Pel < curMVPbits)
              {
                curMVPbits = bitsMVP1Pel;
                curMVPIdx = mvpIdxTemp;
                pu.cu->imv = 1;
              }

              if ((cMv.getHor() % 16 == 0) && (cMv.getVer() % 16 == 0))
              {
                unsigned int bitsMVP4Pel = MAX_UINT;
                Mv mvPred4Pel = currAMVPInfo4Pel.mvCand[mvpIdxTemp];
                m_pcRdCost->setPredictor(mvPred4Pel);
                bitsMVP4Pel = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer(), 4);
                if (bitsMVP4Pel < curMVPbits)
                {
                  curMVPbits = bitsMVP4Pel;
                  curMVPIdx = mvpIdxTemp;
                  pu.cu->imv = 2;
                }
              }
            }
          }
          curMVPbits += bitsOnRefIdx;

          m_cDistParam.cur.buf = refBufStart + (*it).y*refStride + (*it).x;
          Distortion currSad = m_cDistParam.distFunc(m_cDistParam);
          Distortion currCost = currSad + m_pcRdCost->getCost(curMVPbits);

          if (!isPerfectMatch)
          {
            if (pu.cu->slice->getRefPic(eRefPicList, refIdx)->slices[0]->getSliceQp() <= pu.cu->slice->getSliceQp())
            {
              isPerfectMatch = true;
            }
          }

          if (currCost < bestCost)
          {
            bestCost = currCost;
            bestRefPicList = eRefPicList;
            bestRefIndex = refIdx;
            bestMv = cMv;
            bestMVPIndex = curMVPIdx;
            imvBest = pu.cu->imv;
            if (pu.cu->imv == 2)
            {
              bestMvd = cMv - currAMVPInfo4Pel.mvCand[curMVPIdx];
            }
            else if (pu.cu->imv == 1)
            {
              bestMvd = cMv - currAMVPInfoPel.mvCand[curMVPIdx];
            }
            else
            {
              bestMvd = cMv - currAMVPInfoQPel.mvCand[curMVPIdx];
            }
          }
        }
      }
    }
  }
  delete[] hashValue1s;
  delete[] hashValue2s;
  pu.cu->imv = imvBest;
  if (bestMvd == Mv(0, 0))
  {
    pu.cu->imv = 0;
    return false;
  }
  return (bestCost < MAX_INT);
}

bool InterSearch::xHashInterEstimation(PredictionUnit& pu, RefPicList& bestRefPicList, int& bestRefIndex, Mv& bestMv, Mv& bestMvd, int& bestMVPIndex, bool& isPerfectMatch)
{
  int width = pu.cu->lumaSize().width;
  int height = pu.cu->lumaSize().height;
  if (width != height)
  {
    return xRectHashInterEstimation(pu, bestRefPicList, bestRefIndex, bestMv, bestMvd, bestMVPIndex, isPerfectMatch);
  }
  int xPos = pu.cu->lumaPos().x;
  int yPos = pu.cu->lumaPos().y;

  uint32_t hashValue1;
  uint32_t hashValue2;
  Distortion bestCost = UINT64_MAX;

  if (!TComHash::getBlockHashValue((pu.cs->picture->getOrigBuf()), width, height, xPos, yPos, pu.cu->slice->getSPS()->getBitDepths(), hashValue1, hashValue2))
  {
    return false;
  }

  BlockHash currBlockHash;
  currBlockHash.x = xPos;
  currBlockHash.y = yPos;
  currBlockHash.hashValue2 = hashValue2;

  m_pcRdCost->setDistParam(m_cDistParam, pu.cs->getOrgBuf(pu).Y(), 0, 0, m_lumaClpRng.bd, COMPONENT_Y, 0, 1, false);

  int imvBest = 0;

  int numPredDir = pu.cu->slice->isInterP() ? 1 : 2;
  for (int refList = 0; refList < numPredDir; refList++)
  {
    RefPicList eRefPicList = (refList == 0) ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
    int refPicNumber = pu.cu->slice->getNumRefIdx(eRefPicList);


    for (int refIdx = 0; refIdx < refPicNumber; refIdx++)
    {
      int bitsOnRefIdx = 1;
      if (refPicNumber > 1)
      {
        bitsOnRefIdx += refIdx + 1;
        if (refIdx == refPicNumber - 1)
        {
          bitsOnRefIdx--;
        }
      }
      m_numHashMVStoreds[eRefPicList][refIdx] = 0;

      const std::pair<int, int>& scaleRatio = pu.cu->slice->getScalingRatio( eRefPicList, refIdx );
      if( scaleRatio != SCALE_1X )
      {
        continue;
      }

      CHECK( pu.cu->slice->getRefPic( eRefPicList, refIdx )->getHashMap() == nullptr, "Hash table is not initialized" );

      if (refList == 0 || pu.cu->slice->getList1IdxToList0Idx(refIdx) < 0)
      {
        int count = static_cast<int>(pu.cu->slice->getRefPic(eRefPicList, refIdx)->getHashMap()->count(hashValue1));
        if (count == 0)
        {
          continue;
        }

        list<BlockHash> listBlockHash;
        selectMatchesInter(pu.cu->slice->getRefPic(eRefPicList, refIdx)->getHashMap()->getFirstIterator(hashValue1), count, listBlockHash, currBlockHash);
        m_numHashMVStoreds[eRefPicList][refIdx] = (int)listBlockHash.size();
        if (listBlockHash.empty())
        {
          continue;
        }
        AMVPInfo currAMVPInfoPel;
        AMVPInfo currAMVPInfo4Pel;
        pu.cu->imv = 2;
        PU::fillMvpCand(pu, eRefPicList, refIdx, currAMVPInfo4Pel
#if TM_AMVP
                      , this
#endif
        );
        pu.cu->imv = 1;
        PU::fillMvpCand(pu, eRefPicList, refIdx, currAMVPInfoPel
#if TM_AMVP
                      , this
#endif
        );
        AMVPInfo currAMVPInfoQPel;
        pu.cu->imv = 0;
        PU::fillMvpCand(pu, eRefPicList, refIdx, currAMVPInfoQPel
#if TM_AMVP
                      , this
#endif
        );
#if TM_AMVP
        CHECK(currAMVPInfoPel.numCand != currAMVPInfoQPel.numCand, "The number of full-Pel AMVP candidates and that of Q-Pel should be identical");
        CHECK(currAMVPInfoPel.numCand != currAMVPInfo4Pel.numCand, "The number of full-Pel AMVP candidates and that of 4-Pel should be identical");
        const uint8_t amvpNumCand = currAMVPInfoPel.numCand;
        CHECK(currAMVPInfoPel.numCand == 0, "Wrong")
        for (int mvpIdxTemp = 0; mvpIdxTemp < amvpNumCand; mvpIdxTemp++)
#else
        CHECK(currAMVPInfoPel.numCand <= 1, "Wrong")
        for (int mvpIdxTemp = 0; mvpIdxTemp < 2; mvpIdxTemp++)
#endif
        {
          currAMVPInfoQPel.mvCand[mvpIdxTemp].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
          currAMVPInfoPel.mvCand[mvpIdxTemp].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
          currAMVPInfo4Pel.mvCand[mvpIdxTemp].changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
        }

        bool wrap = pu.cu->slice->getRefPic(eRefPicList, refIdx)->isWrapAroundEnabled( pu.cs->pps );
        const Pel* refBufStart = pu.cu->slice->getRefPic(eRefPicList, refIdx)->getRecoBuf(wrap).get(COMPONENT_Y).buf;
        const int refStride = pu.cu->slice->getRefPic(eRefPicList, refIdx)->getRecoBuf(wrap).get(COMPONENT_Y).stride;

        m_cDistParam.cur.stride = refStride;

        m_pcRdCost->selectMotionLambda( );
        m_pcRdCost->setCostScale(0);

        list<BlockHash>::iterator it;
        int countMV = 0;
        for (it = listBlockHash.begin(); it != listBlockHash.end(); ++it)
        {
          int curMVPIdx = 0;
          unsigned int curMVPbits = MAX_UINT;
          Mv cMv((*it).x - currBlockHash.x, (*it).y - currBlockHash.y);
          m_hashMVStoreds[eRefPicList][refIdx][countMV++] = cMv;
          cMv.changePrecision(MV_PRECISION_INT, MV_PRECISION_QUARTER);

#if TM_AMVP
          for (int mvpIdxTemp = 0; mvpIdxTemp < amvpNumCand; mvpIdxTemp++)
#else
          for (int mvpIdxTemp = 0; mvpIdxTemp < 2; mvpIdxTemp++)
#endif
          {
            Mv cMvPredPel = currAMVPInfoQPel.mvCand[mvpIdxTemp];
            m_pcRdCost->setPredictor(cMvPredPel);

            unsigned int tempMVPbits = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer(), 0);

            if (tempMVPbits < curMVPbits)
            {
              curMVPbits = tempMVPbits;
              curMVPIdx = mvpIdxTemp;
              pu.cu->imv = 0;
            }

            if (pu.cu->slice->getSPS()->getAMVREnabledFlag())
            {
              unsigned int bitsMVP1Pel = MAX_UINT;
              Mv mvPred1Pel = currAMVPInfoPel.mvCand[mvpIdxTemp];
              m_pcRdCost->setPredictor(mvPred1Pel);
              bitsMVP1Pel = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer(), 2);
              if (bitsMVP1Pel < curMVPbits)
              {
                curMVPbits = bitsMVP1Pel;
                curMVPIdx = mvpIdxTemp;
                pu.cu->imv = 1;
              }

              if ((cMv.getHor() % 16 == 0) && (cMv.getVer() % 16 == 0))
              {
                unsigned int bitsMVP4Pel = MAX_UINT;
                Mv mvPred4Pel = currAMVPInfo4Pel.mvCand[mvpIdxTemp];
                m_pcRdCost->setPredictor(mvPred4Pel);
                bitsMVP4Pel = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer(), 4);
                if (bitsMVP4Pel < curMVPbits)
                {
                  curMVPbits = bitsMVP4Pel;
                  curMVPIdx = mvpIdxTemp;
                  pu.cu->imv = 2;
                }
              }
            }
          }

          curMVPbits += bitsOnRefIdx;

          m_cDistParam.cur.buf = refBufStart + (*it).y*refStride + (*it).x;
          Distortion currSad = m_cDistParam.distFunc(m_cDistParam);
          Distortion currCost = currSad + m_pcRdCost->getCost(curMVPbits);

          if (!isPerfectMatch)
          {
            if (pu.cu->slice->getRefPic(eRefPicList, refIdx)->slices[0]->getSliceQp() <= pu.cu->slice->getSliceQp())
            {
              isPerfectMatch = true;
            }
          }

          if (currCost < bestCost)
          {
            bestCost = currCost;
            bestRefPicList = eRefPicList;
            bestRefIndex = refIdx;
            bestMv = cMv;
            bestMVPIndex = curMVPIdx;
            imvBest = pu.cu->imv;
            if (pu.cu->imv == 2)
            {
              bestMvd = cMv - currAMVPInfo4Pel.mvCand[curMVPIdx];
            }
            else if (pu.cu->imv == 1)
            {
              bestMvd = cMv - currAMVPInfoPel.mvCand[curMVPIdx];
            }
            else
            {
              bestMvd = cMv - currAMVPInfoQPel.mvCand[curMVPIdx];
            }
          }
        }
      }
    }
  }
  pu.cu->imv = imvBest;
  if (bestMvd == Mv(0, 0))
  {
    pu.cu->imv = 0;
    return false;
  }
  return (bestCost < MAX_INT);
}

bool InterSearch::predInterHashSearch(CodingUnit& cu, Partitioner& partitioner, bool& isPerfectMatch)
{
  Mv       bestMv, bestMvd;
  RefPicList   bestRefPicList;
  int          bestRefIndex;
  int          bestMVPIndex;

  auto &pu = *cu.firstPU;

  Mv cMvZero;
  pu.mv[REF_PIC_LIST_0] = Mv();
  pu.mv[REF_PIC_LIST_1] = Mv();
  pu.mvd[REF_PIC_LIST_0] = cMvZero;
  pu.mvd[REF_PIC_LIST_1] = cMvZero;
  pu.refIdx[REF_PIC_LIST_0] = NOT_VALID;
  pu.refIdx[REF_PIC_LIST_1] = NOT_VALID;
  pu.mvpIdx[REF_PIC_LIST_0] = NOT_VALID;
  pu.mvpIdx[REF_PIC_LIST_1] = NOT_VALID;
  pu.mvpNum[REF_PIC_LIST_0] = NOT_VALID;
  pu.mvpNum[REF_PIC_LIST_1] = NOT_VALID;

  if (xHashInterEstimation(pu, bestRefPicList, bestRefIndex, bestMv, bestMvd, bestMVPIndex, isPerfectMatch))
  {
    pu.interDir = static_cast<int>(bestRefPicList) + 1;
    pu.mv[bestRefPicList] = bestMv;
    pu.mv[bestRefPicList].changePrecision(MV_PRECISION_QUARTER, MV_PRECISION_INTERNAL);

    pu.mvd[bestRefPicList] = bestMvd;
    pu.mvd[bestRefPicList].changePrecision(MV_PRECISION_QUARTER, MV_PRECISION_INTERNAL);
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if JVET_Z0054_BLK_REF_PIC_REORDER
    if (!PU::useRefPairList(pu) && !PU::useRefCombList(pu))
#endif
    if (pu.isMvsdApplicable())
    {
      std::vector<Mv> cMvdDerivedVec;
      Mv cMvPred = pu.mv[bestRefPicList] - pu.mvd[bestRefPicList];
      Mv cMvdKnownAtDecoder = Mv(pu.mvd[bestRefPicList].getAbsHor(), pu.mvd[bestRefPicList].getAbsVer());
      deriveMvdSign(cMvPred, cMvdKnownAtDecoder, pu, bestRefPicList, bestRefIndex, cMvdDerivedVec);
      int idx = deriveMVSDIdxFromMVDTrans(pu.mvd[bestRefPicList], cMvdDerivedVec);
      CHECK(idx == -1, "");
      pu.mvsdIdx[bestRefPicList] = idx;
    }
#endif
    pu.refIdx[bestRefPicList] = bestRefIndex;
    pu.mvpIdx[bestRefPicList] = bestMVPIndex;

#if TM_AMVP
#if JVET_Y0128_NON_CTC
    pu.mvpNum[bestRefPicList] = PU::checkTmEnableCondition(pu.cs->sps, pu.cs->pps, pu.cu->slice->getRefPic(bestRefPicList, bestRefIndex)) ? 1 : 2;
#else
    pu.mvpNum[bestRefPicList] = 1;
#endif
#else
    pu.mvpNum[bestRefPicList] = 2;
#endif
#if JVET_Z0054_BLK_REF_PIC_REORDER
    if (PU::useRefCombList(pu))
    {
      setUniRefIdxLC(pu);
    }
    else if (PU::useRefPairList(pu))
    {
      setBiRefPairIdx(pu);
    }
#endif

    PU::spanMotionInfo(pu);
    PelUnitBuf predBuf = pu.cs->getPredBuf(pu);
    motionCompensation(pu, predBuf, REF_PIC_LIST_X);
    return true;
  }
  else
  {
    return false;
  }

  return true;
}


//! search of the best candidate for inter prediction
#if JVET_X0083_BM_AMVP_MERGE_MODE
void InterSearch::predInterSearch(CodingUnit& cu, Partitioner& partitioner, bool& bdmvrAmMergeNotValid,
    MvField* mvFieldAmListCommon, Mv* mvBufEncAmBDMVR_L0, Mv* mvBufEncAmBDMVR_L1)
#else
void InterSearch::predInterSearch(CodingUnit& cu, Partitioner& partitioner)
#endif
{
  CodingStructure& cs = *cu.cs;

  AMVPInfo     amvp[2];
  Mv           cMvSrchRngLT;
  Mv           cMvSrchRngRB;

  Mv           cMvZero;

  Mv           cMv[2];
  Mv           cMvBi[2];
  Mv           cMvTemp[2][33];
  Mv           cMvHevcTemp[2][33];
  int          iNumPredDir = cs.slice->isInterP() ? 1 : 2;

  Mv           cMvPred[2][33];

  Mv           cMvPredBi[2][33];
  int          aaiMvpIdxBi[2][33];

  int          aaiMvpIdx[2][33];
  int          aaiMvpNum[2][33];

  AMVPInfo     aacAMVPInfo[2][33];

  int          iRefIdx[2]={0,0}; //If un-initialized, may cause SEGV in bi-directional prediction iterative stage.
  int          iRefIdxBi[2] = { -1, -1 };

  uint32_t         uiMbBits[3] = {1, 1, 0};

  uint32_t         uiLastMode = 0;
  uint32_t         uiLastModeTemp = 0;
  int          iRefStart, iRefEnd;

  int          symMode = 0;

  int          bestBiPRefIdxL1 = 0;
  int          bestBiPMvpL1    = 0;
  Distortion   biPDistTemp     = std::numeric_limits<Distortion>::max();

  uint8_t      bcwIdx          = (cu.cs->slice->isInterB() ? cu.BcwIdx : BCW_DEFAULT);
  bool         enforceBcwPred = false;
  MergeCtx     mergeCtx;

  // Loop over Prediction Units
  CHECK(!cu.firstPU, "CU does not contain any PUs");
  uint32_t         puIdx = 0;
  auto &pu = *cu.firstPU;
  WPScalingParam *wp0;
  WPScalingParam *wp1;
  int tryBipred = 0;
  bool checkAffine    = (pu.cu->imv == 0 || pu.cu->slice->getSPS()->getAffineAmvrEnabledFlag()) && pu.cu->imv != IMV_HPEL;
  bool checkNonAffine = pu.cu->imv == 0 || pu.cu->imv == IMV_HPEL || (pu.cu->slice->getSPS()->getAMVREnabledFlag() &&
                                            pu.cu->imv <= (pu.cu->slice->getSPS()->getAMVREnabledFlag() ? IMV_4PEL : 0));
  CodingUnit *bestCU  = pu.cu->cs->bestCS != nullptr ? pu.cu->cs->bestCS->getCU( CHANNEL_TYPE_LUMA ) : nullptr;
  bool trySmvd        = ( bestCU != nullptr && pu.cu->imv == 2 && checkAffine ) ? ( !bestCU->firstPU->mergeFlag && !bestCU->affine ) : true;
  if ( pu.cu->imv && bestCU != nullptr && checkAffine )
  {
    checkAffine = !( bestCU->firstPU->mergeFlag || !bestCU->affine );
  }

#if JVET_X0083_BM_AMVP_MERGE_MODE
  const bool amvpMergeModeFlag = pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1];
  RefPicList refListAmvp       = REF_PIC_LIST_X;
  RefPicList refListMerge      = REF_PIC_LIST_X;
  int candidateRefIdxCount     = 0;
  if (amvpMergeModeFlag)
  {
#if JVET_Y0128_NON_CTC
    if (pu.cu->slice->getUseAmvpMergeMode() == false)
    {
      m_skipPROF = false;
      m_encOnly = false;
      bdmvrAmMergeNotValid = true;
      return;
    }
#endif
    trySmvd = false;
    checkAffine = false;
    refListMerge = pu.amvpMergeModeFlag[0] ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
    refListAmvp = RefPicList(1 - refListMerge);
    getAmvpMergeModeMergeList(pu, mvFieldAmListCommon);
    for (int iRefIdxTemp = 0; iRefIdxTemp < cs.slice->getNumRefIdx(refListAmvp); iRefIdxTemp++)
    {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      if (mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS_MEM].refIdx < 0
          && mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS_MEM + 1].refIdx < 0
          && mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS_MEM + 2].refIdx < 0)
#else
      if (mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS].refIdx < 0 && mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS + 1].refIdx < 0)
#endif
      {
        continue;
      }
      candidateRefIdxCount++;
    }
  }
#if JVET_Y0128_NON_CTC || JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
  if ( amvpMergeModeFlag && !candidateRefIdxCount )
  {
    m_skipPROF = false;
    m_encOnly = false;
    bdmvrAmMergeNotValid = true;
    return;
  }
#endif
#endif
  if ( pu.cu->imv == 2 && checkNonAffine && pu.cu->slice->getSPS()->getAffineAmvrEnabledFlag() )
  {
#if AMVR_ENC_OPT
    checkNonAffine = m_affineMotion.hevcCost[1] < m_affineMotion.hevcCost[0];
#else
    checkNonAffine = m_affineMotion.hevcCost[1] < m_affineMotion.hevcCost[0] * 1.06f;
#endif
  }

#if MULTI_HYP_PRED
  const bool saveMeResultsForMHP = cs.sps->getUseInterMultiHyp()
    && bcwIdx != BCW_DEFAULT
    && (pu.Y().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE && std::min(pu.Y().width, pu.Y().height) >= MULTI_HYP_PRED_RESTRICT_MIN_WH)
      ;
#endif
  {
    if (pu.cu->cs->bestParent != nullptr && pu.cu->cs->bestParent->getCU(CHANNEL_TYPE_LUMA) != nullptr && pu.cu->cs->bestParent->getCU(CHANNEL_TYPE_LUMA)->affine == false)
    {
      m_skipPROF = true;
    }
    m_encOnly = true;
    // motion estimation only evaluates luma component
    m_maxCompIDToPred = MAX_NUM_COMPONENT;
//    m_maxCompIDToPred = COMPONENT_Y;

    CHECK(pu.cu != &cu, "PU is contained in another CU");

    if (cu.cs->sps->getSbTMVPEnabledFlag())
    {
      Size bufSize = g_miScaling.scale(pu.lumaSize());
      mergeCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
    }

    Distortion   uiHevcCost = std::numeric_limits<Distortion>::max();
    Distortion   uiAffineCost = std::numeric_limits<Distortion>::max();
    Distortion   uiCost[2] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max() };
    Distortion   uiCostBi  =   std::numeric_limits<Distortion>::max();
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
    Distortion   uiCostTemp = std::numeric_limits<Distortion>::max();
#else
    Distortion   uiCostTemp;
#endif

    uint32_t         uiBits[3];
    uint32_t         uiBitsTemp;
    Distortion   bestBiPDist = std::numeric_limits<Distortion>::max();

    Distortion   uiCostTempL0[MAX_NUM_REF];
    for (int iNumRef=0; iNumRef < MAX_NUM_REF; iNumRef++)
    {
      uiCostTempL0[iNumRef] = std::numeric_limits<Distortion>::max();
    }
    uint32_t         uiBitsTempL0[MAX_NUM_REF];

    Mv           mvValidList1;
    int          refIdxValidList1 = 0;
    uint32_t         bitsValidList1   = MAX_UINT;
    Distortion   costValidList1   = std::numeric_limits<Distortion>::max();

    PelUnitBuf origBuf = pu.cs->getOrgBuf( pu );

    xGetBlkBits( cs.slice->isInterP(), puIdx, uiLastMode, uiMbBits );

    m_pcRdCost->selectMotionLambda( );

    unsigned imvShift = pu.cu->imv == IMV_HPEL ? 1 : (pu.cu->imv << 1);
    if ( checkNonAffine )
    {
#if JVET_X0083_BM_AMVP_MERGE_MODE
      if (!amvpMergeModeFlag)
      {
#endif
      //  Uni-directional prediction
      for ( int iRefList = 0; iRefList < iNumPredDir; iRefList++ )
      {
        RefPicList  eRefPicList = ( iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0 );
        for (int iRefIdxTemp = 0; iRefIdxTemp < cs.slice->getNumRefIdx(eRefPicList); iRefIdxTemp++)
        {
          uiBitsTemp = uiMbBits[iRefList];
          if ( cs.slice->getNumRefIdx(eRefPicList) > 1 )
          {
            uiBitsTemp += iRefIdxTemp+1;
            if ( iRefIdxTemp == cs.slice->getNumRefIdx(eRefPicList)-1 )
            {
              uiBitsTemp--;
            }
          }
          xEstimateMvPredAMVP( pu, origBuf, eRefPicList, iRefIdxTemp, cMvPred[iRefList][iRefIdxTemp], amvp[eRefPicList], false, &biPDistTemp);

          aaiMvpIdx[iRefList][iRefIdxTemp] = pu.mvpIdx[eRefPicList];
          aaiMvpNum[iRefList][iRefIdxTemp] = pu.mvpNum[eRefPicList];

          if(cs.picHeader->getMvdL1ZeroFlag() && iRefList==1 && biPDistTemp < bestBiPDist)
          {
            bestBiPDist = biPDistTemp;
            bestBiPMvpL1 = aaiMvpIdx[iRefList][iRefIdxTemp];
            bestBiPRefIdxL1 = iRefIdxTemp;
          }

#if TM_AMVP
          uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdx[iRefList][iRefIdxTemp]][aaiMvpNum[iRefList][iRefIdxTemp]];
#else
          uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdx[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];
#endif

          if ( m_pcEncCfg->getFastMEForGenBLowDelayEnabled() && iRefList == 1 )    // list 1
          {
            if ( cs.slice->getList1IdxToList0Idx( iRefIdxTemp ) >= 0 )
            {
              cMvTemp[1][iRefIdxTemp] = cMvTemp[0][cs.slice->getList1IdxToList0Idx( iRefIdxTemp )];
              uiCostTemp = uiCostTempL0[cs.slice->getList1IdxToList0Idx( iRefIdxTemp )];
              /*first subtract the bit-rate part of the cost of the other list*/
              uiCostTemp -= m_pcRdCost->getCost( uiBitsTempL0[cs.slice->getList1IdxToList0Idx( iRefIdxTemp )] );
              /*correct the bit-rate part of the current ref*/
              m_pcRdCost->setPredictor  ( cMvPred[iRefList][iRefIdxTemp] );
              uiBitsTemp += m_pcRdCost->getBitsOfVectorWithPredictor( cMvTemp[1][iRefIdxTemp].getHor(), cMvTemp[1][iRefIdxTemp].getVer(), imvShift + MV_FRACTIONAL_BITS_DIFF );
              /*calculate the correct cost*/
              uiCostTemp += m_pcRdCost->getCost( uiBitsTemp );
            }
            else
            {
              xMotionEstimation( pu, origBuf, eRefPicList, cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, amvp[eRefPicList] );
            }
          }
          else
          {
            xMotionEstimation( pu, origBuf, eRefPicList, cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, amvp[eRefPicList] );
          }
          if( cu.cs->sps->getUseBcw() && cu.BcwIdx == BCW_DEFAULT && cu.cs->slice->isInterB() )
          {
            const bool checkIdentical = true;
            m_uniMotions.setReadMode(checkIdentical, (uint32_t)iRefList, (uint32_t)iRefIdxTemp);
            m_uniMotions.copyFrom(cMvTemp[iRefList][iRefIdxTemp], uiCostTemp - m_pcRdCost->getCost(uiBitsTemp), (uint32_t)iRefList, (uint32_t)iRefIdxTemp);
          }
          xCopyAMVPInfo( &amvp[eRefPicList], &aacAMVPInfo[iRefList][iRefIdxTemp]); // must always be done ( also when AMVP_MODE = AM_NONE )
          xCheckBestMVP( eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], amvp[eRefPicList], uiBitsTemp, uiCostTemp, pu.cu->imv );

          if ( iRefList == 0 )
          {
            uiCostTempL0[iRefIdxTemp] = uiCostTemp;
            uiBitsTempL0[iRefIdxTemp] = uiBitsTemp;
          }
          if ( uiCostTemp < uiCost[iRefList] )
          {
            uiCost[iRefList] = uiCostTemp;
            uiBits[iRefList] = uiBitsTemp; // storing for bi-prediction

            // set motion
            cMv    [iRefList] = cMvTemp[iRefList][iRefIdxTemp];
            iRefIdx[iRefList] = iRefIdxTemp;
          }
#if JVET_Z0054_BLK_REF_PIC_REORDER
          if (cu.cs->sps->getUseARL() && iRefList == 1 && cs.slice->getList1IdxToList0Idx(iRefIdxTemp) >= 0)
          {
            uiCostTemp = MAX_UINT;
          }
#endif

          if ( iRefList == 1 && uiCostTemp < costValidList1 && cs.slice->getList1IdxToList0Idx( iRefIdxTemp ) < 0 )
          {
            costValidList1 = uiCostTemp;
            bitsValidList1 = uiBitsTemp;

            // set motion
            mvValidList1     = cMvTemp[iRefList][iRefIdxTemp];
            refIdxValidList1 = iRefIdxTemp;
          }
        }
      }

      ::memcpy(cMvHevcTemp, cMvTemp, sizeof(cMvTemp));
      if (cu.imv == 0 && (!cu.slice->getSPS()->getUseBcw() || bcwIdx == BCW_DEFAULT))
      {
        insertUniMvCands(pu.Y(), cMvTemp);

        unsigned idx1, idx2, idx3, idx4;
        getAreaIdx(cu.Y(), *cu.slice->getPPS()->pcv, idx1, idx2, idx3, idx4);
#if INTER_LIC
        if (cu.slice->getUseLIC() && cu.LICFlag)
        {
          ::memcpy(&(g_reusedUniMVsLIC[idx1][idx2][idx3][idx4][0][0]), cMvTemp, 2 * 33 * sizeof(Mv));
          g_isReusedUniMVsFilledLIC[idx1][idx2][idx3][idx4] = true;
        }
        else
        {
#endif
        ::memcpy(&(g_reusedUniMVs[idx1][idx2][idx3][idx4][0][0]), cMvTemp, 2 * 33 * sizeof(Mv));
        g_isReusedUniMVsFilled[idx1][idx2][idx3][idx4] = true;
#if INTER_LIC
        }
#endif
      }
#if JVET_X0083_BM_AMVP_MERGE_MODE
      }
#endif
      //  Bi-predictive Motion estimation
      if( ( cs.slice->isInterB() ) && ( PU::isBipredRestriction( pu ) == false )
        && (cu.slice->getCheckLDC() || bcwIdx == BCW_DEFAULT || !m_affineModeSelected || !m_pcEncCfg->getUseBcwFast())
#if INTER_LIC
        && !cu.LICFlag
#endif
        )
      {
        bool doBiPred = true;
        tryBipred = 1;
        cMvBi[0] = cMv[0];
        cMvBi[1] = cMv[1];
        iRefIdxBi[0] = iRefIdx[0];
        iRefIdxBi[1] = iRefIdx[1];

        ::memcpy( cMvPredBi,   cMvPred,   sizeof( cMvPred   ) );
        ::memcpy( aaiMvpIdxBi, aaiMvpIdx, sizeof( aaiMvpIdx ) );

        uint32_t uiMotBits[2];

#if JVET_X0083_BM_AMVP_MERGE_MODE
        if(cs.picHeader->getMvdL1ZeroFlag() && !pu.amvpMergeModeFlag[1])
#else
        if(cs.picHeader->getMvdL1ZeroFlag())
#endif
        {
          xCopyAMVPInfo(&aacAMVPInfo[1][bestBiPRefIdxL1], &amvp[REF_PIC_LIST_1]);
          aaiMvpIdxBi[1][bestBiPRefIdxL1] = bestBiPMvpL1;
          cMvPredBi  [1][bestBiPRefIdxL1] = amvp[REF_PIC_LIST_1].mvCand[bestBiPMvpL1];

          cMvBi    [1] = cMvPredBi[1][bestBiPRefIdxL1];
          iRefIdxBi[1] = bestBiPRefIdxL1;
          pu.mv    [REF_PIC_LIST_1] = cMvBi[1];
          pu.refIdx[REF_PIC_LIST_1] = iRefIdxBi[1];
          pu.mvpIdx[REF_PIC_LIST_1] = bestBiPMvpL1;

          if( m_pcEncCfg->getMCTSEncConstraint() )
          {
            Mv restrictedMv = pu.mv[REF_PIC_LIST_1];
            Area curTileAreaRestricted;
            curTileAreaRestricted = pu.cs->picture->mctsInfo.getTileAreaSubPelRestricted( pu );
            MCTSHelper::clipMvToArea( restrictedMv, pu.cu->Y(), curTileAreaRestricted, *pu.cs->sps );
            // If sub-pel filter samples are not inside of allowed area
            if( restrictedMv != pu.mv[REF_PIC_LIST_1] )
            {
              uiCostBi = std::numeric_limits<Distortion>::max();
              doBiPred = false;
            }
          }
          PelUnitBuf predBufTmp = m_tmpPredStorage[REF_PIC_LIST_1].getBuf( UnitAreaRelative(cu, pu) );
          motionCompensation( pu, predBufTmp, REF_PIC_LIST_1 );

          uiMotBits[0] = uiBits[0] - uiMbBits[0];
          uiMotBits[1] = uiMbBits[1];

          if ( cs.slice->getNumRefIdx(REF_PIC_LIST_1) > 1 )
          {
            uiMotBits[1] += bestBiPRefIdxL1 + 1;
            if ( bestBiPRefIdxL1 == cs.slice->getNumRefIdx(REF_PIC_LIST_1)-1 )
            {
              uiMotBits[1]--;
            }
          }

#if TM_AMVP
          uiMotBits[1] += m_auiMVPIdxCost[aaiMvpIdxBi[1][bestBiPRefIdxL1]][amvp[REF_PIC_LIST_1].numCand];
#else
          uiMotBits[1] += m_auiMVPIdxCost[aaiMvpIdxBi[1][bestBiPRefIdxL1]][AMVP_MAX_NUM_CANDS];
#endif

          uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];

          cMvTemp[1][bestBiPRefIdxL1] = cMvBi[1];
        }
        else
        {
          uiMotBits[0] = uiBits[0] - uiMbBits[0];
          uiMotBits[1] = uiBits[1] - uiMbBits[1];
          uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];
        }

        if( doBiPred )
        {
        // 4-times iteration (default)
        int iNumIter = 4;

        // fast encoder setting: only one iteration
        if ( m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE2 || cs.picHeader->getMvdL1ZeroFlag() )
        {
          iNumIter = 1;
        }
#if JVET_X0083_BM_AMVP_MERGE_MODE
        if (amvpMergeModeFlag)
        {
          iNumIter = 1;
        }
#endif

        enforceBcwPred = (bcwIdx != BCW_DEFAULT);
        for ( int iIter = 0; iIter < iNumIter; iIter++ )
        {
          int         iRefList    = iIter % 2;

#if JVET_X0083_BM_AMVP_MERGE_MODE
          if (amvpMergeModeFlag)
          {
            iRefList = pu.amvpMergeModeFlag[1] ? 0 : 1;
          }
          else
#endif
          if ( m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE2 )
          {
            if( uiCost[0] <= uiCost[1] )
            {
              iRefList = 1;
            }
            else
            {
              iRefList = 0;
            }
            if( bcwIdx != BCW_DEFAULT )
            {
              iRefList = ( abs( getBcwWeight(bcwIdx, REF_PIC_LIST_0 ) ) > abs( getBcwWeight(bcwIdx, REF_PIC_LIST_1 ) ) ? 1 : 0 );
            }
          }
          else if ( iIter == 0 )
          {
            iRefList = 0;
          }
#if JVET_X0083_BM_AMVP_MERGE_MODE
          if ( iIter == 0 && !cs.picHeader->getMvdL1ZeroFlag() && !amvpMergeModeFlag)
#else
          if ( iIter == 0 && !cs.picHeader->getMvdL1ZeroFlag())
#endif
          {
            pu.mv    [1 - iRefList] = cMv    [1 - iRefList];
            pu.refIdx[1 - iRefList] = iRefIdx[1 - iRefList];

            PelUnitBuf predBufTmp = m_tmpPredStorage[1 - iRefList].getBuf( UnitAreaRelative(cu, pu) );
            motionCompensation( pu, predBufTmp, RefPicList(1 - iRefList) );
          }

          RefPicList  eRefPicList = ( iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0 );

          if(cs.picHeader->getMvdL1ZeroFlag())
          {
            iRefList = 0;
            eRefPicList = REF_PIC_LIST_0;
          }

          bool bChanged = false;

          iRefStart = 0;
          iRefEnd   = cs.slice->getNumRefIdx(eRefPicList)-1;
          for (int iRefIdxTemp = iRefStart; iRefIdxTemp <= iRefEnd; iRefIdxTemp++)
          {
#if JVET_X0083_BM_AMVP_MERGE_MODE
          int numberBestMvpIdxLoop = 1;
          int selectedBestMvpIdx = -1;
          Mv selectedBestMv;
          if (amvpMergeModeFlag)
          {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
            if (mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS_MEM].refIdx < 0
                && mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS_MEM + 1].refIdx < 0
                && mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS_MEM + 2].refIdx < 0)
#else
            if (mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS].refIdx < 0 && mvFieldAmListCommon[iRefIdxTemp * AMVP_MAX_NUM_CANDS + 1].refIdx < 0)
#endif
            {
              continue;
            }
            xEstimateMvPredAMVP( pu, origBuf, refListAmvp, iRefIdxTemp, cMvPred[refListAmvp][iRefIdxTemp], amvp[refListAmvp], false, &biPDistTemp, mvFieldAmListCommon);
            xCopyAMVPInfo( &amvp[refListAmvp], &aacAMVPInfo[refListAmvp][iRefIdxTemp]); // must always be done ( also when AMVP_MODE = AM_NONE )
            numberBestMvpIdxLoop = amvp[eRefPicList].numCand;
          }
          for (int bestMvpIdxLoop = 0; bestMvpIdxLoop < numberBestMvpIdxLoop; bestMvpIdxLoop++)
          {
            if (amvpMergeModeFlag)
            {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
              const int mvFieldMergeIdx = iRefIdxTemp * AMVP_MAX_NUM_CANDS_MEM + bestMvpIdxLoop;
#else
              const int mvFieldMergeIdx = iRefIdxTemp * AMVP_MAX_NUM_CANDS + bestMvpIdxLoop;
#endif
              aaiMvpIdxBi[iRefList][iRefIdxTemp] = bestMvpIdxLoop;
              unsigned idx1, idx2, idx3, idx4;
              getAreaIdx(cu.Y(), *cu.slice->getPPS()->pcv, idx1, idx2, idx3, idx4);
              CHECK(g_isReusedUniMVsFilled[idx1][idx2][idx3][idx4] == false, "this is not possible");
              if (g_isReusedUniMVsFilled[idx1][idx2][idx3][idx4])
              {
                cMvTemp[iRefList][iRefIdxTemp] = g_reusedUniMVs[idx1][idx2][idx3][idx4][refListAmvp][iRefIdxTemp];
              }
              else
              {
                cMvTemp[iRefList][iRefIdxTemp] = amvp[eRefPicList].mvCand[bestMvpIdxLoop];
              }
              cMvPredBi[iRefList][iRefIdxTemp] = amvp[eRefPicList].mvCand[bestMvpIdxLoop];
              // set merge dir mv info and MC
              pu.mv[1 - iRefList] = mvFieldAmListCommon[mvFieldMergeIdx].mv;
              pu.refIdx[1 - iRefList] = mvFieldAmListCommon[mvFieldMergeIdx].refIdx;
            }
#endif
            if( m_pcEncCfg->getUseBcwFast() && (bcwIdx != BCW_DEFAULT)
              && (pu.cu->slice->getRefPic(eRefPicList, iRefIdxTemp)->getPOC() == pu.cu->slice->getRefPic(RefPicList(1 - iRefList), pu.refIdx[1 - iRefList])->getPOC())
              && (!pu.cu->imv && pu.cu->slice->getTLayer()>1)
#if INTER_LIC
              && !cu.LICFlag
#endif
              )
            {
              continue;
            }
#if JVET_X0083_BM_AMVP_MERGE_MODE
            if (amvpMergeModeFlag)
            {
              uiBitsTemp = uiMbBits[2];
            }
            else
            {
#endif
#if JVET_Z0054_BLK_REF_PIC_REORDER
              if (cu.cs->sps->getUseARL())
              {
                int refIdxTemp[2];
                refIdxTemp[iRefList] = iRefIdxTemp;
                refIdxTemp[1 - iRefList] = iRefIdxBi[1 - iRefList];
                if (pu.cu->slice->getRefPicPairIdx(refIdxTemp[0], refIdxTemp[1]) < 0)
                {
                  continue;
                }
              }
#endif
            uiBitsTemp = uiMbBits[2] + uiMotBits[1-iRefList];
            uiBitsTemp += ((cs.slice->getSPS()->getUseBcw() == true) ? getWeightIdxBits(bcwIdx) : 0);
#if JVET_X0083_BM_AMVP_MERGE_MODE
            }
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE
            if (( cs.slice->getNumRefIdx(eRefPicList) > 1 ) && !(amvpMergeModeFlag && candidateRefIdxCount <= 1))
#else
            if ( cs.slice->getNumRefIdx(eRefPicList) > 1 )
#endif
            {
              uiBitsTemp += iRefIdxTemp+1;
              if ( iRefIdxTemp == cs.slice->getNumRefIdx(eRefPicList)-1 )
              {
                uiBitsTemp--;
              }
            }
#if TM_AMVP
            uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdxBi[iRefList][iRefIdxTemp]][aacAMVPInfo[iRefList][iRefIdxTemp].numCand];
#else
            uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdxBi[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];
#endif
            if ( cs.slice->getBiDirPred() )
            {
#if JVET_X0083_BM_AMVP_MERGE_MODE
              if (!amvpMergeModeFlag)
#endif
              uiBitsTemp += 1; // add one bit for symmetrical MVD mode
            }
#if MULTI_HYP_PRED
            if (saveMeResultsForMHP)
              uiBitsTemp++; // terminating 0 mh_flag
#endif
            // call ME
#if JVET_X0083_BM_AMVP_MERGE_MODE
            if (amvpMergeModeFlag)
            {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
              if (bestMvpIdxLoop < 2)
              {
                MvField amvpMvField, mergeMvField;
                amvpMvField.setMvField(cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp);
                mergeMvField.setMvField(cMvPredBi[iRefList][iRefIdxTemp].getSymmvdMv(cMvPredBi[iRefList][iRefIdxTemp], pu.mv[1 - iRefList]), pu.refIdx[1 - iRefList]);
                uiCostTemp = xGetSymmetricCost( pu, origBuf, eRefPicList, amvpMvField, mergeMvField, bcwIdx );
                uiCostTemp += m_pcRdCost->getCost( uiBitsTemp );
                cMvTemp[iRefList][iRefIdxTemp] = amvpMvField.mv;
              }
              else
              {
#endif
              PelUnitBuf predBufTmp = m_tmpPredStorage[1 - iRefList].getBuf( UnitAreaRelative(cu, pu) );
              motionCompensation( pu, predBufTmp, RefPicList(1 - iRefList) );
#if MULTI_HYP_PRED
              CHECK(pu.addHypData.empty() == false, "this is not possible");
#endif
              xMotionEstimation ( pu, origBuf, eRefPicList, cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, amvp[eRefPicList], true );
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
              }
#endif
            }
            else
            {
#endif
            xCopyAMVPInfo(&aacAMVPInfo[iRefList][iRefIdxTemp], &amvp[eRefPicList] );
            xMotionEstimation ( pu, origBuf, eRefPicList, cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, amvp[eRefPicList], true );
            xCheckBestMVP( eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPredBi[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], amvp[eRefPicList], uiBitsTemp, uiCostTemp, pu.cu->imv);
#if JVET_X0083_BM_AMVP_MERGE_MODE
            }
#endif
#if MULTI_HYP_PRED
            if (saveMeResultsForMHP)
            {
              // AMVP bi
              MEResult biPredResult;
              biPredResult.cu = cu;
              biPredResult.pu = pu;
              biPredResult.pu.interDir = 3;
              biPredResult.pu.mv[iRefList] = cMvTemp[iRefList][iRefIdxTemp];
              biPredResult.pu.mv[1 - iRefList] = cMvBi[1 - iRefList];
              biPredResult.pu.mv[0].mvCliptoStorageBitDepth();
              biPredResult.pu.mv[1].mvCliptoStorageBitDepth();

              biPredResult.pu.mvd[iRefList] = cMvTemp[iRefList][iRefIdxTemp] - cMvPredBi[iRefList][iRefIdxTemp];
              biPredResult.pu.mvd[1 - iRefList] = cMvBi[1 - iRefList] - cMvPredBi[1 - iRefList][iRefIdxBi[1 - iRefList]];
              biPredResult.pu.refIdx[iRefList] = iRefIdxTemp;
              biPredResult.pu.refIdx[1 - iRefList] = iRefIdxBi[1 - iRefList];
              biPredResult.pu.mvpIdx[iRefList] = aaiMvpIdxBi[iRefList][iRefIdxTemp];
              biPredResult.pu.mvpIdx[1 - iRefList] = aaiMvpIdxBi[1 - iRefList][iRefIdxBi[1 - iRefList]];
              biPredResult.pu.mvpNum[iRefList] = aaiMvpNum[iRefList][iRefIdxTemp];
              biPredResult.pu.mvpNum[1 - iRefList] = aaiMvpNum[1 - iRefList][iRefIdxBi[1 - iRefList]];
              biPredResult.cost = uiCostTemp;
              biPredResult.bits = uiBitsTemp;

              if (!(cu.imv != 0 && biPredResult.pu.mvd[0] == Mv(0, 0) && biPredResult.pu.mvd[1] == Mv(0, 0)))
              {
                cs.m_meResults.push_back(biPredResult);
              }

            }
#endif
            if ( uiCostTemp < uiCostBi )
            {
              bChanged = true;

              cMvBi[iRefList]     = cMvTemp[iRefList][iRefIdxTemp];
              iRefIdxBi[iRefList] = iRefIdxTemp;
#if JVET_X0083_BM_AMVP_MERGE_MODE
              if (amvpMergeModeFlag)
              {
                selectedBestMvpIdx = bestMvpIdxLoop;
                selectedBestMv = cMvTemp[iRefList][iRefIdxTemp];
              }
#endif

              uiCostBi            = uiCostTemp;
#if JVET_X0083_BM_AMVP_MERGE_MODE
              if (amvpMergeModeFlag)
              {
                uiMotBits[iRefList] = uiBitsTemp - uiMbBits[2];
              }
              else
              {
#endif
              uiMotBits[iRefList] = uiBitsTemp - uiMbBits[2] - uiMotBits[1-iRefList];
              uiMotBits[iRefList] -= ((cs.slice->getSPS()->getUseBcw() == true) ? getWeightIdxBits(bcwIdx) : 0);
#if JVET_X0083_BM_AMVP_MERGE_MODE
              }
#endif
              uiBits[2]           = uiBitsTemp;

              if(iNumIter!=1)
              {
                //  Set motion
                pu.mv    [eRefPicList] = cMvBi    [iRefList];
                pu.refIdx[eRefPicList] = iRefIdxBi[iRefList];

                PelUnitBuf predBufTmp = m_tmpPredStorage[iRefList].getBuf( UnitAreaRelative(cu, pu) );
                motionCompensation( pu, predBufTmp, eRefPicList );
              }
            }
#if JVET_X0083_BM_AMVP_MERGE_MODE
          } // for loop-bestMvpIdxLoop

          if (amvpMergeModeFlag && selectedBestMvpIdx >= 0)
          {
            aaiMvpIdxBi[iRefList][iRefIdxTemp] = selectedBestMvpIdx;
            xCopyAMVPInfo(&aacAMVPInfo[iRefList][iRefIdxTemp], &amvp[eRefPicList] );
            cMvTemp[iRefList][iRefIdxTemp] = selectedBestMv;
            cMvPredBi[iRefList][iRefIdxTemp] = amvp[eRefPicList].mvCand[selectedBestMvpIdx];
          }
#endif
          } // for loop-iRefIdxTemp

          if ( !bChanged )
          {
            if ((uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1]) || enforceBcwPred)
            {
              xCopyAMVPInfo(&aacAMVPInfo[0][iRefIdxBi[0]], &amvp[REF_PIC_LIST_0]);
              xCheckBestMVP( REF_PIC_LIST_0, cMvBi[0], cMvPredBi[0][iRefIdxBi[0]], aaiMvpIdxBi[0][iRefIdxBi[0]], amvp[REF_PIC_LIST_0], uiBits[2], uiCostBi, pu.cu->imv);
              if(!cs.picHeader->getMvdL1ZeroFlag())
              {
                xCopyAMVPInfo(&aacAMVPInfo[1][iRefIdxBi[1]], &amvp[REF_PIC_LIST_1]);
                xCheckBestMVP( REF_PIC_LIST_1, cMvBi[1], cMvPredBi[1][iRefIdxBi[1]], aaiMvpIdxBi[1][iRefIdxBi[1]], amvp[REF_PIC_LIST_1], uiBits[2], uiCostBi, pu.cu->imv);
              }
            }
            break;
          }
        } // for loop-iter
        }
        cu.refIdxBi[0] = iRefIdxBi[0];
        cu.refIdxBi[1] = iRefIdxBi[1];

        if ( cs.slice->getBiDirPred() && trySmvd )
        {
          Distortion symCost;
          Mv cMvPredSym[2];
          int mvpIdxSym[2];

          int curRefList = REF_PIC_LIST_0;
          int tarRefList = 1 - curRefList;
          RefPicList eCurRefList = (curRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
          int refIdxCur = cs.slice->getSymRefIdx( curRefList );
          int refIdxTar = cs.slice->getSymRefIdx( tarRefList );
          CHECK (refIdxCur==-1 || refIdxTar==-1, "Uninitialized reference index not allowed");

          if ( aacAMVPInfo[curRefList][refIdxCur].mvCand[0] == aacAMVPInfo[curRefList][refIdxCur].mvCand[1] )
            aacAMVPInfo[curRefList][refIdxCur].numCand = 1;
          if ( aacAMVPInfo[tarRefList][refIdxTar].mvCand[0] == aacAMVPInfo[tarRefList][refIdxTar].mvCand[1] )
            aacAMVPInfo[tarRefList][refIdxTar].numCand = 1;

          MvField cCurMvField, cTarMvField;
          Distortion costStart = std::numeric_limits<Distortion>::max();
          for ( int i = 0; i < aacAMVPInfo[curRefList][refIdxCur].numCand; i++ )
          {
            for ( int j = 0; j < aacAMVPInfo[tarRefList][refIdxTar].numCand; j++ )
            {
              cCurMvField.setMvField( aacAMVPInfo[curRefList][refIdxCur].mvCand[i], refIdxCur );
              cTarMvField.setMvField( aacAMVPInfo[tarRefList][refIdxTar].mvCand[j], refIdxTar );
              Distortion cost = xGetSymmetricCost( pu, origBuf, eCurRefList, cCurMvField, cTarMvField, bcwIdx );
              if ( cost < costStart )
              {
                costStart = cost;
                cMvPredSym[curRefList] = aacAMVPInfo[curRefList][refIdxCur].mvCand[i];
                cMvPredSym[tarRefList] = aacAMVPInfo[tarRefList][refIdxTar].mvCand[j];
                mvpIdxSym[curRefList] = i;
                mvpIdxSym[tarRefList] = j;
              }
            }
          }
          cCurMvField.mv = cMvPredSym[curRefList];
          cTarMvField.mv = cMvPredSym[tarRefList];

          m_pcRdCost->setCostScale(0);
          Mv pred = cMvPredSym[curRefList];
          pred.changeTransPrecInternal2Amvr(pu.cu->imv);
          m_pcRdCost->setPredictor(pred);
          Mv mv = cCurMvField.mv;
          mv.changeTransPrecInternal2Amvr(pu.cu->imv);
          uint32_t bits = m_pcRdCost->getBitsOfVectorWithPredictor(mv.hor, mv.ver, 0);
#if TM_AMVP
          bits += m_auiMVPIdxCost[mvpIdxSym[curRefList]][aacAMVPInfo[curRefList][refIdxCur].numCand];
          bits += m_auiMVPIdxCost[mvpIdxSym[tarRefList]][aacAMVPInfo[tarRefList][refIdxTar].numCand];
#else
          bits += m_auiMVPIdxCost[mvpIdxSym[curRefList]][AMVP_MAX_NUM_CANDS];
          bits += m_auiMVPIdxCost[mvpIdxSym[tarRefList]][AMVP_MAX_NUM_CANDS];
#endif
          costStart += m_pcRdCost->getCost(bits);

          std::vector<Mv> symmvdCands;
          auto smmvdCandsGen = [&](Mv mvCand, bool mvPrecAdj)
          {
            if (mvPrecAdj && pu.cu->imv)
            {
              mvCand.roundTransPrecInternal2Amvr(pu.cu->imv);
            }

            bool toAddMvCand = true;
            for (std::vector<Mv>::iterator pos = symmvdCands.begin(); pos != symmvdCands.end(); pos++)
            {
              if (*pos == mvCand)
              {
                toAddMvCand = false;
                break;
              }
            }

            if (toAddMvCand)
            {
              symmvdCands.push_back(mvCand);
            }
          };

          smmvdCandsGen(cMvHevcTemp[curRefList][refIdxCur], false);
          smmvdCandsGen(cMvTemp[curRefList][refIdxCur], false);
          if (iRefIdxBi[curRefList] == refIdxCur)
          {
            smmvdCandsGen(cMvBi[curRefList], false);
          }
          for (int i = 0; i < m_uniMvListSize; i++)
          {
            if ( symmvdCands.size() >= 5 )
              break;
            BlkUniMvInfo* curMvInfo = m_uniMvList + ((m_uniMvListIdx - 1 - i + m_uniMvListMaxSize) % (m_uniMvListMaxSize));
            smmvdCandsGen(curMvInfo->uniMvs[curRefList][refIdxCur], true);
          }

          for (auto mvStart : symmvdCands)
          {
            bool checked = false; //if it has been checkin in the mvPred.
            for (int i = 0; i < aacAMVPInfo[curRefList][refIdxCur].numCand && !checked; i++)
            {
              checked |= (mvStart == aacAMVPInfo[curRefList][refIdxCur].mvCand[i]);
            }
            if (checked)
            {
              continue;
            }

            Distortion bestCost = costStart;
            symmvdCheckBestMvp(pu, origBuf, mvStart, (RefPicList)curRefList, aacAMVPInfo, bcwIdx, cMvPredSym, mvpIdxSym, costStart);
            if (costStart < bestCost)
            {
              cCurMvField.setMvField(mvStart, refIdxCur);
              cTarMvField.setMvField(mvStart.getSymmvdMv(cMvPredSym[curRefList], cMvPredSym[tarRefList]), refIdxTar);
            }
          }
          Mv startPtMv = cCurMvField.mv;

#if TM_AMVP
          Distortion mvpCost = m_pcRdCost->getCost(m_auiMVPIdxCost[mvpIdxSym[curRefList]][aacAMVPInfo[curRefList][refIdxCur].numCand]
                                                 + m_auiMVPIdxCost[mvpIdxSym[tarRefList]][aacAMVPInfo[tarRefList][refIdxTar].numCand]);
#else
          Distortion mvpCost = m_pcRdCost->getCost(m_auiMVPIdxCost[mvpIdxSym[curRefList]][AMVP_MAX_NUM_CANDS] + m_auiMVPIdxCost[mvpIdxSym[tarRefList]][AMVP_MAX_NUM_CANDS]);
#endif
          symCost = costStart - mvpCost;

          // ME
          xSymmetricMotionEstimation( pu, origBuf, cMvPredSym[curRefList], cMvPredSym[tarRefList], eCurRefList, cCurMvField, cTarMvField, symCost, bcwIdx );

          symCost += mvpCost;

          if (startPtMv != cCurMvField.mv)
          { // if ME change MV, run a final check for best MVP.
            symmvdCheckBestMvp(pu, origBuf, cCurMvField.mv, (RefPicList)curRefList, aacAMVPInfo, bcwIdx, cMvPredSym, mvpIdxSym, symCost, true);
          }

          bits = uiMbBits[2];
          bits += 1; // add one bit for #symmetrical MVD mode
          bits += ((cs.slice->getSPS()->getUseBcw() == true) ? getWeightIdxBits(bcwIdx) : 0);
          symCost += m_pcRdCost->getCost(bits);
          cTarMvField.setMvField(cCurMvField.mv.getSymmvdMv(cMvPredSym[curRefList], cMvPredSym[tarRefList]), refIdxTar);

          if( m_pcEncCfg->getMCTSEncConstraint() )
          {
            if( !( MCTSHelper::checkMvForMCTSConstraint( pu, cCurMvField.mv ) && MCTSHelper::checkMvForMCTSConstraint( pu, cTarMvField.mv ) ) )
              symCost = std::numeric_limits<Distortion>::max();
          }
#if MULTI_HYP_PRED
          if (saveMeResultsForMHP)
          {
            // SMVD
            MEResult biPredResult;
            biPredResult.cu = cu;
            biPredResult.pu = pu;
            biPredResult.pu.interDir = 3;

            biPredResult.cu.smvdMode = 1 + curRefList;

            biPredResult.pu.mv[curRefList] = cCurMvField.mv;
            biPredResult.pu.mv[tarRefList] = cTarMvField.mv;
            biPredResult.pu.mv[curRefList].mvCliptoStorageBitDepth();
            biPredResult.pu.mv[tarRefList].mvCliptoStorageBitDepth();
            biPredResult.pu.mvd[curRefList] = cCurMvField.mv - cMvPredSym[curRefList];
            biPredResult.pu.mvd[tarRefList] = cTarMvField.mv - cMvPredSym[tarRefList];
            biPredResult.pu.refIdx[curRefList] = cCurMvField.refIdx;
            biPredResult.pu.refIdx[tarRefList] = cTarMvField.refIdx;
            biPredResult.pu.mvpIdx[curRefList] = mvpIdxSym[curRefList];
            biPredResult.pu.mvpIdx[tarRefList] = mvpIdxSym[tarRefList];
            biPredResult.pu.mvpNum[curRefList] = aaiMvpNum[curRefList][cCurMvField.refIdx];
            biPredResult.pu.mvpNum[tarRefList] = aaiMvpNum[tarRefList][cTarMvField.refIdx];

            biPredResult.cost = symCost;
            biPredResult.bits = bits;

            if (!(cu.imv != 0 && biPredResult.pu.mvd[0] == Mv(0, 0) && biPredResult.pu.mvd[1] == Mv(0, 0)))
            {
              cs.m_meResults.push_back(biPredResult);
            }

          }
#endif
          // save results
          if ( symCost < uiCostBi )
          {
            uiCostBi = symCost;
            symMode = 1 + curRefList;

            cMvBi[curRefList] = cCurMvField.mv;
            iRefIdxBi[curRefList] = cCurMvField.refIdx;
            aaiMvpIdxBi[curRefList][cCurMvField.refIdx] = mvpIdxSym[curRefList];
            cMvPredBi[curRefList][iRefIdxBi[curRefList]] = cMvPredSym[curRefList];

            cMvBi[tarRefList] = cTarMvField.mv;
            iRefIdxBi[tarRefList] = cTarMvField.refIdx;
            aaiMvpIdxBi[tarRefList][cTarMvField.refIdx] = mvpIdxSym[tarRefList];
            cMvPredBi[tarRefList][iRefIdxBi[tarRefList]] = cMvPredSym[tarRefList];
          }
        }
      } // if (B_SLICE)



      //  Clear Motion Field
    pu.mv    [REF_PIC_LIST_0] = Mv();
    pu.mv    [REF_PIC_LIST_1] = Mv();
    pu.mvd   [REF_PIC_LIST_0] = cMvZero;
    pu.mvd   [REF_PIC_LIST_1] = cMvZero;
    pu.refIdx[REF_PIC_LIST_0] = NOT_VALID;
    pu.refIdx[REF_PIC_LIST_1] = NOT_VALID;
    pu.mvpIdx[REF_PIC_LIST_0] = NOT_VALID;
    pu.mvpIdx[REF_PIC_LIST_1] = NOT_VALID;
    pu.mvpNum[REF_PIC_LIST_0] = NOT_VALID;
    pu.mvpNum[REF_PIC_LIST_1] = NOT_VALID;


    // Set Motion Field

    cMv    [1] = mvValidList1;
    iRefIdx[1] = refIdxValidList1;
    uiBits [1] = bitsValidList1;
    uiCost [1] = costValidList1;
    if (cu.cs->pps->getWPBiPred() == true && tryBipred && (bcwIdx != BCW_DEFAULT))
    {
      CHECK(iRefIdxBi[0]<0, "Invalid picture reference index");
      CHECK(iRefIdxBi[1]<0, "Invalid picture reference index");
      wp0 = cu.cs->slice->getWpScaling(REF_PIC_LIST_0, iRefIdxBi[0]);
      wp1 = cu.cs->slice->getWpScaling(REF_PIC_LIST_1, iRefIdxBi[1]);
      if (WPScalingParam::isWeighted(wp0) || WPScalingParam::isWeighted(wp1))
      {
        uiCostBi = MAX_UINT;
        enforceBcwPred = false;
      }
    }
    if( enforceBcwPred )
    {
      uiCost[0] = uiCost[1] = MAX_UINT;
    }

      uiLastModeTemp = uiLastMode;
#if JVET_X0083_BM_AMVP_MERGE_MODE
      if (amvpMergeModeFlag)
      {
        if (uiCostBi > ((m_amvpOnlyCost * 5) >> 2))
        {
#if JVET_Y0128_NON_CTC || JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
          m_skipPROF = false;
          m_encOnly = false;
#endif
          bdmvrAmMergeNotValid = true;
          return;
        }
        m_amvpOnlyCost = (uiCostBi < m_amvpOnlyCost) ? uiCostBi : m_amvpOnlyCost;
      }
      if (((uiCostBi <= uiCost[0]) && (uiCostBi <= uiCost[1])) || amvpMergeModeFlag)
      {
        uiLastMode = 2;
        if (pu.amvpMergeModeFlag[1])
        {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
          const int mvFieldMergeIdx = iRefIdxBi[0] * AMVP_MAX_NUM_CANDS_MEM + aaiMvpIdxBi[0][iRefIdxBi[0]];
#else
          const int mvFieldMergeIdx = iRefIdxBi[0] * AMVP_MAX_NUM_CANDS + aaiMvpIdxBi[0][iRefIdxBi[0]];
#endif
          pu.mv[REF_PIC_LIST_1] = mvFieldAmListCommon[mvFieldMergeIdx].mv;
          pu.refIdx[REF_PIC_LIST_1] = mvFieldAmListCommon[mvFieldMergeIdx].refIdx;
          pu.mvpIdx[REF_PIC_LIST_1] = 2;
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
          pu.mvd[REF_PIC_LIST_1].setZero();
#endif
        }
        if (pu.amvpMergeModeFlag[0])
        {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
          const int mvFieldMergeIdx = iRefIdxBi[1] * AMVP_MAX_NUM_CANDS_MEM + aaiMvpIdxBi[1][iRefIdxBi[1]];
#else
          const int mvFieldMergeIdx = iRefIdxBi[1] * AMVP_MAX_NUM_CANDS + aaiMvpIdxBi[1][iRefIdxBi[1]];
#endif
          pu.mv[REF_PIC_LIST_0] = mvFieldAmListCommon[mvFieldMergeIdx].mv;
          pu.refIdx[REF_PIC_LIST_0] = mvFieldAmListCommon[mvFieldMergeIdx].refIdx;
          pu.mvpIdx[REF_PIC_LIST_0] = 2;
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
          pu.mvd[REF_PIC_LIST_0].setZero();
#endif
        }
        pu.interDir = 3;
        if (!pu.amvpMergeModeFlag[0])
        {
          pu.mv[REF_PIC_LIST_0] = cMvBi[0];
          pu.refIdx[REF_PIC_LIST_0] = iRefIdxBi[0];
          pu.mvd[REF_PIC_LIST_0] = cMvBi[0] - cMvPredBi[0][iRefIdxBi[0]];
          pu.mvpIdx[REF_PIC_LIST_0] = aaiMvpIdxBi[0][iRefIdxBi[0]];
          pu.mvpNum[REF_PIC_LIST_0] = aaiMvpNum[0][iRefIdxBi[0]];
        }
        if (!pu.amvpMergeModeFlag[1])
        {
          pu.mv[REF_PIC_LIST_1] = cMvBi[1];
          pu.refIdx[REF_PIC_LIST_1] = iRefIdxBi[1];
          pu.mvd[REF_PIC_LIST_1] = cMvBi[1] - cMvPredBi[1][iRefIdxBi[1]];
          pu.mvpIdx[REF_PIC_LIST_1] = aaiMvpIdxBi[1][iRefIdxBi[1]];
          pu.mvpNum[REF_PIC_LIST_1] = aaiMvpNum[1][iRefIdxBi[1]];
        }
        pu.cu->smvdMode = symMode;
      }
#else
      if ( uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1])
      {
        uiLastMode = 2;
        pu.mv    [REF_PIC_LIST_0] = cMvBi[0];
        pu.mv    [REF_PIC_LIST_1] = cMvBi[1];
        pu.mvd   [REF_PIC_LIST_0] = cMvBi[0] - cMvPredBi[0][iRefIdxBi[0]];
        pu.mvd   [REF_PIC_LIST_1] = cMvBi[1] - cMvPredBi[1][iRefIdxBi[1]];
        pu.refIdx[REF_PIC_LIST_0] = iRefIdxBi[0];
        pu.refIdx[REF_PIC_LIST_1] = iRefIdxBi[1];
        pu.mvpIdx[REF_PIC_LIST_0] = aaiMvpIdxBi[0][iRefIdxBi[0]];
        pu.mvpIdx[REF_PIC_LIST_1] = aaiMvpIdxBi[1][iRefIdxBi[1]];
        pu.mvpNum[REF_PIC_LIST_0] = aaiMvpNum[0][iRefIdxBi[0]];
        pu.mvpNum[REF_PIC_LIST_1] = aaiMvpNum[1][iRefIdxBi[1]];
        pu.interDir = 3;

        pu.cu->smvdMode = symMode;
      }
#endif
      else if ( uiCost[0] <= uiCost[1] )
      {
        uiLastMode = 0;
        pu.mv    [REF_PIC_LIST_0] = cMv[0];
        pu.mvd   [REF_PIC_LIST_0] = cMv[0] - cMvPred[0][iRefIdx[0]];
        pu.refIdx[REF_PIC_LIST_0] = iRefIdx[0];
        pu.mvpIdx[REF_PIC_LIST_0] = aaiMvpIdx[0][iRefIdx[0]];
        pu.mvpNum[REF_PIC_LIST_0] = aaiMvpNum[0][iRefIdx[0]];
        pu.interDir = 1;
      }
      else
      {
        uiLastMode = 1;
        pu.mv    [REF_PIC_LIST_1] = cMv[1];
        pu.mvd   [REF_PIC_LIST_1] = cMv[1] - cMvPred[1][iRefIdx[1]];
        pu.refIdx[REF_PIC_LIST_1] = iRefIdx[1];
        pu.mvpIdx[REF_PIC_LIST_1] = aaiMvpIdx[1][iRefIdx[1]];
        pu.mvpNum[REF_PIC_LIST_1] = aaiMvpNum[1][iRefIdx[1]];
        pu.interDir = 2;
      }

      if( bcwIdx != BCW_DEFAULT )
      {
        cu.BcwIdx = BCW_DEFAULT; // Reset to default for the Non-NormalMC modes.
      }

    uiHevcCost = ( uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1] ) ? uiCostBi : ( ( uiCost[0] <= uiCost[1] ) ? uiCost[0] : uiCost[1] );
#if JVET_X0083_BM_AMVP_MERGE_MODE
    if (!amvpMergeModeFlag && (m_amvpOnlyCost > uiHevcCost))
    {
      m_amvpOnlyCost = uiHevcCost;
    }
#endif
    }
#if INTER_RM_SIZE_CONSTRAINTS
    if (cu.Y().width >= 8 && cu.Y().height >= 8 && cu.slice->getSPS()->getUseAffine()
#else
    if (cu.Y().width > 8 && cu.Y().height > 8 && cu.slice->getSPS()->getUseAffine()
#endif
      && checkAffine
      && m_pcEncCfg->getUseAffineAmvp()
      && (bcwIdx == BCW_DEFAULT || m_affineModeSelected || !m_pcEncCfg->getUseBcwFast())
#if JVET_X0083_BM_AMVP_MERGE_MODE
      && !amvpMergeModeFlag
#endif
      )
    {
      m_hevcCost = uiHevcCost;
      // save normal hevc result
      uint32_t uiMRGIndex = pu.mergeIdx;
      bool bMergeFlag = pu.mergeFlag;
      uint32_t uiInterDir = pu.interDir;
      int  iSymMode = cu.smvdMode;

      Mv cMvd[2];
      uint32_t uiMvpIdx[2], uiMvpNum[2];
      uiMvpIdx[0] = pu.mvpIdx[REF_PIC_LIST_0];
      uiMvpIdx[1] = pu.mvpIdx[REF_PIC_LIST_1];
      uiMvpNum[0] = pu.mvpNum[REF_PIC_LIST_0];
      uiMvpNum[1] = pu.mvpNum[REF_PIC_LIST_1];
      cMvd[0]     = pu.mvd[REF_PIC_LIST_0];
      cMvd[1]     = pu.mvd[REF_PIC_LIST_1];

      MvField cHevcMvField[2];
      cHevcMvField[0].setMvField( pu.mv[REF_PIC_LIST_0], pu.refIdx[REF_PIC_LIST_0] );
      cHevcMvField[1].setMvField( pu.mv[REF_PIC_LIST_1], pu.refIdx[REF_PIC_LIST_1] );

      // do affine ME & Merge
      cu.affineType = AFFINEMODEL_4PARAM;
      Mv acMvAffine4Para[2][33][3];
      int refIdx4Para[2] = { -1, -1 };

      xPredAffineInterSearch(pu, origBuf, puIdx, uiLastModeTemp, uiAffineCost, cMvHevcTemp, acMvAffine4Para, refIdx4Para, bcwIdx, enforceBcwPred,
        ((cu.slice->getSPS()->getUseBcw() == true) ? getWeightIdxBits(bcwIdx) : 0));

      if ( pu.cu->imv == 0 )
      {
        storeAffineMotion( pu.mvAffi, pu.refIdx, AFFINEMODEL_4PARAM, bcwIdx );
      }

      if ( cu.slice->getSPS()->getUseAffineType() )
      {
#if AFFINE_ENC_OPT
        if (uiAffineCost < uiHevcCost * 0.95) ///< condition for 6 parameter affine ME
#else
        if (uiAffineCost < uiHevcCost * 1.05) ///< condition for 6 parameter affine ME
#endif
        {
          // save 4 parameter results
          Mv bestMv[2][3], bestMvd[2][3];
          int bestMvpIdx[2], bestMvpNum[2], bestRefIdx[2];
          uint8_t bestInterDir;

          bestInterDir = pu.interDir;
          bestRefIdx[0] = pu.refIdx[0];
          bestRefIdx[1] = pu.refIdx[1];
          bestMvpIdx[0] = pu.mvpIdx[0];
          bestMvpIdx[1] = pu.mvpIdx[1];
          bestMvpNum[0] = pu.mvpNum[0];
          bestMvpNum[1] = pu.mvpNum[1];

          for ( int refList = 0; refList < 2; refList++ )
          {
            bestMv[refList][0] = pu.mvAffi[refList][0];
            bestMv[refList][1] = pu.mvAffi[refList][1];
            bestMv[refList][2] = pu.mvAffi[refList][2];
            bestMvd[refList][0] = pu.mvdAffi[refList][0];
            bestMvd[refList][1] = pu.mvdAffi[refList][1];
            bestMvd[refList][2] = pu.mvdAffi[refList][2];
          }

          refIdx4Para[0] = bestRefIdx[0];
          refIdx4Para[1] = bestRefIdx[1];

          Distortion uiAffine6Cost = std::numeric_limits<Distortion>::max();
          cu.affineType = AFFINEMODEL_6PARAM;
          xPredAffineInterSearch(pu, origBuf, puIdx, uiLastModeTemp, uiAffine6Cost, cMvHevcTemp, acMvAffine4Para, refIdx4Para, bcwIdx, enforceBcwPred,
            ((cu.slice->getSPS()->getUseBcw() == true) ? getWeightIdxBits(bcwIdx) : 0));

          if ( pu.cu->imv == 0 )
          {
            storeAffineMotion( pu.mvAffi, pu.refIdx, AFFINEMODEL_6PARAM, bcwIdx );
          }

          // reset to 4 parameter affine inter mode
          if ( uiAffineCost <= uiAffine6Cost )
          {
            cu.affineType = AFFINEMODEL_4PARAM;
            pu.interDir = bestInterDir;
            pu.refIdx[0] = bestRefIdx[0];
            pu.refIdx[1] = bestRefIdx[1];
            pu.mvpIdx[0] = bestMvpIdx[0];
            pu.mvpIdx[1] = bestMvpIdx[1];
            pu.mvpNum[0] = bestMvpNum[0];
            pu.mvpNum[1] = bestMvpNum[1];
            pu.mv[0].setZero();
            pu.mv[1].setZero();

            for ( int verIdx = 0; verIdx < 3; verIdx++ )
            {
              pu.mvdAffi[REF_PIC_LIST_0][verIdx] = bestMvd[0][verIdx];
              pu.mvdAffi[REF_PIC_LIST_1][verIdx] = bestMvd[1][verIdx];
              pu.mvAffi[REF_PIC_LIST_0][verIdx] = bestMv[0][verIdx];
              pu.mvAffi[REF_PIC_LIST_1][verIdx] = bestMv[1][verIdx];
            }
          }
          else
          {
            uiAffineCost = uiAffine6Cost;
          }
        }

        uiAffineCost += m_pcRdCost->getCost( 1 ); // add one bit for affine_type
      }

      if( uiAffineCost < uiHevcCost )
      {
        if( m_pcEncCfg->getMCTSEncConstraint() && !MCTSHelper::checkMvBufferForMCTSConstraint( pu ) )
        {
          uiAffineCost = std::numeric_limits<Distortion>::max();
        }
      }
      if ( uiHevcCost <= uiAffineCost )
      {
        // set hevc me result
        cu.affine = false;
        pu.mergeFlag = bMergeFlag;
        pu.regularMergeFlag = false;
        pu.mergeIdx = uiMRGIndex;
        pu.interDir = uiInterDir;
        cu.smvdMode = iSymMode;
        pu.mv    [REF_PIC_LIST_0] = cHevcMvField[0].mv;
        pu.refIdx[REF_PIC_LIST_0] = cHevcMvField[0].refIdx;
        pu.mv    [REF_PIC_LIST_1] = cHevcMvField[1].mv;
        pu.refIdx[REF_PIC_LIST_1] = cHevcMvField[1].refIdx;
        pu.mvpIdx[REF_PIC_LIST_0] = uiMvpIdx[0];
        pu.mvpIdx[REF_PIC_LIST_1] = uiMvpIdx[1];
        pu.mvpNum[REF_PIC_LIST_0] = uiMvpNum[0];
        pu.mvpNum[REF_PIC_LIST_1] = uiMvpNum[1];
        pu.mvd[REF_PIC_LIST_0] = cMvd[0];
        pu.mvd[REF_PIC_LIST_1] = cMvd[1];
      }
      else
      {
        cu.smvdMode = 0;
        CHECK( !cu.affine, "Wrong." );
        uiLastMode = uiLastModeTemp;
      }
    }
    if( cu.firstPU->interDir == 3 && !cu.firstPU->mergeFlag )
    {
      if (bcwIdx != BCW_DEFAULT)
      {
        cu.BcwIdx = bcwIdx;
      }
    }
    m_maxCompIDToPred = MAX_NUM_COMPONENT;

#if JVET_X0083_BM_AMVP_MERGE_MODE
    if (amvpMergeModeFlag && PU::checkBDMVRCondition(pu))
    {
      setBdmvrSubPuMvBuf(mvBufEncAmBDMVR_L0, mvBufEncAmBDMVR_L1);
      pu.bdmvrRefine = true;
      // span motion to subPU
      for (int subPuIdx = 0; subPuIdx < MAX_NUM_SUBCU_DMVR; subPuIdx++)
      {
        mvBufEncAmBDMVR_L0[subPuIdx] = pu.mv[0];
        mvBufEncAmBDMVR_L1[subPuIdx] = pu.mv[1];
      }
    }
    if (!pu.bdmvrRefine)
#endif
    {
      PU::spanMotionInfo( pu, mergeCtx );
    }

    m_skipPROF = false;
    m_encOnly = false;
    //  MC
    PelUnitBuf predBuf = pu.cs->getPredBuf(pu);
#if JVET_X0083_BM_AMVP_MERGE_MODE
    if (( bcwIdx == BCW_DEFAULT || !m_affineMotion.affine4ParaAvail || !m_affineMotion.affine6ParaAvail ) && !amvpMergeModeFlag)
#else
    if ( bcwIdx == BCW_DEFAULT || !m_affineMotion.affine4ParaAvail || !m_affineMotion.affine6ParaAvail )
#endif
    {
      m_affineMotion.hevcCost[pu.cu->imv] = uiHevcCost;
    }
#if INTER_LIC
    if (cu.LICFlag)
    {
#if !TM_AMVP
      m_storeBeforeLIC = true;
#endif
      m_predictionBeforeLIC = m_tmpStorageLCU.getBuf(UnitAreaRelative(*pu.cu, pu));
    }
#endif
    motionCompensation( pu, predBuf, REF_PIC_LIST_X );
#if JVET_X0083_BM_AMVP_MERGE_MODE
    if (pu.bdmvrRefine)
    {
      PU::spanMotionInfo( *cu.firstPU, MergeCtx(), mvBufEncAmBDMVR_L0, mvBufEncAmBDMVR_L1, getBdofSubPuMvOffset() );
    }
#endif
#if INTER_LIC && !TM_AMVP
    m_storeBeforeLIC = false;
#endif
    puIdx++;
  }

#if INTER_LIC
#if !TM_AMVP // This LIC optimization must be off; otherwise, enc/dec mismatching will result. Because the cost metrics (MRSAD or SAD) of TM mode is adaptive to LIC flag, refined MVs would change when LIC flag is 1 or 0.
  if (cu.LICFlag && pu.interDir != 10) // xCheckRDCostInterIMV initializes pu.interDir by using 10. When checkAffine and checkNonAffine are both false, pu.interDir remains 10 which should be avoided
  {
    CHECK(pu.interDir != 1 && pu.interDir != 2, "Invalid InterDir for LIC");

    PelUnitBuf predBuf = pu.cs->getPredBuf(pu);
    DistParam distParam;
    m_pcRdCost->setDistParam(distParam, cs.getOrgBuf().Y(), predBuf.Y(), cs.sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, true);
    Distortion distLicOn = distParam.distFunc(distParam);

    m_pcRdCost->setDistParam(distParam, cs.getOrgBuf().Y(), m_predictionBeforeLIC.Y(), cs.sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, true);
    Distortion distLicOff = distParam.distFunc(distParam);
    if (distLicOn >= distLicOff)
    {
      pu.cu->LICFlag = false;
      PU::spanLICFlags(pu, false);
      predBuf.copyFrom(m_predictionBeforeLIC);
    }
  }
#endif
#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if JVET_Z0054_BLK_REF_PIC_REORDER
  if (cu.imv && !CU::hasSubCUNonZeroMVd(cu) && !CU::hasSubCUNonZeroAffineMVd(cu))
  {
    setWpScalingDistParam(-1, REF_PIC_LIST_X, cu.cs->slice);
    return;
  }
  if(!PU::useRefPairList(pu) && !PU::useRefCombList(pu))
#endif
  if (pu.isMvsdApplicable())
  {
    bool bi = pu.interDir == 3;
    if (cu.affine)
    {
      for (uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
      {
        RefPicList eRefPicList = RefPicList(uiRefListIdx);
        Mv absMvd[3];
        absMvd[0] = Mv(pu.mvdAffi[uiRefListIdx][0].getAbsMv());
        absMvd[1] = Mv(pu.mvdAffi[uiRefListIdx][1].getAbsMv());
        absMvd[2] = (cu.affineType == AFFINEMODEL_6PARAM) ? Mv(pu.mvdAffi[uiRefListIdx][2].getAbsMv()) : Mv(0, 0);
        if (pu.cs->slice->getNumRefIdx(eRefPicList) > 0
          && (pu.interDir & (1 << uiRefListIdx)) && (absMvd[0] != Mv(0, 0) || absMvd[1] != Mv(0, 0) || absMvd[2] != Mv(0, 0)) && pu.isMvsdApplicable())
        {
          AffineAMVPInfo affineAMVPInfo;
          PU::fillAffineMvpCand(pu, eRefPicList, pu.refIdx[uiRefListIdx], affineAMVPInfo);
          const unsigned mvpIdx = pu.mvpIdx[eRefPicList];

          std::vector<Mv> cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3;
#if JVET_Z0054_BLK_REF_PIC_REORDER
          deriveMvdSignAffine(affineAMVPInfo.mvCandLT[mvpIdx], affineAMVPInfo.mvCandRT[mvpIdx], affineAMVPInfo.mvCandLB[mvpIdx],
            absMvd, pu, eRefPicList, pu.refIdx[eRefPicList], cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3);
#else
          deriveMvdSignAffine(affineAMVPInfo.mvCandLT[mvpIdx], affineAMVPInfo.mvCandRT[mvpIdx], affineAMVPInfo.mvCandLB[mvpIdx],
            absMvd[0], absMvd[1], absMvd[2], pu, eRefPicList, pu.refIdx[eRefPicList], cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3);
#endif
          int idx = -1;
          idx = deriveMVSDIdxFromMVDAffine(pu, eRefPicList, cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3);
          CHECK(idx == -1, "no match for mvsdIdx at Encoder");
          pu.mvsdIdx[eRefPicList] = idx;
        }
      }
    }
    else
    {
      for (uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
      {
        RefPicList eRefPicList = RefPicList(uiRefListIdx);
        Mv cMvd = pu.mvd[eRefPicList];
        if (pu.cs->slice->getNumRefIdx(eRefPicList) > 0
          && (pu.interDir & (1 << uiRefListIdx))
          && pu.isMvsdApplicable()
          && cMvd.isMvsdApplicable()
          )
        {
          auto aMvPred = bi ? cMvPredBi : cMvPred;
          auto aRefIdx = bi ? iRefIdxBi : iRefIdx;
          auto aMv = bi ? cMvBi : cMv;
          Mv cMvPred2 = aMvPred[uiRefListIdx][aRefIdx[uiRefListIdx]];
          CHECK(cMvd != aMv[uiRefListIdx] - cMvPred2, "");
          int iRefIdx = pu.refIdx[uiRefListIdx];
          Mv cMvdKnownAtDecoder = Mv(cMvd.getAbsHor(), cMvd.getAbsVer());
          std::vector<Mv> cMvdDerivedVec;
          if (cu.smvdMode)
          {
            if (uiRefListIdx == 1)
            {
              cMvd = pu.mvd[REF_PIC_LIST_0];
              CHECK((pu.mvd[REF_PIC_LIST_0].hor != -pu.mvd[REF_PIC_LIST_1].hor) || (pu.mvd[REF_PIC_LIST_0].ver != -pu.mvd[REF_PIC_LIST_1].ver), "not mirrored MVD for SMVD at Enc");
              CHECK(cs.slice->getSymRefIdx(REF_PIC_LIST_0) != pu.refIdx[REF_PIC_LIST_0], "ref Idx for List 0 does not match for SMVD at Enc");
              CHECK(cs.slice->getSymRefIdx(REF_PIC_LIST_1) != pu.refIdx[REF_PIC_LIST_1], "ref Idx for List 1 does not match for SMVD at Enc");

              deriveMvdSignSMVD(aMvPred[0][aRefIdx[0]], aMvPred[1][aRefIdx[1]], cMvdKnownAtDecoder, pu, cMvdDerivedVec);
              int idx = deriveMVSDIdxFromMVDTrans(cMvd, cMvdDerivedVec);
              CHECK(idx == -1, "");
              pu.mvsdIdx[REF_PIC_LIST_0] = idx;
            }
          }
          else
          {
            deriveMvdSign(cMvPred2, cMvdKnownAtDecoder, pu, eRefPicList, iRefIdx, cMvdDerivedVec);
            int idx = deriveMVSDIdxFromMVDTrans(cMvd, cMvdDerivedVec);
            CHECK(idx == -1, "");
            pu.mvsdIdx[eRefPicList] = idx;
          }
        }
      } //loop end for non-affine
    }
  }
#endif
  setWpScalingDistParam( -1, REF_PIC_LIST_X, cu.cs->slice );

  return;
}

uint32_t InterSearch::xCalcAffineMVBits( PredictionUnit& pu, Mv acMvTemp[3], Mv acMvPred[3] )
{
  int mvNum  = pu.cu->affineType ? 3 : 2;
  m_pcRdCost->setCostScale( 0 );
  uint32_t bitsTemp = 0;

  for ( int verIdx = 0; verIdx < mvNum; verIdx++ )
  {
    Mv pred = verIdx == 0 ? acMvPred[verIdx] : acMvPred[verIdx] + acMvTemp[0] - acMvPred[0];
    pred.changeAffinePrecInternal2Amvr(pu.cu->imv);
    m_pcRdCost->setPredictor( pred );
    Mv mv = acMvTemp[verIdx];
    mv.changeAffinePrecInternal2Amvr(pu.cu->imv);

    bitsTemp += m_pcRdCost->getBitsOfVectorWithPredictor( mv.getHor(), mv.getVer(), 0 );
  }

  return bitsTemp;
}

#if MULTI_HYP_PRED
void InterSearch::predInterSearchAdditionalHypothesis(PredictionUnit& pu, const MEResult& x, MEResultVec& out)
{
  const SPS &sps = *pu.cs->sps;
  CHECK(!sps.getUseInterMultiHyp(), "Multi Hyp is not active");
  CHECK(!pu.cs->slice->isInterB(), "Multi Hyp only allowed in B slices");
  CHECK(pu.cu->predMode != MODE_INTER, "Multi Hyp: pu.cu->predMode != MODE_INTER");
  CHECK(pu.addHypData.size() > sps.getMaxNumAddHyps(), "Multi Hyp: too many hypotheseis");
  if( pu.addHypData.size() == sps.getMaxNumAddHyps() )
  {
    return;
  }

  CHECK(!pu.mergeFlag && pu.cu->BcwIdx == BCW_DEFAULT, "!pu.mergeFlag && pu.cu->BcwIdx == BCW_DEFAULT");
  // get first prediction hypothesis
  PelUnitBuf tempPredBuf;
  if (x.predBuf != nullptr)
  {
    tempPredBuf = *x.predBuf;
  }
  else
  {
    tempPredBuf = pu.cs->getPredBuf(pu);
    pu.mvRefine = true;
    motionCompensation(pu, REF_PIC_LIST_X, true, false);
    pu.mvRefine = false;
  }
  const auto &MHRefPics = pu.cs->slice->getMultiHypRefPicList();
  const int iNumMHRefPics = int(MHRefPics.size());
  CHECK(iNumMHRefPics <= 0, "Multi Hyp: iNumMHRefPics <= 0");

  PelUnitBuf origBuf = pu.cs->getOrgBuf(pu);

  const UnitArea unitAreaFromPredBuf(origBuf.chromaFormat, Area(Position(0, 0), origBuf.Y()));
  // NOTE: tempOrigBuf share the same buffer with tempBuf that is used in xAddHypMC.
  PelUnitBuf tempOrigBuf = m_additionalHypothesisStorage.getBuf(unitAreaFromPredBuf);


  MultiHypPredictionData tempMHPredData;

  m_pcRdCost->selectMotionLambda();

  const int numWeights = sps.getNumAddHypWeights();
  unsigned idx1, idx2, idx3, idx4;
  getAreaIdx(pu.Y(), *pu.cs->slice->getPPS()->pcv, idx1, idx2, idx3, idx4);
#if INTER_LIC
  auto savedLICFlag = pu.cu->LICFlag;
#endif
  tempMHPredData.isMrg = true;
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
  uint8_t maxNumMergeCandidates = pu.cs->sps->getMaxNumMHPCand();
  CHECK(maxNumMergeCandidates >= GEO_MAX_NUM_UNI_CANDS, "");
  if (maxNumMergeCandidates > 0)
  {
#else
  {
    uint8_t maxNumMergeCandidates = pu.cs->sps->getMaxNumGeoCand();
    CHECK(maxNumMergeCandidates >= GEO_MAX_NUM_UNI_CANDS, "");
#endif
    DistParam distParam;
    const bool bUseHadamard = !pu.cs->slice->getDisableSATDForRD();
    m_pcRdCost->setDistParam(distParam, origBuf.Y(), tempOrigBuf.Y(), sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, bUseHadamard);

    if (!(pu.addHypData.size() > pu.numMergedAddHyps && m_mhpMrgTempBufSet)) // non 1st addHyp check should already have the MC results stored
    {
      PredictionUnit fakePredData = pu;
      fakePredData.mergeFlag = false;
      fakePredData.mergeType = MRG_TYPE_DEFAULT_N;
      fakePredData.mmvdMergeFlag = false;
      fakePredData.ciipFlag = false;
      fakePredData.addHypData.clear();
      fakePredData.regularMergeFlag = false;
#if TM_MRG
      fakePredData.tmMergeFlag = false;
#endif
#if JVET_X0049_ADAPT_DMVR
      fakePredData.bmMergeFlag = false;
#endif
#if MULTI_PASS_DMVR
      fakePredData.bdmvrRefine = false;
#endif
      if (!m_mhpMrgTempBufSet)
      {
        PU::getGeoMergeCandidates(fakePredData, m_geoMrgCtx);
      }
#if JVET_W0097_GPM_MMVD_TM
      maxNumMergeCandidates = min((int)maxNumMergeCandidates, m_geoMrgCtx.numValidMergeCand);
#endif
      const auto savedAffine = pu.cu->affine;
      const auto savedIMV = pu.cu->imv;
      for (int i = 0; i < maxNumMergeCandidates; i++)
      {
        if (m_mhpMrgTempBufSet // MC results already stored when checking GEO RD cost
#if INTER_LIC
          && (fakePredData.cu->LICFlag == m_geoMrgCtx.LICFlags[i])
#endif
          )
        {
          continue;
        }
        // get prediction for the additional hypothesis
        int refList = m_geoMrgCtx.interDirNeighbours[i] - 1; CHECK(refList != 0 && refList != 1, "");
        fakePredData.interDir = refList + 1;
        fakePredData.mv[refList] = m_geoMrgCtx.mvFieldNeighbours[(i << 1) + refList].mv;
        fakePredData.refIdx[refList] = m_geoMrgCtx.mvFieldNeighbours[(i << 1) + refList].refIdx;
        fakePredData.refIdx[1 - refList] = -1;
        fakePredData.cu->affine = false;
        fakePredData.cu->imv = m_geoMrgCtx.useAltHpelIf[i] ? IMV_HPEL : 0;
        fakePredData.mvRefine = true;
        motionCompensation(fakePredData, m_mhpMrgTempBuf[i], REF_PIC_LIST_X, true, false);
        fakePredData.mvRefine = false;
        // the restore of affine flag and imv flag has to be here
        fakePredData.cu->imv = savedIMV;
        fakePredData.cu->affine = savedAffine;
      }
      setGeoTmpBuffer();
    }
#if JVET_W0097_GPM_MMVD_TM
    else
    {
      maxNumMergeCandidates = min((int)maxNumMergeCandidates, m_geoMrgCtx.numValidMergeCand);
    }
#endif
    for (int i = 0; i < maxNumMergeCandidates; i++)
    {
      int refList = m_geoMrgCtx.interDirNeighbours[i] - 1; CHECK(refList != 0 && refList != 1, "");
      tempMHPredData.mrgIdx = i;
      tempMHPredData.isMrg = true;
      tempMHPredData.refIdx = m_geoMrgCtx.mvFieldNeighbours[(i << 1) + refList].refIdx;
      tempMHPredData.mv = m_geoMrgCtx.mvFieldNeighbours[(i << 1) + refList].mv;
      tempMHPredData.imv = m_geoMrgCtx.useAltHpelIf[i] ? IMV_HPEL : 0;
#if INTER_LIC 
      tempMHPredData.LICFlag = savedLICFlag;
#endif
      tempMHPredData.refList = refList;
      for (tempMHPredData.weightIdx = 0; tempMHPredData.weightIdx < numWeights; ++tempMHPredData.weightIdx)
      {
        tempOrigBuf.copyFrom(tempPredBuf, true);

        tempOrigBuf.addHypothesisAndClip(m_mhpMrgTempBuf[i], g_addHypWeight[tempMHPredData.weightIdx], pu.cs->slice->clpRngs(), true);
        Distortion uiSad = distParam.distFunc(distParam);
        uint32_t uiBits = x.bits + (i + 1);
        if (i == pu.cs->sps->getMaxNumGeoCand() - 1)
        {
          uiBits--;
        }
        uiBits += tempMHPredData.weightIdx + 1;
        if (tempMHPredData.weightIdx == numWeights - 1)
          uiBits--;
        Distortion uiCostTemp = uiSad + m_pcRdCost->getCost(uiBits);
        if (uiCostTemp < x.cost)
        {
          MEResult result;
          result.cu = *pu.cu;
          result.pu = pu;
          CHECK(tempMHPredData.mrgIdx >= maxNumMergeCandidates, "");
          result.pu.addHypData.push_back(tempMHPredData);
          result.cost = uiCostTemp;
          result.bits = uiBits;
          // store MHP MC result for next additonal hypothesis test
          if (pu.addHypData.size() < sps.getMaxNumAddHyps() && m_mhpTempBufCounter < GEO_MAX_TRY_WEIGHTED_SAD)
          {
            result.predBuf = &m_mhpTempBuf[m_mhpTempBufCounter];
            result.predBufIdx = m_mhpTempBufCounter;
            m_mhpTempBufCounter++;
            result.predBuf->copyFrom(tempOrigBuf, true);
          }
          out.push_back(result);
        }
      } // weightIdx
    } // i
  }
  tempMHPredData.isMrg = false;
#if INTER_LIC 
  tempMHPredData.LICFlag = pu.cu->LICFlag;
#endif
  tempMHPredData.imv = pu.cu->imv;
  for (tempMHPredData.weightIdx = 0; tempMHPredData.weightIdx < numWeights; ++tempMHPredData.weightIdx)
  {
    tempOrigBuf.copyFrom(origBuf, true);
    tempOrigBuf.removeHighFreq(tempPredBuf, m_pcEncCfg->getClipForBiPredMeEnabled(), pu.cu->slice->clpRngs(), g_addHypWeight[tempMHPredData.weightIdx]);
    for (tempMHPredData.refIdx = 0; tempMHPredData.refIdx < iNumMHRefPics; ++tempMHPredData.refIdx)
    {
      tempMHPredData.mvpIdx = 0;
      {
        const int iRefPicList = MHRefPics[tempMHPredData.refIdx].refList;
        const int iRefIdxPred = MHRefPics[tempMHPredData.refIdx].refIdx;
        const RefPicList eRefPicList = RefPicList(iRefPicList);
        uint32_t uiBits = x.bits + getAdditionalHypothesisInitialBits(tempMHPredData, numWeights, iNumMHRefPics);
        auto amvpInfo = PU::getMultiHypMVPCands(pu, tempMHPredData);
        Mv cMvPred = amvpInfo.mvCand[tempMHPredData.mvpIdx];
        if ((pu.addHypData.size() + 1 - pu.numMergedAddHyps) < sps.getMaxNumAddHyps())
          uiBits++;
        Mv cMv(0, 0);
        if (g_isReusedUniMVsFilled[idx1][idx2][idx3][idx4])
        {
          cMv = g_reusedUniMVs[idx1][idx2][idx3][idx4][iRefPicList][iRefIdxPred];
          uint32_t bitsDummy = 0;
          Distortion uiCostDummy = 0;
          xCheckBestMVP(eRefPicList, cMv, cMvPred, tempMHPredData.mvpIdx, amvpInfo, bitsDummy, uiCostDummy, pu.cu->imv);
        }
        else
        {
          cMv = cMvPred;
        }
        Distortion uiCostTemp = 0;
#if INTER_LIC
        pu.cu->LICFlag = tempMHPredData.LICFlag;
#endif
        xMotionEstimation(pu, tempOrigBuf, eRefPicList, cMvPred, iRefIdxPred, cMv, tempMHPredData.mvpIdx, uiBits, uiCostTemp, amvpInfo, false, g_addHypWeight[tempMHPredData.weightIdx]);
        xCheckBestMVP(eRefPicList, cMv, cMvPred, tempMHPredData.mvpIdx, amvpInfo, uiBits, uiCostTemp, pu.cu->imv);
#if INTER_LIC
        pu.cu->LICFlag = savedLICFlag;
#endif
        tempMHPredData.mv = cMv;
        tempMHPredData.mv.mvCliptoStorageBitDepth();

        tempMHPredData.mvd = cMv - cMvPred;

        if (uiCostTemp < x.cost)
        {
          MEResult result;
          result.cu = *pu.cu;
          result.pu = pu;
          result.pu.addHypData.push_back(tempMHPredData);
          result.cost = uiCostTemp;
          result.bits = uiBits;
          out.push_back(result);
        }
      }
    }
    }

  // buffer recycling
  if (m_pcEncCfg->getNumMHPCandsToTest() > 4 && x.predBufIdx >= 0 && m_mhpTempBufCounter > x.predBufIdx + 1)
  {
    if (x.predBufIdx < GEO_MAX_TRY_WEIGHTED_SAD - 1)
    {
      ::memcpy(m_mhpTempBuf + x.predBufIdx, m_mhpTempBuf + x.predBufIdx + 1, (m_mhpTempBufCounter - x.predBufIdx - 1) * sizeof(PelUnitBuf));
    }
    m_mhpTempBufCounter--;
  }
}

inline unsigned InterSearch::getAdditionalHypothesisInitialBits(const MultiHypPredictionData& mhData,
  const int iNumWeights,
  const int iNumMHRefPics)
{
  unsigned uiBits = 0;

  // weight idx
  uiBits += mhData.weightIdx + 1;
  if (mhData.weightIdx == iNumWeights - 1)
    uiBits--;

  // AMVP flag
  uiBits++;

  // ref idx
  uiBits += mhData.refIdx + 1;
  if (mhData.refIdx == iNumMHRefPics - 1)
  {
    uiBits--;
  }

  return uiBits;
}
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
void InterSearch::initGeoAngleSelection(PredictionUnit& pu
#if JVET_Y0065_GPM_INTRA
                                      , IntraPrediction* pcIntraPred, const uint8_t (&mpm)[GEO_NUM_PARTITION_MODE][2][GEO_MAX_NUM_INTRA_CANDS]
#endif
)
{
  xAMLGetCurBlkTemplate(pu, pu.lwidth(), pu.lheight());
  memset(&m_gpmacsSplitModeTmSelAvail[0][0][0], 0, sizeof(m_gpmacsSplitModeTmSelAvail));
  memset(&m_gpmPartTplCost[0][0][0][0], -1, sizeof(m_gpmPartTplCost));

  int16_t wIdx = floorLog2((uint32_t)pu.lwidth ()) - GEO_MIN_CU_LOG2;
  int16_t hIdx = floorLog2((uint32_t)pu.lheight()) - GEO_MIN_CU_LOG2;
  m_tplWeightTbl    = m_tplWeightTblDict   [hIdx][wIdx];
  m_tplColWeightTbl = m_tplColWeightTblDict[hIdx][wIdx];

#if JVET_Y0065_GPM_INTRA
  pcIntraPred->clearPrefilledIntraGPMRefTemplate();
  pcIntraPred->prefillIntraGPMReferenceSamples(pu, GEO_MODE_SEL_TM_SIZE, GEO_MODE_SEL_TM_SIZE);
  pcIntraPred->setPrefilledIntraGPMMPMModeAll(mpm);
#endif
}

void InterSearch::setGeoSplitModeToSyntaxTable(PredictionUnit& pu, MergeCtx& mergeCtx0, int mergeCand0, MergeCtx& mergeCtx1, int mergeCand1
#if JVET_Y0065_GPM_INTRA
                                             , IntraPrediction* pcIntraPred
#endif
                                             , int mmvdCand0, int mmvdCand1)
{
#if JVET_Y0065_GPM_INTRA
  bool isIntra[2];
  xRemapMrgIndexAndMmvdIdx(mergeCand0, mergeCand1, mmvdCand0, mmvdCand1, isIntra[0], isIntra[1]);
#endif
  const int idx0 = mmvdCand0 + 1;
  const int idx1 = mmvdCand1 + 1;

  if ((m_gpmacsSplitModeTmSelAvail[idx0][idx1][mergeCand0] & ((uint16_t)1 << mergeCand1)) == 0)
  {
    uint8_t numValidInList = 0;
    uint8_t modeList[GEO_NUM_SIG_PARTMODE];
    selectGeoSplitModes(pu
#if JVET_Y0065_GPM_INTRA
                      , pcIntraPred
#endif
                      , m_gpmPartTplCost[idx0][mergeCand0]
                      , m_gpmPartTplCost[idx1][mergeCand1]
                      , mergeCtx0
                      , mergeCand0
#if JVET_Y0065_GPM_INTRA
                      + (isIntra[0] ? GEO_MAX_NUM_UNI_CANDS : 0)
#endif
                      , mergeCtx1
                      , mergeCand1
#if JVET_Y0065_GPM_INTRA
                      + (isIntra[1] ? GEO_MAX_NUM_UNI_CANDS : 0)
#endif
                      , numValidInList
                      , modeList
#if JVET_W0097_GPM_MMVD_TM
                      , (mmvdCand0 >= GPM_EXT_MMVD_MAX_REFINE_NUM ? -1 : mmvdCand0)
                      , (mmvdCand1 >= GPM_EXT_MMVD_MAX_REFINE_NUM ? -1 : mmvdCand1)
#endif
    );

    xSetGpmModeToSyntaxModeTable(numValidInList, modeList, m_gpmacsSplitModeTmSel[idx0][idx1][mergeCand0][mergeCand1]);
    m_gpmacsSplitModeTmSelAvail[idx0][idx1][mergeCand0] |= ((uint16_t)1 << mergeCand1);
  }
}

#if JVET_W0097_GPM_MMVD_TM && TM_MRG
void InterSearch::setGeoTMSplitModeToSyntaxTable(PredictionUnit& pu, MergeCtx(&mergeCtx)[GEO_NUM_TM_MV_CAND], int mergeCand0, int mergeCand1, int mmvdCand0, int mmvdCand1)
{
  const int idx0 = mmvdCand0 + 1;
  const int idx1 = mmvdCand1 + 1;

  if ((m_gpmacsSplitModeTmSelAvail[idx0][idx1][mergeCand0] & ((uint16_t)1 << mergeCand1)) == 0)
  {
    uint8_t numValidInList = 0;
    uint8_t modeList[GEO_NUM_SIG_PARTMODE];
    selectGeoTMSplitModes(pu
                        , m_gpmPartTplCost[idx0][mergeCand0]
                        , m_gpmPartTplCost[idx1][mergeCand1]
                        , mergeCtx
                        , mergeCand0
                        , mergeCand1
                        , numValidInList
                        , modeList
    );

    xSetGpmModeToSyntaxModeTable(numValidInList, modeList, m_gpmacsSplitModeTmSel[idx0][idx1][mergeCand0][mergeCand1]);
    m_gpmacsSplitModeTmSelAvail[idx0][idx1][mergeCand0] |= ((uint16_t)1 << mergeCand1);
  }
}
#endif

int InterSearch::convertGeoSplitModeToSyntax(int splitDir, int mergeCand0, int mergeCand1, int mmvdCand0, int mmvdCand1)
{
#if JVET_Y0065_GPM_INTRA
  bool isIntra[2];
  xRemapMrgIndexAndMmvdIdx(mergeCand0, mergeCand1, mmvdCand0, mmvdCand1, isIntra[0], isIntra[1]);
#endif
  return m_gpmacsSplitModeTmSel[mmvdCand0 + 1][mmvdCand1 + 1][mergeCand0][mergeCand1][splitDir];
}

bool InterSearch::selectGeoSplitModes(PredictionUnit &pu,
#if JVET_Y0065_GPM_INTRA
                                      IntraPrediction* pcIntraPred,
#endif
                                      uint32_t (&gpmTplCostPart0)[2][GEO_NUM_PARTITION_MODE],
                                      uint32_t (&gpmTplCostPart1)[2][GEO_NUM_PARTITION_MODE],
                                      MergeCtx& mergeCtx0, int mergeCand0, MergeCtx& mergeCtx1, int mergeCand1, uint8_t& numValidInList, uint8_t (&modeList)[GEO_NUM_SIG_PARTMODE], int mmvdCand0, int mmvdCand1)
{
  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    getBestGeoModeList(pu, numValidInList, modeList, nullptr, nullptr, nullptr, nullptr);
    return false;
  }

  if (PU::checkRprRefExistingInGpm(pu, mergeCtx0, mergeCand0, mergeCtx1, mergeCand1))
  {
    bool backupTplValid[2] = {m_bAMLTemplateAvailabe[0], m_bAMLTemplateAvailabe[1]};
    m_bAMLTemplateAvailabe[0] = false;
    m_bAMLTemplateAvailabe[1] = false;
    getBestGeoModeList(pu, numValidInList, modeList, nullptr, nullptr, nullptr, nullptr);
    m_bAMLTemplateAvailabe[0] = backupTplValid[0];
    m_bAMLTemplateAvailabe[1] = backupTplValid[1];
    return false;
  }

  bool fillRefTplPart0 = gpmTplCostPart0[0][0] == std::numeric_limits<uint32_t>::max();
  bool fillRefTplPart1 = gpmTplCostPart1[1][0] == std::numeric_limits<uint32_t>::max();
  Pel* pRefTopPart0    = m_acYuvRefAMLTemplatePart0[0];
  Pel* pRefLeftPart0   = m_acYuvRefAMLTemplatePart0[1];
  Pel* pRefTopPart1    = m_acYuvRefAMLTemplatePart1[0];
  Pel* pRefLeftPart1   = m_acYuvRefAMLTemplatePart1[1];

  // First partition
  if (fillRefTplPart0)
  {
    fillPartGPMRefTemplate<0, false>(pu, mergeCtx0, mergeCand0, mmvdCand0, pRefTopPart0, pRefLeftPart0);
#if JVET_Y0065_GPM_INTRA
    xCollectIntraGeoPartCost<0>(pu, pcIntraPred, mergeCand0, gpmTplCostPart0[0]);
#endif
  }

  // Second
  if (fillRefTplPart1)
  {
    fillPartGPMRefTemplate<1, false>(pu, mergeCtx1, mergeCand1, mmvdCand1, pRefTopPart1, pRefLeftPart1);
#if JVET_Y0065_GPM_INTRA
    xCollectIntraGeoPartCost<1>(pu, pcIntraPred, mergeCand1, gpmTplCostPart1[1]);
#endif
  }

  // Get mode lists
  getBestGeoModeListEncoder(pu, numValidInList, modeList, pRefTopPart0, pRefLeftPart0, pRefTopPart1, pRefLeftPart1, gpmTplCostPart0, gpmTplCostPart1);
  return true;
}

#if JVET_W0097_GPM_MMVD_TM && TM_MRG
bool InterSearch::selectGeoTMSplitModes (PredictionUnit &pu, 
                                         uint32_t (&gpmTplCostPart0)[2][GEO_NUM_PARTITION_MODE],
                                         uint32_t (&gpmTplCostPart1)[2][GEO_NUM_PARTITION_MODE],
                                         MergeCtx (&mergeCtx)[GEO_NUM_TM_MV_CAND], int mergeCand0, int mergeCand1, uint8_t& numValidInList, uint8_t (&modeList)[GEO_NUM_SIG_PARTMODE])
{
  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    getBestGeoModeList(pu, numValidInList, modeList, nullptr, nullptr, nullptr, nullptr);
    return false;
  }

  if (PU::checkRprRefExistingInGpm(pu, mergeCtx[GEO_TM_OFF], mergeCand0, mergeCtx[GEO_TM_OFF], mergeCand1))
  {
    bool backupTplValid[2] = {m_bAMLTemplateAvailabe[0], m_bAMLTemplateAvailabe[1]};
    m_bAMLTemplateAvailabe[0] = false;
    m_bAMLTemplateAvailabe[1] = false;
    getBestGeoModeList(pu, numValidInList, modeList, nullptr, nullptr, nullptr, nullptr);
    m_bAMLTemplateAvailabe[0] = backupTplValid[0];
    m_bAMLTemplateAvailabe[1] = backupTplValid[1];
    return false;
  }

  bool fillRefTplPart0  = gpmTplCostPart0[0][0] == std::numeric_limits<uint32_t>::max();
  bool fillRefTplPart1  = gpmTplCostPart1[1][0] == std::numeric_limits<uint32_t>::max();
  Pel* pRefTopPart0 [GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart0[0], m_acYuvRefAMLTemplatePart0[2], nullptr                      }; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_A]
  Pel* pRefLeftPart0[GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart0[1], m_acYuvRefAMLTemplatePart0[3], nullptr                      }; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_A]
  Pel* pRefTopPart1 [GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart1[0], nullptr,                       m_acYuvRefAMLTemplatePart1[2]}; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_L]
  Pel* pRefLeftPart1[GEO_NUM_TM_MV_CAND] = {nullptr, m_acYuvRefAMLTemplatePart1[1], nullptr,                       m_acYuvRefAMLTemplatePart1[3]}; // For mergeCtx[GEO_TM_SHAPE_AL] and mergeCtx[GEO_TM_SHAPE_L]

  // First partition
  if (fillRefTplPart0)
  {
    fillPartGPMRefTemplate<0, false>(pu, mergeCtx[GEO_TM_SHAPE_AL], mergeCand0, -1, pRefTopPart0[GEO_TM_SHAPE_AL], pRefLeftPart0[GEO_TM_SHAPE_AL]);
    fillPartGPMRefTemplate<0, false>(pu, mergeCtx[GEO_TM_SHAPE_A ], mergeCand0, -1, pRefTopPart0[GEO_TM_SHAPE_A ], pRefLeftPart0[GEO_TM_SHAPE_A ]);
  }

  // Second
  if (fillRefTplPart1)
  {
    fillPartGPMRefTemplate<1, false>(pu, mergeCtx[GEO_TM_SHAPE_AL], mergeCand1, -1, pRefTopPart1[GEO_TM_SHAPE_AL], pRefLeftPart1[GEO_TM_SHAPE_AL]);
    fillPartGPMRefTemplate<1, false>(pu, mergeCtx[GEO_TM_SHAPE_L ], mergeCand1, -1, pRefTopPart1[GEO_TM_SHAPE_L ], pRefLeftPart1[GEO_TM_SHAPE_L ]);
  }

  // Get mode lists
  getBestGeoTMModeListEncoder(pu, numValidInList, modeList, pRefTopPart0, pRefLeftPart0, pRefTopPart1, pRefLeftPart1, gpmTplCostPart0, gpmTplCostPart1);
  return true;
}
#endif

void InterSearch::getBestGeoModeListEncoder(PredictionUnit &pu, uint8_t& numValidInList,
                                            uint8_t(&modeList)[GEO_NUM_SIG_PARTMODE],
                                            Pel* pRefTopPart0, Pel* pRefLeftPart0,
                                            Pel* pRefTopPart1, Pel* pRefLeftPart1,
                                            uint32_t(&gpmTplCostPart0)[2][GEO_NUM_PARTITION_MODE],
                                            uint32_t(&gpmTplCostPart1)[2][GEO_NUM_PARTITION_MODE])
{
  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    getBestGeoModeList(pu, numValidInList, modeList, nullptr, nullptr, nullptr, nullptr);
    return;
  }

  // Check mode
  bool filledRefTplPart0 = gpmTplCostPart0[0][0] == std::numeric_limits<uint32_t>::max();
  bool filledRefTplPart1 = gpmTplCostPart1[1][0] == std::numeric_limits<uint32_t>::max();
  int bitDepth = pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA);

  if (m_bAMLTemplateAvailabe[0])
  {
    SizeType   szPerLine            = pu.lwidth();
    PelUnitBuf pcBufPredCurTop      = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTopPart0 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart0,                szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTopPart1 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart1,                szPerLine, GEO_MODE_SEL_TM_SIZE));

    const int maskStride2[3] = { -(int)szPerLine, (int)szPerLine, -(int)szPerLine }; // template length
    const int maskStride[3] = { GEO_WEIGHT_MASK_SIZE_EXT, GEO_WEIGHT_MASK_SIZE_EXT, -GEO_WEIGHT_MASK_SIZE_EXT }; // mask stride
    const int stepX[3] = { 1, -1, 1 };

    // Cost of partition 0
    if(filledRefTplPart0)
    {
      GetAbsDiffPerSample(pcBufPredRefTopPart0.Y(), pcBufPredCurTop.Y(), pcBufPredRefTopPart0.Y());
      uint32_t fullCostPart0 = (uint32_t)GetSampleSum(pcBufPredRefTopPart0.Y(), bitDepth);

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 0>(splitDir);
        uint32_t tempDist = (uint32_t)Get01MaskedSampleSum(pcBufPredRefTopPart0.Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart0[0][splitDir] = tempDist;
        gpmTplCostPart0[1][splitDir] = fullCostPart0 - tempDist; // pre-calculated
      }
    }

    // Cost of partition 1
    if(filledRefTplPart1)
    {
      GetAbsDiffPerSample(pcBufPredRefTopPart1.Y(), pcBufPredCurTop.Y(), pcBufPredRefTopPart1.Y());
      uint32_t fullCostPart1 = (uint32_t)GetSampleSum(pcBufPredRefTopPart1.Y(), bitDepth);

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 0>(splitDir);
        uint32_t tempDist = (uint32_t)Get01MaskedSampleSum(pcBufPredRefTopPart1.Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart1[0][splitDir] = tempDist;  // pre-calculated
        gpmTplCostPart1[1][splitDir] = fullCostPart1 - tempDist;
      }
    }
  }
  else
  {
    if (filledRefTplPart0)
    {
      memset(gpmTplCostPart0[0], 0, sizeof(uint32_t) * GEO_NUM_PARTITION_MODE);
      memset(gpmTplCostPart0[1], 0, sizeof(uint32_t) * GEO_NUM_PARTITION_MODE);
    }
    if (filledRefTplPart1)
    {
      memset(gpmTplCostPart1[1], 0, sizeof(uint32_t) * GEO_NUM_PARTITION_MODE);
      memset(gpmTplCostPart1[0], 0, sizeof(uint32_t) * GEO_NUM_PARTITION_MODE);
    }
  }

  if (m_bAMLTemplateAvailabe[1])
  {
    SizeType   szPerLine             = pu.lheight();
    PelUnitBuf pcBufPredCurLeft      = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufPredRefLeftPart0 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart0,               szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufPredRefLeftPart1 = PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart1,               szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation

    const int maskStride2[3] = { -(int)szPerLine, -(int)szPerLine, -(int)szPerLine }; // template length
    const int maskStride[3] = { (int)szPerLine, (int)szPerLine, (int)szPerLine }; // mask stride
    const int stepX[3] = { 1, 1, 1 };

    // Cost of partition 0
    if (filledRefTplPart0)
    {
      GetAbsDiffPerSample(pcBufPredRefLeftPart0.Y(), pcBufPredCurLeft.Y(), pcBufPredRefLeftPart0.Y());
      uint32_t fullCostPart0 = (uint32_t)GetSampleSum(pcBufPredRefLeftPart0.Y(), bitDepth);

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 2>(splitDir);
        uint32_t tempDist = (uint32_t)Get01MaskedSampleSum(pcBufPredRefLeftPart0.Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart0[0][splitDir] += tempDist;
        gpmTplCostPart0[1][splitDir] += fullCostPart0 - tempDist; // pre-calculated
      }
    }

    // Cost of partition 1
    if (filledRefTplPart1)
    {
      GetAbsDiffPerSample(pcBufPredRefLeftPart1.Y(), pcBufPredCurLeft.Y(), pcBufPredRefLeftPart1.Y());
      uint32_t fullCostPart1 = (uint32_t)GetSampleSum(pcBufPredRefLeftPart1.Y(), bitDepth);

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 2>(splitDir);
        uint32_t tempDist = (uint32_t)Get01MaskedSampleSum(pcBufPredRefLeftPart1.Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart1[0][splitDir] += tempDist;  // pre-calculated
        gpmTplCostPart1[1][splitDir] += fullCostPart1 - tempDist;
      }
    }
  }

  // Check split mode cost
  uint32_t uiCost[GEO_NUM_PARTITION_MODE];
  for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
  {
    uiCost[splitDir] = gpmTplCostPart0[0][splitDir] + gpmTplCostPart1[1][splitDir];
  }

  // Find best N candidates
  numValidInList = (uint8_t)getIndexMappingTableToSortedArray1D<uint32_t, GEO_NUM_PARTITION_MODE, uint8_t, GEO_NUM_SIG_PARTMODE>(uiCost, modeList);

}

#if JVET_W0097_GPM_MMVD_TM && TM_MRG
void InterSearch::getBestGeoTMModeListEncoder(PredictionUnit &pu, uint8_t& numValidInList,
                                              uint8_t(&modeList)[GEO_NUM_SIG_PARTMODE],
                                              Pel* (&pRefTopPart0)[GEO_NUM_TM_MV_CAND], Pel* (&pRefLeftPart0)[GEO_NUM_TM_MV_CAND],
                                              Pel* (&pRefTopPart1)[GEO_NUM_TM_MV_CAND], Pel* (&pRefLeftPart1)[GEO_NUM_TM_MV_CAND],
                                              uint32_t(&gpmTplCostPart0)[2][GEO_NUM_PARTITION_MODE],
                                              uint32_t(&gpmTplCostPart1)[2][GEO_NUM_PARTITION_MODE])
{
  if (!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1])
  {
    getBestGeoModeList(pu, numValidInList, modeList, nullptr, nullptr, nullptr, nullptr);
    return;
  }

  // Check mode
  bool filledRefTplPart0 = gpmTplCostPart0[0][0] == std::numeric_limits<uint32_t>::max();
  bool filledRefTplPart1 = gpmTplCostPart1[1][0] == std::numeric_limits<uint32_t>::max();
  int bitDepth = pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA);

  if (m_bAMLTemplateAvailabe[0])
  {
    SizeType   szPerLine       = pu.lwidth();
    PelUnitBuf pcBufPredCurTop = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTopPart0[GEO_NUM_TM_MV_CAND] = {PelUnitBuf(), 
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart0[GEO_TM_SHAPE_AL], szPerLine, GEO_MODE_SEL_TM_SIZE)),
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart0[GEO_TM_SHAPE_A ], szPerLine, GEO_MODE_SEL_TM_SIZE)),
                                                           PelUnitBuf()};
    PelUnitBuf pcBufPredRefTopPart1[GEO_NUM_TM_MV_CAND] = {PelUnitBuf(),
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart1[GEO_TM_SHAPE_AL], szPerLine, GEO_MODE_SEL_TM_SIZE)),
                                                           PelUnitBuf(),
                                                           PelUnitBuf(pu.chromaFormat, PelBuf(pRefTopPart1[GEO_TM_SHAPE_L ], szPerLine, GEO_MODE_SEL_TM_SIZE))};

    const int maskStride2[3] = { -(int)szPerLine, (int)szPerLine, -(int)szPerLine }; // template length
    const int maskStride[3] = { GEO_WEIGHT_MASK_SIZE_EXT, GEO_WEIGHT_MASK_SIZE_EXT, -GEO_WEIGHT_MASK_SIZE_EXT }; // mask stride
    const int stepX[3] = { 1, -1, 1 };

    // Cost of partition 0
    if(filledRefTplPart0)
    {
      GetAbsDiffPerSample(pcBufPredRefTopPart0[GEO_TM_SHAPE_AL].Y(), pcBufPredCurTop.Y(), pcBufPredRefTopPart0[GEO_TM_SHAPE_AL].Y());
      GetAbsDiffPerSample(pcBufPredRefTopPart0[GEO_TM_SHAPE_A ].Y(), pcBufPredCurTop.Y(), pcBufPredRefTopPart0[GEO_TM_SHAPE_A ].Y());

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        uint8_t shapeIdx  = g_geoTmShape[0][g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 0>(splitDir);
        uint32_t tempDist = (uint32_t)Get01MaskedSampleSum(pcBufPredRefTopPart0[shapeIdx].Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart0[0][splitDir] = tempDist;
      }
    }

    // Cost of partition 1
    if(filledRefTplPart1)
    {
      GetAbsDiffPerSample(pcBufPredRefTopPart1[GEO_TM_SHAPE_AL].Y(), pcBufPredCurTop.Y(), pcBufPredRefTopPart1[GEO_TM_SHAPE_AL].Y());
      GetAbsDiffPerSample(pcBufPredRefTopPart1[GEO_TM_SHAPE_L ].Y(), pcBufPredCurTop.Y(), pcBufPredRefTopPart1[GEO_TM_SHAPE_L ].Y());

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        uint8_t shapeIdx  = g_geoTmShape[1][g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 0>(splitDir);
        uint32_t tempDist = (uint32_t)Get01InvMaskedSampleSum(pcBufPredRefTopPart1[shapeIdx].Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart1[1][splitDir] = tempDist;
      }
    }
  }
  else
  {
    if (filledRefTplPart0)
    {
      memset(gpmTplCostPart0[0], 0, sizeof(uint32_t) * GEO_NUM_PARTITION_MODE);
    }
    if (filledRefTplPart1)
    {
      memset(gpmTplCostPart1[1], 0, sizeof(uint32_t) * GEO_NUM_PARTITION_MODE);
    }
  }

  if (m_bAMLTemplateAvailabe[1])
  {
    SizeType   szPerLine        = pu.lheight();
    PelUnitBuf pcBufPredCurLeft = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], szPerLine, GEO_MODE_SEL_TM_SIZE)); // reordered to make it 1 row to enable SIMD
    PelUnitBuf pcBufPredRefLeftPart0[GEO_NUM_TM_MV_CAND] = {PelUnitBuf(),
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart0[GEO_TM_SHAPE_AL], szPerLine, GEO_MODE_SEL_TM_SIZE)), // To enable SIMD for cost computation
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart0[GEO_TM_SHAPE_A ], szPerLine, GEO_MODE_SEL_TM_SIZE)), // To enable SIMD for cost computation
                                                            PelUnitBuf()};
    PelUnitBuf pcBufPredRefLeftPart1[GEO_NUM_TM_MV_CAND] = { PelUnitBuf(),
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart1[GEO_TM_SHAPE_AL], szPerLine, GEO_MODE_SEL_TM_SIZE)), // To enable SIMD for cost computation
                                                            PelUnitBuf(),
                                                            PelUnitBuf(pu.chromaFormat, PelBuf(pRefLeftPart1[GEO_TM_SHAPE_L ], szPerLine, GEO_MODE_SEL_TM_SIZE))}; // To enable SIMD for cost computation

    const int maskStride2[3] = { -(int)szPerLine, -(int)szPerLine, -(int)szPerLine }; // template length
    const int maskStride[3] = { (int)szPerLine, (int)szPerLine, (int)szPerLine }; // mask stride
    const int stepX[3] = { 1, 1, 1 };

    // Cost of partition 0
    if (filledRefTplPart0)
    {
      GetAbsDiffPerSample(pcBufPredRefLeftPart0[GEO_TM_SHAPE_AL].Y(), pcBufPredCurLeft.Y(), pcBufPredRefLeftPart0[GEO_TM_SHAPE_AL].Y());
      GetAbsDiffPerSample(pcBufPredRefLeftPart0[GEO_TM_SHAPE_A ].Y(), pcBufPredCurLeft.Y(), pcBufPredRefLeftPart0[GEO_TM_SHAPE_A ].Y());

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        uint8_t shapeIdx  = g_geoTmShape[0][g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 2>(splitDir);
        uint32_t tempDist = (uint32_t)Get01MaskedSampleSum(pcBufPredRefLeftPart0[shapeIdx].Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart0[0][splitDir] += tempDist;
      }
    }

    // Cost of partition 1
    if (filledRefTplPart1)
    {
      GetAbsDiffPerSample(pcBufPredRefLeftPart1[GEO_TM_SHAPE_AL].Y(), pcBufPredCurLeft.Y(), pcBufPredRefLeftPart1[GEO_TM_SHAPE_AL].Y());
      GetAbsDiffPerSample(pcBufPredRefLeftPart1[GEO_TM_SHAPE_L ].Y(), pcBufPredCurLeft.Y(), pcBufPredRefLeftPart1[GEO_TM_SHAPE_L ].Y());

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
      {
        int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
        uint8_t shapeIdx  = g_geoTmShape[1][g_GeoParams[splitDir][0]];
        Pel* mask = getTplWeightTableCU<false, 2>(splitDir);
        uint32_t tempDist = (uint32_t)Get01InvMaskedSampleSum(pcBufPredRefLeftPart1[shapeIdx].Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        gpmTplCostPart1[1][splitDir] += tempDist;
      }
    }
  }

  // Check split mode cost
  uint32_t uiCost[GEO_NUM_PARTITION_MODE];
  for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
  {
    uiCost[splitDir] = gpmTplCostPart0[0][splitDir] + gpmTplCostPart1[1][splitDir];
  }

  // Find best N candidates
  numValidInList = (uint8_t)getIndexMappingTableToSortedArray1D<uint32_t, GEO_NUM_PARTITION_MODE, uint8_t, GEO_NUM_SIG_PARTMODE>(uiCost, modeList);

}
#endif

#if JVET_Y0065_GPM_INTRA
template <uint8_t partIdx>
void InterSearch::xCollectIntraGeoPartCost(PredictionUnit &pu, IntraPrediction* pcIntraPred, int mergeCand, uint32_t(&gpmTplCost)[GEO_NUM_PARTITION_MODE])
{
  if ((!m_bAMLTemplateAvailabe[0] && !m_bAMLTemplateAvailabe[1]) || gpmTplCost[0] != std::numeric_limits<uint32_t>::max() || mergeCand < GEO_MAX_NUM_UNI_CANDS)
  {
    return;
  }

  std::vector<Pel>* LUT = m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && m_pcReshape->getCTUFlag() ? &m_pcReshape->getInvLUT() : nullptr;
  pcIntraPred->fillIntraGPMRefTemplateAll(pu, m_bAMLTemplateAvailabe[0], m_bAMLTemplateAvailabe[1], true, false, false, LUT, (partIdx == 0 ? mergeCand : 0), (partIdx == 1 ? mergeCand : 0));

  int  realCandIdx = mergeCand - GEO_MAX_NUM_UNI_CANDS;
  int  bitDepth    = pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA);
  Pel* pDiffTop    = partIdx == 0 ? m_acYuvRefAMLTemplatePart0[0] : m_acYuvRefAMLTemplatePart1[0];
  Pel* pDiffLeft   = partIdx == 0 ? m_acYuvRefAMLTemplatePart0[1] : m_acYuvRefAMLTemplatePart1[1];

  static_vector<int, GEO_NUM_PARTITION_MODE> intraModeToSplitDirAll[NUM_INTRA_MODE];
  for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; ++splitDir)
  {
    uint8_t intraMode = pcIntraPred->getPrefilledIntraGPMMPMMode(partIdx, splitDir, realCandIdx);
    intraModeToSplitDirAll[intraMode].push_back(splitDir);
  }

  if (m_bAMLTemplateAvailabe[0])
  {
    SizeType   szPerLine        = pu.lwidth();
    PelUnitBuf pcBufPredCurTop  = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[0][0], szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufPredRefTop  = PelUnitBuf(pu.chromaFormat, PelBuf(nullptr,                     szPerLine, GEO_MODE_SEL_TM_SIZE));
    PelUnitBuf pcBufDiffTop     = PelUnitBuf(pu.chromaFormat, PelBuf(pDiffTop,                    szPerLine, GEO_MODE_SEL_TM_SIZE));

    const int maskStride2[3] = { -(int)szPerLine, (int)szPerLine, -(int)szPerLine }; // template length
    const int maskStride[3] = { GEO_WEIGHT_MASK_SIZE_EXT, GEO_WEIGHT_MASK_SIZE_EXT, -GEO_WEIGHT_MASK_SIZE_EXT }; // mask stride
    const int stepX[3] = { 1, -1, 1 };

    for (uint8_t intraMode = 0; intraMode < NUM_INTRA_MODE; ++intraMode)
    {
      static_vector<int, GEO_NUM_PARTITION_MODE>& toSplitDir = intraModeToSplitDirAll[intraMode];
      if (toSplitDir.size() > 0)
      {
        pcBufPredRefTop.Y().buf = pcIntraPred->getPrefilledIntraGPMRefTemplate(intraMode, 0);
        GetAbsDiffPerSample(pcBufDiffTop.Y(), pcBufPredCurTop.Y(), pcBufPredRefTop.Y());

        for (int i = 0; i < toSplitDir.size(); ++i)
        {
          int splitDir = toSplitDir[i];
          int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
          Pel* mask = getTplWeightTableCU<false, 0>(splitDir);
          gpmTplCost[splitDir] = (uint32_t)GetSampleSumFunc(partIdx + 2, pcBufDiffTop.Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        }
      }
    }
  }
  else
  {
    memset(gpmTplCost, 0, sizeof(gpmTplCost));
  }

  if (m_bAMLTemplateAvailabe[1])
  {
    SizeType   szPerLine        = pu.lheight();
    PelUnitBuf pcBufPredCurLeft = PelUnitBuf(pu.chromaFormat, PelBuf(m_acYuvCurAMLTemplate[1][0], szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufPredRefLeft = PelUnitBuf(pu.chromaFormat, PelBuf(nullptr,                     szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation
    PelUnitBuf pcBufDiffLeft    = PelUnitBuf(pu.chromaFormat, PelBuf(pDiffLeft,                   szPerLine, GEO_MODE_SEL_TM_SIZE)); // To enable SIMD for cost computation

    const int maskStride2[3] = { -(int)szPerLine, -(int)szPerLine, -(int)szPerLine }; // template length
    const int maskStride[3] = { (int)szPerLine, (int)szPerLine, (int)szPerLine }; // mask stride
    const int stepX[3] = { 1, 1, 1 };

    for (uint8_t intraMode = 0; intraMode < NUM_INTRA_MODE; ++intraMode)
    {
      static_vector<int, GEO_NUM_PARTITION_MODE>& toSplitDir = intraModeToSplitDirAll[intraMode];
      if (toSplitDir.size() > 0)
      {
        pcBufPredRefLeft.Y().buf = pcIntraPred->getPrefilledIntraGPMRefTemplate(intraMode, 1);
        GetAbsDiffPerSample(pcBufDiffLeft.Y(), pcBufPredCurLeft.Y(), pcBufPredRefLeft.Y());

        for (int i = 0; i < toSplitDir.size(); ++i)
        {
          int splitDir = toSplitDir[i];
          int16_t mirrorIdx = g_angle2mirror[g_GeoParams[splitDir][0]];
          Pel* mask = getTplWeightTableCU<false, 2>(splitDir);
          gpmTplCost[splitDir] += (uint32_t)GetSampleSumFunc(partIdx + 2, pcBufDiffLeft.Y(), bitDepth, mask, stepX[mirrorIdx], maskStride[mirrorIdx], maskStride2[mirrorIdx]);
        }
      }
    }
  }
}
#endif
#endif

// AMVP
#if JVET_X0083_BM_AMVP_MERGE_MODE
void InterSearch::xEstimateMvPredAMVP( PredictionUnit& pu, PelUnitBuf& origBuf, RefPicList eRefPicList, int iRefIdx, Mv& rcMvPred, AMVPInfo& rAMVPInfo, bool bFilled, Distortion* puiDistBiP, MvField* mvFieldAmListCommon )
#else
void InterSearch::xEstimateMvPredAMVP( PredictionUnit& pu, PelUnitBuf& origBuf, RefPicList eRefPicList, int iRefIdx, Mv& rcMvPred, AMVPInfo& rAMVPInfo, bool bFilled, Distortion* puiDistBiP )
#endif
{
  Mv         cBestMv;
  int        iBestIdx   = 0;
  Distortion uiBestCost = std::numeric_limits<Distortion>::max();
  int        i;

  AMVPInfo*  pcAMVPInfo = &rAMVPInfo;

  // Fill the MV Candidates
  if (!bFilled)
  {
#if JVET_X0083_BM_AMVP_MERGE_MODE
    if (pu.amvpMergeModeFlag[1 - eRefPicList] == true)
    {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      const int mvFieldAmvpIdx0 = MAX_NUM_AMVP_CANDS_MAX_REF + iRefIdx * AMVP_MAX_NUM_CANDS_MEM;
      CHECK(mvFieldAmListCommon[mvFieldAmvpIdx0].refIdx != iRefIdx, "this is not possible");
#else
      const int mvFieldAmvpIdx0 = MAX_NUM_AMVP_CANDS_MAX_REF + iRefIdx * AMVP_MAX_NUM_CANDS;
#endif
      pcAMVPInfo->mvCand[0] = mvFieldAmListCommon[mvFieldAmvpIdx0].mv;
      pcAMVPInfo->numCand = 1;
#if !TM_AMVP || JVET_Y0128_NON_CTC || JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      const int mvFieldAmvpIdx1 = mvFieldAmvpIdx0 + 1;
      if (mvFieldAmListCommon[mvFieldAmvpIdx1].refIdx >= 0)
      {
        pcAMVPInfo->mvCand[1] = mvFieldAmListCommon[mvFieldAmvpIdx1].mv;
        pcAMVPInfo->numCand = 2;
      }
#endif
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      const int mvFieldAmvpIdx2 = mvFieldAmvpIdx0 + 2;
      if (mvFieldAmListCommon[mvFieldAmvpIdx2].refIdx >= 0)
      {
        pcAMVPInfo->mvCand[2] = mvFieldAmListCommon[mvFieldAmvpIdx2].mv;
        pcAMVPInfo->numCand = 3;
      }
#endif
      return;
    }
#endif
    PU::fillMvpCand( pu, eRefPicList, iRefIdx, *pcAMVPInfo
#if TM_AMVP
                   , this
#endif
    );
  }
#if INTER_LIC && RPR_ENABLE
  // xPredInterBlk may call PU::checkRprLicCondition()
#if JVET_Y0128_NON_CTC
  pu.interDir = (uint8_t)(eRefPicList + 1);
#endif
  pu.refIdx[eRefPicList]      = iRefIdx;
  pu.refIdx[1 - eRefPicList]  = NOT_VALID;
#endif

  // initialize Mvp index & Mvp
  iBestIdx = 0;
  cBestMv  = pcAMVPInfo->mvCand[0];

  PelUnitBuf predBuf = m_tmpStorageLCU.getBuf( UnitAreaRelative(*pu.cu, pu) );

  //-- Check Minimum Cost.
  for( i = 0 ; i < pcAMVPInfo->numCand; i++)
  {
#if TM_AMVP
    Distortion uiTmpCost = xGetTemplateCost( pu, origBuf, predBuf, pcAMVPInfo->mvCand[i], i, pcAMVPInfo->numCand, eRefPicList, iRefIdx );
#else
    Distortion uiTmpCost = xGetTemplateCost( pu, origBuf, predBuf, pcAMVPInfo->mvCand[i], i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx );
#endif
    if( uiBestCost > uiTmpCost )
    {
      uiBestCost     = uiTmpCost;
      cBestMv        = pcAMVPInfo->mvCand[i];
      iBestIdx       = i;
      (*puiDistBiP)  = uiTmpCost;
    }
  }

  // Setting Best MVP
  rcMvPred = cBestMv;
  pu.mvpIdx[eRefPicList] = iBestIdx;
  pu.mvpNum[eRefPicList] = pcAMVPInfo->numCand;

  return;
}

uint32_t InterSearch::xGetMvpIdxBits(int iIdx, int iNum)
{
  CHECK(iIdx < 0 || iNum < 0 || iIdx >= iNum, "Invalid parameters");

  if (iNum == 1)
  {
    return 0;
  }

  uint32_t uiLength = 1;
  int iTemp = iIdx;
  if ( iTemp == 0 )
  {
    return uiLength;
  }

  bool bCodeLast = ( iNum-1 > iTemp );

  uiLength += (iTemp-1);

  if( bCodeLast )
  {
    uiLength++;
  }

  return uiLength;
}

void InterSearch::xGetBlkBits( bool bPSlice, int iPartIdx, uint32_t uiLastMode, uint32_t uiBlkBit[3])
{
  uiBlkBit[0] = (! bPSlice) ? 3 : 1;
  uiBlkBit[1] = 3;
  uiBlkBit[2] = 5;
}

void InterSearch::xCopyAMVPInfo (AMVPInfo* pSrc, AMVPInfo* pDst)
{
  pDst->numCand = pSrc->numCand;
  for (int i = 0; i < pSrc->numCand; i++)
  {
    pDst->mvCand[i] = pSrc->mvCand[i];
  }
#if TM_AMVP
  pDst->maxSimilarityThreshold = pSrc->maxSimilarityThreshold;
#endif
}

void InterSearch::xCheckBestMVP ( RefPicList eRefPicList, Mv cMv, Mv& rcMvPred, int& riMVPIdx, AMVPInfo& amvpInfo, uint32_t& ruiBits, Distortion& ruiCost, const uint8_t imv )
{
  if ( imv > 0 && imv < 3 )
  {
    return;
  }

  AMVPInfo* pcAMVPInfo = &amvpInfo;

  CHECK(pcAMVPInfo->mvCand[riMVPIdx] != rcMvPred, "Invalid MV prediction candidate");

  if (pcAMVPInfo->numCand < 2)
  {
    return;
  }

  m_pcRdCost->setCostScale ( 0    );

  int iBestMVPIdx = riMVPIdx;

  Mv pred = rcMvPred;
  pred.changeTransPrecInternal2Amvr(imv);
  m_pcRdCost->setPredictor( pred );
  Mv mv = cMv;
  mv.changeTransPrecInternal2Amvr(imv);
  int iOrgMvBits = m_pcRdCost->getBitsOfVectorWithPredictor(mv.getHor(), mv.getVer(), 0);
#if TM_AMVP
  iOrgMvBits += m_auiMVPIdxCost[riMVPIdx][pcAMVPInfo->numCand];
#else
  iOrgMvBits += m_auiMVPIdxCost[riMVPIdx][AMVP_MAX_NUM_CANDS];
#endif
  int iBestMvBits = iOrgMvBits;

  for (int iMVPIdx = 0; iMVPIdx < pcAMVPInfo->numCand; iMVPIdx++)
  {
    if (iMVPIdx == riMVPIdx)
    {
      continue;
    }

    pred = pcAMVPInfo->mvCand[iMVPIdx];
    pred.changeTransPrecInternal2Amvr(imv);
    m_pcRdCost->setPredictor( pred );
    int iMvBits = m_pcRdCost->getBitsOfVectorWithPredictor(mv.getHor(), mv.getVer(), 0);
#if TM_AMVP
    iMvBits += m_auiMVPIdxCost[iMVPIdx][pcAMVPInfo->numCand];
#else
    iMvBits += m_auiMVPIdxCost[iMVPIdx][AMVP_MAX_NUM_CANDS];
#endif

    if (iMvBits < iBestMvBits)
    {
      iBestMvBits = iMvBits;
      iBestMVPIdx = iMVPIdx;
    }
  }

  if (iBestMVPIdx != riMVPIdx)  //if changed
  {
    rcMvPred = pcAMVPInfo->mvCand[iBestMVPIdx];

    riMVPIdx = iBestMVPIdx;
    uint32_t uiOrgBits = ruiBits;
    ruiBits = uiOrgBits - iOrgMvBits + iBestMvBits;
    ruiCost = (ruiCost - m_pcRdCost->getCost( uiOrgBits ))  + m_pcRdCost->getCost( ruiBits );
  }
}


Distortion InterSearch::xGetTemplateCost( const PredictionUnit& pu,
                                          PelUnitBuf& origBuf,
                                          PelUnitBuf& predBuf,
                                          Mv          cMvCand,
                                          int         iMVPIdx,
                                          int         iMVPNum,
                                          RefPicList  eRefPicList,
                                          int         iRefIdx
)
{
  Distortion uiCost = std::numeric_limits<Distortion>::max();

  const Picture* picRef = pu.cu->slice->getRefPic( eRefPicList, iRefIdx );
  clipMv( cMvCand, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );

  // prediction pattern
  const bool bi = pu.cu->slice->testWeightPred() && pu.cu->slice->getSliceType()==P_SLICE
#if INTER_LIC
    && !pu.cu->LICFlag
#endif
    ;


  xPredInterBlk( COMPONENT_Y, pu, picRef, cMvCand, predBuf, bi, pu.cu->slice->clpRng( COMPONENT_Y )
                , false
                , false
                );

  if ( bi )
  {
    xWeightedPredictionUni( pu, predBuf, eRefPicList, predBuf, iRefIdx, m_maxCompIDToPred );
  }

  // calc distortion

  uiCost = m_pcRdCost->getDistPart(origBuf.Y(), predBuf.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, DF_SAD);
  uiCost += m_pcRdCost->getCost( m_auiMVPIdxCost[iMVPIdx][iMVPNum] );

  return uiCost;
}

Distortion InterSearch::xGetAffineTemplateCost( PredictionUnit& pu, PelUnitBuf& origBuf, PelUnitBuf& predBuf, Mv acMvCand[3], int iMVPIdx, int iMVPNum, RefPicList eRefPicList, int iRefIdx )
{
  Distortion uiCost = std::numeric_limits<Distortion>::max();

  const Picture* picRef = pu.cu->slice->getRefPic( eRefPicList, iRefIdx );
#if INTER_LIC && RPR_ENABLE
  // xPredAffineBlk may call PU::checkRprLicCondition()
#if JVET_Y0128_NON_CTC
  pu.interDir = (uint8_t)(eRefPicList + 1);
#endif
  pu.refIdx[eRefPicList]    = iRefIdx;
  pu.refIdx[1-eRefPicList]  = NOT_VALID;
#endif

  // prediction pattern
  const bool bi = pu.cu->slice->testWeightPred() && pu.cu->slice->getSliceType()==P_SLICE
#if INTER_LIC
    && !pu.cu->LICFlag
#endif
    ;
  Mv mv[3];
  memcpy(mv, acMvCand, sizeof(mv));
  m_iRefListIdx = eRefPicList;
#if JVET_Z0136_OOB
  xPredAffineBlk(COMPONENT_Y, pu, picRef, mv, predBuf, bi, pu.cu->slice->clpRng(COMPONENT_Y), eRefPicList);
#else
  xPredAffineBlk(COMPONENT_Y, pu, picRef, mv, predBuf, bi, pu.cu->slice->clpRng(COMPONENT_Y));
#endif
  if( bi )
  {
    xWeightedPredictionUni( pu, predBuf, eRefPicList, predBuf, iRefIdx, m_maxCompIDToPred );
  }

  // calc distortion
  enum DFunc distFunc = (pu.cs->slice->getDisableSATDForRD()) ? DF_SAD : DF_HAD;
  uiCost  = m_pcRdCost->getDistPart( origBuf.Y(), predBuf.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y
    , distFunc
  );
  uiCost += m_pcRdCost->getCost( m_auiMVPIdxCost[iMVPIdx][iMVPNum] );
  DTRACE( g_trace_ctx, D_COMMON, " (%d) affineTemplateCost=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiCost );
  return uiCost;
}

#if MULTI_HYP_PRED
void InterSearch::xMotionEstimation(PredictionUnit& pu, PelUnitBuf& origBuf, RefPicList eRefPicList, Mv& rcMvPred, int iRefIdxPred, Mv& rcMv, int& riMVPIdx, uint32_t& ruiBits, Distortion& ruiCost, const AMVPInfo& amvpInfo, bool bBi, int weight)
#else
void InterSearch::xMotionEstimation(PredictionUnit& pu, PelUnitBuf& origBuf, RefPicList eRefPicList, Mv& rcMvPred, int iRefIdxPred, Mv& rcMv, int& riMVPIdx, uint32_t& ruiBits, Distortion& ruiCost, const AMVPInfo& amvpInfo, bool bBi)
#endif
{
#if MULTI_HYP_PRED
  if (!weight)
#endif
  if( pu.cu->cs->sps->getUseBcw() && pu.cu->BcwIdx != BCW_DEFAULT && !bBi && xReadBufferedUniMv(pu, eRefPicList, iRefIdxPred, rcMvPred, rcMv, ruiBits, ruiCost) )
  {
    return;
  }

  Mv cMvHalf, cMvQter;
  CHECK(eRefPicList >= MAX_NUM_REF_LIST_ADAPT_SR || iRefIdxPred>=int(MAX_IDX_ADAPT_SR), "Invalid reference picture list");
  m_iSearchRange = m_aaiAdaptSR[eRefPicList][iRefIdxPred];
#if MULTI_HYP_PRED 
  if (weight)
  {
    m_iSearchRange = std::min(m_iSearchRange, MULTI_HYP_PRED_SEARCH_RANGE);
  }
#endif

  int    iSrchRng   = (bBi ? m_bipredSearchRange : m_iSearchRange);
  double fWeight    = 1.0;

  PelUnitBuf  origBufTmp = m_tmpStorageLCU.getBuf( UnitAreaRelative(*pu.cu, pu) );
  PelUnitBuf* pBuf       = &origBuf;

  if(bBi) // Bi-predictive ME
  {
#if MULTI_HYP_PRED
    CHECK(weight, "Multi Hyp: bBi");
#endif
    // NOTE: Other buf contains predicted signal from another direction
    PelUnitBuf otherBuf = m_tmpPredStorage[1 - (int)eRefPicList].getBuf( UnitAreaRelative(*pu.cu, pu ));
    origBufTmp.copyFrom(origBuf);
    origBufTmp.removeHighFreq( otherBuf, m_pcEncCfg->getClipForBiPredMeEnabled(), pu.cu->slice->clpRngs()
                              ,getBcwWeight( pu.cu->BcwIdx, eRefPicList )
                              );
    pBuf = &origBufTmp;

    fWeight = xGetMEDistortionWeight( pu.cu->BcwIdx, eRefPicList );
  }
#if MULTI_HYP_PRED
  else if (weight)
  {
    CHECK(bBi, "Multi Hyp: bBi");
    fWeight = fabs(double(weight) / double(1 << MULTI_HYP_PRED_WEIGHT_BITS));
  }
#endif
  m_cDistParam.isBiPred = bBi;
#if INTER_LIC
  m_cDistParam.useMR = pu.cu->LICFlag;
#endif

  //  Search key pattern initialization
  CPelBuf  tmpPattern   = pBuf->Y();
  CPelBuf* pcPatternKey = &tmpPattern;

  m_lumaClpRng = pu.cs->slice->clpRng( COMPONENT_Y );

  bool wrap =  pu.cu->slice->getRefPic(eRefPicList, iRefIdxPred)->isWrapAroundEnabled( pu.cs->pps );
  CPelBuf buf = pu.cu->slice->getRefPic(eRefPicList, iRefIdxPred)->getRecoBuf(pu.blocks[COMPONENT_Y], wrap);

  IntTZSearchStruct cStruct;
  cStruct.pcPatternKey  = pcPatternKey;
  cStruct.iRefStride    = buf.stride;
  cStruct.piRefY        = buf.buf;
  cStruct.imvShift = pu.cu->imv == IMV_HPEL ? 1 : (pu.cu->imv << 1);
  cStruct.useAltHpelIf = pu.cu->imv == IMV_HPEL;
  cStruct.inCtuSearch = false;
  cStruct.zeroMV = false;
  {
    if (m_useCompositeRef && pu.cs->slice->getRefPic(eRefPicList, iRefIdxPred)->longTerm)
    {
      cStruct.inCtuSearch = true;
    }
  }

  auto blkCache = dynamic_cast<CacheBlkInfoCtrl*>( m_modeCtrl );

  bool bQTBTMV  = false;
  bool bQTBTMV2 = false;
  Mv cIntMv;
#if MULTI_HYP_PRED 
  if (!bBi && !weight)
#else
  if (!bBi)
#endif
  {
    bool bValid = blkCache && blkCache->getMv( pu, eRefPicList, iRefIdxPred, cIntMv );
    if( bValid )
    {
      bQTBTMV2 = true;
      cIntMv.changePrecision( MV_PRECISION_INT, MV_PRECISION_INTERNAL);
    }
  }

  Mv predQuarter = rcMvPred;
  predQuarter.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
  m_pcRdCost->setPredictor( predQuarter );

  m_pcRdCost->setCostScale(2);

#if INTER_LIC
  if (pu.cu->LICFlag)
  {
    m_cDistParam.applyWeight = false;
  }
  else
#endif
  {
    setWpScalingDistParam(iRefIdxPred, eRefPicList, pu.cu->slice);
  }
  m_currRefPicList = eRefPicList;
  m_currRefPicIndex = iRefIdxPred;
  m_skipFracME = false;
  //  Do integer search
  if( ( m_motionEstimationSearchMethod == MESEARCH_FULL ) || bBi || bQTBTMV )
  {
    cStruct.subShiftMode = m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE3 ? 2 : 0;
    m_pcRdCost->setDistParam(m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, cStruct.subShiftMode);

    Mv bestInitMv = (bBi ? rcMv : rcMvPred);
    Mv cTmpMv = bestInitMv;

    clipMv( cTmpMv, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
    cTmpMv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
    m_cDistParam.cur.buf = cStruct.piRefY + (cTmpMv.ver * cStruct.iRefStride) + cTmpMv.hor;
    Distortion uiBestSad = m_cDistParam.distFunc(m_cDistParam);
    uiBestSad += m_pcRdCost->getCostOfVectorWithPredictor(cTmpMv.hor, cTmpMv.ver, cStruct.imvShift);
#if JVET_X0083_BM_AMVP_MERGE_MODE
    if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
    {
      cTmpMv = rcMvPred;
      clipMv( cTmpMv, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
      cTmpMv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
      m_cDistParam.cur.buf = cStruct.piRefY + (cTmpMv.ver * cStruct.iRefStride) + cTmpMv.hor;
      Distortion uiSad = m_cDistParam.distFunc(m_cDistParam);
      uiSad += m_pcRdCost->getCostOfVectorWithPredictor(cTmpMv.hor, cTmpMv.ver, cStruct.imvShift);
      if (uiSad < uiBestSad)
      {
        uiBestSad = uiSad;
        bestInitMv = rcMvPred;
        m_cDistParam.maximumDistortionForEarlyExit = uiSad;
      }
    }
#endif

#if AMVR_ENC_OPT
    const MvPrecision tmpIntMvPrec = (pu.cu->imv == IMV_4PEL ? MV_PRECISION_4PEL : MV_PRECISION_INT);
#endif
    for (int i = 0; i < m_uniMvListSize; i++)
    {
      BlkUniMvInfo* curMvInfo = m_uniMvList + ((m_uniMvListIdx - 1 - i + m_uniMvListMaxSize) % (m_uniMvListMaxSize));
#if AMVR_ENC_OPT
      Mv tmpCurMv = curMvInfo->uniMvs[eRefPicList][iRefIdxPred];
      tmpCurMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);
#endif

      int j = 0;
      for (; j < i; j++)
      {
        BlkUniMvInfo *prevMvInfo = m_uniMvList + ((m_uniMvListIdx - 1 - j + m_uniMvListMaxSize) % (m_uniMvListMaxSize));
#if AMVR_ENC_OPT
        Mv tmpPrevMv = prevMvInfo->uniMvs[eRefPicList][iRefIdxPred];
        tmpPrevMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);

        if (tmpCurMv == tmpPrevMv)
#else
        if (curMvInfo->uniMvs[eRefPicList][iRefIdxPred] == prevMvInfo->uniMvs[eRefPicList][iRefIdxPred])
#endif
        {
          break;
        }
      }
      if (j < i)
        continue;

      cTmpMv = curMvInfo->uniMvs[eRefPicList][iRefIdxPred];
      clipMv( cTmpMv, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
#if AMVR_ENC_OPT
      cTmpMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);
      if (tmpIntMvPrec != MV_PRECISION_INT)
      {
        cTmpMv.changePrecision(tmpIntMvPrec, MV_PRECISION_INT);
      }
#else
      cTmpMv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
#endif
      m_cDistParam.cur.buf = cStruct.piRefY + (cTmpMv.ver * cStruct.iRefStride) + cTmpMv.hor;

      Distortion uiSad = m_cDistParam.distFunc(m_cDistParam);
      uiSad += m_pcRdCost->getCostOfVectorWithPredictor(cTmpMv.hor, cTmpMv.ver, cStruct.imvShift);
      if (uiSad < uiBestSad)
      {
        uiBestSad = uiSad;
        bestInitMv = curMvInfo->uniMvs[eRefPicList][iRefIdxPred];
        m_cDistParam.maximumDistortionForEarlyExit = uiSad;
      }
    }

    if( !bQTBTMV )
    {
      xSetSearchRange(pu, bestInitMv, iSrchRng, cStruct.searchRange, cStruct);
    }
    xPatternSearch( cStruct, rcMv, ruiCost);
  }
  else if( bQTBTMV2 )
  {
    rcMv = cIntMv;

    cStruct.subShiftMode = ( !m_pcEncCfg->getRestrictMESampling() && m_pcEncCfg->getMotionEstimationSearchMethod() == MESEARCH_SELECTIVE ) ? 1 :
                            ( m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE3 ) ? 2 : 0;
    xTZSearch(pu, eRefPicList, iRefIdxPred, cStruct, rcMv, ruiCost, NULL, false, true);
  }
  else
  {
    cStruct.subShiftMode = ( !m_pcEncCfg->getRestrictMESampling() && m_pcEncCfg->getMotionEstimationSearchMethod() == MESEARCH_SELECTIVE ) ? 1 :
                            ( m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE3 ) ? 2 : 0;
#if MULTI_HYP_PRED
    if (weight == 0)
#endif
    rcMv = rcMvPred;
    const Mv *pIntegerMv2Nx2NPred = 0;
#if MULTI_HYP_PRED
    const auto savedMEMethod = m_motionEstimationSearchMethod;
    if( weight )
    {
      m_motionEstimationSearchMethod = MESEARCH_DIAMOND_ENHANCED;
    }
#endif
    xPatternSearchFast(pu, eRefPicList, iRefIdxPred, cStruct, rcMv, ruiCost, pIntegerMv2Nx2NPred);
#if MULTI_HYP_PRED
    if( weight )
    {
      m_motionEstimationSearchMethod = savedMEMethod;
    }
    else
    {
#endif
    if( blkCache )
    {
      blkCache->setMv( pu.cs->area, eRefPicList, iRefIdxPred, rcMv );
    }
    else
    {
      m_integerMv2Nx2N[eRefPicList][iRefIdxPred] = rcMv;
    }
#if MULTI_HYP_PRED
    }
#endif
  }
  DTRACE( g_trace_ctx, D_ME, "%d %d %d :MECostFPel<L%d,%d>: %d,%d,%dx%d, %d", DTRACE_GET_COUNTER( g_trace_ctx, D_ME ), pu.cu->slice->getPOC(), 0, ( int ) eRefPicList, ( int ) bBi, pu.Y().x, pu.Y().y, pu.Y().width, pu.Y().height, ruiCost );
  // sub-pel refinement for sub-pel resolution
  if ( pu.cu->imv == 0 || pu.cu->imv == IMV_HPEL )
  {
    if( m_pcEncCfg->getMCTSEncConstraint() )
    {
      Area curTileAreaSubPelRestricted = pu.cs->picture->mctsInfo.getTileAreaSubPelRestricted( pu );
      // Area adjustment, because subpel refinement is going to (x-1;y-1) direction
      curTileAreaSubPelRestricted.x += 1;
      curTileAreaSubPelRestricted.y += 1;
      curTileAreaSubPelRestricted.width -= 1;
      curTileAreaSubPelRestricted.height -= 1;
      if( ! MCTSHelper::checkMvIsNotInRestrictedArea( pu, rcMv, curTileAreaSubPelRestricted, MV_PRECISION_INT ) )
      {
        MCTSHelper::clipMvToArea( rcMv, pu.Y(), curTileAreaSubPelRestricted, *pu.cs->sps, 0 );
      }
    }
    xPatternSearchFracDIF( pu, eRefPicList, iRefIdxPred, cStruct, rcMv, cMvHalf, cMvQter, ruiCost);
    m_pcRdCost->setCostScale( 0 );
    rcMv <<= 2;
    rcMv  += ( cMvHalf <<= 1 );
    rcMv  += cMvQter;
    uint32_t uiMvBits = m_pcRdCost->getBitsOfVectorWithPredictor( rcMv.getHor(), rcMv.getVer(), cStruct.imvShift );
    ruiBits += uiMvBits;
    ruiCost = ( Distortion ) ( floor( fWeight * ( ( double ) ruiCost - ( double ) m_pcRdCost->getCost( uiMvBits ) ) ) + ( double ) m_pcRdCost->getCost( ruiBits ) );
    rcMv.changePrecision(MV_PRECISION_QUARTER, MV_PRECISION_INTERNAL);
  }
  else // integer refinement for integer-pel and 4-pel resolution
  {
    rcMv.changePrecision(MV_PRECISION_INT, MV_PRECISION_INTERNAL);
    xPatternSearchIntRefine( pu, cStruct, rcMv, rcMvPred, riMVPIdx, ruiBits, ruiCost, amvpInfo, fWeight);
  }

#if INTER_LIC
  if (pu.cu->LICFlag)
  {
    PelUnitBuf predTempBuf = m_tmpStorageLCU.getBuf(UnitAreaRelative(*pu.cu, pu));
    const Picture* picRef = pu.cu->slice->getRefPic(eRefPicList, iRefIdxPred);
    xPredInterBlk(COMPONENT_Y, pu, picRef, rcMv, predTempBuf, false, pu.cu->slice->clpRng(COMPONENT_Y), false, false);

    DistParam distParam;
    m_pcRdCost->setDistParam(distParam, origBuf.Y(), predTempBuf.Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, !pu.cs->slice->getDisableSATDForRD());
#if MULTI_HYP_PRED
    ruiCost = (Distortion)floor(fWeight * (double)distParam.distFunc(distParam)) + m_pcRdCost->getCost(ruiBits);
#else
    ruiCost = distParam.distFunc(distParam) + m_pcRdCost->getCost(ruiBits);
#endif
  }
#endif
  DTRACE(g_trace_ctx, D_ME, "   MECost<L%d,%d>: %6d (%d)  MV:%d,%d\n", (int)eRefPicList, (int)bBi, ruiCost, ruiBits, rcMv.getHor() << 2, rcMv.getVer() << 2);
}



void InterSearch::xSetSearchRange ( const PredictionUnit& pu,
                                    const Mv& cMvPred,
                                    const int iSrchRng,
                                    SearchRange& sr
                                  , IntTZSearchStruct& cStruct
)
{
  const int iMvShift = MV_FRACTIONAL_BITS_INTERNAL;
  Mv cFPMvPred = cMvPred;
  clipMv( cFPMvPred, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
  
  Mv mvTL(cFPMvPred.getHor() - (iSrchRng << iMvShift), cFPMvPred.getVer() - (iSrchRng << iMvShift));
  Mv mvBR(cFPMvPred.getHor() + (iSrchRng << iMvShift), cFPMvPred.getVer() + (iSrchRng << iMvShift));

  if (m_pcEncCfg->getMCTSEncConstraint())
  {
    MCTSHelper::clipMvToArea( mvTL, pu.Y(), pu.cs->picture->mctsInfo.getTileArea(), *pu.cs->sps );
    MCTSHelper::clipMvToArea( mvBR, pu.Y(), pu.cs->picture->mctsInfo.getTileArea(), *pu.cs->sps );
  }
  else
  {
    xClipMv( mvTL, pu.cu->lumaPos(),
            pu.cu->lumaSize(),
            *pu.cs->sps
          , *pu.cs->pps
    );
    xClipMv( mvBR, pu.cu->lumaPos(),
            pu.cu->lumaSize(),
            *pu.cs->sps
          , *pu.cs->pps
    );
  }

  mvTL.divideByPowerOf2( iMvShift );
  mvBR.divideByPowerOf2( iMvShift );

  sr.left   = mvTL.hor;
  sr.top    = mvTL.ver;
  sr.right  = mvBR.hor;
  sr.bottom = mvBR.ver;

  if (m_useCompositeRef && cStruct.inCtuSearch)
  {
    Position posRB = pu.Y().bottomRight();
    Position posTL = pu.Y().topLeft();
    const PreCalcValues *pcv = pu.cs->pcv;
    Position posRBinCTU(posRB.x & pcv->maxCUWidthMask, posRB.y & pcv->maxCUHeightMask);
    Position posLTinCTU = Position(posTL.x & pcv->maxCUWidthMask, posTL.y & pcv->maxCUHeightMask).offset(-4, -4);
    if (sr.left < -posLTinCTU.x)
      sr.left = -posLTinCTU.x;
    if (sr.top < -posLTinCTU.y)
      sr.top = -posLTinCTU.y;
    if (sr.right >((int)pcv->maxCUWidth - 4 - posRBinCTU.x))
      sr.right = (int)pcv->maxCUWidth - 4 - posRBinCTU.x;
    if (sr.bottom >((int)pcv->maxCUHeight - 4 - posRBinCTU.y))
      sr.bottom = (int)pcv->maxCUHeight - 4 - posRBinCTU.y;
    if (posLTinCTU.x == -4 || posLTinCTU.y == -4)
    {
      sr.left = sr.right = sr.bottom = sr.top = 0;
      cStruct.zeroMV = 1;
    }
    if (posRBinCTU.x == pcv->maxCUWidthMask || posRBinCTU.y == pcv->maxCUHeightMask)
    {
      sr.left = sr.right = sr.bottom = sr.top = 0;
      cStruct.zeroMV = 1;
    }
  }
}


void InterSearch::xPatternSearch( IntTZSearchStruct&    cStruct,
                                  Mv&            rcMv,
                                  Distortion&    ruiSAD )
{
  Distortion  uiSad;
  Distortion  uiSadBest = std::numeric_limits<Distortion>::max();
  int         iBestX = 0;
  int         iBestY = 0;

  //-- jclee for using the SAD function pointer
  m_pcRdCost->setDistParam( m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, cStruct.subShiftMode );

  const SearchRange& sr = cStruct.searchRange;

  const Pel* piRef = cStruct.piRefY + (sr.top * cStruct.iRefStride);
  for ( int y = sr.top; y <= sr.bottom; y++ )
  {
    for ( int x = sr.left; x <= sr.right; x++ )
    {
      //  find min. distortion position
      m_cDistParam.cur.buf = piRef + x;

      uiSad = m_cDistParam.distFunc( m_cDistParam );

      // motion cost
      uiSad += m_pcRdCost->getCostOfVectorWithPredictor( x, y, cStruct.imvShift );

      if ( uiSad < uiSadBest )
      {
        uiSadBest = uiSad;
        iBestX    = x;
        iBestY    = y;
        m_cDistParam.maximumDistortionForEarlyExit = uiSad;
      }
    }
    piRef += cStruct.iRefStride;
  }
  rcMv.set( iBestX, iBestY );

  cStruct.uiBestSad = uiSadBest; // th for testing
  ruiSAD = uiSadBest - m_pcRdCost->getCostOfVectorWithPredictor( iBestX, iBestY, cStruct.imvShift );
  return;
}


void InterSearch::xPatternSearchFast( const PredictionUnit& pu,
                                      RefPicList            eRefPicList,
                                      int                   iRefIdxPred,
                                      IntTZSearchStruct&    cStruct,
                                      Mv&                   rcMv,
                                      Distortion&           ruiSAD,
                                      const Mv* const       pIntegerMv2Nx2NPred )
{
  switch ( m_motionEstimationSearchMethod )
  {
  case MESEARCH_DIAMOND:
    xTZSearch         ( pu, eRefPicList, iRefIdxPred, cStruct, rcMv, ruiSAD, pIntegerMv2Nx2NPred, false );
    break;

  case MESEARCH_SELECTIVE:
    xTZSearchSelective( pu, eRefPicList, iRefIdxPred, cStruct, rcMv, ruiSAD, pIntegerMv2Nx2NPred );
    break;

  case MESEARCH_DIAMOND_ENHANCED:
    xTZSearch         ( pu, eRefPicList, iRefIdxPred, cStruct, rcMv, ruiSAD, pIntegerMv2Nx2NPred, true );
    break;

  case MESEARCH_FULL: // shouldn't get here.
  default:
    break;
  }
}


void InterSearch::xTZSearch( const PredictionUnit& pu,
                             RefPicList            eRefPicList,
                             int                   iRefIdxPred,
                             IntTZSearchStruct&    cStruct,
                             Mv&                   rcMv,
                             Distortion&           ruiSAD,
                             const Mv* const       pIntegerMv2Nx2NPred,
                             const bool            bExtendedSettings,
                             const bool            bFastSettings)
{
  const bool bUseRasterInFastMode                    = true; //toggle this to further reduce runtime

  const bool bUseAdaptiveRaster                      = bExtendedSettings;
  const int  iRaster                                 = (bFastSettings && bUseRasterInFastMode) ? 8 : 5;
  const bool bTestZeroVector                         = true && !bFastSettings;
  const bool bTestZeroVectorStart                    = bExtendedSettings;
  const bool bTestZeroVectorStop                     = false;
  const bool bFirstSearchDiamond                     = true;  // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const bool bFirstCornersForDiamondDist1            = bExtendedSettings;
  const bool bFirstSearchStop                        = m_pcEncCfg->getFastMEAssumingSmootherMVEnabled();
  const uint32_t uiFirstSearchRounds                     = bFastSettings ? (bUseRasterInFastMode?3:2) : 3;     // first search stop X rounds after best match (must be >=1)
  const bool bEnableRasterSearch                     = bFastSettings ? bUseRasterInFastMode : true;
  const bool bAlwaysRasterSearch                     = bExtendedSettings;  // true: BETTER but factor 2 slower
  const bool bRasterRefinementEnable                 = false; // enable either raster refinement or star refinement
  const bool bRasterRefinementDiamond                = false; // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const bool bRasterRefinementCornersForDiamondDist1 = bExtendedSettings;
  const bool bStarRefinementEnable                   = true;  // enable either star refinement or raster refinement
  const bool bStarRefinementDiamond                  = true;  // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const bool bStarRefinementCornersForDiamondDist1   = bExtendedSettings;
  const bool bStarRefinementStop                     = false || bFastSettings;
  const uint32_t uiStarRefinementRounds                  = 2;  // star refinement stop X rounds after best match (must be >=1)
  const bool bNewZeroNeighbourhoodTest               = bExtendedSettings;

  int iSearchRange = m_iSearchRange;
  if( m_pcEncCfg->getMCTSEncConstraint() )
  {
    MCTSHelper::clipMvToArea( rcMv, pu.Y(), pu.cs->picture->mctsInfo.getTileArea(), *pu.cs->sps );
  }
  else
  {
    clipMv( rcMv, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
  }
  rcMv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
  rcMv.divideByPowerOf2(2);

  // init TZSearchStruct
  cStruct.uiBestSad = std::numeric_limits<Distortion>::max();

  //
  m_cDistParam.maximumDistortionForEarlyExit = cStruct.uiBestSad;
  m_pcRdCost->setDistParam( m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, cStruct.subShiftMode );

  // distortion


  // set rcMv (Median predictor) as start point and as best point
  xTZSearchHelp( cStruct, rcMv.getHor(), rcMv.getVer(), 0, 0 );

  // test whether zero Mv is better start point than Median predictor
  if ( bTestZeroVector )
  {
    if ((rcMv.getHor() != 0 || rcMv.getVer() != 0) &&
      (0 != cStruct.iBestX || 0 != cStruct.iBestY))
    {
      // only test 0-vector if not obviously previously tested.
      xTZSearchHelp( cStruct, 0, 0, 0, 0 );
    }
  }

  SearchRange& sr = cStruct.searchRange;

  if (pIntegerMv2Nx2NPred != 0)
  {
    Mv integerMv2Nx2NPred = *pIntegerMv2Nx2NPred;
    integerMv2Nx2NPred.changePrecision(MV_PRECISION_INT, MV_PRECISION_INTERNAL);
    if( m_pcEncCfg->getMCTSEncConstraint() )
    {
      MCTSHelper::clipMvToArea( integerMv2Nx2NPred, pu.Y(), pu.cs->picture->mctsInfo.getTileArea(), *pu.cs->sps );
    }
    else
    {
      clipMv( integerMv2Nx2NPred, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
    }
    integerMv2Nx2NPred.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
    integerMv2Nx2NPred.divideByPowerOf2(2);

    if ((rcMv != integerMv2Nx2NPred) &&
      (integerMv2Nx2NPred.getHor() != cStruct.iBestX || integerMv2Nx2NPred.getVer() != cStruct.iBestY))
    {
      // only test integerMv2Nx2NPred if not obviously previously tested.
      xTZSearchHelp( cStruct, integerMv2Nx2NPred.getHor(), integerMv2Nx2NPred.getVer(), 0, 0);
    }
  }

#if AMVR_ENC_OPT
  const MvPrecision tmpIntMvPrec = (pu.cu->imv == IMV_4PEL ? MV_PRECISION_4PEL : MV_PRECISION_INT);
#endif
  for (int i = 0; i < m_uniMvListSize; i++)
  {
    BlkUniMvInfo* curMvInfo = m_uniMvList + ((m_uniMvListIdx - 1 - i + m_uniMvListMaxSize) % (m_uniMvListMaxSize));
#if AMVR_ENC_OPT
    Mv tmpCurMv = curMvInfo->uniMvs[eRefPicList][iRefIdxPred];
    tmpCurMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);
#endif

    int j = 0;
    for (; j < i; j++)
    {
      BlkUniMvInfo *prevMvInfo = m_uniMvList + ((m_uniMvListIdx - 1 - j + m_uniMvListMaxSize) % (m_uniMvListMaxSize));
#if AMVR_ENC_OPT
      Mv tmpPrevMv = prevMvInfo->uniMvs[eRefPicList][iRefIdxPred];
      tmpPrevMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);

      if (tmpCurMv == tmpPrevMv)
#else
      if (curMvInfo->uniMvs[eRefPicList][iRefIdxPred] == prevMvInfo->uniMvs[eRefPicList][iRefIdxPred])
#endif
      {
        break;
      }
    }
    if (j < i)
      continue;

    Mv cTmpMv = curMvInfo->uniMvs[eRefPicList][iRefIdxPred];
    clipMv( cTmpMv, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
#if AMVR_ENC_OPT
    cTmpMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);
    if (tmpIntMvPrec != MV_PRECISION_INT)
    {
      cTmpMv.changePrecision(tmpIntMvPrec, MV_PRECISION_INT);
    }
#else
    cTmpMv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
#endif
    m_cDistParam.cur.buf = cStruct.piRefY + (cTmpMv.ver * cStruct.iRefStride) + cTmpMv.hor;

    Distortion uiSad = m_cDistParam.distFunc(m_cDistParam);
    uiSad += m_pcRdCost->getCostOfVectorWithPredictor(cTmpMv.hor, cTmpMv.ver, cStruct.imvShift);
    if (uiSad < cStruct.uiBestSad)
    {
      cStruct.uiBestSad = uiSad;
      cStruct.iBestX = cTmpMv.hor;
      cStruct.iBestY = cTmpMv.ver;
      m_cDistParam.maximumDistortionForEarlyExit = uiSad;
    }
  }

  {
    // set search range
    Mv currBestMv(cStruct.iBestX, cStruct.iBestY );
    currBestMv <<= MV_FRACTIONAL_BITS_INTERNAL;
    xSetSearchRange(pu, currBestMv, m_iSearchRange >> (bFastSettings ? 1 : 0), sr, cStruct);
  }
  if (m_pcEncCfg->getUseHashME() && (m_currRefPicList == 0 || pu.cu->slice->getList1IdxToList0Idx(m_currRefPicIndex) < 0))
  {
    int minSize = min(pu.cu->lumaSize().width, pu.cu->lumaSize().height);
    if (minSize < 128 && minSize >= 4)
    {
      int numberOfOtherMvps = m_numHashMVStoreds[m_currRefPicList][m_currRefPicIndex];
      for (int i = 0; i < numberOfOtherMvps; i++)
      {
        xTZSearchHelp(cStruct, m_hashMVStoreds[m_currRefPicList][m_currRefPicIndex][i].getHor(), m_hashMVStoreds[m_currRefPicList][m_currRefPicIndex][i].getVer(), 0, 0);
      }
      if (numberOfOtherMvps > 0)
      {
        // write out best match
        rcMv.set(cStruct.iBestX, cStruct.iBestY);
        ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCostOfVectorWithPredictor(cStruct.iBestX, cStruct.iBestY, cStruct.imvShift);
        m_skipFracME = true;
        return;
      }
    }
  }

  // start search
  int  iDist = 0;
  int  iStartX = cStruct.iBestX;
  int  iStartY = cStruct.iBestY;

  const bool bBestCandidateZero = (cStruct.iBestX == 0) && (cStruct.iBestY == 0);

  // first search around best position up to now.
  // The following works as a "subsampled/log" window search around the best candidate
  for ( iDist = 1; iDist <= iSearchRange; iDist*=2 )
  {
    if ( bFirstSearchDiamond == 1 )
    {
      xTZ8PointDiamondSearch ( cStruct, iStartX, iStartY, iDist, bFirstCornersForDiamondDist1 );
    }
    else
    {
      xTZ8PointSquareSearch  ( cStruct, iStartX, iStartY, iDist );
    }

    if ( bFirstSearchStop && ( cStruct.uiBestRound >= uiFirstSearchRounds ) ) // stop criterion
    {
      break;
    }
  }

  if (!bNewZeroNeighbourhoodTest)
  {
    // test whether zero Mv is a better start point than Median predictor
    if ( bTestZeroVectorStart && ((cStruct.iBestX != 0) || (cStruct.iBestY != 0)) )
    {
      xTZSearchHelp( cStruct, 0, 0, 0, 0 );
      if ( (cStruct.iBestX == 0) && (cStruct.iBestY == 0) )
      {
        // test its neighborhood
        for ( iDist = 1; iDist <= iSearchRange; iDist*=2 )
        {
          xTZ8PointDiamondSearch( cStruct, 0, 0, iDist, false );
          if ( bTestZeroVectorStop && (cStruct.uiBestRound > 0) ) // stop criterion
          {
            break;
          }
        }
      }
    }
  }
  else
  {
    // Test also zero neighbourhood but with half the range
    // It was reported that the original (above) search scheme using bTestZeroVectorStart did not
    // make sense since one would have already checked the zero candidate earlier
    // and thus the conditions for that test would have not been satisfied
    if (bTestZeroVectorStart == true && bBestCandidateZero != true)
    {
      for ( iDist = 1; iDist <= (iSearchRange >> 1); iDist*=2 )
      {
        xTZ8PointDiamondSearch( cStruct, 0, 0, iDist, false );
        if ( bTestZeroVectorStop && (cStruct.uiBestRound > 2) ) // stop criterion
        {
          break;
        }
      }
    }
  }

  // calculate only 2 missing points instead 8 points if cStruct.uiBestDistance == 1
  if ( cStruct.uiBestDistance == 1 )
  {
    cStruct.uiBestDistance = 0;
    xTZ2PointSearch( cStruct );
  }

  // raster search if distance is too big
  if (bUseAdaptiveRaster)
  {
    int iWindowSize     = iRaster;
    SearchRange localsr = sr;

    if (!(bEnableRasterSearch && ( ((int)(cStruct.uiBestDistance) >= iRaster))))
    {
      iWindowSize ++;
      localsr.left   /= 2;
      localsr.right  /= 2;
      localsr.top    /= 2;
      localsr.bottom /= 2;
    }
    cStruct.uiBestDistance = iWindowSize;
    for ( iStartY = localsr.top; iStartY <= localsr.bottom; iStartY += iWindowSize )
    {
      for ( iStartX = localsr.left; iStartX <= localsr.right; iStartX += iWindowSize )
      {
        xTZSearchHelp( cStruct, iStartX, iStartY, 0, iWindowSize );
      }
    }
  }
  else
  {
    if ( bEnableRasterSearch && ( ((int)(cStruct.uiBestDistance) >= iRaster) || bAlwaysRasterSearch ) )
    {
      cStruct.uiBestDistance = iRaster;
      for ( iStartY = sr.top; iStartY <= sr.bottom; iStartY += iRaster )
      {
        for ( iStartX = sr.left; iStartX <= sr.right; iStartX += iRaster )
        {
          xTZSearchHelp( cStruct, iStartX, iStartY, 0, iRaster );
        }
      }
    }
  }

  // raster refinement

  if ( bRasterRefinementEnable && cStruct.uiBestDistance > 0 )
  {
    while ( cStruct.uiBestDistance > 0 )
    {
      iStartX = cStruct.iBestX;
      iStartY = cStruct.iBestY;
      if ( cStruct.uiBestDistance > 1 )
      {
        iDist = cStruct.uiBestDistance >>= 1;
        if ( bRasterRefinementDiamond == 1 )
        {
          xTZ8PointDiamondSearch ( cStruct, iStartX, iStartY, iDist, bRasterRefinementCornersForDiamondDist1 );
        }
        else
        {
          xTZ8PointSquareSearch  ( cStruct, iStartX, iStartY, iDist );
        }
      }

      // calculate only 2 missing points instead 8 points if cStruct.uiBestDistance == 1
      if ( cStruct.uiBestDistance == 1 )
      {
        cStruct.uiBestDistance = 0;
        if ( cStruct.ucPointNr != 0 )
        {
          xTZ2PointSearch( cStruct );
        }
      }
    }
  }

  // star refinement
  if ( bStarRefinementEnable && cStruct.uiBestDistance > 0 )
  {
    while ( cStruct.uiBestDistance > 0 )
    {
      iStartX = cStruct.iBestX;
      iStartY = cStruct.iBestY;
      cStruct.uiBestDistance = 0;
      cStruct.ucPointNr = 0;
      for ( iDist = 1; iDist < iSearchRange + 1; iDist*=2 )
      {
        if ( bStarRefinementDiamond == 1 )
        {
          xTZ8PointDiamondSearch ( cStruct, iStartX, iStartY, iDist, bStarRefinementCornersForDiamondDist1 );
        }
        else
        {
          xTZ8PointSquareSearch  ( cStruct, iStartX, iStartY, iDist );
        }
        if ( bStarRefinementStop && (cStruct.uiBestRound >= uiStarRefinementRounds) ) // stop criterion
        {
          break;
        }
      }

      // calculate only 2 missing points instead 8 points if cStrukt.uiBestDistance == 1
      if ( cStruct.uiBestDistance == 1 )
      {
        cStruct.uiBestDistance = 0;
        if ( cStruct.ucPointNr != 0 )
        {
          xTZ2PointSearch( cStruct );
        }
      }
    }
  }

  // write out best match
  rcMv.set( cStruct.iBestX, cStruct.iBestY );
  ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCostOfVectorWithPredictor( cStruct.iBestX, cStruct.iBestY, cStruct.imvShift );
}


void InterSearch::xTZSearchSelective( const PredictionUnit& pu,
                                      RefPicList            eRefPicList,
                                      int                   iRefIdxPred,
                                      IntTZSearchStruct&    cStruct,
                                      Mv                    &rcMv,
                                      Distortion            &ruiSAD,
                                      const Mv* const       pIntegerMv2Nx2NPred )
{
  const bool bTestZeroVector          = true;
  const bool bEnableRasterSearch      = true;
  const bool bAlwaysRasterSearch      = false;  // 1: BETTER but factor 15x slower
  const bool bStarRefinementEnable    = true;   // enable either star refinement or raster refinement
  const bool bStarRefinementDiamond   = true;   // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const bool bStarRefinementStop      = false;
  const uint32_t uiStarRefinementRounds   = 2;  // star refinement stop X rounds after best match (must be >=1)
  const int  iSearchRange             = m_iSearchRange;
  const int  iSearchRangeInitial      = m_iSearchRange >> 2;
  const int  uiSearchStep             = 4;
  const int  iMVDistThresh            = 8;

  int   iStartX                 = 0;
  int   iStartY                 = 0;
  int   iDist                   = 0;

  clipMv( rcMv, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
  rcMv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
  rcMv.divideByPowerOf2(2);

  // init TZSearchStruct
  cStruct.uiBestSad = std::numeric_limits<Distortion>::max();
  cStruct.iBestX = 0;
  cStruct.iBestY = 0;

  m_cDistParam.maximumDistortionForEarlyExit = cStruct.uiBestSad;
  m_pcRdCost->setDistParam( m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, cStruct.subShiftMode );


  // set rcMv (Median predictor) as start point and as best point
  xTZSearchHelp( cStruct, rcMv.getHor(), rcMv.getVer(), 0, 0 );

  // test whether zero Mv is better start point than Median predictor
  if ( bTestZeroVector )
  {
    xTZSearchHelp( cStruct, 0, 0, 0, 0 );
  }

  SearchRange& sr = cStruct.searchRange;

  if ( pIntegerMv2Nx2NPred != 0 )
  {
    Mv integerMv2Nx2NPred = *pIntegerMv2Nx2NPred;
    integerMv2Nx2NPred.changePrecision(MV_PRECISION_INT, MV_PRECISION_INTERNAL);
    clipMv( integerMv2Nx2NPred, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
    integerMv2Nx2NPred.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
    integerMv2Nx2NPred.divideByPowerOf2(2);

    xTZSearchHelp( cStruct, integerMv2Nx2NPred.getHor(), integerMv2Nx2NPred.getVer(), 0, 0);

  }

#if AMVR_ENC_OPT
  const MvPrecision tmpIntMvPrec = (pu.cu->imv == IMV_4PEL ? MV_PRECISION_4PEL : MV_PRECISION_INT);
#endif
  for (int i = 0; i < m_uniMvListSize; i++)
  {
    BlkUniMvInfo* curMvInfo = m_uniMvList + ((m_uniMvListIdx - 1 - i + m_uniMvListMaxSize) % (m_uniMvListMaxSize));
#if AMVR_ENC_OPT
    Mv tmpCurMv = curMvInfo->uniMvs[eRefPicList][iRefIdxPred];
    tmpCurMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);
#endif

    int j = 0;
    for (; j < i; j++)
    {
      BlkUniMvInfo *prevMvInfo = m_uniMvList + ((m_uniMvListIdx - 1 - j + m_uniMvListMaxSize) % (m_uniMvListMaxSize));
#if AMVR_ENC_OPT
      Mv tmpPrevMv = prevMvInfo->uniMvs[eRefPicList][iRefIdxPred];
      tmpPrevMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);

      if (tmpCurMv == tmpPrevMv)
#else
      if (curMvInfo->uniMvs[eRefPicList][iRefIdxPred] == prevMvInfo->uniMvs[eRefPicList][iRefIdxPred])
#endif
      {
        break;
      }
    }
    if (j < i)
      continue;

    Mv cTmpMv = curMvInfo->uniMvs[eRefPicList][iRefIdxPred];
    clipMv( cTmpMv, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
#if AMVR_ENC_OPT
    cTmpMv.changePrecision(MV_PRECISION_INTERNAL, tmpIntMvPrec);
    if (tmpIntMvPrec != MV_PRECISION_INT)
    {
      cTmpMv.changePrecision(tmpIntMvPrec, MV_PRECISION_INT);
    }
#else
    cTmpMv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
#endif
    m_cDistParam.cur.buf = cStruct.piRefY + (cTmpMv.ver * cStruct.iRefStride) + cTmpMv.hor;

    Distortion uiSad = m_cDistParam.distFunc(m_cDistParam);
    uiSad += m_pcRdCost->getCostOfVectorWithPredictor(cTmpMv.hor, cTmpMv.ver, cStruct.imvShift);
    if (uiSad < cStruct.uiBestSad)
    {
      cStruct.uiBestSad = uiSad;
      cStruct.iBestX = cTmpMv.hor;
      cStruct.iBestY = cTmpMv.ver;
      m_cDistParam.maximumDistortionForEarlyExit = uiSad;
    }
  }

  {
    // set search range
    Mv currBestMv(cStruct.iBestX, cStruct.iBestY );
    currBestMv <<= 2;
    xSetSearchRange( pu, currBestMv, m_iSearchRange, sr, cStruct );
  }
  if (m_pcEncCfg->getUseHashME() && (m_currRefPicList == 0 || pu.cu->slice->getList1IdxToList0Idx(m_currRefPicIndex) < 0))
  {
    int minSize = min(pu.cu->lumaSize().width, pu.cu->lumaSize().height);
    if (minSize < 128 && minSize >= 4)
    {
      int numberOfOtherMvps = m_numHashMVStoreds[m_currRefPicList][m_currRefPicIndex];
      for (int i = 0; i < numberOfOtherMvps; i++)
      {
        xTZSearchHelp(cStruct, m_hashMVStoreds[m_currRefPicList][m_currRefPicIndex][i].getHor(), m_hashMVStoreds[m_currRefPicList][m_currRefPicIndex][i].getVer(), 0, 0);
      }
      if (numberOfOtherMvps > 0)
      {
        // write out best match
        rcMv.set(cStruct.iBestX, cStruct.iBestY);
        ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCostOfVectorWithPredictor(cStruct.iBestX, cStruct.iBestY, cStruct.imvShift);
        m_skipFracME = true;
        return;
      }
    }
  }

  // Initial search
  int iBestX = cStruct.iBestX;
  int iBestY = cStruct.iBestY;
  int iFirstSrchRngHorLeft    = ((iBestX - iSearchRangeInitial) > sr.left)   ? (iBestX - iSearchRangeInitial) : sr.left;
  int iFirstSrchRngVerTop     = ((iBestY - iSearchRangeInitial) > sr.top)    ? (iBestY - iSearchRangeInitial) : sr.top;
  int iFirstSrchRngHorRight   = ((iBestX + iSearchRangeInitial) < sr.right)  ? (iBestX + iSearchRangeInitial) : sr.right;
  int iFirstSrchRngVerBottom  = ((iBestY + iSearchRangeInitial) < sr.bottom) ? (iBestY + iSearchRangeInitial) : sr.bottom;

  for ( iStartY = iFirstSrchRngVerTop; iStartY <= iFirstSrchRngVerBottom; iStartY += uiSearchStep )
  {
    for ( iStartX = iFirstSrchRngHorLeft; iStartX <= iFirstSrchRngHorRight; iStartX += uiSearchStep )
    {
      xTZSearchHelp( cStruct, iStartX, iStartY, 0, 0 );
      xTZ8PointDiamondSearch ( cStruct, iStartX, iStartY, 1, false );
      xTZ8PointDiamondSearch ( cStruct, iStartX, iStartY, 2, false );
    }
  }

  int iMaxMVDistToPred = (abs(cStruct.iBestX - iBestX) > iMVDistThresh || abs(cStruct.iBestY - iBestY) > iMVDistThresh);

  //full search with early exit if MV is distant from predictors
  if ( bEnableRasterSearch && (iMaxMVDistToPred || bAlwaysRasterSearch) )
  {
    for ( iStartY = sr.top; iStartY <= sr.bottom; iStartY += 1 )
    {
      for ( iStartX = sr.left; iStartX <= sr.right; iStartX += 1 )
      {
        xTZSearchHelp( cStruct, iStartX, iStartY, 0, 1 );
      }
    }
  }
  //Smaller MV, refine around predictor
  else if ( bStarRefinementEnable && cStruct.uiBestDistance > 0 )
  {
    // start refinement
    while ( cStruct.uiBestDistance > 0 )
    {
      iStartX = cStruct.iBestX;
      iStartY = cStruct.iBestY;
      cStruct.uiBestDistance = 0;
      cStruct.ucPointNr = 0;
      for ( iDist = 1; iDist < iSearchRange + 1; iDist*=2 )
      {
        if ( bStarRefinementDiamond == 1 )
        {
          xTZ8PointDiamondSearch ( cStruct, iStartX, iStartY, iDist, false );
        }
        else
        {
          xTZ8PointSquareSearch  ( cStruct, iStartX, iStartY, iDist );
        }
        if ( bStarRefinementStop && (cStruct.uiBestRound >= uiStarRefinementRounds) ) // stop criterion
        {
          break;
        }
      }

      // calculate only 2 missing points instead 8 points if cStrukt.uiBestDistance == 1
      if ( cStruct.uiBestDistance == 1 )
      {
        cStruct.uiBestDistance = 0;
        if ( cStruct.ucPointNr != 0 )
        {
          xTZ2PointSearch( cStruct );
        }
      }
    }
  }

  // write out best match
  rcMv.set( cStruct.iBestX, cStruct.iBestY );
  ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCostOfVectorWithPredictor( cStruct.iBestX, cStruct.iBestY, cStruct.imvShift );
}

void InterSearch::xPatternSearchIntRefine(PredictionUnit& pu, IntTZSearchStruct&  cStruct, Mv& rcMv, Mv& rcMvPred, int& riMVPIdx, uint32_t& ruiBits, Distortion& ruiCost, const AMVPInfo& amvpInfo, double fWeight)
{

  CHECK( pu.cu->imv == 0 || pu.cu->imv == IMV_HPEL , "xPatternSearchIntRefine(): Sub-pel MV used.");
  CHECK( amvpInfo.mvCand[riMVPIdx] != rcMvPred, "xPatternSearchIntRefine(): MvPred issue.");

  const SPS &sps = *pu.cs->sps;
  m_pcRdCost->setDistParam(m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, 0, 1, m_pcEncCfg->getUseHADME() && !pu.cs->slice->getDisableSATDForRD());

  // -> set MV scale for cost calculation to QPEL (0)
  m_pcRdCost->setCostScale ( 0 );

  Distortion  uiDist, uiSATD = 0;
  Distortion  uiBestDist  = std::numeric_limits<Distortion>::max();
  // subtract old MVP costs because costs for all newly tested MVPs are added in here
#if TM_AMVP
  ruiBits -= m_auiMVPIdxCost[riMVPIdx][amvpInfo.numCand];
#else
  ruiBits -= m_auiMVPIdxCost[riMVPIdx][AMVP_MAX_NUM_CANDS];
#endif

  Mv cBestMv = rcMv;
  Mv cBaseMvd[2];
  int iBestBits = 0;
  int iBestMVPIdx = riMVPIdx;
  Mv testPos[9] = { { 0, 0}, { -1, -1},{ -1, 0},{ -1, 1},{ 0, -1},{ 0, 1},{ 1, -1},{ 1, 0},{ 1, 1} };


  cBaseMvd[0] = (rcMv - amvpInfo.mvCand[0]);
  cBaseMvd[1] = (rcMv - amvpInfo.mvCand[1]);
  CHECK( (cBaseMvd[0].getHor() & 0x03) != 0 || (cBaseMvd[0].getVer() & 0x03) != 0 , "xPatternSearchIntRefine(): AMVP cand 0 Mvd issue.");
  CHECK( (cBaseMvd[1].getHor() & 0x03) != 0 || (cBaseMvd[1].getVer() & 0x03) != 0 , "xPatternSearchIntRefine(): AMVP cand 1 Mvd issue.");

  cBaseMvd[0].roundTransPrecInternal2Amvr(pu.cu->imv);
  cBaseMvd[1].roundTransPrecInternal2Amvr(pu.cu->imv);

  // test best integer position and all 8 neighboring positions
  for (int pos = 0; pos < 9; pos ++)
  {
    Mv cTestMv[2];
    // test both AMVP candidates for each position
    for (int iMVPIdx = 0; iMVPIdx < amvpInfo.numCand; iMVPIdx++)
    {
      cTestMv[iMVPIdx] = testPos[pos];
      cTestMv[iMVPIdx].changeTransPrecAmvr2Internal(pu.cu->imv);
      cTestMv[iMVPIdx] += cBaseMvd[iMVPIdx];
      cTestMv[iMVPIdx] += amvpInfo.mvCand[iMVPIdx];

      // MCTS and IMV
      if( m_pcEncCfg->getMCTSEncConstraint() )
      {
        Mv cTestMVRestr = cTestMv[iMVPIdx];
        MCTSHelper::clipMvToArea( cTestMVRestr, pu.cu->Y(), pu.cs->picture->mctsInfo.getTileAreaIntPelRestricted( pu ), *pu.cs->sps );

        if( cTestMVRestr != cTestMv[iMVPIdx] )
        {
          // Skip this IMV pos, cause clipping affects IMV scaling
          continue;
        }
      }
      if ( iMVPIdx == 0 || cTestMv[0] != cTestMv[1])
      {
        Mv cTempMV = cTestMv[iMVPIdx];
        if( !m_pcEncCfg->getMCTSEncConstraint() )
        {
          clipMv( cTempMV, pu.cu->lumaPos(), pu.cu->lumaSize(), sps, *pu.cs->pps );
        }
        m_cDistParam.cur.buf = cStruct.piRefY  + cStruct.iRefStride * (cTempMV.getVer() >>  MV_FRACTIONAL_BITS_INTERNAL) + (cTempMV.getHor() >> MV_FRACTIONAL_BITS_INTERNAL);
        uiDist = uiSATD = (Distortion) (m_cDistParam.distFunc( m_cDistParam ) * fWeight);
      }
      else
      {
        uiDist = uiSATD;
      }

#if TM_AMVP
      int iMvBits = m_auiMVPIdxCost[iMVPIdx][amvpInfo.numCand];
#else
      int iMvBits = m_auiMVPIdxCost[iMVPIdx][AMVP_MAX_NUM_CANDS];
#endif
      Mv pred = amvpInfo.mvCand[iMVPIdx];
      pred.changeTransPrecInternal2Amvr(pu.cu->imv);
      m_pcRdCost->setPredictor( pred );
      Mv mv = cTestMv[iMVPIdx];
      mv.changeTransPrecInternal2Amvr(pu.cu->imv);
      iMvBits += m_pcRdCost->getBitsOfVectorWithPredictor( mv.getHor(), mv.getVer(), 0 );
      uiDist += m_pcRdCost->getCost(iMvBits);

      if (uiDist < uiBestDist)
      {
        uiBestDist = uiDist;
        cBestMv = cTestMv[iMVPIdx];
        iBestMVPIdx = iMVPIdx;
        iBestBits = iMvBits;
      }
    }
  }
  if( uiBestDist == std::numeric_limits<Distortion>::max() )
  {
    ruiCost = std::numeric_limits<Distortion>::max();
    return;
  }

  rcMv = cBestMv;
  rcMvPred = amvpInfo.mvCand[iBestMVPIdx];
  riMVPIdx = iBestMVPIdx;
  m_pcRdCost->setPredictor( rcMvPred );

  ruiBits += iBestBits;
  // taken from JEM 5.0
  // verify since it makes no sence to subtract Lamda*(Rmvd+Rmvpidx) from D+Lamda(Rmvd)
  // this would take the rate for the MVP idx out of the cost calculation
  // however this rate is always 1 so impact is small
  ruiCost = uiBestDist - m_pcRdCost->getCost(iBestBits) + m_pcRdCost->getCost(ruiBits);
  // taken from JEM 5.0
  // verify since it makes no sense to add rate for MVDs twicce

  return;
}

void InterSearch::xPatternSearchFracDIF(
  const PredictionUnit& pu,
  RefPicList            eRefPicList,
  int                   iRefIdx,
  IntTZSearchStruct&    cStruct,
  const Mv&             rcMvInt,
  Mv&                   rcMvHalf,
  Mv&                   rcMvQter,
  Distortion&           ruiCost
)
{

  //  Reference pattern initialization (integer scale)
  int         iOffset    = rcMvInt.getHor() + rcMvInt.getVer() * cStruct.iRefStride;
  CPelBuf cPatternRoi(cStruct.piRefY + iOffset, cStruct.iRefStride, *cStruct.pcPatternKey);
  if (m_skipFracME)
  {
    Mv baseRefMv(0, 0);
    rcMvHalf.setZero();
    m_pcRdCost->setCostScale(0);
    xExtDIFUpSamplingH(&cPatternRoi, cStruct.useAltHpelIf);
    rcMvQter = rcMvInt;   rcMvQter <<= 2;    // for mv-cost
    ruiCost = xPatternRefinement(cStruct.pcPatternKey, baseRefMv, 1, rcMvQter, !pu.cs->slice->getDisableSATDForRD());
    return;
  }


  if (cStruct.imvShift > IMV_FPEL || (m_useCompositeRef && cStruct.zeroMV))
  {
    m_pcRdCost->setDistParam(m_cDistParam, *cStruct.pcPatternKey, cStruct.piRefY + iOffset, cStruct.iRefStride, m_lumaClpRng.bd, COMPONENT_Y, 0, 1, m_pcEncCfg->getUseHADME() && !pu.cs->slice->getDisableSATDForRD());
    ruiCost = m_cDistParam.distFunc( m_cDistParam );
    ruiCost += m_pcRdCost->getCostOfVectorWithPredictor( rcMvInt.getHor(), rcMvInt.getVer(), cStruct.imvShift );
    return;
  }

  //  Half-pel refinement
  m_pcRdCost->setCostScale(1);
  xExtDIFUpSamplingH(&cPatternRoi, cStruct.useAltHpelIf);

  rcMvHalf = rcMvInt;   rcMvHalf <<= 1;    // for mv-cost
  Mv baseRefMv(0, 0);
  ruiCost = xPatternRefinement(cStruct.pcPatternKey, baseRefMv, 2, rcMvHalf, (!pu.cs->slice->getDisableSATDForRD()));

  //  quarter-pel refinement
  if (cStruct.imvShift == IMV_OFF)
  {
  m_pcRdCost->setCostScale( 0 );
  xExtDIFUpSamplingQ ( &cPatternRoi, rcMvHalf);
  baseRefMv = rcMvHalf;
  baseRefMv <<= 1;

  rcMvQter = rcMvInt;    rcMvQter <<= 1;    // for mv-cost
  rcMvQter += rcMvHalf;  rcMvQter <<= 1;
  ruiCost = xPatternRefinement(cStruct.pcPatternKey, baseRefMv, 1, rcMvQter, (!pu.cs->slice->getDisableSATDForRD()));
  }
}

Distortion InterSearch::xGetSymmetricCost( PredictionUnit& pu, PelUnitBuf& origBuf, RefPicList eCurRefPicList, const MvField& cCurMvField, MvField& cTarMvField, int bcwIdx )
{
  Distortion cost = std::numeric_limits<Distortion>::max();
  RefPicList eTarRefPicList = (RefPicList)(1 - (int)eCurRefPicList);

  // get prediction of eCurRefPicList
  PelUnitBuf predBufA = m_tmpPredStorage[eCurRefPicList].getBuf( UnitAreaRelative( *pu.cu, pu ) );
  const Picture* picRefA = pu.cu->slice->getRefPic( eCurRefPicList, cCurMvField.refIdx );
  Mv mvA = cCurMvField.mv;
  clipMv( mvA, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
  if ( (mvA.hor & 15) == 0 && (mvA.ver & 15) == 0 )
  {
    Position offset = pu.blocks[COMPONENT_Y].pos().offset( mvA.getHor() >> 4, mvA.getVer() >> 4 );
    CPelBuf pelBufA = picRefA->getRecoBuf( CompArea( COMPONENT_Y, pu.chromaFormat, offset, pu.blocks[COMPONENT_Y].size() ), false );
    predBufA.bufs[0].buf = const_cast<Pel *>(pelBufA.buf);
    predBufA.bufs[0].stride = pelBufA.stride;
    predBufA.bufs[0].width = pelBufA.width;
    predBufA.bufs[0].height = pelBufA.height;
  }
  else
  {
    xPredInterBlk( COMPONENT_Y, pu, picRefA, mvA, predBufA, false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false );
  }

  // get prediction of eTarRefPicList
  PelUnitBuf predBufB = m_tmpPredStorage[eTarRefPicList].getBuf( UnitAreaRelative( *pu.cu, pu ) );
  const Picture* picRefB = pu.cu->slice->getRefPic( eTarRefPicList, cTarMvField.refIdx );
  Mv mvB = cTarMvField.mv;
  clipMv( mvB, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
  if ( (mvB.hor & 15) == 0 && (mvB.ver & 15) == 0 )
  {
    Position offset = pu.blocks[COMPONENT_Y].pos().offset( mvB.getHor() >> 4, mvB.getVer() >> 4 );
    CPelBuf pelBufB = picRefB->getRecoBuf( CompArea( COMPONENT_Y, pu.chromaFormat, offset, pu.blocks[COMPONENT_Y].size() ), false );
    predBufB.bufs[0].buf = const_cast<Pel *>(pelBufB.buf);
    predBufB.bufs[0].stride = pelBufB.stride;
  }
  else
  {
    xPredInterBlk( COMPONENT_Y, pu, picRefB, mvB, predBufB, false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false );
  }

  PelUnitBuf bufTmp = m_tmpStorageLCU.getBuf( UnitAreaRelative( *pu.cu, pu ) );
  bufTmp.copyFrom( origBuf );
  bufTmp.removeHighFreq( predBufA, m_pcEncCfg->getClipForBiPredMeEnabled(), pu.cu->slice->clpRngs(), getBcwWeight( pu.cu->BcwIdx, eTarRefPicList ) );
  double fWeight = xGetMEDistortionWeight( pu.cu->BcwIdx, eTarRefPicList );

  // calc distortion
  DFunc distFunc = (!pu.cu->slice->getDisableSATDForRD()) ? DF_HAD : DF_SAD;
  cost = (Distortion)floor( fWeight * (double)m_pcRdCost->getDistPart( bufTmp.Y(), predBufB.Y(), pu.cs->sps->getBitDepth( CHANNEL_TYPE_LUMA ), COMPONENT_Y, distFunc ) );
  return(cost);
}

Distortion InterSearch::xSymmeticRefineMvSearch( PredictionUnit &pu, PelUnitBuf& origBuf, Mv& rcMvCurPred, Mv& rcMvTarPred
  , RefPicList eRefPicList, MvField& rCurMvField, MvField& rTarMvField, Distortion uiMinCost, int SearchPattern, int nSearchStepShift, uint32_t uiMaxSearchRounds, int bcwIdx )
{
  const Mv mvSearchOffsetCross[4] = { Mv( 0 , 1 ) , Mv( 1 , 0 ) , Mv( 0 , -1 ) , Mv( -1 ,  0 ) };
  const Mv mvSearchOffsetSquare[8] = { Mv( -1 , 1 ) , Mv( 0 , 1 ) , Mv( 1 ,  1 ) , Mv( 1 ,  0 ) , Mv( 1 , -1 ) , Mv( 0 , -1 ) , Mv( -1 , -1 ) , Mv( -1 , 0 ) };
  const Mv mvSearchOffsetDiamond[8] = { Mv( 0 , 2 ) , Mv( 1 , 1 ) , Mv( 2 ,  0 ) , Mv( 1 , -1 ) , Mv( 0 , -2 ) , Mv( -1 , -1 ) , Mv( -2 ,  0 ) , Mv( -1 , 1 ) };
  const Mv mvSearchOffsetHexagon[6] = { Mv( 2 , 0 ) , Mv( 1 , 2 ) , Mv( -1 ,  2 ) , Mv( -2 ,  0 ) , Mv( -1 , -2 ) , Mv( 1 , -2 ) };

  int nDirectStart = 0, nDirectEnd = 0, nDirectRounding = 0, nDirectMask = 0;
  const Mv * pSearchOffset;
  if ( SearchPattern == 0 )
  {
    nDirectEnd = 3;
    nDirectRounding = 4;
    nDirectMask = 0x03;
    pSearchOffset = mvSearchOffsetCross;
  }
  else if ( SearchPattern == 1 )
  {
    nDirectEnd = 7;
    nDirectRounding = 8;
    nDirectMask = 0x07;
    pSearchOffset = mvSearchOffsetSquare;
  }
  else if ( SearchPattern == 2 )
  {
    nDirectEnd = 7;
    nDirectRounding = 8;
    nDirectMask = 0x07;
    pSearchOffset = mvSearchOffsetDiamond;
  }
  else if ( SearchPattern == 3 )
  {
    nDirectEnd = 5;
    pSearchOffset = mvSearchOffsetHexagon;
  }
  else
  {
    THROW( "Invalid search pattern" );
  }

  int nBestDirect;
  for ( uint32_t uiRound = 0; uiRound < uiMaxSearchRounds; uiRound++ )
  {
    nBestDirect = -1;
    MvField mvCurCenter = rCurMvField;
    for ( int nIdx = nDirectStart; nIdx <= nDirectEnd; nIdx++ )
    {
      int nDirect;
      if ( SearchPattern == 3 )
      {
        nDirect = nIdx < 0 ? nIdx + 6 : nIdx >= 6 ? nIdx - 6 : nIdx;
      }
      else
      {
        nDirect = (nIdx + nDirectRounding) & nDirectMask;
      }

      Mv mvOffset = pSearchOffset[nDirect];
      mvOffset <<= nSearchStepShift;
      MvField mvCand = mvCurCenter, mvPair;
      mvCand.mv += mvOffset;

      if( m_pcEncCfg->getMCTSEncConstraint() )
      {
        if( !( MCTSHelper::checkMvForMCTSConstraint( pu, mvCand.mv ) ) )
          continue; // Skip this this pos
      }
      // get MVD cost
      Mv pred = rcMvCurPred;
      pred.changeTransPrecInternal2Amvr(pu.cu->imv);
      m_pcRdCost->setPredictor( pred );
      m_pcRdCost->setCostScale( 0 );
      Mv mv = mvCand.mv;
      mv.changeTransPrecInternal2Amvr(pu.cu->imv);
      uint32_t uiMvBits = m_pcRdCost->getBitsOfVectorWithPredictor( mv.getHor(), mv.getVer(), 0 );
      Distortion uiCost = m_pcRdCost->getCost( uiMvBits );

      // get MVD pair and set target MV
      mvPair.refIdx = rTarMvField.refIdx;
      mvPair.mv.set( rcMvTarPred.hor - (mvCand.mv.hor - rcMvCurPred.hor), rcMvTarPred.ver - (mvCand.mv.ver - rcMvCurPred.ver) );
      if( m_pcEncCfg->getMCTSEncConstraint() )
      {
        if( !( MCTSHelper::checkMvForMCTSConstraint( pu, mvPair.mv ) ) )
          continue; // Skip this this pos
      }
      uiCost += xGetSymmetricCost( pu, origBuf, eRefPicList, mvCand, mvPair, bcwIdx );
      if ( uiCost < uiMinCost )
      {
        uiMinCost = uiCost;
        rCurMvField = mvCand;
        rTarMvField = mvPair;
        nBestDirect = nDirect;
      }
    }

    if ( nBestDirect == -1 )
    {
      break;
    }
    int nStep = 1;
    if ( SearchPattern == 1 || SearchPattern == 2 )
    {
      nStep = 2 - (nBestDirect & 0x01);
    }
    nDirectStart = nBestDirect - nStep;
    nDirectEnd = nBestDirect + nStep;
  }

  return(uiMinCost);
}


void InterSearch::xSymmetricMotionEstimation( PredictionUnit& pu, PelUnitBuf& origBuf, Mv& rcMvCurPred, Mv& rcMvTarPred, RefPicList eRefPicList, MvField& rCurMvField, MvField& rTarMvField, Distortion& ruiCost, int bcwIdx )
{
  // Refine Search
  int nSearchStepShift = MV_FRACTIONAL_BITS_DIFF;
  int nDiamondRound = 8;
  int nCrossRound = 1;

  nSearchStepShift += pu.cu->imv == IMV_HPEL ? 1 : (pu.cu->imv << 1);
  nDiamondRound >>= pu.cu->imv;

  ruiCost = xSymmeticRefineMvSearch( pu, origBuf, rcMvCurPred, rcMvTarPred, eRefPicList, rCurMvField, rTarMvField, ruiCost, 2, nSearchStepShift, nDiamondRound, bcwIdx );
  ruiCost = xSymmeticRefineMvSearch( pu, origBuf, rcMvCurPred, rcMvTarPred, eRefPicList, rCurMvField, rTarMvField, ruiCost, 0, nSearchStepShift, nCrossRound, bcwIdx );
}

void InterSearch::xPredAffineInterSearch( PredictionUnit&       pu,
                                          PelUnitBuf&           origBuf,
                                          int                   puIdx,
                                          uint32_t&                 lastMode,
                                          Distortion&           affineCost,
                                          Mv                    hevcMv[2][33]
                                        , Mv                    mvAffine4Para[2][33][3]
                                        , int                   refIdx4Para[2]
                                        , uint8_t               bcwIdx
                                        , bool                  enforceBcwPred
                                        , uint32_t              bcwIdxBits
                                         )
{
  const Slice &slice = *pu.cu->slice;

  affineCost = std::numeric_limits<Distortion>::max();

  Mv        cMvZero;
  Mv        aacMv[2][3];
  Mv        cMvBi[2][3];
  Mv        cMvTemp[2][33][3];

  int       iNumPredDir = slice.isInterP() ? 1 : 2;

  int mvNum = 2;
  mvNum = pu.cu->affineType ? 3 : 2;

  // Mvp
  Mv        cMvPred[2][33][3];
  Mv        cMvPredBi[2][33][3];
  int       aaiMvpIdxBi[2][33];
  int       aaiMvpIdx[2][33];
  int       aaiMvpNum[2][33];

  AffineAMVPInfo aacAffineAMVPInfo[2][33];
  AffineAMVPInfo affiAMVPInfoTemp[2];

  int           iRefIdx[2]={0,0}; // If un-initialized, may cause SEGV in bi-directional prediction iterative stage.
  int           iRefIdxBi[2];

  uint32_t          uiMbBits[3] = {1, 1, 0};

  int           iRefStart, iRefEnd;

  int           bestBiPRefIdxL1 = 0;
  int           bestBiPMvpL1 = 0;
  Distortion biPDistTemp = std::numeric_limits<Distortion>::max();

  Distortion    uiCost[2] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max() };
  Distortion    uiCostBi  = std::numeric_limits<Distortion>::max();
  Distortion    uiCostTemp;

  uint32_t          uiBits[3] = { 0 };
  uint32_t          uiBitsTemp;
  Distortion    bestBiPDist = std::numeric_limits<Distortion>::max();

  Distortion    uiCostTempL0[MAX_NUM_REF];
  for (int iNumRef=0; iNumRef < MAX_NUM_REF; iNumRef++)
  {
    uiCostTempL0[iNumRef] = std::numeric_limits<Distortion>::max();
  }
  uint32_t uiBitsTempL0[MAX_NUM_REF];

  Mv            mvValidList1[4];
  int           refIdxValidList1 = 0;
  uint32_t          bitsValidList1 = MAX_UINT;
  Distortion costValidList1 = std::numeric_limits<Distortion>::max();
  Mv            mvHevc[3];
  const bool affineAmvrEnabled = pu.cu->slice->getSPS()->getAffineAmvrEnabledFlag();
  int tryBipred = 0;
  WPScalingParam *wp0;
  WPScalingParam *wp1;
  xGetBlkBits( slice.isInterP(), puIdx, lastMode, uiMbBits);

  pu.cu->affine = true;
  pu.mergeFlag = false;
  pu.regularMergeFlag = false;
  if( bcwIdx != BCW_DEFAULT )
  {
    pu.cu->BcwIdx = bcwIdx;
  }
#if MULTI_HYP_PRED
  const bool saveMeResultsForMHP = pu.cs->sps->getUseInterMultiHyp()
    && pu.cu->imv == 0
    && bcwIdx != BCW_DEFAULT
    && (pu.Y().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE && std::min(pu.Y().width, pu.Y().height) >= MULTI_HYP_PRED_RESTRICT_MIN_WH)
    ;
#endif

  // Uni-directional prediction
  for ( int iRefList = 0; iRefList < iNumPredDir; iRefList++ )
  {
    RefPicList  eRefPicList = ( iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0 );
    pu.interDir = ( iRefList ? 2 : 1 );
    for (int iRefIdxTemp = 0; iRefIdxTemp < slice.getNumRefIdx(eRefPicList); iRefIdxTemp++)
    {
      // Get RefIdx bits
      uiBitsTemp = uiMbBits[iRefList];
      if ( slice.getNumRefIdx(eRefPicList) > 1 )
      {
        uiBitsTemp += iRefIdxTemp+1;
        if ( iRefIdxTemp == slice.getNumRefIdx(eRefPicList)-1 )
        {
          uiBitsTemp--;
        }
      }

      // Do Affine AMVP
      xEstimateAffineAMVP( pu, affiAMVPInfoTemp[eRefPicList], origBuf, eRefPicList, iRefIdxTemp, cMvPred[iRefList][iRefIdxTemp], &biPDistTemp );
      if ( affineAmvrEnabled )
      {
        biPDistTemp += m_pcRdCost->getCost( xCalcAffineMVBits( pu, cMvPred[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp] ) );
      }
      aaiMvpIdx[iRefList][iRefIdxTemp] = pu.mvpIdx[eRefPicList];
      aaiMvpNum[iRefList][iRefIdxTemp] = pu.mvpNum[eRefPicList];;
      if ( pu.cu->affineType == AFFINEMODEL_6PARAM && refIdx4Para[iRefList] != iRefIdxTemp )
      {
        xCopyAffineAMVPInfo( affiAMVPInfoTemp[eRefPicList], aacAffineAMVPInfo[iRefList][iRefIdxTemp] );
        continue;
      }

      // set hevc ME result as start search position when it is best than mvp
      for ( int i=0; i<3; i++ )
      {
        mvHevc[i] = hevcMv[iRefList][iRefIdxTemp];
        mvHevc[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
      }
      PelUnitBuf predBuf = m_tmpStorageLCU.getBuf( UnitAreaRelative(*pu.cu, pu) );

      Distortion uiCandCost = xGetAffineTemplateCost(pu, origBuf, predBuf, mvHevc, aaiMvpIdx[iRefList][iRefIdxTemp],
                                                     AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdxTemp);

      if ( affineAmvrEnabled )
      {
        uiCandCost += m_pcRdCost->getCost( xCalcAffineMVBits( pu, mvHevc, cMvPred[iRefList][iRefIdxTemp] ) );
      }

      //check stored affine motion
      bool affine4Para    = pu.cu->affineType == AFFINEMODEL_4PARAM;
      bool savedParaAvail = pu.cu->imv && ( ( m_affineMotion.affine4ParaRefIdx[iRefList] == iRefIdxTemp && affine4Para && m_affineMotion.affine4ParaAvail ) ||
                                            ( m_affineMotion.affine6ParaRefIdx[iRefList] == iRefIdxTemp && !affine4Para && m_affineMotion.affine6ParaAvail ) );

      if ( savedParaAvail )
      {
        Mv mvFour[3];
        for ( int i = 0; i < mvNum; i++ )
        {
          mvFour[i] = affine4Para ? m_affineMotion.acMvAffine4Para[iRefList][i] : m_affineMotion.acMvAffine6Para[iRefList][i];
          mvFour[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
        }

        Distortion candCostInherit = xGetAffineTemplateCost( pu, origBuf, predBuf, mvFour, aaiMvpIdx[iRefList][iRefIdxTemp], AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdxTemp );
        candCostInherit += m_pcRdCost->getCost( xCalcAffineMVBits( pu, mvFour, cMvPred[iRefList][iRefIdxTemp] ) );

        if ( candCostInherit < uiCandCost )
        {
          uiCandCost = candCostInherit;
          memcpy( mvHevc, mvFour, 3 * sizeof( Mv ) );
        }
      }

      if (pu.cu->affineType == AFFINEMODEL_4PARAM && m_affMVListSize
        && (!pu.cu->cs->sps->getUseBcw() || bcwIdx == BCW_DEFAULT)
        )
      {
        int shift = MAX_CU_DEPTH;
        for (int i = 0; i < m_affMVListSize; i++)
        {
          AffineMVInfo *mvInfo = m_affMVList + ((m_affMVListIdx - i - 1 + m_affMVListMaxSize) % (m_affMVListMaxSize));
          //check;
          int j = 0;
          for (; j < i; j++)
          {
            AffineMVInfo *prevMvInfo = m_affMVList + ((m_affMVListIdx - j - 1 + m_affMVListMaxSize) % (m_affMVListMaxSize));
            if ((mvInfo->affMVs[iRefList][iRefIdxTemp][0] == prevMvInfo->affMVs[iRefList][iRefIdxTemp][0]) &&
              (mvInfo->affMVs[iRefList][iRefIdxTemp][1] == prevMvInfo->affMVs[iRefList][iRefIdxTemp][1])
              && (mvInfo->x == prevMvInfo->x) && (mvInfo->y == prevMvInfo->y)
              && (mvInfo->w == prevMvInfo->w)
              )
            {
              break;
            }
          }
          if (j < i)
            continue;

          Mv mvTmp[3], *nbMv = mvInfo->affMVs[iRefList][iRefIdxTemp];
          int vx, vy;
          int dMvHorX, dMvHorY, dMvVerX, dMvVerY;
          int mvScaleHor = nbMv[0].getHor() << shift;
          int mvScaleVer = nbMv[0].getVer() << shift;
          Mv dMv = nbMv[1] - nbMv[0];
          dMvHorX = dMv.getHor() << (shift - floorLog2(mvInfo->w));
          dMvHorY = dMv.getVer() << (shift - floorLog2(mvInfo->w));
          dMvVerX = -dMvHorY;
          dMvVerY = dMvHorX;
          vx = mvScaleHor + dMvHorX * (pu.Y().x - mvInfo->x) + dMvVerX * (pu.Y().y - mvInfo->y);
          vy = mvScaleVer + dMvHorY * (pu.Y().x - mvInfo->x) + dMvVerY * (pu.Y().y - mvInfo->y);
          roundAffineMv(vx, vy, shift);
          mvTmp[0] = Mv(vx, vy);
          mvTmp[0].clipToStorageBitDepth();
          clipMv( mvTmp[0], pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
          mvTmp[0].roundAffinePrecInternal2Amvr(pu.cu->imv);
          vx = mvScaleHor + dMvHorX * (pu.Y().x + pu.Y().width - mvInfo->x) + dMvVerX * (pu.Y().y - mvInfo->y);
          vy = mvScaleVer + dMvHorY * (pu.Y().x + pu.Y().width - mvInfo->x) + dMvVerY * (pu.Y().y - mvInfo->y);
          roundAffineMv(vx, vy, shift);
          mvTmp[1] = Mv(vx, vy);
          mvTmp[1].clipToStorageBitDepth();
          clipMv( mvTmp[1], pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
          mvTmp[0].roundAffinePrecInternal2Amvr(pu.cu->imv);
          mvTmp[1].roundAffinePrecInternal2Amvr(pu.cu->imv);
          Distortion tmpCost = xGetAffineTemplateCost(pu, origBuf, predBuf, mvTmp, aaiMvpIdx[iRefList][iRefIdxTemp], AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdxTemp);
          if ( affineAmvrEnabled )
          {
            tmpCost += m_pcRdCost->getCost( xCalcAffineMVBits( pu, mvTmp, cMvPred[iRefList][iRefIdxTemp] ) );
          }
          if (tmpCost < uiCandCost)
          {
            uiCandCost = tmpCost;
            std::memcpy(mvHevc, mvTmp, 3 * sizeof(Mv));
          }
        }
      }
      if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
      {
        Mv mvFour[3];
        mvFour[0] = mvAffine4Para[iRefList][iRefIdxTemp][0];
        mvFour[1] = mvAffine4Para[iRefList][iRefIdxTemp][1];
        mvAffine4Para[iRefList][iRefIdxTemp][0].roundAffinePrecInternal2Amvr(pu.cu->imv);
        mvAffine4Para[iRefList][iRefIdxTemp][1].roundAffinePrecInternal2Amvr(pu.cu->imv);

        int shift = MAX_CU_DEPTH;
        int vx2 = (mvFour[0].getHor() << shift) - ((mvFour[1].getVer() - mvFour[0].getVer()) << (shift + floorLog2(pu.lheight()) - floorLog2(pu.lwidth())));
        int vy2 = (mvFour[0].getVer() << shift) + ((mvFour[1].getHor() - mvFour[0].getHor()) << (shift + floorLog2(pu.lheight()) - floorLog2(pu.lwidth())));
        int offset = (1 << (shift - 1));
        vx2 = (vx2 + offset - (vx2 >= 0)) >> shift;
        vy2 = (vy2 + offset - (vy2 >= 0)) >> shift;
        mvFour[2].hor = vx2;
        mvFour[2].ver = vy2;
        mvFour[2].clipToStorageBitDepth();
        mvFour[0].roundAffinePrecInternal2Amvr(pu.cu->imv);
        mvFour[1].roundAffinePrecInternal2Amvr(pu.cu->imv);
        mvFour[2].roundAffinePrecInternal2Amvr(pu.cu->imv);
        Distortion uiCandCostInherit = xGetAffineTemplateCost( pu, origBuf, predBuf, mvFour, aaiMvpIdx[iRefList][iRefIdxTemp], AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdxTemp );
        if ( affineAmvrEnabled )
        {
          uiCandCostInherit += m_pcRdCost->getCost( xCalcAffineMVBits( pu, mvFour, cMvPred[iRefList][iRefIdxTemp] ) );
        }
        if ( uiCandCostInherit < uiCandCost )
        {
          uiCandCost = uiCandCostInherit;
          for ( int i = 0; i < 3; i++ )
          {
            mvHevc[i] = mvFour[i];
          }
        }
      }

      if ( uiCandCost < biPDistTemp )
      {
        ::memcpy( cMvTemp[iRefList][iRefIdxTemp], mvHevc, sizeof(Mv)*3 );
      }
      else
      {
        ::memcpy( cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], sizeof(Mv)*3 );
      }

      // GPB list 1, save the best MvpIdx, RefIdx and Cost
      if ( slice.getPicHeader()->getMvdL1ZeroFlag() && iRefList==1 && biPDistTemp < bestBiPDist )
      {
        bestBiPDist = biPDistTemp;
        bestBiPMvpL1 = aaiMvpIdx[iRefList][iRefIdxTemp];
        bestBiPRefIdxL1 = iRefIdxTemp;
      }

      // Update bits
      uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdx[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];

      if ( m_pcEncCfg->getFastMEForGenBLowDelayEnabled() && iRefList == 1 )   // list 1
      {
        if ( slice.getList1IdxToList0Idx( iRefIdxTemp ) >= 0 && (pu.cu->affineType != AFFINEMODEL_6PARAM || slice.getList1IdxToList0Idx( iRefIdxTemp ) == refIdx4Para[0]) )
        {
          int iList1ToList0Idx = slice.getList1IdxToList0Idx( iRefIdxTemp );
          ::memcpy( cMvTemp[1][iRefIdxTemp], cMvTemp[0][iList1ToList0Idx], sizeof(Mv)*3 );
          uiCostTemp = uiCostTempL0[iList1ToList0Idx];

          uiCostTemp -= m_pcRdCost->getCost( uiBitsTempL0[iList1ToList0Idx] );
          uiBitsTemp += xCalcAffineMVBits( pu, cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp] );
          /*calculate the correct cost*/
          uiCostTemp += m_pcRdCost->getCost( uiBitsTemp );
          DTRACE( g_trace_ctx, D_COMMON, " (%d) uiCostTemp=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiCostTemp );
        }
        else
        {
          xAffineMotionEstimation( pu, origBuf, eRefPicList, cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
                                   , aaiMvpIdx[iRefList][iRefIdxTemp], affiAMVPInfoTemp[eRefPicList]
          );
        }
      }
      else
      {
        xAffineMotionEstimation( pu, origBuf, eRefPicList, cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp
                                 , aaiMvpIdx[iRefList][iRefIdxTemp], affiAMVPInfoTemp[eRefPicList]
        );
      }
      if(pu.cu->cs->sps->getUseBcw() && pu.cu->BcwIdx == BCW_DEFAULT && pu.cu->slice->isInterB())
      {
        m_uniMotions.setReadModeAffine(true, (uint8_t)iRefList, (uint8_t)iRefIdxTemp, pu.cu->affineType);
        m_uniMotions.copyAffineMvFrom(cMvTemp[iRefList][iRefIdxTemp], uiCostTemp - m_pcRdCost->getCost(uiBitsTemp), (uint8_t)iRefList, (uint8_t)iRefIdxTemp, pu.cu->affineType
                                      , aaiMvpIdx[iRefList][iRefIdxTemp]
        );
      }
      // Set best AMVP Index
      xCopyAffineAMVPInfo( affiAMVPInfoTemp[eRefPicList], aacAffineAMVPInfo[iRefList][iRefIdxTemp] );
      if ( pu.cu->imv != 2 || !m_pcEncCfg->getUseAffineAmvrEncOpt() )
      xCheckBestAffineMVP( pu, affiAMVPInfoTemp[eRefPicList], eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp );

      if ( iRefList == 0 )
      {
        uiCostTempL0[iRefIdxTemp] = uiCostTemp;
        uiBitsTempL0[iRefIdxTemp] = uiBitsTemp;
      }
      DTRACE( g_trace_ctx, D_COMMON, " (%d) uiCostTemp=%d, uiCost[iRefList]=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiCostTemp, uiCost[iRefList] );
      if ( uiCostTemp < uiCost[iRefList] )
      {
        uiCost[iRefList] = uiCostTemp;
        uiBits[iRefList] = uiBitsTemp; // storing for bi-prediction

        // set best motion
        ::memcpy( aacMv[iRefList], cMvTemp[iRefList][iRefIdxTemp], sizeof(Mv) * 3 );
        iRefIdx[iRefList] = iRefIdxTemp;
      }
#if JVET_Z0054_BLK_REF_PIC_REORDER
      if (pu.cs->sps->getUseARL() && iRefList == 1 && slice.getList1IdxToList0Idx(iRefIdxTemp) >= 0)
      {
        uiCostTemp = MAX_UINT;
      }
#endif

      if ( iRefList == 1 && uiCostTemp < costValidList1 && slice.getList1IdxToList0Idx( iRefIdxTemp ) < 0 )
      {
        costValidList1 = uiCostTemp;
        bitsValidList1 = uiBitsTemp;

        // set motion
        memcpy( mvValidList1, cMvTemp[iRefList][iRefIdxTemp], sizeof(Mv)*3 );
        refIdxValidList1 = iRefIdxTemp;
      }
    } // End refIdx loop
  } // end Uni-prediction

  if ( pu.cu->affineType == AFFINEMODEL_4PARAM )
  {
    ::memcpy( mvAffine4Para, cMvTemp, sizeof( cMvTemp ) );
    if ( pu.cu->imv == 0 && ( !pu.cu->cs->sps->getUseBcw() || bcwIdx == BCW_DEFAULT ) )
    {
      AffineMVInfo *affMVInfo = m_affMVList + m_affMVListIdx;

      //check;
      int j = 0;
      for (; j < m_affMVListSize; j++)
      {
        AffineMVInfo *prevMvInfo = m_affMVList + ((m_affMVListIdx - j - 1 + m_affMVListMaxSize) % (m_affMVListMaxSize));
        if ((pu.Y().x == prevMvInfo->x) && (pu.Y().y == prevMvInfo->y) && (pu.Y().width == prevMvInfo->w) && (pu.Y().height == prevMvInfo->h))
        {
          break;
        }
      }
      if (j < m_affMVListSize)
        affMVInfo = m_affMVList + ((m_affMVListIdx - j - 1 + m_affMVListMaxSize) % (m_affMVListMaxSize));

      ::memcpy(affMVInfo->affMVs, cMvTemp, sizeof(cMvTemp));

      if (j == m_affMVListSize)
      {
        affMVInfo->x = pu.Y().x;
        affMVInfo->y = pu.Y().y;
        affMVInfo->w = pu.Y().width;
        affMVInfo->h = pu.Y().height;
        m_affMVListSize = std::min(m_affMVListSize + 1, m_affMVListMaxSize);
        m_affMVListIdx = (m_affMVListIdx + 1) % (m_affMVListMaxSize);
      }
    }
  }

  // Bi-directional prediction
  if ( slice.isInterB() && !PU::isBipredRestriction(pu) 
#if AFFINE_ENC_OPT || MULTI_HYP_PRED 
    // In case refIdx4Para[i] is NOT_VALID, uiMotBits[i] would be undefined since list i will not be searched in 6-para model.
    // Therefore, the undefined bits would be stored in MHP candidates.
    && !(pu.cu->affineType == AFFINEMODEL_6PARAM && (refIdx4Para[0] == NOT_VALID || refIdx4Para[1] == NOT_VALID))
#endif
#if INTER_LIC
    && !pu.cu->LICFlag
#endif
    )
  {
    tryBipred = 1;
    pu.interDir = 3;
    m_isBi = true;
    // Set as best list0 and list1
    iRefIdxBi[0] = iRefIdx[0];
    iRefIdxBi[1] = iRefIdx[1];

    ::memcpy( cMvBi,       aacMv,     sizeof(aacMv)     );
    ::memcpy( cMvPredBi,   cMvPred,   sizeof(cMvPred)   );
    ::memcpy( aaiMvpIdxBi, aaiMvpIdx, sizeof(aaiMvpIdx) );

    uint32_t uiMotBits[2];
    bool doBiPred = true;

    if ( slice.getPicHeader()->getMvdL1ZeroFlag() ) // GPB, list 1 only use Mvp
    {
      xCopyAffineAMVPInfo( aacAffineAMVPInfo[1][bestBiPRefIdxL1], affiAMVPInfoTemp[REF_PIC_LIST_1] );
      pu.mvpIdx[REF_PIC_LIST_1] = bestBiPMvpL1;
      aaiMvpIdxBi[1][bestBiPRefIdxL1] = bestBiPMvpL1;

      // Set Mv for list1
      Mv pcMvTemp[3] = { affiAMVPInfoTemp[REF_PIC_LIST_1].mvCandLT[bestBiPMvpL1],
                         affiAMVPInfoTemp[REF_PIC_LIST_1].mvCandRT[bestBiPMvpL1],
                         affiAMVPInfoTemp[REF_PIC_LIST_1].mvCandLB[bestBiPMvpL1] };
      ::memcpy( cMvPredBi[1][bestBiPRefIdxL1], pcMvTemp, sizeof(Mv)*3 );
      ::memcpy( cMvBi[1],                      pcMvTemp, sizeof(Mv)*3 );
      ::memcpy( cMvTemp[1][bestBiPRefIdxL1],   pcMvTemp, sizeof(Mv)*3 );
      iRefIdxBi[1] = bestBiPRefIdxL1;

      if( m_pcEncCfg->getMCTSEncConstraint() )
      {
        Area curTileAreaRestricted;
        curTileAreaRestricted = pu.cs->picture->mctsInfo.getTileAreaSubPelRestricted( pu );
        for( int i = 0; i < mvNum; i++ )
        {
          Mv restrictedMv = pcMvTemp[i];
          MCTSHelper::clipMvToArea( restrictedMv, pu.cu->Y(), curTileAreaRestricted, *pu.cs->sps );

          // If sub-pel filter samples are not inside of allowed area
          if( restrictedMv != pcMvTemp[i] )
          {
            uiCostBi = std::numeric_limits<Distortion>::max();
            doBiPred = false;
          }
        }
      }
      // Get list1 prediction block
      PU::setAllAffineMv( pu, cMvBi[1][0], cMvBi[1][1], cMvBi[1][2], REF_PIC_LIST_1);
      pu.refIdx[REF_PIC_LIST_1] = iRefIdxBi[1];

      PelUnitBuf predBufTmp = m_tmpPredStorage[REF_PIC_LIST_1].getBuf( UnitAreaRelative(*pu.cu, pu) );
      motionCompensation( pu, predBufTmp, REF_PIC_LIST_1 );

      // Update bits
      uiMotBits[0] = uiBits[0] - uiMbBits[0];
      uiMotBits[1] = uiMbBits[1];

      if( slice.getNumRefIdx(REF_PIC_LIST_1) > 1 )
      {
        uiMotBits[1] += bestBiPRefIdxL1+1;
        if( bestBiPRefIdxL1 == slice.getNumRefIdx(REF_PIC_LIST_1)-1 )
        {
          uiMotBits[1]--;
        }
      }
      uiMotBits[1] += m_auiMVPIdxCost[aaiMvpIdxBi[1][bestBiPRefIdxL1]][AMVP_MAX_NUM_CANDS];
      uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];
    }
    else
    {
      uiMotBits[0] = uiBits[0] - uiMbBits[0];
      uiMotBits[1] = uiBits[1] - uiMbBits[1];
      uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];
    }

    if( doBiPred )
    {
    // 4-times iteration (default)
    int iNumIter = 4;
    // fast encoder setting or GPB: only one iteration
    if ( m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE2 || slice.getPicHeader()->getMvdL1ZeroFlag() )
    {
      iNumIter = 1;
    }

    for ( int iIter = 0; iIter < iNumIter; iIter++ )
    {
      // Set RefList
      int iRefList = iIter % 2;
      if ( m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE2 )
      {
        if( uiCost[0] <= uiCost[1] )
        {
          iRefList = 1;
        }
        else
        {
          iRefList = 0;
        }
        if( bcwIdx != BCW_DEFAULT )
        {
          iRefList = ( abs( getBcwWeight( bcwIdx, REF_PIC_LIST_0 ) ) > abs( getBcwWeight( bcwIdx, REF_PIC_LIST_1 ) ) ? 1 : 0 );
        }
      }
      else if ( iIter == 0 )
      {
        iRefList = 0;
      }

      // First iterate, get prediction block of opposite direction
      if( iIter == 0 && !slice.getPicHeader()->getMvdL1ZeroFlag() )
      {
        PU::setAllAffineMv( pu, aacMv[1-iRefList][0], aacMv[1-iRefList][1], aacMv[1-iRefList][2], RefPicList(1-iRefList));
        pu.refIdx[1-iRefList] = iRefIdx[1-iRefList];

        PelUnitBuf predBufTmp = m_tmpPredStorage[1 - iRefList].getBuf( UnitAreaRelative(*pu.cu, pu) );
        motionCompensation( pu, predBufTmp, RefPicList(1 - iRefList) );
      }

      RefPicList eRefPicList = ( iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0 );

      if ( slice.getPicHeader()->getMvdL1ZeroFlag() ) // GPB, fix List 1, search List 0
      {
        iRefList = 0;
        eRefPicList = REF_PIC_LIST_0;
      }

      bool bChanged = false;

      iRefStart = 0;
      iRefEnd   = slice.getNumRefIdx(eRefPicList) - 1;
      for ( int iRefIdxTemp = iRefStart; iRefIdxTemp <= iRefEnd; iRefIdxTemp++ )
      {
#if JVET_Z0054_BLK_REF_PIC_REORDER
        if (pu.cs->sps->getUseARL())
        {
          int refIdxTemp[2];
          refIdxTemp[iRefList] = iRefIdxTemp;
          refIdxTemp[1 - iRefList] = iRefIdxBi[1 - iRefList];
          if (pu.cu->slice->getRefPicPairIdx(refIdxTemp[0], refIdxTemp[1]) < 0)
          {
            continue;
          }
        }
#endif
        if ( pu.cu->affineType == AFFINEMODEL_6PARAM && refIdx4Para[iRefList] != iRefIdxTemp )
        {
          continue;
        }
        if(m_pcEncCfg->getUseBcwFast() && (bcwIdx != BCW_DEFAULT)
          && (pu.cu->slice->getRefPic(eRefPicList, iRefIdxTemp)->getPOC() == pu.cu->slice->getRefPic(RefPicList(1 - iRefList), pu.refIdx[1 - iRefList])->getPOC())
          && (pu.cu->affineType == AFFINEMODEL_4PARAM && pu.cu->slice->getTLayer()>1))
        {
          continue;
        }
        // update bits
        uiBitsTemp = uiMbBits[2] + uiMotBits[1-iRefList];
        uiBitsTemp += ((pu.cu->slice->getSPS()->getUseBcw() == true) ? bcwIdxBits : 0);
        if( slice.getNumRefIdx(eRefPicList) > 1 )
        {
          uiBitsTemp += iRefIdxTemp+1;
          if ( iRefIdxTemp == slice.getNumRefIdx(eRefPicList)-1 )
          {
            uiBitsTemp--;
          }
        }
        uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdxBi[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];
        // call Affine ME
        xAffineMotionEstimation( pu, origBuf, eRefPicList, cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp,
                                 aaiMvpIdxBi[iRefList][iRefIdxTemp], aacAffineAMVPInfo[iRefList][iRefIdxTemp],
          true );
        xCopyAffineAMVPInfo( aacAffineAMVPInfo[iRefList][iRefIdxTemp], affiAMVPInfoTemp[eRefPicList] );
        if ( pu.cu->imv != 2 || !m_pcEncCfg->getUseAffineAmvrEncOpt() )
        xCheckBestAffineMVP( pu, affiAMVPInfoTemp[eRefPicList], eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPredBi[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp );

#if MULTI_HYP_PRED
        if(saveMeResultsForMHP)
        {
          // Affine bi
          MEResult biPredResult;
          biPredResult.cu = *pu.cu;
          biPredResult.cu.smvdMode = 0;
          biPredResult.pu = pu;
          biPredResult.cost = uiCostTemp;
          biPredResult.bits = uiBitsTemp;

          biPredResult.pu.mv[REF_PIC_LIST_0] = Mv();
          biPredResult.pu.mv[REF_PIC_LIST_1] = Mv();
          biPredResult.pu.mvd[REF_PIC_LIST_0] = cMvZero;
          biPredResult.pu.mvd[REF_PIC_LIST_1] = cMvZero;

          for (int i = 0; i < 3; ++i)
          {
            biPredResult.pu.mvAffi[iRefList][i] = cMvTemp[iRefList][iRefIdxTemp][i];
            biPredResult.pu.mvAffi[1 - iRefList][i] = cMvBi[1 - iRefList][i];
            biPredResult.pu.mvAffi[0][i].roundAffinePrecInternal2Amvr(pu.cu->imv);
            biPredResult.pu.mvAffi[1][i].roundAffinePrecInternal2Amvr(pu.cu->imv);
          }

          biPredResult.pu.refIdx[iRefList] = iRefIdxTemp;
          biPredResult.pu.refIdx[1 - iRefList] = iRefIdxBi[1 - iRefList];

          for (int verIdx = 0; verIdx < mvNum; verIdx++)
          {
            biPredResult.pu.mvdAffi[iRefList][verIdx] = cMvTemp[iRefList][iRefIdxTemp][verIdx] - cMvPredBi[iRefList][iRefIdxTemp][verIdx];
            biPredResult.pu.mvdAffi[1 - iRefList][verIdx] = cMvBi[1 - iRefList][verIdx] - cMvPredBi[1 - iRefList][iRefIdxBi[1 - iRefList]][verIdx];
            if (verIdx != 0)
            {
              biPredResult.pu.mvdAffi[0][verIdx] = biPredResult.pu.mvdAffi[0][verIdx] - biPredResult.pu.mvdAffi[0][0];
              biPredResult.pu.mvdAffi[1][verIdx] = biPredResult.pu.mvdAffi[1][verIdx] - biPredResult.pu.mvdAffi[1][0];
            }
          }

          biPredResult.pu.interDir = 3;

          biPredResult.pu.mvpIdx[iRefList] = aaiMvpIdxBi[iRefList][iRefIdxTemp];
          biPredResult.pu.mvpIdx[1 - iRefList] = aaiMvpIdxBi[1 - iRefList][iRefIdxBi[1 - iRefList]];
          biPredResult.pu.mvpNum[iRefList] = aaiMvpNum[iRefList][iRefIdxTemp];
          biPredResult.pu.mvpNum[1 - iRefList] = aaiMvpNum[1 - iRefList][iRefIdxBi[1 - iRefList]];

          pu.cs->m_meResults.push_back(biPredResult);
        }
#endif
        if ( uiCostTemp < uiCostBi )
        {
          bChanged = true;
          ::memcpy( cMvBi[iRefList], cMvTemp[iRefList][iRefIdxTemp], sizeof(Mv)*3 );
          iRefIdxBi[iRefList] = iRefIdxTemp;

          uiCostBi            = uiCostTemp;
          uiMotBits[iRefList] = uiBitsTemp - uiMbBits[2] - uiMotBits[1-iRefList];
          uiMotBits[iRefList] -= ((pu.cu->slice->getSPS()->getUseBcw() == true) ? bcwIdxBits : 0);
          uiBits[2]           = uiBitsTemp;

          if ( iNumIter != 1 ) // MC for next iter
          {
            //  Set motion
            PU::setAllAffineMv( pu, cMvBi[iRefList][0], cMvBi[iRefList][1], cMvBi[iRefList][2], eRefPicList);
            pu.refIdx[eRefPicList] = iRefIdxBi[eRefPicList];
            PelUnitBuf predBufTmp = m_tmpPredStorage[iRefList].getBuf( UnitAreaRelative(*pu.cu, pu) );
            motionCompensation( pu, predBufTmp, eRefPicList );
          }
        }
      } // for loop-iRefIdxTemp

      if ( !bChanged )
      {
        if ((uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1]) || enforceBcwPred)
        {
          xCopyAffineAMVPInfo( aacAffineAMVPInfo[0][iRefIdxBi[0]], affiAMVPInfoTemp[REF_PIC_LIST_0] );
          xCheckBestAffineMVP( pu, affiAMVPInfoTemp[REF_PIC_LIST_0], REF_PIC_LIST_0, cMvBi[0], cMvPredBi[0][iRefIdxBi[0]], aaiMvpIdxBi[0][iRefIdxBi[0]], uiBits[2], uiCostBi );

          if ( !slice.getPicHeader()->getMvdL1ZeroFlag() )
          {
            xCopyAffineAMVPInfo( aacAffineAMVPInfo[1][iRefIdxBi[1]], affiAMVPInfoTemp[REF_PIC_LIST_1] );
            xCheckBestAffineMVP( pu, affiAMVPInfoTemp[REF_PIC_LIST_1], REF_PIC_LIST_1, cMvBi[1], cMvPredBi[1][iRefIdxBi[1]], aaiMvpIdxBi[1][iRefIdxBi[1]], uiBits[2], uiCostBi );
          }
        }
        break;
      }
    } // for loop-iter
    }
    m_isBi = false;
  } // if (B_SLICE)

  pu.mv    [REF_PIC_LIST_0] = Mv();
  pu.mv    [REF_PIC_LIST_1] = Mv();
  pu.mvd   [REF_PIC_LIST_0] = cMvZero;
  pu.mvd   [REF_PIC_LIST_1] = cMvZero;
  pu.refIdx[REF_PIC_LIST_0] = NOT_VALID;
  pu.refIdx[REF_PIC_LIST_1] = NOT_VALID;
  pu.mvpIdx[REF_PIC_LIST_0] = NOT_VALID;
  pu.mvpIdx[REF_PIC_LIST_1] = NOT_VALID;
  pu.mvpNum[REF_PIC_LIST_0] = NOT_VALID;
  pu.mvpNum[REF_PIC_LIST_1] = NOT_VALID;

  for ( int verIdx = 0; verIdx < 3; verIdx++ )
  {
    pu.mvdAffi[REF_PIC_LIST_0][verIdx] = cMvZero;
    pu.mvdAffi[REF_PIC_LIST_1][verIdx] = cMvZero;
  }

  // Set Motion Field
  memcpy( aacMv[1], mvValidList1, sizeof(Mv)*3 );
  iRefIdx[1] = refIdxValidList1;
  uiBits[1]  = bitsValidList1;
  uiCost[1]  = costValidList1;
  if (pu.cs->pps->getWPBiPred() == true && tryBipred && (bcwIdx != BCW_DEFAULT))
  {
    CHECK(iRefIdxBi[0]<0, "Invalid picture reference index");
    CHECK(iRefIdxBi[1]<0, "Invalid picture reference index");
    wp0 = pu.cs->slice->getWpScaling(REF_PIC_LIST_0, iRefIdxBi[0]);
    wp1 = pu.cs->slice->getWpScaling(REF_PIC_LIST_1, iRefIdxBi[1]);

    if (WPScalingParam::isWeighted(wp0) || WPScalingParam::isWeighted(wp1))
    {
      uiCostBi = MAX_UINT;
      enforceBcwPred = false;
    }
  }
  if( enforceBcwPred )
  {
    uiCost[0] = uiCost[1] = MAX_UINT;
  }

  // Affine ME result set
  if ( uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1] ) // Bi
  {
    lastMode = 2;
    affineCost = uiCostBi;
    pu.interDir = 3;

    pu.refIdx[REF_PIC_LIST_0] = iRefIdxBi[0];
    pu.refIdx[REF_PIC_LIST_1] = iRefIdxBi[1];

    for ( int verIdx = 0; verIdx < mvNum; verIdx++ )
    {
      pu.mvAffi[REF_PIC_LIST_0][verIdx] = cMvBi[0][verIdx];
      pu.mvAffi[REF_PIC_LIST_1][verIdx] = cMvBi[1][verIdx];
      pu.mvdAffi[REF_PIC_LIST_0][verIdx] = cMvBi[0][verIdx] - cMvPredBi[0][iRefIdxBi[0]][verIdx];
      pu.mvdAffi[REF_PIC_LIST_1][verIdx] = cMvBi[1][verIdx] - cMvPredBi[1][iRefIdxBi[1]][verIdx];

      if ( verIdx != 0 )
      {
        pu.mvdAffi[0][verIdx] = pu.mvdAffi[0][verIdx] - pu.mvdAffi[0][0];
        pu.mvdAffi[1][verIdx] = pu.mvdAffi[1][verIdx] - pu.mvdAffi[1][0];
      }
    }


    pu.mvpIdx[REF_PIC_LIST_0] = aaiMvpIdxBi[0][iRefIdxBi[0]];
    pu.mvpNum[REF_PIC_LIST_0] = aaiMvpNum[0][iRefIdxBi[0]];
    pu.mvpIdx[REF_PIC_LIST_1] = aaiMvpIdxBi[1][iRefIdxBi[1]];
    pu.mvpNum[REF_PIC_LIST_1] = aaiMvpNum[1][iRefIdxBi[1]];
  }
  else if ( uiCost[0] <= uiCost[1] ) // List 0
  {
    lastMode = 0;
    affineCost = uiCost[0];
    pu.interDir = 1;
    pu.mv[1].setZero();
    pu.refIdx[REF_PIC_LIST_0] = iRefIdx[0];

    for ( int verIdx = 0; verIdx < mvNum; verIdx++ )
    {
      pu.mvAffi[REF_PIC_LIST_0][verIdx] = aacMv[0][verIdx];
      pu.mvdAffi[REF_PIC_LIST_0][verIdx] = aacMv[0][verIdx] - cMvPred[0][iRefIdx[0]][verIdx];
      if ( verIdx != 0 )
      {
        pu.mvdAffi[0][verIdx] = pu.mvdAffi[0][verIdx] - pu.mvdAffi[0][0];
      }
    }

    pu.mvpIdx[REF_PIC_LIST_0] = aaiMvpIdx[0][iRefIdx[0]];
    pu.mvpNum[REF_PIC_LIST_0] = aaiMvpNum[0][iRefIdx[0]];
  }
  else
  {
    lastMode = 1;
    affineCost = uiCost[1];
    pu.interDir = 2;
    pu.mv[0].setZero();
    pu.refIdx[REF_PIC_LIST_1] = iRefIdx[1];

    for ( int verIdx = 0; verIdx < mvNum; verIdx++ )
    {
      pu.mvAffi[REF_PIC_LIST_1][verIdx] = aacMv[1][verIdx];
      pu.mvdAffi[REF_PIC_LIST_1][verIdx] = aacMv[1][verIdx] - cMvPred[1][iRefIdx[1]][verIdx];
      if ( verIdx != 0 )
      {
        pu.mvdAffi[1][verIdx] = pu.mvdAffi[1][verIdx] - pu.mvdAffi[1][0];
      }
    }

    pu.mvpIdx[REF_PIC_LIST_1] = aaiMvpIdx[1][iRefIdx[1]];
    pu.mvpNum[REF_PIC_LIST_1] = aaiMvpNum[1][iRefIdx[1]];
  }
  if( bcwIdx != BCW_DEFAULT )
  {
    pu.cu->BcwIdx = BCW_DEFAULT;
  }
}

// Ax = b, m = {A, b}
#if AFFINE_ENC_OPT
void solveGaussElimination( double( *m )[7], double *x, int num )
{
#define NEARZERO(x) x == 0.

  const int numM1 = num - 1;

  for( int i = 0; i < numM1; i++ )
  {
    // find non-zero diag
    int tempIdx = i;
    if( NEARZERO( m[i][i] ) )
    {
      for( int j = i + 1; j < num; j++ )
      {
        if( !( NEARZERO( m[j][i] ) ) )
        {
          tempIdx = j;
          break;
        }
      }
    }

    // swap line
    if( tempIdx != i )
    {
      swap( m[i], m[tempIdx] );
    }

    double* currRow = m[i];
    const double diagCoeff = currRow[i];

    if( NEARZERO( diagCoeff ) )
    {
      std::memset( x, 0, sizeof( *x ) * num );
      return;
    }

    // eliminate column
    for( int j = i + 1; j < num; j++ )
    {
      double* rowCoeff = m[j];
      const double coeffRatio = rowCoeff[i] / diagCoeff;

      for( int k = i + 1; k <= num; k++ )
      {
        rowCoeff[k] -= currRow[k] * coeffRatio;
      }
    }
  }

  if( NEARZERO( m[numM1][numM1] ) )
  {
    std::memset( x, 0, sizeof( *x ) * num );
    return;
  }

  double* currRow = m[numM1];
  x[numM1] = currRow[num] / currRow[numM1];

  for( int i = num - 2; i >= 0; i-- )
  {
    currRow = m[i];
    const double diagCoeff = currRow[i];

    if( NEARZERO( diagCoeff ) )
    {
      std::memset( x, 0, sizeof( *x ) * num );
      return;
    }

    double temp = 0;
    for( int j = i + 1; j < num; j++ )
    {
      temp += currRow[j] * x[j];
    }
    x[i] = ( currRow[num] - temp ) / diagCoeff;
  }
#undef NEARZERO
}
#else
void solveEqual( double dEqualCoeff[7][7], int iOrder, double *dAffinePara )
{
  for( int k = 0; k < iOrder; k++ )
  {
    dAffinePara[k] = 0.;
  }

  // row echelon
  for( int i = 1; i < iOrder; i++ )
  {
    // find column max
    double temp = fabs( dEqualCoeff[i][i - 1] );
    int tempIdx = i;
    for( int j = i + 1; j < iOrder + 1; j++ )
    {
      if( fabs( dEqualCoeff[j][i - 1] ) > temp )
      {
        temp = fabs( dEqualCoeff[j][i - 1] );
        tempIdx = j;
      }
    }

    // swap line
    if( tempIdx != i )
    {
      for( int j = 0; j < iOrder + 1; j++ )
      {
        dEqualCoeff[0][j] = dEqualCoeff[i][j];
        dEqualCoeff[i][j] = dEqualCoeff[tempIdx][j];
        dEqualCoeff[tempIdx][j] = dEqualCoeff[0][j];
      }
    }

    // elimination first column
    if( dEqualCoeff[i][i - 1] == 0. )
    {
      return;
    }
    for( int j = i + 1; j < iOrder + 1; j++ )
    {
      for( int k = i; k < iOrder + 1; k++ )
      {
        dEqualCoeff[j][k] = dEqualCoeff[j][k] - dEqualCoeff[i][k] * dEqualCoeff[j][i - 1] / dEqualCoeff[i][i - 1];
      }
    }
  }

  if( dEqualCoeff[iOrder][iOrder - 1] == 0. )
  {
    return;
  }
  dAffinePara[iOrder - 1] = dEqualCoeff[iOrder][iOrder] / dEqualCoeff[iOrder][iOrder - 1];
  for( int i = iOrder - 2; i >= 0; i-- )
  {
    if( dEqualCoeff[i + 1][i] == 0. )
    {
      for( int k = 0; k < iOrder; k++ )
      {
        dAffinePara[k] = 0.;
      }
      return;
    }
    double temp = 0;
    for( int j = i + 1; j < iOrder; j++ )
    {
      temp += dEqualCoeff[i + 1][j] * dAffinePara[j];
    }
    dAffinePara[i] = ( dEqualCoeff[i + 1][iOrder] - temp ) / dEqualCoeff[i + 1][i];
  }
}
#endif

void InterSearch::xCheckBestAffineMVP( PredictionUnit &pu, AffineAMVPInfo &affineAMVPInfo, RefPicList eRefPicList, Mv acMv[3], Mv acMvPred[3], int& riMVPIdx, uint32_t& ruiBits, Distortion& ruiCost )
{
  if ( affineAMVPInfo.numCand < 2 )
  {
    return;
  }

  int mvNum = pu.cu->affineType ? 3 : 2;

  m_pcRdCost->selectMotionLambda( );
  m_pcRdCost->setCostScale ( 0 );

  int iBestMVPIdx = riMVPIdx;

  // Get origin MV bits
  Mv tmpPredMv[3];
  int iOrgMvBits = xCalcAffineMVBits( pu, acMv, acMvPred );
  iOrgMvBits += m_auiMVPIdxCost[riMVPIdx][AMVP_MAX_NUM_CANDS];

  int iBestMvBits = iOrgMvBits;
  for (int iMVPIdx = 0; iMVPIdx < affineAMVPInfo.numCand; iMVPIdx++)
  {
    if (iMVPIdx == riMVPIdx)
    {
      continue;
    }
    tmpPredMv[0] = affineAMVPInfo.mvCandLT[iMVPIdx];
    tmpPredMv[1] = affineAMVPInfo.mvCandRT[iMVPIdx];
    if ( mvNum == 3 )
    {
      tmpPredMv[2] = affineAMVPInfo.mvCandLB[iMVPIdx];
    }
    int iMvBits = xCalcAffineMVBits( pu, acMv, tmpPredMv );
    iMvBits += m_auiMVPIdxCost[iMVPIdx][AMVP_MAX_NUM_CANDS];

    if (iMvBits < iBestMvBits)
    {
      iBestMvBits = iMvBits;
      iBestMVPIdx = iMVPIdx;
    }
  }

  if (iBestMVPIdx != riMVPIdx)  // if changed
  {
    acMvPred[0] = affineAMVPInfo.mvCandLT[iBestMVPIdx];
    acMvPred[1] = affineAMVPInfo.mvCandRT[iBestMVPIdx];
    acMvPred[2] = affineAMVPInfo.mvCandLB[iBestMVPIdx];
    riMVPIdx = iBestMVPIdx;
    uint32_t uiOrgBits = ruiBits;
    ruiBits = uiOrgBits - iOrgMvBits + iBestMvBits;
    ruiCost = (ruiCost - m_pcRdCost->getCost( uiOrgBits )) + m_pcRdCost->getCost( ruiBits );
  }
}

void InterSearch::xAffineMotionEstimation( PredictionUnit& pu,
                                           PelUnitBuf&     origBuf,
                                           RefPicList      eRefPicList,
                                           Mv              acMvPred[3],
                                           int             iRefIdxPred,
                                           Mv              acMv[3],
                                           uint32_t&           ruiBits,
                                           Distortion&     ruiCost,
                                           int&            mvpIdx,
                                           const AffineAMVPInfo& aamvpi,
                                           bool            bBi)
{
  if( pu.cu->cs->sps->getUseBcw() && pu.cu->BcwIdx != BCW_DEFAULT && !bBi && xReadBufferedAffineUniMv(pu, eRefPicList, iRefIdxPred, acMvPred, acMv, ruiBits, ruiCost
      , mvpIdx, aamvpi
  ) )
  {
    return;
  }

  uint32_t dirBits = ruiBits - m_auiMVPIdxCost[mvpIdx][aamvpi.numCand];
  int bestMvpIdx   = mvpIdx;
  const int width  = pu.Y().width;
  const int height = pu.Y().height;

  const Picture* refPic = pu.cu->slice->getRefPic(eRefPicList, iRefIdxPred);

  // Set Origin YUV: pcYuv
  PelUnitBuf*   pBuf = &origBuf;
  double        fWeight       = 1.0;

  PelUnitBuf  origBufTmp = m_tmpStorageLCU.getBuf( UnitAreaRelative( *pu.cu, pu ) );
  enum DFunc distFunc = (pu.cs->slice->getDisableSATDForRD()) ? DF_SAD : DF_HAD;
  m_iRefListIdx = eRefPicList;

  // if Bi, set to ( 2 * Org - ListX )
  if ( bBi )
  {
    // NOTE: Other buf contains predicted signal from another direction
    PelUnitBuf otherBuf = m_tmpPredStorage[1 - (int)eRefPicList].getBuf( UnitAreaRelative( *pu.cu, pu ) );
    origBufTmp.copyFrom(origBuf);
    origBufTmp.removeHighFreq(otherBuf, m_pcEncCfg->getClipForBiPredMeEnabled(), pu.cu->slice->clpRngs()
                             ,getBcwWeight(pu.cu->BcwIdx, eRefPicList)
                             );
    pBuf = &origBufTmp;

    fWeight = xGetMEDistortionWeight( pu.cu->BcwIdx, eRefPicList );
  }

  // pred YUV
  PelUnitBuf  predBuf = m_tmpAffiStorage.getBuf( UnitAreaRelative(*pu.cu, pu) );

  // Set start Mv position, use input mv as started search mv
  Mv acMvTemp[3];
  ::memcpy( acMvTemp, acMv, sizeof(Mv)*3 );
  // Set delta mv
  // malloc buffer
  int iParaNum = pu.cu->affineType ? 7 : 5;
  int affineParaNum = iParaNum - 1;
  int mvNum = pu.cu->affineType ? 3 : 2;

  int64_t  i64EqualCoeff[7][7];
  Pel    *piError = m_tmpAffiError;
#if AFFINE_ENC_OPT // using Pel instead of int for better SIMD
  Pel    *pdDerivate[2];
#else
  int    *pdDerivate[2];
#endif
  pdDerivate[0] = m_tmpAffiDeri[0];
  pdDerivate[1] = m_tmpAffiDeri[1];

  Distortion uiCostBest = std::numeric_limits<Distortion>::max();
  uint32_t uiBitsBest = 0;

  // do motion compensation with origin mv
  if( m_pcEncCfg->getMCTSEncConstraint() )
  {
    Area curTileAreaRestricted = pu.cs->picture->mctsInfo.getTileAreaSubPelRestricted( pu );
    MCTSHelper::clipMvToArea( acMvTemp[0], pu.cu->Y(), curTileAreaRestricted, *pu.cs->sps );
    MCTSHelper::clipMvToArea( acMvTemp[1], pu.cu->Y(), curTileAreaRestricted, *pu.cs->sps );
    if( pu.cu->affineType == AFFINEMODEL_6PARAM )
    {
      MCTSHelper::clipMvToArea( acMvTemp[2], pu.cu->Y(), curTileAreaRestricted, *pu.cs->sps );
    }
  }
  else
  {
    clipMv( acMvTemp[0], pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
    clipMv( acMvTemp[1], pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
    if( pu.cu->affineType == AFFINEMODEL_6PARAM )
    {
      clipMv( acMvTemp[2], pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
    }
  }
  acMvTemp[0].roundAffinePrecInternal2Amvr(pu.cu->imv);
  acMvTemp[1].roundAffinePrecInternal2Amvr(pu.cu->imv);
  if (pu.cu->affineType == AFFINEMODEL_6PARAM)
  {
    acMvTemp[2].roundAffinePrecInternal2Amvr(pu.cu->imv);
  }
#if AFFINE_ENC_OPT
  int gStride = width;
#if JVET_Z0136_OOB
  xPredAffineBlk(COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cs->slice->clpRng(COMPONENT_Y), eRefPicList, false, SCALE_1X, true);
#else
  xPredAffineBlk(COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cs->slice->clpRng(COMPONENT_Y), false, SCALE_1X, true);
#endif
#else
#if JVET_Z0136_OOB
  xPredAffineBlk( COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cs->slice->clpRng(COMPONENT_Y), eRefPicList );
#else
  xPredAffineBlk( COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cs->slice->clpRng( COMPONENT_Y ) );
#endif
#endif

  // get error
  uiCostBest = m_pcRdCost->getDistPart(predBuf.Y(), pBuf->Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, distFunc);

  // get cost with mv
  m_pcRdCost->setCostScale(0);
  uiBitsBest = ruiBits;
  if ( pu.cu->imv == 2 && m_pcEncCfg->getUseAffineAmvrEncOpt() )
  {
    uiBitsBest  = dirBits + xDetermineBestMvp( pu, acMvTemp, mvpIdx, aamvpi );
    acMvPred[0] = aamvpi.mvCandLT[mvpIdx];
    acMvPred[1] = aamvpi.mvCandRT[mvpIdx];
    acMvPred[2] = aamvpi.mvCandLB[mvpIdx];
  }
  else
  {
    DTRACE( g_trace_ctx, D_COMMON, " (%d) xx uiBitsBest=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiBitsBest );
    uiBitsBest += xCalcAffineMVBits( pu, acMvTemp, acMvPred );
    DTRACE( g_trace_ctx, D_COMMON, " (%d) yy uiBitsBest=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiBitsBest );
  }
  uiCostBest = (Distortion)( floor( fWeight * (double)uiCostBest ) + (double)m_pcRdCost->getCost( uiBitsBest ) );

  DTRACE( g_trace_ctx, D_COMMON, " (%d) uiBitsBest=%d, uiCostBest=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiBitsBest, uiCostBest );

  ::memcpy( acMv, acMvTemp, sizeof(Mv) * 3 );

  const int bufStride = pBuf->Y().stride;
  const int predBufStride = predBuf.Y().stride;
  Mv prevIterMv[7][3];
  int iIterTime;
  if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
  {
    iIterTime = bBi ? 3 : 4;
  }
  else
  {
    iIterTime = bBi ? 3 : 5;
  }

  if ( !pu.cu->cs->sps->getUseAffineType() )
  {
    iIterTime = bBi ? 5 : 7;
  }
  for ( int iter=0; iter<iIterTime; iter++ )    // iterate loop
  {
    memcpy( prevIterMv[iter], acMvTemp, sizeof( Mv ) * 3 );
    /*********************************************************************************
     *                         use gradient to update mv
     *********************************************************************************/
    // get Error Matrix
    Pel* pOrg  = pBuf->Y().buf;
    Pel* pPred = predBuf.Y().buf;
    Pel* error = piError;

    for ( int j=0; j< height; j++ )
    {
      for ( int i=0; i< width; i++ )
      {
        error[i] = pOrg[i] - pPred[i];
      }
      pOrg  += bufStride;
      pPred += predBufStride;
      error += width;
    }

#if AFFINE_ENC_OPT
    pdDerivate[0] = m_gradX0 + gStride + 1;
    pdDerivate[1] = m_gradY0 + gStride + 1;
#else
    // sobel x direction
    // -1 0 1
    // -2 0 2
    // -1 0 1
    pPred = predBuf.Y().buf;
    m_HorizontalSobelFilter( pPred, predBufStride, pdDerivate[0], width, width, height );

    // sobel y direction
    // -1 -2 -1
    //  0  0  0
    //  1  2  1
    m_VerticalSobelFilter( pPred, predBufStride, pdDerivate[1], width, width, height );
#endif

    // solve delta x and y
    for ( int row = 0; row < iParaNum; row++ )
    {
      memset( &i64EqualCoeff[row][0], 0, iParaNum * sizeof( int64_t ) );
    }

#if AFFINE_ENC_OPT
    // the "6" is the shift number in gradient (canculated in IF_INTERNAL_PREC precision), "-1" is for gradient normalization
    // the input parameter "shift" in is to compensate dI with regard to the gradient
    m_EqualCoeffComputer( piError, width, pdDerivate, gStride, i64EqualCoeff, width, height
      , (pu.cu->affineType == AFFINEMODEL_6PARAM)
      , 6 - 1 - std::max<int>(2, (IF_INTERNAL_PREC - pu.cs->slice->clpRng(COMPONENT_Y).bd))
    );
#else
    m_EqualCoeffComputer( piError, width, pdDerivate, width, i64EqualCoeff, width, height
      , (pu.cu->affineType == AFFINEMODEL_6PARAM)
    );
#endif

    double dAffinePara[6];
    double dDeltaMv[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0,};
    Mv acDeltaMv[3];

#if AFFINE_ENC_OPT
    double pdEqualCoeff[6][7];

    for( int row = 0; row < affineParaNum; row++ )
    {
      double* dCoeff = pdEqualCoeff[row];
      int64_t* iCoeff = i64EqualCoeff[row + 1];

      for( int i = 0; i < iParaNum; i++ )
      {
        dCoeff[i] = ( double ) iCoeff[i];
      }
    }

    solveGaussElimination( pdEqualCoeff, dAffinePara, affineParaNum );
#else
    double pdEqualCoeff[7][7];
    for( int row = 0; row < iParaNum; row++ )
    {
      for( int i = 0; i < iParaNum; i++ )
      {
        pdEqualCoeff[row][i] = ( double ) i64EqualCoeff[row][i];
      }
    }

    solveEqual( pdEqualCoeff, affineParaNum, dAffinePara );
#endif

    // convert to delta mv
    dDeltaMv[0] = dAffinePara[0];
    dDeltaMv[2] = dAffinePara[2];
    if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
    {
      dDeltaMv[1] = dAffinePara[1] * width + dAffinePara[0];
      dDeltaMv[3] = dAffinePara[3] * width + dAffinePara[2];
      dDeltaMv[4] = dAffinePara[4] * height + dAffinePara[0];
      dDeltaMv[5] = dAffinePara[5] * height + dAffinePara[2];
    }
    else
    {
      dDeltaMv[1] = dAffinePara[1] * width + dAffinePara[0];
      dDeltaMv[3] = -dAffinePara[3] * width + dAffinePara[2];
    }

    const int normShiftTab[3] = { MV_PRECISION_QUARTER - MV_PRECISION_INT, MV_PRECISION_SIXTEENTH - MV_PRECISION_INT, MV_PRECISION_QUARTER - MV_PRECISION_INT };
    const int stepShiftTab[3] = { MV_PRECISION_INTERNAL - MV_PRECISION_QUARTER, MV_PRECISION_INTERNAL - MV_PRECISION_SIXTEENTH, MV_PRECISION_INTERNAL - MV_PRECISION_QUARTER };
    const int multiShift = 1 << normShiftTab[pu.cu->imv];
    const int mvShift = stepShiftTab[pu.cu->imv];
    acDeltaMv[0] = Mv( ( int ) ( dDeltaMv[0] * multiShift + SIGN( dDeltaMv[0] ) * 0.5 ) << mvShift, ( int ) ( dDeltaMv[2] * multiShift + SIGN( dDeltaMv[2] ) * 0.5 ) << mvShift );
    acDeltaMv[1] = Mv( ( int ) ( dDeltaMv[1] * multiShift + SIGN( dDeltaMv[1] ) * 0.5 ) << mvShift, ( int ) ( dDeltaMv[3] * multiShift + SIGN( dDeltaMv[3] ) * 0.5 ) << mvShift );
    if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
    {
      acDeltaMv[2] = Mv( ( int ) ( dDeltaMv[4] * multiShift + SIGN( dDeltaMv[4] ) * 0.5 ) << mvShift, ( int ) ( dDeltaMv[5] * multiShift + SIGN( dDeltaMv[5] ) * 0.5 ) << mvShift );
    }
    if ( !m_pcEncCfg->getUseAffineAmvrEncOpt() )
    {
      bool bAllZero = false;
      for ( int i = 0; i < mvNum; i++ )
      {
        Mv deltaMv = acDeltaMv[i];
        if ( pu.cu->imv == 2 )
        {
          deltaMv.roundToPrecision( MV_PRECISION_INTERNAL, MV_PRECISION_HALF );
        }
        if ( deltaMv.getHor() != 0 || deltaMv.getVer() != 0 )
        {
          bAllZero = false;
          break;
        }
        bAllZero = true;
      }

      if ( bAllZero )
        break;
    }
    // do motion compensation with updated mv
    for ( int i = 0; i < mvNum; i++ )
    {
      acMvTemp[i] += acDeltaMv[i];
      acMvTemp[i].hor = Clip3(MV_MIN, MV_MAX, acMvTemp[i].hor );
      acMvTemp[i].ver = Clip3(MV_MIN, MV_MAX, acMvTemp[i].ver );
      acMvTemp[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
      if( m_pcEncCfg->getMCTSEncConstraint() )
      {
        MCTSHelper::clipMvToArea( acMvTemp[i], pu.cu->Y(), pu.cs->picture->mctsInfo.getTileAreaSubPelRestricted( pu ), *pu.cs->sps );
      }
      else
      {
        clipMv( acMvTemp[i], pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
      }
    }

    if ( m_pcEncCfg->getUseAffineAmvrEncOpt() )
    {
      bool identical = false;
      for ( int k = iter; k >= 0; k-- )
      {
        if ( acMvTemp[0] == prevIterMv[k][0] && acMvTemp[1] == prevIterMv[k][1] )
        {
          identical = pu.cu->affineType ? acMvTemp[2] == prevIterMv[k][2] : true;
          if ( identical )
          {
            break;
          }
        }
      }
      if ( identical )
      {
        break;
      }
    }

#if AFFINE_ENC_OPT
#if JVET_Z0136_OOB
    xPredAffineBlk(COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cu->slice->clpRng(COMPONENT_Y), eRefPicList, false, SCALE_1X, true);
#else
    xPredAffineBlk( COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cu->slice->clpRng( COMPONENT_Y ), false, SCALE_1X, true );
#endif
#else
#if JVET_Z0136_OOB
    xPredAffineBlk( COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cu->slice->clpRng(COMPONENT_Y), eRefPicList  );
#else
    xPredAffineBlk( COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cu->slice->clpRng( COMPONENT_Y ) );
#endif
#endif

    // get error
    Distortion uiCostTemp = m_pcRdCost->getDistPart(predBuf.Y(), pBuf->Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, distFunc);
    DTRACE( g_trace_ctx, D_COMMON, " (%d) uiCostTemp=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiCostTemp );

    // get cost with mv
    m_pcRdCost->setCostScale(0);
    uint32_t uiBitsTemp = ruiBits;
    if ( pu.cu->imv == 2 && m_pcEncCfg->getUseAffineAmvrEncOpt() )
    {
      uiBitsTemp  = dirBits + xDetermineBestMvp( pu, acMvTemp, bestMvpIdx, aamvpi );
      acMvPred[0] = aamvpi.mvCandLT[bestMvpIdx];
      acMvPred[1] = aamvpi.mvCandRT[bestMvpIdx];
      acMvPred[2] = aamvpi.mvCandLB[bestMvpIdx];
    }
    else
    {
      uiBitsTemp += xCalcAffineMVBits( pu, acMvTemp, acMvPred );
    }
    uiCostTemp = (Distortion)( floor( fWeight * (double)uiCostTemp ) + (double)m_pcRdCost->getCost( uiBitsTemp ) );

    // store best cost and mv
    if ( uiCostTemp < uiCostBest )
    {
      uiCostBest = uiCostTemp;
      uiBitsBest = uiBitsTemp;
      memcpy( acMv, acMvTemp, sizeof(Mv) * 3 );
      mvpIdx = bestMvpIdx;
    }
  }

  auto checkCPMVRdCost = [&](Mv ctrlPtMv[3])
  {
#if JVET_Z0136_OOB
    xPredAffineBlk(COMPONENT_Y, pu, refPic, ctrlPtMv, predBuf, false, pu.cu->slice->clpRng(COMPONENT_Y), eRefPicList);
#else
    xPredAffineBlk(COMPONENT_Y, pu, refPic, ctrlPtMv, predBuf, false, pu.cu->slice->clpRng(COMPONENT_Y));
#endif
    // get error
    Distortion costTemp = m_pcRdCost->getDistPart(predBuf.Y(), pBuf->Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, distFunc);
    // get cost with mv
    m_pcRdCost->setCostScale(0);
    uint32_t bitsTemp = ruiBits;
    bitsTemp += xCalcAffineMVBits( pu, ctrlPtMv, acMvPred );
    costTemp = (Distortion)(floor(fWeight * (double)costTemp) + (double)m_pcRdCost->getCost(bitsTemp));
    // store best cost and mv
    if (costTemp < uiCostBest)
    {
      uiCostBest = costTemp;
      uiBitsBest = bitsTemp;
      ::memcpy(acMv, ctrlPtMv, sizeof(Mv) * 3);
    }
  };

  const uint32_t mvShiftTable[3] = {MV_PRECISION_INTERNAL - MV_PRECISION_QUARTER, MV_PRECISION_INTERNAL - MV_PRECISION_INTERNAL, MV_PRECISION_INTERNAL - MV_PRECISION_INT};
  const uint32_t mvShift = mvShiftTable[pu.cu->imv];
  if (uiCostBest <= AFFINE_ME_LIST_MVP_TH*m_hevcCost)
  {

    Mv mvPredTmp[3] = { acMvPred[0], acMvPred[1], acMvPred[2] };
    Mv mvME[3];
    ::memcpy(mvME, acMv, sizeof(Mv) * 3);
    Mv dMv = mvME[0] - mvPredTmp[0];

    for (int j = 0; j < mvNum; j++)
    {
      if ((!j && mvME[j] != mvPredTmp[j]) || (j && mvME[j] != (mvPredTmp[j] + dMv)))
      {
        ::memcpy(acMvTemp, mvME, sizeof(Mv) * 3);
        acMvTemp[j] = mvPredTmp[j];

        if (j)
          acMvTemp[j] += dMv;

        checkCPMVRdCost(acMvTemp);
      }
    }

    //keep the rotation/zoom;
    if (mvME[0] != mvPredTmp[0])
    {
      ::memcpy(acMvTemp, mvME, sizeof(Mv) * 3);
      for (int i = 1; i < mvNum; i++)
      {
        acMvTemp[i] -= dMv;
      }
      acMvTemp[0] = mvPredTmp[0];

      checkCPMVRdCost(acMvTemp);
    }

    //keep the translation;
    if (pu.cu->affineType == AFFINEMODEL_6PARAM && mvME[1] != (mvPredTmp[1] + dMv) && mvME[2] != (mvPredTmp[2] + dMv))
    {
      ::memcpy(acMvTemp, mvME, sizeof(Mv) * 3);

      acMvTemp[1] = mvPredTmp[1] + dMv;
      acMvTemp[2] = mvPredTmp[2] + dMv;

      checkCPMVRdCost(acMvTemp);
    }

    // 8 nearest neighbor search
    int testPos[8][2] = { { -1, 0 },{ 0, -1 },{ 0, 1 },{ 1, 0 },{ -1, -1 },{ -1, 1 },{ 1, 1 },{ 1, -1 } };
    const int maxSearchRound = (pu.cu->imv) ? 3 : ((m_pcEncCfg->getUseAffineAmvrEncOpt() && m_pcEncCfg->getIntraPeriod() == (uint32_t)-1) ? 2 : 3);

    for (int rnd = 0; rnd < maxSearchRound; rnd++)
    {
      bool modelChange = false;
      //search the model parameters with finear granularity;
      for (int j = 0; j < mvNum; j++)
      {
        bool loopChange = false;
        for (int iter = 0; iter < 2; iter++)
        {
          if (iter == 1 && !loopChange)
          {
            break;
          }
          Mv centerMv[3];
          memcpy(centerMv, acMv, sizeof(Mv) * 3);
          memcpy(acMvTemp, acMv, sizeof(Mv) * 3);

          for (int i = ((iter == 0) ? 0 : 4); i < ((iter == 0) ? 4 : 8); i++)
          {
            acMvTemp[j].set(centerMv[j].getHor() + (testPos[i][0] << mvShift), centerMv[j].getVer() + (testPos[i][1] << mvShift));
            clipMv( acMvTemp[j], pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
#if JVET_Z0136_OOB
            xPredAffineBlk(COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cu->slice->clpRng(COMPONENT_Y), eRefPicList);
#else
            xPredAffineBlk(COMPONENT_Y, pu, refPic, acMvTemp, predBuf, false, pu.cu->slice->clpRng(COMPONENT_Y));
#endif
            Distortion costTemp = m_pcRdCost->getDistPart(predBuf.Y(), pBuf->Y(), pu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, distFunc);
            uint32_t bitsTemp = ruiBits;
            bitsTemp += xCalcAffineMVBits(pu, acMvTemp, acMvPred);
            costTemp = (Distortion)(floor(fWeight * (double)costTemp) + (double)m_pcRdCost->getCost(bitsTemp));

            if (costTemp < uiCostBest)
            {
              uiCostBest = costTemp;
              uiBitsBest = bitsTemp;
              ::memcpy(acMv, acMvTemp, sizeof(Mv) * 3);
              modelChange = true;
              loopChange = true;
            }
          }
        }
      }

      if (!modelChange)
      {
        break;
      }
    }
  }
  acMvPred[0] = aamvpi.mvCandLT[mvpIdx];
  acMvPred[1] = aamvpi.mvCandRT[mvpIdx];
  acMvPred[2] = aamvpi.mvCandLB[mvpIdx];

  ruiBits = uiBitsBest;
  ruiCost = uiCostBest;
  DTRACE( g_trace_ctx, D_COMMON, " (%d) uiBitsBest=%d, uiCostBest=%d\n", DTRACE_GET_COUNTER(g_trace_ctx,D_COMMON), uiBitsBest, uiCostBest );
}

void InterSearch::xEstimateAffineAMVP( PredictionUnit&  pu,
                                       AffineAMVPInfo&  affineAMVPInfo,
                                       PelUnitBuf&      origBuf,
                                       RefPicList       eRefPicList,
                                       int              iRefIdx,
                                       Mv               acMvPred[3],
                                       Distortion*      puiDistBiP )
{
  Mv         bestMvLT, bestMvRT, bestMvLB;
  int        iBestIdx = 0;
  Distortion uiBestCost = std::numeric_limits<Distortion>::max();

  // Fill the MV Candidates
  PU::fillAffineMvpCand( pu, eRefPicList, iRefIdx, affineAMVPInfo );
  CHECK( affineAMVPInfo.numCand == 0, "Assertion failed." );

  PelUnitBuf predBuf = m_tmpStorageLCU.getBuf( UnitAreaRelative(*pu.cu, pu) );

  // initialize Mvp index & Mvp
  iBestIdx = 0;
  for( int i = 0 ; i < affineAMVPInfo.numCand; i++ )
  {
    Mv mv[3] = { affineAMVPInfo.mvCandLT[i], affineAMVPInfo.mvCandRT[i], affineAMVPInfo.mvCandLB[i] };

    Distortion uiTmpCost = xGetAffineTemplateCost( pu, origBuf, predBuf, mv, i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx );

    if ( uiBestCost > uiTmpCost )
    {
      uiBestCost = uiTmpCost;
      bestMvLT = affineAMVPInfo.mvCandLT[i];
      bestMvRT = affineAMVPInfo.mvCandRT[i];
      bestMvLB = affineAMVPInfo.mvCandLB[i];
      iBestIdx  = i;
      *puiDistBiP = uiTmpCost;
    }
  }

  // Setting Best MVP
  acMvPred[0] = bestMvLT;
  acMvPred[1] = bestMvRT;
  acMvPred[2] = bestMvLB;

  pu.mvpIdx[eRefPicList] = iBestIdx;
  pu.mvpNum[eRefPicList] = affineAMVPInfo.numCand;
  DTRACE( g_trace_ctx, D_COMMON, "#estAffi=%d \n", affineAMVPInfo.numCand );
}

void InterSearch::xCopyAffineAMVPInfo (AffineAMVPInfo& src, AffineAMVPInfo& dst)
{
  dst.numCand = src.numCand;
  DTRACE( g_trace_ctx, D_COMMON, " (%d) #copyAffi=%d \n", DTRACE_GET_COUNTER( g_trace_ctx, D_COMMON ), src.numCand );
  ::memcpy( dst.mvCandLT, src.mvCandLT, sizeof(Mv)*src.numCand );
  ::memcpy( dst.mvCandRT, src.mvCandRT, sizeof(Mv)*src.numCand );
  ::memcpy( dst.mvCandLB, src.mvCandLB, sizeof(Mv)*src.numCand );
}


/**
* \brief Generate half-sample interpolated block
*
* \param pattern Reference picture ROI
* \param biPred    Flag indicating whether block is for biprediction
*/
void InterSearch::xExtDIFUpSamplingH(CPelBuf* pattern, bool useAltHpelIf)
{
  const ClpRng& clpRng = m_lumaClpRng;
  int width      = pattern->width;
  int height     = pattern->height;
  int srcStride  = pattern->stride;

  int intStride = width + 1;
  int dstStride = width + 1;
  Pel *intPtr;
  Pel *dstPtr;
#if IF_12TAP
  int filterSize = NTAPS_LUMA(0);
#else
  int filterSize = NTAPS_LUMA;
#endif 
  int halfFilterSize = (filterSize>>1);
  const Pel *srcPtr = pattern->buf - halfFilterSize*srcStride - 1;

  const ChromaFormat chFmt = m_currChromaFormat;

  m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, m_filteredBlockTmp[0][0], intStride, width + 1, height + filterSize, 0 << MV_FRACTIONAL_BITS_DIFF, false, chFmt, clpRng, 0, false, useAltHpelIf);
  if (!m_skipFracME)
  {
    m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, m_filteredBlockTmp[2][0], intStride, width + 1, height + filterSize, 2 << MV_FRACTIONAL_BITS_DIFF, false, chFmt, clpRng, 0, false, useAltHpelIf);
  }

  intPtr = m_filteredBlockTmp[0][0] + halfFilterSize * intStride + 1;
  dstPtr = m_filteredBlock[0][0][0];
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 0, height + 0, 0 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng, 0, false, useAltHpelIf);
  if (m_skipFracME)
  {
    return;
  }

  intPtr = m_filteredBlockTmp[0][0] + (halfFilterSize - 1) * intStride + 1;
  dstPtr = m_filteredBlock[2][0][0];
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 0, height + 1, 2 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng, 0, false, useAltHpelIf);

  intPtr = m_filteredBlockTmp[2][0] + halfFilterSize * intStride;
  dstPtr = m_filteredBlock[0][2][0];
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 1, height + 0, 0 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng, 0, false, useAltHpelIf);

  intPtr = m_filteredBlockTmp[2][0] + (halfFilterSize - 1) * intStride;
  dstPtr = m_filteredBlock[2][2][0];
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width + 1, height + 1, 2 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng, 0, false, useAltHpelIf);
}





/**
* \brief Generate quarter-sample interpolated blocks
*
* \param pattern    Reference picture ROI
* \param halfPelRef Half-pel mv
* \param biPred     Flag indicating whether block is for biprediction
*/
void InterSearch::xExtDIFUpSamplingQ( CPelBuf* pattern, Mv halfPelRef)
{
  const ClpRng& clpRng = m_lumaClpRng;
  int width      = pattern->width;
  int height     = pattern->height;
  int srcStride  = pattern->stride;

  Pel const* srcPtr;
  int intStride = width + 1;
  int dstStride = width + 1;
  Pel *intPtr;
  Pel *dstPtr;
#if IF_12TAP
  int filterSize = NTAPS_LUMA(0);
#else
  int filterSize = NTAPS_LUMA;
#endif

  int halfFilterSize = (filterSize>>1);

  int extHeight = (halfPelRef.getVer() == 0) ? height + filterSize : height + filterSize-1;

  const ChromaFormat chFmt = m_currChromaFormat;

  // Horizontal filter 1/4
  srcPtr = pattern->buf - halfFilterSize * srcStride - 1;
  intPtr = m_filteredBlockTmp[1][0];
  if (halfPelRef.getVer() > 0)
  {
    srcPtr += srcStride;
  }
  if (halfPelRef.getHor() >= 0)
  {
    srcPtr += 1;
  }
  m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, intPtr, intStride, width, extHeight, 1 << MV_FRACTIONAL_BITS_DIFF, false, chFmt, clpRng);

  // Horizontal filter 3/4
  srcPtr = pattern->buf - halfFilterSize*srcStride - 1;
  intPtr = m_filteredBlockTmp[3][0];
  if (halfPelRef.getVer() > 0)
  {
    srcPtr += srcStride;
  }
  if (halfPelRef.getHor() > 0)
  {
    srcPtr += 1;
  }
  m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, intPtr, intStride, width, extHeight, 3 << MV_FRACTIONAL_BITS_DIFF, false, chFmt, clpRng);

  // Generate @ 1,1
  intPtr = m_filteredBlockTmp[1][0] + (halfFilterSize-1) * intStride;
  dstPtr = m_filteredBlock[1][1][0];
  if (halfPelRef.getVer() == 0)
  {
    intPtr += intStride;
  }
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);

  // Generate @ 3,1
  intPtr = m_filteredBlockTmp[1][0] + (halfFilterSize-1) * intStride;
  dstPtr = m_filteredBlock[3][1][0];
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);

  if (halfPelRef.getVer() != 0)
  {
    // Generate @ 2,1
    intPtr = m_filteredBlockTmp[1][0] + (halfFilterSize - 1) * intStride;
    dstPtr = m_filteredBlock[2][1][0];
    if (halfPelRef.getVer() == 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 2 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);

    // Generate @ 2,3
    intPtr = m_filteredBlockTmp[3][0] + (halfFilterSize - 1) * intStride;
    dstPtr = m_filteredBlock[2][3][0];
    if (halfPelRef.getVer() == 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 2 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);
  }
  else
  {
    // Generate @ 0,1
    intPtr = m_filteredBlockTmp[1][0] + halfFilterSize * intStride;
    dstPtr = m_filteredBlock[0][1][0];
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 0 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);

    // Generate @ 0,3
    intPtr = m_filteredBlockTmp[3][0] + halfFilterSize * intStride;
    dstPtr = m_filteredBlock[0][3][0];
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 0 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);
  }

  if (halfPelRef.getHor() != 0)
  {
    // Generate @ 1,2
    intPtr = m_filteredBlockTmp[2][0] + (halfFilterSize - 1) * intStride;
    dstPtr = m_filteredBlock[1][2][0];
    if (halfPelRef.getHor() > 0)
    {
      intPtr += 1;
    }
    if (halfPelRef.getVer() >= 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);

    // Generate @ 3,2
    intPtr = m_filteredBlockTmp[2][0] + (halfFilterSize - 1) * intStride;
    dstPtr = m_filteredBlock[3][2][0];
    if (halfPelRef.getHor() > 0)
    {
      intPtr += 1;
    }
    if (halfPelRef.getVer() > 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);
  }
  else
  {
    // Generate @ 1,0
    intPtr = m_filteredBlockTmp[0][0] + (halfFilterSize - 1) * intStride + 1;
    dstPtr = m_filteredBlock[1][0][0];
    if (halfPelRef.getVer() >= 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);

    // Generate @ 3,0
    intPtr = m_filteredBlockTmp[0][0] + (halfFilterSize - 1) * intStride + 1;
    dstPtr = m_filteredBlock[3][0][0];
    if (halfPelRef.getVer() > 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);
  }

  // Generate @ 1,3
  intPtr = m_filteredBlockTmp[3][0] + (halfFilterSize - 1) * intStride;
  dstPtr = m_filteredBlock[1][3][0];
  if (halfPelRef.getVer() == 0)
  {
    intPtr += intStride;
  }
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);

  // Generate @ 3,3
  intPtr = m_filteredBlockTmp[3][0] + (halfFilterSize - 1) * intStride;
  dstPtr = m_filteredBlock[3][3][0];
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3 << MV_FRACTIONAL_BITS_DIFF, false, true, chFmt, clpRng);
}





//! set wp tables
void InterSearch::setWpScalingDistParam( int iRefIdx, RefPicList eRefPicListCur, Slice *pcSlice )
{
  if ( iRefIdx<0 )
  {
    m_cDistParam.applyWeight = false;
    return;
  }

  WPScalingParam  *wp0 , *wp1;

  m_cDistParam.applyWeight = ( pcSlice->getSliceType()==P_SLICE && pcSlice->testWeightPred() ) || ( pcSlice->getSliceType()==B_SLICE && pcSlice->testWeightBiPred() ) ;

  if ( !m_cDistParam.applyWeight )
  {
    return;
  }

  int iRefIdx0 = ( eRefPicListCur == REF_PIC_LIST_0 ) ? iRefIdx : (-1);
  int iRefIdx1 = ( eRefPicListCur == REF_PIC_LIST_1 ) ? iRefIdx : (-1);

  getWpScaling( pcSlice, iRefIdx0, iRefIdx1, wp0 , wp1 );

  if ( iRefIdx0 < 0 )
  {
    wp0 = NULL;
  }
  if ( iRefIdx1 < 0 )
  {
    wp1 = NULL;
  }

  m_cDistParam.wpCur  = NULL;

  if ( eRefPicListCur == REF_PIC_LIST_0 )
  {
    m_cDistParam.wpCur = wp0;
  }
  else
  {
    m_cDistParam.wpCur = wp1;
  }
}

void InterSearch::xEncodeInterResidualQT(CodingStructure &cs, Partitioner &partitioner, const ComponentID &compID)
{
  const UnitArea& currArea    = partitioner.currArea();
  const TransformUnit &currTU = *cs.getTU(isLuma(partitioner.chType) ? currArea.lumaPos() : currArea.chromaPos(), partitioner.chType);
  const CodingUnit &cu        = *currTU.cu;
  const unsigned currDepth    = partitioner.currTrDepth;

  const bool bSubdiv          = currDepth != currTU.depth;

  if (compID == MAX_NUM_TBLOCKS)  // we are not processing a channel, instead we always recurse and code the CBFs
  {
    if( partitioner.canSplit( TU_MAX_TR_SPLIT, cs ) )
    {
      CHECK( !bSubdiv, "Not performing the implicit TU split" );
    }
    else if( cu.sbtInfo && partitioner.canSplit( PartSplit( cu.getSbtTuSplit() ), cs ) )
    {
      CHECK( !bSubdiv, "Not performing the implicit TU split - sbt" );
    }
    else
    {
      CHECK( bSubdiv, "transformsplit not supported" );
    }

    CHECK(CU::isIntra(cu), "Inter search provided with intra CU");

    if( cu.chromaFormat != CHROMA_400
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      && (!cu.isSepTree() || isChroma(partitioner.chType))
#else
      && (!CS::isDualITree(cs) || isChroma(partitioner.chType))
#endif
      )
    {
      {
        {
          const bool  chroma_cbf = TU::getCbfAtDepth( currTU, COMPONENT_Cb, currDepth );
          if (!(cu.sbtInfo && (currDepth == 0 || (currDepth == 1 && currTU.noResidual))))
          m_CABACEstimator->cbf_comp( cs, chroma_cbf, currArea.blocks[COMPONENT_Cb], currDepth );
        }
        {
          const bool  chroma_cbf = TU::getCbfAtDepth( currTU, COMPONENT_Cr, currDepth );
          if (!(cu.sbtInfo && (currDepth == 0 || (currDepth == 1 && currTU.noResidual))))
          m_CABACEstimator->cbf_comp( cs, chroma_cbf, currArea.blocks[COMPONENT_Cr], currDepth, TU::getCbfAtDepth( currTU, COMPONENT_Cb, currDepth ) );
        }
      }
    }

    if( !bSubdiv && !( cu.sbtInfo && currTU.noResidual )
      && !isChroma(partitioner.chType)
      )
    {
      m_CABACEstimator->cbf_comp( cs, TU::getCbfAtDepth( currTU, COMPONENT_Y, currDepth ), currArea.Y(), currDepth );
    }
  }

  if (!bSubdiv)
  {
    if (compID != MAX_NUM_TBLOCKS) // we have already coded the CBFs, so now we code coefficients
    {
      if( currArea.blocks[compID].valid() )
      {
        if( compID == COMPONENT_Cr )
        {
          const int cbfMask = ( TU::getCbf( currTU, COMPONENT_Cb ) ? 2 : 0) + ( TU::getCbf( currTU, COMPONENT_Cr ) ? 1 : 0 );
          m_CABACEstimator->joint_cb_cr( currTU, cbfMask );
        }
        if( TU::getCbf( currTU, compID ) )
        {
          m_CABACEstimator->residual_coding( currTU, compID );
        }
      }
    }
  }
  else
  {
    if( compID == MAX_NUM_TBLOCKS || TU::getCbfAtDepth( currTU, compID, currDepth ) )
    {
      if( partitioner.canSplit( TU_MAX_TR_SPLIT, cs ) )
      {
        partitioner.splitCurrArea( TU_MAX_TR_SPLIT, cs );
      }
      else if( cu.sbtInfo && partitioner.canSplit( PartSplit( cu.getSbtTuSplit() ), cs ) )
      {
        partitioner.splitCurrArea( PartSplit( cu.getSbtTuSplit() ), cs );
      }
      else
        THROW( "Implicit TU split not available!" );

      do
      {
        xEncodeInterResidualQT( cs, partitioner, compID );
      } while( partitioner.nextPart( cs ) );

      partitioner.exitCurrSplit();
    }
  }
}

void InterSearch::calcMinDistSbt( CodingStructure &cs, const CodingUnit& cu, const uint8_t sbtAllowed )
{
  if( !sbtAllowed )
  {
    m_estMinDistSbt[NUMBER_SBT_MODE] = 0;
    for( int comp = 0; comp < getNumberValidTBlocks( *cs.pcv ); comp++ )
    {
      const ComponentID compID = ComponentID( comp );
      CPelBuf pred = cs.getPredBuf( compID );
      CPelBuf org  = cs.getOrgBuf( compID );
      m_estMinDistSbt[NUMBER_SBT_MODE] += m_pcRdCost->getDistPart( org, pred, cs.sps->getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
    }
    return;
  }

  //SBT fast algorithm 2.1 : estimate a minimum RD cost of a SBT mode based on the luma distortion of uncoded part and coded part (assuming distorted can be reduced to 1/16);
  //                         if this cost is larger than the best cost, no need to try a specific SBT mode
  int cuWidth  = cu.lwidth();
  int cuHeight = cu.lheight();
  int numPartX = cuWidth  >= 16 ? 4 : ( cuWidth  == 4 ? 1 : 2 );
  int numPartY = cuHeight >= 16 ? 4 : ( cuHeight == 4 ? 1 : 2 );
  Distortion dist[4][4];
  memset( dist, 0, sizeof( Distortion ) * 16 );

  for( uint32_t c = 0; c < getNumberValidTBlocks( *cs.pcv ); c++ )
  {
    const ComponentID compID   = ComponentID( c );
    const CompArea&   compArea = cu.blocks[compID];
    const CPelBuf orgPel  = cs.getOrgBuf( compArea );
    const CPelBuf predPel = cs.getPredBuf( compArea );
    int lengthX = compArea.width / numPartX;
    int lengthY = compArea.height / numPartY;
    int strideOrg  = orgPel.stride;
    int stridePred = predPel.stride;
    uint32_t   uiShift = DISTORTION_PRECISION_ADJUSTMENT( ( *cs.sps.getBitDepth( toChannelType( compID ) ) - 8 ) << 1 );
    Intermediate_Int iTemp;

    //calc distY of 16 sub parts
    for( int j = 0; j < numPartY; j++ )
    {
      for( int i = 0; i < numPartX; i++ )
      {
        int posX = i * lengthX;
        int posY = j * lengthY;
        const Pel* ptrOrg  = orgPel.bufAt( posX, posY );
        const Pel* ptrPred = predPel.bufAt( posX, posY );
        Distortion uiSum = 0;
        for( int n = 0; n < lengthY; n++ )
        {
          for( int m = 0; m < lengthX; m++ )
          {
            iTemp = ptrOrg[m] - ptrPred[m];
            uiSum += Distortion( ( iTemp * iTemp ) >> uiShift );
          }
          ptrOrg += strideOrg;
          ptrPred += stridePred;
        }
        if( isChroma( compID ) )
        {
          uiSum = (Distortion)( uiSum * m_pcRdCost->getChromaWeight() );
        }
        dist[j][i] += uiSum;
      }
    }
  }

  //SSE of a CU
  m_estMinDistSbt[NUMBER_SBT_MODE] = 0;
  for( int j = 0; j < numPartY; j++ )
  {
    for( int i = 0; i < numPartX; i++ )
    {
      m_estMinDistSbt[NUMBER_SBT_MODE] += dist[j][i];
    }
  }
  //init per-mode dist
  for( int i = SBT_VER_H0; i < NUMBER_SBT_MODE; i++ )
  {
    m_estMinDistSbt[i] = std::numeric_limits<uint64_t>::max();
  }

  //SBT fast algorithm 1: not try SBT if the residual is too small to compensate bits for encoding residual info
  uint64_t minNonZeroResiFracBits = 12 << SCALE_BITS;
  if( m_pcRdCost->calcRdCost( 0, m_estMinDistSbt[NUMBER_SBT_MODE] ) < m_pcRdCost->calcRdCost( minNonZeroResiFracBits, 0 ) )
  {
    m_skipSbtAll = true;
    return;
  }

  //derive estimated minDist of SBT = zero-residual part distortion + non-zero residual part distortion / 16
  int shift = 5;
  Distortion distResiPart = 0, distNoResiPart = 0;

  if( CU::targetSbtAllowed( SBT_VER_HALF, sbtAllowed ) )
  {
    int offsetResiPart = 0;
    int offsetNoResiPart = numPartX / 2;
    distResiPart = distNoResiPart = 0;
    assert( numPartX >= 2 );
    for( int j = 0; j < numPartY; j++ )
    {
      for( int i = 0; i < numPartX / 2; i++ )
      {
        distResiPart   += dist[j][i + offsetResiPart];
        distNoResiPart += dist[j][i + offsetNoResiPart];
      }
    }
    m_estMinDistSbt[SBT_VER_H0] = ( distResiPart >> shift ) + distNoResiPart;
    m_estMinDistSbt[SBT_VER_H1] = ( distNoResiPart >> shift ) + distResiPart;
  }

  if( CU::targetSbtAllowed( SBT_HOR_HALF, sbtAllowed ) )
  {
    int offsetResiPart = 0;
    int offsetNoResiPart = numPartY / 2;
    assert( numPartY >= 2 );
    distResiPart = distNoResiPart = 0;
    for( int j = 0; j < numPartY / 2; j++ )
    {
      for( int i = 0; i < numPartX; i++ )
      {
        distResiPart   += dist[j + offsetResiPart][i];
        distNoResiPart += dist[j + offsetNoResiPart][i];
      }
    }
    m_estMinDistSbt[SBT_HOR_H0] = ( distResiPart >> shift ) + distNoResiPart;
    m_estMinDistSbt[SBT_HOR_H1] = ( distNoResiPart >> shift ) + distResiPart;
  }

  if( CU::targetSbtAllowed( SBT_VER_QUAD, sbtAllowed ) )
  {
    assert( numPartX == 4 );
    m_estMinDistSbt[SBT_VER_Q0] = m_estMinDistSbt[SBT_VER_Q1] = 0;
    for( int j = 0; j < numPartY; j++ )
    {
      m_estMinDistSbt[SBT_VER_Q0] += dist[j][0] + ( ( dist[j][1] + dist[j][2] + dist[j][3] ) << shift );
      m_estMinDistSbt[SBT_VER_Q1] += dist[j][3] + ( ( dist[j][0] + dist[j][1] + dist[j][2] ) << shift );
    }
    m_estMinDistSbt[SBT_VER_Q0] = m_estMinDistSbt[SBT_VER_Q0] >> shift;
    m_estMinDistSbt[SBT_VER_Q1] = m_estMinDistSbt[SBT_VER_Q1] >> shift;
  }

  if( CU::targetSbtAllowed( SBT_HOR_QUAD, sbtAllowed ) )
  {
    assert( numPartY == 4 );
    m_estMinDistSbt[SBT_HOR_Q0] = m_estMinDistSbt[SBT_HOR_Q1] = 0;
    for( int i = 0; i < numPartX; i++ )
    {
      m_estMinDistSbt[SBT_HOR_Q0] += dist[0][i] + ( ( dist[1][i] + dist[2][i] + dist[3][i] ) << shift );
      m_estMinDistSbt[SBT_HOR_Q1] += dist[3][i] + ( ( dist[0][i] + dist[1][i] + dist[2][i] ) << shift );
    }
    m_estMinDistSbt[SBT_HOR_Q0] = m_estMinDistSbt[SBT_HOR_Q0] >> shift;
    m_estMinDistSbt[SBT_HOR_Q1] = m_estMinDistSbt[SBT_HOR_Q1] >> shift;
  }

  //SBT fast algorithm 5: try N SBT modes with the lowest distortion
  Distortion temp[NUMBER_SBT_MODE];
  memcpy( temp, m_estMinDistSbt, sizeof( Distortion ) * NUMBER_SBT_MODE );
  memset( m_sbtRdoOrder, 255, NUMBER_SBT_MODE );
  int startIdx = 0, numRDO;
  numRDO = CU::targetSbtAllowed( SBT_VER_HALF, sbtAllowed ) + CU::targetSbtAllowed( SBT_HOR_HALF, sbtAllowed );
  numRDO = std::min( ( numRDO << 1 ), SBT_NUM_RDO );
  for( int i = startIdx; i < startIdx + numRDO; i++ )
  {
    Distortion minDist = std::numeric_limits<uint64_t>::max();
    for( int n = SBT_VER_H0; n <= SBT_HOR_H1; n++ )
    {
      if( temp[n] < minDist )
      {
        minDist = temp[n];
        m_sbtRdoOrder[i] = n;
      }
    }
    temp[m_sbtRdoOrder[i]] = std::numeric_limits<uint64_t>::max();
  }

  startIdx += numRDO;
  numRDO = CU::targetSbtAllowed( SBT_VER_QUAD, sbtAllowed ) + CU::targetSbtAllowed( SBT_HOR_QUAD, sbtAllowed );
  numRDO = std::min( ( numRDO << 1 ), SBT_NUM_RDO );
  for( int i = startIdx; i < startIdx + numRDO; i++ )
  {
    Distortion minDist = std::numeric_limits<uint64_t>::max();
    for( int n = SBT_VER_Q0; n <= SBT_HOR_Q1; n++ )
    {
      if( temp[n] < minDist )
      {
        minDist = temp[n];
        m_sbtRdoOrder[i] = n;
      }
    }
    temp[m_sbtRdoOrder[i]] = std::numeric_limits<uint64_t>::max();
  }
}

uint8_t InterSearch::skipSbtByRDCost( int width, int height, int mtDepth, uint8_t sbtIdx, uint8_t sbtPos, double bestCost, Distortion distSbtOff, double costSbtOff, bool rootCbfSbtOff )
{
  int sbtMode = CU::getSbtMode( sbtIdx, sbtPos );

  //SBT fast algorithm 2.2 : estimate a minimum RD cost of a SBT mode based on the luma distortion of uncoded part and coded part (assuming distorted can be reduced to 1/16);
  //                         if this cost is larger than the best cost, no need to try a specific SBT mode
  if( m_pcRdCost->calcRdCost( 11 << SCALE_BITS, m_estMinDistSbt[sbtMode] ) > bestCost )
  {
    return 0; //early skip type 0
  }

  if( costSbtOff != MAX_DOUBLE )
  {
    if( !rootCbfSbtOff )
    {
      //SBT fast algorithm 3: skip SBT when the residual is too small (estCost is more accurate than fast algorithm 1, counting PU mode bits)
      uint64_t minNonZeroResiFracBits = 10 << SCALE_BITS;
      Distortion distResiPart;
      if( sbtIdx == SBT_VER_HALF || sbtIdx == SBT_HOR_HALF )
      {
        distResiPart = (Distortion)( ( ( m_estMinDistSbt[NUMBER_SBT_MODE] - m_estMinDistSbt[sbtMode] ) * 9 ) >> 4 );
      }
      else
      {
        distResiPart = (Distortion)( ( ( m_estMinDistSbt[NUMBER_SBT_MODE] - m_estMinDistSbt[sbtMode] ) * 3 ) >> 3 );
      }

      double estCost = ( costSbtOff - m_pcRdCost->calcRdCost( 0 << SCALE_BITS, distSbtOff ) ) + m_pcRdCost->calcRdCost( minNonZeroResiFracBits, m_estMinDistSbt[sbtMode] + distResiPart );
      if( estCost > costSbtOff )
      {
        return 1;
      }
      if( estCost > bestCost )
      {
        return 2;
      }
    }
    else
    {
      //SBT fast algorithm 4: skip SBT when an estimated RD cost is larger than the bestCost
      double weight = sbtMode > SBT_HOR_H1 ? 0.4 : 0.6;
      double estCost = ( ( costSbtOff - m_pcRdCost->calcRdCost( 0 << SCALE_BITS, distSbtOff ) ) * weight ) + m_pcRdCost->calcRdCost( 0 << SCALE_BITS, m_estMinDistSbt[sbtMode] );
      if( estCost > bestCost )
      {
        return 3;
      }
    }
  }
  return MAX_UCHAR;
}

void InterSearch::xEstimateInterResidualQT(CodingStructure &cs, Partitioner &partitioner, Distortion *puiZeroDist /*= NULL*/
  , const bool luma, const bool chroma
  , PelUnitBuf* orgResi
)
{
  const UnitArea& currArea = partitioner.currArea();
  const SPS &sps           = *cs.sps;
  m_pcRdCost->setChromaFormat(sps.getChromaFormatIdc());

  const uint32_t numValidComp  = getNumberValidComponents( sps.getChromaFormatIdc() );
  const uint32_t numTBlocks    = getNumberValidTBlocks   ( *cs.pcv );
  const CodingUnit &cu = *cs.getCU(partitioner.chType);
  const unsigned currDepth = partitioner.currTrDepth;
  const bool colorTransFlag = cs.cus[0]->colorTransform;

  bool bCheckFull  = !partitioner.canSplit( TU_MAX_TR_SPLIT, cs );
  if( cu.sbtInfo && partitioner.canSplit( PartSplit( cu.getSbtTuSplit() ), cs ) )
  {
    bCheckFull = false;
  }
  bool bCheckSplit = !bCheckFull;

  // get temporary data
  CodingStructure *csSplit = nullptr;
  CodingStructure *csFull  = nullptr;
  if (bCheckSplit)
  {
    csSplit = &cs;
  }
  else if (bCheckFull)
  {
    csFull = &cs;
  }

  Distortion uiSingleDist         = 0;
  Distortion uiSingleDistComp [3] = { 0, 0, 0 };
  uint64_t   uiSingleFracBits[3] = { 0, 0, 0 };

  const TempCtx ctxStart  ( m_CtxCache, m_CABACEstimator->getCtx() );
  TempCtx       ctxBest   ( m_CtxCache );

  if (bCheckFull)
  {
    TransformUnit &tu = csFull->addTU(CS::getArea(cs, currArea, partitioner.chType), partitioner.chType);
    tu.depth          = currDepth;
    for (int i = 0; i<MAX_NUM_TBLOCKS; i++) tu.mtsIdx[i] = MTS_DCT2_DCT2;
    tu.checkTuNoResidual( partitioner.currPartIdx() );
    Position tuPos = tu.Y();
    tuPos.relativeTo(cu.Y());
    const UnitArea relativeUnitArea(tu.chromaFormat, Area(tuPos, tu.Y().size()));

    const Slice           &slice = *cs.slice;
    if (slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && !(CS::isDualITree(cs) && slice.isIntra() && tu.cu->predMode == MODE_IBC))
    {
#if LMCS_CHROMA_CALC_CU
      const CompArea      &areaY = tu.cu->blocks[COMPONENT_Y];
#else
      const CompArea      &areaY = tu.blocks[COMPONENT_Y];
#endif
      int adj = m_pcReshape->calculateChromaAdjVpduNei(tu, areaY);
      tu.setChromaAdj(adj);
    }

#if JVET_S0234_ACT_CRS_FIX
    PelUnitBuf colorTransResidual = m_colorTransResiBuf[1].getBuf(relativeUnitArea);
    if (colorTransFlag)
    {
      csFull->getResiBuf(currArea).copyFrom(cs.getOrgResiBuf(currArea));
      if (slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[COMPONENT_Cb].width*tu.blocks[COMPONENT_Cr].height > 4)
      {
        csFull->getResiBuf(currArea).bufs[1].scaleSignal(tu.getChromaAdj(), 1, tu.cu->cs->slice->clpRng(COMPONENT_Cb));
        csFull->getResiBuf(currArea).bufs[2].scaleSignal(tu.getChromaAdj(), 1, tu.cu->cs->slice->clpRng(COMPONENT_Cr));
      }
      csFull->getResiBuf(currArea).colorSpaceConvert(colorTransResidual, true, cu.cs->slice->clpRng(COMPONENT_Y));
    }
#endif
    double minCost            [MAX_NUM_TBLOCKS];

    m_CABACEstimator->resetBits();

    memset(m_pTempPel, 0, sizeof(Pel) * tu.Y().area()); // not necessary needed for inside of recursion (only at the beginning)

    for (uint32_t i = 0; i < numTBlocks; i++)
    {
      minCost[i] = MAX_DOUBLE;
    }

    CodingStructure &saveCS = *m_pSaveCS[0];
    saveCS.pcv     = cs.pcv;
    saveCS.picture = cs.picture;
#if JVET_Z0118_GDR
    saveCS.m_pt = cs.m_pt;
#endif
    saveCS.area.repositionTo(currArea);
    saveCS.clearTUs();
    TransformUnit & bestTU = saveCS.addTU(CS::getArea(cs, currArea, partitioner.chType), partitioner.chType);

    for( uint32_t c = 0; c < numTBlocks; c++ )
    {
      const ComponentID compID    = ComponentID(c);
      if (compID == COMPONENT_Y && !luma)
        continue;
      if (compID != COMPONENT_Y && !chroma)
        continue;
      const CompArea&   compArea  = tu.blocks[compID];
      const int channelBitDepth   = sps.getBitDepth(toChannelType(compID));

      if( !tu.blocks[compID].valid() )
      {
        continue;
      }


      const bool tsAllowed  = TU::isTSAllowed(tu, compID) && (isLuma(compID) || (isChroma(compID) && m_pcEncCfg->getUseChromaTS()));
      const bool mtsAllowed = CU::isMTSAllowed( *tu.cu, compID );

      uint8_t nNumTransformCands = 1 + ( tsAllowed ? 1 : 0 ) + ( mtsAllowed ? 4 : 0 ); // DCT + TS + 4 MTS = 6 tests
      std::vector<TrMode> trModes;
#if TU_256
      if(tu.idx != cu.firstTU->idx)
      {
        trModes.push_back( TrMode( cu.firstTU->mtsIdx[compID], true ) );
        nNumTransformCands = 1;
      }
      else
      {
#endif
      if (m_pcEncCfg->getCostMode() == COST_LOSSLESS_CODING && slice.isLossless())
      {
        nNumTransformCands = 0;
      }
      else
      {
      trModes.push_back( TrMode( 0, true ) ); //DCT2
      nNumTransformCands = 1;
      }
      //for a SBT-no-residual TU, the RDO process should be called once, in order to get the RD cost
      if( tsAllowed && !tu.noResidual )
      {
        trModes.push_back( TrMode( 1, true ) );
        nNumTransformCands++;
      }

#if APPLY_SBT_SL_ON_MTS
      //skip MTS if DCT2 is the best
      if( mtsAllowed && ( !tu.cu->slice->getSPS()->getUseSBT() || CU::getSbtIdx( m_histBestSbt ) != SBT_OFF_DCT ) )
#else
      if( mtsAllowed )
#endif
      {
        for( int i = 2; i < 6; i++ )
        {
#if APPLY_SBT_SL_ON_MTS
          //skip the non-best Mts mode
          if( !tu.cu->slice->getSPS()->getUseSBT() || ( m_histBestMtsIdx == MAX_UCHAR || m_histBestMtsIdx == i ) )
          {
#endif
          trModes.push_back( TrMode( i, true ) );
          nNumTransformCands++;
#if APPLY_SBT_SL_ON_MTS
          }
#endif
        }
      }
#if TU_256
      }
#endif

      if (colorTransFlag && (m_pcEncCfg->getCostMode() != COST_LOSSLESS_CODING || !slice.isLossless()))
      {
        m_pcTrQuant->lambdaAdjustColorTrans(true);
#if JVET_S0234_ACT_CRS_FIX
        if (isChroma(compID) && slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[compID].width*tu.blocks[compID].height > 4)
        {
          int cResScaleInv = tu.getChromaAdj();
          m_pcRdCost->lambdaAdjustColorTrans(true, compID, true, &cResScaleInv);
        }
        else
#endif
        m_pcRdCost->lambdaAdjustColorTrans(true, compID);
      }

      const int numTransformCandidates = nNumTransformCands;
      for( int transformMode = 0; transformMode < numTransformCandidates; transformMode++ )
      {
          const bool isFirstMode  = transformMode == 0;
          // copy the original residual into the residual buffer
#if JVET_S0234_ACT_CRS_FIX
          if (colorTransFlag)
          {
            csFull->getResiBuf(compArea).copyFrom(colorTransResidual.bufs[compID]);
          }
          else
#endif
          csFull->getResiBuf(compArea).copyFrom(cs.getOrgResiBuf(compArea));

          m_CABACEstimator->getCtx() = ctxStart;
          m_CABACEstimator->resetBits();

          {
            if (!(m_pcEncCfg->getCostMode() == COST_LOSSLESS_CODING && slice.isLossless()))
            {
            if (bestTU.mtsIdx[compID] == MTS_SKIP && m_pcEncCfg->getUseTransformSkipFast())
            {
              continue;
            }
            if( !trModes[transformMode].second )
            {
              continue;
            }
            }
            tu.mtsIdx[compID] = trModes[transformMode].first;
          }
          QpParam cQP(tu, compID);  // note: uses tu.transformSkip[compID]

#if RDOQ_CHROMA_LAMBDA
          m_pcTrQuant->selectLambda(compID);
#endif
          if (slice.getLmcsEnabledFlag() && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag())
          {
            double cRescale = (double)(1 << CSCALE_FP_PREC) / (double)(tu.getChromaAdj());
            m_pcTrQuant->setLambda(m_pcTrQuant->getLambda() / (cRescale*cRescale));
          }
          if ( sps.getJointCbCrEnabledFlag() && isChroma( compID ) && ( tu.cu->cs->slice->getSliceQp() > 18 ) )
          {
            m_pcTrQuant->setLambda( 1.05 * m_pcTrQuant->getLambda() );
          }

          TCoeff     currAbsSum = 0;
          uint64_t   currCompFracBits = 0;
          Distortion currCompDist = 0;
          double     currCompCost = 0;
          uint64_t   nonCoeffFracBits = 0;
          Distortion nonCoeffDist = 0;
          double     nonCoeffCost = 0;

#if JVET_S0234_ACT_CRS_FIX 
          if (!colorTransFlag && slice.getLmcsEnabledFlag() && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[compID].width*tu.blocks[compID].height > 4)
#else
          if (slice.getLmcsEnabledFlag() && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[compID].width * tu.blocks[compID].height > 4)
#endif
          {
            PelBuf resiBuf = csFull->getResiBuf(compArea);
            resiBuf.scaleSignal(tu.getChromaAdj(), 1, tu.cu->cs->slice->clpRng(compID));
          }
          if( nNumTransformCands > 1 )
          {
            if( transformMode == 0 )
            {
              m_pcTrQuant->transformNxN( tu, compID, cQP, &trModes, m_pcEncCfg->getMTSInterMaxCand() );
              tu.mtsIdx[compID] = trModes[0].first;
            }
            if (!(m_pcEncCfg->getCostMode() == COST_LOSSLESS_CODING && slice.isLossless() && tu.mtsIdx[compID] == 0))
            {
              m_pcTrQuant->transformNxN( tu, compID, cQP, currAbsSum, m_CABACEstimator->getCtx(), true );
            }
          }
          else
          {
            m_pcTrQuant->transformNxN( tu, compID, cQP, currAbsSum, m_CABACEstimator->getCtx() );
          }

          if (isFirstMode || (currAbsSum == 0))
          {
            const CPelBuf zeroBuf(m_pTempPel, compArea);
#if JVET_S0234_ACT_CRS_FIX
            const CPelBuf orgResi = colorTransFlag ? colorTransResidual.bufs[compID] : csFull->getOrgResiBuf(compArea);
#else
            const CPelBuf orgResi = csFull->getOrgResiBuf( compArea );
#endif

            {
              nonCoeffDist = m_pcRdCost->getDistPart( zeroBuf, orgResi, channelBitDepth, compID, DF_SSE ); // initialized with zero residual distortion
            }

            if( !tu.noResidual )
            {
            const bool prevCbf = ( compID == COMPONENT_Cr ? tu.cbf[COMPONENT_Cb] : false );
            m_CABACEstimator->cbf_comp( *csFull, false, compArea, currDepth, prevCbf );

            }

            nonCoeffFracBits = m_CABACEstimator->getEstFracBits();
#if WCG_EXT
            if( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() )
            {
              nonCoeffCost   = m_pcRdCost->calcRdCost(nonCoeffFracBits, nonCoeffDist, false);
            }
            else
#endif
              if (cs.slice->getSPS()->getUseColorTrans())
              {
                nonCoeffCost = m_pcRdCost->calcRdCost(nonCoeffFracBits, nonCoeffDist, false);
              }
              else
              {
                nonCoeffCost = m_pcRdCost->calcRdCost(nonCoeffFracBits, nonCoeffDist);
              }
          }

          if ((puiZeroDist != NULL) && isFirstMode)
          {
            *puiZeroDist += nonCoeffDist; // initialized with zero residual distortion
          }
          if (m_pcEncCfg->getCostMode() == COST_LOSSLESS_CODING && slice.isLossless() && tu.mtsIdx[compID] == 0)
          {
            currAbsSum = 0;
          }

          if (currAbsSum > 0) //if non-zero coefficients are present, a residual needs to be derived for further prediction
          {
            if (isFirstMode)
            {
              m_CABACEstimator->getCtx() = ctxStart;
              m_CABACEstimator->resetBits();
            }

            const bool prevCbf = ( compID == COMPONENT_Cr ? tu.cbf[COMPONENT_Cb] : false );
            m_CABACEstimator->cbf_comp( *csFull, true, compArea, currDepth, prevCbf );
            if( compID == COMPONENT_Cr )
            {
              const int cbfMask = ( tu.cbf[COMPONENT_Cb] ? 2 : 0 ) + 1;
              m_CABACEstimator->joint_cb_cr( tu, cbfMask );
            }

#if SIGN_PREDICTION
            if ( sps.getNumPredSigns() > 0)
            {
#if JVET_Y0141_SIGN_PRED_IMPROVE
              bool doSignPrediction = true;
              if (isLuma(compID) && tu.mtsIdx[COMPONENT_Y] > MTS_SKIP)
              {
                bool signHiding = slice.getSignDataHidingEnabledFlag();
                CoeffCodingContext  cctx(tu, COMPONENT_Y, signHiding);
                int scanPosLast = -1;
                TCoeff* coeff = tu.getCoeffs(compID).buf;
                for (int scanPos = cctx.maxNumCoeff() - 1; scanPos >= 0; scanPos--)
                {
                  unsigned blkPos = cctx.blockPos(scanPos);
                  if (coeff[blkPos])
                  {
                    scanPosLast = scanPos;
                    break;
                  }
                }
                if (scanPosLast < 1)
                {
                  doSignPrediction = false;
                }
              }
              if (doSignPrediction)
              {
#endif
                bool reshapeChroma = slice.getPicHeader()->getLmcsEnabledFlag() && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[compID].width*tu.blocks[compID].height > 4;
#if JVET_Y0065_GPM_INTRA
                if (isLuma(compID) && slice.getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !tu.cu->firstPU->ciipFlag && !tu.cu->firstPU->gpmIntraFlag && !CU::isIBC(*tu.cu))
#else
			    if (isLuma(compID) && slice.getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !tu.cu->firstPU->ciipFlag && !CU::isIBC(*tu.cu))
#endif
                {
#if JVET_Z0118_GDR
                  cs.updateReconMotIPM(tu.blocks[COMPONENT_Y], cs.getPredBuf(tu.blocks[COMPONENT_Y]));
#else
                  cs.picture->getRecoBuf(tu.blocks[COMPONENT_Y]).copyFrom(cs.getPredBuf(tu.blocks[COMPONENT_Y]));
#endif
                  cs.getPredBuf(tu.blocks[compID]).rspSignal(m_pcReshape->getFwdLUT());
                }
                m_pcTrQuant->predCoeffSigns(tu, compID, reshapeChroma);
#if JVET_Y0065_GPM_INTRA
                if (isLuma(compID) && slice.getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !tu.cu->firstPU->ciipFlag && !tu.cu->firstPU->gpmIntraFlag && !CU::isIBC(*tu.cu))
#else
			    if (isLuma(compID) && slice.getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !tu.cu->firstPU->ciipFlag && !CU::isIBC(*tu.cu))
#endif
                {
                  cs.getPredBuf(tu.blocks[COMPONENT_Y]).copyFrom(cs.picture->getRecoBuf(tu.blocks[COMPONENT_Y]));
                }
#if JVET_Y0141_SIGN_PRED_IMPROVE
              }
#endif
            }
#endif

            CUCtx cuCtx;
            cuCtx.isDQPCoded = true;
            cuCtx.isChromaQpAdjCoded = true;
            m_CABACEstimator->residual_coding(tu, compID, &cuCtx);
            m_CABACEstimator->mts_idx(cu, &cuCtx);

            if (compID == COMPONENT_Y && tu.mtsIdx[compID] > MTS_SKIP && !cuCtx.mtsLastScanPos)
            {
              currCompCost = MAX_DOUBLE;
            }
            else
            {

            currCompFracBits = m_CABACEstimator->getEstFracBits();

            PelBuf resiBuf = csFull->getResiBuf(compArea);
#if JVET_S0234_ACT_CRS_FIX
            CPelBuf orgResiBuf = colorTransFlag ? colorTransResidual.bufs[compID] : csFull->getOrgResiBuf(compArea);
#else
            CPelBuf orgResiBuf = csFull->getOrgResiBuf(compArea);
#endif

            m_pcTrQuant->invTransformNxN(tu, compID, resiBuf, cQP);
#if JVET_S0234_ACT_CRS_FIX
            if (!colorTransFlag && slice.getLmcsEnabledFlag() && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[compID].width*tu.blocks[compID].height > 4)
#else
            if (slice.getLmcsEnabledFlag() && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[compID].width * tu.blocks[compID].height > 4)
#endif
            {
              resiBuf.scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(compID));
            }

#if JVET_V0094_BILATERAL_FILTER
            // getCbf() is going to be 1 since currAbsSum > 0 here, according to the if-statement a couple of lines up.
            bool isInter = (cu.predMode == MODE_INTER) ? true : false;
            if (cs.pps->getUseBIF() && isLuma(compID) && (tu.cu->qp > 17) && (128 > std::max(tu.lumaSize().width, tu.lumaSize().height)) && ((isInter == false) || (32 > std::min(tu.lumaSize().width, tu.lumaSize().height))))
            {
              CompArea tmpArea1(COMPONENT_Y, tu.chromaFormat, Position(0, 0), Size(resiBuf.width, resiBuf.height));
              PelBuf tmpRecLuma = m_tmpStorageLCU.getBuf(tmpArea1);
              tmpRecLuma.copyFrom(resiBuf);
              
              const CPelBuf predBuf = csFull->getPredBuf(compArea);
              PelBuf recIPredBuf = csFull->slice->getPic()->getRecoBuf(compArea);
              std::vector<Pel> invLUT;
              m_bilateralFilter->bilateralFilterRDOdiamond5x5(tmpRecLuma, predBuf, tmpRecLuma, tu.cu->qp, recIPredBuf, cs.slice->clpRng(compID), tu, false, false, invLUT);

              currCompDist = m_pcRdCost->getDistPart(orgResiBuf, tmpRecLuma, channelBitDepth, compID, DF_SSE);
            }
            else
            {
#if JVET_X0071_CHROMA_BILATERAL_FILTER
              if(isChroma(compID))
              {
                if (cs.pps->getUseChromaBIF() && isChroma(compID) && (tu.cu->qp > 17))
                {
                  //chroma and bilateral
                  CompArea tmpArea1(compID, tu.chromaFormat, Position(0, 0), Size(resiBuf.width, resiBuf.height));
                  PelBuf tmpRecChroma = m_tmpStorageLCU.getBuf(tmpArea1);
                  tmpRecChroma.copyFrom(resiBuf);

                  const CPelBuf predBuf = csFull->getPredBuf(compArea);
                  PelBuf recIPredBuf = csFull->slice->getPic()->getRecoBuf(compArea);
                  bool isCb = compID == COMPONENT_Cb ? true : false;
                  m_bilateralFilter->bilateralFilterRDOdiamond5x5Chroma(tmpRecChroma, predBuf, tmpRecChroma, tu.cu->qp, recIPredBuf, cs.slice->clpRng(compID), tu, false, isCb);
                  currCompDist = m_pcRdCost->getDistPart(orgResiBuf, tmpRecChroma, channelBitDepth, compID, DF_SSE);
                }
                else
                {   //chroma but not bilateral
                  currCompDist = m_pcRdCost->getDistPart(orgResiBuf, resiBuf, channelBitDepth, compID, DF_SSE);
                }
              }
              else
              { //luma but not bilateral
                currCompDist = m_pcRdCost->getDistPart(orgResiBuf, resiBuf, channelBitDepth, compID, DF_SSE);
              }
#else
              currCompDist = m_pcRdCost->getDistPart(orgResiBuf, resiBuf, channelBitDepth, compID, DF_SSE);
#endif
            }
#else
#if JVET_X0071_CHROMA_BILATERAL_FILTER
            if(isChroma(compID))
            {
              if(cs.pps->getUseChromaBIF() && isChroma(compID) && (tu.cu->qp > 17))
              {
                //chroma and bilateral
                CompArea tmpArea1(compID, tu.chromaFormat, Position(0, 0), Size(resiBuf.width, resiBuf.height));
                PelBuf tmpRecChroma = m_tmpStorageLCU.getBuf(tmpArea1);
                tmpRecChroma.copyFrom(resiBuf);
                const CPelBuf predBuf = csFull->getPredBuf(compArea);
                PelBuf recIPredBuf = csFull->slice->getPic()->getRecoBuf(compArea);
                bool isCb = compID == COMPONENT_Cb ? true : false;
                m_bilateralFilter->bilateralFilterRDOdiamond5x5Chroma(tmpRecChroma, predBuf, tmpRecChroma, tu.cu->qp, recIPredBuf, cs.slice->clpRng(compID), tu, false, isCb);
                currCompDist = m_pcRdCost->getDistPart(orgResiBuf, tmpRecChroma, channelBitDepth, compID, DF_SSE);
              }
              else
              {
                //chroma but not bilateral
                currCompDist = m_pcRdCost->getDistPart(orgResiBuf, resiBuf, channelBitDepth, compID, DF_SSE);
              }
            }
            else
            {
              //luma but not bilateral
              currCompDist = m_pcRdCost->getDistPart(orgResiBuf, resiBuf, channelBitDepth, compID, DF_SSE);
            }
#else
            currCompDist = m_pcRdCost->getDistPart(orgResiBuf, resiBuf, channelBitDepth, compID, DF_SSE);
#endif
#endif

#if WCG_EXT
            currCompCost = m_pcRdCost->calcRdCost(currCompFracBits, currCompDist, false);
#else
            currCompCost = m_pcRdCost->calcRdCost(currCompFracBits, currCompDist);
#endif
            }
          }
          else if( transformMode > 0 )
          {
            currCompCost = MAX_DOUBLE;
          }
          else
          {
            currCompFracBits = nonCoeffFracBits;
            currCompDist     = nonCoeffDist;
            currCompCost     = nonCoeffCost;

            tu.cbf[compID] = 0;
          }

          // evaluate
#if TU_256
          if( isFirstMode || ( currCompCost < minCost[compID] ) || ( transformMode == 1 && currCompCost == minCost[compID] ) )
#else
          if( ( currCompCost < minCost[compID] ) || ( transformMode == 1 && currCompCost == minCost[compID] ) )
#endif
          {
            // copy component
            if (isFirstMode && ((nonCoeffCost < currCompCost) || (currAbsSum == 0))) // check for forced null
            {
              tu.getCoeffs( compID ).fill( 0 );
              csFull->getResiBuf( compArea ).fill( 0 );
              tu.cbf[compID]   = 0;

              currAbsSum       = 0;
              currCompFracBits = nonCoeffFracBits;
              currCompDist     = nonCoeffDist;
              currCompCost     = nonCoeffCost;
            }

            uiSingleDistComp[compID] = currCompDist;
            uiSingleFracBits[compID] = currCompFracBits;
            minCost[compID]          = currCompCost;

              bestTU.copyComponentFrom( tu, compID );
              saveCS.getResiBuf( compArea ).copyFrom( csFull->getResiBuf( compArea ) );
          }
          if( tu.noResidual )
          {
            CHECK( currCompFracBits > 0 || currAbsSum, "currCompFracBits > 0 when tu noResidual" );
          }
      }

        // copy component
        tu.copyComponentFrom( bestTU, compID );
        csFull->getResiBuf( compArea ).copyFrom( saveCS.getResiBuf( compArea ) );
      if (colorTransFlag && (m_pcEncCfg->getCostMode() != COST_LOSSLESS_CODING || !slice.isLossless()))
      {
        m_pcTrQuant->lambdaAdjustColorTrans(false);
        m_pcRdCost->lambdaAdjustColorTrans(false, compID);
      }
#if SIGN_PREDICTION
      if(cs.sps->getNumPredSigns() > 0)
      {
#if JVET_Z0118_GDR
#if JVET_Y0065_GPM_INTRA
        bool lmcsEnable = cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( compID ) && !tu.cu->firstPU->ciipFlag && !tu.cu->firstPU->gpmIntraFlag && !CU::isIBC( *tu.cu );
#else
        bool lmcsEnable = cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( compID ) && !tu.cu->firstPU->ciipFlag && !CU::isIBC( *tu.cu );
#endif
        cs.reconstructPicture(tu.blocks[compID], m_pcReshape->getFwdLUT(), csFull, lmcsEnable);              
#else
        PelBuf picRecoBuff = tu.cs->picture->getRecoBuf( tu.blocks[compID] );

#if JVET_Y0065_GPM_INTRA
        if( cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( compID ) && !tu.cu->firstPU->ciipFlag && !tu.cu->firstPU->gpmIntraFlag && !CU::isIBC( *tu.cu ) )
#else
        if( cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( compID ) && !tu.cu->firstPU->ciipFlag && !CU::isIBC( *tu.cu ) )
#endif
        {
          picRecoBuff.rspSignal( cs.getPredBuf( tu.blocks[compID] ), m_pcReshape->getFwdLUT() );
          picRecoBuff.reconstruct( picRecoBuff, csFull->getResiBuf( tu.blocks[compID] ), tu.cu->cs->slice->clpRng( compID ) );
        }
        else
        {
          picRecoBuff.reconstruct( cs.getPredBuf( tu.blocks[compID] ), csFull->getResiBuf( tu.blocks[compID] ), tu.cu->cs->slice->clpRng( compID ) );
        }
      
#endif
      }
#endif
    } // component loop

    if (colorTransFlag)
    {
      PelUnitBuf     orgResidual = orgResi->subBuf(relativeUnitArea);
      PelUnitBuf     invColorTransResidual = m_colorTransResiBuf[2].getBuf(relativeUnitArea);
      csFull->getResiBuf(currArea).colorSpaceConvert(invColorTransResidual, false, slice.clpRng(COMPONENT_Y));
#if JVET_S0234_ACT_CRS_FIX
      if (slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[COMPONENT_Cb].width*tu.blocks[COMPONENT_Cb].height > 4)
      {
        invColorTransResidual.bufs[1].scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cb));
        invColorTransResidual.bufs[2].scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cr));
      }
#endif

      for (uint32_t c = 0; c < numTBlocks; c++)
      {
        const ComponentID compID = (ComponentID)c;
        uiSingleDistComp[c] = m_pcRdCost->getDistPart(orgResidual.bufs[c], invColorTransResidual.bufs[c], sps.getBitDepth(toChannelType(compID)), compID, DF_SSE);
        minCost[c] = m_pcRdCost->calcRdCost(uiSingleFracBits[c], uiSingleDistComp[c]);
      }
    }

    if ( chroma && isChromaEnabled(tu.chromaFormat) && tu.blocks[COMPONENT_Cb].valid() )
    {
      const CompArea& cbArea = tu.blocks[COMPONENT_Cb];
      const CompArea& crArea = tu.blocks[COMPONENT_Cr];
      bool checkJointCbCr = (sps.getJointCbCrEnabledFlag()) && (!tu.noResidual) && (TU::getCbf(tu, COMPONENT_Cb) || TU::getCbf(tu, COMPONENT_Cr));
      bool checkDCTOnly = (TU::getCbf(tu, COMPONENT_Cb) && tu.mtsIdx[COMPONENT_Cb] == MTS_DCT2_DCT2 && !TU::getCbf(tu, COMPONENT_Cr)) ||
                          (TU::getCbf(tu, COMPONENT_Cr) && tu.mtsIdx[COMPONENT_Cr] == MTS_DCT2_DCT2 && !TU::getCbf(tu, COMPONENT_Cb)) ||
                          (TU::getCbf(tu, COMPONENT_Cb) && tu.mtsIdx[COMPONENT_Cb] == MTS_DCT2_DCT2 && TU::getCbf(tu, COMPONENT_Cr) && tu.mtsIdx[COMPONENT_Cr] == MTS_DCT2_DCT2);

      bool checkTSOnly = (TU::getCbf(tu, COMPONENT_Cb) && tu.mtsIdx[COMPONENT_Cb] == MTS_SKIP && !TU::getCbf(tu, COMPONENT_Cr)) ||
                         (TU::getCbf(tu, COMPONENT_Cr) && tu.mtsIdx[COMPONENT_Cr] == MTS_SKIP && !TU::getCbf(tu, COMPONENT_Cb)) ||
                         (TU::getCbf(tu, COMPONENT_Cb) && tu.mtsIdx[COMPONENT_Cb] == MTS_SKIP && TU::getCbf(tu, COMPONENT_Cr) && tu.mtsIdx[COMPONENT_Cr] == MTS_SKIP);
      const int channelBitDepth = sps.getBitDepth(toChannelType(COMPONENT_Cb));
      bool      reshape         = slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag()
                               && tu.blocks[COMPONENT_Cb].width * tu.blocks[COMPONENT_Cb].height > 4;
      double minCostCbCr = minCost[COMPONENT_Cb] + minCost[COMPONENT_Cr];
      if (colorTransFlag)
      {
        minCostCbCr += minCost[COMPONENT_Y];  // ACT should consider three-component cost
      }

      CompStorage      orgResiCb[4], orgResiCr[4];   // 0:std, 1-3:jointCbCr
      std::vector<int> jointCbfMasksToTest;
      if ( checkJointCbCr )
      {
        orgResiCb[0].create(cbArea);
        orgResiCr[0].create(crArea);
#if JVET_S0234_ACT_CRS_FIX
        if (colorTransFlag)
        {
          orgResiCb[0].copyFrom(colorTransResidual.bufs[1]);
          orgResiCr[0].copyFrom(colorTransResidual.bufs[2]);
        }
        else
        {
#endif
        orgResiCb[0].copyFrom(cs.getOrgResiBuf(cbArea));
        orgResiCr[0].copyFrom(cs.getOrgResiBuf(crArea));
#if JVET_S0234_ACT_CRS_FIX
        }
        if (!colorTransFlag && reshape)
#else
        if (reshape)
#endif
        {
          orgResiCb[0].scaleSignal(tu.getChromaAdj(), 1, tu.cu->cs->slice->clpRng(COMPONENT_Cb));
          orgResiCr[0].scaleSignal(tu.getChromaAdj(), 1, tu.cu->cs->slice->clpRng(COMPONENT_Cr));
        }
        jointCbfMasksToTest = m_pcTrQuant->selectICTCandidates(tu, orgResiCb, orgResiCr);
      }

      for (int cbfMask: jointCbfMasksToTest)
      {
        ComponentID codeCompId = (cbfMask >> 1 ? COMPONENT_Cb : COMPONENT_Cr);
        ComponentID otherCompId = (codeCompId == COMPONENT_Cr ? COMPONENT_Cb : COMPONENT_Cr);
        bool        tsAllowed = TU::isTSAllowed(tu, codeCompId) && (m_pcEncCfg->getUseChromaTS());
        uint8_t     numTransformCands = 1 + (tsAllowed ? 1 : 0); // DCT + TS = 2 tests
        bool        cbfDCT2 = true;

        std::vector<TrMode> trModes;
        if (checkDCTOnly || checkTSOnly)
        {
          numTransformCands = 1;
        }

        if (!checkTSOnly)
        {
          trModes.push_back(TrMode(0, true)); // DCT2
        }
        if (tsAllowed && !checkDCTOnly)
        {
          trModes.push_back(TrMode(1, true));//TS
        }
        for (int modeId = 0; modeId < numTransformCands; modeId++)
        {
          if (modeId && !cbfDCT2)
          {
            continue;
          }
          if (!trModes[modeId].second)
          {
            continue;
          }
        TCoeff     currAbsSum       = 0;
        uint64_t   currCompFracBits = 0;
        Distortion currCompDistCb   = 0;
        Distortion currCompDistCr   = 0;
        double     currCompCost     = 0;

        tu.jointCbCr = (uint8_t) cbfMask;
          // encoder bugfix: initialize mtsIdx for chroma under JointCbCrMode.
        tu.mtsIdx[codeCompId]  = trModes[modeId].first;
        tu.mtsIdx[otherCompId] = MTS_DCT2_DCT2;
        int         codedCbfMask = 0;
        if (colorTransFlag && (m_pcEncCfg->getCostMode() != COST_LOSSLESS_CODING || !slice.isLossless()))
        {
          m_pcTrQuant->lambdaAdjustColorTrans(true);
          m_pcTrQuant->selectLambda(codeCompId);
        }
        else
        {
          m_pcTrQuant->selectLambda(codeCompId);
        }
        // Lambda is loosened for the joint mode with respect to single modes as the same residual is used for both chroma blocks
        const int    absIct = abs( TU::getICTMode(tu) );
        const double lfact  = ( absIct == 1 || absIct == 3 ? 0.8 : 0.5 );
        m_pcTrQuant->setLambda( lfact * m_pcTrQuant->getLambda() );
        if ( checkJointCbCr && (tu.cu->cs->slice->getSliceQp() > 18))
        {
          m_pcTrQuant->setLambda( 1.05 * m_pcTrQuant->getLambda() );
        }

        m_CABACEstimator->getCtx() = ctxStart;
        m_CABACEstimator->resetBits();

        PelBuf cbResi = csFull->getResiBuf(cbArea);
        PelBuf crResi = csFull->getResiBuf(crArea);
        cbResi.copyFrom(orgResiCb[cbfMask]);
        crResi.copyFrom(orgResiCr[cbfMask]);

        if ( reshape )
        {
          double cRescale = (double)(1 << CSCALE_FP_PREC) / (double)(tu.getChromaAdj());
          m_pcTrQuant->setLambda(m_pcTrQuant->getLambda() / (cRescale*cRescale));
        }

        Distortion currCompDistY = MAX_UINT64;
        QpParam qpCbCr(tu, codeCompId);

        tu.getCoeffs(otherCompId).fill(0);   // do we need that?
        TU::setCbfAtDepth(tu, otherCompId, tu.depth, false);

        PelBuf &codeResi   = (codeCompId == COMPONENT_Cr ? crResi : cbResi);
        TCoeff  compAbsSum = 0;
        if (numTransformCands > 1)
        {
          if (modeId == 0)
          {
            m_pcTrQuant->transformNxN(tu, codeCompId, qpCbCr, &trModes, m_pcEncCfg->getMTSInterMaxCand());
            tu.mtsIdx[codeCompId] = trModes[modeId].first;
            tu.mtsIdx[otherCompId] = MTS_DCT2_DCT2;
          }
          m_pcTrQuant->transformNxN(tu, codeCompId, qpCbCr, compAbsSum, m_CABACEstimator->getCtx(), true);
        }
        else
        m_pcTrQuant->transformNxN(tu, codeCompId, qpCbCr, compAbsSum, m_CABACEstimator->getCtx());
        if (compAbsSum > 0)
        {
          m_pcTrQuant->invTransformNxN(tu, codeCompId, codeResi, qpCbCr);
          codedCbfMask += (codeCompId == COMPONENT_Cb ? 2 : 1);
        }
        else
        {
          codeResi.fill(0);
        }

        if (tu.jointCbCr == 3 && codedCbfMask == 2)
        {
          codedCbfMask = 3;
          TU::setCbfAtDepth(tu, COMPONENT_Cr, tu.depth, true);
        }
        if (codedCbfMask && tu.jointCbCr != codedCbfMask)
        {
          codedCbfMask = 0;
        }
        currAbsSum = codedCbfMask;

        if (!tu.mtsIdx[codeCompId])
        {
          cbfDCT2 = (currAbsSum > 0);
        }
        if (currAbsSum > 0)
        {
          m_CABACEstimator->cbf_comp(cs, codedCbfMask >> 1, cbArea, currDepth, false);
          m_CABACEstimator->cbf_comp(cs, codedCbfMask & 1, crArea, currDepth, codedCbfMask >> 1);
          m_CABACEstimator->joint_cb_cr(tu, codedCbfMask);
          if (codedCbfMask >> 1)
            m_CABACEstimator->residual_coding(tu, COMPONENT_Cb);
          if (codedCbfMask & 1)
            m_CABACEstimator->residual_coding(tu, COMPONENT_Cr);
          currCompFracBits = m_CABACEstimator->getEstFracBits();

          m_pcTrQuant->invTransformICT(tu, cbResi, crResi);
#if JVET_S0234_ACT_CRS_FIX
          if (!colorTransFlag && reshape)
#else
          if (reshape)
#endif
          {
            cbResi.scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cb));
            crResi.scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cr));
          }

          if (colorTransFlag)
          {
            PelUnitBuf     orgResidual = orgResi->subBuf(relativeUnitArea);
            PelUnitBuf     invColorTransResidual = m_colorTransResiBuf[2].getBuf(relativeUnitArea);
            csFull->getResiBuf(currArea).colorSpaceConvert(invColorTransResidual, false, slice.clpRng(COMPONENT_Y));
#if JVET_S0234_ACT_CRS_FIX
            if (reshape)
            {
              invColorTransResidual.bufs[1].scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cb));
              invColorTransResidual.bufs[2].scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cr));
            }
#endif

            currCompDistY = m_pcRdCost->getDistPart(orgResidual.bufs[COMPONENT_Y], invColorTransResidual.bufs[COMPONENT_Y], sps.getBitDepth(toChannelType(COMPONENT_Y)), COMPONENT_Y, DF_SSE);
            currCompDistCb = m_pcRdCost->getDistPart(orgResidual.bufs[COMPONENT_Cb], invColorTransResidual.bufs[COMPONENT_Cb], sps.getBitDepth(toChannelType(COMPONENT_Cb)), COMPONENT_Cb, DF_SSE);
            currCompDistCr = m_pcRdCost->getDistPart(orgResidual.bufs[COMPONENT_Cr], invColorTransResidual.bufs[COMPONENT_Cr], sps.getBitDepth(toChannelType(COMPONENT_Cr)), COMPONENT_Cr, DF_SSE);
            currCompCost = m_pcRdCost->calcRdCost(uiSingleFracBits[COMPONENT_Y] + currCompFracBits, currCompDistY + currCompDistCr + currCompDistCb, false);
          }
          else
          {
          currCompDistCb = m_pcRdCost->getDistPart(csFull->getOrgResiBuf(cbArea), cbResi, channelBitDepth, COMPONENT_Cb, DF_SSE);
          currCompDistCr = m_pcRdCost->getDistPart(csFull->getOrgResiBuf(crArea), crResi, channelBitDepth, COMPONENT_Cr, DF_SSE);
#if WCG_EXT
          currCompCost   = m_pcRdCost->calcRdCost(currCompFracBits, currCompDistCr + currCompDistCb, false);
#else
          currCompCost   = m_pcRdCost->calcRdCost(currCompFracBits, currCompDistCr + currCompDistCb);
#endif
          }
        }
        else
          currCompCost = MAX_DOUBLE;

        // evaluate
        if( currCompCost < minCostCbCr )
        {
          uiSingleDistComp[COMPONENT_Cb] = currCompDistCb;
          uiSingleDistComp[COMPONENT_Cr] = currCompDistCr;
          if (colorTransFlag)
          {
            uiSingleDistComp[COMPONENT_Y] = currCompDistY;
          }
          minCostCbCr                    = currCompCost;
          {
            bestTU.copyComponentFrom(tu, COMPONENT_Cb);
            bestTU.copyComponentFrom(tu, COMPONENT_Cr);
            saveCS.getResiBuf(cbArea).copyFrom(csFull->getResiBuf(cbArea));
            saveCS.getResiBuf(crArea).copyFrom(csFull->getResiBuf(crArea));
          }
        }

        if (colorTransFlag && (m_pcEncCfg->getCostMode() != COST_LOSSLESS_CODING || !slice.isLossless()))
        {
          m_pcTrQuant->lambdaAdjustColorTrans(false);
        }
        }
      }
      // copy component
      tu.copyComponentFrom(bestTU, COMPONENT_Cb);
      tu.copyComponentFrom(bestTU, COMPONENT_Cr);
      csFull->getResiBuf(cbArea).copyFrom(saveCS.getResiBuf(cbArea));
      csFull->getResiBuf(crArea).copyFrom(saveCS.getResiBuf(crArea));

#if SIGN_PREDICTION
      if( tu.jointCbCr )
      {
        for( auto i = (int)COMPONENT_Cb; i <= (int)COMPONENT_Cr; ++i)
        {
          ComponentID comp = (ComponentID) i;
#if JVET_Z0118_GDR
#if JVET_Y0065_GPM_INTRA
          bool lmcsEnable = cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( comp ) && !tu.cu->firstPU->ciipFlag && !tu.cu->firstPU->gpmIntraFlag && !CU::isIBC( *tu.cu );
#else         
          bool lmcsEnable = cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( comp ) && !tu.cu->firstPU->ciipFlag && !CU::isIBC( *tu.cu );
#endif
          cs.reconstructPicture(tu.blocks[comp], m_pcReshape->getFwdLUT(), csFull, lmcsEnable);          
#else
          PelBuf picRecoBuff = tu.cs->picture->getRecoBuf( tu.blocks[comp] );

#if JVET_Y0065_GPM_INTRA
          if( cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( comp ) && !tu.cu->firstPU->ciipFlag && !tu.cu->firstPU->gpmIntraFlag && !CU::isIBC( *tu.cu ) )
#else
          if( cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( comp ) && !tu.cu->firstPU->ciipFlag && !CU::isIBC( *tu.cu ) )
#endif
          {
            picRecoBuff.rspSignal( cs.getPredBuf( tu.blocks[comp] ), m_pcReshape->getFwdLUT() );
            picRecoBuff.reconstruct( picRecoBuff, csFull->getResiBuf( tu.blocks[comp] ), tu.cu->cs->slice->clpRng( comp ) );
          }
          else
          {
            picRecoBuff.reconstruct( cs.getPredBuf( tu.blocks[comp] ), csFull->getResiBuf( tu.blocks[comp] ), tu.cu->cs->slice->clpRng( comp ) );
          }
#endif
        }

        if ( sps.getNumPredSigns() > 0)
        {
          bool bJccrWithCr = tu.jointCbCr && !(tu.jointCbCr >> 1);
          ComponentID jccrCompId = bJccrWithCr ? COMPONENT_Cr : COMPONENT_Cb;
          bool reshapeChroma = slice.getPicHeader()->getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[jccrCompId].width*tu.blocks[jccrCompId].height > 4;
          m_pcTrQuant->predCoeffSigns(tu, COMPONENT_Cb, reshapeChroma);
        }
      }
#endif
    }

    m_CABACEstimator->getCtx() = ctxStart;
    m_CABACEstimator->resetBits();
    if( !tu.noResidual )
    {
      static const ComponentID cbf_getComp[MAX_NUM_COMPONENT] = { COMPONENT_Cb, COMPONENT_Cr, COMPONENT_Y };
      for( unsigned c = isChromaEnabled(tu.chromaFormat)?0 : 2; c < MAX_NUM_COMPONENT; c++)
    {
      const ComponentID compID = cbf_getComp[c];
      if (compID == COMPONENT_Y && !luma)
        continue;
      if (compID != COMPONENT_Y && !chroma)
        continue;
      if( tu.blocks[compID].valid() )
      {
        const bool prevCbf = ( compID == COMPONENT_Cr ? TU::getCbfAtDepth( tu, COMPONENT_Cb, currDepth ) : false );
        m_CABACEstimator->cbf_comp( *csFull, TU::getCbfAtDepth( tu, compID, currDepth ), tu.blocks[compID], currDepth, prevCbf );
      }
    }
    }

    for (uint32_t ch = 0; ch < numValidComp; ch++)
    {
      const ComponentID compID = ComponentID(ch);
      if (compID == COMPONENT_Y && !luma)
        continue;
      if (compID != COMPONENT_Y && !chroma)
        continue;
      if (tu.blocks[compID].valid())
      {
        if( compID == COMPONENT_Cr )
        {
          const int cbfMask = ( TU::getCbf( tu, COMPONENT_Cb ) ? 2 : 0 ) + ( TU::getCbf( tu, COMPONENT_Cr ) ? 1 : 0 );
          m_CABACEstimator->joint_cb_cr(tu, cbfMask);
        }
        if( TU::getCbf( tu, compID ) )
        {
          m_CABACEstimator->residual_coding( tu, compID );
        }
        uiSingleDist += uiSingleDistComp[compID];
      }
    }
    if( tu.noResidual )
    {
      CHECK( m_CABACEstimator->getEstFracBits() > 0, "no residual TU's bits shall be 0" );
    }
#if JVET_S0234_ACT_CRS_FIX
    if (colorTransFlag)
    {
      PelUnitBuf resiBuf = csFull->getResiBuf(currArea);
      resiBuf.colorSpaceConvert(resiBuf, false, slice.clpRng(COMPONENT_Y));
      if (slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && tu.blocks[COMPONENT_Cb].width*tu.blocks[COMPONENT_Cb].height > 4)
      {
        resiBuf.bufs[1].scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cb));
        resiBuf.bufs[2].scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(COMPONENT_Cr));
      }
    }
#endif

    csFull->fracBits += m_CABACEstimator->getEstFracBits();
    csFull->dist     += uiSingleDist;
#if WCG_EXT
    if( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() )
    {
      csFull->cost    = m_pcRdCost->calcRdCost(csFull->fracBits, csFull->dist, false);
    }
    else
#endif
    csFull->cost      = m_pcRdCost->calcRdCost(csFull->fracBits, csFull->dist);
  } // check full

  // code sub-blocks
  if( bCheckSplit )
  {
    if( bCheckFull )
    {
      m_CABACEstimator->getCtx() = ctxStart;
    }

    if( partitioner.canSplit( TU_MAX_TR_SPLIT, cs ) )
    {
      partitioner.splitCurrArea( TU_MAX_TR_SPLIT, cs );
    }
    else if( cu.sbtInfo && partitioner.canSplit( PartSplit( cu.getSbtTuSplit() ), cs ) )
    {
      partitioner.splitCurrArea( PartSplit( cu.getSbtTuSplit() ), cs );
    }
    else
      THROW( "Implicit TU split not available!" );

    do
    {
      xEstimateInterResidualQT(*csSplit, partitioner, bCheckFull ? nullptr : puiZeroDist
        , luma, chroma
        , orgResi
      );

      csSplit->cost = m_pcRdCost->calcRdCost( csSplit->fracBits, csSplit->dist );
    } while( partitioner.nextPart( *csSplit ) );

    partitioner.exitCurrSplit();

    unsigned        anyCbfSet   =   0;
    unsigned        compCbf[3]  = { 0, 0, 0 };

    if( !bCheckFull )
    {
      for( auto &currTU : csSplit->traverseTUs( currArea, partitioner.chType ) )
      {
        for( unsigned ch = 0; ch < numTBlocks; ch++ )
        {
          compCbf[ ch ] |= ( TU::getCbfAtDepth( currTU, ComponentID(ch), currDepth + 1 ) ? 1 : 0 );
        }
      }

      {

        for( auto &currTU : csSplit->traverseTUs( currArea, partitioner.chType ) )
        {
          TU::setCbfAtDepth   ( currTU, COMPONENT_Y,  currDepth, compCbf[ COMPONENT_Y  ] );
          if( currArea.chromaFormat != CHROMA_400 )
          {
            TU::setCbfAtDepth ( currTU, COMPONENT_Cb, currDepth, compCbf[ COMPONENT_Cb ] );
            TU::setCbfAtDepth ( currTU, COMPONENT_Cr, currDepth, compCbf[ COMPONENT_Cr ] );
          }
        }

        anyCbfSet    = compCbf[ COMPONENT_Y  ];
        if( currArea.chromaFormat != CHROMA_400 )
        {
          anyCbfSet |= compCbf[ COMPONENT_Cb ];
          anyCbfSet |= compCbf[ COMPONENT_Cr ];
        }
      }

      m_CABACEstimator->getCtx() = ctxStart;
      m_CABACEstimator->resetBits();

      // when compID isn't a channel, code Cbfs:
      xEncodeInterResidualQT( *csSplit, partitioner, MAX_NUM_TBLOCKS );
      for (uint32_t ch = 0; ch < numValidComp; ch++)
      {
        const ComponentID compID = ComponentID(ch);
        if (compID == COMPONENT_Y && !luma)
          continue;
        if (compID != COMPONENT_Y && !chroma)
          continue;
        xEncodeInterResidualQT( *csSplit, partitioner, ComponentID( ch ) );
      }

      csSplit->fracBits = m_CABACEstimator->getEstFracBits();
      csSplit->cost     = m_pcRdCost->calcRdCost(csSplit->fracBits, csSplit->dist);

      if( bCheckFull && anyCbfSet && csSplit->cost < csFull->cost )
      {
        cs.useSubStructure( *csSplit, partitioner.chType, currArea, false, false, false, true, true );
        cs.cost = csSplit->cost;
      }
    }


    if( csSplit && csFull )
    {
      csSplit->releaseIntermediateData();
      csFull ->releaseIntermediateData();
    }
  }
}

void InterSearch::encodeResAndCalcRdInterCU(CodingStructure &cs, Partitioner &partitioner, const bool &skipResidual
  , const bool luma, const bool chroma
)
{
  m_pcRdCost->setChromaFormat(cs.sps->getChromaFormatIdc());

  CodingUnit &cu = *cs.getCU( partitioner.chType );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( cu.predMode == MODE_INTER )
    CHECK( cu.isSepTree(), "CU with Inter mode must be in single tree" );
#endif
  const ChromaFormat format     = cs.area.chromaFormat;;
  const int  numValidComponents = getNumberValidComponents(format);
  const SPS &sps                = *cs.sps;

  bool colorTransAllowed = cs.slice->getSPS()->getUseColorTrans() && luma && chroma;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (cs.slice->getSPS()->getUseColorTrans())
  {
    CHECK(cu.treeType != TREE_D || partitioner.treeType != TREE_D, "localtree should not be applied when adaptive color transform is enabled");
    CHECK(cu.modeType != MODE_TYPE_ALL || partitioner.modeType != MODE_TYPE_ALL, "localtree should not be applied when adaptive color transform is enabled");
  }
#endif
  if( skipResidual ) //  No residual coding : SKIP mode
  {
    cu.skip    = true;
    cu.rootCbf = false;
    cu.colorTransform = false;
    CHECK( cu.sbtInfo != 0, "sbtInfo shall be 0 if CU has no residual" );
    cs.getResiBuf().fill(0);

#if JVET_Y0065_GPM_INTRA
    if( m_pcEncCfg->getLmcs() && ( cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() ) && !cu.firstPU->ciipFlag && !cu.firstPU->gpmIntraFlag && !CU::isIBC( cu ) )
#else
    if( m_pcEncCfg->getLmcs() && ( cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() ) && !cu.firstPU->ciipFlag && !CU::isIBC( cu ) )
#endif
    {
      cs.getRecoBuf().Y().rspSignal( cs.getPredBuf().Y(), m_pcReshape->getFwdLUT() );
      cs.getRecoBuf().Cb().copyFrom( cs.getPredBuf().Cb() );
      cs.getRecoBuf().Cr().copyFrom( cs.getPredBuf().Cr() );
    }
    else
    {
      cs.getRecoBuf().copyFrom( cs.getPredBuf() );
    }

    // add empty TU(s)
    cs.addEmptyTUs( partitioner );
    Distortion distortion = 0;

    for (int comp = 0; comp < numValidComponents; comp++)
    {
      const ComponentID compID = ComponentID(comp);
      if (compID == COMPONENT_Y && !luma)
        continue;
      if (compID != COMPONENT_Y && !chroma)
        continue;
      CPelBuf reco = cs.getRecoBuf (compID);
      CPelBuf org  = cs.getOrgBuf  (compID);
#if WCG_EXT
      if (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() || (
        m_pcEncCfg->getLmcs() && (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())))
      {
        const CPelBuf orgLuma = cs.getOrgBuf( cs.area.blocks[COMPONENT_Y] );
        if (compID == COMPONENT_Y && !(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()))
        {
          const CompArea &areaY = cu.Y();
          CompArea      tmpArea1(COMPONENT_Y, areaY.chromaFormat, Position(0, 0), areaY.size());
          PelBuf tmpRecLuma = m_tmpStorageLCU.getBuf(tmpArea1);
          tmpRecLuma.rspSignal( reco, m_pcReshape->getInvLUT() );
          distortion += m_pcRdCost->getDistPart(org, tmpRecLuma, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
        }
        else
        distortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE_WTD, &orgLuma );
      }
      else
#endif
      distortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
    }

    m_CABACEstimator->resetBits();

    PredictionUnit &pu = *cs.getPU( partitioner.chType );

    m_CABACEstimator->cu_skip_flag  ( cu );
    m_CABACEstimator->merge_data(pu);
#if INTER_LIC
    m_CABACEstimator->cu_lic_flag(cu);
#endif

    cs.dist     = distortion;
    cs.fracBits = m_CABACEstimator->getEstFracBits();
    cs.cost     = m_pcRdCost->calcRdCost(cs.fracBits, cs.dist);

    return;
  }

  //  Residual coding.
  if( luma )
  {
    if( cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() )
    {
#if JVET_Y0065_GPM_INTRA
      if( !cu.firstPU->ciipFlag && !cu.firstPU->gpmIntraFlag && !CU::isIBC( cu ) )
#else
      if( !cu.firstPU->ciipFlag && !CU::isIBC( cu ) )
#endif
      {
        cs.getResiBuf( COMPONENT_Y ).rspSignalAllAndSubtract( cs.getOrgBuf( COMPONENT_Y ), cs.getPredBuf( COMPONENT_Y ), m_pcReshape->getFwdLUT() );
      }
      else
      {
        cs.getResiBuf( COMPONENT_Y ).rspSignalAndSubtract( cs.getOrgBuf( COMPONENT_Y ), cs.getPredBuf( COMPONENT_Y ), m_pcReshape->getFwdLUT() );
      }
    }
    else
    {
      cs.getResiBuf( COMPONENT_Y ).subtract( cs.getOrgBuf( COMPONENT_Y ), cs.getPredBuf( COMPONENT_Y ) );
    }
  }

  if( chroma && isChromaEnabled( cs.pcv->chrFormat ) )
  {
    cs.getResiBuf( COMPONENT_Cb ).subtract( cs.getOrgBuf( COMPONENT_Cb ), cs.getPredBuf( COMPONENT_Cb ) );
    cs.getResiBuf( COMPONENT_Cr ).subtract( cs.getOrgBuf( COMPONENT_Cr ), cs.getPredBuf( COMPONENT_Cr ) );
  }

  const UnitArea curUnitArea = partitioner.currArea();
  CodingStructure &saveCS = *m_pSaveCS[1];
  saveCS.pcv = cs.pcv;
  saveCS.picture = cs.picture;
#if JVET_Z0118_GDR
  saveCS.m_pt = cs.m_pt;
#endif
  saveCS.area.repositionTo(curUnitArea);
  saveCS.clearCUs();
  saveCS.clearPUs();
  saveCS.clearTUs();
  for (const auto &ppcu : cs.cus)
  {
    CodingUnit &pcu = saveCS.addCU(*ppcu, ppcu->chType);
    pcu = *ppcu;
  }
  for (const auto &ppu : cs.pus)
  {
    PredictionUnit &pu = saveCS.addPU(*ppu, ppu->chType);
    pu = *ppu;
  }

#if JVET_S0234_ACT_CRS_FIX
  PelUnitBuf orgResidual;
#else
  PelUnitBuf orgResidual, colorTransResidual;
#endif
  const UnitArea localUnitArea(cs.area.chromaFormat, Area(0, 0, cu.Y().width, cu.Y().height));
  orgResidual = m_colorTransResiBuf[0].getBuf(localUnitArea);
#if !JVET_S0234_ACT_CRS_FIX
  colorTransResidual = m_colorTransResiBuf[1].getBuf(localUnitArea);
#endif
  orgResidual.copyFrom(cs.getResiBuf());
#if !JVET_S0234_ACT_CRS_FIX
  if (colorTransAllowed)
  {
    cs.getResiBuf().colorSpaceConvert(colorTransResidual, true, cu.cs->slice->clpRng(COMPONENT_Y));
  }
#endif

  const TempCtx ctxStart(m_CtxCache, m_CABACEstimator->getCtx());
  int           numAllowedColorSpace = (colorTransAllowed ? 2 : 1);
  Distortion    zeroDistortion = 0;

  double  bestCost = MAX_DOUBLE;
  bool    bestColorTrans = false;
  bool    bestRootCbf = false;
  uint8_t bestsbtInfo = 0;
  uint8_t orgSbtInfo = cu.sbtInfo;
  int     bestIter = 0;

  auto blkCache = dynamic_cast<CacheBlkInfoCtrl*>(m_modeCtrl);
  bool rootCbfFirstColorSpace = true;

  for (int iter = 0; iter < numAllowedColorSpace; iter++)
  {
    if (colorTransAllowed && !m_pcEncCfg->getRGBFormatFlag() && iter)
    {
      continue;
    }
    char colorSpaceOption = blkCache->getSelectColorSpaceOption(cu);
    if (colorTransAllowed)
    {
      if (colorSpaceOption)
      {
        CHECK(colorSpaceOption > 2 || colorSpaceOption < 0, "invalid color space selection option");
        if (colorSpaceOption == 1 && iter)
        {
          continue;
        }
        if (colorSpaceOption == 2 && !iter)
        {
          continue;
        }
      }
    }
    if (!colorSpaceOption)
    {
      if (iter && !rootCbfFirstColorSpace)
      {
        continue;
      }
      if (colorTransAllowed && cs.bestParent && cs.bestParent->tmpColorSpaceCost != MAX_DOUBLE)
      {
        if (cs.bestParent->firstColorSpaceSelected && iter)
        {
          continue;
        }
        if (m_pcEncCfg->getRGBFormatFlag())
        {
          if (!cs.bestParent->firstColorSpaceSelected && !iter)
          {
            continue;
          }
        }
      }
    }
    bool colorTransFlag = (colorTransAllowed && m_pcEncCfg->getRGBFormatFlag()) ? (1 - iter) : iter;
    cu.colorTransform = colorTransFlag;
    cu.sbtInfo = orgSbtInfo;

    m_CABACEstimator->resetBits();
    m_CABACEstimator->getCtx() = ctxStart;
    cs.clearTUs();
    cs.fracBits = 0;
    cs.dist = 0;
    cs.cost = 0;

  if (colorTransFlag)
  {
#if JVET_S0234_ACT_CRS_FIX
    cs.getOrgResiBuf().bufs[0].copyFrom(orgResidual.bufs[0]);
    cs.getOrgResiBuf().bufs[1].copyFrom(orgResidual.bufs[1]);
    cs.getOrgResiBuf().bufs[2].copyFrom(orgResidual.bufs[2]);
#else
    cs.getOrgResiBuf().bufs[0].copyFrom(colorTransResidual.bufs[0]);
    cs.getOrgResiBuf().bufs[1].copyFrom(colorTransResidual.bufs[1]);
    cs.getOrgResiBuf().bufs[2].copyFrom(colorTransResidual.bufs[2]);
#endif

    memset(m_pTempPel, 0, sizeof(Pel) * localUnitArea.blocks[0].area());
    zeroDistortion = 0;
    for (int compIdx = 0; compIdx < 3; compIdx++)
    {
      ComponentID componentID = (ComponentID)compIdx;
      const CPelBuf zeroBuf(m_pTempPel, localUnitArea.blocks[compIdx]);
      zeroDistortion += m_pcRdCost->getDistPart(zeroBuf, orgResidual.bufs[compIdx], sps.getBitDepth(toChannelType(componentID)), componentID, DF_SSE);
    }
    xEstimateInterResidualQT(cs, partitioner, NULL, luma, chroma, &orgResidual);
  }
  else
  {
    zeroDistortion = 0;
  if (luma)
  {
    cs.getOrgResiBuf().bufs[0].copyFrom(orgResidual.bufs[0]);
  }
  if (chroma && isChromaEnabled(cs.pcv->chrFormat))
  {
    cs.getOrgResiBuf().bufs[1].copyFrom(orgResidual.bufs[1]);
    cs.getOrgResiBuf().bufs[2].copyFrom(orgResidual.bufs[2]);
  }
  xEstimateInterResidualQT(cs, partitioner, &zeroDistortion, luma, chroma);
  }
  TransformUnit &firstTU = *cs.getTU( partitioner.chType );

  cu.rootCbf = false;
  m_CABACEstimator->resetBits();
  m_CABACEstimator->rqt_root_cbf( cu );
  const uint64_t  zeroFracBits = m_CABACEstimator->getEstFracBits();
  double zeroCost;
  {
#if WCG_EXT
    if( m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() )
    {
      zeroCost = m_pcRdCost->calcRdCost( zeroFracBits, zeroDistortion, false );
    }
    else
#endif
    zeroCost = m_pcRdCost->calcRdCost( zeroFracBits, zeroDistortion );
  }

  const int  numValidTBlocks   = ::getNumberValidTBlocks( *cs.pcv );
  for (uint32_t i = 0; i < numValidTBlocks; i++)
  {
    cu.rootCbf |= TU::getCbfAtDepth(firstTU, ComponentID(i), 0);
  }

  // -------------------------------------------------------
  // If a block full of 0's is efficient, then just use 0's.
  // The costs at this point do not include header bits.

  if (zeroCost < cs.cost || !cu.rootCbf)
  {
    cs.cost = zeroCost;
    cu.colorTransform = false;
    cu.sbtInfo = 0;
    cu.rootCbf = false;

    cs.clearTUs();

    // add new "empty" TU(s) spanning the whole CU
    cs.addEmptyTUs( partitioner );
  }
  if (!iter)
  {
    rootCbfFirstColorSpace = cu.rootCbf;
  }
  if (cs.cost < bestCost)
  {
    bestIter = iter;
#if !JVET_S0234_ACT_CRS_FIX
    if (cu.rootCbf && cu.colorTransform)
    {
      cs.getResiBuf(curUnitArea).colorSpaceConvert(cs.getResiBuf(curUnitArea), false, cu.cs->slice->clpRng(COMPONENT_Y));
    }
#endif

    if (iter != (numAllowedColorSpace - 1))
    {
      bestCost = cs.cost;
      bestColorTrans = cu.colorTransform;
      bestRootCbf = cu.rootCbf;
      bestsbtInfo = cu.sbtInfo;

      saveCS.clearTUs();
      for (const auto &ptu : cs.tus)
      {
        TransformUnit &tu = saveCS.addTU(*ptu, ptu->chType);
        tu = *ptu;
      }
      saveCS.getResiBuf(curUnitArea).copyFrom(cs.getResiBuf(curUnitArea));
    }
  }
  }

  if (bestIter != (numAllowedColorSpace - 1))
  {
    cu.colorTransform = bestColorTrans;
    cu.rootCbf = bestRootCbf;
    cu.sbtInfo = bestsbtInfo;

    cs.clearTUs();
    for (const auto &ptu : saveCS.tus)
    {
      TransformUnit &tu = cs.addTU(*ptu, ptu->chType);
      tu = *ptu;
    }
    cs.getResiBuf(curUnitArea).copyFrom(saveCS.getResiBuf(curUnitArea));
  }

  // all decisions now made. Fully encode the CU, including the headers:
  m_CABACEstimator->getCtx() = ctxStart;

  uint64_t finalFracBits = xGetSymbolFracBitsInter( cs, partitioner );
  // we've now encoded the CU, and so have a valid bit cost
  if (!cu.rootCbf)
  {
    if (luma)
    {
      cs.getResiBuf().bufs[0].fill(0); // Clear the residual image, if we didn't code it.
    }
    if (chroma && isChromaEnabled(cs.pcv->chrFormat))
    {
      cs.getResiBuf().bufs[1].fill(0); // Clear the residual image, if we didn't code it.
      cs.getResiBuf().bufs[2].fill(0); // Clear the residual image, if we didn't code it.
    }
  }

  if (luma)
  {
    if (cu.rootCbf && cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
    {
#if JVET_Y0065_GPM_INTRA
      if( !cu.firstPU->ciipFlag && !cu.firstPU->gpmIntraFlag && !CU::isIBC( cu ) )
#else
      if( !cu.firstPU->ciipFlag && !CU::isIBC( cu ) )
#endif
      {
        const CompArea &areaY = cu.Y();
        CompArea      tmpArea( COMPONENT_Y, areaY.chromaFormat, Position( 0, 0 ), areaY.size() );
        PelBuf tmpPred = m_tmpStorageLCU.getBuf( tmpArea );
        tmpPred.rspSignal( cs.getPredBuf( COMPONENT_Y ), m_pcReshape->getFwdLUT() );

        cs.getRecoBuf( COMPONENT_Y ).reconstruct( tmpPred, cs.getResiBuf( COMPONENT_Y ), cs.slice->clpRng( COMPONENT_Y ) );
      }
      else
      {
        cs.getRecoBuf( COMPONENT_Y ).reconstruct( cs.getPredBuf( COMPONENT_Y ), cs.getResiBuf( COMPONENT_Y ), cs.slice->clpRng( COMPONENT_Y ) );
      }
    }
    else
    {
      cs.getRecoBuf().bufs[0].reconstruct(cs.getPredBuf().bufs[0], cs.getResiBuf().bufs[0], cs.slice->clpRngs().comp[0]);
#if JVET_Y0065_GPM_INTRA
      if (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !cu.firstPU->ciipFlag && !cu.firstPU->gpmIntraFlag && !CU::isIBC(cu))
#else
      if (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !cu.firstPU->ciipFlag && !CU::isIBC(cu))
#endif
      {
        cs.getRecoBuf().bufs[0].rspSignal(m_pcReshape->getFwdLUT());
      }
    }
  }
  if (chroma && isChromaEnabled(cs.pcv->chrFormat))
  {
    cs.getRecoBuf().bufs[1].reconstruct(cs.getPredBuf().bufs[1], cs.getResiBuf().bufs[1], cs.slice->clpRngs().comp[1]);
    cs.getRecoBuf().bufs[2].reconstruct(cs.getPredBuf().bufs[2], cs.getResiBuf().bufs[2], cs.slice->clpRngs().comp[2]);
  }

  // update with clipped distortion and cost (previously unclipped reconstruction values were used)
  Distortion finalDistortion = 0;

  for (int comp = 0; comp < numValidComponents; comp++)
  {
    const ComponentID compID = ComponentID(comp);
    if (compID == COMPONENT_Y && !luma)
      continue;
    if (compID != COMPONENT_Y && !chroma)
      continue;
    CPelBuf reco = cs.getRecoBuf (compID);
    CPelBuf org  = cs.getOrgBuf  (compID);
#if JVET_V0094_BILATERAL_FILTER
    const CompArea &areaY = cu.Y();
    CompArea      tmpArea1(COMPONENT_Y, areaY.chromaFormat, Position(0, 0), areaY.size());
    PelBuf tmpRecLuma;
    if(isLuma(compID))
    {
      tmpRecLuma = m_tmpStorageLCU.getBuf(tmpArea1);
      tmpRecLuma.copyFrom(reco);
      if(m_pcEncCfg->getLmcs() && (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() ) && !(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()))
      {
        tmpRecLuma.rspSignal(m_pcReshape->getInvLUT());
      }
      
      if(cs.pps->getUseBIF() && isLuma(compID) && (cu.qp > 17))
      {
        for (auto &currTU : CU::traverseTUs(cu))
        {
          Position tuPosInCu = currTU.lumaPos() - cu.lumaPos();
          PelBuf tmpSubBuf = tmpRecLuma.subBuf(tuPosInCu, currTU.lumaSize());
          
          
          bool isInter = (cu.predMode == MODE_INTER) ? true : false;
          if ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
          {
            CompArea compArea = currTU.blocks[compID];
            PelBuf recIPredBuf = cs.slice->getPic()->getRecoBuf(compArea);
            
            // Only reshape surrounding samples if reshaping is on
            if(m_pcEncCfg->getLmcs() && (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() ) && !(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()))
            {
              m_bilateralFilter->bilateralFilterRDOdiamond5x5(tmpSubBuf, tmpSubBuf, tmpSubBuf, currTU.cu->qp, recIPredBuf, cs.slice->clpRng(compID), currTU, true, true, m_pcReshape->getInvLUT());
            }
            else
            {
              std::vector<Pel> invLUT;
              m_bilateralFilter->bilateralFilterRDOdiamond5x5(tmpSubBuf, tmpSubBuf, tmpSubBuf, currTU.cu->qp, recIPredBuf, cs.slice->clpRng(compID), currTU, true, false, invLUT);
            }
          }
        }
      }
    }
#if JVET_X0071_CHROMA_BILATERAL_FILTER
    PelBuf tmpRecChroma;
    if(isChroma(compID))
    {
      bool isCb = compID == COMPONENT_Cb ? true : false;
      const CompArea &areaUV = isCb ? cu.Cb() : cu.Cr();
      CompArea      tmpArea2(isCb ? COMPONENT_Cb : COMPONENT_Cr, areaUV.chromaFormat, Position(0, 0), areaUV.size());
      tmpRecChroma = m_tmpStorageLCU.getBuf(tmpArea2);
      tmpRecChroma.copyFrom(reco);

      if(cs.pps->getUseChromaBIF() && isChroma(compID) && (cu.qp > 17))
      {
        for (auto &currTU : CU::traverseTUs(cu))
        {
          Position tuPosInCu = currTU.chromaPos() - cu.chromaPos();
          PelBuf tmpSubBuf = tmpRecChroma.subBuf(tuPosInCu, currTU.chromaSize());
          bool isInter = (cu.predMode == MODE_INTER) ? true : false;
          if ((TU::getCbf(currTU, isCb ? COMPONENT_Cb : COMPONENT_Cr) || isInter == false))
          {
            CompArea compArea = currTU.blocks[compID];
            PelBuf recIPredBuf = cs.slice->getPic()->getRecoBuf(compArea);
            m_bilateralFilter->bilateralFilterRDOdiamond5x5Chroma(tmpSubBuf, tmpSubBuf, tmpSubBuf, currTU.cu->qp, recIPredBuf, cs.slice->clpRng(compID), currTU, true, isCb);
          }
        }
      }
    }
#endif
#if WCG_EXT
    if (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() || (
      m_pcEncCfg->getLmcs() && (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())))
    {
      const CPelBuf orgLuma = cs.getOrgBuf( cs.area.blocks[COMPONENT_Y] );
      if (compID == COMPONENT_Y )
      {
        //        if(!(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()))
        //        {
        //          tmpRecLuma.rspSignal(m_pcReshape->getInvLUT());
        //        }
        finalDistortion += m_pcRdCost->getDistPart(org, tmpRecLuma, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
      }
      else
#if JVET_X0071_CHROMA_BILATERAL_FILTER
      {
        finalDistortion += m_pcRdCost->getDistPart(org, tmpRecChroma, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
      }
#else
        finalDistortion += m_pcRdCost->getDistPart(org, reco, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
#endif
    }
    else
#endif
    {
      if (compID == COMPONENT_Y )
      {
        finalDistortion += m_pcRdCost->getDistPart( org, tmpRecLuma, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
      }
      else
      {
#if JVET_X0071_CHROMA_BILATERAL_FILTER
        finalDistortion += m_pcRdCost->getDistPart( org, tmpRecChroma, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
#else
        finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
#endif
      }
    }
#else
#if JVET_X0071_CHROMA_BILATERAL_FILTER
    PelBuf tmpRecChroma;
    if(isChroma(compID))
    {
      bool isCb = compID == COMPONENT_Cb ? true : false;
      const CompArea &areaUV = isCb ? cu.Cb() : cu.Cr();
      CompArea      tmpArea2(isCb ? COMPONENT_Cb : COMPONENT_Cr, areaUV.chromaFormat, Position(0, 0), areaUV.size());
      tmpRecChroma = m_tmpStorageLCU.getBuf(tmpArea2);
      tmpRecChroma.copyFrom(reco);
      if(cs.pps->getUseChromaBIF() && isChroma(compID) && (cu.qp > 17))
      {
        for (auto &currTU : CU::traverseTUs(cu))
        {
          Position tuPosInCu = currTU.chromaPos() - cu.chromaPos();
          PelBuf tmpSubBuf = tmpRecChroma.subBuf(tuPosInCu, currTU.chromaSize());
          bool isInter = (cu.predMode == MODE_INTER) ? true : false;
          if ((TU::getCbf(currTU, isCb ? COMPONENT_Cb : COMPONENT_Cr) || isInter == false))
          {
            CompArea compArea = currTU.blocks[compID];
            PelBuf recIPredBuf = cs.slice->getPic()->getRecoBuf(compArea);
            m_bilateralFilter->bilateralFilterRDOdiamond5x5Chroma(tmpSubBuf, tmpSubBuf, tmpSubBuf, currTU.cu->qp, recIPredBuf, cs.slice->clpRng(compID), currTU, true, isCb);
          }
        }
      }
    }
#endif
#if WCG_EXT
    if (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() || (m_pcEncCfg->getLmcs() && (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())))
    {
      const CPelBuf orgLuma = cs.getOrgBuf( cs.area.blocks[COMPONENT_Y] );
      if (compID == COMPONENT_Y && !(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()) )
      {
        const CompArea &areaY = cu.Y();
        CompArea      tmpArea1(COMPONENT_Y, areaY.chromaFormat, Position(0, 0), areaY.size());
        PelBuf tmpRecLuma = m_tmpStorageLCU.getBuf(tmpArea1);
        tmpRecLuma.rspSignal( reco, m_pcReshape->getInvLUT() );
        finalDistortion += m_pcRdCost->getDistPart(org, tmpRecLuma, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
      }
      else
      {
#if JVET_X0071_CHROMA_BILATERAL_FILTER
        if(isChroma(compID))
        {
          finalDistortion += m_pcRdCost->getDistPart(org, tmpRecChroma, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
        }
        else
        {
          finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE_WTD, &orgLuma );
        }
#else
        finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE_WTD, &orgLuma );
#endif
      }
    }
    else
#endif
    {
#if JVET_X0071_CHROMA_BILATERAL_FILTER
      if(isChroma(compID))
      {
        finalDistortion += m_pcRdCost->getDistPart( org, tmpRecChroma, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
      }
      else
      {
        finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
      }
#else
      finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
#endif
    }
#endif
  }

  cs.dist     = finalDistortion;
  cs.fracBits = finalFracBits;
  cs.cost     = m_pcRdCost->calcRdCost(cs.fracBits, cs.dist);
  if (cs.slice->getSPS()->getUseColorTrans())
  {
    if (cs.cost < cs.tmpColorSpaceCost)
    {
      cs.tmpColorSpaceCost = cs.cost;
      if (m_pcEncCfg->getRGBFormatFlag())
      {
        cs.firstColorSpaceSelected = cu.colorTransform || !cu.rootCbf;
      }
      else
      {
        cs.firstColorSpaceSelected = !cu.colorTransform || !cu.rootCbf;
      }
    }
  }

  CHECK(cs.tus.size() == 0, "No TUs present");
}

uint64_t InterSearch::xGetSymbolFracBitsInter(CodingStructure &cs, Partitioner &partitioner)
{
  uint64_t fracBits   = 0;
  CodingUnit &cu    = *cs.getCU( partitioner.chType );

  m_CABACEstimator->resetBits();

#if MULTI_HYP_PRED
  if (cu.firstPU->mergeFlag && !cu.rootCbf && cu.firstPU->numMergedAddHyps == cu.firstPU->addHypData.size())
#else
  if (cu.firstPU->mergeFlag && !cu.rootCbf)
#endif
  {
    cu.skip = true;
    CHECK(cu.colorTransform, "ACT should not be enabled for skip mode");
    m_CABACEstimator->cu_skip_flag  ( cu );
    if (cu.firstPU->ciipFlag)
    {
      // CIIP shouldn't be skip, the upper level function will deal with it, i.e. setting the overall cost to MAX_DOUBLE
    }
    else
    {
      m_CABACEstimator->merge_data(*cu.firstPU);
    }
    fracBits   += m_CABACEstimator->getEstFracBits();
  }
  else
  {
    CHECK( cu.skip, "Skip flag has to be off at this point!" );

    if (cu.Y().valid())
    m_CABACEstimator->cu_skip_flag( cu );
    m_CABACEstimator->pred_mode   ( cu );
    m_CABACEstimator->cu_pred_data( cu );
    CUCtx cuCtx;
    cuCtx.isDQPCoded = true;
    cuCtx.isChromaQpAdjCoded = true;
    m_CABACEstimator->cu_residual ( cu, partitioner, cuCtx );
    fracBits       += m_CABACEstimator->getEstFracBits();
  }

  return fracBits;
}

double InterSearch::xGetMEDistortionWeight(uint8_t bcwIdx, RefPicList eRefPicList)
{
  if( bcwIdx != BCW_DEFAULT )
  {
    return fabs((double)getBcwWeight(bcwIdx, eRefPicList) / (double)g_BcwWeightBase);
  }
  else
  {
    return 0.5;
  }
}
bool InterSearch::xReadBufferedUniMv(PredictionUnit& pu, RefPicList eRefPicList, int32_t iRefIdx, Mv& pcMvPred, Mv& rcMv, uint32_t& ruiBits, Distortion& ruiCost)
{
  if (m_uniMotions.isReadMode((uint32_t)eRefPicList, (uint32_t)iRefIdx))
  {
    m_uniMotions.copyTo(rcMv, ruiCost, (uint32_t)eRefPicList, (uint32_t)iRefIdx);

    Mv pred = pcMvPred;
    pred.changeTransPrecInternal2Amvr(pu.cu->imv);
    m_pcRdCost->setPredictor(pred);
    m_pcRdCost->setCostScale(0);

    Mv mv = rcMv;
    mv.changeTransPrecInternal2Amvr(pu.cu->imv);
    uint32_t mvBits = m_pcRdCost->getBitsOfVectorWithPredictor(mv.getHor(), mv.getVer(), 0);

    ruiBits += mvBits;
    ruiCost += m_pcRdCost->getCost(ruiBits);
    return true;
  }
  return false;
}

bool InterSearch::xReadBufferedAffineUniMv(PredictionUnit& pu, RefPicList eRefPicList, int32_t iRefIdx, Mv acMvPred[3], Mv acMv[3], uint32_t& ruiBits, Distortion& ruiCost
  , int& mvpIdx, const AffineAMVPInfo& aamvpi
)
{
  if (m_uniMotions.isReadModeAffine((uint32_t)eRefPicList, (uint32_t)iRefIdx, pu.cu->affineType))
  {
    m_uniMotions.copyAffineMvTo(acMv, ruiCost, (uint32_t)eRefPicList, (uint32_t)iRefIdx, pu.cu->affineType, mvpIdx);
    m_pcRdCost->setCostScale(0);
    acMvPred[0] = aamvpi.mvCandLT[mvpIdx];
    acMvPred[1] = aamvpi.mvCandRT[mvpIdx];
    acMvPred[2] = aamvpi.mvCandLB[mvpIdx];

    uint32_t mvBits = 0;
    for (int verIdx = 0; verIdx<(pu.cu->affineType ? 3 : 2); verIdx++)
    {
      Mv pred = verIdx ? acMvPred[verIdx] + acMv[0] - acMvPred[0] : acMvPred[verIdx];
      pred.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
      m_pcRdCost->setPredictor(pred);
      Mv mv = acMv[verIdx];
      mv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
      mvBits += m_pcRdCost->getBitsOfVectorWithPredictor(mv.getHor(), mv.getVer(), 0);
    }
    ruiBits += mvBits;
    ruiCost += m_pcRdCost->getCost(ruiBits);
    return true;
  }
  return false;
}
void InterSearch::initWeightIdxBits()
{
  for (int n = 0; n < BCW_NUM; ++n)
  {
    m_estWeightIdxBits[n] = deriveWeightIdxBits(n);
  }
}

void InterSearch::xClipMv( Mv& rcMv, const Position& pos, const struct Size& size, const SPS& sps, const PPS& pps )
{
  int mvShift = MV_FRACTIONAL_BITS_INTERNAL;
  int offset = 8;
  int horMax = ( pps.getPicWidthInLumaSamples() + offset - (int)pos.x - 1 ) << mvShift;
  int horMin = ( -( int ) sps.getMaxCUWidth()   - offset - ( int ) pos.x + 1 ) << mvShift;

  int verMax = ( pps.getPicHeightInLumaSamples() + offset - (int)pos.y - 1 ) << mvShift;
  int verMin = ( -( int ) sps.getMaxCUHeight()   - offset - ( int ) pos.y + 1 ) << mvShift;
  const SubPic &curSubPic = pps.getSubPicFromPos(pos);
  if (curSubPic.getTreatedAsPicFlag() && m_clipMvInSubPic)
  {
    horMax = ((curSubPic.getSubPicRight() + 1)  + offset - (int)pos.x - 1) << mvShift;
    horMin = (-(int)sps.getMaxCUWidth()  - offset - ((int)pos.x - curSubPic.getSubPicLeft()) + 1) << mvShift;

    verMax = ((curSubPic.getSubPicBottom() + 1) + offset -  (int)pos.y - 1) << mvShift;
    verMin = (-(int)sps.getMaxCUHeight() - offset - ((int)pos.y - curSubPic.getSubPicTop()) + 1) << mvShift;
  }
  if( pps.getWrapAroundEnabledFlag() )
  {
    int horMax = ( pps.getPicWidthInLumaSamples() + sps.getMaxCUWidth() - size.width + offset - (int)pos.x - 1 ) << mvShift;
    int horMin = ( -( int ) sps.getMaxCUWidth()                                      - offset - ( int ) pos.x + 1 ) << mvShift;
    rcMv.setHor( std::min( horMax, std::max( horMin, rcMv.getHor() ) ) );
    rcMv.setVer( std::min( verMax, std::max( verMin, rcMv.getVer() ) ) );
    return;
  }

  rcMv.setHor( std::min( horMax, std::max( horMin, rcMv.getHor() ) ) );
  rcMv.setVer( std::min( verMax, std::max( verMin, rcMv.getVer() ) ) );
}

uint32_t InterSearch::xDetermineBestMvp( PredictionUnit& pu, Mv acMvTemp[3], int& mvpIdx, const AffineAMVPInfo& aamvpi )
{
  bool mvpUpdated  = false;
  uint32_t minBits = std::numeric_limits<uint32_t>::max();
  for ( int i = 0; i < aamvpi.numCand; i++ )
  {
    Mv mvPred[3] = { aamvpi.mvCandLT[i], aamvpi.mvCandRT[i], aamvpi.mvCandLB[i] };
    uint32_t candBits = m_auiMVPIdxCost[i][aamvpi.numCand];
    candBits += xCalcAffineMVBits( pu, acMvTemp, mvPred );

    if ( candBits < minBits )
    {
      minBits    = candBits;
      mvpIdx     = i;
      mvpUpdated = true;
    }
  }
  CHECK( !mvpUpdated, "xDetermineBestMvp() error" );
  return minBits;
}

void InterSearch::symmvdCheckBestMvp(
  PredictionUnit& pu,
  PelUnitBuf& origBuf,
  Mv curMv,
  RefPicList curRefList,
  AMVPInfo amvpInfo[2][33],
  int32_t bcwIdx,
  Mv cMvPredSym[2],
  int32_t mvpIdxSym[2],
  Distortion& bestCost,
  bool skip
)
{
  RefPicList tarRefList = (RefPicList)(1 - curRefList);
  int32_t refIdxCur = pu.cu->slice->getSymRefIdx(curRefList);
  int32_t refIdxTar = pu.cu->slice->getSymRefIdx(tarRefList);

  MvField cCurMvField, cTarMvField;
  cCurMvField.setMvField(curMv, refIdxCur);
  AMVPInfo& amvpCur = amvpInfo[curRefList][refIdxCur];
  AMVPInfo& amvpTar = amvpInfo[tarRefList][refIdxTar];
  m_pcRdCost->setCostScale(0);


  // get prediction of eCurRefPicList
  PelUnitBuf predBufA = m_tmpPredStorage[curRefList].getBuf(UnitAreaRelative(*pu.cu, pu));
  const Picture* picRefA = pu.cu->slice->getRefPic(curRefList, cCurMvField.refIdx);
  Mv mvA = cCurMvField.mv;
  clipMv( mvA, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
  if ( (mvA.hor & 15) == 0 && (mvA.ver & 15) == 0 )
  {
    Position offset = pu.blocks[COMPONENT_Y].pos().offset( mvA.getHor() >> 4, mvA.getVer() >> 4 );
    CPelBuf pelBufA = picRefA->getRecoBuf( CompArea( COMPONENT_Y, pu.chromaFormat, offset, pu.blocks[COMPONENT_Y].size() ), false );
    predBufA.bufs[0].buf = const_cast<Pel *>(pelBufA.buf);
    predBufA.bufs[0].stride = pelBufA.stride;
  }
  else
  {
    xPredInterBlk( COMPONENT_Y, pu, picRefA, mvA, predBufA, false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false );
  }
  PelUnitBuf bufTmp = m_tmpStorageLCU.getBuf( UnitAreaRelative( *pu.cu, pu ) );
  bufTmp.copyFrom( origBuf );
  bufTmp.removeHighFreq( predBufA, m_pcEncCfg->getClipForBiPredMeEnabled(), pu.cu->slice->clpRngs(), getBcwWeight( pu.cu->BcwIdx, tarRefList ) );

  double fWeight = xGetMEDistortionWeight( pu.cu->BcwIdx, tarRefList );

  int32_t skipMvpIdx[2];
  skipMvpIdx[0] = skip ? mvpIdxSym[0] : -1;
  skipMvpIdx[1] = skip ? mvpIdxSym[1] : -1;

  for (int i = 0; i < amvpCur.numCand; i++)
  {
    for (int j = 0; j < amvpTar.numCand; j++)
    {
      if (skipMvpIdx[curRefList] == i && skipMvpIdx[tarRefList] == j)
        continue;

      cTarMvField.setMvField(curMv.getSymmvdMv(amvpCur.mvCand[i], amvpTar.mvCand[j]), refIdxTar);

      // get prediction of eTarRefPicList
      PelUnitBuf predBufB = m_tmpPredStorage[tarRefList].getBuf(UnitAreaRelative(*pu.cu, pu));
      const Picture* picRefB = pu.cu->slice->getRefPic(tarRefList, cTarMvField.refIdx);
      Mv mvB = cTarMvField.mv;
      clipMv( mvB, pu.cu->lumaPos(), pu.cu->lumaSize(), *pu.cs->sps, *pu.cs->pps );
      if ( (mvB.hor & 15) == 0 && (mvB.ver & 15) == 0 )
      {
        Position offset = pu.blocks[COMPONENT_Y].pos().offset( mvB.getHor() >> 4, mvB.getVer() >> 4 );
        CPelBuf pelBufB = picRefB->getRecoBuf( CompArea( COMPONENT_Y, pu.chromaFormat, offset, pu.blocks[COMPONENT_Y].size() ), false );
        predBufB.bufs[0].buf = const_cast<Pel *>(pelBufB.buf);
        predBufB.bufs[0].stride = pelBufB.stride;
      }
      else
      {
        xPredInterBlk( COMPONENT_Y, pu, picRefB, mvB, predBufB, false, pu.cu->slice->clpRng( COMPONENT_Y ), false, false );
      }
      // calc distortion
      DFunc distFunc = (!pu.cu->slice->getDisableSATDForRD()) ? DF_HAD : DF_SAD;
      Distortion cost = (Distortion)floor( fWeight * (double)m_pcRdCost->getDistPart( bufTmp.Y(), predBufB.Y(), pu.cs->sps->getBitDepth( CHANNEL_TYPE_LUMA ), COMPONENT_Y, distFunc ) );

      Mv pred = amvpCur.mvCand[i];
      pred.changeTransPrecInternal2Amvr(pu.cu->imv);
      m_pcRdCost->setPredictor(pred);
      Mv mv = curMv;
      mv.changeTransPrecInternal2Amvr(pu.cu->imv);
      uint32_t bits = m_pcRdCost->getBitsOfVectorWithPredictor(mv.hor, mv.ver, 0);
#if TM_AMVP
      bits += m_auiMVPIdxCost[i][amvpCur.numCand];
      bits += m_auiMVPIdxCost[j][amvpTar.numCand];
#else
      bits += m_auiMVPIdxCost[i][AMVP_MAX_NUM_CANDS];
      bits += m_auiMVPIdxCost[j][AMVP_MAX_NUM_CANDS];
#endif
      cost += m_pcRdCost->getCost(bits);
      if (cost < bestCost)
      {
        bestCost = cost;
        cMvPredSym[curRefList] = amvpCur.mvCand[i];
        cMvPredSym[tarRefList] = amvpTar.mvCand[j];
        mvpIdxSym[curRefList] = i;
        mvpIdxSym[tarRefList] = j;
      }
    }
  }
}

uint64_t InterSearch::xCalcPuMeBits(PredictionUnit& pu)
{
  assert(pu.mergeFlag);
  assert(!CU::isIBC(*pu.cu));
  m_CABACEstimator->resetBits();
  m_CABACEstimator->merge_flag(pu);
  if (pu.mergeFlag)
  {
    m_CABACEstimator->merge_data(pu);
#if MULTI_HYP_PRED
    m_CABACEstimator->mh_pred_data(pu);
#endif
  }
#if MULTI_HYP_PRED
  else if (pu.interDir == 3)
  {
    m_CABACEstimator->mh_pred_data(pu);
  }
#endif
  return m_CABACEstimator->getEstFracBits();
}

#if !JVET_Z0084_IBC_TM
bool InterSearch::searchBv(PredictionUnit& pu, int xPos, int yPos, int width, int height, int picWidth, int picHeight, int xBv, int yBv, int ctuSize)
{
  const int ctuSizeLog2 = floorLog2(ctuSize);

  int refRightX = xPos + xBv + width - 1;
  int refBottomY = yPos + yBv + height - 1;

  int refLeftX = xPos + xBv;
  int refTopY = yPos + yBv;

  if ((xPos + xBv) < 0)
  {
    return false;
  }
  if (refRightX >= picWidth)
  {
    return false;
  }

  if ((yPos + yBv) < 0)
  {
    return false;
  }
  if (refBottomY >= picHeight)
  {
    return false;
  }
  if ((xBv + width) > 0 && (yBv + height) > 0)
  {
    return false;
  }

#if !JVET_Z0153_IBC_EXT_REF
  // Don't search the above CTU row
  if (refTopY >> ctuSizeLog2 < yPos >> ctuSizeLog2)
    return false;
#endif

  // Don't search the below CTU row
  if (refBottomY >> ctuSizeLog2 > yPos >> ctuSizeLog2)
  {
    return false;
  }

  unsigned curTileIdx = pu.cs->pps->getTileIdx(pu.lumaPos());
  unsigned refTileIdx = pu.cs->pps->getTileIdx(Position(refLeftX, refTopY));
  if (curTileIdx != refTileIdx)
  {
    return false;
  }
  refTileIdx = pu.cs->pps->getTileIdx(Position(refLeftX, refBottomY));
  if (curTileIdx != refTileIdx)
  {
    return false;
  }
  refTileIdx = pu.cs->pps->getTileIdx(Position(refRightX, refTopY));
  if (curTileIdx != refTileIdx)
  {
    return false;
  }
  refTileIdx = pu.cs->pps->getTileIdx(Position(refRightX, refBottomY));
  if (curTileIdx != refTileIdx)
  {
    return false;
  }

#if JVET_Z0153_IBC_EXT_REF
  if ((refTopY >> ctuSizeLog2) + 2 < (yPos >> ctuSizeLog2))
  {
    return false;
  };
  if (((refTopY >> ctuSizeLog2) == (yPos >> ctuSizeLog2)) && ((refRightX >> ctuSizeLog2) > (xPos >> ctuSizeLog2)))
  {
    return false;
  }
  if (((refTopY >> ctuSizeLog2) + 2 == (yPos >> ctuSizeLog2)) && ((refLeftX >> ctuSizeLog2) + 2 < (xPos >> ctuSizeLog2)))
  {
    return false;
  }
#else
  // in the same CTU line
#if CTU_256
  int numLeftCTUs = ( 1 << ( ( MAX_CU_DEPTH - ctuSizeLog2 ) << 1 ) ) - ( ( ctuSizeLog2 < MAX_CU_DEPTH ) ? 1 : 0 );
#else
  int numLeftCTUs = (1 << ((7 - ctuSizeLog2) << 1)) - ((ctuSizeLog2 < 7) ? 1 : 0);
#endif
  if ((refRightX >> ctuSizeLog2 <= xPos >> ctuSizeLog2) && (refLeftX >> ctuSizeLog2 >= (xPos >> ctuSizeLog2) - numLeftCTUs))
  {

    // in the same CTU, or left CTU
    // if part of ref block is in the left CTU, some area can be referred from the not-yet updated local CTU buffer
#if CTU_256
    if( ( ( refLeftX >> ctuSizeLog2 ) == ( ( xPos >> ctuSizeLog2 ) - 1 ) ) && ( ctuSizeLog2 == MAX_CU_DEPTH ) )
#else
    if (((refLeftX >> ctuSizeLog2) == ((xPos >> ctuSizeLog2) - 1)) && (ctuSizeLog2 == 7))
#endif
    {
      // ref block's collocated block in current CTU
      const Position refPosCol = pu.Y().topLeft().offset(xBv + ctuSize, yBv);
      int offset64x = (refPosCol.x >> (ctuSizeLog2 - 1)) << (ctuSizeLog2 - 1);
      int offset64y = (refPosCol.y >> (ctuSizeLog2 - 1)) << (ctuSizeLog2 - 1);
      const Position refPosCol64x64 = {offset64x, offset64y};
      if (pu.cs->isDecomp(refPosCol64x64, toChannelType(COMPONENT_Y)))
        return false;
      if (refPosCol64x64 == pu.Y().topLeft())
        return false;
    }
  }
  else
    return false;
#endif

  // in the same CTU, or valid area from left CTU. Check if the reference block is already coded
  const Position refPosLT = pu.Y().topLeft().offset(xBv, yBv);
  const Position refPosBR = pu.Y().bottomRight().offset(xBv, yBv);
  const ChannelType      chType = toChannelType(COMPONENT_Y);
  if (!pu.cs->isDecomp(refPosBR, chType))
    return false;
  if (!pu.cs->isDecomp(refPosLT, chType))
    return false;
  return true;
}
#endif

//! \}
