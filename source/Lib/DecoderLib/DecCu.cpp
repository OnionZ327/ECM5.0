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

/** \file     DecCu.cpp
    \brief    CU decoder class
*/

#include "DecCu.h"

#include "CommonLib/InterPrediction.h"
#include "CommonLib/IntraPrediction.h"
#include "CommonLib/Picture.h"
#include "CommonLib/UnitTools.h"

#include "CommonLib/dtrace_buffer.h"

#if RExt__DECODER_DEBUG_TOOL_STATISTICS
#include "CommonLib/CodingStatistics.h"
#endif
#if K0149_BLOCK_STATISTICS
#include "CommonLib/ChromaFormat.h"
#include "CommonLib/dtrace_blockstatistics.h"
#endif

//! \ingroup DecoderLib
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

DecCu::DecCu()
{
  m_tmpStorageLCU = NULL;
}

DecCu::~DecCu()
{
  m_ciipBuffer.destroy();
}

void DecCu::init( TrQuant* pcTrQuant, IntraPrediction* pcIntra, InterPrediction* pcInter)
{
  m_pcTrQuant       = pcTrQuant;
  m_pcIntraPred     = pcIntra;
  m_pcInterPred     = pcInter;
  m_ciipBuffer.destroy();
  m_ciipBuffer.create(CHROMA_420, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)); // TODO: support other color format
}
void DecCu::initDecCuReshaper  (Reshape* pcReshape, ChromaFormat chromaFormatIDC)
{
  m_pcReshape = pcReshape;
  if (m_tmpStorageLCU == NULL)
  {
    m_tmpStorageLCU = new PelStorage;
    m_tmpStorageLCU->create(UnitArea(chromaFormatIDC, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
  }

}
void DecCu::destoryDecCuReshaprBuf()
{
  if (m_tmpStorageLCU)
  {
    m_tmpStorageLCU->destroy();
    delete m_tmpStorageLCU;
    m_tmpStorageLCU = NULL;
  }
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

void DecCu::decompressCtu( CodingStructure& cs, const UnitArea& ctuArea )
{

  const int maxNumChannelType = cs.pcv->chrFormat != CHROMA_400 && CS::isDualITree( cs ) ? 2 : 1;

  if (cs.resetIBCBuffer)
  {
    m_pcInterPred->resetIBCBuffer(cs.pcv->chrFormat, cs.slice->getSPS()->getMaxCUHeight());
    cs.resetIBCBuffer = false;
  }
#if JVET_Z0153_IBC_EXT_REF
  else
  {
    const int ctuSize = cs.slice->getSPS()->getMaxCUHeight();
    m_pcInterPred->resetVPDUforIBC(cs.pcv->chrFormat, ctuSize, ctuSize, ctuArea.Y().x, ctuArea.Y().y);
  }
#endif

#if JVET_Z0118_GDR
  // reset current IBC Buffer only when VB pass through
  if (cs.isInGdrIntervalOrRecoveryPoc())
  {
#if JVET_Z0153_IBC_EXT_REF    
    m_pcInterPred->resetCurIBCBuffer(
      cs.pcv->chrFormat,
      ctuArea.Y(),
      cs.slice->getSPS()->getMaxCUHeight(),        
      1 << (cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) - 1)
    );
#else
    int gdrEndX = cs.picHeader->getGdrEndX();

    if (ctuArea.lx() <= gdrEndX && gdrEndX < ctuArea.lx() + ctuArea.lwidth())
    {
      m_pcInterPred->resetCurIBCBuffer(
        cs.pcv->chrFormat,
        ctuArea.Y(),
        cs.slice->getSPS()->getMaxCUHeight(),
        1 << (cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) - 1)
      );
    }
#endif
  }
#endif

  for( int ch = 0; ch < maxNumChannelType; ch++ )
  {
    const ChannelType chType = ChannelType( ch );
    Position prevTmpPos;
    prevTmpPos.x = -1; prevTmpPos.y = -1;

    for( auto &currCU : cs.traverseCUs( CS::getArea( cs, ctuArea, chType ), chType ) )
    {
#if JVET_Z0054_BLK_REF_PIC_REORDER
      m_pcInterPred->setFillCurTplAboveARMC(false);
      m_pcInterPred->setFillCurTplLeftARMC(false);
#endif
#if JVET_Z0118_GDR
      Slice   *slice = currCU.slice;      
      Picture *refPic;
      
      bool isInGdrInterval  = slice->getPicHeader()->getInGdrInterval();
      bool isRecoveryPocPic = slice->getPicHeader()->getIsGdrRecoveryPocPic();
      bool isNorPicture = !(isInGdrInterval || isRecoveryPocPic) && slice->isInterB();

      if (isNorPicture) 
      {        
        currCU.cs->setReconBuf(PIC_RECONSTRUCTION_0);
        currCU.cs->picture->setCleanDirty(false);

        // 1.01 use only dirty reference picture
        for (int rlist = REF_PIC_LIST_0; rlist < NUM_REF_PIC_LIST_01; rlist++)
        {
          int n = slice->getNumRefIdx((RefPicList)rlist);
          for (int idx = 0; idx < n; idx++)
          {
            Picture *refPic = slice->getReferencePicture((RefPicList)rlist, idx);
            if (refPic)
            {
              // when cur picture is normal picture and ref picture is gdr/recovery picture
              // note: pic.slice.picHeader and pic.cs.picHeader could be different    
              // bool isGdrPic = refPic->cs->picHeader->getInGdrPeriod();
              bool isRefInGdrInterval  = refPic->cs->picHeader->getInGdrInterval();
              bool isRefRecoveryPocPic = refPic->cs->picHeader->getIsGdrRecoveryPocPic();

              if (isRefInGdrInterval || isRefRecoveryPocPic)
              {
                refPic->setCleanDirty(true);                                
              }
              else 
              {
                refPic->setCleanDirty(false);                                
              }
            }          
          }
        }
      }

      // 1.02 use clean reference pictures for some CU when gdr starts
      if (isInGdrInterval || isRecoveryPocPic)
      {
        // 1.01 switch recon based on clean/dirty current area
        bool cleanDirtyFlag;

        bool isCuClean = currCU.Y().valid() ? cs.isClean(currCU.Y().topLeft(), CHANNEL_TYPE_LUMA) : cs.isClean(currCU.Cb().topLeft(), CHANNEL_TYPE_CHROMA);

        if (isCuClean) 
        {
          cleanDirtyFlag = true;
        }
        else 
        {
          cleanDirtyFlag = false;
        }       
        
        currCU.cs->setReconBuf(( cleanDirtyFlag ) ? PIC_RECONSTRUCTION_1 : PIC_RECONSTRUCTION_0);
        currCU.cs->picture->setCleanDirty( cleanDirtyFlag );

        for (int rlist = REF_PIC_LIST_0; rlist < NUM_REF_PIC_LIST_01; rlist++)
        {
          int n = slice->getNumRefIdx((RefPicList)rlist);
          for (int idx = 0; idx < n; idx++)
          {
            refPic = slice->getReferencePicture((RefPicList)rlist, idx);            
            if (refPic)
            {
              refPic->setCleanDirty( cleanDirtyFlag );
            }
          }
        }
      }
#endif

#if !REMOVE_VPDU
      if(currCU.Y().valid())
      {
        const int vSize = cs.slice->getSPS()->getMaxCUHeight() > 64 ? 64 : cs.slice->getSPS()->getMaxCUHeight();
        if((currCU.Y().x % vSize) == 0 && (currCU.Y().y % vSize) == 0)
        {
          for(int x = currCU.Y().x; x < currCU.Y().x + currCU.Y().width; x += vSize)
          {
            for(int y = currCU.Y().y; y < currCU.Y().y + currCU.Y().height; y += vSize)
            {
              m_pcInterPred->resetVPDUforIBC(cs.pcv->chrFormat, cs.slice->getSPS()->getMaxCUHeight(), vSize, x + g_IBCBufferSize / cs.slice->getSPS()->getMaxCUHeight() / 2, y);
            }
          }
        }
      }
#endif
      if (currCU.predMode != MODE_INTRA && currCU.predMode != MODE_PLT && currCU.Y().valid())
      {
        xDeriveCUMV(currCU);
#if K0149_BLOCK_STATISTICS
        if(currCU.geoFlag)
        {
          storeGeoMergeCtx(m_geoMrgCtx);
        }
#endif
      }
      switch( currCU.predMode )
      {
      case MODE_INTER:
      case MODE_IBC:
#if JVET_Y0065_GPM_INTRA
#if ENABLE_DIMD && JVET_W0123_TIMD_FUSION
        if ((cs.slice->getSPS()->getUseDimd() || cs.slice->getSPS()->getUseTimd()) && currCU.geoFlag && currCU.firstPU->gpmIntraFlag)
#elif ENABLE_DIMD
        if (cs.slice->getSPS()->getUseDimd() && currCU.geoFlag && currCU.firstPU->gpmIntraFlag)
#elif JVET_W0123_TIMD_FUSION
        if (cs.slice->getSPS()->getUseTimd() && currCU.geoFlag && currCU.firstPU->gpmIntraFlag)
#endif
        {
          if ((int)(currCU.firstPU->geoMergeIdx0)-GEO_MAX_NUM_UNI_CANDS > 0 || (int)(currCU.firstPU->geoMergeIdx1)-GEO_MAX_NUM_UNI_CANDS > 0) // dimd/timd
          {
            const CompArea &area = currCU.Y();
#if ENABLE_DIMD
#if JVET_W0123_TIMD_FUSION
            if (cs.slice->getSPS()->getUseDimd())
#endif
            {
              IntraPrediction::deriveDimdMode(currCU.cs->picture->getRecoBuf(area), area, currCU);
            }
#endif
#if JVET_W0123_TIMD_FUSION
#if ENABLE_DIMD
            if (cs.slice->getSPS()->getUseTimd())
#endif
            {
              currCU.timdMode = m_pcIntraPred->deriveTimdMode(currCU.cs->picture->getRecoBuf(area), area, currCU);
            }
#endif
          }
        }
#endif
        xReconInter( currCU );
        break;
      case MODE_PLT:
      case MODE_INTRA:
#if ENABLE_DIMD
        if (currCU.dimd)
        {
          PredictionUnit *pu = currCU.firstPU;
          const CompArea &area = currCU.Y();
          IntraPrediction::deriveDimdMode(currCU.cs->picture->getRecoBuf(area), area, currCU);
          pu->intraDir[0] = currCU.dimdMode;
        }
#if JVET_W0123_TIMD_FUSION
        else if (currCU.timd)
        {
          PredictionUnit *pu = currCU.firstPU;
          const CompArea &area = currCU.Y();
#if SECONDARY_MPM
          IntraPrediction::deriveDimdMode(currCU.cs->picture->getRecoBuf(area), area, currCU);
#endif
          currCU.timdMode = m_pcIntraPred->deriveTimdMode(currCU.cs->picture->getRecoBuf(area), area, currCU);
          pu->intraDir[0] = currCU.timdMode;
        }
#endif
        else if (currCU.firstPU->parseLumaMode)
        {
          const CompArea &area = currCU.Y();
          IntraPrediction::deriveDimdMode(currCU.cs->picture->getRecoBuf(area), area, currCU);
        }

        //redo prediction dir derivation
        if (currCU.firstPU->parseLumaMode)
        {
#if SECONDARY_MPM
          uint8_t* mpm_pred = currCU.firstPU->intraMPM;  // mpm_idx / rem_intra_luma_pred_mode
          uint8_t* non_mpm_pred = currCU.firstPU->intraNonMPM;
          PU::getIntraMPMs( *currCU.firstPU, mpm_pred, non_mpm_pred );
#else
          unsigned int mpm_pred[NUM_MOST_PROBABLE_MODES];  // mpm_idx / rem_intra_luma_pred_mode
          PU::getIntraMPMs(*currCU.firstPU, mpm_pred);
#endif
          if (currCU.firstPU->mpmFlag)
          {
            currCU.firstPU->intraDir[0] = mpm_pred[currCU.firstPU->ipred_idx];
          }
          else
          {
#if SECONDARY_MPM
            if (currCU.firstPU->secondMpmFlag)
            {
              currCU.firstPU->intraDir[0] = mpm_pred[currCU.firstPU->ipred_idx];
            }
            else
            {
              currCU.firstPU->intraDir[0] = non_mpm_pred[currCU.firstPU->ipred_idx];
            }
#else
            //postponed sorting of MPMs (only in remaining branch)
            std::sort(mpm_pred, mpm_pred + NUM_MOST_PROBABLE_MODES);
            unsigned ipred_mode = currCU.firstPU->ipred_idx;

            for (uint32_t i = 0; i < NUM_MOST_PROBABLE_MODES; i++)
            {
              ipred_mode += (ipred_mode >= mpm_pred[i]);
            }
            currCU.firstPU->intraDir[0] = ipred_mode;
#endif
          }
        }
        if (currCU.firstPU->parseChromaMode)
        {
          unsigned chromaCandModes[NUM_CHROMA_MODE];
          PU::getIntraChromaCandModes(*currCU.firstPU, chromaCandModes);

          CHECK(currCU.firstPU->candId >= NUM_CHROMA_MODE, "Chroma prediction mode index out of bounds");
          CHECK(PU::isLMCMode(chromaCandModes[currCU.firstPU->candId]), "The intra dir cannot be LM_CHROMA for this path");
          CHECK(chromaCandModes[currCU.firstPU->candId] == DM_CHROMA_IDX, "The intra dir cannot be DM_CHROMA for this path");
#if JVET_Z0050_DIMD_CHROMA_FUSION && ENABLE_DIMD
          CHECK(chromaCandModes[currCU.firstPU->candId] == DIMD_CHROMA_IDX, "The intra dir cannot be DIMD_CHROMA for this path");
#endif

          currCU.firstPU->intraDir[1] = chromaCandModes[currCU.firstPU->candId];
        }
#else
#if JVET_W0123_TIMD_FUSION
        if (currCU.timd)
        {
          PredictionUnit *pu = currCU.firstPU;
          const CompArea &area = currCU.Y();
          currCU.timdMode = m_pcIntraPred->deriveTimdMode(currCU.cs->picture->getRecoBuf(area), area, currCU);
          pu->intraDir[0] = currCU.timdMode;
        }

        //redo prediction dir derivation
        if (currCU.firstPU->parseLumaMode)
        {
#if SECONDARY_MPM
          uint8_t* mpm_pred = currCU.firstPU->intraMPM;  // mpm_idx / rem_intra_luma_pred_mode
          uint8_t* non_mpm_pred = currCU.firstPU->intraNonMPM;
          PU::getIntraMPMs( *currCU.firstPU, mpm_pred, non_mpm_pred );
#else
          unsigned int mpm_pred[NUM_MOST_PROBABLE_MODES];  // mpm_idx / rem_intra_luma_pred_mode
          PU::getIntraMPMs(*currCU.firstPU, mpm_pred);
#endif
          if (currCU.firstPU->mpmFlag)
          {
            currCU.firstPU->intraDir[0] = mpm_pred[currCU.firstPU->ipred_idx];
          }
          else
          {
#if SECONDARY_MPM
            if (currCU.firstPU->secondMpmFlag)
            {
              currCU.firstPU->intraDir[0] = mpm_pred[currCU.firstPU->ipred_idx];
            }
            else
            {
              currCU.firstPU->intraDir[0] = non_mpm_pred[currCU.firstPU->ipred_idx];
            }
#else
            //postponed sorting of MPMs (only in remaining branch)
            std::sort(mpm_pred, mpm_pred + NUM_MOST_PROBABLE_MODES);
            unsigned ipred_mode = currCU.firstPU->ipred_idx;

            for (uint32_t i = 0; i < NUM_MOST_PROBABLE_MODES; i++)
            {
              ipred_mode += (ipred_mode >= mpm_pred[i]);
            }
            currCU.firstPU->intraDir[0] = ipred_mode;
#endif
          }
        }
        if (currCU.firstPU->parseChromaMode)
        {
          unsigned chromaCandModes[NUM_CHROMA_MODE];
          PU::getIntraChromaCandModes(*currCU.firstPU, chromaCandModes);

          CHECK(currCU.firstPU->candId >= NUM_CHROMA_MODE, "Chroma prediction mode index out of bounds");
          CHECK(PU::isLMCMode(chromaCandModes[currCU.firstPU->candId]), "The intra dir cannot be LM_CHROMA for this path");
          CHECK(chromaCandModes[currCU.firstPU->candId] == DM_CHROMA_IDX, "The intra dir cannot be DM_CHROMA for this path");

          currCU.firstPU->intraDir[1] = chromaCandModes[currCU.firstPU->candId];
        }
#endif
#endif
        xReconIntraQT( currCU );
        break;
      default:
        THROW( "Invalid prediction mode" );
        break;
      }

      m_pcInterPred->xFillIBCBuffer(currCU);
#if JVET_Z0118_GDR // decompressCtu
      cs.updateReconMotIPM( currCU ); // decompressCtu : need
#endif

      DTRACE_BLOCK_REC( cs.picture->getRecoBuf( currCU ), currCU, currCU.predMode );
      if (CU::isInter(currCU))
      {
        DTRACE_MOT_FIELD(g_trace_ctx, *currCU.firstPU);
      }
    }
  }
#if K0149_BLOCK_STATISTICS
  getAndStoreBlockStatistics(cs, ctuArea);
#endif
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

void DecCu::xIntraRecBlk( TransformUnit& tu, const ComponentID compID )
{
  if( !tu.blocks[ compID ].valid() )
  {
    return;
  }

        CodingStructure &cs = *tu.cs;
  const CompArea &area      = tu.blocks[compID];
#if SIGN_PREDICTION
  const bool  isJCCR = tu.jointCbCr && isChroma(compID);
  const CompArea &areaCr      = tu.blocks[isJCCR ? COMPONENT_Cr : compID];
  PelBuf piPredCr;
  if(isJCCR)
  {
    piPredCr = cs.getPredBuf( tu.blocks[COMPONENT_Cr] );
  }
#endif

  const ChannelType chType  = toChannelType( compID );

        PelBuf piPred       = cs.getPredBuf( area );

  const PredictionUnit &pu  = *tu.cs->getPU( area.pos(), chType );
#if ENABLE_DIMD
#if JVET_Z0050_DIMD_CHROMA_FUSION && ENABLE_DIMD
  if (pu.intraDir[1] == DIMD_CHROMA_IDX && compID == COMPONENT_Cb)
  {
    CompArea areaCb = pu.Cb();
    CompArea areaCr = pu.Cr();
    CompArea lumaArea = CompArea(COMPONENT_Y, pu.chromaFormat, areaCb.lumaPos(), recalcSize(pu.chromaFormat, CHANNEL_TYPE_CHROMA, CHANNEL_TYPE_LUMA, areaCb.size()));
    IntraPrediction::deriveDimdChromaMode(cs.picture->getRecoBuf(lumaArea), cs.picture->getRecoBuf(areaCb), cs.picture->getRecoBuf(areaCr), lumaArea, areaCb, areaCr, *pu.cu);
  }
#endif
  uint32_t uiChFinalMode = PU::getFinalIntraMode(pu, chType);
#else
  const uint32_t uiChFinalMode = PU::getFinalIntraMode(pu, chType);
#endif
  PelBuf pReco              = cs.getRecoBuf(area);

  //===== init availability pattern =====
  bool predRegDiffFromTB = CU::isPredRegDiffFromTB(*tu.cu, compID);
  bool firstTBInPredReg = CU::isFirstTBInPredReg(*tu.cu, compID, area);
  CompArea areaPredReg(COMPONENT_Y, tu.chromaFormat, area);
#if SIGN_PREDICTION
  if( !isJCCR || compID != COMPONENT_Cr )
  {
#endif
  if (tu.cu->ispMode && isLuma(compID))
  {
    if (predRegDiffFromTB)
    {
      if (firstTBInPredReg)
      {
        CU::adjustPredArea(areaPredReg);
        m_pcIntraPred->initIntraPatternChTypeISP(*tu.cu, areaPredReg, pReco);
      }
    }
    else
    {
      m_pcIntraPred->initIntraPatternChTypeISP(*tu.cu, area, pReco);
    }
  }
  else
  {
    m_pcIntraPred->initIntraPatternChType(*tu.cu, area);
  }

  //===== get prediction signal =====
  if( compID != COMPONENT_Y && PU::isLMCMode( uiChFinalMode ) )
  {
    const PredictionUnit& pu = *tu.cu->firstPU;
    m_pcIntraPred->xGetLumaRecPixels( pu, area );
    m_pcIntraPred->predIntraChromaLM( compID, piPred, pu, area, uiChFinalMode );
  }
  else
  {
#if JVET_V0130_INTRA_TMP
	  if (PU::isTmp(pu, chType))
	  {
		  int foundCandiNum;
#if JVET_W0069_TMP_BOUNDARY
		  RefTemplateType tempType = m_pcIntraPred->getRefTemplateType(*(tu.cu), tu.cu->blocks[COMPONENT_Y]);

      if( tempType != NO_TEMPLATE )
		  {
        m_pcIntraPred->getTargetTemplate(tu.cu, pu.lwidth(), pu.lheight(), tempType);
        m_pcIntraPred->candidateSearchIntra(tu.cu, pu.lwidth(), pu.lheight(), tempType);
        m_pcIntraPred->generateTMPrediction(piPred.buf, piPred.stride, pu.lwidth(), pu.lheight(), foundCandiNum);
		  }
		  else
		  {
			  foundCandiNum = 1;
        m_pcIntraPred->generateTmDcPrediction(piPred.buf, piPred.stride, pu.lwidth(), pu.lheight(), 1 << (tu.cu->cs->sps->getBitDepth(CHANNEL_TYPE_LUMA) - 1));
		  }
#else
      m_pcIntraPred->getTargetTemplate(tu.cu, pu.lwidth(), pu.lheight());
      m_pcIntraPred->candidateSearchIntra(tu.cu, pu.lwidth(), pu.lheight());
      m_pcIntraPred->generateTMPrediction(piPred.buf, piPred.stride, pu.lwidth(), pu.lheight(), foundCandiNum);
#endif
		  assert(foundCandiNum >= 1);
	  }
	  else if (PU::isMIP(pu, chType))
#else
    if( PU::isMIP( pu, chType ) )
#endif
    {
      m_pcIntraPred->initIntraMip( pu, area );
      m_pcIntraPred->predIntraMip( compID, piPred, pu );
    }
    else
    {
      if (predRegDiffFromTB)
      {
        if (firstTBInPredReg)
        {
          PelBuf piPredReg = cs.getPredBuf(areaPredReg);
          m_pcIntraPred->predIntraAng(compID, piPredReg, pu);
        }
      }
      else
        m_pcIntraPred->predIntraAng(compID, piPred, pu);
#if JVET_Z0050_DIMD_CHROMA_FUSION
      if (compID != COMPONENT_Y && pu.isChromaFusion)
      {
        m_pcIntraPred->geneChromaFusionPred(compID, piPred, pu);
      }
#endif
    }
  }
#if SIGN_PREDICTION
  if(isJCCR && compID == COMPONENT_Cb)
  {
    m_pcIntraPred->initIntraPatternChType(*tu.cu, areaCr);
    if( PU::isLMCMode( uiChFinalMode ) )
    {
      const PredictionUnit& pu = *tu.cu->firstPU;
      m_pcIntraPred->xGetLumaRecPixels( pu, areaCr );
      m_pcIntraPred->predIntraChromaLM( COMPONENT_Cr, piPredCr, pu, areaCr, uiChFinalMode );
    }
    else
    {
      if( PU::isMIP( pu, chType ) )
      {
        m_pcIntraPred->initIntraMip( pu, area );
        m_pcIntraPred->predIntraMip( COMPONENT_Cr, piPredCr, pu );
      }
      else
      {
        m_pcIntraPred->predIntraAng(COMPONENT_Cr, piPredCr, pu);
#if JVET_Z0050_DIMD_CHROMA_FUSION
        if (pu.isChromaFusion)
        {
          m_pcIntraPred->geneChromaFusionPred(COMPONENT_Cr, piPredCr, pu);
        }
#endif
      }
    }
  }
  }
#endif
  const Slice           &slice = *cs.slice;
  bool flag = slice.getLmcsEnabledFlag() && (slice.isIntra() || (!slice.isIntra() && m_pcReshape->getCTUFlag()));
  if (flag && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && (compID != COMPONENT_Y) && (tu.cbf[COMPONENT_Cb] || tu.cbf[COMPONENT_Cr]))
  {
#if LMCS_CHROMA_CALC_CU
    const Area area = tu.cu->Y().valid() ? tu.cu->Y() : Area(recalcPosition(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.cu->blocks[tu.chType].pos()), recalcSize(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.cu->blocks[tu.chType].size()));
#else
    const Area area = tu.Y().valid() ? tu.Y() : Area(recalcPosition(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.blocks[tu.chType].pos()), recalcSize(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.blocks[tu.chType].size()));
#endif
    const CompArea &areaY = CompArea(COMPONENT_Y, tu.chromaFormat, area);
    int adj = m_pcReshape->calculateChromaAdjVpduNei(tu, areaY);
    tu.setChromaAdj(adj);
  }
#if SIGN_PREDICTION
  flag = flag && (tu.blocks[compID].width*tu.blocks[compID].height > 4);
#endif
  //===== inverse transform =====
  PelBuf piResi = cs.getResiBuf( area );

  const QpParam cQP( tu, compID );

#if SIGN_PREDICTION
  bool bJccrWithCr = tu.jointCbCr && !(tu.jointCbCr >> 1);
  bool bIsJccr     = tu.jointCbCr && isChroma(compID);
  ComponentID signPredCompID = bIsJccr ? (bJccrWithCr ? COMPONENT_Cr : COMPONENT_Cb): compID;
  bool reshapeChroma = flag && (TU::getCbf(tu, signPredCompID) || tu.jointCbCr) && isChroma(signPredCompID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag();
  m_pcTrQuant->predCoeffSigns(tu, compID, reshapeChroma);
#endif
  if( tu.jointCbCr && isChroma(compID) )
  {
    if( compID == COMPONENT_Cb )
    {
      PelBuf resiCr = cs.getResiBuf( tu.blocks[ COMPONENT_Cr ] );
      if( tu.jointCbCr >> 1 )
      {
        m_pcTrQuant->invTransformNxN( tu, COMPONENT_Cb, piResi, cQP );
      }
      else
      {
        const QpParam qpCr( tu, COMPONENT_Cr );
        m_pcTrQuant->invTransformNxN( tu, COMPONENT_Cr, resiCr, qpCr );
      }
      m_pcTrQuant->invTransformICT( tu, piResi, resiCr );
    }
  }
  else
  if( TU::getCbf( tu, compID ) )
  {
    m_pcTrQuant->invTransformNxN( tu, compID, piResi, cQP );
  }
  else
  {
    piResi.fill( 0 );
  }

  //===== reconstruction =====
#if !SIGN_PREDICTION
  flag = flag && (tu.blocks[compID].width*tu.blocks[compID].height > 4);
#endif
  if (flag && (TU::getCbf(tu, compID) || tu.jointCbCr) && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag())
  {
    piResi.scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(compID));
  }

  if( !tu.cu->ispMode || !isLuma( compID ) )
  {
    cs.setDecomp( area );
  }
  else if( tu.cu->ispMode && isLuma( compID ) && CU::isISPFirst( *tu.cu, tu.blocks[compID], compID ) )
  {
    cs.setDecomp( tu.cu->blocks[compID] );
  }

#if REUSE_CU_RESULTS
  CompArea    tmpArea(COMPONENT_Y, area.chromaFormat, Position(0, 0), area.size());
  PelBuf tmpPred;
#endif
  if (slice.getLmcsEnabledFlag() && (m_pcReshape->getCTUFlag() || slice.isIntra()) && compID == COMPONENT_Y)
  {
#if REUSE_CU_RESULTS
    {
      tmpPred = m_tmpStorageLCU->getBuf(tmpArea);
      tmpPred.copyFrom(piPred);
    }
#endif
  }
#if KEEP_PRED_AND_RESI_SIGNALS
  pReco.reconstruct( piPred, piResi, tu.cu->cs->slice->clpRng( compID ) );
#else
  piPred.reconstruct( piPred, piResi, tu.cu->cs->slice->clpRng( compID ) );
#endif
#if !KEEP_PRED_AND_RESI_SIGNALS
  pReco.copyFrom( piPred );
#endif
  if (slice.getLmcsEnabledFlag() && (m_pcReshape->getCTUFlag() || slice.isIntra()) && compID == COMPONENT_Y)
  {
#if REUSE_CU_RESULTS
    {
      piPred.copyFrom(tmpPred);
    }
#endif
  }
#if REUSE_CU_RESULTS || SIGN_PREDICTION
#if !SIGN_PREDICTION
  if( cs.pcv->isEncoder )
#endif
  {
#if JVET_Z0118_GDR
    cs.updateReconMotIPM(area, pReco);
#else
    cs.picture->getRecoBuf( area ).copyFrom( pReco );
#endif

    cs.picture->getPredBuf(area).copyFrom(piPred);
  }
#endif
}

void DecCu::xIntraRecACTBlk(TransformUnit& tu)
{
  CodingStructure      &cs = *tu.cs;
  const PredictionUnit &pu = *tu.cs->getPU(tu.blocks[COMPONENT_Y], CHANNEL_TYPE_LUMA);
  const Slice          &slice = *cs.slice;

  CHECK(!tu.Y().valid() || !tu.Cb().valid() || !tu.Cr().valid(), "Invalid TU");
  CHECK(&pu != tu.cu->firstPU, "wrong PU fetch");
  CHECK(tu.cu->ispMode, "adaptive color transform cannot be applied to ISP");
  CHECK(pu.intraDir[CHANNEL_TYPE_CHROMA] != DM_CHROMA_IDX, "chroma should use DM mode for adaptive color transform");

  bool flag = slice.getLmcsEnabledFlag() && (slice.isIntra() || (!slice.isIntra() && m_pcReshape->getCTUFlag()));
#if JVET_S0234_ACT_CRS_FIX
  if (flag && slice.getPicHeader()->getLmcsChromaResidualScaleFlag())
#else
  if (flag && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && (tu.cbf[COMPONENT_Cb] || tu.cbf[COMPONENT_Cr]))
#endif
  {
#if LMCS_CHROMA_CALC_CU
    const Area      area = tu.cu->Y().valid() ? tu.cu->Y() : Area(recalcPosition(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.cu->blocks[tu.chType].pos()), recalcSize(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.cu->blocks[tu.chType].size()));
#else
    const Area      area = tu.Y().valid() ? tu.Y() : Area(recalcPosition(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.blocks[tu.chType].pos()), recalcSize(tu.chromaFormat, tu.chType, CHANNEL_TYPE_LUMA, tu.blocks[tu.chType].size()));
#endif
    const CompArea &areaY = CompArea(COMPONENT_Y, tu.chromaFormat, area);
    int            adj = m_pcReshape->calculateChromaAdjVpduNei(tu, areaY);
    tu.setChromaAdj(adj);
  }

  for (int i = 0; i < getNumberValidComponents(tu.chromaFormat); i++)
  {
    ComponentID          compID = (ComponentID)i;
    const CompArea       &area = tu.blocks[compID];
    const ChannelType    chType = toChannelType(compID);

    PelBuf piPred = cs.getPredBuf(area);
    m_pcIntraPred->initIntraPatternChType(*tu.cu, area);
#if JVET_V0130_INTRA_TMP && ! JVET_W0069_TMP_BOUNDARY
	if (PU::isTmp(pu, chType))
	{
		int foundCandiNum;
		const unsigned int uiStride = cs.picture->getRecoBuf(COMPONENT_Y).stride;
    m_pcIntraPred->getTargetTemplate(tu.cu, pu.lwidth(), pu.lheight());
    m_pcIntraPred->candidateSearchIntra(tu.cu, pu.lwidth(), pu.lheight());
    m_pcIntraPred->generateTMPrediction(piPred.buf, uiStride, pu.lwidth(), pu.lheight(), foundCandiNum);
	}
	else if (PU::isMIP(pu, chType))
#else
    if (PU::isMIP(pu, chType))
#endif
    {
      m_pcIntraPred->initIntraMip(pu, area);
      m_pcIntraPred->predIntraMip(compID, piPred, pu);
    }
    else
    {
      m_pcIntraPred->predIntraAng(compID, piPred, pu);
    }

    PelBuf piResi = cs.getResiBuf(area);

    QpParam cQP(tu, compID);

    if (tu.jointCbCr && isChroma(compID))
    {
      if (compID == COMPONENT_Cb)
      {
        PelBuf resiCr = cs.getResiBuf(tu.blocks[COMPONENT_Cr]);
        if (tu.jointCbCr >> 1)
        {
          m_pcTrQuant->invTransformNxN(tu, COMPONENT_Cb, piResi, cQP);
        }
        else
        {
          QpParam qpCr(tu, COMPONENT_Cr);

          m_pcTrQuant->invTransformNxN(tu, COMPONENT_Cr, resiCr, qpCr);
        }
        m_pcTrQuant->invTransformICT(tu, piResi, resiCr);
      }
    }
    else
    {
      if (TU::getCbf(tu, compID))
      {
        m_pcTrQuant->invTransformNxN(tu, compID, piResi, cQP);
      }
      else
      {
        piResi.fill(0);
      }
    }

#if !JVET_S0234_ACT_CRS_FIX
    flag = flag && (tu.blocks[compID].width*tu.blocks[compID].height > 4);
    if (flag && (TU::getCbf(tu, compID) || tu.jointCbCr) && isChroma(compID) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag())
    {
      piResi.scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(compID));
    }
#endif

    cs.setDecomp(area);
  }

  cs.getResiBuf(tu).colorSpaceConvert(cs.getResiBuf(tu), false, tu.cu->cs->slice->clpRng(COMPONENT_Y));

  for (int i = 0; i < getNumberValidComponents(tu.chromaFormat); i++)
  {
    ComponentID          compID = (ComponentID)i;
    const CompArea       &area = tu.blocks[compID];

    PelBuf piPred = cs.getPredBuf(area);
    PelBuf piResi = cs.getResiBuf(area);
    PelBuf piReco = cs.getRecoBuf(area);

    PelBuf tmpPred;
    if (slice.getLmcsEnabledFlag() && (m_pcReshape->getCTUFlag() || slice.isIntra()) && compID == COMPONENT_Y)
    {
      CompArea tmpArea(COMPONENT_Y, area.chromaFormat, Position(0, 0), area.size());
      tmpPred = m_tmpStorageLCU->getBuf(tmpArea);
      tmpPred.copyFrom(piPred);
    }

#if JVET_S0234_ACT_CRS_FIX
    if (flag && isChroma(compID) && (tu.blocks[compID].width*tu.blocks[compID].height > 4) && slice.getPicHeader()->getLmcsChromaResidualScaleFlag())
    {
      piResi.scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(compID));
    }
#endif
    piPred.reconstruct(piPred, piResi, tu.cu->cs->slice->clpRng(compID));
    piReco.copyFrom(piPred);

    if (slice.getLmcsEnabledFlag() && (m_pcReshape->getCTUFlag() || slice.isIntra()) && compID == COMPONENT_Y)
    {
      piPred.copyFrom(tmpPred);
    }

    if (cs.pcv->isEncoder)
    {
#if JVET_Z0118_GDR
      cs.updateReconMotIPM(area, piReco);
#else
      cs.picture->getRecoBuf(area).copyFrom(piReco);
#endif

      cs.picture->getPredBuf(area).copyFrom(piPred);
    }
  }
}

void DecCu::xReconIntraQT( CodingUnit &cu )
{

  if (CU::isPLT(cu))
  {
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    if (CS::isDualITree(*cu.cs))
#else
    if (cu.isSepTree())
#endif
    {
      if (cu.chType == CHANNEL_TYPE_LUMA)
      {
        xReconPLT(cu, COMPONENT_Y, 1);
      }
      if (cu.chromaFormat != CHROMA_400 && (cu.chType == CHANNEL_TYPE_CHROMA))
      {
        xReconPLT(cu, COMPONENT_Cb, 2);
      }
    }
    else
    {
      if( cu.chromaFormat != CHROMA_400 )
      {
        xReconPLT(cu, COMPONENT_Y, 3);
      }
      else
      {
        xReconPLT(cu, COMPONENT_Y, 1);
      }
    }
    return;
  }

  if (cu.colorTransform)
  {
    xIntraRecACTQT(cu);
  }
  else
  {
    const uint32_t numChType = ::getNumberValidChannels(cu.chromaFormat);

    for (uint32_t chType = CHANNEL_TYPE_LUMA; chType < numChType; chType++)
    {
      if (cu.blocks[chType].valid())
      {
        xIntraRecQT(cu, ChannelType(chType));
      }
    }
  }
#if JVET_W0123_TIMD_FUSION
  if (cu.blocks[CHANNEL_TYPE_LUMA].valid())
  {
    PU::spanIpmInfoIntra(*cu.firstPU);
  }
#endif
}

void DecCu::xReconPLT(CodingUnit &cu, ComponentID compBegin, uint32_t numComp)
{
  const SPS&       sps = *(cu.cs->sps);
  TransformUnit&   tu = *cu.firstTU;
  PelBuf    curPLTIdx = tu.getcurPLTIdx(compBegin);

  uint32_t height = cu.block(compBegin).height;
  uint32_t width = cu.block(compBegin).width;

  //recon. pixels
  uint32_t scaleX = getComponentScaleX(COMPONENT_Cb, sps.getChromaFormatIdc());
  uint32_t scaleY = getComponentScaleY(COMPONENT_Cb, sps.getChromaFormatIdc());
  for (uint32_t y = 0; y < height; y++)
  {
    for (uint32_t x = 0; x < width; x++)
    {
      for (uint32_t compID = compBegin; compID < (compBegin + numComp); compID++)
      {
        const int  channelBitDepth = cu.cs->sps->getBitDepth(toChannelType((ComponentID)compID));
        const CompArea &area = cu.blocks[compID];

        PelBuf       picReco   = cu.cs->getRecoBuf(area);
        PLTescapeBuf escapeValue = tu.getescapeValue((ComponentID)compID);
        if (curPLTIdx.at(x, y) == cu.curPLTSize[compBegin])
        {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
          TCoeff value;
#else
          Pel value;
#endif
          QpParam cQP(tu, (ComponentID)compID);
          int qp = cQP.Qp(true);
          int qpRem = qp % 6;
          int qpPer = qp / 6;
          if (compBegin != COMPONENT_Y || compID == COMPONENT_Y)
          {
            int invquantiserRightShift = IQUANT_SHIFT;
            int add = 1 << (invquantiserRightShift - 1);
            value = ((((escapeValue.at(x, y)*g_invQuantScales[0][qpRem]) << qpPer) + add) >> invquantiserRightShift);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
            value = ClipBD<TCoeff>(value, channelBitDepth);
            picReco.at(x, y) = Pel(value);
#else
            value = Pel(ClipBD<int>(value, channelBitDepth));
            picReco.at(x, y) = value;
#endif
          }
          else if (compBegin == COMPONENT_Y && compID != COMPONENT_Y && y % (1 << scaleY) == 0 && x % (1 << scaleX) == 0)
          {
            uint32_t posYC = y >> scaleY;
            uint32_t posXC = x >> scaleX;
            int invquantiserRightShift = IQUANT_SHIFT;
            int add = 1 << (invquantiserRightShift - 1);
            value = ((((escapeValue.at(posXC, posYC)*g_invQuantScales[0][qpRem]) << qpPer) + add) >> invquantiserRightShift);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
            value = ClipBD<TCoeff>(value, channelBitDepth);
            picReco.at(posXC, posYC) = Pel(value);
#else
            value = Pel(ClipBD<int>(value, channelBitDepth));
            picReco.at(posXC, posYC) = value;
#endif

          }
        }
        else
        {
          uint32_t curIdx = curPLTIdx.at(x, y);
          if (compBegin != COMPONENT_Y || compID == COMPONENT_Y)
          {
            picReco.at(x, y) = cu.curPLT[compID][curIdx];
          }
          else if (compBegin == COMPONENT_Y && compID != COMPONENT_Y && y % (1 << scaleY) == 0 && x % (1 << scaleX) == 0)
          {
            uint32_t posYC = y >> scaleY;
            uint32_t posXC = x >> scaleX;
            picReco.at(posXC, posYC) = cu.curPLT[compID][curIdx];
          }
        }
      }
    }
  }
  for (uint32_t compID = compBegin; compID < (compBegin + numComp); compID++)
  {
    const CompArea &area = cu.blocks[compID];
    PelBuf picReco = cu.cs->getRecoBuf(area);
#if JVET_Z0118_GDR
    cu.cs->updateReconMotIPM(area, picReco);
#else
    cu.cs->picture->getRecoBuf(area).copyFrom(picReco);
#endif
    cu.cs->setDecomp(area);
  }
}

/** Function for deriving reconstructed PU/CU chroma samples with QTree structure
* \param pcRecoYuv pointer to reconstructed sample arrays
* \param pcPredYuv pointer to prediction sample arrays
* \param pcResiYuv pointer to residue sample arrays
* \param chType    texture channel type (luma/chroma)
* \param rTu       reference to transform data
*
\ This function derives reconstructed PU/CU chroma samples with QTree recursive structure
*/

void
DecCu::xIntraRecQT(CodingUnit &cu, const ChannelType chType)
{
  for( auto &currTU : CU::traverseTUs( cu ) )
  {
    if( isLuma( chType ) )
    {
      xIntraRecBlk( currTU, COMPONENT_Y );
    }
    else
    {
      const uint32_t numValidComp = getNumberValidComponents( cu.chromaFormat );

      for( uint32_t compID = COMPONENT_Cb; compID < numValidComp; compID++ )
      {
        xIntraRecBlk( currTU, ComponentID( compID ) );
      }
    }
  }
}

void DecCu::xIntraRecACTQT(CodingUnit &cu)
{
  for (auto &currTU : CU::traverseTUs(cu))
  {
    xIntraRecACTBlk(currTU);
  }
}

#if !REMOVE_PCM
/** Function for filling the PCM buffer of a CU using its reconstructed sample array
* \param pCU   pointer to current CU
* \param depth CU Depth
*/
void DecCu::xFillPCMBuffer(CodingUnit &cu)
{
  for( auto &currTU : CU::traverseTUs( cu ) )
  {
    for (const CompArea &area : currTU.blocks)
    {
      if( !area.valid() ) continue;;

      CPelBuf source      = cu.cs->getRecoBuf(area);
       PelBuf destination = currTU.getPcmbuf(area.compID);

      destination.copyFrom(source);
    }
  }
}
#endif

#include "CommonLib/dtrace_buffer.h"

void DecCu::xReconInter(CodingUnit &cu)
{
  if( cu.geoFlag )
  {
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
    m_pcInterPred->motionCompensationGeo( cu, m_geoMrgCtx 
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
                                        , m_geoTmMrgCtx
#endif
#if JVET_Y0065_GPM_INTRA
                                        , m_pcIntraPred, (cu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag()) ? &m_pcReshape->getFwdLUT() : nullptr
#endif
    );
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
    MergeCtx& m_geoTmMrgCtx0 = m_geoTmMrgCtx[g_geoTmShape[0][g_GeoParams[cu.firstPU->geoSplitDir][0]]];
    MergeCtx& m_geoTmMrgCtx1 = m_geoTmMrgCtx[g_geoTmShape[1][g_GeoParams[cu.firstPU->geoSplitDir][0]]];
#endif
#if JVET_W0097_GPM_MMVD_TM
    PU::spanGeoMMVDMotionInfo
#else
	  PU::spanGeoMotionInfo
#endif
                             ( *cu.firstPU, m_geoMrgCtx
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
							                , m_geoTmMrgCtx0, m_geoTmMrgCtx1
#endif
							                , cu.firstPU->geoSplitDir, cu.firstPU->geoMergeIdx0, cu.firstPU->geoMergeIdx1
#if JVET_W0097_GPM_MMVD_TM
#if TM_MRG
                              , cu.firstPU->geoTmFlag0
#endif
							                , cu.firstPU->geoMMVDFlag0, cu.firstPU->geoMMVDIdx0
#if TM_MRG							 
							                , cu.firstPU->geoTmFlag1
#endif
							                , cu.firstPU->geoMMVDFlag1, cu.firstPU->geoMMVDIdx1
#endif
    );
#else
#if JVET_W0097_GPM_MMVD_TM
#if TM_MRG
#if JVET_Y0065_GPM_INTRA
    m_pcInterPred->motionCompensationGeo( cu, m_geoMrgCtx, m_geoTmMrgCtx0, m_geoTmMrgCtx1, m_pcIntraPred, (cu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag()) ? &m_pcReshape->getFwdLUT() : nullptr );
#else
    m_pcInterPred->motionCompensationGeo(cu, m_geoMrgCtx, m_geoTmMrgCtx0, m_geoTmMrgCtx1);
#endif
    PU::spanGeoMMVDMotionInfo(*cu.firstPU, m_geoMrgCtx, m_geoTmMrgCtx0, m_geoTmMrgCtx1, cu.firstPU->geoSplitDir, cu.firstPU->geoMergeIdx0, cu.firstPU->geoMergeIdx1, cu.firstPU->geoTmFlag0, cu.firstPU->geoMMVDFlag0, cu.firstPU->geoMMVDIdx0, cu.firstPU->geoTmFlag1, cu.firstPU->geoMMVDFlag1, cu.firstPU->geoMMVDIdx1);
#else
#if JVET_Y0065_GPM_INTRA
    m_pcInterPred->motionCompensationGeo( cu, m_geoMrgCtx, m_pcIntraPred, (cu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag()) ? &m_pcReshape->getFwdLUT() : nullptr );
#else
    m_pcInterPred->motionCompensationGeo(cu, m_geoMrgCtx);
#endif
    PU::spanGeoMMVDMotionInfo(*cu.firstPU, m_geoMrgCtx, cu.firstPU->geoSplitDir, cu.firstPU->geoMergeIdx0, cu.firstPU->geoMergeIdx1, cu.firstPU->geoMMVDFlag0, cu.firstPU->geoMMVDIdx0, cu.firstPU->geoMMVDFlag1, cu.firstPU->geoMMVDIdx1);
#endif
#else
#if JVET_Y0065_GPM_INTRA
    m_pcInterPred->motionCompensationGeo( cu, m_geoMrgCtx, m_pcIntraPred, (cu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag()) ? &m_pcReshape->getFwdLUT() : nullptr );
#else
    m_pcInterPred->motionCompensationGeo( cu, m_geoMrgCtx );
#endif
    PU::spanGeoMotionInfo( *cu.firstPU, m_geoMrgCtx, cu.firstPU->geoSplitDir, cu.firstPU->geoMergeIdx0, cu.firstPU->geoMergeIdx1 );
#endif
#endif
  }
  else
  {
#if JVET_X0141_CIIP_TIMD_TM && JVET_W0123_TIMD_FUSION
    if (cu.firstPU->ciipFlag
#if CIIP_PDPC
      && !cu.firstPU->ciipPDPC
#endif
      )
    {
      PredictionUnit *pu = cu.firstPU;
      const CompArea &area = cu.Y();
      if (cu.slice->getSPS()->getUseTimd() && (cu.lwidth() * cu.lheight() <= CIIP_MAX_SIZE))
      {
#if SECONDARY_MPM && ENABLE_DIMD
        IntraPrediction::deriveDimdMode(cu.cs->picture->getRecoBuf(area), area, cu);
#endif
        cu.timdMode = m_pcIntraPred->deriveTimdMode(cu.cs->picture->getRecoBuf(area), area, cu);
        pu->intraDir[0] = MAP131TO67(cu.timdMode);
      }
    }
#endif
    m_pcIntraPred->geneIntrainterPred(cu, m_ciipBuffer);

    // inter prediction
    CHECK(CU::isIBC(cu) && cu.firstPU->ciipFlag, "IBC and Ciip cannot be used together");
    CHECK(CU::isIBC(cu) && cu.affine, "IBC and Affine cannot be used together");
    CHECK(CU::isIBC(cu) && cu.geoFlag, "IBC and geo cannot be used together");
    CHECK(CU::isIBC(cu) && cu.firstPU->mmvdMergeFlag, "IBC and MMVD cannot be used together");
    const bool luma = cu.Y().valid();
    const bool chroma = isChromaEnabled(cu.chromaFormat) && cu.Cb().valid();
    if (luma && (chroma || !isChromaEnabled(cu.chromaFormat)))
    {
      m_pcInterPred->motionCompensation(cu);
    }
    else
    {
      m_pcInterPred->motionCompensation(cu, REF_PIC_LIST_0, luma, chroma);
    }

#if MULTI_PASS_DMVR
    if( cu.firstPU->bdmvrRefine )
    {
      PU::spanMotionInfo( *cu.firstPU, MergeCtx(), m_mvBufBDMVR[0], m_mvBufBDMVR[1], m_pcInterPred->getBdofSubPuMvOffset() );
    }
#endif
  }
  if (cu.Y().valid())
  {
    bool isIbcSmallBlk = CU::isIBC(cu) && (cu.lwidth() * cu.lheight() <= 16);
    CU::saveMotionInHMVP( cu, isIbcSmallBlk );
  }

  if (cu.firstPU->ciipFlag)
  {
    const UnitArea localUnitArea( cu.cs->area.chromaFormat, Area( 0, 0, cu.Y().width, cu.Y().height ) );
#if ENABLE_OBMC && JVET_X0090_CIIP_FIX
    cu.isobmcMC = true;
    m_pcInterPred->subBlockOBMC(*cu.firstPU);
    cu.isobmcMC = false;
#endif
    if( cu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() )
    {
      m_pcIntraPred->geneWeightedPred<true>( COMPONENT_Y, cu.cs->getPredBuf( *cu.firstPU ).Y(), *cu.firstPU, cu.cs->getPredBuf( *cu.firstPU ).Y(), m_ciipBuffer.getBuf( localUnitArea.Y() ), m_pcReshape->getFwdLUT().data() );
    }
    else
    {
      m_pcIntraPred->geneWeightedPred<false>( COMPONENT_Y, cu.cs->getPredBuf( *cu.firstPU ).Y(), *cu.firstPU, cu.cs->getPredBuf( *cu.firstPU ).Y(), m_ciipBuffer.getBuf( localUnitArea.Y() ) );
    }

#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    if( isChromaEnabled( cu.chromaFormat ) )
#else
    if( isChromaEnabled( cu.chromaFormat ) && cu.chromaSize().width > 2 )
#endif
    {
      m_pcIntraPred->geneWeightedPred<false>( COMPONENT_Cb, cu.cs->getPredBuf( *cu.firstPU ).Cb(), *cu.firstPU, cu.cs->getPredBuf( *cu.firstPU ).Cb(), m_ciipBuffer.getBuf( localUnitArea.Cb() ) );
      m_pcIntraPred->geneWeightedPred<false>( COMPONENT_Cr, cu.cs->getPredBuf( *cu.firstPU ).Cr(), *cu.firstPU, cu.cs->getPredBuf( *cu.firstPU ).Cr(), m_ciipBuffer.getBuf( localUnitArea.Cr() ) );
    }
  }
#if ENABLE_OBMC
  cu.isobmcMC = true;
#if JVET_X0090_CIIP_FIX
#if JVET_Y0065_GPM_INTRA
  if (!cu.firstPU->ciipFlag && !cu.firstPU->gpmIntraFlag)
#else
  if (!cu.firstPU->ciipFlag)
#endif
  {
    m_pcInterPred->subBlockOBMC(*cu.firstPU);
  }
#else
  m_pcInterPred->subBlockOBMC(*cu.firstPU);
#endif
  cu.isobmcMC = false;
#endif
  DTRACE    ( g_trace_ctx, D_TMP, "pred " );
  DTRACE_CRC( g_trace_ctx, D_TMP, *cu.cs, cu.cs->getPredBuf( cu ), &cu.Y() );

  // inter recon
  xDecodeInterTexture(cu);

  // clip for only non-zero cbf case
  CodingStructure &cs = *cu.cs;

#if !SIGN_PREDICTION
  if (cu.rootCbf)
  {
#if !JVET_S0234_ACT_CRS_FIX
    if (cu.colorTransform)
    {
      cs.getResiBuf(cu).colorSpaceConvert(cs.getResiBuf(cu), false, cu.cs->slice->clpRng(COMPONENT_Y));
    }
#endif
#if REUSE_CU_RESULTS
    const CompArea &area = cu.blocks[COMPONENT_Y];
    CompArea    tmpArea(COMPONENT_Y, area.chromaFormat, Position(0, 0), area.size());
    PelBuf tmpPred;
#endif
    if (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
    {
#if REUSE_CU_RESULTS
      if (cs.pcv->isEncoder)
      {
        tmpPred = m_tmpStorageLCU->getBuf(tmpArea);
        tmpPred.copyFrom(cs.getPredBuf(cu).get(COMPONENT_Y));
      }
#endif
      if (!cu.firstPU->ciipFlag && !CU::isIBC(cu))
        cs.getPredBuf(cu).get(COMPONENT_Y).rspSignal(m_pcReshape->getFwdLUT());
    }
#if KEEP_PRED_AND_RESI_SIGNALS
    cs.getRecoBuf( cu ).reconstruct( cs.getPredBuf( cu ), cs.getResiBuf( cu ), cs.slice->clpRngs() );
#else
    cs.getResiBuf( cu ).reconstruct( cs.getPredBuf( cu ), cs.getResiBuf( cu ), cs.slice->clpRngs() );
    cs.getRecoBuf( cu ).copyFrom   (                      cs.getResiBuf( cu ) );
#endif
    if (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
    {
#if REUSE_CU_RESULTS
      if (cs.pcv->isEncoder)
      {
        cs.getPredBuf(cu).get(COMPONENT_Y).copyFrom(tmpPred);
      }
#endif
    }
  }
  else
  {
    cs.getRecoBuf(cu).copyClip(cs.getPredBuf(cu), cs.slice->clpRngs());
    if (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !cu.firstPU->ciipFlag && !CU::isIBC(cu))
    {
      cs.getRecoBuf(cu).get(COMPONENT_Y).rspSignal(m_pcReshape->getFwdLUT());
    }
  }
#endif

  DTRACE    ( g_trace_ctx, D_TMP, "reco " );
  DTRACE_CRC( g_trace_ctx, D_TMP, *cu.cs, cu.cs->getRecoBuf( cu ), &cu.Y() );

  cs.setDecomp(cu);
}

void DecCu::xDecodeInterTU( TransformUnit & currTU, const ComponentID compID )
{
  if( !currTU.blocks[compID].valid() ) return;

  const CompArea &area = currTU.blocks[compID];

  CodingStructure& cs = *currTU.cs;

  //===== inverse transform =====
  PelBuf resiBuf  = cs.getResiBuf(area);

  QpParam cQP(currTU, compID);

#if SIGN_PREDICTION
  bool bJccrWithCr = currTU.jointCbCr && !(currTU.jointCbCr >> 1);
  bool bIsJccr     = currTU.jointCbCr && isChroma(compID);
  ComponentID signPredCompID = bIsJccr ? (bJccrWithCr ? COMPONENT_Cr : COMPONENT_Cb): compID;

  bool reshapeChroma = currTU.cs->slice->getPicHeader()->getLmcsEnabledFlag() && isChroma(signPredCompID) && (TU::getCbf(currTU, signPredCompID) || currTU.jointCbCr) && currTU.cs->slice->getPicHeader()->getLmcsChromaResidualScaleFlag() && currTU.blocks[signPredCompID].width * currTU.blocks[signPredCompID].height > 4;
#if JVET_Y0065_GPM_INTRA
  if (currTU.cs->slice->getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma(compID) && !currTU.cu->firstPU->ciipFlag && !currTU.cu->firstPU->gpmIntraFlag && !CU::isIBC(*currTU.cu))
#else
  if (currTU.cs->slice->getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma(compID) && !currTU.cu->firstPU->ciipFlag && !CU::isIBC(*currTU.cu))
#endif
  {
#if JVET_Z0118_GDR
    cs.updateReconMotIPM(currTU.blocks[COMPONENT_Y], cs.getPredBuf(currTU.blocks[COMPONENT_Y]));
#else
    cs.picture->getRecoBuf(currTU.blocks[COMPONENT_Y]).copyFrom(cs.getPredBuf(currTU.blocks[COMPONENT_Y]));
#endif
    cs.getPredBuf(currTU.blocks[COMPONENT_Y]).rspSignal(m_pcReshape->getFwdLUT());
  }
  m_pcTrQuant->predCoeffSigns(currTU, compID, reshapeChroma);
#if JVET_Y0065_GPM_INTRA
  if (currTU.cs->slice->getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma(compID) && !currTU.cu->firstPU->ciipFlag && !currTU.cu->firstPU->gpmIntraFlag && !CU::isIBC(*currTU.cu))
#else
  if (currTU.cs->slice->getPicHeader()->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma(compID) && !currTU.cu->firstPU->ciipFlag && !CU::isIBC(*currTU.cu))
#endif
  {
    cs.getPredBuf(currTU.blocks[COMPONENT_Y]).copyFrom(cs.picture->getRecoBuf(currTU.blocks[COMPONENT_Y]));
  }
#endif

  if( currTU.jointCbCr && isChroma(compID) )
  {
    if( compID == COMPONENT_Cb )
    {
      PelBuf resiCr = cs.getResiBuf( currTU.blocks[ COMPONENT_Cr ] );
      if( currTU.jointCbCr >> 1 )
      {
        m_pcTrQuant->invTransformNxN( currTU, COMPONENT_Cb, resiBuf, cQP );
      }
      else
      {
        QpParam qpCr(currTU, COMPONENT_Cr);
        m_pcTrQuant->invTransformNxN( currTU, COMPONENT_Cr, resiCr, qpCr );
      }
      m_pcTrQuant->invTransformICT( currTU, resiBuf, resiCr );
    }
  }
  else
  if( TU::getCbf( currTU, compID ) )
  {
    m_pcTrQuant->invTransformNxN( currTU, compID, resiBuf, cQP );
  }
  else
  {
    resiBuf.fill( 0 );
  }

  //===== reconstruction =====
  const Slice           &slice = *cs.slice;
#if JVET_S0234_ACT_CRS_FIX
  if (!currTU.cu->colorTransform && slice.getLmcsEnabledFlag() && isChroma(compID) && (TU::getCbf(currTU, compID) || currTU.jointCbCr)
#else
  if (slice.getLmcsEnabledFlag() && isChroma(compID) && (TU::getCbf(currTU, compID) || currTU.jointCbCr)
#endif
   && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && currTU.blocks[compID].width * currTU.blocks[compID].height > 4)
  {
    resiBuf.scaleSignal(currTU.getChromaAdj(), 0, currTU.cu->cs->slice->clpRng(compID));
  }

#if SIGN_PREDICTION
  int firstComponent = compID;
  int lastComponent  = compID;

  if( currTU.jointCbCr && isChroma(compID) )
  {
    if( compID == COMPONENT_Cb )
    {
      firstComponent = MAX_INT;
    }
    else
    {
      firstComponent -= 1;
    }
  }
  for( auto i = firstComponent; i <= lastComponent; ++i)
  {
    ComponentID currCompID = (ComponentID) i;
    PelBuf compResiBuf = cs.getResiBuf(currTU.blocks[currCompID]);

#if JVET_Y0065_GPM_INTRA
    if( cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( currCompID ) && !currTU.cu->firstPU->ciipFlag && !currTU.cu->firstPU->gpmIntraFlag && !CU::isIBC( *currTU.cu ) )
#else
    if( cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && isLuma( currCompID ) && !currTU.cu->firstPU->ciipFlag && !CU::isIBC( *currTU.cu ) )
#endif
    {
      PelBuf picRecoBuff = currTU.cs->picture->getRecoBuf( currTU.blocks[currCompID] );
#if JVET_Z0118_GDR
      currTU.cs->rspSignalPicture(currTU.blocks[currCompID], m_pcReshape->getFwdLUT());
#else
      picRecoBuff.rspSignal( cs.getPredBuf( currTU.blocks[currCompID] ), m_pcReshape->getFwdLUT() );
#endif
      currTU.cs->getRecoBuf( currTU.blocks[currCompID] ).reconstruct( picRecoBuff, compResiBuf, currTU.cu->cs->slice->clpRng( currCompID ) );
    }
    else
    {
      currTU.cs->getRecoBuf( currTU.blocks[currCompID] ).reconstruct( cs.getPredBuf( currTU.blocks[currCompID] ), compResiBuf, currTU.cu->cs->slice->clpRng( currCompID ) );
    }
  }
#endif
}

void DecCu::xDecodeInterTexture(CodingUnit &cu)
{
  if( !cu.rootCbf )
  {
#if SIGN_PREDICTION
    CodingStructure &cs = *cu.cs;
    cs.getRecoBuf(cu).copyClip(cs.getPredBuf(cu), cs.slice->clpRngs());
#if JVET_Y0065_GPM_INTRA
    if (cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !cu.firstPU->ciipFlag && !cu.firstPU->gpmIntraFlag && !CU::isIBC(cu))
#else
    if (cs.picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && !cu.firstPU->ciipFlag && !CU::isIBC(cu))
#endif
    {
      cs.getRecoBuf(cu).get(COMPONENT_Y).rspSignal(m_pcReshape->getFwdLUT());
    }
#endif
    return;
  }

  const uint32_t uiNumVaildComp = getNumberValidComponents(cu.chromaFormat);

#if JVET_S0234_ACT_CRS_FIX
  if (cu.colorTransform)
  {
    CodingStructure  &cs = *cu.cs;
    const Slice &slice = *cs.slice;
    for (auto& currTU : CU::traverseTUs(cu))
    {
      for (uint32_t ch = 0; ch < uiNumVaildComp; ch++)
      {
        const ComponentID compID = ComponentID(ch);
        if (slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && (compID == COMPONENT_Y))
        {
          const CompArea &areaY = currTU.blocks[COMPONENT_Y];
          int adj = m_pcReshape->calculateChromaAdjVpduNei(currTU, areaY);
          currTU.setChromaAdj(adj);
        }
        xDecodeInterTU(currTU, compID);
      }

      cs.getResiBuf(currTU).colorSpaceConvert(cs.getResiBuf(currTU), false, cu.cs->slice->clpRng(COMPONENT_Y));
      if (slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && currTU.blocks[COMPONENT_Cb].width * currTU.blocks[COMPONENT_Cb].height > 4)
      {
        cs.getResiBuf(currTU.blocks[COMPONENT_Cb]).scaleSignal(currTU.getChromaAdj(), 0, currTU.cu->cs->slice->clpRng(COMPONENT_Cb));
        cs.getResiBuf(currTU.blocks[COMPONENT_Cr]).scaleSignal(currTU.getChromaAdj(), 0, currTU.cu->cs->slice->clpRng(COMPONENT_Cr));
      }
    }
  }
  else
  {
#endif
#if SIGN_PREDICTION
  for ( auto& currTU : CU::traverseTUs( cu ) )
  {
    for(uint32_t ch = 0; ch < uiNumVaildComp; ch++)
    {
      const ComponentID compID = ComponentID(ch);
#else
  for (uint32_t ch = 0; ch < uiNumVaildComp; ch++)
  {
    const ComponentID compID = ComponentID(ch);

    for( auto& currTU : CU::traverseTUs( cu ) )
    {
#endif
      CodingStructure  &cs = *cu.cs;
      const Slice &slice = *cs.slice;
      if (slice.getLmcsEnabledFlag() && slice.getPicHeader()->getLmcsChromaResidualScaleFlag() && (compID == COMPONENT_Y) && (currTU.cbf[COMPONENT_Cb] || currTU.cbf[COMPONENT_Cr]))
      {
#if LMCS_CHROMA_CALC_CU
        const CompArea &areaY = cu.blocks[COMPONENT_Y];
#else
        const CompArea &areaY = currTU.blocks[COMPONENT_Y];
#endif
        int adj = m_pcReshape->calculateChromaAdjVpduNei(currTU, areaY);
        currTU.setChromaAdj(adj);
    }
      xDecodeInterTU( currTU, compID );
#if SIGN_PREDICTION
      if(cs.pcv->isEncoder && cs.sps->getNumPredSigns() > 0)
      {
        PelBuf picRecoBuff = currTU.cs->picture->getRecoBuf( currTU.blocks[compID] );
        picRecoBuff.copyFrom(currTU.cs->getRecoBuf( currTU.blocks[compID] ));
      }
#endif
    }
  }
#if JVET_S0234_ACT_CRS_FIX
  }
#endif
}

void DecCu::xDeriveCUMV( CodingUnit &cu )
{
  for( auto &pu : CU::traversePUs( cu ) )
  {
    MergeCtx mrgCtx;

#if RExt__DECODER_DEBUG_TOOL_STATISTICS
    if( pu.cu->affine )
    {
      CodingStatistics::IncrementStatisticTool( CodingStatisticsClassType{ STATS__TOOL_AFF, pu.Y().width, pu.Y().height } );
    }
#endif


    if( pu.mergeFlag )
    {
#if MULTI_HYP_PRED
#if REUSE_CU_RESULTS
      CHECK(!cu.cs->pcv->isEncoder && cu.skip && !pu.addHypData.empty(), "Dec: additional hypotheseis signalled in SKIP mode");
      CHECK(cu.skip && pu.addHypData.size() != pu.numMergedAddHyps, "additional hypotheseis signalled in SKIP mode");
      // save signalled add hyp data, to put it at the end after merging
      MultiHypVec addHypData = { std::begin(pu.addHypData) + pu.numMergedAddHyps, std::end(pu.addHypData) };
#else
      CHECK(cu.skip && !pu.addHypData.empty(), "additional hypotheseis signalled in SKIP mode");
      // save signalled add hyp data, to put it at the end after merging
      auto addHypData = std::move(pu.addHypData);
#endif
      pu.addHypData.clear();
      pu.numMergedAddHyps = 0;
#endif
      if (pu.mmvdMergeFlag || pu.cu->mmvdSkip)
      {
        CHECK(pu.ciipFlag == true, "invalid Ciip");
        if (pu.cs->sps->getSbTMVPEnabledFlag())
        {
          Size bufSize = g_miScaling.scale(pu.lumaSize());
          mrgCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
        }

        int   fPosBaseIdx = pu.mmvdMergeIdx / MMVD_MAX_REFINE_NUM;
        PU::getInterMergeCandidates(pu, mrgCtx, 1, fPosBaseIdx + 1);
        PU::getInterMMVDMergeCandidates(pu, mrgCtx,
          pu.mmvdMergeIdx
        );
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
        uint32_t mmvdLUT[MMVD_ADD_NUM];
        uint8_t mmvdIdx = pu.mmvdMergeIdx;
        m_pcInterPred->sortInterMergeMMVDCandidates(pu, mrgCtx, mmvdLUT, mmvdIdx);
        mrgCtx.setMmvdMergeCandiInfo(pu, mmvdIdx, mmvdLUT[mmvdIdx]);
#else
        mrgCtx.setMmvdMergeCandiInfo(pu, pu.mmvdMergeIdx);
#endif

        PU::spanMotionInfo(pu, mrgCtx);
      }
      else
      {
      {
        if( pu.cu->geoFlag )
        {
          PU::getGeoMergeCandidates( pu, m_geoMrgCtx );
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
          if (pu.geoTmFlag0)
          {
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
            MergeCtx& m_geoTmMrgCtx0 = m_geoTmMrgCtx[GEO_TM_SHAPE_AL];
#endif
            m_geoTmMrgCtx0.numValidMergeCand = m_geoMrgCtx.numValidMergeCand;
            m_geoTmMrgCtx0.BcwIdx[pu.geoMergeIdx0] = BCW_DEFAULT;
            m_geoTmMrgCtx0.useAltHpelIf[pu.geoMergeIdx0] = false;
#if INTER_LIC
            m_geoTmMrgCtx0.LICFlags[pu.geoMergeIdx0] = false;
#endif
#if MULTI_HYP_PRED
            m_geoTmMrgCtx0.addHypNeighbours[pu.geoMergeIdx0] = m_geoMrgCtx.addHypNeighbours[pu.geoMergeIdx0];
#endif
            m_geoTmMrgCtx0.interDirNeighbours[pu.geoMergeIdx0] = m_geoMrgCtx.interDirNeighbours[pu.geoMergeIdx0];
            m_geoTmMrgCtx0.mvFieldNeighbours[(pu.geoMergeIdx0 << 1)].mv = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx0 << 1)].mv;
            m_geoTmMrgCtx0.mvFieldNeighbours[(pu.geoMergeIdx0 << 1) + 1].mv = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx0 << 1) + 1].mv;
            m_geoTmMrgCtx0.mvFieldNeighbours[(pu.geoMergeIdx0 << 1)].refIdx = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx0 << 1)].refIdx;
            m_geoTmMrgCtx0.mvFieldNeighbours[(pu.geoMergeIdx0 << 1) + 1].refIdx = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx0 << 1) + 1].refIdx;
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
            memcpy(&m_geoTmMrgCtx[GEO_TM_SHAPE_A], &m_geoTmMrgCtx0, sizeof(m_geoTmMrgCtx0));
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
            for (uint8_t tmType = GEO_TM_SHAPE_AL; tmType < GEO_NUM_TM_MV_CAND; ++tmType)
            {
              if (tmType == GEO_TM_SHAPE_L || (!pu.cs->slice->getSPS()->getUseAltGPMSplitModeCode() && tmType != g_geoTmShape[0][g_GeoParams[pu.geoSplitDir][0]]))
              {
                continue;
              }
              m_geoTmMrgCtx[tmType].setMergeInfo(pu, pu.geoMergeIdx0);
              pu.geoTmType = tmType;
              m_pcInterPred->deriveTMMv(pu);
              m_geoTmMrgCtx[tmType].mvFieldNeighbours[(pu.geoMergeIdx0 << 1)    ].mv.set(pu.mv[0].getHor(), pu.mv[0].getVer());
              m_geoTmMrgCtx[tmType].mvFieldNeighbours[(pu.geoMergeIdx0 << 1) + 1].mv.set(pu.mv[1].getHor(), pu.mv[1].getVer());
            }
#else
            m_geoTmMrgCtx0.setMergeInfo(pu, pu.geoMergeIdx0);
            pu.geoTmType = g_geoTmShape[0][g_GeoParams[pu.geoSplitDir][0]];
            m_pcInterPred->deriveTMMv(pu);
            m_geoTmMrgCtx0.mvFieldNeighbours[(pu.geoMergeIdx0 << 1)].mv.set(pu.mv[0].getHor(), pu.mv[0].getVer());
            m_geoTmMrgCtx0.mvFieldNeighbours[(pu.geoMergeIdx0 << 1) + 1].mv.set(pu.mv[1].getHor(), pu.mv[1].getVer());
#endif
          }
          if (pu.geoTmFlag1)
          {
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
            MergeCtx& m_geoTmMrgCtx1 = m_geoTmMrgCtx[GEO_TM_SHAPE_AL];
#endif
            m_geoTmMrgCtx1.numValidMergeCand = m_geoMrgCtx.numValidMergeCand;
            m_geoTmMrgCtx1.BcwIdx[pu.geoMergeIdx1] = BCW_DEFAULT;
            m_geoTmMrgCtx1.useAltHpelIf[pu.geoMergeIdx1] = false;
#if INTER_LIC
            m_geoTmMrgCtx1.LICFlags[pu.geoMergeIdx1] = false;
#endif
#if MULTI_HYP_PRED
            m_geoTmMrgCtx1.addHypNeighbours[pu.geoMergeIdx1] = m_geoMrgCtx.addHypNeighbours[pu.geoMergeIdx1];
#endif
            m_geoTmMrgCtx1.interDirNeighbours[pu.geoMergeIdx1] = m_geoMrgCtx.interDirNeighbours[pu.geoMergeIdx1];
            m_geoTmMrgCtx1.mvFieldNeighbours[(pu.geoMergeIdx1 << 1)].mv = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx1 << 1)].mv;
            m_geoTmMrgCtx1.mvFieldNeighbours[(pu.geoMergeIdx1 << 1) + 1].mv = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx1 << 1) + 1].mv;
            m_geoTmMrgCtx1.mvFieldNeighbours[(pu.geoMergeIdx1 << 1)].refIdx = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx1 << 1)].refIdx;
            m_geoTmMrgCtx1.mvFieldNeighbours[(pu.geoMergeIdx1 << 1) + 1].refIdx = m_geoMrgCtx.mvFieldNeighbours[(pu.geoMergeIdx1 << 1) + 1].refIdx;
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
            memcpy(&m_geoTmMrgCtx[GEO_TM_SHAPE_L], &m_geoTmMrgCtx1, sizeof(m_geoTmMrgCtx1));
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
            for (uint8_t tmType = GEO_TM_SHAPE_AL; tmType < GEO_NUM_TM_MV_CAND; ++tmType)
            {
              if (tmType == GEO_TM_SHAPE_A || (!pu.cs->slice->getSPS()->getUseAltGPMSplitModeCode() && tmType != g_geoTmShape[1][g_GeoParams[pu.geoSplitDir][0]]))
              {
                continue;
              }
              m_geoTmMrgCtx[tmType].setMergeInfo(pu, pu.geoMergeIdx1);
              pu.geoTmType = tmType;
              m_pcInterPred->deriveTMMv(pu);
              m_geoTmMrgCtx[tmType].mvFieldNeighbours[(pu.geoMergeIdx1 << 1)    ].mv.set(pu.mv[0].getHor(), pu.mv[0].getVer());
              m_geoTmMrgCtx[tmType].mvFieldNeighbours[(pu.geoMergeIdx1 << 1) + 1].mv.set(pu.mv[1].getHor(), pu.mv[1].getVer());
            }
#else
            m_geoTmMrgCtx1.setMergeInfo(pu, pu.geoMergeIdx1);
            pu.geoTmType = g_geoTmShape[1][g_GeoParams[pu.geoSplitDir][0]];
            m_pcInterPred->deriveTMMv(pu);
            m_geoTmMrgCtx1.mvFieldNeighbours[(pu.geoMergeIdx1 << 1)].mv.set(pu.mv[0].getHor(), pu.mv[0].getVer());
            m_geoTmMrgCtx1.mvFieldNeighbours[(pu.geoMergeIdx1 << 1) + 1].mv.set(pu.mv[1].getHor(), pu.mv[1].getVer());
#endif
          }
#endif
        }
        else
        {
        if( pu.cu->affine )
        {
          AffineMergeCtx affineMergeCtx;
          if (pu.cs->sps->getSbTMVPEnabledFlag())
          {
            Size bufSize = g_miScaling.scale( pu.lumaSize() );
            mrgCtx.subPuMvpMiBuf = MotionBuf( m_SubPuMiBuf, bufSize );
            affineMergeCtx.mrgCtx = &mrgCtx;
          }
#if AFFINE_MMVD
#if JVET_W0090_ARMC_TM
          int affMrgIdx = pu.cs->sps->getUseAML() && (((pu.mergeIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE < pu.cs->sps->getMaxNumAffineMergeCand()) || (pu.mergeIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE) == 0) ? pu.mergeIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE * ADAPTIVE_AFFINE_SUB_GROUP_SIZE + ADAPTIVE_AFFINE_SUB_GROUP_SIZE - 1 : pu.mergeIdx;
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
#if JVET_Z0139_NA_AFF
          PU::getAffineMergeCand(pu, affineMergeCtx, (pu.afMmvdFlag ? AF_MMVD_BASE_NUM : affMrgIdx), pu.afMmvdFlag, pu.mergeIdx == 0);
#else
          PU::getAffineMergeCand(pu, affineMergeCtx, (pu.afMmvdFlag ? AF_MMVD_BASE_NUM : affMrgIdx), pu.afMmvdFlag);
#endif
#else
#if JVET_Z0139_NA_AFF
          PU::getAffineMergeCand(pu, affineMergeCtx, (pu.afMmvdFlag ? pu.afMmvdBaseIdx : affMrgIdx), pu.afMmvdFlag, pu.mergeIdx == 0);
#else
          PU::getAffineMergeCand(pu, affineMergeCtx, (pu.afMmvdFlag ? pu.afMmvdBaseIdx : affMrgIdx), pu.afMmvdFlag);
#endif
#endif

#else
          PU::getAffineMergeCand(pu, affineMergeCtx, (pu.afMmvdFlag ? pu.afMmvdBaseIdx : pu.mergeIdx), pu.afMmvdFlag);
#endif

          if (pu.afMmvdFlag)
          {
            pu.mergeIdx = PU::getMergeIdxFromAfMmvdBaseIdx(affineMergeCtx, pu.afMmvdBaseIdx);
            CHECK(pu.mergeIdx >= pu.cu->slice->getPicHeader()->getMaxNumAffineMergeCand(), "Affine MMVD mode doesn't have a valid base candidate!");
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
            uint32_t affMmvdLUT[AF_MMVD_NUM];
            uint32_t tempVal = pu.afMmvdMergeIdx;
            m_pcInterPred->sortAffineMergeCandidates(pu, affineMergeCtx, affMmvdLUT, tempVal );
            pu.afMmvdMergeIdx = tempVal;
            int uiMergeCandTemp = affMmvdLUT[pu.afMmvdMergeIdx];
            int baseIdx = (int)uiMergeCandTemp / AF_MMVD_MAX_REFINE_NUM;
            int stepIdx = (int)uiMergeCandTemp - baseIdx * AF_MMVD_MAX_REFINE_NUM;
            int dirIdx  = stepIdx % AF_MMVD_OFFSET_DIR;
            stepIdx = stepIdx / AF_MMVD_OFFSET_DIR;
            pu.afMmvdDir      = (uint8_t)dirIdx;
            pu.afMmvdStep     = (uint8_t)stepIdx;
#endif
            PU::getAfMmvdMvf(pu, affineMergeCtx, affineMergeCtx.mvFieldNeighbours + (pu.mergeIdx << 1), pu.mergeIdx, pu.afMmvdStep, pu.afMmvdDir);
          }
#if JVET_W0090_ARMC_TM
          else
          {
            if (pu.cs->sps->getUseAML())
#if JVET_Z0139_NA_AFF
            if (affineMergeCtx.numValidMergeCand > 1)
#endif
            {
              m_pcInterPred->adjustAffineMergeCandidates(pu, affineMergeCtx, pu.mergeIdx);
            }
          }
#endif
#else
#if JVET_W0090_ARMC_TM
          int affMrgIdx = pu.cs->sps->getUseAML() && (((pu.mergeIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE + 1)*ADAPTIVE_AFFINE_SUB_GROUP_SIZE < pu.cs->sps->getMaxNumAffineMergeCand()) || (pu.mergeIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE) == 0) ? pu.mergeIdx / ADAPTIVE_AFFINE_SUB_GROUP_SIZE * ADAPTIVE_AFFINE_SUB_GROUP_SIZE + ADAPTIVE_AFFINE_SUB_GROUP_SIZE - 1 : pu.mergeIdx;
#if !JVET_Z0139_NA_AFF
          PU::getAffineMergeCand(pu, affineMergeCtx, affMrgIdx);
#else
          PU::getAffineMergeCand(pu, affineMergeCtx, affMrgIdx, pu.mergeIdx == 0);
          if (affineMergeCtx.numValidMergeCand > 1)
#endif
          if (pu.cs->sps->getUseAML())
          {
            m_pcInterPred->adjustAffineMergeCandidates(pu, affineMergeCtx, pu.mergeIdx);
          }
#else
          PU::getAffineMergeCand( pu, affineMergeCtx, pu.mergeIdx );
#endif
#endif
          pu.interDir = affineMergeCtx.interDirNeighbours[pu.mergeIdx];
          pu.cu->affineType = affineMergeCtx.affineType[pu.mergeIdx];
          pu.cu->BcwIdx = affineMergeCtx.BcwIdx[pu.mergeIdx];
#if INTER_LIC
          pu.cu->LICFlag = affineMergeCtx.LICFlags[pu.mergeIdx];
#endif
          pu.mergeType = affineMergeCtx.mergeType[pu.mergeIdx];
          if ( pu.mergeType == MRG_TYPE_SUBPU_ATMVP )
          {
            pu.refIdx[0] = affineMergeCtx.mvFieldNeighbours[(pu.mergeIdx << 1) + 0][0].refIdx;
            pu.refIdx[1] = affineMergeCtx.mvFieldNeighbours[(pu.mergeIdx << 1) + 1][0].refIdx;
          }
          else
          {
          for( int i = 0; i < 2; ++i )
          {
            if( pu.cs->slice->getNumRefIdx( RefPicList( i ) ) > 0 )
            {
              MvField* mvField = affineMergeCtx.mvFieldNeighbours[(pu.mergeIdx << 1) + i];
              pu.mvpIdx[i] = 0;
              pu.mvpNum[i] = 0;
              pu.mvd[i]    = Mv();
              pu.refIdx[i] = mvField[0].refIdx;
              pu.mvAffi[i][0] = mvField[0].mv;
              pu.mvAffi[i][1] = mvField[1].mv;
              pu.mvAffi[i][2] = mvField[2].mv;
            }
          }
        }
          PU::spanMotionInfo( pu, mrgCtx );
        }
#if JVET_X0141_CIIP_TIMD_TM && TM_MRG
        else if (pu.ciipFlag && pu.tmMergeFlag)
        {
          int storeMrgIdx = pu.mergeIdx;
          pu.tmMergeFlag = false;
          PU::getInterMergeCandidates(pu, mrgCtx, 0, CIIP_TM_MRG_MAX_NUM_CANDS - 1);
          mrgCtx.numValidMergeCand = int(pu.cs->sps->getMaxNumCiipTMMergeCand());
          pu.tmMergeFlag = true;
          for (uint32_t uiMergeCand = 0; uiMergeCand < CIIP_TM_MRG_MAX_NUM_CANDS; uiMergeCand++)
          {
            mrgCtx.setMergeInfo(pu, uiMergeCand);
            m_pcInterPred->deriveTMMv(pu);
            // Store refined motion back to ciipTmMrgCtx
            mrgCtx.interDirNeighbours[uiMergeCand] = pu.interDir;
            mrgCtx.BcwIdx[uiMergeCand] = pu.cu->BcwIdx;  // Bcw may change, because bi may be reduced to uni by deriveTMMv(pu)
            mrgCtx.mvFieldNeighbours[2 * uiMergeCand].setMvField(pu.mv[0], pu.refIdx[0]);
            mrgCtx.mvFieldNeighbours[2 * uiMergeCand + 1].setMvField(pu.mv[1], pu.refIdx[1]);
            if (pu.interDir == 1)
            {
              mrgCtx.mvFieldNeighbours[2 * uiMergeCand + 1].setMvField(Mv(), NOT_VALID);
            }
            if (pu.interDir == 2)
            {
              mrgCtx.mvFieldNeighbours[2 * uiMergeCand].setMvField(Mv(), NOT_VALID);
            }
          }
#if JVET_W0090_ARMC_TM
          if (pu.cs->sps->getUseAML())
          {
             m_pcInterPred->adjustInterMergeCandidates(pu, mrgCtx, CIIP_TM_MRG_MAX_NUM_CANDS - 1);
          }
#endif

          mrgCtx.setMergeInfo(pu, storeMrgIdx);
          pu.bdmvrRefine = false;
          PU::spanMotionInfo(pu, mrgCtx);
        }
#endif
        else
        {
          if (CU::isIBC(*pu.cu))
          {
#if JVET_Y0058_IBC_LIST_MODIFY && JVET_W0090_ARMC_TM
            if (pu.cs->sps->getUseAML())
            {
#if JVET_Z0075_IBC_HMVP_ENLARGE
              uint8_t mrgCandIdx = pu.mergeIdx;
              PU::getIBCMergeCandidates(pu, mrgCtx);
              m_pcInterPred->adjustIBCMergeCandidates(pu, mrgCtx, 0, IBC_MRG_MAX_NUM_CANDS_MEM);
              pu.mergeIdx = mrgCandIdx;
#else
              PU::getIBCMergeCandidates(pu, mrgCtx, (((pu.mergeIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE + 1)*ADAPTIVE_IBC_SUB_GROUP_SIZE < pu.cs->sps->getMaxNumIBCMergeCand()) || (pu.mergeIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE) == 0) ? pu.mergeIdx / ADAPTIVE_IBC_SUB_GROUP_SIZE * ADAPTIVE_IBC_SUB_GROUP_SIZE + ADAPTIVE_IBC_SUB_GROUP_SIZE - 1 : pu.mergeIdx);
              m_pcInterPred->adjustIBCMergeCandidates(pu, mrgCtx, pu.mergeIdx);
#endif
            }
            else
#endif
              PU::getIBCMergeCandidates(pu, mrgCtx, pu.mergeIdx);
          }
          else
#if JVET_X0049_ADAPT_DMVR
            if (pu.bmMergeFlag)
            {
              uint8_t mergeIdx = pu.bmDir == 2 ? pu.mergeIdx - BM_MRG_MAX_NUM_CANDS : pu.mergeIdx;
#if JVET_W0090_ARMC_TM
              if (pu.cs->sps->getUseAML())
              {
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
                int nWidth = pu.lumaSize().width;
                int nHeight = pu.lumaSize().height;
                bool tplAvail = m_pcInterPred->xAMLGetCurBlkTemplate(pu, nWidth, nHeight);

                uint8_t bmDir = pu.bmDir;
                MergeCtx tmvpMergeCandCtx;
                PU::getTmvpBMCand(pu, tmvpMergeCandCtx);
                pu.bmDir = 0;
                if (tplAvail)
                {
                  m_pcInterPred->adjustMergeCandidatesInOneCandidateGroup(pu, tmvpMergeCandCtx, 1, pu.mergeIdx);
                }
                else
                {
                  tmvpMergeCandCtx.numValidMergeCand = std::min(1, tmvpMergeCandCtx.numValidMergeCand);
                }
                MergeCtx namvpMergeCandCtx;
                PU::getNonAdjacentBMCand(pu, namvpMergeCandCtx);
                pu.bmDir = 0;
                if (tplAvail)
                {
                  m_pcInterPred->adjustMergeCandidatesInOneCandidateGroup(pu, namvpMergeCandCtx, 3, pu.mergeIdx);
                }
                else
                {
                  namvpMergeCandCtx.numValidMergeCand = std::min(3, namvpMergeCandCtx.numValidMergeCand);
                }
                PU::getInterBMCandidates(pu, mrgCtx, -1, &tmvpMergeCandCtx, &namvpMergeCandCtx);
#else
                PU::getInterBMCandidates(pu, mrgCtx, pu.cs->sps->getUseAML() && (((mergeIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE < pu.cs->sps->getMaxNumBMMergeCand()) || (mergeIdx / ADAPTIVE_SUB_GROUP_SIZE) == 0) ? mergeIdx / ADAPTIVE_SUB_GROUP_SIZE * ADAPTIVE_SUB_GROUP_SIZE + ADAPTIVE_SUB_GROUP_SIZE - 1 : mergeIdx);
                uint8_t bmDir = pu.bmDir;
#endif
                pu.bmDir = 0;
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
                if (tplAvail)
                {
#if JVET_Z0102_NO_ARMC_FOR_ZERO_CAND 
                  m_pcInterPred->adjustMergeCandidates(pu, mrgCtx, pu.cs->sps->getMaxNumBMMergeCand());
#else
                  m_pcInterPred->adjustMergeCandidatesInOneCandidateGroup(pu, mrgCtx, mergeIdx + 1, mergeIdx);
#endif
                }
#else
                m_pcInterPred->adjustInterMergeCandidates(pu, mrgCtx, mergeIdx);
#endif
                pu.bmDir = bmDir;
              }
              else
#endif
              PU::getInterBMCandidates(pu, mrgCtx, mergeIdx);
            }
            else
#endif
#if JVET_W0090_ARMC_TM
            if (pu.cs->sps->getUseAML())
            {
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
              int nWidth = pu.lumaSize().width;
              int nHeight = pu.lumaSize().height;
              bool tplAvail = m_pcInterPred->xAMLGetCurBlkTemplate(pu, nWidth, nHeight);

              MergeCtx tmvpMergeCandCtx;
              PU::getTmvpMergeCand(pu, tmvpMergeCandCtx);
              if (tplAvail)
              {
                m_pcInterPred->adjustMergeCandidatesInOneCandidateGroup(pu, tmvpMergeCandCtx, 1, pu.mergeIdx);
              }
              else
              {
                tmvpMergeCandCtx.numValidMergeCand = std::min(1, tmvpMergeCandCtx.numValidMergeCand);
              }
              MergeCtx namvpMergeCandCtx;
              PU::getNonAdjacentMergeCand(pu, namvpMergeCandCtx);
              if (tplAvail)
              {
                m_pcInterPred->adjustMergeCandidatesInOneCandidateGroup(pu, namvpMergeCandCtx, 9, pu.mergeIdx);
              }
              else
              {
                namvpMergeCandCtx.numValidMergeCand = std::min(9, namvpMergeCandCtx.numValidMergeCand);
              }

              PU::getInterMergeCandidates(pu, mrgCtx, 0, -1, &tmvpMergeCandCtx, &namvpMergeCandCtx);

              if (tplAvail)
              {
#if JVET_Z0102_NO_ARMC_FOR_ZERO_CAND
                m_pcInterPred->adjustMergeCandidates(pu, mrgCtx, pu.tmMergeFlag ? pu.cs->sps->getMaxNumTMMergeCand() : pu.cs->sps->getMaxNumMergeCand());
#else
                m_pcInterPred->adjustMergeCandidatesInOneCandidateGroup(pu, mrgCtx, pu.mergeIdx + 1, pu.mergeIdx);
#endif
              }
#else
              PU::getInterMergeCandidates(pu, mrgCtx, 0, pu.cs->sps->getUseAML() && (((pu.mergeIdx / ADAPTIVE_SUB_GROUP_SIZE + 1)*ADAPTIVE_SUB_GROUP_SIZE < pu.cs->sps->getMaxNumMergeCand()) || (pu.mergeIdx / ADAPTIVE_SUB_GROUP_SIZE) == 0) ? pu.mergeIdx / ADAPTIVE_SUB_GROUP_SIZE * ADAPTIVE_SUB_GROUP_SIZE + ADAPTIVE_SUB_GROUP_SIZE - 1 : pu.mergeIdx);
              m_pcInterPred->adjustInterMergeCandidates(pu, mrgCtx, pu.mergeIdx);
#endif
            }
            else
            {
              PU::getInterMergeCandidates(pu, mrgCtx, 0, pu.mergeIdx);
            }
#else
            PU::getInterMergeCandidates(pu, mrgCtx, 0, pu.mergeIdx);
#endif
          mrgCtx.setMergeInfo( pu, pu.mergeIdx );
#if TM_MRG && !MULTI_PASS_DMVR
          if (pu.tmMergeFlag)
          {
            m_pcInterPred->deriveTMMv(pu);
#if JVET_Z0084_IBC_TM && TM_MRG
            if (CU::isIBC(*pu.cu))
            {
              pu.bv = pu.mv[0];
              pu.bv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
            }
#endif
          }
#endif
#if MULTI_PASS_DMVR
          if (PU::checkBDMVRCondition(pu))
          {
            m_pcInterPred->setBdmvrSubPuMvBuf(m_mvBufBDMVR[0], m_mvBufBDMVR[1]);
            pu.bdmvrRefine = true;

            CHECK(mrgCtx.numValidMergeCand <= 0, "this is not possible");

#if JVET_X0049_ADAPT_DMVR
#if TM_MRG
            if (!pu.tmMergeFlag && !pu.bmMergeFlag && mrgCtx.xCheckSimilarMotion(pu.mergeIdx, PU::getBDMVRMvdThreshold(pu)))
#else
            if (!pu.bmMergeFlag && mrgCtx.xCheckSimilarMotion(pu.mergeIdx, PU::getBDMVRMvdThreshold(pu)))
#endif
#else
#if TM_MRG
            if( !pu.tmMergeFlag && mrgCtx.xCheckSimilarMotion( pu.mergeIdx, PU::getBDMVRMvdThreshold( pu ) ) )
#else
            if( mrgCtx.xCheckSimilarMotion( pu.mergeIdx, PU::getBDMVRMvdThreshold( pu ) ) )
#endif
#endif
            {
              // span motion to subPU
              for (int subPuIdx = 0; subPuIdx < MAX_NUM_SUBCU_DMVR; subPuIdx++)
              {
                m_mvBufBDMVR[0][subPuIdx] = pu.mv[0];
                m_mvBufBDMVR[1][subPuIdx] = pu.mv[1];
              }
            }
            else
            {
              pu.bdmvrRefine = m_pcInterPred->processBDMVR( pu );
            }
          }
#if TM_MRG
          else
          {
            if (pu.tmMergeFlag)
            {
              m_pcInterPred->deriveTMMv(pu);
#if JVET_Z0084_IBC_TM && TM_MRG
              if (CU::isIBC(*pu.cu))
              {
                pu.bv = pu.mv[0];
                pu.bv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
              }
#endif
            }
          }
#endif
#endif
#if MULTI_HYP_PRED
          CHECK(pu.Y().area() <= MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE && !pu.addHypData.empty(), "There are add hyps in small block")
#endif
#if MULTI_PASS_DMVR
          if( !pu.bdmvrRefine )
          {
            PU::spanMotionInfo( pu, mrgCtx );
          }
#else
          PU::spanMotionInfo( pu, mrgCtx );
#endif
        }
        }
      }
      }
#if MULTI_HYP_PRED
      // put saved additional hypotheseis to the end
#if JVET_Z0118_GDR
      int n = (int) addHypData.capacity();
      int s = (int) pu.addHypData.size();
      int e = (int) (addHypData.end() - addHypData.begin());

      if ((s + e) < n)
      {
        pu.addHypData.insert(pu.addHypData.end(), std::make_move_iterator(addHypData.begin()), std::make_move_iterator(addHypData.end()));
      }
#else
      pu.addHypData.insert(pu.addHypData.end(), std::make_move_iterator(addHypData.begin()), std::make_move_iterator(addHypData.end()));
#endif
#endif
    }
    else
    {
#if JVET_Z0054_BLK_REF_PIC_REORDER
#if REUSE_CU_RESULTS
      if (!cu.cs->pcv->isEncoder)
#endif
      {
        m_pcInterPred->clearTplAmvpBuffer();
        if (pu.cu->affine)
        {
          for (uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
          {
            if (pu.interDir & (1 << uiRefListIdx))
            {
              pu.mvdAffi[uiRefListIdx][0].changeAffinePrecAmvr2Internal(pu.cu->imv);
              pu.mvdAffi[uiRefListIdx][1].changeAffinePrecAmvr2Internal(pu.cu->imv);
              if (cu.affineType == AFFINEMODEL_6PARAM)
              {
                pu.mvdAffi[uiRefListIdx][2].changeAffinePrecAmvr2Internal(pu.cu->imv);
              }
            }
          }
        }
        else if (CU::isIBC(*pu.cu) && pu.interDir == 1)
        {
          pu.mvd[REF_PIC_LIST_0].changeIbcPrecAmvr2Internal(pu.cu->imv);
        }
        else
        {
          for (uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
          {
            if (pu.interDir & (1 << uiRefListIdx))
            {
              pu.mvd[uiRefListIdx].changeTransPrecAmvr2Internal(pu.cu->imv);
            }
          }
        }
        if (PU::useRefPairList(pu))
        {
          m_pcInterPred->setBiRefIdx(pu);
        }
        if (PU::useRefCombList(pu))
        {
          m_pcInterPred->setUniRefListAndIdx(pu);
        }
      }
#endif

#if !JVET_Z0054_BLK_REF_PIC_REORDER
#if REUSE_CU_RESULTS
      if ( cu.imv && !pu.cu->affine && !cu.cs->pcv->isEncoder )
#else
        if (cu.imv && !pu.cu->affine)
#endif
        {
          PU::applyImv(pu, mrgCtx, m_pcInterPred);
        }
        else
#endif
      {
#if JVET_X0083_BM_AMVP_MERGE_MODE
#if REUSE_CU_RESULTS
        if (!cu.cs->pcv->isEncoder && (pu.amvpMergeModeFlag[REF_PIC_LIST_0] || pu.amvpMergeModeFlag[REF_PIC_LIST_1]))
#else
        if (pu.amvpMergeModeFlag[REF_PIC_LIST_0] || pu.amvpMergeModeFlag[REF_PIC_LIST_1])
#endif
        {
          CHECK(pu.interDir != 3, "this is not possible");
          const RefPicList refListMerge = pu.amvpMergeModeFlag[REF_PIC_LIST_0] ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
          const RefPicList refListAmvp = RefPicList(1 - refListMerge);
          int orgRefIdxAMVP = pu.refIdx[refListAmvp];
          int orgInterDir = pu.interDir;
          int orgMvpIdxL0 = pu.mvpIdx[REF_PIC_LIST_0];
          int orgMvpIdxL1 = pu.mvpIdx[REF_PIC_LIST_1];
          Mv orgMvd0 = pu.mvd[0];
          Mv orgMvd1 = pu.mvd[1];
          // this part is to derive the merge info
          m_pcInterPred->getAmvpMergeModeMergeList(pu, m_mvFieldAmListDec, orgRefIdxAMVP);
          // if there was set PU merge info, restore the AMVP information
          pu.mvpIdx[REF_PIC_LIST_0] = orgMvpIdxL0;
          pu.mvpIdx[REF_PIC_LIST_1] = orgMvpIdxL1;
          pu.interDir = orgInterDir;
          pu.mergeFlag = false;
          pu.mvd[0] = orgMvd0;
          pu.mvd[1] = orgMvd1;
          pu.refIdx[refListAmvp] = orgRefIdxAMVP;
        }
#endif

        if( pu.cu->affine )
        {
          for ( uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            RefPicList eRefList = RefPicList( uiRefListIdx );
            if ( pu.cs->slice->getNumRefIdx( eRefList ) > 0 && ( pu.interDir & ( 1 << uiRefListIdx ) ) )
            {
              AffineAMVPInfo affineAMVPInfo;
              PU::fillAffineMvpCand( pu, eRefList, pu.refIdx[eRefList], affineAMVPInfo );

              const unsigned mvpIdx = pu.mvpIdx[eRefList];

              pu.mvpNum[eRefList] = affineAMVPInfo.numCand;

              //    Mv mv[3];
              CHECK( pu.refIdx[eRefList] < 0, "Unexpected negative refIdx." );
#if JVET_Z0054_BLK_REF_PIC_REORDER
              if (!pu.cs->sps->getUseARL())
              {
#else
              if (!cu.cs->pcv->isEncoder)
              {
                pu.mvdAffi[eRefList][0].changeAffinePrecAmvr2Internal(pu.cu->imv);
                pu.mvdAffi[eRefList][1].changeAffinePrecAmvr2Internal(pu.cu->imv);
                if (cu.affineType == AFFINEMODEL_6PARAM)
                {
                  pu.mvdAffi[eRefList][2].changeAffinePrecAmvr2Internal(pu.cu->imv);
                }
              }
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED 
              if (pu.isMvsdApplicable())
              {
                Mv absMvd[3];
                absMvd[0] = Mv(pu.mvdAffi[eRefList][0].getAbsMv());
                absMvd[1] = Mv(pu.mvdAffi[eRefList][1].getAbsMv());
                absMvd[2] = (cu.affineType == AFFINEMODEL_6PARAM) ? Mv(pu.mvdAffi[eRefList][2].getAbsMv()) : Mv(0, 0);
                if ((absMvd[0] != Mv(0, 0) || absMvd[1] != Mv(0, 0) || absMvd[2] != Mv(0, 0)) && pu.isMvsdApplicable())
                {
                  std::vector<Mv> cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3;
#if JVET_Z0054_BLK_REF_PIC_REORDER
                  m_pcInterPred->deriveMvdSignAffine(affineAMVPInfo.mvCandLT[mvpIdx], affineAMVPInfo.mvCandRT[mvpIdx], affineAMVPInfo.mvCandLB[mvpIdx],
                    absMvd, pu, eRefList, pu.refIdx[eRefList], cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3);
#else
                  m_pcInterPred->deriveMvdSignAffine(affineAMVPInfo.mvCandLT[mvpIdx], affineAMVPInfo.mvCandRT[mvpIdx], affineAMVPInfo.mvCandLB[mvpIdx],
                    absMvd[0], absMvd[1], absMvd[2], pu, eRefList, pu.refIdx[eRefList], cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3);
#endif
                  CHECK(pu.mvsdIdx[eRefList] >= cMvdDerivedVec.size(), "");
                  m_pcInterPred->deriveMVDFromMVSDIdxAffine(pu, eRefList, cMvdDerivedVec, cMvdDerivedVec2, cMvdDerivedVec3);
                }
              }
#endif
#if JVET_Z0054_BLK_REF_PIC_REORDER
              }
#endif
              Mv mvLT = affineAMVPInfo.mvCandLT[mvpIdx] + pu.mvdAffi[eRefList][0];
              Mv mvRT = affineAMVPInfo.mvCandRT[mvpIdx] + pu.mvdAffi[eRefList][1];
              mvRT += pu.mvdAffi[eRefList][0];

              Mv mvLB;
              if ( cu.affineType == AFFINEMODEL_6PARAM )
              {
                mvLB = affineAMVPInfo.mvCandLB[mvpIdx] + pu.mvdAffi[eRefList][2];
                mvLB += pu.mvdAffi[eRefList][0];
              }

              pu.mvAffi[eRefList][0] = mvLT;
              pu.mvAffi[eRefList][1] = mvRT;
              pu.mvAffi[eRefList][2] = mvLB;
            }
          }
        }
        else if (CU::isIBC(*pu.cu) && pu.interDir == 1)
        {
          AMVPInfo amvpInfo;
#if JVET_Z0084_IBC_TM && TM_AMVP
          PU::fillIBCMvpCand(pu, amvpInfo, m_pcInterPred);
#else
          PU::fillIBCMvpCand(pu, amvpInfo);
#endif
          pu.mvpNum[REF_PIC_LIST_0] = amvpInfo.numCand;
          Mv mvd = pu.mvd[REF_PIC_LIST_0];
#if !JVET_Z0054_BLK_REF_PIC_REORDER
#if REUSE_CU_RESULTS
          if (!cu.cs->pcv->isEncoder)
#endif
          {
            mvd.changeIbcPrecAmvr2Internal(pu.cu->imv);
          }
#endif
          if (pu.cs->sps->getMaxNumIBCMergeCand() == 1)
          {
            CHECK( pu.mvpIdx[REF_PIC_LIST_0], "mvpIdx for IBC mode should be 0" );
          }
          pu.mv[REF_PIC_LIST_0] = amvpInfo.mvCand[pu.mvpIdx[REF_PIC_LIST_0]] + mvd;
          pu.mv[REF_PIC_LIST_0].mvCliptoStorageBitDepth();
#if JVET_Z0160_IBC_ZERO_PADDING
          pu.bv = pu.mv[0];
          pu.bv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
#endif
        }
        else
        {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
          Mv cMvpL0;
#endif
          for ( uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            RefPicList eRefList = RefPicList( uiRefListIdx );
            if ((pu.cs->slice->getNumRefIdx(eRefList) > 0 || (eRefList == REF_PIC_LIST_0 && CU::isIBC(*pu.cu))) && (pu.interDir & (1 << uiRefListIdx)))
            {
              AMVPInfo amvpInfo;
#if JVET_X0083_BM_AMVP_MERGE_MODE
              if (pu.amvpMergeModeFlag[eRefList] == true)
              {
#if REUSE_CU_RESULTS
                if (!cu.cs->pcv->isEncoder)
                {
#endif
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
                const int mvFieldMergeIdx = pu.refIdx[1 - eRefList] * AMVP_MAX_NUM_CANDS_MEM + pu.mvpIdx[1 - eRefList];
#else
                const int mvFieldMergeIdx = pu.refIdx[1 - eRefList] * AMVP_MAX_NUM_CANDS + pu.mvpIdx[1 - eRefList];
#endif
                pu.mv[eRefList] = m_mvFieldAmListDec[mvFieldMergeIdx].mv;
                pu.refIdx[eRefList] = m_mvFieldAmListDec[mvFieldMergeIdx].refIdx;
#if REUSE_CU_RESULTS
                }
#endif
              }
              else
              {

                if (pu.amvpMergeModeFlag[1 - eRefList] == true)
                {
#if TM_AMVP
#if JVET_Y0128_NON_CTC
                  amvpInfo.numCand = PU::checkTmEnableCondition(pu.cs->sps, pu.cs->pps, pu.cu->slice->getRefPic(eRefList, pu.refIdx[eRefList])) ? 1 : 2;
#else
                  amvpInfo.numCand = 1;
#endif
#else
                  amvpInfo.numCand = AMVP_MAX_NUM_CANDS;
#endif
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
                  amvpInfo.numCand += 1;
#endif
#if REUSE_CU_RESULTS
                  if (cu.cs->pcv->isEncoder)
                  {
                    amvpInfo.mvCand[pu.mvpIdx[eRefList]] = pu.mv[eRefList] - pu.mvd[eRefList];
                  }
                  else
#endif
                  {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
                    const int mvFieldAmvpIdx = MAX_NUM_AMVP_CANDS_MAX_REF + pu.refIdx[eRefList] * AMVP_MAX_NUM_CANDS_MEM + pu.mvpIdx[eRefList];
#else
                    const int mvFieldAmvpIdx = MAX_NUM_AMVP_CANDS_MAX_REF + pu.refIdx[eRefList] * AMVP_MAX_NUM_CANDS + pu.mvpIdx[eRefList];
#endif
                    amvpInfo.mvCand[pu.mvpIdx[eRefList]] = m_mvFieldAmListDec[mvFieldAmvpIdx].mv;
                  }
                }
                else
#endif
              PU::fillMvpCand(pu, eRefList, pu.refIdx[eRefList], amvpInfo
#if TM_AMVP
                            , m_pcInterPred
#endif
              );
              pu.mvpNum [eRefList] = amvpInfo.numCand;
#if JVET_Z0054_BLK_REF_PIC_REORDER
              if(!PU::useRefCombList(pu) && !PU::useRefPairList(pu))
              {
#else
              if (!cu.cs->pcv->isEncoder)
              {
                pu.mvd[eRefList].changeTransPrecAmvr2Internal(pu.cu->imv);
              }
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
              if (pu.isMvsdApplicable() && pu.mvd[eRefList].isMvsdApplicable())
              {
                if (pu.cu->smvdMode)
                {
                  if (eRefList == REF_PIC_LIST_0)
                  {
                    cMvpL0 = amvpInfo.mvCand[pu.mvpIdx[eRefList]];
                  }
                  else
                  {
                    std::vector<Mv> cMvdDerivedVec;
                    Mv cMvdKnownAtDecoder = Mv(pu.mvd[REF_PIC_LIST_0].getAbsHor(), pu.mvd[REF_PIC_LIST_0].getAbsVer());
                    m_pcInterPred->deriveMvdSignSMVD(cMvpL0, amvpInfo.mvCand[pu.mvpIdx[1]], cMvdKnownAtDecoder, pu, cMvdDerivedVec); //pass the absolute Mvd value as MVd may be negative while CU reuse at encoder.
                    CHECK(pu.mvsdIdx[REF_PIC_LIST_0] >= cMvdDerivedVec.size(), "");
                    int mvsdIdx = pu.mvsdIdx[REF_PIC_LIST_0];
                    Mv cMvd = m_pcInterPred->deriveMVDFromMVSDIdxTrans(mvsdIdx, cMvdDerivedVec);
                    CHECK(cMvd == Mv(0, 0), " zero MVD for SMVD!");
                    pu.mvd[REF_PIC_LIST_0] = cMvd;
                    pu.mv[REF_PIC_LIST_0] = cMvpL0 + pu.mvd[REF_PIC_LIST_0];
                    pu.mv[REF_PIC_LIST_0].mvCliptoStorageBitDepth();
                    pu.mvd[REF_PIC_LIST_1].set(-pu.mvd[REF_PIC_LIST_0].hor, -pu.mvd[REF_PIC_LIST_0].ver);
                  }
                }
                else
                {
                  std::vector<Mv> cMvdDerivedVec;
                  Mv cMvdKnownAtDecoder = Mv(pu.mvd[eRefList].getAbsHor(), pu.mvd[eRefList].getAbsVer());
                  m_pcInterPred->deriveMvdSign(amvpInfo.mvCand[pu.mvpIdx[eRefList]], cMvdKnownAtDecoder, pu, eRefList, pu.refIdx[eRefList], cMvdDerivedVec); //pass the absolute Mvd value as MVd may be negative while CU reuse at encoder.
                  CHECK(pu.mvsdIdx[eRefList] >= cMvdDerivedVec.size(), "");
                  int mvsdIdx = pu.mvsdIdx[eRefList];
                  Mv cMvd = m_pcInterPred->deriveMVDFromMVSDIdxTrans(mvsdIdx, cMvdDerivedVec);
                  CHECK(cMvd == Mv(0, 0), " zero MVD!");
                  pu.mvd[eRefList] = cMvd;
                }
              }
#endif
#if JVET_Z0054_BLK_REF_PIC_REORDER
              }
#endif
              pu.mv[eRefList] = amvpInfo.mvCand[pu.mvpIdx[eRefList]] + pu.mvd[eRefList];
              pu.mv[eRefList].mvCliptoStorageBitDepth();
#if JVET_X0083_BM_AMVP_MERGE_MODE
              }
#endif
            }
          }
#if JVET_X0083_BM_AMVP_MERGE_MODE
          if ((pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1]) && PU::checkBDMVRCondition(pu))
          {
            m_pcInterPred->setBdmvrSubPuMvBuf(m_mvBufBDMVR[0], m_mvBufBDMVR[1]);
            pu.bdmvrRefine = true;
            // span motion to subPU
            for (int subPuIdx = 0; subPuIdx < MAX_NUM_SUBCU_DMVR; subPuIdx++)
            {
              m_mvBufBDMVR[0][subPuIdx] = pu.mv[0];
              m_mvBufBDMVR[1][subPuIdx] = pu.mv[1];
            }
          }
#endif
        }
#if JVET_X0083_BM_AMVP_MERGE_MODE
        if (!pu.bdmvrRefine)
#endif
        PU::spanMotionInfo( pu, mrgCtx );
      }
    }
    if( !cu.geoFlag )
    {
      if( g_mctsDecCheckEnabled && !MCTSHelper::checkMvBufferForMCTSConstraint( pu, true ) )
      {
        printf( "DECODER: pu motion vector across tile boundaries (%d,%d,%d,%d)\n", pu.lx(), pu.ly(), pu.lwidth(), pu.lheight() );
      }
    }
    if (CU::isIBC(cu))
    {
      const int cuPelX = pu.Y().x;
      const int cuPelY = pu.Y().y;
      int roiWidth = pu.lwidth();
      int roiHeight = pu.lheight();
      const unsigned int  lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
      int xPred = pu.mv[0].getHor() >> MV_FRACTIONAL_BITS_INTERNAL;
      int yPred = pu.mv[0].getVer() >> MV_FRACTIONAL_BITS_INTERNAL;
#if JVET_Z0118_GDR
      if (cu.cs->slice->getSPS()->getGDREnabledFlag())
      {
        if (cu.cs->isClean(cu)) 
        {
          CHECK(!m_pcInterPred->isLumaBvValid(lcuWidth, cuPelX, cuPelY, roiWidth, roiHeight, xPred, yPred), "invalid block vector for IBC detected.");
        }
      }
      else 
      {
        CHECK(!m_pcInterPred->isLumaBvValid(lcuWidth, cuPelX, cuPelY, roiWidth, roiHeight, xPred, yPred), "invalid block vector for IBC detected.");
      }
#else
      CHECK(!m_pcInterPred->isLumaBvValid(lcuWidth, cuPelX, cuPelY, roiWidth, roiHeight, xPred, yPred), "invalid block vector for IBC detected.");
#endif
    }
#if MULTI_HYP_PRED
    {
      bool derivedMrgList = false;
#if REUSE_CU_RESULTS
      if (!cu.cs->pcv->isEncoder)
#endif
      for (int i = pu.numMergedAddHyps; i < int(pu.addHypData.size()); ++i)
      {
        auto &mhData = pu.addHypData[i];
#if INTER_LIC
        mhData.LICFlag = pu.cu->LICFlag;
#endif
        mhData.imv = pu.cu->imv;
        if (mhData.isMrg)
        {
          if (!derivedMrgList)
          {
            PU::getGeoMergeCandidates(pu, m_geoMrgCtx);
            derivedMrgList = true;
          }
          int refList = m_geoMrgCtx.interDirNeighbours[mhData.mrgIdx] - 1; CHECK(refList != 0 && refList != 1, "");
          mhData.refIdx = m_geoMrgCtx.mvFieldNeighbours[(mhData.mrgIdx << 1) + refList].refIdx;
          mhData.mv = m_geoMrgCtx.mvFieldNeighbours[(mhData.mrgIdx << 1) + refList].mv;
          mhData.refList = refList;
          mhData.imv = m_geoMrgCtx.useAltHpelIf[mhData.mrgIdx] ? IMV_HPEL : 0;
          continue;
        }
        if (pu.cu->affine)
        {
          mhData.mvd.changeAffinePrecAmvr2Internal(pu.cu->imv);
        }
        else
        {
          mhData.mvd.changeTransPrecAmvr2Internal(pu.cu->imv);
        }
        const Mv mvp = PU::getMultiHypMVP(pu, mhData);
        mhData.mv = mvp + mhData.mvd;
        mhData.mv.mvCliptoStorageBitDepth();
      }
    }
#endif
  }
}
//! \}
