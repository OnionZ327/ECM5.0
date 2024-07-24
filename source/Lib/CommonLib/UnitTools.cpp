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

/** \file     UnitTool.cpp
 *  \brief    defines operations for basic units
 */

#include "UnitTools.h"

#include "dtrace_next.h"

#include "Unit.h"
#include "Slice.h"
#include "Picture.h"

#include <utility>
#include <algorithm>

// CS tools


uint64_t CS::getEstBits(const CodingStructure &cs)
{
  return cs.fracBits >> SCALE_BITS;
}

bool CS::isDualITree( const CodingStructure &cs )
{
  return cs.slice->isIntra() && !cs.pcv->ISingleTree;
}

UnitArea CS::getArea( const CodingStructure &cs, const UnitArea &area, const ChannelType chType )
{
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  return isDualITree(cs) ? area.singleChan(chType) : area;
#else
  return isDualITree( cs ) || cs.treeType != TREE_D ? area.singleChan( chType ) : area;
#endif
}
#if !MULTI_PASS_DMVR
void CS::setRefinedMotionField(CodingStructure &cs)
{
  for (CodingUnit *cu : cs.cus)
  {
    for (auto &pu : CU::traversePUs(*cu))
    {
      PredictionUnit subPu = pu;
      int dx, dy, x, y, num = 0;
      dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
      dx = std::min<int>(pu.lumaSize().width, DMVR_SUBCU_WIDTH);
      Position puPos = pu.lumaPos();
      if (PU::checkDMVRCondition(pu))
      {
        for (y = puPos.y; y < (puPos.y + pu.lumaSize().height); y = y + dy)
        {
          for (x = puPos.x; x < (puPos.x + pu.lumaSize().width); x = x + dx)
          {
            subPu.UnitArea::operator=(UnitArea(pu.chromaFormat, Area(x, y, dx, dy)));
            subPu.mv[0] = pu.mv[0];
            subPu.mv[1] = pu.mv[1];
            subPu.mv[REF_PIC_LIST_0] += pu.mvdL0SubPu[num];
            subPu.mv[REF_PIC_LIST_1] -= pu.mvdL0SubPu[num];
            subPu.mv[REF_PIC_LIST_0].clipToStorageBitDepth();
            subPu.mv[REF_PIC_LIST_1].clipToStorageBitDepth();
            pu.mvdL0SubPu[num].setZero();
            num++;
#if JVET_W0123_TIMD_FUSION
            PU::spanMotionInfo2(subPu);
#else
            PU::spanMotionInfo(subPu);
#endif
          }
        }
      }
    }
  }
}
#endif
// CU tools

bool CU::getRprScaling( const SPS* sps, const PPS* curPPS, Picture* refPic, int& xScale, int& yScale )
{
  const Window& curScalingWindow = curPPS->getScalingWindow();
  int curPicWidth = curPPS->getPicWidthInLumaSamples()   - SPS::getWinUnitX( sps->getChromaFormatIdc() ) * (curScalingWindow.getWindowLeftOffset() + curScalingWindow.getWindowRightOffset());
  int curPicHeight = curPPS->getPicHeightInLumaSamples() - SPS::getWinUnitY( sps->getChromaFormatIdc() ) * (curScalingWindow.getWindowTopOffset()  + curScalingWindow.getWindowBottomOffset());

  const Window& refScalingWindow = refPic->getScalingWindow();
  int refPicWidth = refPic->getPicWidthInLumaSamples()   - SPS::getWinUnitX( sps->getChromaFormatIdc() ) * (refScalingWindow.getWindowLeftOffset() + refScalingWindow.getWindowRightOffset());
  int refPicHeight = refPic->getPicHeightInLumaSamples() - SPS::getWinUnitY( sps->getChromaFormatIdc() ) * (refScalingWindow.getWindowTopOffset()  + refScalingWindow.getWindowBottomOffset());

  xScale = ( ( refPicWidth << SCALE_RATIO_BITS ) + ( curPicWidth >> 1 ) ) / curPicWidth;
  yScale = ( ( refPicHeight << SCALE_RATIO_BITS ) + ( curPicHeight >> 1 ) ) / curPicHeight;

  int curSeqMaxPicWidthY = sps->getMaxPicWidthInLumaSamples();                  // pic_width_max_in_luma_samples
  int curSeqMaxPicHeightY = sps->getMaxPicHeightInLumaSamples();                // pic_height_max_in_luma_samples
  int curPicWidthY = curPPS->getPicWidthInLumaSamples();                        // pic_width_in_luma_samples
  int curPicHeightY = curPPS->getPicHeightInLumaSamples();                      // pic_height_in_luma_samples
  int max8MinCbSizeY = std::max((int)8, (1<<sps->getLog2MinCodingBlockSize())); // Max(8, MinCbSizeY)

  CHECK((curPicWidth * curSeqMaxPicWidthY) < refPicWidth * (curPicWidthY - max8MinCbSizeY), "(curPicWidth * curSeqMaxPicWidthY) should be greater than or equal to refPicWidth * (curPicWidthY - max8MinCbSizeY))");
  CHECK((curPicHeight * curSeqMaxPicHeightY) < refPicHeight * (curPicHeightY - max8MinCbSizeY), "(curPicHeight * curSeqMaxPicHeightY) should be greater than or equal to refPicHeight * (curPicHeightY - max8MinCbSizeY))");

  CHECK(curPicWidth * 2 < refPicWidth, "curPicWidth * 2 shall be greater than or equal to refPicWidth");
  CHECK(curPicHeight * 2 < refPicHeight, "curPicHeight * 2 shall be greater than or equal to refPicHeight");
  CHECK(curPicWidth > refPicWidth * 8, "curPicWidth shall be less than or equal to refPicWidth * 8");
  CHECK(curPicHeight > refPicHeight * 8, "curPicHeight shall be less than or equal to refPicHeight * 8");

#if JVET_S0048_SCALING_OFFSET
  int subWidthC = SPS::getWinUnitX(sps->getChromaFormatIdc());
  int subHeightC = SPS::getWinUnitY(sps->getChromaFormatIdc());

  CHECK(subWidthC * curScalingWindow.getWindowLeftOffset() < (-curPicWidthY) * 15, "The value of SubWidthC * pps_scaling_win_left_offset shall be greater than or equal to -pps_pic_width_in_luma_samples * 15");
  CHECK(subWidthC * curScalingWindow.getWindowLeftOffset() >= curPicWidthY, "The value of SubWidthC * pps_scaling_win_left_offset shall be less than pic_width_in_luma_samples");
  CHECK(subWidthC * curScalingWindow.getWindowRightOffset() < (-curPicWidthY) * 15, "The value of SubWidthC * pps_scaling_win_right_offset shall be greater than or equal to -pps_pic_width_in_luma_samples * 15");
  CHECK(subWidthC * curScalingWindow.getWindowRightOffset() >= curPicWidthY, "The value of SubWidthC * pps_scaling_win_right_offset shall be less than pic_width_in_luma_samples");

  CHECK(subHeightC * curScalingWindow.getWindowTopOffset() < (-curPicHeightY) * 15, "The value of SubHeightC * pps_scaling_win_top_offset shall be greater than or equal to -pps_pic_height_in_luma_samples * 15");
  CHECK(subHeightC * curScalingWindow.getWindowTopOffset() >= curPicHeightY, "The value of SubHeightC * pps_scaling_win_top_offset shall be less than pps_pic_height_in_luma_samples");
  CHECK(subHeightC * curScalingWindow.getWindowBottomOffset() < (-curPicHeightY) * 15, "The value of SubHeightC *pps_scaling_win_bottom_offset shall be greater than or equal to -pps_pic_height_in_luma_samples * 15");
  CHECK(subHeightC * curScalingWindow.getWindowBottomOffset() >= curPicHeightY, "The value of SubHeightC *pps_scaling_win_bottom_offset shall be less than pps_pic_height_in_luma_samples");

  CHECK(subWidthC * (curScalingWindow.getWindowLeftOffset() + curScalingWindow.getWindowRightOffset()) < (-curPicWidthY) * 15, "The value of SubWidthC * ( pps_scaling_win_left_offset + pps_scaling_win_right_offset ) shall be greater than or equal to -pps_pic_width_in_luma_samples * 15");
  CHECK(subWidthC * (curScalingWindow.getWindowLeftOffset() + curScalingWindow.getWindowRightOffset()) >= curPicWidthY, "The value of SubWidthC * ( pps_scaling_win_left_offset + pps_scaling_win_right_offset ) shall be less than pic_width_in_luma_samples");
  CHECK(subHeightC * (curScalingWindow.getWindowTopOffset() + curScalingWindow.getWindowBottomOffset()) < (-curPicHeightY) * 15, "The value of SubHeightC * ( pps_scaling_win_top_offset + pps_scaling_win_bottom_offset ) shall be greater than or equal to -pps_pic_height_in_luma_samples * 15");
  CHECK(subHeightC * (curScalingWindow.getWindowTopOffset() + curScalingWindow.getWindowBottomOffset()) >= curPicHeightY, "The value of SubHeightC * ( pps_scaling_win_top_offset + pps_scaling_win_bottom_offset ) shall be less than pic_height_in_luma_samples");
#else
  CHECK(SPS::getWinUnitX(sps->getChromaFormatIdc()) * (abs(curScalingWindow.getWindowLeftOffset()) + abs(curScalingWindow.getWindowRightOffset())) > curPPS->getPicWidthInLumaSamples(), "The value of SubWidthC * ( Abs(pps_scaling_win_left_offset) + Abs(pps_scaling_win_right_offset) ) shall be less than pic_width_in_luma_samples");
  CHECK(SPS::getWinUnitY(sps->getChromaFormatIdc()) * (abs(curScalingWindow.getWindowTopOffset()) + abs(curScalingWindow.getWindowBottomOffset())) > curPPS->getPicHeightInLumaSamples(), "The value of SubHeightC * ( Abs(pps_scaling_win_top_offset) + Abs(pps_scaling_win_bottom_offset) ) shall be less than pic_height_in_luma_samples");
#endif

  return refPic->isRefScaled( curPPS );
}

void CU::checkConformanceILRP(Slice *slice)
{
  const int numRefList = slice->isInterB() ? 2 : 1;

#if JVET_S0258_SUBPIC_CONSTRAINTS
  int currentSubPicIdx = NOT_VALID;

  // derive sub-picture index for the current slice
  for( int subPicIdx = 0; subPicIdx < slice->getPic()->cs->sps->getNumSubPics(); subPicIdx++ )
  {
    if( slice->getPic()->cs->pps->getSubPic( subPicIdx ).getSubPicID() == slice->getSliceSubPicId() )
    {
      currentSubPicIdx = subPicIdx;
      break;
    }
  }

  CHECK( currentSubPicIdx == NOT_VALID, "Sub-picture was not found" );

  if( !slice->getPic()->cs->sps->getSubPicTreatedAsPicFlag( currentSubPicIdx ) )
  {
    return;
  }
#endif

  //constraint 1: The picture referred to by each active entry in RefPicList[ 0 ] or RefPicList[ 1 ] has the same subpicture layout as the current picture
  bool isAllRefSameSubpicLayout = true;
  for (int refList = 0; refList < numRefList; refList++) // loop over l0 and l1
  {
    RefPicList  eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);

    for (int refIdx = 0; refIdx < slice->getNumRefIdx(eRefPicList); refIdx++)
    {
#if JVET_S0258_SUBPIC_CONSTRAINTS
      const Picture* refPic = slice->getRefPic( eRefPicList, refIdx );

      if( refPic->subPictures.size() != slice->getPic()->cs->pps->getNumSubPics() )
#else
      const Picture* refPic = slice->getRefPic(eRefPicList, refIdx)->unscaledPic;

      if (refPic->numSubpics != slice->getPic()->cs->pps->getNumSubPics())
#endif
      {
        isAllRefSameSubpicLayout = false;
        refList = numRefList;
        break;
      }
      else
      {
#if JVET_S0258_SUBPIC_CONSTRAINTS
        for( int i = 0; i < refPic->subPictures.size(); i++ )
        {
          const SubPic& refSubPic = refPic->subPictures[i];
          const SubPic& curSubPic = slice->getPic()->cs->pps->getSubPic( i );

          if( refSubPic.getSubPicWidthInCTUs() != curSubPic.getSubPicWidthInCTUs()
            || refSubPic.getSubPicHeightInCTUs() != curSubPic.getSubPicHeightInCTUs()
            || refSubPic.getSubPicCtuTopLeftX() != curSubPic.getSubPicCtuTopLeftX()
            || refSubPic.getSubPicCtuTopLeftY() != curSubPic.getSubPicCtuTopLeftY()
            || ( refPic->layerId != slice->getPic()->layerId && refSubPic.getSubPicID() != curSubPic.getSubPicID() )
            || refSubPic.getTreatedAsPicFlag() != curSubPic.getTreatedAsPicFlag())
#else
        for (int i = 0; i < refPic->numSubpics; i++)
        {
          if (refPic->subpicWidthInCTUs[i] != slice->getPic()->cs->pps->getSubPic(i).getSubPicWidthInCTUs()
            || refPic->subpicHeightInCTUs[i] != slice->getPic()->cs->pps->getSubPic(i).getSubPicHeightInCTUs()
            || refPic->subpicCtuTopLeftX[i] != slice->getPic()->cs->pps->getSubPic(i).getSubPicCtuTopLeftX()
            || refPic->subpicCtuTopLeftY[i] != slice->getPic()->cs->pps->getSubPic(i).getSubPicCtuTopLeftY())
#endif
          {
            isAllRefSameSubpicLayout = false;
            refIdx = slice->getNumRefIdx(eRefPicList);
            refList = numRefList;
            break;
          }
        }

#if JVET_S0258_SUBPIC_CONSTRAINTS
        // A picture with different sub-picture ID of the collocated sub-picture cannot be used as an active reference picture in the same layer
        if( refPic->layerId == slice->getPic()->layerId )
        {
          isAllRefSameSubpicLayout = isAllRefSameSubpicLayout && refPic->subPictures[currentSubPicIdx].getSubPicID() == slice->getSliceSubPicId();
        }
#endif
      }
    }
  }

  //constraint 2: The picture referred to by each active entry in RefPicList[ 0 ] or RefPicList[ 1 ] is an ILRP for which the value of sps_num_subpics_minus1 is equal to 0
  if (!isAllRefSameSubpicLayout)
  {
    for (int refList = 0; refList < numRefList; refList++) // loop over l0 and l1
    {
      RefPicList  eRefPicList = (refList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
      for (int refIdx = 0; refIdx < slice->getNumRefIdx(eRefPicList); refIdx++)
      {
#if JVET_S0258_SUBPIC_CONSTRAINTS
        const Picture* refPic = slice->getRefPic( eRefPicList, refIdx );
        CHECK( refPic->layerId == slice->getPic()->layerId || refPic->subPictures.size() > 1, "The inter-layer reference shall contain a single subpicture or have same subpicture layout with the current picture" );
#else
        const Picture* refPic = slice->getRefPic(eRefPicList, refIdx)->unscaledPic;
        CHECK(!(refPic->layerId != slice->getPic()->layerId && refPic->numSubpics == 1), "The inter-layer reference shall contain a single subpicture or have same subpicture layout with the current picture");
#endif
      }
    }
  }

  return;
}

bool CU::isIntra(const CodingUnit &cu)
{
  return cu.predMode == MODE_INTRA;
}

bool CU::isInter(const CodingUnit &cu)
{
  return cu.predMode == MODE_INTER;
}

bool CU::isIBC(const CodingUnit &cu)
{
  return cu.predMode == MODE_IBC;
}

bool CU::isPLT(const CodingUnit &cu)
{
  return cu.predMode == MODE_PLT;
}

bool CU::isSameSlice(const CodingUnit& cu, const CodingUnit& cu2)
{
  return cu.slice->getIndependentSliceIdx() == cu2.slice->getIndependentSliceIdx();
}

bool CU::isSameTile(const CodingUnit& cu, const CodingUnit& cu2)
{
  return cu.tileIdx == cu2.tileIdx;
}

bool CU::isSameSliceAndTile(const CodingUnit& cu, const CodingUnit& cu2)
{
  return ( cu.slice->getIndependentSliceIdx() == cu2.slice->getIndependentSliceIdx() ) && ( cu.tileIdx == cu2.tileIdx );
}

bool CU::isSameSubPic(const CodingUnit& cu, const CodingUnit& cu2)
{
  return (cu.slice->getPPS()->getSubPicFromCU(cu).getSubPicIdx() == cu2.slice->getPPS()->getSubPicFromCU(cu2).getSubPicIdx()) ;
}

bool CU::isSameCtu(const CodingUnit& cu, const CodingUnit& cu2)
{
  uint32_t ctuSizeBit = floorLog2(cu.cs->sps->getMaxCUWidth());

  Position pos1Ctu(cu.lumaPos().x  >> ctuSizeBit, cu.lumaPos().y  >> ctuSizeBit);
  Position pos2Ctu(cu2.lumaPos().x >> ctuSizeBit, cu2.lumaPos().y >> ctuSizeBit);

  return pos1Ctu.x == pos2Ctu.x && pos1Ctu.y == pos2Ctu.y;
}

bool CU::isLastSubCUOfCtu( const CodingUnit &cu )
{
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const Area cuAreaY = CS::isDualITree(*cu.cs) ? Area(recalcPosition(cu.chromaFormat, cu.chType, CHANNEL_TYPE_LUMA, cu.blocks[cu.chType].pos()), recalcSize(cu.chromaFormat, cu.chType, CHANNEL_TYPE_LUMA, cu.blocks[cu.chType].size())) : (const Area&)cu.Y();
#else
  const Area cuAreaY = cu.isSepTree() ? Area( recalcPosition( cu.chromaFormat, cu.chType, CHANNEL_TYPE_LUMA, cu.blocks[cu.chType].pos() ), recalcSize( cu.chromaFormat, cu.chType, CHANNEL_TYPE_LUMA, cu.blocks[cu.chType].size() ) ) : (const Area&)cu.Y();
#endif

  return ( ( ( ( cuAreaY.x + cuAreaY.width  ) & cu.cs->pcv->maxCUWidthMask  ) == 0 || cuAreaY.x + cuAreaY.width  == cu.cs->pps->getPicWidthInLumaSamples()  ) &&
           ( ( ( cuAreaY.y + cuAreaY.height ) & cu.cs->pcv->maxCUHeightMask ) == 0 || cuAreaY.y + cuAreaY.height == cu.cs->pps->getPicHeightInLumaSamples() ) );
}

uint32_t CU::getCtuAddr( const CodingUnit &cu )
{
  return getCtuAddr( cu.blocks[cu.chType].lumaPos(), *cu.cs->pcv );
}
#if JVET_V0130_INTRA_TMP
Position CU::getCtuXYAddr(const CodingUnit& cu)
{
	return Position((cu.blocks[cu.chType].lumaPos().x >> cu.cs->pcv->maxCUWidthLog2) << cu.cs->pcv->maxCUWidthLog2, (cu.blocks[cu.chType].lumaPos().y >> cu.cs->pcv->maxCUHeightLog2) << cu.cs->pcv->maxCUHeightLog2);
}
#endif
int CU::predictQP( const CodingUnit& cu, const int prevQP )
{
  const CodingStructure &cs = *cu.cs;

  uint32_t  ctuRsAddr       = getCtuAddr( cu );
  uint32_t  ctuXPosInCtus   = ctuRsAddr % cs.pcv->widthInCtus;
  uint32_t  tileColIdx      = cu.slice->getPPS()->ctuToTileCol( ctuXPosInCtus );
  uint32_t  tileXPosInCtus  = cu.slice->getPPS()->getTileColumnBd( tileColIdx );
  if( ctuXPosInCtus == tileXPosInCtus &&
      !( cu.blocks[cu.chType].x & ( cs.pcv->maxCUWidthMask  >> getChannelTypeScaleX( cu.chType, cu.chromaFormat ) ) ) &&
      !( cu.blocks[cu.chType].y & ( cs.pcv->maxCUHeightMask >> getChannelTypeScaleY( cu.chType, cu.chromaFormat ) ) ) &&
      ( cs.getCU( cu.blocks[cu.chType].pos().offset( 0, -1 ), cu.chType) != NULL ) &&
      CU::isSameSliceAndTile( *cs.getCU( cu.blocks[cu.chType].pos().offset( 0, -1 ), cu.chType), cu ) )
  {
    return ( ( cs.getCU( cu.blocks[cu.chType].pos().offset( 0, -1 ), cu.chType ) )->qp );
  }
  else
  {
    const int a = ( cu.blocks[cu.chType].y & ( cs.pcv->maxCUHeightMask >> getChannelTypeScaleY( cu.chType, cu.chromaFormat ) ) ) ? ( cs.getCU(cu.blocks[cu.chType].pos().offset( 0, -1 ), cu.chType))->qp : prevQP;
    const int b = ( cu.blocks[cu.chType].x & ( cs.pcv->maxCUWidthMask  >> getChannelTypeScaleX( cu.chType, cu.chromaFormat ) ) ) ? ( cs.getCU(cu.blocks[cu.chType].pos().offset( -1, 0 ), cu.chType))->qp : prevQP;

    return ( a + b + 1 ) >> 1;
  }
}


uint32_t CU::getNumPUs( const CodingUnit& cu )
{
  uint32_t cnt = 0;
  PredictionUnit *pu = cu.firstPU;

  do
  {
    cnt++;
  } while( ( pu != cu.lastPU ) && ( pu = pu->next ) );

  return cnt;
}

void CU::addPUs( CodingUnit& cu )
{
  cu.cs->addPU( CS::getArea( *cu.cs, cu, cu.chType ), cu.chType );
}

void CU::saveMotionInHMVP( const CodingUnit& cu, const bool isToBeDone )
{
  const PredictionUnit& pu = *cu.firstPU;
#if JVET_Z0118_GDR
  bool isClean = pu.cs->isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

#if JVET_Z0139_HIST_AFF
  if (cu.affine)
  {
    if (pu.mergeType != MRG_TYPE_DEFAULT_N)
    {
      return;
    }

    AffineMotionInfo addMi[2];
    int addRefIdx[2];

    pu.getAffineMotionInfo(addMi, addRefIdx);
    
#if JVET_Z0118_GDR    
    cu.cs->addAffMiToLut((isClean) ? cu.cs->motionLut.lutAff1 : cu.cs->motionLut.lutAff0, addMi, addRefIdx);
#else
    cu.cs->addAffMiToLut(cu.cs->motionLut.lutAff, addMi, addRefIdx);
#endif

    AffineInheritInfo addAffInherit;
    addAffInherit.basePos = cu.lumaPos();
    addAffInherit.baseMV[0] = MvField(pu.mvAffi[0][0], pu.refIdx[0]);
    addAffInherit.baseMV[1] = MvField(pu.mvAffi[1][0], pu.refIdx[1]);
    addAffInherit.oneSetAffineParametersPattern0 = addMi[0].oneSetAffineParametersPattern;
    addAffInherit.oneSetAffineParametersPattern1 = addMi[1].oneSetAffineParametersPattern;
    if (addAffInherit.oneSetAffineParametersPattern0 == 0)
    {
      addAffInherit.baseMV[0].refIdx = -1;
    }
    if (addAffInherit.oneSetAffineParametersPattern1 == 0)
    {
      addAffInherit.baseMV[1].refIdx = -1;
    }
    if (addAffInherit.baseMV[0].refIdx != -1 || addAffInherit.baseMV[1].refIdx != -1)
    {
#if JVET_Z0118_GDR
      cu.cs->addAffInheritToLut((isClean) ? cu.cs->motionLut.lutAffInherit1 : cu.cs->motionLut.lutAffInherit0, addAffInherit);
#else
      cu.cs->addAffInheritToLut(cu.cs->motionLut.lutAffInherit, addAffInherit);
#endif
    }
    return;
  }
#endif

  if (!cu.geoFlag && !cu.affine && !isToBeDone)
  {
    MotionInfo mi = pu.getMotionInfo();
#if MULTI_HYP_PRED
    mi.addHypData = pu.addHypData;
#endif

    mi.BcwIdx = (mi.interDir == 3) ? cu.BcwIdx : BCW_DEFAULT;

    const unsigned log2ParallelMergeLevel = (pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2);
    const unsigned xBr = pu.cu->Y().width + pu.cu->Y().x;
    const unsigned yBr = pu.cu->Y().height + pu.cu->Y().y;
    bool enableHmvp = ((xBr >> log2ParallelMergeLevel) > (pu.cu->Y().x >> log2ParallelMergeLevel)) && ((yBr >> log2ParallelMergeLevel) > (pu.cu->Y().y >> log2ParallelMergeLevel));
    bool enableInsertion = CU::isIBC(cu) || enableHmvp;
    if (enableInsertion)
    {
#if JVET_Z0075_IBC_HMVP_ENLARGE
#if JVET_Z0118_GDR      
      if (isClean)
      {
        if (CU::isIBC(cu))
        {
          cu.cs->addMiToLutIBC(cu.cs->motionLut.lutIbc1, mi);
        }
        else
        {
          cu.cs->addMiToLut(cu.cs->motionLut.lut1, mi);
        }
      }

      if (CU::isIBC(cu))
      {
        cu.cs->addMiToLutIBC(cu.cs->motionLut.lutIbc0, mi);
      }
      else
      {
        cu.cs->addMiToLut(cu.cs->motionLut.lut0, mi);
      }
#else 
      if (CU::isIBC(cu))
      {
        cu.cs->addMiToLutIBC(cu.cs->motionLut.lutIbc, mi);
      }
      else
      {
        cu.cs->addMiToLut(cu.cs->motionLut.lut, mi);
      }
#endif
#else
#if JVET_Z0118_GDR      
      if (isClean)
      {
        cu.cs->addMiToLut(CU::isIBC(cu) ? cu.cs->motionLut.lutIbc1 : cu.cs->motionLut.lut1, mi);
      }
      cu.cs->addMiToLut(CU::isIBC(cu) ? cu.cs->motionLut.lutIbc0 : cu.cs->motionLut.lut0, mi);      
#else
      cu.cs->addMiToLut(CU::isIBC(cu) ? cu.cs->motionLut.lutIbc : cu.cs->motionLut.lut, mi);
#endif
#endif
    }
  }
}

PartSplit CU::getSplitAtDepth( const CodingUnit& cu, const unsigned depth )
{
  if (depth >= cu.depth)
  {
    return CU_DONT_SPLIT;
  }

  const PartSplit cuSplitType = PartSplit( ( cu.splitSeries >> ( depth * SPLIT_DMULT ) ) & SPLIT_MASK );

  if (cuSplitType == CU_QUAD_SPLIT)
  {
    return CU_QUAD_SPLIT;
  }

  else if (cuSplitType == CU_HORZ_SPLIT)
  {
    return CU_HORZ_SPLIT;
  }

  else if (cuSplitType == CU_VERT_SPLIT)
  {
    return CU_VERT_SPLIT;
  }

  else if (cuSplitType == CU_TRIH_SPLIT)
  {
    return CU_TRIH_SPLIT;
  }
  else if (cuSplitType == CU_TRIV_SPLIT)
  {
    return CU_TRIV_SPLIT;
  }
  else
  {
    THROW("Unknown split mode");
    return CU_QUAD_SPLIT;
  }
}
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
ModeType CU::getModeTypeAtDepth( const CodingUnit& cu, const unsigned depth )
{
  ModeType modeType = ModeType( (cu.modeTypeSeries >> (depth * 3)) & 0x07 );
  CHECK( depth > cu.depth, " depth is wrong" );
  return modeType;
}
#endif


bool CU::divideTuInRows( const CodingUnit &cu )
{
  CHECK( cu.ispMode != HOR_INTRA_SUBPARTITIONS && cu.ispMode != VER_INTRA_SUBPARTITIONS, "Intra Subpartitions type not recognized!" );
  return cu.ispMode == HOR_INTRA_SUBPARTITIONS ? true : false;
}


PartSplit CU::getISPType( const CodingUnit &cu, const ComponentID compID )
{
  if( cu.ispMode && isLuma( compID ) )
  {
    const bool tuIsDividedInRows = CU::divideTuInRows( cu );

    return tuIsDividedInRows ? TU_1D_HORZ_SPLIT : TU_1D_VERT_SPLIT;
  }
  return TU_NO_ISP;
}

bool CU::isISPLast( const CodingUnit &cu, const CompArea &tuArea, const ComponentID compID )
{
  PartSplit partitionType = CU::getISPType( cu, compID );

  Area originalArea = cu.blocks[compID];
  switch( partitionType )
  {
    case TU_1D_HORZ_SPLIT:
      return tuArea.y + tuArea.height == originalArea.y + originalArea.height;
    case TU_1D_VERT_SPLIT:
      return tuArea.x + tuArea.width == originalArea.x + originalArea.width;
    default:
      THROW( "Unknown ISP processing order type!" );
      return false;
  }
}

bool CU::isISPFirst( const CodingUnit &cu, const CompArea &tuArea, const ComponentID compID )
{
  return tuArea == cu.firstTU->blocks[compID];
}

bool CU::canUseISP( const CodingUnit &cu, const ComponentID compID )
{
  const int width     = cu.blocks[compID].width;
  const int height    = cu.blocks[compID].height;
  const int maxTrSize = cu.cs->sps->getMaxTbSize();
  return CU::canUseISP( width, height, maxTrSize );
}

bool CU::canUseISP( const int width, const int height, const int maxTrSize )
{
  bool  notEnoughSamplesToSplit = ( floorLog2(width) + floorLog2(height) <= ( floorLog2(MIN_TB_SIZEY) << 1 ) );
  bool  cuSizeLargerThanMaxTrSize = width > maxTrSize || height > maxTrSize;
  if ( notEnoughSamplesToSplit || cuSizeLargerThanMaxTrSize )
  {
    return false;
  }
  return true;
}

bool CU::canUseLfnstWithISP( const CompArea& cuArea, const ISPType ispSplitType )
{
  if( ispSplitType == NOT_INTRA_SUBPARTITIONS )
  {
    return false;
  }
  Size tuSize = ( ispSplitType == HOR_INTRA_SUBPARTITIONS ) ? Size( cuArea.width, CU::getISPSplitDim( cuArea.width, cuArea.height, TU_1D_HORZ_SPLIT ) ) :
    Size( CU::getISPSplitDim( cuArea.width, cuArea.height, TU_1D_VERT_SPLIT ), cuArea.height );

  if( !( tuSize.width >= MIN_TB_SIZEY && tuSize.height >= MIN_TB_SIZEY ) )
  {
    return false;
  }
  return true;
}

bool CU::canUseLfnstWithISP( const CodingUnit& cu, const ChannelType chType )
{
  CHECK( !isLuma( chType ), "Wrong ISP mode!" );
  return CU::canUseLfnstWithISP( cu.blocks[chType == CHANNEL_TYPE_LUMA ? 0 : 1], (ISPType)cu.ispMode );
}

#if JVET_W0119_LFNST_EXTENSION
Size CU::getLfnstSize( const CodingUnit& cu, const ChannelType chType )
{
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  int chIdx = CS::isDualITree( *cu.cs ) && cu.chType == CHANNEL_TYPE_CHROMA ? 1 : 0;
#else
  int chIdx = cu.isSepTree() && cu.chType == CHANNEL_TYPE_CHROMA ? 1 : 0;
#endif
  Size lfnstSize;
  if( cu.ispMode )
  {
    CHECK( !isLuma( chType ), "Wrong ISP mode!" );
    const CompArea& cuArea = cu.blocks[ 0 ];
    lfnstSize = ( cu.ispMode == HOR_INTRA_SUBPARTITIONS ) ? Size( cuArea.width, CU::getISPSplitDim( cuArea.width, cuArea.height, TU_1D_HORZ_SPLIT ) ) :
                Size( CU::getISPSplitDim( cuArea.width, cuArea.height, TU_1D_VERT_SPLIT ), cuArea.height );
  }
  else
  {
    lfnstSize.width  = cu.blocks[ chIdx ].lumaSize().width;
    lfnstSize.height = cu.blocks[ chIdx ].lumaSize().height;
  }
  return lfnstSize;
}
#endif

uint32_t CU::getISPSplitDim( const int width, const int height, const PartSplit ispType )
{
  bool divideTuInRows = ispType == TU_1D_HORZ_SPLIT;
  uint32_t splitDimensionSize, nonSplitDimensionSize, partitionSize, divShift = 2;

  if( divideTuInRows )
  {
    splitDimensionSize    = height;
    nonSplitDimensionSize = width;
  }
  else
  {
    splitDimensionSize    = width;
    nonSplitDimensionSize = height;
  }

  const int minNumberOfSamplesPerCu = 1 << ( ( floorLog2(MIN_TB_SIZEY) << 1 ) );
  const int factorToMinSamples = nonSplitDimensionSize < minNumberOfSamplesPerCu ? minNumberOfSamplesPerCu >> floorLog2(nonSplitDimensionSize) : 1;
  partitionSize = ( splitDimensionSize >> divShift ) < factorToMinSamples ? factorToMinSamples : ( splitDimensionSize >> divShift );

  CHECK( floorLog2(partitionSize) + floorLog2(nonSplitDimensionSize) < floorLog2(minNumberOfSamplesPerCu), "A partition has less than the minimum amount of samples!" );
  return partitionSize;
}

bool CU::allLumaCBFsAreZero(const CodingUnit& cu)
{
  if (!cu.ispMode)
  {
    return TU::getCbf(*cu.firstTU, COMPONENT_Y) == false;
  }
  else
  {
    int numTotalTUs = cu.ispMode == HOR_INTRA_SUBPARTITIONS ? cu.lheight() >> floorLog2(cu.firstTU->lheight()) : cu.lwidth() >> floorLog2(cu.firstTU->lwidth());
    TransformUnit* tuPtr = cu.firstTU;
    for (int tuIdx = 0; tuIdx < numTotalTUs; tuIdx++)
    {
      if (TU::getCbf(*tuPtr, COMPONENT_Y) == true)
      {
        return false;
      }
      tuPtr = tuPtr->next;
    }
    return true;
  }
}

#if JVET_W0123_TIMD_FUSION
TEMPLATE_TYPE CU::deriveTimdRefType( int iCurX, int iCurY, uint32_t uiCurWidth, uint32_t uiCurHeight, int iTemplateWidth, int iTemplateHeight, int& iRefX, int& iRefY, uint32_t& uiRefWidth, uint32_t& uiRefHeight )
{
  if(iCurX == 0 && iCurY == 0)
  {
    return NO_NEIGHBOR;
  }

  TEMPLATE_TYPE eTempType = NO_NEIGHBOR;
  iRefX = iRefY = -1;
  if(iCurX > 0 && iCurY > 0)
  {
    iRefX       = iCurX - iTemplateWidth;
    iRefY       = iCurY - iTemplateHeight;
    uiRefWidth  = uiCurWidth + iTemplateWidth;
    uiRefHeight = uiCurHeight + iTemplateHeight;
    eTempType   = LEFT_ABOVE_NEIGHBOR;
  }
  else if(iCurX == 0 && iCurY > 0)
  {
    iRefX       = iCurX;
    iRefY       = iCurY - iTemplateHeight;
    uiRefWidth  = uiCurWidth;
    uiRefHeight = uiCurHeight;
    eTempType   = ABOVE_NEIGHBOR;
  }
  else if(iCurX > 0 && iCurY == 0)
  {
    iRefX       = iCurX - iTemplateWidth;
    iRefY       = iCurY;
    uiRefWidth  = uiCurWidth;
    uiRefHeight = uiCurHeight;
    eTempType   = LEFT_NEIGHBOR;
  }
  else
  {
    assert(0);
  }
  return eTempType;
}
#endif


PUTraverser CU::traversePUs( CodingUnit& cu )
{
  return PUTraverser( cu.firstPU, cu.lastPU->next );
}

TUTraverser CU::traverseTUs( CodingUnit& cu )
{
  return TUTraverser( cu.firstTU, cu.lastTU->next );
}

cPUTraverser CU::traversePUs( const CodingUnit& cu )
{
  return cPUTraverser( cu.firstPU, cu.lastPU->next );
}

cTUTraverser CU::traverseTUs( const CodingUnit& cu )
{
  return cTUTraverser( cu.firstTU, cu.lastTU->next );
}

// PU tools

#if SECONDARY_MPM
int PU::getIntraMPMs( const PredictionUnit &pu, uint8_t* mpm, uint8_t* non_mpm, const ChannelType &channelType /*= CHANNEL_TYPE_LUMA*/ )
#else
int PU::getIntraMPMs(const PredictionUnit &pu, unsigned* mpm, const ChannelType &channelType /*= CHANNEL_TYPE_LUMA*/)
#endif
{
#if SECONDARY_MPM
  bool includedMode[NUM_INTRA_MODE];
  memset(includedMode, false, sizeof(includedMode));

  int numValidMPM = 0;
  mpm[numValidMPM++] = PLANAR_IDX;
  includedMode[PLANAR_IDX] = true;
#endif

  const int numMPMs = NUM_MOST_PROBABLE_MODES;
  {
    CHECK(channelType != CHANNEL_TYPE_LUMA, "Not harmonized yet");
    int numCand      = -1;
#if !SECONDARY_MPM
    int leftIntraDir = PLANAR_IDX, aboveIntraDir = PLANAR_IDX;
#endif

    const CompArea &area = pu.block(getFirstComponentOfChannel(channelType));
    const Position posRT = area.topRight();
    const Position posLB = area.bottomLeft();

    // Get intra direction of left PU
#if SECONDARY_MPM
    const PredictionUnit *puLeft = (pu.lheight() >= pu.lwidth())
      ? pu.cs->getPURestricted(posRT.offset(0, -1), pu, channelType)
      : pu.cs->getPURestricted(posLB.offset(-1, 0), pu, channelType);
#else
    const PredictionUnit *puLeft = pu.cs->getPURestricted(posLB.offset(-1, 0), pu, channelType);
#endif
    if (puLeft && CU::isIntra(*puLeft->cu))
    {
#if SECONDARY_MPM
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puLeft->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puLeft)) : PU::getIntraDirLuma(*puLeft);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puLeft);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
#else
      leftIntraDir = PU::getIntraDirLuma( *puLeft );
#if JVET_W0123_TIMD_FUSION
      if (puLeft->cu->timd)
      {
        leftIntraDir = MAP131TO67(leftIntraDir);
      }
#endif
#endif
    }
#if JVET_W0123_TIMD_FUSION && !SECONDARY_MPM
    if (puLeft && CU::isInter(*puLeft->cu))
    {
      leftIntraDir = puLeft->getIpmInfo(posLB.offset(-1, 0));
    }
#endif

    // Get intra direction of above PU
#if SECONDARY_MPM
    const PredictionUnit *puAbove = (pu.lheight() >= pu.lwidth())
      ? pu.cs->getPURestricted(posLB.offset(-1, 0), pu, channelType)
      : pu.cs->getPURestricted(posRT.offset(0, -1), pu, channelType);
#else
    const PredictionUnit *puAbove = pu.cs->getPURestricted(posRT.offset(0, -1), pu, channelType);
#endif
    if (puAbove && CU::isIntra(*puAbove->cu) && CU::isSameCtu(*pu.cu, *puAbove->cu))
    {
#if SECONDARY_MPM
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puAbove->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puAbove)) : PU::getIntraDirLuma(*puAbove);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puAbove);
#endif
      if (!includedMode[mpm[numValidMPM]])
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
#else
      aboveIntraDir = PU::getIntraDirLuma(*puAbove);
#if JVET_W0123_TIMD_FUSION
      if (puAbove->cu->timd)
      {
        aboveIntraDir = MAP131TO67(aboveIntraDir);
      }
#endif
#endif
    }
#if JVET_W0123_TIMD_FUSION && !SECONDARY_MPM
    if (puAbove && CU::isInter(*puAbove->cu))
    {
      aboveIntraDir = puAbove->getIpmInfo(posRT.offset(0, -1));
    }
#endif

#if SECONDARY_MPM
#if JVET_W0123_TIMD_FUSION
    if (puLeft && CU::isInter(*puLeft->cu))
    {
      mpm[numValidMPM] = puLeft->getIpmInfo(pu.lheight() >= pu.lwidth() ? posRT.offset(0, -1) : posLB.offset(-1, 0));
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }
    if (puAbove && CU::isInter(*puAbove->cu))
    {
      mpm[numValidMPM] = puAbove->getIpmInfo(pu.lheight() >= pu.lwidth() ? posLB.offset(-1, 0) : posRT.offset(0, -1));
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }
#endif
    // Get intra direction of below-left PU
    const PredictionUnit *puBelowLeft = pu.cs->getPURestricted(posLB.offset(-1, 1), pu, channelType);
    if (puBelowLeft && CU::isIntra(*puBelowLeft->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puBelowLeft->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puBelowLeft)) : PU::getIntraDirLuma(*puBelowLeft);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puBelowLeft);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }

    // Get intra direction of above-right PU
    const PredictionUnit *puAboveRight = pu.cs->getPURestricted(posRT.offset(1, -1), pu, channelType);
    if (puAboveRight && CU::isIntra(*puAboveRight->cu) && CU::isSameCtu(*pu.cu, *puAboveRight->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puAboveRight->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puAboveRight)) : PU::getIntraDirLuma(*puAboveRight);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puAboveRight);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }

    // Get intra direction of above-left PU
    const Position posTL = area.topLeft();
    const PredictionUnit *puAboveLeft = pu.cs->getPURestricted(posTL.offset(-1, -1), pu, channelType);
    if (puAboveLeft && CU::isIntra(*puAboveLeft->cu) && CU::isSameCtu(*pu.cu, *puAboveLeft->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puAboveLeft->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puAboveLeft)) : PU::getIntraDirLuma(*puAboveLeft);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puAboveLeft);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }
#if JVET_W0123_TIMD_FUSION
    if (puBelowLeft && CU::isInter(*puBelowLeft->cu))
    {
      mpm[numValidMPM] = puBelowLeft->getIpmInfo(posLB.offset(-1, 1));
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }
    if (puAboveRight && CU::isInter(*puAboveRight->cu))
    {
      mpm[numValidMPM] = puAboveRight->getIpmInfo(posRT.offset(1, -1));
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }
    if (puAboveLeft && CU::isInter(*puAboveLeft->cu))
    {
      mpm[numValidMPM] = puAboveLeft->getIpmInfo(posTL.offset(-1, -1));
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }
#endif
#endif

    CHECK(2 >= numMPMs, "Invalid number of most probable modes");

    const int offset = ( int ) NUM_LUMA_MODE - 6;
    const int mod = offset + 3;

    {
#if SECONDARY_MPM
      numCand = numValidMPM;
#else
      mpm[0] = PLANAR_IDX;
      mpm[1] = DC_IDX;
      mpm[2] = VER_IDX;
      mpm[3] = HOR_IDX;
      mpm[4] = VER_IDX - 4;
      mpm[5] = VER_IDX + 4;
#endif
#if ENABLE_DIMD && SECONDARY_MPM
      //adding dimd modes
      if (pu.cu->slice->getSPS()->getUseDimd())
      {
        if (pu.cu->dimdMode != -1)
        {
          mpm[numValidMPM] = pu.cu->dimdMode;
          if( !includedMode[mpm[numValidMPM]] )
          {
            includedMode[mpm[numValidMPM++]] = true;
          }
        }

        if (pu.cu->dimdBlendMode[0] != -1)
        {
          mpm[numValidMPM] = pu.cu->dimdBlendMode[0];
          if( !includedMode[mpm[numValidMPM]] )
          {
            includedMode[mpm[numValidMPM++]] = true;
          }
        }
      }
#endif

#if SECONDARY_MPM
      bool checkDCEnabled = false;

      // Derived modes of mpm[1]
      if (numCand >= 2)
      {
        if (mpm[1] > DC_IDX)
        {
          for (int i = 0; i < 4 && numValidMPM < numMPMs; i++)
          {
            mpm[numValidMPM] = ((mpm[1] + offset - i) % mod) + 2;
            if( !includedMode[mpm[numValidMPM]] )
            {
              includedMode[mpm[numValidMPM++]] = true;
            }

            if( numValidMPM >= numMPMs )
            {
              break;
            }

            mpm[numValidMPM] = ((mpm[1] - 1 + i) % mod) + 2;
            if( !includedMode[mpm[numValidMPM]] )
            {
              includedMode[mpm[numValidMPM++]] = true;
            }
          }
        }
        else if( mpm[1] == DC_IDX )
        {
          checkDCEnabled = true;
        }
      }


      // Derived modes of mpm[2]
      if (numCand >= 3)
      {
        if (mpm[2] > DC_IDX)
        {
          for (int i = 0; i < 4 && numValidMPM < numMPMs; i++)
          {
            mpm[numValidMPM] = ((mpm[2] + offset - i) % mod) + 2;
            if( !includedMode[mpm[numValidMPM]] )
            {
              includedMode[mpm[numValidMPM++]] = true;
            }

            if (numValidMPM >= numMPMs)
              break;

            mpm[numValidMPM] = ((mpm[2] - 1 + i) % mod) + 2;
            if( !includedMode[mpm[numValidMPM]] )
            {
              includedMode[mpm[numValidMPM++]] = true;
            }
          }
        }
        else if( mpm[2] == DC_IDX )
        {
          checkDCEnabled = true;
        }
      }


      // Derived modes of mpm[3]
      if (checkDCEnabled && numCand >= 4 && mpm[3] > DC_IDX)
      {
        for (int i = 0; i < 3 && numValidMPM < numMPMs; i++)
        {
          mpm[numValidMPM] = ((mpm[3] + offset - i) % mod) + 2;
          if( !includedMode[mpm[numValidMPM]] )
          {
            includedMode[mpm[numValidMPM++]] = true;
          }

          if( numValidMPM >= numMPMs )
          {
            break;
          }

          mpm[numValidMPM] = ((mpm[3] - 1 + i) % mod) + 2;
          if( !includedMode[mpm[numValidMPM]] )
          {
            includedMode[mpm[numValidMPM++]] = true;
          }
        }
      }
#else
      if (leftIntraDir == aboveIntraDir)
      {
        numCand = 1;
        if (leftIntraDir > DC_IDX)
        {
          mpm[0] = PLANAR_IDX;
          mpm[1] = leftIntraDir;
          mpm[2] = ((leftIntraDir + offset) % mod) + 2;
          mpm[3] = ((leftIntraDir - 1) % mod) + 2;
          mpm[4] = ((leftIntraDir + offset - 1) % mod) + 2;
          mpm[5] = ( leftIntraDir               % mod) + 2;
        }
      }
      else //L!=A
      {
        numCand = 2;
        int  maxCandModeIdx = mpm[0] > mpm[1] ? 0 : 1;

        if ((leftIntraDir > DC_IDX) && (aboveIntraDir > DC_IDX))
        {
          mpm[0] = PLANAR_IDX;
          mpm[1] = leftIntraDir;
          mpm[2] = aboveIntraDir;
          maxCandModeIdx = mpm[1] > mpm[2] ? 1 : 2;
          int minCandModeIdx = mpm[1] > mpm[2] ? 2 : 1;
          if (mpm[maxCandModeIdx] - mpm[minCandModeIdx] == 1)
          {
            mpm[3] = ((mpm[minCandModeIdx] + offset)     % mod) + 2;
            mpm[4] = ((mpm[maxCandModeIdx] - 1)          % mod) + 2;
            mpm[5] = ((mpm[minCandModeIdx] + offset - 1) % mod) + 2;
          }
          else if (mpm[maxCandModeIdx] - mpm[minCandModeIdx] >= 62)
          {
            mpm[3] = ((mpm[minCandModeIdx] - 1)      % mod) + 2;
            mpm[4] = ((mpm[maxCandModeIdx] + offset) % mod) + 2;
            mpm[5] = ( mpm[minCandModeIdx]           % mod) + 2;
          }
          else if (mpm[maxCandModeIdx] - mpm[minCandModeIdx] == 2)
          {
            mpm[3] = ((mpm[minCandModeIdx] - 1)      % mod) + 2;
            mpm[4] = ((mpm[minCandModeIdx] + offset) % mod) + 2;
            mpm[5] = ((mpm[maxCandModeIdx] - 1)      % mod) + 2;
          }
          else
          {
            mpm[3] = ((mpm[minCandModeIdx] + offset) % mod) + 2;
            mpm[4] = ((mpm[minCandModeIdx] - 1)      % mod) + 2;
            mpm[5] = ((mpm[maxCandModeIdx] + offset) % mod) + 2;
          }
        }
        else if (leftIntraDir + aboveIntraDir >= 2)
        {
          mpm[0] = PLANAR_IDX;
          mpm[1] = (leftIntraDir < aboveIntraDir) ? aboveIntraDir : leftIntraDir;
          maxCandModeIdx = 1;
          mpm[2] = ((mpm[maxCandModeIdx] + offset)     % mod) + 2;
          mpm[3] = ((mpm[maxCandModeIdx] - 1)          % mod) + 2;
          mpm[4] = ((mpm[maxCandModeIdx] + offset - 1) % mod) + 2;
          mpm[5] = ( mpm[maxCandModeIdx]               % mod) + 2;
        }
      }
#endif

#if SECONDARY_MPM
      unsigned mpm_default[numMPMs - 1] = { DC_IDX, VER_IDX, HOR_IDX, VER_IDX - 4, VER_IDX + 4, 14, 22, 42, 58, 10, 26,
                                           38, 62, 6, 30, 34, 66, 2, 48, 52, 16 };
      for (int idx = 0; (idx < numMPMs - 1) && numValidMPM < numMPMs; idx++)
      {
        mpm[numValidMPM] = mpm_default[idx];
        if( !includedMode[mpm[numValidMPM]] )
        {
          includedMode[mpm[numValidMPM++]] = true;
        }
      }

      int numNonMPM = 0;
      for (int idx = 0; idx < NUM_LUMA_MODE; idx++)
      {
        if( !includedMode[idx] )
        {
          non_mpm[numNonMPM++] = idx;
        }
      }
#endif
    }
    for (int i = 0; i < numMPMs; i++)
    {
      CHECK(mpm[i] >= NUM_LUMA_MODE, "Invalid MPM");
    }
    CHECK(numCand == 0, "No candidates found");
    return numCand;
  }
}

#if JVET_Y0065_GPM_INTRA
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
void PU::getGeoIntraMPMs(const PredictionUnit &pu, uint8_t* mpm, uint8_t splitDir, uint8_t shape, bool doInit, bool doInitAL, bool doInitA, bool doInitL)
{
  static uint8_t  partialMPMsAll[GEO_NUM_TM_MV_CAND - 1][GEO_MAX_NUM_INTRA_CANDS]; // [0] for above-left, [1] for above [1] for left
  if (doInit)
  {
    if (doInitAL)
    {
      PU::getGeoIntraMPMs(pu, partialMPMsAll[0], GEO_NUM_PARTITION_MODE, GEO_TM_SHAPE_AL);
    }
    if (doInitA)
    {
      PU::getGeoIntraMPMs(pu, partialMPMsAll[1], GEO_NUM_PARTITION_MODE, GEO_TM_SHAPE_A);
    }
    if (doInitL)
    {
      PU::getGeoIntraMPMs(pu, partialMPMsAll[2], GEO_NUM_PARTITION_MODE, GEO_TM_SHAPE_L);
    }
  }

  const  uint8_t* partialMPMs = partialMPMsAll[shape - 1];
  uint8_t numValidMPM = 1;
  mpm[0] = g_geoAngle2IntraAng[g_GeoParams[splitDir][0]];
  for (int i = 0; i < GEO_MAX_NUM_INTRA_CANDS; ++i)
  {
    if (partialMPMs[i] == NOMODE_IDX)
    {
      break;
    }

    if (partialMPMs[i] != mpm[0])
    {
      mpm[numValidMPM++] = partialMPMs[i];
      if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
      {
        return;
      }
    }
  }

  mpm[numValidMPM] = (mpm[0] > DIA_IDX) ? (mpm[0] - 32) : (mpm[0] + 32);
  for (int i = 1; i < numValidMPM; ++i)
  {
    if (mpm[numValidMPM] == mpm[i])
    {
      --numValidMPM;
      break;
    }
  }
  ++numValidMPM;
  if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
  {
    return;
  }

  mpm[numValidMPM] = PLANAR_IDX;
}
#endif

void PU::getGeoIntraMPMs( const PredictionUnit &pu, uint8_t* mpm, uint8_t splitDir, uint8_t shape )
{
  bool includedMode[NUM_INTRA_MODE];
  memset(includedMode, false, sizeof(includedMode));

  int numValidMPM = 0;
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  bool outputFullMPMs = splitDir < GEO_NUM_PARTITION_MODE;
  if(outputFullMPMs)
  {
#endif
  mpm[numValidMPM++] = g_geoAngle2IntraAng[g_GeoParams[splitDir][0]];
  includedMode[mpm[0]] = true;
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  }
#endif

  CodingUnit* cu = pu.cu;
#if ENABLE_DIMD
  if (cu->slice->getSPS()->getUseDimd())
  {
    if (cu->dimdMode != -1)
    {
      mpm[numValidMPM] = cu->dimdMode;
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
      }
    }
  }
#endif
#if JVET_W0123_TIMD_FUSION
  if (cu->slice->getSPS()->getUseTimd() && cu->timdMode != -1)
  {
    mpm[numValidMPM] = MAP131TO67(cu->timdMode);
    if( !includedMode[mpm[numValidMPM]] )
    {
      includedMode[mpm[numValidMPM++]] = true;
      if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
      {
        return;
      }
    }
  }
#endif
  const CompArea &area = pu.block(COMPONENT_Y);
  const Position posA = area.topRight().offset(0, -1);
  const Position posAR = area.topRight().offset(1, -1);
  const Position posL = area.bottomLeft().offset(-1, 0);
  const Position posBL = area.bottomLeft().offset(-1, 1);
  const Position posAL = area.topLeft().offset(-1, -1);

  if (shape == GEO_TM_SHAPE_L || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puLeft = pu.cs->getPURestricted(posL, pu, CHANNEL_TYPE_LUMA);
    if (puLeft && CU::isIntra(*puLeft->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puLeft->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puLeft)) : PU::getIntraDirLuma(*puLeft);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puLeft);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

  if (shape == GEO_TM_SHAPE_A || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puAbove = pu.cs->getPURestricted(posA, pu, CHANNEL_TYPE_LUMA);
    if (puAbove && CU::isIntra(*puAbove->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puAbove->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puAbove)) : PU::getIntraDirLuma(*puAbove);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puAbove);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

#if JVET_W0123_TIMD_FUSION
  if (shape == GEO_TM_SHAPE_L || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puLeft = pu.cs->getPURestricted(posL, pu, CHANNEL_TYPE_LUMA);
    if (puLeft && CU::isInter(*puLeft->cu))
    {
      mpm[numValidMPM] = puLeft->getIpmInfo(posL);
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

  if (shape == GEO_TM_SHAPE_A || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puAbove = pu.cs->getPURestricted(posA, pu, CHANNEL_TYPE_LUMA);
    if (puAbove && CU::isInter(*puAbove->cu))
    {
      mpm[numValidMPM] = puAbove->getIpmInfo(posA);
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }
#endif

  if (shape == GEO_TM_SHAPE_L || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puBelowLeft = pu.cs->getPURestricted(posBL, pu, CHANNEL_TYPE_LUMA);
    if (puBelowLeft && CU::isIntra(*puBelowLeft->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puBelowLeft->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puBelowLeft)) : PU::getIntraDirLuma(*puBelowLeft);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puBelowLeft);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

  if (shape == GEO_TM_SHAPE_A || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puAboveRight = pu.cs->getPURestricted(posAR, pu, CHANNEL_TYPE_LUMA);
    if (puAboveRight && CU::isIntra(*puAboveRight->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puAboveRight->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puAboveRight)) : PU::getIntraDirLuma(*puAboveRight);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puAboveRight);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

  {
    const PredictionUnit *puAboveLeft = pu.cs->getPURestricted(posAL, pu, CHANNEL_TYPE_LUMA);
    if (puAboveLeft && CU::isIntra(*puAboveLeft->cu))
    {
#if JVET_W0123_TIMD_FUSION
      mpm[numValidMPM] = puAboveLeft->cu->timd ? MAP131TO67(PU::getIntraDirLuma(*puAboveLeft)) : PU::getIntraDirLuma(*puAboveLeft);
#else
      mpm[numValidMPM] = PU::getIntraDirLuma(*puAboveLeft);
#endif
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

#if JVET_W0123_TIMD_FUSION
  if (shape == GEO_TM_SHAPE_L || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puBelowLeft = pu.cs->getPURestricted(posBL, pu, CHANNEL_TYPE_LUMA);
    if (puBelowLeft && CU::isInter(*puBelowLeft->cu))
    {
      mpm[numValidMPM] = puBelowLeft->getIpmInfo(posBL);
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

  if (shape == GEO_TM_SHAPE_A || shape == GEO_TM_SHAPE_AL)
  {
    const PredictionUnit *puAboveRight = pu.cs->getPURestricted(posAR, pu, CHANNEL_TYPE_LUMA);
    if (puAboveRight && CU::isInter(*puAboveRight->cu))
    {
      mpm[numValidMPM] = puAboveRight->getIpmInfo(posAR);
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }

  {
    const PredictionUnit *puAboveLeft = pu.cs->getPURestricted(posAL, pu, CHANNEL_TYPE_LUMA);
    if (puAboveLeft && CU::isInter(*puAboveLeft->cu))
    {
      mpm[numValidMPM] = puAboveLeft->getIpmInfo(posAL);
      if( !includedMode[mpm[numValidMPM]] )
      {
        includedMode[mpm[numValidMPM++]] = true;
        if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
        {
          return;
        }
      }
    }
  }
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  if (!outputFullMPMs)
  {
    for (int i = numValidMPM; i < GEO_MAX_NUM_INTRA_CANDS; ++i)
    {
      mpm[i] = NOMODE_IDX;
    }
    return;
  }
#endif

  mpm[numValidMPM] = (mpm[0] > DIA_IDX) ? (mpm[0]-32) : (mpm[0]+32);
  if( !includedMode[mpm[numValidMPM]] )
  {
    includedMode[mpm[numValidMPM++]] = true;
    if (numValidMPM == GEO_MAX_NUM_INTRA_CANDS)
    {
      return;
    }
  }
  mpm[numValidMPM] = PLANAR_IDX;
}
#endif

bool PU::isMIP(const PredictionUnit &pu, const ChannelType &chType)
{
  if (chType == CHANNEL_TYPE_LUMA)
  {
    // Default case if chType is omitted.
    return pu.cu->mipFlag;
  }
  else
  {
    return isDMChromaMIP(pu) && (pu.intraDir[CHANNEL_TYPE_CHROMA] == DM_CHROMA_IDX);
  }
}
#if JVET_V0130_INTRA_TMP
bool PU::isTmp(const PredictionUnit& pu, const ChannelType& chType)
{
	return (chType == CHANNEL_TYPE_LUMA && pu.cu->tmpFlag);
}
#endif
bool PU::isDMChromaMIP(const PredictionUnit &pu)
{
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  return !pu.cu->isSepTree() && (pu.chromaFormat == CHROMA_444) && getCoLocatedLumaPU(pu).cu->mipFlag;
#else
  return !(CS::isDualITree(*pu.cs)) && (pu.chromaFormat == CHROMA_444) && getCoLocatedLumaPU(pu).cu->mipFlag;
#endif
}

uint32_t PU::getIntraDirLuma( const PredictionUnit &pu )
{
#if JVET_V0130_INTRA_TMP
	if (isMIP(pu) || isTmp(pu))
#else
  if (isMIP(pu))
#endif
  {
    return PLANAR_IDX;
  }
  else
  {
    return pu.intraDir[CHANNEL_TYPE_LUMA];
  }
}

void PU::getIntraChromaCandModes(const PredictionUnit &pu, unsigned modeList[NUM_CHROMA_MODE])
{
  modeList[0] = PLANAR_IDX;
  modeList[1] = VER_IDX;
  modeList[2] = HOR_IDX;
  modeList[3] = DC_IDX;
  modeList[4] = LM_CHROMA_IDX;
#if MMLM
  modeList[5] = MDLM_L_IDX;
  modeList[6] = MDLM_T_IDX;
  modeList[7] = MMLM_CHROMA_IDX;
  modeList[8] = MMLM_L_IDX;
  modeList[9] = MMLM_T_IDX;
#if JVET_Z0050_DIMD_CHROMA_FUSION && ENABLE_DIMD
  modeList[10] = DM_CHROMA_IDX;
  modeList[11] = DIMD_CHROMA_IDX;
#else
  modeList[10] = DM_CHROMA_IDX;
#endif
#else
  modeList[5] = MDLM_L_IDX;
  modeList[6] = MDLM_T_IDX;
#if JVET_Z0050_DIMD_CHROMA_FUSION && ENABLE_DIMD
  modeList[7] = DM_CHROMA_IDX;
  modeList[8] = DIMD_CHROMA_IDX;
#else
  modeList[7] = DM_CHROMA_IDX;
#endif
#endif

  // If Direct Mode is MIP, mode cannot be already in the list.
  if (isDMChromaMIP(pu))
  {
    return;
  }

  const uint32_t lumaMode = getCoLocatedIntraLumaMode(pu);
  for (int i = 0; i < 4; i++)
  {
    if (lumaMode == modeList[i])
    {
      modeList[i] = VDIA_IDX;
      break;
    }
  }
}
  
#if JVET_Z0050_CCLM_SLOPE
  bool PU::hasCclmDeltaFlag(const PredictionUnit &pu, const int mode)
{
  int  chrMode      = mode < 0 ? pu.intraDir[1] : mode;
#if MMLM
  bool hasDeltaFlag = chrMode == LM_CHROMA_IDX || chrMode == MMLM_CHROMA_IDX;
#else
  bool hasDeltaFlag = chrMode == LM_CHROMA_IDX;
#endif
  hasDeltaFlag     &= pu.Cb().width * pu.Cb().height >= 128;

  return hasDeltaFlag;
}
#endif

bool PU::isLMCMode(unsigned mode)
{
#if MMLM
  return (mode >= LM_CHROMA_IDX && mode <= MMLM_T_IDX);
#else
  return (mode >= LM_CHROMA_IDX && mode <= MDLM_T_IDX);
#endif
}
#if MMLM
bool PU::isMultiModeLM(unsigned mode)
{
  return (mode == MMLM_CHROMA_IDX || mode == MMLM_L_IDX || mode == MMLM_T_IDX);
}
#endif
bool PU::isLMCModeEnabled(const PredictionUnit &pu, unsigned mode)
{
#if CCLM_LATENCY_RESTRICTION_RMV
  if (pu.cs->sps->getUseLMChroma() )
#else
  if ( pu.cs->sps->getUseLMChroma() && pu.cu->checkCCLMAllowed() )
#endif
  {
    return true;
  }
  return false;
}

int PU::getLMSymbolList(const PredictionUnit &pu, int *modeList)
{
  int idx = 0;

  modeList[idx++] = LM_CHROMA_IDX;
#if MMLM
  modeList[idx++] = MMLM_CHROMA_IDX;
#endif
  modeList[idx++] = MDLM_L_IDX;
  modeList[idx++] = MDLM_T_IDX;
#if MMLM
  modeList[idx++] = MMLM_L_IDX;
  modeList[idx++] = MMLM_T_IDX;
#endif
  return idx;
}

bool PU::isChromaIntraModeCrossCheckMode( const PredictionUnit &pu )
{
  return !pu.cu->bdpcmModeChroma && pu.intraDir[CHANNEL_TYPE_CHROMA] == DM_CHROMA_IDX;
}

uint32_t PU::getFinalIntraMode( const PredictionUnit &pu, const ChannelType &chType )
{
  uint32_t uiIntraMode = pu.intraDir[chType];

  if( uiIntraMode == DM_CHROMA_IDX && !isLuma( chType ) )
  {
    uiIntraMode = getCoLocatedIntraLumaMode(pu);
  }
  if( pu.chromaFormat == CHROMA_422 && !isLuma( chType ) && uiIntraMode < NUM_LUMA_MODE ) // map directional, planar and dc
  {
    uiIntraMode = g_chroma422IntraAngleMappingTable[uiIntraMode];
  }
#if JVET_Z0050_DIMD_CHROMA_FUSION && ENABLE_DIMD
  if (uiIntraMode == DIMD_CHROMA_IDX && !isLuma(chType))
  {
    uiIntraMode = pu.cu->dimdChromaMode;
  }
#endif
  return uiIntraMode;
}

const PredictionUnit &PU::getCoLocatedLumaPU(const PredictionUnit &pu)
{
  Position              topLeftPos = pu.blocks[pu.chType].lumaPos();
  Position              refPos     = topLeftPos.offset(pu.blocks[pu.chType].lumaSize().width  >> 1,
                                                       pu.blocks[pu.chType].lumaSize().height >> 1);
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const PredictionUnit &lumaPU     = CS::isDualITree(*pu.cs) ? *pu.cs->picture->cs->getPU(refPos, CHANNEL_TYPE_LUMA)
                                                        : *pu.cs->getPU(topLeftPos, CHANNEL_TYPE_LUMA);
#else
  const PredictionUnit &lumaPU     = pu.cu->isSepTree() ? *pu.cs->picture->cs->getPU(refPos, CHANNEL_TYPE_LUMA)
                                                        : *pu.cs->getPU(topLeftPos, CHANNEL_TYPE_LUMA);
#endif
  return lumaPU;
}

uint32_t PU::getCoLocatedIntraLumaMode(const PredictionUnit &pu)
{
#if JVET_W0123_TIMD_FUSION
  if (PU::getCoLocatedLumaPU(pu).cu->timd)
  {
    return MAP131TO67(PU::getIntraDirLuma(PU::getCoLocatedLumaPU(pu)));
  }
#endif
  return PU::getIntraDirLuma(PU::getCoLocatedLumaPU(pu));
}

int PU::getWideAngle( const TransformUnit &tu, const uint32_t dirMode, const ComponentID compID )
{
  //This function returns a wide angle index taking into account that the values 0 and 1 are reserved
  //for Planar and DC respectively, as defined in the Spec. Text.
  if( dirMode < 2 )
  {
    return ( int ) dirMode;
  }

  const CompArea&  area         = tu.cu->ispMode && isLuma(compID) ? tu.cu->blocks[compID] : tu.blocks[ compID ];
  int              width        = area.width;
  int              height       = area.height;
  int              modeShift[ ] = { 0, 6, 10, 12, 14, 15 };
  int              deltaSize    = abs( floorLog2( width ) - floorLog2( height ) );
  int              predMode     = dirMode;

  if( width > height && dirMode < 2 + modeShift[ deltaSize ] )
  {
    predMode += ( VDIA_IDX - 1 );
  }
  else if( height > width && predMode > VDIA_IDX - modeShift[ deltaSize ] )
  {
    predMode -= ( VDIA_IDX + 1 );
  }

  return predMode;
}

#if JVET_W0119_LFNST_EXTENSION
int PU::getLFNSTMatrixDim( int width, int height )
{
  int dimension = ( width == 8 && height == 8 ) ? 16 : ( ( ( width >= 8 ) && ( height >= 8 ) ) ? 32 : 16 );

  if( ( width >= 16 ) && ( height >= 16 ) )
  {
    dimension = L16H;
  }

  return dimension;
}

bool PU::getUseLFNST8( int width, int height )
{
  return ( width >= 8 ) && ( height >= 8 );
}

bool PU::getUseLFNST16( int width, int height )
{
  return ( width >= 16 ) && ( height >= 16 );
}

uint8_t PU::getLFNSTIdx( int intraMode, int mtsMode )
{
  return g_lfnstLut[ intraMode ];
}
#endif

bool PU::addMergeHMVPCand(const CodingStructure &cs, MergeCtx &mrgCtx, const int &mrgCandIdx,
                          const uint32_t maxNumMergeCandMin1, int &cnt, const bool isAvailableA1,
                          const MotionInfo miLeft, const bool isAvailableB1, const MotionInfo miAbove
#if !JVET_Z0075_IBC_HMVP_ENLARGE
                        , const bool ibcFlag, const bool isGt4x4
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE || (JVET_Y0058_IBC_LIST_MODIFY && !JVET_Z0075_IBC_HMVP_ENLARGE)
                        , const PredictionUnit &pu
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE
#if !JVET_Y0128_NON_CTC
                        , const int curPoc
                        , const int amvpPoc
#endif
#endif
#if TM_MRG || (JVET_Z0084_IBC_TM && !JVET_Z0075_IBC_HMVP_ENLARGE)
                        , const uint32_t mvdSimilarityThresh
#endif
  )
{
  const Slice& slice = *cs.slice;
  MotionInfo miNeighbor;

#if JVET_Z0075_IBC_HMVP_ENLARGE
#if JVET_Z0118_GDR  
  bool isClean = cs.isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

#if JVET_Z0118_GDR  
  auto &lut = (isClean) ? cs.motionLut.lut1 : cs.motionLut.lut0;
#else
  auto &lut = cs.motionLut.lut;
#endif
#else
#if JVET_Z0118_GDR  
  auto &lut = ibcFlag ? (isClean ? cs.motionLut.lutIbc1 : cs.motionLut.lutIbc0) : (isClean ? cs.motionLut.lut1 : cs.motionLut.lut0);
#else
  auto &lut = ibcFlag ? cs.motionLut.lutIbc : cs.motionLut.lut;
#endif

#endif

  int num_avai_candInLUT = (int)lut.size();

  for (int mrgIdx = 1; mrgIdx <= num_avai_candInLUT; mrgIdx++)
  {
    miNeighbor = lut[num_avai_candInLUT - mrgIdx];

#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
    bool isValidAmMode = checkIsValidMergeMvCand(pu, miNeighbor.refIdx);
#else
    bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, miNeighbor.refIdx);
#endif
    if (isValidAmMode &&
#if JVET_Z0075_IBC_HMVP_ENLARGE
        ( mrgIdx > 2
#else
        ( mrgIdx > 2 || ((mrgIdx > 1 || !isGt4x4) && ibcFlag)
#endif
          || ((!isAvailableA1 || (miLeft != miNeighbor)) && (!isAvailableB1 || (miAbove != miNeighbor))) )
          )
#else
#if JVET_Z0075_IBC_HMVP_ENLARGE
    if ( mrgIdx > 2
#else
    if ( mrgIdx > 2 || ((mrgIdx > 1 || !isGt4x4) && ibcFlag)
#endif
      || ((!isAvailableA1 || (miLeft != miNeighbor)) && (!isAvailableB1 || (miAbove != miNeighbor))) )
#endif
    {
#if JVET_Y0058_IBC_LIST_MODIFY && JVET_Z0084_IBC_TM && !JVET_Z0075_IBC_HMVP_ENLARGE
      if (!ibcFlag || checkIsIBCCandidateValid(pu, miNeighbor))
      {
#endif
      mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
#if JVET_Z0075_IBC_HMVP_ENLARGE
      mrgCtx.useAltHpelIf      [cnt] = miNeighbor.useAltHpelIf;
#else
      mrgCtx.useAltHpelIf      [cnt] = !ibcFlag && miNeighbor.useAltHpelIf;
#endif
      mrgCtx.BcwIdx            [cnt] = (miNeighbor.interDir == 3) ? miNeighbor.BcwIdx : BCW_DEFAULT;
#if INTER_LIC
      mrgCtx.LICFlags          [cnt] = miNeighbor.usesLIC;
#if !JVET_Z0075_IBC_HMVP_ENLARGE
      if (ibcFlag)
      {
        CHECK(mrgCtx.LICFlags[cnt], "addMergeHMVPCand: LIC is not used with IBC mode")
      }
#endif
#endif

      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
      if (slice.isInterB())
      {
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);
      }
#if MULTI_HYP_PRED
      mrgCtx.addHypNeighbours[cnt] = miNeighbor.addHypData;
#endif
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
      {
        mrgCtx.interDirNeighbours[cnt] = pu.amvpMergeModeFlag[0] ? 1 : 2;
        mrgCtx.mvFieldNeighbours[(cnt << 1) + (pu.amvpMergeModeFlag[0] ? 1 : 0)].setMvField(Mv(), -1);
      }
#endif

#if JVET_Z0084_IBC_TM && !JVET_Z0075_IBC_HMVP_ENLARGE
      if (ibcFlag)
      {
        if (mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh))
        {
          continue;
        }
      }
#if NON_ADJACENT_MRG_CAND
      else
#endif
#endif
#if NON_ADJACENT_MRG_CAND
      if (mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                            , mvdSimilarityThresh
#endif
      ))
      {
        continue;
      }
#endif
#if JVET_Y0058_IBC_LIST_MODIFY && !JVET_Z0075_IBC_HMVP_ENLARGE && !JVET_Z0084_IBC_TM
      if (ibcFlag)
      {
        if (!checkIsIBCCandidateValid(pu, miNeighbor))
        {
          continue;
        }
      }
#endif
      if (mrgCandIdx == cnt)
      {
        return true;
      }
      cnt ++;

      if (cnt  == maxNumMergeCandMin1)
      {
        break;
      }
#if JVET_Y0058_IBC_LIST_MODIFY && JVET_Z0084_IBC_TM && !JVET_Z0075_IBC_HMVP_ENLARGE
      }
#endif
    }
  }

  if (cnt < maxNumMergeCandMin1)
  {
    mrgCtx.useAltHpelIf[cnt] = false;
  }
  return false;
}

#if JVET_Z0075_IBC_HMVP_ENLARGE
bool PU::addIBCMergeHMVPCand(const CodingStructure &cs, MergeCtx &mrgCtx, const int &mrgCandIdx,
  const uint32_t maxNumMergeCandMin1, int &cnt
#if JVET_Y0058_IBC_LIST_MODIFY
  , const PredictionUnit &pu
#endif
#if TM_MRG || JVET_Z0084_IBC_TM
  , const uint32_t mvdSimilarityThresh
#endif
)
{
#if !JVET_Z0084_IBC_TM
  const Slice& slice = *cs.slice;
#endif
  MotionInfo miNeighbor;
#if JVET_Z0118_GDR
  bool isClean = cs.isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

#if JVET_Z0118_GDR
  auto &lut = (isClean) ? cs.motionLut.lutIbc1 : cs.motionLut.lutIbc0;
#else
  auto &lut = cs.motionLut.lutIbc;
#endif
  int num_avai_candInLUT = (int)lut.size();
  int compareNum = cnt;

  for (int mrgIdx = 1; mrgIdx <= num_avai_candInLUT; mrgIdx++)
  {
    miNeighbor = lut[num_avai_candInLUT - mrgIdx];

#if JVET_Y0058_IBC_LIST_MODIFY && JVET_Z0084_IBC_TM
    if (checkIsIBCCandidateValid(pu, miNeighbor))
    {
#endif
      mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
#if !JVET_Z0084_IBC_TM
      mrgCtx.useAltHpelIf      [cnt] = false;
      mrgCtx.BcwIdx            [cnt] = (miNeighbor.interDir == 3) ? miNeighbor.BcwIdx : BCW_DEFAULT;
#if INTER_LIC
      mrgCtx.LICFlags          [cnt] = miNeighbor.usesLIC;
      CHECK(mrgCtx.LICFlags[cnt], "addIBCMergeHMVPCand: LIC is not used with IBC mode")
#endif
#endif

      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
#if !JVET_Z0084_IBC_TM
      if (slice.isInterB())
      {
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);
      }
#if MULTI_HYP_PRED
      mrgCtx.addHypNeighbours[cnt] = miNeighbor.addHypData;
#endif
#endif

#if JVET_Z0084_IBC_TM
      if (mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh, compareNum))
      {
        continue;
      }
#else
      if (mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
        , mvdSimilarityThresh
#endif
        , compareNum
      ))
      {
        continue;
      }
#endif
#if JVET_Y0058_IBC_LIST_MODIFY && !JVET_Z0084_IBC_TM
      if (!checkIsIBCCandidateValid(pu, miNeighbor))
      {
        continue;
      }
#endif
      if (mrgCandIdx == cnt)
      {
        return true;
      }
      cnt++;

      if (cnt == maxNumMergeCandMin1)
      {
        break;
      }
#if JVET_Y0058_IBC_LIST_MODIFY && JVET_Z0084_IBC_TM
    }
#endif
  }

#if !JVET_Z0084_IBC_TM
  if (cnt < maxNumMergeCandMin1)
  {
    mrgCtx.useAltHpelIf[cnt] = false;
  }
#endif
  return false;
}
#endif
void PU::getIBCMergeCandidates(const PredictionUnit &pu, MergeCtx& mrgCtx, const int& mrgCandIdx)
{
  const CodingStructure &cs = *pu.cs;
#if JVET_Z0075_IBC_HMVP_ENLARGE
  const uint32_t maxNumMergeCand = IBC_MRG_MAX_NUM_CANDS_MEM;
#else
  const uint32_t maxNumMergeCand = pu.cs->sps->getMaxNumIBCMergeCand();
#endif
#if JVET_Z0084_IBC_TM
#if TM_MRG
  const uint32_t mvdSimilarityThresh = pu.tmMergeFlag ? PU::getTMMvdThreshold(pu) : 1;
#else
  const uint32_t mvdSimilarityThresh = 1;
#endif
#endif

  for (uint32_t ui = 0; ui < maxNumMergeCand; ++ui)
  {
    mrgCtx.BcwIdx[ui] = BCW_DEFAULT;
    mrgCtx.interDirNeighbours[ui] = 0;
    mrgCtx.mvFieldNeighbours[ui * 2].refIdx = NOT_VALID;
    mrgCtx.mvFieldNeighbours[ui * 2 + 1].refIdx = NOT_VALID;
    mrgCtx.useAltHpelIf[ui] = false;
#if INTER_LIC
    mrgCtx.LICFlags[ui] = false;
#endif
  }

#if JVET_Z0075_IBC_HMVP_ENLARGE
  mrgCtx.numValidMergeCand = pu.cs->sps->getMaxNumIBCMergeCand();
#else 
  mrgCtx.numValidMergeCand = maxNumMergeCand;
#endif
  // compute the location of the current PU

  int cnt = 0;

#if JVET_Y0058_IBC_LIST_MODIFY
  const Position posLT = pu.Y().topLeft();
#endif
  const Position posRT = pu.Y().topRight();
  const Position posLB = pu.Y().bottomLeft();

  MotionInfo miAbove, miLeft, miAboveLeft, miAboveRight, miBelowLeft;

  //left
  const PredictionUnit* puLeft = cs.getPURestricted(posLB.offset(-1, 0), pu, pu.chType);
  bool isGt4x4 = pu.lwidth() * pu.lheight() > 16;
  const bool isAvailableA1 = puLeft && pu.cu != puLeft->cu && CU::isIBC(*puLeft->cu);
  if (isGt4x4 && isAvailableA1)
  {
    miLeft = puLeft->getMotionInfo(posLB.offset(-1, 0));

#if JVET_Y0058_IBC_LIST_MODIFY
    if (checkIsIBCCandidateValid(pu, miLeft))
    {
#endif
    // get Inter Dir
    mrgCtx.interDirNeighbours[cnt] = miLeft.interDir;
    // get Mv from Left
    mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miLeft.mv[0], miLeft.refIdx[0]);
    if (mrgCandIdx == cnt)
    {
      return;
    }
    cnt++;

    // early termination
    if (cnt == maxNumMergeCand)
    {
      return;
    }
#if JVET_Y0058_IBC_LIST_MODIFY
    }
#endif
  }

  // above
  const PredictionUnit *puAbove = cs.getPURestricted(posRT.offset(0, -1), pu, pu.chType);
  bool isAvailableB1 = puAbove && pu.cu != puAbove->cu && CU::isIBC(*puAbove->cu);
  if (isGt4x4 && isAvailableB1)
  {
    miAbove = puAbove->getMotionInfo(posRT.offset(0, -1));

    if (!isAvailableA1 || (miAbove != miLeft))
    {
#if JVET_Y0058_IBC_LIST_MODIFY
      if (checkIsIBCCandidateValid(pu, miAbove))
      {
#endif
      // get Inter Dir
      mrgCtx.interDirNeighbours[cnt] = miAbove.interDir;
      // get Mv from Above
      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miAbove.mv[0], miAbove.refIdx[0]);

#if JVET_Z0084_IBC_TM
      if( !mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh) )
#endif
      {
        if (mrgCandIdx == cnt)
        {
          return;
        }
        cnt++;

        // early termination
        if (cnt == maxNumMergeCand)
        {
          return;
        }
      }
#if JVET_Y0058_IBC_LIST_MODIFY
      }
#endif
    }
  }

#if JVET_Y0058_IBC_LIST_MODIFY
  // above right
  const PredictionUnit *puAboveRight = cs.getPURestricted(posRT.offset(1, -1), pu, pu.chType);
  bool isAvailableB0 = puAboveRight && pu.cu != puAboveRight->cu && CU::isIBC(*puAboveRight->cu);
  if (isGt4x4 && isAvailableB0)
  {
    miAboveRight = puAboveRight->getMotionInfo(posRT.offset(1, -1));

    if ((!isAvailableB1 || (miAbove != miAboveRight)) && (!isAvailableA1 || (miLeft != miAboveRight)))
    {
      if (checkIsIBCCandidateValid(pu, miAboveRight))
      {
        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miAboveRight.interDir;
        // get Mv from Above-right
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miAboveRight.mv[0], miAboveRight.refIdx[0]);

#if JVET_Z0084_IBC_TM
        if( !mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh) )
#endif
        {
          if (mrgCandIdx == cnt)
          {
            return;
          }
          cnt++;

          // early termination
          if (cnt == maxNumMergeCand)
          {
            return;
          }
        }
      }
    }
  }

  //left bottom
  const PredictionUnit *puLeftBottom = cs.getPURestricted(posLB.offset(-1, 1), pu, pu.chType);
  bool isAvailableA0 = puLeftBottom && pu.cu != puLeftBottom->cu && CU::isIBC(*puLeftBottom->cu);
  if (isGt4x4 && isAvailableA0)
  {
    miBelowLeft = puLeftBottom->getMotionInfo(posLB.offset(-1, 1));

    if ((!isAvailableA1 || (miBelowLeft != miLeft)) && (!isAvailableB1 || (miBelowLeft != miAbove)) && (!isAvailableB0 || (miBelowLeft != miAboveRight)))
    {
      if (checkIsIBCCandidateValid(pu, miBelowLeft))
      {
        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miBelowLeft.interDir;
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miBelowLeft.mv[0], miBelowLeft.refIdx[0]);

#if JVET_Z0084_IBC_TM
        if( !mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh) )
#endif
        {
          if (mrgCandIdx == cnt)
          {
            return;
          }
          cnt++;

          // early termination
          if (cnt == maxNumMergeCand)
          {
            return;
          }
        }
      }
    }
  }

  // above left
#if JVET_Z0084_IBC_TM && JVET_Z0075_IBC_HMVP_ENLARGE && JVET_W0090_ARMC_TM
  if ((cnt < 4) && pu.cs->sps->getUseAML() && (mrgCandIdx >= 0)) //Only for AMVP case
#elif !JVET_Z0075_IBC_HMVP_ENLARGE
  if (cnt < 4)
#endif
  {
    const PredictionUnit *puAboveLeft = cs.getPURestricted(posLT.offset(-1, -1), pu, pu.chType);
    bool isAvailableB2 = puAboveLeft && pu.cu != puAboveLeft->cu && CU::isIBC(*puAboveLeft->cu);
    if (isGt4x4 && isAvailableB2)
    {
      miAboveLeft = puAboveLeft->getMotionInfo(posLT.offset(-1, -1));

      if ((!isAvailableA1 || (miLeft != miAboveLeft)) && (!isAvailableB1 || (miAbove != miAboveLeft)) && (!isAvailableA0 || (miBelowLeft != miAboveLeft)) && (!isAvailableB0 || (miAboveRight != miAboveLeft)))
      {
        if (checkIsIBCCandidateValid(pu, miAboveLeft))
        {
          // get Inter Dir
          mrgCtx.interDirNeighbours[cnt] = miAboveLeft.interDir;
          mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miAboveLeft.mv[0], miAboveLeft.refIdx[0]);

#if JVET_Z0084_IBC_TM
          if( !mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh) )
#endif
          {
            if (mrgCandIdx == cnt)
            {
              return;
            }
            cnt++;

            // early termination
            if (cnt == maxNumMergeCand)
            {
              return;
            }
          }
        }
      }
    }
  }
#endif

  if (cnt != maxNumMergeCand)
  {
#if JVET_Z0075_IBC_HMVP_ENLARGE
    bool bFound = addIBCMergeHMVPCand(cs, mrgCtx, mrgCandIdx, maxNumMergeCand, cnt
#if JVET_Y0058_IBC_LIST_MODIFY
                                     , pu
#endif
#if JVET_Z0084_IBC_TM
                                     , mvdSimilarityThresh
#endif
                                     );
#else
    bool bFound = addMergeHMVPCand(cs, mrgCtx, mrgCandIdx, maxNumMergeCand, cnt
#if JVET_Z0084_IBC_TM
                                  , isGt4x4 && isAvailableA1, miLeft, isGt4x4 && isAvailableB1, miAbove
#else
                                  , isAvailableA1, miLeft, isAvailableB1, miAbove
#endif
                                  , true, isGt4x4
#if JVET_X0083_BM_AMVP_MERGE_MODE || JVET_Y0058_IBC_LIST_MODIFY || JVET_Z0118_GDR
                                  , pu
#endif
#if JVET_Z0084_IBC_TM
#if JVET_X0083_BM_AMVP_MERGE_MODE
#if !JVET_Y0128_NON_CTC
                                  , 0, 0
#endif
#endif
                                  , mvdSimilarityThresh
#endif
                                  );
#endif

    if (bFound)
    {
      return;
    }

#if JVET_Z0075_IBC_HMVP_ENLARGE
    if (cnt == maxNumMergeCand)
    {
      return;
    }
#endif
  }

#if JVET_Y0058_IBC_LIST_MODIFY
  // pairwise-average candidates
  if (cnt>1 && cnt <maxNumMergeCand)
  {
    mrgCtx.mvFieldNeighbours[cnt * 2    ].setMvField(Mv(0, 0), NOT_VALID);
    mrgCtx.mvFieldNeighbours[cnt * 2 + 1].setMvField(Mv(0, 0), NOT_VALID);

    const Mv& MvI = mrgCtx.mvFieldNeighbours[0 * 2].mv;
    const Mv& MvJ = mrgCtx.mvFieldNeighbours[1 * 2].mv;
    // average two MVs
    Mv avgMv = MvI;

    avgMv += MvJ;
    roundAffineMv(avgMv.hor, avgMv.ver, 1);
    avgMv.roundToPrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
    MotionInfo miAvg;
    miAvg.mv[0] = avgMv;
    miAvg.refIdx[0] = MAX_NUM_REF;
    if (checkIsIBCCandidateValid(pu, miAvg))
    {
      mrgCtx.mvFieldNeighbours[cnt * 2].setMvField(avgMv, MAX_NUM_REF);
      mrgCtx.interDirNeighbours[cnt] = 1;
#if JVET_Z0084_IBC_TM
      if( !mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh) )
#else
      if (!mrgCtx.xCheckSimilarMotion(cnt))
#endif
      {
        if (mrgCandIdx == cnt)
        {
          return;
        }
        cnt++;

        // early termination
        if (cnt == maxNumMergeCand)
        {
          return;
        }
      }
    }
  }
#endif

#if JVET_Z0160_IBC_ZERO_PADDING
  int  validIdx = 0;
  Mv   validBvCand[IBC_MRG_MAX_NUM_CANDS];
  bool padList = true;

  if (cnt < maxNumMergeCand)
  {
    const int ctbSize     = pu.cs->sps->getCTUSize();
    const int log2CtbSize = floorLog2(ctbSize);
    const int cbX         = pu.lx() & (ctbSize - 1);
    const int cbY         = pu.ly() & (ctbSize - 1);
    const int cbWidth     = pu.lwidth();
    const int cbHeight    = pu.lheight();
#if !JVET_Z0153_IBC_EXT_REF
    const int leftEdge = -ctbSize
                         * (std::min(pu.Y().x >> log2CtbSize, (1 << ((MAX_CU_DEPTH - log2CtbSize) << 1))
                                                                - ((log2CtbSize < MAX_CU_DEPTH) ? 1 : 0)));
    int deltaX = cbX - leftEdge - cbWidth;
    int deltaY = cbY - cbHeight;
#else
    const int deltaX = ctbSize * (pu.lx() >> log2CtbSize) + cbX - cbWidth;
    const int deltaY = ((pu.ly() >> log2CtbSize) == 0)
                         ? (cbY - cbHeight)
                         : ((((pu.ly() >> log2CtbSize) == 1) ? 1 : 2) * ctbSize + cbY - cbHeight);
#endif

    if (deltaX >= 0 && deltaY >= 0)
    {
      validBvCand[0] = Mv(-cbWidth, -cbHeight);
      validBvCand[1] = Mv(-cbWidth, 0);
      validBvCand[2] = Mv(0, -cbHeight);
      validBvCand[3] = Mv(-cbWidth - (deltaX >> 1), -(cbHeight >> 1));
      validBvCand[4] = Mv(-(cbWidth >> 1), -cbHeight - (deltaY >> 1));
      validBvCand[5] = Mv(-cbWidth - (deltaX >> 1), -cbHeight - (deltaY >> 1));
    }
    else if (deltaY >= 0)
    {
      validBvCand[0] = Mv(-cbX, -cbHeight);
      validBvCand[1] = Mv(0, -cbHeight);
      validBvCand[2] = Mv(-(cbX >> 1), -cbHeight - (deltaY >> 1));
      validBvCand[3] = Mv(-cbX, -cbY);
      validBvCand[4] = Mv(0, -cbY);
      validBvCand[5] = Mv(cbWidth, -cbHeight - (deltaY >> 1));
    }
    else if (deltaX >= 0)
    {
      validBvCand[0] = Mv(-cbWidth, -cbY);
      validBvCand[1] = Mv(-cbWidth, 0);
      validBvCand[2] = Mv(-cbWidth - (deltaX >> 1), -(cbY >> 1));
      validBvCand[3] = Mv(-cbX, -cbY);
      validBvCand[4] = Mv(-cbX, 0);
      validBvCand[5] = Mv(-cbWidth - (deltaX >> 1), cbHeight);
    }
    else
    {
      padList = false;
    }
  }

  while (cnt < maxNumMergeCand && padList && validIdx < IBC_MRG_MAX_NUM_CANDS)
  {
#if !JVET_Z0084_IBC_TM
    bool duplicateBVP = false;
#endif
    Mv mvp = validBvCand[validIdx];
    mvp.changePrecision(MV_PRECISION_INT, MV_PRECISION_INTERNAL);

#if !JVET_Z0084_IBC_TM
    for (int i = 0; i < cnt; i++)
    {
      if (mvp == mrgCtx.mvFieldNeighbours[i * 2].mv)
      {
        duplicateBVP = true;
        break;
      }
    }
#endif
    MotionInfo miCand;
    miCand.mv[0] = mvp;
    mrgCtx.mvFieldNeighbours[cnt * 2].setMvField(mvp, MAX_NUM_REF);
    mrgCtx.interDirNeighbours[cnt] = 1;

#if JVET_Z0084_IBC_TM
    if (!mrgCtx.xCheckSimilarIBCMotion(cnt, mvdSimilarityThresh) && checkIsIBCCandidateValid(pu, miCand))
#else
    if (!duplicateBVP && checkIsIBCCandidateValid(pu, miCand))
#endif
    {
      if (mrgCandIdx == cnt)
      {
        return;
      }
      cnt++;
    }
    validIdx++;
  }

#else
  while (cnt < maxNumMergeCand)
  {
    mrgCtx.interDirNeighbours[cnt] = 1;
    mrgCtx.mvFieldNeighbours[cnt * 2].setMvField(Mv(0, 0), MAX_NUM_REF);
    if (mrgCandIdx == cnt)
    {
      return;
    }
    cnt++;
  }
#endif

#if JVET_Z0084_IBC_TM && JVET_Z0075_IBC_HMVP_ENLARGE
  mrgCtx.numValidMergeCand = std::min((int)pu.cs->sps->getMaxNumIBCMergeCand(), cnt);
#elif !JVET_Z0075_IBC_HMVP_ENLARGE
  mrgCtx.numValidMergeCand = cnt;
#endif

#if (JVET_Z0160_IBC_ZERO_PADDING || JVET_Z0084_IBC_TM) && JVET_Z0075_IBC_HMVP_ENLARGE
  // Add a zero motion to be able to stop reordering
  if (cnt < maxNumMergeCand)
  {
    mrgCtx.interDirNeighbours[cnt] = 1;
    mrgCtx.mvFieldNeighbours[cnt * 2].setMvField(Mv(0, 0), MAX_NUM_REF);
  }
#endif
}

#if  JVET_Y0058_IBC_LIST_MODIFY
bool PU::checkIsIBCCandidateValid(const PredictionUnit& pu, const MotionInfo miNeighbor)
{
  Mv bv = miNeighbor.mv[REF_PIC_LIST_0];
  bv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT); // used for only integer resolution
  const int cuPelX = pu.Y().x;
  const int cuPelY = pu.Y().y;
  int roiWidth = pu.lwidth();
  int roiHeight = pu.lheight();
  const int picWidth = pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
  const int picHeight = pu.cs->slice->getPPS()->getPicHeightInLumaSamples();
  const unsigned int  lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
  int xPred = bv.getHor();
  int yPred = bv.getVer();

  if (searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, xPred, yPred, lcuWidth)) // not valid bv derived
  {
    return true;
  }
  else
  {
    return false;
  }
}
#endif

#if JVET_Y0058_IBC_LIST_MODIFY || JVET_Z0084_IBC_TM
bool PU::searchBv(const PredictionUnit& pu, int xPos, int yPos, int width, int height, int picWidth, int picHeight, int xBv, int yBv, int ctuSize)
{
  const int ctuSizeLog2 = floorLog2(ctuSize);

  int refRightX = xPos + xBv + width - 1;
  int refLeftX  = xPos + xBv;

  int refBottomY = yPos + yBv + height - 1;
  int refTopY    = yPos + yBv;

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
  {
    return false;
  }
#endif

  // Don't search the below CTU row
  if (refBottomY >> ctuSizeLog2 > yPos >> ctuSizeLog2)
  {
    return false;
  }

#if JVET_Z0084_IBC_TM
  unsigned curTileIdx = pu.cs->pps->getTileIdx(Position(xPos, yPos));
#else
  unsigned curTileIdx = pu.cs->pps->getTileIdx(pu.lumaPos());
#endif
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
  }
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
#if JVET_Z0084_IBC_TM
      const Position refPosCol(refLeftX + ctuSize,  refTopY);
#else
      const Position refPosCol = pu.Y().topLeft().offset(xBv + ctuSize, yBv);
#endif
      int offset64x = (refPosCol.x >> (ctuSizeLog2 - 1)) << (ctuSizeLog2 - 1);
      int offset64y = (refPosCol.y >> (ctuSizeLog2 - 1)) << (ctuSizeLog2 - 1);
      const Position refPosCol64x64 = {offset64x, offset64y};
      if (pu.cs->isDecomp(refPosCol64x64, toChannelType(COMPONENT_Y)))
      {
        return false;
      }
      if (refPosCol64x64 == pu.Y().topLeft())
      {
        return false;
      }
    }
  }
  else
  {
    return false;
  }
#endif

  // in the same CTU, or valid area from left CTU. Check if the reference block is already coded
#if JVET_Z0084_IBC_TM
  const Position refPosLT(refLeftX,  refTopY);
  const Position refPosBR(refRightX, refBottomY);
#else
  const Position refPosLT = pu.Y().topLeft().offset(xBv, yBv);
  const Position refPosBR = pu.Y().bottomRight().offset(xBv, yBv);
#endif
  const ChannelType      chType = toChannelType(COMPONENT_Y);
  if (!pu.cs->isDecomp(refPosBR, chType))
  {
    return false;
  }
  if (!pu.cs->isDecomp(refPosLT, chType))
  {
    return false;
  }
  return true;
}
#endif

#if MULTI_PASS_DMVR || JVET_W0097_GPM_MMVD_TM
uint32_t PU::getBDMVRMvdThreshold(const PredictionUnit &pu)
{
  uint32_t numPixels = pu.lwidth() * pu.lheight();
  if (numPixels < 64)
  {
    return (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 2;
  }
  else if (numPixels < 256)
  {
    return (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 1;
  }
  else
  {
    return (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 0;
  }
}
#endif

#if TM_MRG
uint32_t PU::getTMMvdThreshold(const PredictionUnit &pu)
{
  uint32_t numPixels = pu.lwidth() * pu.lheight();
  if (numPixels < 64)
  {
    return 1 << MV_FRACTIONAL_BITS_INTERNAL;
  }
  else if (numPixels < 256)
  {
    return 2 << MV_FRACTIONAL_BITS_INTERNAL;
  }
  else
  {
    return 4 << MV_FRACTIONAL_BITS_INTERNAL;
  }
}

int PU::reorderInterMergeCandidates(const PredictionUnit &pu, MergeCtx& mrgCtx, int numCand, uint32_t mvdSimilarityThresh)
{
  Slice& slice = *pu.cu->slice;

  if (numCand > 0 && slice.isInterB() && slice.getRefPOC(REF_PIC_LIST_0, 0) != slice.getRefPOC(REF_PIC_LIST_1, 0))
  {
    int numCandUni = 0;
    int numCandBi  = 0;

    // Promote bi candidates and collect uni candidates
    MergeCtx mrgCtxUni;

    for (int mergeCand = 0; mergeCand < numCand; ++mergeCand)
    {
      if (mrgCtx.interDirNeighbours[mergeCand] == 3)
      {
        mrgCtx.copyRegularMergeCand(numCandBi, mrgCtx, mergeCand);
        ++numCandBi;
      }
      else
      {
        mrgCtxUni.copyRegularMergeCand(numCandUni, mrgCtx, mergeCand);
        mrgCtxUni.convertRegularMergeCandToBi(numCandUni);
        ++numCandUni;
      }
    }

    // Insert uni-to-bi candidates
    for (int uniCandIdx = 0; uniCandIdx < numCandUni; ++uniCandIdx)
    {
      int mergeCand = numCandBi;
      mrgCtx.copyRegularMergeCand(mergeCand, mrgCtxUni, uniCandIdx);

      if (!mrgCtx.xCheckSimilarMotion(mergeCand, mvdSimilarityThresh))
      {
        ++numCandBi;
      }
    }

    return numCandBi;
  }
  else
  {
    return numCand;
  }
}
#endif

void PU::getInterMergeCandidates( const PredictionUnit &pu, MergeCtx& mrgCtx,
                                 int mmvdList,
                                 const int& mrgCandIdx
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
                                , MergeCtx* tmvpMrgCtx
                                , MergeCtx* namvpMrgCtx
#endif
)
{
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
  const CodingStructure &cs  = *pu.cs;
  const Slice &slice         = *pu.cs->slice;
#if TM_MRG
  const uint32_t mvdSimilarityThresh = pu.tmMergeFlag ? PU::getTMMvdThreshold(pu) : 1;
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
  const uint32_t maxNumMergeCand = pu.tmMergeFlag ? (pu.cs->sps->getUseAML()) ? TM_MRG_MAX_NUM_INIT_CANDS : pu.cs->sps->getMaxNumTMMergeCand() : pu.cs->sps->getMaxNumMergeCand();
#else
  const uint32_t maxNumMergeCand     = pu.tmMergeFlag ? pu.cs->sps->getMaxNumTMMergeCand() : pu.cs->sps->getMaxNumMergeCand();
#endif
#else
  const uint32_t maxNumMergeCand = pu.cs->sps->getMaxNumMergeCand();
#endif
  for (uint32_t ui = 0; ui < maxNumMergeCand; ++ui)
  {
    mrgCtx.BcwIdx[ui] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[ui] = false;
#endif
    mrgCtx.interDirNeighbours[ui] = 0;
    mrgCtx.mvFieldNeighbours[(ui << 1)    ].refIdx = NOT_VALID;
    mrgCtx.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtx.useAltHpelIf[ui] = false;
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[ui].clear();
#endif
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
    mrgCtx.candCost[ui] = MAX_UINT64;
#endif
  }

  mrgCtx.numValidMergeCand = maxNumMergeCand;
#if JVET_Z0102_NO_ARMC_FOR_ZERO_CAND
  mrgCtx.numCandToTestEnc = maxNumMergeCand;
#endif
  // compute the location of the current PU
#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
  int amvpMergeCtxMergeDir = -1;
  RefPicList amvpRefList = REF_PIC_LIST_X;
  bool useAmvpMergeMode = false;
#endif
#if JVET_Y0128_NON_CTC
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    mrgCtx.numValidMergeCand = (mrgCandIdx + 1) > maxNumMergeCand ? maxNumMergeCand : (mrgCandIdx + 1);
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
    useAmvpMergeMode = true;
    if (pu.amvpMergeModeFlag[0] == true)
    {
      amvpMergeCtxMergeDir = 1;
      amvpRefList = REF_PIC_LIST_1;
    }
    else
    {
      amvpMergeCtxMergeDir = 2;
      amvpRefList = REF_PIC_LIST_0;
    }
#endif
  }
#else
  const int curPoc = slice.getPOC();
#if !JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
  RefPicList amvpRefList = REF_PIC_LIST_X;
#endif
  RefPicList mergeRefList = REF_PIC_LIST_X;
  int amvpPoc = -1;
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    mergeRefList = pu.amvpMergeModeFlag[0] ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
    amvpRefList = RefPicList(1 - mergeRefList);
    amvpPoc = slice.getRefPOC(amvpRefList, pu.refIdx[amvpRefList]);
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
    useAmvpMergeMode = true;
    if (pu.amvpMergeModeFlag[0] == true)
    {
      amvpMergeCtxMergeDir = 1;
    }
    else
    {
      amvpMergeCtxMergeDir = 2;
    }
#endif
    mrgCtx.numValidMergeCand = (mrgCandIdx + 1) > maxNumMergeCand ? maxNumMergeCand : (mrgCandIdx + 1);
  }
#endif
#endif

  int cnt = 0;

  const Position posLT = pu.Y().topLeft();
  const Position posRT = pu.Y().topRight();
  const Position posLB = pu.Y().bottomLeft();
  MotionInfo miAbove, miLeft, miAboveLeft, miAboveRight, miBelowLeft;

  // above
  const PredictionUnit *puAbove = cs.getPURestricted(posRT.offset(0, -1), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
  bool isAvailableB1 = puAbove && isDiffMER(pu.lumaPos(), posRT.offset(0, -1), plevel) && pu.cu != puAbove->cu && CU::isInter(*puAbove->cu) && puAbove->getMotionInfo(posRT.offset(0, -1)).isInter;
#else
  bool isAvailableB1 = puAbove && isDiffMER(pu.lumaPos(), posRT.offset(0, -1), plevel) && pu.cu != puAbove->cu && CU::isInter(*puAbove->cu);
#endif

  if (isAvailableB1)
  {
    miAbove = puAbove->getMotionInfo(posRT.offset(0, -1));

#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
    bool isValidAmMode = checkIsValidMergeMvCand(pu, miAbove.refIdx);
#else
    bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, miAbove.refIdx);
#endif
    if (isValidAmMode)
    {
#endif
    // get Inter Dir
    mrgCtx.interDirNeighbours[cnt] = miAbove.interDir;
    mrgCtx.useAltHpelIf[cnt] = miAbove.useAltHpelIf;
    // get Mv from Above
    mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puAbove->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[cnt] = miAbove.usesLIC;
#endif
    mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miAbove.mv[0], miAbove.refIdx[0]);

    if (slice.isInterB())
    {
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miAbove.mv[1], miAbove.refIdx[1]);
#if MULTI_HYP_PRED
      mrgCtx.addHypNeighbours[cnt] = puAbove->addHypData;
#endif
    }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
    if (useAmvpMergeMode)
    {
      mrgCtx.interDirNeighbours[cnt] = amvpMergeCtxMergeDir;
      mrgCtx.mvFieldNeighbours[(cnt << 1) + amvpRefList].setMvField(Mv(), -1);
    }
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
    if( !mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                           , mvdSimilarityThresh
#endif
    ) )
    {
#endif
    if (mrgCandIdx == cnt)
    {
#if TM_MRG
      if( !pu.tmMergeFlag )
#endif
      return;
    }

    cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
    }
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE
    } // if (isValidAmMode)
#endif
  }

  // early termination
  if (cnt == maxNumMergeCand)
  {
#if TM_MRG
    if (!pu.tmMergeFlag)
#endif
    return;
  }

  //left
  const PredictionUnit* puLeft = cs.getPURestricted(posLB.offset(-1, 0), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
  const bool isAvailableA1 = puLeft && isDiffMER(pu.lumaPos(), posLB.offset(-1, 0), plevel) && pu.cu != puLeft->cu && CU::isInter(*puLeft->cu) && puLeft->getMotionInfo(posLB.offset(-1, 0)).isInter;
#else
  const bool isAvailableA1 = puLeft && isDiffMER(pu.lumaPos(), posLB.offset(-1, 0), plevel) && pu.cu != puLeft->cu && CU::isInter(*puLeft->cu);
#endif

  if (isAvailableA1)
  {
    miLeft = puLeft->getMotionInfo(posLB.offset(-1, 0));

#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
    bool isValidAmMode = checkIsValidMergeMvCand(pu, miLeft.refIdx);
#else
    bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, miLeft.refIdx);
#endif
    if (isValidAmMode && (!isAvailableB1 || (miAbove != miLeft)))
#else
    if (!isAvailableB1 || (miAbove != miLeft))
#endif
    {
      // get Inter Dir
      mrgCtx.interDirNeighbours[cnt] = miLeft.interDir;
      mrgCtx.useAltHpelIf[cnt] = miLeft.useAltHpelIf;
      mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puLeft->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = miLeft.usesLIC;
#endif
      // get Mv from Left
      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miLeft.mv[0], miLeft.refIdx[0]);

      if (slice.isInterB())
      {
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miLeft.mv[1], miLeft.refIdx[1]);
#if MULTI_HYP_PRED
        mrgCtx.addHypNeighbours[cnt] = puLeft->addHypData;
#endif
      }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      if (useAmvpMergeMode)
      {
        mrgCtx.interDirNeighbours[cnt] = amvpMergeCtxMergeDir;
        mrgCtx.mvFieldNeighbours[(cnt << 1) + amvpRefList].setMvField(Mv(), -1);
      }
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
      if( !mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                             , mvdSimilarityThresh
#endif
      ) )
      {
#endif
      if (mrgCandIdx == cnt)
      {
#if TM_MRG
        if (!pu.tmMergeFlag)
#endif
        return;
      }

      cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
      }
#endif
    }
  }

  // early termination
  if( cnt == maxNumMergeCand )
  {
#if TM_MRG
    if (!pu.tmMergeFlag)
#endif
    return;
  }

  // above right
  const PredictionUnit *puAboveRight = cs.getPURestricted( posRT.offset( 1, -1 ), pu, pu.chType );

#if JVET_Y0065_GPM_INTRA
  bool isAvailableB0 = puAboveRight && isDiffMER( pu.lumaPos(), posRT.offset(1, -1), plevel) && CU::isInter( *puAboveRight->cu ) && puAboveRight->getMotionInfo( posRT.offset( 1, -1 ) ).isInter;
#else
  bool isAvailableB0 = puAboveRight && isDiffMER( pu.lumaPos(), posRT.offset(1, -1), plevel) && CU::isInter( *puAboveRight->cu );
#endif

  if( isAvailableB0 )
  {
    miAboveRight = puAboveRight->getMotionInfo( posRT.offset( 1, -1 ) );

#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
    bool isValidAmMode = checkIsValidMergeMvCand(pu, miAboveRight.refIdx);
#else
    bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, miAboveRight.refIdx);
#endif
    if (isValidAmMode && (!isAvailableB1 || ( miAbove != miAboveRight )))
#else
    if( !isAvailableB1 || ( miAbove != miAboveRight ) )
#endif
    {

      // get Inter Dir
      mrgCtx.interDirNeighbours[cnt] = miAboveRight.interDir;
      mrgCtx.useAltHpelIf[cnt] = miAboveRight.useAltHpelIf;
      // get Mv from Above-right
      mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puAboveRight->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = miAboveRight.usesLIC;
#endif
      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField( miAboveRight.mv[0], miAboveRight.refIdx[0] );

      if( slice.isInterB() )
      {
        mrgCtx.mvFieldNeighbours[( cnt << 1 ) + 1].setMvField( miAboveRight.mv[1], miAboveRight.refIdx[1] );
#if MULTI_HYP_PRED
        mrgCtx.addHypNeighbours[cnt] = puAboveRight->addHypData;
#endif
      }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      if (useAmvpMergeMode)
      {
        mrgCtx.interDirNeighbours[cnt] = amvpMergeCtxMergeDir;
        mrgCtx.mvFieldNeighbours[(cnt << 1) + amvpRefList].setMvField(Mv(), -1);
      }
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
      if( !mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                             , mvdSimilarityThresh
#endif
      ) )
      {
#endif
      if (mrgCandIdx == cnt)
      {
#if TM_MRG
        if (!pu.tmMergeFlag)
#endif
        return;
      }

      cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
      }
#endif
    }
  }
  // early termination
  if( cnt == maxNumMergeCand )
  {
#if TM_MRG
    if (!pu.tmMergeFlag)
#endif
    return;
  }

  //left bottom
  const PredictionUnit *puLeftBottom = cs.getPURestricted( posLB.offset( -1, 1 ), pu, pu.chType );

#if JVET_Y0065_GPM_INTRA
  bool isAvailableA0 = puLeftBottom && isDiffMER( pu.lumaPos(), posLB.offset(-1, 1), plevel) && CU::isInter( *puLeftBottom->cu ) && puLeftBottom->getMotionInfo( posLB.offset( -1, 1 ) ).isInter;
#else
  bool isAvailableA0 = puLeftBottom && isDiffMER( pu.lumaPos(), posLB.offset(-1, 1), plevel) && CU::isInter( *puLeftBottom->cu );
#endif

  if( isAvailableA0 )
  {
    miBelowLeft = puLeftBottom->getMotionInfo( posLB.offset( -1, 1 ) );

#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
    bool isValidAmMode = checkIsValidMergeMvCand(pu, miBelowLeft.refIdx);
#else
    bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, miBelowLeft.refIdx);
#endif
    if (isValidAmMode && (!isAvailableA1 || ( miBelowLeft != miLeft )))
#else
    if( !isAvailableA1 || ( miBelowLeft != miLeft ) )
#endif
    {
      // get Inter Dir
      mrgCtx.interDirNeighbours[cnt] = miBelowLeft.interDir;
      mrgCtx.useAltHpelIf[cnt] = miBelowLeft.useAltHpelIf;
      mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puLeftBottom->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = miBelowLeft.usesLIC;
#endif
      // get Mv from Bottom-Left
      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField( miBelowLeft.mv[0], miBelowLeft.refIdx[0] );

      if( slice.isInterB() )
      {
        mrgCtx.mvFieldNeighbours[( cnt << 1 ) + 1].setMvField( miBelowLeft.mv[1], miBelowLeft.refIdx[1] );
#if MULTI_HYP_PRED
        mrgCtx.addHypNeighbours[cnt] = puLeftBottom->addHypData;
#endif
      }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      if (useAmvpMergeMode)
      {
        mrgCtx.interDirNeighbours[cnt] = amvpMergeCtxMergeDir;
        mrgCtx.mvFieldNeighbours[(cnt << 1) + amvpRefList].setMvField(Mv(), -1);
      }
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
      if( !mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                             , mvdSimilarityThresh
#endif
      ) )
      {
#endif
      if (mrgCandIdx == cnt)
      {
#if TM_MRG
        if (!pu.tmMergeFlag)
#endif
        return;
      }

      cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
      }
#endif
    }
  }
  // early termination
  if( cnt == maxNumMergeCand )
  {
#if TM_MRG
    if (!pu.tmMergeFlag)
#endif
    return;
  }

#if TM_MRG
  if (pu.tmMergeFlag)
  {
    cnt = PU::reorderInterMergeCandidates(pu, mrgCtx, cnt, mvdSimilarityThresh);

    if ((mrgCandIdx >= 0 && mrgCandIdx < cnt) || cnt >= maxNumMergeCand)
    {
      return;
    }
  }
#endif

  // above left
  if ( cnt < 4 )
  {
    const PredictionUnit *puAboveLeft = cs.getPURestricted( posLT.offset( -1, -1 ), pu, pu.chType );

#if JVET_Y0065_GPM_INTRA
    bool isAvailableB2 = puAboveLeft && isDiffMER( pu.lumaPos(), posLT.offset(-1, -1), plevel ) && CU::isInter( *puAboveLeft->cu ) && puAboveLeft->getMotionInfo( posLT.offset( -1, -1 ) ).isInter;
#else
    bool isAvailableB2 = puAboveLeft && isDiffMER( pu.lumaPos(), posLT.offset(-1, -1), plevel ) && CU::isInter( *puAboveLeft->cu );
#endif

    if( isAvailableB2 )
    {
      miAboveLeft = puAboveLeft->getMotionInfo( posLT.offset( -1, -1 ) );

#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
      bool isValidAmMode = checkIsValidMergeMvCand(pu, miAboveLeft.refIdx);
#else
      bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, miAboveLeft.refIdx);
#endif
      if (isValidAmMode && ( !isAvailableA1 || ( miLeft != miAboveLeft ) ) && ( !isAvailableB1 || ( miAbove != miAboveLeft ) ))
#else
      if( ( !isAvailableA1 || ( miLeft != miAboveLeft ) ) && ( !isAvailableB1 || ( miAbove != miAboveLeft ) ) )
#endif
      {

        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miAboveLeft.interDir;
        mrgCtx.useAltHpelIf[cnt] = miAboveLeft.useAltHpelIf;
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puAboveLeft->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
        mrgCtx.LICFlags[cnt] = miAboveLeft.usesLIC;
#endif
        // get Mv from Above-Left
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField( miAboveLeft.mv[0], miAboveLeft.refIdx[0] );

        if( slice.isInterB() )
        {
          mrgCtx.mvFieldNeighbours[( cnt << 1 ) + 1].setMvField( miAboveLeft.mv[1], miAboveLeft.refIdx[1] );
#if MULTI_HYP_PRED
          mrgCtx.addHypNeighbours[cnt] = puAboveLeft->addHypData;
#endif
        }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
        if (useAmvpMergeMode)
        {
          mrgCtx.interDirNeighbours[cnt] = amvpMergeCtxMergeDir;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + amvpRefList].setMvField(Mv(), -1);
        }
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
        if( !mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                               , mvdSimilarityThresh
#endif
        ) )
        {
#endif
        if (mrgCandIdx == cnt)
        {
          return;
        }

        cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
        }
#endif
      }
    }
  }
  // early termination
  if (cnt == maxNumMergeCand)
  {
    return;
  }

#if INTER_RM_SIZE_CONSTRAINTS
  if (slice.getPicHeader()->getEnableTMVPFlag())
#else
  if (slice.getPicHeader()->getEnableTMVPFlag() && (pu.lumaSize().width + pu.lumaSize().height > 12))
#endif
  {
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
    if (tmvpMrgCtx == NULL)
    {
#endif
    //>> MTK colocated-RightBottom
    // offset the pos to be sure to "point" to the same position the uiAbsPartIdx would've pointed to
    Position posRB = pu.Y().bottomRight().offset( -3, -3 );
    const PreCalcValues& pcv = *cs.pcv;

    Position posC0;
    Position posC1 = pu.Y().center();
    bool C0Avail = false;
    bool boundaryCond = ((posRB.x + pcv.minCUWidth) < pcv.lumaWidth) && ((posRB.y + pcv.minCUHeight) < pcv.lumaHeight);
    const SubPic& curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
    if (curSubPic.getTreatedAsPicFlag())
    {
      boundaryCond = ((posRB.x + pcv.minCUWidth) <= curSubPic.getSubPicRight() &&
                      (posRB.y + pcv.minCUHeight) <= curSubPic.getSubPicBottom());
    }
    if (boundaryCond)
    {
      int posYInCtu = posRB.y & pcv.maxCUHeightMask;
      if (posYInCtu + 4 < pcv.maxCUHeight)
      {
        posC0 = posRB.offset(4, 4);
        C0Avail = true;
      }
    }

    Mv        cColMv;
    int       iRefIdx     = 0;
    int       dir         = 0;
    unsigned  uiArrayAddr = cnt;
    bool      bExistMV    = ( C0Avail && getColocatedMVP(pu, REF_PIC_LIST_0, posC0, cColMv, iRefIdx, false
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
      , &iRefIdx
#endif
    ) )
                              || getColocatedMVP( pu, REF_PIC_LIST_0, posC1, cColMv, iRefIdx, false
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
                                , &iRefIdx
#endif
                              );
    if (bExistMV)
    {
      dir     |= 1;
      mrgCtx.mvFieldNeighbours[2 * uiArrayAddr].setMvField(cColMv, iRefIdx);
    }

    if (slice.isInterB())
    {
      bExistMV = ( C0Avail && getColocatedMVP(pu, REF_PIC_LIST_1, posC0, cColMv, iRefIdx, false
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
        , &iRefIdx
#endif
      ) )
                   || getColocatedMVP( pu, REF_PIC_LIST_1, posC1, cColMv, iRefIdx, false
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
                     , &iRefIdx
#endif
                   );
      if (bExistMV)
      {
        dir     |= 2;
        mrgCtx.mvFieldNeighbours[2 * uiArrayAddr + 1].setMvField(cColMv, iRefIdx);
      }
    }

#if JVET_X0083_BM_AMVP_MERGE_MODE
    int8_t tempRefIdx[2] = { mrgCtx.mvFieldNeighbours[2 * uiArrayAddr].refIdx, mrgCtx.mvFieldNeighbours[2 * uiArrayAddr + 1].refIdx };
#if JVET_Y0128_NON_CTC
    bool isValidAmMode = checkIsValidMergeMvCand(pu, tempRefIdx);
#else
    bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, tempRefIdx);
#endif
    if (isValidAmMode && ( dir != 0 ))
#else
    if( dir != 0 )
#endif
    {
      bool addTMvp = true;
      if( addTMvp )
      {
        mrgCtx.interDirNeighbours[uiArrayAddr] = dir;
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
        if (useAmvpMergeMode)
        {
          mrgCtx.interDirNeighbours[cnt] = amvpMergeCtxMergeDir;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + amvpRefList].setMvField(Mv(), -1);
        }
#endif
        mrgCtx.BcwIdx[uiArrayAddr] = BCW_DEFAULT;
#if INTER_LIC
        mrgCtx.LICFlags[uiArrayAddr] = false;
#endif
        mrgCtx.useAltHpelIf[uiArrayAddr] = false;
#if MULTI_HYP_PRED
        mrgCtx.addHypNeighbours[uiArrayAddr].clear();
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
        if( !mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                               , mvdSimilarityThresh
#endif
        ) )
        {
#endif
        if (mrgCandIdx == cnt)
        {
          return;
        }

        cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
        }
#endif
      }
    }
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
    }
    else if (tmvpMrgCtx->numValidMergeCand > 0)
    {
      mrgCtx.BcwIdx[cnt] = tmvpMrgCtx->BcwIdx[0];
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = tmvpMrgCtx->LICFlags[0];
#endif
      mrgCtx.interDirNeighbours[cnt] = tmvpMrgCtx->interDirNeighbours[0];
      mrgCtx.mvFieldNeighbours[cnt << 1] = tmvpMrgCtx->mvFieldNeighbours[0];
      mrgCtx.useAltHpelIf[cnt] = tmvpMrgCtx->useAltHpelIf[0];
      if (slice.isInterB())
      {
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1] = tmvpMrgCtx->mvFieldNeighbours[1];
#if MULTI_HYP_PRED
        mrgCtx.addHypNeighbours[cnt] = tmvpMrgCtx->addHypNeighbours[0];
#endif
      }

#if NON_ADJACENT_MRG_CAND || TM_MRG
      if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
        , mvdSimilarityThresh
#endif
      ))
      {
#endif
        mrgCtx.candCost[cnt] = tmvpMrgCtx->candCost[0];
        if (mrgCandIdx == cnt)
        {
          return;
        }

        cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
      }
#endif
    }
#endif
  }

  // early termination
  if (cnt == maxNumMergeCand)
  {
    return;
  }

#if NON_ADJACENT_MRG_CAND
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
  if (namvpMrgCtx == NULL)
  {
#endif
  MotionInfo miNeighbor;
  int offsetX = 0;
  int offsetY = 0;
  const int iNACANDIDATE_NUM[4] = { 3, 5, 5, 5 };
  const int idxMap[4][5] = { { 0, 1, 4 },{ 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4 } };

  for (int iDistanceIndex = 0; iDistanceIndex < NADISTANCE_LEVEL && cnt < maxNumMergeCand - 1; iDistanceIndex++)
  {
    const int iNADistanceHor = pu.Y().width  * (iDistanceIndex + 1);
    const int iNADistanceVer = pu.Y().height * (iDistanceIndex + 1);

    for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM[iDistanceIndex] && cnt < maxNumMergeCand - 1; NASPIdx++)
    {
      switch (idxMap[iDistanceIndex][NASPIdx])
      {
      case 0:offsetX = -iNADistanceHor - 1; offsetY = pu.Y().height + iNADistanceVer - 1; break;
      case 1:offsetX = pu.Y().width + iNADistanceHor - 1; offsetY = -iNADistanceVer - 1; break;
      case 2:offsetX = pu.Y().width >> 1;   offsetY = -iNADistanceVer - 1;    break;
      case 3:offsetX = -iNADistanceHor - 1; offsetY = pu.Y().height >> 1; break;
      case 4:offsetX = -iNADistanceHor - 1; offsetY = -iNADistanceVer - 1;    break;
      default: printf("error!"); exit(0); break;
      }

      const PredictionUnit *puNonAdjacent = cs.getPURestricted(posLT.offset(offsetX, offsetY), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER( pu.lumaPos(), posLT.offset( offsetX, offsetY ) , plevel) && CU::isInter(*puNonAdjacent->cu) && puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY)).isInter;
#else
      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER( pu.lumaPos(), posLT.offset( offsetX, offsetY ) , plevel) && CU::isInter(*puNonAdjacent->cu);
#endif

      if (isAvailableNonAdjacent)
      {
        miNeighbor = puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY));


#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
        bool isValidAmMode = checkIsValidMergeMvCand(pu, miNeighbor.refIdx);
#else
        bool isValidAmMode = checkIsValidMergeMvCand(cs, pu, curPoc, amvpPoc, miNeighbor.refIdx);
#endif
        if (isValidAmMode)
        {
#endif
        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
        // get Mv from Above-Left
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
        mrgCtx.useAltHpelIf[cnt] = miNeighbor.useAltHpelIf;
        // get Mv from Above-right
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puNonAdjacent->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
        mrgCtx.LICFlags[cnt] = miNeighbor.usesLIC;
#endif
        if (slice.isInterB())
        {
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);
#if MULTI_HYP_PRED
          mrgCtx.addHypNeighbours[cnt] = puNonAdjacent->addHypData;
#endif
        }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
        if (useAmvpMergeMode)
        {
          mrgCtx.interDirNeighbours[cnt] = amvpMergeCtxMergeDir;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + amvpRefList].setMvField(Mv(), -1);
        }
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                               , mvdSimilarityThresh
#endif
        ))
#endif
        {
          if (mrgCandIdx == cnt)
          {
            return;
          }
          cnt++;
        }
#if JVET_X0083_BM_AMVP_MERGE_MODE
        } // if (isValidAmMode)
#endif
      }

    }
  }
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
  }
  else
  {
    for (uint32_t ui = 0; ui < namvpMrgCtx->numValidMergeCand && cnt < maxNumMergeCand - 1; ++ui)
    {
      mrgCtx.BcwIdx[cnt] = namvpMrgCtx->BcwIdx[ui];
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = namvpMrgCtx->LICFlags[ui];
#endif
      mrgCtx.interDirNeighbours[cnt] = namvpMrgCtx->interDirNeighbours[ui];
      mrgCtx.mvFieldNeighbours[cnt << 1] = namvpMrgCtx->mvFieldNeighbours[ui << 1];
      mrgCtx.useAltHpelIf[cnt] = namvpMrgCtx->useAltHpelIf[ui];
      if (slice.isInterB())
      {
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1] = namvpMrgCtx->mvFieldNeighbours[(ui << 1) + 1];
#if MULTI_HYP_PRED
        mrgCtx.addHypNeighbours[cnt] = namvpMrgCtx->addHypNeighbours[ui];
#endif
      }

#if NON_ADJACENT_MRG_CAND || TM_MRG
      if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
        , mvdSimilarityThresh
#endif
      ))
      {
#endif
        mrgCtx.candCost[cnt] = namvpMrgCtx->candCost[ui];
        if (mrgCandIdx == cnt)
        {
          return;
        }

        cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
      }
#endif
    }
  }
#endif
#endif

  int maxNumMergeCandMin1 = maxNumMergeCand - 1;
  if (cnt != maxNumMergeCandMin1)
  {
#if !JVET_Z0075_IBC_HMVP_ENLARGE
    bool isGt4x4 = true;
#endif
    bool bFound = addMergeHMVPCand(cs, mrgCtx, mrgCandIdx, maxNumMergeCandMin1, cnt
      , isAvailableA1, miLeft, isAvailableB1, miAbove
#if !JVET_Z0075_IBC_HMVP_ENLARGE
      , CU::isIBC(*pu.cu)
      , isGt4x4
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE || (JVET_Y0058_IBC_LIST_MODIFY && !JVET_Z0075_IBC_HMVP_ENLARGE) || JVET_Z0118_GDR
      , pu
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE
#if !JVET_Y0128_NON_CTC
      , curPoc
      , amvpPoc
#endif
#endif
#if TM_MRG
      , mvdSimilarityThresh
#endif
      );

    if (bFound)
    {
      return;
    }
  }

  // pairwise-average candidates
  {
    if (cnt > 1 && cnt < maxNumMergeCand)
    {
      mrgCtx.mvFieldNeighbours[cnt * 2].setMvField( Mv( 0, 0 ), NOT_VALID );
      mrgCtx.mvFieldNeighbours[cnt * 2 + 1].setMvField( Mv( 0, 0 ), NOT_VALID );
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = false;
#endif
#if MULTI_HYP_PRED
      mrgCtx.addHypNeighbours[cnt].clear();
#endif
      // calculate average MV for L0 and L1 seperately
      unsigned char interDir = 0;
#if INTER_LIC
      bool averageUsed = false;
#endif

      mrgCtx.useAltHpelIf[cnt] = (mrgCtx.useAltHpelIf[0] == mrgCtx.useAltHpelIf[1]) ? mrgCtx.useAltHpelIf[0] : false;
      for( int refListId = 0; refListId < (slice.isInterB() ? 2 : 1); refListId++ )
      {
        const short refIdxI = mrgCtx.mvFieldNeighbours[0 * 2 + refListId].refIdx;
        const short refIdxJ = mrgCtx.mvFieldNeighbours[1 * 2 + refListId].refIdx;

        // both MVs are invalid, skip
        if( (refIdxI == NOT_VALID) && (refIdxJ == NOT_VALID) )
        {
          continue;
        }

        interDir += 1 << refListId;
        // both MVs are valid, average these two MVs
        if( (refIdxI != NOT_VALID) && (refIdxJ != NOT_VALID) )
        {
          const Mv& MvI = mrgCtx.mvFieldNeighbours[0 * 2 + refListId].mv;
          const Mv& MvJ = mrgCtx.mvFieldNeighbours[1 * 2 + refListId].mv;

          // average two MVs
          Mv avgMv = MvI;
          avgMv += MvJ;
          roundAffineMv(avgMv.hor, avgMv.ver, 1);

          mrgCtx.mvFieldNeighbours[cnt * 2 + refListId].setMvField( avgMv, refIdxI );
#if INTER_LIC
          mrgCtx.LICFlags[cnt] = false;
          averageUsed = true;
#endif
        }
        // only one MV is valid, take the only one MV
        else if( refIdxI != NOT_VALID )
        {
          Mv singleMv = mrgCtx.mvFieldNeighbours[0 * 2 + refListId].mv;
          mrgCtx.mvFieldNeighbours[cnt * 2 + refListId].setMvField( singleMv, refIdxI );
#if INTER_LIC
          if (!averageUsed)
          {
            mrgCtx.LICFlags[cnt] |= mrgCtx.LICFlags[0];
          }
#endif
        }
        else if( refIdxJ != NOT_VALID )
        {
          Mv singleMv = mrgCtx.mvFieldNeighbours[1 * 2 + refListId].mv;
          mrgCtx.mvFieldNeighbours[cnt * 2 + refListId].setMvField( singleMv, refIdxJ );
#if INTER_LIC
          if (!averageUsed)
          {
            mrgCtx.LICFlags[cnt] |= mrgCtx.LICFlags[1];
          }
#endif
        }
      }

      mrgCtx.interDirNeighbours[cnt] = interDir;
      if( interDir > 0 )
      {
#if INTER_LIC
        if (interDir == 3)
        {
          mrgCtx.LICFlags[cnt] = false;
        }
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
        if( !mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
                               , mvdSimilarityThresh
#endif
        ) )
#endif
        cnt++;
      }
    }

    // early termination
    if( cnt == maxNumMergeCand )
    {
      return;
    }
  }

#if JVET_Z0139_HIST_AFF
  if (addMergeHMVPCandFromAffModel(pu, mrgCtx, mrgCandIdx, cnt
#if TM_MRG
    , mvdSimilarityThresh
#endif
  ))
  {
    return;
  }
#endif

#if JVET_Z0139_HIST_AFF
  int uiArrayAddr = cnt;
#else
  uint32_t uiArrayAddr = cnt;
#endif

#if JVET_Z0102_NO_ARMC_FOR_ZERO_CAND 
  mrgCtx.numCandToTestEnc = cnt;
#endif
  int iNumRefIdx = slice.isInterB() ? std::min(slice.getNumRefIdx(REF_PIC_LIST_0), slice.getNumRefIdx(REF_PIC_LIST_1)) : slice.getNumRefIdx(REF_PIC_LIST_0);
#if JVET_X0083_BM_AMVP_MERGE_MODE
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    int mergeRefIdx = 0;
    const int targetMaxNumMergeCand = mrgCtx.numValidMergeCand;
    mrgCtx.numValidMergeCand = uiArrayAddr;
    while ((mrgCtx.numValidMergeCand < targetMaxNumMergeCand) && (mergeRefIdx < iNumRefIdx))
    {
#if JVET_Y0128_NON_CTC
      if (pu.cu->slice->getAmvpMergeModeValidCandPair(pu.amvpMergeModeFlag[0] ? mergeRefIdx  : pu.refIdx[0],
                                                      pu.amvpMergeModeFlag[0] ? pu.refIdx[1] : mergeRefIdx) == false)
#else
      const int mergePoc = slice.getRefPOC(mergeRefList, mergeRefIdx);
      if ((amvpPoc - curPoc) * (mergePoc - curPoc) > 0)
#endif
      {
        mergeRefIdx++;
        continue;
      }
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      mrgCtx.interDirNeighbours [ mrgCtx.numValidMergeCand ] = amvpMergeCtxMergeDir;
#else
      mrgCtx.interDirNeighbours [ mrgCtx.numValidMergeCand ] = 3;
#endif
      mrgCtx.BcwIdx             [ mrgCtx.numValidMergeCand ] = BCW_DEFAULT;
#if INTER_LIC
      mrgCtx.LICFlags           [ mrgCtx.numValidMergeCand ] = false;
#endif
      mrgCtx.useAltHpelIf       [ mrgCtx.numValidMergeCand ] = false;
#if MULTI_HYP_PRED
      mrgCtx.addHypNeighbours   [ mrgCtx.numValidMergeCand ].clear();
#endif
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
      mrgCtx.mvFieldNeighbours  [(mrgCtx.numValidMergeCand << 1) + (1 - amvpRefList)].setMvField(Mv(0, 0), mergeRefIdx);
      mrgCtx.mvFieldNeighbours  [(mrgCtx.numValidMergeCand << 1) + amvpRefList].setMvField(Mv(0, 0), -1);
#else
      mrgCtx.mvFieldNeighbours  [ mrgCtx.numValidMergeCand << 1].setMvField(Mv(0, 0), mergeRefIdx);
      mrgCtx.mvFieldNeighbours  [(mrgCtx.numValidMergeCand << 1) + 1].setMvField(Mv(0, 0), mergeRefIdx);
#endif
      mrgCtx.numValidMergeCand++;
      mergeRefIdx++;
    }
    return;
  }
#endif

  int r = 0;
  int refcnt = 0;
  while (uiArrayAddr < maxNumMergeCand)
  {
    mrgCtx.interDirNeighbours [uiArrayAddr     ] = 1;
    mrgCtx.BcwIdx             [uiArrayAddr     ] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags           [uiArrayAddr     ] = false;
#endif
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[uiArrayAddr].clear();
#endif
    mrgCtx.mvFieldNeighbours  [uiArrayAddr << 1].setMvField(Mv(0, 0), r);
    mrgCtx.useAltHpelIf[uiArrayAddr] = false;

    if (slice.isInterB())
    {
      mrgCtx.interDirNeighbours [ uiArrayAddr          ] = 3;
      mrgCtx.mvFieldNeighbours  [(uiArrayAddr << 1) + 1].setMvField(Mv(0, 0), r);
    }

    uiArrayAddr++;

    if (refcnt == iNumRefIdx - 1)
    {
      r = 0;
    }
    else
    {
      ++r;
      ++refcnt;
    }

#if JVET_Z0139_HIST_AFF
    if (uiArrayAddr < maxNumMergeCand && uiArrayAddr == cnt + 1)
    {
      if (addMergeHMVPCandFromAffModel(pu, mrgCtx, mrgCandIdx, uiArrayAddr
#if TM_MRG
        , mvdSimilarityThresh
#endif
      ))
      {
        return;
      }
      if (uiArrayAddr >= mrgCtx.numValidMergeCand)
      {
        return;
      }
    }
#endif
  }
#if NON_ADJACENT_MRG_CAND
  CHECK(mrgCtx.numValidMergeCand != uiArrayAddr, "not enough number of merge candidates!");
#else
  mrgCtx.numValidMergeCand = uiArrayAddr;
#endif
}

#if JVET_Y0128_NON_CTC
bool  PU::isBiRefScaled(const CodingStructure& cs, const int refIdx0, const int refIdx1)
{
  const bool isResamplingPossible = cs.sps->getRprEnabledFlag();
  if (!isResamplingPossible)
  {
    return false;
  }

  const bool ref0IsScaled = refIdx0 < 0 || refIdx0 >= MAX_NUM_REF
    ? false
    : isResamplingPossible && cs.slice->getRefPic(REF_PIC_LIST_0, refIdx0)->isRefScaled(cs.pps);
  const bool ref1IsScaled = refIdx1 < 0 || refIdx1 >= MAX_NUM_REF
    ? false
    : isResamplingPossible && cs.slice->getRefPic(REF_PIC_LIST_1, refIdx1)->isRefScaled(cs.pps);

  return (ref0IsScaled || ref1IsScaled);
}
#endif

#if JVET_X0049_ADAPT_DMVR
bool PU::isBMMergeFlagCoded(const PredictionUnit& pu)
{

  if (pu.cs->slice->getSPS()->getUseDMVDMode()
#if JVET_Y0128_NON_CTC
    && pu.cs->slice->getUseBM()
#else
    && !pu.cs->slice->getCheckLDC()
#endif
    )
  {
    return (pu.cs->sps->getMaxNumBMMergeCand() > 0);
  }
  return false;
}

bool PU::isBiPredFromDifferentDirEqDistPoc(const PredictionUnit& pu, int refIdx0, int refIdx1)
{
  if (refIdx0 >= 0 && refIdx1 >= 0)
  {
    if (pu.cu->slice->getRefPic(REF_PIC_LIST_0, refIdx0)->longTerm
      || pu.cu->slice->getRefPic(REF_PIC_LIST_1, refIdx1)->longTerm)
    {
      return false;
    }
#if JVET_Y0128_NON_CTC
    if ( PU::isBiRefScaled( *pu.cs, refIdx0, refIdx1 ) )
    {
      return false;
    }
#endif
    const int poc0 = pu.cu->slice->getRefPOC(REF_PIC_LIST_0, refIdx0);
    const int poc1 = pu.cu->slice->getRefPOC(REF_PIC_LIST_1, refIdx1);
    const int poc = pu.cu->slice->getPOC();
    if ((poc - poc0)*(poc - poc1) < 0)
    {
      if (abs(poc - poc0) == abs(poc - poc1))
      {
        return true;
      }
    }
  }
  return false;
}

bool PU::addBMMergeHMVPCand(const CodingStructure &cs, MergeCtx &mrgCtx, const int &mrgCandIdx,
  const uint32_t maxNumMergeCandMin1, int &cnt, const bool isAvailableA1,
  const MotionInfo miLeft, const bool isAvailableB1, const MotionInfo miAbove
#if !JVET_Z0075_IBC_HMVP_ENLARGE
  , const bool ibcFlag, const bool isGt4x4
#endif
#if TM_MRG
  , const uint32_t mvdSimilarityThresh
#endif
)
{
  const Slice& slice = *cs.slice;
  MotionInfo miNeighbor;
#if JVET_Z0118_GDR  
  bool isClean = cs.isClean(cs.area.Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

#if JVET_Z0075_IBC_HMVP_ENLARGE
#if JVET_Z0118_GDR  
  auto &lut = (isClean) ? cs.motionLut.lut1 : cs.motionLut.lut0;
#else
  auto &lut = cs.motionLut.lut;
#endif
#else
#if JVET_Z0118_GDR  
  auto &lut = ibcFlag ? (isClean ? cs.motionLut.lutIbc1 : cs.motionLut.lutIbc0) : (isClean ? cs.motionLut.lut1 : cs.motionLut.lut0);
#else
  auto &lut = ibcFlag ? cs.motionLut.lutIbc : cs.motionLut.lut;
#endif
#endif

  int num_avai_candInLUT = (int)lut.size();

  for (int mrgIdx = 1; mrgIdx <= num_avai_candInLUT; mrgIdx++)
  {
    miNeighbor = lut[num_avai_candInLUT - mrgIdx];

#if JVET_Z0075_IBC_HMVP_ENLARGE
    if ( mrgIdx > 2
#else
    if (mrgIdx > 2 || ((mrgIdx > 1 || !isGt4x4) && ibcFlag)
#endif
      || ((!isAvailableA1 || (miLeft != miNeighbor)) && (!isAvailableB1 || (miAbove != miNeighbor))))
    {
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = miNeighbor.usesLIC;

#if !JVET_Z0075_IBC_HMVP_ENLARGE
      if (ibcFlag)
      {
        CHECK(mrgCtx.LICFlags[cnt], "addMergeHMVPCand: LIC is not used with IBC mode")
      }
#endif
#endif

      bool isBmCand = false;
      int refIdx0 = miNeighbor.refIdx[0];
      int refIdx1 = miNeighbor.refIdx[1];
      if (refIdx0 >= 0 && refIdx1 >= 0 
#if !JVET_Y0089_DMVR_BCW
        && miNeighbor.BcwIdx == BCW_DEFAULT
#endif
#if MULTI_HYP_PRED
        && miNeighbor.addHypData.empty()
#endif
        )
      {
        if (cs.slice->getRefPic(REF_PIC_LIST_0, refIdx0)->longTerm
          || cs.slice->getRefPic(REF_PIC_LIST_1, refIdx1)->longTerm)
        {
          isBmCand = false;
        }
        const int poc0 = cs.slice->getRefPOC(REF_PIC_LIST_0, refIdx0);
        const int poc1 = cs.slice->getRefPOC(REF_PIC_LIST_1, refIdx1);
        const int poc = cs.slice->getPOC();
        if ((poc - poc0)*(poc - poc1) < 0)
        {
          if (abs(poc - poc0) == abs(poc - poc1))
          {
            isBmCand = true;
          }
        }
#if JVET_Y0128_NON_CTC
        if ( isBmCand && PU::isBiRefScaled( cs, refIdx0, refIdx1 ) )
        {
          isBmCand = false;
        }
#endif
      }
      if (!isBmCand)
      {
        continue;
      }
      mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
#if JVET_Z0075_IBC_HMVP_ENLARGE
      mrgCtx.useAltHpelIf      [cnt] = miNeighbor.useAltHpelIf;
#else
      mrgCtx.useAltHpelIf[cnt] = !ibcFlag && miNeighbor.useAltHpelIf;
#endif

      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
      if (slice.isInterB())
      {
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);
#if JVET_Y0089_DMVR_BCW
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? miNeighbor.BcwIdx : BCW_DEFAULT;
#endif
      }
#if NON_ADJACENT_MRG_CAND
      if (mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
        , mvdSimilarityThresh
#endif
      ))
      {
        continue;
      }
#endif
      if (mrgCandIdx == cnt)
      {
        return true;
      }
      cnt++;

      if (cnt == maxNumMergeCandMin1)
      {
        break;
      }
    }
  }

  if (cnt < maxNumMergeCandMin1)
  {
    mrgCtx.useAltHpelIf[cnt] = false;
  }
  return false;
}

void PU::getInterBMCandidates(const PredictionUnit &pu, MergeCtx& mrgCtx,
  const int& mrgCandIdx
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
  , MergeCtx* mvpMrgCtx1
  , MergeCtx* mvpMrgCtx2
#endif
)
{
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
  const CodingStructure &cs = *pu.cs;
  const Slice &slice = *pu.cs->slice;
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
  const uint32_t maxNumMergeCand = (pu.cs->sps->getUseAML()) ? BM_MRG_MAX_NUM_INIT_CANDS : pu.cs->sps->getMaxNumBMMergeCand();
#else
  const uint32_t maxNumMergeCand = pu.cs->sps->getMaxNumBMMergeCand();
#endif
  for (uint32_t ui = 0; ui < maxNumMergeCand; ++ui)
  {
    mrgCtx.BcwIdx[ui] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[ui] = false;
#endif
    mrgCtx.interDirNeighbours[ui] = 0;
    mrgCtx.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtx.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtx.useAltHpelIf[ui] = false;
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[ui].clear();
#endif
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
    mrgCtx.candCost[ui] = MAX_UINT64;
#endif
  }

  mrgCtx.numValidMergeCand = maxNumMergeCand;
  mrgCtx.numCandToTestEnc = maxNumMergeCand;
  // compute the location of the current PU

  int mvThreshod = 1;

  int cnt = 0;

  const Position posLT = pu.Y().topLeft();
  const Position posRT = pu.Y().topRight();
  const Position posLB = pu.Y().bottomLeft();
  MotionInfo miAbove, miLeft, miAboveLeft, miAboveRight, miBelowLeft;

  // above
  const PredictionUnit *puAbove = cs.getPURestricted(posRT.offset(0, -1), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
  bool isAvailableB1 = puAbove && isDiffMER(pu.lumaPos(), posRT.offset(0, -1), plevel) && pu.cu != puAbove->cu && CU::isInter(*puAbove->cu) && puAbove->getMotionInfo(posRT.offset(0, -1)).isInter;
#else
  bool isAvailableB1 = puAbove && isDiffMER(pu.lumaPos(), posRT.offset(0, -1), plevel) && pu.cu != puAbove->cu && CU::isInter(*puAbove->cu);
#endif

  if (isAvailableB1)
  {
    miAbove = puAbove->getMotionInfo(posRT.offset(0, -1));

    if (isBiPredFromDifferentDirEqDistPoc(pu, miAbove.refIdx[0], miAbove.refIdx[1])
      )
    {
      mrgCtx.interDirNeighbours[cnt] = miAbove.interDir;
      mrgCtx.useAltHpelIf[cnt] = miAbove.useAltHpelIf;
      mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miAbove.mv[0], miAbove.refIdx[0]);
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miAbove.mv[1], miAbove.refIdx[1]);
#if JVET_Y0089_DMVR_BCW
      mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puAbove->cu->BcwIdx : BCW_DEFAULT;
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
      if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
      {
#endif
        if (mrgCandIdx == cnt)
        {
          mrgCtx.numValidMergeCand = cnt + 1;
          return;
        }
        cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
      }
#endif
    }
  }

  // early termination
  if (cnt == maxNumMergeCand)
  {
    mrgCtx.numValidMergeCand = cnt;
    return;
  }

  //left
  const PredictionUnit* puLeft = cs.getPURestricted(posLB.offset(-1, 0), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
  const bool isAvailableA1 = puLeft && isDiffMER(pu.lumaPos(), posLB.offset(-1, 0), plevel) && pu.cu != puLeft->cu && CU::isInter(*puLeft->cu) && puLeft->getMotionInfo(posLB.offset(-1, 0)).isInter;
#else
  const bool isAvailableA1 = puLeft && isDiffMER(pu.lumaPos(), posLB.offset(-1, 0), plevel) && pu.cu != puLeft->cu && CU::isInter(*puLeft->cu);
#endif

  if (isAvailableA1)
  {
    miLeft = puLeft->getMotionInfo(posLB.offset(-1, 0));

    if (!isAvailableB1 || (miAbove != miLeft))
    {
      if (isBiPredFromDifferentDirEqDistPoc(pu, miLeft.refIdx[0], miLeft.refIdx[1])
        )
      {
        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miLeft.interDir;
        mrgCtx.useAltHpelIf[cnt] = miLeft.useAltHpelIf;
        // get Mv from Left
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miLeft.mv[0], miLeft.refIdx[0]);
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miLeft.mv[1], miLeft.refIdx[1]);
#if JVET_Y0089_DMVR_BCW
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puLeft->cu->BcwIdx : BCW_DEFAULT;
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
        {
#endif
          if (mrgCandIdx == cnt)
          {
            mrgCtx.numValidMergeCand = cnt + 1;
            return;
          }

          cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
        }
#endif
      }
    }
  }

  // early termination
  if (cnt == maxNumMergeCand)
  {
    mrgCtx.numValidMergeCand = cnt;
    return;
  }

  // above right
  const PredictionUnit *puAboveRight = cs.getPURestricted(posRT.offset(1, -1), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
  bool isAvailableB0 = puAboveRight && isDiffMER(pu.lumaPos(), posRT.offset(1, -1), plevel) && CU::isInter(*puAboveRight->cu) && puAboveRight->getMotionInfo(posRT.offset(1, -1)).isInter;
#else
  bool isAvailableB0 = puAboveRight && isDiffMER(pu.lumaPos(), posRT.offset(1, -1), plevel) && CU::isInter(*puAboveRight->cu);
#endif

  if (isAvailableB0)
  {
    miAboveRight = puAboveRight->getMotionInfo(posRT.offset(1, -1));

    if (!isAvailableB1 || (miAbove != miAboveRight))
    {
      if (isBiPredFromDifferentDirEqDistPoc(pu, miAboveRight.refIdx[0], miAboveRight.refIdx[1])
        )
      {
        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miAboveRight.interDir;
        mrgCtx.useAltHpelIf[cnt] = miAboveRight.useAltHpelIf;
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miAboveRight.mv[0], miAboveRight.refIdx[0]);
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miAboveRight.mv[1], miAboveRight.refIdx[1]);
#if JVET_Y0089_DMVR_BCW
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puAboveRight->cu->BcwIdx : BCW_DEFAULT;
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
        {
#endif
          if (mrgCandIdx == cnt)
          {
            mrgCtx.numValidMergeCand = cnt + 1;
            return;
          }
          cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
        }
#endif
      }
    }
  }
  // early termination
  if (cnt == maxNumMergeCand)
  {
    mrgCtx.numValidMergeCand = cnt;
    return;
  }

  //left bottom
  const PredictionUnit *puLeftBottom = cs.getPURestricted(posLB.offset(-1, 1), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
  bool isAvailableA0 = puLeftBottom && isDiffMER(pu.lumaPos(), posLB.offset(-1, 1), plevel) && CU::isInter(*puLeftBottom->cu) && puLeftBottom->getMotionInfo(posLB.offset(-1, 1)).isInter;
#else
  bool isAvailableA0 = puLeftBottom && isDiffMER(pu.lumaPos(), posLB.offset(-1, 1), plevel) && CU::isInter(*puLeftBottom->cu);
#endif

  if (isAvailableA0)
  {
    miBelowLeft = puLeftBottom->getMotionInfo(posLB.offset(-1, 1));

    if (!isAvailableA1 || (miBelowLeft != miLeft))
    {
      if (isBiPredFromDifferentDirEqDistPoc(pu, miBelowLeft.refIdx[0], miBelowLeft.refIdx[1])
        )
      {
        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miBelowLeft.interDir;
        mrgCtx.useAltHpelIf[cnt] = miBelowLeft.useAltHpelIf;
        // get Mv from Bottom-Left
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miBelowLeft.mv[0], miBelowLeft.refIdx[0]);
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miBelowLeft.mv[1], miBelowLeft.refIdx[1]);
#if JVET_Y0089_DMVR_BCW
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puLeftBottom->cu->BcwIdx : BCW_DEFAULT;
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
        {
#endif
          if (mrgCandIdx == cnt)
          {
#if TM_MRG
            if (!pu.tmMergeFlag)
#endif
              return;
          }
          cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
        }
#endif
      }
    }
  }
  // early termination
  if (cnt == maxNumMergeCand)
  {
    mrgCtx.numValidMergeCand = cnt;
    return;
  }


  const PredictionUnit *puAboveLeft = cs.getPURestricted(posLT.offset(-1, -1), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
  bool isAvailableB2 = puAboveLeft && isDiffMER(pu.lumaPos(), posLT.offset(-1, -1), plevel) && CU::isInter(*puAboveLeft->cu) && puAboveLeft->getMotionInfo(posLT.offset(-1, -1)).isInter;
#else
  bool isAvailableB2 = puAboveLeft && isDiffMER(pu.lumaPos(), posLT.offset(-1, -1), plevel) && CU::isInter(*puAboveLeft->cu);
#endif

  if (isAvailableB2)
  {
    miAboveLeft = puAboveLeft->getMotionInfo(posLT.offset(-1, -1));

    if ((!isAvailableA1 || (miLeft != miAboveLeft)) && (!isAvailableB1 || (miAbove != miAboveLeft)))
    {
      if (isBiPredFromDifferentDirEqDistPoc(pu, miAboveLeft.refIdx[0], miAboveLeft.refIdx[1])
        )
      {

        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miAboveLeft.interDir;
        mrgCtx.useAltHpelIf[cnt] = miAboveLeft.useAltHpelIf;
        // get Mv from Above-Left
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miAboveLeft.mv[0], miAboveLeft.refIdx[0]);
        mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miAboveLeft.mv[1], miAboveLeft.refIdx[1]);
#if JVET_Y0089_DMVR_BCW
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puAboveLeft->cu->BcwIdx : BCW_DEFAULT;
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
        {
#endif
          if (mrgCandIdx == cnt)
          {
            mrgCtx.numValidMergeCand = cnt + 1;
            return;
          }

          cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
        }
#endif
      }
    }
  }
  // early termination
  if (cnt == maxNumMergeCand)
  {
    mrgCtx.numValidMergeCand = cnt;
    return;
  }

#if INTER_RM_SIZE_CONSTRAINTS
  if (slice.getPicHeader()->getEnableTMVPFlag())
#else
  if (slice.getPicHeader()->getEnableTMVPFlag() && (pu.lumaSize().width + pu.lumaSize().height > 12))
#endif
  {
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
    if (mvpMrgCtx1 == NULL)
    {
#endif
    //>> MTK colocated-RightBottom
    // offset the pos to be sure to "point" to the same position the uiAbsPartIdx would've pointed to
    Position posRB = pu.Y().bottomRight().offset(-3, -3);
    const PreCalcValues& pcv = *cs.pcv;

    Position posC0;
    Position posC1 = pu.Y().center();
    bool C0Avail = false;
    bool boundaryCond = ((posRB.x + pcv.minCUWidth) < pcv.lumaWidth) && ((posRB.y + pcv.minCUHeight) < pcv.lumaHeight);
    const SubPic& curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
    if (curSubPic.getTreatedAsPicFlag())
    {
      boundaryCond = ((posRB.x + pcv.minCUWidth) <= curSubPic.getSubPicRight() &&
        (posRB.y + pcv.minCUHeight) <= curSubPic.getSubPicBottom());
    }
    if (boundaryCond)
    {
      int posYInCtu = posRB.y & pcv.maxCUHeightMask;
      if (posYInCtu + 4 < pcv.maxCUHeight)
      {
        posC0 = posRB.offset(4, 4);
        C0Avail = true;
      }
    }

    Mv        cColMv;
    int       iRefIdx = 0;
    int       dir = 0;
    unsigned  uiArrayAddr = cnt;
    bool      bExistMV = (C0Avail && getColocatedMVP(pu, REF_PIC_LIST_0, posC0, cColMv, iRefIdx, false))
      || getColocatedMVP(pu, REF_PIC_LIST_0, posC1, cColMv, iRefIdx, false);
    if (bExistMV)
    {
      dir |= 1;
      mrgCtx.mvFieldNeighbours[2 * uiArrayAddr].setMvField(cColMv, iRefIdx);
    }

    if (slice.isInterB())
    {
      bExistMV = (C0Avail && getColocatedMVP(pu, REF_PIC_LIST_1, posC0, cColMv, iRefIdx, false))
        || getColocatedMVP(pu, REF_PIC_LIST_1, posC1, cColMv, iRefIdx, false);
      if (bExistMV)
      {
        dir |= 2;
        mrgCtx.mvFieldNeighbours[2 * uiArrayAddr + 1].setMvField(cColMv, iRefIdx);
      }
    }

    if (dir != 0)
    {
      bool addTMvp = isBiPredFromDifferentDirEqDistPoc(pu, mrgCtx.mvFieldNeighbours[2 * uiArrayAddr + 0].refIdx, mrgCtx.mvFieldNeighbours[2 * uiArrayAddr + 1].refIdx);
      if (addTMvp)
      {
        mrgCtx.interDirNeighbours[uiArrayAddr] = dir;
        mrgCtx.useAltHpelIf[uiArrayAddr] = false;
#if MULTI_HYP_PRED
        mrgCtx.addHypNeighbours[uiArrayAddr].clear();
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
        {
#endif
          if (mrgCandIdx == cnt)
          {
            mrgCtx.numValidMergeCand = cnt + 1;
            return;
          }

          cnt++;
#if NON_ADJACENT_MRG_CAND || TM_MRG
        }
#endif
      }
    }
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
    }
    else if (mvpMrgCtx1->numValidMergeCand > 0)
    {
      mrgCtx.interDirNeighbours[cnt] = mvpMrgCtx1->interDirNeighbours[0];
      mrgCtx.mvFieldNeighbours[cnt << 1] = mvpMrgCtx1->mvFieldNeighbours[0];
      mrgCtx.useAltHpelIf[cnt] = mvpMrgCtx1->useAltHpelIf[0];
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 1] = mvpMrgCtx1->mvFieldNeighbours[1];
#if MULTI_HYP_PRED
      mrgCtx.addHypNeighbours[cnt] = mvpMrgCtx1->addHypNeighbours[0];
#endif
#if JVET_Y0089_DMVR_BCW
      mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? mvpMrgCtx1->BcwIdx[0] : BCW_DEFAULT;
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
      if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
#endif
      {
        mrgCtx.candCost[cnt] = mvpMrgCtx1->candCost[0];
        if (mrgCandIdx == cnt)
        {
          mrgCtx.numValidMergeCand = cnt + 1;
          return;
        }
        cnt++;
      }
    }
#endif
  }

  // early termination
  if (cnt == maxNumMergeCand)
  {
    mrgCtx.numValidMergeCand = cnt;
    return;
  }

  int maxNumMergeCandMin1 = maxNumMergeCand - 1;
#if NON_ADJACENT_MRG_CAND
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
  if (mvpMrgCtx2 == NULL)
  {
#endif
  MotionInfo miNeighbor;
  int offsetX = 0;
  int offsetY = 0;
  const int iNACANDIDATE_NUM[4] = { 3, 5, 5, 5 };
  const int idxMap[4][5] = { { 0, 1, 4 },{ 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4 } };

  for (int iDistanceIndex = 0; iDistanceIndex < NADISTANCE_LEVEL && cnt < maxNumMergeCandMin1; iDistanceIndex++)
  {
    const int iNADistanceHor = pu.Y().width  * (iDistanceIndex + 1);
    const int iNADistanceVer = pu.Y().height * (iDistanceIndex + 1);

    for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM[iDistanceIndex] && cnt < maxNumMergeCandMin1; NASPIdx++)
    {
      switch (idxMap[iDistanceIndex][NASPIdx])
      {
      case 0:offsetX = -iNADistanceHor - 1; offsetY = pu.Y().height + iNADistanceVer - 1; break;
      case 1:offsetX = pu.Y().width + iNADistanceHor - 1; offsetY = -iNADistanceVer - 1; break;
      case 2:offsetX = pu.Y().width >> 1;   offsetY = -iNADistanceVer - 1;    break;
      case 3:offsetX = -iNADistanceHor - 1; offsetY = pu.Y().height >> 1; break;
      case 4:offsetX = -iNADistanceHor - 1; offsetY = -iNADistanceVer - 1;    break;
      default: printf("error!"); exit(0); break;
      }

      const PredictionUnit *puNonAdjacent = cs.getPURestricted(posLT.offset(offsetX, offsetY), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu) && puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY)).isInter;
#else
      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu);
#endif

      if (isAvailableNonAdjacent)
      {
        miNeighbor = puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY));


        if (isBiPredFromDifferentDirEqDistPoc(pu, miNeighbor.refIdx[0], miNeighbor.refIdx[1])
          )
        {
          // get Inter Dir
          mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
          // get Mv from Above-Left
          mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
          mrgCtx.useAltHpelIf[cnt] = miNeighbor.useAltHpelIf;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);
#if JVET_Y0089_DMVR_BCW
          mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puNonAdjacent->cu->BcwIdx : BCW_DEFAULT;
#endif

#if NON_ADJACENT_MRG_CAND || TM_MRG
          if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
#endif
          {
            if (mrgCandIdx == cnt)
            {
              mrgCtx.numValidMergeCand = cnt + 1;
              return;
            }
            cnt++;
          }
        }
      }
    }
  }
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
  }
  else
  {
    for (int uiNumCand = 0; uiNumCand < mvpMrgCtx2->numValidMergeCand && cnt < maxNumMergeCandMin1; uiNumCand++)
    {
      mrgCtx.interDirNeighbours[cnt] = mvpMrgCtx2->interDirNeighbours[uiNumCand];
      mrgCtx.mvFieldNeighbours[cnt << 1] = mvpMrgCtx2->mvFieldNeighbours[uiNumCand << 1];
      mrgCtx.useAltHpelIf[cnt] = mvpMrgCtx2->useAltHpelIf[uiNumCand];
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 1] = mvpMrgCtx2->mvFieldNeighbours[(uiNumCand << 1) + 1];
#if JVET_Y0089_DMVR_BCW
      mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? mvpMrgCtx2->BcwIdx[uiNumCand] : BCW_DEFAULT;
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
      if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
#endif
      {
        mrgCtx.candCost[cnt] = mvpMrgCtx2->candCost[uiNumCand];
        if (mrgCandIdx == cnt)
        {
          mrgCtx.numValidMergeCand = cnt + 1;
          return;
        }
        cnt++;
      }
    }
  }
#endif
#endif

  if (cnt != maxNumMergeCandMin1)
  {
#if !JVET_Z0075_IBC_HMVP_ENLARGE
    bool isGt4x4 = true;
#endif
    bool bFound = addBMMergeHMVPCand(cs, mrgCtx, mrgCandIdx, maxNumMergeCandMin1, cnt
      , isAvailableA1, miLeft, isAvailableB1, miAbove
#if !JVET_Z0075_IBC_HMVP_ENLARGE
      , CU::isIBC(*pu.cu)
      , isGt4x4
#endif
#if TM_MRG
      , mvThreshod
#endif
    );

    if (bFound)
    {
      return;
    }
  }

  {
    if (cnt > 1 && cnt < maxNumMergeCand)
    {
      mrgCtx.mvFieldNeighbours[cnt * 2].setMvField(Mv(0, 0), NOT_VALID);
      mrgCtx.mvFieldNeighbours[cnt * 2 + 1].setMvField(Mv(0, 0), NOT_VALID);
#if INTER_LIC
      mrgCtx.LICFlags[cnt] = false;
#endif
#if MULTI_HYP_PRED
      mrgCtx.addHypNeighbours[cnt].clear();
#endif
      mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
      // calculate average MV for L0 and L1 seperately
      unsigned char interDir = 0;
      mrgCtx.useAltHpelIf[cnt] = (mrgCtx.useAltHpelIf[0] == mrgCtx.useAltHpelIf[1]) ? mrgCtx.useAltHpelIf[0] : false;
      for (int refListId = 0; refListId < (slice.isInterB() ? 2 : 1); refListId++)
      {
        const short refIdxI = mrgCtx.mvFieldNeighbours[0 * 2 + refListId].refIdx;
        const short refIdxJ = mrgCtx.mvFieldNeighbours[1 * 2 + refListId].refIdx;

        // both MVs are invalid, skip
        if (refIdxI != refIdxJ)
        {
          continue;
        }

        interDir += 1 << refListId;
        const Mv& MvI = mrgCtx.mvFieldNeighbours[0 * 2 + refListId].mv;
        const Mv& MvJ = mrgCtx.mvFieldNeighbours[1 * 2 + refListId].mv;

        // average two MVs
        Mv avgMv = MvI;
        avgMv += MvJ;
        roundAffineMv(avgMv.hor, avgMv.ver, 1);

        mrgCtx.mvFieldNeighbours[cnt * 2 + refListId].setMvField(avgMv, refIdxI);
      }

      mrgCtx.interDirNeighbours[cnt] = interDir;
      if (interDir == 3 && isBiPredFromDifferentDirEqDistPoc(pu, mrgCtx.mvFieldNeighbours[cnt * 2 + 0].refIdx, mrgCtx.mvFieldNeighbours[cnt * 2 + 1].refIdx))
      {
#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
#endif
          cnt++;
      }
    }

    // early termination
    if (cnt == maxNumMergeCand)
    {
      mrgCtx.numValidMergeCand = cnt;
      return;
    }
  }

  mrgCtx.numCandToTestEnc = cnt;

#if JVET_Y0128_NON_CTC
  while(cnt < maxNumMergeCand)
#else
  if (cnt < maxNumMergeCand)
#endif
  {
    mrgCtx.interDirNeighbours[cnt] = 3;
    mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[cnt] = false;
#endif
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[cnt].clear();
#endif
    mrgCtx.useAltHpelIf[cnt] = false;
#if JVET_Y0128_NON_CTC
    mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(Mv(0, 0), pu.cs->slice->getBMDefaultRefIdx(0));
    mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(0, 0), pu.cs->slice->getBMDefaultRefIdx(1));
#else
    mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(Mv(0, 0), 0);
    mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(0, 0), 0);
#endif
#if !JVET_Y0128_NON_CTC
#if NON_ADJACENT_MRG_CAND || TM_MRG
    if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
#endif
#endif
    {
      if (mrgCandIdx == cnt)
      {
        mrgCtx.numValidMergeCand = cnt + 1;
        return;
      }
      cnt++;
    }
  }
  mrgCtx.numValidMergeCand = cnt;
}
#endif

#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
void PU::getTmvpBMCand(const PredictionUnit &pu, MergeCtx& mrgCtx)
{
  const CodingStructure &cs = *pu.cs;
  const Slice &slice = *pu.cs->slice;

  const uint32_t mvThreshod = 1;
  const uint32_t maxNumMergeCand = NUM_TMVP_CANDS;
  for (uint32_t ui = 0; ui < maxNumMergeCand; ++ui)
  {
    mrgCtx.BcwIdx[ui] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[ui] = false;
#endif
    mrgCtx.interDirNeighbours[ui] = 0;
    mrgCtx.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtx.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtx.useAltHpelIf[ui] = false;
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[ui].clear();
#endif
    mrgCtx.candCost[ui] = MAX_UINT64;
  }

  int cnt = 0;

#if INTER_RM_SIZE_CONSTRAINTS
  if (slice.getPicHeader()->getEnableTMVPFlag())
#else
  if (slice.getPicHeader()->getEnableTMVPFlag() && (pu.lumaSize().width + pu.lumaSize().height > 12))
#endif
  {
    //>> MTK colocated-RightBottom
    // offset the pos to be sure to "point" to the same position the uiAbsPartIdx would've pointed to
    const PreCalcValues& pcv = *cs.pcv;
    bool C0Avail;
    bool C1Avail;
    bool boundaryCond;
    const SubPic& curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
    Position posRB = pu.Y().bottomRight().offset(-3, -3);
    Position posCenter = pu.Y().center();
    Position posC0;
    Position posC1;

    int       iRefIdx = 0;
    bool      bExistMV0, bExistMV1;
    Mv        cColMv0, cColMv1;
    int       dir;

    int offsetX0 = 0, offsetX1 = 0, offsetX2 = 0, offsetX3 = pu.Y().width >> 1;
    int offsetY0 = 0, offsetY1 = 0, offsetY2 = 0, offsetY3 = pu.Y().height >> 1;

    const int iNACANDIDATE_NUM[5] = { 2, 2, 2, 2, 2 };
    const int idxMap[5][2] = { { 0, 1 },{ 0, 2 },{ 0, 2 },{ 0, 2 },{ 0, 2 } };
    for (int iDistanceIndex = 0; iDistanceIndex < TMVP_DISTANCE_LEVEL && cnt < maxNumMergeCand; iDistanceIndex++)
    {
      const int iNADistanceHor = pu.Y().width  * iDistanceIndex;
      const int iNADistanceVer = pu.Y().height * iDistanceIndex;

      for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM[iDistanceIndex] && cnt < maxNumMergeCand; NASPIdx++)
      {
        switch (idxMap[iDistanceIndex][NASPIdx])
        {
        case 0: offsetX0 = offsetX2 = 4 + iNADistanceHor; offsetY0 = offsetY2 = 4 + iNADistanceVer; offsetX1 = iNADistanceHor; offsetY1 = iNADistanceVer; break;
        case 1: offsetX0 = 4; offsetY0 = 0; offsetX1 = 0; offsetY1 = 4; break;
        case 2: offsetX0 = offsetX2; offsetY0 = 4 - offsetY3; offsetX1 = 4 - offsetX3; offsetY1 = offsetY2; break;
        default: printf("error!"); exit(0); break;
        }
        C0Avail = false;
        if (curSubPic.getTreatedAsPicFlag())
        {
          boundaryCond = ((posRB.x + offsetX0) <= curSubPic.getSubPicRight() && (posRB.y + offsetY0) <= curSubPic.getSubPicBottom());
        }
        else
        {
          boundaryCond = ((posRB.x + offsetX0) < pcv.lumaWidth) && ((posRB.y + offsetY0) < pcv.lumaHeight);
        }
        if (boundaryCond)
        {
          int posYInCtu = posRB.y & pcv.maxCUHeightMask;
          if (posYInCtu + offsetY0 < pcv.maxCUHeight)
          {
            posC0 = posRB.offset(offsetX0, offsetY0);
            C0Avail = true;
          }
        }

        if (idxMap[iDistanceIndex][NASPIdx] == 0)
        {
          C1Avail = false;
          if (curSubPic.getTreatedAsPicFlag())
          {
            boundaryCond = ((posCenter.x + offsetX1) <= curSubPic.getSubPicRight() && (posCenter.y + offsetY1) <= curSubPic.getSubPicBottom());
          }
          else
          {
            boundaryCond = ((posCenter.x + offsetX1) < pcv.lumaWidth) && ((posCenter.y + offsetY1) < pcv.lumaHeight);
          }
          if (boundaryCond)
          {
            int posYInCtu = posCenter.y & pcv.maxCUHeightMask;
            if (posYInCtu + offsetY1 < pcv.maxCUHeight)
            {
              posC1 = posCenter.offset(offsetX1, offsetY1);
              C1Avail = true;
            }
          }
        }
        else
        {
          C1Avail = false;
          if (curSubPic.getTreatedAsPicFlag())
          {
            boundaryCond = ((posRB.x + offsetX1) <= curSubPic.getSubPicRight() && (posRB.y + offsetY1) <= curSubPic.getSubPicBottom());
          }
          else
          {
            boundaryCond = ((posRB.x + offsetX1) < pcv.lumaWidth) && ((posRB.y + offsetY1) < pcv.lumaHeight);
          }
          if (boundaryCond)
          {
            int posYInCtu = posRB.y & pcv.maxCUHeightMask;
            if (posYInCtu + offsetY1 < pcv.maxCUHeight)
            {
              posC1 = posRB.offset(offsetX1, offsetY1);
              C1Avail = true;
            }
          }
        }

        bExistMV0 = bExistMV1 = false;

        // Candidate with L0 and L1
        dir = 0;
        bExistMV0 = (C0Avail && getColocatedMVP(pu, REF_PIC_LIST_0, posC0, cColMv0, iRefIdx, false))
          || (C1Avail && getColocatedMVP(pu, REF_PIC_LIST_0, posC1, cColMv0, iRefIdx, false));
        if (bExistMV0)
        {
          dir |= 1;
          mrgCtx.mvFieldNeighbours[2 * cnt].setMvField(cColMv0, iRefIdx);
        }
        else
        {
          mrgCtx.mvFieldNeighbours[2 * cnt].setMvField(Mv(), NOT_VALID);
        }
        if (slice.isInterB())
        {
          bExistMV1 = (C0Avail && getColocatedMVP(pu, REF_PIC_LIST_1, posC0, cColMv1, iRefIdx, false))
            || (C1Avail && getColocatedMVP(pu, REF_PIC_LIST_1, posC1, cColMv1, iRefIdx, false));
          if (bExistMV1)
          {
            dir |= 2;
            mrgCtx.mvFieldNeighbours[2 * cnt + 1].setMvField(cColMv1, iRefIdx);
          }
          else
          {
            mrgCtx.mvFieldNeighbours[2 * cnt + 1].setMvField(Mv(), NOT_VALID);
          }
        }
        if (dir != 0)
        {
          bool addTMvp = isBiPredFromDifferentDirEqDistPoc(pu, mrgCtx.mvFieldNeighbours[2 * cnt + 0].refIdx, mrgCtx.mvFieldNeighbours[2 * cnt + 1].refIdx);
          if (addTMvp)
          {
            mrgCtx.interDirNeighbours[cnt] = dir;
            mrgCtx.useAltHpelIf[cnt] = false;
#if MULTI_HYP_PRED
            mrgCtx.addHypNeighbours[cnt].clear();
#endif
            if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod
            ))
            {
              cnt++;
              if (cnt == maxNumMergeCand)
              {
                break;
              }
            }
          }
        }
      }
    }
}

  mrgCtx.numValidMergeCand = cnt;
}
#endif

#if NON_ADJACENT_MRG_CAND && JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
void PU::getNonAdjacentBMCand(const PredictionUnit &pu, MergeCtx& mrgCtx)
{
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
  const CodingStructure &cs = *pu.cs;

  const uint32_t mvThreshod = getBDMVRMvdThreshold(pu);

  const uint32_t maxNumMergeCand = NUM_NON_ADJ_CANDS;

  for (uint32_t ui = 0; ui < maxNumMergeCand; ++ui)
  {
    mrgCtx.BcwIdx[ui] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[ui] = false;
#endif
    mrgCtx.interDirNeighbours[ui] = 0;
    mrgCtx.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtx.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtx.useAltHpelIf[ui] = false;
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[ui].clear();
#endif
    mrgCtx.candCost[ui] = MAX_UINT64;
  }

  int cnt = 0;

  const Position posLT = pu.Y().topLeft();

  MotionInfo miNeighbor;
  int offsetX = 0;
  int offsetY = 0;
  int offsetX0 = 0; int offsetX1 = 0; int offsetX2 = pu.Y().width >> 1;
  int offsetY0 = 0; int offsetY1 = 0; int offsetY2 = pu.Y().height >> 1;

  const int iNACANDIDATE_NUM[7] = { 5, 9, 9, 9, 9, 9, 9 };
  const int idxMap[7][9] = { { 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 } };

  for (int iDistanceIndex = 0; iDistanceIndex < 7 && cnt < maxNumMergeCand; iDistanceIndex++)
  {
    const int iNADistanceHor = pu.Y().width  * (iDistanceIndex + 1);
    const int iNADistanceVer = pu.Y().height * (iDistanceIndex + 1);

    for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM[iDistanceIndex] && cnt < maxNumMergeCand; NASPIdx++)
    {
      switch (idxMap[iDistanceIndex][NASPIdx])
      {
      case 0:offsetX = offsetX0 = -iNADistanceHor - 1;               offsetY = offsetY0 = pu.Y().height + iNADistanceVer - 1; break;
      case 1:offsetX = offsetX1 = pu.Y().width + iNADistanceHor - 1; offsetY = offsetY1 = -iNADistanceVer - 1;                break;
      case 2:offsetX = offsetX2;                                     offsetY = offsetY1;                                      break;
      case 3:offsetX = offsetX0;                                     offsetY = offsetY2;                                      break;
      case 4:offsetX = offsetX0;                                     offsetY = offsetY1;                                      break;
      case 5:offsetX = -1;                                           offsetY = offsetY0;                                      break;
      case 6:offsetX = offsetX1;                                     offsetY = -1;                                            break;
      case 7:offsetX = offsetX0 >> 1;                                offsetY = offsetY0;                                      break;
      case 8:offsetX = offsetX1;                                     offsetY = offsetY1 >> 1;                                 break;
      default: printf("error!"); exit(0); break;
      }

      const PredictionUnit *puNonAdjacent = cs.getPURestricted(posLT.offset(offsetX, offsetY), pu, pu.chType);

      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu);

      if (isAvailableNonAdjacent)
      {
        miNeighbor = puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY));

        if (isBiPredFromDifferentDirEqDistPoc(pu, miNeighbor.refIdx[0], miNeighbor.refIdx[1])
          )
        {
          // get Inter Dir
          mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
          mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
          mrgCtx.useAltHpelIf[cnt] = miNeighbor.useAltHpelIf;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);

          if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
          {
            cnt++;
          }
        }
      }
    }
  }

  const int iNACANDIDATE_NUM2[7] = { 4, 4, 4, 4, 4, 4, 4 };
  const int idxMap2[7][5] = { { 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 } };

  for (int iDistanceIndex = 0; iDistanceIndex < 7 && cnt < maxNumMergeCand; iDistanceIndex++)
  {
    const int iNADistanceHor = pu.Y().width  * (iDistanceIndex + 1);
    const int iNADistanceVer = pu.Y().height * (iDistanceIndex + 1);

    for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM2[iDistanceIndex] && cnt < maxNumMergeCand; NASPIdx++)
    {
      switch (idxMap2[iDistanceIndex][NASPIdx])
      {
      case 0:offsetX = offsetX0 = -iNADistanceHor - 1;                                                         offsetY = offsetY2 + ((pu.Y().height + iNADistanceVer - 1 - offsetY2) >> 1); break;
      case 1:offsetX = offsetX2 + ((pu.Y().width + iNADistanceHor - 1 - offsetX2) >> 1);                       offsetY = offsetY0 = -iNADistanceVer - 1; break;
      case 2:offsetX = offsetX0;                                                                               offsetY = offsetY0 + ((offsetY2 - offsetY0) >> 1); break;
      case 3:offsetX = offsetX0 + ((offsetX2 - offsetX0) >> 1);                                                offsetY = offsetY0; break;
      default: printf("error!"); exit(0); break;
      }

      const PredictionUnit *puNonAdjacent = cs.getPURestricted(posLT.offset(offsetX, offsetY), pu, pu.chType);

      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu);

      if (isAvailableNonAdjacent)
      {
        miNeighbor = puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY));

        if (isBiPredFromDifferentDirEqDistPoc(pu, miNeighbor.refIdx[0], miNeighbor.refIdx[1])
          )
        {
          // get Inter Dir
          mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
          // get Mv from Above-Left
          mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
          mrgCtx.useAltHpelIf[cnt] = miNeighbor.useAltHpelIf;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);

          if (!mrgCtx.xCheckSimilarMotion(cnt, mvThreshod))
          {
            cnt++;
          }
        }
      }
    }
  }

  mrgCtx.numValidMergeCand = cnt;
}
#endif

#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
void PU::getTmvpMergeCand(const PredictionUnit &pu, MergeCtx& mrgCtx)
{
  const CodingStructure &cs = *pu.cs;
  const Slice &slice = *pu.cs->slice;

  const uint32_t mvdSimilarityThresh = 1;
  const uint32_t maxNumMergeCand = NUM_TMVP_CANDS;
  for (uint32_t ui = 0; ui < maxNumMergeCand; ++ui)
  {
    mrgCtx.BcwIdx[ui] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[ui] = false;
#endif
    mrgCtx.interDirNeighbours[ui] = 0;
    mrgCtx.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtx.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtx.useAltHpelIf[ui] = false;
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[ui].clear();
#endif
    mrgCtx.candCost[ui] = MAX_UINT64;
  }

  int cnt = 0;

#if INTER_RM_SIZE_CONSTRAINTS
  if (slice.getPicHeader()->getEnableTMVPFlag())
#else
  if (slice.getPicHeader()->getEnableTMVPFlag() && (pu.lumaSize().width + pu.lumaSize().height > 12))
#endif
  {
    //>> MTK colocated-RightBottom
    // offset the pos to be sure to "point" to the same position the uiAbsPartIdx would've pointed to
    const PreCalcValues& pcv = *cs.pcv;
    bool C0Avail;
    bool C1Avail;
    bool boundaryCond;
    const SubPic& curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
    Position posRB = pu.Y().bottomRight().offset(-3, -3);
    Position posCenter = pu.Y().center();
    Position posC0;
    Position posC1;

    int       iRefIdx = 0;
    bool      bExistMV0, bExistMV1;
    Mv        cColMv0, cColMv1;
    int       dir;

    int offsetX0 = 0, offsetX1 = 0, offsetX2 = 0, offsetX3 = pu.Y().width >> 1;
    int offsetY0 = 0, offsetY1 = 0, offsetY2 = 0, offsetY3 = pu.Y().height >> 1;

    const int iNACANDIDATE_NUM[5] = { 2, 2, 2, 2, 2 };
    const int idxMap[5][2] = { { 0, 1 },{ 0, 2 },{ 0, 2 },{ 0, 2 },{ 0, 2 } };
    for (int iDistanceIndex = 0; iDistanceIndex < TMVP_DISTANCE_LEVEL && cnt < maxNumMergeCand; iDistanceIndex++)
    {
      const int iNADistanceHor = pu.Y().width  * iDistanceIndex;
      const int iNADistanceVer = pu.Y().height * iDistanceIndex;

      for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM[iDistanceIndex] && cnt < maxNumMergeCand; NASPIdx++)
      {
        switch (idxMap[iDistanceIndex][NASPIdx])
        {
        case 0: offsetX0 = offsetX2 = 4 + iNADistanceHor; offsetY0 = offsetY2 = 4 + iNADistanceVer; offsetX1 = iNADistanceHor; offsetY1 = iNADistanceVer; break;
        case 1: offsetX0 = 4; offsetY0 = 0; offsetX1 = 0; offsetY1 = 4; break;
        case 2: offsetX0 = offsetX2; offsetY0 = 4 - offsetY3; offsetX1 = 4 - offsetX3; offsetY1 = offsetY2; break;
        default: printf("error!"); exit(0); break;
        }
        C0Avail = false;
        if (curSubPic.getTreatedAsPicFlag())
        {
          boundaryCond = ((posRB.x + offsetX0) <= curSubPic.getSubPicRight() && (posRB.y + offsetY0) <= curSubPic.getSubPicBottom());
        }
        else
        {
          boundaryCond = ((posRB.x + offsetX0) < pcv.lumaWidth) && ((posRB.y + offsetY0) < pcv.lumaHeight);
        }
        if (boundaryCond)
        {
          int posYInCtu = posRB.y & pcv.maxCUHeightMask;
          if (posYInCtu + offsetY0 < pcv.maxCUHeight)
          {
            posC0 = posRB.offset(offsetX0, offsetY0);
            C0Avail = true;
          }
        }

        if (idxMap[iDistanceIndex][NASPIdx] == 0)
        {
          C1Avail = false;
          if (curSubPic.getTreatedAsPicFlag())
          {
            boundaryCond = ((posCenter.x + offsetX1) <= curSubPic.getSubPicRight() && (posCenter.y + offsetY1) <= curSubPic.getSubPicBottom());
          }
          else
          {
            boundaryCond = ((posCenter.x + offsetX1) < pcv.lumaWidth) && ((posCenter.y + offsetY1) < pcv.lumaHeight);
          }
          if (boundaryCond)
          {
            int posYInCtu = posCenter.y & pcv.maxCUHeightMask;
            if (posYInCtu + offsetY1 < pcv.maxCUHeight)
            {
              posC1 = posCenter.offset(offsetX1, offsetY1);
              C1Avail = true;
            }
          }
        }
        else
        {
          C1Avail = false;
          if (curSubPic.getTreatedAsPicFlag())
          {
            boundaryCond = ((posRB.x + offsetX1) <= curSubPic.getSubPicRight() && (posRB.y + offsetY1) <= curSubPic.getSubPicBottom());
          }
          else
          {
            boundaryCond = ((posRB.x + offsetX1) < pcv.lumaWidth) && ((posRB.y + offsetY1) < pcv.lumaHeight);
          }
          if (boundaryCond)
          {
            int posYInCtu = posRB.y & pcv.maxCUHeightMask;
            if (posYInCtu + offsetY1 < pcv.maxCUHeight)
            {
              posC1 = posRB.offset(offsetX1, offsetY1);
              C1Avail = true;
            }
          }
        }

        bExistMV0 = bExistMV1 = false;

        // Candidate with L0 and L1
        dir = 0;
        int refIdx[NUM_REF_PIC_LIST_01] = { 0, 0 };
        bExistMV0 = (C0Avail && getColocatedMVP(pu, REF_PIC_LIST_0, posC0, cColMv0, iRefIdx, false, &refIdx[REF_PIC_LIST_0]))
          || (C1Avail && getColocatedMVP(pu, REF_PIC_LIST_0, posC1, cColMv0, iRefIdx, false, &refIdx[REF_PIC_LIST_0]));

        if (bExistMV0)
        {
          dir |= 1;
          mrgCtx.mvFieldNeighbours[2 * cnt].setMvField(cColMv0, refIdx[REF_PIC_LIST_0]);
        }
        else
        {
          mrgCtx.mvFieldNeighbours[2 * cnt].setMvField(Mv(), NOT_VALID);
        }
        if (slice.isInterB())
        {
          bExistMV1 = (C0Avail && getColocatedMVP(pu, REF_PIC_LIST_1, posC0, cColMv1, iRefIdx, false, &refIdx[REF_PIC_LIST_1]))
            || (C1Avail && getColocatedMVP(pu, REF_PIC_LIST_1, posC1, cColMv1, iRefIdx, false, &refIdx[REF_PIC_LIST_1]));

          if (bExistMV1)
          {
            dir |= 2;
            mrgCtx.mvFieldNeighbours[2 * cnt + 1].setMvField(cColMv1, refIdx[REF_PIC_LIST_1]);
          }
          else
          {
            mrgCtx.mvFieldNeighbours[2 * cnt + 1].setMvField(Mv(), NOT_VALID);
          }
        }
        if (dir != 0)
        {
          mrgCtx.interDirNeighbours[cnt] = dir;
          mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
          mrgCtx.LICFlags[cnt] = false;
#endif
          mrgCtx.useAltHpelIf[cnt] = false;
#if MULTI_HYP_PRED
          mrgCtx.addHypNeighbours[cnt].clear();
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
          if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
            , mvdSimilarityThresh
#endif
          ))
          {
#endif
            cnt++;
            if (cnt == maxNumMergeCand)
            {
              break;
            }
#if NON_ADJACENT_MRG_CAND || TM_MRG
          }
#endif
        }

        if (!slice.getCheckLDC() && bExistMV0 && bExistMV1)
        {
          // Candidate without L1
          dir = 0;
          dir |= 1;
          mrgCtx.mvFieldNeighbours[2 * cnt].setMvField(cColMv0, refIdx[REF_PIC_LIST_0]);
          mrgCtx.mvFieldNeighbours[2 * cnt + 1].setMvField(Mv(), NOT_VALID);
          if (dir != 0)
          {
            mrgCtx.interDirNeighbours[cnt] = dir;
            mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
            mrgCtx.LICFlags[cnt] = false;
#endif
            mrgCtx.useAltHpelIf[cnt] = false;
#if MULTI_HYP_PRED
            mrgCtx.addHypNeighbours[cnt].clear();
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
            if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
              , mvdSimilarityThresh
#endif
            ))
            {
#endif
              cnt++;
              if (cnt == maxNumMergeCand) break;
#if NON_ADJACENT_MRG_CAND || TM_MRG
            }
#endif
          }

          // Candidate without L0
          dir = 0;
          mrgCtx.mvFieldNeighbours[2 * cnt].setMvField(Mv(), NOT_VALID);
          dir |= 2;
          mrgCtx.mvFieldNeighbours[2 * cnt + 1].setMvField(cColMv1, refIdx[REF_PIC_LIST_1]);
          if (dir != 0)
          {
            mrgCtx.interDirNeighbours[cnt] = dir;
            mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
            mrgCtx.LICFlags[cnt] = false;
#endif
            mrgCtx.useAltHpelIf[cnt] = false;
#if MULTI_HYP_PRED
            mrgCtx.addHypNeighbours[cnt].clear();
#endif
#if NON_ADJACENT_MRG_CAND || TM_MRG
            if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
              , mvdSimilarityThresh
#endif
            ))
            {
#endif
              cnt++;
              if (cnt == maxNumMergeCand)
              {
                break;
              }
#if NON_ADJACENT_MRG_CAND || TM_MRG
            }
#endif
          }
        }
      }
    }
  }

  mrgCtx.numValidMergeCand = cnt;
}
#endif

#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM && NON_ADJACENT_MRG_CAND
void PU::getNonAdjacentMergeCand(const PredictionUnit &pu, MergeCtx& mrgCtx)
{
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
  const CodingStructure &cs = *pu.cs;
  const Slice &slice = *pu.cs->slice;

  const uint32_t mvdSimilarityThresh = getBDMVRMvdThreshold(pu);
  const uint32_t maxNumMergeCand = NUM_NON_ADJ_CANDS;

  for (uint32_t ui = 0; ui < maxNumMergeCand; ++ui)
  {
    mrgCtx.BcwIdx[ui] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[ui] = false;
#endif
    mrgCtx.interDirNeighbours[ui] = 0;
    mrgCtx.mvFieldNeighbours[(ui << 1)].refIdx = NOT_VALID;
    mrgCtx.mvFieldNeighbours[(ui << 1) + 1].refIdx = NOT_VALID;
    mrgCtx.useAltHpelIf[ui] = false;
#if MULTI_HYP_PRED
    mrgCtx.addHypNeighbours[ui].clear();
#endif
    mrgCtx.candCost[ui] = MAX_UINT64;
  }

  int cnt = 0;

  const Position posLT = pu.Y().topLeft();

  MotionInfo miNeighbor;
  int offsetX = 0;
  int offsetY = 0;
  int offsetX0 = 0; int offsetX1 = 0; int offsetX2 = pu.Y().width >> 1;
  int offsetY0 = 0; int offsetY1 = 0; int offsetY2 = pu.Y().height >> 1;

  const int iNACANDIDATE_NUM[7] = { 5, 9, 9, 9, 9, 9, 9 };
  const int idxMap[7][9] = { { 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8 } };

  for (int iDistanceIndex = 0; iDistanceIndex < 7 && cnt < maxNumMergeCand; iDistanceIndex++)
  {
    const int iNADistanceHor = pu.Y().width  * (iDistanceIndex + 1);
    const int iNADistanceVer = pu.Y().height * (iDistanceIndex + 1);

    for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM[iDistanceIndex] && cnt < maxNumMergeCand; NASPIdx++)
    {
      switch (idxMap[iDistanceIndex][NASPIdx])
      {
      case 0:offsetX = offsetX0 = -iNADistanceHor - 1;               offsetY = offsetY0 = pu.Y().height + iNADistanceVer - 1; break;
      case 1:offsetX = offsetX1 = pu.Y().width + iNADistanceHor - 1; offsetY = offsetY1 = -iNADistanceVer - 1;                break;
      case 2:offsetX = offsetX2;                                     offsetY = offsetY1;                                      break;
      case 3:offsetX = offsetX0;                                     offsetY = offsetY2;                                      break;
      case 4:offsetX = offsetX0;                                     offsetY = offsetY1;                                      break;
      case 5:offsetX = -1;                                           offsetY = offsetY0;                                      break;
      case 6:offsetX = offsetX1;                                     offsetY = -1;                                            break;
      case 7:offsetX = offsetX0 >> 1;                                offsetY = offsetY0;                                      break;
      case 8:offsetX = offsetX1;                                     offsetY = offsetY1 >> 1;                                 break;
      default: printf("error!"); exit(0); break;
      }

      const PredictionUnit *puNonAdjacent = cs.getPURestricted(posLT.offset(offsetX, offsetY), pu, pu.chType);

      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu);

      if (isAvailableNonAdjacent)
      {
        miNeighbor = puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY));

        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
        // get Mv from Above-Left
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
        mrgCtx.useAltHpelIf[cnt] = miNeighbor.useAltHpelIf;
        // get Mv from Above-right
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puNonAdjacent->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
        mrgCtx.LICFlags[cnt] = miNeighbor.usesLIC;
#endif
        if (slice.isInterB())
        {
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);
#if MULTI_HYP_PRED
          mrgCtx.addHypNeighbours[cnt] = puNonAdjacent->addHypData;
#endif
        }

#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
          , mvdSimilarityThresh
#endif
        ))
#endif
        {
          cnt++;
        }
      }

    }
  }

  const int iNACANDIDATE_NUM2[7] = { 4, 4, 4, 4, 4, 4, 4 };
  const int idxMap2[7][5] = { { 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 } };

  for (int iDistanceIndex = 0; iDistanceIndex < 7 && cnt < maxNumMergeCand; iDistanceIndex++)
  {
    const int iNADistanceHor = pu.Y().width  * (iDistanceIndex + 1);
    const int iNADistanceVer = pu.Y().height * (iDistanceIndex + 1);

    for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM2[iDistanceIndex] && cnt < maxNumMergeCand; NASPIdx++)
    {
      switch (idxMap2[iDistanceIndex][NASPIdx])
      {
      case 0:offsetX = offsetX0 = -iNADistanceHor - 1;                                                         offsetY = offsetY2 + ((pu.Y().height + iNADistanceVer - 1 - offsetY2) >> 1); break;
      case 1:offsetX = offsetX2 + ((pu.Y().width + iNADistanceHor - 1 - offsetX2) >> 1);                       offsetY = offsetY0 = -iNADistanceVer - 1; break;
      case 2:offsetX = offsetX0;                                                                               offsetY = offsetY0 + ((offsetY2 - offsetY0) >> 1); break;
      case 3:offsetX = offsetX0 + ((offsetX2 - offsetX0) >> 1);                                                offsetY = offsetY0; break;
      default: printf("error!"); exit(0); break;
      }

      const PredictionUnit *puNonAdjacent = cs.getPURestricted(posLT.offset(offsetX, offsetY), pu, pu.chType);

      bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu);

      if (isAvailableNonAdjacent)
      {
        miNeighbor = puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY));

        // get Inter Dir
        mrgCtx.interDirNeighbours[cnt] = miNeighbor.interDir;
        // get Mv from Above-Left
        mrgCtx.mvFieldNeighbours[cnt << 1].setMvField(miNeighbor.mv[0], miNeighbor.refIdx[0]);
        mrgCtx.useAltHpelIf[cnt] = miNeighbor.useAltHpelIf;
        // get Mv from Above-right
        mrgCtx.BcwIdx[cnt] = (mrgCtx.interDirNeighbours[cnt] == 3) ? puNonAdjacent->cu->BcwIdx : BCW_DEFAULT;
#if INTER_LIC
        mrgCtx.LICFlags[cnt] = miNeighbor.usesLIC;
#endif
        if (slice.isInterB())
        {
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(miNeighbor.mv[1], miNeighbor.refIdx[1]);
#if MULTI_HYP_PRED
          mrgCtx.addHypNeighbours[cnt] = puNonAdjacent->addHypData;
#endif
        }

#if NON_ADJACENT_MRG_CAND || TM_MRG
        if (!mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
          , mvdSimilarityThresh
#endif
        ))
#endif
        {
          cnt++;
        }
      }
    }
  }

  mrgCtx.numValidMergeCand = cnt;
}
#endif

bool PU::checkDMVRCondition(const PredictionUnit& pu)
{
#if !MULTI_PASS_DMVR
  if( pu.cs->sps->getUseDMVR() && !pu.cs->picHeader->getDisDmvrFlag() )
  {
    const int refIdx0 = pu.refIdx[REF_PIC_LIST_0];
    const int refIdx1 = pu.refIdx[REF_PIC_LIST_1];

    const WPScalingParam *wp0 = pu.cu->slice->getWpScaling( REF_PIC_LIST_0, refIdx0 );
    const WPScalingParam *wp1 = pu.cu->slice->getWpScaling( REF_PIC_LIST_1, refIdx1 );

    const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();

    const bool ref0IsScaled = refIdx0 < 0 || refIdx0 >= MAX_NUM_REF
      ? false
      : isResamplingPossible && pu.cu->slice->getRefPic( REF_PIC_LIST_0, refIdx0 )->isRefScaled( pu.cs->pps );
    const bool ref1IsScaled = refIdx1 < 0 || refIdx1 >= MAX_NUM_REF
      ? false
      : isResamplingPossible && pu.cu->slice->getRefPic( REF_PIC_LIST_1, refIdx1 )->isRefScaled( pu.cs->pps );

    return pu.mergeFlag && pu.mergeType == MRG_TYPE_DEFAULT_N && !pu.ciipFlag && !pu.cu->affine && !pu.mmvdMergeFlag
      && !pu.cu->mmvdSkip && PU::isBiPredFromDifferentDirEqDistPoc( pu )
#if TM_MRG
      && !pu.tmMergeFlag
#endif
      && ( pu.lheight() >= 8 ) && ( pu.lwidth() >= 8 )
      && ( ( pu.lheight() * pu.lwidth() ) >= 128 )
      && ( pu.cu->BcwIdx == BCW_DEFAULT )
      && !WPScalingParam::isWeighted( wp0 ) && !WPScalingParam::isWeighted( wp1 ) && !ref0IsScaled && !ref1IsScaled;
  }
  else
#endif
  {
    return false;
  }
}

#if MULTI_PASS_DMVR
bool PU::checkBDMVRCondition(const PredictionUnit& pu)
{
  if( pu.cs->sps->getUseDMVDMode() )
  {
    const int refIdx0 = pu.refIdx[REF_PIC_LIST_0];
    const int refIdx1 = pu.refIdx[REF_PIC_LIST_1];

    const WPScalingParam *wp0 = pu.cu->slice->getWpScaling( REF_PIC_LIST_0, refIdx0 );
    const WPScalingParam *wp1 = pu.cu->slice->getWpScaling( REF_PIC_LIST_1, refIdx1 );

    const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();

    const bool ref0IsScaled = refIdx0 < 0 || refIdx0 >= MAX_NUM_REF
      ? false
      : isResamplingPossible && pu.cu->slice->getRefPic( REF_PIC_LIST_0, refIdx0 )->isRefScaled( pu.cs->pps );
    const bool ref1IsScaled = refIdx1 < 0 || refIdx1 >= MAX_NUM_REF
      ? false
      : isResamplingPossible && pu.cu->slice->getRefPic( REF_PIC_LIST_1, refIdx1 )->isRefScaled( pu.cs->pps );

#if JVET_X0083_BM_AMVP_MERGE_MODE
    return ((pu.mergeFlag && pu.mergeType == MRG_TYPE_DEFAULT_N) || (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])) && !pu.ciipFlag && !pu.cu->affine && !pu.mmvdMergeFlag
#else
    return pu.mergeFlag && pu.mergeType == MRG_TYPE_DEFAULT_N && !pu.ciipFlag && !pu.cu->affine && !pu.mmvdMergeFlag
#endif
      && !pu.cu->mmvdSkip
      && PU::isBiPredFromDifferentDirEqDistPoc( pu )
#if !JVET_Y0089_DMVR_BCW
      && ( pu.cu->BcwIdx == BCW_DEFAULT )
#endif
      && !WPScalingParam::isWeighted( wp0 ) && !WPScalingParam::isWeighted( wp1 ) && !ref0IsScaled && !ref1IsScaled;
  }
  else
  {
    return false;
  }
}
#endif

#if INTER_LIC && RPR_ENABLE
bool PU::checkRprLicCondition( const PredictionUnit & pu )
{
  if ( pu.cs->sps->getLicEnabledFlag() )
  {
#if JVET_Y0128_NON_CTC
    if(pu.interDir != 1 && pu.interDir != 2)
#else
    if( pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0 )
#endif
    {
      return false;
    }

    bool  bLicEnabled = true;
#if JVET_Y0128_NON_CTC
    int iRefIdx             = (pu.interDir == 1) ? pu.refIdx[0] : pu.refIdx[1];
    RefPicList  eRefPicList = (pu.interDir == 1) ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
#else
    int iRefIdx             = (pu.refIdx[0] >= 0) ? pu.refIdx[0] : pu.refIdx[1];
    RefPicList  eRefPicList = (pu.refIdx[0] >= 0) ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
#endif

    const bool isResamplingPossible = pu.cs->sps->getRprEnabledFlag();
    if( isResamplingPossible && pu.cu->slice->getRefPic( eRefPicList, iRefIdx )->isRefScaled( pu.cs->pps ) )
    {
      bLicEnabled = false;
    }

    return bLicEnabled;
  }
  else
  {
    return false;
  }
}
#endif

#if JVET_Y0128_NON_CTC
bool PU::checkTmEnableCondition(const SPS* sps, const PPS* pps, const Picture* refPic)
{
#if TM_AMVP || TM_MRG || MULTI_PASS_DMVR
  if (sps->getUseDMVDMode())
  {
    if (sps->getRprEnabledFlag())
    {
      return !refPic->isRefScaled(pps);
    }
    else
    {
      return true;
    }
  }
  else
#endif
  {
    return false;
  }
}
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
bool PU::checkRprRefExistingInGpm(const PredictionUnit& pu, const MergeCtx& geoMrgCtx0, uint8_t candIdx0, const MergeCtx& geoMrgCtx1, uint8_t candIdx1)
{
  if (pu.cs->sps->getRprEnabledFlag())
  {
    auto xCheckUseRprPerPart = [&pu](const MergeCtx& mrgCtx, uint8_t candIdx)
    {
#if JVET_Y0065_GPM_INTRA
      if (candIdx >= GEO_MAX_NUM_UNI_CANDS)
      {
        return false;
      }
#endif
      int refList = mrgCtx.interDirNeighbours[candIdx] - 1;
      int refIdx  = mrgCtx.mvFieldNeighbours[(candIdx << 1) + refList].refIdx;
      return pu.cu->slice->getRefPic((RefPicList)refList, refIdx)->isRefScaled(pu.cs->pps);
    };

    return xCheckUseRprPerPart(geoMrgCtx0, candIdx0) || xCheckUseRprPerPart(geoMrgCtx1, candIdx1);
  }

  return false;
}
#endif

static int xGetDistScaleFactor(const int &iCurrPOC, const int &iCurrRefPOC, const int &iColPOC, const int &iColRefPOC)
{
  int iDiffPocD = iColPOC - iColRefPOC;
  int iDiffPocB = iCurrPOC - iCurrRefPOC;

  if (iDiffPocD == iDiffPocB)
  {
    return 4096;
  }
  else
  {
    int iTDB = Clip3(-128, 127, iDiffPocB);
    int iTDD = Clip3(-128, 127, iDiffPocD);
    int iX = (0x4000 + abs(iTDD / 2)) / iTDD;
    int iScale = Clip3(-4096, 4095, (iTDB * iX + 32) >> 6);
    return iScale;
  }
}

int convertMvFixedToFloat(int32_t val)
{
  int sign  = val >> 31;
  int scale = floorLog2((val ^ sign) | MV_MANTISSA_UPPER_LIMIT) - (MV_MANTISSA_BITCOUNT - 1);

  int exponent;
  int mantissa;
  if (scale >= 0)
  {
    int round = (1 << scale) >> 1;
    int n     = (val + round) >> scale;
    exponent  = scale + ((n ^ sign) >> (MV_MANTISSA_BITCOUNT - 1));
    mantissa  = (n & MV_MANTISSA_UPPER_LIMIT) | (sign << (MV_MANTISSA_BITCOUNT - 1));
  }
  else
  {
    exponent = 0;
    mantissa = val;
  }

  return exponent | (mantissa << MV_EXPONENT_BITCOUNT);
}

int convertMvFloatToFixed(int val)
{
  int exponent = val & MV_EXPONENT_MASK;
  int mantissa = val >> MV_EXPONENT_BITCOUNT;
  return exponent == 0 ? mantissa : (mantissa ^ MV_MANTISSA_LIMIT) << (exponent - 1);
}

int roundMvComp(int x)
{
  return convertMvFloatToFixed(convertMvFixedToFloat(x));
}

int PU::getDistScaleFactor(const int &currPOC, const int &currRefPOC, const int &colPOC, const int &colRefPOC)
{
  return xGetDistScaleFactor(currPOC, currRefPOC, colPOC, colRefPOC);
}

void PU::getInterMMVDMergeCandidates(const PredictionUnit &pu, MergeCtx& mrgCtx, const int& mrgCandIdx)
{
  int refIdxList0, refIdxList1;
  int k;
  int currBaseNum = 0;
  const uint16_t maxNumMergeCand = mrgCtx.numValidMergeCand;

  for( k = 0; k < maxNumMergeCand; k++ )
  {
    refIdxList0 = mrgCtx.mvFieldNeighbours[( k << 1 )].refIdx;
    refIdxList1 = mrgCtx.mvFieldNeighbours[( k << 1 ) + 1].refIdx;

    if( ( refIdxList0 >= 0 ) && ( refIdxList1 >= 0 ) )
    {
      mrgCtx.mmvdBaseMv[currBaseNum][0] = mrgCtx.mvFieldNeighbours[( k << 1 )];
      mrgCtx.mmvdBaseMv[currBaseNum][1] = mrgCtx.mvFieldNeighbours[( k << 1 ) + 1];
    }
    else if( refIdxList0 >= 0 )
    {
      mrgCtx.mmvdBaseMv[currBaseNum][0] = mrgCtx.mvFieldNeighbours[( k << 1 )];
      mrgCtx.mmvdBaseMv[currBaseNum][1] = MvField( Mv( 0, 0 ), -1 );
    }
    else if( refIdxList1 >= 0 )
    {
      mrgCtx.mmvdBaseMv[currBaseNum][0] = MvField( Mv( 0, 0 ), -1 );
      mrgCtx.mmvdBaseMv[currBaseNum][1] = mrgCtx.mvFieldNeighbours[( k << 1 ) + 1];
    }
    mrgCtx.mmvdUseAltHpelIf[currBaseNum] = mrgCtx.useAltHpelIf[k];
    currBaseNum++;

    if( currBaseNum == MMVD_BASE_MV_NUM )
    {
      break;
    }
  }
}
bool PU::getColocatedMVP(const PredictionUnit &pu, const RefPicList &eRefPicList, const Position &_pos, Mv& rcMv, const int &refIdx, bool sbFlag
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
  , int* targetRefIdx
#endif
)
{
  // don't perform MV compression when generally disabled or subPuMvp is used
  const unsigned scale = 4 * std::max<int>(1, 4 * AMVP_DECIMATION_FACTOR / 4);
  const unsigned mask  = ~( scale - 1 );

  const Position pos = Position{ PosType( _pos.x & mask ), PosType( _pos.y & mask ) };

  const Slice &slice = *pu.cs->slice;

  // use coldir.
  const Picture* const pColPic = slice.getRefPic(RefPicList(slice.isInterB() ? 1 - slice.getColFromL0Flag() : 0), slice.getColRefIdx());

  if( !pColPic )
  {
    return false;
  }

  // Check the position of colocated block is within a subpicture
  const SubPic &curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
  if (curSubPic.getTreatedAsPicFlag())
  {
    if (!curSubPic.isContainingPos(pos))
    {
      return false;
    }
  }
  RefPicList eColRefPicList = slice.getCheckLDC() ? eRefPicList : RefPicList(slice.getColFromL0Flag());

  const MotionInfo& mi = pColPic->cs->getMotionInfo( pos );

  if( !mi.isInter )
  {
    return false;
  }
  if (mi.isIBCmot)
  {
    return false;
  }
  if (CU::isIBC(*pu.cu))
  {
    return false;
  }
  int iColRefIdx = mi.refIdx[eColRefPicList];

  if (sbFlag && !slice.getCheckLDC())
  {
    eColRefPicList = eRefPicList;
    iColRefIdx = mi.refIdx[eColRefPicList];
    if (iColRefIdx < 0)
    {
      return false;
    }
  }
  else
  {
    if (iColRefIdx < 0)
    {
      eColRefPicList = RefPicList(1 - eColRefPicList);
      iColRefIdx = mi.refIdx[eColRefPicList];

      if (iColRefIdx < 0)
      {
        return false;
      }
    }
  }

  const Slice *pColSlice = nullptr;

  for( const auto s : pColPic->slices )
  {
    if( s->getIndependentSliceIdx() == mi.sliceIdx )
    {
      pColSlice = s;
      break;
    }
  }

  CHECK( pColSlice == nullptr, "Slice segment not found" );

  const Slice &colSlice = *pColSlice;
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
  if (targetRefIdx != nullptr)
  {
    *targetRefIdx = slice.getImRefIdx(mi.sliceIdx, eColRefPicList, eRefPicList, iColRefIdx);
    if (*targetRefIdx == -1)
    {
      return false;
    }

    const bool bIsCurrRefLongTerm = slice.getRefPic(eRefPicList, *targetRefIdx)->longTerm;
    const bool bIsColRefLongTerm = colSlice.getIsUsedAsLongTerm(eColRefPicList, iColRefIdx);
    if (bIsCurrRefLongTerm != bIsColRefLongTerm)
    {
      return false;
    }

    // Scale the vector.
    Mv cColMv = mi.mv[eColRefPicList];
    cColMv.setHor(roundMvComp(cColMv.getHor()));
    cColMv.setVer(roundMvComp(cColMv.getVer()));

    if (bIsCurrRefLongTerm /*|| bIsColRefLongTerm*/)
    {
      rcMv = cColMv;
      rcMv.clipToStorageBitDepth();
    }
    else
    {
      const int currPOC = slice.getPOC();
      const int colPOC = colSlice.getPOC();
      const int colRefPOC = colSlice.getRefPOC(eColRefPicList, iColRefIdx);
      const int currRefPOC = slice.getRefPic(eRefPicList, *targetRefIdx)->getPOC();
      const int distscale = xGetDistScaleFactor(currPOC, currRefPOC, colPOC, colRefPOC);

      if (distscale == 4096)
      {
        rcMv = cColMv;
        rcMv.clipToStorageBitDepth();
      }
      else
      {
        rcMv = cColMv.scaleMv(distscale);
      }
    }

    return true;
  }
  else
  {
#endif
  const bool bIsCurrRefLongTerm = slice.getRefPic(eRefPicList, refIdx)->longTerm;
  const bool bIsColRefLongTerm  = colSlice.getIsUsedAsLongTerm(eColRefPicList, iColRefIdx);

  if (bIsCurrRefLongTerm != bIsColRefLongTerm)
  {
    return false;
  }


  // Scale the vector.
  Mv cColMv = mi.mv[eColRefPicList];
  cColMv.setHor(roundMvComp(cColMv.getHor()));
  cColMv.setVer(roundMvComp(cColMv.getVer()));

  if (bIsCurrRefLongTerm /*|| bIsColRefLongTerm*/)
  {
    rcMv = cColMv;
    rcMv.clipToStorageBitDepth();
  }
  else
  {
    const int currPOC    = slice.getPOC();
    const int colPOC     = colSlice.getPOC();
    const int colRefPOC  = colSlice.getRefPOC(eColRefPicList, iColRefIdx);
    const int currRefPOC = slice.getRefPic(eRefPicList, refIdx)->getPOC();
    const int distscale  = xGetDistScaleFactor(currPOC, currRefPOC, colPOC, colRefPOC);

    if (distscale == 4096)
    {
      rcMv = cColMv;
      rcMv.clipToStorageBitDepth();
    }
    else
    {
      rcMv = cColMv.scaleMv(distscale);
    }
  }

  return true;
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
  }
#endif
}

bool PU::isDiffMER(const Position &pos1, const Position &pos2, const unsigned plevel)
{
  const unsigned xN = pos1.x;
  const unsigned yN = pos1.y;
  const unsigned xP = pos2.x;
  const unsigned yP = pos2.y;

  if ((xN >> plevel) != (xP >> plevel))
  {
    return true;
  }

  if ((yN >> plevel) != (yP >> plevel))
  {
    return true;
  }

  return false;
}

bool PU::isAddNeighborMv(const Mv& currMv, Mv* neighborMvs, int numNeighborMv)
{
  bool existed = false;
  for (uint32_t cand = 0; cand < numNeighborMv && !existed; cand++)
  {
    if (currMv == neighborMvs[cand])
    {
      existed = true;
    }
  }

  if (!existed)
  {
    return true;
  }
  else
  {
    return false;
  }
}

void PU::getIbcMVPsEncOnly(PredictionUnit &pu, Mv* mvPred, int& nbPred)
{
  const PreCalcValues   &pcv = *pu.cs->pcv;
  const int  cuWidth = pu.blocks[COMPONENT_Y].width;
  const int  cuHeight = pu.blocks[COMPONENT_Y].height;
  const int  log2UnitWidth = floorLog2(pcv.minCUWidth);
  const int  log2UnitHeight = floorLog2(pcv.minCUHeight);
  const int  totalAboveUnits = (cuWidth >> log2UnitWidth) + 1;
  const int  totalLeftUnits = (cuHeight >> log2UnitHeight) + 1;

  nbPred = 0;
  Position posLT = pu.Y().topLeft();

#if JVET_Z0118_GDR  
  bool isClean = pu.cs->isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

  // above-left
  const PredictionUnit *aboveLeftPU = pu.cs->getPURestricted(posLT.offset(-1, -1), pu, CHANNEL_TYPE_LUMA);
  if (aboveLeftPU && CU::isIBC(*aboveLeftPU->cu))
  {
    if (isAddNeighborMv(aboveLeftPU->bv, mvPred, nbPred))
    {
      mvPred[nbPred++] = aboveLeftPU->bv;
    }
  }

  // above neighbors
  for (uint32_t dx = 0; dx < totalAboveUnits && nbPred < IBC_NUM_CANDIDATES; dx++)
  {
    const PredictionUnit* tmpPU = pu.cs->getPURestricted(posLT.offset((dx << log2UnitWidth), -1), pu, CHANNEL_TYPE_LUMA);
    if (tmpPU && CU::isIBC(*tmpPU->cu))
    {
      if (isAddNeighborMv(tmpPU->bv, mvPred, nbPred))
      {
        mvPred[nbPred++] = tmpPU->bv;
      }
    }
  }

  // left neighbors
  for (uint32_t dy = 0; dy < totalLeftUnits && nbPred < IBC_NUM_CANDIDATES; dy++)
  {
    const PredictionUnit* tmpPU = pu.cs->getPURestricted(posLT.offset(-1, (dy << log2UnitHeight)), pu, CHANNEL_TYPE_LUMA);
    if (tmpPU && CU::isIBC(*tmpPU->cu))
    {
      if (isAddNeighborMv(tmpPU->bv, mvPred, nbPred))
      {
        mvPred[nbPred++] = tmpPU->bv;
      }
    }
  }

#if JVET_Z0118_GDR  
  size_t numAvaiCandInLUT = isClean ? pu.cs->motionLut.lutIbc1.size() : pu.cs->motionLut.lutIbc0.size();
#else
  size_t numAvaiCandInLUT = pu.cs->motionLut.lutIbc.size();
#endif

  for (uint32_t cand = 0; cand < numAvaiCandInLUT && nbPred < IBC_NUM_CANDIDATES; cand++)
  {
#if JVET_Z0118_GDR
    MotionInfo neibMi = isClean ? pu.cs->motionLut.lutIbc1[cand] : pu.cs->motionLut.lutIbc0[cand];
#else
    MotionInfo neibMi = pu.cs->motionLut.lutIbc[cand];
#endif
    if (isAddNeighborMv(neibMi.bv, mvPred, nbPred))
    {
      mvPred[nbPred++] = neibMi.bv;
    }
  }

  bool isBvCandDerived[IBC_NUM_CANDIDATES];
  ::memset(isBvCandDerived, false, IBC_NUM_CANDIDATES);

  int curNbPred = nbPred;
  if (curNbPred < IBC_NUM_CANDIDATES)
  {
    do
    {
      curNbPred = nbPred;
      for (uint32_t idx = 0; idx < curNbPred && nbPred < IBC_NUM_CANDIDATES; idx++)
      {
        if (!isBvCandDerived[idx])
        {
          Mv derivedBv;
          if (getDerivedBV(pu, mvPred[idx], derivedBv))
          {
            if (isAddNeighborMv(derivedBv, mvPred, nbPred))
            {
              mvPred[nbPred++] = derivedBv;
            }
          }
          isBvCandDerived[idx] = true;
        }
      }
    } while (nbPred > curNbPred && nbPred < IBC_NUM_CANDIDATES);
  }
}

bool PU::getDerivedBV(PredictionUnit &pu, const Mv& currentMv, Mv& derivedMv)
{
  int   cuPelX = pu.lumaPos().x;
  int   cuPelY = pu.lumaPos().y;
  int rX = cuPelX + currentMv.getHor();
  int rY = cuPelY + currentMv.getVer();
  int offsetX = currentMv.getHor();
  int offsetY = currentMv.getVer();

  if( rX < 0 || rY < 0 || rX >= pu.cs->slice->getPPS()->getPicWidthInLumaSamples() || rY >= pu.cs->slice->getPPS()->getPicHeightInLumaSamples() )
  {
    return false;
  }

  const PredictionUnit *neibRefPU = NULL;
  neibRefPU = pu.cs->getPURestricted(pu.lumaPos().offset(offsetX, offsetY), pu, CHANNEL_TYPE_LUMA);

  bool isIBC = (neibRefPU) ? CU::isIBC(*neibRefPU->cu) : 0;
  if (isIBC)
  {
    derivedMv = neibRefPU->bv;
    derivedMv += currentMv;
  }
  return isIBC;
}

/**
 * Constructs a list of candidates for IBC AMVP (See specification, section "Derivation process for motion vector predictor candidates")
 */
#if JVET_Z0084_IBC_TM && TM_AMVP
void PU::fillIBCMvpCand(PredictionUnit &pu, AMVPInfo &amvpInfo, InterPrediction* pcInter)
#else
void PU::fillIBCMvpCand(PredictionUnit &pu, AMVPInfo &amvpInfo)
#endif
{
  AMVPInfo *pInfo = &amvpInfo;
  pInfo->numCand = 0;

  MergeCtx mergeCtx;
#if JVET_Z0084_IBC_TM && TM_AMVP
  pInfo->maxSimilarityThreshold = (pu.cs->sps->getUseDMVDMode() && pcInter) ? PU::getTMMvdThreshold(pu) : 1;

  pu.tmMergeFlag = true;
#if JVET_Z0075_IBC_HMVP_ENLARGE
  PU::getIBCMergeCandidates(pu, mergeCtx, pu.cs->sps->getMaxNumIBCMergeCand());
#else
  PU::getIBCMergeCandidates(pu, mergeCtx);
#endif
  pu.tmMergeFlag = false;

  int candIdx = 0;
  while ((pInfo->numCand < AMVP_MAX_NUM_CANDS_MEM) && (candIdx < mergeCtx.numValidMergeCand))
  {
    pInfo->mvCand[pInfo->numCand] = mergeCtx.mvFieldNeighbours[candIdx << 1].mv;
    pInfo->mvCand[pInfo->numCand].roundIbcPrecInternal2Amvr(pu.cu->imv);
    if (!pInfo->xCheckSimilarMotion(pInfo->numCand))
    {
      pInfo->numCand++;
    }
    candIdx++;
  }

  if (pu.cs->sps->getUseDMVDMode() && pcInter && pInfo->numCand > 0)
  {
    struct AMVPSort
    {
      Mv AMVPCand;
      Distortion cost;
    };
    AMVPSort temp;
    std::vector<AMVPSort> input;
    const auto CostIncSort = [](const AMVPSort &x, const AMVPSort &y) { return x.cost < y.cost; };
    Distortion tmCost = std::numeric_limits<Distortion>::max();

    for (int candIdx = 0; candIdx < pInfo->numCand; ++candIdx)
    {
      tmCost = pcInter->deriveTMMv(pu, true, std::numeric_limits<Distortion>::max(), REF_PIC_LIST_0, MAX_NUM_REF, TM_MAX_NUM_OF_ITERATIONS, pInfo->mvCand[candIdx]);
      pInfo->mvCand[candIdx].roundIbcPrecInternal2Amvr(pu.cu->imv);
      temp.AMVPCand = pInfo->mvCand[candIdx];
      temp.cost     = tmCost;
      input.push_back(temp);
    }

    stable_sort(input.begin(), input.end(), CostIncSort);
    for (int candIdx = 0; candIdx < pInfo->numCand; ++candIdx)
    {
      pInfo->mvCand[candIdx] = input.at(candIdx).AMVPCand;
    }
  }

  if (pInfo->numCand > AMVP_MAX_NUM_CANDS)
  {
    pInfo->numCand = AMVP_MAX_NUM_CANDS;
  }

  while (pInfo->numCand < AMVP_MAX_NUM_CANDS)
  {
    pInfo->mvCand[pInfo->numCand++] = Mv(0, 0);
  }
#else
  PU::getIBCMergeCandidates(pu, mergeCtx, AMVP_MAX_NUM_CANDS - 1);
  int candIdx = 0;
  while (pInfo->numCand < AMVP_MAX_NUM_CANDS)
  {
    pInfo->mvCand[pInfo->numCand] = mergeCtx.mvFieldNeighbours[(candIdx << 1) + 0].mv;;
    pInfo->numCand++;
    candIdx++;
  }

  for (Mv &mv : pInfo->mvCand)
  {
    mv.roundIbcPrecInternal2Amvr(pu.cu->imv);
  }
#endif
}


/** Constructs a list of candidates for AMVP (See specification, section "Derivation process for motion vector predictor candidates")
* \param uiPartIdx
* \param uiPartAddr
* \param eRefPicList
* \param iRefIdx
* \param pInfo
*/
void PU::fillMvpCand(PredictionUnit &pu, const RefPicList &eRefPicList, const int &refIdx, AMVPInfo &amvpInfo
#if TM_AMVP
                   , InterPrediction* interPred
#endif
)
{
  CodingStructure &cs = *pu.cs;
#if TM_AMVP
#if JVET_Y0128_NON_CTC
  interPred = PU::checkTmEnableCondition(pu.cs->sps, pu.cs->pps, pu.cu->slice->getRefPic(eRefPicList, refIdx)) ? interPred : nullptr;
#else
  interPred = pu.cu->cs->sps->getUseDMVDMode() ? interPred : nullptr;
#endif
#if JVET_Z0054_BLK_REF_PIC_REORDER // AmvpList might be already derived during ref. reordering
  if (interPred != nullptr && interPred->readTplAmvpBuffer(amvpInfo, *pu.cu, eRefPicList, refIdx))
#else
  if (cs.pcv->isEncoder && interPred != nullptr && interPred->readTplAmvpBuffer(amvpInfo, *pu.cu, eRefPicList, refIdx))
#endif
  {
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
    if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
    {
      amvpInfo.numCand = 2;
      amvpInfo.mvCand[1] = amvpInfo.mvCand[0];
    }
#endif
    return;
  }
#endif
  AMVPInfo *pInfo = &amvpInfo;

  pInfo->numCand = 0;
#if TM_AMVP
  pInfo->maxStorageSize         = interPred != nullptr ? REGULAR_AMVP_MAX_NUM_CANDS : AMVP_MAX_NUM_CANDS;
  pInfo->maxSimilarityThreshold = interPred != nullptr ? 1 : 0;
#endif

  if (refIdx < 0)
  {
    return;
  }

  //-- Get Spatial MV
  Position posLT = pu.Y().topLeft();
  Position posRT = pu.Y().topRight();
  Position posLB = pu.Y().bottomLeft();

  {
    bool bAdded = addMVPCandUnscaled( pu, eRefPicList, refIdx, posLB, MD_BELOW_LEFT, *pInfo );

    if( !bAdded )
    {
      bAdded = addMVPCandUnscaled( pu, eRefPicList, refIdx, posLB, MD_LEFT, *pInfo );

    }
  }

  // Above predictor search
  {
    bool bAdded = addMVPCandUnscaled( pu, eRefPicList, refIdx, posRT, MD_ABOVE_RIGHT, *pInfo );

    if( !bAdded )
    {
      bAdded = addMVPCandUnscaled( pu, eRefPicList, refIdx, posRT, MD_ABOVE, *pInfo );

      if( !bAdded )
      {
        addMVPCandUnscaled( pu, eRefPicList, refIdx, posLT, MD_ABOVE_LEFT, *pInfo );
      }
    }
  }

  for( int i = 0; i < pInfo->numCand; i++ )
  {
    pInfo->mvCand[i].roundTransPrecInternal2Amvr(pu.cu->imv);
  }

  if( pInfo->numCand == 2 )
  {
    if( pInfo->mvCand[0] == pInfo->mvCand[1] )
    {
      pInfo->numCand = 1;
    }
  }

#if TM_AMVP
  if (cs.picHeader->getEnableTMVPFlag() && pInfo->numCand < pInfo->maxStorageSize && (pu.lumaSize().width + pu.lumaSize().height > 12))
#else
  if (cs.picHeader->getEnableTMVPFlag() && pInfo->numCand < AMVP_MAX_NUM_CANDS && (pu.lumaSize().width + pu.lumaSize().height > 12))
#endif
  {
    // Get Temporal Motion Predictor
    const int refIdx_Col = refIdx;

    Position posRB = pu.Y().bottomRight().offset(-3, -3);

    const PreCalcValues& pcv = *cs.pcv;

    Position posC0;
    bool C0Avail = false;
    Position posC1 = pu.Y().center();
    Mv cColMv;

    bool boundaryCond = ((posRB.x + pcv.minCUWidth) < pcv.lumaWidth) && ((posRB.y + pcv.minCUHeight) < pcv.lumaHeight);
    const SubPic &curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
    if (curSubPic.getTreatedAsPicFlag())
    {
      boundaryCond = ((posRB.x + pcv.minCUWidth) <= curSubPic.getSubPicRight() &&
                      (posRB.y + pcv.minCUHeight) <= curSubPic.getSubPicBottom());
    }
    if (boundaryCond)
    {
      int posYInCtu = posRB.y & pcv.maxCUHeightMask;
      if (posYInCtu + 4 < pcv.maxCUHeight)
      {
        posC0 = posRB.offset(4, 4);
        C0Avail = true;
      }
    }
    if ( ( C0Avail && getColocatedMVP( pu, eRefPicList, posC0, cColMv, refIdx_Col, false ) ) || getColocatedMVP( pu, eRefPicList, posC1, cColMv, refIdx_Col, false ) )
    {
      cColMv.roundTransPrecInternal2Amvr(pu.cu->imv);
#if TM_AMVP
      pInfo->mvCand[pInfo->numCand] = cColMv;
      if (!pInfo->xCheckSimilarMotion(pInfo->numCand))
      {
        pInfo->numCand++;
      }
#else
      pInfo->mvCand[pInfo->numCand++] = cColMv;
#endif
    }
  }

#if TM_AMVP
  // Non-adjacent candidates
  if (pInfo->numCand < pInfo->maxStorageSize && interPred != nullptr)
  {
    const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
    MotionInfo miNeighbor;
    int offsetX = 0;
    int offsetY = 0;
    const int iNACANDIDATE_NUM[4] = { 3, 5, 5, 5 };
    const int idxMap[4][5] = { { 0, 1, 4 },{ 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4 },{ 0, 1, 2, 3, 4 } };

    for (int iDistanceIndex = 0; iDistanceIndex < NADISTANCE_LEVEL && pInfo->numCand < pInfo->maxStorageSize; iDistanceIndex++)
    {
      const int iNADistanceHor = pu.Y().width  * (iDistanceIndex + 1);
      const int iNADistanceVer = pu.Y().height * (iDistanceIndex + 1);

      for (int NASPIdx = 0; NASPIdx < iNACANDIDATE_NUM[iDistanceIndex] && pInfo->numCand < pInfo->maxStorageSize; NASPIdx++)
      {
        switch (idxMap[iDistanceIndex][NASPIdx])
        {
        case 0: offsetX = -iNADistanceHor - 1;               offsetY = pu.Y().height + iNADistanceVer - 1; break;
        case 1: offsetX = pu.Y().width + iNADistanceHor - 1; offsetY = -iNADistanceVer - 1;                break;
        case 2: offsetX = pu.Y().width >> 1;                 offsetY = -iNADistanceVer - 1;                break;
        case 3: offsetX = -iNADistanceHor - 1;               offsetY = pu.Y().height >> 1;                 break;
        case 4: offsetX = -iNADistanceHor - 1;               offsetY = -iNADistanceVer - 1;                break;
        default: CHECK(true, "Unknown index for non-adjacent AMVP candidate"); break;
        }

        const PredictionUnit *puNonAdjacent = cs.getPURestricted(posLT.offset(offsetX, offsetY), pu, pu.chType);

#if JVET_Y0065_GPM_INTRA
        bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu) && puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY)).isInter;
#else
        bool isAvailableNonAdjacent = puNonAdjacent && isDiffMER(pu.lumaPos(), posLT.offset(offsetX, offsetY), plevel) && CU::isInter(*puNonAdjacent->cu);
#endif

        if (isAvailableNonAdjacent)
        {
          miNeighbor = puNonAdjacent->getMotionInfo(posLT.offset(offsetX, offsetY));

          const int        currRefPOC = cs.slice->getRefPic(eRefPicList, refIdx)->getPOC();
          const RefPicList eRefPicList2nd = (eRefPicList == REF_PIC_LIST_0) ? REF_PIC_LIST_1 : REF_PIC_LIST_0;

          for (int predictorSource = 0; predictorSource < 2; predictorSource++) // examine the indicated reference picture list, then if not available, examine the other list.
          {
            const RefPicList eRefPicListIndex = (predictorSource == 0) ? eRefPicList : eRefPicList2nd;
            const int        neibRefIdx = miNeighbor.refIdx[eRefPicListIndex];

            if (neibRefIdx >= 0 && currRefPOC == cs.slice->getRefPOC(eRefPicListIndex, neibRefIdx))
            {
              miNeighbor.mv[eRefPicListIndex].roundTransPrecInternal2Amvr(pu.cu->imv);
              pInfo->mvCand[(pInfo->numCand)] = miNeighbor.mv[eRefPicListIndex];
              if (!pInfo->xCheckSimilarMotion(pInfo->numCand))
              {
                pInfo->numCand++;
              }
            }
            if (pInfo->numCand > pInfo->maxStorageSize)
            {
              pInfo->numCand = pInfo->maxStorageSize;
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < pInfo->numCand; i++)
  {
    pInfo->mvCand[i].roundTransPrecInternal2Amvr(pu.cu->imv);
  }
#endif

#if TM_AMVP
  if (pInfo->numCand < pInfo->maxStorageSize)
#else
  if (pInfo->numCand < AMVP_MAX_NUM_CANDS)
#endif
  {
    const int currRefPOC = cs.slice->getRefPic(eRefPicList, refIdx)->getPOC();
    addAMVPHMVPCand(pu, eRefPicList, currRefPOC, *pInfo);
  }
#if TM_AMVP
  if (pInfo->numCand > pInfo->maxStorageSize)
  {
    pInfo->numCand = pInfo->maxStorageSize;
  }
#else
  if (pInfo->numCand > AMVP_MAX_NUM_CANDS)
  {
    pInfo->numCand = AMVP_MAX_NUM_CANDS;
  }
#endif

#if TM_AMVP
  if (pu.cu->cs->sps->getUseDMVDMode() && interPred != nullptr && pInfo->numCand > 0)
  {
    struct AMVPSort
    {
      Mv AMVPCand;
      Distortion cost;
    };
    AMVPSort temp;
    std::vector<AMVPSort> input;
    const auto CostIncSort = [](const AMVPSort &x, const AMVPSort &y) { return x.cost < y.cost; };
    Distortion tmCost[REGULAR_AMVP_MAX_NUM_CANDS];
    for (int i = 0; i < REGULAR_AMVP_MAX_NUM_CANDS; i++)
    {
      tmCost[i] = std::numeric_limits<Distortion>::max();
    }

    for (int candIdx = 0; candIdx < pInfo->numCand; ++candIdx)
    {
      tmCost[candIdx] = interPred->deriveTMMv(pu, true, std::numeric_limits<Distortion>::max(), eRefPicList, refIdx, 0, pInfo->mvCand[candIdx]);
      temp.AMVPCand = pInfo->mvCand[candIdx];
      temp.cost = tmCost[candIdx];
      input.push_back(temp);
    }

    stable_sort(input.begin(), input.end(), CostIncSort);
    for (int candIdx = 0; candIdx < pInfo->numCand; ++candIdx)
    {
      pInfo->mvCand[candIdx] = input.at(candIdx).AMVPCand;
      tmCost[candIdx] = input.at(candIdx).cost;
    }

    pInfo->numCand = 1;
    interPred->deriveTMMv(pu, true, tmCost[0], eRefPicList, refIdx, TM_MAX_NUM_OF_ITERATIONS, pInfo->mvCand[0]);
  }

  while (pInfo->numCand < pInfo->maxStorageSize)
#else
  while (pInfo->numCand < AMVP_MAX_NUM_CANDS)
#endif
  {
    pInfo->mvCand[pInfo->numCand] = Mv( 0, 0 );
    pInfo->numCand++;
  }
#if TM_AMVP
  pInfo->numCand = (pu.cu->cs->sps->getUseDMVDMode() && interPred != nullptr ? 1 : pInfo->numCand);
#endif

  for (Mv &mv : pInfo->mvCand)
  {
    mv.roundTransPrecInternal2Amvr(pu.cu->imv);
  }

#if TM_AMVP
#if JVET_Z0054_BLK_REF_PIC_REORDER
  if (interPred != nullptr)
#else
  if (cs.pcv->isEncoder && interPred != nullptr)
#endif
  {
    interPred->writeTplAmvpBuffer(*pInfo, *pu.cu, eRefPicList, refIdx);
  }
#endif
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    if (amvpInfo.numCand == 1)
    {
      amvpInfo.numCand = 2;
      amvpInfo.mvCand[1] = amvpInfo.mvCand[0];
    }
    else
    {
      CHECK(amvpInfo.numCand != 2, "this is not possible");
      amvpInfo.numCand = 3;
      amvpInfo.mvCand[2] = amvpInfo.mvCand[1];
      amvpInfo.mvCand[1] = amvpInfo.mvCand[0];
    }
  }
#endif
}

bool PU::addAffineMVPCandUnscaled( const PredictionUnit &pu, const RefPicList &refPicList, const int &refIdx, const Position &pos, const MvpDir &dir, AffineAMVPInfo &affiAMVPInfo 
#if JVET_Z0139_HIST_AFF
  , int aiNeibeInherited[5]
#endif
)
{
  CodingStructure &cs = *pu.cs;
  const PredictionUnit *neibPU = NULL;
  Position neibPos;

  switch ( dir )
  {
  case MD_LEFT:
    neibPos = pos.offset( -1, 0 );
    break;
  case MD_ABOVE:
    neibPos = pos.offset( 0, -1 );
    break;
  case MD_ABOVE_RIGHT:
    neibPos = pos.offset( 1, -1 );
    break;
  case MD_BELOW_LEFT:
    neibPos = pos.offset( -1, 1 );
    break;
  case MD_ABOVE_LEFT:
    neibPos = pos.offset( -1, -1 );
    break;
  default:
    break;
  }

  neibPU = cs.getPURestricted( neibPos, pu, pu.chType );

  if (neibPU == NULL || !CU::isInter(*neibPU->cu) || !neibPU->cu->affine || neibPU->mergeType != MRG_TYPE_DEFAULT_N)
  {
    return false;
  }

  Mv outputAffineMv[3];
  const MotionInfo& neibMi = neibPU->getMotionInfo( neibPos );

  const int        currRefPOC = cs.slice->getRefPic( refPicList, refIdx )->getPOC();
  const RefPicList refPicList2nd = (refPicList == REF_PIC_LIST_0) ? REF_PIC_LIST_1 : REF_PIC_LIST_0;

  for ( int predictorSource = 0; predictorSource < 2; predictorSource++ ) // examine the indicated reference picture list, then if not available, examine the other list.
  {
    const RefPicList eRefPicListIndex = (predictorSource == 0) ? refPicList : refPicList2nd;
    const int        neibRefIdx = neibMi.refIdx[eRefPicListIndex];

    if ( ((neibPU->interDir & (eRefPicListIndex + 1)) == 0) || pu.cu->slice->getRefPOC( eRefPicListIndex, neibRefIdx ) != currRefPOC )
    {
      continue;
    }

    xInheritedAffineMv( pu, neibPU, eRefPicListIndex, outputAffineMv );
    outputAffineMv[0].roundAffinePrecInternal2Amvr(pu.cu->imv);
    outputAffineMv[1].roundAffinePrecInternal2Amvr(pu.cu->imv);
    affiAMVPInfo.mvCandLT[affiAMVPInfo.numCand] = outputAffineMv[0];
    affiAMVPInfo.mvCandRT[affiAMVPInfo.numCand] = outputAffineMv[1];
    if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
    {
      outputAffineMv[2].roundAffinePrecInternal2Amvr(pu.cu->imv);
      affiAMVPInfo.mvCandLB[affiAMVPInfo.numCand] = outputAffineMv[2];
    }
#if JVET_Z0139_HIST_AFF
    if (!checkLastAffineAMVPCandRedundancy(pu, affiAMVPInfo))
    {
      continue;
    }
#endif
    affiAMVPInfo.numCand++;
#if JVET_Z0139_HIST_AFF
    aiNeibeInherited[dir] = 1;
#endif
    return true;
  }

  return false;
}

void PU::xInheritedAffineMv( const PredictionUnit &pu, const PredictionUnit* puNeighbour, RefPicList eRefPicList, Mv rcMv[3] )
{
  int posNeiX = puNeighbour->Y().pos().x;
  int posNeiY = puNeighbour->Y().pos().y;
  int posCurX = pu.Y().pos().x;
  int posCurY = pu.Y().pos().y;

  int neiW = puNeighbour->Y().width;
  int curW = pu.Y().width;
  int neiH = puNeighbour->Y().height;
  int curH = pu.Y().height;

  Mv mvLT, mvRT, mvLB;
  mvLT = puNeighbour->mvAffi[eRefPicList][0];
  mvRT = puNeighbour->mvAffi[eRefPicList][1];
  mvLB = puNeighbour->mvAffi[eRefPicList][2];

#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  bool isTopCtuBoundary = false;
  if ( (posNeiY + neiH) % pu.cs->sps->getCTUSize() == 0 && (posNeiY + neiH) == posCurY )
  {
    // use bottom-left and bottom-right sub-block MVs for inheritance
    const Position posRB = puNeighbour->Y().bottomRight();
    const Position posLB = puNeighbour->Y().bottomLeft();
    mvLT = puNeighbour->getMotionInfo( posLB ).mv[eRefPicList];
    mvRT = puNeighbour->getMotionInfo( posRB ).mv[eRefPicList];
    posNeiY += neiH;
    isTopCtuBoundary = true;
  }
#endif

  int shift = MAX_CU_DEPTH;
  int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;

  iDMvHorX = (mvRT - mvLT).getHor() << (shift - floorLog2(neiW));
  iDMvHorY = (mvRT - mvLT).getVer() << (shift - floorLog2(neiW));
#if AFFINE_RM_CONSTRAINTS_AND_OPT
  if ( puNeighbour->cu->affineType == AFFINEMODEL_6PARAM )
#else
  if ( puNeighbour->cu->affineType == AFFINEMODEL_6PARAM && !isTopCtuBoundary )
#endif
  {
    iDMvVerX = (mvLB - mvLT).getHor() << (shift - floorLog2(neiH));
    iDMvVerY = (mvLB - mvLT).getVer() << (shift - floorLog2(neiH));
  }
  else
  {
    iDMvVerX = -iDMvHorY;
    iDMvVerY = iDMvHorX;
  }

  int iMvScaleHor = mvLT.getHor() << shift;
  int iMvScaleVer = mvLT.getVer() << shift;
  int horTmp, verTmp;

  // v0
  horTmp = iMvScaleHor + iDMvHorX * (posCurX - posNeiX) + iDMvVerX * (posCurY - posNeiY);
  verTmp = iMvScaleVer + iDMvHorY * (posCurX - posNeiX) + iDMvVerY * (posCurY - posNeiY);
  roundAffineMv( horTmp, verTmp, shift );
  rcMv[0].hor = horTmp;
  rcMv[0].ver = verTmp;
  rcMv[0].clipToStorageBitDepth();

  // v1
  horTmp = iMvScaleHor + iDMvHorX * (posCurX + curW - posNeiX) + iDMvVerX * (posCurY - posNeiY);
  verTmp = iMvScaleVer + iDMvHorY * (posCurX + curW - posNeiX) + iDMvVerY * (posCurY - posNeiY);
  roundAffineMv( horTmp, verTmp, shift );
  rcMv[1].hor = horTmp;
  rcMv[1].ver = verTmp;
  rcMv[1].clipToStorageBitDepth();

  // v2
#if !JVET_Z0139_HIST_AFF
  if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
#endif
  {
    horTmp = iMvScaleHor + iDMvHorX * (posCurX - posNeiX) + iDMvVerX * (posCurY + curH - posNeiY);
    verTmp = iMvScaleVer + iDMvHorY * (posCurX - posNeiX) + iDMvVerY * (posCurY + curH - posNeiY);
    roundAffineMv( horTmp, verTmp, shift );
    rcMv[2].hor = horTmp;
    rcMv[2].ver = verTmp;
    rcMv[2].clipToStorageBitDepth();
  }
}

#if JVET_Z0139_HIST_AFF
void PU::deriveAffineParametersFromMVs(const PredictionUnit& pu, const Mv acMvTemp[3], int* affinePara, EAffineModel affModel)
{
  Mv mvLT = acMvTemp[0];
  Mv mvRT = acMvTemp[1];
  Mv mvLB = acMvTemp[2];
  const int iBit = MAX_CU_DEPTH;
  int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;
  iDMvHorX = (mvRT - mvLT).getHor() << (iBit - g_aucLog2[pu.lwidth()]);
  iDMvHorY = (mvRT - mvLT).getVer() << (iBit - g_aucLog2[pu.lwidth()]);
  if (affModel == AFFINEMODEL_6PARAM)
  {
    iDMvVerX = (mvLB - mvLT).getHor() << (iBit - g_aucLog2[pu.lheight()]);
    iDMvVerY = (mvLB - mvLT).getVer() << (iBit - g_aucLog2[pu.lheight()]);
  }
  else
  {
    iDMvVerX = -iDMvHorY;
    iDMvVerY = iDMvHorX;
  }

  affinePara[0] = iDMvHorX;
  affinePara[1] = iDMvVerX;
  affinePara[2] = iDMvHorY;
  affinePara[3] = iDMvVerY;
}
void PU::deriveMVsFromAffineParameters(const PredictionUnit & pu, Mv rcMv[3], int* affinePara, const Mv & cBaseMv, const Position & cBasePos)
{
  int posCurX = pu.Y().pos().x;
  int posCurY = pu.Y().pos().y;

  int curW = pu.Y().width;
  int curH = pu.Y().height;


  int shift = MAX_CU_DEPTH;
  int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;

  iDMvHorX = affinePara[0];
  iDMvVerX = affinePara[1];
  iDMvHorY = affinePara[2];
  iDMvVerY = affinePara[3];

  int iMvScaleHor = cBaseMv.getHor() << shift;
  int iMvScaleVer = cBaseMv.getVer() << shift;

  int horTmp, verTmp;

  // v0
  horTmp = iMvScaleHor + iDMvHorX * (posCurX - cBasePos.x) + iDMvVerX * (posCurY - cBasePos.y);
  verTmp = iMvScaleVer + iDMvHorY * (posCurX - cBasePos.x) + iDMvVerY * (posCurY - cBasePos.y);
  roundAffineMv(horTmp, verTmp, shift);

  rcMv[0].hor = horTmp;
  rcMv[0].ver = verTmp;


  // v1
  horTmp = iMvScaleHor + iDMvHorX * (posCurX + curW - cBasePos.x) + iDMvVerX * (posCurY - cBasePos.y);
  verTmp = iMvScaleVer + iDMvHorY * (posCurX + curW - cBasePos.x) + iDMvVerY * (posCurY - cBasePos.y);
  roundAffineMv(horTmp, verTmp, shift);

  rcMv[1].hor = horTmp;
  rcMv[1].ver = verTmp;


  // v2

  horTmp = iMvScaleHor + iDMvHorX * (posCurX - cBasePos.x) + iDMvVerX * (posCurY + curH - cBasePos.y);
  verTmp = iMvScaleVer + iDMvHorY * (posCurX - cBasePos.x) + iDMvVerY * (posCurY + curH - cBasePos.y);
  roundAffineMv(horTmp, verTmp, shift);

  rcMv[2].hor = horTmp;
  rcMv[2].ver = verTmp;

}

void PU::deriveCenterMVFromAffineParameters(const PredictionUnit& pu, Mv& rcMv, int* affinePara, const Mv& cBaseMv, const Position& cBasePos)
{

  int posCenterX = pu.Y().pos().x + (pu.Y().width >> 1);
  int posCenterY = pu.Y().pos().y + (pu.Y().height >> 1);

  int shift = MAX_CU_DEPTH;
  int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;

  iDMvHorX = affinePara[0];
  iDMvVerX = affinePara[1];
  iDMvHorY = affinePara[2];
  iDMvVerY = affinePara[3];

  int iMvScaleHor = cBaseMv.getHor() << shift;
  int iMvScaleVer = cBaseMv.getVer() << shift;

  int horTmp, verTmp;

  horTmp = iMvScaleHor + iDMvHorX * (posCenterX - cBasePos.x) + iDMvVerX * (posCenterY - cBasePos.y);
  verTmp = iMvScaleVer + iDMvHorY * (posCenterX - cBasePos.x) + iDMvVerY * (posCenterY - cBasePos.y);
  roundAffineMv(horTmp, verTmp, shift);

  rcMv.hor = horTmp;
  rcMv.ver = verTmp;
}

void PU::storeAffParas(int* affinePara)
{
  const int MinV = -(1 << (AFF_PARA_STORE_BITS - 1));
  const int MaxV = (1 << (AFF_PARA_STORE_BITS - 1)) - 1;
  const int roundingAdd = (1 << AFF_PARA_SHIFT) >> 1;
  affinePara[0] = Clip3(MinV, MaxV, affinePara[0] >= 0 ? ((affinePara[0] + roundingAdd) >> AFF_PARA_SHIFT) : -((-affinePara[0] + roundingAdd) >> AFF_PARA_SHIFT));
  affinePara[1] = Clip3(MinV, MaxV, affinePara[1] >= 0 ? ((affinePara[1] + roundingAdd) >> AFF_PARA_SHIFT) : -((-affinePara[1] + roundingAdd) >> AFF_PARA_SHIFT));
  affinePara[2] = Clip3(MinV, MaxV, affinePara[2] >= 0 ? ((affinePara[2] + roundingAdd) >> AFF_PARA_SHIFT) : -((-affinePara[2] + roundingAdd) >> AFF_PARA_SHIFT));
  affinePara[3] = Clip3(MinV, MaxV, affinePara[3] >= 0 ? ((affinePara[3] + roundingAdd) >> AFF_PARA_SHIFT) : -((-affinePara[3] + roundingAdd) >> AFF_PARA_SHIFT));
}


void PU::xGetAffineMvFromLUT(AffineMotionInfo * affHistInfo, int rParameters[4])
{
  rParameters[0] = affHistInfo->oneSetAffineParameters[0] << AFF_PARA_SHIFT;
  rParameters[1] = affHistInfo->oneSetAffineParameters[1] << AFF_PARA_SHIFT;
  rParameters[2] = affHistInfo->oneSetAffineParameters[2] << AFF_PARA_SHIFT;
  rParameters[3] = affHistInfo->oneSetAffineParameters[3] << AFF_PARA_SHIFT;
}

void PU::xGetAffineMvFromLUT(short affineParameters[4], int rParameters[4])
{
  rParameters[0] = affineParameters[0] << AFF_PARA_SHIFT;
  rParameters[1] = affineParameters[1] << AFF_PARA_SHIFT;
  rParameters[2] = affineParameters[2] << AFF_PARA_SHIFT;
  rParameters[3] = affineParameters[3] << AFF_PARA_SHIFT;
}


// Add affine candidates from LUT
// return true to early terminate

bool PU::checkLastAffineMergeCandRedundancy(const PredictionUnit& pu, AffineMergeCtx& affMrgCtx)
{
  if (affMrgCtx.numValidMergeCand < 1)
  {
    return true;
  }

  const int CPMV_SIMILARITY_THREH = 1;
  const int PARA_SIMILARITY_THREH = 1;

  const int lastIdx = affMrgCtx.numValidMergeCand;
  if (affMrgCtx.mergeType[lastIdx] == MRG_TYPE_SUBPU_ATMVP)
  {
    return true;
  }

  for (int idx = 0; idx < affMrgCtx.numValidMergeCand; idx++)
  {
    if (affMrgCtx.mergeType[idx] == MRG_TYPE_SUBPU_ATMVP)
    {
      continue;
    }
    if (affMrgCtx.affineType[idx] != affMrgCtx.affineType[lastIdx])
    {
      continue;
    }
    CHECK(affMrgCtx.affineType[idx] != AFFINEMODEL_6PARAM && affMrgCtx.affineType[idx] != AFFINEMODEL_4PARAM, "Invalid parameter" );

    if ( affMrgCtx.interDirNeighbours[idx] != affMrgCtx.interDirNeighbours[lastIdx])
    {
      continue;
    }
    if (affMrgCtx.interDirNeighbours[idx] == 3 && affMrgCtx.BcwIdx[idx] != affMrgCtx.BcwIdx[lastIdx])
    {
      continue;
    }
#if INTER_LIC
    if (affMrgCtx.interDirNeighbours[idx] != 3 && affMrgCtx.LICFlags[idx] != affMrgCtx.LICFlags[lastIdx])
    {
      continue;
    }
#endif

    if ((affMrgCtx.interDirNeighbours[lastIdx] & 1) != 0)
    {
      if (affMrgCtx.mvFieldNeighbours[(idx << 1)][0].refIdx != affMrgCtx.mvFieldNeighbours[(lastIdx << 1)][0].refIdx)
      {
        continue;
      }
      Mv acMvTemp[3];
      int affinePara[4], affineParaLast[4];
      acMvTemp[0] = affMrgCtx.mvFieldNeighbours[(idx << 1)][0].mv;
      acMvTemp[1] = affMrgCtx.mvFieldNeighbours[(idx << 1)][1].mv;
      acMvTemp[2] = affMrgCtx.mvFieldNeighbours[(idx << 1)][2].mv;
      deriveAffineParametersFromMVs(pu, acMvTemp,  affinePara, affMrgCtx.affineType[idx]);
      acMvTemp[0] = affMrgCtx.mvFieldNeighbours[(lastIdx << 1)][0].mv;
      acMvTemp[1] = affMrgCtx.mvFieldNeighbours[(lastIdx << 1)][1].mv;
      acMvTemp[2] = affMrgCtx.mvFieldNeighbours[(lastIdx << 1)][2].mv;
      deriveAffineParametersFromMVs(pu, acMvTemp, affineParaLast, affMrgCtx.affineType[idx]);

      if (
        abs(affMrgCtx.mvFieldNeighbours[(idx << 1)][0].mv.getHor() - affMrgCtx.mvFieldNeighbours[(lastIdx << 1)][0].mv.getHor()) > CPMV_SIMILARITY_THREH ||      
        abs(affMrgCtx.mvFieldNeighbours[(idx << 1)][0].mv.getVer() - affMrgCtx.mvFieldNeighbours[(lastIdx << 1)][0].mv.getVer()) > CPMV_SIMILARITY_THREH ||
        abs(affinePara[0] - affineParaLast[0]) > PARA_SIMILARITY_THREH ||
        abs(affinePara[2] - affineParaLast[2]) > PARA_SIMILARITY_THREH        
        )
      {
        continue;
      }


      if (affMrgCtx.affineType[idx] == AFFINEMODEL_6PARAM)
      {
        if( abs(affinePara[1] - affineParaLast[1]) > PARA_SIMILARITY_THREH ||
            abs(affinePara[3] - affineParaLast[3]) > PARA_SIMILARITY_THREH )
        {
          continue;
        }
      }
    }

    if ((affMrgCtx.interDirNeighbours[lastIdx] & 2) != 0)
    {
      if (affMrgCtx.mvFieldNeighbours[(idx << 1) + 1][0].refIdx != affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 1][0].refIdx)
      {
        continue;
      }
      Mv acMvTemp[3];
      int affinePara[4], affineParaLast[4];
      acMvTemp[0] = affMrgCtx.mvFieldNeighbours[(idx << 1) + 1][0].mv;
      acMvTemp[1] = affMrgCtx.mvFieldNeighbours[(idx << 1) + 1][1].mv;
      acMvTemp[2] = affMrgCtx.mvFieldNeighbours[(idx << 1) + 1][2].mv;
      deriveAffineParametersFromMVs(pu, acMvTemp, affinePara, affMrgCtx.affineType[idx]);
      acMvTemp[0] = affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 1][0].mv;
      acMvTemp[1] = affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 1][1].mv;
      acMvTemp[2] = affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 1][2].mv;
      deriveAffineParametersFromMVs(pu, acMvTemp, affineParaLast, affMrgCtx.affineType[idx]);

      if (
        abs(affMrgCtx.mvFieldNeighbours[(idx << 1) + 1][0].mv.getHor() - affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 1][0].mv.getHor()) > CPMV_SIMILARITY_THREH ||
        abs(affMrgCtx.mvFieldNeighbours[(idx << 1) + 1][0].mv.getVer() - affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 1][0].mv.getVer()) > CPMV_SIMILARITY_THREH ||
        abs(affinePara[0] - affineParaLast[0]) > PARA_SIMILARITY_THREH ||
        abs(affinePara[2] - affineParaLast[2]) > PARA_SIMILARITY_THREH)
      {
        continue;
      }

      if (affMrgCtx.affineType[idx] == AFFINEMODEL_6PARAM)
      {
        if (abs(affinePara[1] - affineParaLast[1]) > PARA_SIMILARITY_THREH ||
            abs(affinePara[3] - affineParaLast[3]) > PARA_SIMILARITY_THREH)
        {
          continue;
        }
      }
    }

    for (int mvNum = 0; mvNum < 3; mvNum++)
    {
      affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 0][mvNum].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(lastIdx << 1) + 1][mvNum].setMvField(Mv(), -1);
    }
    affMrgCtx.interDirNeighbours[lastIdx] = 0;
    affMrgCtx.affineType[lastIdx] = AFFINEMODEL_4PARAM;
    affMrgCtx.mergeType[lastIdx] = MRG_TYPE_DEFAULT_N;
    affMrgCtx.BcwIdx[lastIdx] = BCW_DEFAULT;
#if INTER_LIC
    affMrgCtx.LICFlags[lastIdx] = false;
#endif
    return false;
  }
  return true;
}

bool PU::addMergeHMVPCandFromAffModel(const PredictionUnit& pu, MergeCtx& mrgCtx, const int& mrgCandIdx, int& cnt
#if TM_MRG
  , const uint32_t mvdSimilarityThresh
#endif
)
{
  if ((pu.cu->slice->getPicHeader()->getMvdL1ZeroFlag() || pu.cu->slice->getNumRefIdx(REF_PIC_LIST_1) ==0 )
    && abs(pu.cu->slice->getRefPOC( REF_PIC_LIST_0, 0) - pu.cu->slice->getPOC()) == 1)
  {
    return false;
  }

#if JVET_Z0118_GDR
  bool isClean = pu.cs->isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

  const PredictionUnit* npuGroup1[5];
  const PredictionUnit* npuGroup2[5];
  Position posGroup1[5];
  Position posGroup2[5];
  int numGroup1, numGroup2;
  numGroup1 = numGroup2 = 0;

  const int CHECKED_NEI_NUM = 5;

  const Position posLB = pu.Y().bottomLeft();
  const Position posLT = pu.Y().topLeft();
  const Position posRT = pu.Y().topRight();

  Position neiPositions[5] = { posLB.offset(-2, -1), posRT.offset(-1, -2), posRT.offset(3, -2), posLB.offset(-2, 3), posLT.offset(-2, -2) };

  int iTotalAffHMVPCandNum = 0;
  const int MAX_ALLOWED_AFF_HMVP_CAND = 1;

#if JVET_X0083_BM_AMVP_MERGE_MODE && !JVET_Y0128_NON_CTC
  RefPicList mergeRefList = REF_PIC_LIST_X;
  RefPicList amvpRefList = REF_PIC_LIST_X;
  int amvpPoc = -1;
  if (pu.amvpMergeModeFlag[0] || pu.amvpMergeModeFlag[1])
  {
    mergeRefList = pu.amvpMergeModeFlag[0] ? REF_PIC_LIST_0 : REF_PIC_LIST_1;
    amvpRefList = RefPicList(1 - mergeRefList);
    amvpPoc = pu.cu->slice->getRefPOC(amvpRefList, pu.refIdx[amvpRefList]);
  }
#endif
  for (int nei = 0; nei < CHECKED_NEI_NUM; nei++)
  {
    const PredictionUnit* puNei = pu.cs->getPURestricted(neiPositions[nei], pu, pu.chType);
    if (!puNei)
    {
      continue;
    }
    if (puNei->cu->predMode != MODE_INTER)
    {
      continue;
    }
    MotionInfo mvInfo = puNei->getMotionInfo(neiPositions[nei]);
    if (!mvInfo.isInter || mvInfo.interDir <= 0 || mvInfo.interDir > 3)
    {
      continue;
    }
    if (puNei->cu->affine && !(puNei->mergeFlag && puNei->mergeType != MRG_TYPE_DEFAULT_N))
    {
      posGroup1[numGroup1] = neiPositions[nei];
      npuGroup1[numGroup1++] = puNei;
    }
    else
    {
      posGroup2[numGroup2] = neiPositions[nei];
      npuGroup2[numGroup2++] = puNei;
    }
  }
  for (int iAffListIdx = 0; iAffListIdx < MAX_NUM_AFF_HMVP_CANDS && iTotalAffHMVPCandNum < MAX_ALLOWED_AFF_HMVP_CAND; iAffListIdx++)
  {
    for (int i = 0; i < numGroup1 && iTotalAffHMVPCandNum < MAX_ALLOWED_AFF_HMVP_CAND; i++)
    {
      const PredictionUnit* puNei = npuGroup1[i];
      MotionInfo mvInfo = puNei->getMotionInfo(posGroup1[i]);

#if JVET_Z0118_GDR
      if (addOneMergeHMVPCandFromAffModel(pu, mrgCtx, cnt, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, iAffListIdx, mvInfo, posGroup1[i], mvInfo.interDir == 3 ? puNei->cu->BcwIdx : BCW_DEFAULT
#else
      if (addOneMergeHMVPCandFromAffModel(pu, mrgCtx, cnt, pu.cs->motionLut.lutAff, iAffListIdx, mvInfo, posGroup1[i], mvInfo.interDir == 3 ? puNei->cu->BcwIdx : BCW_DEFAULT
#endif
#if INTER_LIC
        , mvInfo.interDir != 3 ? puNei->cu->LICFlag : false
#endif
#if TM_MRG
        , mvdSimilarityThresh
#endif
      ))
      {

#if JVET_X0083_BM_AMVP_MERGE_MODE
        int8_t candRefIdx[2];
        candRefIdx[0] = mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].refIdx;
        candRefIdx[1] = mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].refIdx;
#if JVET_Y0128_NON_CTC
        bool isValidAmMode = checkIsValidMergeMvCand(pu, candRefIdx);
#else
        bool isValidAmMode = checkIsValidMergeMvCand(*pu.cs, pu, pu.cu->slice->getPOC(), amvpPoc, candRefIdx);
#endif
        if (!isValidAmMode)
        {
          mrgCtx.interDirNeighbours[cnt] = 0;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].setMvField(Mv(), -1);
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(), -1);
          mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
          mrgCtx.LICFlags[cnt] = false;
#endif
          continue;
        }
#endif
        if (cnt == mrgCandIdx)
        {
#if TM_MRG
          if (!pu.tmMergeFlag)
#endif
            return true;
        }
        cnt++;
        if (cnt == mrgCtx.numValidMergeCand)
        {
          return true;
        }
        iTotalAffHMVPCandNum++;
      }
    }
    for (int i = 0; i < numGroup2 && iTotalAffHMVPCandNum < MAX_ALLOWED_AFF_HMVP_CAND; i++)
    {
      const PredictionUnit* puNei = npuGroup2[i];
      MotionInfo mvInfo = puNei->getMotionInfo(posGroup2[i]);

#if JVET_Z0118_GDR
      if (addOneMergeHMVPCandFromAffModel(pu, mrgCtx, cnt, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, iAffListIdx, mvInfo, posGroup2[i], mvInfo.interDir == 3 ? puNei->cu->BcwIdx : BCW_DEFAULT
#else
      if (addOneMergeHMVPCandFromAffModel(pu, mrgCtx, cnt, pu.cs->motionLut.lutAff, iAffListIdx, mvInfo, posGroup2[i], mvInfo.interDir == 3 ? puNei->cu->BcwIdx : BCW_DEFAULT
#endif
#if INTER_LIC
        , mvInfo.interDir != 3 ? puNei->cu->LICFlag : false
#endif
#if TM_MRG
        , mvdSimilarityThresh
#endif
      ))
      {
#if JVET_X0083_BM_AMVP_MERGE_MODE
        int8_t candRefIdx[2];
        candRefIdx[0] = mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].refIdx;
        candRefIdx[1] = mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].refIdx;
#if JVET_Y0128_NON_CTC
        bool isValidAmMode = checkIsValidMergeMvCand(pu, candRefIdx);
#else
        bool isValidAmMode = checkIsValidMergeMvCand(*pu.cs, pu, pu.cu->slice->getPOC(), amvpPoc, candRefIdx);
#endif
        if (!isValidAmMode)
        {
          mrgCtx.interDirNeighbours[cnt] = 0;
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].setMvField(Mv(), -1);
          mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(), -1);
          mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
          mrgCtx.LICFlags[cnt] = false;
#endif
          continue;
        }
#endif
        if (cnt == mrgCandIdx)
        {
#if TM_MRG
          if (!pu.tmMergeFlag)
#endif
          {
            return true;
          }
        }
        cnt++;
        if (cnt == mrgCtx.numValidMergeCand)
        {
            return true;
        }
        iTotalAffHMVPCandNum++;
      }
    }
  }
  CHECK(iTotalAffHMVPCandNum > MAX_ALLOWED_AFF_HMVP_CAND, "Invalid number of aff-HMVP based merge candidates");
  return false;
}

bool PU::addOneMergeHMVPCandFromAffModel(const PredictionUnit& pu, MergeCtx& mrgCtx, int& cnt, static_vector<AffineMotionInfo, MAX_NUM_AFF_HMVP_CANDS> *lutAff, int listIdx, const MotionInfo& mvInfo, Position neiPosition
  , int iGBiIdx 
#if INTER_LIC
  ,bool bICflag
#endif
#if TM_MRG
  , const uint32_t mvdSimilarityThresh
#endif
)
{
  Mv cMv[2];
  int HistParameters[2][4];
  mrgCtx.interDirNeighbours[cnt] = 0;
  mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
  mrgCtx.LICFlags[cnt] = false;
#endif
  mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].setMvField(Mv(), -1);
  mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(), -1);

  if ((mvInfo.interDir & 1) != 0)
  {
    CHECK(mvInfo.refIdx[0] == -1, "invalid Refidx");
    int idxInLUT = std::min((int)mvInfo.refIdx[0], MAX_NUM_AFFHMVP_ENTRIES_ONELIST - 1);
    int lutSize = (int) lutAff[idxInLUT].size();
    if (listIdx >= lutSize || lutAff[idxInLUT][lutSize - 1 - listIdx].oneSetAffineParametersPattern == 0)
    {
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].setMvField(Mv(), -1);
    }
    else
    {
      AffineMotionInfo& affHistInfo = lutAff[idxInLUT][lutSize - 1 - listIdx];
      mrgCtx.interDirNeighbours[cnt] |= 1;
      xGetAffineMvFromLUT(&affHistInfo, HistParameters[0]);
      deriveCenterMVFromAffineParameters(pu, cMv[0], HistParameters[0], mvInfo.mv[0], neiPosition);
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].setMvField(cMv[0], mvInfo.refIdx[0]);
    }
  }
  if ((mvInfo.interDir & 2) != 0)
  {
    CHECK(mvInfo.refIdx[1] == -1, "invalid Refidx");
    int idxInLUT = MAX_NUM_AFFHMVP_ENTRIES_ONELIST + std::min((int)mvInfo.refIdx[1], MAX_NUM_AFFHMVP_ENTRIES_ONELIST - 1);
    int lutSize = (int)lutAff[idxInLUT].size();
    if (listIdx >= lutSize || lutAff[idxInLUT][lutSize - 1 - listIdx].oneSetAffineParametersPattern == 0)
    {
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(), -1);
    }
    else
    {
      AffineMotionInfo& affHistInfo = lutAff[idxInLUT][lutSize - 1 - listIdx];
      mrgCtx.interDirNeighbours[cnt] |= 2;
      xGetAffineMvFromLUT(&affHistInfo, HistParameters[1]);
      deriveCenterMVFromAffineParameters(pu, cMv[1], HistParameters[1], mvInfo.mv[1], neiPosition);
      mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(cMv[1], mvInfo.refIdx[1]);

    }
  }
  if (mrgCtx.interDirNeighbours[cnt] == 0)
  {
    mrgCtx.BcwIdx  [cnt] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[cnt] = false;
#endif
    mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].setMvField(Mv(), -1);
    mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(), -1);
    return false;
  }
  mrgCtx.BcwIdx[cnt] = mrgCtx.interDirNeighbours[cnt] == 3 ? iGBiIdx : BCW_DEFAULT;

#if INTER_LIC
  mrgCtx.LICFlags[cnt] = bICflag;
  CHECK(bICflag && mvInfo.interDir == 3, "LIC cannot be used for Bi");
#endif

#if NON_ADJACENT_MRG_CAND
  if (mrgCtx.xCheckSimilarMotion(cnt
#if TM_MRG
    , mvdSimilarityThresh
#endif
  ))
  {
    mrgCtx.interDirNeighbours[cnt] = 0;
    mrgCtx.BcwIdx[cnt] = BCW_DEFAULT;
#if INTER_LIC
    mrgCtx.LICFlags[cnt] = false;
#endif
    mrgCtx.mvFieldNeighbours[(cnt << 1) + 0].setMvField(Mv(), -1);
    mrgCtx.mvFieldNeighbours[(cnt << 1) + 1].setMvField(Mv(), -1);
    return false;
  }
#endif
  return true;
}

bool PU::addOneAffineMergeHMVPCand(const PredictionUnit & pu, AffineMergeCtx & affMrgCtx, static_vector<AffineMotionInfo, MAX_NUM_AFF_HMVP_CANDS>* lutAff, int listIdx, const MotionInfo & mvInfo, Position neiPosition, int iGBiIdx
#if INTER_LIC
                                 ,       bool bICflag
#endif
)
{
  Mv cMv[2][3];
  int aiHistParameters[2][4];
  affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] = 0;
  if ((mvInfo.interDir & 1) != 0)
  {
    CHECK(mvInfo.refIdx[0] == -1, "invalid Refidx");
    int idxInLUT = std::min((int)mvInfo.refIdx[0], MAX_NUM_AFFHMVP_ENTRIES_ONELIST - 1);
    int lutSize = (int) lutAff[idxInLUT].size();
    if (listIdx >= lutSize || lutAff[idxInLUT][lutSize - 1 - listIdx].oneSetAffineParametersPattern == 0)
    {
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][0].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][1].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][2].setMvField(Mv(), -1);
    }
    else
    {
      AffineMotionInfo & affHistInfo = lutAff[idxInLUT][lutSize - 1 - listIdx];
      affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] |= 1;
      xGetAffineMvFromLUT(&affHistInfo, aiHistParameters[0]);
      deriveMVsFromAffineParameters(pu, cMv[0], aiHistParameters[0], mvInfo.mv[0], neiPosition);
      for (int mvNum = 0; mvNum < 3; mvNum++)
      {
        affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][mvNum].setMvField(cMv[0][mvNum], mvInfo.refIdx[0]);
      }
    }
  }
  if ((mvInfo.interDir & 2) != 0)
  {
    CHECK(mvInfo.refIdx[1] == -1, "invalid Refidx");
    int idxInLUT = MAX_NUM_AFFHMVP_ENTRIES_ONELIST + std::min((int)mvInfo.refIdx[1], MAX_NUM_AFFHMVP_ENTRIES_ONELIST - 1);
    int lutSize = (int) lutAff[idxInLUT].size();
    if (listIdx >= lutSize || lutAff[idxInLUT][lutSize - 1 - listIdx].oneSetAffineParametersPattern == 0)
    {
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][0].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][1].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][2].setMvField(Mv(), -1);
    }
    else
    {
      AffineMotionInfo& affHistInfo = lutAff[idxInLUT][lutSize - 1 - listIdx];
      affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] |= 2;
      xGetAffineMvFromLUT(&affHistInfo, aiHistParameters[1]);
      deriveMVsFromAffineParameters(pu, cMv[1], aiHistParameters[1], mvInfo.mv[1], neiPosition);
      for (int mvNum = 0; mvNum < 3; mvNum++)
      {
        affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][mvNum].setMvField(cMv[1][mvNum], mvInfo.refIdx[1]);
      }
    }
  }
  if (affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] == 0)
  {
    return false;
  }
  affMrgCtx.affineType[affMrgCtx.numValidMergeCand] = AFFINEMODEL_6PARAM;

  affMrgCtx.BcwIdx[affMrgCtx.numValidMergeCand] = iGBiIdx;

#if INTER_LIC
  affMrgCtx.LICFlags[affMrgCtx.numValidMergeCand] = bICflag;
  CHECK( bICflag && mvInfo.interDir == 3, "LIC cannot be used for Bi");
#endif


  if (checkLastAffineMergeCandRedundancy(pu, affMrgCtx))
  {
    return true;
  }
  return false;
}



bool PU::addOneInheritedHMVPAffineMergeCand(const PredictionUnit& pu, AffineMergeCtx& affMrgCtx, static_vector<AffineInheritInfo, MAX_NUM_AFF_INHERIT_HMVP_CANDS>& lutAffInherit, int listIdx)
{
  Mv cMv[2][3];
  int aiHistParameters[2][4];

  int lutSize = (int)lutAffInherit.size();

  if (listIdx >= lutSize)
  {
    return false;
  }
  AffineInheritInfo& affHistInfo = lutAffInherit[lutSize - 1 - listIdx];

  affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] = 0;

  affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][0].setMvField(Mv(), -1);
  affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][1].setMvField(Mv(), -1);
  affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][2].setMvField(Mv(), -1);
  affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][0].setMvField(Mv(), -1);
  affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][1].setMvField(Mv(), -1);
  affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][2].setMvField(Mv(), -1);

  if (affHistInfo.oneSetAffineParametersPattern0 != 0 && affHistInfo.baseMV[0].refIdx != -1)
  {
    affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] |= 1;
    xGetAffineMvFromLUT(affHistInfo.oneSetAffineParameters0, aiHistParameters[0]);
    deriveMVsFromAffineParameters(pu, cMv[0], aiHistParameters[0], affHistInfo.baseMV[0].mv, affHistInfo.basePos);
    for (int mvNum = 0; mvNum < 3; mvNum++)
    {
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][mvNum].setMvField(cMv[0][mvNum], affHistInfo.baseMV[0].refIdx);
    }
  }

  if (affHistInfo.oneSetAffineParametersPattern1 != 0 && affHistInfo.baseMV[1].refIdx != -1)
  {
    affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] |= 2;
    xGetAffineMvFromLUT(affHistInfo.oneSetAffineParameters1, aiHistParameters[1]);
    deriveMVsFromAffineParameters(pu, cMv[1], aiHistParameters[1], affHistInfo.baseMV[1].mv, affHistInfo.basePos);
    for (int mvNum = 0; mvNum < 3; mvNum++)
    {
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][mvNum].setMvField(cMv[1][mvNum], affHistInfo.baseMV[1].refIdx);
    }
  }

  if (affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] == 0)
  {
    return false;
  }
  affMrgCtx.affineType[affMrgCtx.numValidMergeCand] = AFFINEMODEL_6PARAM;

  affMrgCtx.BcwIdx[affMrgCtx.numValidMergeCand] = BCW_DEFAULT;

#if INTER_LIC
  affMrgCtx.LICFlags[affMrgCtx.numValidMergeCand] = false;
#endif

  if (checkLastAffineMergeCandRedundancy(pu, affMrgCtx))
  {
    return true;
  }
  return false;
}

bool PU::addSpatialAffineMergeHMVPCand(const PredictionUnit& pu, AffineMergeCtx& affMrgCtx, static_vector<AffineMotionInfo, MAX_NUM_AFF_HMVP_CANDS>* lutAff, int affHMVPIdx, const PredictionUnit* neiPUs[], Position neiPositions[], int iNeiNum, const int mrgCandIdx)
{
  const Slice& slice = *pu.cs->slice;
  const uint32_t maxNumAffineMergeCand = slice.getPicHeader()->getMaxNumAffineMergeCand();

  for (int nei = 0; nei < iNeiNum; nei++)
  {
    const PredictionUnit* puNei = neiPUs[nei];
    CHECK(!neiPUs[nei], "Invalid neighbour PU");
    MotionInfo mvInfo = puNei->getMotionInfo(neiPositions[nei]);
    if (mvInfo.isIBCmot)
    {
      continue;
    }
    if (addOneAffineMergeHMVPCand(pu, affMrgCtx, lutAff, affHMVPIdx, mvInfo, neiPositions[nei], puNei->interDir == 3? puNei->cu->BcwIdx : BCW_DEFAULT
#if INTER_LIC
                                , puNei->cu->LICFlag
#endif
    ))
    {
      if (affMrgCtx.numValidMergeCand == mrgCandIdx) // for decoder 
      {
        return true;
      }

      affMrgCtx.numValidMergeCand++;

      // early termination
      if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
      {
        return true;
      }
    }
  } //for nei
  return false;
}

bool PU::addSpatialAffineAMVPHMVPCand(PredictionUnit & pu, const RefPicList & eRefPicList, const int& refIdx, AffineAMVPInfo & affiAMVPInfo, static_vector<AffineMotionInfo, MAX_NUM_AFF_HMVP_CANDS>* lutAff, int iHMVPlistIdx,
  int neiIdx[], int iNeiNum, int aiNeibeInherited[], bool bFoundOne)
{
  if (affiAMVPInfo.numCand >= AMVP_MAX_NUM_CANDS)
  {
    return false;
  }

  const Slice& slice = *pu.cs->slice;
  const Position posLB = pu.Y().bottomLeft();
  const Position posLT = pu.Y().topLeft();
  const Position posRT = pu.Y().topRight();

  Position neiPositions[5] = { posLB.offset(-2, -1), posRT.offset(-1, -2), posRT.offset(3, -2), posLB.offset(-2, 3), posLT.offset(-2, -2) };


  const int        currRefPOC = pu.cs->slice->getRefPic(eRefPicList, refIdx)->getPOC();

  AffineMotionInfo affHistInfo;
  int idxInLUT = (int)eRefPicList * MAX_NUM_AFFHMVP_ENTRIES_ONELIST + std::min(refIdx, MAX_NUM_AFFHMVP_ENTRIES_ONELIST - 1);
  if (iHMVPlistIdx >= lutAff[idxInLUT].size())
  {
    return false;
  }
  affHistInfo = lutAff[idxInLUT][lutAff[idxInLUT].size() - 1 - iHMVPlistIdx];
  if (affHistInfo.oneSetAffineParametersPattern == 0)
  {
    return false;
  }

  Mv cMv[3];
  int aiHistParameters[4];

  for (int idx = 0; idx < iNeiNum; idx++)
  {
    int nei = neiIdx[idx];
    if (aiNeibeInherited[nei])
    {
      continue;
    }
    aiNeibeInherited[nei] = 1;
    const PredictionUnit* puNei = pu.cs->getPURestricted(neiPositions[nei], pu, pu.chType);
    if (!puNei || !CU::isInter(*pu.cu))
    {
      continue;
    }

    MotionInfo mvInfo = puNei->getMotionInfo(neiPositions[nei]);
    if (mvInfo.isIBCmot)
    {
      continue;
    }
    if (!mvInfo.isInter || mvInfo.interDir <= 0 || mvInfo.interDir > 3)
    {
      continue;
    }
    RefPicList selRefPicListIndex = eRefPicList;
    if (((mvInfo.interDir & (selRefPicListIndex + 1)) == 0) || ((slice.getRefPic(selRefPicListIndex, mvInfo.refIdx[selRefPicListIndex])->getPOC()) != currRefPOC))
    {
      selRefPicListIndex = RefPicList(1 - eRefPicList);
      if (((mvInfo.interDir & (selRefPicListIndex + 1)) == 0) || ((slice.getRefPic(selRefPicListIndex, mvInfo.refIdx[selRefPicListIndex])->getPOC()) != currRefPOC))
      {
        continue;
      }
    }

    xGetAffineMvFromLUT(&affHistInfo, aiHistParameters);
    deriveMVsFromAffineParameters(pu, cMv, aiHistParameters, mvInfo.mv[selRefPicListIndex], neiPositions[nei]);

    if (pu.cu->imv == 0)
    {
      cMv[0].roundToPrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
      cMv[1].roundToPrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
      cMv[2].roundToPrecision(MV_PRECISION_INTERNAL, MV_PRECISION_QUARTER);
    }
    else if (pu.cu->imv == 2)
    {
      cMv[0].roundToPrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
      cMv[1].roundToPrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
      cMv[2].roundToPrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
    }

    affiAMVPInfo.mvCandLT[affiAMVPInfo.numCand] = cMv[0];
    affiAMVPInfo.mvCandRT[affiAMVPInfo.numCand] = cMv[1];
    affiAMVPInfo.mvCandLB[affiAMVPInfo.numCand] = cMv[2];


    if (!checkLastAffineAMVPCandRedundancy(pu, affiAMVPInfo))
    {
      continue;
    }

    affiAMVPInfo.numCand++;

    if (affiAMVPInfo.numCand == AMVP_MAX_NUM_CANDS)
    {
      return true;
    }
    if (bFoundOne)
    {
      return false;
    }
  }

  return false;
}

#endif
#if JVET_Z0139_HIST_AFF || JVET_Z0139_NA_AFF
bool PU::checkLastAffineAMVPCandRedundancy(const PredictionUnit& pu, AffineAMVPInfo& affiAMVPInfo)
{
  if (affiAMVPInfo.numCand < 1)
  {
    return true;
  }

  const int CPMV_SIMILARITY_THREH = 0;
  const int lastIdx = affiAMVPInfo.numCand;

  for (int idx = 0; idx < affiAMVPInfo.numCand; idx++)
  {

    if (
      abs(affiAMVPInfo.mvCandLT[idx].getHor() - affiAMVPInfo.mvCandLT[lastIdx].getHor()) > CPMV_SIMILARITY_THREH ||
      abs(affiAMVPInfo.mvCandLT[idx].getVer() - affiAMVPInfo.mvCandLT[lastIdx].getVer()) > CPMV_SIMILARITY_THREH ||
      abs(affiAMVPInfo.mvCandRT[idx].getHor() - affiAMVPInfo.mvCandRT[lastIdx].getHor()) > CPMV_SIMILARITY_THREH ||
      abs(affiAMVPInfo.mvCandRT[idx].getVer() - affiAMVPInfo.mvCandRT[lastIdx].getVer()) > CPMV_SIMILARITY_THREH)
    {
      continue;
    }
    if (pu.cu->affineType == AFFINEMODEL_6PARAM)
    {
      if (
        abs(affiAMVPInfo.mvCandLB[idx].getHor() - affiAMVPInfo.mvCandLB[lastIdx].getHor()) > CPMV_SIMILARITY_THREH ||
        abs(affiAMVPInfo.mvCandLB[idx].getVer() - affiAMVPInfo.mvCandLB[lastIdx].getVer()) > CPMV_SIMILARITY_THREH)
      {
        continue;
      }
    }
    return false;
  }
  return true;
}
#endif
void PU::fillAffineMvpCand(PredictionUnit &pu, const RefPicList &eRefPicList, const int &refIdx, AffineAMVPInfo &affiAMVPInfo)
{
  affiAMVPInfo.numCand = 0;

  if (refIdx < 0)
  {
    return;
  }

#if JVET_Z0118_GDR
  bool isClean = pu.cs->isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

  // insert inherited affine candidates
  Mv outputAffineMv[3];
  Position posLT = pu.Y().topLeft();
  Position posRT = pu.Y().topRight();
  Position posLB = pu.Y().bottomLeft();

#if JVET_Z0139_HIST_AFF
  int aiNeibeInherited[5];
  memset(aiNeibeInherited, 0, sizeof(aiNeibeInherited));
  Position neiPositions[5] = { posLB.offset(-2, -1), posRT.offset(-1, -2), posRT.offset(3, -2), posLB.offset(-2, 3), posLT.offset(-2, -2) };
  const PredictionUnit * neiPU[5];
  for (int nei = 0; nei < 5; nei++)
  {
    neiPU[nei] = pu.cs->getPURestricted(neiPositions[nei], pu, pu.chType);
  }
  int neiIdx[5] = { 0, 1, 2, 3, 4, };
  int leftNeiIdx[2];
  int leftAffNeiNum = 0;

  int aboveNeiIdx[3];
  int aboveAffNeiNum = 0;

  if (neiPU[3] && CU::isInter(*neiPU[3]->cu) && neiPU[3]->cu->affine && neiPU[3]->mergeType == MRG_TYPE_DEFAULT_N)
  {
    MotionInfo mvInfo = neiPU[3]->getMotionInfo(neiPositions[3]);
    if (mvInfo.isInter && mvInfo.interDir > 1 && mvInfo.interDir <= 3)
    {
      leftNeiIdx[leftAffNeiNum++] = 3;
    }
  }
  if (neiPU[0] && CU::isInter(*neiPU[0]->cu) && neiPU[0]->cu->affine && neiPU[0]->mergeType == MRG_TYPE_DEFAULT_N)
  {
    MotionInfo mvInfo = neiPU[0]->getMotionInfo(neiPositions[0]);
    if (mvInfo.isInter && mvInfo.interDir > 1 && mvInfo.interDir <= 3)
    {
      leftNeiIdx[leftAffNeiNum++] = 0;
    }
  }
  if (neiPU[2] && CU::isInter(*neiPU[2]->cu) && neiPU[2]->cu->affine && neiPU[2]->mergeType == MRG_TYPE_DEFAULT_N)
  {
    MotionInfo mvInfo = neiPU[2]->getMotionInfo(neiPositions[2]);
    if (mvInfo.isInter && mvInfo.interDir > 1 && mvInfo.interDir <= 3)
    {
      aboveNeiIdx[aboveAffNeiNum++] = 2;
    }
  }
  if (neiPU[1] && CU::isInter(*neiPU[1]->cu) && neiPU[1]->cu->affine && neiPU[1]->mergeType == MRG_TYPE_DEFAULT_N)
  {
    MotionInfo mvInfo = neiPU[1]->getMotionInfo(neiPositions[1]);
    if (mvInfo.isInter && mvInfo.interDir > 1 && mvInfo.interDir <= 3)
    {
      aboveNeiIdx[aboveAffNeiNum++] = 1;
    }
  }
  if (neiPU[4] && CU::isInter(*neiPU[4]->cu) && neiPU[4]->cu->affine && neiPU[4]->mergeType == MRG_TYPE_DEFAULT_N)
  {
    MotionInfo mvInfo = neiPU[4]->getMotionInfo(neiPositions[4]);
    if (mvInfo.isInter && mvInfo.interDir > 1 && mvInfo.interDir <= 3)
    {
      aboveNeiIdx[aboveAffNeiNum++] = 4;
    }
  }
  // check left neighbor 
  if (!addAffineMVPCandUnscaled(pu, eRefPicList, refIdx, posLB, MD_BELOW_LEFT, affiAMVPInfo, aiNeibeInherited))
  {
    addAffineMVPCandUnscaled(pu, eRefPicList, refIdx, posLB, MD_LEFT, affiAMVPInfo, aiNeibeInherited);
  }
  leftAffNeiNum = 0;

  // check above neighbor
  if (!addAffineMVPCandUnscaled(pu, eRefPicList, refIdx, posRT, MD_ABOVE_RIGHT, affiAMVPInfo, aiNeibeInherited))
  {
    if (!addAffineMVPCandUnscaled(pu, eRefPicList, refIdx, posRT, MD_ABOVE, affiAMVPInfo, aiNeibeInherited))
    {
      addAffineMVPCandUnscaled(pu, eRefPicList, refIdx, posLT, MD_ABOVE_LEFT, affiAMVPInfo, aiNeibeInherited);
    }
  }
  aboveAffNeiNum = 0;

  for (int affHMVPIdx = 0; affHMVPIdx < 1; affHMVPIdx++)
  {
    if (affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS && leftAffNeiNum > 0)
    {
#if JVET_Z0118_GDR
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, 0, leftNeiIdx, leftAffNeiNum, aiNeibeInherited, true);
#else
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, pu.cs->motionLut.lutAff, 0, leftNeiIdx, leftAffNeiNum, aiNeibeInherited, true);
#endif
    }
    if (affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS && aboveAffNeiNum > 0)
    {
#if JVET_Z0118_GDR
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, 0, aboveNeiIdx, aboveAffNeiNum, aiNeibeInherited, true);
#else
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, pu.cs->motionLut.lutAff, 0, aboveNeiIdx, aboveAffNeiNum, aiNeibeInherited, true);
#endif
    }
  }
#else
  // check left neighbor
  if ( !addAffineMVPCandUnscaled( pu, eRefPicList, refIdx, posLB, MD_BELOW_LEFT, affiAMVPInfo ) )
  {
    addAffineMVPCandUnscaled( pu, eRefPicList, refIdx, posLB, MD_LEFT, affiAMVPInfo );
  }

  // check above neighbor
  if ( !addAffineMVPCandUnscaled( pu, eRefPicList, refIdx, posRT, MD_ABOVE_RIGHT, affiAMVPInfo ) )
  {
    if ( !addAffineMVPCandUnscaled( pu, eRefPicList, refIdx, posRT, MD_ABOVE, affiAMVPInfo ) )
    {
      addAffineMVPCandUnscaled( pu, eRefPicList, refIdx, posLT, MD_ABOVE_LEFT, affiAMVPInfo );
    }
  }
#endif

  if ( affiAMVPInfo.numCand >= AMVP_MAX_NUM_CANDS )
  {
    for (int i = 0; i < affiAMVPInfo.numCand; i++)
    {
      affiAMVPInfo.mvCandLT[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
      affiAMVPInfo.mvCandRT[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
      affiAMVPInfo.mvCandLB[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
    }
    return;
  }

  // insert constructed affine candidates
  int cornerMVPattern = 0;

  //-------------------  V0 (START) -------------------//
  AMVPInfo amvpInfo0;
  amvpInfo0.numCand = 0;

  // A->C: Above Left, Above, Left
  addMVPCandUnscaled( pu, eRefPicList, refIdx, posLT, MD_ABOVE_LEFT, amvpInfo0 );
  if ( amvpInfo0.numCand < 1 )
  {
    addMVPCandUnscaled( pu, eRefPicList, refIdx, posLT, MD_ABOVE, amvpInfo0 );
  }
  if ( amvpInfo0.numCand < 1 )
  {
    addMVPCandUnscaled( pu, eRefPicList, refIdx, posLT, MD_LEFT, amvpInfo0 );
  }
  cornerMVPattern = cornerMVPattern | amvpInfo0.numCand;

  //-------------------  V1 (START) -------------------//
  AMVPInfo amvpInfo1;
  amvpInfo1.numCand = 0;

  // D->E: Above, Above Right
  addMVPCandUnscaled( pu, eRefPicList, refIdx, posRT, MD_ABOVE, amvpInfo1 );
  if ( amvpInfo1.numCand < 1 )
  {
    addMVPCandUnscaled( pu, eRefPicList, refIdx, posRT, MD_ABOVE_RIGHT, amvpInfo1 );
  }
  cornerMVPattern = cornerMVPattern | (amvpInfo1.numCand << 1);

  //-------------------  V2 (START) -------------------//
  AMVPInfo amvpInfo2;
  amvpInfo2.numCand = 0;

  // F->G: Left, Below Left
  addMVPCandUnscaled( pu, eRefPicList, refIdx, posLB, MD_LEFT, amvpInfo2 );
  if ( amvpInfo2.numCand < 1 )
  {
    addMVPCandUnscaled( pu, eRefPicList, refIdx, posLB, MD_BELOW_LEFT, amvpInfo2 );
  }
  cornerMVPattern = cornerMVPattern | (amvpInfo2.numCand << 2);

  outputAffineMv[0] = amvpInfo0.mvCand[0];
  outputAffineMv[1] = amvpInfo1.mvCand[0];
  outputAffineMv[2] = amvpInfo2.mvCand[0];

  outputAffineMv[0].roundAffinePrecInternal2Amvr(pu.cu->imv);
  outputAffineMv[1].roundAffinePrecInternal2Amvr(pu.cu->imv);
  outputAffineMv[2].roundAffinePrecInternal2Amvr(pu.cu->imv);

  if ( cornerMVPattern == 7 || (cornerMVPattern == 3 && pu.cu->affineType == AFFINEMODEL_4PARAM) )
  {
    affiAMVPInfo.mvCandLT[affiAMVPInfo.numCand] = outputAffineMv[0];
    affiAMVPInfo.mvCandRT[affiAMVPInfo.numCand] = outputAffineMv[1];
    affiAMVPInfo.mvCandLB[affiAMVPInfo.numCand] = outputAffineMv[2];

#if JVET_Z0139_HIST_AFF
    if (checkLastAffineAMVPCandRedundancy(pu, affiAMVPInfo))
#endif
    affiAMVPInfo.numCand++;
  }

#if JVET_Z0139_HIST_AFF
  for (int affHMVPIdx = 1; affHMVPIdx < MAX_NUM_AFF_HMVP_CANDS; affHMVPIdx++)
  {
    if (affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS && leftAffNeiNum > 0)
    {
#if JVET_Z0118_GDR
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, 0, leftNeiIdx, leftAffNeiNum, aiNeibeInherited, true);
#else
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, pu.cs->motionLut.lutAff, 0, leftNeiIdx, leftAffNeiNum, aiNeibeInherited, true);
#endif
    }
    if (affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS && aboveAffNeiNum > 0)
    {
#if JVET_Z0118_GDR
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, 0, aboveNeiIdx, aboveAffNeiNum, aiNeibeInherited, true);
#else
      addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, pu.cs->motionLut.lutAff, 0, aboveNeiIdx, aboveAffNeiNum, aiNeibeInherited, true);
#endif
    }
  }
#endif

  if ( affiAMVPInfo.numCand < 2 )
  {
    // check corner MVs
    for ( int i = 2; i >= 0 && affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS; i-- )
    {
      if ( cornerMVPattern & (1 << i) ) // MV i exist
      {
        affiAMVPInfo.mvCandLT[affiAMVPInfo.numCand] = outputAffineMv[i];
        affiAMVPInfo.mvCandRT[affiAMVPInfo.numCand] = outputAffineMv[i];
        affiAMVPInfo.mvCandLB[affiAMVPInfo.numCand] = outputAffineMv[i];
#if JVET_Z0139_HIST_AFF
        if (checkLastAffineAMVPCandRedundancy(pu, affiAMVPInfo))
#endif
          affiAMVPInfo.numCand++;
      }
    }

    // Get Temporal Motion Predictor
    if ( affiAMVPInfo.numCand < 2 && pu.cs->picHeader->getEnableTMVPFlag() )
    {
      const int refIdxCol = refIdx;

      Position posRB = pu.Y().bottomRight().offset( -3, -3 );

      const PreCalcValues& pcv = *pu.cs->pcv;

      Position posC0;
      bool bC0Avail = false;
      Position posC1 = pu.Y().center();
      Mv cColMv;
      bool boundaryCond = ((posRB.x + pcv.minCUWidth) < pcv.lumaWidth) && ((posRB.y + pcv.minCUHeight) < pcv.lumaHeight);
      const SubPic &curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
      if (curSubPic.getTreatedAsPicFlag())
      {
        boundaryCond = ((posRB.x + pcv.minCUWidth) <= curSubPic.getSubPicRight() &&
          (posRB.y + pcv.minCUHeight) <= curSubPic.getSubPicBottom());
      }
      if (boundaryCond)
      {
        int posYInCtu = posRB.y & pcv.maxCUHeightMask;
        if (posYInCtu + 4 < pcv.maxCUHeight)
        {
          posC0 = posRB.offset(4, 4);
          bC0Avail = true;
        }
      }
      if ( (bC0Avail && getColocatedMVP( pu, eRefPicList, posC0, cColMv, refIdxCol, false ) ) || getColocatedMVP( pu, eRefPicList, posC1, cColMv, refIdxCol, false ) )
      {
        cColMv.roundAffinePrecInternal2Amvr(pu.cu->imv);
        affiAMVPInfo.mvCandLT[affiAMVPInfo.numCand] = cColMv;
        affiAMVPInfo.mvCandRT[affiAMVPInfo.numCand] = cColMv;
        affiAMVPInfo.mvCandLB[affiAMVPInfo.numCand] = cColMv;
#if JVET_Z0139_HIST_AFF
        if (checkLastAffineAMVPCandRedundancy(pu, affiAMVPInfo))
#endif
          affiAMVPInfo.numCand++;
      }
    }

#if JVET_Z0139_NA_AFF
    if (affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS)
    {
      addNonAdjCstAffineMVPCandUnscaled(pu, eRefPicList, refIdx, affiAMVPInfo);
    }
#endif

#if JVET_Z0139_HIST_AFF
    if (affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS)
    {
      for (int affHMVPIdx = 0; affHMVPIdx < MAX_NUM_AFF_HMVP_CANDS; affHMVPIdx++)
      {
#if JVET_Z0118_GDR
        addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, affHMVPIdx, neiIdx, 5, aiNeibeInherited, false);
#else
        addSpatialAffineAMVPHMVPCand(pu, eRefPicList, refIdx, affiAMVPInfo, pu.cs->motionLut.lutAff, affHMVPIdx, neiIdx, 5, aiNeibeInherited, false);
#endif
      }
    }
#endif
#if JVET_Z0139_HIST_AFF
    if (affiAMVPInfo.numCand < AMVP_MAX_NUM_CANDS)
#else
    if (affiAMVPInfo.numCand < 2)
#endif
    {
      // add zero MV
      for ( int i = affiAMVPInfo.numCand; i < AMVP_MAX_NUM_CANDS; i++ )
      {
        affiAMVPInfo.mvCandLT[affiAMVPInfo.numCand].setZero();
        affiAMVPInfo.mvCandRT[affiAMVPInfo.numCand].setZero();
        affiAMVPInfo.mvCandLB[affiAMVPInfo.numCand].setZero();
        affiAMVPInfo.numCand++;
      }
    }
  }

  for (int i = 0; i < affiAMVPInfo.numCand; i++)
  {
    affiAMVPInfo.mvCandLT[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
    affiAMVPInfo.mvCandRT[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
    affiAMVPInfo.mvCandLB[i].roundAffinePrecInternal2Amvr(pu.cu->imv);
  }
}

bool PU::addMVPCandUnscaled( const PredictionUnit &pu, const RefPicList &eRefPicList, const int &iRefIdx, const Position &pos, const MvpDir &eDir, AMVPInfo &info )
{
        CodingStructure &cs    = *pu.cs;
  const PredictionUnit *neibPU = NULL;
        Position neibPos;

  switch (eDir)
  {
  case MD_LEFT:
    neibPos = pos.offset( -1,  0 );
    break;
  case MD_ABOVE:
    neibPos = pos.offset(  0, -1 );
    break;
  case MD_ABOVE_RIGHT:
    neibPos = pos.offset(  1, -1 );
    break;
  case MD_BELOW_LEFT:
    neibPos = pos.offset( -1,  1 );
    break;
  case MD_ABOVE_LEFT:
    neibPos = pos.offset( -1, -1 );
    break;
  default:
    break;
  }

  neibPU = cs.getPURestricted( neibPos, pu, pu.chType );

#if JVET_Y0065_GPM_INTRA
  if( neibPU == NULL || !CU::isInter( *neibPU->cu ) || !neibPU->getMotionInfo( neibPos ).isInter)
#else
  if( neibPU == NULL || !CU::isInter( *neibPU->cu ) )
#endif
  {
    return false;
  }

  const MotionInfo& neibMi        = neibPU->getMotionInfo( neibPos );

  const int        currRefPOC     = cs.slice->getRefPic( eRefPicList, iRefIdx )->getPOC();
  const RefPicList eRefPicList2nd = ( eRefPicList == REF_PIC_LIST_0 ) ? REF_PIC_LIST_1 : REF_PIC_LIST_0;

  for( int predictorSource = 0; predictorSource < 2; predictorSource++ ) // examine the indicated reference picture list, then if not available, examine the other list.
  {
    const RefPicList eRefPicListIndex = ( predictorSource == 0 ) ? eRefPicList : eRefPicList2nd;
    const int        neibRefIdx       = neibMi.refIdx[eRefPicListIndex];

    if( neibRefIdx >= 0 && currRefPOC == cs.slice->getRefPOC( eRefPicListIndex, neibRefIdx ) )
    {
#if TM_AMVP
      info.mvCand[info.numCand] = neibMi.mv[eRefPicListIndex];
      if (!info.xCheckSimilarMotion(info.numCand))
      {
        info.numCand++;
        return true;
      }
#else
      info.mvCand[info.numCand++] = neibMi.mv[eRefPicListIndex];
      return true;
#endif
    }
  }

  return false;
}


void PU::addAMVPHMVPCand(const PredictionUnit &pu, const RefPicList eRefPicList, const int currRefPOC, AMVPInfo &info)
{
  const Slice &slice = *(*pu.cs).slice;

  MotionInfo neibMi;
#if JVET_Z0118_GDR  
  bool isClean = pu.cs->isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

#if JVET_Z0075_IBC_HMVP_ENLARGE
#if JVET_Z0118_GDR  
  auto &lut = (isClean) ? pu.cs->motionLut.lut1 : pu.cs->motionLut.lut0;
#else
  auto &lut = pu.cs->motionLut.lut;
#endif
#else
#if JVET_Z0118_GDR  
  auto &lut = CU::isIBC(*pu.cu) ? (isClean ? pu.cs->motionLut.lutIbc1 : pu.cs->motionLut.lutIbc0) : (isClean ? pu.cs->motionLut.lut1 : pu.cs->motionLut.lut0);
#else
  auto &lut = CU::isIBC(*pu.cu) ? pu.cs->motionLut.lutIbc : pu.cs->motionLut.lut;
#endif
#endif

  int num_avai_candInLUT = (int) lut.size();
  int num_allowedCand = std::min(MAX_NUM_HMVP_AVMPCANDS, num_avai_candInLUT);
  const RefPicList eRefPicList2nd = (eRefPicList == REF_PIC_LIST_0) ? REF_PIC_LIST_1 : REF_PIC_LIST_0;

  for (int mrgIdx = 1; mrgIdx <= num_allowedCand; mrgIdx++)
  {
#if TM_AMVP
    if (info.numCand >= info.maxStorageSize)
#else
    if (info.numCand >= AMVP_MAX_NUM_CANDS)
#endif
    {
      return;
    }
    neibMi = lut[mrgIdx - 1];

    for (int predictorSource = 0; predictorSource < 2; predictorSource++)
    {
      const RefPicList eRefPicListIndex = (predictorSource == 0) ? eRefPicList : eRefPicList2nd;
      const int        neibRefIdx = neibMi.refIdx[eRefPicListIndex];

#if JVET_Z0075_IBC_HMVP_ENLARGE
      if (neibRefIdx >= 0 && ((currRefPOC == slice.getRefPOC(eRefPicListIndex, neibRefIdx))))
#else
      if (neibRefIdx >= 0 && (CU::isIBC(*pu.cu) || (currRefPOC == slice.getRefPOC(eRefPicListIndex, neibRefIdx))))
#endif
      {
        Mv pmv = neibMi.mv[eRefPicListIndex];
        pmv.roundTransPrecInternal2Amvr(pu.cu->imv);

#if TM_AMVP
        info.mvCand[info.numCand] = pmv;
        if (!info.xCheckSimilarMotion(info.numCand))
        {
          info.numCand++;
        }

        if( info.numCand >= info.maxStorageSize )
        {
          return;
        }
#else
        info.mvCand[info.numCand++] = pmv;
        if (info.numCand >= AMVP_MAX_NUM_CANDS)
        {
          return;
        }
#endif
      }
    }
  }
}

bool PU::isBipredRestriction(const PredictionUnit &pu)
{
#if !INTER_RM_SIZE_CONSTRAINTS
  if(pu.cu->lumaSize().width == 4 && pu.cu->lumaSize().height ==4 )
  {
    return true;
  }
  /* disable bi-prediction for 4x8/8x4 */
  if ( pu.cu->lumaSize().width + pu.cu->lumaSize().height == 12 )
  {
    return true;
  }
#endif
  return false;
}

#if JVET_Z0139_HIST_AFF
bool PU::getAffineControlPointCand(const PredictionUnit &pu, MotionInfo mi[4], bool isAvailable[4], int verIdx[4], int8_t bcwIdx, int modelIdx, int verNum, AffineMergeCtx& affMrgType)
#else
void PU::getAffineControlPointCand(const PredictionUnit &pu, MotionInfo mi[4], bool isAvailable[4], int verIdx[4], int8_t bcwIdx, int modelIdx, int verNum, AffineMergeCtx& affMrgType)
#endif
{
  int cuW = pu.Y().width;
  int cuH = pu.Y().height;
  int vx, vy;
  int shift = MAX_CU_DEPTH;
  int shiftHtoW = shift + floorLog2(cuW) - floorLog2(cuH);

  // motion info
  Mv cMv[2][4];
  int refIdx[2] = { -1, -1 };
  int dir = 0;
  EAffineModel curType = (verNum == 2) ? AFFINEMODEL_4PARAM : AFFINEMODEL_6PARAM;
#if INTER_LIC
  bool LICFlag = false;
#endif
  if ( verNum == 2 )
  {
    int idx0 = verIdx[0], idx1 = verIdx[1];
    if ( !isAvailable[idx0] || !isAvailable[idx1] )
    {
#if JVET_Z0139_HIST_AFF
      return false;
#else
      return;
#endif
    }

    for ( int l = 0; l < 2; l++ )
    {
      if ( mi[idx0].refIdx[l] >= 0 && mi[idx1].refIdx[l] >= 0 )
      {
        // check same refidx and different mv
        if ( mi[idx0].refIdx[l] == mi[idx1].refIdx[l])
        {
          dir |= (l + 1);
          refIdx[l] = mi[idx0].refIdx[l];
        }
      }
    }
#if INTER_LIC
    LICFlag = mi[idx0].usesLIC || mi[idx1].usesLIC;
#endif
  }
  else if ( verNum == 3 )
  {
    int idx0 = verIdx[0], idx1 = verIdx[1], idx2 = verIdx[2];
    if ( !isAvailable[idx0] || !isAvailable[idx1] || !isAvailable[idx2] )
    {
#if JVET_Z0139_HIST_AFF
      return false;
#else
      return;
#endif
    }

    for ( int l = 0; l < 2; l++ )
    {
      if ( mi[idx0].refIdx[l] >= 0 && mi[idx1].refIdx[l] >= 0 && mi[idx2].refIdx[l] >= 0 )
      {
        // check same refidx and different mv
        if ( mi[idx0].refIdx[l] == mi[idx1].refIdx[l] && mi[idx0].refIdx[l] == mi[idx2].refIdx[l])
        {
          dir |= (l + 1);
          refIdx[l] = mi[idx0].refIdx[l];
        }
      }
    }

#if INTER_LIC
    LICFlag = mi[idx0].usesLIC || mi[idx1].usesLIC || mi[idx2].usesLIC;
#endif
  }

  if ( dir == 0 )
  {
#if JVET_Z0139_HIST_AFF
    return false;
#else
    return;
#endif
  }

  for ( int l = 0; l < 2; l++ )
  {
    if ( dir & (l + 1) )
    {
      for ( int i = 0; i < verNum; i++ )
      {
        cMv[l][verIdx[i]] = mi[verIdx[i]].mv[l];
      }

      // convert to LT, RT[, [LB]]
      switch ( modelIdx )
      {
      case 0: // 0 : LT, RT, LB
        break;

      case 1: // 1 : LT, RT, RB
        cMv[l][2].hor = cMv[l][3].hor + cMv[l][0].hor - cMv[l][1].hor;
        cMv[l][2].ver = cMv[l][3].ver + cMv[l][0].ver - cMv[l][1].ver;
        cMv[l][2].clipToStorageBitDepth();
        break;

      case 2: // 2 : LT, LB, RB
        cMv[l][1].hor = cMv[l][3].hor + cMv[l][0].hor - cMv[l][2].hor;
        cMv[l][1].ver = cMv[l][3].ver + cMv[l][0].ver - cMv[l][2].ver;
        cMv[l][1].clipToStorageBitDepth();
        break;

      case 3: // 3 : RT, LB, RB
        cMv[l][0].hor = cMv[l][1].hor + cMv[l][2].hor - cMv[l][3].hor;
        cMv[l][0].ver = cMv[l][1].ver + cMv[l][2].ver - cMv[l][3].ver;
        cMv[l][0].clipToStorageBitDepth();
        break;

      case 4: // 4 : LT, RT
#if JVET_Z0139_HIST_AFF
      {
        Mv mvLT, mvRT, mvLB;
        mvLT = cMv[l][0];
        mvRT = cMv[l][1];
        int shift = MAX_CU_DEPTH;
        int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;

        iDMvHorX = (mvRT - mvLT).getHor() << (shift - floorLog2(pu.lwidth()));
        iDMvHorY = (mvRT - mvLT).getVer() << (shift - floorLog2(pu.lwidth()));
        iDMvVerX = -iDMvHorY;
        iDMvVerY = iDMvHorX;
        
        int iMvScaleHor = mvLT.getHor() << shift;
        int iMvScaleVer = mvLT.getVer() << shift;
        int horTmp, verTmp;

        horTmp = iMvScaleHor + iDMvVerX * pu.lheight();
        verTmp = iMvScaleVer + iDMvVerY * pu.lheight();
        roundAffineMv(horTmp, verTmp, shift);
        cMv[l][2].hor = horTmp;
        cMv[l][2].ver = verTmp;
        cMv[l][2].clipToStorageBitDepth();
      }
#endif
        break;

      case 5: // 5 : LT, LB
        vx = (cMv[l][0].hor << shift) + ((cMv[l][2].ver - cMv[l][0].ver) << shiftHtoW);
        vy = (cMv[l][0].ver << shift) - ((cMv[l][2].hor - cMv[l][0].hor) << shiftHtoW);
        roundAffineMv( vx, vy, shift );
        cMv[l][1].set( vx, vy );
        cMv[l][1].clipToStorageBitDepth();
        break;

      default:
        CHECK( 1, "Invalid model index!\n" );
        break;
      }
    }
    else
    {
      for ( int i = 0; i < 4; i++ )
      {
        cMv[l][i].hor = 0;
        cMv[l][i].ver = 0;
      }
    }
  }

  for ( int i = 0; i < 3; i++ )
  {
    affMrgType.mvFieldNeighbours[(affMrgType.numValidMergeCand << 1) + 0][i].mv = cMv[0][i];
    affMrgType.mvFieldNeighbours[(affMrgType.numValidMergeCand << 1) + 0][i].refIdx = refIdx[0];

    affMrgType.mvFieldNeighbours[(affMrgType.numValidMergeCand << 1) + 1][i].mv = cMv[1][i];
    affMrgType.mvFieldNeighbours[(affMrgType.numValidMergeCand << 1) + 1][i].refIdx = refIdx[1];
  }
  affMrgType.interDirNeighbours[affMrgType.numValidMergeCand] = dir;
  affMrgType.affineType[affMrgType.numValidMergeCand] = curType;
  affMrgType.BcwIdx[affMrgType.numValidMergeCand] = (dir == 3) ? bcwIdx : BCW_DEFAULT;
#if INTER_LIC
  affMrgType.LICFlags[affMrgType.numValidMergeCand] = LICFlag;
#endif
#if JVET_Z0139_HIST_AFF
  if (!checkLastAffineMergeCandRedundancy(pu, affMrgType))
  {
    return false;
  }
  return true;
#else
  affMrgType.numValidMergeCand++;
  return;
#endif
}

#if JVET_Z0139_NA_AFF
int PU::getMvDiffThresholdByWidthAndHeight(const PredictionUnit &pu, bool width)
{
  uint32_t numPixels = (width ? pu.lwidth() : pu.lheight());
  if (numPixels <= 8)
  {
    return (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 4;
  }
  else if (numPixels <= 32)
  {
    return (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 3;
  }
  else if (numPixels <= 64)
  {
    return (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 2;
  }
  else
  {
    return (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 1;
  }
}

bool PU::addNonAdjCstAffineMVPCandUnscaled(const PredictionUnit &pu, const RefPicList &refPicList, const int &refIdx, AffineAMVPInfo &affiAmvpInfo)
{
  const Position posLT[3] = { pu.Y().topLeft().offset( -1, -1 ), pu.Y().topLeft().offset( 0, -1 ), pu.Y().topLeft().offset( -1, 0 ) };
  const Position posRT[2] = { pu.Y().topRight().offset(0, -1), pu.Y().topRight().offset(1, -1) };
  const Position posLB[2] = { pu.Y().bottomLeft().offset(-1, 0), pu.Y().bottomLeft().offset(-1, 1) };
  const unsigned plevel   = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;

  for (int i = 1; i < (AFF_NON_ADJACENT_DIST + 1); i++)
  {
    MotionInfo miNew[3];
    Position   posNew[3];
    bool       isAvailableNew[3] = { false, false, false };

    const Position        posTRNew     = Position(posRT[0].x, posRT[0].y - (i * ((int) pu.Y().height)));
    const PredictionUnit *puNeighTRNew = pu.cs->getPURestricted(posTRNew, pu, pu.chType);
#if JVET_Y0065_GPM_INTRA
    if (puNeighTRNew && CU::isInter(*puNeighTRNew->cu) && puNeighTRNew->getMotionInfo(posTRNew).isInter && PU::isDiffMER(pu.lumaPos(), posTRNew, plevel))
#else
    if (puNeighTRNew && CU::isInter(*puNeighTRNew->cu) && PU::isDiffMER(pu.lumaPos(), posTRNew, plevel))
#endif
    {
      isAvailableNew[1] = true;
      miNew[1]          = puNeighTRNew->getMotionInfo(posTRNew);
      posNew[1]         = posTRNew;
      CHECK(posTRNew.x < 0 || posTRNew.y < 0, "posTRNew < 0");
    }

    for (int j = 1; j < (AFF_NON_ADJACENT_DIST + 1); j++)
    {
      isAvailableNew[0]                  = false;
      isAvailableNew[2]                  = false;
      const Position        posLBNew     = Position(posLB[0].x - (j * ((int) pu.Y().width)), posLB[0].y);
      const PredictionUnit *puNeighLBNew = pu.cs->getPURestricted(posLBNew, pu, pu.chType);
#if JVET_Y0065_GPM_INTRA
      if (puNeighLBNew && CU::isInter(*puNeighLBNew->cu) && puNeighLBNew->getMotionInfo(posLBNew).isInter && PU::isDiffMER(pu.lumaPos(), posLBNew, plevel))
#else
      if (puNeighLBNew && CU::isInter(*puNeighLBNew->cu) && PU::isDiffMER(pu.lumaPos(), posLBNew, plevel))
#endif
      {
        isAvailableNew[2] = true;
        miNew[2]          = puNeighLBNew->getMotionInfo(posLBNew);
        posNew[2]         = posLBNew;
        CHECK(posLBNew.x < 0 || posLBNew.y < 0, "posLBNew < 0");
      }

      PosType posX = isAvailableNew[2] ? posNew[2].x : (posLT[0].x - (j * ((int) pu.Y().width)));
      PosType posY = isAvailableNew[1] ? posNew[1].y : (posLT[0].y - (i * ((int) pu.Y().height)));
      if (posX < 0)
      {
        posX = isAvailableNew[1] ? posLT[1].x : -1;
      }
      if (posY < 0)
      {
        posY = isAvailableNew[2] ? posLT[2].y : -1;
      }

      const Position        posLTNew     = Position(posX, posY);
      const PredictionUnit *puNeighLTNew = pu.cs->getPURestricted(posLTNew, pu, pu.chType);
#if JVET_Y0065_GPM_INTRA
      if (puNeighLTNew && CU::isInter(*puNeighLTNew->cu) && puNeighLTNew->getMotionInfo(posLTNew).isInter && PU::isDiffMER(pu.lumaPos(), posLTNew, plevel))
#else
      if (puNeighLTNew && CU::isInter(*puNeighLTNew->cu) && PU::isDiffMER(pu.lumaPos(), posLTNew, plevel))
#endif
      {
        isAvailableNew[0] = true;
        miNew[0]          = puNeighLTNew->getMotionInfo(posLTNew);
        posNew[0]         = posLTNew;
        CHECK(posLTNew.x < 0 || posLTNew.y < 0, "posLTNew < 0");
      }

      if (addNonAdjCstAffineMVPConstructedCPMV(pu, miNew, isAvailableNew, posNew, refPicList, refIdx, affiAmvpInfo))
      {
        if (affiAmvpInfo.numCand >= AMVP_MAX_NUM_CANDS)
        {
          return true;
        }
      }
    }
  }

  return false;
}

int PU::getNonAdjAffParaDivFun(int num1, int num2)
{
  int divTable[16] = { 0, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 0 };
  int x = floorLog2(num2);
  int normNum1 = (num2 << 4 >> x) & 15;
  int v = divTable[normNum1] | 8;
  x += (normNum1 != 0);
  int shift = 13 - x;
  int retVal = 0;
  if (shift < 0)
  {
    shift = -shift;
    int add = (1 << (shift - 1));
    retVal = (num1 * v + add) >> shift;
  }
  else
  {
    retVal = (num1 * v) << shift;
  }
  return (retVal>>(16-MAX_CU_DEPTH));
}

bool PU::addNonAdjCstAffineMVPConstructedCPMV( const PredictionUnit &pu, MotionInfo miNew[3], bool isAvaNew[3], Position pos[3], const RefPicList &refPicList, const int &refIdx, AffineAMVPInfo &affiAmvpInfo)
{
  if (!isAvaNew[0] || (!isAvaNew[1] && !isAvaNew[2]))
  {
    return false;
  }
  
  const int        currRefPOC    = pu.cs->slice->getRefPic(refPicList, refIdx)->getPOC();
  const RefPicList refPicList2nd = (refPicList == REF_PIC_LIST_0) ? REF_PIC_LIST_1 : REF_PIC_LIST_0;

  int  shift  = MAX_CU_DEPTH;
  int posNeiX = pos[0].x;
  int posNeiY = pos[0].y;
  int posCurX = pu.Y().pos().x;
  int posCurY = pu.Y().pos().y;

  int curW = pu.Y().width;
  int curH = pu.Y().height;
  int neiW = pos[1].x - pos[0].x;
  int neiH = pos[2].y - pos[0].y;
  if (!isAvaNew[1])
  {
    neiW = pu.Y().topRight().x - pos[0].x;
  }
  if (!isAvaNew[2])
  {
    neiH = pu.Y().bottomLeft().y - pos[0].y;
  }

  bool isConverted = false;

  for (int predictorSource = 0; predictorSource < 2; predictorSource++)
  {
    const RefPicList refListindex     = (predictorSource == 0) ? refPicList : refPicList2nd;
    const int        neibRefIdx = miNew[0].refIdx[refListindex];
    if (neibRefIdx < 0 || pu.cu->slice->getRefPOC(refListindex, neibRefIdx) != currRefPOC)
    {
      continue;
    }
    for (int modelIdx = 0; modelIdx < 3; modelIdx++)
    {
      if ((modelIdx == 1 && !isAvaNew[1]) || (modelIdx == 2 && !isAvaNew[2]))
      {
        continue;
      }
      Mv  outputAffineMv[3];
      Mv  mvLT, mvRT, mvLB;
      int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;
      int horTmp, verTmp;

      mvLT = miNew[0].mv[refListindex];
      mvRT = miNew[1].mv[refListindex];
      mvLB = miNew[2].mv[refListindex];

      int iMvScaleHor = mvLT.getHor() << shift;
      int iMvScaleVer = mvLT.getVer() << shift;
      iDMvHorX        = getNonAdjAffParaDivFun((mvRT - mvLT).getHor(), neiW);
      iDMvHorY        = getNonAdjAffParaDivFun((mvRT - mvLT).getVer(), neiW);
      iDMvVerX        = getNonAdjAffParaDivFun((mvLB - mvLT).getHor(), neiH);
      iDMvVerY        = getNonAdjAffParaDivFun((mvLB - mvLT).getVer(), neiH);

      if (!modelIdx && isAvaNew[0] && isAvaNew[1] && isAvaNew[2]
          && miNew[0].refIdx[refListindex] == miNew[1].refIdx[refListindex] 
          && miNew[0].refIdx[refListindex] == miNew[2].refIdx[refListindex])
      {
      }
      else if (modelIdx == 1 && isAvaNew[0] && isAvaNew[1]
               && miNew[0].refIdx[refListindex] == miNew[1].refIdx[refListindex])
      {
        iDMvVerX = -iDMvHorY;
        iDMvVerY = iDMvHorX;
      }
      else if (modelIdx == 2 && isAvaNew[0] && isAvaNew[2]
               && miNew[0].refIdx[refListindex] == miNew[2].refIdx[refListindex])
      {
        iDMvHorX = iDMvVerY;
        iDMvHorY = -iDMvVerX;
      }
      else
      {
        continue;
      }

      horTmp = iMvScaleHor + iDMvHorX * (posCurX - posNeiX) + iDMvVerX * (posCurY - posNeiY);
      verTmp = iMvScaleVer + iDMvHorY * (posCurX - posNeiX) + iDMvVerY * (posCurY - posNeiY);
      roundAffineMv(horTmp, verTmp, shift);
      outputAffineMv[0].hor = horTmp;
      outputAffineMv[0].ver = verTmp;
      outputAffineMv[0].clipToStorageBitDepth();

      horTmp = iMvScaleHor + iDMvHorX * (posCurX + curW - posNeiX) + iDMvVerX * (posCurY - posNeiY);
      verTmp = iMvScaleVer + iDMvHorY * (posCurX + curW - posNeiX) + iDMvVerY * (posCurY - posNeiY);
      roundAffineMv(horTmp, verTmp, shift);
      outputAffineMv[1].hor = horTmp;
      outputAffineMv[1].ver = verTmp;
      outputAffineMv[1].clipToStorageBitDepth();

      if (pu.cu->affineType == AFFINEMODEL_6PARAM)
      {
        horTmp = iMvScaleHor + iDMvHorX * (posCurX - posNeiX) + iDMvVerX * (posCurY + curH - posNeiY);
        verTmp = iMvScaleVer + iDMvHorY * (posCurX - posNeiX) + iDMvVerY * (posCurY + curH - posNeiY);
        roundAffineMv(horTmp, verTmp, shift);
        outputAffineMv[2].hor = horTmp;
        outputAffineMv[2].ver = verTmp;
        outputAffineMv[2].clipToStorageBitDepth();
      }

      outputAffineMv[0].roundAffinePrecInternal2Amvr(pu.cu->imv);
      outputAffineMv[1].roundAffinePrecInternal2Amvr(pu.cu->imv);
      if (pu.cu->affineType == AFFINEMODEL_6PARAM)
      {
        outputAffineMv[2].roundAffinePrecInternal2Amvr(pu.cu->imv);
      }
      affiAmvpInfo.mvCandLT[affiAmvpInfo.numCand] = outputAffineMv[0];
      affiAmvpInfo.mvCandRT[affiAmvpInfo.numCand] = outputAffineMv[1];
      affiAmvpInfo.mvCandLB[affiAmvpInfo.numCand] = Mv();
      if (pu.cu->affineType == AFFINEMODEL_6PARAM)
      {
        affiAmvpInfo.mvCandLB[affiAmvpInfo.numCand] = outputAffineMv[2];
      }
      if (!checkLastAffineAMVPCandRedundancy(pu, affiAmvpInfo))
      {
        affiAmvpInfo.mvCandLT[affiAmvpInfo.numCand] = Mv();
        affiAmvpInfo.mvCandRT[affiAmvpInfo.numCand] = Mv();
        affiAmvpInfo.mvCandLB[affiAmvpInfo.numCand] = Mv();
        continue;
      }
      affiAmvpInfo.numCand++;
      isConverted = true;
      if (affiAmvpInfo.numCand >= AMVP_MAX_NUM_CANDS)
      {
        return true;
      }
    }
  }

  return isConverted;
}

bool PU::xCPMVSimCheck(const PredictionUnit &pu, AffineMergeCtx &affMrgCtx, Mv curCpmv[2][3], unsigned char curdir, int8_t curRefIdx[2], EAffineModel curType, int bcwIdx, bool LICFlag)
{
  if (affMrgCtx.numValidMergeCand == 0 || affMrgCtx.mergeType[affMrgCtx.numValidMergeCand - 1] != MRG_TYPE_DEFAULT_N)
  {
    return false;
  }
  int iDMvHorX[2], iDMvHorY[2], iDMvVerX[2], iDMvVerY[2], iMvScaleHor[2], iMvScaleVer[2];
  int iDMvHorXCand, iDMvHorYCand, iDMvVerXCand, iDMvVerYCand, iMvScaleHorCand, iMvScaleVerCand;
  Mv  mvLT, mvRT, mvLB;

  const int mvDifThldWidth  = getMvDiffThresholdByWidthAndHeight(pu, true);
  const int mvDifThldHeight = getMvDiffThresholdByWidthAndHeight(pu, false);

  for (int l = 0; l < 2; l++)
  {
    if (curdir & (l + 1))
    {
      mvLT        = curCpmv[l][0];
      mvRT        = curCpmv[l][1];
      mvLB        = curCpmv[l][2];
      iDMvHorX[l] = (mvRT - mvLT).getHor();
      iDMvHorY[l] = (mvRT - mvLT).getVer();
      if (curType == AFFINEMODEL_6PARAM)
      {
        iDMvVerX[l] = (mvLB - mvLT).getHor();
        iDMvVerY[l] = (mvLB - mvLT).getVer();
      }
      else
      {
        iDMvVerX[l] = -iDMvHorY[l];
        iDMvVerY[l] = iDMvHorX[l];
      }

      iMvScaleHor[l] = mvLT.getHor();
      iMvScaleVer[l] = mvLT.getVer();
    }
    else
    {
      iDMvHorX[l]    = iDMvHorY[l] = 0;
      iDMvVerX[l]    = iDMvVerY[l] = 0;
      iMvScaleHor[l] = iMvScaleVer[l] = 0;
    }
  }


  for (uint32_t ui = 0; ui < affMrgCtx.numValidMergeCand; ui++)
  {
    bool isSimilar = true;
    unsigned char candDir    = affMrgCtx.interDirNeighbours[ui];
    int8_t        candRefIdx[2];
    if ((affMrgCtx.mergeType[ui] != MRG_TYPE_DEFAULT_N) || (candDir != curdir))
    {
      continue;
    }
    candRefIdx[0] = affMrgCtx.mvFieldNeighbours[(ui << 1) + 0][0].refIdx;
    candRefIdx[1] = affMrgCtx.mvFieldNeighbours[(ui << 1) + 1][0].refIdx;

    if ((candDir == 3 && (curRefIdx[0] != candRefIdx[0] || curRefIdx[1] != candRefIdx[1]))
      || (candDir == 1 && curRefIdx[0] != candRefIdx[0]) || (candDir == 2 && curRefIdx[1] != candRefIdx[1]))
    {
      continue;
    }
    EAffineModel candAffType = affMrgCtx.affineType[ui];
    for (int l = 0; l < 2; l++) 
    {
      if (curdir & (l + 1)) 
      {
        mvLT = affMrgCtx.mvFieldNeighbours[(ui << 1) + l][0].mv;
        mvRT = affMrgCtx.mvFieldNeighbours[(ui << 1) + l][1].mv;
        mvLB = affMrgCtx.mvFieldNeighbours[(ui << 1) + l][2].mv;
        iDMvHorXCand = (mvRT - mvLT).getHor();
        iDMvHorYCand = (mvRT - mvLT).getVer();
        if (candAffType == AFFINEMODEL_6PARAM)
        {
          iDMvVerXCand = (mvLB - mvLT).getHor();
          iDMvVerYCand = (mvLB - mvLT).getVer();
        }
        else
        {
          iDMvVerXCand = -iDMvHorYCand;
          iDMvVerYCand = iDMvHorXCand;
        }

        iMvScaleHorCand = mvLT.getHor();
        iMvScaleVerCand = mvLT.getVer();
      
        int diffHorX = iDMvHorX[l] - iDMvHorXCand;
        int diffHorY = iDMvHorY[l] - iDMvHorYCand;
        int diffVerX = iDMvVerX[l] - iDMvVerXCand;
        int diffVerY = iDMvVerY[l] - iDMvVerYCand;
        int diffiMvScaleHor = iMvScaleHor[l] - iMvScaleHorCand;
        int diffiMvScaleVer = iMvScaleVer[l] - iMvScaleVerCand;
        
        if (abs(diffHorX) >= mvDifThldWidth || abs(diffHorY) >= mvDifThldWidth || abs(diffVerX) >= mvDifThldHeight
          || abs(diffVerY) >= mvDifThldHeight || abs(diffiMvScaleHor) >= 1
          || abs(diffiMvScaleVer) >= 1)
        {
          isSimilar = false;
        }
        if (!isSimilar)
        {
          break;
        }
      }
    }
    if (isSimilar)
    {
      return true;
    }
  }

  return false;
}

bool PU::addNonAdjAffineConstructedCPMV(const PredictionUnit &pu, MotionInfo miNew[4], bool isAvaNew[4], Position pos[4], int8_t bcwId, AffineMergeCtx &affMrgCtx, int mrgCandIdx)
{
  if (!isAvaNew[0] || (!isAvaNew[1] && !isAvaNew[2]))
  {
    return true;
  }
  int  shift        = MAX_CU_DEPTH;

  int posNeiX = pos[0].x;
  int posNeiY = pos[0].y;
  int posCurX = pu.Y().pos().x;
  int posCurY = pu.Y().pos().y;

  int curW = pu.Y().width;
  int curH = pu.Y().height;
  int neiW = pos[1].x - pos[0].x;
  int neiH = pos[2].y - pos[0].y;
  if (!isAvaNew[1])
  {
    neiW = pu.Y().topRight().x - pos[0].x;
  }
  if (!isAvaNew[2])
  {
    neiH = pu.Y().bottomLeft().y - pos[0].y;
  }

  bool isConverted = false;
  for (int modelIdx = 0; modelIdx < 3; modelIdx++)
  {
    if ((modelIdx == 1 && !isAvaNew[1]) || (modelIdx == 2 && !isAvaNew[2]))
    {
      continue;
    }
    Mv cMv[2][3];
    int8_t refIdx[2] = { -1, -1 };
    int dir = 0;
    EAffineModel curType = AFFINEMODEL_6PARAM;
    bool bLICFlag = false;

    for (int l = 0; l < 2; l++)
    {
      Mv mvLT, mvRT, mvLB;
      int iDMvHorX, iDMvHorY, iDMvVerX, iDMvVerY;
      int horTmp, verTmp;

      mvLT = miNew[0].mv[l];
      mvRT = miNew[1].mv[l];
      mvLB = miNew[2].mv[l];

      int iMvScaleHor = mvLT.getHor() << shift;
      int iMvScaleVer = mvLT.getVer() << shift;
      iDMvHorX        = getNonAdjAffParaDivFun((mvRT - mvLT).getHor(), neiW);
      iDMvHorY        = getNonAdjAffParaDivFun((mvRT - mvLT).getVer(), neiW);
      iDMvVerX        = getNonAdjAffParaDivFun((mvLB - mvLT).getHor(), neiH);
      iDMvVerY        = getNonAdjAffParaDivFun((mvLB - mvLT).getVer(), neiH);

      if (!modelIdx && isAvaNew[0] && isAvaNew[1] && isAvaNew[2] && miNew[0].refIdx[l] >= 0
          && miNew[0].refIdx[l] == miNew[1].refIdx[l] && miNew[0].refIdx[l] == miNew[2].refIdx[l])
      {
#if INTER_LIC
        bLICFlag = bLICFlag || miNew[0].usesLIC || miNew[1].usesLIC || miNew[2].usesLIC;
#endif
      }
      else if (modelIdx == 1 && isAvaNew[0] && isAvaNew[1] && miNew[0].refIdx[l] >= 0
               && miNew[0].refIdx[l] == miNew[1].refIdx[l])
      {
        iDMvVerX = -iDMvHorY;
        iDMvVerY = iDMvHorX;
#if INTER_LIC
        bLICFlag = bLICFlag || miNew[0].usesLIC || miNew[1].usesLIC;
#endif
      }
      else if (modelIdx == 2 && isAvaNew[0] && isAvaNew[2] && miNew[0].refIdx[l] >= 0
               && miNew[0].refIdx[l] == miNew[2].refIdx[l])
      {
        iDMvHorX = iDMvVerY;
        iDMvHorY = -iDMvVerX;
#if INTER_LIC
        bLICFlag = bLICFlag || miNew[0].usesLIC || miNew[2].usesLIC;
#endif
      }
      else
      {
        continue;
      }

      dir |= (l + 1);
      refIdx[l] = miNew[0].refIdx[l];

      horTmp = iMvScaleHor + iDMvHorX * (posCurX - posNeiX) + iDMvVerX * (posCurY - posNeiY);
      verTmp = iMvScaleVer + iDMvHorY * (posCurX - posNeiX) + iDMvVerY * (posCurY - posNeiY);
      roundAffineMv(horTmp, verTmp, shift);
      cMv[l][0].hor = horTmp;
      cMv[l][0].ver = verTmp;
      cMv[l][0].clipToStorageBitDepth();

      horTmp = iMvScaleHor + iDMvHorX * (posCurX + curW - posNeiX) + iDMvVerX * (posCurY - posNeiY);
      verTmp = iMvScaleVer + iDMvHorY * (posCurX + curW - posNeiX) + iDMvVerY * (posCurY - posNeiY);
      roundAffineMv(horTmp, verTmp, shift);
      cMv[l][1].hor = horTmp;
      cMv[l][1].ver = verTmp;
      cMv[l][1].clipToStorageBitDepth();

      {
        horTmp = iMvScaleHor + iDMvHorX * (posCurX - posNeiX) + iDMvVerX * (posCurY + curH - posNeiY);
        verTmp = iMvScaleVer + iDMvHorY * (posCurX - posNeiX) + iDMvVerY * (posCurY + curH - posNeiY);
        roundAffineMv(horTmp, verTmp, shift);
        cMv[l][2].hor = horTmp;
        cMv[l][2].ver = verTmp;
        cMv[l][2].clipToStorageBitDepth();
      }
    }

    if (!dir)
    {
      continue;
    }
    if (xCPMVSimCheck(pu, affMrgCtx, cMv, dir, refIdx, curType, (dir == 3) ? bcwId : BCW_DEFAULT, (dir != 3) ? bLICFlag : false))
    {
      continue;
    }
    isConverted = true;
    for (int i = 0; i < 3; i++)
    {
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][i].mv      = cMv[0][i];
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][i].refIdx = refIdx[0];

      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][i].mv      = cMv[1][i];
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][i].refIdx = refIdx[1];
    }
    affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand]  = dir;
    affMrgCtx.affineType[affMrgCtx.numValidMergeCand]         = curType;
    affMrgCtx.mergeType[affMrgCtx.numValidMergeCand]           = MRG_TYPE_DEFAULT_N;
    affMrgCtx.BcwIdx[affMrgCtx.numValidMergeCand]              = (dir == 3) ? bcwId : BCW_DEFAULT;
#if INTER_LIC
    affMrgCtx.LICFlags[affMrgCtx.numValidMergeCand] = (dir != 3) ? bLICFlag : false;
#endif

    affMrgCtx.numValidMergeCand++;

    if (affMrgCtx.numValidMergeCand != 0 && affMrgCtx.numValidMergeCand - 1 == mrgCandIdx)
    {
      return false;
    }

    if (affMrgCtx.numValidMergeCand == affMrgCtx.maxNumMergeCand)
    {
      return false;
    }
    break;
  }

  return (!isConverted);
}

void PU::getNonAdjCstMergeCand(const PredictionUnit &pu, AffineMergeCtx &affMrgCtx, const int mrgCandIdx, bool isInitialized)
{
  const CodingStructure &cs    = *pu.cs;
  const Slice &          slice = *pu.cs->slice;
  const uint32_t maxNumAffineMergeCand = slice.getPicHeader()->getMaxNumAffineMergeCand();
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
  const Position posLT[3] = { pu.Y().topLeft().offset( -1, -1 ), pu.Y().topLeft().offset( 0, -1 ), pu.Y().topLeft().offset( -1, 0 ) };
  const Position posRT[2] = { pu.Y().topRight().offset(0, -1), pu.Y().topRight().offset(1, -1) };
  const Position posLB[2] = { pu.Y().bottomLeft().offset(-1, 0), pu.Y().bottomLeft().offset(-1, 1) };
  if (!isInitialized)
  {
    for (int i = 0; i < maxNumAffineMergeCand; i++)
    {
      for (int mvNum = 0; mvNum < 3; mvNum++)
      {
        affMrgCtx.mvFieldNeighbours[(i << 1) + 0][mvNum].setMvField(Mv(), -1);
        affMrgCtx.mvFieldNeighbours[(i << 1) + 1][mvNum].setMvField(Mv(), -1);
      }
      affMrgCtx.interDirNeighbours[i] = 0;
      affMrgCtx.affineType[i]         = AFFINEMODEL_4PARAM;
      affMrgCtx.mergeType[i]          = MRG_TYPE_DEFAULT_N;
      affMrgCtx.BcwIdx[i]             = BCW_DEFAULT;
#if INTER_LIC
      affMrgCtx.LICFlags[i] = false;
#endif
    }

    affMrgCtx.numValidMergeCand = 0;
    affMrgCtx.maxNumMergeCand   = maxNumAffineMergeCand;
  }

  for (int i = 1; i < (AFF_NON_ADJACENT_DIST + 1); i++)
  {
    MotionInfo miNew[4];
    Position   posNew[4];
    bool       isAvailableNew[4] = { false, false, false, false };
    int8_t     neighBcwNew[2]    = { BCW_DEFAULT, BCW_DEFAULT };

    const Position        posTRNew     = Position(posRT[0].x, posRT[0].y - (i * ((int) pu.Y().height)));
    const PredictionUnit *puNeighTRNew = cs.getPURestricted(posTRNew, pu, pu.chType);
#if JVET_Y0065_GPM_INTRA
    if (puNeighTRNew && CU::isInter(*puNeighTRNew->cu) && puNeighTRNew->getMotionInfo(posTRNew).isInter && PU::isDiffMER(pu.lumaPos(), posTRNew, plevel))
#else
    if (puNeighTRNew && CU::isInter(*puNeighTRNew->cu) && PU::isDiffMER(pu.lumaPos(), posTRNew, plevel))
#endif
    {
      isAvailableNew[1] = true;
      miNew[1]          = puNeighTRNew->getMotionInfo(posTRNew);
      neighBcwNew[1]    = puNeighTRNew->cu->BcwIdx;
      posNew[1]         = posTRNew;
      CHECK(posTRNew.x < 0 || posTRNew.y < 0, "posTRNew < 0");
    }

    for (int j = 1; j < (AFF_NON_ADJACENT_DIST + 1); j++)
    {
      isAvailableNew[0]                  = false;
      isAvailableNew[2]                  = false;
      const Position        posLBNew     = Position(posLB[0].x - (j * ((int) pu.Y().width)), posLB[0].y);
      const PredictionUnit *puNeighLBNew = cs.getPURestricted(posLBNew, pu, pu.chType);
#if JVET_Y0065_GPM_INTRA
      if (puNeighLBNew && CU::isInter(*puNeighLBNew->cu) && puNeighLBNew->getMotionInfo(posLBNew).isInter && PU::isDiffMER(pu.lumaPos(), posLBNew, plevel))
#else
      if (puNeighLBNew && CU::isInter(*puNeighLBNew->cu) && PU::isDiffMER(pu.lumaPos(), posLBNew, plevel))
#endif
      {
        isAvailableNew[2] = true;
        miNew[2]          = puNeighLBNew->getMotionInfo(posLBNew);
        posNew[2]         = posLBNew;
        CHECK(posLBNew.x < 0 || posLBNew.y < 0, "posLBNew < 0");
      }

      PosType posX = isAvailableNew[2] ? posNew[2].x : (posLT[0].x - (j * ((int) pu.Y().width)));
      PosType posY = isAvailableNew[1] ? posNew[1].y : (posLT[0].y - (i * ((int) pu.Y().height)));
      if (posX < 0)
      {
        posX = isAvailableNew[1] ? posLT[1].x : -1;
      }
      if (posY < 0)
      {
        posY = isAvailableNew[2] ? posLT[2].y : -1;
      }

      const Position        posLTNew     = Position(posX, posY);
      const PredictionUnit *puNeighLTNew = cs.getPURestricted(posLTNew, pu, pu.chType);
#if JVET_Y0065_GPM_INTRA
      if (puNeighLTNew && CU::isInter(*puNeighLTNew->cu) && puNeighLTNew->getMotionInfo(posLTNew).isInter && PU::isDiffMER(pu.lumaPos(), posLTNew, plevel))
#else
      if (puNeighLTNew && CU::isInter(*puNeighLTNew->cu) && PU::isDiffMER(pu.lumaPos(), posLTNew, plevel))
#endif
      {
        isAvailableNew[0] = true;
        miNew[0]          = puNeighLTNew->getMotionInfo(posLTNew);
        neighBcwNew[0]    = puNeighLTNew->cu->BcwIdx;
        posNew[0]         = posLTNew;
        CHECK(posLTNew.x < 0 || posLTNew.y < 0, "posLTNew < 0");
      }

      if (!addNonAdjAffineConstructedCPMV(pu, miNew, isAvailableNew, posNew, neighBcwNew[0], affMrgCtx, mrgCandIdx))
      {
        if (affMrgCtx.numValidMergeCand != 0 && affMrgCtx.numValidMergeCand - 1 == mrgCandIdx)
        {
          return;
        }

        if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
        {
          return;
        }
      }
    }
  }
}
#endif

const int getAvailableAffineNeighboursForLeftPredictor( const PredictionUnit &pu, const PredictionUnit* npu[] 
#if JVET_Z0139_HIST_AFF
  , int neiIdx[]
#endif
)
{
  const Position posLB = pu.Y().bottomLeft();
  int num = 0;
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;

  const PredictionUnit *puLeftBottom = pu.cs->getPURestricted( posLB.offset( -1, 1 ), pu, pu.chType );
  if (puLeftBottom && puLeftBottom->cu->affine && puLeftBottom->mergeType == MRG_TYPE_DEFAULT_N
      && PU::isDiffMER(pu.lumaPos(), posLB.offset(-1, 1), plevel))
  {
#if JVET_Z0139_HIST_AFF
    neiIdx[num] = 3;
#endif
    npu[num++] = puLeftBottom;
    return num;
  }

  const PredictionUnit* puLeft = pu.cs->getPURestricted( posLB.offset( -1, 0 ), pu, pu.chType );
  if (puLeft && puLeft->cu->affine && puLeft->mergeType == MRG_TYPE_DEFAULT_N
      && PU::isDiffMER(pu.lumaPos(), posLB.offset(-1, 0), plevel))
  {
#if JVET_Z0139_HIST_AFF
    neiIdx[num] = 0;
#endif
    npu[num++] = puLeft;
    return num;
  }

  return num;
}

const int getAvailableAffineNeighboursForAbovePredictor( const PredictionUnit &pu, const PredictionUnit* npu[], int numAffNeighLeft 
#if JVET_Z0139_HIST_AFF
  , int neiIdx[]
#endif
)
{
  const Position posLT = pu.Y().topLeft();
  const Position posRT = pu.Y().topRight();
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
  int num = numAffNeighLeft;

  const PredictionUnit* puAboveRight = pu.cs->getPURestricted( posRT.offset( 1, -1 ), pu, pu.chType );
  if (puAboveRight && puAboveRight->cu->affine && puAboveRight->mergeType == MRG_TYPE_DEFAULT_N
      && PU::isDiffMER(pu.lumaPos(), posRT.offset(1, -1), plevel))
  {
#if JVET_Z0139_HIST_AFF
    neiIdx[num] = 2;
#endif
    npu[num++] = puAboveRight;
    return num;
  }

  const PredictionUnit* puAbove = pu.cs->getPURestricted( posRT.offset( 0, -1 ), pu, pu.chType );
  if (puAbove && puAbove->cu->affine && puAbove->mergeType == MRG_TYPE_DEFAULT_N
      && PU::isDiffMER(pu.lumaPos(), posRT.offset(0, -1), plevel))
  {
#if JVET_Z0139_HIST_AFF
    neiIdx[num] = 1;
#endif
    npu[num++] = puAbove;
    return num;
  }

  const PredictionUnit *puAboveLeft = pu.cs->getPURestricted( posLT.offset( -1, -1 ), pu, pu.chType );
  if (puAboveLeft && puAboveLeft->cu->affine && puAboveLeft->mergeType == MRG_TYPE_DEFAULT_N
      && PU::isDiffMER(pu.lumaPos(), posLT.offset(-1, -1), plevel))
  {
#if JVET_Z0139_HIST_AFF
    neiIdx[num] = 4;
#endif
    npu[num++] = puAboveLeft;
    return num;
  }

  return num;
}

void PU::getAffineMergeCand( const PredictionUnit &pu, AffineMergeCtx& affMrgCtx,
#if !AFFINE_MMVD
                             const int mrgCandIdx
#else
                                   int mrgCandIdx, bool isAfMmvd
#endif
#if JVET_Z0139_NA_AFF && JVET_W0090_ARMC_TM
    , bool isZeroCandIdx
#endif
)
{
  const CodingStructure &cs = *pu.cs;
  const Slice &slice = *pu.cs->slice;
  const uint32_t maxNumAffineMergeCand = slice.getPicHeader()->getMaxNumAffineMergeCand();
  const unsigned plevel = pu.cs->sps->getLog2ParallelMergeLevelMinus2() + 2;
#if JVET_Z0118_GDR
  bool isClean = pu.cs->isClean(pu.cu->Y().bottomRight(), CHANNEL_TYPE_LUMA);
#endif

  for ( int i = 0; i < maxNumAffineMergeCand; i++ )
  {
    for ( int mvNum = 0; mvNum < 3; mvNum++ )
    {
      affMrgCtx.mvFieldNeighbours[(i << 1) + 0][mvNum].setMvField( Mv(), -1 );
      affMrgCtx.mvFieldNeighbours[(i << 1) + 1][mvNum].setMvField( Mv(), -1 );
    }
    affMrgCtx.interDirNeighbours[i] = 0;
    affMrgCtx.affineType[i] = AFFINEMODEL_4PARAM;
    affMrgCtx.mergeType[i] = MRG_TYPE_DEFAULT_N;
    affMrgCtx.BcwIdx[i] = BCW_DEFAULT;
#if INTER_LIC
    affMrgCtx.LICFlags[i] = false;
#endif
  }

  affMrgCtx.numValidMergeCand = 0;
  affMrgCtx.maxNumMergeCand = maxNumAffineMergeCand;

  bool sbTmvpEnableFlag = slice.getSPS()->getSbTMVPEnabledFlag()
                          && !(slice.getPOC() == slice.getRefPic(REF_PIC_LIST_0, 0)->getPOC() && slice.isIRAP());
  bool isAvailableSubPu = false;
  if (sbTmvpEnableFlag && slice.getPicHeader()->getEnableTMVPFlag())
  {
    MergeCtx mrgCtx = *affMrgCtx.mrgCtx;
    bool tmpLICFlag = false;

#if JVET_Z0139_NA_AFF && JVET_W0090_ARMC_TM
    if (cs.pcv->isEncoder || isZeroCandIdx)
    {
#endif
    CHECK( mrgCtx.subPuMvpMiBuf.area() == 0 || !mrgCtx.subPuMvpMiBuf.buf, "Buffer not initialized" );
    mrgCtx.subPuMvpMiBuf.fill( MotionInfo() );
#if JVET_Z0139_NA_AFF && JVET_W0090_ARMC_TM
    }
#endif

    int pos = 0;
    // Get spatial MV
    const Position posCurLB = pu.Y().bottomLeft();
    MotionInfo miLeft;

    //left
    const PredictionUnit* puLeft = cs.getPURestricted( posCurLB.offset( -1, 0 ), pu, pu.chType );
#if JVET_Y0065_GPM_INTRA
    const bool isAvailableA1 = puLeft && isDiffMER(pu.lumaPos(), posCurLB.offset(-1, 0), plevel) && pu.cu != puLeft->cu && CU::isInter( *puLeft->cu ) && puLeft->getMotionInfo( posCurLB.offset( -1, 0 ) ).isInter;
#else
    const bool isAvailableA1 = puLeft && isDiffMER(pu.lumaPos(), posCurLB.offset(-1, 0), plevel) && pu.cu != puLeft->cu && CU::isInter( *puLeft->cu );
#endif
    if ( isAvailableA1 )
    {
#if INTER_LIC
      affMrgCtx.mrgCtx->LICFlags[pos] = false;
#endif

      miLeft = puLeft->getMotionInfo( posCurLB.offset( -1, 0 ) );
      // get Inter Dir
      mrgCtx.interDirNeighbours[pos] = miLeft.interDir;

      // get Mv from Left
      mrgCtx.mvFieldNeighbours[pos << 1].setMvField( miLeft.mv[0], miLeft.refIdx[0] );

      if ( slice.isInterB() )
      {
        mrgCtx.mvFieldNeighbours[(pos << 1) + 1].setMvField( miLeft.mv[1], miLeft.refIdx[1] );
      }
      pos++;
    }

    mrgCtx.numValidMergeCand = pos;

#if JVET_Z0139_NA_AFF && JVET_W0090_ARMC_TM
    isAvailableSubPu = getInterMergeSubPuMvpCand(pu, mrgCtx, tmpLICFlag, pos, (!isZeroCandIdx && !cs.pcv->isEncoder) ? 1: 0);
#else
    isAvailableSubPu = getInterMergeSubPuMvpCand(pu, mrgCtx, tmpLICFlag, pos, 0);
#endif
    if ( isAvailableSubPu )
    {
      for ( int mvNum = 0; mvNum < 3; mvNum++ )
      {
        affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][mvNum].setMvField( mrgCtx.mvFieldNeighbours[(pos << 1) + 0].mv, mrgCtx.mvFieldNeighbours[(pos << 1) + 0].refIdx );
        affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][mvNum].setMvField( mrgCtx.mvFieldNeighbours[(pos << 1) + 1].mv, mrgCtx.mvFieldNeighbours[(pos << 1) + 1].refIdx );
      }
      affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] = mrgCtx.interDirNeighbours[pos];

      affMrgCtx.affineType[affMrgCtx.numValidMergeCand] = AFFINE_MODEL_NUM;
      affMrgCtx.mergeType[affMrgCtx.numValidMergeCand] = MRG_TYPE_SUBPU_ATMVP;
#if AFFINE_MMVD
      mrgCandIdx += (isAfMmvd && mrgCandIdx >= 0 ? 1 : 0);
#endif
      if ( affMrgCtx.numValidMergeCand == mrgCandIdx )
      {
        return;
      }

      affMrgCtx.numValidMergeCand++;

#if JVET_Z0139_NA_AFF && JVET_W0090_ARMC_TM
      if (isZeroCandIdx && !cs.pcv->isEncoder)
      {
        return;
      }
#endif
      // early termination
      if ( affMrgCtx.numValidMergeCand == maxNumAffineMergeCand )
      {
        return;
      }
    }
  }

  if ( slice.getSPS()->getUseAffine() )
  {
    ///> Start: inherited affine candidates
#if JVET_Z0139_HIST_AFF
    const int CHECKED_NEI_NUM = 7;
    const PredictionUnit* npu[CHECKED_NEI_NUM];
    const PredictionUnit* npuGroup2[CHECKED_NEI_NUM];
    Position posGroup2[CHECKED_NEI_NUM];
    int numGroup2;
    numGroup2 = 0;
    int neiIdx[CHECKED_NEI_NUM];
    int aiNeibeInherited[CHECKED_NEI_NUM];
    memset(aiNeibeInherited, 0, sizeof(aiNeibeInherited));

    const Position posLB = pu.Y().bottomLeft();
    const Position posLT = pu.Y().topLeft();
    const Position posRT = pu.Y().topRight();

    Position neiPositions[CHECKED_NEI_NUM] = { posLB.offset(-2, -1), posRT.offset(-1, -2), posRT.offset(3, -2), posLB.offset(-2, 3), posLT.offset(-2, -2),
    posLT.offset(-2, 1), posLT.offset(1, -2)};

    int numAffNeighLeft = getAvailableAffineNeighboursForLeftPredictor(pu, npu, neiIdx);
    if (numAffNeighLeft > 0)
    {
      aiNeibeInherited[neiIdx[0]] = 1;
    }
    int numAffNeigh = getAvailableAffineNeighboursForAbovePredictor(pu, npu, numAffNeighLeft, neiIdx);
    if (numAffNeigh > 0)
    {
      aiNeibeInherited[neiIdx[numAffNeigh - 1]] = 1;
    }
    for (int nei = 0; nei < CHECKED_NEI_NUM; nei++)
    {
      const PredictionUnit* puNei = pu.cs->getPURestricted(neiPositions[nei], pu, pu.chType);
      if (aiNeibeInherited[nei])
      {
        continue;
      }
      if (!puNei)
      {
        continue;
      }
      if (puNei->cu->predMode != MODE_INTER)
      {
        continue;
      }
      MotionInfo mvInfo = puNei->getMotionInfo(neiPositions[nei]);
      if (!mvInfo.isInter || mvInfo.interDir <= 0 || mvInfo.interDir > 3)
      {
        continue;
      }
      posGroup2[numGroup2] = neiPositions[nei];
      npuGroup2[numGroup2++] = puNei;
    }
#else
    const PredictionUnit* npu[5];
    int numAffNeighLeft = getAvailableAffineNeighboursForLeftPredictor( pu, npu );
    int numAffNeigh = getAvailableAffineNeighboursForAbovePredictor( pu, npu, numAffNeighLeft );
#endif
    for ( int idx = 0; idx < numAffNeigh; idx++ )
    {
      // derive Mv from Neigh affine PU
      Mv cMv[2][3];
      const PredictionUnit* puNeigh = npu[idx];
      pu.cu->affineType = puNeigh->cu->affineType;
      if ( puNeigh->interDir != 2 )
      {
        xInheritedAffineMv( pu, puNeigh, REF_PIC_LIST_0, cMv[0] );
      }
      if ( slice.isInterB() )
      {
        if ( puNeigh->interDir != 1 )
        {
          xInheritedAffineMv( pu, puNeigh, REF_PIC_LIST_1, cMv[1] );
        }
      }

      for ( int mvNum = 0; mvNum < 3; mvNum++ )
      {
        affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][mvNum].setMvField( cMv[0][mvNum], puNeigh->refIdx[0] );
        affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][mvNum].setMvField( cMv[1][mvNum], puNeigh->refIdx[1] );
      }
      affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] = puNeigh->interDir;
      affMrgCtx.affineType[affMrgCtx.numValidMergeCand] = (EAffineModel)(puNeigh->cu->affineType);
      affMrgCtx.BcwIdx[affMrgCtx.numValidMergeCand] = puNeigh->cu->BcwIdx;
#if INTER_LIC
      CHECK( puNeigh->interDir == 3 && puNeigh->cu->LICFlag, "LIC should not be enabled for affine bi-pred" );
      affMrgCtx.LICFlags[affMrgCtx.numValidMergeCand] = puNeigh->cu->LICFlag;
#endif
#if JVET_Z0139_HIST_AFF
      if (!checkLastAffineMergeCandRedundancy(pu, affMrgCtx))
      {
        continue;
      }
#endif
      if ( affMrgCtx.numValidMergeCand == mrgCandIdx )
      {
        return;
      }

      // early termination
      affMrgCtx.numValidMergeCand++;
      if ( affMrgCtx.numValidMergeCand == maxNumAffineMergeCand )
      {
        return;
      }
    }
#if JVET_Z0139_HIST_AFF
    MotionInfo tmvpInfo;
    tmvpInfo.interDir = 0;
    tmvpInfo.refIdx[0] = tmvpInfo.refIdx[1] = -1;
    bool       isTmvpAvailable = false;
#endif
    ///> End: inherited affine candidates

    ///> Start: Constructed affine candidates
    {
      MotionInfo mi[4];
      bool isAvailable[4] = { false };

      int8_t neighBcw[2] = { BCW_DEFAULT, BCW_DEFAULT };
      // control point: LT B2->B3->A2
      const Position posLT[3] = { pu.Y().topLeft().offset( -1, -1 ), pu.Y().topLeft().offset( 0, -1 ), pu.Y().topLeft().offset( -1, 0 ) };
      for ( int i = 0; i < 3; i++ )
      {
        const Position pos = posLT[i];
        const PredictionUnit* puNeigh = cs.getPURestricted( pos, pu, pu.chType );

#if JVET_Y0065_GPM_INTRA
        if (puNeigh && CU::isInter(*puNeigh->cu) && puNeigh->getMotionInfo( pos ).isInter && PU::isDiffMER(pu.lumaPos(), pos, plevel))
#else
        if (puNeigh && CU::isInter(*puNeigh->cu) && PU::isDiffMER(pu.lumaPos(), pos, plevel))
#endif
        {
          isAvailable[0] = true;
          mi[0] = puNeigh->getMotionInfo( pos );
          neighBcw[0] = puNeigh->cu->BcwIdx;
#if JVET_Z0139_HIST_AFF
          if (puNeigh->interDir != 3) neighBcw[0] = BCW_DEFAULT;
#endif
          break;
        }
      }

      // control point: RT B1->B0
      const Position posRT[2] = { pu.Y().topRight().offset( 0, -1 ), pu.Y().topRight().offset( 1, -1 ) };
      for ( int i = 0; i < 2; i++ )
      {
        const Position pos = posRT[i];
        const PredictionUnit* puNeigh = cs.getPURestricted( pos, pu, pu.chType );

#if JVET_Y0065_GPM_INTRA
        if (puNeigh && CU::isInter(*puNeigh->cu) && puNeigh->getMotionInfo( pos ).isInter && PU::isDiffMER(pu.lumaPos(), pos, plevel))
#else
        if (puNeigh && CU::isInter(*puNeigh->cu) && PU::isDiffMER(pu.lumaPos(), pos, plevel))
#endif
        {
          isAvailable[1] = true;
          mi[1] = puNeigh->getMotionInfo( pos );
          neighBcw[1] = puNeigh->cu->BcwIdx;
#if JVET_Z0139_HIST_AFF
          if (puNeigh->interDir != 3) neighBcw[1] = BCW_DEFAULT;
#endif
          break;
        }
      }

      // control point: LB A1->A0
      const Position posLB[2] = { pu.Y().bottomLeft().offset( -1, 0 ), pu.Y().bottomLeft().offset( -1, 1 ) };
      for ( int i = 0; i < 2; i++ )
      {
        const Position pos = posLB[i];
        const PredictionUnit* puNeigh = cs.getPURestricted( pos, pu, pu.chType );

#if JVET_Y0065_GPM_INTRA
        if (puNeigh && CU::isInter(*puNeigh->cu) && puNeigh->getMotionInfo( pos ).isInter && PU::isDiffMER(pu.lumaPos(), pos, plevel))
#else
        if (puNeigh && CU::isInter(*puNeigh->cu) && PU::isDiffMER(pu.lumaPos(), pos, plevel))
#endif
        {
          isAvailable[2] = true;
          mi[2] = puNeigh->getMotionInfo( pos );
          break;
        }
      }

      // control point: RB
      if ( slice.getPicHeader()->getEnableTMVPFlag() )
      {
        //>> MTK colocated-RightBottom
        // offset the pos to be sure to "point" to the same position the uiAbsPartIdx would've pointed to
        Position posRB = pu.Y().bottomRight().offset( -3, -3 );

        const PreCalcValues& pcv = *cs.pcv;
        Position posC0;
        bool C0Avail = false;

        bool boundaryCond = ((posRB.x + pcv.minCUWidth) < pcv.lumaWidth) && ((posRB.y + pcv.minCUHeight) < pcv.lumaHeight);
        const SubPic &curSubPic = pu.cs->slice->getPPS()->getSubPicFromPos(pu.lumaPos());
        if (curSubPic.getTreatedAsPicFlag())
        {
          boundaryCond = ((posRB.x + pcv.minCUWidth) <= curSubPic.getSubPicRight() &&
            (posRB.y + pcv.minCUHeight) <= curSubPic.getSubPicBottom());
        }
        if (boundaryCond)
        {
          int posYInCtu = posRB.y & pcv.maxCUHeightMask;
          if (posYInCtu + 4 < pcv.maxCUHeight)
          {
            posC0 = posRB.offset(4, 4);
            C0Avail = true;
          }
        }

        Mv        cColMv;
        int       refIdx = 0;
        bool      bExistMV = C0Avail && getColocatedMVP( pu, REF_PIC_LIST_0, posC0, cColMv, refIdx, false 
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
                                                       , &refIdx
#endif
        );
        if ( bExistMV )
        {
          mi[3].mv[0] = cColMv;
          mi[3].refIdx[0] = refIdx;
          mi[3].interDir = 1;
          isAvailable[3] = true;
#if JVET_Z0139_HIST_AFF
          tmvpInfo.mv[0] = cColMv;
          tmvpInfo.refIdx[0] = refIdx;
          tmvpInfo.interDir |= 1;
          isTmvpAvailable = true;
#endif
        }

        if ( slice.isInterB() )
        {
          bExistMV = C0Avail && getColocatedMVP( pu, REF_PIC_LIST_1, posC0, cColMv, refIdx, false
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
                                               , &refIdx
#endif
          );
          if ( bExistMV )
          {
            mi[3].mv[1] = cColMv;
            mi[3].refIdx[1] = refIdx;
            mi[3].interDir |= 2;
            isAvailable[3] = true;
#if JVET_Z0139_HIST_AFF
            tmvpInfo.mv[1] = cColMv;
            tmvpInfo.refIdx[1] = refIdx;
            tmvpInfo.interDir |= 2;
            isTmvpAvailable = true;
#endif
          }
        }
#if INTER_LIC
        mi[3].usesLIC = false;
#endif
      }

      //-------------------  insert model  -------------------//
      int order[6] = { 0, 1, 2, 3, 4, 5 };
      int modelNum = 6;
      int model[6][4] = {
        { 0, 1, 2 },          // 0:  LT, RT, LB
        { 0, 1, 3 },          // 1:  LT, RT, RB
        { 0, 2, 3 },          // 2:  LT, LB, RB
        { 1, 2, 3 },          // 3:  RT, LB, RB
        { 0, 1 },             // 4:  LT, RT
        { 0, 2 },             // 5:  LT, LB
      };

      int verNum[6] = { 3, 3, 3, 3, 2, 2 };
      int startIdx = pu.cs->sps->getUseAffineType() ? 0 : 4;
      for ( int idx = startIdx; idx < modelNum; idx++ )
      {
        int modelIdx = order[idx];
#if JVET_Z0139_HIST_AFF
        if (getAffineControlPointCand(pu, mi, isAvailable, model[modelIdx], ((modelIdx == 3) ? neighBcw[1] : neighBcw[0]), modelIdx, verNum[modelIdx], affMrgCtx))
        {
          if (affMrgCtx.numValidMergeCand == mrgCandIdx)
          {
            return;
          }
          affMrgCtx.numValidMergeCand ++;
          if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
          {
            return;
          }
        }
#else
        getAffineControlPointCand(pu, mi, isAvailable, model[modelIdx], ((modelIdx == 3) ? neighBcw[1] : neighBcw[0]), modelIdx, verNum[modelIdx], affMrgCtx);
        if ( affMrgCtx.numValidMergeCand != 0 && affMrgCtx.numValidMergeCand - 1 == mrgCandIdx )
        {
          return;
        }

        // early termination
        if ( affMrgCtx.numValidMergeCand == maxNumAffineMergeCand )
        {
          return;
        }
#endif
#if JVET_Z0139_HIST_AFF
        if (idx == startIdx)
        {
#if JVET_Z0118_GDR
          if (addSpatialAffineMergeHMVPCand(pu, affMrgCtx, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, 0, npuGroup2, posGroup2, numGroup2, mrgCandIdx))
#else
          if (addSpatialAffineMergeHMVPCand(pu, affMrgCtx, pu.cs->motionLut.lutAff, 0, npuGroup2, posGroup2, numGroup2, mrgCandIdx))
#endif
          {
            return;
          }

#if JVET_Z0118_GDR
          if (addOneInheritedHMVPAffineMergeCand(pu, affMrgCtx, (isClean) ? pu.cs->motionLut.lutAffInherit1 : pu.cs->motionLut.lutAffInherit0, 0))
#else
          if (addOneInheritedHMVPAffineMergeCand(pu, affMrgCtx, pu.cs->motionLut.lutAffInherit, 0))
#endif
          {
            if (affMrgCtx.numValidMergeCand == mrgCandIdx) // for decoder 
            {
              return;
            }

            affMrgCtx.numValidMergeCand++;

            // early termination
            if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
            {
              return;
            }
          }
        }
#endif
#if JVET_Z0139_NA_AFF
        if(!((pu.cu->slice->getPOC() - pu.cu->slice->getRefPOC( REF_PIC_LIST_0, 0)) == 1 && pu.cu->slice->getPicHeader()->getMvdL1ZeroFlag()))
        {
          if (idx == startIdx)
          {
            getNonAdjCstMergeCand(pu, affMrgCtx, mrgCandIdx, true);
            if (affMrgCtx.numValidMergeCand != 0 && affMrgCtx.numValidMergeCand - 1 == mrgCandIdx)
            {
              return;
            }

            if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
            {
              return;
            }
          }
        }
#endif
      }

#if JVET_Z0139_NA_AFF
      if((pu.cu->slice->getPOC() - pu.cu->slice->getRefPOC( REF_PIC_LIST_0, 0)) == 1 && pu.cu->slice->getPicHeader()->getMvdL1ZeroFlag())
      {
        getNonAdjCstMergeCand(pu, affMrgCtx, mrgCandIdx, true);
        if (affMrgCtx.numValidMergeCand != 0 && affMrgCtx.numValidMergeCand - 1 == mrgCandIdx)
        {
          return;
        }

        if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
        {
          return;
        }
      }
#endif
    }
    ///> End: Constructed affine candidates
#if JVET_Z0139_HIST_AFF 
    if (isTmvpAvailable)
    {
      Position posRB = pu.Y().bottomRight().offset(3, 3);
      if (addOneAffineMergeHMVPCand(pu, affMrgCtx
#if JVET_Z0118_GDR
                                 , (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0
#else
                                 , pu.cs->motionLut.lutAff
#endif
                                 , 0, tmvpInfo, posRB, BCW_DEFAULT
#if INTER_LIC
                                 , false
#endif
      ))
      {
        if (affMrgCtx.numValidMergeCand == mrgCandIdx) // for decoder 
        {
          return;
        }

        affMrgCtx.numValidMergeCand++;

        // early termination
        if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
        {
          return;
        }
      }
    }

    for (int iAffListIdx = 1; iAffListIdx < MAX_NUM_AFF_HMVP_CANDS; iAffListIdx++)
    {
#if JVET_Z0118_GDR
      if (addSpatialAffineMergeHMVPCand(pu, affMrgCtx, (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0, iAffListIdx, npuGroup2, posGroup2, numGroup2, mrgCandIdx))
#else
      if (addSpatialAffineMergeHMVPCand(pu, affMrgCtx, pu.cs->motionLut.lutAff, iAffListIdx, npuGroup2, posGroup2, numGroup2, mrgCandIdx))
#endif
      {
        return;
      }
      if (isTmvpAvailable)
      {
        Position posRB = pu.Y().bottomRight().offset(3, 3);
        if (addOneAffineMergeHMVPCand(pu, affMrgCtx
#if JVET_Z0118_GDR
                                    , (isClean) ? pu.cs->motionLut.lutAff1 : pu.cs->motionLut.lutAff0
#else
                                    , pu.cs->motionLut.lutAff
#endif
                                    , iAffListIdx, tmvpInfo, posRB, BCW_DEFAULT
#if INTER_LIC
                                    , false
#endif
        ))
        {
          if (affMrgCtx.numValidMergeCand == mrgCandIdx) // for decoder 
          {
            return;
          }

          affMrgCtx.numValidMergeCand++;

          // early termination
          if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
          {
            return;
          }
        }
      }
    }

    for (int iAffListIdx = 1; iAffListIdx < MAX_NUM_AFF_INHERIT_HMVP_CANDS; iAffListIdx++)
    {
#if JVET_Z0118_GDR
      if (addOneInheritedHMVPAffineMergeCand(pu, affMrgCtx, (isClean) ? pu.cs->motionLut.lutAffInherit1 : pu.cs->motionLut.lutAffInherit0, iAffListIdx))
#else
      if (addOneInheritedHMVPAffineMergeCand(pu, affMrgCtx, pu.cs->motionLut.lutAffInherit, iAffListIdx))
#endif
      {
        if (affMrgCtx.numValidMergeCand == mrgCandIdx) // for decoder 
        {
          return;
        }

        affMrgCtx.numValidMergeCand++;

        // early termination
        if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
        {
          return;
        }
      }
    }
    const int MAX_PAIRWISE_NUM = 9;
    const int preDefinedPairs[MAX_PAIRWISE_NUM][2] = { {0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 3}, { 2, 3}, { 0, 4}, {1, 4}, { 2, 4} };
    int iATMVPoffset = affMrgCtx.mergeType[0] == MRG_TYPE_SUBPU_ATMVP ? 1 : 0;
    int currSize = affMrgCtx.numValidMergeCand;
    for (int pairIdx = 0; pairIdx < MAX_PAIRWISE_NUM; pairIdx++)
    {
      int idx0 = preDefinedPairs[pairIdx][0] + iATMVPoffset;
      int idx1 = preDefinedPairs[pairIdx][1] + iATMVPoffset;

      CHECK(affMrgCtx.mergeType[idx0] == MRG_TYPE_SUBPU_ATMVP || affMrgCtx.mergeType[idx1] == MRG_TYPE_SUBPU_ATMVP, "Invalid Index");

      if (idx0 >= currSize || idx1 >= currSize)
      {
        break;
      }
      if ((affMrgCtx.interDirNeighbours[idx0] & affMrgCtx.interDirNeighbours[idx1]) == 0)
      {
        continue;
      }

      affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] = 0;

      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][0].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][1].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][2].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][0].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][1].setMvField(Mv(), -1);
      affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][2].setMvField(Mv(), -1);

      affMrgCtx.BcwIdx[affMrgCtx.numValidMergeCand] = BCW_DEFAULT;

#if INTER_LIC
      affMrgCtx.LICFlags[affMrgCtx.numValidMergeCand] = false;
#endif

      affMrgCtx.affineType[affMrgCtx.numValidMergeCand] = AFFINEMODEL_6PARAM;

      int numCPMVs = 3;

      if (affMrgCtx.mvFieldNeighbours[(idx0 << 1)][0].refIdx != -1 && affMrgCtx.mvFieldNeighbours[(idx1 << 1)][0].refIdx == affMrgCtx.mvFieldNeighbours[(idx0 << 1)][0].refIdx)
      {
        affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] |= 1;

        for (int mvIdx = 0; mvIdx < numCPMVs; mvIdx++)
        {
          Mv avgMv = affMrgCtx.mvFieldNeighbours[(idx0 << 1)][mvIdx].mv;
          avgMv += affMrgCtx.mvFieldNeighbours[(idx1 << 1)][mvIdx].mv;
          roundAffineMv(avgMv.hor, avgMv.ver, 1);
          affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 0][mvIdx].setMvField(avgMv, affMrgCtx.mvFieldNeighbours[(idx0 << 1)][0].refIdx);
        }
      }
      if (affMrgCtx.mvFieldNeighbours[(idx0 << 1) + 1][0].refIdx != -1 && affMrgCtx.mvFieldNeighbours[(idx1 << 1) + 1][0].refIdx == affMrgCtx.mvFieldNeighbours[(idx0 << 1) + 1][0].refIdx)
      {
        affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] |= 2;

        for (int mvIdx = 0; mvIdx < numCPMVs; mvIdx++)
        {
          Mv avgMv = affMrgCtx.mvFieldNeighbours[(idx0 << 1) + 1][mvIdx].mv;
          avgMv += affMrgCtx.mvFieldNeighbours[(idx1 << 1) + 1][mvIdx].mv;
          roundAffineMv(avgMv.hor, avgMv.ver, 1);
          affMrgCtx.mvFieldNeighbours[(affMrgCtx.numValidMergeCand << 1) + 1][mvIdx].setMvField(avgMv, affMrgCtx.mvFieldNeighbours[(idx0 << 1) + 1][0].refIdx);
        }
      }

      if (affMrgCtx.interDirNeighbours[affMrgCtx.numValidMergeCand] == 0)
      {
        continue;
      }
      if (checkLastAffineMergeCandRedundancy(pu, affMrgCtx))
      {
        if (affMrgCtx.numValidMergeCand == mrgCandIdx) // for decoder 
        {
          return;
        }

        affMrgCtx.numValidMergeCand++;

        // early termination
        if (affMrgCtx.numValidMergeCand == maxNumAffineMergeCand)
        {
          return;
        }
      }
    }
#endif
  }

  ///> zero padding
  int cnt = affMrgCtx.numValidMergeCand;
  while ( cnt < maxNumAffineMergeCand )
  {
    for ( int mvNum = 0; mvNum < 3; mvNum++ )
    {
      affMrgCtx.mvFieldNeighbours[(cnt << 1) + 0][mvNum].setMvField( Mv( 0, 0 ), 0 );
    }
    affMrgCtx.interDirNeighbours[cnt] = 1;

    if ( slice.isInterB() )
    {
      for ( int mvNum = 0; mvNum < 3; mvNum++ )
      {
        affMrgCtx.mvFieldNeighbours[(cnt << 1) + 1][mvNum].setMvField( Mv( 0, 0 ), 0 );
      }
      affMrgCtx.interDirNeighbours[cnt] = 3;
    }
    affMrgCtx.affineType[cnt] = AFFINEMODEL_4PARAM;
    if ( cnt == mrgCandIdx )
    {
      return;
    }
    cnt++;
    affMrgCtx.numValidMergeCand++;
  }
}

#if AFFINE_MMVD
uint8_t PU::getMergeIdxFromAfMmvdBaseIdx(AffineMergeCtx& affMrgCtx, uint16_t afMmvdBaseIdx)
{
  uint8_t mergeIdx = (uint8_t)afMmvdBaseIdx;
  for (int i = 0; i < affMrgCtx.numValidMergeCand; i++)
  {
    if (affMrgCtx.mergeType[i] != MRG_TYPE_DEFAULT_N)
    {
      ++mergeIdx;
    }
    else
    {
      break;
    }
  }

  return mergeIdx;
}

void PU::getAfMmvdMvf(const PredictionUnit& pu, const AffineMergeCtx& affineMergeCtx, MvField mvfMmvd[2][3], const uint16_t afMmvdBase, const uint16_t offsetStep, const uint16_t offsetDir)
{
  CHECK(!pu.cu->affine, "Affine flag is not on for Affine MMVD mode!");
  CHECK(affineMergeCtx.mergeType[afMmvdBase] != MRG_TYPE_DEFAULT_N, "AFF_MMVD base candidate type is not regular Affine!");

  static const int32_t refMvdCands[AF_MMVD_STEP_NUM] = { 1, 2, 4, 8, 16 };
#if JVET_Y0128_NON_CTC
  const int32_t iPicSize = pu.cu->slice->getPic()->lumaSize().area();
  const int32_t mvShift  = iPicSize < 921600 ? 0 : ( iPicSize < 4096000 ? 2 : MV_FRACTIONAL_BITS_INTERNAL - 1); // 921600 = 1280x720, 4096000 = 2560x1600
#else
  static const int32_t iPicSize = pu.cu->slice->getPic()->lumaSize().area();
  static const int32_t mvShift  = iPicSize < 921600 ? 0 : ( iPicSize < 4096000 ? 2 : MV_FRACTIONAL_BITS_INTERNAL - 1); // 921600 = 1280x720, 4096000 = 2560x1600
#endif
  int step = refMvdCands[offsetStep] << mvShift;
  int affineType = affineMergeCtx.affineType[afMmvdBase];

  uint8_t interDir = affineMergeCtx.interDirNeighbours[afMmvdBase];
  int8_t refIdxL0 = affineMergeCtx.mvFieldNeighbours[(afMmvdBase << 1)    ][0].refIdx;
  int8_t refIdxL1 = affineMergeCtx.mvFieldNeighbours[(afMmvdBase << 1) + 1][0].refIdx;

  Mv baseMv[2][3];
  for (int i = 0; i < 3; i++)
  {
    baseMv[0][i] = affineMergeCtx.mvFieldNeighbours[(afMmvdBase << 1)    ][i].mv;
    baseMv[1][i] = affineMergeCtx.mvFieldNeighbours[(afMmvdBase << 1) + 1][i].mv;
    mvfMmvd[0][i].refIdx = refIdxL0;
    mvfMmvd[1][i].refIdx = refIdxL1;
    mvfMmvd[0][i].mv = Mv();
    mvfMmvd[1][i].mv = Mv();
  }

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  const int xDir[] = {1, -1,  0,  0,  1, -1,  1, -1, 2, -2,  2, -2, 1,  1, -1, -1};
  const int yDir[] = {0,  0,  1, -1,  1, -1, -1,  1, 1,  1, -1, -1, 2, -2,  2, -2};
  int offsetX = xDir[offsetDir]  * step;
  int offsetY = yDir[offsetDir] * step;
#else
  int magY = (offsetDir >> 1) & 0x1;
  int sign = (offsetDir & 0x1) ? -1 : 1;
  int offsetX = (1 - magY) * sign * step;
  int offsetY = (    magY) * sign * step;
#endif
  Mv offsetMv(offsetX, offsetY);

  int numCp = (affineType == AFFINEMODEL_4PARAM) ? 2 : 3;
  for (int cpIdx = 0; cpIdx < numCp; cpIdx++)
  {
    if (interDir == 1)
    {
      mvfMmvd[0][cpIdx].mv = baseMv[0][cpIdx] + offsetMv;
    }
    else if (interDir == 2)
    {
      mvfMmvd[1][cpIdx].mv = baseMv[1][cpIdx] + offsetMv;
    }
    else if (interDir == 3)
    {
      int poc_cur = pu.cu->slice->getPOC();
      int poc_l0  = pu.cu->slice->getRefPOC(REF_PIC_LIST_0, refIdxL0);
      int poc_l1  = pu.cu->slice->getRefPOC(REF_PIC_LIST_1, refIdxL1);

      int distL0 = poc_l0 - poc_cur;
      int distL1 = poc_l1 - poc_cur;
      mvfMmvd[0][cpIdx].mv =                                                     baseMv[0][cpIdx] + offsetMv;
      mvfMmvd[1][cpIdx].mv = distL0 * distL1 < 0 ? baseMv[1][cpIdx] - offsetMv : baseMv[1][cpIdx] + offsetMv;
    }
  }
}

int32_t PU::getAfMmvdEstBits(const PredictionUnit &pu)
{
  int baseBits = (AF_MMVD_BASE_NUM == 1 ? 0 : pu.afMmvdBaseIdx + (pu.afMmvdBaseIdx == AF_MMVD_BASE_NUM - 1 ? 0: 1));
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  int stepBits = pu.afMmvdMergeIdx + (pu.afMmvdMergeIdx == AF_MMVD_MAX_REFINE_NUM - 1 ? 0 : 1);
  return stepBits + baseBits;
#else
  int stepBits = pu.afMmvdStep + (pu.afMmvdStep == AF_MMVD_STEP_NUM - 1 ? 0 : 1);
  int dirBits  = gp_sizeIdxInfo->idxFrom(AF_MMVD_OFFSET_DIR);
  
  return stepBits + dirBits + baseBits;
#endif
}
#endif

void PU::setAllAffineMvField( PredictionUnit &pu, MvField *mvField, RefPicList eRefList )
{
  // Set Mv
  Mv mv[3];
  for ( int i = 0; i < 3; i++ )
  {
    mv[i] = mvField[i].mv;
  }
  setAllAffineMv( pu, mv[0], mv[1], mv[2], eRefList );

  // Set RefIdx
  CHECK( mvField[0].refIdx != mvField[1].refIdx || mvField[0].refIdx != mvField[2].refIdx, "Affine mv corners don't have the same refIdx." );
  pu.refIdx[eRefList] = mvField[0].refIdx;
}

void PU::setAllAffineMv(PredictionUnit& pu, Mv affLT, Mv affRT, Mv affLB, RefPicList eRefList, bool clipCPMVs)
{
  int width  = pu.Y().width;
  int shift = MAX_CU_DEPTH;
  const bool isTranslational = (affLT == affRT && affLT == affLB);

  if (clipCPMVs)
  {
    affLT.mvCliptoStorageBitDepth();
    affRT.mvCliptoStorageBitDepth();
    if (pu.cu->affineType == AFFINEMODEL_6PARAM)
    {
      affLB.mvCliptoStorageBitDepth();
    }
  }
  int deltaMvHorX, deltaMvHorY, deltaMvVerX, deltaMvVerY;
  deltaMvHorX = (affRT - affLT).getHor() << (shift - floorLog2(width));
  deltaMvHorY = (affRT - affLT).getVer() << (shift - floorLog2(width));
  int height = pu.Y().height;
  if ( pu.cu->affineType == AFFINEMODEL_6PARAM )
  {
    deltaMvVerX = (affLB - affLT).getHor() << (shift - floorLog2(height));
    deltaMvVerY = (affLB - affLT).getVer() << (shift - floorLog2(height));
  }
  else
  {
    deltaMvVerX = -deltaMvHorY;
    deltaMvVerY = deltaMvHorX;
  }

  int mvScaleHor = affLT.getHor() << shift;
  int mvScaleVer = affLT.getVer() << shift;

  int blockWidth = AFFINE_MIN_BLOCK_SIZE;
  int blockHeight = AFFINE_MIN_BLOCK_SIZE;
  const int halfBW = blockWidth >> 1;
  const int halfBH = blockHeight >> 1;

  MotionBuf mb = pu.getMotionBuf();
  int mvScaleTmpHor, mvScaleTmpVer;
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  const bool subblkMVSpreadOverLimit = InterPrediction::isSubblockVectorSpreadOverLimit( deltaMvHorX, deltaMvHorY, deltaMvVerX, deltaMvVerY, pu.interDir );
#endif

  if (isTranslational)
  {
    mvScaleTmpHor = mvScaleHor;
    mvScaleTmpVer = mvScaleVer;
    roundAffineMv(mvScaleTmpHor, mvScaleTmpVer, shift);
    Mv curMv(mvScaleTmpHor, mvScaleTmpVer);
    curMv.clipToStorageBitDepth();
    MotionInfo *mi = mb.buf;

    for (int y = 0; y < (pu.Y().height >> MIN_CU_LOG2); y++)
    {
      for (int x = 0; x < (pu.Y().width >> MIN_CU_LOG2); x++)
      {
        mi[x].mv[eRefList] = curMv;
      }
      mi += mb.stride;
    }
  }
#if !AFFINE_RM_CONSTRAINTS_AND_OPT
  else if (subblkMVSpreadOverLimit)
  {
      mvScaleTmpHor = mvScaleHor + deltaMvHorX * ( pu.Y().width >> 1 ) + deltaMvVerX * ( pu.Y().height >> 1 );
      mvScaleTmpVer = mvScaleVer + deltaMvHorY * ( pu.Y().width >> 1 ) + deltaMvVerY * ( pu.Y().height >> 1 );
      roundAffineMv(mvScaleTmpHor, mvScaleTmpVer, shift);
      Mv curMv(mvScaleTmpHor, mvScaleTmpVer);
      curMv.clipToStorageBitDepth();
      MotionInfo *mi = mb.buf;
      for (int y = 0; y < (pu.Y().height >> MIN_CU_LOG2); y++)
      {
        for (int x = 0; x < (pu.Y().width >> MIN_CU_LOG2); x++)
        {
          mi[x].mv[eRefList] = curMv;
        }
        mi += mb.stride;
      }
  }
#endif
  else
  {
    int mvScaleHorLine = mvScaleHor + deltaMvHorX * halfBW + deltaMvVerX * halfBH;
    int mvScaleVerLine = mvScaleVer + deltaMvHorY * halfBW + deltaMvVerY * halfBH;
    int deltaMvVerXLine = deltaMvVerX * blockHeight;
    int deltaMvVerYLine = deltaMvVerY * blockHeight;
    int deltaMvHorXBlk = deltaMvHorX * blockWidth;
    int deltaMvHorYBlk = deltaMvHorY * blockWidth;

    MotionInfo *mi = mb.buf;

    for( int h = 0; h < pu.Y().height; h += blockHeight )
    {
      int mvScaleHorBlk = mvScaleHorLine;
      int mvScaleVerBlk = mvScaleVerLine;

      for( int w = 0; w < pu.Y().width; w += blockWidth )
      {
        mvScaleTmpHor = mvScaleHorBlk;
        mvScaleTmpVer = mvScaleVerBlk;

        mvScaleHorBlk += deltaMvHorXBlk;
        mvScaleVerBlk += deltaMvHorYBlk;

        roundAffineMv( mvScaleTmpHor, mvScaleTmpVer, shift );
        Mv curMv( mvScaleTmpHor, mvScaleTmpVer );
        curMv.clipToStorageBitDepth();

        mi[w >> MIN_CU_LOG2].mv[eRefList] = curMv;
      }

      mi += ( blockHeight >> MIN_CU_LOG2 ) * mb.stride;
      mvScaleHorLine += deltaMvVerXLine;
      mvScaleVerLine += deltaMvVerYLine;
    }
  }

  pu.mvAffi[eRefList][0] = affLT;
  pu.mvAffi[eRefList][1] = affRT;
  pu.mvAffi[eRefList][2] = affLB;
}

void clipColPos(int& posX, int& posY, const PredictionUnit& pu)
{
  Position puPos = pu.lumaPos();
  int log2CtuSize = floorLog2(pu.cs->sps->getCTUSize());
  int ctuX = ((puPos.x >> log2CtuSize) << log2CtuSize);
  int ctuY = ((puPos.y >> log2CtuSize) << log2CtuSize);
  int horMax;
  const SubPic &curSubPic = pu.cu->slice->getPPS()->getSubPicFromPos(puPos);
  if (curSubPic.getTreatedAsPicFlag())
  {
    horMax = std::min((int)curSubPic.getSubPicRight(), ctuX + (int)pu.cs->sps->getCTUSize() + 3);
  }
  else
  {
    horMax = std::min((int)pu.cs->pps->getPicWidthInLumaSamples() - 1, ctuX + (int)pu.cs->sps->getCTUSize() + 3);
  }
  int horMin = std::max((int)0, ctuX);
  int verMax = std::min( (int)pu.cs->pps->getPicHeightInLumaSamples() - 1, ctuY + (int)pu.cs->sps->getCTUSize() - 1 );
  int verMin = std::max((int)0, ctuY);

  posX = std::min(horMax, std::max(horMin, posX));
  posY = std::min(verMax, std::max(verMin, posY));
}

bool PU::getInterMergeSubPuMvpCand(const PredictionUnit &pu, MergeCtx& mrgCtx, bool& LICFlag, const int count
  , int mmvdList
)
{
  const Slice   &slice = *pu.cs->slice;
  const unsigned scale = 4 * std::max<int>(1, 4 * AMVP_DECIMATION_FACTOR / 4);
  const unsigned mask = ~(scale - 1);
#if INTER_LIC
  LICFlag = false;
#endif

  const Picture *pColPic = slice.getRefPic(RefPicList(slice.isInterB() ? 1 - slice.getColFromL0Flag() : 0), slice.getColRefIdx());
  Mv cTMv;

  if ( count )
  {
    if ( (mrgCtx.interDirNeighbours[0] & (1 << REF_PIC_LIST_0)) && slice.getRefPic( REF_PIC_LIST_0, mrgCtx.mvFieldNeighbours[REF_PIC_LIST_0].refIdx ) == pColPic )
    {
      cTMv = mrgCtx.mvFieldNeighbours[REF_PIC_LIST_0].mv;
    }
    else if ( slice.isInterB() && (mrgCtx.interDirNeighbours[0] & (1 << REF_PIC_LIST_1)) && slice.getRefPic( REF_PIC_LIST_1, mrgCtx.mvFieldNeighbours[REF_PIC_LIST_1].refIdx ) == pColPic )
    {
      cTMv = mrgCtx.mvFieldNeighbours[REF_PIC_LIST_1].mv;
    }
  }

  ///////////////////////////////////////////////////////////////////////
  ////////          GET Initial Temporal Vector                  ////////
  ///////////////////////////////////////////////////////////////////////

  Mv cTempVector = cTMv;
  bool  tempLICFlag = false;

  // compute the location of the current PU
  Position puPos = pu.lumaPos();
  Size puSize = pu.lumaSize();
  int numPartLine = std::max(puSize.width >> ATMVP_SUB_BLOCK_SIZE, 1u);
  int numPartCol = std::max(puSize.height >> ATMVP_SUB_BLOCK_SIZE, 1u);
  int puHeight = numPartCol == 1 ? puSize.height : 1 << ATMVP_SUB_BLOCK_SIZE;
  int puWidth = numPartLine == 1 ? puSize.width : 1 << ATMVP_SUB_BLOCK_SIZE;

  Mv cColMv;
  int refIdx = 0;
  // use coldir.
  bool     bBSlice = slice.isInterB();

  Position centerPos;

  bool found = false;
  cTempVector = cTMv;

  cTempVector.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
  int tempX = cTempVector.getHor();
  int tempY = cTempVector.getVer();

  centerPos.x = puPos.x + (puSize.width >> 1) + tempX;
  centerPos.y = puPos.y + (puSize.height >> 1) + tempY;

  clipColPos(centerPos.x, centerPos.y, pu);

  centerPos = Position{ PosType(centerPos.x & mask), PosType(centerPos.y & mask) };

  // derivation of center motion parameters from the collocated CU
  const MotionInfo &mi = pColPic->cs->getMotionInfo(centerPos);

  if (mi.isInter && mi.isIBCmot == false)
  {
    mrgCtx.interDirNeighbours[count] = 0;

    for (unsigned currRefListId = 0; currRefListId < (bBSlice ? 2 : 1); currRefListId++)
    {
      RefPicList  currRefPicList = RefPicList(currRefListId);

      if (getColocatedMVP(pu, currRefPicList, centerPos, cColMv, refIdx, true
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
                        , &refIdx
#endif
      ))
      {
        // set as default, for further motion vector field spanning
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
        mrgCtx.mvFieldNeighbours[(count << 1) + currRefListId].setMvField(cColMv, refIdx);
#else
        mrgCtx.mvFieldNeighbours[(count << 1) + currRefListId].setMvField(cColMv, 0);
#endif
        mrgCtx.interDirNeighbours[count] |= (1 << currRefListId);
        LICFlag = tempLICFlag;
        mrgCtx.BcwIdx[count] = BCW_DEFAULT;
        found = true;
      }
      else
      {
        mrgCtx.mvFieldNeighbours[(count << 1) + currRefListId].setMvField(Mv(), NOT_VALID);
        mrgCtx.interDirNeighbours[count] &= ~(1 << currRefListId);
      }
    }
  }

  if (!found)
  {
    return false;
  }
  if (mmvdList != 1)
  {
    int xOff = (puWidth >> 1) + tempX;
    int yOff = (puHeight >> 1) + tempY;

    MotionBuf &mb = mrgCtx.subPuMvpMiBuf;

    const bool isBiPred = isBipredRestriction(pu);

    for (int y = puPos.y; y < puPos.y + puSize.height; y += puHeight)
    {
      for (int x = puPos.x; x < puPos.x + puSize.width; x += puWidth)
      {
        Position colPos{ x + xOff, y + yOff };

        clipColPos(colPos.x, colPos.y, pu);

        colPos = Position{ PosType(colPos.x & mask), PosType(colPos.y & mask) };

        const MotionInfo &colMi = pColPic->cs->getMotionInfo(colPos);

        MotionInfo mi;

        found       = false;
        mi.isInter  = true;
        mi.sliceIdx = slice.getIndependentSliceIdx();
        mi.isIBCmot = false;
        if (colMi.isInter && colMi.isIBCmot == false)
        {
          for (unsigned currRefListId = 0; currRefListId < (bBSlice ? 2 : 1); currRefListId++)
          {
            RefPicList currRefPicList = RefPicList(currRefListId);
            if (getColocatedMVP(pu, currRefPicList, colPos, cColMv, refIdx, true
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
              , &refIdx
#endif
            ))
            {
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING
              mi.refIdx[currRefListId] = refIdx;
#else
              mi.refIdx[currRefListId] = 0;
#endif
              mi.mv[currRefListId]     = cColMv;
              found                    = true;
            }
          }
        }
        if (!found)
        {
          mi.mv[0]     = mrgCtx.mvFieldNeighbours[(count << 1) + 0].mv;
          mi.mv[1]     = mrgCtx.mvFieldNeighbours[(count << 1) + 1].mv;
          mi.refIdx[0] = mrgCtx.mvFieldNeighbours[(count << 1) + 0].refIdx;
          mi.refIdx[1] = mrgCtx.mvFieldNeighbours[(count << 1) + 1].refIdx;
        }

        mi.interDir = (mi.refIdx[0] != -1 ? 1 : 0) + (mi.refIdx[1] != -1 ? 2 : 0);

        if (isBiPred && mi.interDir == 3)
        {
          mi.interDir  = 1;
          mi.mv[1]     = Mv();
          mi.refIdx[1] = NOT_VALID;
        }

        mb.subBuf(g_miScaling.scale(Position{ x, y } - pu.lumaPos()), g_miScaling.scale(Size(puWidth, puHeight)))
          .fill(mi);
      }
    }
  }
  return true;
}

#if MULTI_PASS_DMVR
void PU::spanMotionInfo( PredictionUnit &pu, const MergeCtx &mrgCtx, Mv* bdmvrSubPuMv0, Mv* bdmvrSubPuMv1, Mv* bdofSubPuMvOffset)
#else
void PU::spanMotionInfo( PredictionUnit &pu, const MergeCtx &mrgCtx )
#endif
{
#if !MULTI_PASS_DMVR
  MotionBuf mb = pu.getMotionBuf();
#endif
#if JVET_W0123_TIMD_FUSION
  IpmBuf ib = pu.getIpmBuf();
#endif

  if (!pu.mergeFlag || pu.mergeType == MRG_TYPE_DEFAULT_N || pu.mergeType == MRG_TYPE_IBC)
  {
    MotionInfo mi;

    mi.isInter = !CU::isIntra(*pu.cu);
    mi.isIBCmot = CU::isIBC(*pu.cu);
    mi.sliceIdx = pu.cu->slice->getIndependentSliceIdx();
#if INTER_LIC
    mi.usesLIC = pu.cu->LICFlag;
#endif

    if( mi.isInter )
    {
      mi.interDir = pu.interDir;
      mi.useAltHpelIf = pu.cu->imv == IMV_HPEL;
      for( int i = 0; i < NUM_REF_PIC_LIST_01; i++ )
      {
        mi.mv[i]     = pu.mv[i];
        mi.refIdx[i] = pu.refIdx[i];
      }
      if (mi.isIBCmot)
      {
        mi.bv = pu.bv;
      }
    }

#if MULTI_PASS_DMVR
    if (pu.bdmvrRefine)
    {
      CHECK(bdmvrSubPuMv0 == nullptr, "this is not possible");
      const int dx = std::min<int>(pu.lwidth (), BDOF_SUBPU_DIM);
      const int dy = std::min<int>(pu.lheight(), BDOF_SUBPU_DIM);
      int subPuIdx = 0;
      const int bioSubPuIdxStrideIncr = BDOF_SUBPU_STRIDE - std::max(1, (int)(pu.lwidth() >> BDOF_SUBPU_DIM_LOG2));

      for (int yStart = 0; yStart < pu.lheight(); yStart += dy)
      {
        for (int xStart = 0; xStart < pu.lwidth(); xStart += dx)
        {
          const int bdmvrSubPuIdx = (yStart >> DMVR_SUBCU_HEIGHT_LOG2) * DMVR_SUBPU_STRIDE + (xStart >> DMVR_SUBCU_WIDTH_LOG2);
#if JVET_Y0089_DMVR_BCW
          if (pu.cu->BcwIdx != BCW_DEFAULT)
          {
            mi.mv[0] = bdmvrSubPuMv0[bdmvrSubPuIdx];
            mi.mv[1] = bdmvrSubPuMv1[bdmvrSubPuIdx];
          }
          else
#endif
          {
            mi.mv[0] = bdmvrSubPuMv0[bdmvrSubPuIdx] + bdofSubPuMvOffset[subPuIdx];
            mi.mv[1] = bdmvrSubPuMv1[bdmvrSubPuIdx] - bdofSubPuMvOffset[subPuIdx];
          }

          subPuIdx++;
          MotionBuf mb = pu.cs->getMotionBuf(Area(pu.lx() + xStart, pu.ly() + yStart, dx, dy));
          mb.fill(mi);
        }
        subPuIdx += bioSubPuIdxStrideIncr;
      }
#if JVET_W0123_TIMD_FUSION
      MotionBuf mb = pu.getMotionBuf();
      spanIpmInfoInter(pu, mb, ib);
#endif
      return;
    }
    MotionBuf mb = pu.getMotionBuf();
#endif
    if (pu.cu->affine)
    {
      mi.mv[0].setZero(); // to make sure filling of MV in unused reference list
      mi.mv[1].setZero();
      mb.fill(mi);
      if (pu.refIdx[0] >= 0)
      {
        PU::setAllAffineMv(pu, pu.mvAffi[0][0], pu.mvAffi[0][1], pu.mvAffi[0][2], REF_PIC_LIST_0);
      }
      if (pu.refIdx[1] >= 0)
      {
        PU::setAllAffineMv(pu, pu.mvAffi[1][0], pu.mvAffi[1][1], pu.mvAffi[1][2], REF_PIC_LIST_1);
      }
#if JVET_W0123_TIMD_FUSION
      spanIpmInfoInter(pu, mb, ib);
#endif
    }
    else
    {
      mb.fill(mi);
#if JVET_W0123_TIMD_FUSION
      if (mi.isIBCmot
#if JVET_X0141_CIIP_TIMD_TM
        || pu.ciipFlag
#endif
        )
      {
        ib.fill(PLANAR_IDX);
      }
      else
      {
        spanIpmInfoInter(pu, mb, ib);
      }
#endif
    }
  }
  else if (pu.mergeType == MRG_TYPE_SUBPU_ATMVP)
  {
    CHECK(mrgCtx.subPuMvpMiBuf.area() == 0 || !mrgCtx.subPuMvpMiBuf.buf, "Buffer not initialized");
#if MULTI_PASS_DMVR
    MotionBuf mb = pu.getMotionBuf();
#endif
    mb.copyFrom(mrgCtx.subPuMvpMiBuf);
#if JVET_W0123_TIMD_FUSION
    spanIpmInfoInter(pu, mb, ib);
#endif
  }
}

#if JVET_W0123_TIMD_FUSION
#if MULTI_PASS_DMVR
void PU::spanMotionInfo2( PredictionUnit &pu, const MergeCtx &mrgCtx, Mv* bdmvrSubPuMv0, Mv* bdmvrSubPuMv1, Mv* bdofSubPuMvOffset)
#else
void PU::spanMotionInfo2( PredictionUnit &pu, const MergeCtx &mrgCtx )
#endif
{
#if !MULTI_PASS_DMVR
  MotionBuf mb = pu.getMotionBuf();
#endif

  if (!pu.mergeFlag || pu.mergeType == MRG_TYPE_DEFAULT_N || pu.mergeType == MRG_TYPE_IBC)
  {
    MotionInfo mi;

    mi.isInter = !CU::isIntra(*pu.cu);
    mi.isIBCmot = CU::isIBC(*pu.cu);
    mi.sliceIdx = pu.cu->slice->getIndependentSliceIdx();
#if INTER_LIC
    mi.usesLIC = pu.cu->LICFlag;
#endif

    if( mi.isInter )
    {
      mi.interDir = pu.interDir;
      mi.useAltHpelIf = pu.cu->imv == IMV_HPEL;
      for( int i = 0; i < NUM_REF_PIC_LIST_01; i++ )
      {
        mi.mv[i]     = pu.mv[i];
        mi.refIdx[i] = pu.refIdx[i];
      }
      if (mi.isIBCmot)
      {
        mi.bv = pu.bv;
      }
    }

#if MULTI_PASS_DMVR
    if (pu.bdmvrRefine)
    {
      CHECK(bdmvrSubPuMv0 == nullptr, "this is not possible");
      const int dx = std::min<int>(pu.lwidth (), BDOF_SUBPU_DIM);
      const int dy = std::min<int>(pu.lheight(), BDOF_SUBPU_DIM);
      int subPuIdx = 0;
      const int bioSubPuIdxStrideIncr = BDOF_SUBPU_STRIDE - std::max(1, (int)(pu.lwidth() >> BDOF_SUBPU_DIM_LOG2));

      for (int yStart = 0; yStart < pu.lheight(); yStart += dy)
      {
        for (int xStart = 0; xStart < pu.lwidth(); xStart += dx)
        {
          const int bdmvrSubPuIdx = (yStart >> DMVR_SUBCU_HEIGHT_LOG2) * DMVR_SUBPU_STRIDE + (xStart >> DMVR_SUBCU_WIDTH_LOG2);
          mi.mv[0] = bdmvrSubPuMv0[bdmvrSubPuIdx] + bdofSubPuMvOffset[subPuIdx];
          mi.mv[1] = bdmvrSubPuMv1[bdmvrSubPuIdx] - bdofSubPuMvOffset[subPuIdx];

          subPuIdx++;
          MotionBuf mb = pu.cs->getMotionBuf(Area(pu.lx() + xStart, pu.ly() + yStart, dx, dy));
          mb.fill(mi);
        }
        subPuIdx += bioSubPuIdxStrideIncr;
      }
      return;
    }
    MotionBuf mb = pu.getMotionBuf();
#endif
    if (pu.cu->affine)
    {
      mi.mv[0].setZero(); // to make sure filling of MV in unused reference list
      mi.mv[1].setZero();
      mb.fill(mi);
      if (pu.refIdx[0] >= 0)
      {
        PU::setAllAffineMv(pu, pu.mvAffi[0][0], pu.mvAffi[0][1], pu.mvAffi[0][2], REF_PIC_LIST_0);
      }
      if (pu.refIdx[1] >= 0)
      {
        PU::setAllAffineMv(pu, pu.mvAffi[1][0], pu.mvAffi[1][1], pu.mvAffi[1][2], REF_PIC_LIST_1);
      }
    }
    else
    {
      mb.fill(mi);
    }
  }
  else if (pu.mergeType == MRG_TYPE_SUBPU_ATMVP)
  {
    CHECK(mrgCtx.subPuMvpMiBuf.area() == 0 || !mrgCtx.subPuMvpMiBuf.buf, "Buffer not initialized");
#if MULTI_PASS_DMVR
    MotionBuf mb = pu.getMotionBuf();
#endif
    mb.copyFrom(mrgCtx.subPuMvpMiBuf);
  }
}

void PU::spanIpmInfoIntra( PredictionUnit &pu)
{
  int ipm = PU::getIntraDirLuma(pu);
  if (pu.cu->timd)
  {
    ipm = MAP131TO67(ipm);
  }
  IpmBuf ib = pu.getIpmBuf();
  ib.fill(ipm);
}

#if RPR_ENABLE
void scalePositionInRef( PredictionUnit& pu, const PPS& pps, RefPicList refList, int refIdx, Position& PosY )
{
  const Picture* refPic = pu.cu->slice->getRefPic( refList, refIdx )->unscaledPic;
  const bool scaled = refPic->isRefScaled( &pps );
  if (scaled)
  {
    const SPS* sps = pu.cu->cs->sps;
    ChromaFormat chFmt = sps->getChromaFormatIdc();

    const std::pair<int, int> scalingRatio = pu.cu->slice->getScalingRatio( refList, refIdx );

    int64_t x0Int;
    int64_t y0Int;

    int posX = PosY.x;
    int posY = PosY.y;

    const int posShift = SCALE_RATIO_BITS - 4;

    int shiftHor = MV_FRACTIONAL_BITS_INTERNAL ;
    int shiftVer = MV_FRACTIONAL_BITS_INTERNAL ;
    int offX = 1 << (posShift - shiftHor - 1);
    int offY = 1 << (posShift - shiftVer - 1);

    x0Int = (posX << (4)) * (int64_t)scalingRatio.first ;
    x0Int = SIGN(x0Int) * ((llabs(x0Int) + ((long long)1 << 7)) >> 8) + ((refPic->getScalingWindow().getWindowLeftOffset() * SPS::getWinUnitX(chFmt)) << (posShift));

    y0Int = (posY << (4)) * (int64_t)scalingRatio.second ;
    y0Int = SIGN(y0Int) * ((llabs(y0Int) + ((long long)1 << 7)) >> 8) + ((refPic->getScalingWindow().getWindowTopOffset() * SPS::getWinUnitY(chFmt)) << (posShift));

    int xInt0 = ((int32_t)x0Int + offX) >> posShift;
    int yInt0 = ((int32_t)y0Int + offY) >> posShift;

    PosY.x = std::min( std::max(xInt0, 0), (int)refPic->unscaledPic->getPicWidthInLumaSamples()  - 1 );
    PosY.y = std::min( std::max(yInt0, 0), (int)refPic->unscaledPic->getPicHeightInLumaSamples() - 1 );
  }
  //else
  //  clipColPos( PosY.x, PosY.y, pu );
}
#endif

void PU::spanIpmInfoInter( PredictionUnit &pu, MotionBuf &mb, IpmBuf &ib)
{
  const unsigned scale = 4 * std::max<int>(1, 4 * AMVP_DECIMATION_FACTOR / 4);
  const unsigned mask = ~(scale - 1);
  Mv cMv;
  RefPicList refList;
  int refIdx;
  Position PosY;
  MotionInfo tempMi;
  MotionInfo mi0;
  MotionInfo mi1;
  Position PosY0;
  Position PosY1;
  Mv cMv0;
  Mv cMv1;
#if RPR_ENABLE
  const Picture* pRefPic0;
  const Picture* pRefPic1;
#else
  Picture* pRefPic0;
  Picture* pRefPic1;
#endif
  uint8_t* ii = ib.buf;
  int ibH = mb.height;
  int ibW = mb.width;
  for (int y = 0; y < ibH; y++)
  {
    for (int x = 0; x < ibW; x++)
    {
      uint8_t ipm = PLANAR_IDX;
      tempMi = mb.at(x, y);
      if (tempMi.interDir != 3)
      {
        if (tempMi.interDir != 2)
        {
          cMv = tempMi.mv[0];
          refList = REF_PIC_LIST_0;
          refIdx = tempMi.refIdx[0];
        }
        else
        {
          cMv = tempMi.mv[1];
          refList = REF_PIC_LIST_1;
          refIdx = tempMi.refIdx[1];
        }
        if (refList < 0 || refIdx < 0)
        {
          ipm = PLANAR_IDX;
        }
        else
        {
          cMv.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
          PosY.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv.getHor();
          PosY.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv.getVer();
          clipColPos(PosY.x, PosY.y, pu);
#if RPR_ENABLE
          scalePositionInRef( pu, *pu.cs->pps, refList, refIdx, PosY );
          PosY.x = (PosY.x & mask);
          PosY.y = (PosY.y & mask);
          ipm = pu.cu->slice->getRefPic(refList, refIdx)->unscaledPic->cs->getIpmInfo(PosY);
#else
          PosY.x = (PosY.x & mask);
          PosY.y = (PosY.y & mask);
          ipm = pu.cu->slice->getRefPic(refList, refIdx)->cs->getIpmInfo(PosY);
#endif
        }
      }
      else
      {
#if RPR_ENABLE
        pRefPic0 = pu.cu->slice->getRefPic(REF_PIC_LIST_0, tempMi.refIdx[0])->unscaledPic;
#else
        pRefPic0 = pu.cu->slice->getRefPic(REF_PIC_LIST_0, tempMi.refIdx[0]);
#endif
        cMv0 = tempMi.mv[0];
        cMv0.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
        PosY0.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv0.getHor();
        PosY0.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv0.getVer();
        clipColPos(PosY0.x, PosY0.y, pu);
#if RPR_ENABLE
        scalePositionInRef( pu, *pu.cs->pps, REF_PIC_LIST_0, tempMi.refIdx[0], PosY0 );
#endif
        PosY0.x = (PosY0.x & mask);
        PosY0.y = (PosY0.y & mask);
        mi0 = pRefPic0->cs->getMotionInfo(PosY0);
        int ipm0 = pRefPic0->cs->getIpmInfo(PosY0);
        int pocDiff0 = abs(pRefPic0->getPOC() - pu.cu->slice->getPOC());

#if RPR_ENABLE
        pRefPic1 = pu.cu->slice->getRefPic(REF_PIC_LIST_1, tempMi.refIdx[1])->unscaledPic;
#else
        pRefPic1 = pu.cu->slice->getRefPic(REF_PIC_LIST_1, tempMi.refIdx[1]);
#endif
        cMv1 = tempMi.mv[1];
        cMv1.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
        PosY1.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv1.getHor();
        PosY1.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv1.getVer();
        clipColPos(PosY1.x, PosY1.y, pu);
#if RPR_ENABLE
        scalePositionInRef( pu, *pu.cs->pps, REF_PIC_LIST_1, tempMi.refIdx[1], PosY1 );
#endif
        PosY1.x = (PosY1.x & mask);
        PosY1.y = (PosY1.y & mask);
        mi1 = pRefPic1->cs->getMotionInfo(PosY1);
        int ipm1 = pRefPic1->cs->getIpmInfo(PosY1);
        int pocDiff1 = abs(pRefPic1->getPOC() - pu.cu->slice->getPOC());

        if (!mi0.isInter && mi1.isInter)
        {
          ipm = ipm0;
        }
        else if (!mi1.isInter && mi0.isInter)
        {
          ipm = ipm1;
        }
        else if (ipm0 > DC_IDX && ipm1 <= DC_IDX)
        {
          ipm = ipm0;
        }
        else if (ipm0 <= DC_IDX && ipm1 > DC_IDX)
        {
          ipm = ipm1;
        }
        else if (pocDiff0 < pocDiff1)
        {
          ipm = ipm0;
        }
        else if (pocDiff1 < pocDiff0)
        {
          ipm = ipm1;
        }
        else if (pRefPic0->m_prevQP[0] > pRefPic1->m_prevQP[0])
        {
          ipm = ipm1;
        }
        else
        {
          ipm = ipm0;
        }
      }
      ii[x] = ipm;
    }
    ii += ib.stride;
  }
}
#endif

#if !JVET_Z0054_BLK_REF_PIC_REORDER
void PU::applyImv( PredictionUnit& pu, MergeCtx &mrgCtx, InterPrediction *interPred )
{
  if( !pu.mergeFlag )
  {
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    Mv cMvpL0;
#endif
    if( pu.interDir != 2 /* PRED_L1 */ )
    {
#if !JVET_Z0084_IBC_TM
      pu.mvd[0].changeTransPrecAmvr2Internal(pu.cu->imv);
#endif
      unsigned mvpIdx = pu.mvpIdx[0];
      AMVPInfo amvpInfo;
      if (CU::isIBC(*pu.cu))
      {
#if JVET_Z0084_IBC_TM
        pu.mvd[0].changeIbcPrecAmvr2Internal(pu.cu->imv);
#endif
#if JVET_Z0084_IBC_TM && TM_AMVP
        PU::fillIBCMvpCand(pu, amvpInfo, interPred);
#else
        PU::fillIBCMvpCand(pu, amvpInfo);
#endif
      }
      else
      {
#if JVET_Z0084_IBC_TM
        pu.mvd[0].changeTransPrecAmvr2Internal(pu.cu->imv);
#endif
        PU::fillMvpCand(pu, REF_PIC_LIST_0, pu.refIdx[0], amvpInfo
#if TM_AMVP
                      , interPred
#endif
          );
      }
      pu.mvpNum[0] = amvpInfo.numCand;
      pu.mvpIdx[0] = mvpIdx;
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      if (!CU::isIBC(*pu.cu))
      {
        RefPicList eRefList(REF_PIC_LIST_0);
        if (pu.isMvsdApplicable() && pu.mvd[eRefList].isMvsdApplicable())
        {
          if (pu.cu->smvdMode)
          {
            cMvpL0 = amvpInfo.mvCand[pu.mvpIdx[eRefList]];
          }
          else
          {
            std::vector<Mv> cMvdDerivedVec;
            interPred->deriveMvdSign(amvpInfo.mvCand[mvpIdx], pu.mvd[eRefList], pu, eRefList, pu.refIdx[eRefList], cMvdDerivedVec);
            CHECK(pu.mvsdIdx[eRefList] >= cMvdDerivedVec.size(), "");
            int mvsdIdx = pu.mvsdIdx[eRefList];
            Mv cMvd = interPred->deriveMVDFromMVSDIdxTrans(mvsdIdx, cMvdDerivedVec);
            CHECK(cMvd == Mv(0, 0), " zero MVD!");
            pu.mvd[eRefList] = cMvd;
          }
        }
      }
#endif
      pu.mv    [0] = amvpInfo.mvCand[mvpIdx] + pu.mvd[0];
      pu.mv[0].mvCliptoStorageBitDepth();
#if JVET_Z0160_IBC_ZERO_PADDING
      pu.bv = pu.mv[0];
      pu.bv.changePrecision(MV_PRECISION_INTERNAL, MV_PRECISION_INT);
#endif
    }

    if (pu.interDir != 1 /* PRED_L0 */)
    {
      if( !( pu.cu->cs->picHeader->getMvdL1ZeroFlag() && pu.interDir == 3 ) && pu.cu->imv )/* PRED_BI */
      {
        pu.mvd[1].changeTransPrecAmvr2Internal(pu.cu->imv);
      }
      unsigned mvpIdx = pu.mvpIdx[1];
      AMVPInfo amvpInfo;
      PU::fillMvpCand(pu, REF_PIC_LIST_1, pu.refIdx[1], amvpInfo
#if TM_AMVP
                    , interPred
#endif
      );
      pu.mvpNum[1] = amvpInfo.numCand;
      pu.mvpIdx[1] = mvpIdx;
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      RefPicList eRefList(REF_PIC_LIST_1);
      if (pu.isMvsdApplicable() && pu.mvd[eRefList].isMvsdApplicable())
      {
        if (pu.cu->smvdMode)
        {
          std::vector<Mv> cMvdDerivedVec;
          interPred->deriveMvdSignSMVD(cMvpL0, amvpInfo.mvCand[pu.mvpIdx[1]], pu.mvd[REF_PIC_LIST_0], pu, cMvdDerivedVec);
          CHECK(pu.mvsdIdx[REF_PIC_LIST_0] >= cMvdDerivedVec.size(), "");
          int mvsdIdx = pu.mvsdIdx[REF_PIC_LIST_0];
          Mv cMvd = interPred->deriveMVDFromMVSDIdxTrans(mvsdIdx, cMvdDerivedVec);
          CHECK(cMvd == Mv(0, 0), " zero MVD for SMVD!");
          pu.mvd[REF_PIC_LIST_0] = cMvd;
          pu.mv[REF_PIC_LIST_0] = cMvpL0 + pu.mvd[REF_PIC_LIST_0];
          pu.mv[REF_PIC_LIST_0].mvCliptoStorageBitDepth();
          pu.mvd[REF_PIC_LIST_1].set(-pu.mvd[REF_PIC_LIST_0].hor, -pu.mvd[REF_PIC_LIST_0].ver);
        }
        else
        {
          std::vector<Mv> cMvdDerivedVec;
          interPred->deriveMvdSign(amvpInfo.mvCand[mvpIdx], pu.mvd[eRefList], pu, eRefList, pu.refIdx[eRefList], cMvdDerivedVec);
          CHECK(pu.mvsdIdx[eRefList] >= cMvdDerivedVec.size(), "");
          int mvsdIdx = pu.mvsdIdx[eRefList];
          Mv cMvd = interPred->deriveMVDFromMVSDIdxTrans(mvsdIdx, cMvdDerivedVec );
          CHECK(cMvd == Mv(0, 0), " zero MVD!");
          pu.mvd[eRefList] = cMvd;
        }
      }
#endif
      pu.mv    [1] = amvpInfo.mvCand[mvpIdx] + pu.mvd[1];
      pu.mv[1].mvCliptoStorageBitDepth();
    }
  }
  else
  {
    // this function is never called for merge
    THROW("unexpected");
    PU::getInterMergeCandidates(pu, mrgCtx, 0);

    mrgCtx.setMergeInfo( pu, pu.mergeIdx );
  }

  PU::spanMotionInfo( pu, mrgCtx );
}
#endif

bool PU::isBiPredFromDifferentDirEqDistPoc(const PredictionUnit& pu)
{
  if (pu.refIdx[0] >= 0 && pu.refIdx[1] >= 0)
  {
    if (pu.cu->slice->getRefPic(REF_PIC_LIST_0, pu.refIdx[0])->longTerm
      || pu.cu->slice->getRefPic(REF_PIC_LIST_1, pu.refIdx[1])->longTerm)
    {
      return false;
    }
#if JVET_Y0128_NON_CTC
    if ( PU::isBiRefScaled( *pu.cs, pu.refIdx[0], pu.refIdx[1] ) )
    {
      return false;
    }
#endif
    const int poc0 = pu.cu->slice->getRefPOC(REF_PIC_LIST_0, pu.refIdx[0]);
    const int poc1 = pu.cu->slice->getRefPOC(REF_PIC_LIST_1, pu.refIdx[1]);
    const int poc = pu.cu->slice->getPOC();
    if ((poc - poc0)*(poc - poc1) < 0)
    {
      if (abs(poc - poc0) == abs(poc - poc1))
      {
        return true;
      }
    }
  }
  return false;
}

void PU::restrictBiPredMergeCandsOne(PredictionUnit &pu)
{
  if (PU::isBipredRestriction(pu))
  {
    if (pu.interDir == 3)
    {
      pu.interDir = 1;
      pu.refIdx[1] = -1;
      pu.mv[1] = Mv(0, 0);
      pu.cu->BcwIdx = BCW_DEFAULT;
#if MULTI_HYP_PRED
      pu.addHypData.clear();
      pu.numMergedAddHyps = 0;
#endif
    }
  }
}

#if JVET_W0097_GPM_MMVD_TM
#if TM_MRG
void PU::getGeoMergeCandidates(PredictionUnit &pu, MergeCtx& geoMrgCtx, MergeCtx* mergeCtx)
#else
void PU::getGeoMergeCandidates(const PredictionUnit &pu, MergeCtx& geoMrgCtx, MergeCtx* mergeCtx)
#endif
#else
void PU::getGeoMergeCandidates( const PredictionUnit &pu, MergeCtx& geoMrgCtx )
#endif
{
  MergeCtx tmpMergeCtx;

  const uint32_t maxNumMergeCand = pu.cs->sps->getMaxNumMergeCand();
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
  const uint32_t maxNumGeoMergeCand = std::max(pu.cs->sps->getMaxNumMHPCand(), pu.cs->sps->getMaxNumGeoCand());
#endif
  geoMrgCtx.numValidMergeCand = 0;

#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
  for (int32_t i = 0; i < maxNumGeoMergeCand; i++)
#else
  for (int32_t i = 0; i < GEO_MAX_NUM_UNI_CANDS; i++)
#endif
  {
    geoMrgCtx.BcwIdx[i] = BCW_DEFAULT;
    geoMrgCtx.interDirNeighbours[i] = 0;
    geoMrgCtx.mvFieldNeighbours[(i << 1)].refIdx = NOT_VALID;
    geoMrgCtx.mvFieldNeighbours[(i << 1) + 1].refIdx = NOT_VALID;
    geoMrgCtx.mvFieldNeighbours[(i << 1)].mv = Mv();
    geoMrgCtx.mvFieldNeighbours[(i << 1) + 1].mv = Mv();
    geoMrgCtx.useAltHpelIf[i] = false;
#if INTER_LIC
    geoMrgCtx.LICFlags[i] = false;
#endif
  }
#if JVET_W0097_GPM_MMVD_TM
  if (mergeCtx == NULL)
  {
#if TM_MRG
    const bool tmMergeFlag = pu.tmMergeFlag;
    pu.tmMergeFlag = false;
#endif
#endif
  PU::getInterMergeCandidates(pu, tmpMergeCtx, 0);
#if JVET_W0097_GPM_MMVD_TM
#if TM_MRG
  pu.tmMergeFlag = tmMergeFlag;
#endif
  }
  else
  {
    memcpy(tmpMergeCtx.interDirNeighbours, mergeCtx->interDirNeighbours, maxNumMergeCand * sizeof(unsigned char));
    memcpy(tmpMergeCtx.mvFieldNeighbours, mergeCtx->mvFieldNeighbours, (maxNumMergeCand << 1) * sizeof(MvField));
  }
#endif
  for (int32_t i = 0; i < maxNumMergeCand; i++)
  {
    int parity = i & 1;
    if( tmpMergeCtx.interDirNeighbours[i] & (0x01 + parity) )
    {
      geoMrgCtx.interDirNeighbours[geoMrgCtx.numValidMergeCand] = 1 + parity;
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].mv = Mv(0, 0);
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].mv = tmpMergeCtx.mvFieldNeighbours[(i << 1) + parity].mv;
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].refIdx = -1;
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].refIdx = tmpMergeCtx.mvFieldNeighbours[(i << 1) + parity].refIdx;
#if JVET_W0097_GPM_MMVD_TM
      if (geoMrgCtx.xCheckSimilarMotion(geoMrgCtx.numValidMergeCand, PU::getBDMVRMvdThreshold(pu)))
      {
        continue;
      }
#endif
      geoMrgCtx.numValidMergeCand++;
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
      if (geoMrgCtx.numValidMergeCand == maxNumGeoMergeCand)
#else
      if (geoMrgCtx.numValidMergeCand == GEO_MAX_NUM_UNI_CANDS)
#endif
      {
        return;
      }
      continue;
    }

    if (tmpMergeCtx.interDirNeighbours[i] & (0x02 - parity))
    {
      geoMrgCtx.interDirNeighbours[geoMrgCtx.numValidMergeCand] = 2 - parity;
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].mv = tmpMergeCtx.mvFieldNeighbours[(i << 1) + !parity].mv;
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].mv = Mv(0, 0);
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].refIdx = tmpMergeCtx.mvFieldNeighbours[(i << 1) + !parity].refIdx;
      geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].refIdx = -1;
#if JVET_W0097_GPM_MMVD_TM
      if (geoMrgCtx.xCheckSimilarMotion(geoMrgCtx.numValidMergeCand, PU::getBDMVRMvdThreshold(pu)))
      {
        continue;
      }
#endif
      geoMrgCtx.numValidMergeCand++;
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
      if (geoMrgCtx.numValidMergeCand == maxNumGeoMergeCand)
#else
      if (geoMrgCtx.numValidMergeCand == GEO_MAX_NUM_UNI_CANDS)
#endif
      {
        return;
      }
    }
  }
#if JVET_W0097_GPM_MMVD_TM
  // add more parity based geo candidates, in an opposite parity rule
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
  if (geoMrgCtx.numValidMergeCand < maxNumGeoMergeCand)
#else
  if (geoMrgCtx.numValidMergeCand < pu.cs->sps->getMaxNumGeoCand())
#endif
  {
    for (int32_t i = 0; i < maxNumMergeCand; i++)
    {
      int parity = i & 1;
      if (tmpMergeCtx.interDirNeighbours[i] == 3)
      {
        geoMrgCtx.interDirNeighbours[geoMrgCtx.numValidMergeCand] = 2 - parity;
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].mv = tmpMergeCtx.mvFieldNeighbours[(i << 1) + !parity].mv;
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].mv = Mv(0, 0);
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].refIdx = tmpMergeCtx.mvFieldNeighbours[(i << 1) + !parity].refIdx;
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].refIdx = -1;
        if (geoMrgCtx.xCheckSimilarMotion(geoMrgCtx.numValidMergeCand, PU::getBDMVRMvdThreshold(pu)))
        {
          continue;
        }
        geoMrgCtx.numValidMergeCand++;
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
        if (geoMrgCtx.numValidMergeCand == maxNumGeoMergeCand)
#else
        if (geoMrgCtx.numValidMergeCand == pu.cs->sps->getMaxNumGeoCand())
#endif
        {
          return;
        }
      }
    }
  }

  // add at most two average based geo candidates
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
  if (geoMrgCtx.numValidMergeCand == maxNumGeoMergeCand)
#else
  if (geoMrgCtx.numValidMergeCand == pu.cs->sps->getMaxNumGeoCand())
#endif
  {
    // add one L0 cand by averaging the first two available L0 candidates
    int cnt = 0;
    int firstAvailRefIdx = -1;
    Mv  avgMv;
    avgMv.setZero();
    for (int i = 0; i < geoMrgCtx.numValidMergeCand; i++)
    {
      if (cnt == 2)
      {
        break;
      }
      if (geoMrgCtx.interDirNeighbours[i] == 1)
      {
        avgMv += geoMrgCtx.mvFieldNeighbours[i * 2].mv;
        if (firstAvailRefIdx == -1)
        {
          firstAvailRefIdx = geoMrgCtx.mvFieldNeighbours[i * 2].refIdx;
        }
        cnt++;
      }
    }
    if (cnt == 2)
    {
      roundAffineMv(avgMv.hor, avgMv.ver, 1);
      geoMrgCtx.interDirNeighbours[geoMrgCtx.numValidMergeCand] = 1;
      geoMrgCtx.mvFieldNeighbours[geoMrgCtx.numValidMergeCand * 2].setMvField(avgMv, firstAvailRefIdx);
      geoMrgCtx.mvFieldNeighbours[geoMrgCtx.numValidMergeCand * 2 + 1].setMvField(Mv(0, 0), NOT_VALID);

      if (!geoMrgCtx.xCheckSimilarMotion(geoMrgCtx.numValidMergeCand, PU::getBDMVRMvdThreshold(pu)))
      {
        geoMrgCtx.numValidMergeCand++;
      }
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
      if (geoMrgCtx.numValidMergeCand == maxNumGeoMergeCand)
#else
      if (geoMrgCtx.numValidMergeCand == pu.cs->sps->getMaxNumGeoCand())
#endif
      {
        return;
      }
    }

    // add one L1 cand by averaging the first two available L1 candidates
    cnt = 0;
    firstAvailRefIdx = -1;
    avgMv.setZero();
    for (int i = 0; i < geoMrgCtx.numValidMergeCand; i++)
    {
      if (cnt == 2)
      {
        break;
      }
      if (geoMrgCtx.interDirNeighbours[i] == 2)
      {
        avgMv += geoMrgCtx.mvFieldNeighbours[i * 2 + 1].mv;
        if (firstAvailRefIdx == -1)
        {
          firstAvailRefIdx = geoMrgCtx.mvFieldNeighbours[i * 2 + 1].refIdx;
        }
        cnt++;
      }
    }
    if (cnt == 2)
    {
      roundAffineMv(avgMv.hor, avgMv.ver, 1);
      geoMrgCtx.interDirNeighbours[geoMrgCtx.numValidMergeCand] = 2;
      geoMrgCtx.mvFieldNeighbours[geoMrgCtx.numValidMergeCand * 2 + 1].setMvField(avgMv, firstAvailRefIdx);
      geoMrgCtx.mvFieldNeighbours[geoMrgCtx.numValidMergeCand * 2].setMvField(Mv(0, 0), NOT_VALID);
      if (!geoMrgCtx.xCheckSimilarMotion(geoMrgCtx.numValidMergeCand, PU::getBDMVRMvdThreshold(pu)))
      {
        geoMrgCtx.numValidMergeCand++;
      }
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
      if (geoMrgCtx.numValidMergeCand == maxNumGeoMergeCand)
#else
      if (geoMrgCtx.numValidMergeCand == pu.cs->sps->getMaxNumGeoCand())
#endif
      {
        return;
      }
    }
  }
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
  if (geoMrgCtx.numValidMergeCand < maxNumGeoMergeCand)
#else
  if (geoMrgCtx.numValidMergeCand < pu.cs->sps->getMaxNumGeoCand())
#endif
  {
    const Slice &slice = *pu.cs->slice;
#if JVET_Y0065_GPM_INTRA
    int         iNumRefIdx = pu.cs->slice->isInterP() ? slice.getNumRefIdx(REF_PIC_LIST_0) : std::min(slice.getNumRefIdx(REF_PIC_LIST_0), slice.getNumRefIdx(REF_PIC_LIST_1));
#else
    int         iNumRefIdx = std::min(slice.getNumRefIdx(REF_PIC_LIST_0), slice.getNumRefIdx(REF_PIC_LIST_1));
#endif

    int r = 0;
    int refcnt = 0;

#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
    for (int32_t i = geoMrgCtx.numValidMergeCand; i < maxNumGeoMergeCand; i++)
#else
    for (int32_t i = geoMrgCtx.numValidMergeCand; i < pu.cs->sps->getMaxNumGeoCand(); i++)
#endif
    {
#if JVET_Y0065_GPM_INTRA
      int parity = pu.cs->slice->isInterP() ? 0 : (i & 1);
#else
      int parity = i & 1;
#endif
      if (0x01 + parity)
      {
        geoMrgCtx.interDirNeighbours[geoMrgCtx.numValidMergeCand] = 1 + parity;
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].mv = Mv(0, 0);
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].mv = Mv(0, 0);
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + !parity].refIdx = -1;
        geoMrgCtx.mvFieldNeighbours[(geoMrgCtx.numValidMergeCand << 1) + parity].refIdx = r;

        if (refcnt == iNumRefIdx - 1)
        {
          r = 0;
        }
        else
        {
          ++r;
          ++refcnt;
        }

        geoMrgCtx.numValidMergeCand++;
#if JVET_Z0127_SPS_MHP_MAX_MRG_CAND
        if (geoMrgCtx.numValidMergeCand == maxNumGeoMergeCand)
#else
        if (geoMrgCtx.numValidMergeCand == pu.cs->sps->getMaxNumGeoCand())
#endif
        {
          return;
        }
      }
    }
  }
#endif
}

void PU::spanGeoMotionInfo( PredictionUnit &pu, MergeCtx &geoMrgCtx, const uint8_t splitDir, const uint8_t candIdx0, const uint8_t candIdx1)
{
  pu.geoSplitDir  = splitDir;
  pu.geoMergeIdx0 = candIdx0;
  pu.geoMergeIdx1 = candIdx1;
  MotionBuf mb = pu.getMotionBuf();
#if JVET_W0123_TIMD_FUSION
  IpmBuf ib = pu.getIpmBuf();
#endif

#if JVET_Y0065_GPM_INTRA
  bool isIntra0 = candIdx0 >= GEO_MAX_NUM_UNI_CANDS;
  bool isIntra1 = candIdx1 >= GEO_MAX_NUM_UNI_CANDS;
  uint32_t sliceIdx = pu.cs->slice->getIndependentSliceIdx();
#else
  MotionInfo biMv;
  biMv.isInter  = true;
  biMv.sliceIdx = pu.cs->slice->getIndependentSliceIdx();

  if( geoMrgCtx.interDirNeighbours[candIdx0] == 1 && geoMrgCtx.interDirNeighbours[candIdx1] == 2 )
  {
    biMv.interDir  = 3;
    biMv.mv[0]     = geoMrgCtx.mvFieldNeighbours[ candIdx0 << 1     ].mv;
    biMv.mv[1]     = geoMrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv;
    biMv.refIdx[0] = geoMrgCtx.mvFieldNeighbours[ candIdx0 << 1     ].refIdx;
    biMv.refIdx[1] = geoMrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx;
  }
  else if( geoMrgCtx.interDirNeighbours[candIdx0] == 2 && geoMrgCtx.interDirNeighbours[candIdx1] == 1 )
  {
    biMv.interDir  = 3;
    biMv.mv[0]     = geoMrgCtx.mvFieldNeighbours[ candIdx1 << 1     ].mv;
    biMv.mv[1]     = geoMrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].mv;
    biMv.refIdx[0] = geoMrgCtx.mvFieldNeighbours[ candIdx1 << 1     ].refIdx;
    biMv.refIdx[1] = geoMrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].refIdx;
  }
  else if( geoMrgCtx.interDirNeighbours[candIdx0] == 1 && geoMrgCtx.interDirNeighbours[candIdx1] == 1 )
  {
    biMv.interDir = 1;
    biMv.mv[0] = geoMrgCtx.mvFieldNeighbours[candIdx1 << 1].mv;
    biMv.mv[1] = Mv(0, 0);
    biMv.refIdx[0] = geoMrgCtx.mvFieldNeighbours[candIdx1 << 1].refIdx;
    biMv.refIdx[1] = -1;
  }
  else if( geoMrgCtx.interDirNeighbours[candIdx0] == 2 && geoMrgCtx.interDirNeighbours[candIdx1] == 2 )
  {
    biMv.interDir = 2;
    biMv.mv[0] = Mv(0, 0);
    biMv.mv[1] = geoMrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv;
    biMv.refIdx[0] = -1;
    biMv.refIdx[1] = geoMrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx;
  }
#endif

  int16_t angle = g_GeoParams[splitDir][0];
  int tpmMask = 0;
  int lookUpY = 0, motionIdx = 0;
  bool isFlip = angle >= 13 && angle <= 27;
  int distanceIdx = g_GeoParams[splitDir][1];
  int distanceX = angle;
  int distanceY = (distanceX + (GEO_NUM_ANGLES >> 2)) % GEO_NUM_ANGLES;
  int offsetX = (-(int)pu.lwidth()) >> 1;
  int offsetY = (-(int)pu.lheight()) >> 1;

  if (distanceIdx > 0)
  {
    if( angle % 16 == 8 || ( angle % 16 != 0 && pu.lheight() >= pu.lwidth() ) )
    {
      offsetY += angle < 16 ? ( ( distanceIdx * pu.lheight() ) >> 3 ) : -( int ) ( ( distanceIdx * pu.lheight() ) >> 3 );
    }
    else
    {
      offsetX += angle < 16 ? ( ( distanceIdx * pu.lwidth() ) >> 3 ) : -( int ) ( ( distanceIdx * pu.lwidth() ) >> 3 );
    }
  }

  MotionInfo *motionInfo = mb.buf;

  for (int y = 0; y < mb.height; y++)
  {
    lookUpY = (((4 * y + offsetY) << 1) + 5) * g_Dis[distanceY];
    for (int x = 0; x < mb.width; x++)
    {
      motionIdx = (((4 * x + offsetX) << 1) + 5) * g_Dis[distanceX] + lookUpY;
#if JVET_Y0065_GPM_INTRA
      tpmMask = motionIdx <= 0 ? (1 - isFlip) : isFlip;
      if (tpmMask == 0 && isIntra0)
      {
        motionInfo[x].isInter = false;
        motionInfo[x].interDir = MAX_UCHAR;
        motionInfo[x].refIdx[0] = -1;
        motionInfo[x].refIdx[1] = -1;
        motionInfo[x].mv[0] = Mv();
        motionInfo[x].mv[1] = Mv();
        motionInfo[x].sliceIdx = sliceIdx;
      }
#else
      tpmMask = abs(motionIdx) < 32 ? 2 : (motionIdx <= 0 ? (1 - isFlip) : isFlip);
      if (tpmMask == 2)
      {
        motionInfo[x].isInter = true;
        motionInfo[x].interDir = biMv.interDir;
        motionInfo[x].refIdx[0] = biMv.refIdx[0];
        motionInfo[x].refIdx[1] = biMv.refIdx[1];
        motionInfo[x].mv[0] = biMv.mv[0];
        motionInfo[x].mv[1] = biMv.mv[1];
        motionInfo[x].sliceIdx = biMv.sliceIdx;
      }
#endif
      else if (tpmMask == 0)
      {
        motionInfo[x].isInter = true;
        motionInfo[x].interDir = geoMrgCtx.interDirNeighbours[candIdx0];
        motionInfo[x].refIdx[0] = geoMrgCtx.mvFieldNeighbours[candIdx0 << 1].refIdx;
        motionInfo[x].refIdx[1] = geoMrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].refIdx;
        motionInfo[x].mv[0] = geoMrgCtx.mvFieldNeighbours[candIdx0 << 1].mv;
        motionInfo[x].mv[1] = geoMrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].mv;
#if JVET_Y0065_GPM_INTRA
        motionInfo[x].sliceIdx = sliceIdx;
#else
        motionInfo[x].sliceIdx = biMv.sliceIdx;
#endif
      }
#if JVET_Y0065_GPM_INTRA
      else if (tpmMask == 1 && isIntra1)
      {
        motionInfo[x].isInter = false;
        motionInfo[x].interDir = MAX_UCHAR;
        motionInfo[x].refIdx[0] = -1;
        motionInfo[x].refIdx[1] = -1;
        motionInfo[x].mv[0] = Mv();
        motionInfo[x].mv[1] = Mv();
        motionInfo[x].sliceIdx = sliceIdx;
      }
#endif
      else
      {
        motionInfo[x].isInter = true;
        motionInfo[x].interDir = geoMrgCtx.interDirNeighbours[candIdx1];
        motionInfo[x].refIdx[0] = geoMrgCtx.mvFieldNeighbours[candIdx1 << 1].refIdx;
        motionInfo[x].refIdx[1] = geoMrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx;
        motionInfo[x].mv[0] = geoMrgCtx.mvFieldNeighbours[candIdx1 << 1].mv;
        motionInfo[x].mv[1] = geoMrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv;
#if JVET_Y0065_GPM_INTRA
        motionInfo[x].sliceIdx = sliceIdx;
#else
        motionInfo[x].sliceIdx = biMv.sliceIdx;
#endif
      }
    }
    motionInfo += mb.stride;
  }
#if JVET_W0123_TIMD_FUSION
#if JVET_Y0065_GPM_INTRA
  if (pu.gpmIntraFlag)
  {
    const unsigned scale = 4 * std::max<int>(1, 4 * AMVP_DECIMATION_FACTOR / 4);
    const unsigned mask = ~(scale - 1);
    uint8_t* ii = ib.buf;
    int geoIpm[3];
    geoIpm[0] = isIntra0 ? pu.intraMPM[candIdx0 - GEO_MAX_NUM_UNI_CANDS] : -1;
    geoIpm[1] = isIntra1 ? pu.intraMPM[candIdx1 - GEO_MAX_NUM_UNI_CANDS + GEO_MAX_NUM_INTRA_CANDS] : -1;
    geoIpm[2] = (isIntra0 && isIntra1) ? ((candIdx1 < candIdx0) ? geoIpm[1] : geoIpm[0]) : isIntra0 ? geoIpm[0] : geoIpm[1];

    for (int y = 0; y < mb.height; y++)
    {
      lookUpY = (((4 * y + offsetY) << 1) + 5) * g_Dis[distanceY];
      for (int x = 0; x < mb.width; x++)
      {
        motionIdx = (((4 * x + offsetX) << 1) + 5) * g_Dis[distanceX] + lookUpY;
        tpmMask = motionIdx <= 0 ? (1 - isFlip) : isFlip;
        uint8_t ipm;
        if (geoIpm[tpmMask] >= 0)
        {
          ipm = (uint8_t)geoIpm[tpmMask];
        }
        else
        {
          MotionInfo tempMi = mb.at(x, y);
          if (tempMi.interDir != 3)
          {
            CHECK( tempMi.interDir != 1 && tempMi.interDir != 2, "Uncorrect interDir" );
            int list = tempMi.interDir-1;
            Mv cMv = tempMi.mv[list];
            int refIdx = tempMi.refIdx[list];

            cMv.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
            Position PosY;
            PosY.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv.getHor();
            PosY.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv.getVer();
            clipColPos(PosY.x, PosY.y, pu);
#if JVET_Z0067_RPR_ENABLE
            scalePositionInRef(pu, *pu.cs->pps, RefPicList(list), refIdx, PosY);
            PosY.x = (PosY.x & mask);
            PosY.y = (PosY.y & mask);
            ipm = pu.cu->slice->getRefPic(RefPicList(list), refIdx)->unscaledPic->cs->getIpmInfo(PosY);
#else
            PosY.x = (PosY.x & mask);
            PosY.y = (PosY.y & mask);
            ipm = pu.cu->slice->getRefPic(RefPicList(list), refIdx)->cs->getIpmInfo(PosY);
#endif
          }
          else
          {
#if JVET_Z0067_RPR_ENABLE
            const Picture* pRefPic0 = pu.cu->slice->getRefPic(REF_PIC_LIST_0, tempMi.refIdx[0])->unscaledPic;
#else
            Picture* pRefPic0 = pu.cu->slice->getRefPic(REF_PIC_LIST_0, tempMi.refIdx[0]);
#endif
            Mv cMv0 = tempMi.mv[0];
            cMv0.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
            Position PosY0;
            PosY0.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv0.getHor();
            PosY0.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv0.getVer();
            clipColPos(PosY0.x, PosY0.y, pu);
#if JVET_Z0067_RPR_ENABLE
            scalePositionInRef(pu, *pu.cs->pps, REF_PIC_LIST_0, tempMi.refIdx[0], PosY0);
#endif
            PosY0.x = (PosY0.x & mask);
            PosY0.y = (PosY0.y & mask);
            MotionInfo mi0 = pRefPic0->cs->getMotionInfo(PosY0);
            int ipm0 = pRefPic0->cs->getIpmInfo(PosY0);
            int pocDiff0 = abs(pRefPic0->getPOC() - pu.cu->slice->getPOC());

#if JVET_Z0067_RPR_ENABLE
            const Picture* pRefPic1 = pu.cu->slice->getRefPic(REF_PIC_LIST_1, tempMi.refIdx[1])->unscaledPic;
#else
            Picture* pRefPic1 = pu.cu->slice->getRefPic(REF_PIC_LIST_1, tempMi.refIdx[1]);
#endif
            Mv cMv1 = tempMi.mv[1];
            cMv1.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
            Position PosY1;
            PosY1.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv1.getHor();
            PosY1.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv1.getVer();
            clipColPos(PosY1.x, PosY1.y, pu);
#if JVET_Z0067_RPR_ENABLE
            scalePositionInRef(pu, *pu.cs->pps, REF_PIC_LIST_1, tempMi.refIdx[1], PosY1);
#endif
            PosY1.x = (PosY1.x & mask);
            PosY1.y = (PosY1.y & mask);
            MotionInfo mi1 = pRefPic1->cs->getMotionInfo(PosY1);
            int ipm1 = pRefPic1->cs->getIpmInfo(PosY1);
            int pocDiff1 = abs(pRefPic1->getPOC() - pu.cu->slice->getPOC());

            if (!mi0.isInter && mi1.isInter)
            {
              ipm = ipm0;
            }
            else if (!mi1.isInter && mi0.isInter)
            {
              ipm = ipm1;
            }
            else if (ipm0 > DC_IDX && ipm1 <= DC_IDX)
            {
              ipm = ipm0;
            }
            else if (ipm0 <= DC_IDX && ipm1 > DC_IDX)
            {
              ipm = ipm1;
            }
            else if (pocDiff0 < pocDiff1)
            {
              ipm = ipm0;
            }
            else if (pocDiff1 < pocDiff0)
            {
              ipm = ipm1;
            }
            else if (pRefPic0->m_prevQP[0] > pRefPic1->m_prevQP[0])
            {
              ipm = ipm1;
            }
            else
            {
              ipm = ipm0;
            }
          }
        }
        ii[x] = ipm;
      }
      ii += ib.stride;
    }
  }
  else
#endif
  spanIpmInfoInter(pu, mb, ib);
#endif
}

#if JVET_W0097_GPM_MMVD_TM
#if TM_MRG
void PU::spanGeoMMVDMotionInfo(PredictionUnit &pu, MergeCtx &geoMrgCtx, MergeCtx &geoTmMrgCtx0, MergeCtx &geoTmMrgCtx1, const uint8_t splitDir, const uint8_t mergeIdx0, const uint8_t mergeIdx1, const bool tmFlag0, const bool mmvdFlag0, const uint8_t mmvdIdx0, const bool tmFlag1, const bool mmvdFlag1, const uint8_t mmvdIdx1)
#else
void PU::spanGeoMMVDMotionInfo( PredictionUnit &pu, MergeCtx &geoMrgCtx, const uint8_t splitDir, const uint8_t mergeIdx0, const uint8_t mergeIdx1, const bool mmvdFlag0, const uint8_t mmvdIdx0, const bool mmvdFlag1, const uint8_t mmvdIdx1 )
#endif
{
  pu.geoSplitDir = splitDir;
  pu.geoMergeIdx0 = mergeIdx0;
  pu.geoMergeIdx1 = mergeIdx1;
  pu.geoMMVDFlag0 = mmvdFlag0;
  pu.geoMMVDIdx0 = mmvdIdx0;
  pu.geoMMVDFlag1 = mmvdFlag1;
  pu.geoMMVDIdx1 = mmvdIdx1;

#if TM_MRG
  pu.geoTmFlag0 = tmFlag0;
  pu.geoTmFlag1 = tmFlag1;
  MergeCtx *mergeCtx0 = (tmFlag0 ? &geoTmMrgCtx0 : &geoMrgCtx);
  MergeCtx *mergeCtx1 = (tmFlag1 ? &geoTmMrgCtx1 : &geoMrgCtx);
#else
  MergeCtx *mergeCtx0 = &geoMrgCtx;
  MergeCtx *mergeCtx1 = &geoMrgCtx;
#endif

  const int mvShift = MV_FRACTIONAL_BITS_DIFF;
  const bool extMMVD = pu.cs->picHeader->getGPMMMVDTableFlag();
  const int mmvdCands[8] = { 1 << mvShift , 2 << mvShift , 4 << mvShift , 8 << mvShift , 16 << mvShift , 32 << mvShift,  64 << mvShift , 128 << mvShift };
  const int extMmvdCands[9] = { 1 << mvShift , 2 << mvShift , 4 << mvShift , 8 << mvShift , 12 << mvShift , 16 << mvShift, 24 << mvShift, 32 << mvShift, 64 << mvShift };
  const int* refMvdCands = (extMMVD ? extMmvdCands : mmvdCands);
  Mv mvOffset0[2], mvOffset1[2], deltaMv;

  if (mmvdFlag0)
  {
    int fPosStep = (extMMVD ? (mmvdIdx0 >> 3) : (mmvdIdx0 >> 2));
    int fPosPosition = (extMMVD ? (mmvdIdx0 - (fPosStep << 3)) : (mmvdIdx0 - (fPosStep << 2)));

    if (fPosPosition == 0)
    {
      deltaMv = Mv(refMvdCands[fPosStep], 0);
    }
    else if (fPosPosition == 1)
    {
      deltaMv = Mv(-refMvdCands[fPosStep], 0);
    }
    else if (fPosPosition == 2)
    {
      deltaMv = Mv(0, refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 3)
    {
      deltaMv = Mv(0, -refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 4)
    {
      deltaMv = Mv(refMvdCands[fPosStep], refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 5)
    {
      deltaMv = Mv(refMvdCands[fPosStep], -refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 6)
    {
      deltaMv = Mv(-refMvdCands[fPosStep], refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 7)
    {
      deltaMv = Mv(-refMvdCands[fPosStep], -refMvdCands[fPosStep]);
    }

    if (mergeCtx0->interDirNeighbours[mergeIdx0] == 1)
    {
      mvOffset0[0] = deltaMv;
    }
    else
    {
      mvOffset0[1] = deltaMv;
    }
  }

  if (mmvdFlag1)
  {
    int fPosStep = (extMMVD ? (mmvdIdx1 >> 3) : (mmvdIdx1 >> 2));
    int fPosPosition = (extMMVD ? (mmvdIdx1 - (fPosStep << 3)) : (mmvdIdx1 - (fPosStep << 2)));

    if (fPosPosition == 0)
    {
      deltaMv = Mv(refMvdCands[fPosStep], 0);
    }
    else if (fPosPosition == 1)
    {
      deltaMv = Mv(-refMvdCands[fPosStep], 0);
    }
    else if (fPosPosition == 2)
    {
      deltaMv = Mv(0, refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 3)
    {
      deltaMv = Mv(0, -refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 4)
    {
      deltaMv = Mv(refMvdCands[fPosStep], refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 5)
    {
      deltaMv = Mv(refMvdCands[fPosStep], -refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 6)
    {
      deltaMv = Mv(-refMvdCands[fPosStep], refMvdCands[fPosStep]);
    }
    else if (fPosPosition == 7)
    {
      deltaMv = Mv(-refMvdCands[fPosStep], -refMvdCands[fPosStep]);
    }

    if (mergeCtx1->interDirNeighbours[mergeIdx1] == 1)
    {
      mvOffset1[0] = deltaMv;
    }
    else
    {
      mvOffset1[1] = deltaMv;
    }
  }

  MotionBuf mb = pu.getMotionBuf();

#if JVET_Y0065_GPM_INTRA
  bool isIntra0 = mergeIdx0 >= GEO_MAX_NUM_UNI_CANDS;
  bool isIntra1 = mergeIdx1 >= GEO_MAX_NUM_UNI_CANDS;
  uint32_t sliceIdx = pu.cs->slice->getIndependentSliceIdx();
#else
  MotionInfo biMv;
  biMv.isInter = true;
  biMv.sliceIdx = pu.cs->slice->getIndependentSliceIdx();

  if (mergeCtx0->interDirNeighbours[mergeIdx0] == 1 && mergeCtx1->interDirNeighbours[mergeIdx1] == 2)
  {
    biMv.interDir = 3;
    biMv.mv[0] = mergeCtx0->mvFieldNeighbours[mergeIdx0 << 1].mv + mvOffset0[0];
    biMv.mv[1] = mergeCtx1->mvFieldNeighbours[(mergeIdx1 << 1) + 1].mv + mvOffset1[1];
    biMv.refIdx[0] = mergeCtx0->mvFieldNeighbours[mergeIdx0 << 1].refIdx;
    biMv.refIdx[1] = mergeCtx1->mvFieldNeighbours[(mergeIdx1 << 1) + 1].refIdx;
  }
  else if (mergeCtx0->interDirNeighbours[mergeIdx0] == 2 && mergeCtx1->interDirNeighbours[mergeIdx1] == 1)
  {
    biMv.interDir = 3;
    biMv.mv[0] = mergeCtx1->mvFieldNeighbours[mergeIdx1 << 1].mv + mvOffset1[0];
    biMv.mv[1] = mergeCtx0->mvFieldNeighbours[(mergeIdx0 << 1) + 1].mv + mvOffset0[1];
    biMv.refIdx[0] = mergeCtx1->mvFieldNeighbours[mergeIdx1 << 1].refIdx;
    biMv.refIdx[1] = mergeCtx0->mvFieldNeighbours[(mergeIdx0 << 1) + 1].refIdx;
  }
  else if (mergeCtx0->interDirNeighbours[mergeIdx0] == 1 && mergeCtx1->interDirNeighbours[mergeIdx1] == 1)
  {
    biMv.interDir = 1;
    biMv.mv[0] = mergeCtx1->mvFieldNeighbours[mergeIdx1 << 1].mv + mvOffset1[0];
    biMv.mv[1] = Mv(0, 0);
    biMv.refIdx[0] = mergeCtx1->mvFieldNeighbours[mergeIdx1 << 1].refIdx;
    biMv.refIdx[1] = -1;
  }
  else if (mergeCtx0->interDirNeighbours[mergeIdx0] == 2 && mergeCtx1->interDirNeighbours[mergeIdx1] == 2)
  {
    biMv.interDir = 2;
    biMv.mv[0] = Mv(0, 0);
    biMv.mv[1] = mergeCtx1->mvFieldNeighbours[(mergeIdx1 << 1) + 1].mv + mvOffset1[1];
    biMv.refIdx[0] = -1;
    biMv.refIdx[1] = mergeCtx1->mvFieldNeighbours[(mergeIdx1 << 1) + 1].refIdx;
  }
#endif

  int16_t angle = g_GeoParams[splitDir][0];
  int tpmMask = 0;
  int lookUpY = 0, motionIdx = 0;
  bool isFlip = angle >= 13 && angle <= 27;
  int distanceIdx = g_GeoParams[splitDir][1];
  int distanceX = angle;
  int distanceY = (distanceX + (GEO_NUM_ANGLES >> 2)) % GEO_NUM_ANGLES;
  int offsetX = (-(int)pu.lwidth()) >> 1;
  int offsetY = (-(int)pu.lheight()) >> 1;
  if (distanceIdx > 0)
  {
    if (angle % 16 == 8 || (angle % 16 != 0 && pu.lheight() >= pu.lwidth()))
    {
      offsetY += angle < 16 ? ((distanceIdx * pu.lheight()) >> 3) : -(int)((distanceIdx * pu.lheight()) >> 3);
    }
    else
    {
      offsetX += angle < 16 ? ((distanceIdx * pu.lwidth()) >> 3) : -(int)((distanceIdx * pu.lwidth()) >> 3);
    }
  }
  for (int y = 0; y < mb.height; y++)
  {
    lookUpY = (((4 * y + offsetY) << 1) + 5) * g_Dis[distanceY];
    for (int x = 0; x < mb.width; x++)
    {
      motionIdx = (((4 * x + offsetX) << 1) + 5) * g_Dis[distanceX] + lookUpY;
#if JVET_Y0065_GPM_INTRA
      tpmMask = motionIdx <= 0 ? (1 - isFlip) : isFlip;
      if (tpmMask == 0 && isIntra0)
      {
        mb.at(x, y).isInter = false;
        mb.at(x, y).interDir = MAX_UCHAR;
        mb.at(x, y).refIdx[0] = -1;
        mb.at(x, y).refIdx[1] = -1;
        mb.at(x, y).mv[0] = Mv();
        mb.at(x, y).mv[1] = Mv();
        mb.at(x, y).sliceIdx = sliceIdx;
      }
#else
      tpmMask = abs(motionIdx) < 32 ? 2 : (motionIdx <= 0 ? (1 - isFlip) : isFlip);
      if (tpmMask == 2)
      {
        mb.at(x, y).isInter = true;
        mb.at(x, y).interDir = biMv.interDir;
        mb.at(x, y).refIdx[0] = biMv.refIdx[0];
        mb.at(x, y).refIdx[1] = biMv.refIdx[1];
        mb.at(x, y).mv[0] = biMv.mv[0];
        mb.at(x, y).mv[1] = biMv.mv[1];
        mb.at(x, y).sliceIdx = biMv.sliceIdx;
      }
#endif
      else if (tpmMask == 0)
      {
        mb.at(x, y).isInter = true;
        mb.at(x, y).interDir = mergeCtx0->interDirNeighbours[mergeIdx0];
        mb.at(x, y).refIdx[0] = mergeCtx0->mvFieldNeighbours[mergeIdx0 << 1].refIdx;
        mb.at(x, y).refIdx[1] = mergeCtx0->mvFieldNeighbours[(mergeIdx0 << 1) + 1].refIdx;
        mb.at(x, y).mv[0] = mergeCtx0->mvFieldNeighbours[mergeIdx0 << 1].mv + mvOffset0[0];
        mb.at(x, y).mv[1] = mergeCtx0->mvFieldNeighbours[(mergeIdx0 << 1) + 1].mv + mvOffset0[1];
#if JVET_Y0065_GPM_INTRA
        mb.at(x, y).sliceIdx = sliceIdx;
#else
        mb.at(x, y).sliceIdx = biMv.sliceIdx;
#endif
      }
#if JVET_Y0065_GPM_INTRA
      else if (tpmMask == 1 && isIntra1)
      {
        mb.at(x, y).isInter = false;
        mb.at(x, y).interDir = MAX_UCHAR;
        mb.at(x, y).refIdx[0] = -1;
        mb.at(x, y).refIdx[1] = -1;
        mb.at(x, y).mv[0] = Mv();
        mb.at(x, y).mv[1] = Mv();
        mb.at(x, y).sliceIdx = sliceIdx;
      }
#endif
      else
      {
        mb.at(x, y).isInter = true;
        mb.at(x, y).interDir = mergeCtx1->interDirNeighbours[mergeIdx1];
        mb.at(x, y).refIdx[0] = mergeCtx1->mvFieldNeighbours[mergeIdx1 << 1].refIdx;
        mb.at(x, y).refIdx[1] = mergeCtx1->mvFieldNeighbours[(mergeIdx1 << 1) + 1].refIdx;
        mb.at(x, y).mv[0] = mergeCtx1->mvFieldNeighbours[mergeIdx1 << 1].mv + mvOffset1[0];
        mb.at(x, y).mv[1] = mergeCtx1->mvFieldNeighbours[(mergeIdx1 << 1) + 1].mv + mvOffset1[1];
#if JVET_Y0065_GPM_INTRA
        mb.at(x, y).sliceIdx = sliceIdx;
#else
        mb.at(x, y).sliceIdx = biMv.sliceIdx;
#endif
      }
    }
  }

#if JVET_W0123_TIMD_FUSION
  IpmBuf ib = pu.getIpmBuf();
#if JVET_Y0065_GPM_INTRA
  if (pu.gpmIntraFlag)
  {
    const unsigned scale = 4 * std::max<int>(1, 4 * AMVP_DECIMATION_FACTOR / 4);
    const unsigned mask = ~(scale - 1);
    uint8_t* ii = ib.buf;
    int geoIpm[3];
    geoIpm[0] = isIntra0 ? pu.intraMPM[mergeIdx0 - GEO_MAX_NUM_UNI_CANDS] : -1;
    geoIpm[1] = isIntra1 ? pu.intraMPM[mergeIdx1 - GEO_MAX_NUM_UNI_CANDS + GEO_MAX_NUM_INTRA_CANDS] : -1;
    geoIpm[2] = (isIntra0 && isIntra1) ? ((mergeIdx1 < mergeIdx0) ? geoIpm[1] : geoIpm[0]) : isIntra0 ? geoIpm[0] : geoIpm[1];

    for (int y = 0; y < mb.height; y++)
    {
      lookUpY = (((4 * y + offsetY) << 1) + 5) * g_Dis[distanceY];
      for (int x = 0; x < mb.width; x++)
      {
        motionIdx = (((4 * x + offsetX) << 1) + 5) * g_Dis[distanceX] + lookUpY;
        tpmMask = motionIdx <= 0 ? (1 - isFlip) : isFlip;
        uint8_t ipm;
        if (geoIpm[tpmMask] >= 0)
        {
          ipm = (uint8_t)geoIpm[tpmMask];
        }
        else
        {
          MotionInfo tempMi = mb.at(x, y);
          if (tempMi.interDir != 3)
          {
            CHECK( tempMi.interDir != 1 && tempMi.interDir != 2, "Uncorrect interDir" );
            int list = tempMi.interDir-1;
            Mv cMv = tempMi.mv[list];
            int refIdx = tempMi.refIdx[list];

            cMv.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
            Position PosY;
            PosY.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv.getHor();
            PosY.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv.getVer();
            clipColPos(PosY.x, PosY.y, pu);
#if JVET_Z0067_RPR_ENABLE
            scalePositionInRef(pu, *pu.cs->pps, RefPicList(list), refIdx, PosY);
            PosY.x = (PosY.x & mask);
            PosY.y = (PosY.y & mask);
            ipm = pu.cu->slice->getRefPic(RefPicList(list), refIdx)->unscaledPic->cs->getIpmInfo(PosY);
#else
            PosY.x = (PosY.x & mask);
            PosY.y = (PosY.y & mask);
            ipm = pu.cu->slice->getRefPic(RefPicList(list), refIdx)->cs->getIpmInfo(PosY);
#endif
          }
          else
          {
#if JVET_Z0067_RPR_ENABLE
            const Picture* pRefPic0 = pu.cu->slice->getRefPic(REF_PIC_LIST_0, tempMi.refIdx[0])->unscaledPic;
#else
            Picture* pRefPic0 = pu.cu->slice->getRefPic(REF_PIC_LIST_0, tempMi.refIdx[0]);
#endif
            Mv cMv0 = tempMi.mv[0];
            cMv0.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
            Position PosY0;
            PosY0.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv0.getHor();
            PosY0.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv0.getVer();
            clipColPos(PosY0.x, PosY0.y, pu);
#if JVET_Z0067_RPR_ENABLE
            scalePositionInRef(pu, *pu.cs->pps, REF_PIC_LIST_0, tempMi.refIdx[0], PosY0);
#endif
            PosY0.x = (PosY0.x & mask);
            PosY0.y = (PosY0.y & mask);
            MotionInfo mi0 = pRefPic0->cs->getMotionInfo(PosY0);
            int ipm0 = pRefPic0->cs->getIpmInfo(PosY0);
            int pocDiff0 = abs(pRefPic0->getPOC() - pu.cu->slice->getPOC());

#if JVET_Z0067_RPR_ENABLE
            const Picture* pRefPic1 = pu.cu->slice->getRefPic(REF_PIC_LIST_1, tempMi.refIdx[1])->unscaledPic;
#else
            Picture* pRefPic1 = pu.cu->slice->getRefPic(REF_PIC_LIST_1, tempMi.refIdx[1]);
#endif
            Mv cMv1 = tempMi.mv[1];
            cMv1.changePrecision(MV_PRECISION_SIXTEENTH, MV_PRECISION_INT);
            Position PosY1;
            PosY1.x = pu.Y().x + (x << MIN_CU_LOG2) + cMv1.getHor();
            PosY1.y = pu.Y().y + (y << MIN_CU_LOG2) + cMv1.getVer();
            clipColPos(PosY1.x, PosY1.y, pu);
#if JVET_Z0067_RPR_ENABLE
            scalePositionInRef(pu, *pu.cs->pps, REF_PIC_LIST_1, tempMi.refIdx[1], PosY1);
#endif
            PosY1.x = (PosY1.x & mask);
            PosY1.y = (PosY1.y & mask);
            MotionInfo mi1 = pRefPic1->cs->getMotionInfo(PosY1);
            int ipm1 = pRefPic1->cs->getIpmInfo(PosY1);
            int pocDiff1 = abs(pRefPic1->getPOC() - pu.cu->slice->getPOC());

            if (!mi0.isInter && mi1.isInter)
            {
              ipm = ipm0;
            }
            else if (!mi1.isInter && mi0.isInter)
            {
              ipm = ipm1;
            }
            else if (ipm0 > DC_IDX && ipm1 <= DC_IDX)
            {
              ipm = ipm0;
            }
            else if (ipm0 <= DC_IDX && ipm1 > DC_IDX)
            {
              ipm = ipm1;
            }
            else if (pocDiff0 < pocDiff1)
            {
              ipm = ipm0;
            }
            else if (pocDiff1 < pocDiff0)
            {
              ipm = ipm1;
            }
            else if (pRefPic0->m_prevQP[0] > pRefPic1->m_prevQP[0])
            {
              ipm = ipm1;
            }
            else
            {
              ipm = ipm0;
            }
          }
        }
        ii[x] = ipm;
      }
      ii += ib.stride;
    }
  }
  else
#endif
  spanIpmInfoInter( pu, mb, ib );
#endif
}
#endif

#if JVET_Z0054_BLK_REF_PIC_REORDER
bool PU::useRefCombList(const PredictionUnit &pu)
{
  if (!pu.cs->sps->getUseARL())
  {
    return false;
  }
  if (CU::isIBC(*pu.cu))
  {
    return false;
  }
  if (!pu.mergeFlag && (pu.interDir == 1 || pu.interDir == 2))
  {
    return true;
  }
#if JVET_X0083_BM_AMVP_MERGE_MODE
  if (pu.amvpMergeModeFlag[REF_PIC_LIST_0] || pu.amvpMergeModeFlag[REF_PIC_LIST_1])
  {
    return true;
  }
#endif
  return false;
}

bool PU::useRefPairList(const PredictionUnit &pu)
{
  if (!pu.cs->sps->getUseARL())
  {
    return false;
  }
  if (CU::isIBC(*pu.cu))
  {
    return false;
  }
  if (pu.interDir != 3 || pu.mergeFlag)
  {
    return false;
  }
  if (pu.cu->smvdMode)
  {
    return false;
  }
#if JVET_X0083_BM_AMVP_MERGE_MODE
  if (pu.amvpMergeModeFlag[REF_PIC_LIST_0] || pu.amvpMergeModeFlag[REF_PIC_LIST_1])
  {
    return false;
  }
#endif
  return true;
}
#endif

bool CU::hasSubCUNonZeroMVd( const CodingUnit& cu )
{
  bool bNonZeroMvd = false;

  for( const auto &pu : CU::traversePUs( cu ) )
  {
    if( ( !pu.mergeFlag ) && ( !cu.skip ) )
    {
      if( pu.interDir != 2 /* PRED_L1 */ )
      {
        bNonZeroMvd |= pu.mvd[REF_PIC_LIST_0].getHor() != 0;
        bNonZeroMvd |= pu.mvd[REF_PIC_LIST_0].getVer() != 0;
      }
      if( pu.interDir != 1 /* PRED_L0 */ )
      {
        if( !pu.cu->cs->picHeader->getMvdL1ZeroFlag() || pu.interDir != 3 /* PRED_BI */ )
        {
          bNonZeroMvd |= pu.mvd[REF_PIC_LIST_1].getHor() != 0;
          bNonZeroMvd |= pu.mvd[REF_PIC_LIST_1].getVer() != 0;
        }
      }
    }
  }

  return bNonZeroMvd;
}

bool CU::hasSubCUNonZeroAffineMVd( const CodingUnit& cu )
{
  bool nonZeroAffineMvd = false;

  if ( !cu.affine || cu.firstPU->mergeFlag )
  {
    return false;
  }

  for ( const auto &pu : CU::traversePUs( cu ) )
  {
    if ( ( !pu.mergeFlag ) && ( !cu.skip ) )
    {
      if ( pu.interDir != 2 /* PRED_L1 */ )
      {
        for ( int i = 0; i < ( cu.affineType == AFFINEMODEL_6PARAM ? 3 : 2 ); i++ )
        {
          nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_0][i].getHor() != 0;
          nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_0][i].getVer() != 0;
        }
      }

      if ( pu.interDir != 1 /* PRED_L0 */ )
      {
        if ( !pu.cu->cs->picHeader->getMvdL1ZeroFlag() || pu.interDir != 3 /* PRED_BI */ )
        {
          for ( int i = 0; i < ( cu.affineType == AFFINEMODEL_6PARAM ? 3 : 2 ); i++ )
          {
            nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_1][i].getHor() != 0;
            nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_1][i].getVer() != 0;
          }
        }
      }
    }
  }

  return nonZeroAffineMvd;
}

uint8_t CU::getSbtInfo( uint8_t idx, uint8_t pos )
{
  return ( pos << 4 ) + ( idx << 0 );
}

uint8_t CU::getSbtIdx( const uint8_t sbtInfo )
{
  return ( sbtInfo >> 0 ) & 0xf;
}

uint8_t CU::getSbtPos( const uint8_t sbtInfo )
{
  return ( sbtInfo >> 4 ) & 0x3;
}

uint8_t CU::getSbtMode( uint8_t sbtIdx, uint8_t sbtPos )
{
  uint8_t sbtMode = 0;
  switch( sbtIdx )
  {
  case SBT_VER_HALF: sbtMode = sbtPos + SBT_VER_H0;  break;
  case SBT_HOR_HALF: sbtMode = sbtPos + SBT_HOR_H0;  break;
  case SBT_VER_QUAD: sbtMode = sbtPos + SBT_VER_Q0;  break;
  case SBT_HOR_QUAD: sbtMode = sbtPos + SBT_HOR_Q0;  break;
  default:           assert( 0 );
  }

  assert( sbtMode < NUMBER_SBT_MODE );
  return sbtMode;
}

uint8_t CU::getSbtIdxFromSbtMode( uint8_t sbtMode )
{
  if( sbtMode <= SBT_VER_H1 )
  {
    return SBT_VER_HALF;
  }
  else if( sbtMode <= SBT_HOR_H1 )
  {
    return SBT_HOR_HALF;
  }
  else if( sbtMode <= SBT_VER_Q1 )
  {
    return SBT_VER_QUAD;
  }
  else if( sbtMode <= SBT_HOR_Q1 )
  {
    return SBT_HOR_QUAD;
  }
  else
  {
    assert( 0 );
    return 0;
  }
}

uint8_t CU::getSbtPosFromSbtMode( uint8_t sbtMode )
{
  if( sbtMode <= SBT_VER_H1 )
  {
    return sbtMode - SBT_VER_H0;
  }
  else if( sbtMode <= SBT_HOR_H1 )
  {
    return sbtMode - SBT_HOR_H0;
  }
  else if( sbtMode <= SBT_VER_Q1 )
  {
    return sbtMode - SBT_VER_Q0;
  }
  else if( sbtMode <= SBT_HOR_Q1 )
  {
    return sbtMode - SBT_HOR_Q0;
  }
  else
  {
    assert( 0 );
    return 0;
  }
}

uint8_t CU::targetSbtAllowed( uint8_t sbtIdx, uint8_t sbtAllowed )
{
  uint8_t val = 0;
  switch( sbtIdx )
  {
  case SBT_VER_HALF: val = ( ( sbtAllowed >> SBT_VER_HALF ) & 0x1 ); break;
  case SBT_HOR_HALF: val = ( ( sbtAllowed >> SBT_HOR_HALF ) & 0x1 ); break;
  case SBT_VER_QUAD: val = ( ( sbtAllowed >> SBT_VER_QUAD ) & 0x1 ); break;
  case SBT_HOR_QUAD: val = ( ( sbtAllowed >> SBT_HOR_QUAD ) & 0x1 ); break;
  default:           CHECK( 1, "unknown SBT type" );
  }
  return val;
}

uint8_t CU::numSbtModeRdo( uint8_t sbtAllowed )
{
  uint8_t num = 0;
  uint8_t sum = 0;
  num = targetSbtAllowed( SBT_VER_HALF, sbtAllowed ) + targetSbtAllowed( SBT_HOR_HALF, sbtAllowed );
  sum += std::min( SBT_NUM_RDO, ( num << 1 ) );
  num = targetSbtAllowed( SBT_VER_QUAD, sbtAllowed ) + targetSbtAllowed( SBT_HOR_QUAD, sbtAllowed );
  sum += std::min( SBT_NUM_RDO, ( num << 1 ) );
  return sum;
}

bool CU::isSbtMode( const uint8_t sbtInfo )
{
  uint8_t sbtIdx = getSbtIdx( sbtInfo );
  return sbtIdx >= SBT_VER_HALF && sbtIdx <= SBT_HOR_QUAD;
}

bool CU::isSameSbtSize( const uint8_t sbtInfo1, const uint8_t sbtInfo2 )
{
  uint8_t sbtIdx1 = getSbtIdxFromSbtMode( sbtInfo1 );
  uint8_t sbtIdx2 = getSbtIdxFromSbtMode( sbtInfo2 );
  if( sbtIdx1 == SBT_HOR_HALF || sbtIdx1 == SBT_VER_HALF )
  {
    return sbtIdx2 == SBT_HOR_HALF || sbtIdx2 == SBT_VER_HALF;
  }
  else if( sbtIdx1 == SBT_HOR_QUAD || sbtIdx1 == SBT_VER_QUAD )
  {
    return sbtIdx2 == SBT_HOR_QUAD || sbtIdx2 == SBT_VER_QUAD;
  }
  else
  {
    return false;
  }
}

bool CU::isPredRegDiffFromTB(const CodingUnit &cu, const ComponentID compID)
{
  return (compID == COMPONENT_Y)
         && (cu.ispMode == VER_INTRA_SUBPARTITIONS
             && CU::isMinWidthPredEnabledForBlkSize(cu.blocks[compID].width, cu.blocks[compID].height));
}

bool CU::isMinWidthPredEnabledForBlkSize(const int w, const int h)
{
  return ((w == 8 && h > 4) || w == 4);
}

bool CU::isFirstTBInPredReg(const CodingUnit& cu, const ComponentID compID, const CompArea &area)
{
  return (compID == COMPONENT_Y) && cu.ispMode && ((area.topLeft().x - cu.Y().topLeft().x) % PRED_REG_MIN_WIDTH == 0);
}

void CU::adjustPredArea(CompArea &area)
{
  area.width = std::max<int>(PRED_REG_MIN_WIDTH, area.width);
}

#if MULTI_HYP_PRED
AMVPInfo PU::getMultiHypMVPCandsMerge(PredictionUnit &pu, const RefPicList eRefPicList, const int refIdx)
{
  CHECK(!pu.mergeFlag, "this function is only for merge");

  AMVPInfo amvpInfo;
  PU::fillMvpCand(pu, eRefPicList, refIdx, amvpInfo);

  return amvpInfo;
}

AMVPInfo PU::getMultiHypMVPCandsAMVP(PredictionUnit &pu, const RefPicList eRefPicList, const int refIdx)
{
  CHECK(pu.mergeFlag, "this function is only for AMVP");
  CHECK(pu.interDir != 3, "Multi-Hyp is only allowed for bi-predictive AMVP mode")

  const int refListForAMVP = eRefPicList;
  const int refIdxForAMVP = refIdx;

  AMVPInfo amvpInfo;
  if (pu.cu->affine) // TODO: check performance impact
  {
    AffineAMVPInfo affineAMVPInfo;
    PU::fillAffineMvpCand(pu, RefPicList(refListForAMVP), refIdxForAMVP, affineAMVPInfo);
    amvpInfo.numCand = affineAMVPInfo.numCand;
    memcpy(amvpInfo.mvCand, affineAMVPInfo.mvCandLT, sizeof(amvpInfo.mvCand));
  }
  else
  {
    PU::fillMvpCand(pu, RefPicList(refListForAMVP), refIdxForAMVP, amvpInfo);
  }

  return amvpInfo;
}

AMVPInfo PU::getMultiHypMVPCands(PredictionUnit &pu, const MultiHypPredictionData &mhData)
{
  CHECK(mhData.isMrg, "mhData is merge mode");
  const auto mhRefIdxForAMVPList = mhData.refIdx;
  const auto &MHRefPics = pu.cs->slice->getMultiHypRefPicList();
  CHECK(MHRefPics.empty(), "Multi Hyp: MHRefPics.empty()");
  const auto eRefPicList = RefPicList(MHRefPics[mhRefIdxForAMVPList].refList);
  const int iRefIdx = MHRefPics[mhRefIdxForAMVPList].refIdx;

  return (pu.mergeFlag ? PU::getMultiHypMVPCandsMerge : PU::getMultiHypMVPCandsAMVP)(pu, eRefPicList, iRefIdx);
}

Mv PU::getMultiHypMVP(PredictionUnit &pu, const MultiHypPredictionData &mhData)
{
  const auto amvpInfo = getMultiHypMVPCands(pu, mhData);
  CHECK(mhData.mvpIdx < 0, "Multi Hyp: mhData.mvpIdx < 0");
  CHECK(mhData.mvpIdx >= amvpInfo.numCand, "Multi Hyp: mhData.mvpIdx >= amvpInfo.numCand");
  const Mv mvp = amvpInfo.mvCand[mhData.mvpIdx];
  return mvp;
}
#endif


bool CU::isBcwIdxCoded( const CodingUnit &cu )
{
  if( cu.cs->sps->getUseBcw() == false )
  {
    CHECK(cu.BcwIdx != BCW_DEFAULT, "Error: cu.BcwIdx != BCW_DEFAULT");
    return false;
  }
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  if (cu.skip)
  {
    return false;
  }
#endif
  if (cu.predMode == MODE_IBC)
  {
    return false;
  }

  if( cu.predMode == MODE_INTRA || cu.cs->slice->isInterP() )
  {
    return false;
  }

  if( cu.lwidth() * cu.lheight() < BCW_SIZE_CONSTRAINT )
  {
    return false;
  }

  if( !cu.firstPU->mergeFlag )
  {
    if( cu.firstPU->interDir == 3 )
    {
#if !JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
      const int refIdx0 = cu.firstPU->refIdx[REF_PIC_LIST_0];
      const int refIdx1 = cu.firstPU->refIdx[REF_PIC_LIST_1];
      
      const WPScalingParam *wp0 = cu.cs->slice->getWpScaling(REF_PIC_LIST_0, refIdx0);
      const WPScalingParam *wp1 = cu.cs->slice->getWpScaling(REF_PIC_LIST_1, refIdx1);
      
      return !(WPScalingParam::isWeighted(wp0) || WPScalingParam::isWeighted(wp1));
#else
      return true;
#endif
    }
  }

  return false;
}

uint8_t CU::getValidBcwIdx( const CodingUnit &cu )
{
  if( cu.firstPU->interDir == 3 && !cu.firstPU->mergeFlag )
  {
    return cu.BcwIdx;
  }
  else if( cu.firstPU->interDir == 3 && cu.firstPU->mergeFlag && cu.firstPU->mergeType == MRG_TYPE_DEFAULT_N )
  {
    // This is intended to do nothing here.
  }
  else if( cu.firstPU->mergeFlag && cu.firstPU->mergeType == MRG_TYPE_SUBPU_ATMVP )
  {
    CHECK(cu.BcwIdx != BCW_DEFAULT, " cu.BcwIdx != BCW_DEFAULT ");
  }
  else
  {
    CHECK(cu.BcwIdx != BCW_DEFAULT, " cu.BcwIdx != BCW_DEFAULT ");
  }

  return BCW_DEFAULT;
}

void CU::setBcwIdx( CodingUnit &cu, uint8_t uh )
{
  int8_t uhCnt = 0;

  if( cu.firstPU->interDir == 3 && !cu.firstPU->mergeFlag )
  {
    cu.BcwIdx = uh;
    ++uhCnt;
  }
  else if( cu.firstPU->interDir == 3 && cu.firstPU->mergeFlag && cu.firstPU->mergeType == MRG_TYPE_DEFAULT_N )
  {
    // This is intended to do nothing here.
  }
  else if( cu.firstPU->mergeFlag && cu.firstPU->mergeType == MRG_TYPE_SUBPU_ATMVP )
  {
    cu.BcwIdx = BCW_DEFAULT;
  }
  else
  {
    cu.BcwIdx = BCW_DEFAULT;
  }

  CHECK(uhCnt <= 0, " uhCnt <= 0 ");
}

uint8_t CU::deriveBcwIdx( uint8_t bcwLO, uint8_t bcwL1 )
{
  if( bcwLO == bcwL1 )
  {
    return bcwLO;
  }
  const int8_t w0 = getBcwWeight(bcwLO, REF_PIC_LIST_0);
  const int8_t w1 = getBcwWeight(bcwL1, REF_PIC_LIST_1);
  const int8_t th = g_BcwWeightBase >> 1;
  const int8_t off = 1;

  if( w0 == w1 || (w0 < (th - off) && w1 < (th - off)) || (w0 >(th + off) && w1 >(th + off)) )
  {
    return BCW_DEFAULT;
  }
  else
  {
    if( w0 > w1 )
    {
      return ( w0 >= th ? bcwLO : bcwL1 );
    }
    else
    {
      return ( w1 >= th ? bcwL1 : bcwLO );
    }
  }
}

bool CU::bdpcmAllowed( const CodingUnit& cu, const ComponentID compID )
{
  SizeType transformSkipMaxSize = 1 << cu.cs->sps->getLog2MaxTransformSkipBlockSize();

  bool bdpcmAllowed = cu.cs->sps->getBDPCMEnabledFlag();
       bdpcmAllowed &= CU::isIntra( cu );
       if (isLuma(compID))
       {
         bdpcmAllowed &= (cu.lwidth() <= transformSkipMaxSize && cu.lheight() <= transformSkipMaxSize);
#if ENABLE_DIMD
         bdpcmAllowed &= !cu.dimd;
#endif
       }
       else
       {
         bdpcmAllowed &=
           (cu.chromaSize().width <= transformSkipMaxSize && cu.chromaSize().height <= transformSkipMaxSize)
           && !cu.colorTransform;
       }
  return bdpcmAllowed;
}

bool CU::isMTSAllowed(const CodingUnit &cu, const ComponentID compID)
{
  SizeType tsMaxSize = 1 << cu.cs->sps->getLog2MaxTransformSkipBlockSize();
  const int maxSize  = CU::isIntra( cu ) ? MTS_INTRA_MAX_CU_SIZE : MTS_INTER_MAX_CU_SIZE;
  const int cuWidth  = cu.blocks[0].lumaSize().width;
  const int cuHeight = cu.blocks[0].lumaSize().height;
  bool mtsAllowed    = cu.chType == CHANNEL_TYPE_LUMA && compID == COMPONENT_Y;

  mtsAllowed &= CU::isIntra( cu ) ? cu.cs->sps->getUseIntraMTS() : cu.cs->sps->getUseInterMTS() && CU::isInter( cu );
  mtsAllowed &= cuWidth <= maxSize && cuHeight <= maxSize;
  mtsAllowed &= !cu.ispMode;
  mtsAllowed &= !cu.sbtInfo;
#if JVET_V0130_INTRA_TMP
  mtsAllowed &= !cu.tmpFlag;
#endif
#if JVET_W0123_TIMD_FUSION
  mtsAllowed &= !(cu.timd && cu.firstPU->multiRefIdx);
#endif
  mtsAllowed &= !(cu.bdpcmMode && cuWidth <= tsMaxSize && cuHeight <= tsMaxSize);
  return mtsAllowed;
}

#if ENABLE_OBMC
unsigned int PU::getSameNeigMotion(PredictionUnit &pu, MotionInfo& mi, Position off, int iDir, int& iLength, int iMaxLength)
{

  PredictionUnit* tmpPu = nullptr;
  Position posNeighborMotion = Position(0, 0);
  const Position  posSubBlock(pu.lumaPos().offset(off));

  const unsigned int  uiMinCUW = pu.cs->pcv->minCUWidth;
  if (iDir == 0)
  {
    posNeighborMotion = posSubBlock.offset(0, -1);
  }
  else if (iDir == 1)
  {
    posNeighborMotion = posSubBlock.offset(-1, 0);
  }
  else
  {
    THROW("iDir not supported");
  }
  tmpPu = pu.cs->getPU(posNeighborMotion, pu.chType);

  //if neighbor pu is not available, return.
  if (!tmpPu)
  {
    return 0;
  }

  mi = tmpPu->getMotionInfo(posNeighborMotion);

  iLength = 1;
  if (mi.isInter && !mi.isIBCmot) // inter
  {
    MotionInfo currMi = pu.getMotionInfo(posSubBlock);
    // check similar MV

    bool bIsSimilarMV = PU::identicalMvOBMC(currMi, mi, pu.cs->slice->getCheckLDC());
    for (unsigned int idx = 1; idx < iMaxLength; idx++)
    {

      if (iDir == 0)//top
      {
        posNeighborMotion = posSubBlock.offset(idx*uiMinCUW, -1);
      }
      else//left
      {
        posNeighborMotion = posSubBlock.offset(-1, idx*uiMinCUW);
      }
      tmpPu = pu.cs->getPU(posNeighborMotion, pu.chType);

      MotionInfo miNeigh = tmpPu->getMotionInfo(posNeighborMotion);

      bool       bSameMv = true;

      if (mi.interDir != miNeigh.interDir || miNeigh.isIBCmot)
      {
        bSameMv = false;
      }
      else
      {
        for (unsigned int iRefList = 0; iRefList < 2; iRefList++)
        {
          if (miNeigh.interDir & (1 << iRefList))
          {
            if (!(miNeigh.mv[iRefList] == mi.mv[iRefList] && miNeigh.refIdx[iRefList] == mi.refIdx[iRefList]))
            {
              bSameMv = false; break;
            }
          }
        }
      }
      if (bSameMv)
      {
        iLength++;
      }
      else
      {
        break;
      }
    }
    return 2 + bIsSimilarMV; // 2: not skip; 3: skip
  }
  else             // intra
  {
    for (unsigned int idx = 1; idx < iMaxLength; idx++)
    {
      if (iDir == 0)
      {
        posNeighborMotion = posSubBlock.offset(idx*uiMinCUW, -1);
      }
      else
      {
        posNeighborMotion = posSubBlock.offset(-1, idx*uiMinCUW);
      }
      tmpPu = pu.cs->getPU(posNeighborMotion, pu.chType);

      if (!tmpPu->getMotionInfo(posNeighborMotion).isInter || tmpPu->getMotionInfo(posNeighborMotion).isIBCmot)
      {
        iLength++;
      }
      else
      {
        break;
      }
    }
    return 1;
  }
}

bool PU::identicalMvOBMC(MotionInfo curMI, MotionInfo neighMI, bool bLD)
{
  if (curMI.interDir != neighMI.interDir)
  {
    return false;
  }
  else
  {
    for (int i = 0; i < 2; i++)
    {
      if (curMI.refIdx[i] != -1)
      {
        if (!(curMI.mv[i] == neighMI.mv[i] && curMI.refIdx[i] == neighMI.refIdx[i]))
        {
          return false;
        }
      }
    }
    return true;
  }
}

bool PU::getNeighborMotion(PredictionUnit &pu, MotionInfo& mi, Position off, Size unitSize, int iDir)
{
  PredictionUnit* tmpPu = &pu;
  Position posNeighborMotion = Position(0, 0);
  const Position posSubBlock(pu.lumaPos().offset(off));


  if (iDir == 0) //above
  {
    posNeighborMotion = posSubBlock.offset(0, -1);
  }
  else if (iDir == 1) //left
  {
    posNeighborMotion = posSubBlock.offset(-1, 0);
  }
  else if (iDir == 2) //below
  {
    posNeighborMotion = posSubBlock.offset(0, unitSize.height);
    CHECK(pu.cu != tmpPu->cu, "Got a PU from a different CU when fetching a below PU");
  }
  else if (iDir == 3) //right
  {
    posNeighborMotion = posSubBlock.offset(unitSize.width, 0);
    CHECK(pu.cu != tmpPu->cu, "Got a PU from a different CU when fetching a right PU");
  }

  const bool bNoAdjacentMotion = !tmpPu || CU::isIntra(*tmpPu->cu);

  if (bNoAdjacentMotion)
  {
    return false;
  }

  mi = tmpPu->getMotionInfo(posNeighborMotion);
  const MotionInfo currMotion = pu.getMotionInfo(posSubBlock);

  if (mi.interDir)
  {
    if (mi.interDir != currMotion.interDir)
    {
      return true;
    }
    else
    {
      for (uint32_t iRefList = 0; iRefList < 2; iRefList++)
      {
        if (currMotion.interDir & (1 << iRefList))
        {
          Mv mvd = currMotion.mv[iRefList] - mi.mv[iRefList];

          if (!(mvd.getAbsHor() < 1 && mvd.getAbsVer() < 1 && currMotion.refIdx[iRefList] == mi.refIdx[iRefList]))
          {
            return true;
          }
        }
      }
      return false;
    }
  }
  else
  {
    return false;
  }
}
#endif

// TU tools

bool TU::isNonTransformedResidualRotated(const TransformUnit &tu, const ComponentID &compID)
{
  return tu.cs->sps->getSpsRangeExtension().getTransformSkipRotationEnabledFlag() && tu.blocks[compID].width == 4 && tu.cu->predMode == MODE_INTRA;
}

bool TU::getCbf( const TransformUnit &tu, const ComponentID &compID )
{
  return getCbfAtDepth( tu, compID, tu.depth );
}

bool TU::getCbfAtDepth(const TransformUnit &tu, const ComponentID &compID, const unsigned &depth)
{
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( !tu.blocks[compID].valid() )
  {
    CHECK(tu.cbf[compID] != 0, "cbf must be 0 if the component is not available");
  }
#endif
  return ((tu.cbf[compID] >> depth) & 1) == 1;
}

void TU::setCbfAtDepth(TransformUnit &tu, const ComponentID &compID, const unsigned &depth, const bool &cbf)
{
  // first clear the CBF at the depth
  tu.cbf[compID] &= ~(1  << depth);
  // then set the CBF
  tu.cbf[compID] |= ((cbf ? 1 : 0) << depth);
}

bool TU::isTSAllowed(const TransformUnit &tu, const ComponentID compID)
{
  const int maxSize = tu.cs->sps->getLog2MaxTransformSkipBlockSize();

  bool tsAllowed = tu.cs->sps->getTransformSkipEnabledFlag();
  tsAllowed &= ( !tu.cu->ispMode || !isLuma(compID) );
  SizeType transformSkipMaxSize = 1 << maxSize;
  tsAllowed &= !(tu.cu->bdpcmMode && isLuma(compID));
  tsAllowed &= !(tu.cu->bdpcmModeChroma && isChroma(compID));
  tsAllowed &= tu.blocks[compID].width <= transformSkipMaxSize && tu.blocks[compID].height <= transformSkipMaxSize;
  tsAllowed &= !tu.cu->sbtInfo;

  return tsAllowed;
}


int TU::getICTMode( const TransformUnit& tu, int jointCbCr )
{
  if( jointCbCr < 0 )
  {
    jointCbCr = tu.jointCbCr;
  }
  return g_ictModes[ tu.cs->picHeader->getJointCbCrSignFlag() ][ jointCbCr ];
}


bool TU::needsSqrt2Scale( const TransformUnit &tu, const ComponentID &compID )
{
  const Size &size=tu.blocks[compID];
  const bool isTransformSkip = (tu.mtsIdx[compID] == MTS_SKIP);
  return (!isTransformSkip) && (((floorLog2(size.width) + floorLog2(size.height)) & 1) == 1);
}

bool TU::needsBlockSizeTrafoScale( const TransformUnit &tu, const ComponentID &compID )
{
  return needsSqrt2Scale( tu, compID ) || isNonLog2BlockSize( tu.blocks[compID] );
}

TransformUnit* TU::getPrevTU( const TransformUnit &tu, const ComponentID compID )
{
  TransformUnit* prevTU = tu.prev;

  if( prevTU != nullptr && ( prevTU->cu != tu.cu || !prevTU->blocks[compID].valid() ) )
  {
    prevTU = nullptr;
  }

  return prevTU;
}

bool TU::getPrevTuCbfAtDepth( const TransformUnit &currentTu, const ComponentID compID, const int trDepth )
{
  const TransformUnit* prevTU = getPrevTU( currentTu, compID );
  return ( prevTU != nullptr ) ? TU::getCbfAtDepth( *prevTU, compID, trDepth ) : false;
}


// other tools

uint32_t getCtuAddr( const Position& pos, const PreCalcValues& pcv )
{
  return ( pos.x >> pcv.maxCUWidthLog2 ) + ( pos.y >> pcv.maxCUHeightLog2 ) * pcv.widthInCtus;
}

int getNumModesMip(const Size& block)
{
  switch( getMipSizeId(block) )
  {
  case 0: return 16;
  case 1: return  8;
  case 2: return  6;
  default: THROW( "Invalid mipSizeId" );
  }
}


int getMipSizeId(const Size& block)
{
  if( block.width == 4 && block.height == 4 )
  {
    return 0;
  }
  else if( block.width == 4 || block.height == 4 || (block.width == 8 && block.height == 8) )
  {
    return 1;
  }
  else
  {
    return 2;
  }

}

bool allowLfnstWithMip(const Size& block)
{
  if (block.width >= 16 && block.height >= 16)
  {
    return true;
  }
  return false;
}
#if JVET_V0130_INTRA_TMP
bool allowLfnstWithTmp()
{
	return true;
}
#endif
#if INTER_LIC
bool CU::isLICFlagPresent(const CodingUnit& cu)
{
  if (CU::isIntra(cu) || !cu.slice->getUseLIC())
  {
    return false;
  }

  if (cu.geoFlag || cu.firstPU->mergeFlag || cu.firstPU->interDir == 3 || CU::isIBC(cu) || cu.Y().area() < LIC_MIN_CU_PIXELS)
  {
    return false;
  }

#if JVET_Z0054_BLK_REF_PIC_REORDER
  if (!cu.cs->pcv->isEncoder)
  {
    if (PU::useRefCombList(*cu.firstPU))
    {
      CHECK(cu.firstPU->refIdxLC < 0, "cu.firstPU->refIdxLC < 0");
      return (cu.firstPU->refIdxLC >= cu.cs->slice->getNumNonScaledRefPic() ? false : true);
    }
  }
#endif
#if JVET_Y0128_NON_CTC
  if (!PU::checkRprLicCondition(*cu.firstPU))
  {
    return false;
  }
#endif

  return true;
}

void PU::spanLICFlags(PredictionUnit &pu, const bool LICFlag)
{
  MotionBuf mb = pu.getMotionBuf();
  MotionInfo *motionInfo = mb.buf;

  for (int y = 0; y < mb.height; y++)
  {
    for (int x = 0; x < mb.width; x++)
    {
      motionInfo[x].usesLIC = LICFlag;
    }
    motionInfo += mb.stride;
  }
}
#endif

#if SIGN_PREDICTION
bool TU::getDelayedSignCoding( const TransformUnit &tu, const ComponentID compID )
{
  const uint32_t maxSize = CU::isIntra( *tu.cu ) ? SIGN_PRED_MAX_BS_INTRA : SIGN_PRED_MAX_BS_INTER;

  const uint32_t width = tu.blocks[compID].width;
  const uint32_t height = tu.blocks[compID].height;
  if( tu.cs->sps->getNumPredSigns() <= 0 || tu.mtsIdx[compID] == MTS_SKIP )
  {
    return false;
  }

  if( width < 4 || height < 4 || width > maxSize || height > maxSize )
  {
    return false;
  }
  return true;
}

bool TU::getUseSignPred( const TransformUnit &tu, const ComponentID compID )
{
#if JVET_Y0141_SIGN_PRED_IMPROVE
#if JVET_V0130_INTRA_TMP
  return TU::getDelayedSignCoding( tu, compID ) && (!tu.cu->lfnstIdx || !tu.cu->tmpFlag);
#else
  return TU::getDelayedSignCoding( tu, compID ) && !tu.cu->lfnstIdx;
#endif
#else
  return TU::getDelayedSignCoding( tu, compID ) && !tu.cu->lfnstIdx;
#endif
}

void TU::predBorderResi( const Position blkPos, const CPelBuf &recoBuf, const CPelBuf &predBuf, const ComponentID compID,
  const uint32_t uiWidth, const uint32_t uiHeight, Pel *predResiBorder, const Pel defaultPel )
{
  Pel r1, r2, p0;

  const Pel *pReco1, *pReco2, *pPred0;

  const Pel *pReco = recoBuf.buf;
  const Pel *pPred = predBuf.buf;
  uint32_t strideReco = recoBuf.stride;
  uint32_t stridePred = predBuf.stride;


  if( blkPos.y )
  {
    pReco1 = pReco - strideReco;
    pReco2 = pReco - 2 * strideReco;
    pPred0 = pPred;
    for( int32_t x = 0; x < uiWidth; x++ )
    {
      p0 = pPred0[x];
      r1 = pReco1[x];
      r2 = pReco2[x];
      predResiBorder[uiHeight + x] = ( ( r1 << 1 ) - r2 - p0 );
    }
  }
  else
  {
    for( int32_t x = 0; x < uiWidth; x++ )
    {
      predResiBorder[uiHeight + x] = defaultPel - pPred[x];
    }
  }

  if( blkPos.x )
  {
    pReco1 = pReco - 1;
    pReco2 = pReco - 2;
    pPred0 = pPred;
    for( int32_t y = 0; y < uiHeight; y++ )
    {
      p0 = pPred0[0];
      r1 = pReco1[0];
      r2 = pReco2[0];
      predResiBorder[uiHeight - 1 - y] = ( ( r1 << 1 ) - r2 - p0 );
      pReco1 += strideReco;
      pReco2 += strideReco;
      pPred0 += stridePred;
    }
  }
  else
  {
    for( int32_t y = 0; y < uiHeight; y++ )
    {
      predResiBorder[uiHeight - 1 - y] = defaultPel - pPred[0];
      pPred += stridePred;
    }
  }

}

Position TU::posSignHidingFirstCG( const TransformUnit &tu, ComponentID compID )
{
  Position pos( -1, -1 );
  int width = tu.blocks[compID].width;
  int height = tu.blocks[compID].height;
  if( tu.cu->cs->slice->getSignDataHidingEnabledFlag() && width >= 4 && height >= 4 )
  {
    // Mayte use 4x3 size is enough
    const ScanElement * scan = g_scanOrder[SCAN_GROUPED_4x4][SCAN_DIAG][gp_sizeIdxInfo->idxFrom( width )][gp_sizeIdxInfo->idxFrom( height )];
    const int cgSize = 16;
    int   lastNZPosInCG = -1, firstNZPosInCG = cgSize;
    Position posHiding( -1, -1 );
    const CCoeffBuf &buf = tu.getCoeffs( compID );

    for( auto n = cgSize - 1; n >= 0; --n )
    {
      ScanElement currSE = scan[n];
      Position    currPos = Position( currSE.x, currSE.y );

      if( buf.at( currPos ) )
      {
        lastNZPosInCG = n;
        posHiding = currPos;
        break;
      }
    }

    if( lastNZPosInCG != -1 )
    {
      for( auto n = 0; n < cgSize; n++ )
      {
        ScanElement currSE = scan[n];
        Position    currPos = Position( currSE.x, currSE.y );
        if( buf.at( currPos ) )
        {
          firstNZPosInCG = n;
          break;
        }
      }

    }
    if( firstNZPosInCG - lastNZPosInCG >= SBH_THRESHOLD )
    {
      return posHiding;
    }
  }
  return Position( -1, -1 );
}
#endif

#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0128_NON_CTC
bool PU::checkIsValidMergeMvCand(const PredictionUnit &pu, int8_t mergeRefIdx[ NUM_REF_PIC_LIST_01 ])
{
  if (pu.amvpMergeModeFlag[REF_PIC_LIST_0] || pu.amvpMergeModeFlag[REF_PIC_LIST_1])
  {
    if (mergeRefIdx[pu.amvpMergeModeFlag[REF_PIC_LIST_0] ? 0 : 1] < 0)
    {
      return false;
    }
    return pu.cu->slice->getAmvpMergeModeValidCandPair(pu.amvpMergeModeFlag[REF_PIC_LIST_0] ? mergeRefIdx[0] : pu.refIdx[0],
                                                       pu.amvpMergeModeFlag[REF_PIC_LIST_0] ? pu.refIdx[1]   : mergeRefIdx[1]);
  }
  return true;
}
#else
bool PU::checkIsValidMergeMvCand(const CodingStructure &cs, const PredictionUnit &pu, const int curPoc, const int amvpPoc, int8_t mergeRefIdx[ NUM_REF_PIC_LIST_01 ])
{
  if (pu.amvpMergeModeFlag[REF_PIC_LIST_0])
  {
    if (mergeRefIdx[REF_PIC_LIST_0] < 0)
    {
      return false;
    }
    const int mergePoc = cs.slice->getRefPOC(REF_PIC_LIST_0, mergeRefIdx[REF_PIC_LIST_0]);
    if ((amvpPoc - curPoc) * (mergePoc -curPoc) > 0)
    {
      return false;
    }
  }
  if (pu.amvpMergeModeFlag[REF_PIC_LIST_1])
  {
    if (mergeRefIdx[REF_PIC_LIST_1] < 0)
    {
      return false;
    }
    const int mergePoc = cs.slice->getRefPOC(REF_PIC_LIST_1, mergeRefIdx[REF_PIC_LIST_1]);
    if ((amvpPoc - curPoc) * (mergePoc -curPoc) > 0)
    {
      return false;
    }
  }
  return true;
}
#endif
#endif

#if JVET_Z0135_TEMP_CABAC_WIN_WEIGHT
bool storeContexts( const Slice* slice, const int ctuXPosInCtus, const int ctuYPosInCtus )
{
  if( slice->getSPS()->getTempCabacInitMode() && !slice->isIntra() )
  {
    const PreCalcValues&  pcv = *slice->getPPS()->pcv;
    const int             ctuRsAddr = ctuXPosInCtus + ctuYPosInCtus * pcv.widthInCtus;
    return ctuRsAddr == pcv.sizeInCtus - 1;
  }
  return false;
}
#endif
