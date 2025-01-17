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

/** \file     CABACReader.h
 *  \brief    Reader for low level syntax
 */

#ifndef __CABACREADER__
#define __CABACREADER__

#include "BinDecoder.h"

#include "CommonLib/ContextModelling.h"
#include "CommonLib/MotionInfo.h"
#include "CommonLib/UnitPartitioner.h"


class CABACReader
{
public:
#if JVET_Z0135_TEMP_CABAC_WIN_WEIGHT
  CABACReader( BinDecoderBase& binDecoder, CABACDataStore *cabacDataStore ) : m_BinDecoder( binDecoder ), m_Bitstream( 0 )
  {
    m_CABACDataStore = cabacDataStore;
  }
#else
  CABACReader(BinDecoderBase& binDecoder) : m_BinDecoder(binDecoder), m_Bitstream(0) {}
#endif
  virtual ~CABACReader() {}

public:
  void        initCtxModels             ( Slice&                        slice );
  void        initBitstream             ( InputBitstream*               bitstream )           { m_Bitstream = bitstream; m_BinDecoder.init( m_Bitstream ); }
  const Ctx&  getCtx                    ()                                            const   { return m_BinDecoder.getCtx();  }
  Ctx&        getCtx                    ()                                                    { return m_BinDecoder.getCtx();  }

public:
  // slice segment data (clause 7.3.8.1)
  bool        terminating_bit           ();
  void        remaining_bytes           ( bool                          noTrailingBytesExpected );

  // coding tree unit (clause 7.3.8.2)
  void        coding_tree_unit          ( CodingStructure&              cs,     const UnitArea& area,     int (&qps)[2],   unsigned  ctuRsAddr );

  // sao (clause 7.3.8.3)
  void        sao                       ( CodingStructure&              cs,     unsigned        ctuRsAddr );

#if JVET_V0094_BILATERAL_FILTER
  void        bif                      (CodingStructure&              cs);
  void        bif                      (CodingStructure&              cs, unsigned ctuRsAddr);
#endif
#if JVET_X0071_CHROMA_BILATERAL_FILTER
  void        chromaBifCb              (CodingStructure&              cs);
  void        chromaBifCb              (CodingStructure&              cs, unsigned ctuRsAddr);
  void        chromaBifCr              (CodingStructure&              cs);
  void        chromaBifCr              (CodingStructure&              cs, unsigned ctuRsAddr);
#endif
  
#if JVET_W0066_CCSAO
  void        ccSaoControlIdc           ( CodingStructure &cs, const ComponentID compID, const int curIdx, uint8_t *controlIdc, Position lumaPos, int setNum );
#endif
  
  void        readAlfCtuFilterIndex(CodingStructure&              cs, unsigned        ctuRsAddr);

  void ccAlfFilterControlIdc(CodingStructure &cs, const ComponentID compID, const int curIdx, uint8_t *filterControlIdc,
                             Position lumaPos, int filterCount);

  // coding (quad)tree (clause 7.3.8.4)
  void        coding_tree               ( CodingStructure&              cs,     Partitioner&    pm,       CUCtx& cuCtx, Partitioner* pPartitionerChroma = nullptr, CUCtx* pCuCtxChroma = nullptr);
  PartSplit   split_cu_mode             ( CodingStructure&              cs,     Partitioner&    pm );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS  
  ModeType    mode_constraint           ( CodingStructure&              cs,     Partitioner&    pm,       const PartSplit splitMode );
#endif
  // coding unit (clause 7.3.8.5)
  void        coding_unit               ( CodingUnit&                   cu,     Partitioner&    pm,       CUCtx& cuCtx );
  void        cu_skip_flag              ( CodingUnit&                   cu );
  void        pred_mode                 ( CodingUnit&                   cu );
  void        bdpcm_mode                ( CodingUnit&                   cu,     const ComponentID compID );
  void        cu_pred_data              ( CodingUnit&                   cu );
#if ENABLE_OBMC
  void        obmc_flag                 ( CodingUnit&                   cu );
#endif
  void        cu_bcw_flag               ( CodingUnit&                   cu );
  void        extend_ref_line           (CodingUnit&                     cu);
  void        intra_luma_pred_modes     ( CodingUnit&                   cu );
#if JVET_W0123_TIMD_FUSION
  void        cu_timd_flag              ( CodingUnit&                   cu );
#endif
  void        intra_chroma_pred_modes   ( CodingUnit&                   cu );
  bool        intra_chroma_lmc_mode     ( PredictionUnit&               pu );
  void        intra_chroma_pred_mode    ( PredictionUnit&               pu );
#if ENABLE_DIMD
  void        cu_dimd_flag              (CodingUnit&                   cu);
 #endif
  void        cu_residual               ( CodingUnit&                   cu,     Partitioner&    pm,       CUCtx& cuCtx );
  void        rqt_root_cbf              ( CodingUnit&                   cu );
  void        adaptive_color_transform(CodingUnit&             cu);
  void        sbt_mode                  ( CodingUnit&                   cu );
  void        end_of_ctu                ( CodingUnit&                   cu,     CUCtx&          cuCtx );
#if JVET_V0130_INTRA_TMP
  void        tmp_flag                  ( CodingUnit&                   cu );
#endif
#if JVET_X0049_ADAPT_DMVR
  void        bm_merge_flag             ( PredictionUnit&               pu );
#endif
  void        mip_flag                  ( CodingUnit&                   cu );
  void        mip_pred_modes            ( CodingUnit&                   cu );
  void        mip_pred_mode             ( PredictionUnit&               pu );
  void        cu_palette_info           ( CodingUnit&                   cu,     ComponentID     compBegin, uint32_t numComp, CUCtx& cuCtx );
  void        cuPaletteSubblockInfo     ( CodingUnit&                   cu,     ComponentID     compBegin, uint32_t numComp, int subSetId, uint32_t& prevRunPos, unsigned& prevRunType );
  // prediction unit (clause 7.3.8.6)
  void        prediction_unit           ( PredictionUnit&               pu,     MergeCtx&       mrgCtx );
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  void        mvsd_data                 ( PredictionUnit&               pu );
#endif
  void        merge_flag                ( PredictionUnit&               pu );
  void        merge_data                ( PredictionUnit&               pu );
  void        affine_flag               ( CodingUnit&                   cu );
  void        subblock_merge_flag       ( CodingUnit&                   cu );
  void        merge_idx                 ( PredictionUnit&               pu );
  void        mmvd_merge_idx(PredictionUnit&               pu);
#if AFFINE_MMVD
  void        affine_mmvd_data          ( PredictionUnit&               pu );
#endif
#if TM_MRG
  void        tm_merge_flag             ( PredictionUnit&               pu );
#endif
#if JVET_W0097_GPM_MMVD_TM
  void        geo_mmvd_idx(PredictionUnit&          pu, RefPicList eRefPicList);
  void        geo_merge_idx(PredictionUnit&          pu);
#if JVET_Y0065_GPM_INTRA
  void        geo_merge_idx1            ( PredictionUnit&          pu, bool isIntra0, bool isIntra1);
#else
  void        geo_merge_idx1(PredictionUnit&          pu);
#endif
#endif
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
  void        geoModeIdx                ( PredictionUnit&          pu );
#endif
  void        imv_mode                  ( CodingUnit&                   cu,     MergeCtx&       mrgCtx );
  void        affine_amvr_mode          ( CodingUnit&                   cu,     MergeCtx&       mrgCtx );
  void        inter_pred_idc            ( PredictionUnit&               pu );
  void        ref_idx                   ( PredictionUnit&               pu,     RefPicList      eRefList );
#if JVET_Z0054_BLK_REF_PIC_REORDER
  void        refIdxLC                ( PredictionUnit&               pu );
  void        refPairIdx              ( PredictionUnit&               pu );
#endif
  void        mvp_flag                  ( PredictionUnit&               pu,     RefPicList      eRefList );
  void        Ciip_flag              ( PredictionUnit&               pu );
  void        smvd_mode              ( PredictionUnit&               pu );
#if MULTI_HYP_PRED
  int         ref_idx_mh(const int                     numRef);
  void        mh_pred_data(PredictionUnit&               pu);
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE
  void        amvpMerge_mode         ( PredictionUnit&               pu );
#endif
#if JVET_Z0050_CCLM_SLOPE
  void        cclmDelta             ( PredictionUnit&               pu, int8_t &delta );
  void        cclmDeltaSlope       ( PredictionUnit&               pu );
#endif


  // transform tree (clause 7.3.8.8)
  void        transform_tree            ( CodingStructure&              cs, Partitioner&    pm, CUCtx& cuCtx, const PartSplit ispType = TU_NO_ISP, const int subTuIdx = -1 );
  bool        cbf_comp                  ( CodingStructure&              cs,     const CompArea& area,     unsigned depth, const bool prevCbf = false, const bool useISP = false );

  // mvd coding (clause 7.3.8.9)
  void        mvd_coding                ( Mv &rMvd 
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
    , bool codeSign = true
#endif
  );
#if JVET_Z0131_IBC_BVD_BINARIZATION
  void        bvdCoding                ( Mv &rMvd );
  unsigned    xReadBvdContext(unsigned ctxT, int offset, int param);
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
  void        mvsdIdxFunc(PredictionUnit &pu, RefPicList eRefList);
  void        mvsdAffineIdxFunc(PredictionUnit &pu, RefPicList eRefList);
#endif

  // transform unit (clause 7.3.8.10)
  void        transform_unit            ( TransformUnit&                tu,     CUCtx&          cuCtx, Partitioner& pm,        const int subTuCounter = -1 );
  void        cu_qp_delta               ( CodingUnit&                   cu,     int             predQP, int8_t& qp );
  void        cu_chroma_qp_offset       ( CodingUnit&                   cu );

  // residual coding (clause 7.3.8.11)
  void        residual_coding           ( TransformUnit&                tu,     ComponentID     compID, CUCtx& cuCtx );
  void        ts_flag                   ( TransformUnit&                tu,     ComponentID     compID );
  void        mts_idx                   ( CodingUnit&                   cu,     CUCtx&          cuCtx  );
  void        residual_lfnst_mode       ( CodingUnit&                   cu,     CUCtx&          cuCtx  );
  void        isp_mode                  ( CodingUnit&                   cu );
  int         last_sig_coeff            ( CoeffCodingContext&           cctx,   TransformUnit& tu, ComponentID   compID );
#if TCQ_8STATES
#if SIGN_PREDICTION
  void        residual_coding_subblock  ( CoeffCodingContext&           cctx,   TCoeff*         coeff,   TCoeff*         sign, const uint64_t stateTransTable, int& state );
#else
  void        residual_coding_subblock  ( CoeffCodingContext&           cctx,   TCoeff*         coeff, const uint64_t stateTransTable, int& state );
#endif
#else
#if SIGN_PREDICTION
	void        residual_coding_subblock  ( CoeffCodingContext&           cctx,   TCoeff*         coeff,   TCoeff*         sign, const int stateTransTable, int& state );
#else
	void        residual_coding_subblock  ( CoeffCodingContext&           cctx,   TCoeff*         coeff, const int stateTransTable, int& state );
#endif
#endif
#if SIGN_PREDICTION
  void        parsePredictedSigns ( TransformUnit &tu, ComponentID compID);
#endif
	void        residual_codingTS         ( TransformUnit&                tu,     ComponentID     compID );
  void        residual_coding_subblockTS( CoeffCodingContext&           cctx,   TCoeff*         coeff  );
  void        joint_cb_cr               ( TransformUnit&                tu,     const int cbfMask );


private:
  unsigned    unary_max_symbol          ( unsigned ctxId0, unsigned ctxIdN, unsigned maxSymbol );
  unsigned    unary_max_eqprob          (                                   unsigned maxSymbol );
  unsigned    exp_golomb_eqprob         ( unsigned count );
  unsigned    get_num_bits_read         () { return m_BinDecoder.getNumBitsRead(); }
  unsigned    code_unary_fixed          ( unsigned ctxId, unsigned unary_max, unsigned fixed );

  void        xReadTruncBinCode(uint32_t& symbol, uint32_t maxSymbol);
  void        parseScanRotationModeFlag ( CodingUnit& cu,           ComponentID compBegin );
  void        xDecodePLTPredIndicator   ( CodingUnit& cu,           uint32_t maxPLTSize,   ComponentID compBegin );
  void        xAdjustPLTIndex           ( CodingUnit& cu,           Pel curLevel,          uint32_t idx, PelBuf& paletteIdx, PLTtypeBuf& paletteRunType, int maxSymbol, ComponentID compBegin );
public:
#if INTER_LIC
  void        cu_lic_flag               ( CodingUnit& cu );
#endif

#if JVET_Z0135_TEMP_CABAC_WIN_WEIGHT
  CABACDataStore*         m_CABACDataStore;
#endif

private:
  BinDecoderBase& m_BinDecoder;
  InputBitstream* m_Bitstream;
  ScanElement*    m_scanOrder;
};


class CABACDecoder
{
public:
#if JVET_Z0135_TEMP_CABAC_WIN_WEIGHT
  CABACDecoder()
    : m_CABACReaderStd  ( m_BinDecoderStd, nullptr )
    , m_CABACReader     { &m_CABACReaderStd }
    , m_CABACDataStore  ( nullptr )
  {
    m_CABACDataStore = new CABACDataStore;

    m_CABACReaderStd.m_CABACDataStore = m_CABACDataStore;

    for( int i = 0; i < BPM_NUM - 1; i++ )
    {
      m_CABACReader[i]->m_CABACDataStore = m_CABACDataStore;
    }
  }

  virtual ~CABACDecoder()
  {
    if( m_CABACDataStore )
    {
      delete m_CABACDataStore;
    }
  }
#else
  CABACDecoder()
    : m_CABACReaderStd( m_BinDecoderStd )
    , m_CABACReader{ &m_CABACReaderStd }
  {}
#endif

  CABACReader*                getCABACReader    ( int           id    )       { return m_CABACReader[id]; }

private:
  BinDecoder_Std          m_BinDecoderStd;
  CABACReader             m_CABACReaderStd;
  CABACReader*            m_CABACReader[BPM_NUM-1];

#if JVET_Z0135_TEMP_CABAC_WIN_WEIGHT
  CABACDataStore*         m_CABACDataStore;
#endif
};

#endif
