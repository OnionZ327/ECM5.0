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

/** \file     CommonDef.h
    \brief    Defines version information, constants and small in-line functions
*/

#ifndef __COMMONDEF__
#define __COMMONDEF__

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>


#if _MSC_VER > 1000
// disable "signed and unsigned mismatch"
#pragma warning( disable : 4018 )
// disable bool coercion "performance warning"
#pragma warning( disable : 4800 )
#endif // _MSC_VER > 1000
#include "TypeDef.h"
#include "version.h"

// MS Visual Studio before 2014 does not support required C++11 features
#ifdef _MSC_VER
#if _MSC_VER < 1900
#error "MS Visual Studio version not supported. Please upgrade to Visual Studio 2015 or higher (or use other compilers)"
#endif

#include <intrin.h>
#endif

//! \ingroup CommonLib
//! \{

// ====================================================================================================================
// Platform information
// ====================================================================================================================

#ifdef __clang__
#define NVM_COMPILEDBY  "[clang %d.%d.%d]", __clang_major__, __clang_minor__, __clang_patchlevel__
#ifdef __IA64__
#define NVM_ONARCH    "[on 64-bit] "
#else
#define NVM_ONARCH    "[on 32-bit] "
#endif
#elif __GNUC__
#define NVM_COMPILEDBY  "[GCC %d.%d.%d]", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__
#ifdef __IA64__
#define NVM_ONARCH    "[on 64-bit] "
#else
#define NVM_ONARCH    "[on 32-bit] "
#endif
#endif

#ifdef __INTEL_COMPILER
#define NVM_COMPILEDBY  "[ICC %d]", __INTEL_COMPILER
#elif  _MSC_VER
#define NVM_COMPILEDBY  "[VS %d]", _MSC_VER
#endif

#ifndef NVM_COMPILEDBY
#define NVM_COMPILEDBY "[Unk-CXX]"
#endif

#ifdef _WIN32
#define NVM_ONOS        "[Windows]"
#elif  __linux
#define NVM_ONOS        "[Linux]"
#elif  __CYGWIN__
#define NVM_ONOS        "[Cygwin]"
#elif __APPLE__
#define NVM_ONOS        "[Mac OS X]"
#else
#define NVM_ONOS "[Unk-OS]"
#endif

#define NVM_BITS          "[%d bit] ", (sizeof(void*) == 8 ? 64 : 32) ///< used for checking 64-bit O/S

#ifndef NULL
#define NULL              0
#endif

typedef enum
{
  AFFINEMODEL_4PARAM,
  AFFINEMODEL_6PARAM,
  AFFINE_MODEL_NUM
} EAffineModel;

static const int AFFINE_ME_LIST_SIZE =                             4;
static const int AFFINE_ME_LIST_SIZE_LD =                          3;
static const double AFFINE_ME_LIST_MVP_TH =                        1.0;

#if JVET_V0094_BILATERAL_FILTER || JVET_X0071_CHROMA_BILATERAL_FILTER
static const int32_t NUMBER_PADDED_SAMPLES =  2;
static const int32_t BIF_ROUND_ADD =         32;
static const int32_t BIF_ROUND_SHIFT =        6;
#endif

// ====================================================================================================================
// Common constants
// ====================================================================================================================
static const uint64_t   MAX_UINT64 =                  0xFFFFFFFFFFFFFFFFU;
static const uint32_t   MAX_UINT =                            0xFFFFFFFFU; ///< max. value of unsigned 32-bit integer
static const int    MAX_INT =                              2147483647; ///< max. value of signed 32-bit integer
static const uint8_t  MAX_UCHAR =                                   255;
static const uint8_t  MAX_SCHAR =                                   127;
static const double MAX_DOUBLE =                             1.7e+308; ///< max. value of double-type value

// ====================================================================================================================
// Coding tool configuration
// ====================================================================================================================
// Most of these should not be changed - they resolve the meaning of otherwise magic numbers.

static const int MAX_GOP =                                         64; ///< max. value of hierarchical GOP size
static const int MAX_NUM_REF_PICS =                                29; ///< max. number of pictures used for reference
static const int MAX_NUM_REF =                                     16; ///< max. number of entries in picture reference list
static const int MAX_QP =                                          63;
static const int NOT_VALID =                                       -1;

#if TM_AMVP
static const int REGULAR_AMVP_MAX_NUM_CANDS =                       5; ///< AMVP: advanced motion vector prediction - max number of final candidate for regular inter mode
#endif
static const int AMVP_MAX_NUM_CANDS =                               2; ///< AMVP: advanced motion vector prediction - max number of final candidates
static const int AMVP_MAX_NUM_CANDS_MEM =                           3; ///< AMVP: advanced motion vector prediction - max number of candidates
#if INTER_RM_SIZE_CONSTRAINTS
static const int AMVP_DECIMATION_FACTOR =                           1;
#else
static const int AMVP_DECIMATION_FACTOR =                           2;
#endif
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
static const int NUM_MERGE_CANDS =                                 30; ///< for maximum buffer of merging candidates
static const int NUM_TMVP_CANDS =                                   9; ///< TMVP
#elif JVET_Z0075_IBC_HMVP_ENLARGE
static const int NUM_MERGE_CANDS =                                 20;
#endif
#if NON_ADJACENT_MRG_CAND
static const int MRG_MAX_NUM_CANDS =                                15; ///< MERGE
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
static const int NUM_NON_ADJ_CANDS =                                18; ///< Non-Adj
#endif
static const int LAST_MERGE_IDX_CABAC =                              5;
#else
static const int MRG_MAX_NUM_CANDS =                                6; ///< MERGE
#endif

#if JVET_Z0139_NA_AFF
static const int AFF_NON_ADJACENT_DIST    =                          4;
#endif

#if JVET_Z0139_HIST_AFF || JVET_Z0139_NA_AFF
static const int AFFINE_MRG_MAX_NUM_CANDS =                         15; ///< AFFINE MERGE
#else
static const int AFFINE_MRG_MAX_NUM_CANDS =                         5; ///< AFFINE MERGE
#endif
static const int IBC_MRG_MAX_NUM_CANDS =                            6; ///< IBC MERGE
#if JVET_Z0075_IBC_HMVP_ENLARGE
static const int IBC_MRG_MAX_NUM_CANDS_MEM =                        20; ///< IBC MERGE- max number of candidates 
#endif
#if JVET_X0083_BM_AMVP_MERGE_MODE
#if JVET_Y0129_MVD_SIGNAL_AMVP_MERGE_MODE
static const int MAX_NUM_AMVP_CANDS_MAX_REF =                       MAX_NUM_REF * AMVP_MAX_NUM_CANDS_MEM;
#else
static const int MAX_NUM_AMVP_CANDS_MAX_REF =                       MAX_NUM_REF * AMVP_MAX_NUM_CANDS;
#endif
static const int AMVP_MERGE_MODE_MERGE_LIST_MAX_CANDS =             6;
static const int AMVP_MERGE_MODE_REDUCED_MV_REFINE_SEARCH_ROUND =   8;
#endif

static const int MAX_TLAYER =                                       7; ///< Explicit temporal layer QP offset - max number of temporal layer

static const int ADAPT_SR_SCALE =                                   1; ///< division factor for adaptive search range

static const int MIN_TB_LOG2_SIZEY = 2;
#if TU_256
static const int MAX_TB_LOG2_SIZEY =                                8;
#else
static const int MAX_TB_LOG2_SIZEY =                                6;
#endif

static const int MIN_TB_SIZEY = 1 << MIN_TB_LOG2_SIZEY;
static const int MAX_TB_SIZEY = 1 << MAX_TB_LOG2_SIZEY;

static const int MAX_NUM_PICS_IN_SOP =                           1024;

static const int MAX_NESTING_NUM_OPS =                           1024;
static const int MAX_NESTING_NUM_LAYER =                           64;

static const int MAX_VPS_NUM_HRD_PARAMETERS =                       1;
static const int MAX_VPS_LAYERS =                                  64;
static const int MAX_VPS_SUBLAYERS =                                7;
static const int MAX_NUM_REF_LAYERS =                               7;
static const int MAX_NUM_OLSS =                                   256;
static const int MAX_VPS_OLS_MODE_IDC =                             2;
static const int MAXIMUM_INTRA_FILTERED_WIDTH =                    16;
static const int MAXIMUM_INTRA_FILTERED_HEIGHT =                   16;

#if TU_256
static const int MAX_INTRA_SIZE =                                  128;
static const int MIP_MAX_WIDTH =                                   MAX_INTRA_SIZE;
static const int MIP_MAX_HEIGHT =                                  MAX_INTRA_SIZE;
#else
static const int MIP_MAX_WIDTH =                                   MAX_TB_SIZEY;
static const int MIP_MAX_HEIGHT =                                  MAX_TB_SIZEY;
#endif

#if JVET_W0066_CCSAO
#define          MAX_CCSAO_SET_NUM                                  4
static const int MAX_CCSAO_CAND_POS_Y        =                      9;
static const int MAX_CCSAO_CAND_POS_Y_BITS   =                      4;
static const int MAX_CCSAO_BAND_NUM_Y        =                     16;
static const int MAX_CCSAO_BAND_NUM_Y_BITS   =                      4;
static const int MAX_CCSAO_BAND_NUM_U        =                      4;
static const int MAX_CCSAO_BAND_NUM_U_BITS   =                      2;
static const int MAX_CCSAO_BAND_NUM_V        =                      4;
static const int MAX_CCSAO_BAND_NUM_V_BITS   =                      2;
static const int MAX_CCSAO_CLASS_NUM         =                     64;
static const int MAX_CCSAO_OFFSET_THR        =                     15;
static const int MAX_CCSAO_FILTER_LENGTH     =                      3;
#endif

#if JVET_Y0106_CCSAO_EDGE_CLASSIFIER
static const int N_C                            = 3; /* Num components*/
static const int Y_C                            = 0; /* Y luma */
static const int U_C                            = 1; /* Cb Chroma */
static const int V_C                            = 2; /* Cr Chroma */
static const int CCSAO_EDGE_NUM                 = 16;
static const int CCSAO_EDGE_COMPARE_VALUE       = 2;
static const int CCSAO_EDGE_TYPE                = 4;
static const int CCSAO_QUAN_NUM                 = 16;
static const int CCSAO_EDGE_BAND_NUM_Y          = 4;
static const int CCSAO_EDGE_BAND_NUM_C          = 4;
static const int MAX_CCSAO_BAND_NUM_U_BAND_BITS = 4;
#endif

static const int MAX_NUM_ALF_ALTERNATIVES_CHROMA =                  8;
static const int MAX_NUM_ALF_CLASSES         =                     25;
#if ALF_IMPROVEMENT
static const int MAX_NUM_ALF_ALTERNATIVES_LUMA = 4;
static const int EXT_LENGTH = 2;
static const int MAX_NUM_ALF_LUMA_COEFF      =                     21 + EXT_LENGTH;
static const int MAX_NUM_ALF_CHROMA_COEFF    =                     21;
static const int MAX_ALF_FILTER_LENGTH       =                      9;
#else
static const int MAX_NUM_ALF_LUMA_COEFF      =                     13;
static const int MAX_NUM_ALF_CHROMA_COEFF    =                      7;
static const int MAX_ALF_FILTER_LENGTH       =                      7;
#endif
static const int MAX_NUM_ALF_COEFF           =                     MAX_ALF_FILTER_LENGTH * MAX_ALF_FILTER_LENGTH / 2 + 1;
#if JVET_Z0105_LOOP_FILTER_VIRTUAL_BOUNDARY
static const int       MAX_ALF_PADDING_SIZE           =             8;
#else
static const int       MAX_ALF_PADDING_SIZE           =             4;
#endif
#if JVET_X0071_LONGER_CCALF
#define MAX_NUM_CC_ALF_FILTERS                                      16
static constexpr int MAX_NUM_CC_ALF_CHROMA_COEFF    =               25;
#else
#define MAX_NUM_CC_ALF_FILTERS                                      4
static constexpr int MAX_NUM_CC_ALF_CHROMA_COEFF    =               8;
#endif
static constexpr int CCALF_DYNAMIC_RANGE            =               6;
static constexpr int CCALF_BITS_PER_COEFF_LEVEL     =               3;

static const int ALF_FIXED_FILTER_NUM        =                     64;
static const int ALF_CTB_MAX_NUM_APS         =                      8;
#if ALF_IMPROVEMENT
static const int ALF_ORDER                   =                      4;
static const int NUM_FIXED_FILTER_SETS       =                      2;
#if JVET_X0071_ALF_BAND_CLASSIFIER
static const int ALF_NUM_CLASSIFIER          =                      2;
static const int ALF_CLASSES_NEW             =                     25;
static const int ALF_NUM_CLASSES_CLASSIFIER[ALF_NUM_CLASSIFIER] = { MAX_NUM_ALF_CLASSES,  ALF_CLASSES_NEW };
#endif
#else 
static const int NUM_FIXED_FILTER_SETS       =                     16;
static const int NUM_TOTAL_FILTER_SETS       =                     NUM_FIXED_FILTER_SETS + ALF_CTB_MAX_NUM_APS;
#endif

#if JVET_Z0139_HIST_AFF
static const int MAX_NUM_AFFHMVP_ENTRIES_ONELIST = 5;
static const int MAX_NUM_AFFHMVP_ENTRIES = MAX_NUM_AFFHMVP_ENTRIES_ONELIST * 2;
#endif

#if !BDOF_RM_CONSTRAINTS
static const int MAX_BDOF_APPLICATION_REGION =                     16;
#endif

static const int MAX_CPB_CNT =                                     32; ///< Upper bound of (cpb_cnt_minus1 + 1)
static const int MAX_NUM_LAYER_IDS =                               64;
static const int COEF_REMAIN_BIN_REDUCTION =                        5; ///< indicates the level at which the VLC transitions from Golomb-Rice to TU+EG(k)
#if !JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
static const int COEFF_MIN =                                   -32768;
static const int COEFF_MAX =                                    32767;
#endif
static const int CU_DQP_TU_CMAX =                                   5; ///< max number bins for truncated unary
static const int CU_DQP_EG_k =                                      0; ///< expgolomb order

static const int SBH_THRESHOLD =                                    4; ///< value of the fixed SBH controlling threshold

static const int MAX_NUM_VPS =                                     16;
static const int MAX_NUM_DPS =                                     16;
static const int MAX_NUM_SPS =                                     16;
static const int MAX_NUM_PPS =                                     64;
static const int MAX_NUM_APS =                                     32;  //Currently APS ID has 5 bits
static const int NUM_APS_TYPE_LEN =                                 3;  //Currently APS Type has 3 bits
static const int MAX_NUM_APS_TYPE =                                 8;  //Currently APS Type has 3 bits so the max type is 8

static const int MAX_TILE_COLS =                                   20;  ///< Maximum number of tile columns
#if JVET_S0156_LEVEL_DEFINITION
static const int MAX_TILES =                                      440;  ///< Maximum number of tiles
#else
static const int MAX_TILE_ROWS =                                   22;  ///< Maximum number of tile rows
static const int MAX_TILES =            MAX_TILE_COLS * MAX_TILE_ROWS;  ///< Maximum number of tiles
#endif
static const int MAX_SLICES =                                     600;  ///< Maximum number of slices per picture
#if TU_256
static const int MLS_GRP_NUM =       MAX_TB_SIZEY * MAX_TB_SIZEY >> 4;  ///< Max number of coefficient groups
#else
static const int MLS_GRP_NUM =                                   1024; ///< Max number of coefficient groups, max(16, 256)
#endif

static const int MLS_CG_SIZE =                                      4; ///< Coefficient group size of 4x4; = MLS_CG_LOG2_WIDTH + MLS_CG_LOG2_HEIGHT


static const int RVM_VCEGAM10_M =                                   4;

#if JVET_Y0116_EXTENDED_MRL_LIST
static const int MAX_REF_LINE_IDX =                                 13; //highest refLine offset in the list
static const int MRL_NUM_REF_LINES =                                6; //number of candidates in the array
static const int MULTI_REF_LINE_IDX[6] =               { 0, 1, 3, 5, 7, 12 };
#else
static const int MAX_REF_LINE_IDX =                                 3; //highest refLine offset in the list
static const int MRL_NUM_REF_LINES =                                3; //number of candidates in the array
static const int MULTI_REF_LINE_IDX[4] =               { 0, 1, 2, 0 };
#endif

static const int PRED_REG_MIN_WIDTH =                               4;  // Minimum prediction region width for ISP subblocks

static const int NUM_LUMA_MODE =                                   67; ///< Planar + DC + 65 directional mode (4*16 + 1)
#if MMLM
static const int NUM_LMC_MODE = 1 + 2 + 3; ///< LMC + MDLM_T + MDLM_L + MMLM + MMLM_L + MMLM_T
#else
static const int NUM_LMC_MODE =                                    1 + 2; ///< LMC + MDLM_T + MDLM_L
#endif
static const int NUM_INTRA_MODE = (NUM_LUMA_MODE + NUM_LMC_MODE);

#if TU_256
static const int NUM_EXT_LUMA_MODE =                               30;
#else
static const int NUM_EXT_LUMA_MODE =                               28;
#endif

static const int NUM_DIR =           (((NUM_LUMA_MODE - 3) >> 2) + 1);
static const int PLANAR_IDX =                                       0; ///< index for intra PLANAR mode
static const int DC_IDX =                                           1; ///< index for intra DC     mode
static const int HOR_IDX =                    (1 * (NUM_DIR - 1) + 2); ///< index for intra HORIZONTAL mode
static const int DIA_IDX =                    (2 * (NUM_DIR - 1) + 2); ///< index for intra DIAGONAL   mode
static const int VER_IDX =                    (3 * (NUM_DIR - 1) + 2); ///< index for intra VERTICAL   mode
static const int VDIA_IDX =                   (4 * (NUM_DIR - 1) + 2); ///< index for intra VDIAGONAL  mode
#if JVET_W0123_TIMD_FUSION
static const int BDPCM_IDX =                                      162;
#else
static const int BDPCM_IDX =                  (5 * (NUM_DIR - 1) + 2); ///< index for intra VDIAGONAL  mode
#endif
static const int NOMODE_IDX =                               MAX_UCHAR; ///< indicating uninitialized elements
#if JVET_Z0050_DIMD_CHROMA_FUSION && ENABLE_DIMD
static const int NUM_CHROMA_MODE = (6 + NUM_LMC_MODE); ///< total number of chroma modes
#else
static const int NUM_CHROMA_MODE = (5 + NUM_LMC_MODE); ///< total number of chroma modes
#endif
static const int LM_CHROMA_IDX = NUM_LUMA_MODE; ///< chroma mode index for derived from LM mode
#if ENABLE_DIMD
static const int DIMD_IDX =                                        99; ///< index for intra DIMD mode
#endif
#if JVET_W0123_TIMD_FUSION
static const int TIMD_IDX =                                       199; ///< index for intra TIMD mode
static const int DIMD_MAX_TEMP_SIZE =                               4;
static const int EXT_HOR_IDX =                                     34;
static const int EXT_DIA_IDX =                                     66;
static const int EXT_VER_IDX =                                     98;
static const int EXT_VDIA_IDX =                                   130;
#define MAP131TO67( mode )                 (mode<2?mode:((mode>>1)+1))
#define MAP67TO131( mode )                 (mode<2?mode:((mode<<1)-2))
#endif
#if MMLM
static const int MMLM_CHROMA_IDX = LM_CHROMA_IDX + 1; ///< MDLM_L
static const int MDLM_L_IDX = LM_CHROMA_IDX + 2; ///< MDLM_L
static const int MMLM_L_IDX = LM_CHROMA_IDX + 3; ///< MDLM_T
static const int MDLM_T_IDX = LM_CHROMA_IDX + 4; ///< MDLM_L
static const int MMLM_T_IDX = LM_CHROMA_IDX + 5; ///< MDLM_T
#else
static const int MDLM_L_IDX =                          LM_CHROMA_IDX + 1; ///< MDLM_L
static const int MDLM_T_IDX =                          LM_CHROMA_IDX + 2; ///< MDLM_T
#endif
#if JVET_Z0050_DIMD_CHROMA_FUSION && ENABLE_DIMD
static const int DIMD_CHROMA_IDX =                     NUM_INTRA_MODE; ///< chroma mode index for derived by DIMD method
static const int DM_CHROMA_IDX =                       NUM_INTRA_MODE + 1; ///< chroma mode index for derived from luma intra mode
#else
static const int DM_CHROMA_IDX =                       NUM_INTRA_MODE; ///< chroma mode index for derived from luma intra mode
#endif

#if JVET_Z0131_IBC_BVD_BINARIZATION
static const int BVD_CODING_GOLOMB_ORDER = 1; 
static const int NUM_HOR_BVD_CTX         = 5;
static const int NUM_VER_BVD_CTX         = 5;
static const int HOR_BVD_CTX_OFFSET      = 0;
static const int BVD_IBC_MAX_PREFIX      = 16;
static const int VER_BVD_CTX_OFFSET      = 6;
#endif

#if JVET_W0119_LFNST_EXTENSION
static const int NUM_LFNST_INTRA_MODES   = NUM_LUMA_MODE + NUM_EXT_LUMA_MODE;
static const uint32_t  L16W_ZO           = 96;
static const uint32_t  L16W              = 96;
static const uint32_t  L16H              = 32;
static const uint32_t  L16H_ZO           = 32;
static const uint32_t  L8W               = 64;
static const uint32_t  L8W_ZO            = 64;
static const uint32_t  L8H               = 32;
static const uint32_t  L8H_ZO            = 32;
#endif

#if JVET_Y0142_ADAPT_INTRA_MTS
static const uint32_t  MTS_TH_COEFF[2] =                                { 6, 32};
static const uint32_t  MTS_NCANDS[3] =                                  { 1, 4, 6};
static const uint32_t  NUM_TRAFO_MODES_MTS =                            8;
#else
static const uint32_t  NUM_TRAFO_MODES_MTS =                            6; ///< Max Intra CU size applying EMT, supported values: 8, 16, 32, 64, 128
#endif
#if TU_256
static const uint32_t  MTS_INTRA_MAX_CU_SIZE =                         256; ///< Max Intra CU size applying EMT, supported values: 8, 16, 32, 64, 128
static const uint32_t  MTS_INTER_MAX_CU_SIZE =                         256; ///< Max Inter CU size applying EMT, supported values: 8, 16, 32, 64, 128
#else
static const uint32_t  MTS_INTRA_MAX_CU_SIZE =                         32; ///< Max Intra CU size applying EMT, supported values: 8, 16, 32, 64, 128
static const uint32_t  MTS_INTER_MAX_CU_SIZE =                         32; ///< Max Inter CU size applying EMT, supported values: 8, 16, 32, 64, 128
#endif
#if SECONDARY_MPM
static const int NUM_PRIMARY_MOST_PROBABLE_MODES = 6;
static const int NUM_SECONDARY_MOST_PROBABLE_MODES = 16;
static const int NUM_MOST_PROBABLE_MODES = NUM_PRIMARY_MOST_PROBABLE_MODES + NUM_SECONDARY_MOST_PROBABLE_MODES;
#if SECONDARY_MPM
static const int NUM_NON_MPM_MODES = NUM_LUMA_MODE - NUM_MOST_PROBABLE_MODES;
#endif
#else
static const int NUM_MOST_PROBABLE_MODES = 6;
#endif
static const int LM_SYMBOL_NUM = (1 + NUM_LMC_MODE);

static const int MAX_NUM_MIP_MODE =                                32; ///< maximum number of MIP pred. modes
static const int FAST_UDI_MAX_RDMODE_NUM = (NUM_LUMA_MODE + MAX_NUM_MIP_MODE); ///< maximum number of RD comparison in fast-UDI estimation loop

static const int MAX_LFNST_COEF_NUM =                              16;

static const int LFNST_LAST_SIG_LUMA =                              1;
static const int LFNST_LAST_SIG_CHROMA =                            1;

#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
static const int NUM_LFNST_NUM_PER_SET =                            4;
#else
static const int NUM_LFNST_NUM_PER_SET =                            3;
#endif

static const int LOG2_MAX_NUM_COLUMNS_MINUS1 =                      7;
static const int LOG2_MAX_NUM_ROWS_MINUS1 =                         7;

static const int CABAC_INIT_PRESENT_FLAG =                          1;

static const int MV_FRACTIONAL_BITS_INTERNAL                      = 4;
static const int MV_FRACTIONAL_BITS_SIGNAL                        = 2;
static const int MV_FRACTIONAL_BITS_DIFF = MV_FRACTIONAL_BITS_INTERNAL - MV_FRACTIONAL_BITS_SIGNAL;
static const int LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL = 1 << MV_FRACTIONAL_BITS_SIGNAL;
static const int LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS = 1 << MV_FRACTIONAL_BITS_INTERNAL;
static const int CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS = 1 << (MV_FRACTIONAL_BITS_INTERNAL + 1);

static const int MAX_NUM_SUB_PICS =                         (1 << 16);
static const int MAX_NUM_LONG_TERM_REF_PICS =                      33;
static const int NUM_LONG_TERM_REF_PIC_SPS =                        0;


static const int MAX_QP_OFFSET_LIST_SIZE =                          6; ///< Maximum size of QP offset list is 6 entries
static const int MAX_NUM_CQP_MAPPING_TABLES =                       3; ///< Maximum number of chroma QP mapping tables (Cb, Cr and joint Cb-Cr)
static const int MIN_QP_VALUE_FOR_16_BIT   =                      -48; ////< Minimum value for QP (-6*(bitdepth - 8) ) for bit depth 16 ; actual minimum QP value is bit depth dependent
static const int MAX_NUM_QP_VALUES =    MAX_QP + 1 - MIN_QP_VALUE_FOR_16_BIT; ////< Maximum number of QP values possible - bit depth dependent

// Cost mode support
static const int LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP =      0; ///< QP to use for lossless coding.
static const int LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME =4; ///< QP' to use for mixed_lossy_lossless coding.

static const int RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS =     4;

static const int RExt__PREDICTION_WEIGHTING_ANALYSIS_DC_PRECISION = 0; ///< Additional fixed bit precision used during encoder-side weighting prediction analysis. Currently only used when high_precision_prediction_weighting_flag is set, for backwards compatibility reasons.

static const int MAX_TIMECODE_SEI_SETS =                            3; ///< Maximum number of time sets

#if CTU_256
static const int MAX_CU_DEPTH =                                     8; ///< log2(CTUSize)
#else
static const int MAX_CU_DEPTH =                                     7; ///< log2(CTUSize)
#endif
static const int MAX_CU_SIZE =                        1<<MAX_CU_DEPTH;
static const int MIN_CU_LOG2 =                                      2;
static const int MIN_PU_SIZE =                                      4;
static const int MAX_NUM_PARTS_IN_CTU =                         ( ( MAX_CU_SIZE * MAX_CU_SIZE ) >> ( MIN_CU_LOG2 << 1 ) );
#if !CONVERT_NUM_TU_SPLITS_TO_CFG
static const int MAX_NUM_TUS =                                     16; ///< Maximum number of TUs within one CU. When max TB size is 32x32, up to 16 TUs within one CU (128x128) is supported
#endif
static const int MAX_LOG2_DIFF_CU_TR_SIZE =                         3;
static const int MAX_CU_TILING_PARTITIONS = 1 << ( MAX_LOG2_DIFF_CU_TR_SIZE << 1 );
#if TU_256
static const int JVET_C0024_ZERO_OUT_TH =                          256;
#else
static const int JVET_C0024_ZERO_OUT_TH =                          32;
#endif

static const int MAX_NUM_PART_IDXS_IN_CTU_WIDTH = MAX_CU_SIZE/MIN_PU_SIZE; ///< maximum number of partition indices across the width of a CTU (or height of a CTU)
static const int SCALING_LIST_REM_NUM =                             6;

static const int QUANT_SHIFT =                                     14; ///< Q(4) = 2^14
static const int IQUANT_SHIFT =                                     6;

static constexpr int    SCALE_BITS      = 15;   // Precision for fractional bit estimates
static constexpr double FRAC_BITS_SCALE = 1.0 / (1 << SCALE_BITS);

static constexpr int SCALING_LIST_PRED_MODES = 2;
static const int SCALING_LIST_NUM = MAX_NUM_COMPONENT * SCALING_LIST_PRED_MODES; ///< list number for quantization matrix

static const int SCALING_LIST_START_VALUE =                         8; ///< start value for dpcm mode
static const int MAX_MATRIX_COEF_NUM =                             64; ///< max coefficient number for quantization matrix
static const int MAX_MATRIX_SIZE_NUM =                              8; ///< max size number for quantization matrix
static const int SCALING_LIST_BITS =                                8; ///< bit depth of scaling list entries
static const int LOG2_SCALING_LIST_NEUTRAL_VALUE =                  4; ///< log2 of the value that, when used in a scaling list, has no effect on quantisation
static const int SCALING_LIST_DC =                                 16; ///< default DC value

#if TU_256
static const int LAST_SIGNIFICANT_GROUPS =                         16;
#else
static const int LAST_SIGNIFICANT_GROUPS =                         14;
#endif

static const int AFFINE_MIN_BLOCK_SIZE =                            4; ///< Minimum affine MC block size

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int MMVD_REFINE_STEP =                                 6; ///< max number of distance step
#else
static const int MMVD_REFINE_STEP =                                 8; ///< max number of distance step
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int MMVD_MAX_DIR =                                     16;
#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int MMVD_SIZE_SHIFT =                                   3;
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int MMVD_MAX_REFINE_NUM =                              (MMVD_REFINE_STEP * MMVD_MAX_DIR); ///< max number of candidate from a base candidate
#else
static const int MMVD_MAX_REFINE_NUM =                              (MMVD_REFINE_STEP * 4); ///< max number of candidate from a base candidate
#endif

static const int MMVD_BASE_MV_NUM =                                 2; ///< max number of base candidate
static const int MMVD_ADD_NUM =                                     (MMVD_MAX_REFINE_NUM * MMVD_BASE_MV_NUM);///< total number of mmvd candidate
#if MERGE_ENC_OPT
static const int MMVD_MRG_MAX_RD_NUM =                              20;
#else
static const int MMVD_MRG_MAX_RD_NUM =                              MRG_MAX_NUM_CANDS;
#endif
static const int MMVD_MRG_MAX_RD_BUF_NUM =                          (MMVD_MRG_MAX_RD_NUM + 1);///< increase buffer size by 1
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int LAST_MERGE_MMVD_IDX_CABAC =                        5;
#endif
#if JVET_W0097_GPM_MMVD_TM
static const int GPM_MMVD_REFINE_STEP =                             8;
static const int GPM_MMVD_REFINE_DIRECTION =                        4;
static const int GPM_MMVD_MAX_REFINE_NUM =                          (GPM_MMVD_REFINE_STEP * GPM_MMVD_REFINE_DIRECTION);
static const int GPM_EXT_MMVD_REFINE_STEP =                         9;
static const int GPM_EXT_MMVD_REFINE_DIRECTION =                    8;
static const int GPM_EXT_MMVD_MAX_REFINE_NUM =                      (GPM_EXT_MMVD_REFINE_STEP * GPM_EXT_MMVD_REFINE_DIRECTION);
#endif
static const int MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA =      28;
static const int MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA =    28;

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
static const int BDOF_SUBPU_DIM_LOG2          =                     3;
static const int BDOF_SUBPU_DIM               =                     (1 << BDOF_SUBPU_DIM_LOG2);
static const int BDOF_SUBPU_MAX_NUM           =                     ((MAX_CU_SIZE * MAX_CU_SIZE) >> (BDOF_SUBPU_DIM_LOG2 << 1));
static const int BDOF_SUBPU_STRIDE            =                     (MAX_CU_SIZE >> BDOF_SUBPU_DIM_LOG2);
static const int BDOF_SUBPU_SIZE              =                     (1 << (BDOF_SUBPU_DIM_LOG2 << 1));
static const int BIO_EXTEND_SIZE              =                     3;
#else
static const int BIO_EXTEND_SIZE              =                     1;
#endif
static const int BIO_TEMP_BUFFER_SIZE         =                     (MAX_CU_SIZE + 2 * BIO_EXTEND_SIZE) * (MAX_CU_SIZE + 2 * BIO_EXTEND_SIZE);

static const int PROF_BORDER_EXT_W            =                     1;
static const int PROF_BORDER_EXT_H            =                     1;
static const int BCW_NUM =                                          5; ///< the number of weight options
static const int BCW_DEFAULT =                                      ((uint8_t)(BCW_NUM >> 1)); ///< Default weighting index representing for w=0.5
static const int BCW_SIZE_CONSTRAINT =                            256; ///< disabling Bcw if cu size is smaller than 256
#if NON_ADJACENT_MRG_CAND
static const int MAX_NUM_HMVP_CANDS =                              5; ///< maximum number of HMVP candidates to be stored and used in merge list
#else
static const int MAX_NUM_HMVP_CANDS =                              (MRG_MAX_NUM_CANDS-1); ///< maximum number of HMVP candidates to be stored and used in merge list
#endif
#if JVET_Z0139_HIST_AFF
static const int MAX_NUM_AFF_HMVP_CANDS =                          7; 
#endif
#if JVET_Z0139_HIST_AFF
static const int MAX_NUM_AFF_INHERIT_HMVP_CANDS =                  9; 
#endif
#if JVET_Z0075_IBC_HMVP_ENLARGE
static const int MAX_NUM_HMVP_IBC_CANDS =                          25; ///< maximum number of HMVP candidates to be stored and used in IBC merge list
#endif
static const int MAX_NUM_HMVP_AVMPCANDS =                          4; ///< maximum number of HMVP candidates to be used in AMVP list

static const int ALF_VB_POS_ABOVE_CTUROW_LUMA = 4;
static const int ALF_VB_POS_ABOVE_CTUROW_CHMA = 2;

#if W0038_DB_OPT
static const int MAX_ENCODER_DEBLOCKING_QUALITY_LAYERS =           8 ;
#endif

#if SHARP_LUMA_DELTA_QP
static const uint32_t LUMA_LEVEL_TO_DQP_LUT_MAXSIZE =                1024; ///< max LUT size for QP offset based on luma

#endif
static const int DMVR_SUBCU_WIDTH = 16;
static const int DMVR_SUBCU_HEIGHT = 16;
static const int DMVR_SUBCU_WIDTH_LOG2 = 4;
static const int DMVR_SUBCU_HEIGHT_LOG2 = 4;
static const int DMVR_SUBPU_STRIDE = (MAX_CU_SIZE >> DMVR_SUBCU_WIDTH_LOG2);
static const int MAX_NUM_SUBCU_DMVR = ((MAX_CU_SIZE * MAX_CU_SIZE) >> (DMVR_SUBCU_WIDTH_LOG2 + DMVR_SUBCU_HEIGHT_LOG2));
#if MULTI_PASS_DMVR
static const int MAX_NUM_SUBCU_DMVR_LOG2 = MAX_CU_DEPTH + MAX_CU_DEPTH - DMVR_SUBCU_WIDTH_LOG2 - DMVR_SUBCU_HEIGHT_LOG2;
#endif
static const int DMVR_NUM_ITERATION = 2;

//QTBT high level parameters
//for I slice luma CTB configuration para.
static const int    MAX_BT_DEPTH  =                                 4;      ///<  <=7
                                                                            //for P/B slice CTU config. para.
static const int    MAX_BT_DEPTH_INTER =                            4;      ///< <=7
                                                                            //for I slice chroma CTB configuration para. (in luma samples)
static const int    MAX_BT_DEPTH_C      =                           0;      ///< <=7
static const int    MIN_DUALTREE_CHROMA_WIDTH  =                    4;
static const int    MIN_DUALTREE_CHROMA_SIZE   =                   16;
static const SplitSeries SPLIT_BITS         =                       5;
static const SplitSeries SPLIT_DMULT        =                       5;
static const SplitSeries SPLIT_MASK         =                      31;      ///< = (1 << SPLIT_BITS) - 1

static const int    SKIP_DEPTH =                                    3;
static const int    PICTURE_DISTANCE_TH =                           1;
static const int    FAST_SKIP_DEPTH =                               2;

static const double PBINTRA_RATIO     =                             1.1;
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int    THRES_TRANS =                                  16;
static const int    THRES_AFFINE =                                  4;
#endif
#if !MERGE_ENC_OPT
static const int    NUM_MRG_SATD_CAND =                             4;
#endif
static const double MRG_FAST_RATIO    =                             1.25;
static const int    NUM_AFF_MRG_SATD_CAND =                         2;
#if AFFINE_MMVD
static const int    AF_MMVD_BASE_NUM =                              1;
static const int    AF_MMVD_STEP_NUM =                              5; // number of distance offset
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int    AF_MMVD_OFFSET_DIR =                            8;
#else
static const int    AF_MMVD_OFFSET_DIR =                            4; // 00: (+, 0); 01: (-, 0); 10: (0, +); 11 (0, -);
#endif
#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int    AFFINE_MMVD_SIZE_SHIFT =                        1;
#endif
static const int    AF_MMVD_MAX_REFINE_NUM = AF_MMVD_STEP_NUM * AF_MMVD_OFFSET_DIR; ///< max number of candidate from a base candidate
static const int    AF_MMVD_NUM = AF_MMVD_BASE_NUM * AF_MMVD_MAX_REFINE_NUM;        ///< total number of affine mmvd candidate
#if !MERGE_ENC_OPT
static const int    NUM_AF_MMVD_SATD_CAND = std::min((int)1, MRG_MAX_NUM_CANDS);
#endif
#endif
#if INTER_LIC
static const int    LIC_MIN_CU_PIXELS =                            32; ///< smallest CU size (in terms of number of luma samples) of LIC
static const double LIC_AMVP_SKIP_TH  =                           1.2; ///< Given a IMV mode, LIC is not tested if RD cost of non-LIC IMV AMVP mode is 1.2x worse than the current best RD cost
#endif
#if TM_AMVP || TM_MRG
static const int    TM_TPL_SIZE =                                   4; ///< template size for template matching
static const int    TM_SEARCH_RANGE =                               8; ///< search range (in pixel) of TM refinement
static const int    TM_MIN_CU_SIZE_FOR_ALT_WEIGHTED_COST =          8; ///< minimal CU size (both width and height) to use alternative weight to compute template matching cost
static const int    TM_MAX_NUM_OF_ITERATIONS =                    375; ///< maximum number of refinement loops for entire CU in template mode
static const int    TM_LOG2_BASE_WEIGHT =                           5; ///< base weight of weighted template matching in log2 scale
static const int    TM_DISTANCE_WEIGHTS[][4] = { { 0, 1, 2, 3 }, { 1, 2, 3, 3 } }; ///< far to near
static const int    TM_SPATIAL_WEIGHTS [][4] = { { 2, 2, 2, 2 }, { 0, 1, 1, 2 } }; ///< "left to right for above template" or "top to bottom for left template"
#if TM_MRG
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
static const int    TM_MRG_MAX_NUM_INIT_CANDS =                    10; ///< maximum number of TM merge candidates for ARMC (note: should be at most equal to MRG_MAX_NUM_CANDS)
#endif
static const int    TM_MRG_MAX_NUM_CANDS =                          4; ///< maximum number of TM merge candidates (note: should be at most equal to MRG_MAX_NUM_CANDS)
#if JVET_X0141_CIIP_TIMD_TM
static const int    CIIP_TM_MRG_MAX_NUM_CANDS =                     2; ///< maximum number of CIIP TM merge candidates (note: should be at most equal to CIIP_TM_MRG_MAX_NUM_CANDS)
#endif
#if MERGE_ENC_OPT
static const int    TM_MAX_NUM_SATD_CAND = std::min((int)2, TM_MRG_MAX_NUM_CANDS);
#else
static const int    TM_MAX_NUM_SATD_CAND = std::min((int)4, TM_MRG_MAX_NUM_CANDS);
#endif
#endif
#endif
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
static const int    BM_MRG_MAX_NUM_INIT_CANDS =                    10; ///< maximum number of BM merge candidates (note: should be at most equal to MRG_MAX_NUM_CANDS)
#endif
#if MULTI_PASS_DMVR
static const int    BDMVR_INTME_RANGE =                             8; ///< Bilateral matching search range (n represents for -n pel to n pel, inclusive)
static const int    BDMVR_INTME_STRIDE = (BDMVR_INTME_RANGE << 1) + 1; ///< Bilateral matching search stride
static const int    BDMVR_INTME_AREA   = BDMVR_INTME_STRIDE * BDMVR_INTME_STRIDE; ///< Bilateral matching search area
static const int    BDMVR_INTME_CENTER = BDMVR_INTME_STRIDE * BDMVR_INTME_RANGE + BDMVR_INTME_RANGE; ///< Bilateral matching search area center
static const int    BDMVR_SIMD_IF_FACTOR =                          8; ///< Specify the pixel alignment factor for SIMD IF. (Usually this factor is 8)
static const int    BDMVR_INTME_SQUARE_SEARCH_MAX_NUM_ITERATIONS =  26; ///< maximum number of iterations in BDMVR integer square search at CU level
static const int    BDMVR_INTME_FULL_SEARCH_MAX_NUM_ITERATIONS   =  5; ///< maximum number of iterations in BDMVR integer full search at CU level (up to 5)
#if JVET_X0049_BDMVR_SW_OPT
static const int    BDMVR_BUF_STRIDE = MAX_CU_SIZE + (BDMVR_INTME_RANGE << 1) + (BDMVR_SIMD_IF_FACTOR - 2);
static const int    BDMVR_CENTER_POSITION = BDMVR_INTME_RANGE * BDMVR_BUF_STRIDE + BDMVR_INTME_RANGE;
static const int    BM_MRG_MAX_NUM_CANDS = 6; ///< maximum number of BM merge candidates (note: should be at most equal to MRG_MAX_NUM_CANDS)
static const int    BM_MRG_SUB_PU_INT_MAX_SRCH_ROUND = 3;
#endif
#endif
#if TM_AMVP || TM_MRG || MULTI_PASS_DMVR
static const int    DECODER_SIDE_MV_WEIGHT =                        4; ///< lambda for decoder-side derived MVs
#endif

static const double AMAXBT_TH32 =                                  15.0;
static const double AMAXBT_TH64 =                                  30.0;
#if CTU_256
static const double AMAXBT_TH128 =                                 60.0;
#endif

#if JVET_W0090_ARMC_TM
static const int AML_MERGE_TEMPLATE_SIZE =                         1;
static const int ADAPTIVE_SUB_GROUP_SIZE =                         5;
#if JVET_Z0139_HIST_AFF || JVET_Z0139_NA_AFF
static const int ADAPTIVE_AFFINE_SUB_GROUP_SIZE = AFFINE_MRG_MAX_NUM_CANDS;
#else
static const int ADAPTIVE_AFFINE_SUB_GROUP_SIZE =                  3;
#endif
#if JVET_Y0058_IBC_LIST_MODIFY
static const int ADAPTIVE_IBC_SUB_GROUP_SIZE =                     6;
#endif
#endif

#if JVET_Z0061_TM_OBMC
static const int TM_OBMC_TEMPLATE_SIZE =                           1;
#if (!JVET_W0090_ARMC_TM) && (!JVET_Z0056_GPM_SPLIT_MODE_REORDERING)
static const int AML_MERGE_TEMPLATE_SIZE =                         1;
#endif
#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int ADAPTIVE_SUB_GROUP_SIZE_MMVD =   MMVD_MAX_REFINE_NUM;
#endif

#if JVET_Y0067_ENHANCED_MMVD_MVD_SIGN_PRED
static const int ADAPTIVE_SUB_GROUP_SIZE_MMVD_AFF = AF_MMVD_MAX_REFINE_NUM;
#endif

#if JVET_Y0152_TT_ENC_SPEEDUP
static constexpr int   FAST_METHOD_TT_ENC_SPEEDUP = 0x0001;  ///< Embedding flag, which, if false, de-activates all the following ABT_ENC_SPEEDUP_* modes
static constexpr int FAST_METHOD_HOR_XOR_VER = 0x0002;
static constexpr int FAST_METHOD_ENC_SPEEDUP_BT_BASED = 0x0004;
static constexpr int FAST_METHOD_TT_ENC_SPEEDUP_BSLICE = 0x0008;
static constexpr int FAST_METHOD_TT_ENC_SPEEDUP_ISLICE = 0x0010;
#endif

// need to know for static memory allocation
static const int MAX_DELTA_QP   =                                   7;      ///< maximum supported delta QP value
static const int MAX_TESTED_QPs =   ( 1 + 1 + ( MAX_DELTA_QP << 1 ) );      ///< dqp=0 +- max_delta_qp + lossless mode

static const int COM16_C806_TRANS_PREC =                            0;
#if IF_12TAP
#define NTAPS_LUMA(x)                           ( (x) == 0 ? 12 : 8 )  // 12-tap filter for index 0. 8-tap fitler for other indices.
#else
static const int NTAPS_LUMA               =                         8; ///< Number of taps for luma
#endif
#if JVET_Z0117_CHROMA_IF
static const int NTAPS_CHROMA             =                         6; ///< Number of taps for chroma
static const int NTAPS_CHROMA_RPR         =                         4; ///< Number of taps for chroma RPR
#else
static const int NTAPS_CHROMA             =                         4; ///< Number of taps for chroma
#endif
#if LUMA_ADAPTIVE_DEBLOCKING_FILTER_QP_OFFSET
static const int MAX_LADF_INTERVALS       =                         5; /// max number of luma adaptive deblocking filter qp offset intervals
#endif

static const int NTAPS_BILINEAR           =                         2; ///< Number of taps for bilinear filter

#if INTER_RM_SIZE_CONSTRAINTS
static const int ATMVP_SUB_BLOCK_SIZE =                             2; ///< sub-block size for ATMVP
#else
static const int ATMVP_SUB_BLOCK_SIZE =                             3; ///< sub-block size for ATMVP
#endif

#if NON_ADJACENT_MRG_CAND
static const int GEO_MAX_NUM_UNI_CANDS =                            15;
#else
static const int GEO_MAX_NUM_UNI_CANDS =                            6;
#endif
#if JVET_Y0065_GPM_INTRA
static const int GEO_MAX_NUM_INTRA_CANDS =                          3;
static const int GEO_NUM_INTRA_RDO_BUFFER =                         23;
static const int GEO_NUM_RDO_BUFFER =                               GEO_MAX_NUM_UNI_CANDS + GEO_NUM_INTRA_RDO_BUFFER;
static const int GEO_MAX_NUM_CANDS = (GEO_MAX_NUM_UNI_CANDS+GEO_MAX_NUM_INTRA_CANDS) * ((GEO_MAX_NUM_UNI_CANDS+GEO_MAX_NUM_INTRA_CANDS) - 1);
#else
static const int GEO_MAX_NUM_CANDS = GEO_MAX_NUM_UNI_CANDS * (GEO_MAX_NUM_UNI_CANDS - 1);
#endif
static const int GEO_MIN_CU_LOG2 =                                  3;
static const int GEO_MAX_CU_LOG2 =                                  6;
static const int GEO_MIN_CU_SIZE =               1 << GEO_MIN_CU_LOG2;
static const int GEO_MAX_CU_SIZE =               1 << GEO_MAX_CU_LOG2;
static const int GEO_NUM_CU_SIZE = ( GEO_MAX_CU_LOG2 - GEO_MIN_CU_LOG2 ) + 1;
static const int GEO_NUM_PARTITION_MODE =                          64;
static const int GEO_NUM_ANGLES =                                  32;
static const int GEO_NUM_DISTANCES =                                4;
static const int GEO_NUM_PRESTORED_MASK =                           6;
static const int GEO_WEIGHT_MASK_SIZE = 3 * (GEO_MAX_CU_SIZE >> 3) * 2 + GEO_MAX_CU_SIZE;
#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
#if !JVET_W0090_ARMC_TM
static const int AML_MERGE_TEMPLATE_SIZE =                          1;
#endif
static const int GEO_MODE_SEL_TM_SIZE =       AML_MERGE_TEMPLATE_SIZE;
static const int GEO_TM_ADDED_WEIGHT_MASK_SIZE = GEO_MODE_SEL_TM_SIZE;
static const int GEO_WEIGHT_MASK_SIZE_EXT = GEO_WEIGHT_MASK_SIZE + GEO_TM_ADDED_WEIGHT_MASK_SIZE * 2;
static const int GEO_SPLIT_MODE_RICE_CODE_DIVISOR =                 4;
static const int GEO_MODE_COMPRESSION_RATIO =                       2;
static const int GEO_NUM_SIG_PARTMODE = GEO_NUM_PARTITION_MODE / GEO_MODE_COMPRESSION_RATIO; ///< max number of splitting modes for signaling
static const int GEO_ENC_MMVD_MAX_REFINE_NUM_ADJ = 1 // regular merge(1) 
#if JVET_W0097_GPM_MMVD_TM
                                                 + GPM_EXT_MMVD_MAX_REFINE_NUM + 1 // mmvd(GPM_EXT_MMVD_MAX_REFINE_NUM) + TM(1)
#endif
#if JVET_Y0065_GPM_INTRA
                                                 + 1 // intra(1)
#endif
;
#else
static const int GEO_MV_MASK_SIZE =         GEO_WEIGHT_MASK_SIZE >> 2;
#endif
#if JVET_W0097_GPM_MMVD_TM
static const int GEO_MAX_TRY_WEIGHTED_SAD = 70;
#if TM_MRG
static const int GEO_TM_MAX_NUM_CANDS = GEO_MAX_NUM_UNI_CANDS * (GEO_NUM_TM_MV_CAND - 1);
#endif
#else
static const int GEO_MAX_TRY_WEIGHTED_SAD = 60;
#endif
static const int GEO_MAX_TRY_WEIGHTED_SATD = 8;

#if ENABLE_OBMC
static const unsigned int defaultWeight[2][4] = { {27, 16, 6, 0}, {27, 0, 0, 0} };
#endif

#if JVET_X0141_CIIP_TIMD_TM
static const int CIIP_MAX_SIZE =                                 1024; ///< maximum CU size for using CIIP TIMD
#endif

#if TU_256
static const int SBT_MAX_SIZE =                          MAX_TB_SIZEY; ///< maximum CU size for using SBT
#else
static const int SBT_MAX_SIZE =                                    64; ///< maximum CU size for using SBT
#endif
static const int SBT_NUM_SL =                                      10; ///< maximum number of historical PU decision saved for a CU
static const int SBT_NUM_RDO =                                      2; ///< maximum number of SBT mode tried for a PU
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
static const int NUM_INTER_CU_INFO_SAVE =                           8; ///< maximum number of inter cu information saved for fast algorithm
static const int LDT_MODE_TYPE_INHERIT =                            0; ///< No need to signal mode_constraint_flag, and the modeType of the region is inherited from its parent node
static const int LDT_MODE_TYPE_INFER =                              1; ///< No need to signal mode_constraint_flag, and the modeType of the region is inferred as MODE_TYPE_INTRA
static const int LDT_MODE_TYPE_SIGNAL =                             2; ///< Need to signal mode_constraint_flag, and the modeType of the region is determined by the flag
#endif
static const int IBC_MAX_CAND_SIZE = 16; // max block size for ibc search
static const int IBC_NUM_CANDIDATES = 64; ///< Maximum number of candidates to store/test
static const int CHROMA_REFINEMENT_CANDIDATES = 8; /// 8 candidates BV to choose from
static const int IBC_FAST_METHOD_NOINTRA_IBCCBF0 = 0x01;
static const int IBC_FAST_METHOD_BUFFERBV = 0X02;
static const int IBC_FAST_METHOD_ADAPTIVE_SEARCHRANGE = 0X04;
static constexpr int MV_EXPONENT_BITCOUNT    = 4;
static constexpr int MV_MANTISSA_BITCOUNT    = 6;
static constexpr int MV_MANTISSA_UPPER_LIMIT = ((1 << (MV_MANTISSA_BITCOUNT - 1)) - 1);
static constexpr int MV_MANTISSA_LIMIT       = (1 << (MV_MANTISSA_BITCOUNT - 1));
static constexpr int MV_EXPONENT_MASK        = ((1 << MV_EXPONENT_BITCOUNT) - 1);

static constexpr int MV_BITS =                                   18;
static constexpr int MV_MAX =              (1 << (MV_BITS - 1)) - 1;
static constexpr int MV_MIN =                 -(1 << (MV_BITS - 1));

static const int MVD_MAX =                            (1 << 17) - 1;
static const int MVD_MIN =                               -(1 << 17);

static const int PIC_ANALYZE_CW_BINS =                           32;
static const int PIC_CODE_CW_BINS =                              16;
static const int LMCS_SEG_NUM =                                  32;
static const int FP_PREC =                                       11;
static const int CSCALE_FP_PREC =                                11;
static const int  NEIG_NUM_LOG  =                                 6;
static const int  NEIG_NUM =                      1 << NEIG_NUM_LOG;
static const int LOG2_PALETTE_CG_SIZE =                           4;
static const int RUN_IDX_THRE =                                   4;
static const int MAX_CU_BLKSIZE_PLT =                            64;
static const int NUM_TRELLIS_STATE =                              3;
static const double ENC_CHROMA_WEIGHTING =                      0.8;
static const int MAXPLTPREDSIZE = 63;
static const int MAXPLTSIZE = 31;
static const int MAXPLTPREDSIZE_DUALTREE = 31;
static const int MAXPLTSIZE_DUALTREE = 15;
static const double PLT_CHROMA_WEIGHTING =                      0.8;
static const int PLT_ENCBITDEPTH = 8;
static const int PLT_FAST_RATIO = 100;
#if RExt__DECODER_DEBUG_TOOL_MAX_FRAME_STATS
static const int  EPBIN_WEIGHT_FACTOR =                           4;
#endif
static const int ENC_PPS_ID_RPR =                                 3;
static const int SCALE_RATIO_BITS =                              14;
static const int MAX_SCALING_RATIO =                              2;  // max downsampling ratio for RPR
static const std::pair<int, int> SCALE_1X = std::pair<int, int>( 1 << SCALE_RATIO_BITS, 1 << SCALE_RATIO_BITS );  // scale ratio 1x
static const int DELTA_QP_ACT[4] =                  { -5, 1, 3, 1 };

#if SIGN_PREDICTION
static const int SIGN_PRED_MAX_NUM      = 8;
static const int SIGN_PRED_MAX_BS       = 128; ///< not configurable
static const int SIGN_PRED_MAX_BS_INTRA = 32;
static const int SIGN_PRED_MAX_BS_INTER = 128;
#if JVET_Y0141_SIGN_PRED_IMPROVE
static const int SIGN_PRED_FREQ_RANGE   = 32;///< not configurable
#else
static const int SIGN_PRED_FREQ_RANGE   = 4; ///< not configurable
#endif
static const int SIGN_PRED_SHIFT        = 8; ///< not configurable
static const int SIGN_PRED_OFFSET       = 1 << ( SIGN_PRED_SHIFT - 1 ); ///< not configurable
#endif
#if MULTI_HYP_PRED
static const auto MULTI_HYP_PRED_MAX_CANDS =                     4;
static const auto MULTI_HYP_PRED_NUM_WEIGHTS =                   2;
static const auto MULTI_HYP_PRED_WEIGHT_BITS =                   3;
static const auto MULTI_HYP_PRED_SEARCH_RANGE =                 16;
#endif
#if NON_ADJACENT_MRG_CAND || TM_AMVP
static const auto NADISTANCE_LEVEL =                             4;
#endif
#if JVET_Y0134_TMVP_NAMVP_CAND_REORDERING && JVET_W0090_ARMC_TM
static const auto TMVP_DISTANCE_LEVEL =                          5;
#endif
#if MULTI_HYP_PRED
static const auto MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE =          64; // disable multi-hyp for all blocks <= this number
static const auto MULTI_HYP_PRED_RESTRICT_MIN_WH =               8;
#endif

#if JVET_Z0139_HIST_AFF
static const auto AFF_PARA_STORE_BITS =                         16;
static const auto AFF_PARA_SHIFT =                               0;
#endif

#if JVET_Z0135_TEMP_CABAC_WIN_WEIGHT
static const int TEMP_CABAC_BUFFER_SIZE =                        5;
static const int ADJUSTMENT_RANGE =                              7;
#endif

// ====================================================================================================================
// Macro functions
// ====================================================================================================================

struct ClpRng
{
  int min;
  int max;
  int bd;
  int n;
};

struct ClpRngs
{
  ClpRng comp[MAX_NUM_COMPONENT]; ///< the bit depth as indicated in the SPS
  bool used;
  bool chroma;
};

template <typename T> inline T Clip3 (const T minVal, const T maxVal, const T a) { return std::min<T> (std::max<T> (minVal, a) , maxVal); }  ///< general min/max clip
template <typename T> inline T ClipBD( const T x, const int bitDepth ) { return Clip3( T( 0 ), T( ( 1 << bitDepth ) - 1 ), x ); }
template <typename T> inline T ClipPel (const T a, const ClpRng& clpRng)         { return std::min<T> (std::max<T> (clpRng.min, a) , clpRng.max); }  ///< clip reconstruction

template <typename T> inline void Check3( T minVal, T maxVal, T a)
{
  CHECK( ( a > maxVal ) || ( a < minVal ), "ERROR: Range check " << minVal << " >= " << a << " <= " << maxVal << " failed" );
}  ///< general min/max clip

extern MsgLevel g_verbosity;

#include <stdarg.h>
inline void msg( MsgLevel level, const char* fmt, ... )
{
  if( g_verbosity >= level )
  {
    va_list args;
    va_start( args, fmt );
    vfprintf( level == ERROR ? stderr : stdout, fmt, args );
    va_end( args );
  }
}

template<typename T> bool isPowerOf2( const T val ) { return ( val & ( val - 1 ) ) == 0; }

#define MEMORY_ALIGN_DEF_SIZE       32  // for use with avx2 (256 bit)
#define CACHE_MEM_ALIGN_SIZE      1024

#define ALIGNED_MALLOC              1   ///< use 32-bit aligned malloc/free

#if ALIGNED_MALLOC
#if JVET_J0090_MEMORY_BANDWITH_MEASURE
void *cache_mem_align_malloc(int size, int align_size);
void cache_mem_align_free(void *ptr);
#define xMalloc(type, len)          cache_mem_align_malloc(sizeof(type) * len, CACHE_MEM_ALIGN_SIZE)
#define xFree(ptr)                  cache_mem_align_free(ptr)
#elif     ( _WIN32 && ( _MSC_VER > 1300 ) ) || defined (__MINGW64_VERSION_MAJOR)
#define xMalloc( type, len )        _aligned_malloc( sizeof(type)*(len), MEMORY_ALIGN_DEF_SIZE )
#define xFree( ptr )                _aligned_free  ( ptr )
#elif defined (__MINGW32__)
#define xMalloc( type, len )        __mingw_aligned_malloc( sizeof(type)*(len), MEMORY_ALIGN_DEF_SIZE )
#define xFree( ptr )                __mingw_aligned_free( ptr )
#else
namespace detail {
template<typename T>
T* aligned_malloc(size_t len, size_t alignement) {
  T* p = NULL;
  if( posix_memalign( (void**)&p, alignement, sizeof(T)*(len) ) )
  {
    THROW("posix_memalign failed");
  }
  return p;
}
}
#define xMalloc( type, len )        detail::aligned_malloc<type>( len, MEMORY_ALIGN_DEF_SIZE )
#define xFree( ptr )                free( ptr )
#endif

#else
#define xMalloc( type, len )        malloc   ( sizeof(type)*(len) )
#define xFree( ptr )                free     ( ptr )
#endif //#if ALIGNED_MALLOC

#if defined _MSC_VER
#define ALIGN_DATA(nBytes,v) __declspec(align(nBytes)) v
#else
//#elif defined linux
#define ALIGN_DATA(nBytes,v) v __attribute__ ((aligned (nBytes)))
//#else
//#error unknown platform
#endif

#if defined(__GNUC__) && !defined(__clang__)
#    define GCC_VERSION_AT_LEAST(x,y) (__GNUC__ > x || __GNUC__ == x && __GNUC_MINOR__ >= y)
#else
#    define GCC_VERSION_AT_LEAST(x,y) 0
#endif

#ifdef __clang__
#    define CLANG_VERSION_AT_LEAST(x,y) (__clang_major__ > x || __clang_major__ == x && __clang_minor__ >= y)
#else
#    define CLANG_VERSION_AT_LEAST(x,y) 0
#endif

#ifdef __GNUC__
#    define ALWAYS_INLINE __attribute__((always_inline)) inline
#elif defined _MSC_VER
#    define ALWAYS_INLINE __forceinline
#else
#    define ALWAYS_INLINE
#endif

#if ENABLE_SIMD_OPT

#if defined(__i386__) || defined(i386) || defined(__x86_64__) || defined(_M_X64) || defined (_WIN32) || defined (_MSC_VER)
#define TARGET_SIMD_X86
typedef enum{
  SCALAR = 0,
  SSE41,
  SSE42,
  AVX,
  AVX2,
  AVX512
} X86_VEXT;
#elif defined (__ARM_NEON__)
#define TARGET_SIMD_ARM 1
#else
#error no simd target
#endif

#ifdef TARGET_SIMD_X86
X86_VEXT read_x86_extension_flags(const std::string &extStrId = std::string());
const char* read_x86_extension(const std::string &extStrId);
#endif

#endif //ENABLE_SIMD_OPT

template <typename ValueType> inline ValueType leftShift       (const ValueType value, const int shift) { return (shift >= 0) ? ( value                                  << shift) : ( value                                   >> -shift); }
template <typename ValueType> inline ValueType rightShift      (const ValueType value, const int shift) { return (shift >= 0) ? ( value                                  >> shift) : ( value                                   << -shift); }
template <typename ValueType> inline ValueType leftShift_round (const ValueType value, const int shift) { return (shift >= 0) ? ( value                                  << shift) : ((value + (ValueType(1) << (-shift - 1))) >> -shift); }
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
template <typename ValueType> inline ValueType rightShift_round(const ValueType value, const int shift) { return (shift > 0) ? ((value + (ValueType(1) << (shift - 1))) >> shift) : ( value                                   << -shift); }
#else
template <typename ValueType> inline ValueType rightShift_round(const ValueType value, const int shift) { return (shift >= 0) ? ((value + (ValueType(1) << (shift - 1))) >> shift) : ( value                                   << -shift); }
#endif

static inline int floorLog2(uint32_t x)
{
  if (x == 0)
  {
    // note: ceilLog2() expects -1 as return value
    return -1;
  }
#ifdef __GNUC__
  return 31 - __builtin_clz(x);
#else
#ifdef _MSC_VER
  unsigned long r = 0;
  _BitScanReverse(&r, x);
  return r;
#else
  int result = 0;
  if (x & 0xffff0000)
  {
    x >>= 16;
    result += 16;
  }
  if (x & 0xff00)
  {
    x >>= 8;
    result += 8;
  }
  if (x & 0xf0)
  {
    x >>= 4;
    result += 4;
  }
  if (x & 0xc)
  {
    x >>= 2;
    result += 2;
  }
  if (x & 0x2)
  {
    x >>= 1;
    result += 1;
  }
  return result;
#endif
#endif
}

#if JVET_X0149_TIMD_DIMD_LUT
static inline int floorLog2_uint64(uint64_t x)
{
  if (x == 0)
  {
    // note: ceilLog2() expects -1 as return value
    return -1;
  }
  int result = 0;
  if (x & 0xffffffff00000000)
  {
    x >>= 32;
    result += 32;
  }
  if (x & 0xffff0000)
  {
    x >>= 16;
    result += 16;
  }
  if (x & 0xff00)
  {
    x >>= 8;
    result += 8;
  }
  if (x & 0xf0)
  {
    x >>= 4;
    result += 4;
  }
  if (x & 0xc)
  {
    x >>= 2;
    result += 2;
  }
  if (x & 0x2)
  {
    x >>= 1;
    result += 1;
  }
  return result;
}
#endif

static inline int ceilLog2(uint32_t x)
{
  return (x==0) ? -1 : floorLog2(x - 1) + 1;
}

#if INTER_LIC
static inline int getMSB(unsigned x)
{
  int msb = 0, bits = (sizeof(int) << 3), y = 1;
  while (x > 1u)
  {
    bits >>= 1;
    y = x >> bits;
    if (y)
    {
      x = y;
      msb += bits;
    }
}
  msb += y;
  return msb;
}
#endif

#if JVET_Z0056_GPM_SPLIT_MODE_REORDERING
template <typename InputValueType, uint32_t inputValueArraySize, typename OutputIndexType, uint32_t outputIndexArraySize>
static inline uint32_t getIndexMappingTableToSortedArray1D(InputValueType (&in)[inputValueArraySize], OutputIndexType (&tbl)[outputIndexArraySize])
{
  uint32_t numValidInList = 1;
  InputValueType sortedlist[outputIndexArraySize];
  sortedlist[0] = in[0];
  tbl[0] = (OutputIndexType)0;

  for (int inIdx = 1; inIdx < inputValueArraySize; ++inIdx)
  {
    // Find insertion index position using binary search
    int insertIdx = 0;
    InputValueType* subList = sortedlist;
    uint32_t subListSize = numValidInList;
    while (subListSize > 1)
    {
      int middleIdx = subListSize >> 1;
      if (in[inIdx] < subList[middleIdx])
      {
        subListSize = middleIdx;
      }
      else
      {
        subList     += middleIdx;
        subListSize -= middleIdx;
        insertIdx   += middleIdx;
      }
    }
    insertIdx += (in[inIdx] < subList[0] ? 0 : 1);

    // Perform insertion at found index position
    if (insertIdx < outputIndexArraySize)
    {
      int startIdx = outputIndexArraySize - 1;
      if (numValidInList < outputIndexArraySize)
      {
        startIdx = numValidInList;
        ++numValidInList;
      }

      for (int i = startIdx; i > insertIdx; --i)
      {
        sortedlist[i] = sortedlist[i - 1];
        tbl       [i] = tbl       [i - 1];
      }
      sortedlist[insertIdx] = in[inIdx];
      tbl       [insertIdx] = (OutputIndexType)inIdx;
    }
  }

  return numValidInList;
}
#endif

//CASE-BREAK for breakpoints
#if defined ( _MSC_VER ) && defined ( _DEBUG )
#define _CASE(_x) if(_x)
#define _BREAK while(0);
#define _AREA_AT(_a,_x,_y,_w,_h)  (_a.x==_x && _a.y==_y && _a.width==_w && _a.height==_h)
#define _AREA_CONTAINS(_a,_x,_y)  (_a.contains( Position{ _x, _y} ))
#define _UNIT_AREA_AT(_a,_x,_y,_w,_h) (_a.Y().x==_x && _a.Y().y==_y && _a.Y().width==_w && _a.Y().height==_h)
#else
#define _CASE(...)
#define _BREAK
#define _AREA_AT(...)
#define _AREA_CONTAINS(_a,_x,_y)
#define _UNIT_AREA_AT(_a,_x,_y,_w,_h)
#endif

#if ENABLE_SPLIT_PARALLELISM
#include <omp.h>

#define PARL_PARAM(DEF) , DEF
#define PARL_PARAM0(DEF) DEF
#else
#define PARL_PARAM(DEF)
#define PARL_PARAM0(DEF)
#endif

static const uint32_t CCALF_CANDS_COEFF_NR = 8;
static const int CCALF_SMALL_TAB[CCALF_CANDS_COEFF_NR] = { 0, 1, 2, 4, 8, 16, 32, 64 };

//! \}
#if ALF_IMPROVEMENT
static const int NUM_FIXED_FILTERS       = 512;
static const int ALF_CLASSIFIER_FL       = 5;
static const int NUM_CLASSIFIER          = 2;
static const int NUM_SETS_FIXED_FILTERS  = 8;
static const int NUM_DIR_FIX             = 7;
static const int NUM_ACT_FIX             = 16;
static const int NUM_CLASSES_FIX         = ((NUM_DIR_FIX*(NUM_DIR_FIX + 1))*NUM_ACT_FIX);
static const int MAX_FILTER_LENGTH_FIXED = 13;
static const int FIX_FILTER_NUM_COEFF    = 42;
#endif

#if JVET_V0130_INTRA_TMP
static const int TMP_TEMPLATE_SIZE =            4; // must be multiple of 4 for SIMD
static const int TMP_MAXSIZE_DEPTH =            6; // should be log2(TMP_TEMPLATE_SIZE): keep as 6 to avoid any error
static const int USE_MORE_BLOCKSIZE_DEPTH_MAX = TMP_MAXSIZE_DEPTH - 1;
static const int INIT_THRESHOULD_SHIFTBITS =    2;  ///< (default 2) Early skip threshold for checking distance.
static const int TMP_SEARCH_RANGE_MULT_FACTOR = 5;
#endif

#endif // end of #ifndef  __COMMONDEF__
