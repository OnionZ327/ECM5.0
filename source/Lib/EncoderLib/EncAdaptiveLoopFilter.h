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

/** \file     EncAdaptiveLoopFilter.h
 \brief    estimation part of adaptive loop filter class (header)
 */

#ifndef __ENCADAPTIVELOOPFILTER__
#define __ENCADAPTIVELOOPFILTER__

#include "CommonLib/AdaptiveLoopFilter.h"
#include "CommonLib/ParameterSetManager.h"

#include "CABACWriter.h"
#include "EncCfg.h"

struct AlfCovariance
{
  static constexpr int MaxAlfNumClippingValues = AdaptiveLoopFilter::MaxAlfNumClippingValues;

#if JVET_X0071_LONGER_CCALF
  using TE = double[MAX_NUM_CC_ALF_CHROMA_COEFF][MAX_NUM_CC_ALF_CHROMA_COEFF];
  using Ty = double[MAX_NUM_CC_ALF_CHROMA_COEFF];
#else
  using TE = double[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
  using Ty = double[MAX_NUM_ALF_LUMA_COEFF];
#endif


  using TKE = TE[AdaptiveLoopFilter::MaxAlfNumClippingValues][AdaptiveLoopFilter::MaxAlfNumClippingValues];
  using TKy = Ty[AdaptiveLoopFilter::MaxAlfNumClippingValues];

  int numCoeff;
  int numBins;
  TKy y;
  TKE E;
  double pixAcc;

  AlfCovariance() {}
  ~AlfCovariance() {}

  void create( int size, int num_bins = MaxAlfNumClippingValues )
  {
    numCoeff = size;
    numBins = num_bins;
    std::memset( y, 0, sizeof( y ) );
    std::memset( E, 0, sizeof( E ) );
  }

  void destroy()
  {
  }

  void reset( int num_bins = -1 )
  {
    if ( num_bins > 0 )
      numBins = num_bins;
    pixAcc = 0;
    std::memset( y, 0, sizeof( y ) );
    std::memset( E, 0, sizeof( E ) );
  }

  const AlfCovariance& operator=( const AlfCovariance& src )
  {
    numCoeff = src.numCoeff;
    numBins = src.numBins;
    std::memcpy( E, src.E, sizeof( E ) );
    std::memcpy( y, src.y, sizeof( y ) );
    pixAcc = src.pixAcc;

    return *this;
  }

  void add( const AlfCovariance& lhs, const AlfCovariance& rhs )
  {
    numCoeff = lhs.numCoeff;
    numBins = lhs.numBins;
    for( int b0 = 0; b0 < numBins; b0++ )
    {
      for( int b1 = 0; b1 < numBins; b1++ )
      {
        for( int j = 0; j < numCoeff; j++ )
        {
          for( int i = 0; i < numCoeff; i++ )
          {
            E[b0][b1][j][i] = lhs.E[b0][b1][j][i] + rhs.E[b0][b1][j][i];
          }
        }
      }
    }
    for( int b = 0; b < numBins; b++ )
    {
      for( int j = 0; j < numCoeff; j++ )
      {
        y[b][j] = lhs.y[b][j] + rhs.y[b][j];
      }
    }
    pixAcc = lhs.pixAcc + rhs.pixAcc;
  }

  const AlfCovariance& operator+= ( const AlfCovariance& src )
  {
    for( int b0 = 0; b0 < numBins; b0++ )
    {
      for( int b1 = 0; b1 < numBins; b1++ )
      {
        for( int j = 0; j < numCoeff; j++ )
        {
          for( int i = 0; i < numCoeff; i++ )
          {
            E[b0][b1][j][i] += src.E[b0][b1][j][i];
          }
        }
      }
    }
    for( int b = 0; b < numBins; b++ )
    {
      for( int j = 0; j < numCoeff; j++ )
      {
        y[b][j] += src.y[b][j];
      }
    }
    pixAcc += src.pixAcc;

    return *this;
  }

  const AlfCovariance& operator-= ( const AlfCovariance& src )
  {
    for( int b0 = 0; b0 < numBins; b0++ )
    {
      for( int b1 = 0; b1 < numBins; b1++ )
      {
        for( int j = 0; j < numCoeff; j++ )
        {
          for( int i = 0; i < numCoeff; i++ )
          {
            E[b0][b1][j][i] -= src.E[b0][b1][j][i];
          }
        }
      }
    }
    for( int b = 0; b < numBins; b++ )
    {
      for( int j = 0; j < numCoeff; j++ )
      {
        y[b][j] -= src.y[b][j];
      }
    }
    pixAcc -= src.pixAcc;

    return *this;
  }

  void setEyFromClip(const int* clip, TE _E, Ty _y, int size) const
  {
    for (int k=0; k<size; k++)
    {
      _y[k] = y[clip[k]][k];
      for (int l=0; l<size; l++)
      {
        _E[k][l] = E[clip[k]][clip[l]][k][l];
      }
    }
  }

  double optimizeFilter(const int* clip, double *f, int size) const
  {
    gnsSolveByChol( clip, f, size );
    return calculateError( clip, f );
  }

  double optimizeFilter(const AlfFilterShape& alfShape, int* clip, double *f, bool optimize_clip) const;
  double optimizeFilterClip(const AlfFilterShape& alfShape, int* clip) const
  {
    Ty f;
    return optimizeFilter(alfShape, clip, f, true);
  }

  double calculateError( const int *clip ) const;
  double calculateError( const int *clip, const double *coeff ) const { return calculateError(clip, coeff, numCoeff); }
  double calculateError( const int *clip, const double *coeff, const int numCoeff ) const;
  double calcErrorForCoeffs( const int *clip, const int *coeff, const int numCoeff, const int bitDepth ) const;
  double calcErrorForCcAlfCoeffs(const int16_t *coeff, const int numCoeff, const int bitDepth) const;

  void getClipMax(const AlfFilterShape& alfShape, int *clip_max) const;
  void reduceClipCost(const AlfFilterShape& alfShape, int *clip) const;

  int  gnsSolveByChol( TE LHS, double* rhs, double *x, int numEq ) const;

private:
  // Cholesky decomposition

  int  gnsSolveByChol( const int *clip, double *x, int numEq ) const;
  void gnsBacksubstitution( TE R, double* z, int size, double* A ) const;
  void gnsTransposeBacksubstitution( TE U, double* rhs, double* x, int order ) const;
  int  gnsCholeskyDec( TE inpMatr, TE outMatr, int numEq ) const;
};

class EncAdaptiveLoopFilter : public AdaptiveLoopFilter
{
public:
  inline void           setAlfWSSD(int alfWSSD) { m_alfWSSD = alfWSSD; }
  static std::vector<double>  m_lumaLevelToWeightPLUT;
  inline std::vector<double>& getLumaLevelWeightTable() { return m_lumaLevelToWeightPLUT; }

private:
  int                    m_alfWSSD;
  const EncCfg*          m_encCfg;
#if ALF_IMPROVEMENT
#if JVET_X0071_ALF_BAND_CLASSIFIER
  AlfCovariance*****       m_alfCovariance[MAX_NUM_COMPONENT];          // [compIdx][shapeIdx][ctbAddr][fixedFilterSetIdx][classifierInd][classIdx]
#else
  AlfCovariance****       m_alfCovariance[MAX_NUM_COMPONENT];          // [compIdx][shapeIdx][ctbAddr][fixedFilterSetIdx][classIdx]
#endif
#else
  AlfCovariance***       m_alfCovariance[MAX_NUM_COMPONENT];          // [compIdx][shapeIdx][ctbAddr][classIdx]
#endif
  AlfCovariance**        m_alfCovarianceFrame[MAX_NUM_CHANNEL_TYPE];   // [CHANNEL][shapeIdx][lumaClassIdx/chromaAltIdx]
  uint8_t*               m_ctuEnableFlagTmp[MAX_NUM_COMPONENT];
  uint8_t*               m_ctuEnableFlagTmp2[MAX_NUM_COMPONENT];
  uint8_t*               m_ctuAlternativeTmp[MAX_NUM_COMPONENT];
  AlfCovariance***       m_alfCovarianceCcAlf[2];           // [compIdx-1][shapeIdx][filterIdx][ctbAddr]
  AlfCovariance**        m_alfCovarianceFrameCcAlf[2];      // [compIdx-1][shapeIdx][filterIdx]

#if ALF_IMPROVEMENT
  bool classChanged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES];
  int clipHistory[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
  double errorHistory[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES];
#endif

  //for RDO
  AlfParam               m_alfParamTemp;
  ParameterSetMap<APS>*  m_apsMap;
  AlfCovariance          m_alfCovarianceMerged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES + 2];
  int                    m_alfClipMerged[ALF_NUM_OF_FILTER_TYPES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
  CABACWriter*           m_CABACEstimator;
  CtxCache*              m_CtxCache;
  double                 m_lambda[MAX_NUM_COMPONENT];

  int**                  m_filterCoeffSet; // [lumaClassIdx/chromaAltIdx][coeffIdx]
  int**                  m_filterClippSet; // [lumaClassIdx/chromaAltIdx][coeffIdx]
  int**                  m_diffFilterCoeff;
  short                  m_filterIndices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES];
  unsigned               m_bitsNewFilter[MAX_NUM_CHANNEL_TYPE];
  int&                   m_apsIdStart;
  double                 *m_ctbDistortionUnfilter[MAX_NUM_COMPONENT];
  double                 m_unFiltDistCompnent[MAX_NUM_COMPONENT];
#if ALF_IMPROVEMENT
  double                 *m_distCtbApsLuma[ALF_CTB_MAX_NUM_APS][MAX_NUM_ALF_ALTERNATIVES_LUMA][2];
  double                 *m_distCtbLumaNewFilt[MAX_NUM_ALF_ALTERNATIVES_LUMA][2];
  double                 *m_ctbDistortionFixedFilter[NUM_FIXED_FILTER_SETS][2];
#else
  double                 *m_distCtbApsLuma[ALF_CTB_MAX_NUM_APS];
  double                 *m_distCtbLumaNewFilt;
  double                 *m_ctbDistortionFixedFilter[NUM_FIXED_FILTER_SETS];
  int                    m_clipDefaultEnc[MAX_NUM_ALF_LUMA_COEFF];
#endif

  std::vector<short>     m_alfCtbFilterSetIndexTmp;
  AlfParam               m_alfParamTempNL;
  int                    m_filterTmp[MAX_NUM_ALF_LUMA_COEFF];
  int                    m_clipTmp[MAX_NUM_ALF_LUMA_COEFF];

  int m_apsIdCcAlfStart[2];

  short                  m_bestFilterCoeffSet[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
  bool                   m_bestFilterIdxEnabled[MAX_NUM_CC_ALF_FILTERS];
  uint8_t                m_bestFilterCount;
  uint8_t*               m_trainingCovControl;
  Pel*                   m_bufOrigin;
  PelBuf*                m_buf;
  uint64_t*              m_trainingDistortion[MAX_NUM_CC_ALF_FILTERS];    // for current block size
  uint64_t*              m_lumaSwingGreaterThanThresholdCount;
  uint64_t*              m_chromaSampleCountNearMidPoint;
  uint8_t*               m_filterControl;         // current iterations filter control
  uint8_t*               m_bestFilterControl;     // best saved filter control
  int                    m_reuseApsId[2];
  bool                   m_limitCcAlf;

public:
  EncAdaptiveLoopFilter( int& apsIdStart );
  virtual ~EncAdaptiveLoopFilter() {}

  template<bool alfWSSD>
  void  initDistortion(
#if ALF_IMPROVEMENT
    CodingStructure& cs
#endif
  );
  std::vector<int> getAvaiApsIdsLuma(CodingStructure& cs, int &newApsId);
  void  alfEncoderCtb(CodingStructure& cs, AlfParam& alfParamNewFilters
#if ENABLE_QPA
    , const double lambdaChromaWeight
#endif
  );
  void   alfReconstructor(CodingStructure& cs, const PelUnitBuf& recExtBuf);
  void ALFProcess(CodingStructure& cs, const double *lambdas
#if ENABLE_QPA
    , const double lambdaChromaWeight
#endif
    , Picture* pcPic, uint32_t numSliceSegments
  );

  void getDistApsFilter( CodingStructure& cs, std::vector<int> apsIds );
  void getDistNewFilter( AlfParam& alfParam );
  int getNewCcAlfApsId( CodingStructure &cs, int cIdx );
  void initCABACEstimator( CABACEncoder* cabacEncoder, CtxCache* ctxCache, Slice* pcSlice, ParameterSetMap<APS>* apsMap );
  void create( const EncCfg* encCfg, const int picWidth, const int picHeight, const ChromaFormat chromaFormatIDC, const int maxCUWidth, const int maxCUHeight, const int maxCUDepth, const int inputBitDepth[MAX_NUM_CHANNEL_TYPE], const int internalBitDepth[MAX_NUM_CHANNEL_TYPE] );
  void destroy();
  void setApsIdStart( int i) { m_apsIdStart = i; }
#if ALF_IMPROVEMENT
  static int lengthGolomb(int coeffVal, int k, bool signed_coeff = true);
#endif

private:
  void   alfEncoder( CodingStructure& cs, AlfParam& alfParam, const PelUnitBuf& orgUnitBuf, const PelUnitBuf& recExtBuf, const PelUnitBuf& recBuf, const ChannelType channel
#if ENABLE_QPA
                   , const double lambdaChromaWeight = 0.0
#endif
                   );

  void   copyAlfParam( AlfParam& alfParamDst, AlfParam& alfParamSrc, ChannelType channel );
#if JVET_X0071_ALF_BAND_CLASSIFIER
  double mergeFiltersAndCost( AlfParam& alfParam, AlfFilterShape& alfShape, AlfCovariance* covFrame, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], int& uiCoeffBits, int altIdx, int classifierIdx, int numFiltersLinear );
  void   getFrameStats( ChannelType channel, int iShapeIdx, int altIdx, int fixedFilterSetIdx, int classifierIdx );
#else
  double mergeFiltersAndCost( AlfParam& alfParam, AlfFilterShape& alfShape, AlfCovariance* covFrame, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], int& uiCoeffBits
#if ALF_IMPROVEMENT
    , int altIdx
#endif
  );
  void   getFrameStats( ChannelType channel, int iShapeIdx, int altIdx 
#if ALF_IMPROVEMENT
    , int fixedFilterSetIdx
#endif
  );
#endif
#if ALF_IMPROVEMENT
#if JVET_X0071_ALF_BAND_CLASSIFIER
  void   getFrameStat( AlfCovariance* frameCov, AlfCovariance**** ctbCov, uint8_t* ctbEnableFlags, uint8_t* ctbAltIdx, const int numClasses, int altIdx, int fixedFilterSetIdx, int classifierIdx );
#else
  void   getFrameStat( AlfCovariance* frameCov, AlfCovariance*** ctbCov, uint8_t* ctbEnableFlags, uint8_t* ctbAltIdx, const int numClasses, int altIdx, int fixedFilterSetIdx);
#endif
#else
  void   getFrameStat( AlfCovariance* frameCov, AlfCovariance** ctbCov, uint8_t* ctbEnableFlags, uint8_t* ctbAltIdx, const int numClasses, int altIdx );
#endif
  template<bool alfWSSD>
  void   deriveStatsForFiltering( PelUnitBuf& orgYuv, PelUnitBuf& recYuv, CodingStructure& cs );

#if ALF_IMPROVEMENT
  template<bool m_alfWSSD>
#if JVET_X0071_ALF_BAND_CLASSIFIER
  void   getBlkStats(AlfCovariance* alfCovariace, const AlfFilterShape& shape, AlfClassifier** classifier, const Pel* org, const int orgStride, const Pel* rec, const int recStride, const CompArea& areaDst, const CompArea& area, const ChannelType channel, int fixedFilterSetIdx, int classifierIdx);
#else
  void   getBlkStats(AlfCovariance* alfCovariace, const AlfFilterShape& shape, AlfClassifier** classifier, const Pel* org, const int orgStride, const Pel* rec, const int recStride, const CompArea& areaDst, const CompArea& area, const ChannelType channel, int fixedFilterSetIdx);
#endif
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  void   calcCovariance(Pel ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, Pel ***fixedFitlerResults, Position pos, int fixedFilterSetIdx);
#else  
  void   calcCovariance(int ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, Pel ***fixedFitlerResults, Position pos, int fixedFilterSetIdx);
#endif   
#else
  void   getBlkStats(AlfCovariance* alfCovariace, const AlfFilterShape& shape, AlfClassifier** classifier, Pel* org, const int orgStride, Pel* rec, const int recStride, const CompArea& areaDst, const CompArea& area, const ChannelType channel, int vbCTUHeight, int vbPos);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  void   calcCovariance(Pel ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, int vbDistance);
#else
  void   calcCovariance(int ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, int vbDistance);
#endif
#endif

  template<bool alfWSSD>
  void   deriveStatsForCcAlfFiltering(const PelUnitBuf &orgYuv, const PelUnitBuf &recYuv, const int compIdx,
                                      const int maskStride, const uint8_t filterIdc, CodingStructure &cs);
  template<bool m_alfWSSD>
  void   getBlkStatsCcAlf(AlfCovariance &alfCovariance, const AlfFilterShape &shape, const PelUnitBuf &orgYuv,
                          const PelUnitBuf &recYuv, const UnitArea &areaDst, const UnitArea &area,
                          const ComponentID compID, const int yPos);
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
#if ALF_IMPROVEMENT
  void   calcCovarianceCcAlf( Pel ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel* rec, const int stride, const AlfFilterShape& shape );
#else
  void   calcCovarianceCcAlf(Pel ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel* rec, const int stride, const AlfFilterShape& shape, int vbDistance);
#endif
#else
#if ALF_IMPROVEMENT
  void   calcCovarianceCcAlf( int ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel* rec, const int stride, const AlfFilterShape& shape );
#else
  void   calcCovarianceCcAlf(int ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel* rec, const int stride, const AlfFilterShape& shape, int vbDistance);
#endif
#endif
#if JVET_X0071_ALF_BAND_CLASSIFIER
   void   mergeClasses(const AlfFilterShape& alfShape, AlfCovariance* cov, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], const int numClasses, short filterIndices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES], const int altIdx, int mergedPair[MAX_NUM_ALF_CLASSES][2]);
#else
  void   mergeClasses(const AlfFilterShape& alfShape, AlfCovariance* cov, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], const int numClasses, short filterIndices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES]
#if ALF_IMPROVEMENT
    , const int altIdx
#endif
  );
#endif

  double getFilterCoeffAndCost( CodingStructure& cs, double distUnfilter, ChannelType channel, bool bReCollectStat, int iShapeIdx, int& uiCoeffBits, 
#if ALF_IMPROVEMENT
    int fixedFilterSetIdx,
#endif
    bool onlyFilterCost = false );
#if JVET_X0071_ALF_BAND_CLASSIFIER
  double deriveFilterCoeffs(AlfCovariance* cov, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], AlfFilterShape& alfShape, short* filterIndices, int numFilters, double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2], AlfParam& alfParam, bool nonLinear, int classifierIdx, bool isMaxNum, int mergedPair[MAX_NUM_ALF_CLASSES][2], int mergedCoeff[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], double mergedErr[MAX_NUM_ALF_CLASSES]);
#else
  double deriveFilterCoeffs(AlfCovariance* cov, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], AlfFilterShape& alfShape, short* filterIndices, int numFilters, double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2], AlfParam& alfParam
#if ALF_IMPROVEMENT
  , bool nonLinear
#endif
  );
#endif
  int    deriveFilterCoefficientsPredictionMode( AlfFilterShape& alfShape, int **filterSet, int** filterCoeffDiff, const int numFilters );
  double deriveCoeffQuant( int *filterClipp, int *filterCoeffQuant, const AlfCovariance& cov, const AlfFilterShape& shape, const int bitDepth, const bool optimizeClip );
  double deriveCtbAlfEnableFlags( CodingStructure& cs, const int iShapeIdx, ChannelType channel,
#if ENABLE_QPA
                                  const double chromaWeight,
#endif
                                  const int numClasses, const int numCoeff, double& distUnfilter 
#if ALF_IMPROVEMENT
                                , int fixedFilterSetIdx
#endif
                                );
  void   roundFiltCoeff( int *filterCoeffQuant, double *filterCoeff, const int numCoeff, const int factor );
  void   roundFiltCoeffCCALF(int16_t *filterCoeffQuant, double *filterCoeff, const int numCoeff, const int factor);

  double getDistCoeffForce0( bool* codedVarBins, double errorForce0CoeffTab[MAX_NUM_ALF_CLASSES][2], int* bitsVarBin, int zeroBitsVarBin, const int numFilters);
  int    lengthUvlc( int uiCode );
#if ALF_IMPROVEMENT
#if JVET_X0071_ALF_BAND_CLASSIFIER
  int    getNonFilterCoeffRate( AlfParam& alfParam, int altIdx, int classifierIdx );
#else
  int    getNonFilterCoeffRate( AlfParam& alfParam, int altIdx );
#endif
  int    getCostFilterCoeffForce0( AlfFilterShape& alfShape, int **pDiffQFilterCoeffIntPP, const int numFilters, bool* codedVarBins, int altIdx);
  double getDistForce0( AlfFilterShape& alfShape, const int numFilters, double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2], bool* codedVarBins, int altIdx);
  double getFilteredDistortion( AlfCovariance* cov, const int numClasses, const int numFiltersMinus1, const int numCoeff, int altIdx );
  void initCtuAlternativeLuma( uint8_t* ctuAlts[MAX_NUM_COMPONENT] );
#else
  int    getNonFilterCoeffRate( AlfParam& alfParam );
  int    getCostFilterCoeffForce0( AlfFilterShape& alfShape, int **pDiffQFilterCoeffIntPP, const int numFilters, bool* codedVarBins );
  double getDistForce0( AlfFilterShape& alfShape, const int numFilters, double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2], bool* codedVarBins );
  double getFilteredDistortion( AlfCovariance* cov, const int numClasses, const int numFiltersMinus1, const int numCoeff );
#endif
  int    getCostFilterCoeff( AlfFilterShape& alfShape, int **pDiffQFilterCoeffIntPP, const int numFilters );
  int    getCostFilterClipp( AlfFilterShape& alfShape, int **pDiffQFilterCoeffIntPP, const int numFilters );
  int    lengthFilterCoeffs( AlfFilterShape& alfShape, const int numFilters, int **FilterCoeff );  
  int    getChromaCoeffRate( AlfParam& alfParam, int altIdx );
  double getUnfilteredDistortion( AlfCovariance* cov, ChannelType channel );
  double getUnfilteredDistortion( AlfCovariance* cov, const int numClasses );

  void setEnableFlag( AlfParam& alfSlicePara, ChannelType channel, bool val );
  void setEnableFlag( AlfParam& alfSlicePara, ChannelType channel, uint8_t** ctuFlags );
  void setCtuEnableFlag( uint8_t** ctuFlags, ChannelType channel, uint8_t val );
  void copyCtuEnableFlag( uint8_t** ctuFlagsDst, uint8_t** ctuFlagsSrc, ChannelType channel );
  void initCtuAlternativeChroma( uint8_t* ctuAlts[MAX_NUM_COMPONENT] );
  void setCtuAlternativeChroma( uint8_t* ctuAlts[MAX_NUM_COMPONENT], uint8_t val );
  void copyCtuAlternativeChroma( uint8_t* ctuAltsDst[MAX_NUM_COMPONENT], uint8_t* ctuAltsSrc[MAX_NUM_COMPONENT] );
  int getMaxNumAlternativesChroma( );
  int  getCoeffRateCcAlf(short chromaCoeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF], bool filterEnabled[MAX_NUM_CC_ALF_FILTERS], uint8_t filterCount, ComponentID compID);
  void deriveCcAlfFilterCoeff( ComponentID compID, const PelUnitBuf& recYuv, const PelUnitBuf& recYuvExt, short filterCoeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF], const uint8_t filterIdx );
  void determineControlIdcValues(CodingStructure &cs, const ComponentID compID, const PelBuf *buf, const int ctuWidthC,
                                 const int ctuHeightC, const int picWidthC, const int picHeightC,
                                 double **unfilteredDistortion, uint64_t *trainingDistortion[MAX_NUM_CC_ALF_FILTERS],
                                 uint64_t *lumaSwingGreaterThanThresholdCount,
                                 uint64_t *chromaSampleCountNearMidPoint,
                                 bool reuseFilterCoeff, uint8_t *trainingCovControl, uint8_t *filterControl,
                                 uint64_t &curTotalDistortion, double &curTotalRate,
                                 bool     filterEnabled[MAX_NUM_CC_ALF_FILTERS],
                                 uint8_t  mapFilterIdxToFilterIdc[MAX_NUM_CC_ALF_FILTERS + 1],
                                 uint8_t &ccAlfFilterCount);
  void deriveCcAlfFilter( CodingStructure& cs, ComponentID compID, const PelUnitBuf& orgYuv, const PelUnitBuf& tempDecYuvBuf, const PelUnitBuf& dstYuv );
  std::vector<int> getAvailableCcAlfApsIds(CodingStructure& cs, ComponentID compID);
  void xSetupCcAlfAPS( CodingStructure& cs );
  void countLumaSwingGreaterThanThreshold(const Pel* luma, int lumaStride, int height, int width, int log2BlockWidth, int log2BlockHeight, uint64_t* lumaSwingGreaterThanThresholdCount, int lumaCountStride);
  void countChromaSampleValueNearMidPoint(const Pel* chroma, int chromaStride, int height, int width, int log2BlockWidth, int log2BlockHeight, uint64_t* chromaSampleCountNearMidPoint, int chromaSampleCountNearMidPointStride);
  void getFrameStatsCcalf(ComponentID compIdx, int filterIdc);
  void initDistortionCcalf();
};

#endif
