
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxh.h    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0 $

=========================================================================*/

extern int glbVerboseOutput ;
extern float glbfmDeriMetric[16] ;  

#ifndef zxh_h
#define zxh_h
#include "zxhDllExport.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#define HAS_BUILD_RELEASE
namespace zxh
{
///
/// \namespace zxh
/// \brief: proivde methods/functions
///
/// \ingroup zxh
///

#ifndef zxhuchar
#define zxhuchar unsigned char
#endif
#ifndef zxhushort
#define zxhushort unsigned short int
#endif
#ifndef zxhlfloat
#define zxhlfloat long double
#endif
	
#ifndef float8
#define float8 long double
#endif
#ifndef float4
#define float4 float
#endif

#ifndef ZXH_Floatvtk // is double
#define ZXH_Floatvtk double
#endif


#ifndef ZXH_ImageDimensionMax
#define ZXH_ImageDimensionMax 4
#endif
	#ifndef ImageDimensionMax
#define ImageDimensionMax 4
#endif
#ifndef ZXH_PI
#define ZXH_PI 3.14159f
#endif
#ifndef ZXH_GrayLevelOfSurface
#define ZXH_GrayLevelOfSurface 1
#endif
#ifndef ZXH_Protect_Small_Float
#define ZXH_Protect_Small_Float 100000
#endif
#ifndef ZXH_FloatPrecision
#define ZXH_FloatPrecision 1.0e-20f
#define ZXH_FloatInfinitesimal ZXH_FloatPrecision //1.0e-20f
#define ZXH_FloatEpsilon ZXH_FloatPrecision       //1.0e-20f
#endif
#ifndef ZXH_LFloatPrecision
#define ZXH_LFloatPrecision 1.0e-200
#endif
#ifndef ZXH_Background
#define ZXH_Background 0
#endif
#ifndef ZXH_EMPTY  
#define ZXH_EMPTY -32768
#endif
#ifndef ZXH_Foreground
#define ZXH_Foreground 1 //Now for both char and short
#endif
#ifndef ZXH_DefaultOptimiseStep 
#define ZXH_DefaultOptimiseStep 50
#endif
#ifndef ZXH_BSPLINE
#define ZXH_BSPLINE 1024
#define ZXH_BSPLINE_PLUS 1025
#endif

#ifndef ZXH_InfiniteLargeFloat
#define ZXH_InfiniteLargeFloat 1.0e20f
#endif

#ifndef ZXH_MaxRegistrations
#define ZXH_MaxRegistrations 20
#endif
	
#ifndef ZXH_MaxNumOfLabelOfImage
#define ZXH_MaxNumOfLabelOfImage 16
#endif

	
#ifndef ZXH_SegMultiSequenceNumOfLabelMax
#define ZXH_SegMultiSequenceNumOfLabelMax 8
#endif
#ifndef ZXH_SegMultiSequenceNumOfSeqMax
#define ZXH_SegMultiSequenceNumOfSeqMax 8
#endif

#ifndef ZXH_DefaultAccumulatedDecreaseStepsForStop
#define ZXH_DefaultAccumulatedDecreaseStepsForStop 6 //2
#endif

#ifndef ZXH_PHASE_FLOAT_INT
#define ZXH_PHASE_FLOAT_INT 1000.0
#endif

#ifndef ZXH_Mask
#define ZXH_Mask 0
#endif
#ifndef ZXH_WeightIntRatio
#define ZXH_WeightIntRatio 100.0
#endif

#ifndef ZXH_MATRIX44
#define ZXH_MATRIX44
typedef struct { float m[4][4] ;}  zxhmat44 ;
typedef zxhmat44 mat44 ;
#endif

//#ifndef ZXH_FFD_PRECOMPUTE_SIZE
//#define ZXH_FFD_PRECOMPUTE_SIZE 100
//#endif

#ifndef ZXH_LOCALAFFINE_MACRO
#define ZXH_LocalAffineMaxNumber 256
#define ZXH_LocalAffineExponent 1 //2 //1.5 2, 4
#define ZXH_LocalAffineFacetWeighting 0.5f //1.0/6=0.16667, 1.0/3=0.333333, 1/2, 
// Weight=1+1*( ((1/0.5)^e-1) /0.5 ), if NOT pre-compute weight , then is Distance =0
#define ZXH_LocalAffineWithinRegionWeight ZXH_FloatPrecision
#endif

// verbose output:
	// 0 only output error, and resultant infoa
	// 1 include 0 and success info
	// 2 include 1 and debug info
enum ZXH_ComputeIntensityRange { ZXH_ComputeIntensityRange_NOACTION=0,
								 ZXH_ComputeIntensityRange_NORMAL=1,
								 ZXH_ComputeIntensityRange_MASK=2,
								 ZXH_ComputeIntensityRange_EXCLUDE2PERCENT=3, // 2% excluded but intensity range extended to add two extra bins for each side
								 ZXH_ComputeIntensityRange_MASKEXCLUDE2PERCENT=4,
								 ZXH_ComputeIntensityRange_PRESET=5
							};
enum TypeHistogramCompute {PDFType_ParzenWindow_BSpline=1, PDFType_Histogram=2, PDFType_PV=3, UniformVolumeHistogram} ;
enum TypeMutualInformation { MutualInformation=1, NormalizedMutualInformation=2, MixtureProbabilityModel=3, MinusConditionalEntropyOfTestGivenRef=4} ;
enum TypeSEMIKernelFunction { SEMIGaussian=-1, SEMIBSplineZeroOrder=0, SEMIBSpline=3, SEMIMixtureModel=10};
typedef TypeHistogramCompute     ZXH_TypeHistogramCompute;
typedef TypeMutualInformation    ZXH_TypeMutualInformation;
typedef TypeSEMIKernelFunction   ZXH_TypeSEMIKernelFunction;
#ifndef SpatialMutualInformationSum
#define SpatialMutualInformationSum float(-0.998)
#endif
#ifndef SpatialMutualInformationMagnitude
#define SpatialMutualInformationMagnitude float(-0.997)
#endif
//#ifndef SpatialMutualInformationSearchScope
//#define SpatialMutualInformationSearchScope 1 int(-smisigma * 2)
//#endif

#ifndef ZXH_NV_NOISELEVEL
#define ZXH_NV_NOISELEVEL float(0.02)
#endif

#ifndef BSplineHistogramPaddingBins
#define BSplineHistogramPaddingBins 4 //2*2
#define BSplineHistogramPaddingBinsFrom -1
#define BSplineHistogramPaddingBinsTo 2
#define BSplineHistogramPaddingBinsOffset 2 //BSplineHistogramPaddingBins/2
#endif


/// \return save string to filename 
ZXH_DLL_EXPORT bool SaveString2File(std::string str2save,std::string filename);
///
ZXH_DLL_EXPORT bool SaveString2FileAppend(std::string str2save,std::string filename);

/// \begin{string&stream operation}
/// \return modified string s with spaces trimmed from left
ZXH_DLL_EXPORT void trim_left(std::string & s) ;

/// \return modified string s with spaces trimmed from right
ZXH_DLL_EXPORT void trim_right(std::string& s) ;

/// \return modified string s with spaces trimmed from edges
ZXH_DLL_EXPORT void trim_both(std::string& s) ;

/// \return modified string to upper case
ZXH_DLL_EXPORT void case_upper(std::string& s);

/// \return modified string to lower case
ZXH_DLL_EXPORT void case_lower(std::string& s);
///
ZXH_DLL_EXPORT std::string str_get_lowercase(std::string s);
///
ZXH_DLL_EXPORT std::string str_get_uppercase(std::string s);

///
ZXH_DLL_EXPORT bool is_option( std::string s );
///
ZXH_DLL_EXPORT bool is_digital( std::string s );

/// \return sContent sComents of the strLine
ZXH_DLL_EXPORT void ParseStringLine(std::string& sContent,std::string& sComment, std::string& strLine);

/// \return find a line with searching string
ZXH_DLL_EXPORT bool FindStringFromStream(char*buffer,int buffersize,std::ifstream&ifs,char* substring);

/// \return find next line with sContent.empty()==false
ZXH_DLL_EXPORT bool NextContentLine(std::string& sContent,char*buffer,int buffersize,std::ifstream&ifs);

/// \return int from the stream
ZXH_DLL_EXPORT bool ScanInteger(std::ifstream & ifs,char*buffer,int ibuffersize,int& idata);

/// \return float scan from the stream
ZXH_DLL_EXPORT bool ScanFloat(std::ifstream&ifs,char*buffer,int ibuffersize,float& fdata);

/// \return string scan from stream
ZXH_DLL_EXPORT bool ScanString(std::ifstream&ifs,char*buffer,int ibuffersize,char*pdata);

/// \end{string&stream operation}

/// perform machine dependent byte swapping. Byte swapping is often used when reading or writing binary files.
ZXH_DLL_EXPORT void Swap2Bytes(void* p);
/// \return
ZXH_DLL_EXPORT void Swap4Bytes(void* p);
/// \return
ZXH_DLL_EXPORT void Swap8Bytes(void* p);
/// \return
ZXH_DLL_EXPORT void Swap10Bytes(void* p);
/// \return nbyte:2,4,8
ZXH_DLL_EXPORT void SwapBytes(void* p,unsigned char nbyte);
/// end of swap

/// \return extension file name from a string; and abc.nii.gz => nii.gz
ZXH_DLL_EXPORT std::string GetExtension(const std::string filename);
/// 
ZXH_DLL_EXPORT std::string GetFileExtension(const std::string filename) ;
///
ZXH_DLL_EXPORT std::string GetFileNameNoPath(const std::string filename);
///
ZXH_DLL_EXPORT std::string GetFileNamePath(const std::string filename);
/// no point+ext, e.g.: abc.ext -> abc, abc.nii/abc.nii.gz -> abc
ZXH_DLL_EXPORT std::string GetFileNameNoExtension(const std::string filename);

///\return abs value of float type
ZXH_DLL_EXPORT float absf( const float f ) ;
///\return round value of float type
ZXH_DLL_EXPORT int round( const float f ) ;
///
ZXH_DLL_EXPORT float maxf( const float f1, const float f2 ) ;
///\return min value of float type
ZXH_DLL_EXPORT float minf( const float f1, const float f2 );
///\return echo the application arguments
ZXH_DLL_EXPORT void echo_arguments( int argc, char* argv[] ) ;
///
ZXH_DLL_EXPORT void echo_verbose(int argc, char*argv[] );
//
ZXH_DLL_EXPORT void echo_zxh(int argc, char*argv[] );
///
ZXH_DLL_EXPORT void echo_zxh();
///
ZXH_DLL_EXPORT void NormaliseVector( float v[], int dimension ) ; 
	
///
template<typename T>
ZXH_DLL_EXPORT T MagnitudeOfVector( const T v[], const int dimension )
{
	T mag = 0 ;
	for( int idim=0; idim<dimension; ++idim )
		mag += v[idim]*v[idim] ;
	return (T) sqrt(mag) ;
} 
///
ZXH_DLL_EXPORT void VectorOP_Normalise( float v[], int dimension )	;
///
ZXH_DLL_EXPORT float VectorOP_Magnitude( const float v[], int dimension )	;
///
ZXH_DLL_EXPORT float VectorOP_DotProduct( const float v1[], const float v2[], int dimension );
///
ZXH_DLL_EXPORT void VectorOP_CrossProduct3D( const float A[], const float B[], float *AxB ) ;
///
ZXH_DLL_EXPORT float VectorOP_Distance( const float v1[], const float v2[], int dimension );
///
ZXH_DLL_EXPORT float VectorOP_DistanceSquare( const float v1[], const float v2[], int dimension );
///
ZXH_DLL_EXPORT float VectorOP_DistanceFirstOrderNorm( const float v1[], const float v2[], int dimension );
///v1 dotproduct v2 = cos(angle) mag(v1)*mag(v2)
ZXH_DLL_EXPORT float VectorOP_Cosine( const float v1[], const float v2[], int dimension ); 
///
ZXH_DLL_EXPORT void VectorOP_Substract( const float A[], const float B[], float *A_B, int dimension ) ;
///
ZXH_DLL_EXPORT void VectorOP_Add( const float A[], const float B[], float *AaB, int dimension ) ;
///
ZXH_DLL_EXPORT void VectorOP_Multiple( const float A[], const float f, float *Axf, int dimension ) ;

///
ZXH_DLL_EXPORT float VectorOP_Point2LineDistance( const float pnt[3], const float pL1[3], const float pL2[3], float bWithinLineSegment ) ; 

///
template<typename T>
T** MallacTwoDimArrayWithDefault( T** &p, int num1dim, int num2dim, T defaultvalue )
{
	if( p ) FreeTwoDimArray(p,num1dim);
	p = new T*[num1dim]; 
	for( int ix=0; ix<num1dim; ++ix )
	{
		p[ix] = new T[num2dim]; 
		for( int iy=0; iy<num2dim; ++iy )
			p[ix][iy] = defaultvalue;
	}
	return p; 
}
/// p[num1dim][num2dim]
template<typename T>
T** MallacTwoDimArray( T** &p, int num1dim, int num2dim )
{
	if( p ) FreeTwoDimArray(p,num1dim);
	p = new T*[num1dim]; 
	for( int i1=0; i1<num1dim; ++i1 ) 
		p[i1] = new T[num2dim];   
	return p;
}
template<typename T>
bool FreeTwoDimArray( T** &p, int num1dim )
{ 
	if( p )
	{
		for( int i1=0; i1<num1dim; ++i1 ) 
			if( p[i1] ) delete [] p[i1] ;
		delete [] p ;
		p = 0 ;
	}
	return true;
}
///
ZXH_DLL_EXPORT
template<typename T, typename T2>
bool VectorOP_MeanAndStd( T const * pData, int num, T2 &mean, T2 &standD )
{
    if (pData ==NULL || num<2 )
    {
        return false;
    }
    if( num<10000)
    {
        zxhlfloat lmean = 0.0 ;
        zxhlfloat lsqx = 0.0 ;
        for (int i=0; i< num; ++i)
        {
            lmean += pData[i];
            lsqx += pData[i]*pData[i];
        }
        lmean /= (zxhlfloat) num;

        mean = T2(lmean) ;
        standD = T2(sqrt( (lsqx-num*lmean*lmean)/(zxhlfloat) (num-1) ) );
    }
    else
    {
        zxhlfloat lmean = 0.0, temsum = 0.0 ;
        zxhlfloat lsqx = 0.0 , temsqx = 0.0 ;
        for (int i=0; i< num; ++i)
        {
            temsum += pData[i];
            temsqx += pData[i]*pData[i];
            if( i%10000 == 0 )
            {
                lmean += temsum ; temsum = 0 ;
                lsqx += temsqx ; temsqx = 0 ;
            }
        }
        lmean += temsum ; temsum = 0 ;
        lsqx += temsqx ; temsqx = 0 ;

        lmean /= (zxhlfloat) num;
        mean = T2(lmean) ;
        standD = T2(sqrt( (lsqx-num*lmean*lmean)/(zxhlfloat)(num-1) ) );
    }
    return true;
}

///
ZXH_DLL_EXPORT
template<typename T, typename T2>
bool VectorOP_Sum( T const * pData, int num, T2 &sum )
{
    if (pData ==NULL || num<2 )
    {
        return false;
    }
    if( num<10000)
    {
        zxhlfloat lsum = 0.0 ; 
        for (int i=0; i< num; ++i)
        {
            lsum += pData[i]; 
        } 
		sum=(T2)lsum;
    }
    else
    {
        zxhlfloat lsum = 0.0, temsum = 0.0 ; 
        for (int i=0; i< num; ++i)
        {
            temsum += pData[i]; 
            if( i%10000 == 0 )
            {
                lsum += temsum ; temsum = 0 ; 
            }
        } 
		sum=(T2)lsum;
    }
    return true;
}



///
ZXH_DLL_EXPORT bool FileExist(const std::string filename);

///\return gaussian
ZXH_DLL_EXPORT template <class type> type Gaussian( const type f, const type u, const type sigma )  
{ // 0.39894228 = 1/sqrt(2*pi)
	if( sigma> 0 )
		return (0.39894228/sigma *exp(-0.5*(f-u)*(f-u)/(sigma*sigma)));
	else std::cerr<<"error: sigma=0 in Gaussian function\n";
	return sigma ; 
};
///\return Ngaussian(f) = Gaussian(f,0,1)/Gaussian(0,0,1)     sqrt(2*pi)
/// NGaussian(0)=1;  NGaussian(1)=0.6065; Ngaussian(2)=0.1353; NGaussian(3)=0.0111
ZXH_DLL_EXPORT template <class type> type NGaussian( const type f )
{
	return (exp(-0.5*f*f ));
}; 
///\return  
ZXH_DLL_EXPORT template <class type> type NGaussianTruncate( const type f, const float fTruncateTimesOfSigma )
{
	if( abs(f) >= fTruncateTimesOfSigma && fTruncateTimesOfSigma>0 ) return (exp(-0.5*fTruncateTimesOfSigma*fTruncateTimesOfSigma )) ;
	return (exp(-0.5*f*f ));
}; 
///\return gaussian
ZXH_DLL_EXPORT template <class type> type GaussianTruncate( type f, const type u, const type sigma, const float fTruncateTimesOfSigma )  
{ // 0.39894228 = 1/sqrt(2*pi)
	if( sigma> 0 )
	{
		type value = abs((f-u)/sigma) ; 
		if(  value< fTruncateTimesOfSigma || fTruncateTimesOfSigma<=0 )  
			return (0.39894228/sigma *exp(-0.5*value*value));
		else 
			return (0.39894228/sigma *exp(-0.5*fTruncateTimesOfSigma*fTruncateTimesOfSigma ));
	}
	else std::cerr<<"error: sigma=0 in Gaussian function\n";
	return sigma ; 
};
///\return gaussian derivative
ZXH_DLL_EXPORT template <class type> type GaussianTruncateDerivative( type f, const type u, const type sigma, const float fTruncateTimesOfSigma )  
{ // 0.39894228 = 1/sqrt(2*pi)
	if( sigma> 0 )
	{
		type value = (f-u)/sigma ; 
		if(  abs(value)< fTruncateTimesOfSigma || fTruncateTimesOfSigma<=0 )  
			return (0.39894228/sigma *exp(-0.5*value*value))  * (-value / sigma) ;
		else 
			return (0.39894228/sigma *exp(-0.5*fTruncateTimesOfSigma*fTruncateTimesOfSigma ))  * (-fTruncateTimesOfSigma / sigma);
	}
	else std::cerr<<"error: sigma=0 in Gaussian function\n";
	return sigma ; 
};
///
ZXH_DLL_EXPORT
template <class type>
bool InclusiveBetween( type value, type from, type to )
{
	if( value>=from && value <= to )
		return true ;
	return false ;
};
///\return abs value of  the type
ZXH_DLL_EXPORT
template <class TYPE>
TYPE abs( TYPE f )
{
	if( f>=0 ) return f ;
	return -f ;
};
///\return  
ZXH_DLL_EXPORT
template <class TYPE>
void ExchangeTwoValues( TYPE &f1, TYPE &f2 )
{
	TYPE tmp = f2 ; 
	f2 = f1 ; 
	f1 = tmp ;  
}
ZXH_DLL_EXPORT
template <class TYPE>
void CorrectFromTo( TYPE *From, TYPE *To, int idim )
{
	for( int id=0; id<idim; ++id ) 
		if( From[id] > To[id] )
			ExchangeTwoValues( From[id], To[id] ) ;  
}
///\return distance squared of two coordinates
ZXH_DLL_EXPORT
template <class TYPE>
TYPE SqaredDistance( TYPE *co1, TYPE *co2, int Dimension )
{
	TYPE dis = 0;
	for(int i=0;i<Dimension;++i)
		dis += (co1[i]-co2[i])*(co1[i]-co2[i]);
	return dis;
}
///
ZXH_DLL_EXPORT void MatrixIdentity(float *pMatDes, const int iMatrixLength ) ;
///
ZXH_DLL_EXPORT void SetMatrix(float *pMatDes, const float* pMatSrc1, const int iMatrixLength ) ;
///
ZXH_DLL_EXPORT void MultiplyMatrix(float *pMatDes, const float* pMatSrc1, const float* pMatSrc2, const int iMatrixLength ) ;
///
ZXH_DLL_EXPORT void MultiplyMatrixVector(float *pVecDes, const float* pMatSrc, const float* pVecSrc, const int iMatrixLength ) ;
///
ZXH_DLL_EXPORT float DeterminentMatrix3D(const float *pMat, const int iMatrixRowLength ) ;
/// assume the arrays are both 3*iMatrixRowLength
ZXH_DLL_EXPORT void InvertMatrix3D(float* pInvMat, const float* pForMat, const int iMatrixRowLength);
/* matrix op */
///
/// assume the arrays are both 3*iMatrixRowLength
ZXH_DLL_EXPORT void MatrixOP_InvertMatrix3D(float* pInvMat, const float* pForMat, const int iMatrixRowLength);
///
ZXH_DLL_EXPORT void MatrixOP_InvertMatrix4D3D( zxhmat44 &inv, const zxhmat44 &forward ) ;
/// 3x4 matrix
ZXH_DLL_EXPORT void MatrixOP_InvertMatrix4D3D(float* pInvMat, const float* pForMat, const int iMatrixRowLength);
///
ZXH_DLL_EXPORT void MatrixOP_SetIdentity(float *pMatDes, const int iMatrixLength )  ;
///
ZXH_DLL_EXPORT void MatrixOP_SetMatrix(float *pMatDes, const float* pMatSrc1, const int iMatrixLength )  ;
///
ZXH_DLL_EXPORT void MatrixOP_Multiply(float *pMatDes, const float* pMatSrc1, const float* pMatSrc2, const int iMatrixLength ) ; ;
///
ZXH_DLL_EXPORT void MatrixOP_MultiplyVector(float *pVecDes, const float* pMatSrc, const float* pVecSrc, const int iMatrixLength )   ; 
///
ZXH_DLL_EXPORT float MatrixOP_DeterminentMatrix3D(const float *pMat, const int iMatrixRowLength )  ;
 
/// \return BSpline
//template <class TYPE>
//TYPE BSpline ( TYPE u )
ZXH_DLL_EXPORT float BSpline(float u) ;
/// \return  
ZXH_DLL_EXPORT float BSplineDerivative(float u) ;
///
ZXH_DLL_EXPORT float BSplineSecondDerivative ( float u ) ;
/// \return BSpline
ZXH_DLL_EXPORT float BSplinei(int iOrd,float u) ;
///
ZXH_DLL_EXPORT float DerivativeOfBSplinei(int iOrd, float u) ;
///
ZXH_DLL_EXPORT float SecondDerivativeOfBSplinei(int iOrd,float u) ;
 

//************ Macro ***************//
#define MIN(f1,f2) ((f1)>=(f2)?(f2):(f1))
#define MAX(f1,f2) ((f1)<=(f2)?(f2):(f1))
#define ZXH_MIN(f1,f2) ((f1)>=(f2)?(f2):(f1))
#define ZXH_MAX(f1,f2) ((f1)<=(f2)?(f2):(f1))
#define ZXH_ABS( f ) ((f)>0?(f):-(f))
#define ZXH_ROUND( f ) ((f)>=0?int((f)+0.5):-int(-(f)+0.5)) 
#define ZXH_InclusiveBetween( value, from, to )  ((value)>=(from) && (value) <= (to) ? true:false)
 
/// ***** macro for metric type conversion *****   
#define ISMITYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_MI_")!=NULL)

#define ISSEMITYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_SEMI_")!=NULL)

#define ISCONTAINERTYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_CONTAINER_")!=NULL)

#define ISCONSTRAINTSHAPETYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_CONSTRAINT_SHAPE_")!=NULL)

#define ISPHASETYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_PHASE_")!=NULL)

#define ISANALYTICALTYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_ANALYTICAL_")!=NULL)

#define ISCONSTDISCRETEPATHTYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_CONSTDISCRETEPATH_")!=NULL)

#define ISMICNSTTYPE(pMetric) \
(strstr(pMetric->GetMetricType().c_str(),"_MI_")!=NULL&&strstr(pMetric->GetMetricType().c_str(),"_CNST_")!=NULL)

#define ISRMSTYPE(pMetric)\
(strstr(pMetric->GetMetricType().c_str(),"_RMS_")!=NULL)

#define ISCCTYPE(pMetric)\
(strstr(pMetric->GetMetricType().c_str(),"_CC_")!=NULL)

#define ISNVTYPE(pMetric)\
(strstr(pMetric->GetMetricType().c_str(),"_NV_")!=NULL)

#define ISFFDMETRIC(pMetric)\
(strstr(pMetric->GetMetricType().c_str(),"_FFD_")!=NULL)

#define ISFFD2METRIC(pMetric)\
(strstr(pMetric->GetMetricType().c_str(),"_FFD2_")!=NULL)

#define ISLOCALAFFINEMETRIC(pMetric)\
(strstr(pMetric->GetMetricType().c_str(),"_LAFFS_")!=NULL)

#define SetFFDMetric2Metric(FFDMetricObject,MetricObject) \
{	if(ISMITYPE(MetricObject)) FFDMetricObject = ((zxhMetricMIFFD*) MetricObject); \
	if(ISPHASETYPE(MetricObject))FFDMetricObject = ((zxhMetricPhaseFFD*) MetricObject); \
	if(ISRMSTYPE(MetricObject))FFDMetricObject = ((zxhMetricRMS*) MetricObject); \
	if(ISNVTYPE(MetricObject)) FFDMetricObject = ((zxhMetricNV*) MetricObject); \
	if(ISCCTYPE(MetricObject)) FFDMetricObject = ((zxhMetricCCFFD2*) MetricObject); \
};

//if(ISPHASETYPE(MetricObject)&&ISANALYTICALTYPE ) LAMetricObject = (zxhMetricPhaseLocalAffines*) MetricObject ; 

#define SetLocalAffinesMetric2Metric(LAMetricObject,MetricObject) \
{	if( ISMITYPE(MetricObject)&&ISANALYTICALTYPE(MetricObject) ) LAMetricObject = (zxhMetricMILocalAffineAnalytical*) MetricObject ; \
	else if( ISCCTYPE(MetricObject)&&ISANALYTICALTYPE(MetricObject) ) LAMetricObject = (zxhMetricCCLocalAffineAnalytical*) MetricObject ; \
	else if( ISMITYPE(MetricObject)&& ! ISANALYTICALTYPE(MetricObject) ) LAMetricObject = (zxhMetricMILocalAffine*) MetricObject ; \
	else if( ISCCTYPE(MetricObject)&& ! ISANALYTICALTYPE(MetricObject) ) LAMetricObject = (zxhMetricCCLocalAffine*) MetricObject ; \
	else if(ISPHASETYPE(MetricObject)) LAMetricObject = (zxhMetricPhaseLocalAffines*) MetricObject ; \
	else if(ISMITYPE(MetricObject)) LAMetricObject = (zxhMetricMILocalAffines*) MetricObject ; \
	else if(ISCONSTRAINTSHAPETYPE(MetricObject)) LAMetricObject = (zxhMetricConstraintLocalAffineShape*) MetricObject ; } ;
 
/// ***** macro for metric type conversion end *****   
#ifdef __WIN32__
	#define ZXH_IS_WINDOWS 
#else
	#ifdef _WIN32
		#define ZXH_IS_WINDOWS 
	#else
		#ifdef _MSC_VER
			#define ZXH_IS_WINDOWS 
		#else
#define ZXH_IS_NOT_WINDOWS 
		#endif
	#endif
#endif

#ifdef ZXH_IS_WINDOWS
#include<Windows.h>
#else
#include <unistd.h>
#endif

ZXH_DLL_EXPORT bool IsWindowsOS() ;
///

ZXH_DLL_EXPORT size_t GetSizeOfMemoryMB();

/// a line, three float delimit by tab for coronary/ostiia class in the future
ZXH_DLL_EXPORT bool ReadOstiumCoordinateFromTxt( std::string filename, float * wco ) ; 
/// append a line, three float delimit by tab
ZXH_DLL_EXPORT bool AppendSaveOstiumCoordinateToTxt( float * wco, std::string filename ) ; 


}//end of zxh
//#ifdef __cpluscplus
//	}
//#endif

#endif //zxh_h

