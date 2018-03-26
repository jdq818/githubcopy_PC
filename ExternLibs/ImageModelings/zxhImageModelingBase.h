
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhImageModelingBase.h    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0 $

=========================================================================*/


#ifndef zxhImageModelingBaseT_H
#define zxhImageModelingBaseT_H
#include "zxhImageData.h"
#include "zxhTransformBase.h"

namespace zxh
{
#ifndef ImageModeling_XYZ_Type
#define ImageModeling_XYZ_Type
	enum ImageModeling_XYZ
	{
		ImageModeling_X=0,
		ImageModeling_Y=1,
		ImageModeling_Z=2,
		ImageModeling_T=3,
		ImageModeling_5D=4
	};
#endif
}

///
/// \class zxhImageModelingBaseT
/// \brief  image modeling base class
///		 use for getting pixel, gradient,
///		 different claas stands for different interpolation
///		 Here: pixelvalue and gradient type should be float type
///  for registration image modeling, the pixeltype is default ZXH_PixelTypeDefault(short)
///  the template only used in interpolate the bspline ffd
///  implemented derived sub-class, Nearest, Linear, BSpline(1D,2D,3D,3D+1D)
///
///	\ingroup zxhImageModelings
///
//class zxhTransformBase ;
template<class PixelType=ZXH_PixelTypeDefault,class DerivativeType=float>
class zxhImageModelingBaseT
{
public:
	/// typedef
	typedef float	PixelFloatType;
	/// typedef
	typedef float	DerivativeFloatType;

	/// constructor
	zxhImageModelingBaseT(void)	{m_pImage=0;};

	///\return
	~zxhImageModelingBaseT(void){};
	///
	virtual void SetImage(zxhImageDataT<PixelType>* pImage)					{m_pImage=pImage; };

	///\return
	virtual zxhImageDataT<PixelType>* GetImage()	const						{return m_pImage;};

	/* ************** necessary consider re-implementation in derived classes ******************* */
	;
	///
	virtual zxhImageModelingBaseT<PixelType,DerivativeType>*	CloneTo(zxhImageModelingBaseT<PixelType,DerivativeType>* &pRet)const
	{
		if(pRet==0)
			pRet=new zxhImageModelingBaseT<PixelType,DerivativeType>;
		pRet->SetImage(m_pImage);
		return pRet;
	}; 
	///
	virtual inline PixelFloatType	GetPixelFloatValueWithCheckClosest(float fx,float fy,float fz=0,float ft=0)const
	{ 	return this->m_pImage->GetPixelGreyscaleClosest( zxh::round(fx), zxh::round(fy), zxh::round(fz), zxh::round(ft)) ; };
	/// without boundary check, assume [] \in 0~size-1, diff from different modelling
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy)const
	{	return m_pImage->GetPixelGreyscale( zxh::round(fx ), zxh::round(fy ) ) ;};
	/// without boundary check
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy,float fz)const
	{	return m_pImage->GetPixelGreyscale( zxh::round(fx ), zxh::round(fy ), zxh::round(fz ) ); };
	/// without boundary check
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy,float fz,float ft)const
	{	return m_pImage->GetPixelGreyscale( zxh::round(fx ), zxh::round(fy ), zxh::round(fz ), zxh::round(ft ) ); };  ;
	///
	virtual void GetPixelGradientByWorld(DerivativeType pwGrad[ZXH_ImageDimensionMax],float wx,float wy=0,float wz=0,float wt=0) const;
	/// Image Gradient of nearest modeling
	/// 2015-7-10 zxh revise: gradient should be dI/dWx, where Wx is world coordinate, but not image coordinate
	virtual void GetPixelGradient(DerivativeType pwGrad[ZXH_ImageDimensionMax],float x,float y=0,float z=0,float t=0)const
	{
		if( m_pImage == 0 )
			return ;
		float w[] = {x,y,z,t} ;
		m_pImage->ImageToWorld( w ) ;
		return GetPixelGradientByWorld( pwGrad, w[0],w[1],w[2],w[3] ) ;
	};
 
	/* ************** necessary consider re-implementation in derived classes -- end ******************* */
	
	///  closest
	virtual inline PixelFloatType	GetPixelFloatValueWithCheck(float fx,float fy,float fz=0,float ft=0)const
	{	return GetPixelFloatValueWithCheckClosest(fx,fy,fz,ft) ; };

	///
	///\return
	virtual DerivativeType GetPixelGradientMagnitude( float x, float y=0, float z=0, float t=0 )const
	{
		DerivativeType pGrad[ZXH_ImageDimensionMax] ;
		this->GetPixelGradient( pGrad, x, y, z, t ) ;
		return zxh::MagnitudeOfVector( pGrad, m_pImage->GetDimension() ) ;
	};
	///
	virtual DerivativeType GetPixelGradientMagnitudeMax() const;
	///
	virtual DerivativeType GetPixelGradientMagnitudeMin() const;

	/// partial differential on physical coordinate using image grid index
	virtual DerivativeType	GetPixelPhysicalPartialDiff(zxh::ImageModeling_XYZ dim,float x,float y=0,float z=0,float t=0) const;

	//virtual DerivativeType	GetPixel2OrderDerivative(ImageModeling_XYZ firstDim,ImageModeling_XYZ secondDim,float x,float y,float z=0)=0;

	/* * Image operation * */

	/// resize image
	virtual bool ResizeImage(const int size[], const float spacing[] );
	/// resize using spacing
	virtual bool ResizeImage(const  float spacing[] ) ;
	///
	virtual bool ResizeImage(const  zxhImageInfo *imageinfo );
	/// using closest for outliers
	virtual bool ResizeImageUsingClosestForOutliers(const  zxhImageInfo *imageinfo );
	/// fill background ZXH_Background
	virtual bool FowardTransformImageUsingInverse(	zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const;
	/// use closest points for outliers
	virtual bool FowardTransformImageUsingInverseClosest(	zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const;
	/// fill background ZXH_Background
	virtual bool BackwardTransformImage(zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const;
	///
	virtual bool BackwardTransformImageMirrorOutlier(zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const;
	///
	virtual bool BackwardTransformImageClosest(zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const;
	/// 1) if m_pImage is 2D, then make it 3D by adding two slices with spacing sp3d ;
	//virtual bool BackwardTransformImage2DTo3DBg( zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform, float sp3d, bool nFillBg, PixelType vBg ) ;
	;
	 
	/// \return backward transform using transform list
	/// \para 
	///  bBackground = false, donot use background value but use closest pixel points to fill
	virtual bool BackwardTransformImageByList(zxhImageDataT<PixelType>*pSave, 
												  const zxhTransformBase*const*pList,int nTransform, // const pList makes the function hard to use
												  bool bBackground=true, PixelType background=ZXH_Background) const; 

	/* *********************** none virtual method, do NOT inheritate them ************** */
	;
	/// interpolation pixel value, with outlier using mirror reflect and then GetPixelFloatValueWithoutCheck
	PixelFloatType	GetPixelFloatValueWithCheckMirror(float fx,float fy,float fz=0,float ft=0)const
	{
		int nx,ny,nz,nt ;
		m_pImage->GetImageSize( nx,ny,nz,nt ) ;
		if(fx<0)
			return this->GetPixelFloatValueWithCheckMirror( -1*fx, fy, fz, ft ) ;
		if(fx>static_cast<float>(nx-1))
			return this->GetPixelFloatValueWithCheckMirror( static_cast<float>(nx*2-2)-fx, fy, fz, ft ) ;
		if(fy<0)
			return this->GetPixelFloatValueWithCheckMirror( fx, -1*fy, fz, ft ) ;
		if(fy>static_cast<float>(ny-1))
			return this->GetPixelFloatValueWithCheckMirror( fx, static_cast<float>(ny*2-2)-fy, fz, ft ) ;
		if(fz<0)
			return this->GetPixelFloatValueWithCheckMirror( fx, fy, -1*fz, ft ) ;
		if(fz>static_cast<float>(nz-1))
			return this->GetPixelFloatValueWithCheckMirror( fx, fy, static_cast<float>(nz*2-2)-fz, ft ) ;
		if(ft<0)
			return this->GetPixelFloatValueWithCheckMirror( fx, fy, fz, -1*ft ) ;
		if(ft>static_cast<float>(nt-1))
			return this->GetPixelFloatValueWithCheckMirror( fx, fy, fz, static_cast<float>(nt*2-2)-ft ) ;

		return this->GetPixelFloatValueWithoutCheck( fx, fy, fz, ft ) ;
	};
	///
	PixelFloatType	GetPixelFloatValueWithCheckBk(PixelFloatType bk, float fx,float fy,float fz=0,float ft=0)const
	{   
		if( m_pImage->InsideImageSafe(fx,fy,fz,ft) ) 
			return this->GetPixelFloatValueWithoutCheck( fx, fy, fz, ft ) ;
		return bk ;
	};
	/// get the pixel by world coordinate
	PixelFloatType	GetPixelFloatValueWithCheckClosestByWorld(float fx,float fy,float fz=0,float ft=0)const
	{
		return GetPixelFloatValueWithCheckByWorld(fx,fy,fz,ft) ;
	}
	/// get the pixel by world coordinate
	PixelFloatType	GetPixelFloatValueWithCheckByWorld(float fx,float fy,float fz=0,float ft=0)const
	{
		if( m_pImage == 0 )
			return 0;
		float fv[] = {fx,fy,fz,ft} ;
		m_pImage->WorldToImage( fv ) ;
		return this->GetPixelFloatValueWithCheck( fv[0],fv[1],fv[2],fv[3] );
	};
	/// get the pixel by physical coordinate
	PixelFloatType	GetPixelFloatValueWithCheckMirrorByWorld(float fx,float fy,float fz=0,float ft=0)const
	{
		if( m_pImage == 0 )
			return 0;
		float fv[] = {fx,fy,fz,ft} ;
		m_pImage->WorldToImage( fv ) ;
		return this->GetPixelFloatValueWithCheckMirror( fv[0],fv[1],fv[2],fv[3] );
	};
	/// get the pixel by physical coordinate
	PixelFloatType	GetPixelFloatValueWithCheckBkByWorld(PixelFloatType bk, float fx,float fy,float fz=0,float ft=0)const
	{
		if( m_pImage == 0 )
			return 0;
		float fv[] = {fx,fy,fz,ft} ;
		m_pImage->WorldToImage( fv ) ;
		return this->GetPixelFloatValueWithCheckBk(bk, fv[0],fv[1],fv[2],fv[3] );
	};
	/* *********************** none virtual method, do NOT inheritate them ************** */

protected:
	/// \return
	zxhImageDataT<PixelType>	*	m_pImage;
};

template<class PixelType,class DerivativeType>
DerivativeType
zxhImageModelingBaseT<PixelType, DerivativeType>::GetPixelGradientMagnitudeMax()const
{
	DerivativeType max = 0, mag;
	DerivativeType gd[ZXH_ImageDimensionMax]={0,0,0,0} ;
	int size[] = {1,1,1,1} ;
	m_pImage->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	int dimension = m_pImage->GetDimension() ;

	for( float it = 0 ; it <= size[3]-1 ; ++it )
	for( float iz = 0 ; iz <= size[2]-1 ; ++iz )
	for( float iy = 0 ; iy <= size[1]-1 ; ++iy )
	for( float ix = 0 ; ix <= size[0]-1 ; ++ix )
	{
		this->GetPixelGradient( gd, ix, iy, iz, it ) ;
		mag = zxh::MagnitudeOfVector( gd, dimension ) ;
		if( max < mag )
			max = mag ;
	}
	return max ;
}

template<class PixelType,class DerivativeType>
DerivativeType
zxhImageModelingBaseT<PixelType, DerivativeType>::GetPixelGradientMagnitudeMin()const
{
	DerivativeType min , mag;
	DerivativeType gd[ZXH_ImageDimensionMax]={0,0,0,0} ;
	int size[] = {1,1,1,1} ;
	m_pImage->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	int dimension = m_pImage->GetDimension() ;

	this->GetPixelGradient( gd, 0,0,0,0 ) ;
	min = zxh::MagnitudeOfVector( gd, dimension ) ;

	for( float it = 0 ; it <= size[3]-1 ; ++it )
	for( float iz = 0 ; iz <= size[2]-1 ; ++iz )
	for( float iy = 0 ; iy <= size[1]-1 ; ++iy )
	for( float ix = 0 ; ix <= size[0]-1 ; ++ix )
	{
		this->GetPixelGradient( gd, ix, iy, iz, it ) ;
		mag = zxh::MagnitudeOfVector( gd, dimension ) ;
		if( min > mag )
			min = mag ;
	}
	return min ;
}

template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::ResizeImage(const  zxhImageInfo *imageinfo )
{
	if( imageinfo == 0 ) return false ;
	if( imageinfo->SameDimSizeSpacingAs( this->m_pImage->GetImageInfo() ) == true &&
		imageinfo->SameOrientationAs( this->m_pImage->GetImageInfo() ) == true )
		return true ;
	zxhImageDataT<PixelType> img ;
	img.NewImage( imageinfo->Dimension, imageinfo->Size, imageinfo->Spacing, imageinfo ) ;
	this->BackwardTransformImage( &img, 0 ) ;
	this->m_pImage->NewImage( imageinfo->Dimension, imageinfo->Size, imageinfo->Spacing, imageinfo ) ;
	this->m_pImage->MemCopyImageDataFrom( &img ) ;
	return true ;
}

template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::ResizeImageUsingClosestForOutliers(const  zxhImageInfo *imageinfo )
{
	if( imageinfo == 0 ) return false ;
	if( imageinfo->SameDimSizeSpacingAs( this->m_pImage->GetImageInfo() ) == true &&
		imageinfo->SameOrientationAs( this->m_pImage->GetImageInfo() ) == true )
		return true ;
	zxhImageDataT<PixelType> img ;
	img.NewImage( imageinfo->Dimension, imageinfo->Size, imageinfo->Spacing, imageinfo ) ;
	this->BackwardTransformImageClosest( &img, 0 ) ;
	this->m_pImage->NewImage( imageinfo->Dimension, imageinfo->Size, imageinfo->Spacing, imageinfo ) ;
	this->m_pImage->MemCopyImageDataFrom( &img ) ;
	return true ;
}


template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::ResizeImage(const  float spacing[] )
{
	if( this->m_pImage == 0 || this->m_pImage->GetNumberOfPixels()<= 1 )
		return false ;
	int size[] = {1,1,1,1} ;
	this->m_pImage->GetImageInfo()->GetSizeUsingSpacing( spacing, size ) ;
	return this->ResizeImage( size, spacing ) ;
};
template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::ResizeImage(const int size[], const float spacing[] )
{
	if( this->m_pImage == 0 || this->m_pImage->GetNumberOfPixels()<= 1 )
		return false ;
	bool same = true ;

	int imagesize[] = {1,1,1,1}, imagedimension = this->m_pImage->GetDimension() ;
	float imagespacing[] = {1,1,1,1} ;

	this->m_pImage->GetImageSize( imagesize[0], imagesize[1], imagesize[2], imagesize[3] ) ;
	this->m_pImage->GetImageSpacing( imagespacing[0], imagespacing[1], imagespacing[2], imagespacing[3] ) ;

	for( int idim = 0 ; idim < imagedimension; ++idim )
		if( imagesize[idim] != size[idim] ||
			imagespacing[idim] != spacing[idim] )
			same=false ;
	if( same ) return true ;
	zxhImageInfo imageinfo ;
	this->m_pImage->GetImageInfo( & imageinfo ) ;

	for( int i=0; i<ZXH_ImageDimensionMax; ++i )
	{
		imageinfo.Spacing[i] = spacing[i] ;
		imageinfo.Size[i] = size[i] ;
	}
	imageinfo.UpdateOrientationInfo( -2 ) ;

	zxhImageDataT<PixelType> copy ;
	copy.NewImage( m_pImage->GetDimension(), size, spacing, &imageinfo );
	copy.SetImageInfo( &imageinfo ) ; // to copy the image name info etc.

 	PixelType value ;

	switch(	m_pImage->GetDataType() )
	{
	case GIPL_CHAR:
	case GIPL_U_CHAR:
	case GIPL_SHORT:
	case GIPL_U_SHORT:
	case GIPL_U_INT:
	case GIPL_INT:
		//char uchar short zxhushort int uint
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float w[] = { float(ix), float(iy), float(iz), float(it) } ;
			copy.ImageToWorld( w ) ;
			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckByWorld( w[0], w[1], w[2], w[3] ) +0.5);
			copy.SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}
		break;
	case GIPL_FLOAT:
	case GIPL_DOUBLE:
	default:
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float w[] = { float(ix), float(iy), float(iz), float(it) } ;
			copy.ImageToWorld( w ) ;
			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckByWorld( w[0], w[1], w[2], w[3] ) );
			copy.SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}

		break;
	};
 	copy.CloneTo( this->m_pImage ) ;
	return true ;
};

template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::FowardTransformImageUsingInverseClosest(	
	zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const
{
	if( pSave==0 || pSave->IsEmpty() )
		return false ;
	int size[] = {1,1,1,1} ;
	pSave->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	PixelType value ;
	PixelType background = ZXH_Background ;
	PixelType vmin = this->m_pImage->GetPixelGreyscaleMin() ;
	if( background > vmin )
		background = vmin ;

	switch(	m_pImage->GetDataType() )
	{
	case GIPL_CHAR:
	case GIPL_U_CHAR:
	case GIPL_SHORT:
	case GIPL_U_SHORT:
	case GIPL_U_INT:
	case GIPL_INT:
		//char uchar short zxhushort int uint
		for( int ix=0; ix<size[0]; ++ix )
		for( int iy=0; iy<size[1]; ++iy )
		for( int iz=0; iz<size[2]; ++iz )
		for( int it=0; it<size[3]; ++it )
		{
			float wto[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wto ) ;
			float wfrom[] = { wto[0],wto[1],wto[2],wto[3] } ;
			if( pTransform )
				pTransform->InverseTransformWorld( wto, wfrom );

			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckClosestByWorld( wfrom[0], wfrom[1], wfrom[2], wfrom[3] ) +0.5);

			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}
		break;
	case GIPL_FLOAT:
	case GIPL_DOUBLE:
	default:
		for( int ix=0; ix<size[0]; ++ix )
		for( int iy=0; iy<size[1]; ++iy )
		for( int iz=0; iz<size[2]; ++iz )
		for( int it=0; it<size[3]; ++it )
		{
			float wto[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wto ) ;
			float wfrom[] = { wto[0],wto[1],wto[2],wto[3] } ;
			if( pTransform )
				pTransform->InverseTransformWorld( wto, wfrom );

			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckClosestByWorld( wfrom[0], wfrom[1], wfrom[2], wfrom[3] ) );
			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}

		break;
	};


	return true ;
}

template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::FowardTransformImageUsingInverse(	
	zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const
{
	if( pSave==0 || pSave->IsEmpty() )
		return false ;
	int size[] = {1,1,1,1} ;
	pSave->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	PixelType value ;
	PixelType background = ZXH_Background ;
	PixelType vmin = this->m_pImage->GetPixelGreyscaleMin() ;
	if( background > vmin )
		background = vmin ;

	switch(	m_pImage->GetDataType() )
	{
	case GIPL_CHAR:
	case GIPL_U_CHAR:
	case GIPL_SHORT:
	case GIPL_U_SHORT:
	case GIPL_U_INT:
	case GIPL_INT:
		//char uchar short zxhushort int uint
		for( int ix=0; ix<size[0]; ++ix )
		for( int iy=0; iy<size[1]; ++iy )
		for( int iz=0; iz<size[2]; ++iz )
		for( int it=0; it<size[3]; ++it )
		{
			float wto[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wto ) ;
			float wfrom[] = { wto[0],wto[1],wto[2],wto[3] } ;
			if( pTransform )
				pTransform->InverseTransformWorld( wto, wfrom );

			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckBkByWorld(background, wfrom[0], wfrom[1], wfrom[2], wfrom[3] ) +0.5);

			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}
		break;
	case GIPL_FLOAT:
	case GIPL_DOUBLE:
	default:
		for( int ix=0; ix<size[0]; ++ix )
		for( int iy=0; iy<size[1]; ++iy )
		for( int iz=0; iz<size[2]; ++iz )
		for( int it=0; it<size[3]; ++it )
		{
			float wto[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wto ) ;
			float wfrom[] = { wto[0],wto[1],wto[2],wto[3] } ;
			if( pTransform )
				pTransform->InverseTransformWorld( wto, wfrom );

			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckBkByWorld(background, wfrom[0], wfrom[1], wfrom[2], wfrom[3] ) );
			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}

		break;
	};


	return true ;
}

 
/// \return backward transform using transform list
/// \para 
/// bBackground = false, donot use background value but use closest pixel points to fill
template<class PixelType,class DerivativeType>
bool zxhImageModelingBaseT<PixelType, DerivativeType>::BackwardTransformImageByList(
								zxhImageDataT<PixelType>*pSave,
								const zxhTransformBase*const*pList,int nTransform, // const pList makes the function hard to use
								bool bBackground, PixelType background ) const
{
	if( pSave==0 || pSave->IsEmpty() )
		return false ;
	zxhImageDataT<PixelType> *pImg = this->m_pImage ; 
	const zxhImageModelingBaseT<PixelType,DerivativeType>* pModel = this ; 
	const zxhImageDataT<PixelType> * pTempImage = this->m_pImage ;  
	const int * size = pSave->GetImageSize() ;
	int pimgsize[4] = {1,1,1,1} ; 
	pImg->GetImageSize( pimgsize[0], pimgsize[1], pimgsize[2], pimgsize[3] ) ; 
	int iDimension = this->m_pImage->GetDimension() ; 
	zxhushort datatype = pSave->GetDataType();
	for(int it=0;it<size[3];++it)
	for(int iz=0;iz<size[2];++iz)
	for(int iy=0;iy<size[1];++iy)
	for(int ix=0;ix<size[0];++ix)
	{
		float from[ZXH_ImageDimensionMax]={ix,iy,iz,it};
		pSave->ImageToWorld(from);
		float to[ZXH_ImageDimensionMax]={from[0],from[1],from[2],from[3]};
		zxh::TransformWorld2WorldByCollection(pList,nTransform,pSave->GetDimension(), from,to);
		pImg->WorldToImage(to);
		
		if( pImg->InsideImageSafe( to[0],to[1],to[2],to[3] ) == false )
		{ 
			if( bBackground == true )
				pSave->SetPixelByGreyscale( ix,iy,iz,it, background ) ;
			else 
				pSave->SetPixelByGreyscale(ix,iy,iz,it, 
										pImg->GetPixelGreyscaleClosest(
										zxh::round(to[0]),zxh::round(to[1]),zxh::round(to[2]),zxh::round(to[3])) ) ;
			continue ;
		}

		switch( datatype )
		{
		case GIPL_CHAR:
		case GIPL_U_CHAR:
		case GIPL_SHORT:
		case GIPL_U_SHORT:
		case GIPL_U_INT:
		case GIPL_INT:
			pSave->SetPixelByGreyscale(ix,iy,iz,it, PixelType(zxh::round(pModel->GetPixelFloatValueWithCheck(to[0],to[1],to[2],to[3])))); 
			break;
		case GIPL_FLOAT:
		case GIPL_DOUBLE:
		case GIPL_C_FLOAT:
		case GIPL_C_DOUBLE:
		default:
			pSave->SetPixelByGreyscale(ix,iy,iz,it, PixelType( pModel->GetPixelFloatValueWithCheck(to[0],to[1],to[2],to[3]) ) );  
			break;
		}
	}  
	return true ;
}
/// fill background ZXH_Background
template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::BackwardTransformImage(
	zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const
{
	if( pSave==0 || pSave->IsEmpty() )
		return false ;
	int size[] = {1,1,1,1} ;
	pSave->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	PixelType value ;
	PixelType background = ZXH_Background ;
	PixelType vmin = this->m_pImage->GetPixelGreyscaleMin() ;
	if( background > vmin )
		background = vmin ;

	switch(	m_pImage->GetDataType() )
	{
	case GIPL_CHAR:
	case GIPL_U_CHAR:
	case GIPL_SHORT:
	case GIPL_U_SHORT:
	case GIPL_U_INT:
	case GIPL_INT:
		//char uchar short zxhushort int uint
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float wfrom[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wfrom ) ;
			float wto[] = { wfrom[0],wfrom[1],wfrom[2],wfrom[3] } ;
			if( pTransform )
				pTransform->TransformWorldToWorld(wfrom, wto );
			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckBkByWorld(background, wto[0], wto[1], wto[2], wto[3] )+0.5);

			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}
		break;
	case GIPL_FLOAT:
	case GIPL_DOUBLE:
	default:
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float wfrom[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wfrom ) ;
			float wto[] = { wfrom[0],wfrom[1],wfrom[2],wfrom[3] } ;
			if( pTransform )
				pTransform->TransformWorldToWorld(wfrom, wto );
			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckBkByWorld(background, wto[0], wto[1], wto[2], wto[3] ) );
			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}

		break;
	};

	return true ;
}

template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::BackwardTransformImageClosest(
	zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const
{
	if( pSave==0 || pSave->IsEmpty() )
		return false ;
	int size[] = {1,1,1,1} ;
	pSave->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	PixelType value ;

	switch(	m_pImage->GetDataType() )
	{
	case GIPL_CHAR:
	case GIPL_U_CHAR:
	case GIPL_SHORT:
	case GIPL_U_SHORT:
	case GIPL_U_INT:
	case GIPL_INT:
		//char uchar short zxhushort int uint
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float wfrom[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wfrom ) ;
			float wto[] = { wfrom[0],wfrom[1],wfrom[2],wfrom[3] } ;
			if( pTransform )
				pTransform->TransformWorldToWorld( wfrom, wto );

			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckClosestByWorld(wto[0], wto[1], wto[2], wto[3] )+0.5);

			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}
		break;
	case GIPL_FLOAT:
	case GIPL_DOUBLE:
	default:
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float wfrom[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wfrom ) ;
			float wto[] = { wfrom[0],wfrom[1],wfrom[2],wfrom[3] } ;
			if( pTransform )
				pTransform->TransformWorldToWorld( wfrom, wto );
			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckClosestByWorld(wto[0], wto[1], wto[2], wto[3] ) );
			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}

		break;
	};

	return true ;
}
template<class PixelType,class DerivativeType>
bool
zxhImageModelingBaseT<PixelType, DerivativeType>::BackwardTransformImageMirrorOutlier(
	zxhImageDataT<PixelType>* pSave, zxhTransformBase* pTransform )const
{
	if( pSave==0 || pSave->IsEmpty() )
		return false ;
	int size[] = {1,1,1,1} ;
	pSave->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	PixelType value ;

	switch(	m_pImage->GetDataType() )
	{
	case GIPL_CHAR:
	case GIPL_U_CHAR:
	case GIPL_SHORT:
	case GIPL_U_SHORT:
	case GIPL_U_INT:
	case GIPL_INT:
		//char uchar short zxhushort int uint
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float wfrom[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wfrom ) ;
			float wto[] = { wfrom[0],wfrom[1],wfrom[2],wfrom[3] } ;
			if( pTransform )
				pTransform->TransformWorldToWorld( wfrom, wto );

			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckMirrorByWorld(wto[0], wto[1], wto[2], wto[3] )+0.5);

			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}
		break;
	case GIPL_FLOAT:
	case GIPL_DOUBLE:
	default:
		for( int it=0; it<size[3]; ++it )
		for( int iz=0; iz<size[2]; ++iz )
		for( int iy=0; iy<size[1]; ++iy )
		for( int ix=0; ix<size[0]; ++ix )
		{
			float wfrom[] = { float(ix), float(iy), float(iz), float(it) } ;
			pSave->ImageToWorld( wfrom ) ;
			float wto[] = { wfrom[0],wfrom[1],wfrom[2],wfrom[3] } ;
			if( pTransform )
				pTransform->TransformWorldToWorld( wfrom, wto );
			value = static_cast<PixelType>( this->GetPixelFloatValueWithCheckMirrorByWorld(wto[0], wto[1], wto[2], wto[3] ) );
			pSave->SetPixelByGreyscale( ix, iy, iz, it, value ) ;
		}

		break;
	};

	return true ;
}
 

template<class PixelType,class DerivativeType>
void
zxhImageModelingBaseT<PixelType, DerivativeType>::GetPixelGradientByWorld(
					DerivativeType pwGrad[ZXH_ImageDimensionMax],float wx,float wy,float wz,float wt)const
{
	if( this->m_pImage == 0 )
	{
		std::cerr<<"error: have not set image before getting gradient\n" ;
		return ;
	}
	const float *spacing = this->m_pImage->GetImageSpacing() ;
	float pixel[]={wx,wy,wz,wt};
	for(int idim=0;idim<this->m_pImage->GetDimension();++idim)
	{
		if( spacing[idim] <= 0 )
		{
			std::cerr<<"error: get gradient of dimension "<<idim<<" which has "<<spacing[idim]<<" mm spacing\n" ;
			return ;
		}
		pixel[idim] += 0.5*spacing[idim];
		pwGrad[idim] = static_cast<DerivativeType>(this->GetPixelFloatValueWithCheckByWorld(pixel[0],pixel[1],pixel[2],pixel[3]));
		pixel[idim] -= spacing[idim] ;
		pwGrad[idim] -= static_cast<DerivativeType>(this->GetPixelFloatValueWithCheckByWorld(pixel[0],pixel[1],pixel[2],pixel[3]));
		pixel[idim] += 0.5*spacing[idim] ;
		pwGrad[idim] /= (spacing[idim]) ;
	}
};


/// partial differential on physical coordinate using image grid index
template<class PixelType,class DerivativeType>
DerivativeType
zxhImageModelingBaseT<PixelType,DerivativeType>::GetPixelPhysicalPartialDiff( zxh::ImageModeling_XYZ dim, float x, float y, float z, float t )const
{
	if( this->m_pImage == 0 )
	{
		std::cerr<<"error: have not set image before getting partial differential\n" ;
		return 0 ;
	}
	const float *spacing = this->m_pImage->GetImageSpacing() ;
	if( spacing[dim] <= 0 )
	{
		std::cerr<<"error: get differential of dimension "<<dim<<" which has "<<spacing[dim]<<" mm spacing\n" ;
		return 0 ;
	}
	float c[] = {x,y,z,t} ;

	DerivativeType g=0;
	c[dim] += 0.5 ;
	g = static_cast<DerivativeType>(this->GetPixelFloatValueWithCheck(c[0],c[1],c[2],c[3]));
	c[dim] -= 1.0 ;
	g -= static_cast<DerivativeType>(this->GetPixelFloatValueWithCheck(c[0],c[1],c[2],c[3]));
	c[dim] += 0.5 ;
	g /= (spacing[dim]) ;
	return g ;
};

typedef zxhImageModelingBaseT<ZXH_PixelTypeDefault,float> zxhImageModelingBase;
#define zxhImageModelingNearestT zxhImageModelingBaseT  //<PixelType,DerivativeType>
typedef zxhImageModelingNearestT<ZXH_PixelTypeDefault,float> zxhImageModelingNearest;
#endif


