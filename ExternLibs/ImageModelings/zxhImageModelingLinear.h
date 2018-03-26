
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhImageModelingLinear.h    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0, 2.1 $

=========================================================================*/


#ifndef _ZXHIMAGEMODELINGLinear_H_20051013
#define _ZXHIMAGEMODELINGLinear_H_20051013

#include "zxhImageModelingBase.h"
#include <math.h>

///
/// \class  zxhImageModelingLinearT
/// \brief image modeling by Linear interpolation
///      I(x) = (1-w)*(1-v)*A00 + w*v*A11 + (1-w)*v*A01 + w*(1-v)*A10
/// \ingroup zxhImageModelings
///


template<class PixelType=ZXH_PixelTypeDefault,class DerivativeType=float>
class zxhImageModelingLinearT :
	public zxhImageModelingBaseT<PixelType,DerivativeType>
{
public:
	/// typedef
	typedef typename zxhImageModelingBaseT<PixelType,DerivativeType>::PixelFloatType PixelFloatType;
	/// constructor
	zxhImageModelingLinearT(void){};

	///\return
	~zxhImageModelingLinearT(void){};

	///
	virtual zxhImageModelingBaseT<PixelType,DerivativeType>*	CloneTo(zxhImageModelingBaseT<PixelType,DerivativeType>* &pRet)const
	{
		if(pRet==0) 
			pRet=new zxhImageModelingLinearT<PixelType,DerivativeType>;
		pRet->SetImage(this->m_pImage);
		return pRet;
	};

	/// inline interpolate without dimension check, without boundary check, assume []\in [0~size-1) if int(fx) >= size-1, then int(fx)=size-2, u=1 
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy) const;
	/// inline interpolate without dimension check, without boundary check
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy,float fz) const;
	/// inline interpolate without dimension check, without boundary check
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy,float fz,float ft) const;

	///  closest: GetPixelFloatValueWithCheck(float fx,float fy,float fz=0,float ft=0)
	virtual inline PixelFloatType	GetPixelFloatValueWithCheckClosest(float fx,float fy,float fz=0,float ft=0) const;

};
/// inline
template<class PixelType,class DerivativeType>
inline typename zxhImageModelingLinearT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingLinearT<PixelType,DerivativeType>::GetPixelFloatValueWithoutCheck(float fx,float fy)const
{ 
	int nx,ny,nz,nt ;
	this->m_pImage->GetImageSize( nx,ny,nz,nt ) ;
	if( fx<0 ) fx = 0 ; 
	if( fy<0 ) fy = 0 ; 
	int ix	= static_cast<int>(fx),
		iy	= static_cast<int>(fy);

	float timx	= fx-static_cast<float>(ix);
	float timy	= fy-static_cast<float>(iy);
 
	// ----------- 2013-1 to place GetPixelGreyscaleClose as GetPixelGreyscale 
	int ixa1 = ix+1, iya1 = iy+1 ; 
	if( ixa1>=nx ) ixa1=ix; 
	if( iya1>=ny ) iya1=iy; 

	float op1, op4 ;

	op1	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy )+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy );
	op4	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1 )+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1 );
	return ((1.0f-timy)*op1+timy*op4) ;
};

template<class PixelType,class DerivativeType>
inline typename zxhImageModelingLinearT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingLinearT<PixelType,DerivativeType>::GetPixelFloatValueWithoutCheck(float fx,float fy,float fz)const
{
	if( this->m_pImage->GetDimension()==2 )
		return this->GetPixelFloatValueWithoutCheck(fx,fy); 
	
	if( fx<0 ) fx = 0 ; 
	if( fy<0 ) fy = 0 ; 
	if( fz<0 ) fz = 0 ; 
	int ix	= static_cast<int>(fx),
			iy	= static_cast<int>(fy),
			iz	= static_cast<int>(fz) ;

	float timx	= fx-static_cast<float>(ix);
	float timy	= fy-static_cast<float>(iy);
	float timz	= fz-static_cast<float>(iz);
 
//	if( fx== float(size[0]-1) )   ------------- NOT needed anymore because we use GetPixelGreyscale below
//	{ timx = 1 ; ix = (size[0]-2); }
//	if( fy== float(size[1]-1) )
//	{ timy = 1 ; iy	= (size[1]-2); }
//	if( fz== float(size[2]-1) )
//	{ timz = 1 ; iz	= (size[2]-2); }
	
	// ----------- 2013-1 to place GetPixelGreyscaleClose as GetPixelGreyscale 
	int nx,ny,nz,nt ;
	this->m_pImage->GetImageSize( nx,ny,nz,nt ) ;
	int ixa1 = ix+1, iya1 = iy+1, iza1 = iz+1 ; 
	if( ixa1>=nx ) ixa1=ix; 
	if( iya1>=ny ) iya1=iy;
	if( iza1>=nz ) iza1=iz; 


	float op1,op2,op3,op4,op11,op12 ;

	op1	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iz)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iz);
	op4	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iz)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iz);
	op11= (1.0f-timy)*op1+timy*op4;

	op2	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iza1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iza1);
	op3	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iza1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iza1);
	op12= (1.0f-timy)*op2+timy*op3;

	return ((1.0f-timz)*op11+timz*op12) ;
};

template<class PixelType,class DerivativeType>
inline typename zxhImageModelingLinearT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingLinearT<PixelType,DerivativeType>::GetPixelFloatValueWithoutCheck(float fx,float fy,float fz,float ft)const
{
	if( this->m_pImage->GetDimension()==3 )
		return this->GetPixelFloatValueWithoutCheck(fx,fy,fz);
	if( this->m_pImage->GetDimension()==2 )
		return this->GetPixelFloatValueWithoutCheck(fx,fy);
	
	if( fx<0 ) fx = 0 ; 
	if( fy<0 ) fy = 0 ; 
	if( fz<0 ) fz = 0 ; 
	if( ft<0 ) ft = 0 ; 
	int ix	= static_cast<int>(fx), 
			iy	= static_cast<int>(fy),
			iz	= static_cast<int>(fz),
			it	= static_cast<int>(ft);

	float timx	= fx-static_cast<float>(ix);
	float timy	= fy-static_cast<float>(iy);
	float timz	= fz-static_cast<float>(iz);
	float timt	= ft-static_cast<float>(it);
  
	 
	int nx,ny,nz,nt ;
	this->m_pImage->GetImageSize( nx,ny,nz,nt ) ;
	int ixa1 = ix+1, iya1 = iy+1, iza1 = iz+1, ita1 = it+1 ; 
	if( ixa1>=nx ) ixa1=ix; 
	if( iya1>=ny ) iya1=iy;
	if( iza1>=nz ) iza1=iz;
	if( ita1>=nt ) ita1=it;

	float op1,op2,op3,op4,op11,op12,op1t,op2t;

	op1	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iz,it)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iz,it);
	op4	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iz,it)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iz,it);
	op11= (1.0f-timy)*op1+timy*op4;

	op2	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iza1,it)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iza1,it);
	op3	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iza1,it)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iza1,it);
	op12= (1.0f-timy)*op2+timy*op3;
	op1t= (1.0f-timz)*op11+timz*op12;

	op1	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iz,ita1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iz,ita1);
	op2	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iza1,ita1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iza1,ita1);
	op3	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iza1,ita1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iza1,ita1);
	op4	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iz,ita1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iz,ita1);
	op11= (1.0f-timy)*op1+timy*op4;
	op12= (1.0f-timy)*op2+timy*op3;
	op2t= (1.0f-timz)*op11+timz*op12;

	return ((1.0f-timt)*op1t+timt*op2t);
};


template<class PixelType,class DerivativeType>
inline typename zxhImageModelingLinearT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingLinearT<PixelType,DerivativeType>::GetPixelFloatValueWithCheckClosest(float fx,float fy,float fz,float ft)const
{
	m_pImage->GetClosestPixelCoordinateWithinImage( fx,fy,fz,ft ) ;
	switch(this->m_pImage->GetDimension())
	{
	case 3: 
	{
		int nx,ny,nz,nt ;
		this->m_pImage->GetImageSize( nx,ny,nz,nt ) ; 
		int ix	= static_cast<int>(fx),
		iy	= static_cast<int>(fy),
		iz	= static_cast<int>(fz) ;
		int ixa1 = ix+1, iya1 = iy+1, iza1 = iz+1 ;  
		if( ixa1>=nx ) ixa1=ix; 
		if( iya1>=ny ) iya1=iy;
		if( iza1>=nz ) iza1=iz; 

		float timx	= fx-static_cast<float>(ix);
		float timy	= fy-static_cast<float>(iy);
		float timz	= fz-static_cast<float>(iz);  
		float op1,op2,op3,op4,op11,op12 ;

		op1	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iz)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iz);
		op4	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iz)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iz);
		op11= (1.0f-timy)*op1+timy*op4;

		op2	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iy,iza1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iy,iza1);
		op3	= (1.0f-timx)*(float)this->m_pImage->GetPixelGreyscale(ix,iya1,iza1)+(timx)*(float)this->m_pImage->GetPixelGreyscale(ixa1,iya1,iza1);
		op12= (1.0f-timy)*op2+timy*op3;

		return ((1.0f-timz)*op11+timz*op12) ;
	}
	break;
	case 2: return this->GetPixelFloatValueWithoutCheck(fx,fy); break;
	case 4: return this->GetPixelFloatValueWithoutCheck(fx,fy,fz,ft);break;
	}
	return this->GetPixelFloatValueWithoutCheck(fx,fy,fz,ft);
};


typedef zxhImageModelingLinearT<ZXH_PixelTypeDefault,float> zxhImageModelingLinear;
#endif

