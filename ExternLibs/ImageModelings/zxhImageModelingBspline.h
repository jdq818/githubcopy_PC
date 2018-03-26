/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhImageModelingBspline.h    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0 $

=========================================================================*/


#ifndef _ZXHIMAGEMODELINGBSPLINE_H
#define _ZXHIMAGEMODELINGBSPLINE_H

#include "zxhImageModelingBase.h"
#include <math.h>
#ifndef DBL_EPSILON
#define DBL_EPSILON ZXH_FloatPrecision
#endif
///
/// \class  zxhImageModelingBsplineT
/// \brief image modeling by bspline intepolation, 3d validated
///       call ComputeCoefficient() before use it
/// \ingroup zxhImageModelings
///


template<class PixelType=ZXH_PixelTypeDefault,class DerivativeType=float>
class zxhImageModelingBsplineT :
	public zxhImageModelingBaseT<PixelType,DerivativeType>
{
public:
	/// typedef
	typedef typename zxhImageModelingBaseT<PixelType,DerivativeType>::PixelFloatType PixelFloatType;
	/// constructor, default cubic bspline
	zxhImageModelingBsplineT(void){m_iBSplineDegree=3;};

	///\return
	~zxhImageModelingBsplineT(void){};

	///
	virtual zxhImageModelingBaseT<PixelType,DerivativeType>*	CloneTo(zxhImageModelingBaseT<PixelType,DerivativeType>* &pRet)const
	{
		if(pRet==0) 
			pRet=new zxhImageModelingBsplineT<PixelType,DerivativeType>;
		pRet->SetImage(this->m_pImage);
 		zxhImageDataT<float> * pco = &((zxhImageModelingBsplineT<PixelType,DerivativeType>*)pRet)->m_imgInterpolationCoefficient;
		this->m_imgInterpolationCoefficient.CloneTo( pco ) ;
		((zxhImageModelingBsplineT<PixelType,DerivativeType>*)pRet)->m_iBSplineDegree = m_iBSplineDegree ;
		((zxhImageModelingBsplineT<PixelType,DerivativeType>*)pRet)->m_ImageMin = m_ImageMin ;
		((zxhImageModelingBsplineT<PixelType,DerivativeType>*)pRet)->m_ImageMax = m_ImageMax ;
		return pRet;
	};

	/// using BS interpolation, first step is to compute coefficient
	bool ComputeCoefficient();

	/// inline interpolate without dimension check, without boundary check, assume []\in0~size-1
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy) const;
	/// inline interpolate without dimension check, without boundary check
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy,float fz) const;
	/// inline interpolate without dimension check, without boundary check
	virtual inline PixelFloatType	GetPixelFloatValueWithoutCheck(float fx,float fy,float fz,float ft) const;

	///  closest: GetPixelFloatValueWithCheck(float fx,float fy,float fz=0,float ft=0)
	virtual inline PixelFloatType	GetPixelFloatValueWithCheckClosest(float fx,float fy,float fz=0,float ft=0) const
	{
		if( this->m_pImage->GetDimension()==3 )
			return this->GetPixelValueWithinImage3D(fx,fy,fz) ;
		if( this->m_pImage->GetDimension()==2 )
			return this->GetPixelValueWithinImage2D(fx,fy) ; 

		std::cerr<<"error: zxhtodo, 4D bspline interpolation has not implemented yet \n" ;
		return 0; // zxhtodo
	};

	/// 	virtual void GetPixelPhysicalGradientOnImageGridCoordinate(DerivativeType pGrad[ZXH_ImageDimensionMax],float x,float y=0,float z=0,float t=0)const;

private:
	/// assume fxfyfz be within image size
	PixelFloatType GetPixelValueWithinImage3D(float fx,float fy, float fz) const;
	/// assume fxfy be within image size
	PixelFloatType GetPixelValueWithinImage2D(float fx,float fy) const;
	/// Initialize anti-causal coefficients
	float InitialAntiCausalCoefficient(float c[], int DataLength, float z) const
	{
	  /* this initialization corresponds to mirror boundaries */
	  return((z / (z * z - 1.0)) * (z * c[DataLength - 2] + c[DataLength - 1]));
	}
	/// Initialize causal coefficients
	float InitialCausalCoefficient(float c[], int DataLength, float z, float Tolerance) const
	{
	  float Sum, zn, z2n, iz;
	  int n, Horizon;
	  /* this initialization corresponds to mirror boundaries */
	  Horizon = DataLength;
	  if (Tolerance > 0.0)
	    Horizon = (int)ceil(log(Tolerance) / log(fabs(z)));
	  if (Horizon < DataLength)
	  {
	    /* accelerated loop */
	    zn = z; Sum = c[0];
	    for (n = 1; n < Horizon; n++)
	    {
	      Sum += zn * c[n]; zn *= z;
	    }
	    return(Sum);
	  } else
	  {
	    /* full loop */
	    zn = z; iz = 1.0 / z;
	    z2n = pow(z, (float)(DataLength - 1));
	    Sum = c[0] + z2n * c[DataLength - 1];
	    z2n *= z2n * iz;
	    for (n = 1; n <= DataLength - 2; n++)
	    {
	      Sum += (zn + z2n) * c[n];
	      zn *= z; z2n *= iz;
	    }
	    return(Sum / (1.0 - zn * zn));
	  }
	};

	/// Convert voxel values to B-spline coefficients
	void ConvertToInterpolationCoefficients(float *c, int DataLength, float *z, int NbPoles, float Tolerance) const
	{
	  float Lambda = 1.0; int n, k;
	  /* special case required by mirror boundaries */
	  if (DataLength == 1)
	    return;
	  /* compute the overall gain */
	  for (k = 0; k < NbPoles; k++)
	    Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
	  /* apply the gain */
	  for (n = 0; n < DataLength; n++)
	    c[n] *= Lambda;
	  /* loop over all poles */
	  for (k = 0; k < NbPoles; k++)
	  {
	    /* causal initialization */
	    c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
	    /* causal recursion */
	    for (n = 1; n < DataLength; n++)
	      c[n] += z[k] * c[n - 1];
	    /* anticausal initialization */
	    c[DataLength - 1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
	    /* anticausal recursion */
	    for (n = DataLength - 2; 0 <= n; n--)
	      c[n] = z[k] * (c[n + 1] - c[n]);
	  }
	};
	///
	zxhImageDataT<float> m_imgInterpolationCoefficient;
	///
	int m_iBSplineDegree ;
	///
	PixelType m_ImageMin ;
	///
	PixelType m_ImageMax ;

};

template<class PixelType,class DerivativeType>
bool zxhImageModelingBsplineT<PixelType,DerivativeType>::ComputeCoefficient()
{
	if( this->m_pImage==0 || this->m_pImage->IsEmpty() )
		return false ;
	int imagesize[] = {1,1,1,1}, imagedimension = this->m_pImage->GetDimension() ;
	float imagespacing[] = {1,1,1,1} ;
	this->m_pImage->GetImageSize( imagesize[0], imagesize[1], imagesize[2], imagesize[3] ) ;
	this->m_pImage->GetImageSpacing( imagespacing[0], imagespacing[1], imagespacing[2], imagespacing[3] ) ;

	if( imagedimension==3 && imagesize[2]>1 )
	{ //  500x500x350 about 8 seconds
	  this->m_ImageMin = this->m_pImage->GetPixelGreyscaleMin() ;
	  this->m_ImageMax = this->m_pImage->GetPixelGreyscaleMax() ;
	  zxhImageInfo imageinfo ;
	  this->m_pImage->GetImageInfo( &imageinfo ) ;

	  this->m_imgInterpolationCoefficient.NewImage(imagedimension, imagesize, imagespacing, &imageinfo);
	  float *data; float Pole[2];
	  int NbPoles; int x, y, z;

	  // recover the poles from a lookup table
	  switch (m_iBSplineDegree)
	  {
	  case 2:
		NbPoles = 1;
		Pole[0] = sqrt(8.0) - 3.0;
		break;
	  case 3:
		NbPoles = 1;
		Pole[0] = sqrt(3.0) - 2.0;
		break;
	  case 4:
		NbPoles = 2;
		Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
		Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
		break;
	  case 5:
		NbPoles = 2;
		Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
		  - 13.0 / 2.0;
		Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
		  - 13.0 / 2.0;
		break;
	  default:
		std::cerr << "error: Invalid spline degree\n" ;
		exit(1);
	  }

	  // convert the image samples into interpolation coefficients

	  int _x=imagesize[0], _y=imagesize[1], _z=imagesize[2] ;
	  // in-place separable process, along x
	  data = new float[_x];
	  for (z = 0; z < imagesize[2]; z++)
	  {
		for (y = 0; y < imagesize[1]; y++)
		{
		  for (x = 0; x < imagesize[0]; x++)
		  {
			data[x] = this->m_pImage->GetPixelGreyscale(x, y, z);
		  }
		  this->ConvertToInterpolationCoefficients(data, imagesize[0], Pole, NbPoles, DBL_EPSILON);
		  for (x = 0; x < imagesize[0]; x++)
		  {
			this->m_imgInterpolationCoefficient.SetPixelByGreyscale(x, y, z, 0, data[x]);
		  }
		}
	  }
	  delete []data;

	  //in-place separable process, along y
	  data = new float[_y];
	  for (z = 0; z < _z; z++)
	  {
		for (x = 0; x < _x; x++)
		{
		  for (y = 0; y < _y; y++)
		  {
			data[y] = this->m_imgInterpolationCoefficient.GetPixelGreyscale(x, y, z);
		  }
		  ConvertToInterpolationCoefficients(data, _y, Pole, NbPoles, DBL_EPSILON);
		  for (y = 0; y < _y; y++)
		  {
			this->m_imgInterpolationCoefficient.SetPixelByGreyscale(x, y, z, 0, data[y]);
		  }
		}
	  }
	  delete []data;

	  // in-place separable process, along z
	  data = new float[_z];
	  for (y = 0; y < _y; y++)
	  {
		for (x = 0; x < _x; x++)
		{
		  for (z = 0; z < _z; z++)
		  {
			data[z] = this->m_imgInterpolationCoefficient.GetPixelGreyscale(x, y, z);
		  }
		  ConvertToInterpolationCoefficients(data, _z, Pole, NbPoles, DBL_EPSILON);
		  for (z = 0; z < _z; z++)
		  {
			 this->m_imgInterpolationCoefficient.SetPixelByGreyscale(x, y, z, 0, data[z]);
		  }
		}
	  }
	  delete []data;
      return true ;
	}
	if( imagedimension==2 || imagesize[2]==1 )
	{
	  this->m_ImageMin = this->m_pImage->GetPixelGreyscaleMin() ;
	  this->m_ImageMax = this->m_pImage->GetPixelGreyscaleMax() ;
	  zxhImageInfo imageinfo ;
	  this->m_pImage->GetImageInfo( &imageinfo ) ;

	  this->m_imgInterpolationCoefficient.NewImage(imagedimension, imagesize, imagespacing, &imageinfo);

	  float *data, Pole[2];
	  int NbPoles, x, y;

	  /* recover the poles from a lookup table */
	  switch (m_iBSplineDegree) {
	  case 2:
	    NbPoles = 1;
	    Pole[0] = sqrt(8.0) - 3.0;
	    break;
	  case 3:
	    NbPoles = 1;
	    Pole[0] = sqrt(3.0) - 2.0;
	    break;
	  case 4:
	    NbPoles = 2;
	    Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
	    Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
	    break;
	  case 5:
	    NbPoles = 2;
	    Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
	      - 13.0 / 2.0;
	    Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
	      - 13.0 / 2.0;
	    break;
	  default:
	    std::cerr << "error: Invalid spline degree\n";
	    exit(1);
	  }

	  /* convert the image samples into interpolation coefficients */
	  int _x=imagesize[0], _y=imagesize[1];

	  /* in-place separable process, along x */
	  data = new float[_x];
	  for (y = 0; y < _y; y++)
	  {
	    for (x = 0; x < _x; x++)
	      data[x] = this->m_pImage->GetPixelGreyscale(x, y, 0);
	    ConvertToInterpolationCoefficients(data, _x, Pole, NbPoles, DBL_EPSILON);
	    for (x = 0; x < _x; x++)
	    	this->m_imgInterpolationCoefficient.SetPixelByGreyscale(x, y, 0,0, data[x]);
	  }
	  delete []data;

	  /* in-place separable process, along y */
	  data = new float[_y];
	  for (x = 0; x < _x; x++)
	  {
	    for (y = 0; y < _y; y++)
	      data[y] = this->m_imgInterpolationCoefficient.GetPixelGreyscale(x, y, 0);
	    ConvertToInterpolationCoefficients(data, _y, Pole, NbPoles, DBL_EPSILON);
	    for (y = 0; y < _y; y++)
	    	this->m_imgInterpolationCoefficient.SetPixelByGreyscale(x, y, 0,0, data[y]);
	  }
	  delete []data;
	  return true ;
	}

	return false ;
}
template<class PixelType,class DerivativeType>
inline typename zxhImageModelingBsplineT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingBsplineT<PixelType,DerivativeType>::GetPixelValueWithinImage2D(float fx,float fy) const
{
	if( this->m_imgInterpolationCoefficient.IsEmpty() )
	{
		std::cerr<<"error: wrong usage of bspline interpolation modeling, should call compute coefficient first\n";
		exit(-1);
		return 0;
	}
	if( fx<0 ) fx = 0 ; 
	if( fy<0 ) fy = 0 ; 
	if( fx>this->m_pImage->GetImageSize()[0]-1 ) fx = this->m_pImage->GetImageSize()[0]-1 ; 
	if( fy>this->m_pImage->GetImageSize()[1]-1 ) fy = this->m_pImage->GetImageSize()[1]-1 ; 

	  int i, j, m, xIndex[6], yIndex[6];
	  float xWeight[6], yWeight[6], value, w, w2, w4, t, t0, t1;

	  /* compute the interpolation indexes */
	  if (m_iBSplineDegree & 1)
	  {
	    i = (int)floor(fx) - this->m_iBSplineDegree / 2;
	    j = (int)floor(fy) - this->m_iBSplineDegree / 2;
	    for (m = 0; m <= this->m_iBSplineDegree; m++)
	    {
	      xIndex[m] = i++;
	      yIndex[m] = j++;
	    }
	  } else
	  {
	    i = (int)floor(fx + 0.5) - this->m_iBSplineDegree / 2;
	    j = (int)floor(fy + 0.5) - this->m_iBSplineDegree / 2;
	    for (m = 0; m <= this->m_iBSplineDegree; m++)
	    {
	      xIndex[m] = i++;
	      yIndex[m] = j++;
	    }
	  }

	  /* compute the interpolation weights */
	  switch (m_iBSplineDegree)
	  {
	  case 2:
	    /* x */
	    w = fx - (float)xIndex[1];
	    xWeight[1] = 3.0 / 4.0 - w * w;
	    xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
	    xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
	    /* y */
	    w = fy - (float)yIndex[1];
	    yWeight[1] = 3.0 / 4.0 - w * w;
	    yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
	    yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
	    break;
	  case 3:
	    /* x */
	    w = fx - (float)xIndex[1];
	    xWeight[3] = (1.0 / 6.0) * w * w * w;
	    xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
	    xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
	    xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
	    /* y */
	    w = fy - (float)yIndex[1];
	    yWeight[3] = (1.0 / 6.0) * w * w * w;
	    yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
	    yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
	    yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
	    break;
	  case 4:
	    /* x */
	    w = fx - (float)xIndex[2];
	    w2 = w * w;
	    t = (1.0 / 6.0) * w2;
	    xWeight[0] = 1.0 / 2.0 - w;
	    xWeight[0] *= xWeight[0];
	    xWeight[0] *= (1.0 / 24.0) * xWeight[0];
	    t0 = w * (t - 11.0 / 24.0);
	    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
	    xWeight[1] = t1 + t0;
	    xWeight[3] = t1 - t0;
	    xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
	    xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
	    /* y */
	    w = fy - (float)yIndex[2];
	    w2 = w * w;
	    t = (1.0 / 6.0) * w2;
	    yWeight[0] = 1.0 / 2.0 - w;
	    yWeight[0] *= yWeight[0];
	    yWeight[0] *= (1.0 / 24.0) * yWeight[0];
	    t0 = w * (t - 11.0 / 24.0);
	    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
	    yWeight[1] = t1 + t0;
	    yWeight[3] = t1 - t0;
	    yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
	    yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
	    break;
	  case 5:
	    /* x */
	    w = fx - (float)xIndex[2];
	    w2 = w * w;
	    xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
	    w2 -= w;
	    w4 = w2 * w2;
	    w -= 1.0 / 2.0;
	    t = w2 * (w2 - 3.0);
	    xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
	    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
	    t1 = (-1.0 / 12.0) * w * (t + 4.0);
	    xWeight[2] = t0 + t1;
	    xWeight[3] = t0 - t1;
	    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
	    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
	    xWeight[1] = t0 + t1;
	    xWeight[4] = t0 - t1;
	    /* y */
	    w = fy - (float)yIndex[2];
	    w2 = w * w;
	    yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
	    w2 -= w;
	    w4 = w2 * w2;
	    w -= 1.0 / 2.0;
	    t = w2 * (w2 - 3.0);
	    yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
	    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
	    t1 = (-1.0 / 12.0) * w * (t + 4.0);
	    yWeight[2] = t0 + t1;
	    yWeight[3] = t0 - t1;
	    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
	    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
	    yWeight[1] = t0 + t1;
	    yWeight[4] = t0 - t1;
	    break;
	  default:
	    std::cerr<<"error: Invalid spline degree\n";
	    return(0.0);
	  }

	  /* perform interpolation */
	  value = 0.0;
	  for (j = 0; j <= this->m_iBSplineDegree; j++)
	  {
	    for (i = 0; i <= this->m_iBSplineDegree; i++)
	      value += xWeight[i] * yWeight[j] * this->m_imgInterpolationCoefficient.GetPixelGreyscale(xIndex[i], yIndex[j]);
	  }

	  //if (this->_clamped)
	  if (value > this->m_ImageMax) return PixelFloatType(this->m_ImageMax);
	  if (value < this->m_ImageMin) return PixelFloatType(this->m_ImageMin);

	  return(value);

}

template<class PixelType,class DerivativeType>
		inline typename zxhImageModelingBsplineT<PixelType,DerivativeType>::PixelFloatType
		zxhImageModelingBsplineT<PixelType,DerivativeType>::GetPixelValueWithinImage3D(float fx,float fy,float fz) const
{ 
	if( this->m_imgInterpolationCoefficient.IsEmpty() )
	{
		std::cerr<<"error: wrong usage of bspline interpolation modeling, should call compute coefficient first\n";
		exit(-1);
		return 0;
	}
	if( fx<0 ) fx = 0 ; 
	if( fy<0 ) fy = 0 ; 
	if( fz<0 ) fz = 0 ; 
	if( fx>this->m_pImage->GetImageSize()[0]-1 ) fx = this->m_pImage->GetImageSize()[0]-1 ; 
	if( fy>this->m_pImage->GetImageSize()[1]-1 ) fy = this->m_pImage->GetImageSize()[1]-1 ; 
	if( fz>this->m_pImage->GetImageSize()[2]-1 ) fz = this->m_pImage->GetImageSize()[2]-1 ; 

  int i, j, k, m;
  int xIndex[6], yIndex[6], zIndex[6];
  PixelFloatType xWeight[6], yWeight[6], zWeight[6];
  PixelFloatType value, w, w2, w4, t, t0, t1;

  // compute the interpolation indexes
  if (m_iBSplineDegree & 1)
  {
    i = (int)floor(fx) - this->m_iBSplineDegree / 2;
    j = (int)floor(fy) - this->m_iBSplineDegree / 2;
    k = (int)floor(fz) - this->m_iBSplineDegree / 2;
    for (m = 0; m <= this->m_iBSplineDegree; m++)
    {
      xIndex[m] = i++;
      yIndex[m] = j++;
      zIndex[m] = k++;
    }
  } else
  {
    i = (int)floor(fx + 0.5) - this->m_iBSplineDegree / 2;
    j = (int)floor(fy + 0.5) - this->m_iBSplineDegree / 2;
    k = (int)floor(fz + 0.5) - this->m_iBSplineDegree / 2;
    for (m = 0; m <= this->m_iBSplineDegree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
      zIndex[m] = k++;
    }
  }

  // compute the interpolation weights
  switch (m_iBSplineDegree)
  {
  case 2:
    // fx
    w = fx - (float)xIndex[1];
    xWeight[1] = 3.0 / 4.0 - w * w;
    xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
    xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
    // fy
    w = fy - (float)yIndex[1];
    yWeight[1] = 3.0 / 4.0 - w * w;
    yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
    yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
    // fz
    w = fz - (float)zIndex[1];
    zWeight[1] = 3.0 / 4.0 - w * w;
    zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
    zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
    break;
  case 3:
    // fx
    w = fx - (float)xIndex[1];
    xWeight[3] = (1.0 / 6.0) * w * w * w;
    xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
    xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
    xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
    // fy
    w = fy - (float)yIndex[1];
    yWeight[3] = (1.0 / 6.0) * w * w * w;
    yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
    yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
    yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
    // fz
    w = fz - (float)zIndex[1];
    zWeight[3] = (1.0 / 6.0) * w * w * w;
    zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3];
    zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
    zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
    break;
  case 4:
    // fx
    w = fx - (float)xIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    xWeight[0] = 1.0 / 2.0 - w;
    xWeight[0] *= xWeight[0];
    xWeight[0] *= (1.0 / 24.0) * xWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    xWeight[1] = t1 + t0;
    xWeight[3] = t1 - t0;
    xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
    xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
    // fy
    w = fy - (float)yIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    yWeight[0] = 1.0 / 2.0 - w;
    yWeight[0] *= yWeight[0];
    yWeight[0] *= (1.0 / 24.0) * yWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    yWeight[1] = t1 + t0;
    yWeight[3] = t1 - t0;
    yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
    yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
    // fz
    w = fz - (float)zIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    zWeight[0] = 1.0 / 2.0 - w;
    zWeight[0] *= zWeight[0];
    zWeight[0] *= (1.0 / 24.0) * zWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    zWeight[1] = t1 + t0;
    zWeight[3] = t1 - t0;
    zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
    zWeight[2] = 1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4];
    break;
  case 5:
    // fx
    w = fx - (float)xIndex[2];
    w2 = w * w;
    xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    xWeight[2] = t0 + t1;
    xWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    xWeight[1] = t0 + t1;
    xWeight[4] = t0 - t1;
    //  fy
    w = fy - (float)yIndex[2];
    w2 = w * w;
    yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    yWeight[2] = t0 + t1;
    yWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    yWeight[1] = t0 + t1;
    yWeight[4] = t0 - t1;
    // fz
    w = fz - (float)zIndex[2];
    w2 = w * w;
    zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    zWeight[2] = t0 + t1;
    zWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    zWeight[1] = t0 + t1;
    zWeight[4] = t0 - t1;
    break;
  default:
    printf("Invalid spline degree\n");
    return(0.0);
  }

  // perform interpolation
  value = 0.0;
  for (k = 0; k <= this->m_iBSplineDegree; k++)
  {
    for (j = 0; j <= this->m_iBSplineDegree; j++)
    {
      for (i = 0; i <= this->m_iBSplineDegree; i++)
      {

    	PixelType pixelvalue = 0 ;
    	if( this->m_imgInterpolationCoefficient.InsideImage( xIndex[i], yIndex[j], zIndex[k], 0) )
    	  pixelvalue = this->m_imgInterpolationCoefficient.GetPixelGreyscale(xIndex[i], yIndex[j], zIndex[k]);
    	else
    	  pixelvalue = this->m_imgInterpolationCoefficient.GetPixelGreyscaleClosest(xIndex[i], yIndex[j], zIndex[k],0) ;
    	value += xWeight[i] * yWeight[j] * zWeight[k] * pixelvalue ;
      }
    }
  }

  //if (this->_clamped)
  if (value > this->m_ImageMax) return PixelFloatType(this->m_ImageMax);
  if (value < this->m_ImageMin) return PixelFloatType(this->m_ImageMin);

  return(value);
};
/// inline
template<class PixelType,class DerivativeType>
inline typename zxhImageModelingBsplineT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingBsplineT<PixelType,DerivativeType>::GetPixelFloatValueWithoutCheck(float fx,float fy)const
{
	return this->GetPixelValueWithinImage2D(fx,fy) ;
	return 0; //zxhtodo 
};

template<class PixelType,class DerivativeType>
inline typename zxhImageModelingBsplineT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingBsplineT<PixelType,DerivativeType>::GetPixelFloatValueWithoutCheck(float fx,float fy,float fz)const
{
	if( this->m_pImage->GetDimension()==2 )
		return this->GetPixelValueWithinImage2D(fx,fy) ;
	return this->GetPixelValueWithinImage3D(fx,fy,fz) ;
};

template<class PixelType,class DerivativeType>
inline typename zxhImageModelingBsplineT<PixelType,DerivativeType>::PixelFloatType
zxhImageModelingBsplineT<PixelType,DerivativeType>::GetPixelFloatValueWithoutCheck(float fx,float fy,float fz,float ft)const
{
	if( this->m_pImage->GetDimension()==3 )
		return this->GetPixelValueWithinImage3D(fx,fy,fz) ;
	if( this->m_pImage->GetDimension()==2 )
		return this->GetPixelValueWithinImage2D(fx,fy) ; 

	std::cerr<<"error: zxhtodo, 4D bspline interpolation has not implemented yet \n" ;
	return 0; // zxhtodo
};
  
typedef zxhImageModelingBsplineT<ZXH_PixelTypeDefault,float> zxhImageModelingBspline;

#endif

