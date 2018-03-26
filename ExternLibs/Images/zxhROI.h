/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 ZHUANG, Xia Hai
  Date:      From 2010-11
  Version:	 V 2.1

=========================================================================*/

#ifndef zxhROI_h
#define zxhROI_h
#include "zxh.h"

///
/// \class zxhROI
/// \brief:  ROI will be from <= to (automatically corrected in the left-hand coordinate when setting)
/// \ingroup zxhImageDataT
/// 
class zxhROI
{
public:
	///
	zxhROI();

	///
	virtual ~zxhROI();

	/// new an image and clone
	virtual zxhROI * CloneTo(zxhROI * & pRet) const;
	
	///     
	virtual void SetRoiFromTo( const float from[], const float to[], int dim ) ; 
	///
	virtual void SetRoiFromTo( const int from[], const int to[], int dim ) ; 
	///
	virtual void SetRoi3DFromTo( int xfrom, int yfrom, int zfrom, int xto, int yto, int zto ) ; 
	///
	virtual void SetRoi3DFromTo( float xfrom, float yfrom, float zfrom, float xto, float yto, float zto ) ; 
	///
	virtual void SetRoiCenterRadius( const float center[], float radius, int dim ) ; 
	///
	virtual void SetRoiCenterSideLength( const float center[], const float length[], int dim ) ; 

	/// \return whether the ROI has been set
	virtual bool GetRoiFromTo( float from[], float to[] ) const ; 
	///
	virtual bool GetRoiFromTo( int from[], int to[] ) const ; 
	///
	virtual bool GetRoi3DFromTo( int& xfrom, int& yfrom, int& zfrom, int& xto, int& yto, int& zto  ) const ; 
	///
	virtual bool GetRoi3DFromTo( float& xfrom, float& yfrom, float& zfrom, float& xto, float& yto, float& zto  ) const ; 

protected:
	///
	int m_iDimension ; 
	///
	int m_aiFrom[ZXH_ImageDimensionMax] ;
	///
	int m_aiTo[ZXH_ImageDimensionMax] ;
	///
	float m_afFrom[ZXH_ImageDimensionMax] ;
	///
	float m_afTo[ZXH_ImageDimensionMax] ;
	/// 

	/// \return whether need to perform the correction
	/// correct ROI to be from <= to
	bool CorrectFromTo() ; 
} ;

 
#endif //zxhROI_h

