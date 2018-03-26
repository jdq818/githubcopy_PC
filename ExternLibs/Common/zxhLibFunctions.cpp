
 /*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhLibFunctions.cpp    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0 $

=========================================================================*/

#include "zxhLibFunctions.h"

#ifdef HAS_VTK
#include "vtkImageData.h"
#include "vtkImageThreshold.h"
#include "vtkImageEuclideanDistance.h"
#include "vtkImageProgressIterator.h"
#endif

namespace zxh
{
bool	vtkEuclideanDistanceTransform(
					  zxhImageData		* pImage,
					  zxhImageDataT<float>* pDistance,
					  int	iForegroundFrom, int	iForegroundTo, // inclusive
					  bool	OnForeground ,
					  bool	ReturnSquare )
{
#ifdef HAS_VTK
	vtkImageData *pvtkImage = vtkImageData::New() ;
	ConvertToVtkData( pImage, pvtkImage, GIPL_SHORT, (short)1 ) ;
	vtkImageThreshold *Threshold = vtkImageThreshold::New();
	Threshold->SetInput( pvtkImage );
	Threshold->ThresholdBetween( iForegroundFrom, iForegroundTo );
	Threshold->SetOutputScalarTypeToShort();
	if(OnForeground)
	{
		Threshold->SetInValue( 1 );
		Threshold->SetOutValue( 0 );
	}
	else
	{
		Threshold->SetInValue( 0 );
		Threshold->SetOutValue( 1 );
	}
	Threshold->Update();

	vtkImageEuclideanDistance * DistanceTransform = vtkImageEuclideanDistance::New();
	DistanceTransform->SetInput(Threshold->GetOutput());
	DistanceTransform->ConsiderAnisotropyOn();
	DistanceTransform->Update();
	DistanceTransform->GetOutput();
	//vtkImageData * ImageSave = vtkImageData::New() ;
	//ImageSave->SetScalarTypeToFloat() ;
	//ImageSave->SetDimensions( DistanceTransform->GetOutput()->GetDimensions() ) ;
	//ImageSave->SetExtent( DistanceTransform->GetOutput()->GetExtent() ) ;
	//ImageSave->SetNumberOfScalarComponents( DistanceTransform->GetOutput()->GetNumberOfScalarComponents() ) ;
 //   //	ImageSave->SetOrigin( ImageReader->GetOutput()->GetOrigin() ) ;
	//zxhlfloat origin[4] = {0,0,0,0} ;
	//ImageSave->SetOrigin( origin ) ; // because ImageReader hasn't implemented properly the origin reader
	//ImageSave->SetSpacing( DistanceTransform->GetOutput()->GetSpacing() ) ;
	//ImageSave->Update() ;
	VtkDataConvertTo( DistanceTransform->GetOutput(), (double)0.0f, pDistance, pImage->GetDimension() );
	zxhImageInfo imageinfo ;
	pImage->GetImageInfo( &imageinfo ) ;
	pDistance->SetImageOrientationInfo( &imageinfo ) ;
	if( ReturnSquare == false )
	{
		long int lsize = pDistance->GetNumberOfPixels();
		for( int ipos=0; ipos<lsize; ++ipos )
			pDistance->SetImageData(ipos, sqrt(pDistance->GetImageData(ipos)) ) ;
	}
	//ImageSave->Delete();
	DistanceTransform->Delete();
	Threshold->Delete();
	pvtkImage->Delete();
	return true ;
#else
	std::cerr<<"error: build without vtk support, can not compute distance transform\n" ; 
	return false ;
#endif
}


void	ComputeNormalVectorUsingDistanceTransform( zxhImageData * pImageRegion,
												  zxhImageDataT<float>*pNV )
{
	zxhImageDataT<float>* pNVminus = 0 ;
	pNV->CloneTo( pNVminus );
	vtkEuclideanDistanceTransform( pImageRegion, pNV, 1, 32000, false, false ) ;
	vtkEuclideanDistanceTransform( pImageRegion, pNVminus, 1, 32000, true, false ) ;

	long int lsize = pNV->GetNumberOfPixels() ;

	for( long int ip=0; ip<lsize; ++ip )
		pNV->SetImageData( ip, pNV->GetImageData(ip)-pNVminus->GetImageData(ip) ) ;
	delete pNVminus;
};

}//end namespace zxh


