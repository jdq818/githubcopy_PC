
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhLibFunctions.h    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0 $

=========================================================================*/
#ifndef zxhLibFunctions_h
#define zxhLibFunctions_h

#include <iostream>
#include <fstream>
#include <math.h>

#ifdef HAS_VTK
#include "vtkImageData.h"
#include "vtkImageThreshold.h"
#include "vtkImageEuclideanDistance.h"
#include "vtkImageProgressIterator.h"
#endif

#include "zxh.h"
#include "zxhImageData.h"
#include "zxhImageGipl.h"
//#include "zxhImageModelingBase.h"

namespace zxh
{
#ifdef HAS_VTK

/// vtktype should be the consistent as vtkDataType
template <class pixeltype, class vtktype>
bool	ConvertToVtkData( zxhImageDataT<pixeltype> * pInImage,
						 vtkImageData * pOutImage, int vtkDataType, vtktype type)
{
	if( pOutImage == 0 )
		pOutImage = vtkImageData::New() ;
	int size[] = {1,1,1,1} ;
	float sp[] = {1,1,1,1} ;
	pInImage->GetImageSize( size[0], size[1], size[2], size[3] );
	pInImage->GetImageSpacing( sp[0], sp[1], sp[2], sp[3] ) ;
	/// double only for vtk
	double spacing[ImageDimensionMax] = { sp[0], sp[1], sp[2], sp[3] };
	pOutImage->SetDimensions( size );
	int ext[ImageDimensionMax*2] = { 0,size[0]-1, 0,size[1]-1, 0,size[2]-1, 0,size[3]-1 } ;
	pOutImage->SetExtent( ext );
	double org[]={0,0,0,0} ;
	pOutImage->SetOrigin( org );
	pOutImage->SetSpacing( spacing ) ;
	pOutImage->SetNumberOfScalarComponents( 1 ) ;

	switch ( vtkDataType )
	{
	case GIPL_CHAR:
		pOutImage->SetScalarTypeToChar();
		break;
	case GIPL_U_CHAR:
		pOutImage->SetScalarTypeToUnsignedChar();
		break;
	case GIPL_SHORT:
		pOutImage->SetScalarTypeToShort();
		break;
	case GIPL_U_SHORT:
		pOutImage->SetScalarTypeToUnsignedShort();
		break;
	case GIPL_FLOAT:
		pOutImage->SetScalarTypeToFloat();
		break;
	case GIPL_U_INT:
		pOutImage->SetScalarTypeToUnsignedInt();
		break;
	case GIPL_INT:
		pOutImage->SetScalarTypeToInt();
		break;
	case GIPL_DOUBLE:
		pOutImage->SetScalarTypeToDouble();
		break;
	default:
		std::cerr<<"warning: unknown data type set to float\n " ;
		pOutImage->SetScalarTypeToFloat();
		break;
	}
	pOutImage->Update() ;
	for(int t = ext[6]; t <= ext[7]; t++)
	for(int z = ext[4]; z <= ext[5]; z++)
    {
		for(int y = ext[2]; y <= ext[3]; y++)
		{
			for(int x = ext[0]; x <= ext[1]; x++)
			{
				int coo[]={x,y,z,t};
				vtktype* pFpixel = static_cast<vtktype*>( pOutImage->GetScalarPointer( coo ) ) ;
				*pFpixel = static_cast<vtktype>( pInImage->GetPixelGreyscale( x,y,z,t ) ) ;
			}
		}
    }
	return true ;
};
/// vtktype, e.g. (double)0.0, should be the consistent as vtkDataType
template <class pixeltype, class vtktype>
bool	VtkDataConvertTo( vtkImageData * pInImage, vtktype type,
						  zxhImageDataT<pixeltype>* pOutImage, int dimension )
{
	if( pOutImage == 0 )
		pOutImage = new zxhImageDataT<pixeltype>() ;
	int vtksize[ImageDimensionMax] = {1,1,1,1} ;
	double vtkspacing[ImageDimensionMax] = {1,1,1,1} ;
	pInImage->GetDimensions( vtksize );
	pInImage->GetSpacing( vtkspacing ) ;

	int size[ImageDimensionMax] = {vtksize[0], vtksize[1], vtksize[2], vtksize[3] } ;
	float spacing[ImageDimensionMax] = {vtkspacing[0], vtkspacing[1], vtkspacing[2], vtkspacing[3] } ;
	zxhImageInfo imageinfo ;
	pOutImage->NewImage( dimension, size, spacing, &imageinfo ) ;

	int *ext = pInImage->GetExtent();
	//for(int t = ext[6]; t <= ext[7]; t++)
	int t=0; //only for 3D
	for(int z = ext[4]; z <= ext[5]; z++)
    {
		for(int y = ext[2]; y <= ext[3]; y++)
		{
			for(int x = ext[0]; x <= ext[1]; x++)
			{
				int coo[]={x,y,z,t} ;
				vtktype* pFpixel = static_cast<vtktype*>( pInImage->GetScalarPointer( coo ) ) ;
				pOutImage->SetPixelByGreyscale( x,y,z,t, static_cast<pixeltype>( *pFpixel ) ) ;
			}
		}
    }
	return true ;
};

#endif
bool	vtkEuclideanDistanceTransform(
					  zxhImageData		* pImage,
					  zxhImageDataT<float>* pDistance,
					  int	iForegroundFrom, int	iForegroundTo, // inclusive
					  bool	OnForeground = true,
					  bool	ReturnSquare = false ) ;


void	ComputeNormalVectorUsingDistanceTransform( zxhImageData * pImageRegion,
												  zxhImageDataT<float>*pNV );


}//endof namespace zxh

#endif //zxhLibFunctions_h



