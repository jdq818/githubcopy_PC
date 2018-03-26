/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhImageGipl.h    $
  Language:  C++
  Date:      $Date: From  2010-11 $
  Version:   $Revision: 2.0, V2.1 $

=========================================================================*/

#ifndef zxhImageNifti_h
#define zxhImageNifti_h

#include "zxhImageData.h"
#include "nifti1_io.h"
#include <string>

///
//namespace
//{
///
/// \class zxhImageNiftiT
/// \brief image reader
///		   read gipl image, need to swap according to magic number
///        ADD: read analyze format (.hdr,.img),y-dim index should be mirror flip in order to match gipl
/// \ingroup zxhImageData
///

template<class PixelType=ZXH_PixelTypeDefault>
class zxhImageNiftiT
{
public:
	/// constructor
	zxhImageNiftiT(){};

	/// \return
	virtual ~zxhImageNiftiT(){};
 	///
	static bool OpenImage( zxhImageDataT<PixelType> * pImage, const std::string sFilename );
	///
	///OpenImageSafe to open arbitrary type image into specific type image
	static bool OpenImageSafe(zxhImageDataT<PixelType> * pImage, const std::string sFilename ); //, const ZXH_EnumerateDataType datatype);
	/// only read image header
	/// \return whether recognized gipl file
	static bool ReadHeader( zxhImageInfo &imageinfo, const std::string strFilename, bool showinfo=false);
	 

	/// save giple image, call SetDataType before calling this
	static	bool	SaveImage(const zxhImageDataT<PixelType>*pImg, const zxhushort type, const std::string strFilename ) ;

protected:
	static bool	SetImageInfoFromNIFTI( zxhImageInfo &imageinfo, nifti_image* pNimg ) ;
	static bool SetNIFTIFromImageInfo( nifti_image* pNimg, zxhImageInfo & imageinfo ) ;
	static std::string GetFilenameCheckGz(const std::string strFilename);
	static long int GetValueAsLongIntFromOriginalNiiStruct(const nifti_image * pNimg, const int ipos);
	static float8 GetValueAsDoubleFromOriginalNiiStruct(const nifti_image * pNimg, const int ipos);

};
/////////////////////////

template<class PixelType>
bool zxhImageNiftiT<PixelType>::SetNIFTIFromImageInfo( nifti_image* nim, zxhImageInfo & imageinfo )
{
	if( nim==0 ) return false ;
	// these not very useful, nor stored info for initialization of nifti_image
	nim->scl_slope = imageinfo.RescaleSlope ;
	nim->scl_inter = imageinfo.RescaleIntercept ;
	nim->intent_code = 0;
	nim->intent_p1 = nim->intent_p2 = nim->intent_p3 = 0;
	nim->data = NULL; // only available from images

	// size and spacing
	nim->datatype = zxhImageInfo::Gipl_to_analyze_type( imageinfo.DataType ) ; // Will be NIFTI_TYPE_UINT8 | NIFTI_TYPE_INT16 | NIFTI_TYPE_FLOAT32
	nim->nx = imageinfo.Size[0] ;
	nim->ny = imageinfo.Size[1] ;
	nim->nz = imageinfo.Size[2] ;
	nim->nt = imageinfo.Size[3] ;
	nim->dx = imageinfo.Spacing[0] ;
	nim->dy = imageinfo.Spacing[1] ;
	nim->dz = imageinfo.Spacing[2] ;
	nim->dt = imageinfo.Spacing[3] ;
	nim->nifti_type = 1 ; // should set the magic in nhdr
	nim->ndim = imageinfo.Dimension>3? imageinfo.Dimension: 3 ;
	// Default values for 3D
	nim->nu   = 1;
	nim->nv   = 1;
	nim->nw   = 1;
	nim->nbyper = imageinfo.ByteOfDataType ;
	nim->nvox	= nim->nx * nim->ny * nim->nz * nim->nt * nim->nu * nim->nv * nim->nw;

	// oreientation
	nim->sform_code = 0 ;
	nim->qform_code = 1 ; // ---- always save using quatern parameters for the orientation info
	zxh::SetMatrix( &nim->qto_xyz.m[0][0], &imageinfo.ImageToWorldMatrix[0][0], 4 ) ;
	nifti_mat44_to_quatern( nim->qto_xyz, &(nim->quatern_b), &(nim->quatern_c), &(nim->quatern_d),
			&(nim->qoffset_x), &(nim->qoffset_y), &(nim->qoffset_z), &(nim->dx), &(nim->dy), &(nim->dz), &(nim->qfac));

	// Set the units
	nim->xyz_units  = NIFTI_UNITS_MM;
	nim->time_units = NIFTI_UNITS_SEC;
	nim->toffset = 0;

	// set du/dv/dw(pixdim[5-7]), cal_min/max, slice_code/start/end to default 0
	nim->du = nim->dv = nim->dw = 0; 
	nim->cal_max = nim->cal_min = 0;
	nim->slice_code = nim->slice_start = nim->slice_end = 0; 
	return true ;
}
template<class PixelType>
bool zxhImageNiftiT<PixelType>::SetImageInfoFromNIFTI( zxhImageInfo &imageinfo, nifti_image*pNimg )
{
	if( pNimg ==0 ) return false ;

	if( pNimg->scl_slope != 0 ) 
		imageinfo.RescaleSlope = pNimg->scl_slope ;
	else imageinfo.RescaleSlope = 1 ;
	imageinfo.RescaleIntercept = pNimg->scl_inter ;

	imageinfo.FileName = pNimg->fname ;
	imageinfo.ImageFormat = "nii" ;
	imageinfo.Dimension = pNimg->dim[0] ;

	for( int i=0; i<4; ++i )
	{
		imageinfo.Spacing[i] = zxh::abs( pNimg->pixdim[i+1] ); // force to positive, LR info is in sform
		if( imageinfo.Spacing[i]==0 )
			imageinfo.Spacing[i] = 1 ;
	}
	for( int i=0; i<3; ++i )
	{
		imageinfo.Size[i] = pNimg->dim[i+1] ; // index from 1, dim[0] is Dimension
	}
	if (pNimg->dim[0] >= 4)
	{
		imageinfo.Dimension = 4 ; 
		imageinfo.Size[3] = pNimg->dim[4];
		for( int i=4; i<pNimg->dim[0]; ++i ) imageinfo.Size[3] *= pNimg->dim[i+1] ;
	}
	else
	{
		imageinfo.Size[3] = 1 ; 
	}
	if(imageinfo.Size[3]==1 && imageinfo.Dimension==4 ) imageinfo.Dimension = 3 ;
	// Check which coordinate system to use
	mat44 mat_44 ;
	int omethod = 3 ;
	if (pNimg->qform_code > 0)  // quatern
	{	mat_44 = pNimg->qto_xyz; omethod=-2; }
	else if (pNimg->sform_code > 0) // matrix
		mat_44     = pNimg->sto_xyz;
	else
	{
		zxh::MatrixIdentity( &mat_44.m[0][0], 4 ) ;
		mat_44.m[0][0] = -imageinfo.Spacing[0] ;  // always left-hand coordinate
		mat_44.m[1][1] = imageinfo.Spacing[1] ;
		mat_44.m[2][2] = imageinfo.Spacing[2] ;
		mat_44.m[0][3] = (imageinfo.Size[0]-1)/2.0 ;  // origin becomes zeros
		mat_44.m[1][3] = -(imageinfo.Size[1]-1)/2.0 ;
		mat_44.m[2][3] = -(imageinfo.Size[2]-1)/2.0 ;
		omethod=3 ;
	}
	for( int r=0; r<4; ++r )
		for( int c=0; c<4; ++c )
			imageinfo.ImageToWorldMatrix[r][c] = mat_44.m[r][c] ;
	imageinfo.UpdateOrientationInfo( 3 ) ; // using matrix
	imageinfo.OrientationMethod = omethod ;

	//  check left or right hand coordinate
	/*float det = zxh::DeterminentMatrix3D( &imageinfo.ImageToWorldMatrix[0][0], 4 ) ;
	if( glbVerboseOutput>0 )
	{
		if (det < 0.0)
			std::cout<<"success: \t nifti image is left hand coordinate\n" ;
		else std::cout<<"success: \t nifti image is right ahnd coordinate\n" ; //NIFTI_RADIOLOGICAL; else is  NIFTI_NEUROLOGICAL;
	}*/
	imageinfo.DataType = zxhImageInfo::Analyze_to_gipl_type( pNimg->datatype ) ; // all change to GIPL_
	switch (pNimg->datatype) //definition in nifti1.h DT_????? or NIFTI_TYPE_???
	{
		case NIFTI_TYPE_UINT8:
			imageinfo.ByteOfDataType = 1;
			break;
		case NIFTI_TYPE_INT8:
			imageinfo.ByteOfDataType = 1;
			break;
		case NIFTI_TYPE_INT16:
			imageinfo.ByteOfDataType = 2;
			break;
		case NIFTI_TYPE_UINT16:
			imageinfo.ByteOfDataType = 2;
			break;
		case NIFTI_TYPE_FLOAT32:
			imageinfo.ByteOfDataType = 4;
			break;
		case NIFTI_TYPE_INT32:
			imageinfo.ByteOfDataType = 4;
			break;
		case NIFTI_TYPE_UINT32:
			imageinfo.ByteOfDataType = 4;
			break;
		case NIFTI_TYPE_FLOAT64:
			imageinfo.ByteOfDataType = 8;
			break;
		default:
			std::cerr<<"error: open image:"<<imageinfo.FileName<<" failed because un-recognized data type !\n";
			return false ;
	}
	return true ;
}
template<class PixelType>
bool zxhImageNiftiT<PixelType>::OpenImage( zxhImageDataT<PixelType> * pImage, const std::string sFilename )
{
	std::string sFilenameExist = GetFilenameCheckGz(sFilename);
	nifti_image * pNimg = nifti_image_read(sFilenameExist.c_str(), 1); 
	if( pNimg == NULL )
	{
		std::cerr<<"error: Open nifti image "<<sFilename<<" failed\n" ;
		return false ;
	}
	zxhImageInfo imageinfo ;
	bool bhead = SetImageInfoFromNIFTI( imageinfo, pNimg ) ;
	if( imageinfo.ByteOfDataType != sizeof( PixelType ) )
	{
		//std::cerr<<"error: open image "<<imageinfo.FileName<<" failed because unmatched data type from image file and the template reader !\n";
		return false ;
	}
	if( bhead == false )
	{
		if( pNimg != 0 ) nifti_image_free(pNimg) ; pNimg = 0 ;
		return false ;
	} 
	pImage->NewImage( (PixelType*) pNimg->data, &imageinfo ) ; 
	pImage->SetImageInfo( &imageinfo ) ; //clone image info 
	char * pData = (char*) pNimg->data ;
	pNimg->data = NULL ; 
	long int inum = pImage->GetNumberOfPixels() ;
	// copy in pixel values, then delete pNimg 

	// rescale slope and intercept
	float slope = imageinfo.RescaleSlope ;
	float inter = imageinfo.RescaleIntercept ;
	switch(imageinfo.DataType)
	{
		case GIPL_CHAR:
		case GIPL_U_CHAR:
		case GIPL_SHORT:
		case GIPL_U_SHORT:
		case GIPL_U_INT:
		case GIPL_INT:
		{
			if(  slope!=1 || inter!=0)
				for( long int ipos = 0 ; ipos < inum ; ++ipos )
				{
					PixelType * ppix = (PixelType*) ( pData+ imageinfo.ByteOfDataType*ipos ) ;
					pImage->SetImageData( ipos, (PixelType) zxh::round((*ppix)*slope+inter) ) ; // problem in negative figure using int(+0.5)
				}
			/*else
				for( long int ipos = 0 ; ipos < inum ; ++ipos )
				{
					PixelType * ppix = (PixelType*) ( pData+ imageinfo.ByteOfDataType*ipos ) ;
					pImage->SetImageData( ipos, (PixelType) (*ppix) ) ;
				}*/
			break;
		}
		case GIPL_DOUBLE:
		case GIPL_FLOAT:
		default:
		{
			if(  slope!=1 || inter!=0)
				for( long int ipos = 0 ; ipos < inum ; ++ipos )
				{
					PixelType * ppix = (PixelType*) ( pData+ imageinfo.ByteOfDataType*ipos ) ;
					pImage->SetImageData( ipos, (PixelType) ((*ppix)*slope+inter) ) ;
				}
			/*else
				for( long int ipos = 0 ; ipos < inum ; ++ipos )
				{
					PixelType * ppix = (PixelType*) ( pData+ imageinfo.ByteOfDataType*ipos ) ;
					pImage->SetImageData( ipos, (PixelType) (*ppix) ) ;
				}*/
		};break;
	}


	if( pNimg != 0 ) nifti_image_free(pNimg); pNimg = 0 ;
	if( glbVerboseOutput>0 ) std::cout<<"success: open image "<<sFilenameExist<<" !\n";
	return true;
};

template<class PixelType>
bool zxhImageNiftiT<PixelType>::OpenImageSafe(zxhImageDataT<PixelType> * pImage, const std::string sFilename ) // 2017-07-19 update , const ZXH_EnumerateDataType datatype)
{
	std::string sFilenameExist = GetFilenameCheckGz(sFilename);
	nifti_image * pNimg = nifti_image_read(sFilenameExist.c_str(), 1);
	if (pNimg == NULL)
	{
		std::cerr << "error: Open nifti image " << sFilename << " failed\n";
		return false;
	}
	zxhImageInfo imageinfo;
	bool bhead = SetImageInfoFromNIFTI(imageinfo, pNimg); 
	imageinfo.DataType = (zxhushort) pImage->GetImageInfo()->DataType;
	imageinfo.ByteOfDataType = sizeof(PixelType);  
	if (bhead == false)
	{
		if (pNimg != 0) nifti_image_free(pNimg); pNimg = 0;
		return false;
	}
	pImage->NewImage(imageinfo.Dimension, imageinfo.Size, imageinfo.Spacing, &imageinfo);
	pImage->SetImageInfo(&imageinfo); //clone image info   
	long int inum = pImage->GetNumberOfPixels();
	// copy in pixel values, then delete pNimg 

	// rescale slope and intercept
	float slope = imageinfo.RescaleSlope;
	float inter = imageinfo.RescaleIntercept;
	switch (imageinfo.DataType)
	{
	case GIPL_CHAR:
	case GIPL_U_CHAR:
	case GIPL_SHORT:
	case GIPL_U_SHORT:
	case GIPL_U_INT:
	case GIPL_INT:
	{ 
		for (long int ipos = 0; ipos < inum; ++ipos)
		{
			long int niivalue = GetValueAsLongIntFromOriginalNiiStruct(pNimg, ipos);  
			pImage->SetImageData(ipos, (PixelType)zxh::round((niivalue)*slope + inter)); // problem in negative figure using int(+0.5)
		} 
		break;
	}
	case GIPL_DOUBLE:
	case GIPL_FLOAT:
	default:
	{ 
		for (long int ipos = 0; ipos < inum; ++ipos)
		{
			float8 niivalue = GetValueAsDoubleFromOriginalNiiStruct(pNimg, ipos);
			pImage->SetImageData(ipos, (PixelType)((niivalue)*slope + inter));
		} 
	}; break;
	}


	if (pNimg != 0) nifti_image_free(pNimg); pNimg = 0;
	if (glbVerboseOutput>0) std::cout << "success: open image " << sFilenameExist << " !\n";
	return true;
};
////
template<class PixelType>
bool zxhImageNiftiT<PixelType>::ReadHeader( zxhImageInfo &imageinfo, const std::string strFilename, bool showinfo)
{
	std::string sFilenameExist = GetFilenameCheckGz( strFilename ) ;
	nifti_image * pNimg = nifti_image_read( sFilenameExist.c_str(), 0) ;
	if( pNimg == NULL )
	{
		std::cerr<<"error: open nifti image header "<<strFilename<<" failed\n" ;
		return false ;
	}
	bool bhead = SetImageInfoFromNIFTI( imageinfo, pNimg ) ;
	if( pNimg != 0 ) nifti_image_free(pNimg) ; pNimg = 0 ;
	if( bhead == false )
	{
		std::cerr<<"error: setting imageinfo from nifti failed\n" ;
		return false ;
	}

	if( glbVerboseOutput>0 && showinfo==true)
	{
		std::cout<< imageinfo.GetPrintString() ;
		std::cout<<"success: check header over !\n";
	}
	return true;
};
 
template<class PixelType>
bool zxhImageNiftiT<PixelType>::SaveImage(const zxhImageDataT<PixelType>*pImg,const zxhushort type, const std::string strFilename)
{
	if( pImg==0 ) return false ;
	zxhImageInfo imageinfo ;
	nifti_image niimg ;
	pImg->GetImageInfo( & imageinfo ) ;
	imageinfo.DataType = type ;
	imageinfo.FileName = strFilename ;

	SetNIFTIFromImageInfo( &niimg, imageinfo ) ;

	struct nifti_1_header niihdr = nifti_convert_nim2nhdr( &niimg ) ;
	nifti_image *niisave = nifti_convert_nhdr2nim(niihdr, strFilename.c_str()); // This sets fname and iname
	niisave->iname_offset = 352;                           			// reference from irtk Some nifti versions lose this on the way!

	/// Set data pointer in nifti image struct
	long int inum = pImg->GetNumberOfPixels() ;
	PixelType *pdata = new PixelType[ inum ] ;

	float slope = imageinfo.RescaleSlope ;
	float inter = imageinfo.RescaleIntercept ;
	switch(imageinfo.DataType) // for rescale slope
	{
		case GIPL_CHAR:
		case GIPL_U_CHAR:
		case GIPL_SHORT:
		case GIPL_U_SHORT:
		case GIPL_U_INT:
		case GIPL_INT:
		{
			if(  (slope!=1 || inter!=0 ) && slope!=0 )
			{
				for( long int ipos = 0 ; ipos<inum ; ++ipos )
					pdata[ipos] = static_cast<PixelType>( zxh::round((pImg->GetImageData( ipos )-inter)/slope) ) ; // problem in negative figure using int(+0.5)
			}
			else
			{
				for( long int ipos = 0 ; ipos<inum ; ++ipos )
					pdata[ipos] = pImg->GetImageData( ipos ) ;
			}
			break;
		}
	case GIPL_DOUBLE:
	case GIPL_FLOAT:
	default:
		{
			if(  (slope!=1 || inter!=0 ) && slope!=0 )
			{
				for( long int ipos = 0 ; ipos<inum ; ++ipos )
					pdata[ipos] = static_cast<PixelType>( (pImg->GetImageData( ipos )-inter)/slope ) ;
			}
			else
			{
				for( long int ipos = 0 ; ipos<inum ; ++ipos )
					pdata[ipos] = pImg->GetImageData( ipos ) ;
			}
			break;
		}
	}

	niisave->data = pdata;
	//delete [] pdata ; will be delete in nifti_image_free(niisave);

	nifti_image_write(niisave ) ;
	/*znzFile fp = nifti_image_write_hdr2_img2(niisave, 1, "wb", NULL, NULL);
	if (fp){
		if (g_opts.debug > 2) fprintf(stderr, "-d niw: done with znzFile\n");
		free(fp);
	}
	if (g_opts.debug > 1) fprintf(stderr, "-d nifti_image_write: done\n");*/

	nifti_image_free( niisave ) ;
	if( glbVerboseOutput>0 ) std::cout<<"success: Save image "<<strFilename<<"\n";
	return true ;
}

template<class PixelType>
long int zxhImageNiftiT<PixelType>::GetValueAsLongIntFromOriginalNiiStruct(const nifti_image * pNimg, const int ipos)
{ 
	switch (pNimg->datatype) //definition in nifti1.h DT_????? or NIFTI_TYPE_???
	{
	case NIFTI_TYPE_UINT8:
	{unsigned char *pdata = (unsigned char *)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_INT8:
	{char *pdata = (char *)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_INT16:
	{short int *pdata = (short int *)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_UINT16:
	{unsigned short int *pdata = (unsigned short int *)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_FLOAT32:
	{float *pdata = (float *)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_INT32:
	{int *pdata = (int *)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_UINT32:
	{unsigned int *pdata = (unsigned int *)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_FLOAT64:
	{float8 *pdata = (float8*)pNimg->data;
	return (long int)*(pdata + ipos); }
		break;
	default: 
		return 0;
	}
}
template<class PixelType>
float8 zxhImageNiftiT<PixelType>::GetValueAsDoubleFromOriginalNiiStruct(const nifti_image * pNimg, const int ipos)
{ 
	switch (pNimg->datatype) //definition in nifti1.h DT_????? or NIFTI_TYPE_???
	{
	case NIFTI_TYPE_UINT8:
	{unsigned char *pdata = (unsigned char *)pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_INT8:
	{char *pdata = (char *)pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_INT16:
	{short int *pdata = (short int *)pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_UINT16:
	{unsigned short int *pdata = (unsigned short int *)pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_FLOAT32:
	{float *pdata = (float *)pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_INT32:
	{int *pdata = (int *)pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_UINT32:
	{unsigned int *pdata = (unsigned int *)pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	case NIFTI_TYPE_FLOAT64:
	{float8 *pdata = (float8* )pNimg->data;
	return (float8)*(pdata + ipos); }
		break;
	default:
		return 0;
	}
}

template<class PixelType>
std::string  zxhImageNiftiT<PixelType>::GetFilenameCheckGz( const std::string strFilename )
{
	if( ! zxh::FileExist( strFilename ) )
	{
		std::string sret = strFilename;
		if( zxh::GetExtension(sret).compare("nii.gz") == 0 )
		{
			sret = sret.substr( 0, strFilename.length()-3 ) ; // .gz
			if( zxh::FileExist( sret ) )
				return sret ;
		}
		else if( zxh::GetExtension(sret).compare("nii") == 0 && zxh::FileExist( sret+".gz" ) )
			return sret+".gz" ;
	}
	return strFilename ;
}

typedef zxhImageNiftiT<ZXH_PixelTypeDefault> zxhImageNifti;

//}//end of namespace
#endif //zxhImageNifti_h


