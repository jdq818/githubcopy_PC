/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhImageGipl.h    $
  Language:  C++
  Date:      $Date: From  2006-11 $
  Version:   $Revision: 2.0 $

=========================================================================*/

#ifndef zxhImageGipl_h
#define zxhImageGipl_h

#include "zxhImageData.h"
#include "zxhImageNifti.h"
#include <string>

///
//namespace
//{
///
/// \class zxhImageGiplT
/// \brief image reader
///		   read gipl image, need to swap according to magic number
///        Orientation issue:
///           read: will consider header[106,110,114] and header[122,126,130] as x-axis and y-axis if header[138,142,146]==0 (method2), else header[106-154] is image to world matrix(method3)
///           write: always use method2 unless det(imagetoworldmatrix)<0 then use method3
///        ADD: read analyze format (.hdr,.img),y-dim index should be mirror flip in order to match gipl
///	       2012-4-24 zxh add: if image has rescale slop/intersept info
///                           when set to original intensity when read in, set to rescale value when save out
/// \ingroup zxhImageData
///
//typedef ZXH_PixelTypeDefault GiplType;
 
template<class PixelType=ZXH_PixelTypeDefault>
class zxhImageGiplT //:	public zxhImageDataT<PixelType>
{
public:
	/// constructor
	zxhImageGiplT(){};

	/// \return
	virtual ~zxhImageGiplT(){};
 	///
	/// 2012-4-24 zxh add: if image has rescale slop/intersept info, set to original intensity
	static bool OpenImage( zxhImageDataT<PixelType> * pImage, const std::string sFilename );

	/// only read image header
	/// \return whether recognized gipl file
	static bool ReadHeader( zxhImageInfo &imageinfo, const std::string strFilename, bool showinfo=false);
	 

	/// save giple image, call SetDataType before calling this
	/// 2012-4-24 zxh add: if image has rescale slop/intersept info, set to rescale intensity
	static	bool	SaveImage(const zxhImageDataT<PixelType>*pImg, const zxhushort type, const std::string strFilename ) ;
	///

protected:
	/// Read Analyze image format, filename should be .hdr or .img; no orientation
	static int ReadAnalyze(const std::string filename, zxhImageDataT<PixelType>*pImg, bool read_image) ;
};
/////////////////////////

template<class PixelType>
bool zxhImageGiplT<PixelType>::OpenImage( zxhImageDataT<PixelType> * pImage, const std::string sFilename )
{
	if( strcmp( zxh::GetExtension(sFilename).c_str(), "hdr")==0 || // image is analyze format
		strcmp( zxh::GetExtension(sFilename).c_str(), "img")==0 )
	{
		int iret=ReadAnalyze(sFilename, pImage, true ) ;
		if (iret<0 ) return false ;
		return true ;
	}

	// read header imageinfo
	zxhImageInfo imageinfo ;
	if( zxhImageGiplT<PixelType>::ReadHeader( imageinfo, sFilename, false ) ==false )
		return false;
	pImage->NewImage( imageinfo.Dimension, &imageinfo.Size[0], &imageinfo.Spacing[0], &imageinfo ) ;
	pImage->SetImageInfo( &imageinfo ) ;

	char bufferheader[GIPL_HEADERSIZE];
	bool bswap;
	std::ifstream ifs;
	ifs.open(sFilename.c_str(),std::ios::in|std::ios::binary);
	if(ifs.fail())
	{
		std::cerr<<"error: Open Image:"<<sFilename<<" Failed!\n";
		return false;
	}
	ifs.seekg(0,std::ios_base::beg);
	ifs.read(bufferheader,GIPL_HEADERSIZE);
	unsigned int magic	= *(unsigned int *)(bufferheader+252);
	if(magic==GIPL_MAGIC1||magic==GIPL_MAGIC2)
		bswap=false;
	else
	{
		zxh::Swap4Bytes((void*)&magic);
		if(magic==GIPL_MAGIC1||magic==GIPL_MAGIC2)
			bswap=true;
		else
		{
			std::cerr<<"error: Open Image:"<<sFilename<<" Failed due to wrong header -- magic "
				<<magic <<"\n";
			return false;
		}
	}

	short int imageorient = *(short int *)&bufferheader[186];
	if( bswap ) zxh::Swap2Bytes( (void*)&imageorient ) ;
	pImage->SetImageOrient( imageorient ) ;
	int extensionflag = *(int*)&bufferheader[244] ;
	if( bswap ) zxh::Swap4Bytes( (void*)&extensionflag ) ;
	pImage->SetExtensionFlag(extensionflag) ;

	// read data
	void* bufferdata=0;
	char* pbufferdata=0;
	unsigned short typesize;
	long int iVolume4D = imageinfo.Size[0]*imageinfo.Size[1]*imageinfo.Size[2]*imageinfo.Size[3] ;
	// store data
	switch(imageinfo.DataType)
	{
	case GIPL_CHAR:
		bufferdata=new signed char[iVolume4D];typesize=sizeof(signed char);
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_CHAR\n";
		break;
	case GIPL_U_CHAR:bufferdata=new unsigned char[iVolume4D];typesize=sizeof(unsigned char);
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_CHAR\n";
		break;
	case GIPL_SHORT:bufferdata=new signed short[iVolume4D];typesize=sizeof(signed short);
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_SHORT\n";
		break;
	case GIPL_U_SHORT:bufferdata=new unsigned short[iVolume4D];typesize=sizeof(unsigned short);
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_SHORT\n";
		break;
	case GIPL_FLOAT:bufferdata=new float[iVolume4D];typesize=sizeof(float);
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_FLOAT\n";
		break;
	case GIPL_U_INT: bufferdata = new unsigned int[iVolume4D]; typesize = sizeof( unsigned int ) ;
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_INT\n";
		break;
	case GIPL_INT: bufferdata = new int[iVolume4D]; typesize = sizeof( int ) ;
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_INT\n";
		break;
	case GIPL_DOUBLE: bufferdata = new float8[iVolume4D];typesize=sizeof(float8);
		//if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_DOUBLE \n";
		break;
	default:
		pImage->ReleaseMem() ;
		if(bufferdata)
			delete [] (PixelType*) bufferdata;
		std::cerr<<"error: Open Image:"<<sFilename<<" failed because unrecoginised data type !\n";
		ifs.close();
		return false;
		break;
	};
	if( typesize != sizeof( PixelType ) )
	{
		pImage->ReleaseMem() ;
		if(bufferdata)
			delete [] (PixelType*) bufferdata;
		std::cerr<<"error: Open Image:"<<sFilename<<" failed because unmatched data type from image file and the gipltemplate reader !\n";
		ifs.close();
		return false;
	}

	pbufferdata=(char*)(bufferdata);
	ifs.seekg(GIPL_HEADERSIZE);
	{
		if(ifs.eof()==true||ifs.fail()==true)
		{
			pImage->ReleaseMem() ;
			if(bufferdata)
				delete [] (PixelType*) bufferdata;
			std::cerr<<"error: Open Image:"<<sFilename<<" Failed because the image raw data size inconsistent with Header!\n";
			ifs.close();
			return false;
		}
		ifs.read( pbufferdata, iVolume4D*typesize ) ;
		if(ifs.eof()==true||ifs.fail()==true)
		{
			std::cerr<<"error: read gipl file error when reading graylevel data\n";
			pImage->ReleaseMem() ;
			if(bufferdata)
				delete [] (PixelType*) bufferdata;
			ifs.close();
			return false;
		}
		if(bswap&&typesize>1)
		{
			for(long int ipos=0; ipos<iVolume4D;++ipos)
				zxh::SwapBytes((void*)(pbufferdata+ipos*typesize),static_cast<unsigned char>(typesize));
		}
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
					for(long int ipos=0;ipos<iVolume4D;++ipos) 
						pImage->SetImageData( ipos, static_cast<PixelType>(  // problem in negative figure using int(+0.5)
													zxh::round(static_cast<PixelType>( *((PixelType*)(pbufferdata+ipos*typesize))) *slope +inter )) ) ;
				else
					for(long int ipos=0;ipos<iVolume4D;++ipos)
						pImage->SetImageData( ipos, static_cast<PixelType>(
													*((PixelType*)(pbufferdata+ipos*typesize))) );
			}; break;
		case GIPL_DOUBLE:
		case GIPL_FLOAT:
			{
				if(  slope!=1 || inter!=0)
					for(long int ipos=0;ipos<iVolume4D;++ipos)
						pImage->SetImageData( ipos, static_cast<PixelType>(
													static_cast<PixelType>( *((PixelType*)(pbufferdata+ipos*typesize))) *slope +inter ) ) ;
				else
					for(long int ipos=0;ipos<iVolume4D;++ipos)
						pImage->SetImageData( ipos, static_cast<PixelType>(
													*((PixelType*)(pbufferdata+ipos*typesize))) );
			}; break;
		default:
			{
				pImage->ReleaseMem() ;
				if(bufferdata)
					delete [] (PixelType*) bufferdata;
				std::cerr<<"error: Open Image "<<sFilename<<" Failed because data more than 2 bytes!\n";
				ifs.close();
				return false;
			}
			break;
		}
	}
	ifs.close();
	delete [] (PixelType*) bufferdata;
	if( glbVerboseOutput>0 ) std::cout<<"success: open image "<<sFilename<<" !\n";
	return true;
};
////
template<class PixelType>
bool zxhImageGiplT<PixelType>::ReadHeader( zxhImageInfo &imageinfo, const std::string strFilename, bool showinfo)
{
		// is analyze image formate ?
	if( strcmp( zxh::GetExtension(strFilename).c_str(), "hdr")==0 ||
		strcmp( zxh::GetExtension(strFilename).c_str(), "img")==0 )
	{
		zxhImageDataT<PixelType> image;
		int iret=ReadAnalyze(strFilename, &image, false ) ;
		if (iret<0 ) return false ;
		if( glbVerboseOutput>0 && showinfo==true)
		{
			image.GetImageInfo( &imageinfo ) ;
			std::cout<< imageinfo.GetPrintString() ;
			std::cout<<"success: check header over !\n";
		}
		return true ;
	}

	float pixsize[] = {1,1,1,1} ;
	unsigned short imagesize[4]={1,1,1,1};
	unsigned short datatype=0;//data type only concern char, uchar, short, zxhushort
	float origin[] = {0,0,0,0} ;
	float giplorm[12] ;
	float xaxis[]={0,0,0,0}, yaxis[]={0,0,0,0}, zaxis[]={0,0,0,0};

	char bufferheader[GIPL_HEADERSIZE];
	bool bswap;
	std::ifstream ifs;
	ifs.open(strFilename.c_str() , std::ios::in|std::ios::binary);
	if(ifs.fail())
	{
		std::cerr<<"error: open file "<<strFilename<<" failed!\n";
		return false;
	}
	ifs.seekg(0,std::ios_base::beg);
	ifs.read(bufferheader,GIPL_HEADERSIZE);
	ifs.close();
	// check header of magic number to decide whether swap

	unsigned int magic, magic1;
	magic=magic1= *(unsigned int *)(bufferheader+252);
	if(magic==GIPL_MAGIC1||magic==GIPL_MAGIC2)
		bswap=false;
	else
	{
		zxh::Swap4Bytes((void*)&magic1);
		if(magic1==GIPL_MAGIC1||magic1==GIPL_MAGIC2)
			bswap=true;
		else
		{
			std::cerr<<"error: Open Image:"<<strFilename<<" Failed due to wrong header -- magic "
				<<magic <<"\n";
			return false;
		}
	}
	//
	for(int i=0;i<4;++i)
	{
		imagesize[i]	= *(unsigned short*)&(bufferheader[i*2]);
		pixsize[i]		= *(float*)&(bufferheader[10+i*4]);
		if(bswap)
		{
			zxh::Swap2Bytes((void*)&imagesize[i]);
			zxh::Swap4Bytes((void*)&pixsize[i]);
		}
		if( imagesize[i] <=0 )
			imagesize[i] = 1 ;
		if( pixsize[i] <=0 )
			pixsize[i] = 1 ;
	}
	datatype=*(unsigned short*)&(bufferheader[8]);
	if(bswap)
		zxh::Swap2Bytes((void*)&datatype);

	for(int idim=0;idim<4;++idim)
	{
		float8 lf8 = *(float8*) &(bufferheader[204+idim*8]) ;
		if( bswap ) zxh::SwapBytes((void*)&lf8, sizeof(float8));
		origin[idim] = (float) lf8 ;
	}
	for( int i=0; i<12; ++i )
	{
		giplorm[i] = *(float*) &(bufferheader[106+i*4]) ;
		if( bswap ) zxh::SwapBytes((void*) &giplorm[i], sizeof(float)) ;
	}

	// set image info
	imageinfo.FileName = strFilename ;
	imageinfo.ImageFormat = "gipl" ;
	imageinfo.DataType = datatype ;
	imageinfo.ByteOfDataType = sizeof(PixelType) ;
	int Dimension = 1;
	if(imagesize[3]>1)
		Dimension=4;
	else if(imagesize[3]==1&&imagesize[2]>1)
		Dimension=3;
	else if(imagesize[3]==1&&imagesize[2]==1)
		Dimension=2;
	for(int idim=0;idim<Dimension;++idim)
	{
		imageinfo.Spacing[idim] = pixsize[idim] ;
		imageinfo.Size[idim] = imagesize[idim] ;
		imageinfo.Origin[idim] = origin[idim] ;
	}

	// set orientation info
	if( giplorm[0]==0 &&giplorm[1]==0 &&giplorm[2]==0
	    &&giplorm[3]==0 &&giplorm[4]==0 &&giplorm[5]==0 )
	{
		imageinfo.UpdateOrientationInfo(1);
		imageinfo.OrientationMethod = 1 ;
	}
	else if( giplorm[8]==0 &&giplorm[9]==0 &&giplorm[10]==0 ) // orientation use method2
	{
		imageinfo.OrientationMethod = 2 ;
		for( int idim=0; idim<3; ++idim )
		{
			xaxis[idim] = giplorm[idim] ;
			yaxis[idim] = giplorm[idim+4] ;
		}
		// Construct the z-axis using a right handed (`neurological') coordinate
		zaxis[0] = xaxis[1]*yaxis[2] - xaxis[2]*yaxis[1];
		zaxis[1] = xaxis[2]*yaxis[0] - xaxis[0]*yaxis[2];
		zaxis[2] = xaxis[0]*yaxis[1] - xaxis[1]*yaxis[0];
		for(int i=0; i<3; ++i )
		{
			imageinfo.OrientationRotationMatrix[i][0] = xaxis[i] ;
			imageinfo.OrientationRotationMatrix[i][1] = yaxis[i] ;
			imageinfo.OrientationRotationMatrix[i][2] = zaxis[i] ;
		}
		imageinfo.UpdateOrientationInfo(2);
	}
	else // orientation use method3
	{
		imageinfo.OrientationMethod = 3 ;
		float *pm = & imageinfo.ImageToWorldMatrix[0][0] ;
		for( int ip=0; ip<12; ++ip )
			pm[ip] = giplorm[ip] ;
		pm[12]=pm[13]=pm[14]=0;
		pm[15]=1;
		imageinfo.UpdateOrientationInfo(3);
	}
	imageinfo.Dimension = 3;
	if(imagesize[3]>1)
		imageinfo.Dimension=4;

	if( glbVerboseOutput>0 && showinfo==true)
	{
		// store data
		switch(datatype)
		{
		case GIPL_CHAR:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_CHAR\n";
			break;
		case GIPL_U_CHAR:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_CHAR\n";
			break;
		case GIPL_SHORT:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_SHORT\n";
			break;
		case GIPL_U_SHORT:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_SHORT\n";
			break;
		case GIPL_FLOAT:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_FLOAT\n";
			break;
		case GIPL_U_INT:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_INT\n";
			break;
		case GIPL_INT:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_INT\n";
			break;
		case GIPL_DOUBLE:
			if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_DOUBLE (long double 64bits) \n";
			break;
		default:
			std::cerr<<"error: unrecognized data type !\n";
			return false;
			break;
		};
		std::cout<< imageinfo.GetPrintString() ;
		std::cout<<"success: check header over !\n";
	}
	return true;
};
 
template<class PixelType>
bool zxhImageGiplT<PixelType>::SaveImage(const zxhImageDataT<PixelType>*pImg,const zxhushort type, const std::string strFilename)
{
	if(type>65||type<=0)
	{
		std::cerr<<"set idatatype before save image\n";
		return false;
	}
	int typesize ;
	switch(type)
	{
	case GIPL_CHAR: typesize = sizeof( char ) ;
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_CHAR\n";
		break;
	case GIPL_U_CHAR: typesize = sizeof( unsigned char ) ;
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_CHAR\n";
		break;
	case GIPL_SHORT: typesize = sizeof( short int ) ;
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_SHORT\n";
		break;
	case GIPL_U_SHORT: typesize = sizeof( unsigned short int ) ;
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_SHORT\n";
		break;
	case GIPL_FLOAT: typesize = sizeof( float );
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_FLOAT\n";
		break;
	case GIPL_U_INT: typesize = sizeof( unsigned int ) ;
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_U_INT\n";
		break;
	case GIPL_INT: typesize = sizeof( int ) ;
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_INT\n";
		break;
	case GIPL_DOUBLE: typesize = sizeof( float8 ) ;
		if( glbVerboseOutput>0 ) std::cout<<"success: image data type GIPL_DOUBLE (long double 64bits)\n";
		break;
	default:
		std::cerr<<"Error: unrecognized image data type\n";
		return false;
		break;
	};
	if( typesize != sizeof( PixelType ) )
	{
		std::cerr<<"error: unmatched data type from image and saving type !\n";
		return false;
	}

	char bufferheader[GIPL_HEADERSIZE];
	memset((void*)bufferheader,0,GIPL_HEADERSIZE);
	std::ofstream ofs;
	ofs.open(strFilename.c_str(),std::ios::out|std::ios::binary);
	if(ofs.fail())
	{
		std::cerr<<"error: Open File "<<strFilename<<" as gipl file name Failed!\n";
		return false;
	}
	unsigned short swap2;
	float swap4;
	float8 swap8;
	long iVolume4D = pImg->GetNumberOfPixels() ;
	if( glbVerboseOutput>=2 )
		std::cout<<"debug: save image volume"<<iVolume4D ;
	int size[] = {1, 1, 1, 1} ;
	pImg->GetImageSize( size[0], size[1], size[2], size[3] ) ;
	float spacing[] = {1, 1, 1, 1} ;
	pImg->GetImageSpacing( spacing[0], spacing[1], spacing[2], spacing[3] ) ;
	for(int i=0;i<4;++i)
	{
		swap2 = size[i];
		zxh::SwapBytes((void*)&swap2,sizeof(unsigned short));
		*((unsigned short*)&(bufferheader[i*2]))=swap2;
		swap4 = spacing[i];
		zxh::SwapBytes((void*)&swap4,sizeof(float));
		*((float*)&(bufferheader[10+i*4]))=swap4;
	}

	swap2=type;
	zxh::SwapBytes((void*)&swap2,sizeof(unsigned short));
	*(unsigned short*)&(bufferheader[8])=swap2;//datatype

	zxhImageInfo imageinfo ;
	pImg->GetImageInfo(&imageinfo) ;
	// orientation  and origin
	*(short int*)&bufferheader[186] = pImg->GetImageOrient() ;
	zxh::Swap2Bytes( (void*)&bufferheader[186] ) ;

	if( zxh::DeterminentMatrix3D( &imageinfo.ImageToWorldMatrix[0][0], 4)>0 )// this is method2 to use orientation info,
	{
		float orm_t[] = {imageinfo.OrientationRotationMatrix[0][0], imageinfo.OrientationRotationMatrix[1][0], imageinfo.OrientationRotationMatrix[2][0], 0,
					 imageinfo.OrientationRotationMatrix[0][1], imageinfo.OrientationRotationMatrix[1][1], imageinfo.OrientationRotationMatrix[2][1], 0,
					0,0,0,0} ;
		for(int ipos=0;ipos<ZXH_ImageDimensionMax*3;++ipos)
		{
			swap4= orm_t[ipos];
			zxh::SwapBytes((void*)&swap4,sizeof(float));
			((float*)(bufferheader+106))[ipos]=swap4;
		}
	}
	else // method3 to save orientation
	{
		float *pm= &imageinfo.ImageToWorldMatrix[0][0] ;
		for( int ip=0; ip<12; ++ip )
		{
			float swap4 = pm[ip];
			zxh::SwapBytes((void*)&swap4,sizeof(float));
			((float*)(bufferheader+106))[ip] = swap4;
		}
	}
	// origin may not be useful anymore in method3
	for(int idim=0;idim<ZXH_ImageDimensionMax;++idim)
	{
		swap8=static_cast<float8>(imageinfo.Origin[idim]);
		zxh::SwapBytes((void*)&swap8,sizeof(float8));
		((float8*)(bufferheader+204))[idim]=swap8;
	}

	*(int*)&bufferheader[244] = pImg->GetExtensionFlag();
	zxh::Swap4Bytes( (void*)&bufferheader[244] ) ;

	*(unsigned int *)(bufferheader+252) = (unsigned int)GIPL_MAGIC_END;//need to swap
	float8 max = pImg->GetPixelGreyscaleMax() ;
	float8 min = pImg->GetPixelGreyscaleMin() ;
	PixelType * pswap=new PixelType[iVolume4D];
	if( glbVerboseOutput>=2 )
		std::cout<<"debug: starting check swap "<< iVolume4D;

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
			if(  (slope!=1&&slope!=0) || inter!=0)
			{
				for(long int ipos=0;ipos< iVolume4D;++ipos)
					pswap[ipos]= static_cast<PixelType>( zxh::round((pImg->GetImageData(ipos)-inter)/slope) ) ;
			}
			else
			{
				for(long int ipos=0;ipos< iVolume4D;++ipos)
					pswap[ipos]=pImg->GetImageData(ipos);
			}
			break;
		}
	case GIPL_DOUBLE:
	case GIPL_FLOAT:
	//case NIFTI_TYPE_FLOAT128: 
	default:
		{
			if(  (slope!=1&&slope!=0) || inter!=0)
			{
				for(long int ipos=0;ipos< iVolume4D;++ipos)
					pswap[ipos]=(pImg->GetImageData(ipos)-inter)/slope;
			}
			else
			{
				for(long int ipos=0;ipos< iVolume4D;++ipos)
					pswap[ipos]=pImg->GetImageData(ipos);
			}
			break;
		}
	}
	if( sizeof( PixelType ) > 1 )
		for(long int ipos=0;ipos< iVolume4D;++ipos)
			zxh::SwapBytes((void*)(&pswap[ipos]),sizeof(PixelType));

	zxh::SwapBytes((void*)&min,sizeof(float8));
	zxh::SwapBytes((void*)&max,sizeof(float8));
	((float8*)(bufferheader+188))[0]=min;
	((float8*)(bufferheader+196))[0]=max;
	if( glbVerboseOutput>=2 )
		std::cout<<"debug: starting writing to file\n"  ;
	ofs.write(bufferheader,GIPL_HEADERSIZE);
	ofs.write((char*)pswap, iVolume4D*sizeof(PixelType));
	delete [] pswap;
	ofs.close();
	if( glbVerboseOutput>0 ) std::cout<<"success: Save image "<<strFilename<<"\n";

	return true;
}

template<class PixelType>
int zxhImageGiplT<PixelType>::ReadAnalyze(const std::string filename, zxhImageDataT<PixelType>*pImg, bool read_image)
{
	char header[348];
	unsigned int magic0, magic1;

	// Open stream
	std::string strHeader = zxh::GetFileNameNoExtension( filename )+".hdr";
	std::ifstream infile(strHeader.c_str(), std::ios::in | std::ios::binary);

	if (infile.fail())
	{
		std::cerr << "error: couldn't open " << strHeader <<" for read \n";
		return -1;
	}
	// Read header into data-block
	infile.read((char *)(&header[0]), 348);
	infile.close();

	// Use magic number to figure out whether byte-swapping is required
	bool bswap = false ;
	int size[] = {1,1,1,1} ;
	float spacing[] = {1,1,1,1} ;
	magic0 = magic1 = *(unsigned int *)&header[0];
	short int imageorient;
	float origin[4];

	zxh::Swap4Bytes((void*)&magic1);

	if (magic0 == ANALYZEMAGIC)
		bswap=false ;//std::cerr << "error: no byte-swapping required " << endl;
	if (magic1 == ANALYZEMAGIC)
		bswap=true ; //if (BCgverbose >= 2) cerr << " BCreadAnalyze:\tbyte-swapping for read " << endl;
	else
	{
		std::cerr <<"error: magic number didn't match ANALYZE image ("<< magic0 << ", " << magic1 << ")\n";
		return -1;
	}

	// Get image info and close header file
		// begin header
	short idatatype ;
	if (bswap)
	{
		for (int i = 0; i < ZXH_ImageDimensionMax; ++i)
		{
			short int reso = *(short *)&header[42 + i * sizeof(short)];
			zxh::Swap2Bytes((void*)&reso) ;
			size[i] = reso;
			spacing[i] = *(float *)&header[80 + i * sizeof(float)];
			zxh::Swap4Bytes((void*)&spacing[i]) ;

			// origin
			float8 orig= *(float8*)&header[204+i*sizeof(float8)] ;
			zxh::Swap8Bytes( (void*)&orig ) ;
			origin[i] = orig;
		}
		idatatype = *(short *)&header[70];
		zxh::Swap2Bytes( (void*) &idatatype ) ;

		// number of bypes of data type
		short bitpix = *(short *)&header[72];
		zxh::Swap2Bytes( (void*) &bitpix ) ;
	}
	else
	{
		for (int i = 0; i < ZXH_ImageDimensionMax; ++i)
		{
			size[i] = *(short *)&header[42 + i * sizeof(short)];
			spacing[i] = *(float *)&header[80 + i * sizeof(float)];
			// origin
			origin[i]= *(float8*)&header[204+i*sizeof(float8)] ;
		}
		idatatype = *(short *)&header[70];
		short bitpix = *(short *)&header[72];
	}
	// orientation, origin, analyze_orient
	imageorient = header[186] ;

	for (int i = 0; i < ZXH_ImageDimensionMax; ++i) // Make sure no dimensions are zero
	{
		if(size[i]<1) size[i]=1;
		if(spacing[i]<=0) spacing[i]=1;
	}

	// Initialize number of in-plane pixels and total Pixels
	long int iVolume4D = size[0]*size[1] *size[2]*size[3];
	int dimension = 0 ;
	if(size[3]>1)
		dimension=4;
	else if(size[3]==1&&size[2]>1)
		dimension=3;
	else if(size[3]==1&&size[2]==1)
		dimension=2;
	short iDataType = zxhImageInfo::Analyze_to_gipl_type(idatatype);
	if ( iDataType == IT_NONE)
	{
		std::cerr << "error: unknown image type from Analyze header \n";
		return -1 ;
	}
	pImg->SetDataType( iDataType ) ;

	// end header

	zxhImageInfo imageinfo ;
	imageinfo.FileName = filename ;
	imageinfo.ImageFormat = "hdr" ;
	imageinfo.DataType = iDataType ;
	imageinfo.ByteOfDataType = sizeof(PixelType) ;
	imageinfo.OrientationMethod = 1 ;
	int Dimension = 1;
	if(size[3]>1)
		Dimension=4;
	else if(size[3]==1&&size[2]>1)
		Dimension=3;
	else if(size[3]==1&&size[2]==1)
		Dimension=2;
	for(int idim=0;idim<Dimension;++idim)
	{
		imageinfo.Spacing[idim] = spacing[idim] ;
		imageinfo.Size[idim] = size[idim] ;
		imageinfo.Origin[idim] = origin[idim] ;
	}

	imageinfo.UpdateOrientationInfo(1);
	imageinfo.Dimension = 3;
	if(size[3]>1)
		imageinfo.Dimension=4;

	pImg->SetImageOrientationInfo(&imageinfo);
	if( read_image==false ) return 1 ;

	pImg->NewImage( dimension, size, spacing, &imageinfo ) ;
	// Read image if requested
	int typesize;
	void* bufferdata=0;
	char* pbufferdata=0;
	switch(iDataType)
	{
	case GIPL_CHAR:typesize=sizeof(signed char);
		break;
	case GIPL_U_CHAR:typesize=sizeof(unsigned char); //bufferdata=new unsigned char[iVolume4D];
		break;
	case GIPL_SHORT:typesize=sizeof(signed short); //bufferdata=new signed short[iVolume4D];
		break;
	case GIPL_U_SHORT:typesize=sizeof(unsigned short); //bufferdata=new unsigned short[iVolume4D];
		break;
	case GIPL_FLOAT:typesize=sizeof(float);  //bufferdata=new float[iVolume4D];
		break;
	case GIPL_U_INT:typesize = sizeof( unsigned int ) ;// bufferdata = new unsigned int[iVolume4D];
		break;
	case GIPL_INT:typesize = sizeof( int ) ;// bufferdata = new int[iVolume4D];
		break;
	case GIPL_DOUBLE:typesize=sizeof(float8);  
	//case NIFTI_TYPE_FLOAT128:typesize=sizeof(long double); // bufferdata = new double[iVolume4D];
		break;
	default:
		std::cerr<<"error: open image:"<<filename<<" failed because un-recognized data type !\n";
		infile.close();
		return -1;
		break;
	};
	if( typesize != sizeof( PixelType ) )
	{
		std::cerr<<"error: Open Image:"<<filename<<" failed because unmatched data type from image file and the template reader !\n";
		infile.close();
		return -1;
	}

	long int ibytes = sizeof(PixelType) * iVolume4D ;
	//if (pi->agz == BCFALSE) // Read uncompressed using ifstream
	std::string strFile = zxh::GetFileNameNoExtension( filename ) + ".img" ;
	infile.open(strFile.c_str(), std::ios::in | std::ios::binary);

	if (infile.fail())
	{
		std::cerr << "error: couldn't open " << strFile	<< " for read \n";
		return -1 ;
	}

	bufferdata = new char[ibytes] ;
	infile.read((char*)bufferdata, ibytes);
	infile.close();

	/*else if (pi->agz == BCTRUE)     // Read compressed using pipe
	{
#ifdef WIN32
		cerr << " BCreadAnalyze:\terror - cannot open compressed file " <<
			filename << ".gz for read on Win32" << endl;
		exit(1);
#else
		pipename = new char[strlen(headername) + 30];

		sprintf(pipename, "gunzip -f -c %s.img.gz", filename);
		FILE *fp = popen(pipename, "r");

		if (fp == (FILE *) 0)
		{
			cerr << " BCreadAnalyze:\terror - couldn't open " << filename
				<< ".img.gz for read " << endl;
			exit(1);
		}
		fread((void *)pi->image, sizeof(char), nbytes, fp);
		pclose(fp);
		delete[]pipename;
#endif
	}*/

	if(bswap&&typesize>1)
	{
		pbufferdata = (char*)bufferdata ;
		for(long int ipos=0; ipos<iVolume4D;++ipos)
		{
			zxh::SwapBytes((void*)(pbufferdata+ipos*typesize),static_cast<unsigned char>(typesize));
			pImg->SetImageData( ipos, (PixelType) * (PixelType*)(pbufferdata+ipos*typesize) ) ;
		}
	}
	// y-dim mirror flip
	for(int it=0;it<size[3]; ++it)
	for(int iz=0;iz<size[2]; ++iz)
	for(int ix=0;ix<size[0]; ++ix)
	{
		for(int iy=0; iy<size[1]/2; ++iy)
		{
			int indy = size[1]-1-iy ;
			PixelType tmpleft =pImg->GetPixelGreyscale(ix,iy,iz,it) ;
			PixelType tmpright=pImg->GetPixelGreyscale(ix,indy,iz,it) ;
			pImg->SetPixelByGreyscale(ix,iy,iz,it, tmpright);
			pImg->SetPixelByGreyscale(ix,indy,iz,it, tmpleft);
		}
	}

	// there is no info in Analyze image format to set to zxhImageInfo.RescaleSlope/RescaleIntercept
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
				for(long int ipos=0;ipos<iVolume4D;++ipos) // problem in negative figure using int(+0.5)
					pImg->SetImageData( ipos, static_cast<PixelType>( zxh::round(pImg->GetImageData(ipos)*slope +inter)) ) ; 
			break;
		}
	case GIPL_DOUBLE: // 64bit, long double
	case GIPL_FLOAT:  // 4 bit, float/double
	//case NIFTI_TYPE_FLOAT128:
	default:
		{
			if(  slope!=1 || inter!=0)
				for(long int ipos=0;ipos<iVolume4D;++ipos)
					pImg->SetImageData( ipos, static_cast<PixelType>( pImg->GetImageData(ipos)*slope +inter ) ) ;
			break;
		}
	}

	// Flip image in y-direction for compatibility
	// BCflipy(pi);
	return 0;
}

typedef zxhImageGiplT<ZXH_PixelTypeDefault> zxhImageGipl;

///////////////////////move to zxhRegistration.h later ---------------------------///////////////////////////////////////////////////////
namespace zxh{

bool ReadHeader(zxhImageInfo & imageinfo, std::string strFilename, bool showinfo ) ;

template <class PixelType>
bool OpenImage(zxhImageDataT<PixelType>*pImg, std::string strFilename)
{
	if( strcmp( zxh::GetExtension(strFilename).c_str(), "gipl" ) == 0 )
		return zxhImageGiplT<PixelType>::OpenImage( pImg, strFilename ) ;
	else if( strcmp( zxh::GetExtension(strFilename).c_str(), "nii" ) == 0 )
		return zxhImageNiftiT<PixelType>::OpenImage( pImg, strFilename ) ;
	else if( strcmp( zxh::GetExtension(strFilename).c_str(), "nii.gz" ) == 0 &&
			 strstr( strFilename.c_str(), "nii.gz" )!=NULL )
		return zxhImageNiftiT<PixelType>::OpenImage( pImg, strFilename ) ;
	return false ;
}

///OpenImageSafe to open arbitrary type image into specific type image
template <class PixelType>
bool OpenImageSafe(zxhImageDataT<PixelType>*pImg, std::string strFilename )
{
	if (strcmp(zxh::GetExtension(strFilename).c_str(), "gipl") == 0)
		return zxhImageGiplT<PixelType>::OpenImage(pImg, strFilename);
	else if ((strcmp(zxh::GetExtension(strFilename).c_str(), "nii") == 0) ||
		(strcmp(zxh::GetExtension(strFilename).c_str(), "nii.gz") == 0 && strstr(strFilename.c_str(), "nii.gz") != NULL))
	{
		if (zxhImageNiftiT<PixelType>::OpenImage(pImg, strFilename) == true)
			return true;
		else return zxhImageNiftiT<PixelType>::OpenImageSafe(pImg, strFilename);
	}
	return false;
}
template <class PixelType>
bool SaveImage(const zxhImageDataT<PixelType>*pImg, zxhushort type, std::string strFilename)
{
	if( pImg==0 )
	{
		std::cerr<<"error: try to save zero pointer image\n" ;
		return false ;
	}
	std::string ext = zxh::GetExtension(strFilename) ;
	zxh::case_lower( ext ) ;
	if( strcmp( ext.c_str(), "gipl" ) == 0 )
		return zxhImageGiplT<PixelType>::SaveImage( pImg, type, strFilename ) ;
	else if( strcmp( ext.c_str(), "nii" ) == 0 )
		return zxhImageNiftiT<PixelType>::SaveImage( pImg, type, strFilename ) ;
	else if( strcmp( ext.c_str(), "nii.gz" ) == 0 &&
			 strstr( strFilename.c_str(), "nii.gz" )!=NULL )
		return zxhImageNiftiT<PixelType>::SaveImage( pImg, type, strFilename ) ;
	// default save as nii
	if( strFilename.length()>1 && strFilename.at(strFilename.length()-1) == '.')
		return zxhImageNiftiT<PixelType>::SaveImage( pImg, type, strFilename + ZXHDefaultImageType ) ;
	else
		return zxhImageNiftiT<PixelType>::SaveImage( pImg, type, strFilename + "." + ZXHDefaultImageType ) ;
}
template <class PixelType>
bool SaveImage(const zxhImageDataT<PixelType>*pImg, std::string strFilename)
{
	zxhushort type = pImg->GetDataType() ;
	if( type <= 0 || type >= 256)
		switch( sizeof(PixelType) )
		{
		case 1:	type = GIPL_CHAR ;break;
		case 2:	type = GIPL_SHORT ;break;
		case 4:	type = GIPL_FLOAT ;break;
		case 8:	type = GIPL_DOUBLE ;break;
		case 16:	type = NIFTI_TYPE_FLOAT128 ;break;
		default: return false ;
		}
	return SaveImage( pImg, type, strFilename ) ;
}
}
//}//end of namespace
#endif //zxhImageGipl_h


