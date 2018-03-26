
#ifndef zxhImageGipl_cpp
#define zxhImageGipl_cpp

#include "zxhImageGipl.h"
#include <iostream>

namespace zxh{

bool ReadHeader(zxhImageInfo & imageinfo, std::string strFilename, bool showinfo )
{
	if( strcmp( zxh::GetExtension(strFilename).c_str(), "gipl" ) == 0 )
		return zxhImageGipl::ReadHeader( imageinfo, strFilename, showinfo ) ;
	else if( strcmp( zxh::GetExtension(strFilename).c_str(), "nii" ) == 0 )
		return zxhImageNifti::ReadHeader( imageinfo, strFilename, showinfo ) ;
	else if( strcmp( zxh::GetExtension(strFilename).c_str(), "nii.gz" ) == 0 ||
			 (strcmp( zxh::GetExtension(strFilename).c_str(), "gz" ) == 0 && strstr( strFilename.c_str(), "nii.gz" )!=NULL ) )
		return zxhImageNifti::ReadHeader( imageinfo, strFilename, showinfo ) ;
	return false ;
};
}//end of namespace

#endif
