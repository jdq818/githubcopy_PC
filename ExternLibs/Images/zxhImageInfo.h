
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 ZHUANG, Xia Hai
  Date:      From 2010-11
  Version:	 V 2.1

=========================================================================*/

#ifndef zxhImageInfo_h
#define zxhImageInfo_h
#include "zxh.h"

#include <iostream>
#include <sstream>
///
/// \class zxhImageInfo
/// \brief: record image info ; when copy or remove orientation info, always use method(-2); 
///          Image coordinate and world coordinate conversion only for 3D 
///          4th-D World = Image*spacing[3]+Origin[3]
/// \ingroup zxhImageDataT
///
#ifndef _ZXH_DATATYPE_DEFINITION
#define _ZXH_DATATYPE_DEFINITION

// GIPL magic number
#define GIPL_MAGIC1         719555000u
#define GIPL_MAGIC2        4026526128u
#define GIPL_MAGIC_EXT            815u

 #define GIPL_MAGIC_END  3096044330  // need to swap
 #define GIPL_MAGIC2_END 2968125423

// GIPL header size
#define GIPL_HEADERSIZE   256

// GIPL filter types
#define GIPL_BINARY       1
#define GIPL_CHAR         7
#define GIPL_U_CHAR       8
#define GIPL_SHORT        15
#define GIPL_U_SHORT      16
#define GIPL_U_INT        31
#define GIPL_INT          32
#define GIPL_FLOAT        64
#define GIPL_DOUBLE       65
#define GIPL_C_SHORT      144
#define GIPL_C_INT        160
#define GIPL_C_FLOAT      192
#define GIPL_C_DOUBLE     193
#define GIPL_NONE		  1024 
enum ZXH_EnumerateDataType { ZXH_DataTypeBinary = GIPL_BINARY, ZXH_DataTypeChar=GIPL_CHAR, ZXH_DataTypeUChar=GIPL_U_CHAR,
	ZXH_DataTypeShort=GIPL_SHORT, ZXH_DataTypeUShort=GIPL_U_SHORT, ZXH_DataTypeInt=GIPL_INT, ZXH_DataTypeUInt=GIPL_U_INT,
    ZXH_DataTypeFloat=GIPL_FLOAT, ZXH_DataTypeDouble=GIPL_DOUBLE, 
	ZXH_DataTypeCShort=GIPL_C_SHORT, ZXH_DataTypeCInt=GIPL_C_INT, ZXH_DataTypeCFloat=GIPL_C_FLOAT, ZXH_DataTypeCDouble=GIPL_C_DOUBLE,
	ZXH_DataTypeNone=GIPL_NONE};

#define ANALYZEMAGIC 348
// GIPL image types
#define IT_NONE 0
#define IT_BINARY 1
#define IT_CHAR 7
#define IT_UCHAR 8
#define IT_SHORT 15
#define IT_USHORT 16
#define IT_UINT 31
#define IT_INT 32
#define IT_FLOAT 64
#define IT_DOUBLE 65
// ANALYZE image types
#define DT_NONE  0
#define DT_UNKNOWN 0
#define DT_BINARY 1
#define DT_UNSIGNED_CHAR 2
#define DT_SIGNED_SHORT 4
#define DT_SIGNED_INT 8
#define DT_FLOAT 16
#define DT_COMPLEX 32
#define DT_DOUBLE 64
#define DT_RGB 128
#define DT_ALL 255
#endif //_ZXH_DATATYPE_DEFINITION


#ifndef NIFTI_TYPE_UINT8
#define NIFTI_TYPE_UINT8           2
									   /*! signed short. */
#define NIFTI_TYPE_INT16           4
									   /*! signed int. */
#define NIFTI_TYPE_INT32           8
									   /*! 32 bit float. */
#define NIFTI_TYPE_FLOAT32        16
									   /*! 64 bit complex = 2 32 bit floats. */
#define NIFTI_TYPE_COMPLEX64      32
									   /*! 64 bit float = floag8. */
#define NIFTI_TYPE_FLOAT64        64
									   /*! 3 8 bit bytes. */
#define NIFTI_TYPE_RGB24         128
									   /*! signed char. */
#define NIFTI_TYPE_INT8          256
									   /*! unsigned short. */
#define NIFTI_TYPE_UINT16        512
									   /*! unsigned int. */
#define NIFTI_TYPE_UINT32        768
									   /*! signed long long. */
#define NIFTI_TYPE_INT64        1024
									   /*! unsigned long long. */
#define NIFTI_TYPE_UINT64       1280
									   /*! 128 bit float = long-long- double. */
#define NIFTI_TYPE_FLOAT128     1536
									   /*! 128 bit complex = 2 64 bit floats. */
#define NIFTI_TYPE_COMPLEX128   1792
									   /*! 256 bit complex = 2 128 bit floats */
#define NIFTI_TYPE_COMPLEX256   2048
#endif

class zxhImageInfo
{
public:
	///
	zxhImageInfo();

	///
	virtual ~zxhImageInfo();

	/// new an image and clone
	virtual zxhImageInfo * CloneTo(zxhImageInfo * & pRet) const;

	// different image type reader need to re-implement this function
	//virtual bool	CopyImageInfoFrom(const zxhImageInfo *pSource)

	/// remove orientation or set orientation info to identity
	virtual void	RemoveOrientationInfo() ;
	/// set orientation correspondingly using source info without guarantee consistent, need to call update afterwards
	virtual bool	CopyOrientationInfoFrom(const zxhImageInfo*pSource) ;
	///
	virtual bool	CopyIntensityRescaleInfoFrom( const zxhImageInfo*pSource) ;
	/// set orientation using quatenion paramters
	virtual bool	SetOrientationQuaternFrom(const zxhImageInfo*pSource)
	{	bool b = CopyOrientationInfoFrom(pSource);UpdateOrientationInfo(-2); return b ;};
	/// update orientation method 1,2,-2,3 using given method
	///   always assume spacing set, method1 will update all the others, while the others do not update method1
	///   set OrientationMethod at the end
	virtual void	UpdateOrientationInfo(int iOrientatationMethod) ;
	///
	virtual bool	UpdateImageInfoUsingNewSpacing( const float * spacing ) ; 
	///
	virtual bool	UpdateImageInfoByExtendRoi( const int *addRoiFrom, const int *addRoiTo ) ;

	/// whether two images have same image and voxel dimensions
	virtual bool	SameDimSizeSpacingAs( const zxhImageInfo * pTest ) const;
	/// whether two images have same image and voxel dimensions
	virtual bool	SameOrientationAs( const zxhImageInfo * pTest ) const;

	///
	virtual void	GetExtent( float e[] ) const ;
	/// return the start and end of extent in world coordinates
	virtual void 	GetExtent( float worldfrom[], float worldto[] ) const ;
	///
	virtual float	GetVolumeOfPixel() const ;
	///
	virtual int GetNumberOfPixels()	const	{ return Size[0]* Size[1]* Size[2]* Size[3]; };

	/// get new size if change to different extent
	virtual void GetSizeUsingExtent(const float e[ZXH_ImageDimensionMax], int s[ZXH_ImageDimensionMax]) const ;
	/// get new size if change to different spacing
	virtual void GetSizeUsingSpacing(const float sp[ZXH_ImageDimensionMax], int sz[ZXH_ImageDimensionMax]) const ;
	///
	virtual std::string	GetPrintString() const;
	///
	virtual std::string	GetPrintStringSimpleQuuatern() const ;  
	/// 
	virtual bool SetImageInfoSimpleQuaternFromStream( std::ifstream& ifs) ;
	///
	void WorldToImage( float fv[] ) const ;
	///
	void ImageToWorld( float fv[] ) const ;
	///
	void TemporalImageToWorld( float&f ) const ;
	///
	void TemporalWorldToImage( float&f ) const ;
	///
	void ImageToPhysical( float fv[] ) const ;
	///
	void PhysicalToImage( float fv[] ) const ;
	///
	void ProjectPhysVectorToWorldCoordinate( float &fvx, float &fvy, float &fvz ) const  ;
	///
	void ProjectWorldVectorToPhysCoordinate( float &fvx, float &fvy, float &fvz ) const ;

	/// this image should be 2D
	bool Project3DWorldVectorTo2DPlane( float *fc ) const ;

	///	image grid indices to imagedata index
	inline int		GridToIndex(int ix, int iy, int iz=0, int it=0 ) const;

	/// image data index to image grid indices
	inline void		IndexToGrid(const int index, int *ix, int *iy, int *iz=0, int*it=0 ) const;
	///
	inline void		IndexToWorld(const int index, float *w ) const;
	
	///
	bool	InsideImageWithSliceThickness( const float pvox[] ) const ; 
		
	///
	bool	InsideImage(const float fx, const float fy, const float fz, const float ft) const
	{ 
		float fgrid[] = {fx,fy,fz,ft} ;
		for(int idim=0;idim< Dimension;++idim)
			if( fgrid[idim]<0 || fgrid[idim]> Size[idim]-1 )
			return false;
		return true;
	};
	/// \return whether the physical coordinate is inside the image
	bool	InsideImageWorld(const float fx, const float fy, const float fz, const float ft) const
	{
		float fgrid[] = {fx,fy,fz,ft} ;
		this->WorldToImage( fgrid ) ;
		for(int idim=0;idim< Dimension;++idim)
			if( fgrid[idim]<0 || fgrid[idim]> Size[idim]-1 )
			return false;
		return true;
	} ;
	///analyze_to_gipl_type: Convert from Analyze/nifti to gipl image types
	static short Analyze_to_gipl_type(short analyze_type)
	{
		short gipl_type = IT_NONE;

		switch ( analyze_type )
		{
		case DT_BINARY:			gipl_type = GIPL_BINARY ; break;

		case NIFTI_TYPE_UINT8:	gipl_type = GIPL_U_CHAR ; break; //NIFTI_TYPE_UINT8
		case NIFTI_TYPE_INT8:		gipl_type = GIPL_CHAR ; break;

		case NIFTI_TYPE_INT16:	gipl_type = GIPL_SHORT ; break ;
		case NIFTI_TYPE_UINT16:	gipl_type = GIPL_U_SHORT ; break ;

		case NIFTI_TYPE_INT32:	gipl_type = GIPL_INT ; break ;
		case NIFTI_TYPE_UINT32:	gipl_type = GIPL_U_INT ; break ;

		case NIFTI_TYPE_FLOAT32:	gipl_type = GIPL_FLOAT ; break ;
		case NIFTI_TYPE_FLOAT64:	gipl_type = GIPL_DOUBLE ; break ;

		default: gipl_type = IT_NONE ; break;
		}

		return gipl_type;
	}
	///
	static short Gipl_to_analyze_type(short gipl_type)
	{
		short analyze_type = DT_NONE;

		switch ( gipl_type )
		{
		case GIPL_BINARY:	analyze_type = DT_BINARY ; break;

		case GIPL_U_CHAR:	analyze_type = NIFTI_TYPE_UINT8 ; break; //NIFTI_TYPE_UINT8
		case GIPL_CHAR:		analyze_type = NIFTI_TYPE_INT8 ; break;

		case GIPL_SHORT:		analyze_type = NIFTI_TYPE_INT16 ; break ;
		case GIPL_U_SHORT:	analyze_type = NIFTI_TYPE_UINT16 ; break ;

		case GIPL_INT:		analyze_type = NIFTI_TYPE_INT32 ; break ;
		case GIPL_U_INT:		analyze_type = NIFTI_TYPE_UINT32 ; break ;

		case GIPL_FLOAT:		analyze_type = NIFTI_TYPE_FLOAT32 ; break ;
		case GIPL_DOUBLE:	analyze_type = NIFTI_TYPE_FLOAT64 ; break ;

		default: gipl_type = DT_NONE ; break;
		}

		return analyze_type;
	}
	
	/// m[4][4]
	virtual void GetMatrixImageToPhysical( float *m ) const ; 
	/// m[4][4]
	virtual void GetMatrixPhysicalToWorld( float *m ) const ; 

	/// note that for eroade 1 pixel, the kernel should be only ONE pixel, as it will remove the current pixel
	/// return num of pixels in the offset array, which is a cube; subx/y/z is the sampling interval 
	int	ConstructIndexOffsetCube( int radiusx, int radiusy, int radiusz, int subx, int suby, int subz, int * offset ) const; 
	/// return num of pixels in the offset array, which is a sphere with radius in mm unit 
	int	ConstructIndexOffsetSphere( float fPhysRadiusMM, int* offset ) const; 
	/// return num of pixels in the offset array, which is a 2D circle with radius in mm unit 
	int	ConstructIndexOffset2DCirc( float fPhysRadiusMM, int* offset ) const; 
	/// return num of pixels (7, first one is the center) in the offset array, which is a sphere with radius in 1 pixel 
	int	ConstructIndexOffsetNeighbour( int* offset ) const; 
	/// return num of pixels (5, first one is the center) in the offset array, which is a sphere with radius in 1 pixel 
	int	ConstructIndexOffsetNeighbour2D( int* offset ) const; 


	///
	void GetOrthogonalOrientationWorldCoordinateCornersSpacingSize( float Corner0World[], float Corner1World[], float SpacingWorld[], int WorldImageSize[] ) const;
	/// recompute Spacing, Size and oreientation, preset Dimension and Quatern four parameters before using this
	void UpdateOrthogonalImageQuaternInfoUsingWorldInfo( const float* Corner0World, const float* Corner1World, const float* SpacingWorld, const int* WorldImageSize ) ;

public:
	///
	std::string		FileName;
	///
	int				Dimension;
	/// image format is related to datatype
	std::string		ImageFormat;
	/// data type, definition is based on GIPL_??? ONLY
	zxhushort		DataType;
	///
	int				ByteOfDataType ;

	/// define orientation method
	int				OrientationMethod ;

	// -------- method 1 -------- xyz = ijk *spacing
	/// image spacing info default 1, [0-2] unit is mm, [3]is temporal dimension and unit is second
	float			Spacing[ZXH_ImageDimensionMax];
	/// image size default 1, [3] is num. of phases
	int				Size[ZXH_ImageDimensionMax];

	// -------- method 2 -------- for gipl
	/// gipl origin
	float			Origin[ZXH_ImageDimensionMax];
	/// orientation rotation matrix (ORM) ;
	///		ORM[0-2][0] is x-axis, ORM[0-2][1] is y-axis, ORM[0-2][2] is z-axis
	/// 	notice that x-axis -> gipl.header[106,110,114]; y-axis -> gipl.header[122,126,130]
	/// 	ImageToWorldMatrix = [identity origin] x [ORM 0] x [scaling matrix] x [identity -(image_size-1)/2]
	float			OrientationRotationMatrix[3][3];

	// -------- method (-2) -------- for nii (most reliable method)
	/// Quaternion factor, sign of spacingZ step (< 0 is negative, improper left-hand coordinate; >= 0 is positive, proper right-hand coordinate)
	float		QuaternFactor ;
	/// Quaternion b parameter
	float		QuaternB ;
	/// Quaternion c parameter
	float 		QuaternC ;
	/// Quaternion d parameter
	float		QuaternD ;
	/// Quaternion x,y,z shift
	float		Qoffsetxyz[3] ;

	// -------- method 3 -------- [x y z ]^T = imagetoworldmatrix x [ i j k ]^T
	float			ImageToWorldMatrix[4][4] ;
	float			WorldToImageMatrix[4][4] ;

	/// intensity rescale
	/// original intensity value = image intensity * rescaleslope + rescaleintercept
	/// currently this is not  used in image post processing
	float		RescaleSlope ;
	float		RescaleIntercept ;
} ;


///	image grid indices to imagedata index
inline int	zxhImageInfo::GridToIndex(int ix, int iy, int iz, int it ) const
{
	return ( ((it*Size[2]+iz)*Size[1]+iy)*Size[0]+ix ) ;
}

/// image data index to image grid indices
inline void	zxhImageInfo::IndexToGrid(const int index, int *ix, int *iy, int *iz, int* it ) const
{
	if( ix!=0 ) *ix = index%(Size[0]*Size[1]*Size[2])%(Size[1]*Size[0])%Size[0] ;
	if( iy!=0 ) *iy = index%(Size[0]*Size[1]*Size[2])%(Size[1]*Size[0])/Size[0] ;
	if( iz!=0 ) *iz = index%(Size[0]*Size[1]*Size[2])/(Size[1]*Size[0]) ;
	if( it!=0 ) *it = index/(Size[0]*Size[1]*Size[2]) ;
}
///
inline void	zxhImageInfo::IndexToWorld(const int index, float *w ) const
{ 
	int ix,iy,iz,it ;
	ix = index%(Size[0]*Size[1]*Size[2])%(Size[1]*Size[0])%Size[0] ;
	iy = index%(Size[0]*Size[1]*Size[2])%(Size[1]*Size[0])/Size[0] ;
	iz = index%(Size[0]*Size[1]*Size[2])/(Size[1]*Size[0]) ;
	it = index/(Size[0]*Size[1]*Size[2]) ;
	w[0]=ix;w[1]=iy;w[2]=iz;w[3]=it;
	this->ImageToWorld( w ) ; 
}

#endif //zxhImageInfo_h

