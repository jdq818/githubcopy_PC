
/*=========================================================================

Program:   ZXH Registration Software
Author:	 ZHUANG, Xia Hai
Date:      From 2004-01
Version:   1.0, 2.0, 2.1
V 2.1: change images to be with a matrix which transforms the image to world coordinate

=========================================================================*/

#ifndef zxhImageData_h
#define zxhImageData_h
//#include "stdafx.h" //////////////////////////////////////// for win32
#include "zxhImageInfo.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>

#ifndef ZXH_PixelTypeDefault
#define ZXH_PixelTypeDefault short
#endif 
#ifndef ZXHPixelTypeDefault
#define ZXHPixelTypeDefault ZXH_PixelTypeDefault
#endif
#ifndef ZXHDefaultImageType
#define ZXHDefaultImageType "nii.gz"
#endif

/*template <class IndexType, class PixelType>
struct zxhPoint{
IndexType x;
IndexType y;
IndexType z;
IndexType t;
PixelType value;
void zxhPoint():x(0),y(0),z(0),t(0),value(0) {} ;
};*/

///
/// \class zxhImageDataT
/// \brief: interface for the zxhImageDataT class 4D image
/// zxhImageDataT include datavalue stored in signed short []
/// can not read image from file
/// default pixel valut type is short when reading real image
/// float/struct type image only used in interpolate spline ffd
/// Do not have origin and orientation matrix YET
/// \ingroup zxhImageDataT
///

template <class PixelType = ZXH_PixelTypeDefault>
class zxhImageDataT
{
public:

	///constract
	zxhImageDataT();

	// new a image with data store
	//zxhImageDataT(int idim,int size[],float spacing[]);
	// new a image with data store
	//virtual bool NewImage(int idim,int size[],float spacing[]);

	/// new a image without data store
	virtual bool NewImage(int idim, const int size[], const float spacing[], const zxhImageInfo*pImageInfo);
	///
	virtual bool NewImage(const zxhImageInfo*pImageInfo)
	{
		return NewImage(pImageInfo->Dimension, pImageInfo->Size, pImageInfo->Spacing, pImageInfo);
	};
	/// new a image with data store
	virtual bool NewImage(PixelType*pData, const zxhImageInfo*pImageInfo);
	/// \return
	virtual void	SetDataType(zxhushort i)	{ this->m_ImageInfo.DataType = i; };
	/// \return
	virtual zxhushort	GetDataType()const		{ return this->m_ImageInfo.DataType; };

	/// \return
	virtual ~zxhImageDataT();

	/// new an image and clone
	virtual zxhImageDataT<PixelType>* CloneTo(zxhImageDataT<PixelType>* & pRet)const;
	/// 
	virtual bool CloneFrom(const zxhImageDataT<PixelType>* pSrc);

	/// copy image data, assuming same NumberOfPixels
	void MemCopyImageDataFrom(zxhImageDataT<PixelType>* src)
	{
		memcpy((void*)m_pImageData, src->GetImageData(), GetNumberOfPixels()*sizeof(PixelType));
	};



	/// only for 2D/3D/4D--- without ---   boundary check
	virtual inline bool   SetPixelByGreyscale(int x, int y, int z, int t, PixelType val);

	///
	virtual inline PixelType GetPixelGreyscaleClosest(int x, int y, int z, int t) const;
	virtual inline PixelType GetPixelGreyscaleClosest(int x, int y, int z) const;
	virtual inline PixelType GetPixelGreyscaleClosest(int x, int y) const { return GetPixelGreyscaleClosest(x, y, 0); };
	///
	virtual PixelType GetPixelGreyscaleMirror(int x, int y, int z, int t) const;

	/// only for 2D respectively in order to accelerate image,-- without ---   boundary check
	virtual inline PixelType GetPixelGreyscale(int x, int y) const
	{
		return (*(m_pImageData + x + y*m_ImageInfo.Size[0]));
	};

	/// only for 3D respectively in order to accelerate image,-- without ---   boundary check
	virtual inline PixelType GetPixelGreyscale(int x, int y, int z) const
	{
		return (*(m_pImageData + z*m_iResolution + y*m_ImageInfo.Size[0] + x));
	}
	/// only for 4D respectively in order to accelerate image,-- without ---   boundary check
	virtual inline PixelType GetPixelGreyscale(int x, int y, int z, int t) const
	{
		return (*(m_pImageData + t*m_iVolume + z*m_iResolution + y*m_ImageInfo.Size[0] + x));
	}
	/// only for 4D respectively in order to accelerate image,-- without ---   boundary check
	virtual inline PixelType GetPixelGreyscale(float x, float y = 0, float z = 0, float t = 0) const
	{
		return this->GetPixelGreyscale(zxh::round(x), zxh::round(y), zxh::round(z), zxh::round(t));
	}

	///
	virtual inline PixelType GetPixelGreyscale(int iVector[ZXH_ImageDimensionMax]) const;

	/// index < number of pixels
	virtual PixelType GetImageData(int index)	const  	{ return m_pImageData[index]; };
	/// only used to check whether have data
	virtual const PixelType * GetImageData()	const 	{ return m_pImageData; };
	/// index < number of pixels
	virtual void SetImageData(int index, PixelType v) 	{ m_pImageData[index] = v; };
	///
	virtual PixelType GetPixelGreyscaleMax() const;
	///
	virtual PixelType GetPixelGreyscaleMin() const;
	/// dimension
	virtual int 	GetDimension()		const	{ return m_ImageInfo.Dimension; };
	///
	virtual void	SetDimension(int i)			{ m_ImageInfo.Dimension = i; }

	/// get image info different from 4D: size
	virtual void GetImageSize(int &sx, int &sy, int &sz, int &st) const
	{
		sx = m_ImageInfo.Size[0]; sy = m_ImageInfo.Size[1]; sz = m_ImageInfo.Size[2]; st = m_ImageInfo.Size[3];
	};
	///
	const int * GetImageSize()		const  	{ return &m_ImageInfo.Size[0]; };
	///
	virtual long int GetNumberOfPixels()	const 	{ if (m_pImageData != 0) return m_iVolume*m_ImageInfo.Size[3]; else return 0; };
	/// get image info different from 4D:spacing
	virtual void GetImageSpacing(float&dx, float&dy, float&dz, float&dt) const
	{
		dx = m_ImageInfo.Spacing[0]; dy = m_ImageInfo.Spacing[1]; dz = m_ImageInfo.Spacing[2]; dt = m_ImageInfo.Spacing[3];
	};
	///
	const float * GetImageSpacing()	const 	{ return &m_ImageInfo.Spacing[0]; };

	/// shearing in imagetoworld matrix is ignored
	virtual void GetImageExtent(float e[]) const
	{
		m_ImageInfo.GetExtent(e);
	}
	/// return the start and end of extent in world coordinates
	virtual void GetImageExtent(float worldfrom[], float worldto[]) const
	{
		m_ImageInfo.GetExtent(worldfrom, worldto);
	}

	// open image, only implemented in derived class
	//virtual bool OpenImage(const std::string Filename) ;

	/// Release memory of m_pImageData
	virtual void ReleaseMem();
	///
	virtual bool IsEmpty()const			{ return (m_pImageData == 0); };

	/// change image grid in to world coordinate
	virtual inline void ImageGridToWorldCoordinate(float fv[ZXH_ImageDimensionMax]) const
	{
		ImageToWorld(fv);
	};
	/// change world coordinate in to image grid
	virtual inline void WorldCoordinateToImageGrid(float fv[ZXH_ImageDimensionMax]) const
	{
		WorldToImage(fv);
	};

	/// change image grid in to world coordinate
	virtual inline void ImageToWorld(float fv[ZXH_ImageDimensionMax]) const { m_ImageInfo.ImageToWorld(fv); };
	virtual inline void ImageToWorld(float &fx, float&fy, float&fz, float&ft) const
	{
		float w[] = { fx, fy, fz, ft }; m_ImageInfo.ImageToWorld(w);
		fx = w[0]; fy = w[1]; fz = w[2]; ft = w[3];
	}
	;
	/// change world coordinate in to image grid
	virtual inline void WorldToImage(float fv[ZXH_ImageDimensionMax]) const		{ m_ImageInfo.WorldToImage(fv); };
	virtual inline void WorldToImage(float &fx, float&fy, float&fz, float&ft) const
	{
		float w[] = { fx, fy, fz, ft }; m_ImageInfo.WorldToImage(w);
		fx = w[0]; fy = w[1]; fz = w[2]; ft = w[3];
	};

	/// change extent to image size for storage in memory
	virtual void ExtentToImageSize(const float e[ZXH_ImageDimensionMax], int s[ZXH_ImageDimensionMax]) const
	{
		return m_ImageInfo.GetSizeUsingExtent(e, s);
	};
	///
	virtual void GetImageSizeUsingExtent(const float e[ZXH_ImageDimensionMax], int s[ZXH_ImageDimensionMax]) const
	{
		return m_ImageInfo.GetSizeUsingExtent(e, s);
	};
	/// get new size if change to different spacing
	virtual void GetImageSizeUsingSpacing(const float sp[ZXH_ImageDimensionMax], int sz[ZXH_ImageDimensionMax]) const
	{
		return m_ImageInfo.GetSizeUsingSpacing(sp, sz);
	};

	/// \return whether the image grid coordinate is inside the image
	virtual bool	InsideImage(float fx, float fy, float fz, float ft) const
	{
		if (m_pImageData == 0) return false;
		if (fx<0 || fx>(m_ImageInfo.Size[0]) - 1 ||
			fy<0 || fy>(m_ImageInfo.Size[1]) - 1 ||
			fz<0 || fz>(m_ImageInfo.Size[2]) - 1 ||
			ft<0 || ft>(m_ImageInfo.Size[3]) - 1)
			return false;
		return true;
	}
	virtual bool	InsideImageSafe(float fx, float fy, float fz, float ft) const
	{
		if (m_pImageData == 0) return false;
		// float f = m_ImageInfo.Size[0]-1-ZXH_FloatPrecision, then f = m_ImageInfo.Size[0]-1, due to limited figures of float 
		//if( zxh::absf( fx ) < ZXH_FloatPrecision ) 
		//	fx = 0 ; 
		//if( zxh::absf( fx - (m_ImageInfo.Size[0]-1)) < ZXH_FloatPrecision  )
		//	fx = (m_ImageInfo.Size[0]-1) ; 
		//if( zxh::absf( fy ) < ZXH_FloatPrecision ) 
		//	fy = 0 ; 
		//if( zxh::absf( fy - (m_ImageInfo.Size[1]-1)) < ZXH_FloatPrecision  )
		//	fy = (m_ImageInfo.Size[1]-1) ;
		//if( zxh::absf( fz ) < ZXH_FloatPrecision ) 
		//	fz = 0 ; 
		//if( zxh::absf( fz - (m_ImageInfo.Size[2]-1)) < ZXH_FloatPrecision  )
		//	fz = (m_ImageInfo.Size[2]-1) ;
		//if( zxh::absf( ft ) < ZXH_FloatPrecision ) 
		//	ft = 0 ; 
		//if( zxh::absf( ft - (m_ImageInfo.Size[3]-1)) < ZXH_FloatPrecision  )
		//	ft = (m_ImageInfo.Size[3]-1) ;
		float precision = 0.001;
		if (m_ImageInfo.Dimension == 3 && (
			(fx<-precision*m_ImageInfo.Spacing[0] || fx - precision*m_ImageInfo.Spacing[0]>(m_ImageInfo.Size[0] - 1)) ||
			(fy<-precision*m_ImageInfo.Spacing[1] || fy - precision*m_ImageInfo.Spacing[1]>(m_ImageInfo.Size[1] - 1)) ||
			(fz<-precision*m_ImageInfo.Spacing[2] || fz - precision*m_ImageInfo.Spacing[2]>(m_ImageInfo.Size[2] - 1))
			))
			return false;
		else if (m_ImageInfo.Dimension == 2 && (
			(fx<-precision*m_ImageInfo.Spacing[0] || fx - precision*m_ImageInfo.Spacing[0]>(m_ImageInfo.Size[0] - 1)) ||
			(fy<-precision*m_ImageInfo.Spacing[1] || fy - precision*m_ImageInfo.Spacing[1]>(m_ImageInfo.Size[1] - 1))
			))
			return false;
		else if (m_ImageInfo.Dimension == 4 && (
			(fx<-precision*m_ImageInfo.Spacing[0] || fx - precision*m_ImageInfo.Spacing[0]>(m_ImageInfo.Size[0] - 1)) ||
			(fy<-precision*m_ImageInfo.Spacing[1] || fy - precision*m_ImageInfo.Spacing[1]>(m_ImageInfo.Size[1] - 1)) ||
			(fz<-precision*m_ImageInfo.Spacing[2] || fz - precision*m_ImageInfo.Spacing[2]>(m_ImageInfo.Size[2] - 1)) ||
			(ft<-precision*m_ImageInfo.Spacing[3] || ft - precision*m_ImageInfo.Spacing[3]>(m_ImageInfo.Size[3] - 1))
			))
			return false;
		return true;
	}
	/// \return whether the image grid coordinate is inside the image
	virtual bool	InsideImage(int fx, int fy, int fz, int ft) const
	{
		if (m_pImageData == 0) return false;
		if (fx<0 || fx>(m_ImageInfo.Size[0]) - 1 ||
			fy<0 || fy>(m_ImageInfo.Size[1]) - 1 ||
			fz<0 || fz>(m_ImageInfo.Size[2]) - 1 ||
			ft<0 || ft>(m_ImageInfo.Size[3]) - 1)
			return false;
		return true;
	}
	/// 
	virtual void	GetClosestPixelCoordinateWithinImage(float&fx, float&fy, float&fz, float&ft) const
	{
		if (fx<0) fx = 0; // float f = m_ImageInfo.Size[0]-1-ZXH_FloatPrecision, then f = m_ImageInfo.Size[0]-1, due to limited figures of float 
		if (fy<0) fy = 0;
		if (fz<0) fz = 0;
		if (ft<0) ft = 0;
		if (fx>m_ImageInfo.Size[0] - 1) fx = m_ImageInfo.Size[0] - 1;
		if (fy>m_ImageInfo.Size[1] - 1) fy = m_ImageInfo.Size[1] - 1;
		if (fz>m_ImageInfo.Size[2] - 1) fz = m_ImageInfo.Size[2] - 1;
		if (ft>m_ImageInfo.Size[3] - 1) ft = m_ImageInfo.Size[3] - 1;
	}
	/// \return whether the image grid coordinate is inside the image
	virtual bool	InsideImage(float * fgrid) const
	{
		for (int idim = 0; idim<m_ImageInfo.Dimension; ++idim)
			if (fgrid[idim]<0 || fgrid[idim]> m_ImageInfo.Size[idim] - 1)
				return false;
		return true;
	}
	/// \return whether the physical coordinate is inside the image
	virtual bool	InsideImageWorld(const float fx, const float fy, const float fz, const float ft) const
	{
		float fv[] = { fx, fy, fz, ft };
		this->WorldToImage(fv);
		return this->InsideImage(fv);
	};
	///
	virtual bool	InsideImageWorld(const float *f) const
	{
		float fv[] = { f[0], f[1], f[2], f[3] };
		this->WorldToImage(fv);
		return this->InsideImage(fv);
	};

	/// \return whether computing grey level Centre on image grid
	virtual bool	ComputeGreyCentre(bool forcecompute = true);
	/// image grid coordinate
	virtual void GetGreyCentre(float &gcx, float &gcy, float &gcz, float &gct) const
	{
		gcx = m_afGreyLevelCentre[0]; gcy = m_afGreyLevelCentre[1]; gcz = m_afGreyLevelCentre[2]; gct = m_afGreyLevelCentre[3];
	};

	/// copy the image orientation image and set to corresponding
	/// typerefpoint = 0 bottom corner (0,0,0), 1 center, 2 upper corner
	virtual void SetImageOrientationInfo(const zxhImageInfo *pSource, int typerefpoint = 0);
	///
	virtual bool CopyIntensityRescaleInfoFrom(const zxhImageInfo *pSource)
	{
		return this->m_ImageInfo.CopyIntensityRescaleInfoFrom(pSource);
	}
	/// get image info object for getting or setting info
	virtual void GetImageOrientationInfo(zxhImageInfo&ret) const
	{
		ret.CopyOrientationInfoFrom(&m_ImageInfo);
	};

	/// only test spacing and image size, and dimension
	virtual bool	ImageNvoxelDimensionSame(const zxhImageDataT<PixelType>* pTest) const;

	///
	virtual void SetExtensionFlag(int i) 			{ m_iExtensionFlag = i; };
	///
	virtual int GetExtensionFlag(void) 	const 	{ return m_iExtensionFlag; };
	///
	virtual short int GetImageOrient() 	const	{ return m_iImageOrient; };
	///
	virtual void SetImageOrient(short int c) 		{ m_iImageOrient = c; };

	///
	bool	Is2DSlice()	const						{ return (m_ImageInfo.Size[2] == 1); };
	/// bFillBg==false, closest ; bFillBg==true, vbg; add slices on the FROM corner (afrom) to TO corner (ato), equivalent to -addroi afrom ato
	bool	ExtendRoi(const int addroifrom[], const int addroito[], bool bFillBg, PixelType vBg, float twoD2threeDsp);

	///
	virtual void	SetImageFilename(std::string filename)		{ m_ImageInfo.FileName = filename; };
	///
	virtual std::string GetImageFilename()	const				{ return m_ImageInfo.FileName; };

	///
	const zxhImageInfo * GetImageInfo() const					{ return &m_ImageInfo; };
	/// clone the image info to p (pointer to an object, not null)
	virtual void 	GetImageInfo(zxhImageInfo*p)const		{ if (p == 0){ std::cerr << "error: get image info\n"; return; }m_ImageInfo.CloneTo(p); };
	/// will clone the imageinfo -- be careful
	/// if only to set the orientation, please use Get/SetImageOrientationInfo
	virtual void	SetImageInfo(zxhImageInfo*p) 				{ if (p == 0) return; zxhImageInfo *pii = &m_ImageInfo; p->CloneTo(pii); };
	///
	virtual void	SetImageTimeSpacing(float f)				{ m_ImageInfo.Spacing[3] = f; };
	///if f[id]==0, use original
	virtual bool	SetImageSpacing(float f[])
	{
		if (f == 0) return false;
		for (int id = 0; id<ZXH_ImageDimensionMax; ++id)
			if (zxh::abs(f[id])>ZXH_FloatPrecision)
				m_ImageInfo.Spacing[id] = f[id];
		return true;
	};
protected:
	/// pointer to image data array
	PixelType	*m_pImageData;

	/// Assistant info of the image for access acceleration, 2D slice
	unsigned long				m_iResolution; //size[0]*size[1]

	/// 3D volume
	unsigned long				m_iVolume;		//size[0]*size[1]*size[2]

	/// grey level Centre, image coordinate
	float			m_afGreyLevelCentre[ZXH_ImageDimensionMax];

	///

	///
	zxhImageInfo	m_ImageInfo;
	// data type   ----------- ImageInfo


	/// for gipl
	short int			m_iImageOrient;
	/// four bytes
	int             m_iExtensionFlag;

};
////////////////////////////////////////////////// for zxhImageDataT class ///////////////////////////////////////
template<class PixelType>
zxhImageDataT<PixelType>::zxhImageDataT()
{
	m_ImageInfo.Dimension = 1;
	m_pImageData = 0;
	m_iResolution = 1;
	m_iVolume = 1;
	for (int i = 0; i<ZXH_ImageDimensionMax; i++)
	{
		m_iImageOrient = 0;
		m_afGreyLevelCentre[i] = 0;
	}
	m_ImageInfo.ByteOfDataType = sizeof(PixelType);
	if (m_ImageInfo.ByteOfDataType == 1)
		m_ImageInfo.DataType = GIPL_CHAR;
	if (m_ImageInfo.ByteOfDataType == 2)
		m_ImageInfo.DataType = GIPL_SHORT;
	if (m_ImageInfo.ByteOfDataType == 4)
		m_ImageInfo.DataType = GIPL_FLOAT;
	m_iExtensionFlag = 815;
	m_afGreyLevelCentre[0] = ZXH_InfiniteLargeFloat;
};

template<class PixelType> bool
zxhImageDataT<PixelType>::ExtendRoi(const int afrom[], const int ato[], bool bFillBg, PixelType vBg, float twoD2threeDsp)
{
	zxhImageInfo newimageinfo;
	this->GetImageInfo(&newimageinfo);

	if (this->Is2DSlice())
	{
		newimageinfo.Spacing[2] = twoD2threeDsp;
		newimageinfo.UpdateOrientationInfo(-2);
		m_ImageInfo.Spacing[2] = twoD2threeDsp;
		m_ImageInfo.UpdateOrientationInfo(-2);
	}
	zxhImageDataT<PixelType> imgRoi;
	int size[] = { 1, 1, 1, 1 };
	this->GetImageSize(size[0], size[1], size[2], size[3]);
	for (int idim = 0; idim<4; ++idim)
	{
		size[idim] = size[idim] + afrom[idim] + ato[idim];
		if (size[idim] < 1)
		{
			std::cerr << "warning: size of dimension " << idim << " will be less than 1 if ExtendROI, hence not action be applied\n";
			return false;
		}
	}
	int idimension = 4;
	if (size[3] == 1) idimension = 3;
	if (size[3] == 1 && size[2] == 1) idimension = 2;

	float leftcorner[] = { -afrom[0], -afrom[1], -afrom[2], -afrom[3] };

	this->ImageToWorld(leftcorner);
	for (int id = 0; id<3; ++id)
		newimageinfo.ImageToWorldMatrix[id][3] = leftcorner[id];
	newimageinfo.UpdateOrientationInfo(3);

	if (glbVerboseOutput >= 2)
		std::cout << "debug: \t resize " << size[0] << "x" << size[1] << "x" << size[2] << "x" << size[3] << "\n";

	imgRoi.NewImage(idimension, size, &newimageinfo.Spacing[0], &newimageinfo);

	for (int it = 0; it<size[3]; ++it)
		for (int iz = 0; iz<size[2]; ++iz)
			for (int iy = 0; iy<size[1]; ++iy)
				for (int ix = 0; ix<size[0]; ++ix)
				{
					if (this->InsideImage(ix - afrom[0], iy - afrom[1], iz - afrom[2], it - afrom[3]))
					{
						imgRoi.SetPixelByGreyscale(ix, iy, iz, it, this->GetPixelGreyscale(ix - afrom[0], iy - afrom[1], iz - afrom[2], it - afrom[3]));
					}
					else if (bFillBg == false) // closet point
					{
						imgRoi.SetPixelByGreyscale(ix, iy, iz, it, this->GetPixelGreyscaleClosest(ix - afrom[0], iy - afrom[1], iz - afrom[2], it - afrom[3]));
					}
					else imgRoi.SetPixelByGreyscale(ix, iy, iz, it, vBg);
				}
	this->NewImage(idimension, size, &newimageinfo.Spacing[0], &newimageinfo);
	for (int ip = 0; ip<imgRoi.GetNumberOfPixels(); ++ip)
		this->SetImageData(ip, imgRoi.GetImageData(ip));
	return true;
}

;
/// new a image with data store, and update the matrices using quatern parameters
template<class PixelType>
bool zxhImageDataT<PixelType>::NewImage(int idim, const int size[], const float spacing[], const zxhImageInfo*pImageInfo)
{
	int nPixel = 1;
	for (int id = 0; id<idim; ++id)
		nPixel *= size[id];
	if (GetNumberOfPixels() != nPixel ||
		size[0] * size[1] * size[2] * size[3] == 1)
	{
		if (m_pImageData != 0)
			delete[]m_pImageData;
		m_pImageData = new PixelType[nPixel];
	}
	for (int i = 0; i<ZXH_ImageDimensionMax; i++)
	{
		m_ImageInfo.Spacing[i] = 1.0f;
		m_ImageInfo.Size[i] = 1;
		m_iImageOrient = 0;
		m_afGreyLevelCentre[i] = 0;
	}
	m_afGreyLevelCentre[0] = ZXH_InfiniteLargeFloat;
	m_iExtensionFlag = 815;

	m_ImageInfo.Dimension = idim;
	for (int i = 0; i<m_ImageInfo.Dimension; i++)
	{
		m_ImageInfo.Spacing[i] = spacing[i];
		m_ImageInfo.Size[i] = size[i];
	}
	m_iResolution = m_ImageInfo.Size[0] * m_ImageInfo.Size[1];
	m_iVolume = m_iResolution*m_ImageInfo.Size[2];
	memset((void*)m_pImageData, 0, sizeof(PixelType)*m_iVolume*m_ImageInfo.Size[3]);
	this->m_ImageInfo.CopyOrientationInfoFrom(pImageInfo);
	this->m_ImageInfo.UpdateOrientationInfo(-2);
	//this->m_ImageInfo.CopyIntensityRescaleInfoFrom( pImageInfo ) ; -- cause problem
	return true;
};
template<class PixelType>
bool zxhImageDataT<PixelType>::NewImage(PixelType*pData, const zxhImageInfo*pImageInfo)
{
	if (m_pImageData != 0)
		delete[]m_pImageData;
	m_pImageData = pData;
	for (int i = 0; i<ZXH_ImageDimensionMax; i++)
	{
		m_ImageInfo.Spacing[i] = 1.0f;
		m_ImageInfo.Size[i] = 1;
		m_iImageOrient = 0;
		m_afGreyLevelCentre[i] = 0;
	}
	m_afGreyLevelCentre[0] = ZXH_InfiniteLargeFloat;
	m_iExtensionFlag = 815;

	m_ImageInfo.Dimension = pImageInfo->Dimension;
	for (int i = 0; i<m_ImageInfo.Dimension; i++)
	{
		m_ImageInfo.Spacing[i] = pImageInfo->Spacing[i];
		m_ImageInfo.Size[i] = pImageInfo->Size[i];
	}
	m_iResolution = m_ImageInfo.Size[0] * m_ImageInfo.Size[1];
	m_iVolume = m_iResolution*m_ImageInfo.Size[2];
	this->m_ImageInfo.CopyOrientationInfoFrom(pImageInfo);
	this->m_ImageInfo.UpdateOrientationInfo(-2);
	//this->m_ImageInfo.CopyIntensityRescaleInfoFrom( pImageInfo ) ; -- cause problem
	return true;
};


/// only for 2D/3D/4D--- without ---   boundary check
template<class PixelType>
inline bool zxhImageDataT<PixelType>::SetPixelByGreyscale(int x, int y, int z, int t, PixelType val)
{
	if (m_pImageData)
	{
		long int ipos = x + y*m_ImageInfo.Size[0] + z*m_iResolution + t*m_iVolume;
		*(m_pImageData + ipos) = val;
		return true;
	}
	std::cerr << "error: zxhImageData does not have m_pImageData for SetPixelByGreyscale \n";
	exit(1);
	return false;
}

///
template<class PixelType>
PixelType zxhImageDataT<PixelType>::GetPixelGreyscaleClosest(int x, int y, int z, int t) const
{
	if (x<0) x = 0;
	if (y<0) y = 0;
	if (z<0) z = 0;
	if (t<0) t = 0;
	if (x>m_ImageInfo.Size[0] - 1) x = m_ImageInfo.Size[0] - 1;
	if (y>m_ImageInfo.Size[1] - 1) y = m_ImageInfo.Size[1] - 1;
	if (z>m_ImageInfo.Size[2] - 1) z = m_ImageInfo.Size[2] - 1;
	if (t>m_ImageInfo.Size[3] - 1) t = m_ImageInfo.Size[3] - 1;
	return (*(m_pImageData + t*m_iVolume + z*m_iResolution + y*m_ImageInfo.Size[0] + x));
};
///
template<class PixelType>
PixelType zxhImageDataT<PixelType>::GetPixelGreyscaleClosest(int x, int y, int z) const
{
	if (x<0) x = 0;
	if (y<0) y = 0;
	if (z<0) z = 0;
	if (x>m_ImageInfo.Size[0] - 1) x = m_ImageInfo.Size[0] - 1;
	if (y>m_ImageInfo.Size[1] - 1) y = m_ImageInfo.Size[1] - 1;
	if (z>m_ImageInfo.Size[2] - 1) z = m_ImageInfo.Size[2] - 1;
	return (*(m_pImageData + z*m_iResolution + y*m_ImageInfo.Size[0] + x));
};
///
template<class PixelType>
PixelType zxhImageDataT<PixelType>::GetPixelGreyscaleMirror(int x, int y, int z, int t)const
{
	if (x<0) x = -x;
	if (y<0) y = -y;
	if (z<0) z = -z;
	if (t<0) t = -t;
	if (x>m_ImageInfo.Size[0] - 1) x = 2 * m_ImageInfo.Size[0] - 2 - x;
	if (y>m_ImageInfo.Size[1] - 1) y = 2 * m_ImageInfo.Size[1] - 2 - y;
	if (z>m_ImageInfo.Size[2] - 1) z = 2 * m_ImageInfo.Size[2] - 2 - z;
	if (t>m_ImageInfo.Size[3] - 1) t = 2 * m_ImageInfo.Size[3] - 2 - t;
	if (InsideImage(x, y, z, t))
		return (*(m_pImageData + t*m_iVolume + z*m_iResolution + y*m_ImageInfo.Size[0] + x));
	else return GetPixelGreyscaleMirror(x, y, z, t);
};

template<class PixelType>
inline PixelType zxhImageDataT<PixelType>::GetPixelGreyscale(int iVector[ZXH_ImageDimensionMax])const
{
	switch (this->GetDimension())
	{
	case 2:return this->GetPixelGreyscale(iVector[0], iVector[1]);
		break;
	case 3:return this->GetPixelGreyscale(iVector[0], iVector[1], iVector[2]);
		break;
	case 4:return this->GetPixelGreyscale(iVector[0], iVector[1], iVector[2], iVector[3]);
		break;
	default:
		std::cerr << "dimension error:" << this->GetDimension() << "\n";
		break;
	}
	return 0;
}

/// copy the image orientation image and set to corresponding
/// typerefpoint = 0 bottom corner (0,0,0), 1 center, 2 upper corner
template<class PixelType>
void zxhImageDataT<PixelType>::SetImageOrientationInfo(const zxhImageInfo *pSource, int typerefpoint)
{
	if (pSource == 0) return;
	this->m_ImageInfo.SetOrientationQuaternFrom(pSource);
	float refpointnew[] = { 0, 0, 0, 0 }, refpointsrc[] = { 0, 0, 0, 0 };
	if (typerefpoint == 1)
	{
		for (int i = 0; i<pSource->Dimension; ++i)
		{
			refpointnew[i] = (m_ImageInfo.Size[i] - 1.0) / 2.0;
			refpointsrc[i] = (pSource->Size[i] - 1.0) / 2.0;
		}
	}
	if (typerefpoint == 2)
	{
		for (int i = 0; i<pSource->Dimension; ++i)
		{
			refpointnew[i] = (m_ImageInfo.Size[i] - 1.0);
			refpointsrc[i] = (pSource->Size[i] - 1.0);
		}
	}
	m_ImageInfo.ImageToWorld(refpointnew);
	pSource->ImageToWorld(refpointsrc);

	for (int i = 0; i<pSource->Dimension; ++i)
	{
		m_ImageInfo.ImageToWorldMatrix[i][3] += (refpointsrc[i] - refpointnew[i]);
	}
	m_ImageInfo.UpdateOrientationInfo(3);
	return;
}
/// \return whether computing grey level Centre on image grid
template<class PixelType>
bool zxhImageDataT<PixelType>::ComputeGreyCentre(bool forcecompute)
{
	if (m_afGreyLevelCentre[0] != ZXH_InfiniteLargeFloat && forcecompute == false)
		return false;
	zxhlfloat sumofgrey = 0;
	zxhlfloat sumofgreyxyz[ZXH_ImageDimensionMax] = { 0, 0, 0, 0 };
	int ico[ZXH_ImageDimensionMax] = { 0, 0, 0, 0 };
	for (ico[3] = 0; ico[3]<m_ImageInfo.Size[3]; ++ico[3])
		for (ico[2] = 0; ico[2]<m_ImageInfo.Size[2]; ++ico[2])
			for (ico[1] = 0; ico[1]<m_ImageInfo.Size[1]; ++ico[1])
				for (ico[0] = 0; ico[0]<m_ImageInfo.Size[0]; ++ico[0])
				{
					for (int idim = 0; idim<m_ImageInfo.Dimension; ++idim)
					{
						sumofgreyxyz[idim] += this->GetPixelGreyscale(ico[0], ico[1], ico[2], ico[3])*ico[idim];
					}
					sumofgrey += this->GetPixelGreyscale(ico[0], ico[1], ico[2], ico[3]);
				}
	if (sumofgrey == 0)
	{
		if (glbVerboseOutput>0)
			std::cerr << "warning: sum of image " << m_ImageInfo.FileName << "grey level is zero\n";
		for (int idim = 0; idim<m_ImageInfo.Dimension; ++idim)
		{
			m_afGreyLevelCentre[idim] = 0; // m_ImageInfo.Size[idim]*m_ImageInfo.Spacing[idim]/2.0f;, hence set to the centre of the image extent
		}
		return true;
	}
	for (int idim = 0; idim<m_ImageInfo.Dimension; ++idim)
	{
		m_afGreyLevelCentre[idim] = static_cast<float>(sumofgreyxyz[idim] / (sumofgrey));
	}
	return true;
}

///
template<class PixelType>
bool zxhImageDataT<PixelType>::ImageNvoxelDimensionSame(const zxhImageDataT<PixelType>* pTest)const
{
	if (pTest == 0) return false;
	return m_ImageInfo.SameDimSizeSpacingAs(pTest->GetImageInfo());
}
template<class PixelType>
PixelType zxhImageDataT<PixelType>::GetPixelGreyscaleMax()const
{
	PixelType max = this->m_pImageData[0];
	for (unsigned long int ipos = 1; ipos < this->m_iVolume*this->m_ImageInfo.Size[3]; ++ipos)
		if (max < this->m_pImageData[ipos])
			max = this->m_pImageData[ipos];
	return max;
}
template<class PixelType>
PixelType zxhImageDataT<PixelType>::GetPixelGreyscaleMin()const
{
	PixelType min = this->m_pImageData[0];
	for (long int ipos = 1; ipos < (long int)(this->m_iVolume)*this->m_ImageInfo.Size[3]; ++ipos)
		if (min > this->m_pImageData[ipos])
			min = this->m_pImageData[ipos];
	return min;
}
template<class PixelType>
zxhImageDataT<PixelType>::~zxhImageDataT()
{
	if (m_pImageData)
		delete[] m_pImageData;
};


template<class PixelType>
bool
zxhImageDataT<PixelType>::CloneFrom(const zxhImageDataT<PixelType>* pSrc)
{
	if (pSrc == 0) return false;
	zxhImageDataT<PixelType> * p = this;
	pSrc->CloneTo(p);
	return true;
}

template<class PixelType>
zxhImageDataT<PixelType>*
zxhImageDataT<PixelType>::CloneTo(zxhImageDataT<PixelType>* &pRet)const
{
	if (pRet == 0)pRet = new zxhImageDataT();

	if (m_pImageData)
	{
		long int size = this->GetNumberOfPixels();
		if (pRet->GetNumberOfPixels() != size)
		{
			if (pRet->m_pImageData)
				delete[] pRet->m_pImageData;
			pRet->m_pImageData = 0;
		}
		if (pRet->m_pImageData == 0)
			pRet->m_pImageData = new PixelType[size];
		memcpy(pRet->m_pImageData, m_pImageData, sizeof(PixelType)*size);
	}

	zxhImageInfo * pImageInfo = &pRet->m_ImageInfo;
	this->m_ImageInfo.CloneTo(pImageInfo);
	pRet->m_iResolution = m_iResolution;
	pRet->m_iVolume = m_iVolume;
	for (int i = 0; i<ZXH_ImageDimensionMax; i++)
		pRet->m_afGreyLevelCentre[i] = m_afGreyLevelCentre[i];

	// orientation  and origin
	pRet->m_iImageOrient = m_iImageOrient;
	pRet->m_iExtensionFlag = m_iExtensionFlag;

	return pRet;
};

template<class PixelType>
void zxhImageDataT<PixelType>::ReleaseMem()
{
	if (m_pImageData)
	{
		delete[]m_pImageData;
		m_pImageData = 0;
	}
	m_iResolution = m_iVolume = m_ImageInfo.Size[0] = m_ImageInfo.Size[1] = m_ImageInfo.Size[2] = m_ImageInfo.Size[3] = 0;
};


typedef zxhImageDataT<ZXH_PixelTypeDefault> zxhImageData;
typedef zxhImageDataT<float> zxhImageDataF;


#endif //zxhImageData_h

