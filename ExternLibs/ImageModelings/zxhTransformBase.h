
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhTransformBase.h    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0 $

=========================================================================*/

#ifndef zxhTransformBase_H
#define zxhTransformBase_H

//#include "zxhMetricBase.h"
#include "zxhImageData.h"
#include <string.h>
#include <stdio.h>

//class zxhTransformGradientBase
//{
//public:
//	/// construction
//	zxhTransformGradientBase(){};
//	/// deconstruct
//	virtual ~zxhTransformGradientBase(){};
//};

///
/// \class zxhTransformBase
/// \brief
/// This class is used as tranformation class
/// Transformation
/// This class is an abstract class to represent the transformation, or gradient type of Metrics
/// 在Transform中，只是针对imageTest，因为所谓的transform是对test image
/// 进行transform，其实不涉及到别的image。只有在把test image transform完
/// 后的image 映射到reference image时才与之有关系。而test与ref image的关系
/// 是通过test原来的原点坐标（就是世界坐标）与ref image的原点坐标（也是世界坐标）
/// 是一样的。
/// 所以 T(x) = x' 为image 内部坐标x 变换到 内部坐标 x'
/// AdjustTestToWorldicalx'*spacing or x*spacing 为世界坐标 smm
/// smm 映射到ref image 上的坐标为 xref=smm/spacing_ref
/// AdjustToRefGrid: x'(x)*spacing/spacing_ref = xref
/// 1. setdimension, setimage before reading transform from file in order to keep consistency
/// 2. transform type:AFF,FFD,FFDMI,FLD
///
/// 2010-sep-10, xiahai add: GetTransformationDerivativeByParameter
///
/// \ingroup zxhTransforms
///

class zxhTransformBase //: public zxhTransformGradientBase
{
public:
	/// construction
	zxhTransformBase();
	/// deconstruct
	virtual ~zxhTransformBase();

	///  will set dimension of transformation same to image if the two image's dimension same
	virtual bool	SetImage(zxhImageData*pTest, zxhImageData*pRef);

	///\return image
	virtual zxhImageData* GetTestImage()const						{ return m_pImageTest; };
	/// \return image
	virtual zxhImageData* GetRefImage()const						{ return m_pImageRef; };

	/// dimension
	virtual int	GetDimension()	const						{ return m_iDimension; };

	/// 专门为把Transform当作导数类型用在梯度下降算法中
	virtual bool	SetDimension(int i)						{ m_iDimension = i; return true; };

	/// return image type string in upper case, currently have, used for extension of transform files
	virtual std::string GetTransformType()const				{ return "TFM"; };//AFF,FFD,FFDMI,FLD

	///
	virtual bool SameTransformType(const zxhTransformBase*p) const
	{
		return (strcmp(this->GetTransformType().c_str(), p->GetTransformType().c_str()) == 0);
	}

	/// if pRet==0, new one object in method of clone
	virtual zxhTransformBase*	CloneTo(zxhTransformBase*&pRet) const;

	// abtract function

	/// return whether is outlier point given reference image
	/// fVectorFrom is the grid coordinate of test image and
	/// fVectorTo is a return value of grid coordinate of reference image
	/// if the test/ref image is not set then the coordinate is regarded as world coordinate
	/// hence they can be seen as with image spacing 1mm
	virtual bool	TransformPointTo(const float fVectorFrom[ZXH_ImageDimensionMax], float fVectorTo[ZXH_ImageDimensionMax]) const
	{
		this->TransformPointToWorld(fVectorFrom, fVectorTo);
		this->AdjustWorldToRefImage(fVectorTo);
		if (m_pImageRef) return !m_pImageRef->InsideImage(fVectorTo);
		return false;
	};
	/// \return the world coordinate after transform
	virtual void	TransformPointToWorld(const float fVectorFrom[ZXH_ImageDimensionMax], float fVectorToWorld[ZXH_ImageDimensionMax]) const
	{
		float fv[] = { fVectorFrom[0], fVectorFrom[1], fVectorFrom[2], fVectorFrom[3] };
		this->AdjustTestImageToWorld(fv);
		this->TransformWorldToWorld(fv, fVectorToWorld);
	};
	///
	virtual void	TransformWorldToWorld(const float fVectorFromWorld[ZXH_ImageDimensionMax], float fVectorToWorld[ZXH_ImageDimensionMax])const
	{
		return;
	};

	/// inverse transform, input fVectorToWorld, output fVectorFromWorld
	virtual void	InverseTransformWorld(const float fVectorToWorld[ZXH_ImageDimensionMax], float fVectorFromWorld[ZXH_ImageDimensionMax]) const
	{
		std::cerr << "error: inverse transform world has not be implemented.\n";
	};

	/// \return xJacobian[0][0..3] dX/dx.dy.dz.dt
	//virtual bool	GetJacobianMatrix(float fVector[ZXH_ImageDimensionMax],float* xJacobian[ZXH_ImageDimensionMax]){return false;};
	/// pJacobian[i*4+0..3] = dX/dx.dy.dz.dt
	virtual bool	GetJacobianMatrix(const float fVector[ZXH_ImageDimensionMax], float pJacobian[ZXH_ImageDimensionMax*ZXH_ImageDimensionMax]) const;
	///  \return if dimension==3, 
	/// xJacobian = [ dY1/dx1 dY1/dx2 dY1/dx3 0 
	///               dY2/dx1 dY2/dx2 dY2/dx3 0
	///               dY2/dx1 dY2/dx2 dY2/dx3 0
	///                  0       0       0    1 ]
	virtual bool	GetJacobianMatrixWorld(const float fWorld[ZXH_ImageDimensionMax], float pJacobian[ZXH_ImageDimensionMax*ZXH_ImageDimensionMax])const;

	///
	virtual void	TransformIntensityGradientVectorTo(float fVectorFrom[ZXH_ImageDimensionMax], float fVectorTo[ZXH_ImageDimensionMax], float* fWorldFrom) const;

	/// inverse dC/dx1=dC/dy1*dy1/dx1 + dC/dy2*dy2/dx1 + dC/dy3*dy3/dx1
	virtual void	InverseTransformIntensityGradientVector(float fVectorRef[ZXH_ImageDimensionMax], float fInverseVectorTest[ZXH_ImageDimensionMax], float* fWorldTest) const;


	virtual	bool	SetTransformPara(float);
	/// similar to SetTransformPara, unless calling SetParameterValueWithoutUpdate, which can be different in affine without update matrix
	virtual	bool	SetTransformParaWithoutUpdate(float);
	/// operators
	/// A+B
	virtual bool	GetAdd(const zxhTransformBase*, const zxhTransformBase*);
	/// +
	virtual bool	Add(const zxhTransformBase*);
	/// A-B
	virtual bool	GetSubtract(const zxhTransformBase*, const zxhTransformBase*);
	/// -
	virtual bool	Subtract(const zxhTransformBase*);
	/// A*f
	virtual bool	GetMultiply(const zxhTransformBase*, float);
	/// *
	virtual bool	Multiply(float);
	/// A*B
	virtual bool	GetMultiplyByPara(const zxhTransformBase*, const zxhTransformBase*);
	/// *B
	virtual bool	MultiplyByPara(const zxhTransformBase*);

	/// transformation regarded as a vector
	virtual float	GetMagnitudeAsVector() const;


	/// for optimization step length normalization
	/* *** validate results on FFDs:
	/* 1) [X] using assume each control point has own support local regions to avoid turbulent changes--> reduce accuracy
	/* 2) [X] using a global maximal mag from all CPs for normalization, a big turbulent changes dominant the optimization length
	*    --> reduce robustness in registration where large deformation fields are needed
	* ****************************************/
	virtual float	GetMaxLocalMagnitudeAndNormalisation(bool bNormalisationByLocalMagnitude)
	{
		float fmag = GetMagnitudeAsVector();
		if (fmag > 0) this->Multiply(1 / fmag);
		return fmag;
	};
	/*virtual void	SetNormalisationByLocalMagnitude( bool b )
	{	m_bNormalisationByLocalMagnitude	= b ; }
	virtual bool	GetNormalisationByLocalMagnitude()
	{	return m_bNormalisationByLocalMagnitude ; }*/
	;
	///
	virtual float	GetMaxAbsValueFromParameters() const;

	/// \return
	virtual float	GetPointMutiplyAsVector(const zxhTransformBase*) const;

	/// set transformation as a identity transform
	virtual bool	SetTransformIdentity() { return false; };
	/// read transformation from stream, donot care the comments
	virtual bool	SetTransformFromStream(std::ifstream & ifs);
	/// from file
	virtual bool	SetTransformFromFile(const char* pFile);
	///
	virtual std::string GetTransformFileName() 	const	{ return m_strFileName; };

	/// description of the transformation type(not gradient), can be saved to disk as a file
	/// and read by SetTransformFromFile
	virtual	std::string GetPrintString() const;

	/// save to disc, if sFileNameWithoutExt is with Ext and Ext same as TransformType, then just save if
	virtual bool SaveTransform2Disc(std::string sFileNameWithoutExt) const;

	/// make the transformation one to one, currently only for FFD
	/// return the portion of degrees needs corrections
	virtual float Guarantee1To1(void)	{ return -1; };
	/// set
	virtual void SetGuarantee1To1Effected(bool b = true)				{ m_bEffect1to1 = b; };
	/// for registration optimization
	virtual bool	AdvanceStep(const zxhTransformBase*p, float f = 1)
	{
		zxhTransformBase *pf = 0; p->CloneTo(pf); pf->GetMultiply(p, f);
		bool b = this->Add(pf); delete pf; return b;
	};

	/// adjust the coordinate of the transformed point to the coordinate system
	/// of the reference image, considering the different spacing
	/// this is much more complicated than statement here, consider interpolation
	/// is based on image lattice points
	inline void AdjustTestImageToWorld(float afTest[]) const
	{
		if (m_pImageTest)
			m_pImageTest->ImageToWorld(afTest);
		else
		{
			std::cerr << "error: no set transform images to AdjustTestToworld\n"; exit(1);
		};
	}
	/// ajust the world coodinate in reference image's world coordinate to reference image grid
	inline void AdjustWorldToRefImage(float afRef[]) const
	{
		if (m_pImageRef)
			m_pImageRef->WorldToImage(afRef);
		else
		{
			std::cerr << "error: no set transform images to AdjustToRefGrid\n"; exit(1);
		};
	};

	/// Get partial differential of transformation [Y1 Y2 Y3]=T(x), with respect to the transformation parameter Ti with index
	/// return whether success, results (e.g. 3D, [d(Y1)/d(Ti) d(Y2)/d(Ti) d(Y3)/d(Ti)]) returned in g[4]
	//virtual bool	GetTransformationDerivativeByParameter(float g[4], int index, float Worldx[]);
	///
	virtual int		GetNoParameters()const
	{
		return-1;
	};
	///  
	virtual float	GetParameters(int index)const{ return 0; };
	/// set transformation parameter value,
	/// if the transformation model need pre-computation from parameters for efficiency such as affine,
	/// then re-implement this function to call the pre-computation method after set the value
	virtual bool	SetParameterValue(int index, float f)
	{
		return false;
	};
	/// default same as SetParameterValue, but for affine, here no update matrix operation will be taken
	virtual bool	SetParameterValueWithoutUpdate( int index, float f ) {return SetParameterValue(index,f);} ;
	/// get transforms DOF
	virtual int GetDegreeOfFreedom() const	{return GetNoParameters() ; };
	/// get parameter value >= f , dof
	virtual int GetDOFWithValue( float f ) const	
	{ 
		int idof=0 , no=GetNoParameters(); 
		for( int index=0; index<no; ++index )
			if( GetParameters(index) >= f ) idof++ ;
		return idof ;	
	};
	///
	virtual int GetActiveDof()const				{ return GetDOFWithValue(ZXH_FloatPrecision) ;};
	 
protected:
	/// constant value
	const static unsigned int static_buffer_size_of_base=1024;

	/// dimension of transformation, currently only same as image
	int				m_iDimension;

	/// state of whether guarantee the transformation as smooth as one to one
	bool			m_bEffect1to1;

	/// not transformation property
	zxhImageData*		m_pImageTest;
	/// ref image
	zxhImageData*		m_pImageRef;
	///
	std::string			m_strFileName ;
	
}  ;

namespace zxh{
/// \return transformation collection to transform from world to world coordinates
bool TransformWorld2WorldByCollection(const zxhTransformBase*const* pListTransform, int nTransform, // const pListTransform makes it difficult to use
		 int dimension,const float fWorldFrom[ZXH_ImageDimensionMax],float fWorldTo[ZXH_ImageDimensionMax]);

/// \return whether affine
bool IsAffTransform(zxhTransformBase *pTrans);
///
bool IsAffineMatrix(zxhTransformBase *pTrans);
/// \return whether FFD based transform
bool IsFFDTransform(zxhTransformBase*pTrans);
/// \return whether FFD based transform
bool IsFFDTransform(std::string type);
///
enum ZXH_TransformType{TransformTypeRigid=0, TransformTypeScale=1, TransformTypeAffine=2, TransformTypeFFD=3,
TransformTypeSkew=5,
TransformTypeLocalAffine=10};

}
#endif //
