
#include "zxhTransformBase.h"
#include <fstream>
#include <iostream>

zxhTransformBase::zxhTransformBase():
	m_bEffect1to1( true ),
	m_iDimension( 0 ),
	m_pImageTest( 0 ),
	m_pImageRef( 0 )
//	,m_bNormalisationByLocalMagnitude(false)
{
	m_strFileName = "";
}

zxhTransformBase::~zxhTransformBase()
{

}

zxhTransformBase* zxhTransformBase::CloneTo(zxhTransformBase*&pRet) const
{
	if( pRet==0 ) return 0;
	pRet->m_pImageTest	= m_pImageTest;
	pRet->m_pImageRef	= m_pImageRef;
	pRet->m_iDimension	= m_iDimension;
	pRet->m_bEffect1to1	= m_bEffect1to1;
	pRet->m_strFileName = m_strFileName ;
//	pRet->m_bNormalisationByLocalMagnitude = m_bNormalisationByLocalMagnitude ;
	return pRet;
};

bool	zxhTransformBase::SetImage(zxhImageData*pTest,zxhImageData*pRef)
{
	if(pTest==0||pRef==0)
		return false;
	m_pImageTest	= pTest;
	m_pImageRef		= pRef;
	if(this->GetDimension()<1||this->GetDimension()>ZXH_ImageDimensionMax)//not set transform
	{
		if(pTest->GetDimension()==pRef->GetDimension())
		{	this->SetDimension( pTest->GetDimension() );}
		else
		{std::cerr<<"error: test image and ref image should have same dimension\n"; exit(1);}
	}
	return true;
};

std::string zxhTransformBase::GetPrintString() const
{
	std::string s;
	char buffer[static_buffer_size_of_base];
	memset((void*)buffer,0,static_buffer_size_of_base);
	sprintf(buffer,"#dimension#\n%i\n#one to one#\n%i\n",
		m_iDimension,int(m_bEffect1to1));
	s=buffer;
	return s;
};
/// save to disc, if sFileNameWithoutExt is with Ext and Ext same as TransformType, then just save if
bool	zxhTransformBase::SaveTransform2Disc( std::string sFileNameWithoutExt )const
{
	std::string strTransformType = this->GetTransformType() ;
	if( strcmp( zxh::GetExtension( sFileNameWithoutExt ).c_str(), strTransformType.c_str() ) != 0 )
		sFileNameWithoutExt += "." + strTransformType ;
	std::ofstream ofs ;
	ofs.open( sFileNameWithoutExt.c_str() ) ;
	if(  ofs.fail() )
	{
		std::cerr<< "error: save to file name "<< sFileNameWithoutExt <<" fail! " ;
		return false ;
	}
	std::string result = this->GetPrintString() ;
	ofs.write( result.c_str(), result.length() ) ;
	ofs.close();
	return true ;
}

bool	zxhTransformBase::SetTransformFromFile(const char* pFile)
{
	m_strFileName = pFile ;
	std::ifstream ifs;
	bool bsuc=true;
	ifs.open(pFile,std::ios_base::in);
	if(ifs.fail()==true)
	{
		std::cerr<<"error: fail in opening transformation file "<<pFile<<"!\n";
		return false;
	}
	ifs.seekg(0,std::ios_base::beg);
	bsuc=SetTransformFromStream(ifs);
	if( bsuc == false )
	{
		m_strFileName = "";
		std::cerr<<"error: failed in reading file stream of "<<pFile<<"\n";
	}
	ifs.close();
	return bsuc;
};

bool	zxhTransformBase::SetTransformFromStream(std::ifstream & ifs)
{
	bool bsuc=true;
	char buffer[static_buffer_size_of_base] ;
	memset((void*)buffer,0,static_buffer_size_of_base);
	//dimension
	int idim;
	bsuc = bsuc&zxh::ScanInteger(ifs,buffer,static_buffer_size_of_base,idim);
	SetDimension( idim ) ;
	//differential step, property of gradient object
	//bsuc = bsuc&zxh::ScanFloat(ifs,buffer,static_buffer_size_of_base,m_fDifferentialStep);
	//one to one
	int ionetoone;
	bsuc = bsuc&zxh::ScanInteger(ifs,buffer,static_buffer_size_of_base,ionetoone);
	m_bEffect1to1=static_cast<bool>(ionetoone);
	return bsuc;
}
bool	zxhTransformBase::GetJacobianMatrix(const float fVector[ZXH_ImageDimensionMax],float pJacobian[ZXH_ImageDimensionMax*ZXH_ImageDimensionMax]) const
{
	// pJacobian[i*4+0..3] = dX/dx.dy.dz.dt
	float World[] = {fVector[0],fVector[1],fVector[2],fVector[3]} ;
	this->AdjustTestImageToWorld( World ) ;
	return GetJacobianMatrixWorld( World, pJacobian );
}

/// Get partial differential of transformation [Y1 Y2 Y3]=T(x), with respect to the transformation parameter Ti with index
/// return whether success, results (e.g. 3D, [d(Y1)/d(Ti) d(Y2)/d(Ti) d(Y3)/d(Ti)]) returned in g[4]
/*bool	zxhTransformBase::GetTransformationDerivativeByParameter(float g[], int index, float Worldx[])
{
	float Ti = this->GetParameters(index) ;

	float Worldtof[] = { Worldx[0],Worldx[1],Worldx[2],Worldx[3] };
	float Worldtob[] = { Worldx[0],Worldx[1],Worldx[2],Worldx[3] };
	float step = 0.5 ;

	float currvalue = pTi[0] ;

	//forward
	this->SetParameterValue( index, currvalue+step ) ;
	this->TransformWorldToWorld( Worldx, Worldtof ) ;

	// backward
	this->SetParameterValue( index, currvalue-step ) ;
	this->TransformWorldToWorld( Worldx, Worldtob ) ;

	// return and result
	this->SetParameterValue( index, currvalue ) ;

	for( int jd=0; jd<4; ++jd )
		g[jd] = Worldtof[jd] -Worldtob[jd] ;
	return true ;
}*/
bool	zxhTransformBase::GetJacobianMatrixWorld(const float fWorld[ZXH_ImageDimensionMax],float pJacobian[ZXH_ImageDimensionMax*ZXH_ImageDimensionMax])const
{
	float Worldtof[] = { fWorld[0],fWorld[1],fWorld[2],fWorld[3] };
	float Worldtob[] = { fWorld[0],fWorld[1],fWorld[2],fWorld[3] };
	float step = 0.5 ;
	for( int ip=0; ip<ZXH_ImageDimensionMax*ZXH_ImageDimensionMax; ++ip )
			pJacobian[ip] = 0;

	pJacobian[0]    = pJacobian[ZXH_ImageDimensionMax+1] = pJacobian[2*ZXH_ImageDimensionMax+2] = pJacobian[3*ZXH_ImageDimensionMax+3] = 1 ;//identity

	float wco[] = { fWorld[0], fWorld[1], fWorld[2], fWorld[3] } ; 									
	for( int idx=0; idx<m_iDimension; ++idx )
	{
		// forward
		wco[idx] += step ;
		this->TransformWorldToWorld( wco, Worldtof ) ;

		// backward
		wco[idx] -= step*2 ;
		this->TransformWorldToWorld( wco, Worldtob ) ;

		// return and result
		wco[idx] += step ;
		for( int jd=0; jd<m_iDimension; ++jd )
				pJacobian[jd*ZXH_ImageDimensionMax+ idx] = Worldtof[jd]-Worldtob[jd] ; // $d(X'Y'Z')^T / d x $
	}
	return true ;
}

void	zxhTransformBase::TransformIntensityGradientVectorTo(float fVectorFrom[ZXH_ImageDimensionMax],float fVectorTo[ZXH_ImageDimensionMax],float* fWorldFrom ) const
{
	float forward[]={0,0,0,0}, backward[]={0,0,0,0}, transforward[]={0,0,0,0}, transbackward[]={0,0,0,0}, mag=0,magj=0,mag1=0; 
	fVectorTo[0]=fVectorTo[1]=fVectorTo[2]=fVectorTo[3]=0 ;
	mag=zxh::MagnitudeOfVector( fVectorFrom, m_iDimension ) ;
	if( mag==0 ) return ;

	for( int id=0 ; id<m_iDimension; ++id )
	{	
		forward[id]  = fWorldFrom[id] + 0.5*fVectorFrom[id]/mag ;
		backward[id] = fWorldFrom[id] - 0.5*fVectorFrom[id]/mag ;
	}
	this->TransformWorldToWorld( forward, transforward ) ;
	this->TransformWorldToWorld( backward,transbackward) ;

	for( int id=0; id<m_iDimension; ++id )
		fVectorTo[id] = transforward[id] - transbackward[id] ; 

		// one to compute is to use 
	mag1 = zxh::MagnitudeOfVector( fVectorTo, m_iDimension ) ;
	if( mag1!=0 )  
	for( int id=0; id<m_iDimension; ++id )
		fVectorTo[id] = fVectorTo[id]*mag/(mag1*mag1) ; // mag/mag1 * 1/mag1 
	
	return ;

		// the other using G(y) = Jac(inv T) G(x) and assuming |G_ref(y)| = |G_test(x)|/Jacob(Tx)|
	/*zxh::NormaliseVector( fVectorTo[id], m_iDimension ) ;
	
	float jacobian[]={0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 } ;
	GetJacobianMatrixWorld( fWorldTest, jacobian ) ;
	magj = zxh::DeterminentMatrix3D( jacobian, 4 ) ;
	if( magj!=0 ) 
	for( int id=0; id<m_iDimension; ++id )
		fVectorTo[id] = fVectorTo[id]*mag/magj ; */
}

	// inverse dC/dx1=dC/dy1*dy1/dx1 + dC/dy2*dy2/dx1 + dC/dy3*dy3/dx1
void	zxhTransformBase::InverseTransformIntensityGradientVector( float fVectorRef[ZXH_ImageDimensionMax],float fInverseVectorTest[ZXH_ImageDimensionMax],float* fWorldTest ) const
{
	float jacobian[]={0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 } ;
	GetJacobianMatrixWorld( fWorldTest, jacobian ) ;
	fInverseVectorTest[0]=fInverseVectorTest[1]=fInverseVectorTest[2]=fInverseVectorTest[3]=0;
	for( int id=0; id<m_iDimension; ++id )
		fInverseVectorTest[id] = fVectorRef[0]*jacobian[id] + fVectorRef[1]*jacobian[ZXH_ImageDimensionMax+id] + fVectorRef[2]*jacobian[ZXH_ImageDimensionMax*2+id] ; 
};
/* ************************** OPERATORS ************************************** */
bool	zxhTransformBase::GetAdd(const zxhTransformBase* p1,const zxhTransformBase* p2)
{
	int no = this->GetNoParameters();
	bool b=true ;
	for(int idx=0;idx<no;++idx)
	{
		b=b&this->SetParameterValue( idx, p1->GetParameters(idx) + p2->GetParameters(idx) );
	}
	return b;
};
bool	zxhTransformBase::Add(const zxhTransformBase* p1)
{
	int no = this->GetNoParameters();
	bool b=true ;
	for(int idx=0;idx<no;++idx)
	{
		b=b&this->SetParameterValue( idx, this->GetParameters(idx) + p1->GetParameters(idx) ) ;
	}
	return b;
};
bool	zxhTransformBase::GetSubtract(const zxhTransformBase*p1,const zxhTransformBase*p2)
{
	int no = this->GetNoParameters();
	bool b=true ;
	for(int idx=0;idx<no;++idx)
	{
		b=b&this->SetParameterValue( idx, p1->GetParameters(idx) - p2->GetParameters(idx) );
	}
	return b;
};
bool	zxhTransformBase::Subtract(const zxhTransformBase*p1)
{
	int no = this->GetNoParameters();
	bool b=true ;
	for(int idx=0;idx<no;++idx)
	{		
		b=b&this->SetParameterValue( idx, this->GetParameters(idx) - p1->GetParameters(idx) );
	}
	return b;
};
bool	zxhTransformBase::GetMultiply(const zxhTransformBase* p1,float f)
{
	int no = this->GetNoParameters();
	bool b=true ;
	for(int idx=0;idx<no;++idx)
	{
		b=b&this->SetParameterValue( idx, p1->GetParameters(idx) * f );
	}
	return b;
};
bool	zxhTransformBase::Multiply(float f)
{
	int l=this->GetNoParameters();
	for(int idx=0;idx<l;++idx)
		this->SetParameterValue( idx, GetParameters(idx)*f ) ;
	return true;
};
bool	zxhTransformBase::GetMultiplyByPara(const zxhTransformBase*p1,const zxhTransformBase*p2)
{
	int no = this->GetNoParameters();
	for(int idx=0;idx<no;++idx)
	{
		this->SetParameterValue( idx, p1->GetParameters(idx) * p2->GetParameters(idx) );
	}
	return true;
};
bool	zxhTransformBase::MultiplyByPara(const zxhTransformBase* p1)
{
	int no = this->GetNoParameters();
	for(int idx=0;idx<no;++idx)
	{
		this->SetParameterValue( idx, GetParameters(idx) * p1->GetParameters(idx) ) ;
	}
	return true;
}
bool	zxhTransformBase::SetTransformPara( float f )
{
	int no = this->GetNoParameters();
	for(int idx=0;idx<no;++idx)
	{
		this->SetParameterValue( idx, f ) ;
	}
	return true ;
}
bool	zxhTransformBase::SetTransformParaWithoutUpdate( float f )
{
	int no = this->GetNoParameters();
	for(int idx=0;idx<no;++idx)
	{
		this->SetParameterValueWithoutUpdate( idx, f ) ;
	}
	return true ;
}
float	zxhTransformBase::GetMagnitudeAsVector( ) const
{
	float fmag=0; 
	int no = this->GetNoParameters();
	for(int idx=0;idx<no;++idx)
	{
		fmag += this->GetParameters(idx)*this->GetParameters(idx);
	}
	return sqrt(fmag) ;
} 
float	zxhTransformBase::GetMaxAbsValueFromParameters( ) const
{
	float fmax=0; 
	int no = this->GetNoParameters();
	for(int idx=0;idx<no;++idx)
	{
		fmax = zxh::maxf( zxh::absf(this->GetParameters(idx)), fmax );
	}
	return fmax ;
}   
float	zxhTransformBase::GetPointMutiplyAsVector( const zxhTransformBase*p )const
{
	float fprod=0 ;
	int no = this->GetNoParameters();
	for(int idx=0;idx<no;++idx)
	{
		fprod += this->GetParameters(idx)*p->GetParameters(idx);
	}
	return fprod ;
}
 
/* ************************* OPERATORS end *********************************** */
namespace zxh{
/// \return transformation collection to transform from Worldical to Worldical coordinates
bool TransformWorld2WorldByCollection(const zxhTransformBase*const *pListTransform, int nTransform,// const pListTransform makes it difficult to use
									int dimension,const float fWorldFrom[ZXH_ImageDimensionMax],float fWorldTo[ZXH_ImageDimensionMax])
{ 
	for( int idim = 0 ; idim < dimension; ++idim )
		fWorldTo[idim] = fWorldFrom[idim] ;
	if(pListTransform==0||nTransform<1) return false;
	float fMedium[ZXH_ImageDimensionMax] ;
	for(int inum=0;inum<nTransform;++inum)
	{
		if( pListTransform[inum] == 0 )
			continue ;
		fMedium[0]=fWorldTo[0];fMedium[1]=fWorldTo[1];
		fMedium[2]=fWorldTo[2];fMedium[3]=fWorldTo[3];
		pListTransform[inum]->TransformWorldToWorld(fMedium, fWorldTo);
	}
	return true;
};

bool IsAffTransform(zxhTransformBase *pTrans)
{
	if( pTrans==0 ) return false ;
	std::string transformtype=pTrans->GetTransformType();
	zxh::case_upper(transformtype);
	if(strstr(transformtype.c_str(),"AFF")!=NULL)
		return true;
	else return false;
}
bool IsAffineMatrix(zxhTransformBase *pTrans)
{
	if( pTrans==0 ) return false ;
	std::string transformtype=pTrans->GetTransformType();
	zxh::case_upper(transformtype);
	if(strstr(transformtype.c_str(),"AFF")!=NULL||strstr(transformtype.c_str(),"MTX")!=NULL)
		return true;
	else return false;
}
bool IsFFDTransform(std::string type)
{
	zxh::case_upper(type);
	std::string transformtype=zxh::GetExtension(type);
	if( strcmp(type.c_str(),"FFD")==0 ||
		strcmp(type.c_str(),"GRID")==0 ||
		strcmp(type.c_str(),"FBS")==0 ||
		strcmp(transformtype.c_str(),"FFD")==0 ||
		strcmp(transformtype.c_str(),"GRID")==0 ||
		strcmp(transformtype.c_str(),"FBS")==0  )
		return true;
	return false;
}
bool IsFFDTransform(zxhTransformBase*pTrans)
{
	if( pTrans==0 ) return false ;
	std::string transformtype=pTrans->GetTransformType();
	return IsFFDTransform(transformtype);
}


}

