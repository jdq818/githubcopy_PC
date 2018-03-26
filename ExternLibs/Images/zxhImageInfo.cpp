
#ifndef zxhImageInfo_cpp
#define zxhImageInfo_cpp
#include "zxhImageInfo.h"

#ifdef HAVE_NIFTI
#include "nifti1_io.h"
#endif

zxhImageInfo::zxhImageInfo()
{
	FileName = "" ;			//unknown
	Dimension	= 0 ;		//un-set

	OrientationMethod = 0 ; //unknown
	ByteOfDataType = 0	;	//unknown
	DataType=0 ;			//unknown
	ImageFormat="";			//unknown

	for(int i=0;i<ZXH_ImageDimensionMax;i++)
	{
		Spacing[i]	= 1.0f;	//default
		Size[i]		= 1;	//default
		Origin[i]	= 0;	//unknown
	}
	Spacing[3] = 1.0f / 60.0f; //  unit second
	//method2 (for gipl)
	OrientationRotationMatrix[0][0] = 1 ;
	OrientationRotationMatrix[0][1] = 0 ;
	OrientationRotationMatrix[0][2] = 0 ;
	OrientationRotationMatrix[1][0] = 0 ;
	OrientationRotationMatrix[1][1] = 1 ;
	OrientationRotationMatrix[1][2] = 0 ;
	OrientationRotationMatrix[2][0] = 0 ;
	OrientationRotationMatrix[2][1] = 0 ;
	OrientationRotationMatrix[2][2] = 1 ;

	//method(-2) (for nii)
	QuaternFactor	= 1 ; // default using right-hand coordinate
	QuaternB		= 0 ;
	QuaternC 		= 0 ;
	QuaternD 		= 0 ;
	Qoffsetxyz[0] 	= 0 ;
	Qoffsetxyz[1] 	= 0 ;
	Qoffsetxyz[2] 	= 0 ;

	//method3 (can be used for all)
	for( int r=0; r<4; ++r )
		for( int c=0; c<4; ++c )
		{
			if( r==c )
				ImageToWorldMatrix[r][c] = 1 ;
			else
				ImageToWorldMatrix[r][c] = 0 ;
			WorldToImageMatrix[r][c] = ImageToWorldMatrix[r][c] ;
		}
	//rescale
	RescaleSlope = 1 ;
	RescaleIntercept = 0 ;
}


zxhImageInfo::~zxhImageInfo(){};

zxhImageInfo * zxhImageInfo::CloneTo(zxhImageInfo * & pRet) const
{
	if( pRet==0 ) pRet = new zxhImageInfo();
	pRet->FileName 		= this->FileName ;
	pRet->ImageFormat	= this->ImageFormat ;
	pRet->ByteOfDataType 	= this->ByteOfDataType ;
	pRet->DataType 		= this->DataType ;
	pRet->ImageFormat	= this->ImageFormat ;
	pRet->Dimension		= this->Dimension ;

	pRet->OrientationMethod = this->OrientationMethod ;

	for(int i=0;i<ZXH_ImageDimensionMax;i++)
	{
		pRet->Spacing[i]	= this->Spacing[i];
		pRet->Size[i]		= this->Size[i];
		pRet->Origin[i]	= this->Origin[i];
	}

	for( int i=0; i<3; ++i )
	{
		pRet->OrientationRotationMatrix[0][i] = this->OrientationRotationMatrix[0][i] ;
		pRet->OrientationRotationMatrix[1][i] = this->OrientationRotationMatrix[1][i] ;
		pRet->OrientationRotationMatrix[2][i] = this->OrientationRotationMatrix[2][i] ;
	}

	pRet->QuaternFactor	= this->QuaternFactor ;
	pRet->QuaternB		= this->QuaternB ;
	pRet->QuaternC 		= this->QuaternC ;
	pRet->QuaternD 		= this->QuaternD ;
	pRet->Qoffsetxyz[0] 	= this->Qoffsetxyz[0] ;
	pRet->Qoffsetxyz[1] 	= this->Qoffsetxyz[1] ;
	pRet->Qoffsetxyz[2] 	= this->Qoffsetxyz[2] ;

	for( int r=0; r<4; ++r )
		for( int c=0; c<4; ++c )
		{
			pRet->ImageToWorldMatrix[r][c] = this->ImageToWorldMatrix[r][c] ;
			pRet->WorldToImageMatrix[r][c] = this->WorldToImageMatrix[r][c] ;
		}

	pRet->RescaleIntercept = this->RescaleIntercept ;
	pRet->RescaleSlope = this->RescaleSlope ;
	return pRet ;
}
std::string zxhImageInfo::GetPrintString() const
{

	std::string ret="filename:" ;
	char buffer[1024];
	ret += FileName ;
	sprintf( buffer, "\nproblem of dimension: %d\n", Dimension ) ;
	ret += buffer ;
	sprintf( buffer, "data type:%d ; byte of type: %d ; slope/intercept: %f/%f\n", DataType, ByteOfDataType, RescaleSlope, RescaleIntercept ) ;
	ret += buffer ;
	float det = zxh::DeterminentMatrix3D( &ImageToWorldMatrix[0][0], 4 ) ;
	sprintf( buffer, "orientation method: %d; %s \n", OrientationMethod, det<0?"Left hand coordinate":"Right hand coordinate");
	ret += buffer ;
	sprintf( buffer, "spacing: (%.3f,%.3f,%.3f,%.3f); size: (%d,%d,%d,%d); origin: (%.3f,%.3f,%.3f,%.3f)\n",
			Spacing[0],Spacing[1],Spacing[2],Spacing[3],Size[0],Size[1],Size[2],Size[3],Origin[0],Origin[1],Origin[2],Origin[3] );
	ret += buffer ;

	//method2 (for gipl)
	sprintf( buffer, "orientation rotation matrix:\n(%f\t%f\t%f)\n(%f\t%f\t%f)\n(%f\t%f\t%f)\n",
			OrientationRotationMatrix[0][0], OrientationRotationMatrix[0][1], OrientationRotationMatrix[0][2],
			OrientationRotationMatrix[1][0], OrientationRotationMatrix[1][1], OrientationRotationMatrix[1][2],
			OrientationRotationMatrix[2][0], OrientationRotationMatrix[2][1], OrientationRotationMatrix[2][2] ) ;
	ret += buffer ;

	//method(-2) (for nii)
	sprintf( buffer, "quatern parameters (qf, qb, qc, qd, qoffset[3]): %f; %f,%f,%f; %f,%f,%f. \n",
			QuaternFactor, QuaternB, QuaternC, QuaternD, Qoffsetxyz[0], Qoffsetxyz[1], Qoffsetxyz[2] );
	ret += buffer ;

	//method3 (can be used for all)
	sprintf( buffer, "image grid to world coordinate matrix:\n(%f\t%f\t%f\t%f)\n(%f\t%f\t%f\t%f)\n(%f\t%f\t%f\t%f)\n",
			ImageToWorldMatrix[0][0], ImageToWorldMatrix[0][1], ImageToWorldMatrix[0][2], ImageToWorldMatrix[0][3],
			ImageToWorldMatrix[1][0], ImageToWorldMatrix[1][1], ImageToWorldMatrix[1][2], ImageToWorldMatrix[1][3],
			ImageToWorldMatrix[2][0], ImageToWorldMatrix[2][1], ImageToWorldMatrix[2][2], ImageToWorldMatrix[2][3] ) ;
	ret += buffer ;

	sprintf( buffer, "world coordinate to image grid matrix:\n(%f\t%f\t%f\t%f)\n(%f\t%f\t%f\t%f)\n(%f\t%f\t%f\t%f)\n",
			WorldToImageMatrix[0][0], WorldToImageMatrix[0][1], WorldToImageMatrix[0][2], WorldToImageMatrix[0][3],
			WorldToImageMatrix[1][0], WorldToImageMatrix[1][1], WorldToImageMatrix[1][2], WorldToImageMatrix[1][3],
			WorldToImageMatrix[2][0], WorldToImageMatrix[2][1], WorldToImageMatrix[2][2], WorldToImageMatrix[2][3] ) ;
	ret += buffer ;

	return ret ;
}


std::string	zxhImageInfo::GetPrintStringSimpleQuuatern() const
{
	char buffer[1024];
	sprintf(buffer,"#quaternion qf,qb,qc,qd,qx,qy,qz#\n%f\t %f\t%f\t%f\t %f\t%f\t%f\n",
		QuaternFactor, QuaternB, QuaternC, QuaternD, Qoffsetxyz[0], Qoffsetxyz[1], Qoffsetxyz[2] ); 
	std::string str = buffer;
	return str; 
}
bool zxhImageInfo::SetImageInfoSimpleQuaternFromStream( std::ifstream& ifs )  
{
 // else read in quaternion orientation info:
	 //      m_ImageInfo.QuaternFactor, m_ImageInfo.QuaternB, m_ImageInfo.QuaternC, m_ImageInfo.QuaternD, m_ImageInfo.Qoffsetxyz[0], m_ImageInfo.Qoffsetxyz[1], m_ImageInfo.Qoffsetxyz[2]
	char buffer[1024];
	float qp[7] = {0,0,0,0,0,0,0} ;
	int iqp = 0 ;
	std::string sLine,sContent,sComment ;
	while( iqp<7 && ifs.eof()==false&&ifs.fail()==false )
	{
		do{
			ifs.getline(buffer,1024);
			sLine=buffer;
			zxh::ParseStringLine(sContent,sComment,sLine);
			zxh::trim_both(sContent);
		}while( sContent.empty() && ifs.eof()==false&&ifs.fail()==false );
		std::istringstream strstrm(sContent.c_str());
		for(int i=iqp;i<7;++i)
		{	strstrm>>qp[i] ; ++iqp ; }
	}
	if( iqp == 7 )
	{
		QuaternFactor = qp[0] ; QuaternB = qp[1] ;
		QuaternC = qp[2] ; 		QuaternD = qp[3] ;
		for( int i=0; i<3; ++i )	
			Qoffsetxyz[i] = qp[i+4] ;
		UpdateOrientationInfo(-2);
	} 
	return iqp==7 ;
}

void zxhImageInfo::GetExtent( float e[] ) const
{
	for(int id=0; id<this->Dimension; ++id) e[id] = Spacing[id]*(Size[id]-1);
}
/// return the start and end of extent in world coordinates
void zxhImageInfo::GetExtent( float worldfrom[], float worldto[] ) const
{
	for( int i=0; i<ZXH_ImageDimensionMax; ++i )
	{
		worldfrom[i] = 0;
		worldto[i] = static_cast<float>(Size[i]-1 );
	}
	this->ImageToWorld( worldfrom ) ;
	this->ImageToWorld( worldto ) ;
	for( int i=0; i<ZXH_ImageDimensionMax; ++i )
		if( worldfrom[i] > worldto[i] )
		{
			float a =worldto[i] ;
			worldto[i] = worldfrom[i] ;
			worldfrom[i] = a ;
		}
};
float zxhImageInfo::GetVolumeOfPixel( ) const
{
	if( Dimension>=3 ) 
		return Spacing[0]*Spacing[1]*Spacing[2] ; 
	else return Spacing[0]*Spacing[1] ; 
}
bool zxhImageInfo::SameOrientationAs( const zxhImageInfo * pTest ) const
{
	if( pTest==0 ) return false ;
	if( this->Qoffsetxyz[0] != pTest->Qoffsetxyz[0] ||
 		this->Qoffsetxyz[1] != pTest->Qoffsetxyz[1] ||
 		this->Qoffsetxyz[2] != pTest->Qoffsetxyz[2] ||
		this->QuaternB != pTest->QuaternB ||
		this->QuaternC != pTest->QuaternC ||
		this->QuaternD != pTest->QuaternD ||
		this->QuaternFactor != pTest->QuaternFactor )
		return false ;
	return true ;
}
bool zxhImageInfo::SameDimSizeSpacingAs( const zxhImageInfo * pTest ) const
{
	if( pTest==0 ) return false ;
	if( Dimension != pTest->Dimension ) return false ;
	for(int i=0;i<Dimension;i++)
	{
		if( Spacing[i] != pTest->Spacing[i] ) return false ;
		if( Size[i] != pTest->Size[i] ) return false ;
	}
	return true ;
}
void zxhImageInfo::GetSizeUsingExtent(const float e[ZXH_ImageDimensionMax], int s[ZXH_ImageDimensionMax]) const
{
	for(int idim=0;idim<Dimension;++idim)
		s[idim] = int(ceil( e[idim]/Spacing[idim]+1 ) );
};
/// get new size if change to different spacing
void zxhImageInfo::GetSizeUsingSpacing(const float sp[ZXH_ImageDimensionMax], int sz[ZXH_ImageDimensionMax]) const
{
	for(int idim=0;idim<Dimension;++idim)
		sz[idim] = int(ceil( float(Size[idim]-1)*Spacing[idim]/sp[idim]+1 ) );
};


///
void zxhImageInfo::TemporalImageToWorld( float&f ) const 
{	f = f*Spacing[3] +Origin[3] ; }; 
///
void zxhImageInfo::TemporalWorldToImage( float&f ) const 
{	f = (f-Origin[3])/Spacing[3] ; }; 
void zxhImageInfo::WorldToImage( float fv[] ) const
{
	const float *p = &WorldToImageMatrix[0][0] ;
	float x = p[0]*fv[0] + p[1]*fv[1] + p[2]*fv[2] + p[3] ;
	float y = p[4]*fv[0] + p[5]*fv[1] + p[6]*fv[2] + p[7] ;
	float z = p[8]*fv[0] + p[9]*fv[1] + p[10]*fv[2] + p[11] ;
	fv[0] = x ; fv[1] = y ; fv[2] = z ;
	if(Dimension>3&&Spacing[3]!=0) 
		fv[3]=(fv[3]-Origin[3])/Spacing[3];
} ;
///
void zxhImageInfo::ImageToWorld( float fv[] ) const
{
	const float *p = &ImageToWorldMatrix[0][0] ;
	float x = p[0]*fv[0] + p[1]*fv[1] + p[2]*fv[2] + p[3] ;
	float y = p[4]*fv[0] + p[5]*fv[1] + p[6]*fv[2] + p[7] ;
	float z = p[8]*fv[0] + p[9]*fv[1] + p[10]*fv[2] + p[11] ;
	fv[0] = x ; fv[1] = y ; fv[2] = z ; 
	if(Dimension>3&&Spacing[3]!=0) 
		fv[3] = fv[3]*Spacing[3] +Origin[3] ;
} ;
///
void zxhImageInfo::ImageToPhysical( float fv[] ) const
{	
	for(int id=0;id<Dimension;++id) 
		fv[id] = fv[id] *Spacing[id]; 
}; 
///
void zxhImageInfo::PhysicalToImage( float fv[] ) const
{	
	for(int id=0;id<Dimension;++id) 
		fv[id] = fv[id] /Spacing[id]; 
}; 
///
void zxhImageInfo::ProjectPhysVectorToWorldCoordinate( float &fvx, float &fvy, float &fvz ) const
{
	float fv[] = {fvx,fvy,fvz} ;
	fvx = OrientationRotationMatrix[0][0]*fv[0] + OrientationRotationMatrix[0][1]*fv[1] + OrientationRotationMatrix[0][2]*fv[2] ;
	fvy = OrientationRotationMatrix[1][0]*fv[0] + OrientationRotationMatrix[1][1]*fv[1] + OrientationRotationMatrix[1][2]*fv[2] ;
	fvz = OrientationRotationMatrix[2][0]*fv[0] + OrientationRotationMatrix[2][1]*fv[1] + OrientationRotationMatrix[2][2]*fv[2] ;
} ;
///
void zxhImageInfo::ProjectWorldVectorToPhysCoordinate( float &fvx, float &fvy, float &fvz ) const
{
	float fv[] = {fvx,fvy,fvz} ;
	fvx = (WorldToImageMatrix[0][0]*fv[0] + WorldToImageMatrix[0][1]*fv[1] + WorldToImageMatrix[0][2]*fv[2])*Spacing[0] ;
	fvy = (WorldToImageMatrix[1][0]*fv[0] + WorldToImageMatrix[1][1]*fv[1] + WorldToImageMatrix[1][2]*fv[2])*Spacing[1] ;
	fvz = (WorldToImageMatrix[2][0]*fv[0] + WorldToImageMatrix[2][1]*fv[1] + WorldToImageMatrix[2][2]*fv[2])*Spacing[2] ;
} ;

bool zxhImageInfo::Project3DWorldVectorTo2DPlane( float *fc ) const
{
	// voxel[0,0,0] -> world is ImageToWorldMatrix[0-2][3]
	float fs[] = { fc[0]+ImageToWorldMatrix[0][3],
				   fc[1]+ImageToWorldMatrix[1][3],
				   fc[2]+ImageToWorldMatrix[2][3], 0 } ;
	WorldToImage( fs ) ;
	fs[2] = fs[3] = 0 ;
	ImageToWorld( fs ) ;
	fc[0] = fs[0] - ImageToWorldMatrix[0][3] ;
	fc[1] = fs[1] - ImageToWorldMatrix[1][3] ;
	fc[2] = fs[2] - ImageToWorldMatrix[2][3] ;

	return true ;

	/* tested same result ^_^ 2011-04-08
	std::cout<<"original ("<<fc[0]<<","<<fc[1]<<","<<fc[2]<<") changed to -- " ;
	std::cout<<"("<<fc[0]<<","<<fc[1]<<","<<fc[2]<<")" ;

	// test
	{
	///   gradient fc, project to 2D slice pImageSlice
	///   Origin_W   --->  {Origin_W  +  fc_W}
	///               |(W2I)
	///              \|/
	///   Origin_I   --->  {Origin_I  +  (_Ix,_Iy,_Iz)=(fc_I-w2i[0-2][3]) }    {noted:Origin_I=(0,0,0)}
	///               |(I2W){Origin_I +  | (_Ix, _Iy, 0)=proj }
	///              \|/                \|/
	///   Origin_W   --->  {Origin_W  +  (proj_W - i2w[0-2][3]) }
	/// SO.... fc_W_on2dplane =
	// o_Image={0,0,0,0}, o_world={ImageToWorldMatrix[0-2][3]};
	float fs[] = {fc[0], fc[1], fc[2], 0} ;
	WorldToImage( fs ) ;
	fs[0] = fs[0] - WorldToImageMatrix[0][3] ;
	fs[1] = fs[1] - WorldToImageMatrix[1][3] ;
	fs[2] = fs[3] = 0 ;
	ImageToWorld( fs ) ;
	fc[0] = fs[0] - ImageToWorldMatrix[0][3] ;
	fc[1] = fs[1] - ImageToWorldMatrix[1][3] ;
	fc[2] = fs[2] - ImageToWorldMatrix[2][3] ;


	std::cout<<" = ("<<fc[0]<<","<<fc[1]<<","<<fc[2]<<")\n" ;
	} */

}


bool	zxhImageInfo::InsideImageWithSliceThickness( const float pvox[] ) const 
{
	for( int id=0; id<Dimension; ++id )
		if( pvox[id] < -0.5 || pvox[id] > Size[id]-0.5 ) 
			return false ;
	return true ;
}



/// remove orientation or set orientation info to identity
void	zxhImageInfo::RemoveOrientationInfo()
{
	QuaternFactor	= 1 ; // default using right-hand coordinate
	QuaternB		= 0 ;
	QuaternC 		= 0 ;
	QuaternD 		= 0 ;
	Qoffsetxyz[0] 	= 0 ;
	Qoffsetxyz[1] 	= 0 ;
	Qoffsetxyz[2] 	= 0 ;
	UpdateOrientationInfo(-2) ;
	this->OrientationMethod = 0 ;
}
bool	zxhImageInfo::CopyIntensityRescaleInfoFrom( const zxhImageInfo*pSource)
{
	if( pSource == 0 )
	{
		std::cerr<<"error: image info source is null in \n" ;
		return false ;
	}
	RescaleSlope = pSource->RescaleSlope ;
	RescaleIntercept = pSource->RescaleIntercept ;
	return true ;
}

/// set orientation correspondingly using source info without guarantee consistent, need to call update afterwards
bool	zxhImageInfo::CopyOrientationInfoFrom(const zxhImageInfo*pSource)
{
	if( pSource == 0 )
	{
		std::cerr<<"error: image info source is null in \n" ;
		return false ;
	}
	if( OrientationMethod==0 ) //un-set
		OrientationMethod = pSource->OrientationMethod ;

	if( Dimension<1||Dimension>4 )
		Dimension = pSource->Dimension ;
	for(int i=0;i<ZXH_ImageDimensionMax;i++)
	{
		//Spacing[i]	= pSource->Spacing[i];
		//S/ize[i]		= pSource->Size[i];
		Origin[i]	= pSource->Origin[i];
	}

	for( int i=0; i<3; ++i )
	{
		OrientationRotationMatrix[0][i] = pSource->OrientationRotationMatrix[0][i] ;
		OrientationRotationMatrix[1][i] = pSource->OrientationRotationMatrix[1][i] ;
		OrientationRotationMatrix[2][i] = pSource->OrientationRotationMatrix[2][i] ;
	}

	QuaternFactor	= pSource->QuaternFactor ;
	QuaternB		= pSource->QuaternB ;
	QuaternC 		= pSource->QuaternC ;
	QuaternD 		= pSource->QuaternD ;
	Qoffsetxyz[0] 	= pSource->Qoffsetxyz[0] ;
	Qoffsetxyz[1] 	= pSource->Qoffsetxyz[1] ;
	Qoffsetxyz[2] 	= pSource->Qoffsetxyz[2] ;

	for( int r=0; r<4; ++r )
		for( int c=0; c<4; ++c )
		{
			ImageToWorldMatrix[r][c] = pSource->ImageToWorldMatrix[r][c] ;
			WorldToImageMatrix[r][c] = pSource->WorldToImageMatrix[r][c] ;
		}
	//  UpdateOrientationInfo(-2) ; // assuming all orientation info has been updated to be consistent
	return true ;
}

bool zxhImageInfo::UpdateImageInfoUsingNewSpacing( const float * spacing ) 
{
	if( spacing == 0 ) return false ;
	int size[] = {1,1,1,1} ; 
	GetSizeUsingSpacing( spacing, size ) ; 
	for( int i=0; i<ZXH_ImageDimensionMax; ++i ) 
	{
		Spacing[i] = spacing[i] ; 
		Size[i] = size[i] ; 
	}
	UpdateOrientationInfo(-2) ; 
	return true ;
}
bool zxhImageInfo::UpdateImageInfoByExtendRoi(  const int *addRoiFrom, const int *addRoiTo ) 
{    
	int newsize[]={1,1,1,1}; 
	for(int idim=0;idim<4;++idim)
	{
		newsize[idim] = Size[idim] + addRoiFrom[idim] + addRoiTo[idim] ;
		if( newsize[idim] < 1 )
		{
			std::cerr<<"warning: size of dimension "<<idim<<" will be less than 1 if extend ROI, hence not action be applied\n";
			return false ;
		}
	}
	int idimension = 4 ;
	if( newsize[3]==1 ) idimension = 3 ;
	if( newsize[3]==1 && newsize[2]==1 ) idimension = 2 ; 
	float leftcorner[] = { -addRoiFrom[0], -addRoiFrom[1], -addRoiFrom[2], -addRoiFrom[3] } ; 
	this->ImageToWorld( leftcorner ) ;

	for( int id=0; id<3; ++id )
		this->ImageToWorldMatrix[ id ][3] = leftcorner[id] ;
	this->UpdateOrientationInfo( 3 ) ; 
	this->Dimension = idimension ; 
	for(int idim=0;idim<4;++idim)
		this->Size[idim] = newsize[idim] ;
	
	return true ;
}

/// update orientation method 1,2,-2,3 using given method
///   always assume spacing set, method1 will update all the others, while the others do not update method1
///   set OrientationMethod at the end
void	zxhImageInfo::UpdateOrientationInfo(int iOrientatationMethod)
{
	switch (iOrientatationMethod) {
		case 1:
		{
			zxh::MatrixIdentity( &ImageToWorldMatrix[0][0], 4 ) ;
			ImageToWorldMatrix[0][0] = Spacing[0] ;
			ImageToWorldMatrix[1][1] = Spacing[1] ;
			ImageToWorldMatrix[2][2] = Spacing[2] ;
			UpdateOrientationInfo( 3 ) ;
		}
			break;
		case 2:
		{//ImageToWorldMatrix = [identity origin] x [ORM 0] x [scaling matrix] x [identity -(image_size-1)/2]
			// update method3
			float ORM[4][4] ;
			zxh::MatrixIdentity( &ORM[0][0], 4 ) ;
			for( int i=0; i<3; ++i )
			for( int j=0; j<3; ++j )
				ORM[i][j] = OrientationRotationMatrix[i][j] ;

		   	float idenORI[4][4];
			zxh::MatrixIdentity( &idenORI[0][0], 4 ) ;
 			idenORI[0][3] = Origin[0] ;
 			idenORI[1][3] = Origin[1] ;
 			idenORI[2][3] = Origin[2] ;

			float idenSize[4][4] ;
			zxh::MatrixIdentity( &idenSize[0][0], 4 ) ;
			idenSize[0][3] = - float(Size[0] - 1) / 2.0;
			idenSize[1][3] = - float(Size[1] - 1) / 2.0;
			idenSize[2][3] = - float(Size[2] - 1) / 2.0;

			float scalingM[4][4] ;
			zxh::MatrixIdentity( &scalingM[0][0], 4 ) ;
			scalingM[0][0] = Spacing[0];
			scalingM[1][1] = Spacing[1];
			scalingM[2][2] = Spacing[2];

			float temp1[4][4], temp2[4][4] ;
			zxh::MultiplyMatrix( &temp1[0][0], &idenORI[0][0], &ORM[0][0], 4 ) ;
			zxh::MultiplyMatrix( &temp2[0][0], &scalingM[0][0], &idenSize[0][0], 4 ) ;
			zxh::MultiplyMatrix( &ImageToWorldMatrix[0][0], &temp1[0][0], &temp2[0][0], 4 ) ;

			mat44 R ;
			for( int i=0;i<4; ++i ) for(int j=0;j<4;++j) R.m[i][j]=ImageToWorldMatrix[i][j] ;

				// update WorldToImageMatrix
			mat44 invR = nifti_mat44_inverse( R ) ;
			for( int i=0;i<4; ++i ) for(int j=0;j<4;++j) WorldToImageMatrix[i][j] = invR.m[i][j] ;

			// update method(-2): quatern
			float spx,spy,spz ;
			nifti_mat44_to_quatern( R,
					 &QuaternB, &QuaternC, &QuaternD, &Qoffsetxyz[0], &Qoffsetxyz[1], &Qoffsetxyz[2],
					 &spx,&spy,&spz, &QuaternFactor );

		}
			break;
		case -2: //nifti_quatern_to_mat44 and nifti_mat44_to_quatern)
		{
			// update method3
			mat44 R = nifti_quatern_to_mat44( QuaternB, QuaternC, QuaternD, Qoffsetxyz[0], Qoffsetxyz[1], Qoffsetxyz[2],
					Spacing[0], Spacing[1], Spacing[2], QuaternFactor ) ;
			zxh::SetMatrix( &ImageToWorldMatrix[0][0], &R.m[0][0], 4 ) ;
				// update WorldToImageMatrix
			mat44 invR = nifti_mat44_inverse( R ) ;
			for( int i=0;i<4; ++i ) for(int j=0;j<4;++j) WorldToImageMatrix[i][j] = invR.m[i][j] ;

			// update method2
			for (int i = 0; i < 3; i++)
			{
			    OrientationRotationMatrix[i][0] = ImageToWorldMatrix[i][0]/Spacing[0] ; // x-axis
			    OrientationRotationMatrix[i][1] = ImageToWorldMatrix[i][1]/Spacing[1] ; // y-axis
			    OrientationRotationMatrix[i][2] = ImageToWorldMatrix[i][2]/Spacing[2] ; // z-axis
			}
				// set origin
			mat44 identcen, invscaleM, invORM, temp1, temp2, temp3 ;
			zxh::MatrixIdentity( &identcen.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &invscaleM.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &invORM.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temp1.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temp2.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temp3.m[0][0], 4 ) ;
			for( int r=0; r<3; r++ )
			{
				identcen.m[r][3] = float(Size[r]-1)/2.0 ;
				invscaleM.m[r][r]= 1/Spacing[r] ;
			}
			for( int r=0; r<3; r++ )
				for( int c=0; c<3; ++c )
					invORM.m[r][c] = OrientationRotationMatrix[r][c] ;
			invORM = nifti_mat44_inverse( invORM ) ;
			zxh::MultiplyMatrix( &temp1.m[0][0], &ImageToWorldMatrix[0][0], &identcen.m[0][0], 4 ) ;
			zxh::MultiplyMatrix( &temp2.m[0][0], &temp1.m[0][0], &invscaleM.m[0][0], 4 ) ;
			zxh::MultiplyMatrix( &temp3.m[0][0], &temp2.m[0][0], &invORM.m[0][0], 4 ) ;
			for( int r=0; r<3; r++ )
			 	Origin[r] = temp3.m[r][3] ;

		}
			break;
		case 3: // only use ImageToWorldMatrix to set method2 and method(-2), assuming spacing is correctly set
		{
			// method2: origin and OrientationRotationMatrix
				//  ImageToWorldMatrix = [identity origin] x [ORM 0] x [scaling matrix] x [identity -(image_size-1)/2] =>
				//  [identity origin] = ImageToWorldMatrix x [identity (image_size-1)/2] x [scaling matrix]^{-1} x [ORM 0]^{-1}
			/*zxh::mat44 wtoiMatrix, itowMatrix, temMatrixImageCenter, temICSMORM ;
			zxh::MatrixIdentity( &wtoiMatrix.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &itowMatrix.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temMatrixImageCenter.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temICSMORM.m[0][0], 4 ) ;

			for(int c = 0; c < 4; c++)
				for(int r = 0; r < 4; r++)
					itowMatrix.m[r][c] = ImageToWorldMatrix[r][c] ;
			for(int j=0; j<3; ++j )
				temMatrixImageCenter.m[j][3] = float(Size[j] - 1) / 2.0;
			zxh::MatrixOP_InvertMatrix4D3D( wtoiMatrix, itowMatrix ) ;
			zxh::MatrixOP_MultiplyMatrix( &temICSMORM.m[0][0], &temMatrixImageCenter.m[0][0], &wtoiMatrix.m[0][0], 4 ) ;
			zxh::MatrixOP_MultiplyMatrix( &temMatrixImageCenter.m[0][0], &itowMatrix.m[0][0], &temICSMORM.m[0][0], 4 ) ;

			for( int i=0; i<3; ++i )
				Origin[i]  = temMatrixImageCenter.m[i][3] + ImageToWorldMatrix[i][3]; */

				// set orientation matrix
			for (int i = 0; i < 3; i++)
			{
			    OrientationRotationMatrix[i][0] = ImageToWorldMatrix[i][0]/Spacing[0] ; // x-axis
			    OrientationRotationMatrix[i][1] = ImageToWorldMatrix[i][1]/Spacing[1] ; // y-axis
			    OrientationRotationMatrix[i][2] = ImageToWorldMatrix[i][2]/Spacing[2] ; // z-axis
			}
				// set origin
			mat44 identcen, invscaleM, invORM, temp1, temp2, temp3 ;
			zxh::MatrixIdentity( &identcen.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &invscaleM.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &invORM.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temp1.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temp2.m[0][0], 4 ) ;
			zxh::MatrixIdentity( &temp3.m[0][0], 4 ) ;
			for( int r=0; r<3; r++ )
			{
				identcen.m[r][3] = float(Size[r]-1)/2.0 ;
				invscaleM.m[r][r]= 1/Spacing[r] ;
			}
			for( int r=0; r<3; r++ )
				for( int c=0; c<3; ++c )
					invORM.m[r][c] = OrientationRotationMatrix[r][c] ;
			invORM = nifti_mat44_inverse( invORM ) ;
			zxh::MultiplyMatrix( &temp1.m[0][0], &ImageToWorldMatrix[0][0], &identcen.m[0][0], 4 ) ;
			zxh::MultiplyMatrix( &temp2.m[0][0], &temp1.m[0][0], &invscaleM.m[0][0], 4 ) ;
			zxh::MultiplyMatrix( &temp3.m[0][0], &temp2.m[0][0], &invORM.m[0][0], 4 ) ;
			for( int r=0; r<3; r++ )
			 	Origin[r] = temp3.m[r][3] ;

			// method(-2): quatern
			mat44 R ; float computex, computey, computez;
			for( int i=0;i<4; ++i ) for(int j=0;j<4;++j) R.m[i][j]=ImageToWorldMatrix[i][j] ;
			nifti_mat44_to_quatern( R,
                     &QuaternB, &QuaternC, &QuaternD, &Qoffsetxyz[0], &Qoffsetxyz[1], &Qoffsetxyz[2],
                     &computex, &computey, &computez, &QuaternFactor );

			// update WorldToImageMatrix
			mat44 invR = nifti_mat44_inverse( R ) ;
			for( int i=0;i<4; ++i ) for(int j=0;j<4;++j) WorldToImageMatrix[i][j] = invR.m[i][j] ;
		}
			break;
		default:
			std::cerr<<"error: undefined orientation method\n" ;
			break;
	}
	this->OrientationMethod = iOrientatationMethod ;
}



/// m[4][4]
void zxhImageInfo::GetMatrixImageToPhysical( float *m ) const 
{
	for( int i=0;i<16; ++i ) m[i] = 0 ; 
	for( int id=0; id<Dimension&&id<3; ++id ) 
		m[id*4+id] = Spacing[id] ; 
	m[15]=1;
}
/// m[4][4] 
void zxhImageInfo::GetMatrixPhysicalToWorld( float *m ) const  
{
	float minusMip[16] ;
	for( int i=0;i<16; ++i ) m[i] = minusMip[i] = 0 ; 
	m[15]=minusMip[15] = 1 ;
	
	// M(i->w)=M(p->w) * M(i->p)
	// M(p->w)=M(i->w) * M^{-1}(i->p)
	for( int id=0; id<Dimension&&id<3; ++id ) 
		minusMip[id*4+id] = 1/Spacing[id] ; 
	zxh::MultiplyMatrix( m, &ImageToWorldMatrix[0][0], minusMip, 4 ) ;  
}

 

int	zxhImageInfo::ConstructIndexOffsetCube( int radiusx, int radiusy, int radiusz, int subx, int suby, int subz, int * offset ) const
{ 
	if( offset==0 ) return -1 ; 
	int num=0;
	for( int sez=-radiusz; sez<= radiusz; ++sez ) // search volume
	for( int sey=-radiusy; sey<= radiusy; ++sey )
	for( int sex=-radiusx; sex<= radiusx; ++sex ) 
	{ 
		offset[num++] = this->GridToIndex( sex*subx,sey*suby,sez*subz ) ;  
	}
	return num ; 
}
/// return num of pixels in the offset array, which is a sphere with radius in mm unit 
int	zxhImageInfo::ConstructIndexOffsetSphere( float fPhysRadiusMM, int* offset ) const
{
	if( offset==0 ) return -1 ; 
	int radiusx = fPhysRadiusMM/Spacing[0]; 
	int radiusy = fPhysRadiusMM/Spacing[1]; 
	int radiusz = fPhysRadiusMM/Spacing[2];
	int num=0;
	for( int sez=-radiusz; sez<= radiusz; ++sez ) // search volume
	for( int sey=-radiusy; sey<= radiusy; ++sey )
	for( int sex=-radiusx; sex<= radiusx; ++sex ) 
	{ 
		if( (sex*Spacing[0])*(sex*Spacing[0])+
			(sey*Spacing[1])*(sey*Spacing[1])+
			(sez*Spacing[2])*(sez*Spacing[2]) <= fPhysRadiusMM*fPhysRadiusMM )
			offset[num++] = this->GridToIndex( sex,sey,sez ) ;  
	}
	return num ; 
} 
/// return num of pixels in the offset array, which is a 2D circle with radius in mm unit 
int	zxhImageInfo::ConstructIndexOffset2DCirc( float fPhysRadiusMM, int* offset ) const
{
	if( offset==0 ) return -1 ; 
	int radiusx = fPhysRadiusMM/Spacing[0]; 
	int radiusy = fPhysRadiusMM/Spacing[1];  
	int num=0; 
	for( int sey=-radiusy; sey<= radiusy; ++sey )
	for( int sex=-radiusx; sex<= radiusx; ++sex ) 
	{ 
		if( (sex*Spacing[0])*(sex*Spacing[0])+
			(sey*Spacing[1])*(sey*Spacing[1]) <= fPhysRadiusMM*fPhysRadiusMM )
			offset[num++] = this->GridToIndex( sex,sey,0 ) ;  
	}
	return num ; 
} ;

/// return num of pixels in the offset array, which is a sphere with radius in 1 pixel 
int	zxhImageInfo::ConstructIndexOffsetNeighbour( int* offset ) const
{
	if( offset==0 ) return -1 ; 
	int nei[7][3] = { 0,0,0,   0,0,1, 0,0,-1,   0,1,0, 0,-1,0,   1,0,0, -1,0,0 } ; 
	int num=0; 
	for( int i=0; i<7; ++i )
		offset[num++] = this->GridToIndex( nei[i][0], nei[i][1], nei[i][2] ) ; 
	return num ; 
} ;
	
/// return num of pixels in the offset array, which is a sphere with radius in 1 pixel 
int	zxhImageInfo::ConstructIndexOffsetNeighbour2D( int* offset ) const
{
	if( offset==0 ) return -1 ; 
	int nei[5][3] = { 0,0,0,   0,1,0, 0,-1,0,   1,0,0, -1,0,0 } ; 
	int num=0; 
	for( int i=0; i<5; ++i )
		offset[num++] = this->GridToIndex( nei[i][0], nei[i][1], nei[i][2] ) ; 
	return num ; 
} ;

void zxhImageInfo::GetOrthogonalOrientationWorldCoordinateCornersSpacingSize( float Corner0World[], float Corner1World[], float SpacingWorld[], int WorldImageSize[] ) const
{
	for( int id=0; id<Dimension; ++id )
	{
		float from[ZXH_ImageDimensionMax]={0}, to[ZXH_ImageDimensionMax]={0} ;
		to[id] = Size[id]-1 ; 
		this->ImageToWorld(from) ; 
		this->ImageToWorld(to);
		float spacingforworld = 0 ; 
		int jDimOfWorld = 0 ; 
		for( int jspc=0; jspc<Dimension; ++jspc )
		{
			float f = zxh::absf(from[jspc]-to[jspc]) ;
			if( f>spacingforworld ) 
			{
				spacingforworld = f ; 
				jDimOfWorld = jspc;
			}
		}
		SpacingWorld[jDimOfWorld] = Spacing[id] ;
		WorldImageSize[jDimOfWorld] = Size[id] ;
		if( from[id]>to[id] )
		{
			Corner1World[jDimOfWorld] = from[id] ; 
			Corner0World[jDimOfWorld] = to[id] ; 
		}
		else
		{
			Corner1World[jDimOfWorld] = to[id] ; 
			Corner0World[jDimOfWorld] = from[id] ; 
		} 
	}
}

/// recompute Spacing, Size and 
void zxhImageInfo::UpdateOrthogonalImageQuaternInfoUsingWorldInfo( const float* Corner0World, const float* Corner1World, const float* SpacingWorld, const int* WorldImageSize )
{
	for( int id=0; id<Dimension; ++id )
	{
		float	from[ZXH_ImageDimensionMax] = {Corner0World[0],Corner0World[1],Corner0World[2],Corner0World[3]}, 
				to[ZXH_ImageDimensionMax] = {Corner0World[0],Corner0World[1],Corner0World[2],Corner0World[3]} ;
		to[id] = Corner1World[id] ; 
		this->WorldToImage(from) ; 
		this->WorldToImage(to) ;
		float lengthofimage = 0 ; 
		int jDimOfImage = 0 ; 
		for( int jspc=0; jspc<Dimension; ++jspc )
		{
			float f = zxh::absf(from[jspc]-to[jspc]) ;
			if( f>lengthofimage ) 
			{
				lengthofimage = f ; 
				jDimOfImage = jspc;
			}
		}
		Spacing[jDimOfImage] = SpacingWorld[id] ;
		Size[jDimOfImage] = WorldImageSize[id] ;
	}
	Qoffsetxyz[0] = Qoffsetxyz[1] = Qoffsetxyz[2] = 0 ; 
	this->UpdateOrientationInfo( -2 ) ; 
	float	afCenterCorr[] = { 0.5*(Size[0]-1), 0.5*(Size[1]-1), 0.5*(Size[2]-1), 0},
		afCenterOrig[] = { 0.5*(Corner0World[0]+Corner1World[0]), 0.5*(Corner0World[1]+Corner1World[1]), 0.5*(Corner0World[2]+Corner1World[2]), 0 } ; 
	this->ImageToWorld( afCenterCorr ) ; 
	for( int id=0; id<Dimension&&id<3; ++id )
		Qoffsetxyz[id] = afCenterOrig[id] - afCenterCorr[id] ; 

}




#endif //zxhImageInfo_cpp


