
#include "zxh.h"
#include <math.h>
#include <iostream>

//#ifdef __cplusplus
//	extern "C"
//	{
//#endif
int glbVerboseOutput = 1 ;
float glbfmDeriMetric[16] ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} ;
namespace zxh
{

	
/// \return save string to filename
/* bool SaveString2File(std::string str2save,char*pfilename)
{
	std::ofstream ofs;
	ofs.open(pfilename);
	ofs.write(str2save.c_str(),str2save.length());
	ofs.close();
}*/
ZXH_DLL_EXPORT bool SaveString2File(std::string str2save,std::string filename)
{ 
	std::ofstream ofs;
	if(  ofs.fail() )
	{
		std::cerr<< "error: open file name "<< filename <<" fail! " ; 
		return false ;
	} 
	ofs.open(filename.c_str());
	ofs.write(str2save.c_str(),str2save.length());
	ofs.close();
	return true ;
}

ZXH_DLL_EXPORT bool SaveString2FileAppend(std::string str2save,std::string filename)
{ 
	std::ofstream ofs;
	if(  ofs.fail() )
	{
		std::cerr<< "error: open file name "<< filename <<" fail! " ; 
		return false ;
	} 
	ofs.open(filename.c_str(), std::ios_base::app);
	ofs<< str2save ;
	ofs.close(); 
	return true ;
}
/// \return modified string s with spaces trimmed from left
ZXH_DLL_EXPORT void trim_left(std::string & s)
{
	int pos(0);
	int size((int)s.size());
	for(;pos<size&&(s[pos]==' '||s[pos]=='\t'||s[pos]=='\n'||s[pos]=='\r'||s[pos]=='\v');++pos);
	s=s.substr(pos,size-pos);
}

/// \return modified string s with spaces trimmed from right
ZXH_DLL_EXPORT void trim_right(std::string& s)
{
	//int size(s.size());
	int pos((int)s.size()-1);
	for ( ; pos>=0&&(s[pos]==' '||s[pos]=='\t'||s[pos]=='\n'||s[pos]=='\r'||s[pos]=='\v');--pos );
	s=s.substr(0,pos+1);
}

/// \return modified string s with spaces trimmed from edges
ZXH_DLL_EXPORT void trim_both(std::string& s)
{
	trim_right(s);
	trim_left(s);
}

/// \return modified string to upper case
ZXH_DLL_EXPORT void case_upper(std::string& s)
{
	int pos(0);
	int size((int)s.size());
	for(;pos<size;++pos)
		if(s[pos]>='a'&&s[pos]<='z')s[pos]-=32;
}

/// \return modified string to lower case
ZXH_DLL_EXPORT void case_lower(std::string& s)
{
	int pos(0);
	int size((int)s.size());
	for(;pos<size;++pos)
		if(s[pos]>='A'&&s[pos]<='Z')s[pos]+=32;
} 
///
ZXH_DLL_EXPORT std::string str_get_lowercase(std::string s)
{
	case_lower(s) ; 
	return s;
}
///
ZXH_DLL_EXPORT std::string str_get_uppercase(std::string s)
{
	case_upper(s) ; 
	return s;
}

///
ZXH_DLL_EXPORT bool is_option( std::string s )
{
	trim_both( s );
	if( s.length()>1 && s.c_str()[0]=='-' &&
		(s.c_str()[1]>'9' || s.c_str()[1]<'0') && s.c_str()[1]!='.')
		return true ;
	return false ; 
}
/// not exactly digital
ZXH_DLL_EXPORT bool is_digital( std::string s )
{
	float f=atof(s.c_str() );
	if( f!=0.0 && f!=HUGE_VAL && f!=-HUGE_VAL) return true ;
	trim_both( s );
	if( s.length()>0 )
	{
		if(s[0] =='+' ) s[0]='0';
		if(s[0] =='-' ) s[0]='0';
		if(s.find_first_of('.')!=std::string::npos)
		{s[s.find_first_of('.')]='0';}
		for( int ipos=0;ipos<s.length();++ipos )
			if( s.c_str()[ipos]!='0' )
				return false ;
	}else return false ;
	return true ;
}
/// \return sContent sComents of the strLine
ZXH_DLL_EXPORT void ParseStringLine(std::string& sContent,std::string& sComment, std::string& strLine)
{
	// extract the subline of the # % ? // comments and others
	int iSize	= (int)strLine.length();
	sComment="";
	int iPos = (int)strLine.find("#",0);
	if(iPos==-1) iPos=(int)strLine.find("%",0);
	if(iPos==-1) iPos=(int)strLine.find("?",0);
	//if(iPos==-1) iPos=(int)strLine.find("//",0);
	if(iPos!=-1) sComment=strLine.substr(iPos,iSize-iPos);
	sContent = strLine;
	if(iPos!=-1) sContent=strLine.substr(0,iPos);
}

/// \return find a line with searching string
ZXH_DLL_EXPORT bool FindStringFromStream(char*buffer,int buffersize,std::ifstream&ifs,char* substring)
{
	std::string s;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,buffersize);
		s=buffer;
		if(int(s.find(substring,0))!=-1)
			return true;
	}
	return false;
}

/// \return find next line with sContent.empty()==false
ZXH_DLL_EXPORT bool NextContentLine(std::string& sContent,char*buffer,int buffersize,std::ifstream&ifs)
{
	std::string sLine,sComment;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,buffersize);
		sLine=buffer;
		ParseStringLine(sContent,sComment,sLine);
		trim_both(sContent);
		if(sContent.empty()==false) return true;
	}
	return true;
}

/// \return int from the stream
ZXH_DLL_EXPORT bool ScanInteger(std::ifstream & ifs,char*buffer,int ibuffersize,int& idata)
{
	std::string sComment,sContent,sLine;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,ibuffersize);
		sLine=buffer;
		ParseStringLine(sContent,sComment,sLine);
		trim_both(sContent);
		if(sContent.length()!=0)
		{
			std::istringstream istr;
			istr.str(sContent);
			istr>>idata;
			return true;
		};
	}
	return false;
}

/// \return float scan from the stream
ZXH_DLL_EXPORT bool ScanFloat(std::ifstream&ifs,char*buffer,int ibuffersize,float& fdata)
{
	std::string sComment,sContent,sLine;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,ibuffersize);
		sLine=buffer;
		ParseStringLine(sContent,sComment,sLine);
		trim_both(sContent);
		if(sContent.length()!=0)
		{
			std::istringstream istr;
			istr.str(sContent);
			istr>>fdata;
			return true;
		};
	}
	return false;
}

/// \return string scan from stream, stop before space ect.
ZXH_DLL_EXPORT bool ScanString(std::ifstream&ifs,char*buffer,int ibuffersize,char*pdata)
{
	std::string sComment,sContent,sLine;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,ibuffersize);
		sLine=buffer;
		ParseStringLine(sContent,sComment,sLine);
		trim_both(sContent);
		if(sContent.length()!=0)
		{
			std::istringstream istr;
			istr.str(sContent);
			istr>>pdata;
			return true;
		};
	}
	return false;
}

/// swap data
ZXH_DLL_EXPORT void Swap2Bytes(void* p)
{
    char one_byte;
    char* data = static_cast<char*>(p);
    one_byte = data[0]; data[0] = data[1]; data[1] = one_byte;
};
ZXH_DLL_EXPORT void Swap4Bytes(void* p)
{
    char one_byte;
    char* data = static_cast<char*>(p);
    one_byte = data[0]; data[0] = data[3]; data[3] = one_byte;
    one_byte = data[1]; data[1] = data[2]; data[2] = one_byte;
};
ZXH_DLL_EXPORT void Swap8Bytes(void* p)
{
    char one_byte;
    char* data = static_cast<char*>(p);
    one_byte = data[0]; data[0] = data[7]; data[7] = one_byte;
    one_byte = data[1]; data[1] = data[6]; data[6] = one_byte;
    one_byte = data[2]; data[2] = data[5]; data[5] = one_byte;
    one_byte = data[3]; data[3] = data[4]; data[4] = one_byte;
};
ZXH_DLL_EXPORT void Swap10Bytes(void* p)
{
    char one_byte;
    char* data = static_cast<char*>(p);
    one_byte = data[0]; data[0] = data[9]; data[9] = one_byte;
    one_byte = data[1]; data[1] = data[8]; data[8] = one_byte;
    one_byte = data[2]; data[2] = data[7]; data[7] = one_byte;
    one_byte = data[3]; data[3] = data[6]; data[6] = one_byte;
    one_byte = data[4]; data[4] = data[5]; data[5] = one_byte;
};
ZXH_DLL_EXPORT void SwapBytes(void* p,unsigned char nbyte)
{
	switch(nbyte)
	{
	case 2:Swap2Bytes(p);break;
	case 4:Swap4Bytes(p);break;
	case 8:Swap8Bytes(p);break;
	case 10:Swap10Bytes(p);break;
	default:
		return ;//throw "unreconised swap bypes" ;
	}
}
///

ZXH_DLL_EXPORT float VectorOP_Point2LineDistance( float pnt[3],float pL1[3], float pL2[3], float bWithinLineSegment ) 
{
	float * pShort = &pL1[0], * pLong = &pL2[0] ; 
	float ABs = VectorOP_Distance( pnt, pShort, 3 ) ; 
	float ABl = VectorOP_Distance( pnt, pLong, 3 ) ; 
	if( ABs > ABl )  // swap
	{
		pShort = &pL2[0] ; 
		pLong  = &pL1[0] ; 
		float t=ABl ; 
		ABl = ABs ; 
		ABs = t ; 
	}
	float vSA[] = {pnt[0]-pShort[0], pnt[1]-pShort[1], pnt[2]-pShort[2]};
	float vSL[] = {pLong[0]-pShort[0], pLong[1]-pShort[1], pLong[2]-pShort[2]};
	float consine = VectorOP_Cosine( vSA, vSL, 3 ) ; 
	float vDis = sqrt(1-consine*consine)*ABs ; 
	if( consine>=0 || bWithinLineSegment==false ) 
		return vDis ; 
	return ABl ; // withinline segment
}

///
ZXH_DLL_EXPORT bool FileExist(const std::string filename)
{
  std::ifstream ifile(filename.c_str());
  if( ifile.good() )
  {
	  ifile.close();
	  return true ;
  }
  return false;
}
/// end swap data
ZXH_DLL_EXPORT std::string GetFileNameNoExtension(const std::string filename)
{
	int ipos = int(filename.find_last_of('.')) ;
	if( ipos != -1 )
	{
		std::string noext=filename.substr(0,ipos);
		std::string ext=filename.substr(ipos+1,filename.length()-ipos-1);
		if( strcmp( ext.c_str(), "gz" ) == 0 )
		{
			ipos = int(noext.find_last_of('.')) ;
			if( ipos != -1 )
			{
				ext=noext.substr(ipos+1,noext.length()-ipos-1);
				if( strcmp( ext.c_str(), "nii" ) == 0 )
					noext=noext.substr(0,ipos);
			}
		}
		return noext;
	}
	else return filename;
} ;
ZXH_DLL_EXPORT std::string GetFileExtension(const std::string filename)
{ return GetExtension(filename);};
ZXH_DLL_EXPORT std::string GetExtension(const std::string filename)
{
	if(int(filename.find_last_of('.'))!=-1)
	{
		std::string ext=filename.substr(filename.find_last_of('.')+1,filename.length()-filename.find_last_of('.')-1);
		if( strcmp( ext.c_str(), "gz" ) == 0 )
		{
			std::string newfilename = filename.substr(0, filename.find_last_of('.'));
			int ipos = int(newfilename.find_last_of('.')) ; 
			if( ipos != -1 && 
				strcmp( newfilename.substr(ipos+1,newfilename.length()-ipos-1).c_str(),"nii") ==0 )
				ext  ="nii.gz";
		}
		return ext;
	}
	else return "";
} ;

ZXH_DLL_EXPORT std::string GetFileNameNoPath(const std::string filename)
{
	if(int(filename.find_last_of('/'))!=-1)
	{
		std::string file=filename.substr(filename.find_last_of('/')+1,filename.length()-filename.find_last_of('/')-1);
		return file;
	}
	else
		if(int(filename.find_last_of('\\'))!=-1)
		{
			std::string file=filename.substr(filename.find_last_of('\\')+1,filename.length()-filename.find_last_of('\\')-1);
			return file;
		}
		else return "";
} ;

ZXH_DLL_EXPORT std::string GetFileNamePath(const std::string filename)
{
	if(int(filename.find_last_of('/'))!=-1)
	{
		std::string path=filename.substr(0,filename.find_last_of('/'));
		return path;
	}
	else
		if(int(filename.find_last_of('\\'))!=-1)
		{
			std::string path=filename.substr(0,filename.find_last_of('\\'));
			return path;
		}
		else return "";
} ;
  
///\return abs value of float type
ZXH_DLL_EXPORT float absf( const float f )
{
	if( f >= 0 ) return f ;
	return -f ;
};

///\return round value of float type
ZXH_DLL_EXPORT int round( const float f )
{
	if( f>=0 ) return int(f+0.5) ;
	return -int(-f+0.5) ;
}
///\return max value of float type
ZXH_DLL_EXPORT float maxf( const float f1, const float f2 )
{	return (f1>f2? f1: f2) ;	};
///\return min value of float type
ZXH_DLL_EXPORT float minf( const float f1, const float f2 )
{	return (f1>f2? f2: f1) ;	};
///\return echo the application arguments
ZXH_DLL_EXPORT void echo_arguments( int argc, char* argv[] )
{
	if( glbVerboseOutput<=0 ) return ;
	for( int iarg=0 ; iarg<argc ; ++iarg )
		std::cout<< argv[iarg] << " " ;
	std::cout<< "\n" ;
};
///
ZXH_DLL_EXPORT void echo_verbose(int argc, char*argv[] )
{
	for( int iarg = 1 ; iarg < argc ; ++iarg )
	{
		if( strcmp( argv[iarg], "-v" ) == 0 )
		{
			if( iarg < argc-1 && !zxh::is_option(argv[iarg+1]) )
			{
				glbVerboseOutput = atoi( argv[ ++iarg ] ) ;
			}
			else glbVerboseOutput = 1 ;
		}
	}
}
ZXH_DLL_EXPORT void echo_zxh(int argc, char*argv[] )
{
	echo_verbose( argc, argv );
	echo_zxh();
#ifdef HAS_BUILD_RELEASE
	if( glbVerboseOutput>=4140 &&  glbVerboseOutput<=4143 )
		glbVerboseOutput-=4140 ; 
	else glbVerboseOutput = -1 ;
#endif
}

ZXH_DLL_EXPORT void echo_zxh()
{
	if( glbVerboseOutput==0 ) 
	{		
		std::cout<<"* ***************************************************************** *\n"; 
		std::cout<<"*      zxhproj, version 2.2 (c) 2004-2016 ZHUANG, XiaHai code.      *\n";
		std::cout<<"*          (zhuangxiahai@163.com). All rights reserved.             *\n";
		std::cout<<"* ***************************************************************** *\n";  
		return ;
	}
#ifdef HAS_BUILD_RELEASE
	std::cout<<"* ***************************************************************** *\n";
	std::cout<<"*  This is testing code and not free from bugs. The software is for *\n";
	std::cout<<"*  research purpose ONLY, and should NOT be used in any clinical    *\n";
	std::cout<<"*  related situations. Any usage is entirely at your own risk.      *\n";
	std::cout<<"*  zxhproj, version 2.2 (c) 2004-2016 ZHUANG, XiaHai Release Build. *\n";
	std::cout<<"*          (zhuangxiahai@163.com). All rights reserved.             *\n";
	std::cout<<"* ***************************************************************** *\n"; 
#else 
	std::cout<<"* ***************************************************************** *\n"; 
	std::cout<<"*      zxhproj, version 2.2 (c) 2004-2016 ZHUANG, XiaHai code.      *\n";
	std::cout<<"*          (zhuangxiahai@163.com). All rights reserved.             *\n";
	std::cout<<"* ***************************************************************** *\n";  
#endif
}
///
ZXH_DLL_EXPORT void VectorOP_Normalise( float v[], int dimension )		
{NormaliseVector(v,dimension);};
///
ZXH_DLL_EXPORT void NormaliseVector( float v[], int dimension )
{
	float mag = 0 ;
	for( int idim=0; idim<dimension; ++idim )
		mag += v[idim]*v[idim] ;
	mag = sqrt(mag) ;
	if( mag>0 )
		for( int idim=0; idim<dimension; ++idim )
			v[idim] /= mag ;
}
ZXH_DLL_EXPORT float VectorOP_Magnitude( const float v[], int dimension )		
{return sqrt(VectorOP_DotProduct( v,v,dimension )) ; } 
ZXH_DLL_EXPORT float VectorOP_DotProduct( const float v1[], const float v2[], const int dimension )
{
	float dot = 0 ; 
	for( int idim=0; idim<dimension; ++idim )
		dot += v1[idim]*v2[idim] ;
	return dot ;
}
///
ZXH_DLL_EXPORT void VectorOP_CrossProduct3D( const float A[], const float B[], float *AxB ) 
{
  AxB[0]=A[1]*B[2]-A[2]*B[1];
  AxB[1]=A[2]*B[0]-A[0]*B[2];
  AxB[2]=A[0]*B[1]-A[1]*B[0];
}
ZXH_DLL_EXPORT float VectorOP_Distance( const float v1[], const float v2[], int dimension )
{ 
	float dis = 0 ; 
	for( int idim=0; idim<dimension; ++idim )
		dis += (v1[idim]-v2[idim])*(v1[idim]-v2[idim]) ;
	return sqrt(dis) ;
}
ZXH_DLL_EXPORT float VectorOP_DistanceSquare( const float v1[], const float v2[], int dimension )
{ 
	float dis2 = 0 ; 
	for( int idim=0; idim<dimension; ++idim )
		dis2 += (v1[idim]-v2[idim])*(v1[idim]-v2[idim]) ;
	return dis2 ;
}
ZXH_DLL_EXPORT float VectorOP_DistanceFirstOrderNorm( const float v1[], const float v2[], int dimension )
{ 
	float norm = 0 ; 
	for( int idim=0; idim<dimension; ++idim )
		norm += zxh::abs(v1[idim]-v2[idim]) ;
	return norm ;
}
ZXH_DLL_EXPORT float VectorOP_Cosine( const float v1[], const float v2[], int dimension )
{
	float dot = 0, mag1=0, mag2=0 ; // cos(angle) = ( v1 dotproduct v2 )/( mag(v1)*mag(v2) )
	for( int idim=0; idim<dimension; ++idim )
	{
		dot += v1[idim] * v2[idim] ; 
		mag1+= v1[idim] * v1[idim] ; 
		mag2+= v2[idim] * v2[idim] ; 
	}	
	if( mag1*mag2 == 0 ) return 0 ; 
	return dot/sqrt(mag1*mag2) ;
}
///
ZXH_DLL_EXPORT void VectorOP_Substract( const float A[], const float B[], float *A_B, int dimension ) 
{
	for( int idim=0; idim<dimension; ++idim )
		A_B[idim] = A[idim]-B[idim] ; 
}
///
ZXH_DLL_EXPORT void VectorOP_Add( const float A[], const float B[], float *AaB, int dimension )
{
	for( int idim=0; idim<dimension; ++idim )
		AaB[idim] = A[idim]+B[idim] ; 
}
///
ZXH_DLL_EXPORT void VectorOP_Multiple( const float A[], const float f, float *Axf, int dimension ) 
{
	for( int idim=0; idim<dimension; ++idim )
		Axf[idim] = A[idim]*f ; 
}


ZXH_DLL_EXPORT void MatrixIdentity(float *pMatDes, const int iMatrixLength )
{
	float *des=0;
	for( int c=0; c<iMatrixLength; ++c ) //column
	{
		for( int l=0; l<iMatrixLength; ++l ) //row/line
		{
			des = &pMatDes[l*iMatrixLength+c];
			if( c==l )
				*des = 1;
			else
				*des = 0 ;
		}
	}
}

ZXH_DLL_EXPORT void SetMatrix(float *pMatDes, const float* pMatSrc1, const int iMatrixLength )
{
	for( int c=0; c<iMatrixLength; ++c ) //column
	for( int l=0; l<iMatrixLength; ++l ) //row/line
		pMatDes[l*iMatrixLength+c] = pMatSrc1[l*iMatrixLength+c] ;
}
ZXH_DLL_EXPORT void MultiplyMatrix(float *pJacDes, const float* pJacSrc1, const float* pJacSrc2, const int iMatrixLength )
{
	float *des=0;
	for( int c=0; c<iMatrixLength; ++c ) //column
	for( int l=0; l<iMatrixLength; ++l ) //row/line
	{
		des = &pJacDes[l*iMatrixLength+c];
		*des = 0 ;
		for( int i=0; i<iMatrixLength; ++i )
			*des += pJacSrc1[l*iMatrixLength+i] * pJacSrc2[i*iMatrixLength+c] ;
	}
}

ZXH_DLL_EXPORT void MultiplyMatrixVector(float *pVecDes, const float* pMatSrc, const float* pVecSrc, const int iMatrixLength )
{
	float *pdes = new float[iMatrixLength] ;
	for( int r=0 ; r< iMatrixLength; ++r )
	{
		pdes[r] = 0 ;
		for( int c=0 ; c< iMatrixLength; ++c )
			pdes[r] += pMatSrc[r*iMatrixLength+c] * pVecSrc[c] ;
	}
	for( int r=0 ; r< iMatrixLength; ++r )
		pVecDes[r] = pdes[r] ;
	delete []pdes ;
}
ZXH_DLL_EXPORT float DeterminentMatrix3D(const float *pMat, const int iMatrixRowLength )
{
	float ret = 0 ;
	ret += pMat[0*iMatrixRowLength+0]*pMat[1*iMatrixRowLength+1]*pMat[2*iMatrixRowLength+2] ;
	ret += pMat[0*iMatrixRowLength+1]*pMat[1*iMatrixRowLength+2]*pMat[2*iMatrixRowLength+0] ;
	ret += pMat[0*iMatrixRowLength+2]*pMat[1*iMatrixRowLength+0]*pMat[2*iMatrixRowLength+1] ;
	ret -= pMat[2*iMatrixRowLength+0]*pMat[1*iMatrixRowLength+1]*pMat[0*iMatrixRowLength+2] ;
	ret -= pMat[1*iMatrixRowLength+0]*pMat[0*iMatrixRowLength+1]*pMat[2*iMatrixRowLength+2] ;
	ret -= pMat[0*iMatrixRowLength+0]*pMat[2*iMatrixRowLength+1]*pMat[1*iMatrixRowLength+2] ;
	return ret ;
}
ZXH_DLL_EXPORT void InvertMatrix3D(float* pInvMat, const float* pForMat, const int iMatrixRowLength)
{ //http://mathworld.wolfram.com/MatrixInverse.html
	float Idet = 1/DeterminentMatrix3D( pForMat, iMatrixRowLength ) ;
	pInvMat[0*iMatrixRowLength+0] = Idet* (pForMat[1*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+2]-pForMat[1*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+1]) ; //
	pInvMat[0*iMatrixRowLength+1] = Idet* (pForMat[0*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+1]-pForMat[0*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+2]) ;
	pInvMat[0*iMatrixRowLength+2] = Idet* (pForMat[0*iMatrixRowLength+1]*pForMat[1*iMatrixRowLength+2]-pForMat[0*iMatrixRowLength+2]*pForMat[1*iMatrixRowLength+1]) ;
	pInvMat[1*iMatrixRowLength+0] = Idet* (pForMat[1*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+0]-pForMat[1*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+2]) ; //
	pInvMat[1*iMatrixRowLength+1] = Idet* (pForMat[0*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+2]-pForMat[0*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+0]) ;
	pInvMat[1*iMatrixRowLength+2] = Idet* (pForMat[0*iMatrixRowLength+2]*pForMat[1*iMatrixRowLength+0]-pForMat[0*iMatrixRowLength+0]*pForMat[1*iMatrixRowLength+2]) ;
	pInvMat[2*iMatrixRowLength+0] = Idet* (pForMat[1*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+1]-pForMat[1*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+0]) ; //
	pInvMat[2*iMatrixRowLength+1] = Idet* (pForMat[0*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+0]-pForMat[0*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+1]) ;
	pInvMat[2*iMatrixRowLength+2] = Idet* (pForMat[0*iMatrixRowLength+0]*pForMat[1*iMatrixRowLength+1]-pForMat[0*iMatrixRowLength+1]*pForMat[1*iMatrixRowLength+0]) ;
}
/// matrix op
ZXH_DLL_EXPORT void MatrixOP_SetIdentity(float *pMatDes, const int iMatrixLength )
{
	float *des=0;
	for( int c=0; c<iMatrixLength; ++c ) //column
	{
		for( int l=0; l<iMatrixLength; ++l ) //row/line
		{
			des = &pMatDes[l*iMatrixLength+c];
			if( c==l )
				*des = 1;
			else
				*des = 0 ;
		}
	}
}

ZXH_DLL_EXPORT void MatrixOP_SetMatrix(float *pMatDes, const float* pMatSrc1, const int iMatrixLength )
{
	for( int c=0; c<iMatrixLength; ++c ) //column
	for( int l=0; l<iMatrixLength; ++l ) //row/line
		pMatDes[l*iMatrixLength+c] = pMatSrc1[l*iMatrixLength+c] ;
}
ZXH_DLL_EXPORT void MatrixOP_Multiply(float *pJacDes, const float* pJacSrc1, const float* pJacSrc2, const int iMatrixLength )
{
	float *des=0;
	for( int c=0; c<iMatrixLength; ++c ) //column
	for( int l=0; l<iMatrixLength; ++l ) //row/line
	{
		des = &pJacDes[l*iMatrixLength+c];
		*des = 0 ;
		for( int i=0; i<iMatrixLength; ++i )
			*des += pJacSrc1[l*iMatrixLength+i] * pJacSrc2[i*iMatrixLength+c] ;
	}
}

ZXH_DLL_EXPORT void MatrixOP_MultiplyVector(float *pVecDes, const float* pMatSrc, const float* pVecSrc, const int iMatrixLength )
{
	float *pdes = new float[iMatrixLength] ;
	for( int r=0 ; r< iMatrixLength; ++r )
	{
		pdes[r] = 0 ;
		for( int c=0 ; c< iMatrixLength; ++c )
			pdes[r] += pMatSrc[r*iMatrixLength+c] * pVecSrc[c] ;
	}
	for( int r=0 ; r< iMatrixLength; ++r )
		pVecDes[r] = pdes[r] ;
	delete []pdes ;
}
ZXH_DLL_EXPORT float MatrixOP_DeterminentMatrix3D(const float *pMat, const int iMatrixRowLength )
{
	float ret = 0 ;
	ret += pMat[0*iMatrixRowLength+0]*pMat[1*iMatrixRowLength+1]*pMat[2*iMatrixRowLength+2] ;
	ret += pMat[0*iMatrixRowLength+1]*pMat[1*iMatrixRowLength+2]*pMat[2*iMatrixRowLength+0] ;
	ret += pMat[0*iMatrixRowLength+2]*pMat[1*iMatrixRowLength+0]*pMat[2*iMatrixRowLength+1] ;
	ret -= pMat[2*iMatrixRowLength+0]*pMat[1*iMatrixRowLength+1]*pMat[0*iMatrixRowLength+2] ;
	ret -= pMat[1*iMatrixRowLength+0]*pMat[0*iMatrixRowLength+1]*pMat[2*iMatrixRowLength+2] ;
	ret -= pMat[0*iMatrixRowLength+0]*pMat[2*iMatrixRowLength+1]*pMat[1*iMatrixRowLength+2] ;
	return ret ;
}
/// assume the arrays are both 3*3=9
ZXH_DLL_EXPORT void MatrixOP_InvertMatrix3D(float* pInvMat, const float* pForMat, const int iMatrixRowLength)
{ //http://mathworld.wolfram.com/MatrixInverse.html
	float Idet = 1/DeterminentMatrix3D( pForMat, iMatrixRowLength ) ;
	pInvMat[0*iMatrixRowLength+0] = Idet* (pForMat[1*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+2]-pForMat[1*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+1]) ; //
	pInvMat[0*iMatrixRowLength+1] = Idet* (pForMat[0*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+1]-pForMat[0*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+2]) ;
	pInvMat[0*iMatrixRowLength+2] = Idet* (pForMat[0*iMatrixRowLength+1]*pForMat[1*iMatrixRowLength+2]-pForMat[0*iMatrixRowLength+2]*pForMat[1*iMatrixRowLength+1]) ;
	pInvMat[1*iMatrixRowLength+0] = Idet* (pForMat[1*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+0]-pForMat[1*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+2]) ; //
	pInvMat[1*iMatrixRowLength+1] = Idet* (pForMat[0*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+2]-pForMat[0*iMatrixRowLength+2]*pForMat[2*iMatrixRowLength+0]) ;
	pInvMat[1*iMatrixRowLength+2] = Idet* (pForMat[0*iMatrixRowLength+2]*pForMat[1*iMatrixRowLength+0]-pForMat[0*iMatrixRowLength+0]*pForMat[1*iMatrixRowLength+2]) ;
	pInvMat[2*iMatrixRowLength+0] = Idet* (pForMat[1*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+1]-pForMat[1*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+0]) ; //
	pInvMat[2*iMatrixRowLength+1] = Idet* (pForMat[0*iMatrixRowLength+1]*pForMat[2*iMatrixRowLength+0]-pForMat[0*iMatrixRowLength+0]*pForMat[2*iMatrixRowLength+1]) ;
	pInvMat[2*iMatrixRowLength+2] = Idet* (pForMat[0*iMatrixRowLength+0]*pForMat[1*iMatrixRowLength+1]-pForMat[0*iMatrixRowLength+1]*pForMat[1*iMatrixRowLength+0]) ;
}
// 3x4 matrix
ZXH_DLL_EXPORT void MatrixOP_InvertMatrix4D3D(float* pInvMat, const float* pForMat, const int iMatrixRowLength)
{
	zxhmat44 inv, forward ; 
	for( int row=0; row<3 ; row++ )
		for( int col=0; col<4; col++ ) 
			forward.m[row][col] = pForMat[row*iMatrixRowLength+col] ;  
	forward.m[3][3] = 1 ; 
	MatrixOP_InvertMatrix4D3D( inv, forward ) ; 
	for( int row=0; row<3 ; row++ )
		for( int col=0; col<4; col++ ) 
			pInvMat[row*iMatrixRowLength+col] = inv.m[row][col] ;
}
ZXH_DLL_EXPORT void MatrixOP_InvertMatrix4D3D( zxhmat44 &inv, const zxhmat44 &forward )
{
   float r11,r12,r13,r21,r22,r23,r31,r32,r33,v1,v2,v3 , deti ;
                                                       /*  INPUT MATRIX IS:  */
   r11 = forward.m[0][0]; r12 = forward.m[0][1]; r13 = forward.m[0][2];  /* [ r11 r12 r13 v1 ] */
   r21 = forward.m[1][0]; r22 = forward.m[1][1]; r23 = forward.m[1][2];  /* [ r21 r22 r23 v2 ] */
   r31 = forward.m[2][0]; r32 = forward.m[2][1]; r33 = forward.m[2][2];  /* [ r31 r32 r33 v3 ] */
   v1  = forward.m[0][3]; v2  = forward.m[1][3]; v3  = forward.m[2][3];  /* [  0   0   0   1 ] */

   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13 ;

   if( deti != 0.0l ) deti = 1.0l / deti ;

   inv.m[0][0] = deti*( r22*r33-r32*r23) ;
   inv.m[0][1] = deti*(-r12*r33+r32*r13) ;
   inv.m[0][2] = deti*( r12*r23-r22*r13) ;
   inv.m[0][3] = deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3
                     -r22*v1*r33-r32*r13*v2+r32*v1*r23) ;

   inv.m[1][0] = deti*(-r21*r33+r31*r23) ;
   inv.m[1][1] = deti*( r11*r33-r31*r13) ;
   inv.m[1][2] = deti*(-r11*r23+r21*r13) ;
   inv.m[1][3] = deti*( r11*r23*v3-r11*v2*r33-r21*r13*v3
                     +r21*v1*r33+r31*r13*v2-r31*v1*r23) ;

   inv.m[2][0] = deti*( r21*r32-r31*r22) ;
   inv.m[2][1] = deti*(-r11*r32+r31*r12) ;
   inv.m[2][2] = deti*( r11*r22-r21*r12) ;
   inv.m[2][3] = deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3
                     -r21*r32*v1-r31*r12*v2+r31*r22*v1) ;

   inv.m[3][0] = inv.m[3][1] = inv.m[3][2] = 0.0l ;
   inv.m[3][3] = (deti == 0.0l) ? 0.0l : 1.0l ; /* failure flag if deti == 0 */
}
/////////--------
ZXH_DLL_EXPORT float BSpline(float u)
{
	if( u<0 ) u = -u ;
	if( u<=1 ) return 0.5*u*u*u-u*u+2.0/3.0;
	if( u<=2 )
	{ u = 2-u ;return u*u*u/6.0;
	}
	return 0 ;
} ;
/// \return
ZXH_DLL_EXPORT float BSplineDerivative(float u)
{
	float minus=1;
	if( u<0 ) u = -u,minus=-1 ;
	if( u<=1 ) return  (1.5*u*u-2.0*u)*minus;
	if( u<=2 )
	{
		u = 2-u ;
		return (-0.5*u*u*minus) ;
	}
	return 0 ;
}
///
ZXH_DLL_EXPORT float BSplineSecondDerivative ( float u )
{
	if( u<0 ) u = -u ;
	if( u<=1 ) return  (3*u-2);
	if( u<=2 ) return (2-u) ;
	return 0 ;
} ;
/// \return BSpline
ZXH_DLL_EXPORT float BSplinei(int iOrd,float u)
{
	float v;
	switch(iOrd)
	{
	case -1:
		v=1.0f-u;
		return v*v*v/6.0;
	case 0:
		v=u*u;
		return 0.5*u*v-v+2.0/3.0;
	case 1:
		v=u*u;
		return (-3.0*u*v+3.0*v+3.0*u+1)/6.0;
	case 2:return u*u*u/6.0;
	}
	return 0;
} ;
ZXH_DLL_EXPORT float DerivativeOfBSplinei(int iOrd, float u)
{ // -1,0,1,2
	float v;
	switch(iOrd)
	{
	case -1:
		v=1.0-u;
		return v*v*(-0.5);
	case 0:return (3.0*u-4.0)*u*0.5;
	case 1:return (u+0.5-1.5*u*u);
	case 2:return u*u*0.5;
	}
	return 0;
};
ZXH_DLL_EXPORT float SecondDerivativeOfBSplinei(int iOrd,float u)
{ // -1,0,1,2
	switch(iOrd)
	{
	case -1:
		return 1.0f-u;
	case 0:
		return (3.0f*u-2.0f);
	case 1:return (1.0f-3.0f*u);
	case 2:return u;
	}
	return 0;
};
 
ZXH_DLL_EXPORT bool IsWindowsOS()
{
#ifndef ZXH_IS_NOT_WINDOWS
	return true ;
#else
	return false ;
#endif 
};
///

ZXH_DLL_EXPORT size_t GetSizeOfMemoryMB()
{ 
#ifndef ZXH_IS_NOT_WINDOWS 
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	return status.ullTotalPhys/(1024*1024); 
#else //ZXH_IS_NOT_WINDOWS 
#ifdef _SC_PHYS_PAGES
	long pages = sysconf(_SC_PHYS_PAGES); //_SC_AVPHYS_PAGES  available
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size/(1024*1024); 
#else
	return 1024*4 ; // 4G zxhtodo for Mac OSX
#endif
#endif
};
///	
ZXH_DLL_EXPORT bool ReadOstiumCoordinateFromTxt( std::string filename, float * wco ) 
{
	std::ifstream instream(filename.c_str());
	if (instream.fail()==true)
	{
		std::cerr<<"error: Open ostium file "<<filename<<" Failed!\n";
		return false;
	}
	instream.seekg(0,std::ios_base::beg);  
	char buffer[1024];
	std::string sLine,sContent,sComment;
	do{
			instream.getline(buffer,1024);
			sLine=buffer;
			zxh::ParseStringLine(sContent,sComment,sLine);
			zxh::trim_both(sContent);
	}while(sContent.empty()&&instream.eof()==false&&instream.fail()==false);
	if(instream.eof()==true||instream.fail()==true)
		return false;
	std::istringstream strstream(sContent.c_str());
	for(int c=0;c<3;++c)//center over
	{
		strstream>>wco[c];
	}  
	instream.close();
	return true ;   
}
///
ZXH_DLL_EXPORT bool AppendSaveOstiumCoordinateToTxt( float * wco, std::string filename ) 
{    
	char buffer[1024];
	sprintf(buffer, "%f \t %f \t %f\n", wco[0],wco[1],wco[2]) ;
	return SaveString2FileAppend( buffer, filename ) ;  
}


}// end of zxh
//#ifdef __cplusplus
//}
//#endif

