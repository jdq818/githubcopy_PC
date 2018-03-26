
/*=========================================================================

  Program:   class zxhMatrixOP
  Author:	 ZHUANG, Xiahai
  Module:    $RCSfle: zxhMatrixOP.h    $
  Language:  C++
  Date:      $Date: From  2012-10 $
  Version:   $Revision: 2.0$

  Update log:  

=========================================================================*/


#ifndef ZXHMATRIXOP_H
#define ZXHMATRIXOP_H

#include "zxh.h"  
#include "zxhMatrix.h"
#ifndef ZXH_MaxIterationStepForComputeSquareRootOfMatrix
#define ZXH_MaxIterationStepForComputeSquareRootOfMatrix 1024
#endif
 
///
/// \class zxhMatrixOP
/// \brief 
/// \ingroup zxhMatrixOP
/// 

#include <iostream>
#include <string>                                               // for STL string
#include <list>
#include <vector>
#include <cmath>
  
template <class DataType=DataTypeDefault>
class zxhMatrixOPT: public zxhMatrixT<DataType>
{
public:
	///a void matrix
	zxhMatrixOPT( ); 
	///Destructor function
	virtual ~zxhMatrixOPT( ); 
 
	/// return whether square root been successfully computed, both pSquareRootMatrix and pMatrix are iMatrixLength*iMatrixLength dimension
	///	1) Babylonian iteration method Xk+1 = 1/2(Xk+AXk^{-1}); X0=Identity ; X^2 = A
	/// 2) Denman¨CBeavers iteration Xk+1=1/2(Xk+Zk^{-1}), Zk+1=1/2(Zk+Xk^{-1}); X0=A,Z0=I; X=A^{1/2},Z=A^{-1/2}
	/// 3) Jordan decomposition (not implemented yet, zxhtodo todo) 
	static bool ComputeSquareRootOfMatrix3D( zxhMatrixT<DataType>*pSquareRoot, const zxhMatrixT<DataType>*pSourceMatrix, float TOL=1.0e-5 ) ;
	/// 
	static bool SetSquareMatrixToIdentity(zxhMatrixT<DataType>*pMatrix) 
	{
		if( pMatrix==0||pMatrix->IsSqrMatrix() == false )
		{
			std::cerr<<"error: SetSquareMatrixToIdentity fail\n" ;
			return false ;
		}
		int nM = pMatrix->GetColNum() ;
		pMatrix->SetAllValue( DataType(0) );
		for( int ic=0; ic<nM; ++ic ) 
			pMatrix->SetValue( DataType(1), ic,ic ) ; 
		return true ;
	}
	/// magnitude or normal of data
	static float ComputeMagnitudeOfData( zxhMatrixT<DataType>*pMatrix )  
	{
		if( pMatrix==0 ) return -1 ; 
		float mag=0;
		for( int ip=0; ip<pMatrix->GetColNum()*pMatrix->GetRowNum(); ++ip )
			mag += pMatrix->GetData()[ip]*pMatrix->GetData()[ip];
		return sqrt(mag);
	}
	/// copy from zxh.h
	static float ComputeDeterminentMatrix3D( const zxhMatrixT<DataType>*pMatrix )  ;
	///
	static bool ComputeInvertOfMatrix3D( zxhMatrixT<DataType>*pInvMatrix, const zxhMatrixT<DataType>*pSourceMatrix )  ;   

};
 
/// return whether square root been successfully computed, both pSquareRootMatrix and pMatrix are iMatrixLength*iMatrixLength dimension
///	1) Babylonian iteration method Xk+1 = 1/2(Xk+AXk^{-1}); X0=Identity ; X^2 = A
/// 2) Denman¨CBeavers iteration Xk+1=1/2(Xk+Zk^{-1}), Zk+1=1/2(Zk+Xk^{-1}); X0=A,Z0=I; X=A^{1/2},Z=A^{-1/2}
/// 3) Jordan decomposition (not implemented yet, zxhtodo todo)
template <class DataType>
bool	zxhMatrixOPT<DataType>::ComputeSquareRootOfMatrix3D( zxhMatrixT<DataType>*pSquareRoot, const zxhMatrixT<DataType>*pSourceMatrix, float TOL)
{ 
	if( pSourceMatrix==0||pSquareRoot==0||pSourceMatrix->GetColNum()!=3||pSourceMatrix->GetRowNum()!=3 ) 
		return false ; 
	int nM = 3 ;  
	zxhMatrixT<DataType> Xk(nM, nM ) ; 
	zxhMatrixT<DataType> invXk(nM, nM ) ; 
	zxhMatrixT<DataType> Xkplus1(nM, nM ) ; 

	float error = 0 ;
	// 1) Babylonian iteration method Xk+1 = 1/2(Xk+A*Xk^{-1}); X0=Identity ; X^2 = A
	if( SetSquareMatrixToIdentity( &Xkplus1 ) == false )
		return false ;
	for( int iit=0; iit<ZXH_MaxIterationStepForComputeSquareRootOfMatrix; ++iit )
	{
		Xk = Xkplus1 ;  
		if( ComputeInvertOfMatrix3D( &invXk, &Xk ) == false )
		{
			break ; // use method 2
			//std::cerr<<"error: ComputeInvertOfMatrix3D fail and hence ComputeSquareRootOfMatrix3D failed\n" ;
		}
		Xkplus1 = (Xk + *pSourceMatrix*invXk)*(DataType)0.5 ; 

		// test whether converge to less than TOL difference 
		zxhMatrixT<DataType> diff ;
		diff = Xkplus1*Xkplus1 - *pSourceMatrix ; 
		error = ComputeMagnitudeOfData( &diff ) ;
		if( error < TOL )
		{
			*pSquareRoot = Xkplus1 ; 
			return true ;
		}
	} 
	// 2) Denman¨CBeavers iteration Xk+1=1/2(Xk+Zk^{-1}), Zk+1=1/2(Zk+Xk^{-1}); X0=A,Z0=I; X=A^{1/2},Z=A^{-1/2}
	zxhMatrixT<DataType> Zk(nM, nM ) ; 
	zxhMatrixT<DataType> invZk(nM, nM ) ; 
	zxhMatrixT<DataType> Zkplus1(nM, nM ) ; 
	Xkplus1 = *pSourceMatrix ; 
	if( SetSquareMatrixToIdentity( &Zkplus1 ) == false ) return false ;
	for( int iit=0; iit<ZXH_MaxIterationStepForComputeSquareRootOfMatrix; ++iit )
	{
		Xk = Xkplus1 ; 
		Zk = Zkplus1 ;   
		if( ComputeInvertOfMatrix3D( &invXk, &Xk ) == false )
		{
			break ; // use method 3
			//std::cerr<<"error: ComputeInvertOfMatrix3D fail and hence ComputeSquareRootOfMatrix3D failed\n" ;
		}
		if( ComputeInvertOfMatrix3D( &invZk, &Zk ) == false )
		{
			break ; // use method 3
			//std::cerr<<"error: ComputeInvertOfMatrix3D fail and hence ComputeSquareRootOfMatrix3D failed\n" ;
		}
		
		Xkplus1 = ( Xk + invZk ) * (DataType)0.5 ;
		Zkplus1 = ( Zk + invXk ) * (DataType)0.5 ;
		// test whether converge to less than TOL difference 
		zxhMatrixT<DataType> diff ;
		diff = Xkplus1*Xkplus1 - *pSourceMatrix ; 

		error = ComputeMagnitudeOfData( &diff ) ;
		if( error < TOL )
		{
			*pSquareRoot = Xkplus1 ; 
			return true ;
		}
	} 

	// 3) Jordan decomposition (not implemented yet, zxhtodo todo)
	std::cerr<<"error: square root affine transformation fail due to the estimated error "<<error<<" > "<<TOL<<" of tollerance\n" ;
	return false ;
}

template <class DataType>
float	zxhMatrixOPT<DataType>::ComputeDeterminentMatrix3D( const zxhMatrixT<DataType>*pMatrix )  
{
	if( pMatrix==0 || pMatrix->GetColNum()!=3 || pMatrix->GetRowNum()!=3 ) 
		return 0 ;  
	const DataType * pMat = pMatrix->GetData() ;  
	float ret = 0 ; 
	ret += pMat[0*3+0]*pMat[1*3+1]*pMat[2*3+2] ;
	ret += pMat[0*3+1]*pMat[1*3+2]*pMat[2*3+0] ;
	ret += pMat[0*3+2]*pMat[1*3+0]*pMat[2*3+1] ;
	ret -= pMat[2*3+0]*pMat[1*3+1]*pMat[0*3+2] ;
	ret -= pMat[1*3+0]*pMat[0*3+1]*pMat[2*3+2] ;
	ret -= pMat[0*3+0]*pMat[2*3+1]*pMat[1*3+2] ;
	return ret ;
}

template <class DataType>
bool	zxhMatrixOPT<DataType>::ComputeInvertOfMatrix3D( zxhMatrixT<DataType>*pInvMatrix, const zxhMatrixT<DataType>*pSourceMatrix ) 
{  
	if( pSourceMatrix==0||pInvMatrix==0||pSourceMatrix->GetColNum()!=3||pSourceMatrix->GetRowNum()!=3 ) 
		return false ; 
	float Idet = ComputeDeterminentMatrix3D( pSourceMatrix ) ;
	if( zxh::abs(Idet) < ZXH_FloatPrecision )
	{
		std::cerr<<"error: can not invert a matrix with zero determinant, ComputeInvertOfMatrix3D failed\n" ;
		return false ;
	}
	Idet = 1/Idet ;  
	DataType pInvMat[9];
	const DataType* pForMat = pSourceMatrix->GetData() ;  
	pInvMat[0*3+0] = Idet* (pForMat[1*3+1]*pForMat[2*3+2]-pForMat[1*3+2]*pForMat[2*3+1]) ; //
	pInvMat[0*3+1] = Idet* (pForMat[0*3+2]*pForMat[2*3+1]-pForMat[0*3+1]*pForMat[2*3+2]) ;
	pInvMat[0*3+2] = Idet* (pForMat[0*3+1]*pForMat[1*3+2]-pForMat[0*3+2]*pForMat[1*3+1]) ;
	pInvMat[1*3+0] = Idet* (pForMat[1*3+2]*pForMat[2*3+0]-pForMat[1*3+0]*pForMat[2*3+2]) ; //
	pInvMat[1*3+1] = Idet* (pForMat[0*3+0]*pForMat[2*3+2]-pForMat[0*3+2]*pForMat[2*3+0]) ;
	pInvMat[1*3+2] = Idet* (pForMat[0*3+2]*pForMat[1*3+0]-pForMat[0*3+0]*pForMat[1*3+2]) ;
	pInvMat[2*3+0] = Idet* (pForMat[1*3+0]*pForMat[2*3+1]-pForMat[1*3+1]*pForMat[2*3+0]) ; //
	pInvMat[2*3+1] = Idet* (pForMat[0*3+1]*pForMat[2*3+0]-pForMat[0*3+0]*pForMat[2*3+1]) ;
	pInvMat[2*3+2] = Idet* (pForMat[0*3+0]*pForMat[1*3+1]-pForMat[0*3+1]*pForMat[1*3+0]) ; 
	pInvMatrix->MallocNewData( pInvMat, 3, 3 ) ; 
	return true ;
} ;   
 

#endif



