
/*=========================================================================

  Program:   class zxhPoint
  Author:	 Zhuang, Xiahai
  Module:    $RCSfle: zxhPoint.h    $
  Language:  C++
  Date:      $Date: From  2012-08 $
  Version:   $Revision: 1.0$
			 $Revision: 1.1$

  Update log: 2012-8-30 some new member functions were included

=========================================================================*/

#ifndef zxhPoint_H
#define zxhPoint_H

#include <iostream>
#include <string>                                               // for STL string
#include <list>
#include <vector>
#include <cmath>
#include "zxh.h"


using namespace std;

///
/// \class zxhPointT
/// \brief
/// zxhPoint.h: interface for the 3D Point construction, Set Value, Get Value, calculation, etc.
/// we can use this class to define a 3d Point or a 3d vector in the x-y-z axis.
/// example:
/// if the coordinate of the point or vector is (a,b,c), we can define as follows:
/// zxhPoint point(a,b,c). 
/// \ingroup zxhPointT
///
 

template <class DataType=float>
class zxhPointT
{
public:
	///a new point in the 3d space
	zxhPointT() {x=0.0; y=0.0; z=0.0;};
	///a copy of point a 
	zxhPointT(const zxhPointT &a) {x = a.x; y = a.y ; z = a.z;};
	///a new point with coordinate (xx,yy,zz)
	zxhPointT(DataType xx ,DataType yy, DataType zz ) {x = xx; y = yy ; z = zz;};
	///
	virtual ~zxhPointT() { };
	///
	zxhPointT<DataType> operator+ (const zxhPointT &P);
	///
	zxhPointT<DataType> operator- (const zxhPointT &P);
	///
	zxhPointT<DataType> operator+ (DataType m);
	///
	zxhPointT<DataType> operator* (DataType m);
	///
	zxhPointT<DataType> operator/ (DataType m);
	/// 
	bool operator== (const zxhPointT &P);

	///set the point as (a,b,c)
	void Set(DataType a, DataType b , DataType c ) {x =a ; y = b ;z = c;};
	///calculate the norm of the point, sqrt(x*x+y*y+z*z)
	DataType Norm() {return sqrt(x*x + y*y +z*z);};
	///dot multiply 
	zxhPointT<DataType> DotMul(const zxhPointT &P);
	///cross multiply 
	zxhPointT<DataType> CrossMul(const zxhPointT &P);
	///vector product with dot multiply
	DataType VecProduct(const zxhPointT &P);
	///display the point on the console
	void Display() {std::cout<<x<<'\t'<<y<<'\t'<<z<<'\n';}


	//2012-8-30
	///normalize
	void Normalize();
	///update the members, x,y,z, cross multiply
	void MulCross(zxhPointT p1, zxhPointT p2); 
	/// 
	zxhPointT<DataType> operator+= (const zxhPointT &P);
	///
	zxhPointT<DataType> operator-= (const zxhPointT &P);
	///
	zxhPointT<DataType> operator*= (DataType m);
	///
	zxhPointT<DataType> operator/= (DataType m);
	/// 

	DataType x, y, z;
};




template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator+( const zxhPointT &P )
{
	zxhPointT<DataType> temp;
	temp.x=x+P.x;
	temp.y=y+P.y;
	temp.z=z+P.z;
	return temp;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator-( const zxhPointT &P )
{
	zxhPointT<DataType> temp;
	temp.x=x-P.x;
	temp.y=y-P.y;
	temp.z=z-P.z;
	return temp;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator+( DataType m )
{
	zxhPointT<DataType> temp;
	temp.x=x+m;
	temp.y=y+m;
	temp.z=z+m;
	return temp;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator*( DataType m )
{
	zxhPointT<DataType> temp;
	temp.x=x*m;
	temp.y=y*m;
	temp.z=z*m;
	return temp;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator/( DataType m )
{
	if (m==0.0)
	{
		cerr<<"the m should not be zero"<<endl;
		exit(1);
	}
	zxhPointT<DataType> temp;
	temp.x=x/m;
	temp.y=y/m;
	temp.z=z/m;
	return temp;
}


template <class DataType>
bool zxhPointT<DataType>::operator==( const zxhPointT &P )
{
	if ((zxh::abs(x-P.x)==0)&&(zxh::abs(y-P.y)==0)&&(zxh::abs(z-P.z)==0))
	{
		return true;
	}
	return false;	
}
 

template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::DotMul( const zxhPointT &P )
{
	zxhPointT<DataType> temp;
	temp.x=x*P.x;
	temp.y=y*P.y;
	temp.z=z*P.z;
	return temp;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::CrossMul( const zxhPointT &P )
{
	zxhPointT<DataType> temp;
	temp.x=y*P.z-z*P.y;
	temp.y=z*P.x-x*P.z;
	temp.z=x*P.y-y*P.x;
	return temp;
}


template <class DataType>
DataType zxhPointT<DataType>::VecProduct( const zxhPointT &P )
{
	zxhPointT<DataType> temp;
	temp=this->DotMul(P);
	return (temp.x+temp.y+temp.z);
}


//2012-8-30

template <class DataType>
void zxhPointT<DataType>::Normalize()
{
	DataType nor=Norm();
	if (nor!=0)
	{
		x/=nor;
		y/=nor;
		z/=nor;
	}	
}


template <class DataType>
void zxhPointT<DataType>::MulCross( zxhPointT p1, zxhPointT p2 )
{
	x = p1.y*p2.z - p1.z*p2.y ; 
	y = p1.z*p2.x - p1.x*p2.z; 
	z = p1.x*p2.y - p1.y*p2.x; 
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator+=( const zxhPointT &P )
{
	x+=P.x; 
	y+=P.y; 
	z+=P.z; 
	return *this;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator-=( const zxhPointT &P )
{
	x-=P.x; 
	y-=P.y; 
	z-=P.z; 
	return *this;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator*=(DataType m)
{
	x*=m; 
	y*=m; 
	z*=m; 
	return *this;
}


template <class DataType>
zxhPointT<DataType> zxhPointT<DataType>::operator/=(DataType m)
{
	x/=m; 
	y/=m; 
	z/=m; 
	return *this;
}
 
typedef zxhPointT<float> zxhPointF;


#endif
