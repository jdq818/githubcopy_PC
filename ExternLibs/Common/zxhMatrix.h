
/*=========================================================================

  Program:   class zxhMatrix
  Author:	 ZHUANG, Xiahai & SHAO, Yi
  Module:    $RCSfle: zxhMatrix.h    $
  Language:  C++
  Date:      $Date: From  2012-08 $
  Version:   $Revision: 1.0$
			 $Revision: 2.0$

  Update log: template <class DataType=DataTypeDefault> class zxhMatrixT,
				DataTypeDefault=float
			2012-9-4 modified with ZXH's comments. 
					 add boundary judgement and memory judgement.
					 marked with //zxh add.
			2012-10-7 zxh add some functions: MallocNewData() GetData()

=========================================================================*/


#ifndef MIIMATRIX_H
#define MIIMATRIX_H

///
/// \class zxhMatrix
/// \brief
/// zxhMatrix.h: interface for the Matrix construction, Set Value, Get Value, calculation, etc.
/// example: zxhMatrix a(3,4);
/// located in the memory as follows:
///	a={ a(0,0), a(0,1), a(0,2), a(0,3),
///		a(1,0), a(1,1), a(1,2), a(1,3),
///		a(2,0), a(2,1), a(2,2), a(2,3),};
/// here m_RowNum=3, m_ColNum=4
/// if we want to get or set the value a(1,2), we could use GetValue(1,2) or SetValue(value,1,2)
/// that's to say, the first number and second number correspond with row and col, respectively.
/// \ingroup zxhMatrix
///


#include <iostream>
#include <string>                                               // for STL string
#include <list>
#include <vector>
#include <cmath>
#include "zxh.h"
 


template <class DataType=float>
class zxhMatrixT
{
public:
	///a void matrix
	zxhMatrixT( );
	///a copy of matrix mtx
	zxhMatrixT(const zxhMatrixT& mtx );
	///a new RowNum*ColNum zero matrix
	zxhMatrixT(int RowNum, int ColNum);
	///a new RowNum*ColNum matrix with data store
	zxhMatrixT(const DataType* data, int RowNum, int ColNum);
	///Destructor function
	virtual ~zxhMatrixT( );
	///get the number of row
	int GetRowNum() const {return m_RowNum;} 
	///get the number of column
	int GetColNum() const {return m_ColNum;}
	///set all the elements of matrix with data
	void SetAllValue(DataType data);
	///set the element (row,col) with data, m:0-m_RowNum-1, n:0-m_ColNum-1
	void SetValue(DataType data, int row, int col);
	///get the value of element (row,col), m:0-m_RowNum-1, n:0-m_ColNum-1
	DataType GetValue(int row, int col) const;
	///multiply the matrix with a double type number
	zxhMatrixT<DataType> MatrixMultNum(DataType num);
	///dot multiply
	zxhMatrixT<DataType> MatrixDotMultiply(zxhMatrixT &m) const;
	///matrix transpose
	zxhMatrixT<DataType> MatrixTranspose() const;
	/// 
	zxhMatrixT<DataType> operator= (const zxhMatrixT &m);
	///
	zxhMatrixT<DataType> operator+ (const zxhMatrixT &m) const;
	///
	zxhMatrixT<DataType> operator- (const zxhMatrixT &m) const;
	///
	zxhMatrixT<DataType> operator* (const zxhMatrixT &m) const;
	///
	zxhMatrixT<DataType> operator* (DataType num) const;
	///determine whether the matrix is a square matrix
	bool IsSqrMatrix();
	///determine whether the matrix is symmetric
	bool IsSymmetric();
	///find the maximum positive of the matrix
	bool FindPosMax(int &row, int &col);
	///copy the m_data[] of Matrix
	void CopyData(DataType Copy[]) const;
	///display the matrix on the console 
	void Display();


	///calculate the EigenvalueVector of a real symmetric matrix 
	// zxh add: it is difficult to know what is a, and what is v, and eps, jt, 
	//  better to provide comments or use a more meaningful name for the variables
	static int EigenvalueVectorRealSymmetryJacobi(zxhMatrixT& a,  
		zxhMatrixT& v, DataType eps, int jt);

	
	///
	virtual bool MallocNewData(int RowNum, int ColNum)
	{ 
		if ( (RowNum<=0) || (ColNum <=0) )
			return false ; 
		else
		{ 
	 		if( RowNum*ColNum!=m_RowNum*m_ColNum )
			{
				if( m_data!=0 )
					delete [] m_data ; 
				m_data=new DataType[RowNum*ColNum];
			}
			m_RowNum=RowNum;
			m_ColNum=ColNum;
			for(int i=0;i<RowNum*ColNum;i++)
				m_data[i]=0;
		}
		return true ;
	}
	///
	virtual bool MallocNewData(const DataType* data, int RowNum, int ColNum)
	{ 
		if ( (RowNum<=0) || (ColNum <=0) )
			return false ; 
		else
		{ 
	 		if( RowNum*ColNum!=m_RowNum*m_ColNum )
			{
				if( m_data!=0 )
					delete [] m_data ; 
				m_data=new DataType[RowNum*ColNum];
			}
			m_RowNum=RowNum;
			m_ColNum=ColNum;
			for(int i=0;i<RowNum*ColNum;i++)
				m_data[i]=data[i];
		}
		return true ;
	}
	///
	virtual const DataType* GetData() const	{return m_data;};
protected:
	int m_RowNum;
	int m_ColNum;
	DataType* m_data;


};



template <class DataType>
zxhMatrixT<DataType>::zxhMatrixT()
{
	m_RowNum=0;
	m_ColNum=0;

	m_data=0;
}

template <class DataType>
zxhMatrixT<DataType>::zxhMatrixT(const zxhMatrixT& mtx )
{
	
	// zxh add 
	if ( (mtx.m_RowNum<=0) || (mtx.m_ColNum <=0) )
	{		
		zxhMatrixT();
	}
	else
	{
		m_RowNum=mtx.m_RowNum;
		m_ColNum=mtx.m_ColNum;
		m_data=new DataType[m_RowNum*m_ColNum];
		for (int i=0;i<m_RowNum*m_ColNum;i++)
		{
			m_data[i]=mtx.m_data[i];
		}
	}
	
}

template <class DataType>
zxhMatrixT<DataType>::zxhMatrixT(int RowNum, int ColNum)
{
	// zxh add: 
	//if (0==RowNum && 0==ColNum)
	if ( (RowNum<=0) || (ColNum <=0) )
	{		
		zxhMatrixT();
	}
	else
	{
		// zxh add
// 		if( m_data!=0 )
// 			delete [] m_data ; 

		m_data=new DataType[RowNum*ColNum];
		m_RowNum=RowNum;
		m_ColNum=ColNum;
		for(int i=0;i<RowNum*ColNum;i++)
			m_data[i]=0;
	}

}

template <class DataType>
zxhMatrixT<DataType>::zxhMatrixT(const DataType* data, int RowNum, int ColNum)
{
	// zxh add: 
	//if (0==RowNum && 0==ColNum)
	if ( (RowNum<=0) || (ColNum <=0) )
	{
		zxhMatrixT();
	}
	else
	{
		// zxh add
// 		if( m_data!=0 )
// 			delete [] m_data ; 

		m_data=new DataType[RowNum*ColNum];
		m_RowNum=RowNum;
		m_ColNum=ColNum;
		for(int i=0;i<RowNum*ColNum;i++)
			m_data[i]=data[i];
	}

}

template <class DataType>
zxhMatrixT<DataType>::~zxhMatrixT()
{
	// zxh add
	if( (m_data!=0) && (m_RowNum>0) && (m_ColNum>0) )
		delete [] m_data ; 

	//delete [] m_data;
}

template <class DataType>
void zxhMatrixT<DataType>::SetAllValue( DataType data )
{
	for(int i=0;i<m_RowNum*m_ColNum;i++)
		m_data[i]=data;
}

template <class DataType>
void zxhMatrixT<DataType>::SetValue( DataType data, int row, int col)
{
	m_data[row*m_ColNum+col]=data;
}

template <class DataType>
DataType zxhMatrixT<DataType>::GetValue( int row, int col) const
{
	return m_data[row*m_ColNum+col];
}

template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::MatrixMultNum( DataType num )
{
	if (m_RowNum<=0 || m_ColNum<=0)
	{
		return *this;
	} 
	else
	{
		for(int i=0;i<m_RowNum*m_ColNum;i++)
			m_data[i]*=num;
		return *this;
	}
	
}


template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::MatrixDotMultiply( zxhMatrixT &m ) const
{
	if (m_RowNum!=m.m_RowNum || m_ColNum!=m.m_ColNum)
	{
		std::cerr<<"Dimension not matched";
		abort();
	}
	zxhMatrixT temp(m_RowNum,m_ColNum);
	for (int i=0;i<m_RowNum*m_ColNum;i++)
	{
		temp.m_data[i]=m_data[i]*m.m_data[i];
	}
	return temp; 
	// zxh add, this is correct compared to the previous function, but because you didn't have operator =, this is useless
}

template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::MatrixTranspose() const
{
	zxhMatrixT temp(m_ColNum,m_RowNum);
	for (int i=0;i<m_RowNum;i++)
	{
		for (int j=0;j<m_ColNum;j++)
		{
			temp.m_data[j*m_RowNum+i]=m_data[i*m_ColNum+j];
		}
	}
	return temp; // zxh add, same problem as previous
}


template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::operator=( const zxhMatrixT &m )
{
	if (this==&m)
	{
		return *this;
	} 
	else
	{
		if( (m_RowNum>0) && (m_ColNum>0) )
			delete [] m_data ; 
		
		if ( (m.m_RowNum<=0) || (m.m_ColNum <=0) )
		{		
			zxhMatrixT();
		}
		else
		{
			m_RowNum=m.m_RowNum;
			m_ColNum=m.m_ColNum;
			m_data=new DataType[m_RowNum*m_ColNum];
			for (int i=0;i<m_RowNum*m_ColNum;i++)
			{
				m_data[i]=m.m_data[i];
			}
		}
	}
	return *this;
}


template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::operator+(const zxhMatrixT &m ) const
{
	if ( m_RowNum!=m.m_RowNum || m_ColNum!=m.m_ColNum )
	{
		std::cerr<<"Dimension not matched";
		abort();
	}
	zxhMatrixT temp(m_RowNum,m_ColNum);
	for (int i=0;i<m_RowNum*m_ColNum;i++)
	{
		temp.m_data[i]=m_data[i]+m.m_data[i];
	}
	return temp;// zxh add, same problem as previous
}

template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::operator-(const zxhMatrixT &m ) const
{
	if ( m_RowNum!=m.m_RowNum || m_ColNum!=m.m_ColNum )
	{
		std::cerr<<"Dimension not matched";
		abort();
	}
	zxhMatrixT temp(m_RowNum,m_ColNum);
	for (int i=0;i<m_RowNum*m_ColNum;i++)
	{
		temp.m_data[i]=m_data[i]-m.m_data[i];
	}
	return temp;// zxh add, same problem as previous
}

template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::operator*(const zxhMatrixT &m ) const
{
	if(m_ColNum!=m.m_RowNum)
	{
		std::cerr<<"Dimension not matched";
		abort();
	}
	zxhMatrixT temp=zxhMatrixT(m_RowNum,m.m_ColNum);
	for(int i=0;i<temp.m_RowNum;i++)
	{
		for(int j=0;j<temp.m_ColNum;j++)
		{
			DataType tmp=0;
			for(int k=0;k<m_ColNum;k++)
			{
				tmp+=m_data[i*m_ColNum+k]*m.m_data[k*m.m_ColNum+j];
			}
			temp.m_data[i*temp.m_ColNum+j]=tmp;
		}
	}
	return temp;// zxh add, same problem as previous
}

template <class DataType>
zxhMatrixT<DataType> zxhMatrixT<DataType>::operator*( DataType num ) const
{
	zxhMatrixT temp(m_RowNum,m_ColNum);
	for(int i=0;i<m_RowNum*m_ColNum;i++)
		temp.m_data[i]=m_data[i]*num;
	return temp;// zxh add, same problem as previous
}

template <class DataType>
bool zxhMatrixT<DataType>::IsSqrMatrix()
{
	if (m_RowNum==m_ColNum)
	{
		return true;
	}
	return false;
}

template <class DataType>
bool zxhMatrixT<DataType>::IsSymmetric()
{
	if (this->IsSqrMatrix())
	{
		int N=m_RowNum;
		for (int i=0;i<N;i++)
		{
			for (int j=i;j<N;j++)
			{
				if (zxh::abs(m_data[i*m_ColNum+j]-m_data[j*m_ColNum+i])>0 )
				{
					return false;
				}
			}
		}
		return true;
	}
	return false;
}

template <class DataType>
bool zxhMatrixT<DataType>::FindPosMax( int &row, int &col )
{
	DataType max=0;
	int count=0;
	for (int i=0;i<m_RowNum;i++)
	{
		for (int j=0;j<m_ColNum;j++)
		{
			if (zxh::abs(m_data[i*m_ColNum+j])>=max)
			{
				max=zxh::abs(m_data[i*m_ColNum+j]);
				row=i;
				col=j;
				count++;
			}
		}
	}
	if (count<=0)
	{
		return false;
	}
	return true;
}


template <class DataType>
void zxhMatrixT<DataType>::CopyData( DataType Copy[] ) const
{
	// zxh add need to do a 
	if( Copy!=0)
	{
		for (int i=0;i<m_RowNum*m_ColNum;i++)
		{
			Copy[i]=m_data[i];
		}
	}
		
}

template <class DataType>
void zxhMatrixT<DataType>::Display()
{
	if (0==m_RowNum || 0==m_ColNum)
	{
		std::cout<<"This matrix is empty!";
	} 
	else
	{
		for (int i=0;i<m_RowNum;i++)
		{
			for (int j=0;j<m_ColNum;j++)
			{
				std::cout<<m_data[i*m_ColNum+j]<<'\t';
			}
			std::cout<<'\n';
		}
	}
	std::cout<<'\n';
}

template <class DataType>
int zxhMatrixT<DataType>::EigenvalueVectorRealSymmetryJacobi( zxhMatrixT& a, zxhMatrixT& v, DataType eps, int jt )
{
	int i,j, p, q, l(1), mrank;
	DataType fm,cn,sn,omega,x,y,d;

	if(!a.IsSymmetric())	//不是对称阵
		return(0);					

	mrank = a.GetRowNum();		// 矩阵阶数

	for(i=0; i<mrank; i++)
	{
		v.SetValue(1.0,i,i);
		for(j=0; j<mrank; j++)
			if(i!=j) v.SetValue(0.0,i,j);
	}
	while(1)
	{ 
		fm=0.0;
		for(i=1; i<mrank; i++)
		{
			for(j=0; j<i; j++)
			{
				d=zxh::abs(a.GetValue(i,j));
				if((i!=j)&&(d>fm))
				{
					fm=d; 
					p=i;
					q=j;
				}
			}
		}
		if(fm<eps)  return(1);
		if(l>jt)  return(-1);
		l=l+1;
		x=-a.GetValue(p,q);
		y=(a.GetValue(q,q)-a.GetValue(p,p))/2.0;
		omega=x/sqrt(x*x+y*y);
		if(y<0.0) omega=-omega;
		sn=1.0+sqrt(1.0-omega*omega);
		sn=omega/sqrt(2.0*sn);
		cn=sqrt(1.0-sn*sn);
		fm=a.GetValue(p,p);
		a.SetValue(fm*cn*cn+a.GetValue(q,q)*sn*sn+a.GetValue(p,q)*omega,p,p);
		a.SetValue(fm*sn*sn+a.GetValue(q,q)*cn*cn-a.GetValue(p,q)*omega,q,q);
		a.SetValue(0.0,p,q); 
		a.SetValue(0.0,q,p);
		for(j=0; j<mrank; j++)
		{
			if((j!=p)&&(j!=q))
			{
				fm=a.GetValue(p,j);
				a.SetValue(fm*cn+a.GetValue(q,j)*sn,p,j);
				a.SetValue(-fm*sn+a.GetValue(q,j)*cn,q,j);
			}
		}
		for(i=0; i<mrank; i++)
		{
			if((i!=p)&&(i!=q))
			{ 
				fm=a.GetValue(i,p);
				a.SetValue(fm*cn+a.GetValue(i,q)*sn,i,p);
				a.SetValue(-fm*sn+a.GetValue(i,q)*cn,i,q);
			}
		}
		for(i=0; i<mrank; i++)
		{
			fm=v.GetValue(i,p);
			v.SetValue(fm*cn+v.GetValue(i,q)*sn,i,p);
			v.SetValue(-fm*sn+v.GetValue(i,q)*cn,i,q);
		}
	}
	return(1);
};



typedef zxhMatrixT<float> zxhMatrixF;
typedef zxhMatrixT<float> zxhMatrix;


#endif



