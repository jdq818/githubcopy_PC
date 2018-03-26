 #include "zxhROI.h"
 

zxhROI::zxhROI()
{ 
	m_iDimension = -1 ;
	for( int i=0; i<ZXH_ImageDimensionMax; ++i )
	{
		m_aiFrom[i] = m_aiTo[i] = 0 ;
		m_afFrom[i] = m_afTo[i] = 0.0f ;
	}
}


zxhROI::~zxhROI(){};

zxhROI * zxhROI::CloneTo(zxhROI * & pRet) const
{
	if( pRet==0 ) pRet = new zxhROI(); 
	pRet->m_iDimension = this->m_iDimension ; 
	for( int i=0; i<ZXH_ImageDimensionMax; ++i )
	{
		pRet->m_aiFrom[i] = this->m_aiFrom[i] ;
		pRet->m_aiTo[i] = this->m_aiTo[i] ;
		pRet->m_afFrom[i] = this->m_afFrom[i] ; 
		pRet->m_afTo[i] = this->m_afTo[i] ;
	}
	return pRet ;
}  
bool zxhROI::CorrectFromTo() 
{
	bool b=false ; 
	
	for( int i=0; i<ZXH_ImageDimensionMax; ++i )
	{
		if( m_aiFrom[i] > m_aiTo[i] )
		{
			b = true ;
			zxh::ExchangeTwoValues( m_aiFrom[i], m_aiTo[i] ) ; 
		};
		if( m_afFrom[i] > m_afTo[i] )
		{
			b = true ;
			zxh::ExchangeTwoValues( m_afFrom[i], m_afTo[i] ) ; 
		};
	}
	return b ; 
}

	///     
void zxhROI::SetRoiFromTo( const float from[], const float to[], int dim ) 
{
	m_iDimension = dim ;
	for( int i=0; i<dim; ++i )
	{
		m_afFrom[i] = from[i] ;
		m_aiFrom[i] = zxh::round( from[i] ) ; 
		m_afTo[i] = to[i] ; 
		m_aiTo[i] = zxh::round( to[i] ) ; 
	}
	CorrectFromTo() ;
}
	///
void zxhROI::SetRoiFromTo( const int from[], const int to[], int dim ) 
{
	m_iDimension = dim ;
	for( int i=0; i<dim; ++i )
	{
		m_afFrom[i] = from[i] ;
		m_aiFrom[i] = from[i] ; 
		m_afTo[i] = to[i] ; 
		m_aiTo[i] = to[i] ; 
	}
	CorrectFromTo() ;
} 
	///
void zxhROI::SetRoi3DFromTo( int xfrom, int yfrom, int zfrom, int xto, int yto, int zto ) 
{
	int from[] = {xfrom, yfrom, zfrom}, to[] = {xto, yto, zto} ; 
	m_iDimension = 3 ;
	for( int i=0; i<3; ++i )
	{
		m_afFrom[i] = from[i] ;
		m_aiFrom[i] = from[i] ; 
		m_afTo[i] = to[i] ; 
		m_aiTo[i] = to[i] ; 
	}
	CorrectFromTo() ;
}
	///
void zxhROI::SetRoi3DFromTo( float xfrom, float yfrom, float zfrom, float xto, float yto, float zto ) 
{
	float from[] = {xfrom, yfrom, zfrom}, to[] = {xto, yto, zto} ; 
	m_iDimension = 3 ;
	for( int i=0; i<3; ++i )
	{
		m_afFrom[i] = from[i] ;
		m_aiFrom[i] = zxh::round( from[i] ) ; 
		m_afTo[i] = to[i] ; 
		m_aiTo[i] = zxh::round( to[i] ) ; 
	}
	CorrectFromTo() ;
} 
	///
void zxhROI::SetRoiCenterRadius( const float center[], float radius, int dim )
{
	m_iDimension = dim ;
	for( int i=0; i<dim; ++i )
	{
		m_afFrom[i] = center[i]-radius ;
		m_aiFrom[i] = zxh::round( m_afFrom[i] ) ; 
		m_afTo[i] = center[i]+radius ; 
		m_aiTo[i] = zxh::round( m_afTo[i] ) ; 
	}
	CorrectFromTo() ;
} 
	///
void zxhROI::SetRoiCenterSideLength( const float center[], const float length[], int dim ) 
{
	m_iDimension = dim ;
	for( int i=0; i<dim; ++i )
	{
		m_afFrom[i] = center[i]-length[i] ;
		m_aiFrom[i] = zxh::round( m_afFrom[i] ) ; 
		m_afTo[i] = center[i]+length[i] ; 
		m_aiTo[i] = zxh::round( m_afTo[i] ) ; 
	}
	CorrectFromTo() ;
} 

	/// \return whether the ROI has been set
bool zxhROI::GetRoiFromTo( float from[], float to[] ) const 
{
	if( m_iDimension<=0 )
		return false ; 
	for( int i=0; i<m_iDimension; ++i )
	{
		from[i] = m_afFrom[i] ;
		to[i] = m_afTo[i] ;
	}
	return true ;
}
	///
bool zxhROI::GetRoiFromTo( int from[], int to[] ) const 
{
	if( m_iDimension<=0 )
		return false ; 
	for( int i=0; i<m_iDimension; ++i )
	{
		from[i] = m_aiFrom[i] ;
		to[i] = m_aiTo[i] ;
	}
	return true ;
}
	///
bool zxhROI::GetRoi3DFromTo( int& xfrom, int& yfrom, int& zfrom, int& xto, int& yto, int& zto  ) const 
{
	if( m_iDimension<=0 )
		return false ; 
	xfrom = m_aiFrom[0] ; 
	yfrom = m_aiFrom[1] ; 
	zfrom = m_aiFrom[2] ; 
	xto = m_aiTo[0] ; 
	yto = m_aiTo[1] ; 
	zto = m_aiTo[2] ;
	return true ;
} 
	///
bool zxhROI::GetRoi3DFromTo( float& xfrom, float& yfrom, float& zfrom, float& xto, float& yto, float& zto  ) const  
{
	if( m_iDimension<=0 )
		return false ; 
	xfrom = m_afFrom[0] ; 
	yfrom = m_afFrom[1] ; 
	zfrom = m_afFrom[2] ; 
	xto = m_afTo[0] ; 
	yto = m_afTo[1] ; 
	zto = m_afTo[2] ;
	return true ;
} 


