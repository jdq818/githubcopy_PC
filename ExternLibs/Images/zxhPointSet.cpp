#include "zxhPointSet.h"
#include "zxhImageGipl.h"

zxhPointSet::zxhPointSet()
{ 
	m_iDimension = 3 ; 
	m_iNumOfPoints = 0 ;
	for( int id=0; id<ZXH_ImageDimensionMax; ++id )
		m_pListPoints[id] = 0 ;
	m_iReadInType = 0 ; 
}


zxhPointSet::~zxhPointSet()
{
	if( m_iNumOfPoints>0 )
	{
		for( int id=0; id<ZXH_ImageDimensionMax; ++id )
			if( m_pListPoints[id]!=0 )
				delete [] m_pListPoints[id] ;
	}
};

zxhPointSet * zxhPointSet::CloneTo(zxhPointSet * & pRet) const
{
	if( pRet==0 ) pRet = new zxhPointSet(); 
	pRet->m_iDimension = this->m_iDimension ; 
	pRet->m_iNumOfPoints = this->m_iNumOfPoints ;
	for( int id=0; id<ZXH_ImageDimensionMax; ++id )
	{
		if( m_pListPoints[id]!=0 )
		{
			pRet->m_pListPoints[id] = new ZXHPOINTSET_COORDINATETYPE [m_iNumOfPoints];
			memcpy((void*)pRet->m_pListPoints[id], (void*)m_pListPoints[id], m_iNumOfPoints*sizeof(ZXHPOINTSET_COORDINATETYPE));
		}
	}
	pRet->m_iReadInType = this->m_iReadInType ;
	return pRet ;
} 



void zxhPointSet::ReleaseMemoryOfPointSet() 
{ 
	if( m_iNumOfPoints>0 )
	{
		for( int id=0; id<ZXH_ImageDimensionMax; ++id )
			if( m_pListPoints[id]!=0 )
			{
				delete [] m_pListPoints[id] ;
				m_pListPoints[id] = 0 ;
			}
		m_iNumOfPoints = 0 ;
	}
}

bool zxhPointSet::AllocateMemoryForPointSet( int num, int dim ) 
{ 
	if( num<1 || dim<1 ) 
		return false ;
	if( m_iNumOfPoints!=num || m_iDimension!=dim ) 
	{
		for( int id=0; id<ZXH_ImageDimensionMax; ++id )
			if( m_pListPoints[id]!=0 )
			{
				delete [] m_pListPoints[id] ;
				m_pListPoints[id] = 0 ;
			}
		for( int id=0; id<dim; ++id ) 
		{
			m_pListPoints[id] = new ZXHPOINTSET_COORDINATETYPE [num];
		}
	}
	m_iNumOfPoints = num ;
	m_iDimension = dim ;
	for( int id=0; id<m_iDimension; ++id ) 
		memset((void*)m_pListPoints[id], 0, m_iNumOfPoints*sizeof(ZXHPOINTSET_COORDINATETYPE));
	return true ;
}
 

	///
bool zxhPointSet::SetPointValues( int index, zxhPointF &p )
{
	if( m_iDimension!=3 ) return false ;
	if( PointAccessValid(index) == false ) return false ;
	m_pListPoints[0][index] = p.x;  
	m_pListPoints[1][index] = p.y;  
	m_pListPoints[2][index] = p.y;  
	return true ;
} 
bool zxhPointSet::SetPointValues( int index, const float p[] )
{
	if( PointAccessValid(index) == false ) return false ;
	for( int id=0; id<m_iDimension; ++id )
		m_pListPoints[id][index] = p[id] ;  
	return true ;
}
bool zxhPointSet::SetPointValues( int index, float px, float py, float pz, float pt )
{
	if( PointAccessValid(index) == false ) return false ;
	float p[] = {px,py,pz,pt} ; 
	for( int id=0; id<m_iDimension; ++id )
		m_pListPoints[id][index] = p[id] ;  
	return true ;
}
	///
bool zxhPointSet::GetPointValues( int index, float p[] ) const
{
	if( PointAccessValid(index) == false ) return false ;
	for( int id=0; id<m_iDimension; ++id )
		p[id] = m_pListPoints[id][index] ;  
	return true ;
}
bool zxhPointSet::GetPointValues( int index, int p[] ) const
{
	if( PointAccessValid(index) == false ) return false ;
	for( int id=0; id<m_iDimension; ++id )
		p[id] = zxh::round(m_pListPoints[id][index]) ;  
	return true ;
}
bool zxhPointSet::GetPointValues( int index,  float &px, float &py, float &pz ) const
{
	if( m_iDimension!=3 ) return false ;
	if( PointAccessValid(index) == false ) return false ;
	px = m_pListPoints[0][index] ;
	py = m_pListPoints[1][index] ;
	pz = m_pListPoints[2][index] ;
	return true ;
}
bool zxhPointSet::SetPointSetFromFile( std::string s )  
{
	m_iReadInType = 0 ; 
	std::string ext = zxh::GetExtension( s ) ; 
	zxh::case_lower( ext ) ; 
	if( strcmp( ext.c_str(), "txt" )==0 )
	{
		return SetPointSetFromTxt( s ) ; 
	}
	else if( strcmp( ext.c_str(), "vtk" )==0 )
	{
		if( SetPointSetFromVtkPolyDataFile( s ) == true )
			return true ;
		else return SetPointSetFromVtkUnstructuredGridFile( s ) ;  
	}
	if( strcmp( ext.c_str(), "pst" )==0 )
	{
		std::cerr<<"error: have not implement PST type yet!\n" ; 
		return false ;
	}
	// default, nii.gz type
	return SetPointSetFromNiiImage(s) ;
}

	///
bool zxhPointSet::SavePointSetToFile( std::string s )  const 
{
	if( m_iNumOfPoints<1 || m_iDimension<1 ) return false ;
	std::string ext = zxh::GetExtension( s ) ; 
	zxh::case_lower( ext ) ; 
	if( strcmp( ext.c_str(), "txt" )==0 )
	{
		return SavePointSetToTxt( s ) ; 
	}
	else if( strcmp( ext.c_str(), "vtk" )==0 )
	{
		switch( m_iReadInType )
		{
		case 1: return SavePointSetToVtkPolyDataFile( s ) ; break ;
		case 2: return SavePointSetToVtkUnstructuredGridFile( s ) ; break ;
		default : return SavePointSetToVtkPolyDataFile( s ) ; break ; //default polydata
		}
	}
	if( strcmp( ext.c_str(), "pst" )==0 )
	{
		std::cerr<<"error: have not implement PST type yet!\n" ; 
		return false ;
	}
	// default, nii.gz type
	return SavePointSetToNiiImage(s) ;
}
 
bool zxhPointSet::SetPointSetFromVtkPolyDataFile( std::string s ) 
{
#ifdef HAS_VTK 
	vtkPolyDataReader * reader = vtkPolyDataReader::New() ; 
	reader->SetFileName( s.c_str() ) ; 
	reader->Update(); 
	//unsigned long ecode = reader->GetErrorCode()  ; 
	vtkPolyData * pPolyData = reader->GetOutput(); 
	pPolyData->Register(pPolyData);
	reader->Delete();
		
	if( pPolyData == 0 ) 
		return false;
	AllocateMemoryForPointSet(  pPolyData->GetNumberOfPoints(), 3 ) ;  
	if( m_iNumOfPoints<=0 ) return false ;

	for( int ip=0; ip<m_iNumOfPoints; ++ip )
		for( int id=0; id<m_iDimension; ++id )
			m_pListPoints[id][ip] = pPolyData->GetPoint(ip)[id] ; 
		 
	pPolyData->Delete();
	m_iReadInType = 1 ;
	return true ;
#else 
	return false ;
#endif
}

bool zxhPointSet::SavePointSetToVtkPolyDataFile( std::string s ) const
{ 
#ifdef HAS_VTK 
	if( m_iDimension > 3 )
	{
		std::cerr<<"error: vtkPolyData does not take points with dimension more than 3\n";
		return false ;
	} 
	
	// point 
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();
	//vtkPoints* iPoints = vtkPoints::New() ;  
	for( int ip=0; ip<m_iNumOfPoints; ++ip )
	{ 
		float point[3] = {0,0,0} ; 
		for( int id=0; id<m_iDimension; ++id )
			point[id] = float(m_pListPoints[id][ip]);
		iPoints->InsertNextPoint( point[0],point[1],point[2] ) ; 
	} 

	/// old vtk (<5.0)
	/*
	vtkPolyData *polyData = vtkPolyData::New();
	polyData->SetPoints(iPoints);
	polyData->Update() ;

	vtkPolyDataWriter * writer = vtkPolyDataWriter::New() ; 
	writer->SetFileName( s.c_str() ) ; 
	writer->SetInput( polyData ) ;
	//writer->SetFileTypeToASCII();*/

	/// vtk 6.3, 2017-01-21 zxh
	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(iPoints);

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(s.c_str());
	writer->SetInputData(polyData);
   
	int iRet = writer->Write() ;
	return  (iRet==1);
#else 
	return false ;
#endif
}
 
	///
bool zxhPointSet::SetPointSetFromVtkUnstructuredGridFile( std::string s ) 
{ 
#ifdef HAS_VTK 
	vtkSmartPointer<vtkUnstructuredGridReader> iVtkReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	iVtkReader->SetFileName( s.c_str() );
	iVtkReader->Update();
	//unsigned long ecode = iVtkReader->GetErrorCode()  ; 
	vtkSmartPointer<vtkUnstructuredGrid> iGridRead = iVtkReader->GetOutput();
	if( iGridRead==0 ) return false ; 
	  
	AllocateMemoryForPointSet(  iGridRead->GetMaxCellSize(), 3 ) ; 
	if( m_iNumOfPoints<=0 ) return false ;

	for( int ip=0; ip<m_iNumOfPoints; ++ip )
	{
		ZXH_Floatvtk pnt[3] ; // for vtk default type d ouble
		iGridRead->GetPoint(ip, pnt); 
		m_pListPoints[0][ip] = pnt[0] ; 
		m_pListPoints[1][ip] = pnt[1] ; 
		m_pListPoints[2][ip] = pnt[2] ; 
		 
	}
	m_iReadInType = 2 ;
	return true ; 
#else 
	return false ;
#endif
}
	///
bool zxhPointSet::SavePointSetToVtkUnstructuredGridFile( std::string s )  const 
{
#ifdef HAS_VTK 
	if( m_iDimension!=3 || m_iNumOfPoints<1 )
		return false ;

	vtkSmartPointer<vtkUnstructuredGridWriter> iVtkWriter = vtkUnstructuredGridWriter::New();

	// point 
	vtkSmartPointer<vtkPoints> iPoints = vtkPoints::New() ; 
	for (int ip = 0; ip < m_iNumOfPoints; ip++)
	{ 
		iPoints->InsertNextPoint(m_pListPoints[0][ip], m_pListPoints[1][ip], m_pListPoints[2][ip] ) ;  
	}
	//structure line
	vtkSmartPointer<vtkPolyLine> iLine = vtkSmartPointer<vtkPolyLine>::New();
	iLine->GetPointIds()->SetNumberOfIds(m_iNumOfPoints);
	for (int ip = 0; ip < m_iNumOfPoints; ip++)
	{
		iLine->GetPointIds()->SetId(ip, ip);
	}
	//grid 
	vtkSmartPointer<vtkUnstructuredGrid> iGrid = vtkUnstructuredGrid::New();  
	iGrid->Allocate(1, 1);	
	iGrid->InsertNextCell(iLine->GetCellType(), iLine->GetPointIds());
	iGrid->SetPoints(iPoints); 

	iVtkWriter->SetInputData(iGrid);

	iVtkWriter->SetFileName(s.c_str()); 
	return (iVtkWriter->Write()==1);
#else 
	return false ;
#endif
} 

/// 
bool zxhPointSet::SetPointSetFromNiiImage( std::string s )  
{
	zxhImageDataF img ; 
	if( zxh::OpenImage( &img, s ) == false )
		return false ;
	const int * size = img.GetImageSize(); 
	 
	AllocateMemoryForPointSet( size[0], size[1] ) ; 
	for( int idy=0; idy<size[1]; ++idy )
	for( int ipx=0; ipx<size[0]; ipx++ )
	{
		m_pListPoints[idy][ipx] = img.GetPixelGreyscale( ipx, idy, 0, 0 ) ; 
	}
	m_iReadInType = 4 ;
	return true ;
}
	///
bool zxhPointSet::SavePointSetToNiiImage( std::string s )  const 
{ 
	zxhImageDataF img ;  
	zxhImageInfo imageinfo ;  
	int size[] = { m_iNumOfPoints, m_iDimension, 1, 1} ; 
	float spacing[] = {1,1,1,1} ; 
	img.NewImage( 3, size, spacing, &imageinfo ) ;   
	for( int idy=0; idy<size[1]; ++idy )
	for( int ipx=0; ipx<size[0]; ipx++ )
	{
		img.SetPixelByGreyscale( ipx, idy, 0, 0,  m_pListPoints[idy][ipx] ) ;  
	}
	return zxh::SaveImage( &img, s ) ;
}
	/// 
bool zxhPointSet::SetPointSetFromTxt( std::string s ) 
{ 
	std::ifstream instream(s.c_str());
	if (instream.fail()==true)
	{
		std::cerr<<"error: Open point set txt file "<<s<<" Failed!\n";
		return false;
	}
	instream.seekg(0,std::ios_base::beg);  
	char buffer[1024];
	std::string sLine,sContent,sComment;
	std::vector<float> vcox, vcoy, vcoz ;  
	float wco[3] ;
	while( instream.eof() == false )
	{
		do{
				instream.getline(buffer,1024);
				sLine=buffer;
				zxh::ParseStringLine(sContent,sComment,sLine);
				zxh::trim_both(sContent);
		}while(sContent.empty()&&instream.eof()==false&&instream.fail()==false);
		if(instream.eof()==true||instream.fail()==true)
			break ;
		std::istringstream strstream(sContent.c_str());
		for(int c=0;c<3;++c) 
			strstream>>wco[c]; 
		vcox.push_back( wco[0] ) ; 
		vcoy.push_back( wco[1] ) ; 
		vcoz.push_back( wco[2] ) ; 
	}
	instream.close();
	 
	if( vcox.size()<1 )
		return false ; 

	AllocateMemoryForPointSet( vcox.size(), 3 ) ;
	 
	for( int ip=0; ip<m_iNumOfPoints; ++ip )
	{ 
		m_pListPoints[0][ip] = vcox.at( ip ) ; 
		m_pListPoints[1][ip] = vcoy.at( ip ) ; 
		m_pListPoints[2][ip] = vcoz.at( ip ) ; 		 
	}
	m_iReadInType = 3 ; 
	return true ;  
}
	/// 3 float for each line, for 3D point set 
bool zxhPointSet::SavePointSetToTxt( std::string s )  const 
{
	if( m_iDimension!=3 )
		return false ;

	std::ofstream ofs ;
	ofs.open( s.c_str(), std::ofstream::out);
	if(  ofs.fail() )
	{
		std::cerr<< "error: open file name "<< s <<" fail! " ;
		return false ;
	} 
	char buffer[1024];
	if( m_iNumOfPoints>1 )
	{
		sprintf(buffer, "#number of 3D points: %d #\n", m_iNumOfPoints ) ;
		ofs<< buffer ;
	}
	for( int ip=0; ip<m_iNumOfPoints; ip++ )
	{ 
		sprintf(buffer, "%f \t %f \t %f\n", m_pListPoints[0][ip], m_pListPoints[1][ip], m_pListPoints[2][ip] ) ;
		ofs<< buffer ;
	}
	ofs.flush();
	ofs.close();
	return true ;
	 
}


