/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 ZHUANG, Xia Hai
  Date:      From 2010-11
  Version:	 V 2.1

=========================================================================*/

#ifndef zxhPointSet_h
#define zxhPointSet_h
#ifdef HAS_VTK 
	#include "vtkPolyData.h"
	#include "vtkPolyDataReader.h"  
	#include "vtkPolyDataWriter.h"
	#include "vtkSmartPointer.h"
	#include "vtkPoints.h"
	#include "vtkUnstructuredGridReader.h"
	#include "vtkUnstructuredGridWriter.h"
	#include "vtkUnstructuredGrid.h"
	#include "vtkSmartPointer.h"
	#include "vtkPolyLine.h" 
#endif

#include "zxh.h"
#define ZXHPOINTSET_COORDINATETYPE float
#include "zxhPoint.h"

///
/// \class zxhPointSet
/// \brief:    default save/read type is nii.gz
/// \ingroup zxhImageDataT
/// 
class zxhPointSet
{
public:
	///
	zxhPointSet();

	///
	virtual ~zxhPointSet();

	/// new an image and clone
	virtual zxhPointSet * CloneTo(zxhPointSet * & pRet) const;
	///
	virtual void SetDimension( int i)							{m_iDimension=i;};
	///
	virtual int GetDimension( void ) const						{return m_iDimension;};
	///
	virtual int GetNumOfPoints( void ) const					{return m_iNumOfPoints;};

	///
	virtual void ReleaseMemoryOfPointSet() ; 
	///
	virtual bool AllocateMemoryForPointSet( int num, int dim ) ; 

	/// 
	virtual bool SetPointValues( int index, const float p[] );
	virtual bool SetPointValues( int index, float px, float py, float pz=0, float pt=0 );
	/// 3D
	virtual bool SetPointValues( int index, zxhPointF &p ); 
	///
	virtual bool GetPointValues( int index, float p[] ) const;
	virtual bool GetPointValues( int index, int p[] ) const;
	virtual bool GetPointValues( int index,  float &px, float &py, float &pz ) const;

	
	/// 
	virtual bool SetPointSetFromFile( std::string s )  ;

	///
	virtual bool SavePointSetToFile( std::string s )  const;

	/// 
	virtual bool SetPointSetFromNiiImage( std::string s )  ; 
	///
	virtual bool SavePointSetToNiiImage( std::string s )  const ; 
	/// 
	virtual bool SetPointSetFromTxt( std::string s )  ; 
	/// 3 float for each line, for 3D point set 
	virtual bool SavePointSetToTxt( std::string s )  const ; 

	/// s.pst            zxhtodo
	//virtual bool SetPointSetFromPST( std::string s )  ; 
	/// dimension ; num_of_points ; a point per line
	//virtual bool SavePointSetToPST( std::string s )  const ; 
	 
	/* vtk type:	STRUCTURED_POINTS
					STRUCTURED_GRID
					UNSTRUCTURED_GRID
					POLYDATA
					RECTILINEAR_GRID
					FIELD  */

	///
	virtual bool SetPointSetFromVtkPolyDataFile( std::string s )  ; 
	///
	virtual bool SavePointSetToVtkPolyDataFile( std::string s )  const ; 
	///
	virtual bool SetPointSetFromVtkUnstructuredGridFile( std::string s )  ; 
	///
	virtual bool SavePointSetToVtkUnstructuredGridFile( std::string s )  const ;  

protected:
	/// default 3
	int m_iDimension ;  
	/// 
	int m_iNumOfPoints ;
	///
	ZXHPOINTSET_COORDINATETYPE *m_pListPoints[ZXH_ImageDimensionMax] ;
	/// 0: unset, 1 vtk_polydata, 2 vtk_unstructuredgrid, 3 txt, 4 nii
	int m_iReadInType ; 
	///
	bool PointAccessValid( int index ) const 
	{
		if( m_iDimension<1 || m_iNumOfPoints<1 || index<0 || index>m_iNumOfPoints-1 )
			return false ; 
		return true ;
	}
} ;

 
#endif //zxhPointSet_h

