

 /*=========================================================================

  Program:   ZXH CoronaryArteryExtraction Software
  Author:	 Dengqiang Jia
  Module:    $RCSfle: zxhcaeDMP.cpp    $
  Language:  C++
  Date:      $Date: From  2011-03 $
  Version:   $Revision: 2.2.1 $

=========================================================================*/
/// \brief
/// Spatially encoded mutual information + free-form deformation registration
/// Linear
/// gradient descent
/// save all registration info into zxhRegistrationStruct StructRegistration which set to optimiser
/// For concatenated transforms all (for FFD and LocallyAffines in optimization regularisation step)
///                                 First update zxhRegistrationStruct::m_pConcatenatedTransformsByRegridding(spacing 1mm)
///                                 then transform all images in ref_space using zxhRegistrationStruct.m_CopyRefXXXOrg
///                                 finally set current transform to identity
/// For preset transformation, unless same spacing FFD for the first Reg, otherwise would be treated as concatenation
///
 
//
//void Help()
//{
//	std::cout<<" An simple example for registration of images, target.nii.gz and source.nii.gz, result save as res: \n" ;
//	std::cout<<" zxhsemi0 -test target.nii.gz -ref source.nii.gz -o result0 -ffd 20 20 20 -bending 0.001\n";
//	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -pre 0 result0.FFD -ffd 20 20 20 -ffd 10 10 10 -Reg 2 -sub 2 2 2 -sub 1 1 1 -bending 0.001\n";
//	std::cout<<"OR \n";
//	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -hierarchy 3 -bending 0.0031\n\n";
//
//	std::cout<<"  <-test/target img.>     (test or target image)\n" ;
//	std::cout<<"  <-ref/source img.>      (reference or source image)\n" ;
//	std::cout<<"  <-o savename>           (string for saving transformed -ref/-source image, file names prefix)\n" ;
//	std::cout<<"  <-maskt img.>           (mask image on test image, use -maskr on ref image) \n" ;
//	std::cout<<"  USE -ffd: zxhsemi0 to fast init and get the .FFD for setting -pre, and then set -ffd\n" ;
//	std::cout<<"  <-ffd fx fy fz>         (spacing of free form deformation, FFD, multi-input related to -Reg)\n" ;
//	std::cout<<"  <-pre 0 s>              (pre set transformation field)\n";
//	std::cout<<"  <-sub fx fy fz [ft]>    ([3 3 3], sampling spacing; positive: millimeters interval; negative: pixels interval)\n";
//	std::cout<<"  <-Reg i>                ([1] number of multi-level registrations)\n" ;
//	std::cout<<"  OR USE -hierarchy, simple and not need to set -ffd,-sub,-Reg:\n" ; 
//	std::cout<<"  <-hierarchy n>          ([3] number of multi-level FFD registration, where\n";
//	std::cout<<"                           the first level of -ffd spacing is one forth of test image extent, and halve in subsequence level\n" ; 
//	std::cout<<"                           the final level of -sub sampling is pixel size of the test image\n" ; 	 
//	std::cout<<"\n";
//	std::cout<<"  <-bending f..f>         ([0.001]weighting for bending energy, recommend f=0.001~0.01)\n" ;
//	std::cout<<"  OPTIONS of spatially encoding scheme\n"  ;// similarity computation, default normalized mutual information\n" ; 
//	std::cout<<"  <-semiradius f...f>     (radius of local region, default set to twice ffd spacing)\n";
//	//std::cout<<"  <-semiwidth f...f>      (or -semisize, width/size of local region, default 80mm)\n";
//	//std::cout<<"  <-semikernel s>         ([0], -1='Gaussian', 0='ZeroBSpline', 3='BSpline')\n";
//
//	std::cout<<"\n" ;
//
//} 
//void HELP()
//{ 
//	Help() ;
//	std::cout<<"------------------------------------------------------\n" ;
//	std::cout<<"  OPTIONS of gradient optimization computation; use setting in previous -Reg when un-set \n" ; 
//}  
//int main(int argc, char* argv[])
//{
//	zxh::echo_zxh(argc,argv);
//	if( argc == 1 )
//	{
//		std::cout<<"  zxhsemi [options]        \n";
//		Help();
//		return 0 ;
//	}
//	if( argc == 2 && strcmp( argv[1],"-H" )==0 )
//	{
//		std::cout<<"  zxhsemi [options]        \n";
//		HELP();
//		return -1 ;
//	}
//	if( glbVerboseOutput>0 )
//	{
//		std::cout<<"\n * \n * zxhsemi, version of 2011-03  \n * \n";
//		zxh::echo_arguments( argc, argv ) ;
//	} 
//	return zxhsemi_main(argc,argv);
//} 
//int zxhsemi_main(int argc, char* argv[])
//{}
//


#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"


#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"


#include "vtkUnstructuredGrid.h"

// for read
#include "vtkUnstructuredGridReader.h"

// for write
#include "vtkUnstructuredGridWriter.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include<stdlib.h>

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"

#define M_PI 3.14159265358979323846
using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
}PointCordTypeDef;

typedef struct
{
	float fMean;
	float fStdVar;
	float fMaxErr;
}ErrTypeDef;

void ReadVtk(char *chFileName, vector<PointCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find VTK-file!" << endl;
		return;
	}
	if (!PointCord.empty())
	{
		PointCord.clear();
	}

	vtkSmartPointer<vtkUnstructuredGridReader> iVtkReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	iVtkReader->SetFileName( chFileName );
	iVtkReader->Update();

	vtkSmartPointer<vtkUnstructuredGrid> iGridRead = iVtkReader->GetOutput();

	int nPointNum = iGridRead->GetMaxCellSize();

	double dCord[3];
	PointCordTypeDef strctTempPoint;

	for (int i = 0; i < nPointNum; i++)
	{
		iGridRead->GetPoint(i, dCord);
		strctTempPoint.x = dCord[0];
		strctTempPoint.y = dCord[1];
		strctTempPoint.z = dCord[2];
		PointCord.push_back(strctTempPoint);
	}
}

void WriteVtk(vector< PointCordTypeDef > PointCord, char* chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	/*	int nPointNum = PointCord.size();*/

	float fImgPixel[3];
	int nPointNum = PointCord.size();

	for (int i = 0; i < nPointNum; i++)
	{
		iPoints->InsertNextPoint(PointCord[i].x, PointCord[i].y, PointCord[i].z);
	}	

	vtkSmartPointer<vtkPolyLine> iLine = vtkSmartPointer<vtkPolyLine>::New();
	iLine->GetPointIds()->SetNumberOfIds(nPointNum);
	for (int i = 0; i < nPointNum; i++)
	{
		iLine->GetPointIds()->SetId(i, i);
	}

	vtkSmartPointer<vtkUnstructuredGrid> iGrid = vtkUnstructuredGrid::New();
	iGrid->Allocate(1, 1);	
	iGrid->InsertNextCell(iLine->GetCellType(), iLine->GetPointIds());
	iGrid->SetPoints(iPoints);

	vtkSmartPointer<vtkUnstructuredGridWriter> iVtkWriter = vtkUnstructuredGridWriter::New();
	iVtkWriter->SetInput(iGrid);
	iVtkWriter->SetFileName(chFileName);
	iVtkWriter->Write();
}
float ArcDist(int StartPosi,int EndPosi,vector<PointCordTypeDef> vPathPoints)
{
	float fFPoint[3]={0},fBPoint[3]={0};
	float fLength=0;
	if (StartPosi=EndPosi)fLength=0;
	else
	{
	for (int j=StartPosi;j<EndPosi;j++)
	{
		fFPoint[0]=vPathPoints[j].x;
		fFPoint[1]=vPathPoints[j].y;
		fFPoint[2]=vPathPoints[j].z;
		fBPoint[0]=vPathPoints[j+1].x;
		fBPoint[1]=vPathPoints[j+1].y;
		fBPoint[2]=vPathPoints[j+1].z;
		fLength=fLength+zxh::VectorOP_Distance(fFPoint,fFPoint,3);
	}
	}
	return fLength;
}

int main(int argc, char *argv[])
{
	/*if( argc < 2 )
	{
		cerr << "Usage: " << endl;
		cerr << "miiFindMaxPath	vtkFile(*)(.vtk)	FileNumber" << endl;
		return -1;
	}
*/
	// get file path and header
	char *chFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset00/vessel1/reference.vtk";
	char *chMovedFilePath = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset00/vessel1";
	int SimuNUM=30;

	
	vector<PointCordTypeDef> vPathPoints,vPathMovedPoints;
	// read vtk file
	fstream DetectFile;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
		}
	DetectFile.close();
	delete[] chFileName;
	float fFPoint[3]={0},fBPoint[3]={0};
	float fLength=0,fLengthDeform=0;
	
	PointCordTypeDef StartPoint,DeformPoint,EndPoint,TempPoint,fRandVector[3];
	StartPoint=vPathPoints[0];
	EndPoint=vPathPoints[vPathPoints.size()-1];
	int DeformPointPosi;
	fLength=ArcDist(0,vPathPoints.size()-1,vPathPoints);
	for (int j=0;j<=vPathPoints.size()-1;j++)
	{
		fLengthDeform=ArcDist(0,j,vPathPoints);
		if (fLengthDeform>=0.5*fLength)
		{
			DeformPoint=vPathPoints[j];
			DeformPointPosi=j;
			break;
		}
	}
	for (int i=0;i<SimuNUM;i++)
	{
		float fFi,fTheta,R;
		for (int i=0;i<3;i++)
		{
			fFi=2*M_PI*rand();
			fTheta=M_PI*rand();
			R=5*rand();
			fRandVector[i].x=R*sin(fTheta)*cos(fFi);
			fRandVector[i].y=R*sin(fTheta)*sin(fFi);
			fRandVector[i].z=R*cos(fTheta);
		}
		float StartDeformPosi=0,EndDeformPosi=DeformPointPosi-1;
		for(int i=0;i<2;i++)
		{
			float fFLength=0,fBLength=0;

			float fMVec[3]={0};

			for(int j=StartDeformPosi;j<=EndDeformPosi;j++)
			{
				fFLength=ArcDist(0,j,vPathPoints);
				fBLength=fLength-fFLength;
				fMVec[0]=fBLength/fLength*fRandVector[i].x+fFLength/fLength*fRandVector[i+1].x;
				fMVec[1]=fBLength/fLength*fRandVector[i].y+fFLength/fLength*fRandVector[i+1].y;
				fMVec[2]=fBLength/fLength*fRandVector[i].z+fFLength/fLength*fRandVector[i+1].z;
				TempPoint.x =vPathPoints[j].x+fMVec[0];
				TempPoint.y =vPathPoints[j].y+fMVec[1];
				TempPoint.z =vPathPoints[j].z+fMVec[2];
				vPathMovedPoints.push_back(TempPoint);

			}
			StartDeformPosi=DeformPointPosi,EndDeformPosi=vPathPoints.size()-1;
		}
		// save max-path
		char chTemp[25];
		_itoa_s(i, chTemp, 10);
		int nFileLen = strlen(chMovedFilePath) + strlen("/Simulation")+strlen(chTemp) + strlen(".vtk") + 1;
		char *chFileName = (char*)malloc(nFileLen);
		strcpy(chFileName, chMovedFilePath);
		strcat(chFileName, "/Simulation");
		strcat(chFileName, chTemp);
		strcat(chFileName, ".vtk");
		WriteVtk(vPathMovedPoints, chFileName);
	}
	
	cout << " Simulation Down!" << endl;

	return 1;
}

 