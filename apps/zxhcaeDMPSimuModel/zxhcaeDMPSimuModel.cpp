

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
void Help()
{
	
	//char *chFileName = argv[1] ; //"F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset00/vessel0/reference.vtk";
	//char *chMovedVecFilePath = argv[2] ; //"F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset00/vessel0/Simu_reference/20mmDeform/5PointsDeform";
	//int SimuNUM = atoi(argv[3]) ; // 30;//simulate number.
	//const int VecNUM = atoi(argv[4]) ; // 5;//number of deforming points.
	//const int RAD = atoi(argv[5]) ; // =20;

	std::cout<<" An simple example for siumulate reference line, target.vtk and source.vtk, result save as res: \n" ;
	std::cout<<" zxhcaeDMPSimulation centerlinefilename deformedcenterlinefilenamePre numberofsimulation num_point radius\n";
	std::cout<<"OR \n";
	//std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -hierarchy 3 -bending 0.0031\n\n";

	//std::cout<<"  <-test/target img.>     (test or target image)\n" ;
	//std::cout<<"  <-ref/source img.>      (reference or source image)\n" ;
	//std::cout<<"  <-o savename>           (string for saving transformed -ref/-source image, file names prefix)\n" ;
	//std::cout<<"  <-maskt img.>           (mask image on test image, use -maskr on ref image) \n" ;
	//std::cout<<"  USE -ffd: zxhsemi0 to fast init and get the .FFD for setting -pre, and then set -ffd\n" ;
	//std::cout<<"  <-ffd fx fy fz>         (spacing of free form deformation, FFD, multi-input related to -Reg)\n" ;
	//std::cout<<"  <-pre 0 s>              (pre set transformation field)\n";
	//std::cout<<"  <-sub fx fy fz [ft]>    ([3 3 3], sampling spacing; positive: millimeters interval; negative: pixels interval)\n";
	//std::cout<<"  <-Reg i>                ([1] number of multi-level registrations)\n" ;
	//std::cout<<"  OR USE -hierarchy, simple and not need to set -ffd,-sub,-Reg:\n" ; 
	//std::cout<<"  <-hierarchy n>          ([3] number of multi-level FFD registration, where\n";
	//std::cout<<"                           the first level of -ffd spacing is one forth of test image extent, and halve in subsequence level\n" ; 
	//std::cout<<"                           the final level of -sub sampling is pixel size of the test image\n" ; 	 
	//std::cout<<"\n";
	//std::cout<<"  <-bending f..f>         ([0.001]weighting for bending energy, recommend f=0.001~0.01)\n" ;
	//std::cout<<"  OPTIONS of spatially encoding scheme\n"  ;// similarity computation, default normalized mutual information\n" ; 
	//std::cout<<"  <-semiradius f...f>     (radius of local region, default set to twice ffd spacing)\n";
	////std::cout<<"  <-semiwidth f...f>      (or -semisize, width/size of local region, default 80mm)\n";
	////std::cout<<"  <-semikernel s>         ([0], -1='Gaussian', 0='ZeroBSpline', 3='BSpline')\n";

	std::cout<<"\n" ;

} 
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

void WriteModCA2Txt(char *chFileName,const vector<double>m_sMolPontIntsFinl)
	{
	ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = m_sMolPontIntsFinl.size();
	float fmodints;	
	for (int i = 0; i < nPointNum; i++)
	{
	 fmodints= m_sMolPontIntsFinl[i];
     WriteFileTxt <<right<<fixed<<setfill('0')<<setprecision(4)<< fmodints<<'\n';

	}	
	
	}
bool bInsideImg(int *to,const int ImgNewRawSize[4])
{ bool binside=true;
	for(int i=0;i<3;i++)
	{
		if ((int)to[i]<0)
			{
				to[i]=0;binside=false;
		}
		if ((int)to[i]>ImgNewRawSize[i]-1)
		{
			to[i]=ImgNewRawSize[i]-1;binside=false;
		}
	}
	return binside;
}
short SetModPontInts(PointCordTypeDef TempPoint,const int nImW[4],const short *sImgData,const zxhImageInfo *pImgInfo)//based on one point
{
	
	float fTempP[3]={0};
	fTempP[0]=TempPoint.x;
	fTempP[1]=TempPoint.y;
	fTempP[2]=TempPoint.z;
	pImgInfo->WorldToImage(fTempP);
	int itoTempP[3]={0};
	itoTempP[0]=(int)(fTempP[0]+0.5);
	itoTempP[1]=(int)(fTempP[1]+0.5);
	itoTempP[2]=(int)(fTempP[2]+0.5);
	bInsideImg(itoTempP,nImW);
	
	int m=itoTempP[2] * nImW[1] * nImW[0] + itoTempP[1] * nImW[0] + itoTempP[0];
	return sImgData[itoTempP[2] * nImW[1] * nImW[0] + itoTempP[1] * nImW[0] + itoTempP[0]];
	

}
double ModifiedSetModPontInts(PointCordTypeDef TempPoint,const int nImW[4],const short *sImgData,const zxhImageInfo *pImgInfo)//based on neighbour points.
{
	

	int nSearchRange = 8;
	int nCount = 0;
	double nSum = 0;
	float fCord[3];
	fCord[0] = TempPoint.x;
	fCord[1] = TempPoint.y;
	fCord[2] = TempPoint.z;
	pImgInfo->WorldToImage(fCord);
	int nMax[3]={0};
	int nMin[3]={0};
	nMin[0] = (int)(fCord[0]+0.5) - nSearchRange;
	nMax[0] = (int)(fCord[0]+0.5) + nSearchRange;

	nMin[1]= (int)(fCord[1]+0.5) - nSearchRange;
	nMax[1] = (int)(fCord[1]+0.5) + nSearchRange;

	nMin[2] = (int)(fCord[2]+0.5) - nSearchRange;
	nMax[2] = (int)(fCord[2]+0.5) + nSearchRange;
	bInsideImg(nMin,nImW);
	bInsideImg(nMax,nImW);
	// calculate mean
	for (int iz = nMin[2]; iz < nMax[2]; iz++)
	{
		for (int iy =nMin[1]; iy < nMax[1]; iy++)
		{
			for (int ix = nMin[0]; ix < nMax[0]; ix++)
			{
				if (sImgData[iz * nImW[1] * nImW[0] + iy * nImW[0] + ix] > 10)
				{
					nCount++;
					nSum += (double)sImgData[iz * nImW[1] * nImW[0] + iy *  nImW[0] + ix];
				}
			}
		}
	}

	return nSum / (nCount + 0.00000000001);

}
bool ResetModModPontInts(const vector<double>m_sMolPontInts,vector<double> &m_sMolPontIntsFinal,const float m_meandis_of_coronarymodel)
{
	if (m_sMolPontIntsFinal.size()>0)m_sMolPontIntsFinal.clear();
	float fsearchrange=5;
	int isearchrange=(int)(fsearchrange/m_meandis_of_coronarymodel+0.5);
	for(int i=0;i<m_sMolPontInts.size();i++)
	{ 
		int iF=i-isearchrange;
		int iB=i+isearchrange;
		if (iF<0)iF=0;
		if (iB>m_sMolPontInts.size()-1)iB=m_sMolPontInts.size()-1;
		long temp=0;
	    int iCount=0;
		for (int j=iF;j<=iB;j++)
		{
			temp=temp+m_sMolPontInts[j];
		   iCount++;
		}
		temp=temp/iCount;
		m_sMolPontIntsFinal.push_back(temp);
	}
	return true;
}
float ArcDist(int StartPosi,int EndPosi,vector<PointCordTypeDef> vPathPoints)
{
	float fFPoint[3]={0},fBPoint[3]={0};
	float fLength=0;
	if (StartPosi==EndPosi)
		fLength=0;
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
			fLength=fLength+zxh::VectorOP_Distance(fFPoint,fBPoint,3);
		}
	}
	return fLength;
}
float CalcMeandist(const vector<PointCordTypeDef> vPathPoints)//Calculate the mean distance density of model line.
{
	float fF[3]={0},fB[3]={0};
	float Dist=0;
	for(int i=1;i<vPathPoints.size();i++)
	{
		fF[0]=vPathPoints[i-1].x;
		fF[1]=vPathPoints[i-1].y;
		fF[2]=vPathPoints[i-1].z;
		fB[0]=vPathPoints[i].x;
		fB[1]=vPathPoints[i].y;
		fB[2]=vPathPoints[i].z;
		Dist=Dist+zxh::VectorOP_Distance(fF, fB, 3 ) ;  
	}
	return Dist/(vPathPoints.size()-1);
}
bool GetTangentOnPosij(int j,vector<PointCordTypeDef> vPathPoints,int NUM,PointCordTypeDef &tempPointVec)
{
	PointCordTypeDef FjPoint=vPathPoints[j-NUM];
	PointCordTypeDef BjPoint=vPathPoints[j+NUM];
	float Vector[3]={0,0,0};
	Vector[0]=BjPoint.x-FjPoint.x;
	Vector[1]=BjPoint.y-FjPoint.y;
	Vector[2]=BjPoint.z-FjPoint.z;
	zxh::VectorOP_Normalise(Vector,3);
	tempPointVec.x=Vector[0];
	tempPointVec.y=Vector[1];
	tempPointVec.z=Vector[2];
	return true;
}
int main(int argc, char *argv[])
{
	
	//**Simu-Reference**//
	/*if( argc < 7 )
	{
		cerr << "Usage: " << endl;
		cerr << "miiFindMaxPath	vtkFile(*)(.vtk)	FileNumber" << endl;
		return -1;
	} */
	// get file path and header
	//char *chModelFileName = argv[1] ; //"F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/reference.vtk";
	//char *chModelRawFileName=argv[2];//
	//char *chMovedVecFilePath = argv[3] ; //"F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/20mmDeform/5PointsDeform";
	//int SimuNUM = atoi(argv[4]) ; // 30;//simulate number.
	//int VecNUM = atoi(argv[5]) ; // 5;//number of deforming points.
	//int RAD = atoi(argv[6]) ; // =20
	//char *chModelFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset03/vessel2/reference.vtk";
	////char *chModelRawFileName="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";
	//char *chModelRawFileName="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image03.nii";
	//char *chMovedVecFilePath = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/SimuModel/dataset03/vessel2/Simu_reference";
 //   int SimuNUM = 30;//simulate number.
	//int VecNUM=5;//number of deforming points.
	//int RAD_20 =20;
	//int RAD_15 =15;
	//int RAD_5=5;
	//int RAD_2=2;
	//float DIST=1;
	//float fFPoint[3]={0},fBPoint[3]={0};
	//float TempVector[3]={0,0,0};
	//float fLength=0,fLengthDeform=0,fSectionLength=0;
	//int *DeformPointPosi = new int[VecNUM] ,DPOrder=0;
	//
	//vector<PointCordTypeDef> vPathPoints,vPathMovedPoints,vTangentVector;
	//vector<double> m_sMolPontInts;
	//double dTempInts;
	//// read vtk file
	//fstream DetectFile;
	//DetectFile.open(chModelFileName,ios::in);
	//if(DetectFile)
	//	{
	//		ReadVtk(chModelFileName, vPathPoints);
	//	}
	//DetectFile.close();

	//float m_meandis_of_coronarymodel=CalcMeandist(vPathPoints);//Calculate the mean distance density of model line.
	//int PNum=(int)DIST/m_meandis_of_coronarymodel;
	//zxhImageDataT<short>imgReadMolRaw; //Change by JDQ
	//if( zxh::OpenImage( &imgReadMolRaw, chModelRawFileName ) == false )
	//{
	//	std::cerr << "Raw image(nifti-file) is not found!"; 
	//	return -1;
	//}
	//const short *sImData = imgReadMolRaw.GetImageData();
	//int nImW[4];
	//imgReadMolRaw.GetImageSize(nImW[0], nImW[1], nImW[2], nImW[3]);
	//PointCordTypeDef StartPoint,BStartPoint,EndPoint,FEndPoint,TempPoint;
	//PointCordTypeDef *DeformPoint = new PointCordTypeDef[VecNUM];
	//PointCordTypeDef *fRandVector=new PointCordTypeDef[VecNUM];
	//StartPoint=vPathPoints[0];
	////get the tangent vector of startpoint
	//StartPoint=vPathPoints[0];
	//BStartPoint=vPathPoints[PNum];
	//float TangentStartPointVector[3]={0,0,0};
	//TangentStartPointVector[0]=BStartPoint.x-StartPoint.x;
	//TangentStartPointVector[1]=BStartPoint.y-StartPoint.y;
	//TangentStartPointVector[2]=BStartPoint.z-StartPoint.z;
	//zxh::VectorOP_Normalise(TangentStartPointVector,3);
	//PointCordTypeDef tempPoint;
	//tempPoint.x=TangentStartPointVector[0];
	//tempPoint.y=TangentStartPointVector[1];
	//tempPoint.z=TangentStartPointVector[2];
	//vTangentVector.push_back(tempPoint);	
	////calculate the whole arclength of the line
	//fLength=ArcDist(0,vPathPoints.size()-1,vPathPoints);
	////find the Deform Position of the Deform Points;
	//DeformPointPosi[0]=0;
	//DeformPoint[0]=vPathPoints[0];
	//DeformPointPosi[VecNUM-1]=vPathPoints.size()-1;
	//DeformPoint[VecNUM-1]=vPathPoints[vPathPoints.size()-1];
	//for (int j=0;j<=vPathPoints.size()-1;j++)
	//{
	//	fLengthDeform=ArcDist(DeformPointPosi[DPOrder],j,vPathPoints);
	//	if (fLengthDeform>=(float)fLength/(VecNUM-1))
	//	{
	//		DPOrder++;
	//		DeformPoint[DPOrder]=vPathPoints[j];
	//		DeformPointPosi[DPOrder]=j;
	//		PointCordTypeDef tempPointVec1;
	//		GetTangentOnPosij(j,vPathPoints,PNum,tempPointVec1);
	//		vTangentVector.push_back(tempPointVec1);	
	//		if (DPOrder>=VecNUM-2) break;
	//	}
	//}
	////get the tangent vector of endpoint
	//EndPoint=vPathPoints[vPathPoints.size()-1];
	//FEndPoint=vPathPoints[vPathPoints.size()-1-PNum];
	//float TangentEndPointVector[3]={0,0,0};
	//TangentEndPointVector[0]=EndPoint.x-FEndPoint.x;
	//TangentEndPointVector[1]=EndPoint.y-FEndPoint.y;
	//TangentEndPointVector[2]=EndPoint.z-FEndPoint.z;
	//zxh::VectorOP_Normalise(TangentEndPointVector,3);
	//tempPoint.x=TangentEndPointVector[0];
	//tempPoint.y=TangentEndPointVector[1];
	//tempPoint.z=TangentEndPointVector[2];
	//vTangentVector.push_back(tempPoint);	

	////get the simulation 
	//for (int i=0;i<SimuNUM;i++)
	//{
	//	if (vPathMovedPoints.size()>0)vPathMovedPoints.clear();
	//	if (m_sMolPontInts.size()>0)m_sMolPontInts.clear();
	//	float fFi=0,fTheta=0,R=0;
	//	srand((unsigned)time(NULL));
	//	for (int i=0;i<VecNUM;i++)//get the random defomation of deformpoints.
	//	{   
	//		
	//		fFi=(float)M_PI* (rand()%360);
	//		fTheta=(float)M_PI*(rand()%180);
	//		
	//		if (i==VecNUM-1)
	//		{
	//			R=(float)(rand()%RAD_15);
	//			if(R==0)R=1;
	//			TempVector[0]=R*sin(fTheta)*cos(fFi);
	//			TempVector[1]=R*sin(fTheta)*sin(fFi);
	//			TempVector[2]=R*cos(fTheta);
	//			
	//			if(zxh::VectorOP_Cosine( TempVector, TangentEndPointVector, 3 )<0)
	//			{
	//				float TempVectorNorm[3];//y axsis
	//				float TempVectorDire[3];//x axis
	//				zxh::VectorOP_Normalise(TempVector,3);
	//				zxh::VectorOP_CrossProduct3D(TempVector,TangentEndPointVector,TempVectorDire);
	//				zxh::VectorOP_CrossProduct3D(TempVectorDire,TangentEndPointVector,TempVectorNorm);//calculate y axis
	//				zxh::VectorOP_Normalise(TempVectorNorm,3);
	//				float cosinevalue=zxh::VectorOP_Cosine( TempVector, TangentEndPointVector, 3 );
	//				float sinvalue=sqrt(1-cosinevalue*cosinevalue);
	//				zxh::VectorOP_Multiple(TangentEndPointVector,-R*cosinevalue,TempVectorDire,3);
	//				zxh::VectorOP_Multiple(TempVectorNorm,R*sinvalue,TempVectorNorm,3);
	//				zxh::VectorOP_Add(TempVectorDire,TempVectorDire,TempVector,3);//TempVector is the final vector while cosiinevalue<0
	//			}
	//			fRandVector[i].x=TempVector[0];
	//			fRandVector[i].y=TempVector[1];
	//			fRandVector[i].z=TempVector[2];

	//		}
	//		else if (i>0)
	//		{
	//			R=(float)RAD_20*(rand() % 21)/20.0;
	//			if(R==0)R=1;
	//			TempVector[0]=R*sin(fTheta)*cos(fFi);
	//			TempVector[1]=R*sin(fTheta)*sin(fFi);
	//			TempVector[2]=R*cos(fTheta);
	//			zxh::VectorOP_Normalise(TempVector,3);
	//			float TempVec[3];
	//			float TempVectorNorm[3];//y axsis
	//			float TempVectorDire[3];//x axis
	//			TempVec[0]=vTangentVector[i].x;
	//			TempVec[1]=vTangentVector[i].y;
	//			TempVec[2]=vTangentVector[i].z;
	//			zxh::VectorOP_Normalise(TempVector,3);
	//			zxh::VectorOP_CrossProduct3D(TempVector,TempVec,TempVectorDire);
	//			zxh::VectorOP_CrossProduct3D(TempVectorDire,TempVec,TempVectorNorm);//calculate y axis
	//			zxh::VectorOP_Normalise(TempVectorNorm,3);//calculate the vertical vector of tangenvector
	//			fRandVector[i].x=R*TempVectorNorm[0];
	//			fRandVector[i].y=R*TempVectorNorm[1];
	//			fRandVector[i].z=R*TempVectorNorm[2];
	//		}
	//		else //first model point deforms 2mm.
	//		{
	//			R=(float)RAD_2*(rand() % 20)/19.0;
	//			fRandVector[i].x=R*sin(fTheta)*cos(fFi);
	//			fRandVector[i].y=R*sin(fTheta)*sin(fFi);
	//			fRandVector[i].z=R*cos(fTheta);
	//		}
	//	}
	//	int StartDeformPosi=DeformPointPosi[0],EndDeformPosi=DeformPointPosi[1]-1;
	//	//deform the points
	//	for(int i=0;i<VecNUM-1;i++)
	//	{
	//		float fFLength=0,fBLength=0;

	//		float fMVec[3]={0};
	//		fSectionLength=ArcDist(StartDeformPosi,EndDeformPosi+1,vPathPoints);
	//		for(int j=StartDeformPosi;j<=EndDeformPosi;j++)
	//		{
	//			fFLength=ArcDist(StartDeformPosi,j,vPathPoints);
	//			fBLength=fSectionLength-fFLength;
	//			fMVec[0]=fBLength/fSectionLength*fRandVector[i].x+fFLength/fSectionLength*fRandVector[i+1].x;
	//			fMVec[1]=fBLength/fSectionLength*fRandVector[i].y+fFLength/fSectionLength*fRandVector[i+1].y;
	//			fMVec[2]=fBLength/fSectionLength*fRandVector[i].z+fFLength/fSectionLength*fRandVector[i+1].z;
	//			dTempInts=SetModPontInts(vPathPoints[j],nImW,sImData,imgReadMolRaw.GetImageInfo());
	//			m_sMolPontInts.push_back(dTempInts);
	//			TempPoint.x =vPathPoints[j].x+fMVec[0];
	//			TempPoint.y =vPathPoints[j].y+fMVec[1];
	//			TempPoint.z =vPathPoints[j].z+fMVec[2];
	//			vPathMovedPoints.push_back(TempPoint);

	//		}
 //           StartDeformPosi=EndDeformPosi+1;
	//		if ((i+2)>(VecNUM-1)) break;
	//		else EndDeformPosi=DeformPointPosi[i+2]-1;
	//	}
	//	dTempInts=SetModPontInts(EndPoint,nImW,sImData,imgReadMolRaw.GetImageInfo());
	//	m_sMolPontInts.push_back(dTempInts);;
	//	TempPoint.x =EndPoint.x+fRandVector[VecNUM-1].x;
	//    TempPoint.y =EndPoint.y+fRandVector[VecNUM-1].y;
	//	TempPoint.z =EndPoint.z+fRandVector[VecNUM-1].z;
	//	vPathMovedPoints.push_back(TempPoint);
	//	/*vector<double>m_sMolPontIntsFinl;
	//	ResetModModPontInts(m_sMolPontInts,m_sMolPontIntsFinl,m_meandis_of_coronarymodel);*/
	//	// save simulation files
	//	char chTemp[25];
	//	_itoa_s(i, chTemp, 10);
	//	int nFileLen = strlen(chMovedVecFilePath) + strlen("/Simulation")+strlen(chTemp) + strlen(".vtk") + 1;
	//	char *chFileName = (char*)malloc(nFileLen);
	//	strcpy(chFileName, chMovedVecFilePath);
	//	strcat(chFileName, "/Simulation");
	//	strcat(chFileName, chTemp);
	//	strcat(chFileName, ".vtk");
	//	WriteVtk(vPathMovedPoints, chFileName);
	//	int nFileLen1 = strlen(chMovedVecFilePath) + strlen("/Simuintensity")+strlen(chTemp) + strlen(".txt") + 1;
	//	char *chFileName1 = (char*)malloc(nFileLen1);
	//	strcpy(chFileName1, chMovedVecFilePath);
	//	strcat(chFileName1, "/Simuintensity");
	//	strcat(chFileName1, chTemp);
	//	strcat(chFileName1, ".txt");
	//	WriteModCA2Txt(chFileName1,m_sMolPontInts);
	//	//output the deforming information
	//	string strMCFilename = string(chMovedVecFilePath);
	//	strMCFilename = strMCFilename + "/" + string("SimuModelDf.txt");
	//	ofstream WriteFileMC;
	//	WriteFileMC.open(strMCFilename,ios::app);
	//	for (int i=0;i<VecNUM;i++)
	//	{
	//	WriteFileMC << chTemp<<"-"<<i<<":(" << fRandVector[i].x << ", " << fRandVector[i].y << ", " << fRandVector[i].z << ")" << endl;
	//	}
	//	WriteFileMC.close();
	//}
	//
	//delete []DeformPointPosi ; 
	//delete []DeformPoint;
	//delete []fRandVector;
	//cout << " Simulation Down!" << endl;

	//return 1;

	//**Simu-Reference**//

	//**Simu-Meancenterline

//Release
	if( argc < 2 )
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeDMPSimuModel	ModelFileNamevtkFile(*)(.vtk)	MovedVecFilePath" << endl;
		return -1;
	} 	
	char *chModelFileName = argv[1] ; //"F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image00_meancenterline_v0.vtk";
	char *chMovedVecFilePath=argv[2];//"F:/Coronary_0/Coronary_Niessen/Simu_meancenterline/image00/vessel0";

//Release
	// get file path and header
	//char *chModelFileName = argv[1] ; //"F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/reference.vtk";
	//char *chModelRawFileName=argv[2];//
	//char *chMovedVecFilePath = argv[3] ; //"F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/20mmDeform/5PointsDeform";
	//int SimuNUM = atoi(argv[4]) ; // 30;//simulate number.
	//int VecNUM = atoi(argv[5]) ; // 5;//number of deforming points.
	//int RAD = atoi(argv[6]) ; // =20
	
	



//char *chModelFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset04/vessel3/reference.vtk";
//char *chMovedVecFilePath = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/SimuModel/dataset04/vessel3/Simu_reference2";
    int SimuNUM =20;//simulate number.
	int VecNUM=5;//number of deforming points.
	int RAD_20 =20;
	int RAD_15 =15;
	int RAD_5=5;
	int RAD_2=2;
	float DIST=1;
	float fFPoint[3]={0},fBPoint[3]={0};
	float TempVector[3]={0,0,0};
	float fLength=0,fLengthDeform=0,fSectionLength=0;
	int *DeformPointPosi = new int[VecNUM] ,DPOrder=0;
	
	vector<PointCordTypeDef> vPathPoints,vPathMovedPoints,vTangentVector;
	vector<double> m_sMolPontInts;
	double dTempInts;
	// read vtk file
	fstream DetectFile;
	DetectFile.open(chModelFileName,ios::in);
	if(DetectFile)
		{
			ReadVtk(chModelFileName, vPathPoints);
		}
	DetectFile.close();

	float m_meandis_of_coronarymodel=CalcMeandist(vPathPoints);//Calculate the mean distance density of model line.
	int PNum=(int)DIST/m_meandis_of_coronarymodel;
	zxhImageDataT<short>imgReadMolRaw; //Change by JDQ
	int nImW[4];
	imgReadMolRaw.GetImageSize(nImW[0], nImW[1], nImW[2], nImW[3]);
	PointCordTypeDef StartPoint,BStartPoint,EndPoint,FEndPoint,TempPoint;
	PointCordTypeDef *DeformPoint = new PointCordTypeDef[VecNUM];
	PointCordTypeDef *fRandVector=new PointCordTypeDef[VecNUM];
	StartPoint=vPathPoints[0];
	//get the tangent vector of startpoint
	StartPoint=vPathPoints[0];
	BStartPoint=vPathPoints[PNum];
	float TangentStartPointVector[3]={0,0,0};
	TangentStartPointVector[0]=BStartPoint.x-StartPoint.x;
	TangentStartPointVector[1]=BStartPoint.y-StartPoint.y;
	TangentStartPointVector[2]=BStartPoint.z-StartPoint.z;
	zxh::VectorOP_Normalise(TangentStartPointVector,3);
	PointCordTypeDef tempPoint;
	tempPoint.x=TangentStartPointVector[0];
	tempPoint.y=TangentStartPointVector[1];
	tempPoint.z=TangentStartPointVector[2];
	vTangentVector.push_back(tempPoint);	
	//calculate the whole arclength of the line
	fLength=ArcDist(0,vPathPoints.size()-1,vPathPoints);
	//find the Deform Position of the Deform Points;
	DeformPointPosi[0]=0;
	DeformPoint[0]=vPathPoints[0];
	DeformPointPosi[VecNUM-1]=vPathPoints.size()-1;
	DeformPoint[VecNUM-1]=vPathPoints[vPathPoints.size()-1];
	for (int j=0;j<=vPathPoints.size()-1;j++)
	{
		fLengthDeform=ArcDist(DeformPointPosi[DPOrder],j,vPathPoints);
		if (fLengthDeform>=(float)fLength/(VecNUM-1))
		{
			DPOrder++;
			DeformPoint[DPOrder]=vPathPoints[j];
			DeformPointPosi[DPOrder]=j;
			PointCordTypeDef tempPointVec1;
			GetTangentOnPosij(j,vPathPoints,PNum,tempPointVec1);
			vTangentVector.push_back(tempPointVec1);	
			if (DPOrder>=VecNUM-2) break;
		}
	}
	//get the tangent vector of endpoint
	EndPoint=vPathPoints[vPathPoints.size()-1];
	FEndPoint=vPathPoints[vPathPoints.size()-1-PNum];
	float TangentEndPointVector[3]={0,0,0};
	TangentEndPointVector[0]=EndPoint.x-FEndPoint.x;
	TangentEndPointVector[1]=EndPoint.y-FEndPoint.y;
	TangentEndPointVector[2]=EndPoint.z-FEndPoint.z;
	zxh::VectorOP_Normalise(TangentEndPointVector,3);
	tempPoint.x=TangentEndPointVector[0];
	tempPoint.y=TangentEndPointVector[1];
	tempPoint.z=TangentEndPointVector[2];
	vTangentVector.push_back(tempPoint);	

	//get the simulation 
	for (int i=0;i<SimuNUM;i++)
	{
		if (vPathMovedPoints.size()>0)vPathMovedPoints.clear();
		if (m_sMolPontInts.size()>0)m_sMolPontInts.clear();
		float fFi=0,fTheta=0,R=0;
		srand((unsigned)time(NULL));
		for (int i=0;i<VecNUM;i++)//get the random defomation of deformpoints.
		{   
			
			fFi=(float)M_PI* (rand()%360);
			fTheta=(float)M_PI*(rand()%180);
			
			if (i==VecNUM-1)
			{
				R=(float)(rand()%RAD_15);//R=(float)RAD_2*(rand() % 20)/19.0;
				if(R==0)R=1;
				TempVector[0]=R*sin(fTheta)*cos(fFi);
				TempVector[1]=R*sin(fTheta)*sin(fFi);
				TempVector[2]=R*cos(fTheta);
				
				if(zxh::VectorOP_Cosine( TempVector, TangentEndPointVector, 3 )<0)
				{
					float TempVectorNorm[3];//y axsis
					float TempVectorDire[3];//x axis
					zxh::VectorOP_Normalise(TempVector,3);
					zxh::VectorOP_CrossProduct3D(TempVector,TangentEndPointVector,TempVectorDire);
					zxh::VectorOP_CrossProduct3D(TempVectorDire,TangentEndPointVector,TempVectorNorm);//calculate y axis
					zxh::VectorOP_Normalise(TempVectorNorm,3);
					float cosinevalue=zxh::VectorOP_Cosine( TempVector, TangentEndPointVector, 3 );
					float sinvalue=sqrt(1-cosinevalue*cosinevalue);
					zxh::VectorOP_Multiple(TangentEndPointVector,-R*cosinevalue,TempVectorDire,3);
					zxh::VectorOP_Multiple(TempVectorNorm,R*sinvalue,TempVectorNorm,3);
					zxh::VectorOP_Add(TempVectorDire,TempVectorDire,TempVector,3);//TempVector is the final vector while cosiinevalue<0
				}
				fRandVector[i].x=TempVector[0];
				fRandVector[i].y=TempVector[1];
				fRandVector[i].z=TempVector[2];

			}
			else if (i>0)
			{
				R=(float)(rand()%RAD_15);//R=(float)RAD_2*(rand() % 20)/19.0;
				if(R==0)R=1;
				TempVector[0]=R*sin(fTheta)*cos(fFi);
				TempVector[1]=R*sin(fTheta)*sin(fFi);
				TempVector[2]=R*cos(fTheta);
				zxh::VectorOP_Normalise(TempVector,3);
				float TempVec[3];
				float TempVectorNorm[3];//y axsis
				float TempVectorDire[3];//x axis
				TempVec[0]=vTangentVector[i].x;
				TempVec[1]=vTangentVector[i].y;
				TempVec[2]=vTangentVector[i].z;
				zxh::VectorOP_Normalise(TempVector,3);
				zxh::VectorOP_CrossProduct3D(TempVector,TempVec,TempVectorDire);
				zxh::VectorOP_CrossProduct3D(TempVectorDire,TempVec,TempVectorNorm);//calculate y axis
				zxh::VectorOP_Normalise(TempVectorNorm,3);//calculate the vertical vector of tangenvector
				fRandVector[i].x=R*TempVectorNorm[0];
				fRandVector[i].y=R*TempVectorNorm[1];
				fRandVector[i].z=R*TempVectorNorm[2];
			}
			else //first model point deforms 2mm.
			{
				R=(float)RAD_2*(rand() % 20)/19.0;
				fRandVector[i].x=R*sin(fTheta)*cos(fFi);
				fRandVector[i].y=R*sin(fTheta)*sin(fFi);
				fRandVector[i].z=R*cos(fTheta);
			}
		}
		int StartDeformPosi=DeformPointPosi[0],EndDeformPosi=DeformPointPosi[1]-1;
		//deform the points
		for(int i=0;i<VecNUM-1;i++)
		{
			float fFLength=0,fBLength=0;

			float fMVec[3]={0};
			fSectionLength=ArcDist(StartDeformPosi,EndDeformPosi+1,vPathPoints);
			for(int j=StartDeformPosi;j<=EndDeformPosi;j++)
			{
				fFLength=ArcDist(StartDeformPosi,j,vPathPoints);
				fBLength=fSectionLength-fFLength;
				fMVec[0]=fBLength/fSectionLength*fRandVector[i].x+fFLength/fSectionLength*fRandVector[i+1].x;
				fMVec[1]=fBLength/fSectionLength*fRandVector[i].y+fFLength/fSectionLength*fRandVector[i+1].y;
				fMVec[2]=fBLength/fSectionLength*fRandVector[i].z+fFLength/fSectionLength*fRandVector[i+1].z;
				TempPoint.x =vPathPoints[j].x+fMVec[0];
				TempPoint.y =vPathPoints[j].y+fMVec[1];
				TempPoint.z =vPathPoints[j].z+fMVec[2];
				vPathMovedPoints.push_back(TempPoint);

			}
            StartDeformPosi=EndDeformPosi+1;
			if ((i+2)>(VecNUM-1)) break;
			else EndDeformPosi=DeformPointPosi[i+2]-1;
		}
		TempPoint.x =EndPoint.x+fRandVector[VecNUM-1].x;
	    TempPoint.y =EndPoint.y+fRandVector[VecNUM-1].y;
		TempPoint.z =EndPoint.z+fRandVector[VecNUM-1].z;
		vPathMovedPoints.push_back(TempPoint);
		/*vector<double>m_sMolPontIntsFinl;
		ResetModModPontInts(m_sMolPontInts,m_sMolPontIntsFinl,m_meandis_of_coronarymodel);*/
		// save simulation files
		char chTemp[25];
		_itoa_s(i, chTemp, 10);
		int nFileLen = strlen(chMovedVecFilePath) + strlen("/Simulation")+strlen(chTemp) + strlen(".vtk") + 1;
		char *chFileName = (char*)malloc(nFileLen);
		strcpy(chFileName, chMovedVecFilePath);
		strcat(chFileName, "/Simulation");
		strcat(chFileName, chTemp);
		strcat(chFileName, ".vtk");
		WriteVtk(vPathMovedPoints, chFileName);
		//output the deforming information
		string strMCFilename = string(chMovedVecFilePath);
		strMCFilename = strMCFilename + "/" + string("SimuModelDf.txt");
		ofstream WriteFileMC;
		WriteFileMC.open(strMCFilename,ios::app);
		for (int i=0;i<VecNUM;i++)
		{
		WriteFileMC << chTemp<<"-"<<i<<":(" << fRandVector[i].x << ", " << fRandVector[i].y << ", " << fRandVector[i].z << ")" << endl;
		}
		WriteFileMC.close();
	}
	
	delete []DeformPointPosi ; 
	delete []DeformPoint;
	delete []fRandVector;
	cout << " Simulation Down!" << endl;

	return 1;

	//**Simu-Centerline
}

 