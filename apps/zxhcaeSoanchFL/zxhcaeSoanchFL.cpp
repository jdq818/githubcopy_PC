

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
#include<stdio.h>

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"
#include "miiMinHeap.h"

#define TAB_CHAR	9
#define M_PI 3.14159265358979323846
#define SPHERE_RADIUS 4
using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
}PointCordTypeDef;

typedef struct
{
	int num;
	float dist;
}LineInfo;
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
void WriteTxt(vector< PointCordTypeDef > PointCord, char* chFileName)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = PointCord.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = PointCord[i].x;
		fImgPixel[1] = PointCord[i].y;
		fImgPixel[2] = PointCord[i].z;
		WriteFileTxt <<right<<fixed<<setfill('0')<<setprecision(4) << -fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] <<'\n';

	}	
}
float CalcPathLen(vector< PointCordTypeDef > &PointCord)
{
	float fRefLineLen=0;
	for(int i=0;i<PointCord.size()-1;i++)
	{
		float fF[3]={PointCord[i+1].x,PointCord[i+1].y,PointCord[i+1].z};
		float fi[3]={PointCord[i].x,PointCord[i].y,PointCord[i].z};
		float fTdist=zxh::VectorOP_Distance(fF,fi,3);
		fRefLineLen=fRefLineLen+fTdist;

	}
	return fRefLineLen;
}
char *GetPathEXT(char *chFileName)
{  char path_buffer[_MAX_PATH];  
char drive[_MAX_DRIVE];  
char dir[_MAX_DIR];  
char fname[_MAX_FNAME];  
char ext[_MAX_EXT];  

_splitpath( chFileName, drive, dir, fname, ext );  

return ext;
}
bool ReadRadIoTxt(char *filename,float fRad[3],float fIo[3])
{
	string strNum;
	int nStart = 0, nEnd = 0;
	miiCNode<double, float> strctTempPoint;
	ifstream iFileIn(filename);
	if(!iFileIn)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << filename << endl;
		return false; //exit(1); // terminate with error
	}
	string strLine;
	if (!getline(iFileIn, strLine))
	{
		cout << "The format of Mean and StdDev file is wrong!" << endl;
		return false;
	}

	ifstream ifs(filename);

	char *buffer=new char[255];int ibuffersize=255;
	float s;
	std::string sComment,sContent,sLine;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,ibuffersize);
		sLine=buffer;
		zxh::ParseStringLine(sContent,sComment,sLine);
		zxh::trim_both(sContent);

		if(sContent.length()!=0)
		{
			std::istringstream istr;
			istr.str(sContent);

			istr>>fRad[0]>>fIo[0]>>fRad[1]>>fIo[1]>>fRad[2]>>fIo[2]>>fRad[3]>>fIo[3];
			//cout<<fRad[0]<<fIo[0]<<fRad[1]<<fIo[1]<<fRad[2]<<fIo[2]<<fRad[3]<<fIo[3];
		}
	}

	return true;
}
bool ReadPointTxt(char *filename,vector< miiCNode<double, float> > &cl)
{
	string strNum;
	int nStart = 0, nEnd = 0;
	miiCNode<double, float> strctTempPoint;
	ifstream iFileIn(filename);
	if(!iFileIn)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << filename << endl;
		return false; //exit(1); // terminate with error
	}
	string strLine;
	if (!getline(iFileIn, strLine))
	{
		cout << "The format of Mean and StdDev file is wrong!" << endl;
		return false;
	}
	//// read 1st number
	//nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	//strNum.assign(strLine, nStart, nEnd - nStart);
	//strctTempPoint.x =atof(strNum.c_str());
	////fMean=1024*fMean/1896;//Add by JDQ
	//// read 2nd number
	//nStart = nEnd + 1;
	//nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	//strNum.assign(strLine, nStart, nEnd - nStart);
	//strctTempPoint.y =atof(strNum.c_str());
	//// read 3rd number
	//nStart = nEnd + 1;
	//nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	//strNum.assign(strLine, nStart, nEnd - nStart);
	//strctTempPoint.z =atof(strNum.c_str());
	//cl.push_back(strctTempPoint);

	ifstream ifs(filename);

	char *buffer=new char[100];int ibuffersize=100;
	float s;
	std::string sComment,sContent,sLine;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,ibuffersize);
		sLine=buffer;
		zxh::ParseStringLine(sContent,sComment,sLine);
		zxh::trim_both(sContent);

		if(sContent.length()!=0)
		{
			std::istringstream istr;
			istr.str(sContent);
			float a,b,c;
			istr>>a>>b>>c;
			strctTempPoint.x=-a;
			strctTempPoint.y=-b;
			strctTempPoint.z=c;
		};
	}
	cl.push_back(strctTempPoint);
	return true;
}
bool SoanchMaxL(string strMCLineFolder,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM)
{

	int nMaxPointNum = 0, nMaxPointIdx = -1;
	//read the results as MCLine.vtk from different model

	for (int j=0;j<nModinFileMaxNUM;j++)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
			if (vPathPoints.size() > nMaxPointNum)
			{
				nMaxPointNum = vPathPoints.size();
				nMaxPointIdx = j;
			}			
		}

		DetectFile.close();
		delete[] chFileName;		
	}

	if (nMaxPointIdx < 0)
	{
		cout << "Cannot find any file!" << endl;
		return -1;
	}
	//read points into vector


	char *chFileName=new char[1024];
	fstream DetectFile;
	string strtempstr;
	stringstream sstr;
	string strj;
	sstr.clear();
	sstr<<nMaxPointIdx;
	sstr>>strj;
	strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr.clear();
	sstr<<strtempstr;
	sstr>>chFileName;
	ReadVtk(chFileName, vPathPoints);
	delete[] chFileName;
}
int FindPosiInPathF(vector<PointCordTypeDef> &vPathPointsF,miiCNode<double, float> &REFfourPoints)
{
	float fmindist=10000;
	int minindex=-1;
	float fRefP[3]={REFfourPoints.x,REFfourPoints.y,REFfourPoints.z};
	for (int j=0;j<vPathPointsF.size();j++)
	{
		float fcLP[3]={vPathPointsF[j].x,vPathPointsF[j].y,vPathPointsF[j].z};

		float tempdist=zxh::VectorOP_Distance(fRefP,fcLP,3);
		if(tempdist<fmindist)
		{
			fmindist=tempdist;
			minindex=j;
		}
	}
	return minindex;
}

bool CalcSATFPN(vector<PointCordTypeDef> &vPathPointsF,vector<PointCordTypeDef> &vPathPointsT, miiCNode<double, float> &REFfourPoints,float fRad,double &sumD,float &avgD,int &T,int &F)
{
	//calculate TPM FP
	int nPosi=FindPosiInPathF(vPathPointsF,REFfourPoints);
	float flen=fRad*4;
	float fLineLen=CalcPathLen(vPathPointsF);
	float finten=fLineLen/(vPathPointsF.size()-1);
	int nRange=flen/finten;
	int nS=max(nPosi-nRange,0);
	int nE=(int)zxh::minf(nPosi+nRange,vPathPointsF.size()-1);
	for(int j=nS;j<nE;j++)
	{
		float fRefLP[3]={vPathPointsF[j].x,vPathPointsF[j].y,vPathPointsF[j].z};
		float mindist=100000;
		int minindex=0;
		for (int j=0;j<vPathPointsT.size();j++)
		{
			float fcLP[3]={vPathPointsT[j].x,vPathPointsT[j].y,vPathPointsT[j].z};

			float tempdist=zxh::VectorOP_Distance(fRefLP,fcLP,3);
			if(tempdist<mindist)
			{
				mindist=tempdist;
				minindex=j;


			}
		}
		if(mindist<fRad)
		{
			T++;
		}
		else
		{
			F++;
		}
		sumD=sumD+mindist;


	}
	avgD=sumD/(nE-nS);

	return true;
}
float CalcDistFromBSLineToLine(vector<PointCordTypeDef>vResResPathPoints,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef > vResRefPathPoints,int nModinFileMaxNUM,float fRad[3],float fIo[3],int RNUM)
{
	float distmm=0;

	float fPointWorld[3]={0,0,0};

	int MNUM=4-RNUM;
	//refline length
	float fRefLineLen=CalcPathLen(vResRefPathPoints);
	float finten=fRefLineLen/(vResRefPathPoints.size()-1);
	//method line length
	float fMLineLen=CalcPathLen(vResResPathPoints);
	float fMLinten=fMLineLen/(vResResPathPoints.size()-1);
	int nForRefPontNUM=vREFfourPoints.size()-MNUM;
	vector<double> vsumResD,vsumRefD,vavgResD,vavgRefD;
	vector<int> vTPR,vFN,vTPM,vFP;
	vsumResD.clear();
	vsumRefD.clear();
	vTPR.clear();
	vFN.clear();
	vTPM.clear();
	vFP.clear();
	//calculate the reference point 
	for(int i=0;i<nForRefPontNUM;i++)
	{
		float fPathPointWorld[3]={0,0,0};
		fPointWorld[0]=vREFfourPoints[i].x;
		fPointWorld[1]=vREFfourPoints[i].y;
		fPointWorld[2]=vREFfourPoints[i].z;
		double sumRefD=0;float avgRefD=0;
		int TPR=0;int FN=0;
		double sumResD=0;float avgResD=0;
		int TPM=0;int FP=0;
		CalcSATFPN(vResRefPathPoints,vResResPathPoints,vREFfourPoints[i],fRad[i],sumRefD,avgRefD,TPR,FN);
		vsumRefD.push_back(sumRefD);
		vavgRefD.push_back(avgRefD);
		vTPR.push_back(TPR);
		vFN.push_back(FN);
		//calculate TPR FN
		CalcSATFPN(vResResPathPoints,vResRefPathPoints,vREFfourPoints[i],fRad[i],sumResD,avgResD,TPM,FP);
		vsumResD.push_back(sumResD);
		vavgResD.push_back(avgResD);
		vTPM.push_back(TPM);
		vFP.push_back(FP);

	}
	//add the T P F N
	int T, N;
	for(int i=0;i<vTPR.size()-1;i++)
	{
		T=T+vTPR[i];
	}
	for(int i=0;i<vTPM.size()-1;i++)
	{
		T=T+vTPM[i];
	}
	for(int i=0;i<vFN.size()-1;i++)
	{
		N=N+vFN[i];
	}
	for(int i=0;i<vFP.size()-1;i++)
	{
		N=N+vFP[i];
	}

	float fSumavDist=0;
	for(int i=0;i<vavgResD.size()-1;i++)
	{
		fSumavDist=fSumavDist+vavgResD[i];
	}
	for(int i=0;i<vavgRefD.size()-1;i++)
	{
		fSumavDist=fSumavDist+vavgRefD[i];
	}	

	double fSumDist=0;
	for(int i=0;i<vsumResD.size()-1;i++)
	{
		fSumDist=fSumDist+vsumResD[i];
	}
	for(int i=0;i<vsumRefD.size()-1;i++)
	{
		fSumDist=fSumDist+vsumRefD[i];
	}	
	float favgmindist=fSumDist/(T+N);
	return favgmindist;
}
float CalcMinDistFromPoinToLine(vector<PointCordTypeDef>&vPathPoints,vector< miiCNode<double, float> > vREFfourPoints,int RNUM)
{
	float distmm=0;
	float sumdistmm=0;
	float fPointWorld[3]={0,0,0};
	float fPathPointWorld[3]={0,0,0};
	int MNUM=4-RNUM;
	for(int i=0;i<vREFfourPoints.size()-MNUM;i++)
	{
		fPointWorld[0]=vREFfourPoints[i].x;
		fPointWorld[1]=vREFfourPoints[i].y;
		fPointWorld[2]=vREFfourPoints[i].z;
		float mindist=100000;
		int minindex=0;
		for (int j=0;j<vPathPoints.size();j++)
		{
			fPathPointWorld[0]=vPathPoints[j].x;
			fPathPointWorld[1]=vPathPoints[j].y;
			fPathPointWorld[2]=vPathPoints[j].z;
			float tempdist=zxh::VectorOP_Distance(fPointWorld,fPathPointWorld,3);
			if(tempdist<mindist)
			{
				mindist=tempdist;
				minindex=j;
				int n=1;

			}
		}
		//cout<<"F " <<mindist<<endl;
		sumdistmm=sumdistmm+mindist;
	}
	return sumdistmm;
}
bool SoanchMaxD(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM)
{
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;
	int MCLineMinNUM4;
	for (int j=0;j<nModinFileMaxNUM;j++)
	{
		vector<PointCordTypeDef> vPoints;
		if (vPathPoints.size()!=0)vPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPoints);

			float Min4DistSum=CalcMinDistFromPoinToLine(vPoints,vREFfourPoints,4);
			
            // cout<<vREFfourPoints[0].x<<" "<<vREFfourPoints[0].y<<" "<<vREFfourPoints[0].z<<endl;
			float Min3DistSum=CalcMinDistFromPoinToLine(vPoints,vREFfourPoints,3);
			float Min2DistSum=CalcMinDistFromPoinToLine(vPoints,vREFfourPoints,2);
			float Min1DistSum=CalcMinDistFromPoinToLine(vPoints,vREFfourPoints,1);
			if(Min4DistSum<fMin4Distmm)
			{
				fMin4Distmm=Min4DistSum;
				MCLineMinNUM4=j;
				//cout<<j<<fMin4Distmm<<endl;
			}
		}

		DetectFile.close();

		delete[] chFileName;	

	}


	// read max-path
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	stringstream sstr1;
	string strmacMCL;
	string strj;
	sstr1.clear();
	sstr1<<MCLineMinNUM4;
	sstr1>>strj;
	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr1.clear();
	sstr1<<strmacMCL;
	sstr1>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
	}

	DetectFile.close();
	fMin1Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
	fMin2Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
	fMin3Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
	fMin4Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
	cout<<"MCLine"<<MCLineMinNUM4<<" is the best branch;"<<endl;
	cout<<"Its Sum Distanche is (From proximal to distal): "<<fMin1Distmm<<" mm"<<","<<fMin2Distmm<<" mm"<<","<<fMin3Distmm<<" mm"<<","<<fMin4Distmm<<"mm."<<endl;
	//output to
	string strBestNUMFilename = strMCLineFolder + "/" + string("SMD.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << MCLineMinNUM4 <<" "<<setprecision(4)<< fMin1Distmm<<" "<<setprecision(4)<<fMin2Distmm<<" "<<setprecision(4)<<fMin3Distmm<<" "<<setprecision(4)<<fMin4Distmm<<endl;	
	DetectFile.close();
	return true;
}
int nSoanchMaxD(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM)
{
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;
	int MCLineMinNUM4;
	for (int j=0;j<nModinFileMaxNUM;j++)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);

			float Min4DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
			float Min3DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
			float Min2DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
			float Min1DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
			if(Min4DistSum<fMin4Distmm)
			{
				fMin4Distmm=Min4DistSum;
				MCLineMinNUM4=j;
			}
			if(Min3DistSum<fMin3Distmm)
			{
				fMin3Distmm=Min3DistSum;
			}
			if(Min2DistSum<fMin2Distmm)
			{
				fMin2Distmm=Min2DistSum;
			}
			if(Min1DistSum<fMin1Distmm)
			{
				fMin1Distmm=Min1DistSum;
			}
		}

		DetectFile.close();

		delete[] chFileName;	
	}

	return MCLineMinNUM4;
	// read max-path
	//if (vPathPoints.size()!=0)vPathPoints.clear();
	//	fstream DetectFile;
	//	/*char chIdx[25];
	//	_itoa_s(i, chIdx, 10);*/
	//	char *chFileName=new char[1024];
	//	stringstream sstr1;
	//	string strmacMCL;
	//	string strj;
	//	sstr1.clear();
	//	sstr1<<MCLineMinNUM4;
	//	sstr1>>strj;
	//	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	//	sstr1.clear();
	//	sstr1<<strmacMCL;
	//	sstr1>>chFileName;
	//	DetectFile.open(chFileName,ios::in);
	//	if(DetectFile)
	//	{
	//		ReadVtk(chFileName, vPathPoints);
	//	}

	//	DetectFile.close();
	// cout<<"MCLine"<<MCLineMinNUM4<<" is the best branch;"<<endl;
	// cout<<"Its Sum Distanche is (From proximal to distal): "<<fMin1Distmm<<" mm"<<","<<fMin2Distmm<<" mm"<<","<<fMin3Distmm<<" mm"<<","<<fMin4Distmm<<"mm."<<endl;
	//output to
	/*string strBestNUMFilename = strMCLineFolder + "/" + string("SMD.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << MCLineMinNUM4 <<" "<<setprecision(4)<< fMin1Distmm<<" "<<setprecision(4)<<fMin2Distmm<<" "<<setprecision(4)<<fMin3Distmm<<" "<<setprecision(4)<<fMin4Distmm<<endl;	
	DetectFile.close();
	return true;*/
}
int nSoanchMaxD5R(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM,float fRad[3],float fIo[3])
{
	float fMin4Distmm5=100000;
	float fMin3Distmm5=100000;
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;

	int MCLineMinNUM=-1;
	int MCLineMinNUM2=-1;
	int MCLineMinNUM3=-1;
	int MCLineMinNUM4=-1;
	vector<LineInfo> vMD5NUM4;
	vector<LineInfo> vMD5NUM3;
	float fSUM4R=fRad[0]+fRad[1]+fRad[2]+fRad[3];
	float fSUM3R=fRad[0]+fRad[1]+fRad[2];
	float fSUM2R=fRad[0]+fRad[1];
	float fSUM1R=fRad[0];
	int REFNUM=-1;
	for (int j=nModinFileMaxNUM;j>0;j--)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
			float Min4DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
			float Min3DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
			float Min2DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
			float Min1DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
			if(Min4DistSum<fSUM4R)//store the number of line with the distance less than 5mm
			{
				LineInfo tmpLine;
				tmpLine.num=j;
				tmpLine.dist=Min4DistSum;
				vMD5NUM4.push_back(tmpLine);
				fMin4Distmm5=Min4DistSum;
				MCLineMinNUM4=j;
			}
			if(Min3DistSum<fSUM3R)
			{
				LineInfo tmpLine;
				tmpLine.num=j;
				tmpLine.dist=Min3DistSum;
				vMD5NUM3.push_back(tmpLine);
				fMin3Distmm5=Min3DistSum;
				MCLineMinNUM3=j;
			}

			if(Min1DistSum<fMin1Distmm)
			{
				fMin1Distmm=Min1DistSum;
			}
			if(Min2DistSum<fMin2Distmm)
			{
				fMin2Distmm=Min2DistSum;
				MCLineMinNUM2=j;
			}
			if(Min3DistSum<fMin3Distmm)
			{
				fMin3Distmm=Min3DistSum;
			}	
			if(Min4DistSum<fMin4Distmm)
			{
				fMin4Distmm=Min4DistSum;

			}

		}

		DetectFile.close();

		delete[] chFileName;	
	}
	//select the branch
	if(MCLineMinNUM4==-1)
	{
		if(MCLineMinNUM3==-1)
		{
			return -1;
		}
		else
		{
			float fMin3Dist=100000;
			for(int k=0;k<vMD5NUM3.size();k++)
			{
				if(vMD5NUM3[k].dist<fMin3Dist)
				{
					fMin3Dist=vMD5NUM3[k].dist;
					MCLineMinNUM=vMD5NUM3[k].num;

				}
			}
			return MCLineMinNUM;
		}
	}
	else
	{
		float fMin4Dist=100000;
		float fMin3Dist=100000;

		for(int k=0;k<vMD5NUM4.size();k++)
		{
			if(vMD5NUM4[k].dist<fMin4Dist)
			{
				fMin4Dist=vMD5NUM4[k].dist;
				MCLineMinNUM=vMD5NUM4[k].num;
			}

		}
		return MCLineMinNUM;
	}


}
bool SoanchSDSMD(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM,float fRad[3],float fIo[3])
{
	int nMinNum=-1;
	int nSDMinNum=nSoanchMaxD5R(strMCLineFolder,vREFfourPoints,vPathPoints,nModinFileMaxNUM,fRad,fIo);
	if (nSDMinNum>0)
	{
		nMinNum=nSDMinNum;
	}
	else
	{
		nMinNum	=nSoanchMaxD(strMCLineFolder,vREFfourPoints,vPathPoints,nModinFileMaxNUM);
	}
	if (nMinNum<0)
	{
		cout<<"error: can not find the best branch"<<endl;;
	}

	// read max-path
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	stringstream sstr1;
	string strmacMCL;
	string strj;
	sstr1.clear();
	sstr1<<nMinNum;
	sstr1>>strj;
	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr1.clear();
	sstr1<<strmacMCL;
	sstr1>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
	}
	float fMin1Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
	float fMin2Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
	float fMin3Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
	float fMin4Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
	DetectFile.close();
	string strBestNUMFilename = strMCLineFolder + "/" + string("SDSMD5.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	cout<<"MCLine"<<nMinNum<<" is the best branch;"<<endl;
	cout<<"Its Sum Distanche is (From proximal to distal): "<<fMin1Distmm<<" mm"<<","<<fMin2Distmm<<" mm"<<","<<fMin3Distmm<<" mm"<<","<<fMin4Distmm<<"mm."<<endl;
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << " " << nMinNum<<" "<<setprecision(4)<<fMin1Distmm<<" "<<setprecision(4)<<fMin2Distmm<<" "<<setprecision(4)<<fMin3Distmm<<" "<<setprecision(4)<<fMin4Distmm<<endl;
	return true;
}


bool CalcDistFrPtoD(string strMCLineFolder,char *chResultNameHDr,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints)
{
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;


	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	string strtempstr;
	stringstream sstr;
	strtempstr=strMCLineFolder;
	sstr<<strtempstr;
	sstr>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
		fMin4Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
		fMin3Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
		fMin2Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
		fMin1Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
	}

	DetectFile.close();

	delete[] chFileName;	


	cout<<"Its Sum Distanche is (From proximal to distal): "<<fMin1Distmm<<" mm"<<","<<fMin2Distmm<<" mm"<<","<<fMin3Distmm<<" mm"<<","<<fMin4Distmm<<"mm."<<endl;
	//output to
	string strResultNameHDr;
	strResultNameHDr=chResultNameHDr;
	string strBestNUMFilename = strResultNameHDr + "/" + string("SMDI5_PDDist.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) <<setprecision(4)<< fMin1Distmm<<" "<<setprecision(4)<<fMin2Distmm<<" "<<setprecision(4)<<fMin3Distmm<<" "<<setprecision(4)<<fMin4Distmm<<endl;	
	DetectFile.close();
	return true;
}
bool SoanchMaxD5(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM)
{
	float fMin4Distmm5=100000;
	float fMin3Distmm5=100000;
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;

	int MCLineMinNUM=-1;
	int MCLineMinNUM2=-1;
	int MCLineMinNUM3=-1;
	int MCLineMinNUM4=-1;
	vector<LineInfo> vMD5NUM4;
	vector<LineInfo> vMD5NUM3;
	int REFNUM=-1;
	for (int j=0;j<nModinFileMaxNUM;j++)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
			float Min4DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
			float Min3DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
			float Min2DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
			float Min1DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
			if(Min4DistSum<5)//store the number of line with the distance less than 5mm
			{
				LineInfo tmpLine;
				tmpLine.num=j;
				tmpLine.dist=Min4DistSum;
				vMD5NUM4.push_back(tmpLine);
				fMin4Distmm5=Min4DistSum;
				MCLineMinNUM4=j;
			}
			if(Min3DistSum<5)
			{
				LineInfo tmpLine;
				tmpLine.num=j;
				tmpLine.dist=Min3DistSum;
				vMD5NUM3.push_back(tmpLine);
				fMin3Distmm5=Min3DistSum;
				MCLineMinNUM3=j;
			}

			if(Min1DistSum<fMin1Distmm)
			{
				fMin1Distmm=Min1DistSum;
			}
			if(Min2DistSum<fMin2Distmm)
			{
				fMin2Distmm=Min2DistSum;
				MCLineMinNUM2=j;
			}
			if(Min3DistSum<fMin3Distmm)
			{
				fMin3Distmm=Min3DistSum;
			}	
			if(Min4DistSum<fMin4Distmm)
			{
				fMin4Distmm=Min4DistSum;

			}

		}

		DetectFile.close();

		delete[] chFileName;	
	}
	//select the branch
	if(MCLineMinNUM4=-1)
	{
		if(MCLineMinNUM3=-1)
		{
			MCLineMinNUM=MCLineMinNUM2;
			REFNUM=2;
		}
		else
		{
			float fMin3Dist=100000;
			for(int k=0;k<vMD5NUM3.size();k++)
			{
				if(vMD5NUM3[k].dist<fMin3Dist)
				{
					fMin3Dist=vMD5NUM3[k].dist;
					MCLineMinNUM=vMD5NUM3[k].num;
				}
			}
			REFNUM=3;
		}
	}
	else
	{
		float fMin4Dist=100000;
		float fMin3Dist=100000;

		for(int k=0;k<vMD5NUM4.size();k++)
		{
			if(vMD5NUM4[k].dist<fMin4Dist)
			{
				fMin4Dist=vMD5NUM4[k].dist;
				MCLineMinNUM=vMD5NUM4[k].num;
			}
		}
		REFNUM=4;
	}
	// read max-path
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	stringstream sstr1;
	string strmacMCL;
	string strj;
	sstr1.clear();
	sstr1<<MCLineMinNUM;
	sstr1>>strj;
	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr1.clear();
	sstr1<<strmacMCL;
	sstr1>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
	}
	float MinDistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,REFNUM);
	DetectFile.close();
	string strBestNUMFilename = strMCLineFolder + "/" + string("SMD5.txt");
	ofstream WriteFileBD(strBestNUMFilename);

	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << REFNUM << " " << MCLineMinNUM<<" "<<setprecision(4) <<MinDistSum<<" "<<setprecision(4)<<fMin4Distmm<<" "<<setprecision(4)<<fMin3Distmm<<" "<<setprecision(4)<<fMin2Distmm<<" "<<setprecision(4)<<fMin1Distmm<<endl;	
	return true;
}
bool SoanchMaxD5R(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM,float fRad[3],float fIo[3])
{
	float fMin4Distmm5=100000;
	float fMin3Distmm5=100000;
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;

	int MCLineMinNUM=-1;
	int MCLineMinNUM2=-1;
	int MCLineMinNUM3=-1;
	int MCLineMinNUM4=-1;
	vector<LineInfo> vMD5NUM4;
	vector<LineInfo> vMD5NUM3;
	float fSUM4R=fRad[0]+fRad[1]+fRad[2]+fRad[3];
	float fSUM3R=fRad[0]+fRad[1]+fRad[2];
	float fSUM2R=fRad[0]+fRad[1];
	float fSUM1R=fRad[0];
	int REFNUM=-1;
	for (int j=nModinFileMaxNUM;j>0;j--)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
			float Min4DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
			float Min3DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
			float Min2DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
			float Min1DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
			if(Min4DistSum<fSUM4R)//store the number of line with the distance less than 5mm
			{
				LineInfo tmpLine;
				tmpLine.num=j;
				tmpLine.dist=Min4DistSum;
				vMD5NUM4.push_back(tmpLine);
				fMin4Distmm5=Min4DistSum;
				MCLineMinNUM4=j;
			}
			if(Min3DistSum<fSUM3R)
			{
				LineInfo tmpLine;
				tmpLine.num=j;
				tmpLine.dist=Min3DistSum;
				vMD5NUM3.push_back(tmpLine);
				fMin3Distmm5=Min3DistSum;
				MCLineMinNUM3=j;
			}

			if(Min1DistSum<fMin1Distmm)
			{
				fMin1Distmm=Min1DistSum;
			}
			if(Min2DistSum<fMin2Distmm)
			{
				fMin2Distmm=Min2DistSum;
				MCLineMinNUM2=j;
			}
			if(Min3DistSum<fMin3Distmm)
			{
				fMin3Distmm=Min3DistSum;
			}	
			if(Min4DistSum<fMin4Distmm)
			{
				fMin4Distmm=Min4DistSum;

			}

		}

		DetectFile.close();

		delete[] chFileName;	
	}
	//select the branch
	if(MCLineMinNUM4==-1)
	{
		if(MCLineMinNUM3==-1)
		{
			MCLineMinNUM=MCLineMinNUM2;
			REFNUM=2;
		}
		else
		{
			float fMin3Dist=100000;
			for(int k=0;k<vMD5NUM3.size();k++)
			{
				if(vMD5NUM3[k].dist<fMin3Dist)
				{
					fMin3Dist=vMD5NUM3[k].dist;
					MCLineMinNUM=vMD5NUM3[k].num;
				}
			}
			REFNUM=3;
		}
	}
	else
	{
		float fMin4Dist=100000;
		float fMin3Dist=100000;

		for(int k=0;k<vMD5NUM4.size();k++)
		{
			if(vMD5NUM4[k].dist<fMin4Dist)
			{
				fMin4Dist=vMD5NUM4[k].dist;
				MCLineMinNUM=vMD5NUM4[k].num;
			}
		}
		REFNUM=4;
	}
	// read max-path
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	stringstream sstr1;
	string strmacMCL;
	string strj;
	sstr1.clear();
	sstr1<<MCLineMinNUM;
	sstr1>>strj;
	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr1.clear();
	sstr1<<strmacMCL;
	sstr1>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
	}
	float MinDistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,REFNUM);
	float Min4DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
			float Min3DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
			float Min2DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
			float Min1DistSum=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
	DetectFile.close();
	string strBestNUMFilename = strMCLineFolder + "/" + string("SMD5.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	cout<<"MCLine"<<MCLineMinNUM<<" is the best branch;"<<endl;
	cout<<"Its Sum Distanche is (From proximal to distal): "<<Min1DistSum<<" mm"<<","<<Min2DistSum<<" mm"<<","<<Min3DistSum<<" mm"<<","<<Min4DistSum<<"mm."<<endl;
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << " " << MCLineMinNUM<<" "<<setprecision(4)<<Min1DistSum<<" "<<setprecision(4)<<Min2DistSum<<" "<<setprecision(4)<<Min3DistSum<<" "<<setprecision(4)<<Min4DistSum<<endl;	
	return true;
}

bool SoanchMaxDI5(string strMCLineFolder,vector<PointCordTypeDef > vRefLinePointsWorld,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM)
{
	float fMinSumDistFromModel=100000000;int NMinDistNUM=-1;
	for (int k=0;k<nModinFileMaxNUM;k++)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strk;
		sstr.clear();
		sstr<<k;
		sstr>>strk;
		strtempstr=strMCLineFolder+"/MCLine"+strk+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
		}
		//calculate the mean dist between MCLine and Reference bessel curve
		float fSumDist=0;
		for(int i=0;i<vRefLinePointsWorld.size();i++)
		{
			float fMinDist=10000;
			float fRefPont[3]={vRefLinePointsWorld[i].x,vRefLinePointsWorld[i].y,vRefLinePointsWorld[i].z};

			for(int j=0;j<vPathPoints.size();j++)
			{
				float fPatPont[3]={vPathPoints[j].x,vPathPoints[j].y,vPathPoints[j].z};
				float fDist=zxh::VectorOP_Distance(fRefPont,fPatPont,3);
				if(fDist<fMinDist)
				{
					fMinDist=fDist;
				}
			}//for j
			fSumDist=fSumDist+fMinDist;

		}//for i
		if(fSumDist<fMinSumDistFromModel)
		{
			fMinSumDistFromModel=fSumDist;
			NMinDistNUM=k;
		}

	}//for k
	float fMinAvDist=fMinSumDistFromModel/(vRefLinePointsWorld.size()-1);
	cout<<"The minimal average distance is:"<<fMinAvDist<<endl;
	// read max-path
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	stringstream sstr1;
	string strmacMCL;
	string strj;
	sstr1.clear();
	sstr1<<NMinDistNUM;
	sstr1>>strj;
	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr1.clear();
	sstr1<<strmacMCL;
	sstr1>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
	}
	DetectFile.close();
	string strBestNUMFilename = strMCLineFolder + "/" + string("SMDI5.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << NMinDistNUM << " " <<setprecision(4) <<fMinAvDist<<endl;	
	return true;
}
bool resapEvaLinToRefNum(const vector<PointCordTypeDef> & vOriPoints, float fRefsampling,vector<PointCordTypeDef> &vNewPoints)
{
	vNewPoints.clear();
	float curLen = 0.0;
	float sampleNr = 0;
	int curIdx = 0;

	while (curIdx < int(vOriPoints.size())-1)
	{
		// resample up to nextIdx
		float dF[3]={vOriPoints[curIdx+1].x,vOriPoints[curIdx+1].y,vOriPoints[curIdx+1].z};
		float di[3]={vOriPoints[curIdx].x,vOriPoints[curIdx].y,vOriPoints[curIdx].z};
		float nextLen = curLen + zxh::VectorOP_Distance(dF,di,3);
		while (sampleNr*fRefsampling < nextLen)
		{
			// linearly interpolate between curIdx at curLen, and curIdx at nextLen
			float dist = sampleNr * fRefsampling;
			float a = (nextLen-dist) / (nextLen - curLen); // weight for curIdx
			float b = (dist - curLen) / (nextLen - curLen);// weight for curIdx+1
			PointCordTypeDef PCT;
			PCT.x=a*vOriPoints[curIdx].x+b*vOriPoints[curIdx+1].x;
			PCT.y=a*vOriPoints[curIdx].y+b*vOriPoints[curIdx+1].y;
			PCT.z=a*vOriPoints[curIdx].z+b*vOriPoints[curIdx+1].z;
			vNewPoints.push_back(PCT);
			++sampleNr;
		}
		curLen = nextLen;
		++curIdx;
	}

	return true; 
}
bool SoanchMaxDIR5(string strMCLineFolder,vector<PointCordTypeDef > vRefLinePointsWorld,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM)//with resampling
{
	float fMinSumDistFromModel=100000000;int NMinDistNUM=-1;
	vector<PointCordTypeDef> vResRefPathPoints,vResResPathPoints;
	float fRefLineLen=CalcPathLen(vRefLinePointsWorld);	
	float fResampleSize=10000;
	float fRefsampling=fRefLineLen/fResampleSize;
	vResRefPathPoints.clear();
	resapEvaLinToRefNum( vRefLinePointsWorld, fRefsampling,vResRefPathPoints);
	for (int k=30;k>=0;k--)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strk;
		sstr.clear();
		sstr<<k;
		sstr>>strk;
		strtempstr=strMCLineFolder+"/MCLine"+strk+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		vector<PointCordTypeDef> vNewPoints;
		DetectFile.open(chFileName,ios::in);
		//read cl line
		if(DetectFile)
		{
			ReadVtk(chFileName, vNewPoints);

			//resample the bessel curve to evaluated centerline
			float fMLineLen=CalcPathLen(vNewPoints);
			if(fMLineLen==0)continue;
			float fclsampling=fMLineLen/(vNewPoints.size()-1);
			//resample cl line
			vResResPathPoints.clear();
			if ( fRefsampling/fclsampling < 0.999 || fRefsampling/fclsampling > 1.001 )
			{
				resapEvaLinToRefNum( vNewPoints, fRefsampling,vResResPathPoints);
				//refline length
			}
			//calculate the mean dist between MCLine and Reference bessel curve
			float fSumDist=0;
			for(int i=0;i<vResRefPathPoints.size();i++)
			{
				float fMinDist=10000;
				float fRefPont[3]={vResRefPathPoints[i].x,vResRefPathPoints[i].y,vResRefPathPoints[i].z};

				for(int j=0;j<vResResPathPoints.size();j++)
				{
					float fPatPont[3]={vResResPathPoints[j].x,vResResPathPoints[j].y,vResResPathPoints[j].z};
					float fDist=zxh::VectorOP_Distance(fRefPont,fPatPont,3);
					if(fDist<fMinDist)
					{
						fMinDist=fDist;
					}
				}//for j
				fSumDist=fSumDist+fMinDist;

			}//for i
			if(fSumDist<fMinSumDistFromModel)
			{
				fMinSumDistFromModel=fSumDist;
				NMinDistNUM=k;
			}

		}//for k
		DetectFile.close();
	}
	float fMinAvDist=fMinSumDistFromModel/(vResRefPathPoints.size()-1);
	cout<<"The minimal average distance is:"<<fMinAvDist<<endl;
	// read max-path
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	stringstream sstr1;
	string strmacMCL;
	string strj;
	sstr1.clear();
	sstr1<<NMinDistNUM;
	sstr1>>strj;
	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr1.clear();
	sstr1<<strmacMCL;
	sstr1>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
	}
	DetectFile.close();

	string strBestNUMFilename = strMCLineFolder + "/" + string("SMDRI5.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << NMinDistNUM << " " <<setprecision(4) <<fMinAvDist<<endl;	
	return true;
}
float MinDist(miiCNode<double, float>mdfREFPoints,vector<PointCordTypeDef> &vPathPoints)
{
	float fPathPointWorld[3]={0,0,0};
	float fRefPointWorld[3]={mdfREFPoints.x,mdfREFPoints.y,mdfREFPoints.z};
	float mindist=10000;
	int minindex=100;
	for (int j=0;j<vPathPoints.size();j++)
	{
		fPathPointWorld[0]=vPathPoints[j].x;
		fPathPointWorld[1]=vPathPoints[j].y;
		fPathPointWorld[2]=vPathPoints[j].z;
		float tempdist=zxh::VectorOP_Distance(fRefPointWorld,fPathPointWorld,3);
		if(tempdist<mindist)
		{
			mindist=tempdist;
			minindex=j;
		}
	}
	return mindist;
}

bool SoanchNestAorB(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM)
{
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;
	int MCLineMinNUM;
	int ABNlabel=-1;
	for (int j=0;j<nModinFileMaxNUM;j++)
	{

		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
			float Min3Dist=MinDist(vREFfourPoints[2],vPathPoints);
			float Min2Dist=MinDist(vREFfourPoints[1],vPathPoints);
			float Min4Dist=MinDist(vREFfourPoints[3],vPathPoints);
			float Min1Dist=MinDist(vREFfourPoints[0],vPathPoints);
			if(Min3Dist<5)
			{
				fMin3Distmm=Min3Dist;
				MCLineMinNUM=j;
				ABNlabel=0;

				fMin1Distmm=Min1Dist;
				fMin2Distmm=Min2Dist;
				fMin4Distmm=Min4Dist;
				break;
			}
			else if (Min2Dist<5)
			{
				fMin2Distmm=Min2Dist;
				MCLineMinNUM=j;
				ABNlabel=1;

				fMin1Distmm=Min1Dist;
				fMin3Distmm=Min3Dist;
				fMin4Distmm=Min4Dist;
				break;
			}
			else
			{
				fMin2Distmm=Min2Dist;
				fMin1Distmm=Min1Dist;
				fMin3Distmm=Min3Dist;
				fMin4Distmm=Min4Dist;
				break;
			}

		}

		DetectFile.close();

		delete[] chFileName;	
	}

	if(ABNlabel>=0)
	{
		// read max-path
		if (vPathPoints.size()!=0)vPathPoints.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		stringstream sstr1;
		string strmacMCL;
		string strj;
		sstr1.clear();
		sstr1<<MCLineMinNUM;
		sstr1>>strj;
		strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr1.clear();
		sstr1<<strmacMCL;
		sstr1>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
		}

		DetectFile.close();
		cout<<"MCLine"<<MCLineMinNUM<<" is the best branch;"<<endl;
		if(ABNlabel==0)
		{
			cout<<"It is the closest branch from Point B;"<<endl;

		}
		else
		{
			cout<<"It is the closest branch from Point A;"<<endl;
		}

		//output to
		string strBestNUMFilename = strMCLineFolder + "/" + string("SNAB.txt");
		ofstream WriteFileBD(strBestNUMFilename);
		WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << MCLineMinNUM<<" "<<setprecision(4)<< fMin1Distmm<<" "<<setprecision(4)<<fMin2Distmm<<" "<<setprecision(4)<<fMin3Distmm<<" "<<setprecision(4)<<fMin4Distmm<<endl;	
		DetectFile.close();
	}

	else
	{
		cout<<"The best vessel is not selected;"<<endl;
	}
	cout<<"The  Distances are (From proximal to distal): "<<fMin1Distmm<<" mm"<<","<<fMin2Distmm<<" mm"<<","<<fMin3Distmm<<" mm"<<","<<fMin4Distmm<<"mm."<<endl;
	return true;
}
bool SoanchMaxDI5Cube(string strMCLineFolder,vector< miiCNode<double, float> > vREFfourPoints,vector<PointCordTypeDef > vRefLinePointsWorld, vector<PointCordTypeDef> &vPathPoints,int nModinFileMaxNUM,float fRad[3],float fIo[3])
{
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;
	int MCLineMinNUM4;
	if (vPathPoints.size()!=0)vPathPoints.clear();
	for (int j=nModinFileMaxNUM;j>0;j--)
	{


		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		char *chFileName=new char[1024];
		string strtempstr;
		stringstream sstr;
		string strj;
		sstr.clear();
		sstr<<j;
		sstr>>strj;
		strtempstr=strMCLineFolder+"/MCLine"+strj+".vtk";
		sstr.clear();
		sstr<<strtempstr;
		sstr>>chFileName;	
		vector<PointCordTypeDef> vNewPoints;
		DetectFile.open(chFileName,ios::in);

		if(DetectFile)
		{

			vector<PointCordTypeDef> vResRefPathPoints,vResResPathPoints;
			float fRefLineLen=CalcPathLen(vRefLinePointsWorld);	
			float fResampleSize=10000;
			float fRefsampling=fRefLineLen/fResampleSize;
			vResRefPathPoints.clear();
			resapEvaLinToRefNum( vRefLinePointsWorld, fRefsampling,vResRefPathPoints);
			//read cl line
			ReadVtk(chFileName, vNewPoints);

			//resample the bessel curve to evaluated centerline
			float fMLineLen=CalcPathLen(vNewPoints);	
			float fclsampling=fMLineLen/(vNewPoints.size()-1);
			//resample cl line
			vResResPathPoints.clear();
			if ( fRefsampling/fclsampling < 0.999 || fRefsampling/fclsampling > 1.001 )
			{
				resapEvaLinToRefNum( vNewPoints, fRefsampling,vResResPathPoints);
				//refline length
			}
			//Calculate mean distance and TPR¡¡FN TPM FP

			float Min4DistSum=CalcDistFromBSLineToLine(vResResPathPoints, vREFfourPoints,vResRefPathPoints,nModinFileMaxNUM,fRad, fIo,4);

			if(Min4DistSum<fMin4Distmm)
			{
				fMin4Distmm=Min4DistSum;
				MCLineMinNUM4=j;
			}
		}

		DetectFile.close();

		delete[] chFileName;	
	}


	// read max-path
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	char *chFileName=new char[1024];
	stringstream sstr1;
	string strmacMCL;
	string strj;
	sstr1.clear();
	sstr1<<MCLineMinNUM4;
	sstr1>>strj;
	strmacMCL=strMCLineFolder+"/MCLine"+strj+".vtk";
	sstr1.clear();
	sstr1<<strmacMCL;
	sstr1>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
	}

	DetectFile.close();
	//fMin1Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,1);
	//fMin2Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,2);
	//fMin3Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,3);
	//fMin4Distmm=CalcMinDistFromPoinToLine(vPathPoints,vREFfourPoints,4);
	cout<<"MCLine"<<MCLineMinNUM4<<" is the best branch;"<<endl;
	cout<<"SMDIL: Sum average Distanche is : "<<fMin4Distmm<<"mm."<<endl;
	//output to
	string strBestNUMFilename = strMCLineFolder + "/" + string("SMDIL.txt");
	ofstream WriteFileBD(strBestNUMFilename);
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << MCLineMinNUM4 <<" "<<setprecision(4)<< fMin4Distmm<<endl;	
	DetectFile.close();
	return true;
}


int main(int argc, char *argv[])
{
	//SSL means Soanch-MaxL
	//SSD means Soanch-MinD
	string strMCLineFolder ="";
	char *chResultNameHDr ="";
	char *chRefPFolder ="";
	string RorN="";
	char *REFRFileName="";
	char *chRefLName ="";
	char *chRadUo="";

	vector< miiCNode<double, float> > vREFfourPoints;
	float fRad[4]={0,0,0,0};
	float fIo[4]={0,0,0,0};
	//cout<<argc<<endl;
	//release***
	if( argc == 7 )//SMDI5,SMDIL
	{
		strMCLineFolder =string(argv[1]);
		chResultNameHDr =argv[2];
		chRefPFolder =argv[3];
		chRefLName =argv[4];
		chRadUo=argv[5];
	    RorN=string(argv[6]);
	   // define variables 
	   int iREFRFileLen = strlen(chRefPFolder) + strlen("S.txt");
		char *chLetter[4]={"S.txt","A.txt","B.txt","E.txt"};
		////read the four points selected from reference (pointS,pointA,pointB,pointE)
		for (int i=0;i<4;i++)
		{	
			REFRFileName = (char*)malloc(iREFRFileLen);
			strcpy(REFRFileName, chRefPFolder);
			strcat(REFRFileName, chLetter[i]);
			ifstream str(REFRFileName);
			ReadPointTxt(REFRFileName, vREFfourPoints);
		}
		//order the A B points
		float fPointS[3]={vREFfourPoints[0].x,vREFfourPoints[0].y,vREFfourPoints[0].z};
		float fPointA[3]={vREFfourPoints[1].x,vREFfourPoints[1].y,vREFfourPoints[1].z};
		float fPointB[3]={vREFfourPoints[2].x,vREFfourPoints[2].y,vREFfourPoints[2].z};
		float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
		float fSBDist=zxh::VectorOP_Distance(fPointS,fPointB,3);
		if(fSADist>fSBDist)
		{
			miiCNode<double, float> PTemp=vREFfourPoints[1];
			vREFfourPoints[1]=vREFfourPoints[2];
			vREFfourPoints[2]=PTemp;
		}	
		//SoanchMaxD();
	}
	else if( argc == 6 )//SMD, SD5RSMD,SMD5R
	{
		strMCLineFolder =string(argv[1]);
		chResultNameHDr =argv[2];
		chRefPFolder =argv[3];
		chRadUo=argv[4];
		RorN=string(argv[5]);
		int iREFRFileLen = strlen(chRefPFolder) + strlen("S.txt");

		// define variables 
		char *chLetter[4]={"S.txt","A.txt","B.txt","E.txt"};
		////read the four points selected from reference (pointS,pointA,pointB,pointE)
		for (int i=0;i<4;i++)
		{	
			REFRFileName = (char*)malloc(iREFRFileLen);
			strcpy(REFRFileName, chRefPFolder);
			strcat(REFRFileName, chLetter[i]);
			ifstream str(REFRFileName);
			ReadPointTxt(REFRFileName, vREFfourPoints);
		}
		//order the A B points
		float fPointS[3]={vREFfourPoints[0].x,vREFfourPoints[0].y,vREFfourPoints[0].z};
		float fPointA[3]={vREFfourPoints[1].x,vREFfourPoints[1].y,vREFfourPoints[1].z};
		float fPointB[3]={vREFfourPoints[2].x,vREFfourPoints[2].y,vREFfourPoints[2].z};
		float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
		float fSBDist=zxh::VectorOP_Distance(fPointS,fPointB,3);
		if(fSADist>fSBDist)
		{
			miiCNode<double, float> PTemp=vREFfourPoints[1];
			vREFfourPoints[1]=vREFfourPoints[2];
			vREFfourPoints[2]=PTemp;
		}	
		//SoanchMaxD();
	}
	else if( argc == 5 )//SMD
	{
		strMCLineFolder =string(argv[1]);
		chResultNameHDr =argv[2];
		chRefPFolder =argv[3];
		RorN=string(argv[4]);

		int iREFRFileLen = strlen(chRefPFolder) + strlen("S.txt");

		// define variables 
		char *chLetter[4]={"S.txt","A.txt","B.txt","E.txt"};
		////read the four points selected from reference (pointS,pointA,pointB,pointE)



		for (int i=0;i<4;i++)
		{	
			REFRFileName = (char*)malloc(iREFRFileLen);
			strcpy(REFRFileName, chRefPFolder);
			strcat(REFRFileName, chLetter[i]);
			ifstream str(REFRFileName);
			ReadPointTxt(REFRFileName, vREFfourPoints);
		}
		//order the A B points
		float fPointS[3]={vREFfourPoints[0].x,vREFfourPoints[0].y,vREFfourPoints[0].z};
		float fPointA[3]={vREFfourPoints[1].x,vREFfourPoints[1].y,vREFfourPoints[1].z};
		float fPointB[3]={vREFfourPoints[2].x,vREFfourPoints[2].y,vREFfourPoints[2].z};
		float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
		float fSBDist=zxh::VectorOP_Distance(fPointS,fPointB,3);
		if(fSADist>fSBDist)
		{
			miiCNode<double, float> PTemp=vREFfourPoints[1];
			vREFfourPoints[1]=vREFfourPoints[2];
			vREFfourPoints[2]=PTemp;
		}	

	}
	else
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeSoanch	MCLineFolder ResultName (reference folder SMD needed) -SML(SMD)" << endl;
		return -1;
	}
	//**release

	//**Debug
	/*strMCLineFolder="F:/Coronary_0/trainningdataZXHCAEDMP/DirHighResultsDFM/mod01_to_unseen00_results/meanimg01v0model";
	chResultNameHDr ="F:/Coronary_0/trainningdataZXHCAEDMP/DirHighResultsDFM/mod01_to_unseen00_results/meanimg01v0model/MaxMDI5CLine";
	chRefPFolder="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel0/point";
	chRefLName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel0/BCFromRefPoint1.vtk";
	RorN="-SMDI5";*/
	 /*strMCLineFolder ="F:/Coronary_0/ZXHJDQCAEDMP/DHResDFM_VVS_VSPSC6_BPD_T/mod00_to_unseen01_results/meanimg00v0model";
	chResultNameHDr ="F:/Coronary_0/ZXHJDQCAEDMP/DHResDFM_VVS_VSPSC6_BPD_T/mod00_to_unseen01_results/meanimg00v0model/SMDLine";
	chRefPFolder ="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset01/vessel0/point";
	chRadUo="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/RadIo/vessel0.txt";
	RorN="-SMDIL";
	chRefLName="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset01/vessel0/reference.vtk";*/
	
	//int iREFRFileLen = strlen(chRefPFolder) + strlen("S.txt");
	// define variables 
	//char *chLetter[4]={"S.txt","A.txt","B.txt","E.txt"};
	//read the four points selected from reference (pointS,pointA,pointB,pointE)
	//for (int i=0;i<4;i++)
	//{	
	//	REFRFileName = (char*)malloc(iREFRFileLen);
	//	strcpy(REFRFileName, chRefPFolder);
	//	strcat(REFRFileName, chLetter[i]);
	//	ifstream str(REFRFileName);
	//	ReadPointTxt(REFRFileName, vREFfourPoints);
	//}
	//order the A B points
		/*float fPointS[3]={vREFfourPoints[0].x,vREFfourPoints[0].y,vREFfourPoints[0].z};
	float fPointA[3]={vREFfourPoints[1].x,vREFfourPoints[1].y,vREFfourPoints[1].z};
	float fPointB[3]={vREFfourPoints[2].x,vREFfourPoints[2].y,vREFfourPoints[2].z};
	float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
	float fSBDist=zxh::VectorOP_Distance(fPointS,fPointB,3);
	if(fSADist>fSBDist)
	{
	miiCNode<double, float> PTemp=vREFfourPoints[1];
	vREFfourPoints[1]=vREFfourPoints[2];
	vREFfourPoints[2]=PTemp;
	}*/
	//read the radius
	ReadRadIoTxt(chRadUo, fRad,fIo);
	int nModinFileMaxNUM=30;
	//**Debug
	vector<PointCordTypeDef>vPathPoints;//define the final result points

	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;
	int MCLineMinNUM;
	if (RorN=="-SML")//
	{
		cout<<"Search for the Maximal length line."<<endl;
		SoanchMaxL(strMCLineFolder,vPathPoints,nModinFileMaxNUM);
	}
	if (RorN=="-SMD")//
	{

		cout<<"Search for the Minimal Distance line using SMD."<<endl;
		SoanchMaxD(strMCLineFolder,vREFfourPoints,vPathPoints,nModinFileMaxNUM);
	}
	//Soanch-MinD5mm
	if (RorN=="-SMD5R")//strick conditions
	{

		cout<<"Search for the Minimal Distance line using SD5R."<<endl;
			//read the radius and io
		ReadRadIoTxt(chRadUo, fRad,fIo);
		SoanchMaxD5R(strMCLineFolder,vREFfourPoints,vPathPoints,nModinFileMaxNUM,fRad,fIo);

	}
	vector<PointCordTypeDef>vRefLinePointsWorld;
	if (RorN=="-SMDIR5")//
	{

		if (vRefLinePointsWorld.size()!=0)vRefLinePointsWorld.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		DetectFile.open(chRefLName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chRefLName, vRefLinePointsWorld);
		}
		DetectFile.close();
		cout<<"Search for the best branch using SMDIR5."<<endl;
		SoanchMaxDIR5(strMCLineFolder,vRefLinePointsWorld,vPathPoints,nModinFileMaxNUM);
	}
	if (RorN=="-CalCPDDist")//
	{
		CalcDistFrPtoD(strMCLineFolder,chResultNameHDr,vREFfourPoints,vPathPoints);
		exit(0);
	}
	if (RorN=="-SNAorB")//
	{
		cout<<"Search for the minimal distance less than 5mm line from A or B."<<endl;
		SoanchNestAorB(strMCLineFolder,vREFfourPoints,vPathPoints,nModinFileMaxNUM);
	}
	if (RorN=="-SD5RSMD")//
	{
		cout<<"Select the best branch using SD5R first and then SMD."<<endl;
		SoanchSDSMD(strMCLineFolder,vREFfourPoints,vPathPoints,nModinFileMaxNUM,fRad,fIo);
	}
	if (RorN=="-SMDIL")//
	{
		if (vRefLinePointsWorld.size()!=0)vRefLinePointsWorld.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		DetectFile.open(chRefLName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chRefLName, vRefLinePointsWorld);
		}
		DetectFile.close();

		cout<<"Search for the best branch using SMDIL."<<endl;
		SoanchMaxDI5Cube(strMCLineFolder,vREFfourPoints,vRefLinePointsWorld,vPathPoints,nModinFileMaxNUM,fRad,fIo);
	}
	int nFileLen = strlen(chResultNameHDr) + strlen(".vtk")+1 ;
	char *chFileNamevtk = (char*)malloc(nFileLen);
	strcpy(chFileNamevtk, chResultNameHDr);
	strcat(chFileNamevtk, ".vtk");
	WriteVtk(vPathPoints, chFileNamevtk);
	free(chFileNamevtk);
	char *chFileNametxt = (char*)malloc(nFileLen);
	strcpy(chFileNametxt, chResultNameHDr);
	strcat(chFileNametxt, ".txt");
	WriteTxt(vPathPoints, chFileNametxt);
	free(chFileNametxt);
	cout << RorN<<" "<<"has already completed!" << endl;


}

