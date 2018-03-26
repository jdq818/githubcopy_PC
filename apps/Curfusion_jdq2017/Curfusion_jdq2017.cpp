

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
#define _CRT_SECURE_NO_WARNINGS
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
#include <functional>
#include <numeric>
#include <array>  

#include <string>
#include <vector>
#include <math.h>
#include<stdlib.h>
#include <stdio.h>

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"


#include <algorithm>
#include "jdq2017util.h"
#include "jdqdijkstra.h"
#include "jdqPoint.h"

#include <time.h>

#include <Eigen/Dense>  
using namespace Eigen;  
using namespace Eigen::internal;  
using namespace Eigen::Architecture; 

#define M_PI 3.14159265358979323846
// change to true for debug output
const bool debugOutput = false;
// Change to true to output all the measures on one line
// this line will start with the dataset and vessel number
const bool outputOneLine = true;


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
	float x;
	float y;
	float z;
	vector<int>Lnum;
}PointCLTypeDef;
typedef struct
{
	int x;
	int y;
	int z;
	short val;
}PointImgTypeDef;
typedef struct
{
	int cunt;
	PointImgTypeDef pont;
	vector<PointImgTypeDef> vpont;
}PointLinkDef;
typedef struct
{
	float fMean;
	float fStdVar;
	float fMaxErr;
}ErrTypeDef;

typedef struct
{
	std::vector<jdq2017::point3D> vLine;
	int nLine;
}vLinesDef;

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
void WriteCA2Txt(vector<PointCordTypeDef> vPoints,char *chFileName)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = vPoints.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = vPoints[i].x;
		fImgPixel[1] = vPoints[i].y;
		fImgPixel[2] = vPoints[i].z;
		WriteFileTxt << right<<fixed<<setprecision(4) <<fImgPixel[0] << " " << fImgPixel[1] << " " << fImgPixel[2] << '\n';
	}	

}
void WriteCA2Txt_jdq(std::vector<jdq2017::point3D> vPoints,char *chFileName)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = vPoints.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = vPoints[i]._x;
		fImgPixel[1] = vPoints[i]._y;
		fImgPixel[2] = vPoints[i]._z;
		WriteFileTxt << right<<fixed<<setprecision(4) <<fImgPixel[0] << " " << fImgPixel[1] << " " << fImgPixel[2] << '\n';
	}	

}
void WriteCA2Txt_Skip(vector<PointCordTypeDef> vPoints,char *chFileName)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = vPoints.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i+=6)
	{
		fImgPixel[0] = vPoints[i].x;
		fImgPixel[1] = vPoints[i].y;
		fImgPixel[2] = vPoints[i].z;
		WriteFileTxt << right<<fixed<<setprecision(4) <<fImgPixel[0] << " " << fImgPixel[1] << " " << fImgPixel[2] << '\n';
	}	

}
void ReadCA2Txt(vector<PointCordTypeDef> &vPoints,char *chFileName)
{
	ifstream ReadFileTxt;
	ReadFileTxt.open(chFileName);
	if(!ReadFileTxt)
	{
		cout << "Reading failed！" << endl;
	}
	float fImgPixel[3];	
	while(!ReadFileTxt.eof())
	{
		ReadFileTxt>>fixed>>setprecision(4)>>fImgPixel[0] >>fImgPixel[1]>>fImgPixel[2];
		PointCordTypeDef Ptmp;
		Ptmp.x=fImgPixel[0];
		Ptmp.y=fImgPixel[1];
		Ptmp.z=fImgPixel[2];
		vPoints.push_back(Ptmp);
	}
	ReadFileTxt.close();
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
void init_img(zxhImageData &imgCountmap)
{
	int ImgSize[4]={1};
	imgCountmap.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	for(int it=0;it<ImgSize[3];++it)
		for(int iz=0;iz<ImgSize[2];++iz)
			for(int iy=0;iy<ImgSize[1];++iy)
				for(int ix=0;ix<ImgSize[0];++ix)
				{
					imgCountmap.SetPixelByGreyscale(ix,iy,iz,it,0);
				}
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
float fCalcTwoPoiDist(const PointCordTypeDef pPoint1,const PointCordTypeDef pPoint2)
{
	float fF[3]={0},fB[3]={0};
	fF[0]=pPoint1.x;
	fF[1]=pPoint1.y;
	fF[2]=pPoint1.z;
	fB[0]=pPoint2.x;
	fB[1]=pPoint2.y;
	fB[2]=pPoint2.z;
	return zxh::VectorOP_Distance(fF, fB, 3 ) ;  
}
float CalcTotaldist(const vector<PointCordTypeDef> vPathPoints)//Calculate the mean distance density of model line.
{
	
	float Dist=0;
	for(int i=1;i<vPathPoints.size();i++)
	{
		Dist=Dist+fCalcTwoPoiDist(vPathPoints[i-1],vPathPoints[i]);
		
	}
	return Dist;
}
void ResampleLine(const vector<PointCordTypeDef> vPoints,float ResampleDen,vector<PointCordTypeDef> &vResPoints)
{
	//Check for valid samplingdistance
	if(ResampleDen==0.0)
	{
		vResPoints=vPoints;
		return;
	}
	//Resampled result path is stored in _resultPath
	vResPoints.clear();
	double curLen = 0.0;
	double sampleNr = 0.0;
	int curIdx = 0;
	//temporary point
	PointCordTypeDef ptempPoint;
	while (curIdx < int(vPoints.size())-1)
	{
		// resample up to nextIdx
		double nextLen = curLen + fCalcTwoPoiDist(vPoints[curIdx+1],vPoints[curIdx]);

		while (sampleNr*ResampleDen < nextLen)
		{
			// linearly interpolate between curIdx at curLen, and curIdx at nextLen
			double dist = sampleNr * ResampleDen;
			double a = (nextLen-dist) / (nextLen - curLen); // weight for curIdx
			double b = (dist - curLen) / (nextLen - curLen);// weight for curIdx+1
			ptempPoint.x=a*vPoints[curIdx].x + b*vPoints[curIdx+1].x;
			ptempPoint.y=a*vPoints[curIdx].y + b*vPoints[curIdx+1].y;
			ptempPoint.z=a*vPoints[curIdx].z + b*vPoints[curIdx+1].z;
			vResPoints.push_back(ptempPoint);
			++sampleNr;
		}
		curLen = nextLen;
		++curIdx;
	}

	// add last sample
	vResPoints.push_back(vPoints.back());
}
template <class PixelType>
void SetResMolPoiInt(vector<PointCordTypeDef> vResPoints,zxhImageDataT<PixelType>&zxhImg,vector<short> &sMolPontInts)
{
	float fCord[3];
	int iCord[3];
	int iImgWX,iImgWY,iImgWZ,iImgWT;
	sMolPontInts.clear();
	zxhImg.GetImageSize(iImgWX, iImgWY, iImgWZ, iImgWT);
	for(int i=0;i<vResPoints.size();i++)
	{
		fCord[0]=vResPoints[i].x;
		fCord[1]=vResPoints[i].y;
		fCord[2]=vResPoints[i].z;
		zxhImg.GetImageInfo()->WorldToImage(fCord);
		iCord[0]=int(fCord[0]+0.5);
		iCord[1]=int(fCord[1]+0.5);
		iCord[2]=int(fCord[2]+0.5);
		if(iCord[0]<0)iCord[0]=0;
		if(iCord[0]>=iImgWX)iCord[0]=iImgWX-1;
		if(iCord[1]<0)iCord[1]=0;
		if(iCord[1]>=iImgWY)iCord[1]=iImgWY-1;
		if(iCord[2]<0)iCord[2]=0;
		if(iCord[2]>=iImgWZ)iCord[2]=iImgWZ-1;
		sMolPontInts.push_back(zxhImg.GetPixelGreyscale(iCord[0],iCord[1],iCord[2],0 ));
	}

}
bool GenerateCTAPointsForANewCenterLineWithNewSize(
 const vector<PointCordTypeDef> & vOriPoints, // 这里不用指针或引用，会导致原来的整个vector被copy一份到这个函数中来
 int iNewSize,vector<PointCordTypeDef> &vNewPoints)
{
    int iSizeOri=vOriPoints.size();
    if( iNewSize<2 || iSizeOri<2 )
        return false ; 
    vNewPoints.empty(); //保证这个vector是空的
    for( int i=0; i<iNewSize; ++i )
    {
        float fProportion = float(i)/float(iNewSize-1); //0 to 1
        float fIndexOfOri = fProportion * float(iSizeOri-1) ; 
        int iIndexOfOriFrom = int(fIndexOfOri) ;
        int iIndexOfOriTo=iIndexOfOriFrom+1 ; 
        if( iIndexOfOriTo > iSizeOri-1 )
            iIndexOfOriTo = iSizeOri-1 ; 
        float fPropForTo = fIndexOfOri-iIndexOfOriFrom;
        float fPropForFrom = 1-fPropForTo;  
        PointCordTypeDef ptempPoint;        
        ptempPoint.x=fPropForFrom*vOriPoints[iIndexOfOriFrom].x + fPropForTo*vOriPoints[fPropForTo].x;
        ptempPoint.y=fPropForFrom*vOriPoints[iIndexOfOriFrom].y + fPropForTo*vOriPoints[fPropForTo].y;
        ptempPoint.z=fPropForFrom*vOriPoints[iIndexOfOriFrom].z + fPropForTo*vOriPoints[fPropForTo].z;
        vNewPoints.push_back(ptempPoint);
    }
    return true ; //根据返回值来判断你的这个函数是否正确被执行了
}

bool GetLenOfCurve(vector<PointCordTypeDef> vRefCurPontsWorld,float &fWhLen)
{
	fWhLen=0;
	for(int i=0;i<vRefCurPontsWorld.size()-1;i++)
	{
		float fCurV1[3]={vRefCurPontsWorld[i].x,vRefCurPontsWorld[i].y,vRefCurPontsWorld[i].y};
		float fCurV2[3]={vRefCurPontsWorld[i+1].x,vRefCurPontsWorld[i+1].y,vRefCurPontsWorld[i+1].y};
		float fTemp=zxh::VectorOP_Distance(fCurV1,fCurV2,3);
		fWhLen=fWhLen+fTemp;
	}
	return true;
}
bool ResampleCurve(vector<PointCordTypeDef> path,double samplingDistance,vector<PointCordTypeDef> &vc1)
{
	if(samplingDistance==0)
		vc1=path;
	vc1.clear();
	double curLen = 0.0;
	double sampleNr = 0;
	int curIdx = 0;
	while (curIdx < int(path.size())-1)
	{
		// resample up to nextIdx
		float dNextPont[3]={path[curIdx+1].x,path[curIdx+1].y,path[curIdx+1].z};
		float dCurPont[3]={path[curIdx].x,path[curIdx].y,path[curIdx].z};
		float dLen=zxh::VectorOP_Distance(dNextPont,dCurPont,3);
		double nextLen = curLen + dLen;
		while (sampleNr*samplingDistance < nextLen)
		{
			// linearly interpolate between curIdx at curLen, and curIdx at nextLen
			double dist = sampleNr * samplingDistance;
			double a = (nextLen-dist) / (nextLen - curLen); // weight for curIdx
			double b = (dist - curLen) / (nextLen - curLen);// weight for curIdx+1
			PointCordTypeDef TempP;
			TempP.x=a*path[curIdx].x+b*path[curIdx+1].x;
			TempP.y=a*path[curIdx].y+b*path[curIdx+1].y;
			TempP.z=a*path[curIdx].z+b*path[curIdx+1].z;
			vc1.push_back(TempP);;
			++sampleNr;
		}
		curLen = nextLen;
		++curIdx;
	}

	return true;
}
void CorrectImagePos(PointCordTypeDef &diPoint,zxhImageDataT<short> &imgReadRaw)
{
	int nImWX, nImWY, nImWZ, nImWT;
	imgReadRaw.GetImageSize(nImWX, nImWY, nImWZ, nImWT);
	if (diPoint.x >= nImWX)
		diPoint.x = nImWX - 1;
	if (diPoint.y >= nImWY)
		diPoint.y = nImWY - 1;
	if (diPoint.z >= nImWZ)
		diPoint.z = nImWZ - 1;
	if (diPoint.x <= 0)
		diPoint.x = 0;
	if (diPoint.y <= 0)
		diPoint.y = 0;
	if (diPoint.z <= 0)
		diPoint.z = 0;
};
bool BoundaryCorrect(int *PointPos,int ImgNewVslsSize[4])
{
	for (int i=0;i<3;i++)
	{
		PointPos[i]=zxh::maxf(0,PointPos[i]);
		PointPos[i]=zxh::minf(ImgNewVslsSize[i]-1,PointPos[i]);
	}
	return true;
}
bool TransAllPontsIntoUnorga(vector<vLinesDef>&vline,vector<PointCordTypeDef> &vUnorgaPointsWorld)
{
	int nvline=vline.size();
	if(!vUnorgaPointsWorld.empty())
		vUnorgaPointsWorld.clear();

	for(int i=0;i<nvline;i++)
	{
		int npont=vline[i].vLine.size();
		for(int j=0;j<npont;j++)
		{
			PointCordTypeDef PonDeftmp;
			PonDeftmp.x=(vline[i].vLine)[j]._x;
			PonDeftmp.y=(vline[i].vLine)[j]._y;
			PonDeftmp.z=(vline[i].vLine)[j]._z;
			vUnorgaPointsWorld.push_back(PonDeftmp);
		}
	}
	return true;
}
bool Collet1(PointCordTypeDef PCent,float H,vector<PointCordTypeDef> &vUnorgaPointsWorld,vector<PointCordTypeDef> &vLocalPointsInBallWorld)
{
	float fCent[3]={PCent.x,PCent.y,PCent.z};
	int nsize=vUnorgaPointsWorld.size();
	if(!vLocalPointsInBallWorld.empty())
		vLocalPointsInBallWorld.clear();
	for(int i=0;i<nsize;i++)
	{
		float fTemp[3]={vUnorgaPointsWorld[i].x,vUnorgaPointsWorld[i].y,vUnorgaPointsWorld[i].z};
		float fdist=zxh::VectorOP_Distance(fCent,fTemp,3);
		if(fdist<=H)
		{
			vLocalPointsInBallWorld.push_back(vUnorgaPointsWorld[i]);
		}
	}

	return true;
}

bool StorUnique(vector<PointCordTypeDef> &vUnorgaPointsWorld)
{
	bool bunique=true;
	for(int i=0;i<vUnorgaPointsWorld.size();i++)
	{
		float fi[3]={vUnorgaPointsWorld[i].x,vUnorgaPointsWorld[i].y,vUnorgaPointsWorld[i].z};
		for(int j=i+1;j<vUnorgaPointsWorld.size();j++)
		{
			float fj[3]={vUnorgaPointsWorld[j].x,vUnorgaPointsWorld[j].y,vUnorgaPointsWorld[j].z};
			float dist=zxh::VectorOP_Distance(fi,fj,3);
			if(dist<0.00001)
			{
				vUnorgaPointsWorld.erase(vUnorgaPointsWorld.begin()+j);
				j--;
			}
		}
	}
	return true;
}
bool MapCurvPontsToImage(vector<PointCordTypeDef> &vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls,int SearchRange[4])//map line into original image in a range
{
	int ImgNewVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	float ImgNewVslsSpacing[4]={0};
	float fminpixdist=1*sqrt(ImgNewVslsSpacing[0]*ImgNewVslsSpacing[0]+ImgNewVslsSpacing[1]*ImgNewVslsSpacing[1]+ImgNewVslsSpacing[2]*ImgNewVslsSpacing[2]);
	vector<PointCordTypeDef> vPathPointsWorldMAPT;
	vPathPointsWorldMAPT.clear();
	for (int i=0;i<vPathPointsWorld.size();i++)//map the vetor points to the image
	{
		float PointPosWorld[ZXH_ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointCordTypeDef PointMAPT;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		imgReadNewVsls.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=zxh::round(PointPosWorld[1]);
		PointPos[2]=zxh::round(PointPosWorld[2]);
		BoundaryCorrect(PointPos,ImgNewVslsSize);
		PointMAPT.x=PointPos[0];
		PointMAPT.y=PointPos[1];
		PointMAPT.z=PointPos[2];
		vPathPointsWorldMAPT.push_back(PointMAPT);
		imgReadNewVsls.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);

	}

	for (int mapNUM=1;mapNUM<vPathPointsWorldMAPT.size()-1;mapNUM++)//插值
	{
		PointCordTypeDef PointMAPTS;//start point in every step
		PointCordTypeDef PointMAPTE;//End point in every step
		PointMAPTS=vPathPointsWorldMAPT[mapNUM-1];
		PointMAPTE=vPathPointsWorldMAPT[mapNUM];
		for(int jz=PointMAPTS.z-8;jz<PointMAPTS.z+8;++jz)
			for(int jy=PointMAPTS.y-8;jy<PointMAPTS.y+8;++jy)
				for(int jx=PointMAPTS.x-8;jx<PointMAPTS.x+8;++jx)
				{
					float fL[4]={0};
					float fM[4]={0};
					float fN[4]={0};
					fL[0]=PointMAPTE.x-PointMAPTS.x;
					fL[1]=PointMAPTE.y-PointMAPTS.y;
					fL[2]=PointMAPTE.z-PointMAPTS.z;
					fM[0]=jx-PointMAPTS.x;
					fM[1]=jy-PointMAPTS.y;
					fM[2]=jz-PointMAPTS.z;
					fN[0]=jx-PointMAPTE.x;
					fN[1]=jy-PointMAPTE.y;
					fN[2]=jz-PointMAPTE.z;
					float fcosineLM=zxh::VectorOP_Cosine(fL,fM,3);
					float fcosineLN=-zxh::VectorOP_Cosine(fL,fN,3);
					float fLengthM=zxh::MagnitudeOfVector(fM,3);
					float fDist=fLengthM*sqrt(1-fcosineLM*fcosineLM);
					if(fcosineLM>0.9&&fcosineLN>0.9&&fDist<=fminpixdist)
					{
						imgReadNewVsls.SetPixelByGreyscale(jx,jy,jz,0,ZXH_Foreground);		
					}
				}
	}
	return true;
}
bool MapCurvPontsToTNImage(vector<PointCordTypeDef> &vPathPointsWorld,zxhImageDataT<short>&imgReadNewRaw,int SearchRange[4])//map the curve to a totally new image
{
	int ImgNewSize[4]={1};
	imgReadNewRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					imgReadNewRaw.SetPixelByGreyscale(ix,iy,iz,it,0);
				}
				MapCurvPontsToImage(vPathPointsWorld,imgReadNewRaw,SearchRange);
				return true;
}

bool MapCurvPontsToTarImage(vector<PointCordTypeDef> &vPathPointsWorld,zxhImageDataT<short>&imgReadNewRaw,int SearchRange[4])//map the curve to a target image
{
				MapCurvPontsToImage(vPathPointsWorld,imgReadNewRaw,SearchRange);
				return true;
}

template <typename Point>
bool jdqpontToPonDef(const std::vector<Point> &c1,vector<PointCordTypeDef> &vPathPointsWorld)
{
	if(!vPathPointsWorld.empty())
		vPathPointsWorld.clear();
	for (typename std::vector<Point>::const_iterator i = c1.begin()+1;
		i != c1.end();  ++i)
	{
		PointCordTypeDef PontWorld;
		PontWorld.x=(*i)._x;
		PontWorld.y=(*i)._y;
		PontWorld.z=(*i)._z;
		vPathPointsWorld.push_back(PontWorld);	
	}
	return true;
}

bool PlaneFitting(PointCordTypeDef &PfirPont,float H,vector<PointCordTypeDef> &vPathPointsWorld,float fabc[4])
{
	//z=a*x+b*y+c;
	int n=vPathPointsWorld.size();
	Array<double, 1, 1000> Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mz(0,i)=0;	
		Mw(0,i)=0;	
	}
	//tranform the point position to matrix
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=vPathPointsWorld[i].x;
		My(0,i)=vPathPointsWorld[i].y;	
		Mz(0,i)=vPathPointsWorld[i].z;	
	}
		//tranform the point position to matrix
	float fmfirsyz[3]={PfirPont.x,PfirPont.y,PfirPont.z};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=vPathPointsWorld[i].x;
		My(0,i)=vPathPointsWorld[i].y;	
		Mz(0,i)=vPathPointsWorld[i].z;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxyz[3]={vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z};
		float fr2=(zxh::VectorOP_Distance(fmfirsyz,fmxyz,3))*(zxh::VectorOP_Distance(fmfirsyz,fmxyz,3));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=1;
	}

	Matrix<double, 3, 3>MA,MA_inv;
	Matrix<double, 3, 1>Mb;
	Matrix<double, 3, 1>Mabc;
	//Ax=b
	//generate A matrx
	MA(0,0)=(Mw*Mx*Mx).sum();
	MA(1,1)=(Mw*My*My).sum();
	MA(2,2)=Mw.sum();

	MA(0,1)=(Mw*Mx*My).sum();
	MA(1,0)=(Mw*Mx*My).sum();

	MA(0,2)=Mx.sum();
	MA(2,0)=Mx.sum();

	MA(1,2)=My.sum();
	MA(2,1)=My.sum();
	//generate b matrx
	Mb(0,0)=(Mw*Mx*Mz).sum();
	Mb(1,0)=(Mw*My*Mz).sum();
	Mb(2,0)=(Mw*Mz).sum();
	float fmadet=MA.determinant();


	////gerate A_inverse
	//MA_inv=MA.inverse();

	////result of abc
	//Mabc=MA_inv*Mb;
	Mabc = MA.ldlt().solve(Mb);
	fabc[0]=Mabc(0,0);
	fabc[1]=Mabc(1,0);
	fabc[2]=-1;
	fabc[3]=Mabc(2,0);
	//

	Matrix<double, 3, 1>xx;
	xx=MA*Mabc-Mb;
	////Mean square deviation
	float fmsd=0;
	for(int i=0;i<n;i++)
	{
		float fcurpont[3]={vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z};	
		float f_dist_Curvefit=fabsf(fabc[0]*fcurpont[0]+fabc[1]*fcurpont[1]+fabc[2]*fcurpont[2]+fabc[3])/sqrt(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
		fmsd=fmsd+f_dist_Curvefit;
	}
	return true;
}
bool ProjPontsToPlaen(vector<PointCordTypeDef> & vLocalPointsInBallWorld,float fabc[4],PointCordTypeDef Pcent, Matrix<float,1,3> &MPro_Pcent,vector<PointCordTypeDef> &Pro_vLocalPointsInBallWorld)
{
	if(!Pro_vLocalPointsInBallWorld.empty())
		Pro_vLocalPointsInBallWorld.clear();
	int n=vLocalPointsInBallWorld.size();
	//find the center point
	int ncentid=0;
	for(int i=0;i<n;i++)
	{
		if(Pcent.x==vLocalPointsInBallWorld[i].x&&Pcent.y==vLocalPointsInBallWorld[i].y&&Pcent.z==vLocalPointsInBallWorld[i].z)
			ncentid=i;
	}
	//other points
	for(int i=0;i<n;i++)
	{
		float fxyz[3]={vLocalPointsInBallWorld[i].x,vLocalPointsInBallWorld[i].y,vLocalPointsInBallWorld[i].z};
		float ft=(-fabc[0]*fxyz[0]-fabc[1]*fxyz[1]-fabc[2]*fxyz[2]-fabc[3])/(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
		PointCordTypeDef TempPont;
		TempPont.x=ft*fabc[0]+fxyz[0];
		TempPont.y=ft*fabc[1]+fxyz[1];
		TempPont.z=ft*fabc[2]+fxyz[2];
		Pro_vLocalPointsInBallWorld.push_back(TempPont);
	}
	MPro_Pcent(0,0)=Pro_vLocalPointsInBallWorld[ncentid].x;
	MPro_Pcent(0,1)=Pro_vLocalPointsInBallWorld[ncentid].y;
	MPro_Pcent(0,2)=Pro_vLocalPointsInBallWorld[ncentid].z;
	//Mean square deviation
	float fmsd=0;
	float fmsd1=0;
	for(int i=0;i<n;i++)
	{
		float fcurpont[3]={Pro_vLocalPointsInBallWorld[i].x,Pro_vLocalPointsInBallWorld[i].y,Pro_vLocalPointsInBallWorld[i].z};	
		float f_dist_Curvefit=fabsf(fabc[0]*fcurpont[0]+fabc[1]*fcurpont[1]+fabc[2]*fcurpont[2]+fabc[3])/sqrt(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
		float fcurpont_befpro[3]={vLocalPointsInBallWorld[i].x,vLocalPointsInBallWorld[i].y,vLocalPointsInBallWorld[i].z};	
		float f_dist_Curvefit1=zxh::VectorOP_Distance(fcurpont,fcurpont_befpro,3);
		fmsd1=fmsd1+f_dist_Curvefit1;
		fmsd=fmsd+f_dist_Curvefit;
	}
	return true;
}
int BelotoQuadr(float fabc[2],char* chplane)
{
	if(strcmp(chplane,"yz")==0)
	{
		if(fabc[0]>=0&&fabc[1]>=0)
		{
			return 1;
		}
		if(fabc[0]<0&&fabc[1]>0)
		{
			return 2;
		}
		if(fabc[0]<0&&fabc[1]<0)
		{
			return 3;
		}
		if(fabc[0]>0&&fabc[1]<0)
		{
			return 4;
		}
	}
	if(strcmp(chplane,"xz")==0)
	{
		if(fabc[0]>0&&fabc[1]>0)
		{
			return 2;
		}
		if(fabc[0]<=0&&fabc[1]>=0)
		{
			return 1;
		}
		if(fabc[0]<0&&fabc[1]<0)
		{
			return 4;
		}
		if(fabc[0]>0&&fabc[1]<0)
		{
			return 3;
		}
	}
	return -1;
}
bool CalcRotMatrp(int nxx,float fprj_v_xzp[3],float fglozaxis[3],float &fcostheta,float &fsintheta)
{
	switch(nxx)
	{
	case 1:
		{
			fcostheta=zxh::VectorOP_Cosine(fglozaxis,fprj_v_xzp,3);
			fsintheta=sqrt(1-fcostheta*fcostheta);
			break;
			return true;
		}
	case 2:
		{
			fcostheta=zxh::VectorOP_Cosine(fglozaxis,fprj_v_xzp,3);
			fsintheta=-sqrt(1-fcostheta*fcostheta);
			break;
			return true;
		}
	case 3:
		{
			fcostheta=zxh::VectorOP_Cosine(fglozaxis,fprj_v_xzp,3);
			fsintheta=-sqrt(1-fcostheta*fcostheta);
			break;
			return true;
		}
	case 4:
		{
			fcostheta=zxh::VectorOP_Cosine(fglozaxis,fprj_v_xzp,3);
			fsintheta=sqrt(1-fcostheta*fcostheta);
			break;
			return true;
		}
	}


	return true;
}
bool CalcRotMatrp1(int nxx,float fprj_v_xzp[3],float fglozaxis[3],char* chplane,Matrix<float,3,3> &Mrx)
{
	float fcostheta=zxh::VectorOP_Cosine(fglozaxis,fprj_v_xzp,3);
	float fsintheta=sqrt(1-fcostheta*fcostheta);
	Matrix<float,3,3> Mrx_inv;
	if(strcmp(chplane,"yz")==0)
	{
		Mrx << 1,0,0,  
			0,fcostheta,fsintheta,  
			0,-fsintheta,fcostheta;
	}
	if(strcmp(chplane,"xz")==0)
	{
		Mrx << fcostheta,0,-fsintheta,  
			0,1,0,  
			fsintheta,0,fcostheta;
	}
	if(nxx==2||nxx==3)
	{
		Mrx_inv=Mrx.inverse();
		Mrx=Mrx_inv;
	}
	

	return true;
}
bool CalcRotMatr(float fglozaxis[3],float fprj_v_p[3],char* chplane,Matrix<float,3,3> &Mrx)
{
	int nxx1=0;
	float fcostheta=1000;
	float fsintheta=1000;
	if(strcmp(chplane,"yz")==0)
	{
		float fyzp[2]={fprj_v_p[1],fprj_v_p[2]};
		nxx1=BelotoQuadr(fyzp,chplane);
		//CalcRotMatrp(nxx1,fprj_v_p,fglozaxis,fcostheta,fsintheta);
		CalcRotMatrp(nxx1,fprj_v_p,fglozaxis,fcostheta,fsintheta);
		Mrx << 1,0,0,  
			0,fcostheta,fsintheta,  
			0,-fsintheta,fcostheta;
	}
		if(strcmp(chplane,"xz")==0)
	{
		float fxzp[2]={fprj_v_p[0],fprj_v_p[2]};
		nxx1=BelotoQuadr(fxzp,chplane);
		CalcRotMatrp(nxx1,fprj_v_p,fglozaxis,fcostheta,fsintheta);
		Mrx << fcostheta,0,-fsintheta,  
			0,1,0,  
			fsintheta,0,fcostheta;
	}
		return true;
}
bool CalcRotMatr1(float fglozaxis[3],float fprj_v_p[3],char* chplane,Matrix<float,3,3> &Mrx)
{
	int nxx1=0;
	if(strcmp(chplane,"yz")==0)
	{
		float fyzp[2]={fprj_v_p[1],fprj_v_p[2]};
		nxx1=BelotoQuadr(fyzp,chplane);
		//CalcRotMatrp(nxx1,fprj_v_p,fglozaxis,fcostheta,fsintheta);
		CalcRotMatrp1(nxx1,fprj_v_p,fglozaxis,chplane,Mrx);
		return true;
	}
	if(strcmp(chplane,"xz")==0)
	{
		float fxzp[2]={fprj_v_p[0],fprj_v_p[2]};
		nxx1=BelotoQuadr(fxzp,chplane);
		//CalcRotMatrp(nxx1,fprj_v_p,fglozaxis,fcostheta,fsintheta);
		CalcRotMatrp1(nxx1,fprj_v_p,fglozaxis,chplane,Mrx);
		return true;
	}
		/*if(strcmp(chplane,"xz")==0)
	{
		float fxzp[2]={fprj_v_p[0],fprj_v_p[2]};
		nxx1=BelotoQuadr(fxzp,chplane);
		CalcRotMatrp(nxx1,fprj_v_p,fglozaxis,fcostheta,fsintheta);
		Mrx << fcostheta,0,-fsintheta,  
			0,1,0,  
			fsintheta,0,fcostheta;
	}*/
		return true;
}
bool RotaPontsToLoca1(vector<PointCordTypeDef> &vPathPointsWorld,float fabc[3],Matrix<float,1,3> &MfirPont,Matrix<float,1,2> &MPro_Pcent_2D,vector<PointCordTypeDef> &vPathRotedPointsLocWorld,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2)
{
	// the target of this function is to transform the vector v=(a, b,-1) to vt=(0, 0, 1)
	//step1: rotate v around x-axis to xz plane as v1 ; the angle alpha=arccos<v_prj_to_yz, z>, where v_prj_to_yz is the projection of v to yz plane
	// judge which quadrant
	//a*x+b*y+c*z+d=0;
	float fprj_v_yzp[3]={0,fabc[1],fabc[2]};
	float fglozaxis[3]={0,0,1};
	
	CalcRotMatr1(fglozaxis,fprj_v_yzp,"yz",Mrx1);
	//rotate around x axis
	MatrixXf Mabc(1,3),MabcT(1,3);
	Mabc<<fabc[0],fabc[1],fabc[2];
	MabcT=Mabc*Mrx1;
	//cout<<Mabc<<endl;
	//cout<<MabcT<<endl;
	//step2: rotate v1 around y-axis to z-axis
	float fprj_v_xzp[3]={MabcT(0,0),0,MabcT(0,2)};
	CalcRotMatr1(fglozaxis,fprj_v_xzp,"xz",Mrx2);
	MabcT=MabcT*Mrx2;
	
	//cout<<Mabc<<endl;
	//cout<<MabcT<<endl;

	
	//cout<<MabcT<<endl;
	//find the center point 3D
	int ncentid=0;
	int n=vPathPointsWorld.size();
	for(int i=0;i<n;i++)
	{
		if(MfirPont(0,0)==vPathPointsWorld[i].x&&MfirPont(0,1)==vPathPointsWorld[i].y&&MfirPont(0,2)==vPathPointsWorld[i].z)
			ncentid=i;
	}
	//rotate the points in the vector
	if(!vPathRotedPointsLocWorld.empty())vPathRotedPointsLocWorld.clear();

	for(int i=0;i<vPathPointsWorld.size();i++)
	{
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z;
		MnewPont=(MoriPont-MfirPont)*Mrx1*Mrx2;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		PnewPont.z=MnewPont(0,2);
		vPathRotedPointsLocWorld.push_back(PnewPont);
		//cout<<MnewPont<<endl;
	}
	MPro_Pcent_2D(0,0)=vPathRotedPointsLocWorld[ncentid].x;
	MPro_Pcent_2D(0,1)=vPathRotedPointsLocWorld[ncentid].y;

	//cout<<MfifthPontT<<endl;

	return true;
}
bool RotaPontsToLoca(vector<PointCordTypeDef> &vPathPointsWorld,float fabc[3],vector<PointCordTypeDef> &vPathRotedPointsLocWorld,Matrix<float,1,3> &MfirPont,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2)
{
	// the target of this function is to transform the vector v=(a, b,-1) to vt=(0, 0, 1)
	//step1: rotate v around x-axis to xz plane as v1 ; the angle alpha=arccos<v_prj_to_yz, z>, where v_prj_to_yz is the projection of v to yz plane
	// judge which quadrant

	float fprj_v_yzp[3]={0,fabc[1],-1};
	float fglozaxis[3]={0,0,1};
	
	CalcRotMatr(fglozaxis,fprj_v_yzp,"yz",Mrx1);
	//rotate around x axis
	MatrixXf Mabc(1,3),MabcT(1,3);
	Mabc<<fabc[0],fabc[1],-1;
	MabcT=Mabc*Mrx1;
	//step2: rotate v1 around y-axis to z-axis
	float fprj_v_xzp[3]={MabcT(0,0),0,MabcT(0,2)};
	CalcRotMatr(fglozaxis,fprj_v_xzp,"xz",Mrx2);
	MabcT=MabcT*Mrx2;
	//cout<<MabcT<<endl;
	//rotate the points in the vector
	if(!vPathRotedPointsLocWorld.empty())vPathRotedPointsLocWorld.clear();

	for(int i=0;i<vPathPointsWorld.size();i++)
	{
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z;
		MnewPont=(MoriPont-MfirPont)*Mrx1*Mrx2;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		PnewPont.z=MnewPont(0,2);
		vPathRotedPointsLocWorld.push_back(PnewPont);
		//cout<<MnewPont<<endl;
	}

	//cout<<MfifthPontT<<endl;

	return true;
}
int test_pseudoinverse(Matrix<double, 2, 2>&MA,Matrix<double, 2, 2>&MA_inv)  
{  
    //std::vector<std::vector<float>> vec{ { 0.68f, 0.597f },  
    //              { -0.211f, 0.823f },  
    //              { 0.566f, -0.605f } };  
    //const int rows{ 3 }, cols{ 2 };  
  
   /* std::vector<std::vector<float>> vec;
	vec={ { 0.68f, 0.597f, -0.211f },  { 0.823f, 0.566f, -0.605f } };  */
	std::vector<std::vector<float>> vec;
	vector<float>vvector1,vvector2;
	vvector1.push_back(MA(0,0));
	vvector1.push_back(MA(0,1));
	vvector2.push_back(MA(1,0));
	vvector2.push_back(MA(1,1));
	vec.push_back(vvector1);
	vec.push_back(vvector2);
    const int rows={ 2 }, cols={ 2 };  
  
    std::vector<float> vec_;  
    for (int i = 0; i < rows; ++i) {  
        vec_.insert(vec_.begin() + i * cols, vec[i].begin(), vec[i].end());  
    }  
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> m(vec_.data(), rows, cols);  
  
    fprintf(stderr, "source matrix:\n");  
    std::cout << m << std::endl;  
  
    fprintf(stderr, "\nEigen implement pseudoinverse:\n");  
    auto svd = m.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);  
  
    const auto &singularValues = svd.singularValues();  
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> singularValuesInv(m.cols(), m.rows());  
    singularValuesInv.setZero();  
    double  pinvtoler = 1.e-6; // choose your tolerance wisely  
    for (unsigned int i = 0; i < singularValues.size(); ++i) {  
        if (singularValues(i) > pinvtoler)  
            singularValuesInv(i, i) = 1.0f / singularValues(i);  
        else  
            singularValuesInv(i, i) = 0.f;  
    }  
  
    Eigen::MatrixXf pinvmat = svd.matrixV() * singularValuesInv * svd.matrixU().transpose(); 
	MA_inv(0,0)=pinvmat(0,0);
	MA_inv(1,1)=pinvmat(1,1);
	MA_inv(0,1)=pinvmat(0,1);
	MA_inv(1,0)=pinvmat(1,0);
    std::cout << pinvmat << std::endl;  
	std::cout << MA_inv << std::endl;  
  
    return 0;  
}  
//x=b is not consider
bool LocalRegLine(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &Pro_RotedvLocalPointsInBallWorld,float fab[3])
{
	//a  + b x - cy = 0
	//y=a+bx;
	int n=Pro_RotedvLocalPointsInBallWorld.size();
	Array<double, 1, 1000> Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mz(0,i)=0;	
		Mw(0,i)=0;
	}
	//tranform the point position to matrix
	float fmfirsyz[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=Pro_RotedvLocalPointsInBallWorld[i].x;
		My(0,i)=Pro_RotedvLocalPointsInBallWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxyz[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2))*(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2));
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=fwi;
	}
	Matrix<double, 2, 2>MA,MA_inv;
	Matrix<double, 2, 1>Mb;
	Matrix<double, 2, 1>Mab;
	
	//Ax=b
	//generate A matrx
	MA(0,0)=Mw.sum();
	MA(1,1)=((Mx.square())*Mw).sum();
	
	MA(1,0)=(Mx*Mw).sum();
	MA(0,1)=(Mx*Mw).sum();
	//generate Mb matrx
	Mb(0,0)=(My*Mw).sum();
	Mb(1,0)=(My*Mx*Mw).sum();
	//det of Ma
	float fMA_det=MA.determinant();
	if (fMA_det==0)
	{
	test_pseudoinverse(MA,MA_inv);
	}
	//result of abc
	Mab=MA_inv*Mb;
	//Mab(0,0)=a/c
	//Mab(1,0)=b/c
	fab[0]=Mab(0,0);
	fab[1]=Mab(1,0);
	fab[2]=1;
	return true;
}

bool LocalRegLine_Gener(Matrix<float,1,3> &MfirPont,float H,vector<PointCordTypeDef> &Pro_RotedvLocalPointsInBallWorld,float fabc[3])
{
	//a x + b y + c = 0
	int n=Pro_RotedvLocalPointsInBallWorld.size();
	if (n<2)
	{
		fabc[0] = 0;
       		fabc[1] = 0;
       		fabc[2] = 0;
         return false;
	}
	Array<double, 1, 1000> Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mw(0,i)=0;
	}
	//tranform the point position to matrix
	float fmfirsyz[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=Pro_RotedvLocalPointsInBallWorld[i].x;
		My(0,i)=Pro_RotedvLocalPointsInBallWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxyz[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2))*(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=fwi;
	}
	float fx_mean=(Mx*Mw).sum()/n;
	float fy_mean=(My*Mw).sum()/n;
	float fw_mean=Mw.sum()/n;
	Array<double, 1, 1000> Mx_mean,My_mean,Mw_mean;
	for(int i=0;i<1000;i++)
	{
		Mx_mean(0,i)=fx_mean;
		My_mean(0,i)=fy_mean;	
		Mw_mean(0,i)=fw_mean;
	}
	float fDxx=0,fDxy=0,fDyy=0;
	fDxx=fw_mean*(Mw*Mx*Mx).sum()-n*fx_mean*fx_mean;
	fDyy=fw_mean*(Mw*My*My).sum()-n*fy_mean*fy_mean;
	fDxy=fw_mean*(Mw*Mx*My).sum()-n*fx_mean*fy_mean;
	double lambda = ( (fDxx + fDyy) - sqrt( (fDxx - fDyy) * (fDxx - fDyy) + 4 * fDxy * fDxy) ) / 2.0;
     double den = sqrt( fDxy * fDxy + (n*fw_mean*lambda - fDxx) * (n*fw_mean*lambda - fDxx) );

	 if(fabs(den) < 1e-5)
	 {
		 if( fabs(fDxx / fDyy - 1) < 1e-5) 
		 {
			 return false;
		 }
		 else if (fy_mean==0)
		 {
			  fabc[0] = 1;
			  fabc[1] = 0;
			  fabc[2] = -fx_mean/fw_mean;
		 }
		  else if (fx_mean==0)
		 {
			  fabc[0] = 0;
			  fabc[1] = 1;
			  fabc[2] = -fy_mean/fw_mean;
		 }
	 }
	 else
	 {
		 fabc[0]=fDxy / den;
		 fabc[1]=(n*fw_mean*lambda - fDxx) / den;
		 fabc[2]= (- fabc[0] * fx_mean - fabc[1] * fy_mean)/fw_mean;
	 }
	return true;
}
bool LocalRegLine1(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &Pro_RotedvLocalPointsInBallWorld)
{
	//a x + b y + c = 0
	int n=Pro_RotedvLocalPointsInBallWorld.size();
	Array<double, 1, 1000> Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mw(0,i)=0;
	}
	//tranform the point position to matrix
	float fmfirsyz[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=Pro_RotedvLocalPointsInBallWorld[i].x;
		My(0,i)=Pro_RotedvLocalPointsInBallWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxyz[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2))*(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=1;
	}
	//coefficient matrix
	Matrix3f MA;
	MA(0,0)=Mw.sum();
	MA(1,1)=(Mw*Mx*Mx).sum();
	MA(2,2)=-(Mw*My*My).sum();

	MA(0,1)=(Mw*Mx).sum();
	MA(0,2)=-(Mw*My).sum();

	MA(1,0)=(Mw*Mx).sum();
	MA(1,2)=-(Mw*Mx*My).sum();

	MA(2,0)=(Mw*My).sum();
	MA(2,1)=(Mw*My*Mx).sum();
	float fdet=MA.determinant();
	FullPivLU<Matrix3f> lu_decomp(MA);
	int nrank=lu_decomp.rank();
    //cout << "The rank of A is " << lu_decomp.rank() << endl;
	// nrank<3 means the line is x=0 or y=0
	//then rotate the points with 45°
	float fnewcostheta=cos(M_PI/4);
	float fnewsintheta=sin(M_PI/4);
	if (nrank<3)
	{
		//find the center point 3D
		int ncentid=0;
		int n=Pro_RotedvLocalPointsInBallWorld.size();
		for(int i=0;i<n;i++)
		{
			if(MfirPont(0,0)==Pro_RotedvLocalPointsInBallWorld[i].x&&MfirPont(0,1)==Pro_RotedvLocalPointsInBallWorld[i].y)
				ncentid=i;
		}
		Matrix<float,2,2>Mrot;
		Mrot<<fnewcostheta,fnewsintheta,
			-fnewsintheta,fnewcostheta;
		for(int i=0;i<Pro_RotedvLocalPointsInBallWorld.size();i++)
		{
			MatrixXf MoriPont(1,2),MnewPont(1,2);
			MoriPont<<Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y;
			MnewPont=(MoriPont)*Mrot;
			Pro_RotedvLocalPointsInBallWorld[i].x=MnewPont(0,0);
			Pro_RotedvLocalPointsInBallWorld[i].y=MnewPont(0,1);
		}
		MfirPont<<Pro_RotedvLocalPointsInBallWorld[ncentid].x,Pro_RotedvLocalPointsInBallWorld[ncentid].y;
		return false;
	}
	return true;
}
Matrix3f Eijexchange(int i,int j)
{
	MatrixXf ME(3,3);
	ME= MatrixXf::Identity(3,3);
	RowVector3f rvtemp=ME.row(i);
	ME.row(i)=ME.row(j);
	ME.row(j)=rvtemp;
	return ME;
}
Matrix4f Eijexchange4(int i,int j)
{
	Matrix4f ME;
	ME= MatrixXf::Identity(4,4);
	RowVector4f rvtemp=ME.row(i);
	ME.row(i)=ME.row(j);
	ME.row(j)=rvtemp;
	return ME;
}
Matrix3f Eijexchange3(int i,int j)
{
	Matrix3f ME;
	ME= MatrixXf::Identity(3,3);
	RowVector3f rvtemp=ME.row(i);
	ME.row(i)=ME.row(j);
	ME.row(j)=rvtemp;
	return ME;
}
Matrix3f Eijmul(int i,float k)
{
	MatrixXf ME(3,3);
	ME= MatrixXf::Identity(3,3);
	ME.row(i)=k*ME.row(i);
	return ME;
}
Matrix4f Eijmul4(int i,float k)
{
	MatrixXf ME(4,4);
	ME= MatrixXf::Identity(4,4);
	ME.row(i)=k*ME.row(i);
	return ME;
}
Matrix3f Eijmul3(int i,float k)
{
	MatrixXf ME(3,3);
	ME= MatrixXf::Identity(3,3);
	ME.row(i)=k*ME.row(i);
	return ME;
}
Matrix3f Eijmulex(int i,int j,float k)//j*k+i
{
   MatrixXf ME(3,3);
	ME= MatrixXf::Identity(3,3);
	RowVector3f rvtemp=ME.row(i);
	ME.row(i)=ME.row(j)*k+rvtemp;
	//cout<<ME<<endl;
	return ME;
}
Matrix4f Eijmulex4(int i,int j,float k)//j*k+i
{
   Matrix4f ME;
	ME= MatrixXf::Identity(4,4);
	RowVector4f rvtemp=ME.row(i);
	ME.row(i)=ME.row(j)*k+rvtemp;
	//cout<<ME<<endl;
	return ME;
}
Matrix3f Eijmulex3(int i,int j,float k)//j*k+i
{
   Matrix3f ME;
	ME= MatrixXf::Identity(3,3);
	RowVector3f rvtemp=ME.row(i);
	ME.row(i)=ME.row(j)*k+rvtemp;
	//cout<<ME<<endl;
	return ME;
}
bool ETrans(Matrix3f &MA)
{
int maxRowE;
    double temp; 
	for(int j=0;j<=2;j++)//j=col
	{
		maxRowE=j;
		for(int i=j;i<=2;i++)//i=row
			if(fabs(MA(i,j))>fabs(MA(maxRowE,j)))
				maxRowE = i;
		if(maxRowE!=j)//find the max element
		{
			MA=Eijexchange(j,maxRowE)*MA;
		}
		if(MA(j,j)!=0)
		{
			//eleminate the non-zero element
			for(int i=j+1;i<=2;i++)
			{
				if(MA(i,j)==0)continue;
				temp =MA(i,j)/MA(j,j);
				//cout<<MA<<endl;
				MA=Eijmulex(i,j,-temp)*MA;
			}
		}
	}
	 return true;
}
bool ETrans4(Matrix4f &MA)
{
	int maxRowE;
	double temp; 
	for(int j=0;j<=3;j++)//j=col
	{
		maxRowE=j;
		for(int i=j;i<=3;i++)//i=row
		{
			if(fabsf(MA(i,j))>fabsf(MA(maxRowE,j)))
				maxRowE = i;
		}
		if(maxRowE!=j)//find the max element
		{
			MA=Eijexchange4(j,maxRowE)*MA;
		}
		if(MA(j,j)!=0)
		{
			//eleminate the non-zero element
			for(int i=j+1;i<=3;i++)
			{
				if(MA(i,j)==0)continue;
				temp =MA(i,j)/MA(j,j);
				//cout<<MA<<endl;
				MA=Eijmulex4(i,j,-temp)*MA;
			}
		}
	}
	 return true;
}
bool ETrans4_rref(Matrix4f &MA)
{
	int maxRowE;
	double temp; 
	for(int j=0;j<=3;j++)//j=col
	{
		maxRowE=j;
		for(int i=j;i<=3;i++)//i=row
		{
			if(fabsf(MA(i,j))>fabsf(MA(maxRowE,j)))
				maxRowE = i;
		}
		if(maxRowE!=j)//find the max element
		{
			MA=Eijexchange4(j,maxRowE)*MA;
		}
		if(MA(j,j)!=0)
		{
			//eleminate the non-zero element
			for(int i=j+1;i<=3;i++)
			{
				if(MA(i,j)==0)continue;
				temp =MA(i,j)/MA(j,j);
				//cout<<MA<<endl;
				MA=Eijmulex4(i,j,-temp)*MA;
			}
		}
		if(MA(j,j)!=0)
		{
		float k=1/MA(j,j);
		MA=Eijmul4(j,k)*MA;
		}
			if(MA(j,j)!=0)
		{
			//eleminate the non-zero element
			for(int i=0;i<j;i++)
			{
				if(MA(i,j)==0)continue;
				temp =MA(i,j)/MA(j,j);
				//cout<<MA<<endl;
				MA=Eijmulex4(i,j,-temp)*MA;
			}
		}
	}
	 return true;
}
bool ETrans3_rref(Matrix3f &MA)
{
	int maxRowE;
	double temp; 
	for(int j=0;j<=2;j++)//j=col
	{
		maxRowE=j;
		for(int i=j;i<=2;i++)//i=row
		{
			if(fabsf(MA(i,j))>fabsf(MA(maxRowE,j)))
				maxRowE = i;
		}
		if(maxRowE!=j)//find the max element
		{
			MA=Eijexchange3(j,maxRowE)*MA;
		}
		if(MA(j,j)!=0)
		{
			//eleminate the non-zero element
			for(int i=j+1;i<=2;i++)
			{
				if(MA(i,j)==0)continue;
				temp =MA(i,j)/MA(j,j);
				//cout<<MA<<endl;
				MA=Eijmulex3(i,j,-temp)*MA;
			}
		}
		if(MA(j,j)!=0)
		{
		float k=1/MA(j,j);
		MA=Eijmul3(j,k)*MA;
		}
			if(MA(j,j)!=0)
		{
			//eleminate the non-zero element
			for(int i=0;i<j;i++)
			{
				if(MA(i,j)==0)continue;
				temp =MA(i,j)/MA(j,j);
				//cout<<MA<<endl;
				MA=Eijmulex3(i,j,-temp)*MA;
			}
		}
	}
	 return true;
}
bool ETrans3(Matrix3f &MA)
{
	int maxRowE;
	double temp; 
	for(int j=0;j<=2;j++)//j=col
	{
		maxRowE=j;
		for(int i=j;i<=2;i++)//i=row
		{
			if(fabsf(MA(i,j))>fabsf(MA(maxRowE,j)))
				maxRowE = i;
		}
		if(maxRowE!=j)//find the max element
		{
			MA=Eijexchange3(j,maxRowE)*MA;
		}
		if(MA(j,j)!=0)
		{
			//eleminate the non-zero element
			for(int i=j+1;i<=2;i++)
			{
				if(MA(i,j)==0)continue;
				temp =MA(i,j)/MA(j,j);
				//cout<<MA<<endl;
				MA=Eijmulex3(i,j,-temp)*MA;
			}
		}
	}
	 return true;
}
bool LocalRegLine2(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &Pro_RotedvLocalPointsInBallWorld,float fabc[3])
{
	//a  + b x - cy = 0
	int n=Pro_RotedvLocalPointsInBallWorld.size();
	Array<double, 1, 1000> Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mw(0,i)=0;
	}
	//tranform the point position to matrix
	float fmfirsyz[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=Pro_RotedvLocalPointsInBallWorld[i].x;
		My(0,i)=Pro_RotedvLocalPointsInBallWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxyz[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2))*(zxh::VectorOP_Distance(fmfirsyz,fmxyz,2));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=1;
	}
	//coefficient matrix
	Matrix3f MA;
	MA(0,0)=Mw.sum();
	MA(1,1)=(Mw*Mx*Mx).sum();
	MA(2,2)=-(Mw*My*My).sum();

	MA(0,1)=(Mw*Mx).sum();
	MA(0,2)=-(Mw*My).sum();

	MA(1,0)=(Mw*Mx).sum();
	MA(1,2)=-(Mw*Mx*My).sum();

	MA(2,0)=(Mw*My).sum();
	MA(2,1)=(Mw*My*Mx).sum();
	float fdet=MA.determinant();
	FullPivLU<Matrix3f> lu_decomp(MA);
	int nrank=lu_decomp.rank();
	//b
	Matrix<float,3,1>Mabc;
	Mabc<<0,0,0;
	if(nrank==2)//elementary transformation of matrices
	{
		ETrans(MA);
		for(int i=0;i<3;i++)
		{
			if(MA(i,i)==0)
			{
				Mabc(i,0)=1;
			}

		}
		fabc[0]=Mabc(0,0);
		fabc[1]=Mabc(1,0);
		fabc[2]=Mabc(2,0);
		
	}
	else if(nrank==1)//elementary transformation of matrices
	{
		ETrans(MA);
		for(int i=0;i<3;i++)
		{
			if(fabsf(MA(i,i))<10e-5)
			{
				Mabc(i,0)=1;
			}

		}
		fabc[0]=Mabc(0,0);
		fabc[1]=Mabc(1,0);
		fabc[2]=Mabc(2,0);
	}
	else
	{
		LocalRegLine(MfirPont,H,Pro_RotedvLocalPointsInBallWorld,fabc);
	}


	
	return true;
}

 bool RotLocaRegreLineForRo(int i,float fabc[3],vector<PointCordTypeDef> &Pro_ForRovLocalPointsInBallWorld,vector<PointCordTypeDef> &Pro_ForRovRotaLocalPointsInBallWorld)
 {
	 //angle between the regression line and the positive x axis
	 //a+bx+cy=0
	 //rotate the regression line to the positive x axis

	 //cout the number of 0 values
	 int nXNUM=3;
	 for (int i=0;i<3;i++)
	 {
		 if (fabs(fabc[i])==0)nXNUM--;
	 }	
	 float fdrl[2]={0,0};
	if (nXNUM==1)
	{
		if (fabs(fabc[0])!=0)
		{
		cout<<" Wrong regression line;"<<endl;
			return false;
		}
		if (fabs(fabc[1])!=0)
		{
			fdrl[0]=0;
			fdrl[1]=1;
		}
		if (fabs(fabc[2])!=0)
		{
			fdrl[0]=1;
			fdrl[1]=0;
		}
	}
	else 
	{
		if(fabs(fabc[2])==0)
		{
			fdrl[0]=0;
			fdrl[1]=1;
		}
		if(fabs(fabc[2])!=0)
		{
			fdrl[0]=1;
			fdrl[1]=-fabc[1]/fabc[2];
		}
	}
	 float fxp[2]={1,0};
	 float fcostheta=zxh::VectorOP_Cosine(fdrl,fxp,2);
	 float fsintheta=sqrt(1-fcostheta*fcostheta);
	  Matrix<float,2,2> Mrxp,Mrxp_inv;
	  Mrxp<<fcostheta,fsintheta,
		  -fsintheta,fcostheta;
	  if (fdrl[1]>0)
	  {
		  Mrxp_inv=Mrxp.inverse();
		  Mrxp=Mrxp_inv;
	  }
	  //x1=x*Mrot
	 //generate rotation matrix
	  //and then rotate the regression line with 45°
	 Matrix<float,2,2> Mrot;
	 float fnewcostheta=cos(M_PI/4);
	 float fnewsintheta=sin(M_PI/4);
	 Mrot<<fnewcostheta,fnewsintheta,
		 -fnewsintheta,fnewcostheta;
	
	 if(!Pro_ForRovRotaLocalPointsInBallWorld.empty())Pro_ForRovRotaLocalPointsInBallWorld.clear();
	 //vector<PointCordTypeDef>Pro_v2D_x;
	for(int i=0;i<Pro_ForRovLocalPointsInBallWorld.size();i++)
	{
	    MatrixXf MoriPont(1,2),MnewPont(1,2);//MnewPontx(1,2);
		MoriPont<<Pro_ForRovLocalPointsInBallWorld[i].x,Pro_ForRovLocalPointsInBallWorld[i].y;
		MnewPont=(MoriPont)*Mrxp*Mrot;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		Pro_ForRovRotaLocalPointsInBallWorld.push_back(PnewPont);
		//
	/*	MnewPontx=(MoriPont)*Mrxp;
		PointCordTypeDef PnewPontx;
		PnewPontx.x=MnewPontx(0,0);
		PnewPontx.y=MnewPontx(0,1);
		Pro_v2D_x.push_back(PnewPontx);*/
	}

	//Pro_Localx_Points_2D.txt contains the projected points after plane fitting
	//char *Pro_Localx_Points_2D="F:/Coronary_0/code/Resutsfusion/Pro_Localx_Points_2D.txt";
	//if(i==320)
	//{
	//	WriteCA2Txt(Pro_v2D_x,Pro_Localx_Points_2D);
	//}
	 return true;
 }
  bool RotLocaRegreLineForQu1(int i,float fabc[3],vector<PointCordTypeDef> &Pro_ForRovLocalPointsInBallWorld,Matrix<float,2,2> &Mrot,vector<PointCordTypeDef> &Pro_ForRovRotaLocalPointsInBallWorld)
 {
	 //angle between the regression line and the positive x axis
	 //a+bx+cy=0
	 //rotate the regression line to the positive x axis

	 //cout the number of 0 values
	 int nXNUM=3;
	 for (int i=0;i<3;i++)
	 {
		 if (fabs(fabc[i])==0)nXNUM--;
	 }	
	 float fdrl[2]={0,0};
	if (nXNUM==1)
	{
		if (fabs(fabc[0])!=0)
		{
		cout<<" Wrong regression line;"<<endl;
			return false;
		}
		if (fabs(fabc[1])!=0)
		{
			fdrl[0]=0;
			fdrl[1]=1;
		}
		if (fabs(fabc[2])!=0)
		{
			fdrl[0]=1;
			fdrl[1]=0;
		}
	}
	else 
	{
		if(fabs(fabc[2])==0)
		{
			fdrl[0]=0;
			fdrl[1]=1;
		}
		if(fabs(fabc[2])!=0)
		{
			fdrl[0]=1;
			fdrl[1]=-fabc[1]/fabc[2];
		}
	}
	 float fxp[2]={1,0};
	 float fcostheta=zxh::VectorOP_Cosine(fdrl,fxp,2);
	 float fsintheta=sqrt(1-fcostheta*fcostheta);
	  Matrix<float,2,2> Mrxp,Mrxp_inv;
	  Mrxp<<fcostheta,fsintheta,
		  -fsintheta,fcostheta;
		 
	  if (fdrl[1]>0)
	  {
		  Mrxp_inv=Mrxp.inverse();
		  Mrxp=Mrxp_inv;
	  } 
	  Mrot=Mrxp;
	  //x1=x*Mrot
	 //generate rotation matrix
	  //and then rotate the regression line with 45°
	
	
	 if(!Pro_ForRovRotaLocalPointsInBallWorld.empty())Pro_ForRovRotaLocalPointsInBallWorld.clear();
	 //vector<PointCordTypeDef>Pro_v2D_x;
	for(int i=0;i<Pro_ForRovLocalPointsInBallWorld.size();i++)
	{
	    MatrixXf MoriPont(1,2),MnewPont(1,2);//MnewPontx(1,2);
		MoriPont<<Pro_ForRovLocalPointsInBallWorld[i].x,Pro_ForRovLocalPointsInBallWorld[i].y;
		MnewPont=(MoriPont)*Mrxp;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		Pro_ForRovRotaLocalPointsInBallWorld.push_back(PnewPont);
		//
	/*	MnewPontx=(MoriPont)*Mrxp;
		PointCordTypeDef PnewPontx;
		PnewPontx.x=MnewPontx(0,0);
		PnewPontx.y=MnewPontx(0,1);
		Pro_v2D_x.push_back(PnewPontx);*/
	}

	//Pro_Localx_Points_2D.txt contains the projected points after plane fitting
	//char *Pro_Localx_Points_2D="F:/Coronary_0/code/Resutsfusion/Pro_Localx_Points_2D.txt";
	//if(i==320)
	//{
	//	WriteCA2Txt(Pro_v2D_x,Pro_Localx_Points_2D);
	//}
	 return true;
 }
 bool RotLocaRegreLineForQu(float fabc[3],vector<PointCordTypeDef> &Pro_ForRovLocalPointsInBallWorld,Matrix<float,2,2> &Mrot,vector<PointCordTypeDef> &Pro_ForRovRotaLocalPointsInBallWorld)
 {
	 //angle between the regression line and the positive x axis
	 //a+bx-cy=0
	 float fdrl[2]={0,0};
	 if(fabc[2]<0.00001)
	 {
		 fdrl[0]=0;
		 fdrl[1]=1;
	 }
	 else
	 {
	   fdrl[0]=1;
	   fdrl[1]=fabc[1]/fabc[2];
	 }
	 float fxp[2]={1,0};
	 float fcostheta=zxh::VectorOP_Cosine(fdrl,fxp,2);
	 float fsintheta=sqrt(1-fcostheta*fcostheta);
	 if(fdrl[1]<0)
	 {
	 float fsintheta=-sqrt(1-fcostheta*fcostheta);
	 }
	 //generate rotation matrix
	 Matrix<float,2,2> Mr;
	 Mr<<fcostheta,fsintheta,
		 -fsintheta,fcostheta;
	 if(fsintheta>=0)
         Mrot=Mr.inverse();
	 else
		  Mrot=Mr;
	 if(!Pro_ForRovRotaLocalPointsInBallWorld.empty())Pro_ForRovRotaLocalPointsInBallWorld.clear();

	for(int i=0;i<Pro_ForRovLocalPointsInBallWorld.size();i++)
	{
	    MatrixXf MoriPont(1,2),MnewPont(1,2);
		MoriPont<<Pro_ForRovLocalPointsInBallWorld[i].x,Pro_ForRovLocalPointsInBallWorld[i].y;
		MnewPont=(MoriPont)*Mrot;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		Pro_ForRovRotaLocalPointsInBallWorld.push_back(PnewPont);
	}
	 return true;
 }
float CalcRp(vector<PointCordTypeDef> &vPathRotedPointsWorld)
{
	//ro=cov(X, Y)/delta_x*deta_y;
	float fsumx=0;
	float fsumy=0;
	float fcov_xy=0;
	float fdelta_x2=0;
	float fdelta_y2=0;
    for(int i=0;i<vPathRotedPointsWorld.size();i++)
	{
		fsumx=fsumx+vPathRotedPointsWorld[i].x;
		fsumy=fsumy+vPathRotedPointsWorld[i].y;
	}
	float fmeanx=fsumx/vPathRotedPointsWorld.size();
	float fmeany=fsumy/vPathRotedPointsWorld.size();
	 for(int i=0;i<vPathRotedPointsWorld.size();i++)
	{
		fcov_xy=fcov_xy+(vPathRotedPointsWorld[i].x-fmeanx)*(vPathRotedPointsWorld[i].y-fmeany);
		fdelta_x2=fdelta_x2+(vPathRotedPointsWorld[i].x-fmeanx)*(vPathRotedPointsWorld[i].x-fmeanx);
		fdelta_y2=fdelta_y2+(vPathRotedPointsWorld[i].y-fmeany)*(vPathRotedPointsWorld[i].y-fmeany);
	}
	 float feff=fcov_xy/(sqrt(fdelta_x2)*sqrt(fdelta_y2));
	 return feff;
}
bool CalcCoeMMatrix(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &vPathRotedPointsWorld,Matrix<float,3,3> &Mcoe,Matrix<float,3,1> &Mb)
{
	//y=a+bx+cx^2;
		int m=vPathRotedPointsWorld.size();
	Array<float, 1, 1000> Ax,Ay,Aw;
		//initialization
	for(int i=0;i<1000;i++)
	{
		Ax(0,i)=0;
		Ay(0,i)=0;	
		Aw(0,i)=0;	
	}
	//tranform the point position to matrix
		float fmfirsy[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<m;i++)
	{
		Ax(0,i)=vPathRotedPointsWorld[i].x;
		Ay(0,i)=vPathRotedPointsWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxy[2]={vPathRotedPointsWorld[i].x,vPathRotedPointsWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsy,fmxy,2))*(zxh::VectorOP_Distance(fmfirsy,fmxy,2));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		//float fwi=exp(-pow(fr2,2)/pow(H,2));
		Aw(0,i)=fwi;

	}
	//Mcoe*x=Mb
	//generate Mcoe matrx
	Mcoe(0,0)=Aw.sum();
	Mcoe(1,1)=(Ax*Ax*Aw).sum();
	Mcoe(2,2)=(Ax.square().square()*Aw).sum();

	Mcoe(0,1)=(Ax*Aw).sum();
	Mcoe(1,0)=(Ax*Aw).sum();

	Mcoe(0,2)=(Ax.square()*Aw).sum();
	Mcoe(2,0)=(Ax.square()*Aw).sum();

	Mcoe(1,2)=(Ax.square()*Ax*Aw).sum();
	Mcoe(2,1)=(Ax.square()*Ax*Aw).sum();
	float fdet=Mcoe.determinant();
	//generate Mb matrx
	Mb(0,0)=(Ay*Aw).sum();
	Mb(1,0)=(Ax*Ay*Aw).sum();
	Mb(2,0)=((Ax.square())*Ay*Aw).sum();
	return true;

}
float CurFit2D(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &vPathRotedPointsWorld,float fabc[3])
{
	//Mcoe*x=Mb
	//y=a+bx+cx^2
	Matrix<float,3,3> Mcoe,Mcoe_inv;
	Matrix<float,3,1> Mb,Mabc;

	CalcCoeMMatrix(MfirPont,H,vPathRotedPointsWorld,Mcoe,Mb);
	//gerate A_inverse
	//Mcoe_inv=Mcoe.inverse();
	//result of abc
	//Mabc=Mcoe_inv*Mb;
	Mabc = Mcoe.ldlt().solve(Mb);
	fabc[0]=Mabc(0,0);
	fabc[1]=Mabc(1,0);
	fabc[2]=Mabc(2,0);
		//Mean square deviation
	float fmsd=0;
	int n=vPathRotedPointsWorld.size();
	for(int i=0;i<n;i++)
	{
		float fcurpont[3]={vPathRotedPointsWorld[i].x,vPathRotedPointsWorld[i].y,vPathRotedPointsWorld[i].z};	
		float fdelta2=(abs(fabc[0]+fabc[1]*fcurpont[0]+fabc[2]*fcurpont[0]*fcurpont[0]-fcurpont[1]))*(abs(fabc[0]+fabc[1]*fcurpont[0]+fabc[2]*fcurpont[0]*fcurpont[0]-fcurpont[1]));
		fmsd=fmsd+fdelta2;
	}
	if(fmsd>0.1)
	{
		int xx=0;
	}
	//the displacement of the first new point
		
	return fabc[0];
}
bool BackRptaPontsToGlo(float ffirstpointdist,Matrix<float,2,2> &Mrot,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,Matrix<float,1,3> &MfirPont,PointCordTypeDef &PNewfirPont)
{
	MatrixXf MFirNewPontLoca2D(1,2),MFirNewPontLoca(1,3),MFirNewPontGlo(1,3);
	MFirNewPontLoca2D<<0,ffirstpointdist;//2D plane;
	Matrix<float,2,2> Mrot_inv=Mrot.inverse();
	MFirNewPontLoca2D=MFirNewPontLoca2D*Mrot_inv;//2D plane
	MFirNewPontLoca<<MFirNewPontLoca2D(0,0),MFirNewPontLoca2D(0,1),0;

	MFirNewPontGlo=MFirNewPontLoca*(Mrx2.inverse());
	//cout<<MFirNewPontGlo<<endl;
	MFirNewPontGlo=MFirNewPontGlo*(Mrx1.inverse());
	//cout<<MFirNewPontGlo<<endl;
	MFirNewPontGlo=MFirNewPontGlo+MfirPont;
	//cout<<MFirNewPontGlo<<endl;
	//dist
	float fdist=sqrt(((MFirNewPontGlo-MfirPont)*(MFirNewPontGlo-MfirPont).transpose())(0,0));
	PNewfirPont.x=MFirNewPontGlo(0,0);
	PNewfirPont.y=MFirNewPontGlo(0,1);
	PNewfirPont.z=MFirNewPontGlo(0,2);
	return true;
}

double vectordoubleSum(vector<double>::iterator first,vector<double>::size_type size)
{
    double sum=0.0;
    for(int ix=0;ix!=size;++ix){
        sum+=*first++;    
    } 
    return sum;
}

void Matrix_T(double *K,int m,int n,double *KT)//返回矩阵K的转置KT.k[m][n]

{

      int i,j,a,b;

      for (i=0,a=0;i<m;i++)

      {

            for (j=0,b=0;j<n;j++)

            {

                  KT[b+i]=K[a+j];

                  b+=m;

            }

            a+=n;

      }

}

void Matrix_Mul(double *Mul1,int Mul1_m,double *Mul2,int Mul2_n,int nm,double *Mul)

{
      //Mul1[Mul1_m][nm]*Mul2[nm][Mul2_n]=Mul即矩阵的乘法

      int i,j,k,a,b,c,d;

      for (i=0,a=0,c=0;i<Mul1_m;i++)

      {

            for (j=0;j<Mul2_n;j++)

            {

                  b=a+j;

                  Mul[b]=0;

                  for (k=0,d=0;k<nm;k++)

                  {

                        Mul[b]+=Mul1[c+k]*Mul2[d+j];

                        d+=Mul2_n;

                  }

            }

            c+=nm;

            a+=Mul2_n;

      }
}
bool Matrix_LU(double *K,int n,double *L,double *U)//对方阵K进行LU分解.分解失败返回False.成功返回True以及分解得到的L与U

{

      int i,j,a,b,c,d;

      double temp;

      for (i=0,a=0;i<n;i++)

      {

            for (j=0;j<n;j++)

            {

                  L[a+j]=U[a+j]=0;

            }

            U[a+i]=1;

            a+=n;

      }

      for (j=0,d=0;j<n;j++)

      {

            for (i=j,b=d;i<n;i++)

            {

                  temp=0;

                  a=0,c=j;

                  while (a<j)

                  {

                        temp+=L[b+a]*U[c];

                        c+=n;

                        a++;

                  }

                  L[b+j]=K[b+j]-temp;

                  b+=n;

            }

            i=j+1;

            while(i<n)

            {

                  temp=0;

                  a=0,c=i;

                  while(a<j)

                  {

                        temp+=L[d+a]*U[c];

                        a++;

                        c+=n;

                  }

                  if (L[d+j]==0)

                  {

                        return false;

                  }

                  U[d+i]=(K[d+i]-temp)/L[d+j];

                  i++;

            }

            d+=n;

      }

      return true;

}
bool Matrix_Inv(double *K,int n,double *InvK)//采用LU分解方法求方阵K的逆InvK,K[n][n]

{

      if (1==n)

      {

            if (K[0]==0)

            {

                  return false;

            }

            else

            {

                  InvK[0]=1/K[0];

            }

      }

      else if (n<1)

      {

            return false;

      }

      else

      {

            int i,j,a,b;

            double *L,*U,*d,*x,*e,temp;

            a=n*n;

            L=new double[a];

            U=new double[a];

            if (Matrix_LU(K,n,L,U))

            {

                  d=new double[n];

                  x=new double[n];

                  e=new double[n];

                  for (i=0;i<n;i++)

                  {

                        x[i]=d[i]=0;

                  }

                  for (i=0;i<n;i++)

                  {

                        for(j=0;j<n;j++)

                        {

                              e[j]=0;

                        }

                        e[i]=1;

                        j=0;

                        b=0;

                        while(j<n)

                        {

                              temp=0;

                              a=0;

                              while(a<j)

                              {

                                    temp+=d[a]*L[b+a];

                                    a++;

                              }

                              d[j]=e[j]-temp;

                              d[j]/=L[b+j];

                              j++;

                              b+=n;

                        }

                        j=n-1;

                        b-=n;

                        while(j>-1)

                        {

                              temp=0;

                              a=j+1;

                              while(a<n)

                              {

                                    temp+=U[b+a]*x[a];

                                    a++;

                              }

                              x[j]=d[j]-temp;

                              x[j]/=U[b+j];

                              j--;

                              b-=n;

                        }

                        for(j=0,b=i;j<n;j++)

                        {

                              InvK[b]=x[j];

                              b+=n;

                        }

                  }

                  delete []d;

                  delete []x;

                  delete []e;

            }

            else

            {

                  delete []L;

                  delete []U;

                  return false;

            }

            delete []L;

            delete []U;

      }

      return true;

}
bool Matrix_Solve(double *K,double *B,int m,int n,double *x)//Kx=B求解x。K[m][n]。其结果返回最小二乘解,B[m][1]

{

      double *KT,*Kmul,*Kb,*Kinv;

      int i;

      i=n*n;

      KT=new double[m*n];

      Kmul=new double[i];

      Kinv=new double[i];

      Kb=new double[n];

      Matrix_T(K,m,n,KT);

      Matrix_Mul(KT,n,K,n,m,Kmul);

      Matrix_Mul(KT,n,B,1,m,Kb);

      if (Matrix_Inv(Kmul,n,Kinv))

      {

            Matrix_Mul(Kinv,n,Kb,1,n,x);

            delete []KT;

            delete []Kmul;

            delete []Kinv;

            delete []Kb;

            return true;

      }

      else

      {

            delete []KT;

            delete []Kmul;

            delete []Kinv;

            delete []Kb;

            return false;

      }

}
bool Matrix_Solve_Swithch(vector<float> vfk,vector<float> vfb,vector<float> &vx)//Kx=B求解x。K[m][n]。其结果返回最小二乘解,B[m][1]

{
	int m=vfb.size();
	switch(m)
	{
	case 3:
		{
			Matrix<double, 3, 1>Mx,Mb;
			Mb<<vfb[0],vfb[1],vfb[2];
			Matrix<double, 3, 3>MA;
			MA<<vfk[0],vfk[1],vfk[2],
				vfk[3],vfk[4],vfk[5],
				vfk[6],vfk[7],vfk[8];
			cout<<MA<<endl;
			Mx = MA.ldlt().solve(Mb);

			for(int i=0;i<m;i++)
			{
				vx.push_back(Mx[i]);
			}
			break;
			return true;
		}
		case 2:
		{
				Matrix<double, 2, 1>Mx,Mb;
			Mb<<vfb[0],vfb[1];
			Matrix<double, 2, 2>MA;
			MA<<vfk[0],vfk[1],
				vfk[2],vfk[3];
			Mx = MA.ldlt().solve(Mb);
			for(int i=0;i<m;i++)
			{
				vx.push_back(Mx[i]);
			}
			break;
			return true;
		}
			case 1:
		{
			Matrix<double, 1, 1>Mx,Mb;
			Mb<<vfb[0];
			Matrix<double, 1, 1>MA;
			MA<<vfk[0];
			Mx = MA.ldlt().solve(Mb);
			for(int i=0;i<m;i++)
			{
				vx.push_back(Mx[i]);
			}
			break;
			return true;
		}
	}
	return false;
}
bool Solve_Eqs_Null(double *K,Matrix4f MA,int n,int nr,int np,vector<int>vp,vector<int>vnp)
{
		for(int i=0;i<n;i++)
		for(int j=0;j<np;j++)
		{
			K[i*np+j]=0;
		}
	
		if (n>nr)
		{
			//no p

			for(int i=0;i<vnp.size();i++)
			{
				int ni=vnp[i];

				for(int j=0;j<np;j++)
				{
					if(j==i)
					{
						K[ni*np+j]=1;
					}
					else
					{
						K[ni*np+j]=0;
					}
				}


			}	

			for(int i=0;i<vp.size();i++)
			{
				int ni=vp[i];
				for(int j=0;j<np;j++)
				{
						K[ni*np+j]=-MA(i,nr+j);
				}
			}

			/*
			}
			for(int i=0;i<n;i++)
			{
			for(int j=0;j<np;j++)
			{
			float xxx=K[i*n+j];
			}
			}*/
		}
		return true;
		
}
bool Solve_Eqs_Null3(double *K,Matrix3f MA,int n,int nr,int np,vector<int>vp,vector<int>vnp)
{
		for(int i=0;i<n;i++)
		for(int j=0;j<np;j++)
		{
			K[i*np+j]=0;
		}
	
		if (n>nr)
		{
			//no p

			for(int i=0;i<vnp.size();i++)
			{
				int ni=vnp[i];

				for(int j=0;j<np;j++)
				{
					if(j==i)
					{
						K[ni*np+j]=1;
					}
					else
					{
						K[ni*np+j]=0;
					}
				}


			}	

			for(int i=0;i<vp.size();i++)
			{
				int ni=vp[i];
				for(int j=0;j<np;j++)
				{
						K[ni*np+j]=-MA(i,nr+j);
				}
			}

			/*
			}
			for(int i=0;i<n;i++)
			{
			for(int j=0;j<np;j++)
			{
			float xxx=K[i*n+j];
			}
			}*/
		}
		return true;
		
}
bool PlaneFitting2(PointCordTypeDef &PfirPont,float H,vector<PointCordTypeDef> &vPathPointsWorld,float fabc[4])
{

	//a*x+b*y+c*z+d=0;
	int n=vPathPointsWorld.size();
	Array<double, 1, 1000>Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mz(0,i)=0;	
		Mw(0,i)=0;	
	}
	//tranform the point position to matrix
	float fmfirsyz[3]={PfirPont.x,PfirPont.y,PfirPont.z};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=vPathPointsWorld[i].x;
		My(0,i)=vPathPointsWorld[i].y;	
		Mz(0,i)=vPathPointsWorld[i].z;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxyz[3]={vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z};
		float fr2=(zxh::VectorOP_Distance(fmfirsyz,fmxyz,3))*(zxh::VectorOP_Distance(fmfirsyz,fmxyz,3));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=1;
	}
	Matrix4f MA,MA_copy,MA_inv;
	//Ax=b
	//generate A matrx
	MA(0,0)=(Mw*Mx*Mx).sum();
	MA(1,1)=(Mw*My*My).sum();
	MA(2,2)=(Mw*Mz*Mz).sum();
	MA(3,3)=Mw.sum();

	MA(0,1)=(Mw*Mx*My).sum();
	MA(1,0)=(Mw*Mx*My).sum();

	MA(0,2)=(Mw*Mx*Mz).sum();
	MA(2,0)=(Mw*Mx*Mz).sum();

	MA(0,3)=(Mw*Mx).sum();
	MA(3,0)=(Mw*Mx).sum();

	MA(1,2)=(Mw*My*Mz).sum();
	MA(2,1)=(Mw*My*Mz).sum();

	MA(1,3)=(Mw*My).sum();
	MA(3,1)=(Mw*My).sum();

	MA(2,3)=(Mw*Mz).sum();
	MA(3,2)=(Mw*Mz).sum();
	MA_copy=MA;
	//save to txt
//	std::string fileName = "F:/Coronary_0/code/AboutMinimalPathForExtractArtery/coronary_extraction/package_arteryextraction/package_arteryextraction/ResFusion/MA.txt" ;
//    std::ofstream outfile( fileName.c_str() ) ; // file name and the operation type. 
//
//int i, j ;
//for( i=0; i<4; i++ ){
//    for( j=0; j<4; j++ )
//          outfile<<setiosflags(ios::fixed)<<setprecision(10)<<MA(i,j) << " " ;        
//     outfile << std::endl ;       // 
//}
//outfile.close() ;
///
	//EigenSolver<Matrix4f> es(MA);
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	float fmadet=MA.determinant();
	FullPivLU<Matrix4f> lu_decomp(MA);
	int nrank=lu_decomp.rank();
	//b
	Matrix<float,4,1>Mabc;
	Mabc<<0,0,0,0;
	int nxindex[4]={-1,-1,-1,-1};

	//elementary transformation of matrices
	
		int nXnum=4;
		ETrans4(MA);
		for(int i=0;i<4;i++)
		{
			float fxxx=fabsf(MA(i,i));
			if(fabsf(MA(i,i))==0)
			{
				Mabc(i,0)=1;
				nxindex[i]=1;
				nXnum--;
			}

		}
		if(nXnum!=4)
		{
			//solve the remain x
			int m=nXnum;
			int n=nXnum;
			double *fB,*K,*x;
			fB=new double[m*1];
			K=new double[m*n];
			x=new double[m*1];
			int numbi=0;
			for(int i=0;i<m;i++)
			{
				fB[i]=0;
				x[i]=0;
			}
			for(int i=0;i<m;i++)
				for(int j=0;j<n;j++)
				{
					K[i*n+j]=0;
				}


				for(int i=0;i<4;i++)
				{
					if(nxindex[i]>0)continue;
					vector<double>vB;//calculate B
					for(int j=i+1;j<4;j++)
					{
						if(nxindex[j]<0)continue;
						float fb=-MA(i,j)*Mabc(j,0);
						vB.push_back(fb);
					}
					fB[numbi]=vectordoubleSum(vB.begin(),vB.size());
					int numbj=numbi;
					for(int j=i;j<4;j++)//calculate K matrix
					{
						if(nxindex[j]>0)continue;
						if(numbj>=numbi)
						{
							K[n*numbi+numbj]=MA(i,j);
						}
						numbj++;

					}
					numbi++;
				}
				vector<float>vfk,vfb;
				for(int i=0;i<m;i++)
					for(int j=0;j<n;j++)
					{
						vfk.push_back(K[i*n+j]);
					}

					for(int i=0;i<m;i++)
					{
						vfb.push_back(fB[i]);
					}
					//solve Kx=B;m=3,2,1

					//Matrix_Solve(K,fB,m,n,x);
					vector<float>vfx;
					Matrix_Solve_Swithch(vfk,vfb,vfx);
					
					for(int i=0;i<m;i++)
					{
						vfx.push_back(x[i]);
					}

					for(int i=0;i<4;i++)
					{

						if(nxindex[i]<0)
						{
							Mabc(i,0)=vfx.front();
							vfx.erase( vfx.begin( ) );
						}

					}

					fabc[0]=Mabc(0,0);
					fabc[1]=Mabc(1,0);
					fabc[2]=Mabc(2,0);
					fabc[3]=Mabc(3,0);
					Matrix<float,4,1>xxx;
					xxx<<0,0,0,0;
					xxx=MA_copy*Mabc;
					int xxxx=0;
		}
		else
		{
			//a*x+b*y+c*z+d=0;
			PlaneFitting(PfirPont,H,vPathPointsWorld,fabc);
		}
		//Mean square deviation
		//a*x+b*y+c*z+d=0;
		float fmsd=0;
		float a=fabc[0];
		float b=fabc[1];
		float c=fabc[2];
		float d=fabc[3];
		for(int i=0;i<n;i++)
		{
			float fcurpont[3]={vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z};	
			float f_dist_Curvefit=fabsf(fabc[0]*fcurpont[0]+fabc[1]*fcurpont[1]+fabc[2]*fcurpont[2]+fabc[3])/sqrt(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
			fmsd=fmsd+f_dist_Curvefit;
		}
		if(fmsd>0.1)
		{
			int xx=0;
		}

	return true;
}
bool PlaneFitting3(PointCordTypeDef &PfirPont,float H,vector<PointCordTypeDef> &vPathPointsWorld,float fabc[4])
{

	//a*x+b*y+c*z+d=0;
	int n=vPathPointsWorld.size();
	Array<double, 1, 1000>Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mz(0,i)=0;	
		Mw(0,i)=0;	
	}
	//tranform the point position to matrix
	float fmfirsyz[3]={PfirPont.x,PfirPont.y,PfirPont.z};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=vPathPointsWorld[i].x;
		My(0,i)=vPathPointsWorld[i].y;	
		Mz(0,i)=vPathPointsWorld[i].z;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxyz[3]={vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z};
		float fr2=(zxh::VectorOP_Distance(fmfirsyz,fmxyz,3))*(zxh::VectorOP_Distance(fmfirsyz,fmxyz,3));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=1;
	}
	Matrix4f MA,MA_copy,MA_inv;
	//Ax=b
	//generate A matrx
	MA(0,0)=(Mw*Mx*Mx).sum();
	MA(1,1)=(Mw*My*My).sum();
	MA(2,2)=(Mw*Mz*Mz).sum();
	MA(3,3)=Mw.sum();

	MA(0,1)=(Mw*Mx*My).sum();
	MA(1,0)=(Mw*Mx*My).sum();

	MA(0,2)=(Mw*Mx*Mz).sum();
	MA(2,0)=(Mw*Mx*Mz).sum();

	MA(0,3)=(Mw*Mx).sum();
	MA(3,0)=(Mw*Mx).sum();

	MA(1,2)=(Mw*My*Mz).sum();
	MA(2,1)=(Mw*My*Mz).sum();

	MA(1,3)=(Mw*My).sum();
	MA(3,1)=(Mw*My).sum();

	MA(2,3)=(Mw*Mz).sum();
	MA(3,2)=(Mw*Mz).sum();
	MA_copy=MA;
	//save to txt
//	std::string fileName = "F:/Coronary_0/code/AboutMinimalPathForExtractArtery/coronary_extraction/package_arteryextraction/package_arteryextraction/ResFusion/MA.txt" ;
//    std::ofstream outfile( fileName.c_str() ) ; // file name and the operation type. 
//
//int i, j ;
//for( i=0; i<4; i++ ){
//    for( j=0; j<4; j++ )
//          outfile<<setiosflags(ios::fixed)<<setprecision(10)<<MA(i,j) << " " ;        
//     outfile << std::endl ;       // 
//}
//outfile.close() ;
///
	//EigenSolver<Matrix4f> es(MA);
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	float fmadet=MA.determinant();
	FullPivLU<Matrix4f> lu_decomp(MA);
	int nrank=lu_decomp.rank();
	//b
	Matrix<float,4,1>Mabc;
	Mabc<<0,0,0,0;
	int nxindex[4]={-1,-1,-1,-1};

	//elementary transformation of matrices
	
		int nXnum=4;
		ETrans4_rref(MA);
		for(int i=0;i<4;i++)
		{
			float fxxx=fabsf(MA(i,i));
			if(fabsf(MA(i,i))==0)
			{
				nxindex[i]=1;
				nXnum--;
			}

		}
		if(nXnum!=4)
		{
			//solve the remain x
			int m=MA.rows();
			int n=MA.cols();
			int np=n-nXnum;
			double *K;
			K=new double[n*np];


			vector<int> vp;
			vector<int> vnp;
			for(int i=0;i<n;i++)
			{
				int nc=nxindex[i];
				if(nc<0)
				{
					vp.push_back(i);
				}
				else
				{
					vnp.push_back(i);
				}
			}
			Solve_Eqs_Null(K,MA,n,nXnum,np,vp,vnp);
			for(int i=0;i<n;i++)
			{
				float sumKj=0;
				for(int j=0;j<np;j++)
				{
					sumKj=sumKj+K[np*i+j];
				}
				Mabc(i,0)=sumKj;
			}
			Matrix<float,4,1>xxx;
					xxx<<0,0,0,0;
					xxx=MA_copy*Mabc;
					int xxxx=0;
			fabc[0]=Mabc(0,0);
			fabc[1]=Mabc(1,0);
			fabc[2]=Mabc(2,0);
			fabc[3]=Mabc(3,0);
		}

		else
		{
			//a*x+b*y+c*z+d=0;
			PlaneFitting(PfirPont,H,vPathPointsWorld,fabc);
		}
		//Mean square deviation
		//a*x+b*y+c*z+d=0;
		float fmsd=0;
		float a=fabc[0];
		float b=fabc[1];
		float c=fabc[2];
		float d=fabc[3];
		for(int i=0;i<n;i++)
		{
			float fcurpont[3]={vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z};	
			float f_dist_Curvefit=fabsf(fabc[0]*fcurpont[0]+fabc[1]*fcurpont[1]+fabc[2]*fcurpont[2]+fabc[3])/sqrt(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
			fmsd=fmsd+f_dist_Curvefit;
		}
		if(fmsd>0.1)
		{
			int xx=0;
		}

	return true;
}
bool LocalRegLine4(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &Pro_RotedvLocalPointsInBallWorld,float fab[3])
{
	//a  + b x + cy = 0
	//y=a+bx;
	int n=Pro_RotedvLocalPointsInBallWorld.size();
	Array<double, 1, 1000> Mx,My,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mw(0,i)=0;
	}
	//tranform the point position to matrix
	float fmfirsy[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=Pro_RotedvLocalPointsInBallWorld[i].x;
		My(0,i)=Pro_RotedvLocalPointsInBallWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxy[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsy,fmxy,2))*(zxh::VectorOP_Distance(fmfirsy,fmxy,2));
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=fwi;
	}
	Matrix<double, 2, 2>MA,MA_inv;
	Matrix<double, 2, 1>Mb;
	Matrix<double, 2, 1>Mab;
	
	//Ax=b
	//generate A matrx
	MA(0,0)=Mw.sum();
	MA(1,1)=((Mx.square())*Mw).sum();
	
	MA(1,0)=(Mx*Mw).sum();
	MA(0,1)=(Mx*Mw).sum();
	//generate Mb matrx
	Mb(0,0)=(My*Mw).sum();
	Mb(1,0)=(My*Mx*Mw).sum();
	//det of Ma
	//float fMA_det=MA.determinant();
	//if (fMA_det==0)
	//{
	//test_pseudoinverse(MA,MA_inv);
	//}
	//result of abc
	//Mab=MA_inv*Mb;
	//Mab(0,0)=a/c
	//Mab(1,0)=b/c

	Mab = MA.ldlt().solve(Mb);
	fab[0]=Mab(0,0);
	fab[1]=Mab(1,0);
	fab[2]=-1;
	return true;
}
bool LocalRegLine3(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &Pro_RotedvLocalPointsInBallWorld,float fabc[3])
{
	//a  + b x + cy = 0
	int n=Pro_RotedvLocalPointsInBallWorld.size();
	Array<double, 1, 1000> Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mw(0,i)=0;
	}
	//tranform the point position to matrix
	float fmfirsy[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=Pro_RotedvLocalPointsInBallWorld[i].x;
		My(0,i)=Pro_RotedvLocalPointsInBallWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxy[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsy,fmxy,2))*(zxh::VectorOP_Distance(fmfirsy,fmxy,2));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=fwi;
	}
	//coefficient matrix
	Matrix3f MA,MA_copy;
	MA(0,0)=Mw.sum();
	MA(1,1)=(Mw*Mx*Mx).sum();
	MA(2,2)=(Mw*My*My).sum();

	MA(0,1)=(Mw*Mx).sum();
	MA(0,2)=(Mw*My).sum();

	MA(1,0)=(Mw*Mx).sum();
	MA(1,2)=(Mw*Mx*My).sum();

	MA(2,0)=(Mw*My).sum();
	MA(2,1)=(Mw*My*Mx).sum();
	MA_copy=MA;
	float fdet=MA.determinant();
	//b
	Matrix<float,3,1>Mabc;
	Mabc<<0,0,0;
	int nxindex[3]={-1,-1,-1};

	//elementary transformation of matrices

	int nXnum=3;
	ETrans3(MA);
	for(int i=0;i<3;i++)
	{
		float fxxx=fabsf(MA(i,i));
		if(fabsf(MA(i,i))==0)
		{
			Mabc(i,0)=1;
			nxindex[i]=1;
			nXnum--;
		}
	}
	if (nXnum==1&&nxindex[0]<0)
	{
		cout<<"Can not regress the line;"<<endl;
			return false;
	}
	if(nXnum!=3)
	{
		//solve the remain x
		int m=nXnum;
		int n=nXnum;
		double *fB,*K,*x;
		fB=new double[m*1];
		K=new double[m*n];
		x=new double[m*1];
		int numbi=0;
		for(int i=0;i<m;i++)
		{
			fB[i]=0;
			x[i]=0;
		}
		for(int i=0;i<m;i++)
			for(int j=0;j<n;j++)
			{
				K[i*n+j]=0;
			}
			for(int i=0;i<3;i++)
			{
				if(nxindex[i]>0)continue;
				vector<double>vB;//calculate B
				for(int j=i+1;j<3;j++)
				{
					if(nxindex[j]<0)continue;
					float fb=-MA(i,j)*Mabc(j,0);
					vB.push_back(fb);
				}
				fB[numbi]=vectordoubleSum(vB.begin(),vB.size());
				int numbj=numbi;
				for(int j=i;j<3;j++)//calculate K matrix
				{
					if(nxindex[j]>0)continue;
					if(numbj>=numbi)
					{
						K[n*numbi+numbj]=MA(i,j);
					}
					numbj++;

				}
				numbi++;
			}
			vector<float>vfk,vfb;
			for(int i=0;i<m;i++)
				for(int j=0;j<n;j++)
				{
					vfk.push_back(K[i*n+j]);
				}

				for(int i=0;i<m;i++)
				{
					vfb.push_back(fB[i]);
				}

				vector<float>vfx;
				Matrix_Solve_Swithch(vfk,vfb,vfx);
				for(int i=0;i<m;i++)
				{
					vfx.push_back(x[i]);
				}
				delete []K;
				delete []fB;
				delete []x;



				for(int i=0;i<3;i++)
				{

					if(nxindex[i]<0)
					{
						Mabc(i,0)=vfx.front();
						vfx.erase(vfx.begin( ));
					}

				}

				fabc[0]=Mabc(0,0);
				fabc[1]=Mabc(1,0);
				fabc[2]=Mabc(2,0);
				Matrix<float,3,1>xxx;
				xxx<<0,0,0;
				xxx=MA_copy*Mabc;
	}
	else
	{
		LocalRegLine4(MfirPont,H,Pro_RotedvLocalPointsInBallWorld,fabc);
	}
	//Mean square deviation
	//a  + b x + cy = 0
	float fmsd=0;
	for(int i=0;i<n;i++)
	{
		float fcurpont[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};	
		float f_dist_Curvefit=fabsf(fabc[0]+fabc[1]*fcurpont[0]+fabc[2]*fcurpont[1])/sqrt(fabc[1]*fabc[1]+fabc[2]*fabc[2]);
		fmsd=fmsd+f_dist_Curvefit;
	}
	if(fmsd>0.1)
	{
		int xx=0;
	}
	return true;
}

bool LocalRegLine5(Matrix<float,1,2> &MfirPont,float H,vector<PointCordTypeDef> &Pro_RotedvLocalPointsInBallWorld,float fabc[3])
{
	//a  + b x + cy = 0
	int n=Pro_RotedvLocalPointsInBallWorld.size();
	Array<double, 1, 1000> Mx,My,Mz,Mw;
	//initialization
	for(int i=0;i<1000;i++)
	{
		Mx(0,i)=0;
		My(0,i)=0;	
		Mw(0,i)=0;
	}
	//tranform the point position to matrix
	float fmfirsy[2]={MfirPont(0,0),MfirPont(0,1)};
	for(int i=0;i<n;i++)
	{
		Mx(0,i)=Pro_RotedvLocalPointsInBallWorld[i].x;
		My(0,i)=Pro_RotedvLocalPointsInBallWorld[i].y;	
		//calculate the weight
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		float fmxy[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};
		float fr2=(zxh::VectorOP_Distance(fmfirsy,fmxy,2))*(zxh::VectorOP_Distance(fmfirsy,fmxy,2));
		//w_i=(2*r^3/H^3-3*r^2/H^2+1)
		//float fwi=2*pow(fr2,3)/pow(H,3)-3*pow(fr2,2)/pow(H,2)+1;
		//w_i=exp(-pow(r,2)/pow(H,2))
		float fwi=exp(-pow(fr2,2)/pow(H,2));
		Mw(0,i)=fwi;
	}
	//coefficient matrix
	Matrix3f MA,MA_copy;
	MA(0,0)=Mw.sum();
	MA(1,1)=(Mw*Mx*Mx).sum();
	MA(2,2)=(Mw*My*My).sum();

	MA(0,1)=(Mw*Mx).sum();
	MA(0,2)=(Mw*My).sum();

	MA(1,0)=(Mw*Mx).sum();
	MA(1,2)=(Mw*Mx*My).sum();

	MA(2,0)=(Mw*My).sum();
	MA(2,1)=(Mw*My*Mx).sum();
	MA_copy=MA;
	float fdet=MA.determinant();
	//b
	Matrix<float,3,1>Mabc;
	Mabc<<0,0,0;
	int nxindex[3]={-1,-1,-1};

	//elementary transformation of matrices

	int nXnum=3;
	ETrans3_rref(MA);
	for(int i=0;i<3;i++)
	{
		float fxxx=fabsf(MA(i,i));
		if(fabsf(MA(i,i))==0)
		{
			Mabc(i,0)=1;
			nxindex[i]=1;
			nXnum--;
		}
	}
	if (nXnum==1&&nxindex[0]<0)
	{
		cout<<"Can not regress the line;"<<endl;
			return false;
	}
	if(nXnum!=3)
	{
	//solve the remain x
			int m=MA.rows();
			int n=MA.cols();
			int np=n-nXnum;
			double *K;
			K=new double[n*np];


			vector<int> vp;
			vector<int> vnp;
			for(int i=0;i<n;i++)
			{
				int nc=nxindex[i];
				if(nc<0)
				{
					vp.push_back(i);
				}
				else
				{
					vnp.push_back(i);
				}
			}
			Solve_Eqs_Null3(K,MA,n,nXnum,np,vp,vnp);
			for(int i=0;i<n;i++)
			{
				float sumKj=0;
				for(int j=0;j<np;j++)
				{
					sumKj=sumKj+K[np*i+j];
				}
				Mabc(i,0)=sumKj;
			}
			Matrix<float,3,1>xxx;
					xxx<<0,0,0;
					xxx=MA_copy*Mabc;
					int xxxx=0;
			fabc[0]=Mabc(0,0);
			fabc[1]=Mabc(1,0);
			fabc[2]=Mabc(2,0);
	}
	else
	{
		LocalRegLine4(MfirPont,H,Pro_RotedvLocalPointsInBallWorld,fabc);
	}
	//Mean square deviation
	//a  + b x + cy = 0
	float fmsd=0;
	for(int i=0;i<n;i++)
	{
		float fcurpont[2]={Pro_RotedvLocalPointsInBallWorld[i].x,Pro_RotedvLocalPointsInBallWorld[i].y};	
		float f_dist_Curvefit=fabsf(fabc[0]+fabc[1]*fcurpont[0]+fabc[2]*fcurpont[1])/sqrt(fabc[1]*fabc[1]+fabc[2]*fabc[2]);
		fmsd=fmsd+f_dist_Curvefit;
	}
	if(fmsd>0.1)
	{
		int xx=0;
	}
	return true;
}
int FindPosi(PointCordTypeDef PCent,std::vector<jdq2017::point3D> &path)
{
	float fpxyz[3]={PCent.x,PCent.y,PCent.z};
	float fmindist=100000;
	int curIdx=0;
	for (int i=0;i< int(path.size())-1;i++)
	{
		float fxyz[3]={path[i]._x,path[i]._y,path[i]._z};
		float fpxyzdist=zxh::VectorOP_Distance(fxyz,fpxyz,3);
		if (fpxyzdist<fmindist)
		{
			fmindist=fpxyzdist;
			curIdx=i;
		}
		if(fxyz[0]==fpxyz[0]&&fxyz[1]==fpxyz[1]&&fxyz[2]==fpxyz[2])
		{
			return i;
		}
	} 
	if(fmindist<0.01)
	{
		return curIdx;
	}
	return -1;
  
}
bool FindPosi_Y_ALL(PointCordTypeDef PCent,vector<vLinesDef> &vnLine,vector<pair<int,int>>&vLP)
{
	float fpxyz[3]={PCent.x,PCent.y,PCent.z};
	float fmindist=100000;
	int curIdx=0;
	if(!vLP.empty())vLP.clear();
	for(int i=0;i<vnLine.size();i++)
	{
		vLinesDef vLtmp;
		vLtmp=vnLine[i];

		for (int j=0;j<(vLtmp.vLine).size() ;j++)
		{ 
			std::vector<jdq2017::point3D> vPtmp;
			vPtmp=vLtmp.vLine;
			float fxyz[3]={vPtmp[j]._x,vPtmp[j]._y,vPtmp[j]._z};
			if(fxyz[0]==fpxyz[0]&&fxyz[1]==fpxyz[1]&&fxyz[2]==fpxyz[2])
			{
				pair<int,int>pij;
				pij=make_pair(i,j);
				vLP.push_back(pij);
				return true;
			}
		} 
	}
	return false;
  
}
int FindPosi_Y(PointCordTypeDef PCent,std::vector<jdq2017::point3D> &path)
{
	float fpxyz[3]={PCent.x,PCent.y,PCent.z};
	float fmindist=100000;
	int curIdx=0;
	for (int i=0;i< int(path.size())-1;i++)
	{
		float fxyz[3]={path[i]._x,path[i]._y,path[i]._z};
		if(fxyz[0]==fpxyz[0]&&fxyz[1]==fpxyz[1]&&fxyz[2]==fpxyz[2])
		{
			return i;
		}
	} 
	return -1;
  
}
bool RotaPontsToLoca12(vector<PointCordTypeDef> &vPathPointsWorld,float fabc[3],Matrix<float,1,3> &MfirPont,vector<PointCordTypeDef> &vPathRotedPointsLocWorld,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,Matrix<float,3,3> &Mrx3)
{
	// the target of this function is to transform the vector v=(a, b,-1) to vt=(0, 0, 1)
	//step1: rotate v around x-axis to xz plane as v1 ; the angle alpha=arccos<v_prj_to_yz, z>, where v_prj_to_yz is the projection of v to yz plane
	// judge which quadrant
	//a*x+b*y+c*z+d=0;
	float fprj_v_yzp[3]={0,fabc[1],fabc[2]};
	float fglozaxis[3]={0,0,1};
	
	CalcRotMatr1(fglozaxis,fprj_v_yzp,"yz",Mrx1);
	//rotate around x axis
	MatrixXf Mabc(1,3),MabcT(1,3),MabcTT(1,3);
	Mabc<<fabc[0],fabc[1],fabc[2];
	MabcT=Mabc*Mrx1;
	//cout<<Mabc<<endl;
	//cout<<MabcT<<endl;
	//step2: rotate v1 around y-axis to z-axis
	float fprj_v_xzp[3]={MabcT(0,0),0,MabcT(0,2)};
	CalcRotMatr1(fglozaxis,fprj_v_xzp,"xz",Mrx2);
	MabcT=MabcT*Mrx2;
	////step3: rotate v1 around y-axis to z-axis
	//cout<<MabcT<<endl;
	float fcostheta3=cos(M_PI/2);
	float fsintheta3=sin(M_PI/2);
	Mrx3 << fcostheta3,0,-fsintheta3,  
			0,1,0,  
			fsintheta3,0,fcostheta3;
	MabcTT=Mabc*Mrx1*Mrx2*Mrx3;


	//rotate the points in the vector
	if(!vPathRotedPointsLocWorld.empty())vPathRotedPointsLocWorld.clear();

	for(int i=0;i<vPathPointsWorld.size();i++)
	{
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z;
		MnewPont=(MoriPont-MfirPont)*Mrx1*Mrx2*Mrx3;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		PnewPont.z=MnewPont(0,2);
		vPathRotedPointsLocWorld.push_back(PnewPont);
		//cout<<MnewPont<<endl;
	}
	

	//cout<<MfifthPontT<<endl;

	return true;
}
bool RotaPontsToLoca12_ALL(int k,vector<PointCordTypeDef> &vPathPointsWorld,float fabc[3],Matrix<float,1,3> &MfirPont,vector<PointCordTypeDef> &vPathRotedPointsLocWorld,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,Matrix<float,3,3> &Mrx3)
{
	// the target of this function is to transform the vector v=(a, b,-1) to vt=(0, 0, 1)
	//step1: rotate v around x-axis to xz plane as v1 ; the angle alpha=arccos<v_prj_to_yz, z>, where v_prj_to_yz is the projection of v to yz plane
	// judge which quadrant
	//a*x+b*y+c*z+d=0;
	float fprj_v_yzp[3]={0,fabc[1],fabc[2]};
	float fglozaxis[3]={0,0,1};
	
	CalcRotMatr1(fglozaxis,fprj_v_yzp,"yz",Mrx1);
	//rotate around x axis
	MatrixXf Mabc(1,3),MabcT(1,3),MabcTT(1,3);
	Mabc<<fabc[0],fabc[1],fabc[2];
	MabcT=Mabc*Mrx1;
	//cout<<Mabc<<endl;
	//cout<<MabcT<<endl;
	//step2: rotate v1 around y-axis to z-axis
	float fprj_v_xzp[3]={MabcT(0,0),0,MabcT(0,2)};
	CalcRotMatr1(fglozaxis,fprj_v_xzp,"xz",Mrx2);
	MabcT=MabcT*Mrx2;
	////step3: rotate v1 around y-axis to z-axis
	//cout<<MabcT<<endl;
	float fcostheta3=cos(M_PI/2);
	float fsintheta3=sin(M_PI/2);
	Mrx3 << fcostheta3,0,-fsintheta3,  
			0,1,0,  
			fsintheta3,0,fcostheta3;
	MabcTT=Mabc*Mrx1*Mrx2*Mrx3;

	//rotate the points in the vector
	if(!vPathRotedPointsLocWorld.empty())vPathRotedPointsLocWorld.clear();

	for(int i=0;i<vPathPointsWorld.size();i++)
	{
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<vPathPointsWorld[i].x,vPathPointsWorld[i].y,vPathPointsWorld[i].z;
		MnewPont=(MoriPont-MfirPont)*Mrx1*Mrx2*Mrx3;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		PnewPont.z=MnewPont(0,2);
		vPathRotedPointsLocWorld.push_back(PnewPont);
		//cout<<MnewPont<<endl;
	}
	

	//cout<<MfifthPontT<<endl;

	return true;
}
bool Collet2(PointCordTypeDef PCent,float H,vector<PointCordTypeDef> &vUnorgaPointsWorld,float fRo,int nINUM,std::vector<jdq2017::point3D> ref,std::vector<jdq2017::point3D> cl,vector<PointCordTypeDef> &vLocalPointsInBallWorld)
{
	//find the number of PCent in ref
	int nRefNum=FindPosi_Y(PCent,ref);
	int nSamNUM=500;
	double refSampling = pathLength(ref)/nSamNUM;
	float fv[3]={0,0,0};

	if (nRefNum>=0&&nRefNum<ref.size()-1)//
	{
		int nRefPNum=H/refSampling;
		int nRefNext=nRefNum+nRefPNum;
		
		nRefNext=zxh::minf(ref.size()-1,nRefNext);
		fv[0]=ref[nRefNext]._x-ref[nRefNum]._x;
		fv[1]=ref[nRefNext]._y-ref[nRefNum]._y;
		fv[2]=ref[nRefNext]._z-ref[nRefNum]._z;
	}
	if (nRefNum==ref.size()-1)//generate a ellipsoid
	{
		int nRefNext=ref.size()-1;
		int nRefPNum=H/refSampling;
		int nRefNum=nRefNext-nRefPNum;
		
		fv[0]=ref[nRefNext]._x-ref[nRefNum]._x;
		fv[1]=ref[nRefNext]._y-ref[nRefNum]._y;
		fv[2]=ref[nRefNext]._z-ref[nRefNum]._z;
	}
	//find the number of PCent in cl
	int nCLNum=FindPosi_Y(PCent,cl);
	double clSampling = pathLength(cl)/nSamNUM;
		if (nCLNum>=0&&nCLNum<cl.size()-1)//generate a ellipsoid
	{
		int nRefPNum=H/clSampling;
		int nCLNext=nCLNum+nRefPNum;
		nCLNext=zxh::minf(cl.size()-1,nCLNext);
		fv[0]=cl[nCLNext]._x-cl[nCLNum]._x;
		fv[1]=cl[nCLNext]._y-cl[nCLNum]._y;
		fv[2]=cl[nCLNext]._z-cl[nCLNum]._z;
	}
		if (nCLNum==cl.size()-1)//generate a ellipsoid
	{
		int nCLNext=cl.size()-1;
		int nRefPNum=H/refSampling;
		int nCLNum=nCLNext-nRefPNum;
		
		fv[0]=cl[nCLNext]._x-cl[nCLNum]._x;
		fv[1]=cl[nCLNext]._y-cl[nCLNum]._y;
		fv[2]=cl[nCLNext]._z-cl[nCLNum]._z;
	}
	
	//rotate the fv vector to the x-positive
	vector<PointCordTypeDef> RotedvPoints,vUnorgaPointsWorld_Elli;
	vUnorgaPointsWorld_Elli.assign(vUnorgaPointsWorld.begin(),vUnorgaPointsWorld.end());
	Matrix<float,3,3> Mrx1,Mrx2,Mrx3;
	if(!RotedvPoints.empty())RotedvPoints.clear();
	Matrix<float,1,3> MPCent;
	MPCent<<PCent.x,PCent.y,PCent.z;
	RotaPontsToLoca12(vUnorgaPointsWorld_Elli,fv,MPCent,RotedvPoints,Mrx1,Mrx2,Mrx3);	
	//Ellipe
	float falph=1+nINUM*0.1;
	float fd=zxh::minf(falph*(1-fRo*fRo),4*H*H-0.5);
	float fb=(4*H*H-fd*fd)/4/H;
	float fa=2*H-fb;
	//Within the Ellipe
	for(int i=0;i<RotedvPoints.size();i++)
	{
		float fxyz[3]={RotedvPoints[i].x,RotedvPoints[i].y,RotedvPoints[i].z};
		float fmxyz=(fxyz[0]*fxyz[0]/fa/fa)+(fxyz[1]*fxyz[1]+fxyz[2]*fxyz[2])/fb/fb;
		if(fmxyz<=1)
		{

			Matrix<float,3,3> Mrx1_inv,Mrx2_inv,Mrx3_inv;
			MatrixXf MrotPont(1,3),MnewPont(1,3);
		    MrotPont<<RotedvPoints[i].x,RotedvPoints[i].y,RotedvPoints[i].z;
			Mrx1_inv=Mrx1.inverse();
			Mrx2_inv=Mrx2.inverse();
			Mrx3_inv=Mrx3.inverse();
			MnewPont=MrotPont*Mrx3_inv*Mrx2_inv*Mrx1_inv+MPCent;
			//MnewPont=(MoriPont-MPCent)*Mrx1*Mrx2*Mrx3;
			PointCordTypeDef PnewPont;
			PnewPont.x=MnewPont(0,0);
			PnewPont.y=MnewPont(0,1);
			PnewPont.z=MnewPont(0,2);
			vLocalPointsInBallWorld.push_back(PnewPont);
		}
	}

	return true;
}
bool Collet2_ALL(int k,float clSampling,PointCordTypeDef PCent,float H,vector<PointCordTypeDef> &vUnorgaPointsWorld,float fRo,int nINUM,vector<vLinesDef> &vnLine,vector<PointCordTypeDef> &vLocalPointsInBallWorld)
{
	//find the number of PCent in all lines
	vector<pair<int,int>>vLP;
	FindPosi_Y_ALL(PCent,vnLine,vLP);
	float fv[3]={0,0,0};
	for(int i=0;i<1;i++)
	{
		int nPNum=vLP[i].second;
		int nLNUM=vLP[i].first;
		   std::vector<jdq2017::point3D> ref;
			ref=vnLine[nLNUM].vLine;
		int nRefSiz=ref.size();
		
		
		if (nPNum>=0&&nPNum<nRefSiz-1)//
		{
			int nPPNum=H/clSampling;
			int nPNext=nPNum+nPPNum;
			nPNext=zxh::minf(nRefSiz-1,nPNext);
			
			fv[0]=ref[nPNext]._x-ref[nPNum]._x;
			fv[1]=ref[nPNext]._y-ref[nPNum]._y;
			fv[2]=ref[nPNext]._z-ref[nPNum]._z;
		}
		if (nPNum==ref.size()-1)//generate a ellipsoid
		{
			int nPNext=ref.size()-1;
			int nRefPNum=H/clSampling;
			int nPNum=nPNext-nRefPNum;

			fv[0]=ref[nPNext]._x-ref[nPNum]._x;
			fv[1]=ref[nPNext]._y-ref[nPNum]._y;
			fv[2]=ref[nPNext]._z-ref[nPNum]._z;
		}
	}
		//rotate the fv vector to the x-positive
	vector<PointCordTypeDef> RotedvPoints,vUnorgaPointsWorld_Elli;
	vUnorgaPointsWorld_Elli.assign(vUnorgaPointsWorld.begin(),vUnorgaPointsWorld.end());
	Matrix<float,3,3> Mrx1,Mrx2,Mrx3;
	if(!RotedvPoints.empty())RotedvPoints.clear();
	Matrix<float,1,3> MPCent;
	MPCent<<PCent.x,PCent.y,PCent.z;
	RotaPontsToLoca12_ALL(k,vUnorgaPointsWorld_Elli,fv,MPCent,RotedvPoints,Mrx1,Mrx2,Mrx3);	
	//Ellipe
	float falph=H+0.1*nINUM;
	falph=zxh::minf(falph,4*H);
	float fd=zxh::minf(falph*(1-fRo*fRo),sqrt(4*H*H-4*H*0.1));
	float fb=(4*H*H-fd*fd)/4/H;
	float fa=2*H-fb;
	//Within the Ellipe
	for(int i=0;i<RotedvPoints.size();i++)
	{
		float fxyz[3]={RotedvPoints[i].x,RotedvPoints[i].y,RotedvPoints[i].z};
		float fmxyz=(fxyz[0]*fxyz[0]/fa/fa)+(fxyz[1]*fxyz[1]+fxyz[2]*fxyz[2])/fb/fb;
		if(fmxyz<=1)
		{

			Matrix<float,3,3> Mrx1_inv,Mrx2_inv,Mrx3_inv;
			MatrixXf MrotPont(1,3),MnewPont(1,3);
		    MrotPont<<RotedvPoints[i].x,RotedvPoints[i].y,RotedvPoints[i].z;
			Mrx1_inv=Mrx1.inverse();
			Mrx2_inv=Mrx2.inverse();
			Mrx3_inv=Mrx3.inverse();
			MnewPont=MrotPont*Mrx3_inv*Mrx2_inv*Mrx1_inv+MPCent;
			//MnewPont=(MoriPont-MPCent)*Mrx1*Mrx2*Mrx3;
			PointCordTypeDef PnewPont;
			PnewPont.x=MnewPont(0,0);
			PnewPont.y=MnewPont(0,1);
			PnewPont.z=MnewPont(0,2);
			vLocalPointsInBallWorld.push_back(PnewPont);
		}
	}
	if(vLocalPointsInBallWorld.size()<2)
	{
		return false;
	}
	return true;
}
float DeterIfBranPont(PointCordTypeDef PCent,std::vector<jdq2017::point3D> &ref,std::vector<jdq2017::point3D> &cl)
{
	double ddist=2;
	int nSamNUM=500;
	//ref
	int nRefNum=FindPosi(PCent,ref);
	if (nRefNum<0) return false;
	double refSampling = pathLength(ref)/nSamNUM;
	int nRefNextNum=nRefNum+ddist/refSampling;
	//cl
		int nclNum=FindPosi(PCent,cl);
	if (nclNum<0) return false;
	double clSampling = pathLength(cl)/nSamNUM;
	int nclNextNum=nclNum+ddist/clSampling;
	//calculate the cosine angle
	float frefvector[3]={ref[nRefNextNum]._x-ref[nRefNum]._x,ref[nRefNextNum]._y-ref[nRefNum]._y,ref[nRefNextNum]._z-ref[nRefNum]._z};
		float fclvector[3]={cl[nclNextNum]._x-cl[nclNum]._x,cl[nclNextNum]._y-cl[nclNum]._y,cl[nclNextNum]._z-cl[nclNum]._z};
	float fcos=zxh::VectorOP_Cosine(frefvector,fclvector,3);
		return fcos;
}
bool CalcDirVec(PointCordTypeDef PPont,std::vector<jdq2017::point3D> &path,float fDirVec[3])
{
	int nRefNum=FindPosi(PPont,path);//Pcent is the current point in the unorganised point set
	PointCordTypeDef RefPcent,LRefPcent,NRefPcent;//RefPcent is the point in the ref point set which is nearest to the ref Pcent;
	                                    //LRefPcent is the last point of RefPcent
	RefPcent.x=path[nRefNum]._x;
	RefPcent.y=path[nRefNum]._y;
	RefPcent.z=path[nRefNum]._z;
	
	LRefPcent.x=path[nRefNum-1]._x;
	LRefPcent.y=path[nRefNum-1]._y;
	LRefPcent.z=path[nRefNum-1]._z;
	
	NRefPcent.x=path[nRefNum+1]._x;
	NRefPcent.y=path[nRefNum+1]._y;
	NRefPcent.z=path[nRefNum+1]._z;
	
	fDirVec[0]=NRefPcent.x-LRefPcent.x;
	fDirVec[1]=NRefPcent.y-LRefPcent.y;
	fDirVec[2]=NRefPcent.z-LRefPcent.z;
	return true;
}
bool DeterIfBranPont1(PointCordTypeDef PCent,std::vector<jdq2017::point3D> &ref,std::vector<jdq2017::point3D> &cl)
{
	double ddist=2;
	int nSamNUM=500;
	//ref
	int nRefNum=FindPosi(PCent,ref);//Pcent is the current point in the unorganised point set
	PointCordTypeDef RefPcent,LRefPcent;//RefPcent is the point in the ref point set which is nearest to the ref Pcent;
	                                    //LRefPcent is the last point of RefPcent
	//ref current
	RefPcent.x=ref[nRefNum]._x;
	RefPcent.y=ref[nRefNum]._y;
	RefPcent.z=ref[nRefNum]._z;
	float fRefDirVec[3]={0,0,0};
	CalcDirVec(RefPcent,ref,fRefDirVec);
	//ref last
	
	LRefPcent.x=ref[nRefNum-1]._x;
	LRefPcent.y=ref[nRefNum-1]._y;
	LRefPcent.z=ref[nRefNum-1]._z;
	float fLRefDirVec[3]={0,0,0};
	CalcDirVec(LRefPcent,ref,fLRefDirVec);

	//
	//cl
	int nClNum=FindPosi(PCent,cl);//Pcent is the current point in the unorganised point set
	PointCordTypeDef ClPcent,LClPcent;//RefPcent is the point in the ref point set which is nearest to the ref Pcent;
	                                    //LRefPcent is the last point of RefPcent
	//cl current
	ClPcent.x=cl[nClNum]._x;
	ClPcent.y=cl[nClNum]._y;
	ClPcent.z=cl[nClNum]._z;
	float fClDirVec[3]={0,0,0};
	CalcDirVec(ClPcent,cl,fClDirVec);
	//cl last
	LClPcent.x=cl[nClNum-1]._x;
	LClPcent.y=cl[nClNum-1]._y;
	LClPcent.z=cl[nClNum-1]._z;
	float fLClDirVec[3]={0,0,0};
	CalcDirVec(LRefPcent,cl,fLClDirVec);
	//calculate the angle between the two vectors
	float fcosCur=zxh::VectorOP_Cosine(fRefDirVec,fLRefDirVec,3);
	float fcosLas=zxh::VectorOP_Cosine(fClDirVec,fLClDirVec,3);

	return true;
}
bool AdEvePontsBySec(float clSampling,vector<PointCordTypeDef>&vUnorgaPointsWorld,std::vector<jdq2017::point3D> ref,std::vector<jdq2017::point3D> cl)
{
//set every point in vPathPointsWorld as the origin
	float fro=1;	
	float H=6*clSampling;
	int i=0;
	//intialiation with a start point and start direction


	bool bRoM1=false;
	int nINUM=0;
	while (i<vUnorgaPointsWorld.size())
	{

		if(i==503)
		{
			int x=0;
		}
		//collet points from vUnorgaPointsWorld within a ball of radius H
		vector<PointCordTypeDef> vLocalPointsInBallWorld,vLocalPointsInBallWorld_ori;
		PointCordTypeDef PCent=vUnorgaPointsWorld[i];
		if(bRoM1)
		{
			Collet2(PCent,H,vUnorgaPointsWorld,fro,nINUM,ref,cl,vLocalPointsInBallWorld);
		}
		else
		{
			Collet1(PCent,H,vUnorgaPointsWorld,vLocalPointsInBallWorld);
		}
		vLocalPointsInBallWorld_ori.assign(vLocalPointsInBallWorld.begin(),vLocalPointsInBallWorld.end());
		StorUnique(vLocalPointsInBallWorld);	
		//fit a plane
		float fabc[4]={0,0,0,0};
		//			//z=a*x+b*y+c;
		PlaneFitting2(PCent,H,vLocalPointsInBallWorld,fabc);
		//PlaneFitting(vLocalPointsInBallWorld,fabc);
		//project the points to the plane; points are 3D
		vector<PointCordTypeDef> Pro_vLocalPointsInBallWorld;
		Matrix<float,1,3> MPro_Pcent;//3D
		ProjPontsToPlaen(vLocalPointsInBallWorld,fabc,PCent,MPro_Pcent,Pro_vLocalPointsInBallWorld);
		vector<PointCordTypeDef> Pro_vLocalPointsInBallWorld_ori;
		Pro_vLocalPointsInBallWorld_ori.assign(Pro_vLocalPointsInBallWorld.begin(),Pro_vLocalPointsInBallWorld.end());
		//transform the points to local axis as 2D points
		vector<PointCordTypeDef>Pro_RotedvLocalPointsInBallWorld;
		Matrix<float,3,3> Mrx1,Mrx2;
		if(!Pro_RotedvLocalPointsInBallWorld.empty())Pro_RotedvLocalPointsInBallWorld.clear();
		Matrix<float,1,2> MPro_Pcent_2D;//2D on  the plane
		RotaPontsToLoca1(Pro_vLocalPointsInBallWorld,fabc,MPro_Pcent,MPro_Pcent_2D,Pro_RotedvLocalPointsInBallWorld,Mrx1,Mrx2);	
		//Pro_RotedvLocalPointsInBallWorld are the 2D points set
		float fab[3]={0,0,0};
		//linear regression
		//if det|A|=0,false;
		//bool bxyisnot0=LocalRegLine1(MPro_Pcent_2D,H,Pro_RotedvLocalPointsInBallWorld);
		//LocalRegLine_Gener(MfirPont_Local,H,Pro_RotedvLocalPointsInBallWorld,fab);
		//LocalRegLine1(MPro_Pcent_2D,H,Pro_RotedvLocalPointsInBallWorld);
		LocalRegLine3(MPro_Pcent_2D,H,Pro_RotedvLocalPointsInBallWorld,fab);
		vector<PointCordTypeDef> Pro_ForRovRotaLocalPointsInBallWorld;
		//rotate the local regression line with pi/4-theta

		float fxx=1.5;
		float fs=fxx/50;
		vector<PointCordTypeDef> vPonfxy;
		for (int kx=-50;kx<=50;kx++)
		{
			float fx=kx*fs;
			float fy=-(fab[0]+fab[1]*fx)/(fab[2]+0.00001);
			PointCordTypeDef ptem;
			ptem.x=fx;
			ptem.y=fy;
			vPonfxy.push_back(ptem);
		}

		RotLocaRegreLineForRo(i,fab,Pro_RotedvLocalPointsInBallWorld,Pro_ForRovRotaLocalPointsInBallWorld);
		//calculate the correlation coefficient
		fro=CalcRp(Pro_ForRovRotaLocalPointsInBallWorld);

		//2D curve fitting
		//y=ax2+bx+c
		if (fabs(fro)>0.7)
		{
			float fabc[3]={0,0,0};
			Matrix<float,2,2> Mrot;
			vector<PointCordTypeDef> Pro_ForQuvRotaLocalPointsInBallWorld;
			RotLocaRegreLineForQu1(i,fab,Pro_RotedvLocalPointsInBallWorld,Mrot,Pro_ForQuvRotaLocalPointsInBallWorld);
			//calculate the displcement between the new first point (new_x,new_y, new_z)and the originial point (0,0,0)
			float ffirstpointdist=CurFit2D(MPro_Pcent_2D,H,Pro_ForQuvRotaLocalPointsInBallWorld,fabc);
			//transform the local origin point to the global axis
			PointCordTypeDef PNewfirPont;//3D
			BackRptaPontsToGlo(ffirstpointdist,Mrot,Mrx1,Mrx2,MPro_Pcent,PNewfirPont);
			//update the point sets
			vUnorgaPointsWorld[i]=PNewfirPont;
			i++;
			H=6*clSampling;
			nINUM=0;
			bRoM1=false;
		}
		else
		{	//Points_3D.txt contains the points before plane fitting
			//char *Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D.txt";
			//WriteCA2Txt(vLocalPointsInBallWorld,Points_3D_Filename);
			////Pro_Points_3D.txt contains the projected points after plane fitting
			/*char *Pro_Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Pro_Points_3D.txt";
			WriteCA2Txt(Pro_RotedvLocalPointsInBallWorld,Pro_Points_3D_Filename);*/
			//Pro_Local_Points_2D.txt contains the projected points after plane fitting
			//char *Pro_Local_Points_2D="F:/Coronary_0/code/Resutsfusion/Pro_Local_Points_2D.txt";
			//WriteCA2Txt(Pro_RotedvLocalPointsInBallWorld,Pro_Local_Points_2D);
			//Points_2D.txt contains the points for calculate fro
			//char *Points_2D_Filename="F:/Coronary_0/code/Resutsfusion/Points_2D.txt";
			//WriteCA2Txt(Pro_ForRovRotaLocalPointsInBallWorld,Points_2D_Filename);
			//Determine if it is a bifurcation

			//float fcos=DeterIfBranPont1(PCent,ref,cl);
			//if(fcos<0.6)//it is a bifurcation point
			bRoM1=true;
			nINUM++;
			char *Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D.txt";
			WriteCA2Txt(vUnorgaPointsWorld,Points_3D_Filename);
			if(i==503)
			{
			char *Points_2D_Filename="F:/Coronary_0/code/Resutsfusion/Points_2D.txt";
			WriteCA2Txt(Pro_ForRovRotaLocalPointsInBallWorld,Points_2D_Filename);
			}

		}
	}
	return true;
}
bool AdEvePontsBySec_ALL(float clSampling,vector<PointCordTypeDef>&vUnorgaPointsWorld,vector<vLinesDef> &vnLine)
{
//set every point in vPathPointsWorld as the origin
	float fro=1;	
	float H=3*clSampling;
	int i=0;
	//intialiation with a start point and start direction


	bool bRoM1=false;
	int nINUM=0;
	while (i<vUnorgaPointsWorld.size())
	{

		//collet points from vUnorgaPointsWorld within a ball of radius H
		if(i==972||i==988||i==989)
		{
			int x=0;
		}
		vector<PointCordTypeDef> vLocalPointsInBallWorld,vLocalPointsInBallWorld_ori;
		PointCordTypeDef PCent=vUnorgaPointsWorld[i];
		if(bRoM1)
		{
			Collet2_ALL(i,clSampling,PCent,H,vUnorgaPointsWorld,fro,nINUM,vnLine,vLocalPointsInBallWorld);
		}
		else
		{
			Collet1(PCent,H,vUnorgaPointsWorld,vLocalPointsInBallWorld);
		}
		vLocalPointsInBallWorld_ori.assign(vLocalPointsInBallWorld.begin(),vLocalPointsInBallWorld.end());
		StorUnique(vLocalPointsInBallWorld);	
		//fit a plane
		float fabc[4]={0,0,0,0};
		//			//z=a*x+b*y+c;
		PlaneFitting3(PCent,H,vLocalPointsInBallWorld,fabc);
		//PlaneFitting(vLocalPointsInBallWorld,fabc);
		//project the points to the plane; points are 3D
		vector<PointCordTypeDef> Pro_vLocalPointsInBallWorld;
		Matrix<float,1,3> MPro_Pcent;//3D
		ProjPontsToPlaen(vLocalPointsInBallWorld,fabc,PCent,MPro_Pcent,Pro_vLocalPointsInBallWorld);
		vector<PointCordTypeDef> Pro_vLocalPointsInBallWorld_ori;
		Pro_vLocalPointsInBallWorld_ori.assign(Pro_vLocalPointsInBallWorld.begin(),Pro_vLocalPointsInBallWorld.end());
		//transform the points to local axis as 2D points
		vector<PointCordTypeDef>Pro_RotedvLocalPointsInBallWorld;
		Matrix<float,3,3> Mrx1,Mrx2;
		if(!Pro_RotedvLocalPointsInBallWorld.empty())Pro_RotedvLocalPointsInBallWorld.clear();
		Matrix<float,1,2> MPro_Pcent_2D;//2D on  the plane
		RotaPontsToLoca1(Pro_vLocalPointsInBallWorld,fabc,MPro_Pcent,MPro_Pcent_2D,Pro_RotedvLocalPointsInBallWorld,Mrx1,Mrx2);	
		//Pro_RotedvLocalPointsInBallWorld are the 2D points set
		float fab[3]={0,0,0};
		//linear regression
		//if det|A|=0,false;
		//bool bxyisnot0=LocalRegLine1(MPro_Pcent_2D,H,Pro_RotedvLocalPointsInBallWorld);
		//LocalRegLine_Gener(MfirPont_Local,H,Pro_RotedvLocalPointsInBallWorld,fab);
		//LocalRegLine1(MPro_Pcent_2D,H,Pro_RotedvLocalPointsInBallWorld);
		LocalRegLine5(MPro_Pcent_2D,H,Pro_RotedvLocalPointsInBallWorld,fab);
		vector<PointCordTypeDef> Pro_ForRovRotaLocalPointsInBallWorld;
		//rotate the local regression line with pi/4-theta

		float fxx=1.5;
		float fs=fxx/50;
		vector<PointCordTypeDef> vPonfxy;
		for (int kx=-50;kx<=50;kx++)
		{
			float fx=kx*fs;
			float fy=-(fab[0]+fab[1]*fx)/(fab[2]+0.00001);
			PointCordTypeDef ptem;
			ptem.x=fx;
			ptem.y=fy;
			vPonfxy.push_back(ptem);
		}

		RotLocaRegreLineForRo(i,fab,Pro_RotedvLocalPointsInBallWorld,Pro_ForRovRotaLocalPointsInBallWorld);
		//calculate the correlation coefficient
		fro=CalcRp(Pro_ForRovRotaLocalPointsInBallWorld);

		//2D curve fitting
		//y=ax2+bx+c
		if (fabs(fro)>0.7)
		{
			float fabc[3]={0,0,0};
			Matrix<float,2,2> Mrot;
			vector<PointCordTypeDef> Pro_ForQuvRotaLocalPointsInBallWorld;
			RotLocaRegreLineForQu1(i,fab,Pro_RotedvLocalPointsInBallWorld,Mrot,Pro_ForQuvRotaLocalPointsInBallWorld);
			//calculate the displcement between the new first point (new_x,new_y, new_z)and the originial point (0,0,0)
			float ffirstpointdist=CurFit2D(MPro_Pcent_2D,H,Pro_ForQuvRotaLocalPointsInBallWorld,fabc);
			//transform the local origin point to the global axis
			PointCordTypeDef PNewfirPont;//3D
			BackRptaPontsToGlo(ffirstpointdist,Mrot,Mrx1,Mrx2,MPro_Pcent,PNewfirPont);
			//update the point sets
			vUnorgaPointsWorld[i]=PNewfirPont;
			i++;
			H=6*clSampling;
			nINUM=0;
			bRoM1=false;
		}
		else
		{	//Points_3D.txt contains the points before plane fitting
			//char *Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D.txt";
			//WriteCA2Txt(vLocalPointsInBallWorld,Points_3D_Filename);
			////Pro_Points_3D.txt contains the projected points after plane fitting
			/*char *Pro_Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Pro_Points_3D.txt";
			WriteCA2Txt(Pro_RotedvLocalPointsInBallWorld,Pro_Points_3D_Filename);*/
			//Pro_Local_Points_2D.txt contains the projected points after plane fitting
			//char *Pro_Local_Points_2D="F:/Coronary_0/code/Resutsfusion/Pro_Local_Points_2D.txt";
			//WriteCA2Txt(Pro_RotedvLocalPointsInBallWorld,Pro_Local_Points_2D);
			//Points_2D.txt contains the points for calculate fro
			//char *Points_2D_Filename="F:/Coronary_0/code/Resutsfusion/Points_2D.txt";
			//WriteCA2Txt(Pro_ForRovRotaLocalPointsInBallWorld,Points_2D_Filename);
			//Determine if it is a bifurcation

			//float fcos=DeterIfBranPont1(PCent,ref,cl);
			//if(fcos<0.6)//it is a bifurcation point
			bRoM1=true;
			nINUM++;
		
		    if(i==972||i==988||i==989)
			{
				char *Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D_Local.txt";
			WriteCA2Txt(vLocalPointsInBallWorld,Points_3D_Filename);
			char *Points_2D_Filename="F:/Coronary_0/code/Resutsfusion/Points_2D.txt";
			WriteCA2Txt(Pro_ForRovRotaLocalPointsInBallWorld,Points_2D_Filename);
			char *Bef_Points_2D_Filename="F:/Coronary_0/code/Resutsfusion/Bef_Points_2D.txt";
			WriteCA2Txt(Pro_RotedvLocalPointsInBallWorld,Bef_Points_2D_Filename);
			char *Reg_Points_2D_Filename="F:/Coronary_0/code/Resutsfusion/Reg_Points_2D.txt";
			WriteCA2Txt(vPonfxy,Reg_Points_2D_Filename);
			int x=0;
			}

		}
	}
	return true;
}
bool ReadFilesInAfolder(char *chCurvesfolder,vector<vLinesDef> &svnLine)
{
    FILE *fp = NULL;
    FILE *fq = NULL;
    char filename_body[10];
	if(!svnLine.empty())
		svnLine.clear();
    for(int i = 0; i <= 30; i++)
    {
		char filename[128];
		strcpy(filename,chCurvesfolder);
		strcat(filename,"CL");

		memset(filename_body, 0, sizeof(filename_body));
		_itoa(i, filename_body, 10);
		strcat(filename,filename_body);

		char filename_tail[5] = ".vtk";
		strcat(filename,filename_tail);

		if((fp=fopen(filename, "r")) == NULL)
		{
			cout<<"open filename fail..."<<endl;
			return false;
		}
		std::vector<jdq2017::point3D> line;
		if (! jdq2017::readCenterlinevtk(filename, line) )//read the reference line
		{
			return false;
		}	
		vLinesDef tmp;
		tmp.vLine=line;
		tmp.nLine=i;
		svnLine.push_back(tmp);
		fclose(fp);
    }
	return true;
}
bool ResamEveryCurve(vector<vLinesDef> &vnLine,double clSampling)
{
	vLinesDef vLDTemp;
	std::vector<jdq2017::point3D> vLtemp;
	for (int i=0;i<vnLine.size();i++)
	{
	
		vLtemp=vnLine[i].vLine;	
		if(i==1)
		{
			char *Points_3D_ori_Filename="F:/Coronary_0/code/Resutsfusion/v1_ori.txt";
	//WriteCA2Txt_Skip(vUnorgaPointsWorld_ori,Points_3D_ori_Filename);
	WriteCA2Txt_jdq(vLtemp,Points_3D_ori_Filename);
		}
		jdq2017::ResamplePaths<jdq2017::point3D> resamplerr;
		resamplerr(vLtemp, clSampling);
		vLtemp = resamplerr.resultPath();
			if(i==1)
		{
			char *Points_3D_ori_Filename="F:/Coronary_0/code/Resutsfusion/v1_res.txt";
	//WriteCA2Txt_Skip(vUnorgaPointsWorld_ori,Points_3D_ori_Filename);
	WriteCA2Txt_jdq(vLtemp,Points_3D_ori_Filename);
		}
		vnLine[i].vLine=vLtemp;
	}
	return true;
}

bool MapModelPointsToImage(zxhImageDataT<short>&imgReadNewRaw,vector<PointCordTypeDef> vPathPointsWorld)//map line into original image in a range
{
	int ImgNewSize[4]={1};
	imgReadNewRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);

		for(int it=0;it<ImgNewSize[3];++it)
			for(int iz=0;iz<ImgNewSize[2];++iz)
				for(int iy=0;iy<ImgNewSize[1];++iy)
					for(int ix=0;ix<ImgNewSize[0];++ix)
					{
						imgReadNewRaw.SetPixelByGreyscale(ix,iy,iz,it,0);
					}
	vector<PointCordTypeDef> vPathPointsWorldMAPT;
	vPathPointsWorldMAPT.clear();
	for (int i=0;i<vPathPointsWorld.size();i++)//map the vetor points to the image
	{
		float PointPosWorld[ZXH_ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointCordTypeDef PointMAPT;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		imgReadNewRaw.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=zxh::round(PointPosWorld[1]);
		PointPos[2]=zxh::round(PointPosWorld[2]);
		BoundaryCorrect(PointPos,ImgNewSize);
		PointMAPT.x=PointPos[0];
		PointMAPT.y=PointPos[1];
		PointMAPT.z=PointPos[2];
		vPathPointsWorldMAPT.push_back(PointMAPT);
		if(PointPos[0]==212&&PointPos[1]==184&&PointPos[2]==60)
		{
			int xxx=0;
		}
		if(PointPos[0]==212&&PointPos[1]==185&&PointPos[2]==60)
		{
			int xxx=0;
		}
		if(PointPos[0]==213&&PointPos[1]==185&&PointPos[2]==60)
		{
			int xxx=0;
		}
		imgReadNewRaw.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);

	}
	return true;
}
bool MapImgPointsToImage(zxhImageDataT<short>&imgReadNewRaw,vector<PointImgTypeDef> vPathPointsWorld)//map line into original image in a range
{
	int ImgNewSize[4]={1};
	imgReadNewRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);

		for(int it=0;it<ImgNewSize[3];++it)
			for(int iz=0;iz<ImgNewSize[2];++iz)
				for(int iy=0;iy<ImgNewSize[1];++iy)
					for(int ix=0;ix<ImgNewSize[0];++ix)
					{
						imgReadNewRaw.SetPixelByGreyscale(ix,iy,iz,it,0);
					}
	for (int i=0;i<vPathPointsWorld.size();i++)
	{
		PointImgTypeDef tmp;
		tmp=vPathPointsWorld[i];
		imgReadNewRaw.SetPixelByGreyscale(tmp.x,tmp.y,tmp.z,0,ZXH_Foreground);

	}
	return true;
}
bool Points_Init(zxhImageDataT<short>&imgReadNewRaw,zxhImageDataT<short>&imgCountmap)//map line into original image in a range
{
	int gNbr[26][3] = { 
		{-1, 0, 0}, \
		{-1, -1, 0}, \
		{-1, 1, 0}, \
		{-1, 0, -1}, \
		{-1, -1, -1}, \
		{-1, -1, 1}, \
		{-1, 1, -1}, \
		{-1, 1, 1}, \
		{-1, 0, 1}, \
		{1, 0, 0}, \
		{1, -1, 0}, \
		{1, 1, 0}, \
		{1, 0, -1}, \
		{1, -1, -1}, \
		{1, -1, 1}, \
		{1, 1, -1}, \
		{1, 1, 1}, \
		{1, 0, 1}, \
		{ 0,-1, 0}, \
		{ 0, 1, 0}, \
		{ 0, 0,-1}, \
		{ 0, 0, 1},\
		{ 0, -1,-1}, \
		{ 0, -1,1}, \
		{ 0, 1,-1}, \
		{ 0, 1,1}, \
	};
	/*
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };
	*/
	//store the seed points.
	vector<PointCordTypeDef> vSeedsPonts;
	int ImgNewSize[4]={1};
	imgReadNewRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short sint=imgReadNewRaw.GetPixelGreyscale(ix,iy,iz,it);
					if(sint==ZXH_Foreground)
					{
						PointCordTypeDef pseed;
						pseed.x=ix;
						pseed.y=iy;
						pseed.z=iz;
						vSeedsPonts.push_back(pseed);
					}

				}
				//cout the neibour points

				for (int i=0;i<vSeedsPonts.size();i++)//map the vetor points to the image
				{
					PointCordTypeDef pseed;
					pseed=vSeedsPonts[i];
					for (int i = 0; i < 26; i++)
					{
						int nx = pseed.x + gNbr[i][0];
						int ny = pseed.y + gNbr[i][1];
						int nz = pseed.z + gNbr[i][2];
						short sint=imgReadNewRaw.GetPixelGreyscale(nx,ny,nz,0);
						if(nx>=0&&nx<ImgNewSize[0]&&ny>=0&&ny<ImgNewSize[1]&&nz>=0&&nz<ImgNewSize[2]&&sint==ZXH_Foreground)
						{
							int npx=(int)pseed.x;
							int npy=(int)pseed.y;
							int npz=(int)pseed.z;
							short sintc=imgCountmap.GetPixelGreyscale(npx,npy,npz,0);
							sintc=sintc+1;
							imgCountmap.SetPixelByGreyscale(npx,npy,npz,0,sintc);
						}
					}
				}
char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_Cont.nii.gz";
	string chFileName2(chResultName);
	zxh::SaveImage(&imgCountmap,chFileName2.c_str());
	return true;
}
bool Calc_NumOfNei(zxhImageDataT<short>&imgReadNewRaw,zxhImageDataT<short>&imgCountmap,vector<PointImgTypeDef>&vSeedsPonts)//map line into original image in a range
{
	int gNbr[26][3] = { 
		{-1, 0, 0}, \
		{-1, -1, 0}, \
		{-1, 1, 0}, \
		{-1, 0, -1}, \
		{-1, -1, -1}, \
		{-1, -1, 1}, \
		{-1, 1, -1}, \
		{-1, 1, 1}, \
		{-1, 0, 1}, \
		{1, 0, 0}, \
		{1, -1, 0}, \
		{1, 1, 0}, \
		{1, 0, -1}, \
		{1, -1, -1}, \
		{1, -1, 1}, \
		{1, 1, -1}, \
		{1, 1, 1}, \
		{1, 0, 1}, \
		{ 0,-1, 0}, \
		{ 0, 1, 0}, \
		{ 0, 0,-1}, \
		{ 0, 0, 1},\
		{ 0, -1,-1}, \
		{ 0, -1,1}, \
		{ 0, 1,-1}, \
		{ 0, 1,1}, \
	};
	/*
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };
	*/
	//store the seed points.
				//cout the neibour points
	int ImgNewSize[4]={1};
	imgReadNewRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
				for (int i=0;i<vSeedsPonts.size();i++)//map the vetor points to the image
				{
					PointImgTypeDef pseed;
					pseed=vSeedsPonts[i];
					for (int i = 0; i < 26; i++)
					{
						int nx = pseed.x + gNbr[i][0];
						int ny = pseed.y + gNbr[i][1];
						int nz = pseed.z + gNbr[i][2];
						short sint=imgReadNewRaw.GetPixelGreyscale(nx,ny,nz,0);
						if(nx>=0&&nx<ImgNewSize[0]&&ny>=0&&ny<ImgNewSize[1]&&nz>=0&&nz<ImgNewSize[2]&&sint==ZXH_Foreground)
						{
							int npx=(int)pseed.x;
							int npy=(int)pseed.y;
							int npz=(int)pseed.z;
							short sintc=imgCountmap.GetPixelGreyscale(npx,npy,npz,0);
							sintc=sintc+1;
							imgCountmap.SetPixelByGreyscale(npx,npy,npz,0,sintc);
						}
					}
				}
			
char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_Cont.nii.gz";
	string chFileName2(chResultName);	
	fstream _file;
     _file.open(chFileName2,ios::in);
	  if(_file)
     {
		remove(chResultName);
	  }
	zxh::SaveImage(&imgCountmap,chFileName2.c_str());
	return true;
}
bool Points_Init1(zxhImageDataT<short>&imgReadNewRaw,zxhImageDataT<short>&imgCountmap,vector<PointImgTypeDef>&vSeedsPonts)//map line into original image in a range
{
	int gNbr[26][3] = { 
		{-1, 0, 0}, \
		{-1, -1, 0}, \
		{-1, 1, 0}, \
		{-1, 0, -1}, \
		{-1, -1, -1}, \
		{-1, -1, 1}, \
		{-1, 1, -1}, \
		{-1, 1, 1}, \
		{-1, 0, 1}, \
		{1, 0, 0}, \
		{1, -1, 0}, \
		{1, 1, 0}, \
		{1, 0, -1}, \
		{1, -1, -1}, \
		{1, -1, 1}, \
		{1, 1, -1}, \
		{1, 1, 1}, \
		{1, 0, 1}, \
		{ 0,-1, 0}, \
		{ 0, 1, 0}, \
		{ 0, 0,-1}, \
		{ 0, 0, 1},\
		{ 0, -1,-1}, \
		{ 0, -1,1}, \
		{ 0, 1,-1}, \
		{ 0, 1,1}, \
	};
	/*
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };
	*/
	//store the seed points.
	
	int ImgNewSize[4]={1};
	imgReadNewRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short sint=imgReadNewRaw.GetPixelGreyscale(ix,iy,iz,it);
					if(sint==ZXH_Foreground)
					{
						PointImgTypeDef pseed;
						pseed.x=ix;
						pseed.y=iy;
						pseed.z=iz;
						pseed.val=0;
						vSeedsPonts.push_back(pseed);
					}

				}
	return true;
}
bool Coutmap_vec(zxhImageDataT<short>&imgCountmap,std::vector<std::vector<std::vector<int> > > &BW,std::vector<std::vector<std::vector<int> > > &BW1)
{

		int ImgSize[4]={1};
	imgCountmap.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
		//init
	for(int z=0;z<ImgSize[2];z++)  
    {  
        for (int x=0;x<ImgSize[0];x++)  
        {  
            for (int y=0;y<ImgSize[1];y++)  
            {  
                BW[z][x][y]=0;  
                
            }  
        }  
    }  
	BW1.assign(BW.begin(),BW.end());
	for(int it=0;it<ImgSize[3];++it)
		for(int iz=0;iz<ImgSize[2];++iz)
			for(int ix=0;ix<ImgSize[0];++ix)
				for(int iy=0;iy<ImgSize[1];++iy)
				{
					short scon=imgCountmap.GetPixelGreyscale(ix,iy,iz,it);
					if(scon>=5||scon==1)
					{
						BW[iz][ix][iy]=1;  
					}
				}
	return true;
}
bool Coutmap_vec1(zxhImageDataT<short>&imgCountmap,zxhImageDataT<short>&imgCountConNeimap,std::vector<std::vector<std::vector<int> > > &BW,std::vector<std::vector<std::vector<int> > > &BW1)
{

		int ImgSize[4]={1};
	imgCountmap.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
		//init
	for(int z=0;z<ImgSize[2];z++)  
    {  
        for (int x=0;x<ImgSize[0];x++)  
        {  
            for (int y=0;y<ImgSize[1];y++)  
            {  
                BW[z][x][y]=0;  
                
            }  
        }  
    }  
	BW1.assign(BW.begin(),BW.end());
	zxhImageDataT<short> imgCountmapAndNei;
	imgCountmapAndNei.NewImage(imgCountmap.GetImageInfo());
	for(int it=0;it<ImgSize[3];++it)
		for(int iz=0;iz<ImgSize[2];++iz)
			for(int ix=0;ix<ImgSize[0];++ix)
				for(int iy=0;iy<ImgSize[1];++iy)
				{
					short scon=imgCountmap.GetPixelGreyscale(ix,iy,iz,it);
					short sconnei=imgCountConNeimap.GetPixelGreyscale(ix,iy,iz,it);
					if(scon>=5||scon==1||sconnei==1)
					{
						BW[iz][ix][iy]=1;  
						imgCountmapAndNei.SetPixelByGreyscale(ix,iy,iz,it,1);
					}
				}
				char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_CountAndNei.nii.gz";
				string chFileName2(chResultName);
				zxh::SaveImage(&imgCountmapAndNei,chFileName2.c_str());
	return true;
}
bool FILLRUNVECTORS(std::vector<std::vector<std::vector<int> > > A, int& NumberOfRuns, vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun)
{
	for(int k = 0; k < A.size(); k++)
	{
		std::vector<std::vector<int> >  Ahight;
		Ahight=A[k];

		for (int i = 0; i < Ahight.size(); i++)
		{
			std::vector<int> VrowData=Ahight[i];
			int ik=k*Ahight.size()+i;
			if (VrowData[0] == 1)
			{
				NumberOfRuns++;
				stRun.push_back(0);
				rowRun.push_back(ik);
			}
			for (int j = 1; j < VrowData.size(); j++)
			{
				if (VrowData[j - 1] == 0 && VrowData[j] == 1)
				{
					NumberOfRuns++;
					stRun.push_back(j);
					rowRun.push_back(ik);
				}
				else if (VrowData[j - 1] == 1 && VrowData[j] == 0)
				{
					enRun.push_back(j - 1);
				}
			}
			if (VrowData[VrowData.size() - 1] == 1)
			{
				enRun.push_back(VrowData.size() - 1);
			}
		}
	}

	return true;
}
bool Calc3Hrun(int i,vector<int>& rowRun,int HcurRowIdx,int curRowIdx,int &HfirstRunOnPre,int &HlastRunOnPre)
{
	int Hnumrunbe=0;
	int Hnumrunone=0;

	for(int m=0;m<rowRun.size();m++)
	{
		if(rowRun[m]>=HcurRowIdx&&rowRun[m]<=curRowIdx-1)
		{
			Hnumrunbe++;
		}
		if(rowRun[m]==HcurRowIdx)
		{
			Hnumrunone++;
		}
	}
	if(Hnumrunone==0)
	{
		HfirstRunOnPre= 0;
		HlastRunOnPre = -1;
		return true;
	}
	else
	{

		HlastRunOnPre= i - Hnumrunbe+Hnumrunone-1;
		HfirstRunOnPre = i - Hnumrunbe;
		return true;
	}
	return true;
}
bool CalcLables(int i,int HfirstRunOnPre,int HlastRunOnPre,vector<int>& stRun, vector<int>& enRun,vector<int>& runLabels, vector<pair<int, int>>& equivalences, int offset)
{
	for (int k = HfirstRunOnPre; k <= HlastRunOnPre; k++)
		{
			if (stRun[i] <= enRun[k] + offset && enRun[i] >= stRun[k] - offset)
			{
				if (runLabels[i] == 0) // never marked
					runLabels[i] = runLabels[k];
				else if (runLabels[i] != runLabels[k])// have marked             
					equivalences.push_back(make_pair(runLabels[i], runLabels[k])); // 保存等价对
			}
		}
	return true;
}
bool FIRSTPASS_1(int imgSize[4],vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun, int &NumberOfRuns,
    vector<int>& runLabels, vector<pair<int, int>>& equivalences, int offset)
{
	runLabels.assign(NumberOfRuns, 0);
	int idxLabel = 1;

	int curRowIdx = 0;
	int LcurRowIdx = 0;
	int firstRunOnCur = 0;
	int firstRunOnPre = 0;
	int lastRunOnPre = -1;

	int HcurRowIdx[3] = {0,0,0};
	int HfirstRunOnPre[4] = {0,0,0};
	int HlastRunOnPre[4] ={-1,-1,-1,-1};
	for (int i = 0; i < NumberOfRuns; i++)
	{
		//calculate the start and then end position of the runsl
		if (rowRun[i] != curRowIdx)
		{
			curRowIdx = rowRun[i];
			LcurRowIdx=curRowIdx-1;
			int numrunbe=0;
			int numrunone=0;

			for(int m=0;m<rowRun.size();m++)
			{

				if(rowRun[m]==LcurRowIdx)
				{
					numrunone++;
				}
			}
			if(curRowIdx%imgSize[0]==0||numrunone==0)
			{
				firstRunOnPre = 0;
				lastRunOnPre = -1;
			}
			else
			{
				firstRunOnPre = i - numrunone;
				lastRunOnPre = i -1;
			}
			HfirstRunOnPre[0] = firstRunOnPre;
			HlastRunOnPre[0] =lastRunOnPre;

			//max(rowRun[i]-3,0);
			//calculate the number of runs between current run and the front run.
			HcurRowIdx[0]=max(rowRun[i]-imgSize[0],-1);
			HcurRowIdx[1]=max(rowRun[i]-imgSize[0]-1,-1);
			HcurRowIdx[2]=max(rowRun[i]-imgSize[0]+1,-1);
			// the front row
			int CRow=HcurRowIdx[0];
			int CFRP=0;
			int CLRP=-1;
			Calc3Hrun(i,rowRun,CRow,curRowIdx,CFRP,CLRP);
			HfirstRunOnPre[1] = CFRP;
			HlastRunOnPre[1] =CLRP;
			//the up front row
			CRow=HcurRowIdx[1];
			if(CRow%imgSize[0]==imgSize[0]-1)
			{
				CFRP=0;
				CLRP=-1;
			}
			else
			{
				Calc3Hrun(i,rowRun,CRow,curRowIdx,CFRP,CLRP);
			}
			HfirstRunOnPre[2] = CFRP;
			HlastRunOnPre[2] =CLRP;
			//the down front row
			CRow=HcurRowIdx[2];
			if(CRow%imgSize[0]==0)
			{
				CFRP=0;
				CLRP=-1;
			}
			else
			{
				Calc3Hrun(i,rowRun,CRow,curRowIdx,CFRP,CLRP);
			}
			HfirstRunOnPre[3] = CFRP;
			HlastRunOnPre[3] =CLRP;
			

		}
		for (int h=0;h<4;h++)
		{
			int HCFRP=HfirstRunOnPre[h];
			int HCLRP=HlastRunOnPre[h];
			CalcLables(i,HCFRP,HCLRP,stRun,enRun,runLabels,equivalences, offset);
		}
		if (runLabels[i] == 0) // 没有与前一列的任何run重合
		{
			runLabels[i] = idxLabel++;
		}

	}
	return true;
}
	bool FIRSTPASS_2(int imgSize[4],vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun, int &NumberOfRuns,
    vector<int>& runLabels, vector<pair<int, int>>& equivalences, int offset)
{
	runLabels.assign(NumberOfRuns, 0);
	int idxLabel = 1;

	int curRowIdx = 0;
	int LcurRowIdx = 0;
	int firstRunOnCur = 0;
	int firstRunOnPre = 0;
	int lastRunOnPre = -1;

	int HcurRowIdx[3] = {0,0,0};
	int HfirstRunOnPre[4] = {0,0,0};
	int HlastRunOnPre[4] ={-1,-1,-1,-1};
	for (int i = 0; i < NumberOfRuns; i++)
	{
		//calculate the start and then end position of the runsl
		if (rowRun[i] != curRowIdx)
		{
			curRowIdx = rowRun[i];
			LcurRowIdx=curRowIdx-1;
			int numrunbe=0;
			int numrunone=0;

			for(int m=0;m<rowRun.size();m++)
			{

				if(rowRun[m]==LcurRowIdx)
				{
					numrunone++;
				}
			}
			if(curRowIdx%imgSize[0]==0||numrunone==0)
			{
				firstRunOnPre = 0;
				lastRunOnPre = -1;
			}
			else
			{
				firstRunOnPre = i - numrunone;
				lastRunOnPre = i -1;
			}
			HfirstRunOnPre[0] = firstRunOnPre;
			HlastRunOnPre[0] =lastRunOnPre;

			//max(rowRun[i]-3,0);
			//calculate the number of runs between current run and the front run.
			HcurRowIdx[0]=max(rowRun[i]-imgSize[0],-1);
			HcurRowIdx[1]=max(rowRun[i]-imgSize[0]-1,-1);
			HcurRowIdx[2]=max(rowRun[i]-imgSize[0]+1,-1);
			// the front row
			int CRow=HcurRowIdx[0];
			int CFRP=0;
			int CLRP=-1;
			Calc3Hrun(i,rowRun,CRow,curRowIdx,CFRP,CLRP);
			HfirstRunOnPre[1] = CFRP;
			HlastRunOnPre[1] =CLRP;
			//the up front row
			CRow=HcurRowIdx[1];
			if(CRow%imgSize[0]==imgSize[0]-1)
			{
				CFRP=0;
				CLRP=-1;
			}
			else
			{
				Calc3Hrun(i,rowRun,CRow,curRowIdx,CFRP,CLRP);
			}
			HfirstRunOnPre[2] = 0;
			HlastRunOnPre[2] =-1;
			//the down front row
			CRow=HcurRowIdx[2];
			if(CRow%imgSize[0]==0)
			{
				CFRP=0;
				CLRP=-1;
			}
			else
			{
				Calc3Hrun(i,rowRun,CRow,curRowIdx,CFRP,CLRP);
			}
			HfirstRunOnPre[3] = 0;
			HlastRunOnPre[3] =-1;
			

		}
		for (int h=0;h<4;h++)
		{
			int HCFRP=HfirstRunOnPre[h];
			int HCLRP=HlastRunOnPre[h];
			CalcLables(i,HCFRP,HCLRP,stRun,enRun,runLabels,equivalences, offset);
		}
		if (runLabels[i] == 0) // 没有与前一列的任何run重合
		{
			runLabels[i] = idxLabel++;
		}

	}
	return true;
}
bool FIRSTPASS(int imgSize[4],vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun, int &NumberOfRuns,
    vector<int>& runLabels, vector<pair<int, int>>& equivalences, int offset)
{
	runLabels.assign(NumberOfRuns, 0);

	int idxLabel = 1;

	int curRowIdx = 0;
	int LcurRowIdx = 0;
	int firstRunOnCur = 0;
	int firstRunOnPre = 0;
	int lastRunOnPre = -1;

	int HcurRowIdx = 0;
	int HfirstRunOnPre = 0;
	int HlastRunOnPre = -1;

	for (int i = 0; i < NumberOfRuns; i++)
	{

		//
		if (rowRun[i] != curRowIdx)
		{
			curRowIdx = rowRun[i];
			LcurRowIdx=curRowIdx-1;
			int numrunbe=0;
			int numrunone=0;

			for(int m=0;m<rowRun.size();m++)
			{
	
				if(rowRun[m]==LcurRowIdx)
				{
					numrunone++;
				}
			}
			if(curRowIdx%imgSize[0]==0||numrunone==0)
			{
				firstRunOnPre = 0;
				lastRunOnPre = -1;
			}
			else
			{
				firstRunOnPre = i - numrunone;
				lastRunOnPre = i -1;
			}

			//max(rowRun[i]-3,0);
			//calculate the number of runs between current run and the front run.
			int Hnumrunbe=0;
			int Hnumrunone=0;
			HcurRowIdx=max(rowRun[i]-imgSize[0],-1);
			for(int m=0;m<rowRun.size();m++)
			{
				if(rowRun[m]>=HcurRowIdx&&rowRun[m]<=curRowIdx-1)
				{
					Hnumrunbe++;
				}
				if(rowRun[m]==HcurRowIdx)
				{
					Hnumrunone++;
				}
			}
			if(Hnumrunone==0)
			{
				HfirstRunOnPre = 0;
				HlastRunOnPre = -1;
			}
			else
			{

				HlastRunOnPre = i - Hnumrunbe+Hnumrunone-1;
				HfirstRunOnPre = i - Hnumrunbe;
			}
		}
		for (int j = firstRunOnPre; j <= lastRunOnPre; j++)
		{
			if (stRun[i] <= enRun[j] + offset && enRun[i] >= stRun[j] - offset)
			{
				if (runLabels[i] == 0) // never marked
					runLabels[i] = runLabels[j];
				else if (runLabels[i] != runLabels[j])// have marked             
					equivalences.push_back(make_pair(runLabels[i], runLabels[j])); // 保存等价对
			}
		}
		for (int k = HfirstRunOnPre; k <= HlastRunOnPre; k++)
		{
			if (stRun[i] <= enRun[k] + offset && enRun[i] >= stRun[k] - offset)
			{
				if (runLabels[i] == 0) // never marked
					runLabels[i] = runLabels[k];
				else if (runLabels[i] != runLabels[k])// have marked             
					equivalences.push_back(make_pair(runLabels[i], runLabels[k])); // 保存等价对
			}
		}
		if (runLabels[i] == 0) // 没有与前一列的任何run重合
		{
			runLabels[i] = idxLabel++;
		}

	}
	return true;
}
bool replaceSameLabel(vector<int>& runLabels, vector<pair<int, int>>&
    equivalence)
{
    int maxLabel = *max_element(runLabels.begin(), runLabels.end());
    vector<vector<bool>> eqTab(maxLabel, vector<bool>(maxLabel, false));
    vector<pair<int, int>>::iterator vecPairIt = equivalence.begin();
    while (vecPairIt != equivalence.end())
    {
        eqTab[vecPairIt->first - 1][vecPairIt->second - 1] = true;
        eqTab[vecPairIt->second - 1][vecPairIt->first - 1] = true;
        vecPairIt++;
    }
    vector<int> labelFlag(maxLabel, 0);
    vector<vector<int>> equaList;
    vector<int> tempList;
    //cout << maxLabel << endl;
    for (int i = 1; i <= maxLabel; i++)
    {
        if (labelFlag[i - 1])
        {
            continue;
        }
        labelFlag[i - 1] = equaList.size() + 1;
        tempList.push_back(i);
        for (vector<int>::size_type j = 0; j < tempList.size(); j++)
        {
            for (vector<bool>::size_type k = 0; k != eqTab[tempList[j] - 1].size(); k++)
            {
                if (eqTab[tempList[j] - 1][k] && !labelFlag[k])
                {
                    tempList.push_back(k + 1);
                    labelFlag[k] = equaList.size() + 1;
                }
            }
        }
        equaList.push_back(tempList);
        tempList.clear();
    }
    //cout << equaList.size() << endl;
    for (vector<int>::size_type i = 0; i != runLabels.size(); i++)
    {
        runLabels[i] = labelFlag[runLabels[i] - 1];
    }
	return true;
}
bool BACKLABELED(int imgSize[4],std::vector<std::vector<std::vector<int> > > A,std::vector<std::vector<std::vector<int> > > &B,vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun, int &NumberOfRuns,
    vector<int>& runLabels,zxhImageDataT<short>&imgLabel)
{
	
	for(int i=0;i<NumberOfRuns;i++)
	{
		int nrow=rowRun[i];
		int stcol=stRun[i];
		int encol=enRun[i];
		int cols=imgSize[1];
		int iz=nrow/imgSize[0];
        int ix=nrow%imgSize[0];
		int zz=B.size();
		std::vector<std::vector<int> > Bx=B[iz];
		int xx=Bx.size();
		std::vector<int> By=Bx[ix];
		int yy=By.size();
		for(int j=0;j<cols;j++)
		{
			if(j>=stcol&&j<=encol)
			{
				B[iz][ix][j]=runLabels[i];
				imgLabel.SetPixelByGreyscale(ix,j,iz,0,runLabels[i]);
			}
			int k=j;
		}
	}
	char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_Lab.nii.gz";
	string chFileName2(chResultName);
	zxh::SaveImage(&imgLabel,chFileName2.c_str());

	return true;
}
bool BACKLABELEDW(int imgSize[4],std::vector<std::vector<std::vector<int> > > A,std::vector<std::vector<std::vector<int> > > &B,vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun, int &NumberOfRuns,
    vector<int>& runLabels)
{
	
	for(int i=0;i<NumberOfRuns;i++)
	{
		int nrow=rowRun[i];
		int stcol=stRun[i];
		int encol=enRun[i];
		int cols=imgSize[1];
		int iz=nrow/imgSize[0];
        int ix=nrow%imgSize[0];
		int zz=B.size();
		std::vector<std::vector<int> > Bx=B[iz];
		int xx=Bx.size();
		std::vector<int> By=Bx[ix];
		int yy=By.size();
		for(int j=0;j<cols;j++)
		{
			if(j>=stcol&&j<=encol)
			{
				B[iz][ix][j]=runLabels[i];
			}
		}
	}

	return true;
}
void BWLABEL3_test(int imgSize[3],std::vector<std::vector<std::vector<int> > > A,std::vector<std::vector<std::vector<int> > > &B,int &NumberOfLabs)
{

	vector<int> stRun;
	vector<int> enRun;
	vector<int> rowRun;
	int NumberOfRuns=0;
	vector<int> runLabels;
		//
	FILLRUNVECTORS(A,NumberOfRuns,stRun,enRun,rowRun);
	//

	vector<pair<int, int>> equivalences;
	int offset=1;
	FIRSTPASS_1(imgSize,stRun,enRun,rowRun,NumberOfRuns,runLabels, equivalences, offset);
	replaceSameLabel(runLabels,equivalences);
	std::vector<int>::iterator biggest = std::max_element(std::begin(runLabels), std::end(runLabels));  
	NumberOfLabs=*biggest;
}
void BWLABEL3(int imgSize[4],std::vector<std::vector<std::vector<int> > > A,std::vector<std::vector<std::vector<int> > > &B,int &NumberOfLabs,zxhImageDataT<short>&imgLabel)
{

	vector<int> stRun;
	vector<int> enRun;
	vector<int> rowRun;
	int NumberOfRuns=0;
	vector<int> runLabels;
		//
	FILLRUNVECTORS(A,NumberOfRuns,stRun,enRun,rowRun);
	//

	vector<pair<int, int>> equivalences;
	int offset=1;
	FIRSTPASS_1(imgSize,stRun,enRun,rowRun,NumberOfRuns,runLabels, equivalences, offset);
	replaceSameLabel(runLabels,equivalences);
	std::vector<int>::iterator biggest = std::max_element(std::begin(runLabels), std::end(runLabels));  
	NumberOfLabs=*biggest;
	BACKLABELED(imgSize,A,B,stRun, enRun, rowRun, NumberOfRuns,runLabels,imgLabel);
}
void SelCen(int ImgSize[4],std::vector<std::vector<std::vector<int> > > &BW1,std::vector<std::vector<std::vector<int> > > &BW2,int &NumberOfLabs,zxhImageDataT<short>&imgCenLabel)
{

	for (int i = 1; i <=NumberOfLabs; i++)
	{
		vector<int>vcandi;
		for(int z=0;z<ImgSize[2];z++)  
		{  
			for (int x=0;x<ImgSize[0];x++)  
			{  
				for (int y=0;y<ImgSize[1];y++)  
				{  
					int nlab=BW1[z][x][y];  
					if (nlab==i)
					{
						//iz* nImWY * nImWX + iy* nImWX + ix
						int ntemlab=z*ImgSize[0]*ImgSize[1]+x*ImgSize[1]+y;
						vcandi.push_back(ntemlab);
					}
				}  
			}  
		} 
		
		sort(vcandi.begin(),vcandi.end());
		int ncnumb=vcandi.size()/2;
		int iz=vcandi[ncnumb]/(ImgSize[0]*ImgSize[1]);
		int ix=(vcandi[ncnumb]%(ImgSize[0]*ImgSize[1]))/ImgSize[1];
		int iy=(vcandi[ncnumb]%(ImgSize[0]*ImgSize[1]))%ImgSize[1];
		BW2[iz][ix][iy]=i;
		short sm=i;
		if(ix<0||ix>ImgSize[0]-1||iy<0||iy>ImgSize[1]-1||iz<0||iz>ImgSize[2]-1)
		{
			int xxxx=0;
		}
		imgCenLabel.SetPixelByGreyscale(ix,iy,iz,0,sm);
		
	}
	char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_CentLab.nii.gz";
	string chFileName2(chResultName);
	zxh::SaveImage(&imgCenLabel,chFileName2.c_str());
}

bool fillRunVectors(Matrix<int, 4, 14> A, int& NumberOfRuns, vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun)
{
    for (int i = 0; i < A.rows(); i++)
    {
		Matrix<int, 1, 14> MrowData=A.row(i);

        if (MrowData[0] == 1)
        {
            NumberOfRuns++;
            stRun.push_back(0);
            rowRun.push_back(i);
        }
        for (int j = 1; j < A.cols(); j++)
        {
            if (MrowData[j - 1] == 0 && MrowData[j] == 1)
            {
                NumberOfRuns++;
                stRun.push_back(j);
                rowRun.push_back(i);
            }
            else if (MrowData[j - 1] == 1 && MrowData[j] == 0)
            {
                enRun.push_back(j - 1);
            }
        }
       if (MrowData[A.cols() - 1] == 1)
        {
            enRun.push_back(A.cols() - 1);
        }
    }
	return true;
}
bool firstPass(vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun, int NumberOfRuns,
    vector<int>& runLabels, vector<pair<int, int>>& equivalences, int offset)
{
    runLabels.assign(NumberOfRuns, 0);
    int idxLabel = 1;
    int curRowIdx = 0;
    int firstRunOnCur = 0;
    int firstRunOnPre = 0;
    int lastRunOnPre = -1;
    for (int i = 0; i < NumberOfRuns; i++)
    {
        if (rowRun[i] != curRowIdx)
        {
            curRowIdx = rowRun[i];
            firstRunOnPre = firstRunOnCur;
            lastRunOnPre = i - 1;
            firstRunOnCur = i;

        }
        for (int j = firstRunOnPre; j <= lastRunOnPre; j++)
        {
            if (stRun[i] <= enRun[j] + offset && enRun[i] >= stRun[j] - offset)
            {
                if (runLabels[i] == 0) // never marked
                    runLabels[i] = runLabels[j];
                else if (runLabels[i] != runLabels[j])// have marked             
                    equivalences.push_back(make_pair(runLabels[i], runLabels[j])); // 保存等价对
            }
        }
        if (runLabels[i] == 0) // 没有与前一列的任何run重合
                                          {
            runLabels[i] = idxLabel++;
        }

    }
	return true;
}

bool Backlabeled(Matrix<int, 4, 14> A,Matrix<int, 4, 14> &B,vector<int>& stRun, vector<int>& enRun, vector<int>& rowRun, int NumberOfRuns,
    vector<int>& runLabels)
{
	int rows=A.rows();
	int cols=A.cols();
	B.setZero(rows,cols);
	for(int i=0;i<NumberOfRuns;i++)
	{
		int nrow=rowRun[i];
		int stcol=stRun[i];
		int encol=enRun[i];
		for(int j=0;j<cols;j++)
		{
			if(j>=stcol&&j<=encol)
			{
				B(nrow,j)=runLabels[i];
			}
			int k=j;
		}
	}
	return true;
}

void bwlabel(Matrix<int, 4, 14> A)
{
	int NumberOfRuns=0;
	vector<int> stRun;
	vector<int> enRun;
		vector<int> rowRun;
	fillRunVectors(A,NumberOfRuns,stRun,enRun,rowRun);
	//
	vector<int> runLabels;
	vector<pair<int, int>> equivalences;
	int offset=1;
	firstPass(stRun,enRun, rowRun,NumberOfRuns,runLabels, equivalences, offset);
	//
	replaceSameLabel(runLabels,equivalences);
	//
	Matrix<int, 4, 14> B;
	Backlabeled(A,B,stRun, enRun, rowRun, NumberOfRuns,runLabels);
	//cout<<B<<endl;
	int x=0;
}

bool cmp( PointImgTypeDef a,  PointImgTypeDef b) 
{ 
    if(a.val < b.val)
    { 
        return true; 
    } 
    return false; 
} 
bool GetBrPontsInOrder(zxhImageDataT<short> &imgRotCent,vector <PointImgTypeDef> &vCenPont)
{
	int ImgSize[4]={0,0,0,0};
	imgRotCent.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	for(int it=0;it<ImgSize[3];++it)
		for(int iz=0;iz<ImgSize[2];++iz)
			for(int iy=0;iy<ImgSize[1];++iy)
				for(int ix=0;ix<ImgSize[0];++ix)
				{
					short shval=imgRotCent.GetPixelGreyscale(ix,iy,iz,it);
				
					if(shval!=0)
					{	
						PointImgTypeDef pTmpCenPont;
						pTmpCenPont.x=ix;
						pTmpCenPont.y=iy;
						pTmpCenPont.z=iz;
						pTmpCenPont.val=shval;
						vCenPont.push_back(pTmpCenPont);
					}
				}
				sort(vCenPont.begin(),vCenPont.end(),cmp);
	return true;
}
bool BWLABEL3W2(int imgSize[3],std::vector<std::vector<std::vector<int> > > A,std::vector<std::vector<std::vector<int> > > &B,int &NumberOfLabs)
{
	vector<int> stRun;
	vector<int> enRun;
	vector<int> rowRun;
	int NumberOfRuns=0;
	vector<int> runLabels;
		//
	FILLRUNVECTORS(A,NumberOfRuns,stRun,enRun,rowRun);
	//
	if(NumberOfRuns==0)
	{
		NumberOfLabs=0;
		return true;
	}
	vector<pair<int, int>> equivalences;
	int offset=0;
	FIRSTPASS_2(imgSize,stRun,enRun,rowRun,NumberOfRuns,runLabels, equivalences, offset);
	replaceSameLabel(runLabels,equivalences);
	std::vector<int>::iterator biggest = std::max_element(std::begin(runLabels), std::end(runLabels));  
	NumberOfLabs=*biggest;
	BACKLABELEDW(imgSize,A,B,stRun, enRun, rowRun, NumberOfRuns,runLabels);
	return true;
}
bool BWLABEL3W(int imgSize[3],std::vector<std::vector<std::vector<int> > > A,std::vector<std::vector<std::vector<int> > > &B,int &NumberOfLabs)
{
	vector<int> stRun;
	vector<int> enRun;
	vector<int> rowRun;
	int NumberOfRuns=0;
	vector<int> runLabels;
		//
	FILLRUNVECTORS(A,NumberOfRuns,stRun,enRun,rowRun);
	//

	vector<pair<int, int>> equivalences;
	int offset=1;
	FIRSTPASS_1(imgSize,stRun,enRun,rowRun,NumberOfRuns,runLabels, equivalences, offset);
	replaceSameLabel(runLabels,equivalences);
	std::vector<int>::iterator biggest = std::max_element(std::begin(runLabels), std::end(runLabels));  
	NumberOfLabs=*biggest;
	BACKLABELEDW(imgSize,A,B,stRun, enRun, rowRun, NumberOfRuns,runLabels);
	return true;
}
bool Reglabel(int nSOfRe,std::vector<std::vector<std::vector<int> > > Region,int &NumberOfLabs)
	{
			int ImgSize[3]={nSOfRe,nSOfRe,nSOfRe};
			

		return true;
	}
bool Gen_region(int ImgSize[3],std::vector<std::vector<std::vector<int> > >&Region,std::vector<std::vector<std::vector<int> > >Region_Lab)
{
	std::vector<std::vector<std::vector<int> > > vmask;
	vmask.assign(Region_Lab.begin(),Region_Lab.end());
	int xxx=ImgSize[2]/2;
	int nlab=Region_Lab[ImgSize[2]/2][ImgSize[0]/2][ImgSize[1]/2];
	for(int z=0;z<ImgSize[2];z++)  
	{  
		for (int x=0;x<ImgSize[0];x++)  
		{  
			for (int y=0;y<ImgSize[1];y++)  
			{  
				int nval=Region_Lab[z][x][y];
				if(nval==nlab)
				{
					vmask[z][x][y]=nval;
				}
				else
				{
					vmask[z][x][y]=0;
				}
			}
		}
	}
	std::vector<std::vector<std::vector<int> > > vRegion;
	vRegion.assign(Region_Lab.begin(),Region_Lab.end());
	for(int z=0;z<ImgSize[2];z++)   
		for (int x=0;x<ImgSize[0];x++)  
			for (int y=0;y<ImgSize[1];y++)  
			{ 
                vRegion[z][x][y]=Region[z][x][y]*vmask[z][x][y];
			}
			Region.assign(vRegion.begin(),vRegion.end());
	return true;
}
bool Reglabel1(int nSOfRe,std::vector<std::vector<std::vector<int> > > &Region,std::vector<std::vector<std::vector<int> > > &Region_Lab,int &NumberOfLabs)
{
	int ImgSize[3]={nSOfRe,nSOfRe,nSOfRe};

	BWLABEL3W(ImgSize,Region,Region_Lab,NumberOfLabs);

	Gen_region(ImgSize,Region,Region_Lab);
	return true;
}
bool Reglabel2(int nSOfRe,std::vector<std::vector<std::vector<int> > > &Region,std::vector<std::vector<std::vector<int> > > &Region_Lab,int &NumberOfLabs)
{
	int ImgSize[3]={nSOfRe,nSOfRe,nSOfRe};

	BWLABEL3W2(ImgSize,Region,Region_Lab,NumberOfLabs);
	if(NumberOfLabs==0)
		return true;
	//Gen_region(ImgSize,Region,Region_Lab);
	return true;
}
bool Point_anglevec(zxhImageDataT<short> &imgRot,PointImgTypeDef PontImg,int R)
{
	int ImgSize[4]={0,0,0,0};
	imgRot.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	int nSOfRe=2*R+1;
	if(PontImg.x>R&&PontImg.y>R&&PontImg.z>R&&(PontImg.x+R)<=ImgSize[0]&&(PontImg.y+R)<=ImgSize[1]&&(PontImg.z+R)<=ImgSize[2])
	{
		std::vector<std::vector<std::vector<int> > > Region(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		for(int z=0;z<nSOfRe;z++)  
		{  
			for (int x=0;x<nSOfRe;x++)  
			{  
				for (int y=0;y<nSOfRe;y++)  
				{  
					int imgposi[4]={PontImg.x-R+x,PontImg.y-R+y,PontImg.z-R+z,0};
					if(z==0||y==0||x==0||x==nSOfRe-1||y==nSOfRe-1||z==nSOfRe-1)
					{
						Region[z][x][y]=imgRot.GetPixelGreyscale(imgposi[0],imgposi[1],imgposi[2],imgposi[3]);
					}
					else
					{
						Region[z][x][y]=0;
					}
					//cout<<Region[z][x][y];
				}
				//cout<<endl;
			}
		}
		int NumberOfLabs=0;
		
Reglabel(nSOfRe,Region,NumberOfLabs);
	}

	return true;
}
bool Anglevc(int nSOfRe,std::vector<std::vector<std::vector<int> > > Region,std::vector<std::vector<std::vector<int> > >&Region_Angvec)
{
	for(int z=0;z<nSOfRe;z++)  
		for (int x=0;x<nSOfRe;x++)  
			for (int y=0;y<nSOfRe;y++)  
			{  
				if(z==0||y==0||x==0||x==nSOfRe-1||y==nSOfRe-1||z==nSOfRe-1)
				{
					Region_Angvec[z][x][y]=Region[z][x][y];
					if(Region_Angvec[z][x][y]!=0)
						Region_Angvec[z][x][y]=1;
				}
				else
				{
					Region_Angvec[z][x][y]=0;
				}
			}
			return true;
}
bool Point_anglevec1(zxhImageDataT<short> &imgRot,PointImgTypeDef PontImg,int R,int &NumberOfAngvecLabs)
{
	int ImgSize[4]={0,0,0,0};
	imgRot.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	int nSOfRe=2*R+1;
	if(PontImg.x>R&&PontImg.y>R&&PontImg.z>R&&(PontImg.x+R)<=ImgSize[0]&&(PontImg.y+R)<=ImgSize[1]&&(PontImg.z+R)<=ImgSize[2])
	{
		std::vector<std::vector<std::vector<int> > > Region(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		for(int z=0;z<nSOfRe;z++)  
		{  
			for (int x=0;x<nSOfRe;x++)  
			{  
				for (int y=0;y<nSOfRe;y++)  
				{  
					
					int imgposi[4]={PontImg.x-R+x,PontImg.y-R+y,PontImg.z-R+z,0};
					Region[z][x][y]=imgRot.GetPixelGreyscale(imgposi[0],imgposi[1],imgposi[2],imgposi[3]);
				}

			}
		}
		int NumberOfLabs=0;
		std::vector<std::vector<std::vector<int> > > Region_Lab(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		Reglabel1(nSOfRe,Region,Region_Lab,NumberOfLabs);
		std::vector<std::vector<std::vector<int> > > Region_Angvec(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		Anglevc(nSOfRe,Region,Region_Angvec);
		//
		
		int RegSize[3]={nSOfRe,nSOfRe,nSOfRe};
		std::vector<std::vector<std::vector<int> > > Region_Angvec_Lab(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		BWLABEL3W(RegSize,Region_Angvec,Region_Angvec_Lab,NumberOfAngvecLabs);

	}

	return true;
}
bool Point_anglevec2(zxhImageDataT<short> &imgRot,PointImgTypeDef PontImg,int R,int &NumberOfAngvecLabs)
{
	int ImgSize[4]={0,0,0,0};
	imgRot.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	int nSOfRe=2*R+1;
	if(PontImg.x>R&&PontImg.y>R&&PontImg.z>R&&(PontImg.x+R)<=ImgSize[0]&&(PontImg.y+R)<=ImgSize[1]&&(PontImg.z+R)<=ImgSize[2])
	{
		std::vector<std::vector<std::vector<int> > > Region(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		for(int z=0;z<nSOfRe;z++)  
		{  
			for (int x=0;x<nSOfRe;x++)  
			{  
				for (int y=0;y<nSOfRe;y++)  
				{  

					int imgposi[4]={PontImg.x-R+x,PontImg.y-R+y,PontImg.z-R+z,0};
					if(x==R&&y==R&&z==R)
					{
						Region[z][x][y]=0;
					}
					else
					{
						Region[z][x][y]=imgRot.GetPixelGreyscale(imgposi[0],imgposi[1],imgposi[2],imgposi[3]);
					}
					}

			}
		}
		std::vector<std::vector<std::vector<int> > > Region_Lab(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		//Reglabel2(nSOfRe,Region,Region_Lab,NumberOfLabs);
	
	//std::vector<std::vector<std::vector<int> > > Region_Angvec(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		//Anglevc(nSOfRe,Region,Region_Angvec);
		//
		
		int RegSize[3]={nSOfRe,nSOfRe,nSOfRe};
		std::vector<std::vector<std::vector<int> > > Region_Angvec_Lab(nSOfRe,vector<vector<int> >(nSOfRe,vector<int>(nSOfRe,0)));  
		BWLABEL3W2(RegSize,Region,Region_Lab,NumberOfAngvecLabs);

	}

	return true;
}
bool Calc_ConRegOfNei(zxhImageDataT<short>&imgRot,zxhImageDataT<short>&imgCountConNeimap,vector<PointImgTypeDef>&vSeedsPonts)
{
	int ImgNewSize[4]={1};
	imgRot.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for (int i=0;i<vSeedsPonts.size();i++)//map the vetor points to the image
	{
		PointImgTypeDef pseed;
		pseed=vSeedsPonts[i];
		int R=1;
		int AngNum=1;
		if(pseed.x==160&&pseed.y==83&&pseed.z==109)
		{
			int x=0;
		}
		Point_anglevec2(imgRot,pseed,R,AngNum);
		imgCountConNeimap.SetPixelByGreyscale(pseed.x,pseed.y,pseed.z,0,AngNum);
	}
	char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_Nei_Con.nii.gz";
	string chFileName2(chResultName);
	zxh::SaveImage(&imgCountConNeimap,chFileName2.c_str());

	return true;
}
bool Point_select(zxhImageDataT<short> &imgRot,vector <PointImgTypeDef> &vCenPont,int R,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vptroot)
{
	for (int i=0;i<vCenPont.size();i++)
	{
		PointImgTypeDef PontImg;
		PontImg=vCenPont[i];
		int AngNum=0;
		if(i==40)
		{
			int xxxxx=0;
		}
		Point_anglevec1(imgRot,PontImg,R,AngNum);
		if(AngNum==3)
		{
			PointImgTypeDef tmpPont;
			tmpPont=PontImg;
			tmpPont.val=vpbitu.size()+1;
			vpbitu.push_back(tmpPont);
		}
		if(AngNum==1)
		{
			PointImgTypeDef tmpPont;
			tmpPont=PontImg;
			tmpPont.val=vptroot.size()+1;
			vptroot.push_back(tmpPont);
		}
	}
	return true;
}
bool Point_show(zxhImageDataT<short>&imgRot_cormarg,vector<PointImgTypeDef> &vpbitu,int R)
{

	 for(int i=0;i<vpbitu.size();i++)
	 {
		 PointImgTypeDef tmpPont;
		 tmpPont=vpbitu[i];
		 imgRot_cormarg.SetPixelByGreyscale(tmpPont.x,tmpPont.y,tmpPont.z,0,tmpPont.val);
		 for (int nz = tmpPont.z-R; nz <=tmpPont.z+R; nz++)
			 for (int nx = tmpPont.x-R; nx <=tmpPont.x+R; nx++)
				 for (int ny = tmpPont.y-R; ny <=tmpPont.y+R; ny++)
				 {
					 if(nx==tmpPont.x-R||nx==tmpPont.x+R||ny==tmpPont.y-R||ny==tmpPont.y+R||nz==tmpPont.z-R||nz==tmpPont.z+R)
					 {
						 imgRot_cormarg.SetPixelByGreyscale(nx,ny,nz,0,tmpPont.val);
					 }
				 }
	 }
	// char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_CentMargLab.nii.gz";
	// string chFileName2(chResultName);
	// zxh::SaveImage(&imgRot_cormarg,chFileName2.c_str());
	 return true;
}
bool Point_show_biro(zxhImageDataT<short>&imgRot_cormarg,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,int &R)
{

	 for(int i=0;i<vpbitu.size();i++)
	 {
		 PointImgTypeDef tmpPont;
		 tmpPont=vpbitu[i];
		 imgRot_cormarg.SetPixelByGreyscale(tmpPont.x,tmpPont.y,tmpPont.z,0,tmpPont.val);
		 for (int nz = tmpPont.z-R; nz <=tmpPont.z+R; nz++)
			 for (int nx = tmpPont.x-R; nx <=tmpPont.x+R; nx++)
				 for (int ny = tmpPont.y-R; ny <=tmpPont.y+R; ny++)
				 {
					 if(nx==tmpPont.x-R||nx==tmpPont.x+R||ny==tmpPont.y-R||ny==tmpPont.y+R||nz==tmpPont.z-R||nz==tmpPont.z+R)
					 {
						 imgRot_cormarg.SetPixelByGreyscale(nx,ny,nz,0,tmpPont.val);
					 }
				 }
	 }
	  for(int i=0;i<vproot.size();i++)
	 {
		 PointImgTypeDef tmpPont;
		 tmpPont=vproot[i];
		 int roval=vpbitu.size()+tmpPont.val;
		 imgRot_cormarg.SetPixelByGreyscale(tmpPont.x,tmpPont.y,tmpPont.z,0,roval);
		 for (int nz = tmpPont.z-R; nz <=tmpPont.z+R; nz++)
			 for (int nx = tmpPont.x-R; nx <=tmpPont.x+R; nx++)
				 for (int ny = tmpPont.y-R; ny <=tmpPont.y+R; ny++)
				 {
					 if(nx==tmpPont.x-R||nx==tmpPont.x+R||ny==tmpPont.y-R||ny==tmpPont.y+R||nz==tmpPont.z-R||nz==tmpPont.z+R)
					 {
						 imgRot_cormarg.SetPixelByGreyscale(nx,ny,nz,0,roval);
					 }
				 }
	 }
	 char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_CentMargLab1.nii.gz";
	 string chFileName2(chResultName);
	 zxh::SaveImage(&imgRot_cormarg,chFileName2.c_str());
	 return true;
}
bool Point_show_ro(zxhImageDataT<short>&imgRot_cormarg,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,int &R)
{

	  for(int i=0;i<vproot.size();i++)
	 {
		 PointImgTypeDef tmpPont;
		 tmpPont=vproot[i];
		 int roval=vpbitu.size()+tmpPont.val;
		 imgRot_cormarg.SetPixelByGreyscale(tmpPont.x,tmpPont.y,tmpPont.z,0,roval);
		 for (int nz = tmpPont.z-R; nz <=tmpPont.z+R; nz++)
			 for (int nx = tmpPont.x-R; nx <=tmpPont.x+R; nx++)
				 for (int ny = tmpPont.y-R; ny <=tmpPont.y+R; ny++)
				 {
					 if(nx==tmpPont.x-R||nx==tmpPont.x+R||ny==tmpPont.y-R||ny==tmpPont.y+R||nz==tmpPont.z-R||nz==tmpPont.z+R)
					 {
						 imgRot_cormarg.SetPixelByGreyscale(nx,ny,nz,0,roval);
					 }
				 }
	 }
	 char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_RootMargLab1.nii.gz";
	 string chFileName2(chResultName);
	 zxh::SaveImage(&imgRot_cormarg,chFileName2.c_str());
	 return true;
}

bool SetZeroinmap(zxhImageDataT<short> &imgRot_cormarg,PointImgTypeDef PontSeed)
{
	int ImgSize[4]={0,0,0,0};
	imgRot_cormarg.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	for(int it=0;it<ImgSize[3];++it)
		for(int iz=0;iz<ImgSize[2];++iz)
			for(int iy=0;iy<ImgSize[1];++iy)
				for(int ix=0;ix<ImgSize[0];++ix)
				{
					short sint=imgRot_cormarg.GetPixelGreyscale(ix,iy,iz,it);
					if(sint==PontSeed.val)
					{
						imgRot_cormarg.SetPixelByGreyscale(ix,iy,iz,it,0);
					}
				}
	return true;
}
bool SetZero(zxhImageDataT<short> &imgRot,PointImgTypeDef PontSeed)
{
	imgRot.SetPixelByGreyscale(PontSeed.x,PontSeed.y,PontSeed.z,0,0);
	return true;
}
bool FindCentPont(int nintmarg,vector<PointImgTypeDef> &vpbitu,PointImgTypeDef &CentPont)
{
	for(int i=0;i<vpbitu.size();i++)
	{
		PointImgTypeDef tmp;
		tmp=vpbitu[i];
		if(nintmarg==tmp.val)
		{
			CentPont=tmp;
			break;
		}
	}

	return true;
}
bool FindCentPont_biro(int nintmarg,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,PointImgTypeDef &CentPont)
{
	CentPont.val=-1;
	for(int i=0;i<vpbitu.size();i++)
	{
		PointImgTypeDef tmp;
		tmp=vpbitu[i];
		if(nintmarg==tmp.val)
		{
			CentPont=tmp;
			break;
		}
	}
	if(CentPont.val==-1)
	{
		for(int i=0;i<vproot.size();i++)
	{
		PointImgTypeDef tmp;
		tmp=vproot[i];
		if(nintmarg==vpbitu.size()+tmp.val)
		{
			CentPont=tmp;
			CentPont.val=vpbitu.size()+tmp.val;
			break;
		}
	}
	}

	return true;
}
bool FindCentPont_ro(int nintmarg,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,vector<PointLinkDef> &vbitulink,PointLinkDef &tmplink,PointImgTypeDef &CentPont)
{
	
		for(int i=0;i<vpbitu.size();i++)
	{
		PointImgTypeDef tmp;
		tmp=vpbitu[i];
		if(nintmarg==tmp.val)
		{
			return false;
			
		}
	}

		for(int i=0;i<vproot.size();i++)
	{
		PointImgTypeDef tmp;
		tmp=vproot[i];
		if(nintmarg==tmp.val)
		{
			CentPont=tmp;
			CentPont.val=vpbitu.size()+tmp.val;
			break;
		}
	}
	

	return true;
}
bool CopyImg(zxhImageDataT<short> &imgRot,zxhImageDataT<short> &imgRot_copy)
{
	int ImgSize[4]={0,0,0,0};
	imgRot_copy.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	for(int it=0;it<ImgSize[3];++it)
		for(int iz=0;iz<ImgSize[2];++iz)
			for(int iy=0;iy<ImgSize[1];++iy)
				for(int ix=0;ix<ImgSize[0];++ix)
				{
					short sint=imgRot.GetPixelGreyscale(ix,iy,iz,it);
					imgRot_copy.SetPixelByGreyscale(ix,iy,iz,it,sint);
				}
				return true;
}
bool Point_neighbor(zxhImageDataT<short> &imgRot,vector<PointImgTypeDef> &vpbitu,PointImgTypeDef PontSeed,zxhImageDataT<short>&imgRot_cormarg,int NeighborNum,int R,PointLinkDef &tmplink)
{
		int gNbr[26][3] = { 
		{-1, 0, 0}, \
		{-1, -1, 0}, \
		{-1, 1, 0}, \
		{-1, 0, -1}, \
		{-1, -1, -1}, \
		{-1, -1, 1}, \
		{-1, 1, -1}, \
		{-1, 1, 1}, \
		{-1, 0, 1}, \
		{1, 0, 0}, \
		{1, -1, 0}, \
		{1, 1, 0}, \
		{1, 0, -1}, \
		{1, -1, -1}, \
		{1, -1, 1}, \
		{1, 1, -1}, \
		{1, 1, 1}, \
		{1, 0, 1}, \
		{ 0,-1, 0}, \
		{ 0, 1, 0}, \
		{ 0, 0,-1}, \
		{ 0, 0, 1},\
		{ 0, -1,-1}, \
		{ 0, -1,1}, \
		{ 0, 1,-1}, \
		{ 0, 1,1}, \
	};
		SetZeroinmap(imgRot_cormarg,PontSeed);
		SetZero(imgRot,PontSeed);
		vector<PointImgTypeDef> vmap;
		vmap.push_back(PontSeed);
		int ImgSize[4]={0,0,0,0};
		imgRot_cormarg.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
		int ncout=0;
		vector<PointImgTypeDef> vlink;
		tmplink.pont=PontSeed;
		while(!vmap.empty()&&ncout<NeighborNum)
		{
			PontSeed=vmap[0];
			for (int i = 0; i < 26; i++)
			{
				int nx = PontSeed.x + gNbr[i][0];
				int ny = PontSeed.y + gNbr[i][1];
				int nz = PontSeed.z + gNbr[i][2];
				int nint=imgRot.GetPixelGreyscale(nx,ny,nz,0);
				PointImgTypeDef tmp;
				tmp.x=nx;
				tmp.y=ny;
				tmp.z=nz;
				tmp.val=nint;
				if(nx>=0&&nx<ImgSize[0]&&ny>=0&&ny<ImgSize[1]&&nz>=0&&nz<ImgSize[2]&&nint!=0)
				{

					SetZero(imgRot,tmp);
					int nintmarg=imgRot_cormarg.GetPixelGreyscale(nx,ny,nz,0);
					if(nintmarg!=0)
					{
						PointImgTypeDef CentPont;
						FindCentPont(nintmarg,vpbitu,CentPont);
						vlink.push_back(CentPont);
						ncout++;
						if(ncout>=NeighborNum)
							break;
						SetZeroinmap(imgRot_cormarg,CentPont);
						for (int i = 0; i < 26; i++)
						{
							int nx1 = tmp.x + gNbr[i][0];
							int ny1 = tmp.y + gNbr[i][1];
							int nz1 = tmp.z + gNbr[i][2];
							if(nx1>0&&nx1<ImgSize[0]&&ny1>0&&ny1<ImgSize[1]&&nz1>0&&nz1<ImgSize[2])
							{
								PointImgTypeDef tmp1;
								tmp1.x=nx1;
								tmp1.y=ny1;
								tmp1.z=nz1;
								SetZero(imgRot,tmp1);
							}
						}
					}
					else
					{
						vmap.push_back(tmp);
					}
				}
			}
			vector<PointImgTypeDef>::iterator k = vmap.begin();
			vmap.erase(k);
		}
		tmplink.cunt=ncout;
		tmplink.vpont=vlink;
	return true;
}
bool Point_neighbor_biro(zxhImageDataT<short> &imgRot,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,PointImgTypeDef PontSeed,zxhImageDataT<short>&imgRot_cormarg,int NeighborNum,int R,PointLinkDef &tmplink)
{
		int gNbr[26][3] = { 
		{-1, 0, 0}, \
		{-1, -1, 0}, \
		{-1, 1, 0}, \
		{-1, 0, -1}, \
		{-1, -1, -1}, \
		{-1, -1, 1}, \
		{-1, 1, -1}, \
		{-1, 1, 1}, \
		{-1, 0, 1}, \
		{1, 0, 0}, \
		{1, -1, 0}, \
		{1, 1, 0}, \
		{1, 0, -1}, \
		{1, -1, -1}, \
		{1, -1, 1}, \
		{1, 1, -1}, \
		{1, 1, 1}, \
		{1, 0, 1}, \
		{ 0,-1, 0}, \
		{ 0, 1, 0}, \
		{ 0, 0,-1}, \
		{ 0, 0, 1},\
		{ 0, -1,-1}, \
		{ 0, -1,1}, \
		{ 0, 1,-1}, \
		{ 0, 1,1}, \
	};
		SetZeroinmap(imgRot_cormarg,PontSeed);
		SetZero(imgRot,PontSeed);
		vector<PointImgTypeDef> vmap;
		vmap.push_back(PontSeed);
		int ImgSize[4]={0,0,0,0};
		imgRot_cormarg.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
		int ncout=0;
		vector<PointImgTypeDef> vlink;
		tmplink.pont=PontSeed;
		while(!vmap.empty()&&ncout<NeighborNum)
		{
			PontSeed=vmap[0];
			for (int i = 0; i < 26; i++)
			{
				int nx = PontSeed.x + gNbr[i][0];
				int ny = PontSeed.y + gNbr[i][1];
				int nz = PontSeed.z + gNbr[i][2];
				int nint=imgRot.GetPixelGreyscale(nx,ny,nz,0);
				PointImgTypeDef tmp;
				tmp.x=nx;
				tmp.y=ny;
				tmp.z=nz;
				tmp.val=nint;
				if(nx>=0&&nx<ImgSize[0]&&ny>=0&&ny<ImgSize[1]&&nz>=0&&nz<ImgSize[2]&&nint!=0)
				{

					SetZero(imgRot,tmp);
					int nintmarg=imgRot_cormarg.GetPixelGreyscale(nx,ny,nz,0);
					if(nintmarg!=0)
					{
						PointImgTypeDef CentPont;
						FindCentPont_biro(nintmarg,vpbitu,vproot,CentPont);
						vlink.push_back(CentPont);
						ncout++;
						if(ncout==NeighborNum)
							break;
						SetZeroinmap(imgRot_cormarg,CentPont);
						for (int i = 0; i < 26; i++)
						{
							int nx1 = tmp.x + gNbr[i][0];
							int ny1 = tmp.y + gNbr[i][1];
							int nz1 = tmp.z + gNbr[i][2];
							if(nx1>0&&nx1<ImgSize[0]&&ny1>0&&ny1<ImgSize[1]&&nz1>0&&nz1<ImgSize[2])
							{
								PointImgTypeDef tmp1;
								tmp1.x=nx1;
								tmp1.y=ny1;
								tmp1.z=nz1;
								SetZero(imgRot,tmp1);
							}
						}
					}
					else
					{
						vmap.push_back(tmp);
					}
				}
			}
			vector<PointImgTypeDef>::iterator k = vmap.begin();
			vmap.erase(k);
		}
		tmplink.cunt=ncout;
		tmplink.vpont=vlink;
	return true;
}
bool Point_neighbor_bifurroot(zxhImageDataT<short> &imgRot,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,PointImgTypeDef PontSeed,zxhImageDataT<short>&imgRot_cormarg,int NeighborNum,int R,vector<PointLinkDef> &vbitulink,PointLinkDef &tmplink)
{
		int gNbr[26][3] = { 
		{-1, 0, 0}, \
		{-1, -1, 0}, \
		{-1, 1, 0}, \
		{-1, 0, -1}, \
		{-1, -1, -1}, \
		{-1, -1, 1}, \
		{-1, 1, -1}, \
		{-1, 1, 1}, \
		{-1, 0, 1}, \
		{1, 0, 0}, \
		{1, -1, 0}, \
		{1, 1, 0}, \
		{1, 0, -1}, \
		{1, -1, -1}, \
		{1, -1, 1}, \
		{1, 1, -1}, \
		{1, 1, 1}, \
		{1, 0, 1}, \
		{ 0,-1, 0}, \
		{ 0, 1, 0}, \
		{ 0, 0,-1}, \
		{ 0, 0, 1},\
		{ 0, -1,-1}, \
		{ 0, -1,1}, \
		{ 0, 1,-1}, \
		{ 0, 1,1}, \
	};
		SetZeroinmap(imgRot_cormarg,PontSeed);
		SetZero(imgRot,PontSeed);
		vector<PointImgTypeDef> vmap;
		vmap.push_back(PontSeed);
		int ImgSize[4]={0,0,0,0};
		imgRot_cormarg.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
		int ncout=0;
		vector<PointImgTypeDef> vlink;
		tmplink.pont=PontSeed;
		while(!vmap.empty()&&ncout<NeighborNum)
		{
			PontSeed=vmap[0];
			for (int i = 0; i < 26; i++)
			{
				int nx = PontSeed.x + gNbr[i][0];
				int ny = PontSeed.y + gNbr[i][1];
				int nz = PontSeed.z + gNbr[i][2];
				int nint=imgRot.GetPixelGreyscale(nx,ny,nz,0);
				PointImgTypeDef tmp;
				tmp.x=nx;
				tmp.y=ny;
				tmp.z=nz;
				tmp.val=nint;
				if(nx>=0&&nx<ImgSize[0]&&ny>=0&&ny<ImgSize[1]&&nz>=0&&nz<ImgSize[2]&&nint!=0)
				{

					SetZero(imgRot,tmp);
					int nintmarg=imgRot_cormarg.GetPixelGreyscale(nx,ny,nz,0);
					if(nintmarg!=0)
					{
						PointImgTypeDef CentPont;
						FindCentPont_ro(nintmarg,vpbitu,vproot,vbitulink,tmplink,CentPont);
						vlink.push_back(CentPont);
						ncout++;
						if(ncout>=NeighborNum)
							break;
						SetZeroinmap(imgRot_cormarg,CentPont);
						for (int i = 0; i < 26; i++)
						{
							int nx1 = tmp.x + gNbr[i][0];
							int ny1 = tmp.y + gNbr[i][1];
							int nz1 = tmp.z + gNbr[i][2];
							if(nx1>0&&nx1<ImgSize[0]&&ny1>0&&ny1<ImgSize[1]&&nz1>0&&nz1<ImgSize[2])
							{
								PointImgTypeDef tmp1;
								tmp1.x=nx1;
								tmp1.y=ny1;
								tmp1.z=nz1;
								SetZero(imgRot,tmp1);
							}
						}
					}
					else
					{
						vmap.push_back(tmp);
					}
				}
			}
			vector<PointImgTypeDef>::iterator k = vmap.begin();
			vmap.erase(k);
		}
		tmplink.cunt=ncout;
		tmplink.vpont=vlink;
	return true;
}
bool Point_link(zxhImageDataT<short> &imgRot,zxhImageDataT<short>&imgRotCent,vector<PointImgTypeDef> &vpbitu,int R,vector<PointLinkDef> &vbitulink)
{
	zxhImageDataT<short> imgRot_cormarg;
	imgRot_cormarg.NewImage(imgRot.GetImageInfo());
	Point_show(imgRot_cormarg,vpbitu,R);
	int NeighborNum=3;
	zxhImageDataT<short> imgRot_copy,imgRot_cormarg_copy;
		imgRot_copy.NewImage(imgRot.GetImageInfo());
		imgRot_cormarg_copy.NewImage(imgRot_cormarg.GetImageInfo());
	for(int i=0;i<vpbitu.size();i++)
	{
		PointImgTypeDef PontSeed;
		PontSeed=vpbitu[i];
		PointLinkDef tmplink;
		CopyImg(imgRot,imgRot_copy);
		CopyImg(imgRot_cormarg,imgRot_cormarg_copy);
		Point_neighbor(imgRot,vpbitu,PontSeed,imgRot_cormarg,NeighborNum,R,tmplink);
		vbitulink.push_back(tmplink);
		CopyImg(imgRot_copy,imgRot);
		CopyImg(imgRot_cormarg_copy,imgRot_cormarg);
		
	}
	return true;
}

bool Point_link_biro(zxhImageDataT<short> &imgRot,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,int &R,vector<PointLinkDef> &vbitulink)
{
	zxhImageDataT<short> imgRot_cormarg;
	imgRot_cormarg.NewImage(imgRot.GetImageInfo());
	Point_show_biro(imgRot_cormarg,vpbitu,vproot,R);
	int NeighborNum=3;
	zxhImageDataT<short> imgRot_copy,imgRot_cormarg_copy;
		imgRot_copy.NewImage(imgRot.GetImageInfo());
		imgRot_cormarg_copy.NewImage(imgRot_cormarg.GetImageInfo());
	for(int i=0;i<vpbitu.size();i++)
	{
		PointImgTypeDef PontSeed;
		PontSeed=vpbitu[i];
		PointLinkDef tmplink;
		CopyImg(imgRot,imgRot_copy);
		CopyImg(imgRot_cormarg,imgRot_cormarg_copy);
		Point_neighbor_biro(imgRot,vpbitu,vproot,PontSeed,imgRot_cormarg,NeighborNum,R,tmplink);
		vbitulink.push_back(tmplink);
		CopyImg(imgRot_copy,imgRot);
		CopyImg(imgRot_cormarg_copy,imgRot_cormarg);
		
	}
	return true;
}
bool Point_link_ro(zxhImageDataT<short> &imgRot,vector<PointImgTypeDef> &vpbitu,vector<PointImgTypeDef> &vproot,int &R,vector<PointLinkDef> &vbitulink,vector<PointLinkDef> &vbiturootlink)
{
	zxhImageDataT<short> imgRot_cormarg;
	imgRot_cormarg.NewImage(imgRot.GetImageInfo());
	Point_show_ro(imgRot_cormarg,vpbitu,vproot,R);
	int NeighborNum=3;
	zxhImageDataT<short> imgRot_copy,imgRot_cormarg_copy;
		imgRot_copy.NewImage(imgRot.GetImageInfo());
		imgRot_cormarg_copy.NewImage(imgRot_cormarg.GetImageInfo());
	for(int i=0;i<vpbitu.size();i++)
	{
		PointImgTypeDef PontSeed;
		PontSeed=vpbitu[i];
		PointLinkDef tmplink;
		CopyImg(imgRot,imgRot_copy);
		CopyImg(imgRot_cormarg,imgRot_cormarg_copy);
		Point_neighbor_bifurroot(imgRot,vpbitu,vproot,PontSeed,imgRot_cormarg,NeighborNum,R,vbitulink,tmplink);
		vbiturootlink.push_back(tmplink);
		CopyImg(imgRot_copy,imgRot);
		CopyImg(imgRot_cormarg_copy,imgRot_cormarg);
		
	}
	return true;
}


bool BifurDec(zxhImageDataT<short>&imgReadNewRaw)
{
	
	zxhImageDataT<short>imgCountmap,imgCountConNeimap;
	imgCountmap.NewImage( imgReadNewRaw.GetImageInfo() );
	init_img(imgCountmap);
	vector<PointImgTypeDef> vSeedsPonts;
	Points_Init1(imgReadNewRaw,imgCountmap,vSeedsPonts);
	//count the number of neighbor points
	imgCountConNeimap.NewImage(imgReadNewRaw.GetImageInfo() );
	init_img(imgCountConNeimap);
	Calc_NumOfNei(imgReadNewRaw,imgCountmap,vSeedsPonts);
	//calculate the connected region of the neighbor poins
	Calc_ConRegOfNei(imgReadNewRaw,imgCountConNeimap,vSeedsPonts);
	int ImgSize[4]={1};
	imgCountmap.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	std::vector<std::vector<std::vector<int> > > BW(ImgSize[2],vector<vector<int> >(ImgSize[0],vector<int>(ImgSize[1],0)));  
	std::vector<std::vector<std::vector<int> > > BW1(ImgSize[2],vector<vector<int> >(ImgSize[0],vector<int>(ImgSize[1],0)));  
	int zz=BW1.size();
	int zz1=BW.size();
	//mark the candidate points
	Coutmap_vec1(imgCountmap,imgCountConNeimap,BW,BW1);
	zxhImageDataT<short>imgLabel;
	imgLabel.NewImage( imgReadNewRaw.GetImageInfo() );
		int NumberOfLabs=0;
	BWLABEL3(ImgSize,BW,BW1,NumberOfLabs,imgLabel);
	//select the center point
	std::vector<std::vector<std::vector<int> > > BW2(ImgSize[2],vector<vector<int> >(ImgSize[0],vector<int>(ImgSize[1],0)));  
	zxhImageDataT<short> imgCenLabel;
	imgCenLabel.NewImage( imgReadNewRaw.GetImageInfo() );
	SelCen(ImgSize,BW1,BW2,NumberOfLabs,imgCenLabel);
	return true;
}
int main(int argc, char *argv[])
{
	//if( argc < 4 )
	//{
	//	cerr << "Usage: " << endl;
	//	cerr << "zxhcaeDMPModelIntenGen	RawImage(.nii)	ModelRef(.vtk) MapModelFile(.vtk) ResultPath " << endl;
	//	return -1;
	//} 
	//
	//string strFileNameRaw =string(argv[1]);
	//char *chModelFileName =argv[2];
	//char *chMapModelRawFileName=argv[3];
	//char *chMapMolInteFileName=argv[4];
	
	//read and resample the curve
	//--..--..--..--..
     string strfilenameraw =  "J:/JDQ/CCTA_CAR/RCAA_32/training/dataset00/CAE_ME_L.nii.gz";
	char *chRefCurvefilename ="J:/JDQ/CCTA_CAR/RCAA_32/training/dataset00/CL0.vtk";
	char *chTarCurvefilename="J:/JDQ/CCTA_CAR/RCAA_32/training/dataset00/CL1.vtk";


	
	char *chCurvesfolder="F:/Coronary_0/code/Resutsfusion/DS_For_ResFu/dataset00/";
	vector<vLinesDef> vnLine;
	ReadFilesInAfolder(chCurvesfolder,vnLine);
	

	float flength=0,flengthdeform=0,fsectionlength=0;
	float ResampleNUM=9999;
	//read the raw image
	zxhImageDataT<short> imgReadRaw;
	if( zxh::OpenImage( &imgReadRaw, strfilenameraw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"; 
		return -1;
	}
	std::vector<jdq2017::point3D> ref;
	//std::vector<jdq2017::point3D> cl;
	if (! jdq2017::readCenterlinevtk(chRefCurvefilename, ref) )//read the reference line
	{
		if (!outputOneLine) 
		{
			std::cerr << "Error in reading input data" << std::endl;
			return 1;
		} 

	}


	//calculate the distance step
	int nSamNUM=500;
	double clSampling = pathLength(ref)/nSamNUM;

	ResamEveryCurve(vnLine,clSampling);
	
	//Transform all the points of lines into unorganized points
	vector<PointCordTypeDef> vUnorgaPointsWorld;
	TransAllPontsIntoUnorga(vnLine,vUnorgaPointsWorld);
	//Store the unique points
	vector<PointCordTypeDef> vUnorgaPointsWorld_ori;
	StorUnique(vUnorgaPointsWorld);	
	//first fusion
	vUnorgaPointsWorld_ori.assign(vUnorgaPointsWorld.begin(),vUnorgaPointsWorld.end());
	
	char *Points_3D_ori_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D_ori_n.txt";
	//WriteCA2Txt_Skip(vUnorgaPointsWorld_ori,Points_3D_ori_Filename);
	WriteCA2Txt(vUnorgaPointsWorld_ori,Points_3D_ori_Filename);
	//adjust every points by section
	AdEvePontsBySec_ALL(clSampling,vUnorgaPointsWorld,vnLine);
	
	char *Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D_new.txt";
	//WriteCA2Txt_Skip(vUnorgaPointsWorld,Points_3D_Filename);
	WriteCA2Txt(vUnorgaPointsWorld,Points_3D_Filename);
	//second fusion
	//char *Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D_new.txt";
	//vector<PointCordTypeDef> vAdpoints;
	//ReadCA2Txt(vAdpoints,Points_3D_Filename);
	//AdEvePontsBySec_ALL(clSampling,vUnorgaPointsWorld, vnLine);

	//--..--..--..--..
	//----------------------
	

	//
	
    //---...----...----....
	//get the central points
	//string strfilenameraw =  "J:/JDQ/CCTA_CAR/RCAA_32/training/dataset00/CAE_ME_L.nii.gz";
	////read the raw image
	//zxhImageDataT<short> imgReadRaw;
	//if( zxh::OpenImage( &imgReadRaw, strfilenameraw ) == false )
	//{
	//	std::cerr << "Raw image(nifti-file) is not found!"; 
	//	return -1;
	//}
	//char *Points_3D_Filename="F:/Coronary_0/code/Resutsfusion/Points_3D_new.txt";
	//vector<PointCordTypeDef> vAdpoints;
	//ReadCA2Txt(vAdpoints,Points_3D_Filename);
	//zxhImageDataT<short>imgReadNewRaw;
	//imgReadNewRaw.NewImage( imgReadRaw.GetImageInfo() );
	//MapModelPointsToImage(imgReadNewRaw,vAdpoints);
	////store the points as image
	//char *chResultName="F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot.nii.gz";
	//string chFileName2(chResultName);
	//zxh::SaveImage(&imgReadNewRaw,chFileName2.c_str());


	//vector<PointCordTypeDef> vBifurpoints;
	//BifurDec(imgReadNewRaw);
	//---...----...----....
	
	



//----....----....----....----....----....----....
//get the link of the bifurcation points
//string strfilenameraw_rot =  "F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot.nii.gz";
//zxhImageDataT<short> imgRot;
//if( zxh::OpenImage( &imgRot, strfilenameraw_rot ) == false )
//{
//	std::cerr << "Raw image(nifti-file) is not found!"; 
//	return -1;
//}
//
//string strfilenameraw_rotCent =  "F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_CentLab.nii.gz";
//zxhImageDataT<short> imgRotCent;
//if( zxh::OpenImage( &imgRotCent, strfilenameraw_rotCent ) == false )
//{
//	std::cerr << "Raw image(nifti-file) is not found!"; 
//	return -1;
//}
////get the branch points in order
//
//vector <PointImgTypeDef> vCenPont;
//GetBrPontsInOrder(imgRotCent,vCenPont);
//int R=3;
//vector<PointImgTypeDef> vpbitu;
//	vector<PointImgTypeDef> vptroot;
//Point_select(imgRot,vCenPont,R,vpbitu,vptroot);
//vector<PointLinkDef> vbitulink;
//Point_link(imgRot,imgRotCent,vpbitu,R,vbitulink);
////save bifurcation points as img
//zxhImageDataT<short>imgVPbituRaw,imgVProotRaw;
//imgVPbituRaw.NewImage( imgRot.GetImageInfo() );
//MapImgPointsToImage(imgVPbituRaw,vpbitu);
//char *chResultName1="F:/Coronary_0/code/Resutsfusion/CAE_VPbitu.nii.gz";
//string chFileName1(chResultName1);
//zxh::SaveImage(&imgVPbituRaw,chFileName1.c_str());
////
//imgVProotRaw.NewImage( imgRot.GetImageInfo() );
//MapImgPointsToImage(imgVProotRaw,vptroot);
//char *chResultName2="F:/Coronary_0/code/Resutsfusion/CAE_VProot.nii.gz";
//string chFileName2(chResultName2);
//zxh::SaveImage(&imgVProotRaw,chFileName2.c_str());

//vector <PointImgTypeDef> vCenPont;
//GetBrPontsInOrder(imgRotCent,vCenPont);
//int R=3;
//	vector<PointImgTypeDef> vpbitu;
//	vector<PointImgTypeDef> vptroot;
//Point_select(imgRot,vCenPont,R,vpbitu,vptroot);
//
//Point_link(imgRot,R,vpbitu,vptroot);
//----....----....----....----....----....----....
//-----.....-----.....-----.....-----.....-----.....-----.....-----.....
//get the link of the bifurcation points and root points
//string strfilenameraw_rot =  "F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot.nii.gz";
//zxhImageDataT<short> imgRot;
//if( zxh::OpenImage( &imgRot, strfilenameraw_rot ) == false )
//{
//	std::cerr << "Raw image(nifti-file) is not found!"; 
//	return -1;
//}
//
//string strfilenameraw_rotCent =  "F:/Coronary_0/code/Resutsfusion/CAE_ME_L_Rot_CentLab.nii.gz";
//zxhImageDataT<short> imgRotCent;
//if( zxh::OpenImage( &imgRotCent, strfilenameraw_rotCent ) == false )
//{
//	std::cerr << "Raw image(nifti-file) is not found!"; 
//	return -1;
//}
////get the branch points in order
//
//vector <PointImgTypeDef> vCenPont;
//GetBrPontsInOrder(imgRotCent,vCenPont);
//int R=3;
//vector<PointImgTypeDef> vpbitu;
//	vector<PointImgTypeDef> vptroot;
//Point_select(imgRot,vCenPont,R,vpbitu,vptroot);
//vector<PointLinkDef> vbitulink;
//Point_link(imgRot,imgRotCent,vpbitu,R,vbitulink);
//vector<PointLinkDef> vbiturootlink;
//Point_link_ro(imgRot,vpbitu,vptroot,R,vbitulink,vbiturootlink);


//-----.....-----.....-----.....-----.....-----.....-----.....-----.....

/*
//testing BWLABEL3

	int NZ=1;//hight
int NX=4;//row
int NY=14;//col
std::vector<std::vector<std::vector<int> > > grdarray(NZ,vector<vector<int> >(NX,vector<int>(NY,0)));  
    for(int z=0;z<NZ;z++)  
    {  
        for (int x=0;x<NX;x++)  
        {  
            for (int y=0;y<NY;y++)  
            {  
                grdarray[z][x][y]=0;  
                
            }  
        }  
    }  
	grdarray[0][0][1]=1; 
	grdarray[0][0][3]=1; 
	grdarray[0][0][4]=1; 
	grdarray[0][0][8]=1; 
	grdarray[0][0][9]=1; 
	grdarray[0][0][10]=1; 
	grdarray[0][0][13]=1; 

	grdarray[0][1][2]=1; 
	grdarray[0][1][6]=1; 


	grdarray[0][2][0]=1; 
	grdarray[0][2][3]=1; 
	grdarray[0][2][4]=1; 
	grdarray[0][2][6]=1; 
	grdarray[0][2][7]=1; 
	grdarray[0][2][13]=1; 

	grdarray[0][3][0]=1; 
	grdarray[0][3][13]=1; 
	  for(int z=0;z<NZ;z++)  
    {  
        for (int x=0;x<NX;x++)  
        {  
            for (int y=0;y<NY;y++)  
            {  
                cout<<grdarray[z][x][y];  
                
            }  
			cout<<endl; 
        }  
    }  
BWLABEL3(grdarray);
//testing BWLABEL3
*/

				

//create 3D vector
//int NZ=3;//hight
//int NX=3;//row
//int NY=3;//col
//std::vector<std::vector<std::vector<int> > > grdarray(NZ,vector<vector<int> >(NX,vector<int>(NY,0)));  
//    for(int z=0;z<NZ;z++)  
//    {  
//        for (int y=0;y<NX;y++)  
//        {  
//            for (int x=0;x<NY;x++)  
//            {  
//                grdarray[z][x][y]=0;  
//                
//            }  
//        }  
//    }  
//	grdarray[0][0][0]=1; 
//	grdarray[0][0][1]=1; 
//	grdarray[0][1][1]=1; 
//	grdarray[0][2][0]=1; 
//
//	grdarray[1][1][1]=1; 
//	grdarray[1][2][1]=1; 
//
//	grdarray[2][0][1]=1; 
//	grdarray[2][0][2]=1; 
//	grdarray[2][2][2]=1; 
//	  for(int z=0;z<NZ;z++)  
//    {  
//        for (int x=0;x<NX;x++)  
//        {  
//            for (int y=0;y<NY;y++)  
//            {  
//                cout<<grdarray[z][x][y];  
//                
//            }  
//			cout<<endl; 
//        }  
//		cout<<endl; 
//    }  
//	  int imgSize[3]={3,3,3};
//	  std::vector<std::vector<std::vector<int> > > B(NZ,vector<vector<int> >(NX,vector<int>(NY,0)));  
//	  int numofLabs=0;
//	  BWLABEL3_test(imgSize,grdarray,B,numofLabs);
//	  
//	return 1;
//
//---solve equs of no positive
/*
Matrix4f MA;
//MA<<1,1,1,1,
//0,2,1,1,
//0,0,1,1,
//1,1,1,1;
MA<<1,1,1,1,
	0,1,0.5,0.5,
	1,1,1,1,
	1,1,1,1;
cout<<MA<<endl;
Matrix<float,4,1>Mabc;
Mabc<<0,0,0,0;
int nxindex[4]={-1,-1,-1,-1};
int nXnum=4;
ETrans4_rref(MA);
cout<<MA<<endl;
for(int i=0;i<4;i++)
{
	float fxxx=fabsf(MA(i,i));
	if(fabsf(MA(i,i))==0)
	{
		Mabc(i,0)=1;
		nxindex[i]=1;
		nXnum--;
	}
}
int m=MA.rows();
int n=MA.cols();
int np=n-nXnum;

double *K;
K=new double[n*np];


vector<int> vp;
vector<int> vnp;
for(int i=0;i<n;i++)
{
	int nc=nxindex[i];
	if(nc<0)
	{
		vp.push_back(i);
	}
	else
	{
		vnp.push_back(i);
	}

}
Solve_Eqs_Null(K,MA,n,nXnum,np,vp,vnp);
// cout<<"xxx"<<endl;
//for(int i=0;i<n;i++)
//{
//	for(int j=0;j<np;j++)
//	{
//		float xxx=K[np*i+j];
//		cout<<K[np*i+j]<<" ";
//	}
//	cout<<endl;
//}
//---solve equs of no positive
*/
}

 