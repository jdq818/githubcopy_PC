

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


#include <algorithm>
#include "jdq2017util.h"
#include "jdqdijkstra.h"
#include "jdqPoint.h"

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
 const vector<PointCordTypeDef> & vOriPoints, // ���ﲻ��ָ������ã��ᵼ��ԭ��������vector��copyһ�ݵ������������
 int iNewSize,vector<PointCordTypeDef> &vNewPoints)
{
    int iSizeOri=vOriPoints.size();
    if( iNewSize<2 || iSizeOri<2 )
        return false ; 
    vNewPoints.empty(); //��֤���vector�ǿյ�
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
    return true ; //���ݷ���ֵ���ж������������Ƿ���ȷ��ִ����
}
void WriteCAIntTxt(char *chFileName,vector<short> sMolPontInts)
{
	ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = sMolPontInts.size();
	for (int i = 0; i < nPointNum; i++)
	{
		WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4) <<sMolPontInts[i] <<'\n';
	}
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
	string strfilenameraw =  "K:/JDQ/CCTA_CAR/RCAA_32/training/dataset00/CAE_ME_L.nii.gz";
	char *chRefCurvefilename ="K:/JDQ/CCTA_CAR/RCAA_32/training/dataset00/CL0.vtk";
	char *chTarCurvefilename="K:/JDQ/CCTA_CAR/RCAA_32/training/dataset00/CL1.vtk";
	float flength=0,flengthdeform=0,fsectionlength=0;
	float ResampleNUM=9999;
	vector<PointCordTypeDef> vRefCurPontsWorld,vTarCurPontsWorld,vref,vc1;
	std::vector<jdq2017::point3D> cl;
	// read reference curve 
	fstream DetectFile1;
	DetectFile1.open(chRefCurvefilename,ios::in);
	if(DetectFile1)
		{
			ReadVtk(chRefCurvefilename, vRefCurPontsWorld);
		}
	DetectFile1.close();
	//Get the length of the whole curve
	float fWhLengh;
	GetLenOfCurve(vRefCurPontsWorld,fWhLengh);
     double samlingDistance=fWhLengh/1000;
	//resample the reference curve
	vref.clear();
	ResampleCurve(vRefCurPontsWorld,samlingDistance,vref);

// read target curve 
	DetectFile1.open(chTarCurvefilename,ios::in);
	if(DetectFile1)
		{
			ReadVtk(chTarCurvefilename, vTarCurPontsWorld);
		}
	DetectFile1.close();
	//resample the target curve
	vc1.clear();
	ResampleCurve(vTarCurPontsWorld,samlingDistance,vc1);
	//corespondence between vref and vc1
	//correspondence();
	return 1;
}

 