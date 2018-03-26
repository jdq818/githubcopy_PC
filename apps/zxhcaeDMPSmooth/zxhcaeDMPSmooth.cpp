#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"


#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"


#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
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
	float fMean;
	float fStdVar;
	float fMaxErr;
}ErrTypeDef;

bool ReadVtk(char *chFileName, vector< miiCNode<double, float> > &vPointCord)
{
	if (chFileName == NULL)
	{
		cerr << "Cannot find VTK file!" << endl;
		return false;
	}
	//	char* chFileName = "reference.vtk";
	vtkSmartPointer<vtkUnstructuredGridReader> iVtkReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	iVtkReader->SetFileName( chFileName );
	iVtkReader->Update();

	vtkSmartPointer<vtkUnstructuredGrid> iVtkGridRead = iVtkReader->GetOutput();

	int nPointNum = iVtkGridRead->GetMaxCellSize();

	double fCord[3];
	miiCNode<double, float> strctTempPoint;

	for (int i = 0; i < nPointNum; i++)
	{
		iVtkGridRead->GetPoint(i, fCord);

		strctTempPoint.x = fCord[0];
		strctTempPoint.y = fCord[1];
		strctTempPoint.z = fCord[2];
		vPointCord.push_back(strctTempPoint);
	}
	return true;
}
bool ReadTxt(char *chFileName,vector<short>&sMolPontInts)
{
	 ifstream inMolIntdata(chFileName);
	    double data;string strLine;
   bool bsizeeque=true;
	if(!inMolIntdata)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << chFileName << endl;
		return false; //exit(1); // terminate with error
	}
	else
	{
		
	while(getline(inMolIntdata,strLine,'\n'))
        { 
	 data = atof(strLine.c_str());
     sMolPontInts.push_back(data);
         }
	return true;

   }
}

void WriteVtk(vector<miiCNode<double,float>> vdfPointCord, char* chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();
	int nPointNum = vdfPointCord.size();
	for (int i = 0; i < nPointNum; i++)
	{
		
		iPoints->InsertNextPoint(vdfPointCord[i].x, vdfPointCord[i].y, vdfPointCord[i].z);
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

void WriteVtkPD(vector< PointCordTypeDef > PointCord, char* chFileName)
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
void WriteTxt(vector<miiCNode<double,float>> vdfPointCord, char* chFileName)
{
    ofstream WriteFileTxt(chFileName);
	int nPointNum = vdfPointCord.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
	 fImgPixel[0] = vdfPointCord[i].x;
	 fImgPixel[1] = vdfPointCord[i].y;
	 fImgPixel[2] = vdfPointCord[i].z;
     WriteFileTxt <<right<<fixed<<setfill('0')<<setprecision(4) << -fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] <<'\n'; 

	}	
}
void WriteCAIntTxt(vector<short> sMolPontInts,char* chFileName)
{
ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = sMolPontInts.size();
	for (int i = 0; i < nPointNum; i++)
	{
     WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4) <<sMolPontInts[i] <<"\n";

	}	
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
	// read 1st number
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	strctTempPoint.x =atof(strNum.c_str());
	//fMean=1024*fMean/1896;//Add by JDQ
	// read 2nd number
	nStart = nEnd + 1;
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	strctTempPoint.y =atof(strNum.c_str());
	// read 3rd number
	nStart = nEnd + 1;
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	strctTempPoint.z =atof(strNum.c_str());
	cl.push_back(strctTempPoint);
    return true;
}
float CalcMinDistFromPoinToLine(vector<PointCordTypeDef>vPathPoints,vector< miiCNode<double, float> > vREFfourPoints)
{
	float distmm=0;
	float sumdistmm=0;
	float fPointWorld[3]={0,0,0};
	float fPathPointWorld[3]={0,0,0};
	
	for(int i=0;i<vREFfourPoints.size();i++)
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
		sumdistmm=sumdistmm+mindist;
	}
	return sumdistmm;
}
float CalcDistmm2(miiCNode<double, float> nPixPoint1,miiCNode<double, float> nPixPoint2)
{  //wrong, delete this function zxh noted
	float Distmm2 = (nPixPoint1.x - nPixPoint2.x)*(nPixPoint1.x -nPixPoint2.x) + 
			(nPixPoint1.y -nPixPoint2.y)*(nPixPoint1.y - nPixPoint2.y) + 
			(nPixPoint1.z - nPixPoint2.z)*(nPixPoint1.z - nPixPoint2.z);
	return Distmm2;
}
float fCalcMinPathLength(vector<miiCNode<double,float>> vSgmtMinPathWorld){
	miiCNode<double,float>fSgmFPointWorld,fSgmBPointWorld;
	float BackTrackDistmm=0;

	for (int j=1;j<vSgmtMinPathWorld.size();j++)
	{
	fSgmBPointWorld=vSgmtMinPathWorld[j];
	fSgmFPointWorld=vSgmtMinPathWorld[j-1];
	BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
	}

return BackTrackDistmm;
}
bool LengthenForwardMCline(vector<miiCNode<double,float>>vPointWorld,vector<PointCordTypeDef>&vLengthenPointWorld)
{
	float BackTrackDistmm=fCalcMinPathLength(vPointWorld);
	float meandis_of_coronaryvessel=BackTrackDistmm/(vPointWorld.size() - 1);
	float extendLENGTHmm=5;
	float CalCLENGTHmm=3;
	int xtendNUM=CalCLENGTHmm/meandis_of_coronaryvessel;
	float BPointWorld[3]={0,0,0};
	float SPointWorld[3]={0,0,0};
	float tempVec[3]={0,0,0};
	PointCordTypeDef XPointWorld;
	BPointWorld[0]=vPointWorld[xtendNUM].x;
	BPointWorld[1]=vPointWorld[xtendNUM].y;
	BPointWorld[2]=vPointWorld[xtendNUM].z;
	SPointWorld[0]=vPointWorld[0].x;
	SPointWorld[1]=vPointWorld[0].y;
	SPointWorld[2]=vPointWorld[0].z;
	tempVec[0]=SPointWorld[0]-BPointWorld[0];
	tempVec[1]=SPointWorld[1]-BPointWorld[1];
	tempVec[2]=SPointWorld[2]-BPointWorld[2];
	zxh::VectorOP_Normalise(tempVec,3);
	XPointWorld.x=SPointWorld[0]+extendLENGTHmm*tempVec[0];
	XPointWorld.y=SPointWorld[1]+extendLENGTHmm*tempVec[1];
	XPointWorld.z=SPointWorld[2]+extendLENGTHmm*tempVec[2];
	vLengthenPointWorld.push_back(XPointWorld);
	for (int j=0;j<vPointWorld.size();j++)
	{
	XPointWorld.x=vPointWorld[j].x;
	XPointWorld.y=vPointWorld[j].y;
	XPointWorld.z=vPointWorld[j].z;
	vLengthenPointWorld.push_back(XPointWorld);
	}
	return true;
}
void CorrectImagePos(miiCNode<double, int> *diPoint,zxhImageDataT<short> imgReadRaw)
{
	int nImWX, nImWY, nImWZ, nImWT;
	imgReadRaw.GetImageSize(nImWX, nImWY, nImWZ, nImWT);
	if (diPoint->x >= nImWX)
		diPoint->x = nImWX - 1;
	if (diPoint->y >= nImWY)
		diPoint->y = nImWY - 1;
	if (diPoint->z >= nImWZ)
		diPoint->z = nImWZ - 1;
	if (diPoint->x <= 0)
		diPoint->x = 0;
	if (diPoint->y <= 0)
		diPoint->y = 0;
	if (diPoint->z <= 0)
		diPoint->z = 0;
};
miiCNode<double, float> ImageTransToWorldPoint(miiCNode<double, int>diPoint,zxhImageDataT<short> imgReadRaw )
{
	const zxhImageInfo *m_pBaseImgInfo;
	m_pBaseImgInfo=imgReadRaw.GetImageInfo();
	float fCord[3];
	CorrectImagePos(&diPoint,imgReadRaw);
	fCord[0] = (float)(diPoint.x);
	fCord[1] = (float)(diPoint.y); 
	fCord[2] = (float)(diPoint.z);
	m_pBaseImgInfo->ImageToWorld(fCord);
	miiCNode<double, float> dfPoint;
	dfPoint.x= fCord[0];
	dfPoint.y = fCord[1];
	dfPoint.z= fCord[2];
	dfPoint.val=diPoint.val;
	return dfPoint;

};
miiCNode<double, int> WorldTransToImagePoint(miiCNode<double, float>dfPointWorld,zxhImageDataT<short> imgReadRaw )
{
	const zxhImageInfo *m_pBaseImgInfo;
	m_pBaseImgInfo=imgReadRaw.GetImageInfo();
	float fCord[3];
	fCord[0] = dfPointWorld.x;
	fCord[1] = dfPointWorld.y; 
	fCord[2] = dfPointWorld.z;
	m_pBaseImgInfo->WorldToImage(fCord);

	miiCNode<double, int> diPoint;
	diPoint.x= int(fCord[0]);
	diPoint.y =int(fCord[1]);
	diPoint.z= int(fCord[2]);
	diPoint.val=dfPointWorld.val;
	CorrectImagePos(&diPoint,imgReadRaw);

	return diPoint;

}

bool SMoothPath(vector<miiCNode<double,float>> vPointWorld,vector<miiCNode<double,float>> &vSMOPointWorld)
{
	float BackTrackDistmm=fCalcMinPathLength(vPointWorld);//calculate the total length
	float m_meandis_of_coronaryvessel=BackTrackDistmm/(vPointWorld.size() - 1);

	float fStepDist=BackTrackDistmm/100;
	
	if (vSMOPointWorld.size() > 0)
	{
		vSMOPointWorld.clear();
	}

	vSMOPointWorld.push_back(vPointWorld[0]);

	int idistNum_2mm=2/m_meandis_of_coronaryvessel;
	for(int i=1;i<101;i++)//100 seg point
	{
		int SegNUM=0;
		float dist=0;
		for(int j=0;j<vPointWorld.size();j++)
		{
			dist=dist+sqrt(CalcDistmm2(vPointWorld[j],vPointWorld[j+1]));
			if(dist>=i*fStepDist)
			{
				SegNUM=j;
				break;
			}
		}

		int SegStart=SegNUM-idistNum_2mm;
		int SegEnd=SegNUM+idistNum_2mm;
		if (SegStart<=0)SegStart=0;
		if (SegEnd>=vPointWorld.size()-1)SegEnd=vPointWorld.size()-1;
		float fSegSum[3]={0};
		int n=0;
		for (int k=SegStart;k<=SegEnd;k++)
		{
			fSegSum[0]=fSegSum[0]+vPointWorld[k].x;
			fSegSum[1]=fSegSum[1]+vPointWorld[k].y;
			fSegSum[2]=fSegSum[2]+vPointWorld[k].z;
			n++;
		}
		miiCNode<double,float> mtempPointWorld;
		miiCNode<double,int> mtempPoint;
		mtempPointWorld.x=fSegSum[0]/n;
		mtempPointWorld.y=fSegSum[1]/n;
		mtempPointWorld.z=fSegSum[2]/n;
		vSMOPointWorld.push_back(mtempPointWorld);
		/*float nCord[3];
		float fCord[3];
		nCord[0] =mtempPointWorld.x; 
		nCord[1] =mtempPointWorld.y; 
		nCord[2] =mtempPointWorld.z; 
		m_pBaseImgInfo->WorldToImage(nCord);	
		mtempPoint.x = int(nCord[0]+0.5);
		mtempPoint.y = int(nCord[1]+0.5);
		mtempPoint.z =int(nCord[2]+0.5);
		CorrectImagePos(&mtempPoint,imgReadRaw);*/
		
	}
	return true;
}
int main(int argc, char *argv[])

{
	
	if( argc < 3 )
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeDMPSmooth	vtkFile(*)(.vtk)	rawimag  resultvtk resulttxt" << endl;
		return -1;
	}
	char *MCLinFileName =argv[1];//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResultsSM/unseen00_resultsADE/vessel0/MCLine.vtk";//"F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/mod05_to_unseen01_results/meanimg05v2model/MCLine.vtk";
	char *chResultPathNameVTK=argv[2];//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResultsSM/unseen00_resultsADE/vessel0/MCLine_SMO.vtk";//"F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/mod05_to_unseen00_results/meanimg05v2model_lengthen/MCLine_EXT";
	char *chResultPathNameTXT=argv[3];//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResultsSM/unseen00_resultsADE/vessel0/MCLine_SMO.txt";//"F:/Coronary_0/t

	//char *MCLinFileName ="F:/Coronary_0/ZXHJDQCAEDMP/DHResDFM712_VVS/mod01_to_unseen00_results/meanimg01v0model/SMLLine.vtk";//"F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/mod05_to_unseen01_results/meanimg05v2model/MCLine.vtk";
	//char *chResultPathNameVTK="F:/Coronary_0/ZXHJDQCAEDMP/DHResDFM712_VVS/mod01_to_unseen00_results/meanimg01v0model/MCLineSML_SMO.vtk";//"F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/mod05_to_unseen00_results/meanimg05v2model_lengthen/MCLine_EXT";
	//char *chResultPathNameTXT="F:/Coronary_0/ZXHJDQCAEDMP/DHResDFM712_VVS/mod01_to_unseen00_results/meanimg01v0model/MCLineSML_SMO.txt";//"F:/Coronary_0/t


	//read the results as MCLine.vtk from different model
	vector<miiCNode<double,float>> vPathPointsWorld;
	if (vPathPointsWorld.size()!=0)vPathPointsWorld.clear();
	fstream DetectFile;
	/*char chIdx[25];
	_itoa_s(i, chIdx, 10);*/
	DetectFile.open(MCLinFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(MCLinFileName, vPathPointsWorld);
	}
	DetectFile.close();
	vector<miiCNode<double,float>>vSMOPointWorld;
	SMoothPath(vPathPointsWorld,vSMOPointWorld);

	vector<PointCordTypeDef>vPointWorld;
	PointCordTypeDef ptemp;
	for (int i=0;i<vSMOPointWorld.size();i++)
	{
		ptemp.x=vSMOPointWorld[i].x;
		ptemp.y=vSMOPointWorld[i].y;
		ptemp.z=vSMOPointWorld[i].z;
		vPointWorld.push_back(ptemp);
	}
	// save smoothen-path

	WriteTxt(vSMOPointWorld, chResultPathNameTXT);
	int nFileLen = strlen(chResultPathNameVTK)+1 ;
	char *chResultNameVTK;
	chResultNameVTK= (char*)malloc(nFileLen);
	strcpy(chResultNameVTK, chResultPathNameVTK);

	WriteVtkPD(vPointWorld, chResultNameVTK);
	free(chResultNameVTK);
	cout << "Smoothing Down!" << endl;


}