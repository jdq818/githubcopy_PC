
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

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"
#include "miiMinHeap.h"

#define TAB_CHAR	9
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
void WriteLenTxt(char* chFileName,float Od,float Nd)
{
    ofstream WriteFileTxt(chFileName);
     WriteFileTxt <<right<<fixed<<setfill('0')<<setprecision(4) << Od << " " << Nd<<'\n'; 
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
	float BackTrackDistmm=0.0;
	float fdensityofminmalpath=0.0;
	float CERTAINLENGTH=2.0;
	int BackPintPosi=0;
	for (int j=1;j<vSgmtMinPathWorld.size();j++)
	{
		fSgmBPointWorld=vSgmtMinPathWorld[j];
		fSgmFPointWorld=vSgmtMinPathWorld[j-1];
		BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
	}
	fdensityofminmalpath=BackTrackDistmm/(vSgmtMinPathWorld.size()-1+0.001);
	BackTrackDistmm=0;
	BackPintPosi=(int)(CERTAINLENGTH/fdensityofminmalpath+0.5);
	int i=0;
	for (int j=0;(j+BackPintPosi)<vSgmtMinPathWorld.size();)
	{

		fSgmBPointWorld=vSgmtMinPathWorld[j+BackPintPosi];
		fSgmFPointWorld=vSgmtMinPathWorld[j];
		BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
		j=j+BackPintPosi;
	}
	return BackTrackDistmm;
}



float fCalcMinPathLengthP(vector<PointCordTypeDef> vSgmtMinPathWorld){
	
	float BackTrackDistmm=0;

	for (int j=1;j<vSgmtMinPathWorld.size();j++)
	{
		float fB[3]={vSgmtMinPathWorld[j].x,vSgmtMinPathWorld[j].y,vSgmtMinPathWorld[j].z};
	float fF[3]={vSgmtMinPathWorld[j-1].x,vSgmtMinPathWorld[j-1].y,vSgmtMinPathWorld[j-1].z};
	BackTrackDistmm=BackTrackDistmm+zxh::VectorOP_Distance(fB,fF,3);
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
bool ShortenForwardMCline(vector<miiCNode<double,float>>vPointWorld,vector<PointCordTypeDef>&vShortenPointWorld,float fPers)
{
	float BackTrackDistmm=fCalcMinPathLength(vPointWorld);
	if (BackTrackDistmm==0) 
	{
		cout<<"The length of the line is zero"<<endl;
		return false;
	}
	float meandis_of_coronaryvessel=BackTrackDistmm/(vPointWorld.size() - 1);
	float fLen=fPers*BackTrackDistmm;
	int nNumber=fLen/meandis_of_coronaryvessel;
	cout<<nNumber<<endl;
	cout<<vPointWorld.size() - 1<<endl;
	for (int j=0;j<nNumber;j++)
	{
	PointCordTypeDef XPointWorld;
	XPointWorld.x=vPointWorld[j].x;
	XPointWorld.y=vPointWorld[j].y;
	XPointWorld.z=vPointWorld[j].z;
	vShortenPointWorld.push_back(XPointWorld);
	}
	return true;
}
bool ShortenForwardMClineP(vector<PointCordTypeDef> &vPointWorld,vector<PointCordTypeDef>&vShortenPointWorld,float fPers)
{
	float BackTrackDistmm=fCalcMinPathLengthP(vPointWorld);
	if (BackTrackDistmm==0) 
	{
		cout<<"The length of the line is zero"<<endl;
		return false;
	}
	float meandis_of_coronaryvessel=BackTrackDistmm/(vPointWorld.size() - 1);
	float fLen=fPers*BackTrackDistmm;
	int nNumber=fLen/meandis_of_coronaryvessel;
	cout<<nNumber<<endl;
	cout<<vPointWorld.size() - 1<<endl;
	int Cur=0;
	float fCurLen=0;
	while (fCurLen<fLen)
	{
	PointCordTypeDef XPointWorld;
	XPointWorld.x=vPointWorld[Cur].x;
	XPointWorld.y=vPointWorld[Cur].y;
	XPointWorld.z=vPointWorld[Cur].z;
	vShortenPointWorld.push_back(XPointWorld);
		Cur++;
	    fCurLen=0;
		for(int i=0;i<Cur;i++)
		{
			float fB[3]={vPointWorld[i].x,vPointWorld[i].y,vPointWorld[i].z};
			float fF[3]={vPointWorld[i+1].x,vPointWorld[i+1].y,vPointWorld[i+1].z};
			fCurLen=fCurLen+zxh::VectorOP_Distance(fB,fF,3);
		}
	
	}
	return true;
}
int main(int argc, char *argv[])

{

	if( argc < 5 )
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeDMPShorten	vtkFile(*)(.vtk)	result-path" << endl;
		return -1;
	}
	char *MCLinFileName = argv[1];//"F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/mod05_to_unseen01_results/meanimg05v2model/MCLine.vtk";
	char *chResultPathNameVTK=argv[2];//"F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/mod05_to_unseen00_results/meanimg05v2model_lengthen/MCLine_EXT";
	char *chResultPathNameTXT=argv[3];//"F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/mod05_to_unseen00_results/meanimg05v2model_lengthen/MCLine_EXT";
	char *cheResultPathLenNameTxt=argv[4];
	float fPers=atof(argv[5]);
	//char *MCLinFileName ="F:/Coronary_0/testingdataZXHCAEDMP/HighResoResults/unseen26/vessel3/MCLineF_SMO.vtk";
	//char *chResultPathNameVTK="F:/Coronary_0/testingdataZXHCAEDMP/HighResoResults/unseen26/vessel3/MCLineF_SMO_EXT.vtk";
	//char *chResultPathNameTXT="F:/Coronary_0/testingdataZXHCAEDMP/HighResoResults/unseen26/vessel3/MCLineF_SMO_EXT.txt";

	//read the results as MCLine.vtk from different model
	vector<PointCordTypeDef>vPathPointsWorld,vShortenPointWorld;
		if (vPathPointsWorld.size()!=0)vPathPointsWorld.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		DetectFile.open(MCLinFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(MCLinFileName, vPathPointsWorld);
		}
		int ne=vPathPointsWorld.size()-1;
		int n4e=ne/4;
		float fdist=0;
		for(int i=0;i<4;i++)
		{
			int numB=n4e*i;
			int numF=n4e*(i+1);
		float fs[3]={vPathPointsWorld[numB].x,vPathPointsWorld[numB].y,vPathPointsWorld[numB].z};
		float fe[3]={vPathPointsWorld[numF].x,vPathPointsWorld[numF].y,vPathPointsWorld[numF].z};
		 fdist=fdist+zxh::VectorOP_Distance(fs,fe,3);
		}
		DetectFile.close();
	//vector<miiCNode<double,float>>vPointWorld;
	//miiCNode<double,float> ptemp;
	//for (int i=0;i<vPathPointsWorld.size();i++)
	//{
	//	ptemp.x=vPathPointsWorld[i].x;
	//	ptemp.y=vPathPointsWorld[i].y;
	//	ptemp.z=vPathPointsWorld[i].z;
	//	vPointWorld.push_back(ptemp);
	//}
		ShortenForwardMClineP(vPathPointsWorld,vShortenPointWorld,fPers);

		float foldlength=fCalcMinPathLengthP(vPathPointsWorld);
		float fnewlength=fCalcMinPathLengthP(vShortenPointWorld);
		// save lenthen-path
		cout<<"Old length: "<<foldlength<<"mm"<<endl;
		cout<<"New length: "<<fnewlength<<"mm"<<endl;
		//cout<<"Old size: "<<vPathPointsWorld.size()<<"mm"<<endl;
		cout<<"Shorten: "<<1-fnewlength/foldlength<<endl;
		WriteVtk(vShortenPointWorld, chResultPathNameVTK);
		WriteTxt(vShortenPointWorld, chResultPathNameTXT);
		WriteLenTxt(cheResultPathLenNameTxt,foldlength,fnewlength);
		cout << "Shorten Down!" << endl;

	return 1;
}