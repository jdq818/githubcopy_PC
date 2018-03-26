
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
     WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4)<< -fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] << "\n ";

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
float CalcMinDistFromPoinToLine4(vector<PointCordTypeDef>vPathPoints,vector< miiCNode<double, float> > vREFfourPoints)
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
float CalcMinDistFromPoinToLine3(vector<PointCordTypeDef>vPathPoints,vector< miiCNode<double, float> > vREFfourPoints)
{
	float distmm=0;
	float sumdistmm=0;
	float fPointWorld[3]={0,0,0};
	float fPathPointWorld[3]={0,0,0};
	
	for(int i=0;i<vREFfourPoints.size()-1;i++)
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
float CalcMinDistFromPoinToLine2(vector<PointCordTypeDef>vPathPoints,vector< miiCNode<double, float> > vREFfourPoints)
{
	float distmm=0;
	float sumdistmm=0;
	float fPointWorld[3]={0,0,0};
	float fPathPointWorld[3]={0,0,0};
	
	for(int i=0;i<vREFfourPoints.size()-2;i++)
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
float CalcMinDistFromPoinToLine1(vector<PointCordTypeDef>vPathPoints,vector< miiCNode<double, float> > vREFfourPoints)
{
	float distmm=0;
	float sumdistmm=0;
	float fPointWorld[3]={0,0,0};
	float fPathPointWorld[3]={0,0,0};
	
	for(int i=0;i<vREFfourPoints.size()-3;i++)
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
int main(int argc, char *argv[])

{
	//**.exe**
	/*if( argc < 2 )
	{
		cerr << "Usage: " << endl;
		cerr << "miiFindMaxPath	vtkFile(*)(.vtk)	FileNumber" << endl;
		return -1;
	}*/

	// get file path and header
	//char *chFileNameHrd = argv[1];

	// get out file name
	//char nFileNum = atoi(argv[2]);

	//**.exe**
/*
	 get file path and header
	char *chFileNameHrd = "E:/UI_Projects/lib/coronary/training/results/Exp2/dataset00/model01/vessel0/ModelLineCrct";

	 get out file name
	char nFileNum = 8;
*/
	//**.exe**
	if( argc < 3 )
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeDMPSelctBestBr	path0 path1 path2 path3 path4 path5 path6 refpath resultpath" << endl;
		return -1;
	}


	//**for debug
	string strFileName=string(argv[1]);
	char *REFRFileNameHrd =argv[2];//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset07/vessel1/point";
	char *chResultPath =argv[3];//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/unseen07_results_new/vessel1";


	//for debug**
	int iREFRFileLen = strlen(REFRFileNameHrd) + strlen("S.txt");
	int MCLineMinNUM4,ModelFileNUM4,MCLineMinNUM3,ModelFileNUM3,MCLineMinNUM2,ModelFileNUM2,MCLineMinNUM1,ModelFileNUM1;
	// define variables 
	char *chLetter[4]={"S.txt","A.txt","B.txt","E.txt"};
	////read the four points selected from reference (pointS,pointA,pointB,pointE)

		//read the four points selected from reference (pointS,pointA,pointB,pointE)
	char *REFRFileName;
	vector< miiCNode<double, float> > vREFfourPoints;
	for (int i=0;i<4;i++)
	{	
		REFRFileName = (char*)malloc(iREFRFileLen);
		strcpy(REFRFileName, REFRFileNameHrd);
		strcat(REFRFileName, chLetter[i]);
		ifstream str(REFRFileName);
		ReadPointTxt(REFRFileName, vREFfourPoints);
	}	
	//read the results as MCLine.vtk from different model
	vector<PointCordTypeDef>vPathPoints;
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
	float fMin44Distmm=100000;
	float fMin33Distmm=100000;
	float fMin22Distmm=100000;
	float fMin11Distmm=100000;
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;
	if (vPathPoints.size()!=0)vPathPoints.clear();
	fstream DetectFile;
	char *chFileName=new char[1024];
	string strtempstr;
	stringstream sstr;
	strtempstr=strFileName;
	sstr.clear();
	sstr<<strtempstr;
	sstr>>chFileName;
	DetectFile.open(chFileName,ios::in);
	if(DetectFile)
	{
		ReadVtk(chFileName, vPathPoints);
		float Min4DistSum=CalcMinDistFromPoinToLine4(vPathPoints,vREFfourPoints);
		float Min3DistSum=CalcMinDistFromPoinToLine3(vPathPoints,vREFfourPoints);
		float Min2DistSum=CalcMinDistFromPoinToLine2(vPathPoints,vREFfourPoints);
		float Min1DistSum=CalcMinDistFromPoinToLine1(vPathPoints,vREFfourPoints);
		std::cout<<Min4DistSum<<"\n";
		std::cout<<Min3DistSum<<"\n";
		std::cout<<Min2DistSum<<"\n";
		std::cout<<Min1DistSum<<"\n";
		if(Min4DistSum<fMin4Distmm)
		{
			fMin4Distmm=Min4DistSum;

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


	fMin44Distmm=fMin4Distmm;
	fMin33Distmm=fMin3Distmm;
	fMin22Distmm=fMin2Distmm;
	fMin11Distmm=fMin1Distmm;
		string strBestNUMFilename = string(chResultPath);
		strBestNUMFilename = strBestNUMFilename + "/" + string("DistFromREFPont.txt");
		ofstream WriteFileBF(strBestNUMFilename);
		WriteFileBF <<right<<fixed<<setfill('0')<<setprecision(4) <<fMin44Distmm<<" "<<setprecision(4)<<fMin33Distmm<<" "<<setprecision(4)<<fMin22Distmm<<" "<<setprecision(4)<<fMin11Distmm<<'\n';
		WriteFileBF.close();
		cout << "Calculation Down" << endl;

	return 1;
}