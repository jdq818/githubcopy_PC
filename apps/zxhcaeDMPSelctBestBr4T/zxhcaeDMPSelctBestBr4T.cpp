
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
	if( argc < 10 )
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeDMPSelctBestBr	path0 path1 path2 path3 path4 path5 path6 refpath resultpath" << endl;
		return -1;
	}


	//**for debug
	string *strFileNameHrd;
	strFileNameHrd=new string[8];
	strFileNameHrd[0]=string(argv[1]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod00_to_unseen07_new_results/meanimg00v1model";
	strFileNameHrd[1]=string(argv[2]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod01_to_unseen07_new_results/meanimg01v1model";
	strFileNameHrd[2]=string(argv[3]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod02_to_unseen07_new_results/meanimg02v1model";
	strFileNameHrd[3]=string(argv[4]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod03_to_unseen07_new_results/meanimg03v1model";
	strFileNameHrd[4]=string(argv[5]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod04_to_unseen07_new_results/meanimg04v1model";
	strFileNameHrd[5]=string(argv[6]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod05_to_unseen07_new_results/meanimg05v1model";
	strFileNameHrd[6]=string(argv[7]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod06_to_unseen07_new_results/meanimg06v1model";
	strFileNameHrd[7]=string(argv[8]);//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod06_to_unseen07_new_results/meanimg06v1model";
	char *REFRFileNameHrd =argv[9];//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset07/vessel1/point";
	char *chResultPath =argv[10];//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/unseen07_results_new/vessel1";


	//string *strFileNameHrd;
	//strFileNameHrd=new string[8];
	//strFileNameHrd[0]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod00_to_unseen08_results/meanimg00v0model";
	//strFileNameHrd[1]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod01_to_unseen08_results/meanimg01v0model";
	//strFileNameHrd[2]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod02_to_unseen08_results/meanimg02v0model";
	//strFileNameHrd[3]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod03_to_unseen08_results/meanimg03v0model";
	//strFileNameHrd[4]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod04_to_unseen08_results/meanimg04v0model";
	//strFileNameHrd[5]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod05_to_unseen08_results/meanimg05v0model";
	//strFileNameHrd[6]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod06_to_unseen08_results/meanimg06v0model";
	//strFileNameHrd[7]="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/mod07_to_unseen08_results/meanimg07v0model";
	//char *REFRFileNameHrd ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/testing/dataset08/vessel0/point";
	//char *chResultPath ="F:/Coronary_0/testingdataZXHCAEDMP/LowResoResults_LR/unseen08_results/vessel0";

	//for debug**
	int iREFRFileLen = strlen(REFRFileNameHrd) + strlen("S.txt");
	int nModNum = 8;
	int nModinFileMaxNUM=30;

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
	//Define every best output file
	string strEveryBestDistFilename = string(chResultPath);
	strEveryBestDistFilename = strEveryBestDistFilename + "/" + string("BestDist.txt");
	ofstream WriteFileBD(strEveryBestDistFilename);
	
	
	//calculate the minimal dist case in whole range and for every model
	float fMin4Distmm=100000;
	float fMin3Distmm=100000;
	float fMin2Distmm=100000;
	float fMin1Distmm=100000;
	for (int i = 0; i < nModNum; i++)
	{
	float fMin44Distmm=100000;
	float fMin33Distmm=100000;
	float fMin22Distmm=100000;
	float fMin11Distmm=100000;
	int ModelFileNUM44,ModelFileNUM33,ModelFileNUM22,ModelFileNUM11;
		for (int j=0;j<nModinFileMaxNUM;j++)
		{

			if (vPathPoints.size()!=0)vPathPoints.clear();
			fstream DetectFile;
			char *chFileName=new char[1024];
			string strtempstr;
			stringstream sstr;
			string strj;
			sstr.clear();
			sstr<<j;
			sstr>>strj;
			strtempstr=strFileNameHrd[i]+"/MCLine"+strj+".vtk";
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

				//whole range
				if(Min4DistSum<fMin4Distmm)
				{
					fMin4Distmm=Min4DistSum;
					ModelFileNUM4=i;
					MCLineMinNUM4=j;

				}
				if(Min3DistSum<fMin3Distmm)
				{
					fMin3Distmm=Min3DistSum;

					ModelFileNUM3=i;
					MCLineMinNUM3=j;

				}
				if(Min2DistSum<fMin2Distmm)
				{
					fMin2Distmm=Min2DistSum;

					ModelFileNUM2=i;
					MCLineMinNUM2=j;

				}
					if(Min1DistSum<fMin1Distmm)
				{
					fMin1Distmm=Min1DistSum;

					ModelFileNUM1=i;
					MCLineMinNUM1=j;

				}

			//model file

				if(Min4DistSum<fMin44Distmm)
				{
					fMin44Distmm=Min4DistSum;
					ModelFileNUM44=j;

				}
				if(Min3DistSum<fMin33Distmm)
				{
					fMin33Distmm=Min3DistSum;
					ModelFileNUM33=j;

				}
				if(Min2DistSum<fMin22Distmm)
				{
					fMin22Distmm=Min2DistSum;

					ModelFileNUM22=j;

				}
					if(Min1DistSum<fMin11Distmm)
				{
					fMin11Distmm=Min1DistSum;

				
					ModelFileNUM11=j;

				}
			
			}
			DetectFile.close();
			delete[] chFileName;	

		}//for file
	WriteFileBD <<right<<fixed<<setfill('0')<<setprecision(0) << i << " " << ModelFileNUM44<<" "<<setprecision(4) <<fMin44Distmm<<" "<<setprecision(4)<<fMin33Distmm<<" "<<setprecision(4)<<fMin22Distmm<<" "<<setprecision(4)<<fMin11Distmm<<endl;

	}//for model
	WriteFileBD.close();
	//output the minimal dist in whole range
	string strBestNUMFilename = string(chResultPath);
	strBestNUMFilename = strBestNUMFilename + "/" + string("BestNUM.txt");
	ofstream WriteFileBF(strBestNUMFilename);
	WriteFileBF <<right<<fixed<<setfill('0')<<setprecision(0) << ModelFileNUM4<< " " << MCLineMinNUM4<<" "<<setprecision(4) <<fMin4Distmm<<" "<<setprecision(4)<<fMin3Distmm<<" "<<setprecision(4)<<fMin2Distmm<<" "<<setprecision(4)<<fMin1Distmm<<'\n';
	WriteFileBF.close();


	// read best-path
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
		strmacMCL=strFileNameHrd[ModelFileNUM4]+"/MCLine"+strj+".vtk";
		sstr1.clear();
		sstr1<<strmacMCL;
		sstr1>>chFileName;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
		}

		DetectFile.close();
		// save max-path
		char chTemp1[25];
		_itoa_s(ModelFileNUM4, chTemp1, 10);
		char chTemp2[25];
		_itoa_s(MCLineMinNUM4, chTemp2, 10);

		int nFileLen = strlen(chResultPath) + strlen("/MCLine.vtk")+1;
		chFileName = (char*)malloc(nFileLen);
		strcpy(chFileName, chResultPath);
		strcat(chFileName, "/MCLine.vtk");
		WriteVtk(vPathPoints, chFileName);
		free(chFileName);

		chFileName= (char*)malloc(nFileLen);
		strcpy(chFileName, chResultPath);
		strcat(chFileName, "/MCLine.txt");
		WriteTxt(vPathPoints, chFileName);
		free(chFileName);


		
		cout << "Found best-path!" << endl;

	return 1;
}