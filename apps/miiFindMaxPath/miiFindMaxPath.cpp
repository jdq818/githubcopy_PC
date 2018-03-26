
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
void WriteCAIntTxt(vector<short> sMolPontInts,char* chFileName)
{
ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = sMolPontInts.size();
	for (int i = 0; i < nPointNum; i++)
	{
     WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4) <<sMolPontInts[i] <<"\n";

	}	
}

int main(int argc, char *argv[])

{
	//**.exe**
	if( argc < 2 )
	{
		cerr << "Usage: " << endl;
		cerr << "miiFindMaxPath	vtkFile(*)(.vtk)	FileNumber" << endl;
		return -1;
	}

	// get file path and header
	char *chFileNameHrd = argv[1];

	// get out file name
	char nFileNum = atoi(argv[2]);

	//**.exe**

	 //get file path and header
	//char *chFileNameHrd = "F:/Coronary_0/trainningdataZXHCAEDMP/DirHighResultsDFM/mod01_to_unseen00_results/meanimg01v0model";

	 //get out file name
	//char nFileNum = 30;

	/*char *chFileNameHrd="F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu0/MCLine";
	char nFileNum = 16;*/
	// define variables 
	vector<PointCordTypeDef> vPathPoints;
	vector<short> sMolPontInts;
	int nFileLen = strlen(chFileNameHrd) + strlen("xx.vtk") + 1;
	int nFileLen1 = strlen(chFileNameHrd) + strlen("/_intensity")+strlen("xx.txt") + 1;
	char *chFileName;
	fstream DetectFile;
	int nMaxPointNum = 0, nMaxPointIdx = -1;
	for (int i = 0; i < nFileNum; i++)
	{
		chFileName = (char*)malloc(nFileLen);
		char chIdx[25];
		_itoa_s(i, chIdx, 10);
		strcpy(chFileName, chFileNameHrd);
		strcat(chFileName, chIdx);
		strcat(chFileName, ".vtk");

		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPoints);
			if (vPathPoints.size() > nMaxPointNum)
			{
				nMaxPointNum = vPathPoints.size();
				nMaxPointIdx = i;
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

	// read max-path
	chFileName = (char*)malloc(nFileLen);
	char chIdx[25];
	_itoa_s(nMaxPointIdx, chIdx, 10);
	strcpy(chFileName, chFileNameHrd);
	strcat(chFileName, chIdx);
	strcat(chFileName, ".vtk");
	ReadVtk(chFileName, vPathPoints);
	delete[] chFileName;
	chFileName = (char*)malloc(nFileLen1);
	strcpy(chFileName, chFileNameHrd);
	strcat(chFileName,"_intensity");
	strcat(chFileName, chIdx);
	strcat(chFileName, ".txt");
    ReadTxt(chFileName,sMolPontInts);
    delete[] chFileName;
	// save max-path
	chFileName = (char*)malloc(nFileLen);
	strcpy(chFileName, chFileNameHrd);
	strcat(chFileName, ".vtk");
	WriteVtk(vPathPoints, chFileName);
	delete[] chFileName;

	chFileName= (char*)malloc(nFileLen);
	strcpy(chFileName, chFileNameHrd);
	strcat(chFileName, ".txt");
	WriteTxt(vPathPoints, chFileName);
	delete[] chFileName;

	
	
	int nLen2 = strlen(chFileNameHrd) + strlen("_intensity")+ strlen(".txt")+1;
	chFileName= (char*)malloc(nLen2);
	strcpy(chFileName, chFileNameHrd);
	strcat(chFileName, "_intensity");
	strcat(chFileName, ".txt");
	WriteCAIntTxt(sMolPontInts,chFileName);
	delete[] chFileName;
	cout << "Found max-path!" << endl;

	return 1;
}