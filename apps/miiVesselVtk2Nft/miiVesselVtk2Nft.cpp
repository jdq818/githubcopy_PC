
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCamera.h"
#include "vtkPolyVertex.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkCellLinks.h"
#include "vtkAssembly.h"

#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"

#include "vtkLookupTable.h"
#include "vtkColorTransferFunction.h"

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

#define SPACE_CHAR	32

typedef struct
{
	float x;
	float y;
	float z;
}PointCordTypeDef;

void readVtk(char * chFileName, vector<PointCordTypeDef> &PointCord)
{
//	char* chFileName = "reference.vtk";
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

bool World2Image(const zxhImageInfo *pImageInfo, vector<PointCordTypeDef> &WorldPoints, \
	vector<PointCordTypeDef> &ImPoints)
{
	if (WorldPoints.size() == 0)
	{
		cerr << "The vector is null!" << endl;
		return false;
	}

	float fTempCord[3];
	PointCordTypeDef strctTempPoint;

	int nPointNum = WorldPoints.size();

	// convert coordinate
	for (int i = 0; i < nPointNum; i++)
	{
		fTempCord[0] = WorldPoints[i].x;
		fTempCord[1] = WorldPoints[i].y;
		fTempCord[2] = WorldPoints[i].z;

		pImageInfo->WorldToImage(fTempCord);
//		pImageInfo->ImageToPhysical(fTempCord);

		strctTempPoint.x = fTempCord[0] + 0.5;
		strctTempPoint.y = fTempCord[1] + 0.5;
		strctTempPoint.z = fTempCord[2] + 0.5;

		ImPoints.push_back(strctTempPoint);
	}

	return true;
}

bool SaveImage(const zxhImageInfo * pImageInfo, short *sImgData, string strFileName)
{
	if (pImageInfo == NULL || sImgData == NULL)
	{
		std::cerr << "The pointer of image is null!\n";
		return false;
	}

	//define the variables
	short m_pixel;

	// 1 new image 
	zxhImageDataT<short> imgTest1; 

	imgTest1.NewImage( pImageInfo->Dimension, pImageInfo->Size, pImageInfo->Spacing, pImageInfo ) ; 
	const int *image_size = pImageInfo->Size ;
	for( int i = 0; i < image_size[2]; i++ )
		for( int j = 0; j < image_size[1]; j++ )
			for( int k = 0; k < image_size[0]; k++ )
			{
				//get a pixel's value
				m_pixel = sImgData[i * image_size[1] * image_size[0] + j * image_size[0] + k];

				imgTest1.SetPixelByGreyscale(k, j, i, 0, m_pixel);
			}

			zxh::SaveImage(&imgTest1, strFileName); 

			return true;
}

int main(int argc, char *argv[])
{
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << "miiVesselVtk2Nft InputNiftiFileName   InputVtkFileName	OutputNiftiFileName	LabelStart&EndPoints" << std::endl;
		return -1;
	}

	// get "nifti" file name
	string strInNiiFileName = string(argv[1]);

	// get "vtk" file name
	char *chReadVtkFileName = argv[2];

	string strOutNiiFileName = string(argv[3]);

	int bStartEndFlg = 1;
	if (argc == 5)
	{
		bStartEndFlg = atoi(argv[4]);
	}

	// read a "nifti" file
	zxhImageDataT<short> imgRead; 
	if( zxh::OpenImage( &imgRead, strInNiiFileName ) == false )
	{
		std::cerr << "nifti-file is not found!"; 
		return -1;
	}

	fstream DetectFile;
	DetectFile.open(argv[1],ios::in);
	if(!DetectFile)
	{
		cout << "Cannot find nifti file!" << endl;
		DetectFile.close();
		return -1;
	}
	DetectFile.close();

	DetectFile.open(argv[2],ios::in);
	if(!DetectFile)
	{
		cout << "Cannot find vtk file!" << endl;
		DetectFile.close();
		return -1;
	}
	DetectFile.close();

	// define variables 
	vector<PointCordTypeDef> WorldPoints, ImgPoints;

	// read "*.vtk" file
	readVtk(chReadVtkFileName, WorldPoints);

	// convert coordinate
	World2Image(imgRead.GetImageInfo(), WorldPoints, ImgPoints);

	//get the image size
	int nImWX, nImWY, nImWZ, nImWT;
	imgRead.GetImageSize(nImWX, nImWY, nImWZ, nImWT);

	short *sImgData = new short[nImWX * nImWY * nImWZ];

	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	{
		sImgData[i] = 0;
	}

	for (int i = 0; i < ImgPoints.size(); i++)
	{
		sImgData[(int)ImgPoints[i].z * nImWY * nImWX + (int)ImgPoints[i].y * nImWX + (int)ImgPoints[i].x] = 1000;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// label the starting and end points
	if (bStartEndFlg)
	{
		for (int i = 0; i < 6; i++)
		{
			int nx_1 = (int)ImgPoints[0].x + gNbr[i][0];
			int ny_1 = (int)ImgPoints[0].y + gNbr[i][1];
			int nz_1 = (int)ImgPoints[0].z + gNbr[i][2];

			int nx_2 = (int)ImgPoints[ImgPoints.size()-1].x + gNbr[i][0];
			int ny_2 = (int)ImgPoints[ImgPoints.size()-1].y + gNbr[i][1];
			int nz_2 = (int)ImgPoints[ImgPoints.size()-1].z + gNbr[i][2];

			sImgData[nz_1 * nImWY * nImWX + ny_1 * nImWX + nx_1] = 2000;

			sImgData[nz_2 * nImWY * nImWX + ny_2 * nImWX + nx_2] = 2000;
		}
	}

	// save "*.nii" file
	SaveImage(imgRead.GetImageInfo(), sImgData, strOutNiiFileName);

	delete[] sImgData;
	return 1;
}