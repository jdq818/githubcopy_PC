

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
	vector<int>vLnum;
}PointCLTypeDef;
typedef struct
{
	float x;
	float y;
	float z;
	vector<pair<int,int>>vLPnum;
}PointLPTypeDef;
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
char *GetPathEXT(char *chFileName)
{  char path_buffer[_MAX_PATH];  
   char drive[_MAX_DRIVE];  
   char dir[_MAX_DIR];  
   char fname[_MAX_FNAME];  
   char ext[_MAX_EXT];  
  
   _splitpath( chFileName, drive, dir, fname, ext );  

   return ext;
}
bool txt2vtk(std::vector<jdq2017::point3D> &ref)
{

	for(int i=0;i<ref.size();i++)
	{
		ref[i]._x=-1*ref[i]._x;
		ref[i]._y=-1*ref[i]._y;
	}
	return true;

}
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


void Write2Vtk_jdq(std::vector<jdq2017::point3D> PointCord, char* chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	int nPointNum = PointCord.size();

	for (int i = 0; i < nPointNum; i++)
	{
		iPoints->InsertNextPoint(-1*PointCord[i]._x, -1*PointCord[i]._y, PointCord[i]._z);
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


void Write2Txt_jdq(std::vector<jdq2017::point3D> vPoints,char *chFileName)
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












int main(int argc, char *argv[])
{
	if( argc < 4 )
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeDMPModelIntenGen	RawImage(.nii)	ModelRef(.vtk) MapModelFile(.vtk) ResultPath " << endl;
		return -1;
	} 
	//
	char *chRefCurvefilename =argv[1];
	char *chTarCurvefilename=argv[2];
	char *chResamfilename=argv[3];
	
	//read and resample the curve
	//--..--..--..--..
	
	/*char *chRefCurvefilename ="E:/work_jdq/rebuttal_media18/rebuttal/Exp/LowRes_0307/dataset03/vessel1/reference.txt";
	char *chTarCurvefilename="E:/work_jdq/rebuttal_media18/rebuttal/Exp/LowRes_0307/dataset03/vessel1/SMDMCLine.vtk";
		char *chResamfilename="E:/work_jdq/rebuttal_media18/rebuttal/Exp/LowRes_0307/dataset03/vessel1/SMDMCLine_resa";*/
		

	std::vector<jdq2017::point3D> ref,cl;
	//std::vector<jdq2017::point3D> cl;
	char* chext=GetPathEXT(chRefCurvefilename);
	if(strcmp(chext,".vtk")==0)// read vtk file
	{
		if (! jdq2017::readCenterlinevtk(chRefCurvefilename, ref) )//read the reference line
		{
			if (!outputOneLine) 
			{
				std::cerr << "Error in reading input data" << std::endl;
				return 1;
			} 

		}
		txt2vtk(ref);
		std::cerr << "Type of Ref line is .vtk" << std::endl;
	}
	if(strcmp(chext,".txt")==0)// read vtk file
	{
		if (! jdq2017::readCenterline(chRefCurvefilename, ref) )
		{


			return 0;
		}
		std::cerr << "Type of Ref line is .txt" << std::endl;
	}
	char* chext_tar=GetPathEXT(chTarCurvefilename);
	if(strcmp(chext_tar,".vtk")==0)// read vtk file
	{
		if (! jdq2017::readCenterlinevtk(chTarCurvefilename, cl) )
		{
			if (!outputOneLine) 
			{
				std::cerr << "Error in reading input data" << std::endl;
				return 1;
			} 
		}
		txt2vtk(cl);
		std::cerr << "Type of Tsr line is .vtk" << std::endl;
	}
	if(strcmp(chext_tar,".txt")==0)// read vtk file
	{
		if (! jdq2017::readCenterline(chTarCurvefilename, cl) )
		{


			return 0;
		}
		std::cerr << "Type of Ref line is .txt" << std::endl;
	}
	  // ** Resample centerline to same sampling as reference standard
  // Calculate reference sampling distance
  double refSampling = pathLength(ref)/(ref.size()-1);
  if (debugOutput)
  {
    std::cout << "Reference sampling distance = " << refSampling << std::endl;
  }
  // Calculate method sampling distance
  double clSampling = pathLength(cl)/(cl.size()-1);
  if (debugOutput)
  {
    std::cout << "Centerline sampling distance = " << clSampling << std::endl;
  }
  
  //Check if the sampling distances differ more than 0.1%
  if ( refSampling/clSampling < 0.999 || refSampling/clSampling > 1.001 )
  {
    if (debugOutput) std::cout << "Resampling centerline ..." << std::endl;
	//std::cout<<pathLength(cl)<<std::endl;
    //Resample method centerline
    jdq2017::ResamplePaths<jdq2017::point3D> resampler;
    resampler(cl, refSampling);
    cl = resampler.resultPath();
	//std::cout<<pathLength(cl)<<std::endl;
    if (debugOutput)
    {
      std::cout << "Centerline resampled: " << std::endl;
      for (int i = 0; i < int(cl.size()); ++i)
        std::cout << "Point: " << cl[i]._x << "," << cl[i]._y << "," << cl[i]._z << std::endl;
    }
  }



		int nLen0 = strlen(chResamfilename)+ strlen(".vtk") + 1;
		char *chFileName0 = (char *)malloc(nLen0);
		strcpy(chFileName0, chResamfilename);
		strcat(chFileName0, ".vtk");
			Write2Vtk_jdq(cl,chFileName0);
		free(chFileName0);

			int nLen1= strlen(chResamfilename)+ strlen(".txt") + 1;
		char *chFileName1 = (char *)malloc(nLen1);
		strcpy(chFileName1, chResamfilename);
		strcat(chFileName1, ".txt");
			Write2Txt_jdq(cl,chFileName1);
		free(chFileName0);
	//WriteCA2Txt_Skip(vUnorgaPointsWorld_ori,Points_3D_ori_Filename);

	
	
}

 