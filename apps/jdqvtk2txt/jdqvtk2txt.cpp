

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
#include<iostream>
using namespace std;
void Help()
{
	std::cout<<" An simple example for extraction centerlines  of CTA images, target.nii.gz and source.nii.gz, result save as res: \n" ;
	std::cout<<" zxhsemi0 -test target.nii.gz -ref source.nii.gz -o result0 -ffd 20 20 20 -bending 0.001\n";
	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -pre 0 result0.FFD -ffd 20 20 20 -ffd 10 10 10 -Reg 2 -sub 2 2 2 -sub 1 1 1 -bending 0.001\n";
	std::cout<<"OR \n";
	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -hierarchy 3 -bending 0.0031\n\n";

	std::cout<<"  <-test/target img.>     (test or target image)\n" ;
	std::cout<<"  <-ref/source img.>      (reference or source image)\n" ;
	std::cout<<"  <-o savename>           (string for saving transformed -ref/-source image, file names prefix)\n" ;
	std::cout<<"  <-maskt img.>           (mask image on test image, use -maskr on ref image) \n" ;
	std::cout<<"  USE -ffd: zxhsemi0 to fast init and get the .FFD for setting -pre, and then set -ffd\n" ;
	std::cout<<"  <-ffd fx fy fz>         (spacing of free form deformation, FFD, multi-input related to -Reg)\n" ;
	std::cout<<"  <-pre 0 s>              (pre set transformation field)\n";
	std::cout<<"  <-sub fx fy fz [ft]>    ([3 3 3], sampling spacing; positive: millimeters interval; negative: pixels interval)\n";
	std::cout<<"  <-Reg i>                ([1] number of multi-level registrations)\n" ;
	std::cout<<"  OR USE -hierarchy, simple and not need to set -ffd,-sub,-Reg:\n" ; 
	std::cout<<"  <-hierarchy n>          ([3] number of multi-level FFD registration, where\n";
	std::cout<<"                           the first level of -ffd spacing is one forth of test image extent, and halve in subsequence level\n" ; 
	std::cout<<"                           the final level of -sub sampling is pixel size of the test image\n" ; 	 
	std::cout<<"\n";
	std::cout<<"  <-bending f..f>         ([0.001]weighting for bending energy, recommend f=0.001~0.01)\n" ;
} 
void HELP()
{ 
	Help() ;
	std::cout<<"------------------------------------------------------\n" ;
	std::cout<<"  OPTIONS of gradient optimization computation; use setting in previous -Reg when un-set \n" ; 
}  





#include "miiMinPathModel.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkPolyLine.h"

#include <time.h>

#define TAB_CHAR	9
#define M_PI 3.14159265358979323846

int GetTime() ;


typedef struct
{
	double x;
	double y;
	double z;
}PointCordTypeDef;


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
void WriteCA2Txt(const char *chFileName,const vector<miiCNode<double,float> >&vMinPathWorld)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = vMinPathWorld.size();
	float fImgPixel[3];	
	if (WriteFileTxt.is_open())   
	{  

		for (int i = 0; i < nPointNum; i++)
		{
			fImgPixel[0] = vMinPathWorld[i].x;
			fImgPixel[1] = vMinPathWorld[i].y;
			fImgPixel[2] = vMinPathWorld[i].z;
			WriteFileTxt << right<<fixed<<setprecision(4) <<-fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] << '\n';

		}	
	}
	WriteFileTxt.close();
}

int main(int argc, char *argv[])//Use Resolution Mask to extraction
{
	

if( argc <2 )
{
	std::cerr << "Usage: " << std::endl;
	std::cerr << "jdqvtk2txt	vtkfile  resulttxtfile" << std::endl;
	return -1;

}

char *chvtkFileName=argv[1];
char *chtxtFileName=argv[2];
	/*char *chvtkFileName="F:/Coronary_0/ZXHJDQCAEDMP/trainningLLDMPInt/dataset00/model01/vessel0/vessel0.vtk";
	char *chtxtFileName="F:/Coronary_0/ZXHJDQCAEDMP/trainningLLDMPInt/dataset00/vessel0/MCLine0.txt";*/
	// read model points from "*.vtk" file
	vector< miiCNode<double, float> > vtkPoints;
	ReadVtk(chvtkFileName,vtkPoints);
	// save CA to txt file
	WriteCA2Txt(chtxtFileName,vtkPoints);
	cout << "Transmit down!\n" << endl;
	return 1;	

}



 