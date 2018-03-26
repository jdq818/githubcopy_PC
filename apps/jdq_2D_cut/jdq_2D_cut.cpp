

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
//void Help()
//{
//	std::cout<<" An simple example for registration of images, target.nii.gz and source.nii.gz, result save as res: \n" ;
//	std::cout<<" zxhsemi0 -test target.nii.gz -ref source.nii.gz -o result0 -ffd 20 20 20 -bending 0.001\n";
//	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -pre 0 result0.FFD -ffd 20 20 20 -ffd 10 10 10 -Reg 2 -sub 2 2 2 -sub 1 1 1 -bending 0.001\n";
//	std::cout<<"OR \n";
//	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -hierarchy 3 -bending 0.0031\n\n";
//
//	std::cout<<"  <-test/target img.>     (test or target image)\n" ;
//	std::cout<<"  <-ref/source img.>      (reference or source image)\n" ;
//	std::cout<<"  <-o savename>           (string for saving transformed -ref/-source image, file names prefix)\n" ;
//	std::cout<<"  <-maskt img.>           (mask image on test image, use -maskr on ref image) \n" ;
//	std::cout<<"  USE -ffd: zxhsemi0 to fast init and get the .FFD for setting -pre, and then set -ffd\n" ;
//	std::cout<<"  <-ffd fx fy fz>         (spacing of free form deformation, FFD, multi-input related to -Reg)\n" ;
//	std::cout<<"  <-pre 0 s>              (pre set transformation field)\n";
//	std::cout<<"  <-sub fx fy fz [ft]>    ([3 3 3], sampling spacing; positive: millimeters interval; negative: pixels interval)\n";
//	std::cout<<"  <-Reg i>                ([1] number of multi-level registrations)\n" ;
//	std::cout<<"  OR USE -hierarchy, simple and not need to set -ffd,-sub,-Reg:\n" ; 
//	std::cout<<"  <-hierarchy n>          ([3] number of multi-level FFD registration, where\n";
//	std::cout<<"                           the first level of -ffd spacing is one forth of test image extent, and halve in subsequence level\n" ; 
//	std::cout<<"                           the final level of -sub sampling is pixel size of the test image\n" ; 	 
//	std::cout<<"\n";
//	std::cout<<"  <-bending f..f>         ([0.001]weighting for bending energy, recommend f=0.001~0.01)\n" ;
//	std::cout<<"  OPTIONS of spatially encoding scheme\n"  ;// similarity computation, default normalized mutual information\n" ; 
//	std::cout<<"  <-semiradius f...f>     (radius of local region, default set to twice ffd spacing)\n";
//	//std::cout<<"  <-semiwidth f...f>      (or -semisize, width/size of local region, default 80mm)\n";
//	//std::cout<<"  <-semikernel s>         ([0], -1='Gaussian', 0='ZeroBSpline', 3='BSpline')\n";
//
//	std::cout<<"\n" ;
//
//} 
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


#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"


#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"

#include "jdq2017util.h"
#include "jdqdijkstra.h"
#include "jdqPoint.h"

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
#include<stdio.h>

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"
#include "zxhImageModelingLinear.h"

#include <Eigen/Dense>  
using namespace Eigen;  
using namespace Eigen::internal;  
using namespace Eigen::Architecture; 

// change to true for debug output
const bool debugOutput = false;
// Change to true to output all the measures on one line
// this line will start with the dataset and vessel number
const bool outputOneLine = true;
#define TAB_CHAR	9
#define M_PI 3.14159265358979323846
#define SPHERE_RADIUS 4
using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
	short val;
}PointCordTypeDef;

typedef struct
{
	int x;
	int y;
	int z;
	short val;
}PointImgTypeDef;


char *GetPathEXT(char *chFileName)
{  char path_buffer[_MAX_PATH];  
   char drive[_MAX_DRIVE];  
   char dir[_MAX_DIR];  
   char fname[_MAX_FNAME];  
   char ext[_MAX_EXT];  
   _splitpath( chFileName, drive, dir, fname, ext );  
   return ext;
}
void init_img(zxhImageDataT<short> &imgCountmap)
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
		iPoints->InsertNextPoint(PointCord[i]._x, PointCord[i]._y, PointCord[i]._z);
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
		fImgPixel[0] =-1*vPoints[i]._x;
		fImgPixel[1] =-1*vPoints[i]._y;
		fImgPixel[2] = vPoints[i]._z;
		WriteFileTxt << right<<fixed<<setprecision(4) <<fImgPixel[0] << " " << fImgPixel[1] << " " << fImgPixel[2] << '\n';
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

bool BoundaryCorrect(int *PointPos,int ImgNewVslsSize[4])
{
	for (int i=0;i<3;i++)
	{
	PointPos[i]=zxh::maxf(0,PointPos[i]);
	PointPos[i]=zxh::minf(ImgNewVslsSize[i]-1,PointPos[i]);
	}
	return true;
}


bool ReadPointTxt(char *filename,vector< PointCordTypeDef> &cl)
{
	string strNum;
	int nStart = 0, nEnd = 0;
	PointCordTypeDef strctTempPoint;
	if(cl.size()>0)
		cl.clear();
	ifstream iFileIn(filename,ios::in);

	if(!iFileIn)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << filename << endl;
		return false; //exit(1); // terminate with error
	}
	else
	{
		while(!iFileIn.eof())
		{
			iFileIn>>strctTempPoint.x>>strctTempPoint.y>>strctTempPoint.z;
			cl.push_back(strctTempPoint);
		}
	}
	
	return true;
}
bool Calc_divec(int i,float fdivec[3],std::vector<jdq2017::point3D> &ref)
{
double clSampling = pathLength(ref)/ref.size();
int nFBsize=1/clSampling;
int nFPosi=zxh::minf(i+nFBsize,ref.size()-1);
int nBPosi=zxh::maxf(i-nFBsize,0);
jdq2017::point3D curPontW=ref[i];
jdq2017::point3D FPontW=ref[nFPosi];
jdq2017::point3D BPontW=ref[nBPosi];
float fvec1[3]={curPontW._x-BPontW._x,curPontW._y-BPontW._y,curPontW._z-BPontW._z};
float fvec2[3]={FPontW._x-curPontW._x,FPontW._y-curPontW._y,FPontW._z-curPontW._z};
fdivec[0]=0.5*(fvec1[0]+fvec2[0]);
fdivec[1]=0.5*(fvec1[1]+fvec2[1]);
fdivec[2]=0.5*(fvec1[2]+fvec2[2]);
zxh::VectorOP_Normalise(fdivec,3);
return true;
}
bool Get_related_pix_z(zxhImageDataT<short>&imgReadRaws,float fnormdi[3],int ncurPontI[3],vector<PointImgTypeDef>&v_propont)
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	float fcenpontW[3]={ncurPontI[0],ncurPontI[1],ncurPontI[2]};
	imgReadRaws.GetImageInfo()->ImageToWorld(fcenpontW);
     float res=0.1;
	 zxh::VectorOP_Normalise(fnormdi,3);
	 int nindex[3]={-1,-1,-1};
	 //mark
		short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
		for(int z=0;z<nimgSize[2];z++)  
    {  
        for (int y=0;y<nimgSize[1];y++)  
        {  
            for (int x=0;x<nimgSize[0];x++)  
            {  
                BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
                
            }  
        }  
    }  
	for(int i=-100;i<100;i++)
	{
		float flen=res*i;
		float fmove[3]={flen*fnormdi[0],flen*fnormdi[1],flen*fnormdi[2]};
		float fnewpontW[3]={fcenpontW[0]+ fmove[0],fcenpontW[1]+ fmove[1],fcenpontW[2]+ fmove[2]};
		int nnewpoint[3]={0,0,0};
		imgReadRaws.GetImageInfo()->WorldToImage(fnewpontW);
		nnewpoint[0]= (int)(fnewpontW[0]+0.5);
		nnewpoint[1]= (int)(fnewpontW[1]+0.5);
		nnewpoint[2]= (int)(fnewpontW[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
		short smar=BW[nnewpoint[2]* nimgSize[1] * nimgSize[0] + nnewpoint[1]* nimgSize[0] + nnewpoint[0]];
		if( smar==0)
		{
			PointImgTypeDef tmppropont;
			tmppropont.x=nnewpoint[0];
			tmppropont.y=nnewpoint[1];
			tmppropont.z=nnewpoint[2];
			v_propont.push_back(tmppropont);
			BW[nnewpoint[2]* nimgSize[1] * nimgSize[0] + nnewpoint[1]* nimgSize[0] + nnewpoint[0]]=1;
		}
	}
	delete[]BW;
	return true;
}
bool Get_related_pix_y(zxhImageDataT<short>&imgReadRaws,float fnormdi[3],int ncurPontI[3],vector<PointImgTypeDef>&v_propont)
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	float fcenpontW[3]={ncurPontI[0],ncurPontI[1],ncurPontI[2]};
	imgReadRaws.GetImageInfo()->ImageToWorld(fcenpontW);
     float res=0.1;
	 zxh::VectorOP_Normalise(fnormdi,3);
	 int nindex[3]={-1,-1,-1};
	 //mark
	 std::vector<std::vector<std::vector<int> > > BW(nimgSize[2],vector<vector<int> >(nimgSize[1],vector<int>(nimgSize[0],0))); 
	for(int z=0;z<nimgSize[2];z++)  
    {  
        for (int y=0;y<nimgSize[1];y++)  
        {  
            for (int x=0;x<nimgSize[0];x++)  
            {  
                BW[z][y][x]=0;  
                
            }  
        }  
    }  
	for(int i=-100;i<100;i++)
	{
		float flen=res*i;
		float fmove[3]={flen*fnormdi[0],flen*fnormdi[1],flen*fnormdi[2]};
		float fnewpontW[3]={fcenpontW[0]+ fmove[0],fcenpontW[1]+ fmove[1],fcenpontW[2]+ fmove[2]};
		int nnewpoint[3]={0,0,0};
		imgReadRaws.GetImageInfo()->WorldToImage(fnewpontW);
		nnewpoint[0]= (int)(fnewpontW[0]+0.5);
		nnewpoint[1]= (int)(fnewpontW[1]+0.5);
		nnewpoint[2]= (int)(fnewpontW[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
		if(BW[nnewpoint[2]][nnewpoint[1]][nnewpoint[0]]==0)
		{
			PointImgTypeDef tmppropont;
			tmppropont.x=nnewpoint[0];
			tmppropont.y=nnewpoint[1];
			tmppropont.z=nnewpoint[2];
			v_propont.push_back(tmppropont);
			BW[nnewpoint[2]][nnewpoint[1]][nnewpoint[0]]=1;
		}
	}
	return true;
}
bool Get_related_pix_x(zxhImageDataT<short>&imgReadRaws,float fnormdi[3],int ncurPontI[3],vector<PointImgTypeDef>&v_propont)
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	float fcenpontW[3]={ncurPontI[0],ncurPontI[1],ncurPontI[2]};
	imgReadRaws.GetImageInfo()->ImageToWorld(fcenpontW);
     float res=0.1;
	 zxh::VectorOP_Normalise(fnormdi,3);
	 int nindex[3]={-1,-1,-1};
	 //mark
	 std::vector<std::vector<std::vector<int> > > BW(nimgSize[2],vector<vector<int> >(nimgSize[1],vector<int>(nimgSize[0],0))); 
	for(int z=0;z<nimgSize[2];z++)  
    {  
        for (int y=0;y<nimgSize[1];y++)  
        {  
            for (int x=0;x<nimgSize[0];x++)  
            {  
                BW[z][y][x]=0;  
                
            }  
        }  
    }  
	for(int i=-100;i<100;i++)
	{
		float flen=res*i;
		float fmove[3]={flen*fnormdi[0],flen*fnormdi[1],flen*fnormdi[2]};
		float fnewpontW[3]={fcenpontW[0]+ fmove[0],fcenpontW[1]+ fmove[1],fcenpontW[2]+ fmove[2]};
		int nnewpoint[3]={0,0,0};
		imgReadRaws.GetImageInfo()->WorldToImage(fnewpontW);
		nnewpoint[0]= (int)(fnewpontW[0]+0.5);
		nnewpoint[1]= (int)(fnewpontW[1]+0.5);
		nnewpoint[2]= (int)(fnewpontW[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
		if(BW[nnewpoint[2]][nnewpoint[1]][nnewpoint[0]]==0)
		{
			PointImgTypeDef tmppropont;
			tmppropont.x=nnewpoint[0];
			tmppropont.y=nnewpoint[1];
			tmppropont.z=nnewpoint[2];
			v_propont.push_back(tmppropont);
			BW[nnewpoint[2]][nnewpoint[1]][nnewpoint[0]]=1;
		}
	}
	return true;
}
bool Cutimg_byLines_2D(zxhImageDataT<short>&imgReadRaws,std::vector<jdq2017::point3D> &ref)
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	
	//built mark
	std::vector<std::vector<std::vector<int> > > BW(nimgSize[2],vector<vector<int> >(nimgSize[1],vector<int>(nimgSize[0],0))); 
	for(int z=0;z<nimgSize[2];z++)  
    {  
        for (int y=0;y<nimgSize[1];y++)  
        {  
            for (int x=0;x<nimgSize[0];x++)  
            {  
                BW[z][y][x]=0;  
                
            }  
        }  
    }  
	float fzaxi[3]={0,0,1};
	float fyaxi[3]={0,1,0};
	float fxaxi[3]={1,0,0};
	//visit each node
	vector<vector<PointImgTypeDef>> vv_proponts;
	for (int i=0;i<ref.size();i++)
	{
		jdq2017::point3D curPontW=ref[i];
		float fcurPontI[3]={curPontW._x,curPontW._y,curPontW._z};
		imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI);
		int ncurPontI[3];
		ncurPontI[0]= (int)(fcurPontI[0]+0.5);
		ncurPontI[1]= (int)(fcurPontI[1]+0.5);
		ncurPontI[2]= (int)(fcurPontI[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
		vector<PointImgTypeDef> v_propont;
		if(BW[ncurPontI[2]][ncurPontI[1]][ncurPontI[0]]==0)
		{
			float fdivec[3]={0,0,0};
			Calc_divec(i,fdivec,ref);
			float fcosthetaz=zxh::VectorOP_Cosine(fdivec,fzaxi,3);
			float fcosthetay=zxh::VectorOP_Cosine(fdivec,fyaxi,3);
			float fcosthetax=zxh::VectorOP_Cosine(fdivec,fxaxi,3);
			float fnormdi[3]={0,0,0};
			if((zxh::absf(fcosthetay)>=zxh::absf(fcosthetaz))&&(zxh::absf(fcosthetax)>=zxh::absf(fcosthetaz)))//we see from z-axis
			{
				zxh::VectorOP_CrossProduct3D(fdivec,fzaxi,fnormdi);

	           Get_related_pix_z(imgReadRaws,fnormdi,ncurPontI,v_propont);
				//y=f(x)
			}
			else if((zxh::absf(fcosthetay)>=zxh::absf(fcosthetax))&&(zxh::absf(fcosthetaz)>=zxh::absf(fcosthetax)))//we see from x-axis
			{
				zxh::VectorOP_CrossProduct3D(fdivec,fxaxi,fnormdi);
				Get_related_pix_z(imgReadRaws,fnormdi,ncurPontI,v_propont);
				
			}
				else if((zxh::absf(fcosthetaz)>=zxh::absf(fcosthetay))&&(zxh::absf(fcosthetax)>=zxh::absf(fcosthetay)))//we see from y-axis
			{
				zxh::VectorOP_CrossProduct3D(fdivec,fyaxi,fnormdi);
				Get_related_pix_z(imgReadRaws,fnormdi,ncurPontI,v_propont);
				
			}

			vv_proponts.push_back(v_propont);

			BW[ncurPontI[2]][ncurPontI[1]][ncurPontI[0]]=1;
			if(vv_proponts.size()==300)
			{
				zxhImageDataT<short>Outputimg;
				Outputimg.NewImage(imgReadRaws.GetImageInfo());
				for (int k=0;k<vv_proponts.size();k++)
				{
					vector<PointImgTypeDef> vtmp=vv_proponts[k];
					for (int h=0;h<vtmp.size();h++)
					{
						PointImgTypeDef tmppont=vtmp[h];
						short inten=imgReadRaws.GetPixelGreyscale(tmppont.x,tmppont.y,tmppont.z);
						Outputimg.SetPixelByGreyscale(k,h,0,0,inten);

					}
				}
				char *chResultName="E:/work_jdq/rebuttal_media18/rebuttal/Exp/2D_cut/img_2D.nii.gz";
				string chFileName2(chResultName);
				zxh::SaveImage(&Outputimg,chFileName2.c_str());
			}
		}
	}
	return true;

}
bool Segint_Proj_image(zxhImageDataT<short>&imgReadRaws,int ncurPontI[3],zxhImageDataT<short>&Outputimg)
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	int ixyz[3]={0,0,0};
	for(int iz=-25;iz<25;iz++)
		for(int iy=-25;iy<25;iy++)
		{
			int iiz=ncurPontI[2]+iz;
			int iiy=ncurPontI[1]+iy;
			short maxintx=0;

			for(int ix=-10;ix<10;ix++)
			{
				int iix=ncurPontI[0]+ix;
				 ixyz[0]=iix;
				  ixyz[1]=iiy;
				   ixyz[2]=iiz;
				BoundaryCorrect(ixyz,nimgSize);
				short scurint=imgReadRaws.GetPixelGreyscale(ixyz[0],ixyz[1],ixyz[2]);
				if(scurint>maxintx)
				{
					maxintx=scurint;
				}
			}
			short spreint=Outputimg.GetPixelGreyscale(0,ixyz[1],ixyz[2],0);
			if(spreint>maxintx)
				Outputimg.SetPixelByGreyscale(0,ixyz[1],ixyz[2],0,spreint);
			else
				Outputimg.SetPixelByGreyscale(0,ixyz[1],ixyz[2],0,maxintx);
		}

	
	

	return true;
}

	bool Gen_Int_inWin(float fdivec[3],zxhImageDataT<short>&imgReadRaws,int ncurPontI[3],zxhImageDataT<short>&Outputimg)
{
	zxh::VectorOP_Normalise(fdivec,3);
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
		// get the image data
	const short *sImData = imgReadRaws.GetImageData();
	int winsize=25;
	float fcenpontW[3]={ncurPontI[0],ncurPontI[1],ncurPontI[2]};
	imgReadRaws.GetImageInfo()->ImageToWorld(fcenpontW);
	int nzmin=zxh::maxf(ncurPontI[0]-30,0);
	int nzmax=zxh::minf(ncurPontI[0]+30,nimgSize[2]);
	int nymin=zxh::maxf(ncurPontI[1]-30,0);
	int nymax=zxh::minf(ncurPontI[1]+30,nimgSize[1]);
	int nxmin=zxh::maxf(ncurPontI[2]-30,0);
	int nxmax=zxh::minf(ncurPontI[2]+30,nimgSize[0]);
		for(int z=nzmin;z<=nzmax;z++)   
        for(int y=nymin;y<=nymax;y++)   
           for(int x=nxmin;x<=nxmax;x++)   
            {  
				float ftmppont[3]={x,y,z};
				imgReadRaws.GetImageInfo()->ImageToWorld(ftmppont);
				
				float fdist=zxh::VectorOP_Distance(ftmppont,fcenpontW,3);
				float fa=fdivec[0];
				float fb=fdivec[1];
				float fc=fdivec[2];
				float fd=-1*fdivec[0]*ftmppont[0];
				float fdistplane=zxh::absf(fdivec[0]*(ftmppont[0]-fcenpontW[0])+fdivec[1]*(ftmppont[1]-fcenpontW[1])+fdivec[2]*(ftmppont[2]-fcenpontW[2]));
			}
	

	return true;
}
bool Loca_pro_2D(zxhImageDataT<short>&imgReadRaws,std::vector<jdq2017::point3D> &ref)
{
	zxhImageDataT<short>Outputimg;
	Outputimg.NewImage(imgReadRaws.GetImageInfo());
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	
	//built mark
short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
		for(int z=0;z<nimgSize[2];z++)  
    {  
        for (int y=0;y<nimgSize[1];y++)  
        {  
            for (int x=0;x<nimgSize[0];x++)  
            {  
                BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
				Outputimg.SetPixelByGreyscale(x,y,z,0,0);
            }  
        }  
    }  
	float fzaxi[3]={0,0,1};
	float fyaxi[3]={0,1,0};
	float fxaxi[3]={1,0,0};
	//visit each node
	vector<vector<PointImgTypeDef>> vv_proponts;
	for (int i=0;i<ref.size();i++)
	{
		jdq2017::point3D curPontW=ref[i];
		float fcurPontI[3]={curPontW._x,curPontW._y,curPontW._z};
		imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI);
		int ncurPontI[3];
		ncurPontI[0]= (int)(fcurPontI[0]+0.5);
		ncurPontI[1]= (int)(fcurPontI[1]+0.5);
		ncurPontI[2]= (int)(fcurPontI[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
			short smar=BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]];
		if( smar==0)
		{
		//set the intensity of the projected image
		Segint_Proj_image(imgReadRaws,ncurPontI,Outputimg);
		BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]]=1;
		}

	}
	char *chResultName="E:/work_jdq/rebuttal_media18/rebuttal/Exp/2D_cut/img_2D.nii.gz";
				string chFileName2(chResultName);
				zxh::SaveImage(&Outputimg,chFileName2.c_str());
	return true;
}
bool Gen_sca_ponts(vector<PointCordTypeDef>&vscaponts)
{
	float res=0.25;
	int npnum=150;
	for(int i=-1*npnum;i<=npnum;i++)
		for(int j=-1*npnum;j<=npnum;j++)
		{
			float fx=res*i;
			float fy=res*j;
			PointCordTypeDef temp;
			temp.x=fx;
			temp.y=fy;
			temp.z=0;
			vscaponts.push_back(temp);
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
		
		return true;
}
bool Rotaz_to_vec(float fabc[3],Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2)
{
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
	return true;
}
bool Rotpon_With_Matr(vector<PointCordTypeDef> &vscaponts, Matrix<float,1,3>&MPro_Pcent,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,vector<PointCordTypeDef> &vsc_rot)
{
	for(int i=0;i<vscaponts.size();i++)
	{
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<vscaponts[i].x,vscaponts[i].y,vscaponts[i].z;
		MnewPont=MoriPont*Mrx2.inverse()*Mrx1.inverse();
		PointCordTypeDef PnewPont;

		PnewPont.x=MnewPont(0,0)+MPro_Pcent(0,0);
		PnewPont.y=MnewPont(0,1)+MPro_Pcent(0,1);
		PnewPont.z=MnewPont(0,2)+MPro_Pcent(0,2);
		vsc_rot.push_back(PnewPont);
		//cout<<MnewPont<<endl;
	}
	return true;
}
bool WorldtoImage(zxhImageDataT<short>&imgReadRaws,PointCordTypeDef &curPontW,PointImgTypeDef &curPontI)
{
	

	float fcurPontI[3]={curPontW.x,curPontW.y,curPontW.z};
	imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI);
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	int ncurPontI[3];
	ncurPontI[0]= (int)(fcurPontI[0]+0.5);
	ncurPontI[1]= (int)(fcurPontI[1]+0.5);
	ncurPontI[2]= (int)(fcurPontI[2]+0.5);
	BoundaryCorrect(ncurPontI,nimgSize);
	curPontI.x=ncurPontI[0];
	curPontI.y=ncurPontI[1];
	curPontI.z=ncurPontI[2];
	return true;
}
bool ImagetoWorld(zxhImageDataT<short>&imgReadRaws,PointImgTypeDef &curPontI,PointCordTypeDef &curPontW)
{
	

	float fcurPontI[3]={curPontI.x,curPontI.y,curPontI.z};
	imgReadRaws.GetImageInfo()->ImageToWorld(fcurPontI);
	curPontW.x=fcurPontI[0];
	curPontW.y=fcurPontI[1];
	curPontW.z=fcurPontI[2];
	return true;
}
bool Rotpon_With_Matr12(zxhImageDataT<short>&imgReadRaws,vector<PointCordTypeDef> &vscaponts, Matrix<float,1,3>&MPro_Pcent,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,vector<PointCordTypeDef> &vsc_rot,float fmovvec[3])
{
	
	for(int i=0;i<vscaponts.size();i++)
	{
		PointCordTypeDef curPontW;
		curPontW.x=vscaponts[i].x;
		curPontW.y=vscaponts[i].y;
		curPontW.z=vscaponts[i].z;
	
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<curPontW.x,curPontW.y,curPontW.z;
		MnewPont=MoriPont-MPro_Pcent;
		PointCordTypeDef PnewPont;
		PnewPont.x=(MnewPont*Mrx1*Mrx2)(0,0);
		PnewPont.y=(MnewPont*Mrx1*Mrx2)(0,1);
		PnewPont.z=(MnewPont*Mrx1*Mrx2)(0,2);
		PnewPont.val=vscaponts[i].val;
		vsc_rot.push_back(PnewPont);
		
		//cout<<MnewPont<<endl;
	}
	return true;
}
bool Rot_ponts_vec_zaxis(vector<PointCordTypeDef> &vscaponts, Matrix<float,1,3>&MPro_Pcent,float fdivec[3],vector<PointCordTypeDef> &vsc_rot,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2)
{
	

	//rotate the vector to z-axis
	Rotaz_to_vec(fdivec,Mrx1,Mrx2);	

	//rotate the points to the local axis
	Rotpon_With_Matr(vscaponts,MPro_Pcent,Mrx1,Mrx2,vsc_rot);
	return true;
}
//bool WorldtoImage(zxhImageDataT<short>&imgReadRaws,PointCordTypeDef &curPontW,PointImgTypeDef &curPontI)
//{
//	zxhImageModelingLinear mod ;
//	mod.SetImage£¨&imgReadRaws);
//
//	
//	short inten = mod.GetPixelFloatValueWithCheckByWorld(curPontW.x,curPontW.y,curPontW.z) ; 
//	curPontI.val=inten;
//	/*
//	float fcurPontI[3]={curPontW.x,curPontW.y,curPontW.z};
//	imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI);
//	int nimgSize[4];
//	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
//	int ncurPontI[3];
//	ncurPontI[0]= (int)(fcurPontI[0]+0.5);
//	ncurPontI[1]= (int)(fcurPontI[1]+0.5);
//	ncurPontI[2]= (int)(fcurPontI[2]+0.5);
//	BoundaryCorrect(ncurPontI,nimgSize);
//	short inten=imgReadRaws.GetPixelGreyscale(ncurPontI[0],ncurPontI[1],ncurPontI[2],0);
//	curPontI.x=ncurPontI[0];
//	curPontI.y=ncurPontI[1];
//	curPontI.z=ncurPontI[2];
//	curPontI.val=inten;
//	curPontW.val=inten;*/
//	return true;
//}

bool RotBack(zxhImageDataT<short>&imgReadRaws,vector<PointCordTypeDef>&vsc_rot, Matrix<float,1,3>&MPro_Pcent,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,vector<PointCordTypeDef> &vsc_rotbac_int)
{

	for (int i=0;i<vsc_rot.size();i++)
	{
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<vsc_rot[i].x,vsc_rot[i].y,vsc_rot[i].z;
		MoriPont(0,0)=MoriPont(0,0)-MPro_Pcent(0,0);
		MoriPont(0,1)=MoriPont(0,1)-MPro_Pcent(0,1);
		MoriPont(0,2)=MoriPont(0,2)-MPro_Pcent(0,2);

		MnewPont=MoriPont*Mrx1*Mrx2;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		PnewPont.z=MnewPont(0,2);
		PnewPont.val=vsc_rot[i].val;
		vsc_rotbac_int.push_back(PnewPont);
	}
	return true;
}
bool RotBack_corr(vector<PointCordTypeDef>&vsc_rot, Matrix<float,1,3>&MPro_Pcent,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,vector<PointCordTypeDef> &vsc_rotbac_int,jdq2017::point3D newPontM,jdq2017::point3D &newPontM_rob)
{

	for (int i=0;i<vsc_rot.size();i++)
	{
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<vsc_rot[i].x,vsc_rot[i].y,vsc_rot[i].z;
		MoriPont(0,0)=MoriPont(0,0)-MPro_Pcent(0,0);
		MoriPont(0,1)=MoriPont(0,1)-MPro_Pcent(0,1);
		MoriPont(0,2)=MoriPont(0,2)-MPro_Pcent(0,2);

		MnewPont=MoriPont*Mrx1*Mrx2;
		PointCordTypeDef PnewPont;
		PnewPont.x=MnewPont(0,0);
		PnewPont.y=MnewPont(0,1);
		PnewPont.z=MnewPont(0,2);
		PnewPont.val=vsc_rot[i].val;
		vsc_rotbac_int.push_back(PnewPont);
	}
	MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<newPontM._x,newPontM._y,newPontM._z;
		MoriPont(0,0)=MoriPont(0,0)-MPro_Pcent(0,0);
		MoriPont(0,1)=MoriPont(0,1)-MPro_Pcent(0,1);
		MoriPont(0,2)=MoriPont(0,2)-MPro_Pcent(0,2);

		MnewPont=MoriPont*Mrx1*Mrx2;
		newPontM_rob._x=MnewPont(0,0);
		newPontM_rob._y=MnewPont(0,1);
		newPontM_rob._z=MnewPont(0,2);
	return true;
}
bool Get_int_img(zxhImageDataT<short>&imgReadRaws,vector<PointCordTypeDef>&vsc_rot)
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	/*short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
	for(int z=0;z<nimgSize[2];z++)  
	{  
		for (int y=0;y<nimgSize[1];y++)  
		{  
			for (int x=0;x<nimgSize[0];x++)  
			{  
				BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
			}  
		}  
	} */ 
	for (int i=0;i<vsc_rot.size();i++)
	{
		PointCordTypeDef curPontW;
		curPontW=vsc_rot[i];
		PointImgTypeDef curPontI;
		WorldtoImage(imgReadRaws,curPontW,curPontI);
		vsc_rot[i].val=curPontW.val;
		//short smar=BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x];
	
		//BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x]=1;

	}
	return true;
}
	bool Tran_metp_img(zxhImageDataT<short>&imgReadRaws,jdq2017::point3D newPontM_rob,PointImgTypeDef &newPontM_img,float fmovvec[3])
{
	
		PointCordTypeDef curPontW;
		curPontW.x=newPontM_rob._x+fmovvec[0];
		curPontW.y=newPontM_rob._y+fmovvec[1];
		curPontW.z=newPontM_rob._z+fmovvec[2];
		WorldtoImage(imgReadRaws,curPontW,newPontM_img);
		return true;
}
bool Tans_img(zxhImageDataT<short>&imgReadRaws,vector<PointCordTypeDef>&vsc_rotbac_int,vector<PointImgTypeDef>&vsc_img,float fmovvec[3])
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
	for(int z=0;z<nimgSize[2];z++)  
	{  
		for (int y=0;y<nimgSize[1];y++)  
		{  
			for (int x=0;x<nimgSize[0];x++)  
			{  
				BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
			}  
		}  
	} 

	for (int i=0;i<vsc_rotbac_int.size();i++)
	{
		PointCordTypeDef curPontW;
	    curPontW.x=vsc_rotbac_int[i].x+fmovvec[0];
		curPontW.y=vsc_rotbac_int[i].y+fmovvec[1];
		curPontW.z=vsc_rotbac_int[i].z+fmovvec[2];
		PointImgTypeDef curPontI;
		WorldtoImage(imgReadRaws,curPontW,curPontI);
		short smar=BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x];
		if(smar==0)
		{
			curPontI.val=vsc_rotbac_int[i].val;
			vsc_img.push_back(curPontI);
			BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x]=1;
		}

	}
	delete[]BW;
	return true;
}
bool CMP_inWin(zxhImageDataT<short>&imgReadRaws,std::vector<jdq2017::point3D> &ref,zxhImageDataT<short>&Outputimg)
{
	
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	//calculate the distance from the center
	float norigin[3]={0,0,0};
	float nimgcen[3]={0.5*nimgSize[0],0.5*nimgSize[1],0};
	imgReadRaws.GetImageInfo()->ImageToWorld(norigin);
	imgReadRaws.GetImageInfo()->ImageToWorld(nimgcen);
	float fmovvec[3];
	zxh::VectorOP_Substract(nimgcen,norigin,fmovvec,3);
	
	//built mark
	short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
	for(int z=0;z<nimgSize[2];z++)  
	{  
		for (int y=0;y<nimgSize[1];y++)  
		{  
			for (int x=0;x<nimgSize[0];x++)  
			{  
				BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
				Outputimg.SetPixelByGreyscale(x,y,z,0,0);
			}  
		}  
	}  
	float fzaxi[3]={0,0,1};
	float fyaxi[3]={0,1,0};
	float fxaxi[3]={1,0,0};
	//visit each node
	vector<vector<PointImgTypeDef>> vv_proponts;
	for (int i=0;i<0.05*ref.size();i++)
	{

		jdq2017::point3D curPontW=ref[i];
		float fcurPontI[3]={curPontW._x,curPontW._y,curPontW._z};
		imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI);
		int ncurPontI[3];
		ncurPontI[0]= (int)(fcurPontI[0]+0.5);
		ncurPontI[1]= (int)(fcurPontI[1]+0.5);
		ncurPontI[2]= (int)(fcurPontI[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
		short smar=BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]];
		if( smar==0)
		{
			//set the intensity of the projected image
			float fdivec[3]={0,0,0};
			Calc_divec(i,fdivec,ref);
			vector<PointCordTypeDef> vscaponts;
			//generate the scatter points in a plane
			Gen_sca_ponts(vscaponts);
			//center point
			float fcenpontW[3]={ncurPontI[0],ncurPontI[1],ncurPontI[2]};
			imgReadRaws.GetImageInfo()->ImageToWorld(fcenpontW);
			Matrix<float,1,3>MPro_Pcent;
			MPro_Pcent<<fcenpontW[0],fcenpontW[1],fcenpontW[2];
			vector<PointImgTypeDef> vsc_img;
			vector<PointCordTypeDef> vsc_rot;
			Matrix<float,3,3> Mrx1;
			Matrix<float,3,3> Mrx2;
			Rot_ponts_vec_zaxis(vscaponts,MPro_Pcent,fdivec,vsc_rot,Mrx1,Mrx2);
			//transform to the img cor
			Get_int_img(imgReadRaws,vsc_rot);
			vector<PointCordTypeDef> vsc_rotbac_int;
			//rot back to the origin
			RotBack(imgReadRaws,vsc_rot,MPro_Pcent,Mrx1,Mrx2,vsc_rotbac_int);
			//trans to img
			Tans_img(imgReadRaws,vsc_rotbac_int,vsc_img,fmovvec);
			vv_proponts.push_back(vsc_img);
			BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]]=1;
		}

	}
	for(int h=0;h<vv_proponts.size();h++)
	{
		vector<PointImgTypeDef> vtmp;
		vtmp=vv_proponts[h];
		for(int k=0;k<vtmp.size();k++)
		{
			PointImgTypeDef tmpp;
			tmpp=vtmp[k];
			Outputimg.SetPixelByGreyscale(tmpp.x,tmpp.y,h,0,tmpp.val);
		}
	}
	delete[] BW;


	return true;
}
int Find_corr(int k,const std::vector<jdq2017::Connection> &connections)
{
	int nmet=0;
	for ( std::vector<jdq2017::Connection>::const_iterator i = connections.begin();
		i != connections.end();
		++i )
	{
		int nref=i->first;
		if(k==nref)
		{
			nmet=i->second;
			break;
		}
	}

	return nmet;
}
bool Prj_to_plane(jdq2017::point3D curPontM,float fdivec[3],Matrix<float,1,3> MPro_Pcent,jdq2017::point3D &newPontM)
{
	//ax+by+cz+D=0
	float fxyz[3]={curPontM._x,curPontM._y,curPontM._z};
	float fD=-fdivec[0]*fxyz[0]-fdivec[1]*fxyz[1]-fdivec[2]*fxyz[2];
	float fabc[4]={fdivec[0],fdivec[1],fdivec[2],fD};
		float ft=(-fabc[0]*fxyz[0]-fabc[1]*fxyz[1]-fabc[2]*fxyz[2]-fabc[3])/(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
		newPontM._x=ft*fabc[0]+fxyz[0];
		newPontM._y=ft*fabc[1]+fxyz[1];
		newPontM._z=ft*fabc[2]+fxyz[2];
	return true;
}
bool Prj_to_plane_new(float fcurP[3],float fdivec[3],float fcenP[3],float fnewP[3])
{
	//ax+by+cz+D=0
	float fxyz[3]={fcurP[0],fcurP[1],fcurP[2]};
	float fD=-fdivec[0]*fcenP[0]-fdivec[1]*fcenP[1]-fdivec[2]*fcenP[2];
	float fabc[4]={fdivec[0],fdivec[1],fdivec[2],fD};
		float ft=(-fabc[0]*fxyz[0]-fabc[1]*fxyz[1]-fabc[2]*fxyz[2]-fabc[3])/(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
		fnewP[0]=ft*fabc[0]+fxyz[0];
		fnewP[1]=ft*fabc[1]+fxyz[1];
		fnewP[2]=ft*fabc[2]+fxyz[2];
	return true;
}
bool CMP_inWin_corr(zxhImageDataT<short>&imgReadRaws,std::vector<jdq2017::point3D> &ref,std::vector<jdq2017::point3D> &met,const std::vector<jdq2017::Connection> &connections,zxhImageDataT<short>&Outputimg,zxhImageDataT<short>&Outputimg_met)
{
	
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	//calculate the distance from the center
	float norigin[3]={0,0,0};
	float nimgcen[3]={0.5*nimgSize[0],0.5*nimgSize[1],0};
	imgReadRaws.GetImageInfo()->ImageToWorld(norigin);
	imgReadRaws.GetImageInfo()->ImageToWorld(nimgcen);
	float fmovvec[3];
	zxh::VectorOP_Substract(nimgcen,norigin,fmovvec,3);
	
	//built mark
	short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
	for(int z=0;z<nimgSize[2];z++)  
	{  
		for (int y=0;y<nimgSize[1];y++)  
		{  
			for (int x=0;x<nimgSize[0];x++)  
			{  
				BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
				Outputimg.SetPixelByGreyscale(x,y,z,0,0);
			}  
		}  
	}  
	float fzaxi[3]={0,0,1};
	float fyaxi[3]={0,1,0};
	float fxaxi[3]={1,0,0};
	//visit each node
	vector<vector<PointImgTypeDef>> vv_proponts;
	vector<PointImgTypeDef> v_metpoints;
	for (int i=0;i<0.8*ref.size();i++)
	{
		
		jdq2017::point3D curPontW=ref[i];
		float fcurPontI[3]={curPontW._x,curPontW._y,curPontW._z};
		imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI);
		int ncurPontI[3];
		ncurPontI[0]= (int)(fcurPontI[0]+0.5);
		ncurPontI[1]= (int)(fcurPontI[1]+0.5);
		ncurPontI[2]= (int)(fcurPontI[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
		short smar=BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]];
		if( smar==0)
		{
			//find the connected points in method line
		int nmet=Find_corr(i,connections);
			jdq2017::point3D curPontM=met[nmet];
			float fcurPontI_met[3]={curPontM._x,curPontM._y,curPontM._z};
		imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI_met);
		
		int ncurPontI_met[3];
		ncurPontI_met[0]= (int)(fcurPontI_met[0]+0.5);
		ncurPontI_met[1]= (int)(fcurPontI_met[1]+0.5);
		ncurPontI_met[2]= (int)(fcurPontI_met[2]+0.5);
		BoundaryCorrect(ncurPontI_met,nimgSize);
			//set the intensity of the projected image
			float fdivec[3]={0,0,0};
			Calc_divec(i,fdivec,ref);
			vector<PointCordTypeDef> vscaponts;
			//generate the scatter points in a plane
			Gen_sca_ponts(vscaponts);
			//center point
			float fcenpontW[3]={ncurPontI[0],ncurPontI[1],ncurPontI[2]};
			imgReadRaws.GetImageInfo()->ImageToWorld(fcenpontW);
			Matrix<float,1,3>MPro_Pcent;
			MPro_Pcent<<fcenpontW[0],fcenpontW[1],fcenpontW[2];
			vector<PointImgTypeDef> vsc_img;
			vector<PointCordTypeDef> vsc_rot;
			Matrix<float,3,3> Mrx1;
			Matrix<float,3,3> Mrx2;
			Rot_ponts_vec_zaxis(vscaponts,MPro_Pcent,fdivec,vsc_rot,Mrx1,Mrx2);
			//transform to the img cor
			Get_int_img(imgReadRaws,vsc_rot);
			//project the corresponding method point to the plane
			jdq2017::point3D newPontM,newPontM_rob;
			Prj_to_plane(curPontM,fdivec,MPro_Pcent,newPontM);
			vector<PointCordTypeDef> vsc_rotbac_int;
			//rot back to the origin
			RotBack_corr(vsc_rot,MPro_Pcent,Mrx1,Mrx2,vsc_rotbac_int,newPontM,newPontM_rob);
			//trans to img
			Tans_img(imgReadRaws,vsc_rotbac_int,vsc_img,fmovvec);
			//trans the method point to img
			PointImgTypeDef newPontM_img;
			Tran_metp_img(imgReadRaws,newPontM_rob,newPontM_img,fmovvec);
			//store the slice points
			vv_proponts.push_back(vsc_img);
			v_metpoints.push_back(newPontM_img);
			BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]]=1;
		}

	}
	for(int h=0;h<vv_proponts.size();h++)
	{
		vector<PointImgTypeDef> vtmp;
		vtmp=vv_proponts[h];
		for(int k=0;k<vtmp.size();k++)
		{
			PointImgTypeDef tmpp;
			tmpp=vtmp[k];
			Outputimg.SetPixelByGreyscale(tmpp.x,tmpp.y,h,0,tmpp.val);
		}
	}
	for(int k=0;k<v_metpoints.size();k++)
		{
			PointImgTypeDef tmpp;
			tmpp=v_metpoints[k];
			Outputimg_met.SetPixelByGreyscale(tmpp.x,tmpp.y,k,0,1000);
		}
	delete[] BW;


	return true;
}
bool Get_dirvec_byPrj(float fdivec[3],jdq2017::point3D curPontW,float fdivec1[3])
{
	float fcenP[3]={curPontW._x,curPontW._y,curPontW._z};
	float fcurP[3]={1,0,0};
	float fori[3]={0,0,0};
	float fnewP[3]={0,0,0};
	Prj_to_plane_new( fcurP, fdivec, fcenP, fnewP);
		float fnewori[3]={0,0,0};
	Prj_to_plane_new( fori, fdivec, fcenP, fnewori);
	fdivec1[0]=fnewP[0]-fnewori[0];
	fdivec1[1]=fnewP[1]-fnewori[1];
	fdivec1[2]=fnewP[2]-fnewori[2];
	return true;
}
bool Gen_sca_ponts_inPla(zxhImageDataT<short>&imgReadRaws,float fdivec[3],float fdivec1[3],float fdivec2[3],jdq2017::point3D curPontW,vector<PointCordTypeDef> &vscaponts,Matrix<float,3,3> &Mrx1,
			Matrix<float,3,3> &Mrx2)
{
	float res=0.2;
	int npnum=50;
	for(int i=-1*npnum;i<=npnum;i++)
		for(int j=-1*npnum;j<=npnum;j++)
		{
			float fx[3]={res*i*fdivec1[0],res*i*fdivec1[1],res*i*fdivec1[2]};
			float fy[3]={res*j*fdivec2[0],res*j*fdivec2[1],res*j*fdivec2[2]};
			float fv[3]={0,0,0};
			zxh::VectorOP_Substract(fy,fx,fv,3);
			PointCordTypeDef temp;
			temp.x=curPontW._x+fv[0];
			temp.y=curPontW._y+fv[1];
			temp.z=curPontW._z+fv[2];
			vscaponts.push_back(temp);

	
	/*	MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<temp.x,temp.y,temp.z;
			Matrix<float,1,3>MPro_Pcent;
	MPro_Pcent<<curPontW._x,curPontW._y,curPontW._z;
		MnewPont=MoriPont-MPro_Pcent;
		PointCordTypeDef PnewPont;
		PnewPont.x=(MnewPont*Mrx1*Mrx2)(0,0);
		PnewPont.y=(MnewPont*Mrx1*Mrx2)(0,1);
		PnewPont.z=(MnewPont*Mrx1*Mrx2)(0,2);
			float fcosi=zxh::VectorOP_Cosine(fv,fdivec,3);
				PointCordTypeDef curPontW;
	curPontW.x=PnewPont.x;
		curPontW.y=PnewPont.y;
		curPontW.z=PnewPont.z;
		PointImgTypeDef curPontI;
		WorldtoImage(imgReadRaws,curPontW,curPontI);

			if(curPontI.x==13&&curPontI.y==14)
				int xxx=3;
			int xxx=0;*/

		}

			
		
	return true;
}
bool Get_int(zxhImageModelingLinear &mod,zxhImageDataT<short>&imgReadRaws,vector<PointCordTypeDef>&vscaponts)
{
	//built mark

	mod.SetImage(&imgReadRaws);
	for (int i=0;i<vscaponts.size();i++)
	{
		PointCordTypeDef curPontW=vscaponts[i];
		short inten = mod.GetPixelFloatValueWithCheckByWorld(curPontW.x,curPontW.y,curPontW.z) ; 
		vscaponts[i].val=inten;
	}
	return true;
}
bool Rot_ponts_to_zaxis(zxhImageDataT<short>&imgReadRaws,vector<PointCordTypeDef>&vscaponts,jdq2017::point3D curPontW,vector<PointCordTypeDef>&vscaponts_rot,Matrix<float,3,3> &Mrx1,
			Matrix<float,3,3> &Mrx2,float fmovvec[3])
{
	
	Matrix<float,1,3>MPro_Pcent;
	MPro_Pcent<<curPontW._x,curPontW._y,curPontW._z;
	//rotate the points to the local axis
	Rotpon_With_Matr12(imgReadRaws,vscaponts,MPro_Pcent,Mrx1,Mrx2,vscaponts_rot, fmovvec);

	return true;
}

bool Get_ponts_inPlan(int i,std::vector<jdq2017::point3D> &ref,jdq2017::point3D &curPontW,zxhImageDataT<short>&imgReadRaws,float fmovvec[3],float fdivec[3],Matrix<float,3,3> &Mrx1,Matrix<float,3,3>&Mrx2,vector<vector<PointImgTypeDef>> &vv_proponts)
{//get the tangent direction
			
			Calc_divec(i,fdivec,ref);
			//get the second direction of the plane
			float fdivec1[3]={0,0,0};
			Get_dirvec_byPrj(fdivec,curPontW,fdivec1);
			//get the third direction of the plane
			float fdivec2[3]={0,0,0};
			zxh::VectorOP_CrossProduct3D(fdivec,fdivec1,fdivec2);
			zxh::VectorOP_Normalise(fdivec2,3);
			zxh::VectorOP_Normalise(fdivec1,3);
		
        	//rotate the vector to z-axis
	         Rotaz_to_vec(fdivec,Mrx1,Mrx2);	
			//generate the scatter points in this local  plane
			vector<PointCordTypeDef> vscaponts,vscaponts_rot;
			Gen_sca_ponts_inPla(imgReadRaws,fdivec,fdivec1,fdivec2,curPontW,vscaponts,Mrx1,Mrx2);
			//get the intensity of each points
			zxhImageModelingLinear mod ;
			vector<PointImgTypeDef>vscapontI;
			Get_int(mod,imgReadRaws,vscaponts);
			//rotate the points to the world cor
			Rot_ponts_to_zaxis(imgReadRaws,vscaponts,curPontW,vscaponts_rot,Mrx1,Mrx2,fmovvec);
			//trans to the image points
			Tans_img(imgReadRaws,vscaponts_rot,vscapontI,fmovvec);
			//store the point
			//store the slice points
			vv_proponts.push_back(vscapontI);

	return true;
}
bool Rot_to_z_axis(float fnewP[3],jdq2017::point3D &curPontM,Matrix<float,1,3>MPro_Pcent,Matrix<float,3,3> &Mrx1,Matrix<float,3,3> &Mrx2,PointCordTypeDef &PnewPont)
{
		PointCordTypeDef curPontW;
		MatrixXf MoriPont(1,3),MnewPont(1,3);
		MoriPont<<curPontM._x,curPontM._y,curPontM._z;
		MnewPont=MoriPont-MPro_Pcent;
		PnewPont.x=(MnewPont*Mrx1*Mrx2)(0,0);
		PnewPont.y=(MnewPont*Mrx1*Mrx2)(0,1);
		PnewPont.z=(MnewPont*Mrx1*Mrx2)(0,2);
		PnewPont.val=1;
return true;
}
bool Tran_Cord_img(zxhImageDataT<short>&imgReadRaws,float fmovvec[3],PointCordTypeDef PnewPont,PointImgTypeDef &PnewImgPont)
{

	PointCordTypeDef curPontW;
	    curPontW.x=PnewPont.x+fmovvec[0];
		curPontW.y=PnewPont.y+fmovvec[1];
		curPontW.z=PnewPont.z+fmovvec[2];
		WorldtoImage(imgReadRaws,curPontW,PnewImgPont);	
	return true;
}
bool Get_met_ponts(jdq2017::point3D&curPontM,jdq2017::point3D &curPontW,zxhImageDataT<short>&imgReadRaws,float fmovvec[3],float fdivec[3],Matrix<float,3,3> &Mrx1,Matrix<float,3,3>&Mrx2,std::vector<jdq2017::point3D>&v_Cmetpoints,vector<PointImgTypeDef>&v_metpoints)
{
	//prject the method point to the norm plane
			float fcenP[3]={curPontW._x,curPontW._y,curPontW._z};
	float fcurP[3]={curPontM._x,curPontM._y,curPontM._z};
	float fnewP[3]={0,0,0};
	Prj_to_plane_new( fcurP, fdivec, fcenP, fnewP);
	//rotate the new projected point to the z-axis
		PointCordTypeDef PnewPont;
		Matrix<float,1,3>MPro_Pcent;
	MPro_Pcent<<curPontW._x,curPontW._y,curPontW._z;
	Rot_to_z_axis(fnewP,curPontM,MPro_Pcent,Mrx1,Mrx2,PnewPont);
	jdq2017::point3D PnewjdqPont;
 PnewjdqPont._x=PnewPont.x;
 PnewjdqPont._y=PnewPont.y;
  PnewjdqPont._y=PnewPont.y;
	v_Cmetpoints.push_back(PnewjdqPont);
	////trans Cord to Img
	PointImgTypeDef PnewImgPont;
	Tran_Cord_img(imgReadRaws,fmovvec,PnewPont,PnewImgPont);
	v_metpoints.push_back(PnewImgPont);

	return true;
}
bool Smooth(std::vector<jdq2017::point3D> v_CmetPonts,std::vector<jdq2017::point3D>&v_CmetPonts_SMO)
{

	float BackTrackDistmm=pathLength(v_CmetPonts);//calculate the total length
	float m_meandis_of_coronaryvessel=BackTrackDistmm/(v_CmetPonts.size() - 1);

	float fStepDist=BackTrackDistmm/100;
	
	if (v_CmetPonts_SMO.size() > 0)
	{
		v_CmetPonts_SMO.clear();
	}

	v_CmetPonts_SMO.push_back(v_CmetPonts[0]);

	int idistNum_2mm=2/m_meandis_of_coronaryvessel;
	for(int i=1;i<101;i++)//100 seg point
	{
		int SegNUM=0;
		float dist=0;
		for(int j=0;j<v_CmetPonts.size()-1;j++)
		{
			 double len = (v_CmetPonts[j]-v_CmetPonts[j+1]).length();
			dist=dist+len;
			if(dist>=i*fStepDist)
			{
				SegNUM=j;
				break;
			}
		}

		int SegStart=SegNUM-idistNum_2mm;
		int SegEnd=SegNUM+idistNum_2mm;
		if (SegStart<=0)SegStart=0;
		if (SegEnd>=v_CmetPonts.size()-1)SegEnd=v_CmetPonts.size()-1;
		float fSegSum[3]={0};
		int n=0;
		for (int k=SegStart;k<=SegEnd;k++)
		{
			fSegSum[0]=fSegSum[0]+v_CmetPonts[k]._x;
			fSegSum[1]=fSegSum[1]+v_CmetPonts[k]._y;
			fSegSum[2]=fSegSum[2]+v_CmetPonts[k]._z;
			n++;
		}
		jdq2017::point3D mtempPointWorld;
		mtempPointWorld._x=fSegSum[0]/n;
		mtempPointWorld._y=fSegSum[1]/n;
		mtempPointWorld._z=fSegSum[2]/n;
		v_CmetPonts_SMO.push_back(mtempPointWorld);
	}

	return true;

}
bool Tans_img_new(zxhImageDataT<short>&imgReadRaws,std::vector<jdq2017::point3D>&vsc_rotbac_int,vector<PointImgTypeDef>&vsc_img)
{
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
	for(int z=0;z<nimgSize[2];z++)  
	{  
		for (int y=0;y<nimgSize[1];y++)  
		{  
			for (int x=0;x<nimgSize[0];x++)  
			{  
				BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
			}  
		}  
	} 

	for (int i=0;i<vsc_rotbac_int.size();i++)
	{
		PointCordTypeDef curPontW;
	    curPontW.x=vsc_rotbac_int[i]._x;
		curPontW.y=vsc_rotbac_int[i]._y;
		curPontW.z=vsc_rotbac_int[i]._z;
		PointImgTypeDef curPontI;
		WorldtoImage(imgReadRaws,curPontW,curPontI);
		short smar=BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x];
		if(smar==0)
		{
			vsc_img.push_back(curPontI);
			BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x]=1;
		}

	}
	delete[]BW;
	return true;
}
bool resampleL(std::vector<jdq2017::point3D> ref,std::vector<jdq2017::point3D>&met)
{
	double refSampling = pathLength(ref)/(ref.size()-1);
	if (debugOutput)
	{
		std::cout << "Reference sampling distance = " << refSampling << std::endl;
	}
	// Calculate method sampling distance
	double clSampling = pathLength(met)/(met.size()-1);
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
		resampler(met, refSampling);
		met = resampler.resultPath();
	}
	return true;
}
bool CMP_inWin_corr_new(zxhImageDataT<short>&imgReadRaws,std::vector<jdq2017::point3D> &ref,std::vector<jdq2017::point3D> &met,const std::vector<jdq2017::Connection> &connections,vector<vector<PointImgTypeDef>>&vv_proponts,vector<PointImgTypeDef>&v_metpoints,std::vector<jdq2017::point3D>&v_CmetPonts)
{
	
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	//calculate the distance from the center
	float norigin[3]={0,0,0};
	float nimgcen[3]={0.5*nimgSize[0],0.5*nimgSize[1],0};
	imgReadRaws.GetImageInfo()->ImageToWorld(norigin);
	imgReadRaws.GetImageInfo()->ImageToWorld(nimgcen);
	float fmovvec[3];
	zxh::VectorOP_Substract(nimgcen,norigin,fmovvec,3);
	
	//built mark
	short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
	for(int z=0;z<nimgSize[2];z++)  
	{  
		for (int y=0;y<nimgSize[1];y++)  
		{  
			for (int x=0;x<nimgSize[0];x++)  
			{  
				BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
			}  
		}  
	}  
	float fzaxi[3]={0,0,1};
	float fyaxi[3]={0,1,0};
	float fxaxi[3]={1,0,0};
	//visit each node

	for (int i=0;i<ref.size();i++)
	{
		
		jdq2017::point3D curPontW=ref[i];
		float fcurPontI[3]={curPontW._x,curPontW._y,curPontW._z};
		imgReadRaws.GetImageInfo()->WorldToImage(fcurPontI);
		int ncurPontI[3];
		ncurPontI[0]= (int)(fcurPontI[0]+0.5);
		ncurPontI[1]= (int)(fcurPontI[1]+0.5);
		ncurPontI[2]= (int)(fcurPontI[2]+0.5);
		BoundaryCorrect(ncurPontI,nimgSize);
		short smar=BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]];
		if( smar==0)
		{
			//find the connected points in method line
			int nmet=Find_corr(i,connections);
			jdq2017::point3D curPontM=met[nmet];
			//Get the points in norm plane
			float fdivec[3]={0,0,0};
				Matrix<float,3,3> Mrx1;
			Matrix<float,3,3> Mrx2;
			Get_ponts_inPlan(i,ref,curPontW,imgReadRaws,fmovvec,fdivec,Mrx1,Mrx2,vv_proponts);
			//Get the method points
			Get_met_ponts(curPontM,curPontW,imgReadRaws,fmovvec,fdivec,Mrx1,Mrx2,v_CmetPonts,v_metpoints);
			//mark as visited
			BW[ncurPontI[2]* nimgSize[1] * nimgSize[0] + ncurPontI[1]* nimgSize[0] + ncurPontI[0]]=1;
		}

	}

	//for(int h=0;h<vv_proponts.size();h++)
	//{
	//	vector<PointImgTypeDef> vtmp;
	//	vtmp=vv_proponts[h];
	//	for(int k=0;k<vtmp.size();k++)
	//	{
	//		PointImgTypeDef tmpp;
	//		tmpp=vtmp[k];
	//		//Outputimg.SetPixelByGreyscale(tmpp.x,tmpp.y,h,0,tmpp.val);
	//	}
	//}
	//for(int k=0;k<v_metpoints.size();k++)
	//	{
	//		PointImgTypeDef tmpp, tmpCMP;
	//		tmpp=v_metpoints[k];
	//		//Outputimg_met.SetPixelByGreyscale(tmpp.x,tmpp.y,k,0,1000);
	//		tmpCMP.x=tmpp.x;
	//		tmpCMP.y=tmpp.y;
	//		tmpCMP.z=k;
	//		PointCordTypeDef tmpCord;
	//		ImagetoWorld(Outputimg_met,tmpp,tmpCord);
	//		jdq2017::point3D tmpjdqCord;
	//		tmpjdqCord._x=tmpCord.x;
	//		tmpjdqCord._y=tmpCord.y;
	//		tmpjdqCord._z=k;
	//		//v_CmetPonts_CMP.push_back(tmpjdqCord);
	//	}

		/*int nLen3 = strlen(chResultFilePath)+ strlen("CutLine.vtk") + 1;
		char *chFileName3 = (char *)malloc(nLen3);
		strcpy(chFileName3, chResultFilePath);
		strcat(chFileName3, "CutLine.vtk");
			Write2Vtk_jdq(v_CmetPonts_CMP,chFileName3);
		free(chFileName0);

			int nLen1= strlen(chResamfilename)+ strlen("CutLine.txt") + 1;
		char *chFileName1 = (char *)malloc(nLen1);
		strcpy(chFileName1, chResamfilename);
		strcat(chFileName1, ".txt");
			Write2Txt_jdq(v_CmetPonts_CMP,chFileName1);
		free(chFileName0);*/

	return true;
}
int Calc_num_imgpont(zxhImageDataT<short>&imgReadRaws,std::vector<jdq2017::point3D> ref)
{
int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	short *BW = new short[nimgSize[2] * nimgSize[1] * nimgSize[0]]; 
	for(int z=0;z<nimgSize[2];z++)  
	{  
		for (int y=0;y<nimgSize[1];y++)  
		{  
			for (int x=0;x<nimgSize[0];x++)  
			{  
				BW[z* nimgSize[1] * nimgSize[0] + y* nimgSize[0] + x]=0;  
			}  
		}  
	} 
	int num=0;
	for (int i=0;i<ref.size();i++)
	{
		PointCordTypeDef curPontW;
	    curPontW.x=ref[i]._x;
		curPontW.y=ref[i]._y;
		curPontW.z=ref[i]._z;
		PointImgTypeDef curPontI;
		WorldtoImage(imgReadRaws,curPontW,curPontI);
		short smar=BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x];
		if(smar==0)
		{
			num++;
			BW[curPontI.z* nimgSize[1] * nimgSize[0] + curPontI.y* nimgSize[0]+curPontI.x]=1;
		}

	}
	delete[]BW;
	return num;
}
bool Set_Ponts(zxhImageDataT<short>&imgReadRaws,vector<vector<PointImgTypeDef>>&vv_proponts,vector<PointImgTypeDef>&v_metpoints,std::vector<jdq2017::point3D>&v_CmetPonts_CMP,zxhImageDataT<short>&Outputimg,zxhImageDataT<short>&Outputimg_met)
{
	init_img(Outputimg);

	for(int h=0;h<vv_proponts.size();h++)
	{
		vector<PointImgTypeDef> vtmp;
		vtmp=vv_proponts[h];
		for(int k=0;k<vtmp.size();k++)
		{
			PointImgTypeDef tmpp;
			tmpp=vtmp[k];
			Outputimg.SetPixelByGreyscale(tmpp.x,tmpp.y,h,0,tmpp.val);
		}
	}
	init_img(Outputimg_met);
	for(int k=0;k<v_metpoints.size();k++)
		{
			PointImgTypeDef tmpp, tmpCMP;
			tmpp=v_metpoints[k];
			Outputimg_met.SetPixelByGreyscale(tmpp.x,tmpp.y,k,0,1000);
			PointImgTypeDef curPontI;
			curPontI.x=tmpp.x;
			curPontI.y=tmpp.y;
			curPontI.z=k;
			PointCordTypeDef curPointW;
			ImagetoWorld(Outputimg_met,curPontI,curPointW);
			jdq2017::point3D jdqpnt;
			jdqpnt._x=curPointW.x;
				jdqpnt._y=curPointW.y;
					jdqpnt._z=curPointW.z;
			v_CmetPonts_CMP.push_back(jdqpnt);
		}
	return true;
}

int main(int argc, char *argv[])
{
	//if( argc < 4 )
	//{
	//	cerr << "Usage: " << endl;
	//	cerr << "zxhcaeDMPMToNewImg.cpp	imageRaw(.nii)	imageResoRaw(.nii) MCLine(.vtk) ResultHigh-ResolutionName(.nii.gz) ResultLow-ResolutionName(.nii.gz)" << endl;
	//	return -1;
	//}
	//string strFileNameRaw =string(argv[1]);//  "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	//char *chFileName =argv[2];//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
	//char *chResultName =argv[3];// "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/MCLine.nii.gz";
	//string RorN=string(argv[4]);

	//string strFileNameRaw ="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	//char *chFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
	//char *chResultName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/MCLine.nii.gz";
 //   string RorN="N";

	//string strFileNameRaw ="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	//char *chFileName ="F:/Coronary_0/Exp6_detect_ostia/Ostia0007To00007/result_dataset00_rightostium.txt";
	//char *chResultName ="F:/Coronary_0/Exp6_detect_ostia/MapModelToImg/MP_img00.nii.gz";
 //   string RorN="-R";
//***Training Data as model****//
	
	string strFileNameRaw =  "E:/work_jdq/rebuttal_media18/rebuttal/Exp/LowRes_0307/dataset03/image03_resave.nii.gz";
	char *chRefCurvefilename = "E:/work_jdq/rebuttal_media18/rebuttal/Exp/LowRes_0307/dataset03/vessel1/reference.vtk";
	char *chMetCurvefilename = "E:/work_jdq/rebuttal_media18/rebuttal/Exp/LowRes_0307/dataset03/vessel1/SMD_resm_SMO.txt";
	char *chResultFilePath = "E:/work_jdq/rebuttal_media18/rebuttal/Exp/LowRes_0307/dataset03/vessel1/";
	std::vector<jdq2017::point3D> ref,met;
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
		
		std::cerr << "Type of Ref line is .vtk" << std::endl;
	}
	if(strcmp(chext,".txt")==0)// read vtk file
	{
		if (! jdq2017::readCenterline(chRefCurvefilename, ref) )
		{


			return 0;
		}
		txt2vtk(ref);
		std::cerr << "Type of Ref line is .txt" << std::endl;
	}
	char* chext_tar=GetPathEXT(chMetCurvefilename);
	if(strcmp(chext_tar,".vtk")==0)// read vtk file
	{
		if (! jdq2017::readCenterlinevtk(chMetCurvefilename, met) )
		{
			if (!outputOneLine) 
			{
				std::cerr << "Error in reading input data" << std::endl;
				return 1;
			} 
		}
		std::cerr << "Type of met line is .vtk" << std::endl;
	}
	if(strcmp(chext_tar,".txt")==0)// read vtk file
	{
		if (! jdq2017::readCenterline(chMetCurvefilename, met) )
		{


			return 0;
		}
		std::cerr << "Type of met line is .txt" << std::endl;
		txt2vtk(met);
	}
	// ** Resample centerline to same sampling as reference standard
  // Calculate reference sampling distance
  double refSampling = pathLength(ref)/(ref.size()-1);
  if (debugOutput)
  {
    std::cout << "Reference sampling distance = " << refSampling << std::endl;
  }
  // Calculate method sampling distance
  double clSampling = pathLength(met)/(met.size()-1);
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
    resampler(met, refSampling);
    met = resampler.resultPath();

  }

  

  //The input is correct
  //*******
  
  // Determine correspondence between reference and method centerline
  jdq2017::Dijkstra<jdq2017::point3D> correspondence;
 correspondence(ref, met);
  //Get correspondences as a vector of connections
  //The connections are stored as pairs
  //Connection.first is an index in the reference centerline
  //Connection.second is an index in the method centerline
  const std::vector<jdq2017::Connection> &connections(correspondence.pairing());


	zxhImageDataT<short> imgReadRaws;//Change by JDQ
	if( zxh::OpenImage( &imgReadRaws, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"; 
		return -1;
	}
	//Cutimg_byLines_2D(imgReadRaws,ref);
	//Loca_pro_2D(imgReadRaws,ref);
	int nimgSize[4];
	imgReadRaws.GetImageSize(nimgSize[0], nimgSize[1], nimgSize[2], nimgSize[3]);
	zxhImageDataT<short>Outputimg,Outputimg_met,Outputimg_met_SMO;
	//Outputimg.NewImage(imgReadRaws.GetImageInfo());
	int num=Calc_num_imgpont(imgReadRaws,ref);
	int newsize[4]={nimgSize[0],nimgSize[1],num,1};
	float newspacing[3]={0.5,0.5,0.5};
	Outputimg.NewImage(3,newsize,newspacing,Outputimg.GetImageInfo());
	Outputimg_met.NewImage(3,newsize,newspacing,Outputimg.GetImageInfo());
	Outputimg_met_SMO.NewImage(3,newsize,newspacing,Outputimg.GetImageInfo());
	//CMP_inWin(imgReadRaws,ref,Outputimg);
		vector<vector<PointImgTypeDef>> vv_proponts;
	vector<PointImgTypeDef> v_metpoints;
	std::vector<jdq2017::point3D> v_CmetPonts;
	CMP_inWin_corr_new(imgReadRaws,ref,met,connections,vv_proponts,v_metpoints,v_CmetPonts);
	//
	std::vector<jdq2017::point3D>v_CmetPonts_CMP;
	Set_Ponts(imgReadRaws,vv_proponts,v_metpoints,v_CmetPonts_CMP,Outputimg,Outputimg_met);
	//the new image has a different 
	//save the results
	int nLen0 = strlen(chResultFilePath)+ strlen("img_2Dcut.nii.gz") + 1;
		char *chFileName0 = (char *)malloc(nLen0);
		strcpy(chFileName0, chResultFilePath);
		strcat(chFileName0, "img_2Dcut.nii.gz");	
		string strFileName0(chFileName0);
	zxh::SaveImage(&Outputimg,strFileName0.c_str());
		free(chFileName0);

		int nLen1 = strlen(chResultFilePath)+ strlen("img_2Dcut_met.nii.gz") + 1;
		char *chFileName1 = (char *)malloc(nLen1);
		strcpy(chFileName1, chResultFilePath);
		strcat(chFileName1, "img_2Dcut_met.nii.gz");	
		string strFileName1(chFileName1);
	zxh::SaveImage(&Outputimg_met,strFileName1.c_str());
		free(chFileName1);

		int nLen2 = strlen(chResultFilePath)+ strlen("img_2Dcut_met_C.vtk") + 1;
		char *chFileName2 = (char *)malloc(nLen2);
		strcpy(chFileName2, chResultFilePath);
		strcat(chFileName2, "img_2Dcut_met_C.vtk");	
	Write2Vtk_jdq(v_CmetPonts_CMP,chFileName2);
		free(chFileName2);

			int nLen3 = strlen(chResultFilePath)+ strlen("img_2Dcut_met_C.txt") + 1;
		char *chFileName3 = (char *)malloc(nLen3);
		strcpy(chFileName3, chResultFilePath);
		strcat(chFileName3, "img_2Dcut_met_C.txt");	
	Write2Txt_jdq(v_CmetPonts_CMP,chFileName3);
		free(chFileName3);
		
	cout << "Line has already mapped to image!" << endl;


}

 