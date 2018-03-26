

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



#include "miiMinPathModel.h"

#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkPolyLine.h"
#include "miiMinPath.h"
#include <time.h>

#define TAB_CHAR	9

int GetTime() ;
using namespace std;
typedef struct
{
	double x;
	double y;
	double z;
}PointCordTypeDef;

bool ReadMeanStdDevFromTxt(char *chFileName, double &fMean, double &fStdDev)
{
	string strNum;
	int nStart = 0, nEnd = 0;

	ifstream iFileIn(chFileName);
	if(!iFileIn)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << chFileName << endl;
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
	fMean = atof(strNum.c_str());

	// read 2nd number
	nStart = nEnd + 1;
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	fStdDev = atof(strNum.c_str());


	return true;
}
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

void WriteVtk(vector< miiCNode<double, float> > vPointCord, char *chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

/*	int nPointNum = PointCord.size();*/

	int nPointNum = vPointCord.size();

	for (int i = 0; i < nPointNum; i++)
	{
		iPoints->InsertNextPoint(vPointCord[i].x, vPointCord[i].y, vPointCord[i].z);
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

	
	/*void WriteCA2Txt(vector< miiCNode<double, float> > vPointCord, char *chFileName)
	{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = vPointCord.size();

	for (int i = 0; i < nPointNum; i++)
	{
		
     WriteFileTxt << vPointCord[i].x << " " << vPointCord[i].y << " " << vPointCord[i].z << "/n ";

	}	
	
	}
	*/
	

void CalcUnseenMeanStdDev(double gAortaMean[2], double gAortaStdDev[2], double fAtlsCAMean, double fAtlsCAStdDev,
	double &fUnseenCAMean, double &fUnseenCAStdDev)
{
	double fAtlsAortaMean = gAortaMean[0];
	double fUnseenAortaMean = gAortaMean[1];
	

	double fMeanDiff = fAtlsAortaMean - fAtlsCAMean;
	fUnseenCAMean = fUnseenAortaMean - fMeanDiff;

	double fAtlsAortaStdDev = gAortaStdDev[0];
	double fUnseenAortaStdDev = gAortaStdDev[1];

	fUnseenCAStdDev = fAtlsCAStdDev * fUnseenAortaStdDev / (fAtlsAortaStdDev + 0.00000001);	
}
void ModifiedCalcUnseenMeanStdDev(double gAortaMean[2], double gAortaStdDev[2], double fAtlsCAMean,
	double &fUnseenCAMean, double &fUnseenCAStdDev)
{
	double fAtlsAortaMean = gAortaMean[0];
	double fUnseenAortaMean = gAortaMean[1];
	

	double fMeanDiff = fAtlsAortaMean - fAtlsCAMean;
	fUnseenCAMean = fUnseenAortaMean - fMeanDiff;

	double fAtlsAortaStdDev = gAortaStdDev[0];
	double fUnseenAortaStdDev = gAortaStdDev[1];

	fUnseenCAStdDev = fUnseenAortaStdDev;	
}

bool CalcAtlsCAMeanStdDev(string strFileNameAtls, char *strFileNameCA, double &fAtlsCAMean, double &fAtlsCAStdDev)//Calculate atlas meanstdev based on neighbour points
{
	//read a "nifti" file
	zxhImageDataT<short> imgReadAtls; 
	if( zxh::OpenImage( &imgReadAtls, strFileNameAtls ) == false )
	{
		std::cerr << " Raw image(nifti-file) is not found!"; 
		return false;
	}

	vector< miiCNode<double, float> > vAtlsCAPoints;
	ReadVtk(strFileNameCA, vAtlsCAPoints);

	// get the image size
	int nImgWX, nImgWY, nImgWZ, nImgWT;
	imgReadAtls.GetImageSize(nImgWX, nImgWY, nImgWZ, nImgWT);

	// get the image data
	const short *sImgData = imgReadAtls.GetImageData();

	int nSearchRange = 8;
	int nCount = 0;
	double nSum = 0;
	float fCord[3];
	fCord[0] = vAtlsCAPoints[0].x;
	fCord[1] = vAtlsCAPoints[0].y;
	fCord[2] = vAtlsCAPoints[0].z;
	imgReadAtls.GetImageInfo()->WorldToImage(fCord);
	int nXmin = (int)(fCord[0]+0.5) - nSearchRange;
	int nXmax = (int)(fCord[0]+0.5) + nSearchRange;

	int nYmin = (int)(fCord[1]+0.5) - nSearchRange;
	int nYmax = (int)(fCord[1]+0.5) + nSearchRange;

	int nZmin = (int)(fCord[2]+0.5) - nSearchRange;
	int nZmax = (int)(fCord[2]+0.5) + nSearchRange;
	
	// calculate mean
	for (int iz = nZmin; iz < nZmax; iz++)
	{
		for (int iy = nYmin; iy < nYmax; iy++)
		{
			for (int ix = nXmin; ix < nXmax; ix++)
			{
				if (sImgData[iz * nImgWY * nImgWX + iy * nImgWX + ix] > 10)
				{
					nCount++;
					nSum += (double)sImgData[iz * nImgWY * nImgWX + iy * nImgWX + ix];
				}
			}
		}
	}

	fAtlsCAMean = nSum / (nCount + 0.00000000001);
	nSum = 0;

	// calculate deviation
	for (int iz = nZmin; iz < nZmax; iz++)
	{
		for (int iy = nYmin; iy < nYmax; iy++)
		{
			for (int ix = nXmin; ix < nXmax; ix++)
			{
				if (sImgData[iz * nImgWY * nImgWX + iy * nImgWX + ix] > 10)
				{
					double temp = sImgData[iz * nImgWY * nImgWX + iy * nImgWX + ix] - fAtlsCAMean;
					temp *= temp;
					nSum += temp;
				}
			}
		}
	}

	fAtlsCAStdDev = sqrt(nSum / (nCount + 1));

	imgReadAtls.ReleaseMem();

	return true;
}
bool ModifiedCalcAtlsCAMeanStdDev(string strFileNameAtls, char *strFileNameCA, double &fAtlsCAMean)//Calculate atlas meanstdev based on one points
{
	//read a "nifti" file
	zxhImageDataT<short> imgReadAtls; 
	if( zxh::OpenImage( &imgReadAtls, strFileNameAtls ) == false )
	{
		std::cerr << " Raw image(nifti-file) is not found!"; 
		return false;
	}

	vector< miiCNode<double, float> > vAtlsCAPoints;
	ReadVtk(strFileNameCA, vAtlsCAPoints);

	// get the image size
	int nImgWX, nImgWY, nImgWZ, nImgWT;
	imgReadAtls.GetImageSize(nImgWX, nImgWY, nImgWZ, nImgWT);

	// get the image data
	const short *sImgData = imgReadAtls.GetImageData();

	
	

	float fCord[3];
	fCord[0] = vAtlsCAPoints[0].x;
	fCord[1] = vAtlsCAPoints[0].y;
	fCord[2] = vAtlsCAPoints[0].z;
	imgReadAtls.GetImageInfo()->WorldToImage(fCord);
	int iz=(int)(fCord[2]+0.5);
    int iy=(int)(fCord[1]+0.5);
	int ix=(int)(fCord[0]+0.5);
		
	fAtlsCAMean =(double)sImgData[iz * nImgWY * nImgWX + iy * nImgWX + ix];
	

	imgReadAtls.ReleaseMem();

	return true;
}
bool SetMolPointsIntensity(char *chModelMapCAFileName,string strFileNameAtls,vector< miiCNode<double, float> > &vModelPointsWithIntensity)
{
	
	zxhImageDataT<short> imgReadAtls; 
	if( zxh::OpenImage( &imgReadAtls, strFileNameAtls ) == false )
	{
		std::cerr << " Raw image(nifti-file) is not found!"; 
		return false;
	}

	
	ReadVtk(chModelMapCAFileName, vModelPointsWithIntensity);

	// get the image size
	int nImgWX, nImgWY, nImgWZ, nImgWT;
	imgReadAtls.GetImageSize(nImgWX, nImgWY, nImgWZ, nImgWT);

	// get the image data
	const short *sImgData = imgReadAtls.GetImageData();

	int nSearchRange = 8;
	for (int i=0;i<vModelPointsWithIntensity.size();i++)
	{
	int nCount = 0;
	double nSum = 0;
	float fCord[3];
	fCord[0] =vModelPointsWithIntensity[i].x;
	fCord[1] =vModelPointsWithIntensity[i].y;
	fCord[2] =vModelPointsWithIntensity[i].z;
	imgReadAtls.GetImageInfo()->WorldToImage(fCord);
	int nXmin = (int)(fCord[0]+0.5) - nSearchRange;
	int nXmax = (int)(fCord[0]+0.5) + nSearchRange;
	
	int nYmin = (int)(fCord[1]+0.5) - nSearchRange;
	int nYmax = (int)(fCord[1]+0.5) + nSearchRange;

	int nZmin = (int)(fCord[2]+0.5) - nSearchRange;
	int nZmax = (int)(fCord[2]+0.5) + nSearchRange;
	if( nXmin<0)nXmin=0;
	if( nXmax>nImgWX)nXmin=nImgWX;
	if( nYmin<0)nYmin=0;
	if( nYmax>nImgWY)nYmin=nImgWY;
	if( nZmin<0)nZmin=0;
	if( nZmax>nImgWZ)nZmin=nImgWZ;

	// calculate mean
	for (int iz = nZmin; iz < nZmax; iz++)
	{
		for (int iy = nYmin; iy < nYmax; iy++)
		{
			for (int ix = nXmin; ix < nXmax; ix++)
			{
				if (sImgData[iz * nImgWY * nImgWX + iy * nImgWX + ix] > 10)
				{
					nCount++;
					nSum += (double)sImgData[iz * nImgWY * nImgWX + iy * nImgWX + ix];
				}
			}
		}
	}

	vModelPointsWithIntensity[i].val= nSum / (nCount + 0.00000000001);
	nSum = 0;
		
	}
	imgReadAtls.ReleaseMem();

	return true;
}
bool VesselnessInit(const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, double *m_sVes,double *m_sNormVes)
{
	if (sVslsData == NULL || m_sVes == NULL ||m_sNormVes==NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}
		float MaxVsls=0;
		int m_nImgWX=pVslsInfo->Size[0];
		int m_nImgWY=pVslsInfo->Size[1];
		int m_nImgWZ=pVslsInfo->Size[2];
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				MaxVsls=zxh::maxf(sVslsData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix],MaxVsls);
			}
		}
	}
	if (MaxVsls==0)MaxVsls=1;
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				float fCord[3];
				fCord[0] = (float)ix; fCord[1] = (float)iy; fCord[2] = (float)iz;
				pImgInfo->ImageToWorld(fCord);
				pVslsInfo->WorldToImage(fCord);
				int nCord[3];
				nCord[0] = (int)(fCord[0]+0.5);
				nCord[1] = (int)(fCord[1]+0.5);
				nCord[2] = (int)(fCord[2]+0.5);
				if (nCord[0] >= m_nImgWX)
					nCord[0] = m_nImgWX - 1;
				if (nCord[0] <= 0)
					nCord[0] = 0;
				if (nCord[1] >= m_nImgWY)
					nCord[1] = m_nImgWY - 1;
				if (nCord[1] <= 0)
					nCord[1] = 0;
				if (nCord[2] >= m_nImgWZ)
					nCord[2] = m_nImgWZ - 1;
				if (nCord[2] <= 0)
					nCord[2] = 0;
				m_sVes[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]=(float)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]]);
				m_sNormVes[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]=(float)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]])/MaxVsls;
			}
		}
	}

	return true;
}
bool bInsideImg(int *to,int ImgNewRawSize[4])
{ bool binside=true;
	for(int i=0;i<4;i++)
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
bool bInsideImgNotChange(int to[4],int ImgNewRawSize[4])
{ bool binside=true;
	for(int i=0;i<4;i++)
	{
		if ((int)to[i]<0)
			{
				binside=false;
		}
		if ((int)to[i]>ImgNewRawSize[i]-1)
		{
			binside=false;
		}
	}
	return binside;
}
bool ResoResaTraverseIntensity(const  short *sImgData, const  short *sWinImgData,const zxhImageInfo *pImgInfo,float m_fMean,float m_fStdDev,const short *sVslsData, const short *sWinVslsData,const zxhImageInfo *pVslsInfo,int SearchiRange[4],zxhImageDataT<short> *imgReadNewRaw,zxhImageDataT<short> *imgReadNewVsls)
{       int nImgWX=pVslsInfo->Size[0];
		int nImgWY=pVslsInfo->Size[1];
		int nImgWZ=pVslsInfo->Size[2];
	double *m_sVes=new double[nImgWX * nImgWY * nImgWZ] ;//add by JDQ
	double *m_sNormVes=new double[nImgWX * nImgWY * nImgWZ] ;//add by JDQ
	//initialize the vesselness
	VesselnessInit(pImgInfo, sVslsData, pVslsInfo, m_sVes,m_sNormVes);//Add by JDQ
	int ImgNewRawSize[4]={1};
	int ImgReadRawSize[4]={pImgInfo->Size[0],pImgInfo->Size[1],pImgInfo->Size[2],pImgInfo->Size[3]};

	imgReadNewRaw->GetImageSize(ImgNewRawSize[0],ImgNewRawSize[1],ImgNewRawSize[2],ImgNewRawSize[3]);
	
	for(int it=0;it<ImgNewRawSize[3];++it)
	for(int iz=0;iz<ImgNewRawSize[2];++iz)
	for(int iy=0;iy<ImgNewRawSize[1];++iy)
	for(int ix=0;ix<ImgNewRawSize[0];++ix)
	{
			
		float from[ImageDimensionMax]={ix,iy,iz,it};
		imgReadNewRaw->GetImageInfo()->ImageToWorld(from);
		float toImg[ImageDimensionMax]={from[0],from[1],from[2],from[3]};
		pImgInfo->WorldToImage(toImg);
		int itoImg[ImageDimensionMax]={(int)toImg[0],toImg[1],toImg[2],toImg[3]};
		bInsideImg(itoImg,ImgReadRawSize);
		int locali[4]={0};
		int newix,newiy,newiz,newit;
		double LOCALMIN=ZXH_InfiniteLargeFloat;double dNewImgVal=0;double dNewVslVal=0;
		if(ix==30)
		{int m=ix;}
		for(locali[2]=itoImg[2]-SearchiRange[2];locali[2]<=itoImg[2]+SearchiRange[2];locali[2]++)
        for(locali[1]=itoImg[1]-SearchiRange[1];locali[1]<=itoImg[1]+SearchiRange[1];locali[1]++)
		for(locali[0]=itoImg[0]-SearchiRange[0];locali[0]<=itoImg[0]+SearchiRange[0];locali[0]++)
		{
			if(!bInsideImgNotChange(locali,ImgReadRawSize)) continue;
			
			double dSims=0;
	        double dImgVal=0;
			
			double dvmultiplys=0;
			
				dImgVal = (double)sWinImgData[locali[2] * nImgWY * nImgWX + locali[1] * nImgWX +locali[0]];
				if ( dImgVal < m_fMean)
				{
					dSims= - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSims = exp(dSims);
				}
				else
					dSims= 1;
		dvmultiplys=1/((m_sNormVes[locali[2] * nImgWY * nImgWX + locali[1] * nImgWX + locali[0]])*dSims+0.001);
		if (dvmultiplys<LOCALMIN)
		{
			LOCALMIN=dvmultiplys;
			newix=locali[0];
			newiy=locali[1];
			newiz=locali[2];
			newit=locali[3];
			dNewImgVal=sImgData[newiz * nImgWY * nImgWX + newiy* nImgWX +newix];
			dNewVslVal=sVslsData[newiz * nImgWY * nImgWX + newiy* nImgWX +newix];
		}

		  }
		imgReadNewRaw->SetPixelByGreyscale(ix,iy,iz,it,dNewImgVal);//Normalized intensity
		imgReadNewVsls->SetPixelByGreyscale(ix,iy,iz,it,dNewVslVal);//Normalized intensity
		//cout<<"x:"<<ix<<"y:"<<iy<<"z:"<<iz<<"Intensity:"<<dNewImgVal<<"vesselness:"<<dNewVslVal<<"\n";
	}
	delete[]m_sVes;
	delete[]m_sNormVes;
	
	return true;
}

//int main(int argc, char *argv[])
//{
//
//	/*string strFileNameVsls = string(argv[1]);
//	string strFileNameRaw = string(argv[2]);
//	
//	char *chModelMapCAFileName = argv[3];
//	char *chResultPath = argv[4];
//
//	char *chRCAFileName = argv[5];
//	char *chLADFileName = argv[6];	
//
//	char *chAtlsMeanStdDevFileName = argv[7];
//	char *chUnseenMeanStdDevFileName = argv[8];
//
//	string strFileNameAtls = string(argv[9]);
//	char *chFileNameAtlsCA = argv[10];
//
//	
//
//	int nMaxSgmtNum = atoi(argv[11]);*/
//	
//	//image data,vesselness data
//	string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image02_vesselness.nii";
//	string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";
//	//string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image02_vesselness.nii";
//	//string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";
//
//	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/PotentialPlus/newepsilon_0.01";
//	char *chResultPath = "F:/Coronary_0/MultiResolution_results";
//	//intensity settings data
//	char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_mean_std.txt";
//	char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_to_atlas02_mean_std.txt";
//	string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";
//	char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/reference.vtk";
//	/*char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_mean_std.txt";
//	char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_to_atlas02_mean_std.txt";
//	string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
//	char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset00/vessel0/reference.vtk";*/
//	//
//	double fAtlsCAMean, fAtlsCAStdDev;
//	double gAortaMean[2], gAortaStdDev[2];
//	double fUnseenCAMean, fUnseenCAStdDev; //
//	
//	ReadMeanStdDevFromTxt(chAtlsMeanStdDevFileName, gAortaMean[0], gAortaStdDev[0]);
//	ReadMeanStdDevFromTxt(chUnseenMeanStdDevFileName, gAortaMean[1], gAortaStdDev[1]);
////	gAortaMean[0] --- atlas
////	gAortaMean[1] --- unseen image
////	gAortaStdDev[0] --- atlas
////	gAortaStdDev[1] --- unseen image
//
//	CalcAtlsCAMeanStdDev(strFileNameAtls, chFileNameAtlsCA, fAtlsCAMean, fAtlsCAStdDev);//calculate the mean and stddev of the model image
//	CalcUnseenMeanStdDev(gAortaMean, gAortaStdDev, fAtlsCAMean, fAtlsCAStdDev, fUnseenCAMean, fUnseenCAStdDev);//calculate the mean and stddev of the useen image
//	
//	// for Model with the correction
//	string strMCFilename = string(chResultPath);
//	strMCFilename = strMCFilename + "/" + string("MethodMC.txt");
//	ofstream WriteFileMC(strMCFilename);
//
//	bool bRunResult;
//	bool resave=false,resample=true;
//	//read a "nifti" file (imgdata and vesselness data)
//	zxhImageDataT<short> imgReadVsls, imgReadRaw;//Change by JDQ
//	zxhImageDataT<short> *imgReadNewRaw=new zxhImageDataT<short>;
//    zxhImageDataT<short> *imgReadNewVsls=new zxhImageDataT<short>;
//	if( zxh::OpenImage( &imgReadRaw, strFileNameRaw ) == false )
//	{
//		std::cerr << "Raw image(nifti-file) is not found!"; 
//		return -1;
//	}
//
//	if( zxh::OpenImage( &imgReadVsls, strFileNameVsls ) == false )
//	{
//		std::cerr << "Vesselness image(nifti-file) is not found!"; 
//		return -1;
//	}
//
//	
//	// get the image data
//	const short *sImData = imgReadRaw.GetImageData();
//	// get the image spacing
//	float fNewImgSpacing[]={1,1,1,1};//Add by JDQ
//	float fOldImgSpacing[]={1,1,1,1};
//    float fImgExtend[ImageDimensionMax]={1,1,1,1};
//	float div[ImageDimensionMax] = {1,1,1,1};
//
//	if( glbVerboseOutput>0 ) std::cout<<"success: size "<< imgReadRaw.GetImageSize()[0] <<" " <<imgReadRaw.GetImageSize()[1] <<" "<< imgReadRaw.GetImageSize()[2] <<" "<< imgReadRaw.GetImageSize()[3] <<"\n";
//	if( glbVerboseOutput>0 ) std::cout<<"success: spacing "<< imgReadRaw.GetImageSpacing()[0] <<" "<< imgReadRaw.GetImageSpacing()[1] <<" "<< imgReadRaw.GetImageSpacing()[2] <<" "<< imgReadRaw.GetImageSpacing()[3] <<"\n";
//	imgReadRaw.GetImageSpacing(fOldImgSpacing[0],fOldImgSpacing[1],fOldImgSpacing[2],fOldImgSpacing[3] );//Add by JDQ
//	//imgReadRaw.GetImageExtent( fImgExtend ) ;
//	float fdiv[4]={40,40,10,1};
//	
//	/*	if (resample){
//	          for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
//				div[idim] =fdiv[idim];
//		}*/
//		/*else if(strcmp(argv[iarg],"-spacing")==0&& iarg<argc-3)
//		{
//			for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
//				fImgSpacing[idim] = static_cast<float>(atof(argv[++iarg]));
//		}
//		else if(strcmp(argv[iarg],"-extent")==0&& iarg<argc-3)
//		{
//			for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
//				fImgExtend[idim] = static_cast<float>(atof(argv[++iarg]));*/
//	int size[ImageDimensionMax]={1,1,1,1};
//	if (resample){
//		imgReadRaw.GetImageInfo()->GetSizeUsingSpacing(fNewImgSpacing,size);
//		}
//	int SearchRange[4]={2};
//	for(int i=0;i<4;i++)
//	{
//		SearchRange[i]=int(fNewImgSpacing[i]/fOldImgSpacing[i]+0.5);
//	}
//	/*for( int i=0; i<imgReadRaw.GetDimension(); ++i )
//		if( zxh::abs(fImgSpacing[i]) < ZXH_FloatPrecision )
//			fImgSpacing[i] = imgReadRaw.GetImageSpacing()[i]; 
//		else if( fImgSpacing[i] < 0 )
//			fImgSpacing[i] = -fImgSpacing[i] * imgReadRaw.GetImageSpacing()[i] ; 
//*/
//
//	/*if( glbVerboseOutput>0 ) std::cout<<"success: subsample divide image grid size by div\n";
//
//		for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
//			std::cout<<div[idim]<<" ";
//		std::cout<<"... ... \n";
//		int size[ImageDimensionMax]={1,1,1,1};
//		for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
//		{
//			fImgSpacing[idim]=imgReadRaw.GetImageSpacing()[idim]*div[idim];
//			size[idim]=static_cast<int>((imgReadRaw.GetImageSize()[idim]-1)/div[idim])+1;
//		}*/
//		imgReadNewRaw->NewImage(imgReadRaw.GetDimension(),size,fNewImgSpacing, imgReadRaw.GetImageInfo());
//		imgReadNewVsls->NewImage(imgReadVsls.GetDimension(),size,fNewImgSpacing, imgReadVsls.GetImageInfo());
//	// get the vesselness data 
//	const short *sVslsData = imgReadVsls.GetImageData();	
//	//get the image size
//	int nImWX, nImWY, nImWZ, nImWT;
//	imgReadRaw.GetImageSize(nImWX, nImWY, nImWZ, nImWT);
//	//get the new resample image size
//	int nReImWX, nReImWY, nReImWZ, nReImWT;
//	imgReadNewRaw->GetImageSize(nReImWX, nReImWY, nReImWZ, nReImWT);
//	// normalize
//	short *sNormImg = new short[nImWX * nImWY * nImWZ];
//	short *sNormVsls = new short[nImWX * nImWY * nImWZ];
//	short sRawMin = 0, sRawMax = 0;
//	short sMolRawMin = 0, sMolRawMax = 0;
//	short sVslsMin = 0, sVslsMax = 50;
//
//	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
//	{
//		if (sVslsData[i] > sVslsMax)
//			sNormVsls[i] = sVslsMax;
//		else
//			sNormVsls[i] = sVslsData[i];
//
//		sNormImg[i] = sImData[i];
//	
//		if (sNormImg[i] < sRawMin)
//			sNormImg[i] = sRawMin; //sRawMin = sNormImg[i];
//		if (sNormImg[i] > sRawMax)
//			sRawMax = sNormImg[i];
//		if (sNormVsls[i] < sVslsMin)
//			sVslsMin = sNormVsls[i];
//		if (sNormVsls[i] > sVslsMax)
//			sVslsMax = sNormVsls[i];
//	}	
//
//	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
//	{
//		sNormImg[i] = 1024 * ((double)(sNormImg[i] - sRawMin) / (sRawMax - sRawMin));
//		sNormVsls[i] = 1024 * ((double)(sNormVsls[i] - sVslsMin) / (sVslsMax - sVslsMin));		
//	}
//	
//	fUnseenCAMean = 1024 * (fUnseenCAMean / (sRawMax - sRawMin));
//	fUnseenCAStdDev = 1024 * (fUnseenCAStdDev / (sRawMax - sRawMin));
//	
//	// create an instance for Multi-Resolution
//	
//	
//	//miiMinPathModel *ResolutionResample = new miiMinPathModel(nImWX, nImWY, nImWZ,fImgSpacing,\
//	//	imgReadRaw.GetImageInfo(),chResultPath, true);//creat a class point
//	
//	ResoResaTraverseIntensity(sNormImg, imgReadRaw.GetImageInfo(),fUnseenCAMean, fUnseenCAStdDev, \
//		sNormVsls, imgReadVsls.GetImageInfo(),SearchRange,imgReadNewRaw,imgReadNewVsls);
//
//		//int nLen1 = strlen(chResultPath) + strlen("/Resmol_image02")+strlen(".nii.gz");
//		//char *chFileName1 = (char *)malloc(nLen1);
//		//strcpy(chFileName1, chResultPath);
//		//strcat(chFileName1, "/Resmol_image02");
//		//strcat(chFileName1, ".nii.gz");
//		string chFileName1(chResultPath);
//		chFileName1 +=  "/Resmol_image02";
//		chFileName1 += ".nii.gz";
//		const zxhImageDataT<short> *img=imgReadNewRaw;
//		zxh::SaveImage(img,chFileName1.c_str());
//	
//		string chFileName2(chResultPath);
//		chFileName2 +=  "/Remoeb_image02_vesselness";
//		chFileName2 += ".nii.gz";
//		const zxhImageDataT<short> *vsl=imgReadNewVsls;
//		zxh::SaveImage(vsl,chFileName2.c_str());
//		
//		
//		
//	delete[] sNormImg;
//	delete[] sNormVsls;
//	delete[] imgReadNewRaw;
//	delete[] imgReadNewVsls;
//}
int main(int argc, char *argv[])//without intensity nomalization
{

	/*string strFileNameVsls = string(argv[1]);
	string strFileNameRaw = string(argv[2]);
	
	char *chModelMapCAFileName = argv[3];
	char *chResultPath = argv[4];

	char *chRCAFileName = argv[5];
	char *chLADFileName = argv[6];	

	char *chAtlsMeanStdDevFileName = argv[7];
	char *chUnseenMeanStdDevFileName = argv[8];

	string strFileNameAtls = string(argv[9]);
	char *chFileNameAtlsCA = argv[10];

	

	int nMaxSgmtNum = atoi(argv[11]);*/
	
	//image data,vesselness data
	string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image02_vesselness.nii";
	string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";
	//string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image02_vesselness.nii";
	//string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";

	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/PotentialPlus/newepsilon_0.01";
	char *chResultPath = "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity";
	//intensity settings data
	char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_mean_std.txt";
	char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_to_atlas02_mean_std.txt";
	string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";
	char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/reference.vtk";
	/*char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_mean_std.txt";
	char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_to_atlas02_mean_std.txt";
	string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset00/vessel0/reference.vtk";*/
	//
	double fAtlsCAMean, fAtlsCAStdDev;
	double gAortaMean[2], gAortaStdDev[2];
	double fUnseenCAMean, fUnseenCAStdDev; //
	
	ReadMeanStdDevFromTxt(chAtlsMeanStdDevFileName, gAortaMean[0], gAortaStdDev[0]);
	ReadMeanStdDevFromTxt(chUnseenMeanStdDevFileName, gAortaMean[1], gAortaStdDev[1]);
//	gAortaMean[0] --- atlas
//	gAortaMean[1] --- unseen image
//	gAortaStdDev[0] --- atlas
//	gAortaStdDev[1] --- unseen image
	
	CalcAtlsCAMeanStdDev(strFileNameAtls, chFileNameAtlsCA, fAtlsCAMean, fAtlsCAStdDev);//calculate the mean and stddev of the model image
	//CalcUnseenMeanStdDev(gAortaMean, gAortaStdDev, fAtlsCAMean, fAtlsCAStdDev, fUnseenCAMean, fUnseenCAStdDev);//calculate the mean and stddev of the useen image
	//ModifiedCalcAtlsCAMeanStdDev(strFileNameAtls, chFileNameAtlsCA, fAtlsCAMean);//calculate the mean and stddev of the model image
	ModifiedCalcUnseenMeanStdDev(gAortaMean, gAortaStdDev, fAtlsCAMean, fUnseenCAMean, fUnseenCAStdDev);//calculate the mean and stddev of the useen image
	// for Model with the correction
	string strMCFilename = string(chResultPath);
	strMCFilename = strMCFilename + "/" + string("MethodMC.txt");
	ofstream WriteFileMC(strMCFilename);

	bool bRunResult;
	bool resave=false,resample=true;
	//read a "nifti" file (imgdata and vesselness data)
	zxhImageDataT<short> imgReadVsls, imgReadRaw;//Change by JDQ
	zxhImageDataT<short> *imgReadNewRaw=new zxhImageDataT<short>;
    zxhImageDataT<short> *imgReadNewVsls=new zxhImageDataT<short>;
	if( zxh::OpenImage( &imgReadRaw, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"; 
		return -1;
	}

	if( zxh::OpenImage( &imgReadVsls, strFileNameVsls ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}

	
	// get the image data
	const short *sImData = imgReadRaw.GetImageData();
	// get the image spacing
	float fNewImgSpacing[]={1,1,1,1};//Add by JDQ
	float fOldImgSpacing[]={1,1,1,1};
    float fImgExtend[ImageDimensionMax]={1,1,1,1};
	float div[ImageDimensionMax] = {1,1,1,1};

	if( glbVerboseOutput>0 ) std::cout<<"success: size "<< imgReadRaw.GetImageSize()[0] <<" " <<imgReadRaw.GetImageSize()[1] <<" "<< imgReadRaw.GetImageSize()[2] <<" "<< imgReadRaw.GetImageSize()[3] <<"\n";
	if( glbVerboseOutput>0 ) std::cout<<"success: spacing "<< imgReadRaw.GetImageSpacing()[0] <<" "<< imgReadRaw.GetImageSpacing()[1] <<" "<< imgReadRaw.GetImageSpacing()[2] <<" "<< imgReadRaw.GetImageSpacing()[3] <<"\n";
	imgReadRaw.GetImageSpacing(fOldImgSpacing[0],fOldImgSpacing[1],fOldImgSpacing[2],fOldImgSpacing[3] );//Add by JDQ
	//imgReadRaw.GetImageExtent( fImgExtend ) ;
	float fdiv[4]={40,40,10,1};
	
	/*	if (resample){
	          for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
				div[idim] =fdiv[idim];
		}*/
		/*else if(strcmp(argv[iarg],"-spacing")==0&& iarg<argc-3)
		{
			for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
				fImgSpacing[idim] = static_cast<float>(atof(argv[++iarg]));
		}
		else if(strcmp(argv[iarg],"-extent")==0&& iarg<argc-3)
		{
			for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
				fImgExtend[idim] = static_cast<float>(atof(argv[++iarg]));*/
	int size[ImageDimensionMax]={1,1,1,1};
	if (resample){
		imgReadRaw.GetImageInfo()->GetSizeUsingSpacing(fNewImgSpacing,size);
		}
	int SearchRange[4]={2};
	for(int i=0;i<4;i++)
	{
		SearchRange[i]=int(fNewImgSpacing[i]/fOldImgSpacing[i]+0.5);
	}
	/*for( int i=0; i<imgReadRaw.GetDimension(); ++i )
		if( zxh::abs(fImgSpacing[i]) < ZXH_FloatPrecision )
			fImgSpacing[i] = imgReadRaw.GetImageSpacing()[i]; 
		else if( fImgSpacing[i] < 0 )
			fImgSpacing[i] = -fImgSpacing[i] * imgReadRaw.GetImageSpacing()[i] ; 
*/

	/*if( glbVerboseOutput>0 ) std::cout<<"success: subsample divide image grid size by div\n";

		for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
			std::cout<<div[idim]<<" ";
		std::cout<<"... ... \n";
		int size[ImageDimensionMax]={1,1,1,1};
		for(int idim=0;idim<imgReadRaw.GetDimension();++idim)
		{
			fImgSpacing[idim]=imgReadRaw.GetImageSpacing()[idim]*div[idim];
			size[idim]=static_cast<int>((imgReadRaw.GetImageSize()[idim]-1)/div[idim])+1;
		}*/
		imgReadNewRaw->NewImage(imgReadRaw.GetDimension(),size,fNewImgSpacing, imgReadRaw.GetImageInfo());
		imgReadNewVsls->NewImage(imgReadVsls.GetDimension(),size,fNewImgSpacing, imgReadVsls.GetImageInfo());
	// get the vesselness data 
	const short *sVslsData = imgReadVsls.GetImageData();	
	//get the image size
	int nImWX, nImWY, nImWZ, nImWT;
	imgReadRaw.GetImageSize(nImWX, nImWY, nImWZ, nImWT);
	//get the new resample image size
	int nReImWX, nReImWY, nReImWZ, nReImWT;
	imgReadNewRaw->GetImageSize(nReImWX, nReImWY, nReImWZ, nReImWT);
	
	/*for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	{
		if (sImData[i] > sRawMax)sRawMax=sImData[i];
		if (sImData[i] < sRawMin)sRawMin=sImData[i];
		if (sVslsData[i] > sVslsMax)sVslsMax=sImData[i];
		if (sVslsData[i] < sVslsMin)sVslsMin=sImData[i];

	}*/
	 //windowing
	short *sWinImg = new short[nImWX * nImWY * nImWZ];
	short *sWinVsls = new short[nImWX * nImWY * nImWZ];
	short sRawMin = 0, sRawMax = 0;
	short sMolRawMin = 0, sMolRawMax = 0;
	short sVslsMin = 0, sVslsMax = 200;

	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	{
		if (sVslsData[i] > sVslsMax)
			sWinVsls[i] = sVslsMax;
		else
			sWinVsls[i] = sVslsData[i];
	
		if (sImData[i] < sRawMin)
			sWinImg[i] = sRawMin; //sRawMin = sNormImg[i];
		else
		sWinImg[i] = sImData[i];
	}	
	
	// create an instance for Multi-Resolution
	
	
	//miiMinPathModel *ResolutionResample = new miiMinPathModel(nImWX, nImWY, nImWZ,fImgSpacing,\
	//	imgReadRaw.GetImageInfo(),chResultPath, true);//creat a class point
	
	ResoResaTraverseIntensity(sImData,sWinImg, imgReadRaw.GetImageInfo(),fUnseenCAMean, fUnseenCAStdDev, \
		sVslsData,sWinVsls, imgReadVsls.GetImageInfo(),SearchRange,imgReadNewRaw,imgReadNewVsls);

		//int nLen1 = strlen(chResultPath) + strlen("/Resmol_image02")+strlen(".nii.gz");
		//char *chFileName1 = (char *)malloc(nLen1);
		//strcpy(chFileName1, chResultPath);
		//strcat(chFileName1, "/Resmol_image02");
		//strcat(chFileName1, ".nii.gz");
		string chFileName1(chResultPath);
		chFileName1 +=  "/Resmol_image02";
		chFileName1 += ".nii.gz";
		const zxhImageDataT<short> *img=imgReadNewRaw;
		zxh::SaveImage(img,chFileName1.c_str());
	
		string chFileName2(chResultPath);
		chFileName2 +=  "/Remoeb_image02_vesselness";
		chFileName2 += ".nii.gz";
		const zxhImageDataT<short> *vsl=imgReadNewVsls;
		zxh::SaveImage(vsl,chFileName2.c_str());
		
		
    delete[]sWinImg;
	delete[]sWinVsls;
	
	delete[] imgReadNewRaw;
	delete[] imgReadNewVsls;
}
int GetTime()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}


 