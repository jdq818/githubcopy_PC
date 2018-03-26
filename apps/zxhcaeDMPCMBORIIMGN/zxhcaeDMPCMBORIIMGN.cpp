

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

bool ReadModelPoiIntenFromTxt(char *chFileName,vector< double > &vModelPointsIntensity)
{
	if(vModelPointsIntensity.size()>0)vModelPointsIntensity.clear ();
	ifstream inMolIntdata(chFileName);
	double data;string strLine="";
	bool bsizeeque=true;
	if(inMolIntdata)
	{

		while(getline(inMolIntdata,strLine,'\n'))
		{ 

			data = atof(strLine.c_str());
			vModelPointsIntensity.push_back(data);

		}
	}
	else
	{
		cout << "Unable to open txt-file: " << endl;
		cout << chFileName << endl;
		return false; //exit(1); // terminate with error
	}
	return true;
}
bool CalcUnseenMeanStd(double gAortaMean[2],double gAortaStdDev[2],double fAtlsCAMean,double &fUnseenCAMean,double &fUnseenCAStdDev)
{
	fUnseenCAMean=gAortaMean[1];//fUnseenCAMean=fAtlsCAMean-gAortaMean[0]+gAortaMean[1];
	fUnseenCAStdDev=gAortaStdDev[1];
	return true;
}
bool CMBMaskWithHighRaw(zxhImageDataT<short> &imgReadVsls,zxhImageDataT<short> &imgReadMask,zxhImageDataT<short> &imgCMBVsls)
{
		//get the new resample image size
	int nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT;
	imgCMBVsls.GetImageSize(nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT);
	short sInteThrehold=2;
	for (int it=0;it<nCMBImWT;it++)
		for (int iz=0;iz<nCMBImWZ;iz++)
			for (int iy=0;iy<nCMBImWY;iy++)
				for (int ix=0;ix<nCMBImWX;ix++)
				{
					
					short sgvsls=imgReadVsls.GetPixelGreyscale(ix,iy,iz,it);
					short sgmask=imgReadMask.GetPixelGreyscale(ix,iy,iz,it);
						if(sgvsls<=sInteThrehold&&sgmask==ZXH_Foreground)
					{
						imgCMBVsls.SetPixelByGreyscale(ix,iy,iz,it,sInteThrehold);
					}
					else
					{
						imgCMBVsls.SetPixelByGreyscale(ix,iy,iz,it,sgvsls);
					}
				}
	return true;
}
bool CMBLVslsWithHighRaw1(zxhImageDataT<short> &imgReadVsls,zxhImageDataT<short> &imgResoUpRespVsls,zxhImageDataT<short> &imgCMBVsls)
{
		//get the new resample image size
	int nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT;
	imgCMBVsls.GetImageSize(nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT);
	short sInteThrehold=2;
	for (int it=0;it<nCMBImWT;it++)
		for (int iz=0;iz<nCMBImWZ;iz++)
			for (int iy=0;iy<nCMBImWY;iy++)
				for (int ix=0;ix<nCMBImWX;ix++)
				{
					
					short sgvsls=imgReadVsls.GetPixelGreyscale(ix,iy,iz,it);
					short sgLvsls=imgResoUpRespVsls.GetPixelGreyscale(ix,iy,iz,it);
						if(sgvsls<=sInteThrehold)
					{
						imgCMBVsls.SetPixelByGreyscale(ix,iy,iz,it,sgLvsls);
					}
					else
					{
						imgCMBVsls.SetPixelByGreyscale(ix,iy,iz,it,sgvsls);
					}
				}
	return true;
}
bool CMBLVslsWithHighRaw2(zxhImageDataT<short> &imgReadVsls,zxhImageDataT<short> &imgResoUpRespVsls,zxhImageDataT<short> &imgCMBVsls)
{
		//get the new resample image size
	int nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT;
	imgCMBVsls.GetImageSize(nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT);
	short sInteThrehold=2;
	float w=0.1;
	for (int it=0;it<nCMBImWT;it++)
		for (int iz=0;iz<nCMBImWZ;iz++)
			for (int iy=0;iy<nCMBImWY;iy++)
				for (int ix=0;ix<nCMBImWX;ix++)
				{
					
					short sgvsls=imgReadVsls.GetPixelGreyscale(ix,iy,iz,it);
					short sgLvsls=imgResoUpRespVsls.GetPixelGreyscale(ix,iy,iz,it);
						if(sgvsls<=sInteThrehold)
					{
						short sgNvsls=zxh::minf(sInteThrehold, sgvsls+w*sgLvsls);
						imgCMBVsls.SetPixelByGreyscale(ix,iy,iz,it,sgNvsls);
					}
					else
					{
						imgCMBVsls.SetPixelByGreyscale(ix,iy,iz,it,sgvsls);
					}
				}
	return true;
}
bool CMBLVslsWithHighRaw2N(zxhImageDataT<short> &imgReadVsls,zxhImageDataT<short> &imgReadVslsxyz,zxhImageDataT<short> &imgResoUpRespVslsxyz,zxhImageDataT<short> &imgCMBVslsxyz)
{
		//get the new resample image size
	int nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT;
	imgCMBVslsxyz.GetImageSize(nCMBImWX, nCMBImWY, nCMBImWZ, nCMBImWT);
	float fInteThrehold=2;
	float fInteThreholdxyz=1;
	float w=0.1;
	for (int it=0;it<nCMBImWT;it++)
		for (int iz=0;iz<nCMBImWZ;iz++)
			for (int iy=0;iy<nCMBImWY;iy++)
				for (int ix=0;ix<nCMBImWX;ix++)
				{
					
					float fgvsls=imgReadVsls.GetPixelGreyscale(ix,iy,iz,it);
					float fgvslsxyz=imgReadVslsxyz.GetPixelGreyscale(ix,iy,iz,it);
					short sgvslsxyz=(short)fgvslsxyz;
					float fgLvslsxyz=imgResoUpRespVslsxyz.GetPixelGreyscale(ix,iy,iz,it);
					if(fgvsls<=fInteThrehold)
					{
						float fgvslsNxyz=fgvslsxyz+w*fgLvslsxyz;
						float fabsgNvslsxyz=zxh::minf(fInteThreholdxyz, abs(fgvslsxyz+w*fgLvslsxyz));

						if(fgvslsNxyz<0)
						{
							short sgNvslsxyz=(short)-1*fabsgNvslsxyz;
							imgCMBVslsxyz.SetPixelByGreyscale(ix,iy,iz,it,sgNvslsxyz);
						}
						else
						{
							short sgNvslsxyz=(short)fabsgNvslsxyz;
							imgCMBVslsxyz.SetPixelByGreyscale(ix,iy,iz,it,sgNvslsxyz);
						}
					}
					else
					{
						imgCMBVslsxyz.SetPixelByGreyscale(ix,iy,iz,it,sgvslsxyz);
					}
				}
	return true;
}
int main(int argc, char *argv[])//without intensity nomalization
{

	
	int p=0;
	if( argc < 15 )
	{
		cerr << "usage: " << endl;
		cerr << "zxhcaeDMPCMBORIIMG	vessleness.nii ResoUpRespVsls.nii ResoMAPToHighMask vslsresult vsls1result vsls2result " << endl;
		return -1;
	} 
	
	string strFileNameVsls =string(argv[1]);
	string strFileNameResoUpRespVsls =string(argv[2]);//best low-resolution image is upsampled as higvslsFromL..._vesselness.nii.gz image
	string chResoMAPToHighMask =string(argv[3]);
	string strFileNameVslsx=string(argv[4]);
	string strFileNameResoUpRespVslsx =string(argv[5]);
	string strFileNameVslsy=string(argv[6]);
	string strFileNameResoUpRespVslsy =string(argv[7]);
	string strFileNameVslsz=string(argv[8]);
	string strFileNameResoUpRespVslsz =string(argv[9]);
	string strCMBVslsName =string(argv[10]);//high-resolution image incorporating MCLine.vtk
	string strCMBVsls1Name =string(argv[11]);//high-resolution image incorporating low-resolution image v0=vL(v0<=Threshold)
	string strCMBVsls2Name =string(argv[12]);//high-resolution image incorporating low-resolution image v0=v0+w*vL(v0<=Threshold)
	string strCMBVslsx2Name =string(argv[13]);
	string strCMBVslsy2Name =string(argv[14]);
	string strCMBVslsz2Name =string(argv[15]);
	//**command
	//image data,vesselness data
	//string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image00_vesselness.nii.gz";
	//string strFileNameResoUpRespVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLungFromBLToHN/higvslsFromL00_v0_vesselness.nii.gz";//best low-resolution image is upsampled as higvslsFromL..._vesselness.nii.gz image
	//string chResoMAPToHighMask ="F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResultsN/unseen00_results/vessel0/MCLine_SMO.nii.gz";
	//string strFileNameVslsx="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung_U/vesselness_u1/mol_image00_vesselness_new_x.nii.gz";
	//string strFileNameResoUpRespVslsx = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLungFromBLToHN/higvslsFromL00_v0_vesselness_X.nii.gz";
	//string strFileNameVslsy="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung_U/vesselness_u1/mol_image00_vesselness_new_y.nii.gz";
	//string strFileNameResoUpRespVslsy = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLungFromBLToHN/higvslsFromL00_v0_vesselness_Y.nii.gz";
	//string strFileNameVslsz="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung_U/vesselness_u1/mol_image00_vesselness_new_z.nii.gz";
	//string strFileNameResoUpRespVslsz = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLungFromBLToHN/higvslsFromL00_v0_vesselness_Z.nii.gz";
	//string strCMBVslsName ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImgN/CMB_image00_v0_vesselness.nii.gz";//high-resolution image incorporating MCLine.vtk
	//string strCMBVsls1Name ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImgN/CMB_image00_v0_vesselness1.nii.gz";//high-resolution image incorporating low-resolution image v0=vL(v0<=Threshold)
	//string strCMBVsls2Name ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImgN/CMB_image00_v0_vesselness2.nii.gz";//high-resolution image incorporating low-resolution image v0=v0+w*vL(v0<=Threshold)
	//string strCMBVslsx2Name ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImgN/CMB_image00_v0_vesselness2_x.nii.gz";
	//string strCMBVslsy2Name ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImgN/CMB_image00_v0_vesselness2_y.nii.gz";
	//string strCMBVslsz2Name ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImgN/CMB_image00_v0_vesselness2_z.nii.gz";
	//read a "nifti" file (imgdata and vesselness data)
	zxhImageDataT<short> imgReadVsls,imgResoUpRespVsls,imgReadVslsx,imgResoUpRespVslsx,imgReadVslsy,imgResoUpRespVslsy,imgReadVslsz,imgResoUpRespVslsz,imgReadMask,imgCMBVsls,imgCMBVsls1,imgCMBVsls2,imgCMBVsls2x,imgCMBVsls2y,imgCMBVsls2z;//Change by JDQ

	if( zxh::OpenImage( &imgReadVsls, strFileNameVsls ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgResoUpRespVsls, strFileNameResoUpRespVsls ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgReadVslsx, strFileNameVslsx ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgResoUpRespVslsx, strFileNameResoUpRespVslsx ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgReadVslsy, strFileNameVslsy ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgResoUpRespVslsy, strFileNameResoUpRespVslsy ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgReadVslsz, strFileNameVslsz ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgResoUpRespVslsz, strFileNameResoUpRespVslsz ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgReadMask, chResoMAPToHighMask ) == false )
	{
		std::cerr << "ResoMAPToHighMask image(nifti-file) is not found!"; 
		return -1;
	}


	// Creat combined new raw image and vesselness
	imgCMBVsls.NewImage(imgReadVsls.GetImageInfo());
	imgCMBVsls1.NewImage(imgReadVsls.GetImageInfo());
	imgCMBVsls2.NewImage(imgReadVsls.GetImageInfo());
	imgCMBVsls2x.NewImage(imgReadVsls.GetImageInfo());
	imgCMBVsls2y.NewImage(imgReadVsls.GetImageInfo());
	imgCMBVsls2z.NewImage(imgReadVsls.GetImageInfo());
	CMBLVslsWithHighRaw1(imgReadVsls,imgResoUpRespVsls,imgCMBVsls1);
	CMBLVslsWithHighRaw2(imgReadVsls,imgResoUpRespVsls,imgCMBVsls2);
	CMBMaskWithHighRaw(imgReadVsls,imgReadMask,imgCMBVsls);
	CMBLVslsWithHighRaw2N(imgReadVsls,imgReadVslsx,imgResoUpRespVslsx,imgCMBVsls2x);
	CMBLVslsWithHighRaw2N(imgReadVsls,imgReadVslsy,imgResoUpRespVslsy,imgCMBVsls2y);
	CMBLVslsWithHighRaw2N(imgReadVsls,imgReadVslsz,imgResoUpRespVslsz,imgCMBVsls2z);
	/*const zxhImageDataT<short> *img=&imgCMBRaw;
	zxh::SaveImage(img,strCMBRawName);*/
	const zxhImageDataT<short> *vsl=&imgCMBVsls;
	zxh::SaveImage(vsl,strCMBVslsName);
	const zxhImageDataT<short> *vsl1=&imgCMBVsls1;
	zxh::SaveImage(vsl1,strCMBVsls1Name);
	const zxhImageDataT<short> *vsl2=&imgCMBVsls2;
	zxh::SaveImage(vsl2,strCMBVsls2Name);
	const zxhImageDataT<short> *vsl2x=&imgCMBVsls2x;
	zxh::SaveImage(vsl2x,strCMBVslsx2Name);
	const zxhImageDataT<short> *vsl2y=&imgCMBVsls2y;
	zxh::SaveImage(vsl2y,strCMBVslsy2Name);
	const zxhImageDataT<short> *vsl2z=&imgCMBVsls2z;
	zxh::SaveImage(vsl2z,strCMBVslsz2Name);
}

int GetTime()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}


 