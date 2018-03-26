

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
	short sInteThrehold=12;
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
	short sInteThrehold=12;
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
	short sInteThrehold=12;
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
bool PROVslsimgWithLine(zxhImageDataT<short> &imgReadVsls,vector<miiCNode<double,float>> vPointWorld,zxhImageDataT<short> &imgPROVsls)
{
	int nPROImWX, nPROImWY, nPROImWZ, nPROImWT;
	imgReadVsls.GetImageSize(nPROImWX, nPROImWY, nPROImWZ, nPROImWT);
	const zxhImageInfo *m_pBaseImgInfo;
	m_pBaseImgInfo=imgReadVsls.GetImageInfo();
	short sgvsls=0;
	for (int it=0;it<nPROImWT;it++)
		for (int iy=0;iy<nPROImWY;iy++)
			for (int ix=0;ix<nPROImWX;ix++)
			{
				short sMaxInten=0;
				for (int iz=0;iz<nPROImWZ;iz++)
				{
					sgvsls=imgReadVsls.GetPixelGreyscale(ix,iy,iz,it);

					if(sgvsls>=sMaxInten)
					{
						sMaxInten=sgvsls;
					}
				}
				imgPROVsls.SetPixelByGreyscale(ix,iy,0,it,sMaxInten);
			}
			//map the line points to the plane
			for (int vi=0;vi<vPointWorld.size();vi++)
			{
				float fcord[3]={0,0,0};
				int ncord[3]={0,0,0};
				fcord[0]=vPointWorld[vi].x;
				fcord[1]=vPointWorld[vi].y;
				fcord[2]=vPointWorld[vi].z;
				m_pBaseImgInfo->WorldToImage(fcord);
				ncord[0]= (int)(fcord[0]+0.5);
				ncord[1]= (int)(fcord[1]+0.5);
				ncord[2]= (int)(fcord[2]+0.5);
				if (ncord[0] >= nPROImWX)
					ncord[0]= nPROImWX - 1;
				if (ncord[0] < 0)
					ncord[0]=0;
				if (ncord[1] >= nPROImWY)
					ncord[1] = nPROImWY - 1;
				if (ncord[1] < 0)
					ncord[1]=0;
				if (ncord[2] >= nPROImWZ)
					ncord[2] = nPROImWZ - 1;
				if (ncord[2] < 0)
					ncord[2]=0;
				imgPROVsls.SetPixelByGreyscale(ncord[0],ncord[1],0,0,500);
			}

	return true;
}
bool PROVslsimg(zxhImageDataT<short> &imgReadVsls,zxhImageDataT<short> &imgPROVsls)
{
	int nPROImWX, nPROImWY, nPROImWZ, nPROImWT;
	imgReadVsls.GetImageSize(nPROImWX, nPROImWY, nPROImWZ, nPROImWT);
	const zxhImageInfo *m_pBaseImgInfo;
	m_pBaseImgInfo=imgReadVsls.GetImageInfo();
	short sgvsls=0;
	for (int it=0;it<nPROImWT;it++)
		for (int iy=0;iy<nPROImWY;iy++)
			for (int ix=0;ix<nPROImWX;ix++)
			{
				short sMaxInten=0;
				for (int iz=0;iz<nPROImWZ;iz++)
				{
					sgvsls=imgReadVsls.GetPixelGreyscale(ix,iy,iz,it);

					if(sgvsls>=sMaxInten)
					{
						sMaxInten=sgvsls;
					}
				}
				imgPROVsls.SetPixelByGreyscale(ix,iy,0,it,sMaxInten);
			}
			

	return true;
}
int main(int argc, char *argv[])//without intensity nomalization
{

	string strFileNameVsls="";
	char *MCLinFileName="";
	string strPROResulVsls1Name="";
	string strLorN="";
	int p=0;
	if( argc < 3 )
	{
		cerr << "usage: " << endl;
		cerr << "zxhcaeDMPProjection 3Dimg.nii.gz Line.vtk 2Dpro.nii.gz" << endl;
		return -1;
	} 
	if (argc==5)
	{
		strFileNameVsls =string(argv[1]);;
		MCLinFileName =argv[2];
		strPROResulVsls1Name =string(argv[3]);
		strLorN=string(argv[4]);
	}
	if (argc==4)
	{
		strFileNameVsls =string(argv[1]);;
		strPROResulVsls1Name =string(argv[2]);
		strLorN=string(argv[3]);
	}
	//string strFileNameVsls =string(argv[1]);//"F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image01_vesselness.nii.gz";
	//string strFileNameResoUpRespVsls = string(argv[2]);//"F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image01.nii.gz";	
	//string strResoMAPToHighMask =argv[3];//"F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/unseen01_resultsADE/vessel0/MCLine_SMOW.nii.gz";
	//string strCMBVslsName =string(argv[4]);//"F:/Coronary_0/Coronary_Niessen/CoronaryCMBImg/CMB_image01_v0.nii.gz";
	//string strCMBVsls1Name =string(argv[5]);//"F:/Coronary_0/Coronary_Niessen/CoronaryCMBImg/CMB_image01_v0_vesselness.nii.gz";
	//string strCMBVsls2Name =string(argv[6]);

	/*string strFileNameVsls ="F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image04_vesselness.nii.gz";
	char *MCLinFileName ="F:/Coronary_0/trainningdataZXHCAEDMP/HighResoResults/unseen04_ADE/vessel2/MCLine_EXT_SMO.vtk";
	string strPROResulVsls1Name ="F:/Coronary_0/Coronary_Niessen/VesselnessProjection/Pro_moeb_image04_vesselness.nii.gz";*/


	/*string strFileNameVsls ="F:/Coronary_0/trainningdataZXHCAEDMP/Occlusion/MPREF_OCC_DI4_Gau.nii.gz";
	char *MCLinFileName ="F:/Coronary_0/trainningdataZXHCAEDMP/Occlusion/reference.vtk";
	string strPROResulVsls1Name ="F:/Coronary_0/trainningdataZXHCAEDMP/Occlusion/Pro_MPREF_OCC_DI4_Gau.nii.gz";*/
	//**command
	//image data,vesselness data
	//string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image00_vesselness.nii.gz";
	//string strFileNameResoUpRespVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLungFromBLToH/higvslsFromL00_v0_vesselness.nii.gz";//best low-resolution image is upsampled as higvslsFromL..._vesselness.nii.gz image
	//string chResoMAPToHighMask ="F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults_CMB/unseen00_resultsADE/vessel0/MCLine_SMOW.nii.gz";
	//string strCMBVslsName ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImg/CMB_image00_v0_vesselness.nii.gz";//high-resolution image incorporating MCLine.vtk
	//string strCMBVsls1Name ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImg/CMB_image00_v0_vesselness1.nii.gz";//high-resolution image incorporating low-resolution image v0=vL(v0<=Threshold)
	//string strCMBVsls2Name ="F:/Coronary_0/Coronary_Niessen/CoronaryCMBImg/CMB_image00_v0_vesselness2.nii.gz";//high-resolution image incorporating low-resolution image v0=v0+w*vL(v0<=Threshold)
	//read a "nifti" file (imgdata and vesselness data)
	zxhImageDataT<short> imgReadVsls,imgPROVsls;//Change by JDQ

	if( zxh::OpenImage( &imgReadVsls, strFileNameVsls ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}


	// Creat project new raw image and vesselness
	imgPROVsls.NewImage(imgReadVsls.GetImageInfo());
	//project the vesselness image to a new one
	if(strLorN=="-L")
	{
		vector<PointCordTypeDef>vPathPointsWorld,vSMOPointWorld;
		if (vPathPointsWorld.size()!=0)vPathPointsWorld.clear();
		fstream DetectFile;
		/*char chIdx[25];
		_itoa_s(i, chIdx, 10);*/
		DetectFile.open(MCLinFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(MCLinFileName, vPathPointsWorld);
		}

		DetectFile.close();
		vector<miiCNode<double,float>>vPointWorld;
		miiCNode<double,float> ptemp;
		for (int i=0;i<vPathPointsWorld.size();i++)
		{
			ptemp.x=vPathPointsWorld[i].x;
			ptemp.y=vPathPointsWorld[i].y;
			ptemp.z=vPathPointsWorld[i].z;
			vPointWorld.push_back(ptemp);
		}
		PROVslsimgWithLine(imgReadVsls,vPointWorld,imgPROVsls);
		/*const zxhImageDataT<short> *img=&imgCMBRaw;
		zxh::SaveImage(img,strCMBRawName);*/
		const zxhImageDataT<short> *vsl=&imgPROVsls;
		zxh::SaveImage(vsl,strPROResulVsls1Name);
	}
	if(strLorN=="-N")
	{
		PROVslsimg(imgReadVsls,imgPROVsls);
		/*const zxhImageDataT<short> *img=&imgCMBRaw;
		zxh::SaveImage(img,strCMBRawName);*/
		const zxhImageDataT<short> *vsl=&imgPROVsls;
		zxh::SaveImage(vsl,strPROResulVsls1Name);
	}
}

int GetTime()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}


 