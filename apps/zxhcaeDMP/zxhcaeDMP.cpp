

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

#include <time.h>

#define TAB_CHAR	9

int GetTime() ;


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
	//fMean=1024*fMean/1896;//Add by JDQ
	// read 2nd number
	nStart = nEnd + 1;
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	fStdDev = atof(strNum.c_str());
	//fStdDev=1024*fStdDev/1896;//Add by JDQ

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
bool ReadandSetMolMeanFromTxt(char *chModelMapCAFileName,char *chFileName,vector< miiCNode<double, float> > &vModelPointsWithIntensity)
{
   if(vModelPointsWithIntensity.size()>0)vModelPointsWithIntensity.clear ();
	vector<double>vModelIntensity;
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
inMolIntdata>>data; 
vModelIntensity.push_back(data);
}
ReadVtk(chModelMapCAFileName, vModelPointsWithIntensity);
if(vModelPointsWithIntensity.size()==vModelIntensity.size())
{
		 
	for (int i=0;i<vModelPointsWithIntensity.size();i++)
{
	vModelPointsWithIntensity[i].val=vModelIntensity[i];
}
}
else
	bsizeeque=false;

	}
inMolIntdata.close();
return bsizeeque;
}


void WriteVtk(vector< miiCNode<double, float> > vPointCord, char *chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

/*	int nPointNum = PointCord.size();*/

	float fImgPixel[3];
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
bool ModifiedCalcUnseenMeanStdDev(double gAortaMean[2],double gAortaStdDev[2],double fAtlsCAMean,double &fUnseenCAMean,double &fUnseenCAStdDev)//calculate the first mean std of unseen image.
{
	fUnseenCAMean=fAtlsCAMean-gAortaMean[0]+gAortaMean[1];
	fUnseenCAStdDev=gAortaStdDev[1];
	return true;
}

bool CalcAtlsCAMeanStdDev(string strFileNameAtls, char *strFileNameCA, double &fAtlsCAMean, double &fAtlsCAStdDev)
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
	fCord[0] = vAtlsCAPoints[10].x;
	fCord[1] = vAtlsCAPoints[10].y;
	fCord[2] = vAtlsCAPoints[10].z;
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
bool ModifiedCalcAtlsCAMeanStdDev(vector< miiCNode<double, float> > vModelPointsWithIntensity, double &fAtlsCAMean)
{
	fAtlsCAMean=vModelPointsWithIntensity[0].val;
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

int main(int argc, char *argv[])
{/**/
	
	/*
	if( argc < 12 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << "miiMinPathModelApp	vesselness	image	atlas_vesselx	results	atlas_vessel0	\
			atlas_vessel1	atlas_mean	unseen_mean	atlas_image	atlas_CA	segment_num" << std::endl;
		return -1;
	}
	//set file path and name

	string strFileNameVsls = "E:/UI_Projects/lib/maps/CCTA_0004_80_vesselness.nii.gz"; 
	string strFileNameRaw = "E:/UI_Projects/lib/maps/CCTA_0004_80_vesselness.nii.gz";

	char *chModelMapCAFileName = "E:/UI_Projects/lib/maps/reg/reg_mod_to_ccta_vessel_vtk/reg_mod00_to_ccta0004_80_vessel0.vtk";
	char *chResultPath = "E:/UI_Projects/lib/maps/results/CCTA_0004";

	char *chRCAFileName = "E:/UI_Projects/lib/maps/reg/reg_mod_to_ccta_vessel_vtk/reg_mod00_to_ccta0004_80_vessel0.vtk";
	char *chLADFileName = "E:/UI_Projects/lib/maps/reg/reg_mod_to_ccta_vessel_vtk/reg_mod00_to_ccta0004_80_vessel1.vtk";	

	char *chAtlsMeanStdDevFileName = "E:/UI_Projects/lib/maps/reg/reg_mod_to_ccta_enhanced_mean_std_txt/PriorModel_dataset00_enhanced_mean_std.txt";
	char *chUnseenMeanStdDevFileName = "E:/UI_Projects/lib/maps/reg/reg_mod_to_ccta_enhanced_mean_std_txt/reg_mod00_to_ccta0004_80_enhanced_mean_std.txt";

	string strFileNameAtls = "E:/UI_Projects/lib/coronary/training/DataMaskOutLung/mol_image00.nii.gz";
	char *chFileNameAtlsCA = "E:/UI_Projects/lib/coronary/training/dataset00/vessel0/reference.vtk";

	int nMaxSgmtDist = 100;

	int nMaxSgmtNum = 5;*/
	

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
	
	//string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/image02_vesselness.nii"; 
	/*string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/moeb_image02_vesselness.nii";
	string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";*/
	string strFileNameVsls = "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/Remoeb_image02_vesselness.nii";
	string strFileNameRaw =  "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/Resmol_image02.nii";
	//string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/image02.nii";
	//string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset01/image01_res1mm.nii";
	
	//char *chModelMapCAFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image02_vessel0.vtk";//model points(JDQ)
	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen01_results";
	char *chModelMapCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/SetintensiyByOnePoint/Simulation0.vtk";//model points(JDQ)
	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMaskSiumAsModel/ModifiedMethedwithpreviousmeanstd";
	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMaskSiumAsModel/ModifiedMethedwithintensityupdating";
	
    char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing";
	
	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/withintensityupdated";
	//char *chRCAFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image02_vessel0.vtk";
	//char *chLADFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image02_vessel1.vtk";	
	char *chRCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/vessel0/reference.vtk";
	char *chLADFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/vessel1/reference.vtk";	


	char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_mean_std.txt";
	char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_to_atlas02_mean_std.txt";
	char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/SetintensiyByOnePoint/Simuintensity0.txt";
	//string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image02.nii";
	//string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	//char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/mol_first2other20mm5PontDeformModified/Simulation0.vtk";
	//char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/reference.vtk";

	//int nMaxSgmtDist = 25;
	
	int nMaxSgmtNum =17;
	
	// define variables for Mean and Std-Deviation
	

	double fAtlsCAMean, fAtlsCAStdDev;
	double gAortaMean[2], gAortaStdDev[2];
	double fUnseenCAMean, fUnseenCAStdDev; //
	double dDiffUnseenToAtlasAortaMean;
	double dRatioUnseenToAtlasAortaStdDev;
	
	ReadMeanStdDevFromTxt(chAtlsMeanStdDevFileName, gAortaMean[0], gAortaStdDev[0]);
	ReadMeanStdDevFromTxt(chUnseenMeanStdDevFileName, gAortaMean[1], gAortaStdDev[1]);
//	gAortaMean[0] --- atlas
//	gAortaMean[1] --- unseen image
//	gAortaStdDev[0] --- atlas
//	gAortaStdDev[1] --- unseen image

//read model points intensity
	
	vector< miiCNode<double, float> > vModelPointsWithIntensity;

	if(!ReadandSetMolMeanFromTxt(chModelMapCAFileName,chMolPointsIntensityFileName,vModelPointsWithIntensity))//modified method to get the model-points intensity based on input intensity file.
	{
		std::cerr << "size of model points and model points intensity are not the same !"; 
		return -1;
	}
	//CalcAtlsCAMeanStdDev(strFileNameAtls, chFileNameAtlsCA, fAtlsCAMean, fAtlsCAStdDev);//calculate the mean and stddev of the model image
	//CalcUnseenMeanStdDev(gAortaMean, gAortaStdDev, fAtlsCAMean, fAtlsCAStdDev, fUnseenCAMean, fUnseenCAStdDev);//calculate the mean and stddev of the useen image
	ModifiedCalcAtlsCAMeanStdDev(vModelPointsWithIntensity, fAtlsCAMean);//set the first model mean intensity.
	ModifiedCalcUnseenMeanStdDev(gAortaMean,gAortaStdDev,fAtlsCAMean,fUnseenCAMean,fUnseenCAStdDev);//calculate the first mean std of unseen image.
	// for Model with the correction
	string strMCFilename = string(chResultPath);
	strMCFilename = strMCFilename + "/" + string("MethodMC.txt");
	ofstream WriteFileMC(strMCFilename);

	bool bRunResult;

	//read a "nifti" file
	zxhImageDataT<short> imgReadVsls, imgReadRaw,imgReadMolRaw; //Change by JDQ
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
	float nImgSpacing[4];//Add by JDQ
	imgReadRaw.GetImageSpacing(nImgSpacing[0],nImgSpacing[1],nImgSpacing[2],nImgSpacing[3] );//Add by JDQ
	// get the vesselness data 
	const short *sVslsData = imgReadVsls.GetImageData();	
	
	//get the model image data
	const short *sMolData = imgReadMolRaw.GetImageData();

	//get the image size
	int nImWX, nImWY, nImWZ, nImWT;
	imgReadRaw.GetImageSize(nImWX, nImWY, nImWZ, nImWT);
	
	// normalize
	short *sWinImg = new short[nImWX * nImWY * nImWZ];
	short *sWinVsls = new short[nImWX * nImWY * nImWZ];
	

	short sRawMin = 0, sRawMax = 0;
	short sMolRawMin = 0, sMolRawMax = 0;
	short sVslsMin = 0, sVslsMax = 200;
	short sVslsMax1=0;short sVslsMin1=10000;
	short sImgMax1=0;short sImgMin1=1000;
	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	{
		if (sVslsData[i] > sVslsMax1)
			sVslsMax1=sVslsData[i] ;
		if (sVslsData[i] < sVslsMin1)
			sVslsMin1=sVslsData[i] ;
	
		if (sImData[i]> sImgMax1)
			sImgMax1=sImData[i]; //sRawMin = sNormImg[i];
		if (sImData[i] < sImgMin1)
			sImgMin1=sImData[i] ;
	}	
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

	//short sVslsMax1=0;short sVslsMin1=10000;
	//short sImgMax1=0;short sImgMin1=1000;
	//for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	//{
	//	if (sWinVsls[i] > sVslsMax1)
	//		sVslsMax1=sWinVsls[i] ;
	//	if (sWinVsls[i] < sVslsMin1)
	//		sVslsMin1=sWinVsls[i] ;
	//
	//	if (sWinImg[i]> sImgMax1)
	//		sImgMax1=sWinImg[i]; //sRawMin = sNormImg[i];
	//	if (sWinImg[i] < sImgMin1)
	//		sImgMin1=sWinImg[i] ;
	//}	

	
	
	// change the orientation for test!!!
/*
	short *sTempData = new short[nImWX * nImWY * nImWZ];
//	zxhImageData imgTempData ;
//	imgTempData.NewImage(imgReadRaw.GetDimension(),imgReadRaw.GetImageSize(), imgReadRaw.GetImageSpacing(), imgReadRaw.GetImageInfo() ) ; 
	float fCord[3] = {0};
	int nCord[3] = {0} ;
	for (int iz = 0; iz < nImWZ; iz++)
	{
		for (int iy = 0; iy < nImWY; iy++)
		{
			for (int ix = 0; ix < nImWX; ix++)
			{
				fCord[0] = ix; fCord[1] = iy; fCord[2] = iz;
				imgReadRaw.GetImageInfo()->ImageToWorld(fCord);
				imgReadVsls.GetImageInfo()->WorldToImage(fCord);
				nCord[0] = (int)(fCord[0] + 0.5);
				nCord[1] = (int)(fCord[1] + 0.5);
				nCord[2] = (int)(fCord[2] + 0.5);

//				imgReadRaw.GetImageInfo()->GridToIndex()
//				imgTempData.SetPixelByGreyscale( ix,iy,iz,0, imgReadVsls.GetPixelGreyscaleClosest(nCord[0],nCord[1],nCord[2],0) ) ; 
	
				sTempData[iz * nImWY * nImWX + iy * nImWX + ix] = \
					sNormVsls[nCord[2] * nImWY * nImWX + nCord[1] * nImWX + nCord[0]];
			}
		}
	}
	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	{
//		sNormVsls[i] = sTempData[i];
	}
	delete[] sTempData;
*/
	// read model points from "*.vtk" file
	vector< miiCNode<double, float> > vModelPoints;
	ReadVtk(chModelMapCAFileName, vModelPoints);
	//SetMolPointsIntensity(chModelMapCAFileName,strFileNameAtls,vModelPointsWithIntensity);//set the model-points intensity based on neighbour points intensiy.
	// read the mapped RCA and LAD for calculating the starting point at the model
	vector< miiCNode<double, float> > vRCAPoints;
	ReadVtk(chRCAFileName, vRCAPoints);
	vector< miiCNode<double, float> > vLADPoints;
	ReadVtk(chLADFileName, vLADPoints);

	float nStartPoint[3];
	float nStartPointWithIntisity[5];
	float nEndPoint[3];
	// get starting point
	nStartPoint[0] = (vRCAPoints[0].x + vLADPoints[0].x) / 2;
	nStartPoint[1] = (vRCAPoints[0].y + vLADPoints[0].y) / 2;
	nStartPoint[2] = (vRCAPoints[0].z + vLADPoints[0].z) / 2;
	nStartPointWithIntisity[0]=nStartPoint[0];
	nStartPointWithIntisity[1]=nStartPoint[1];
	nStartPointWithIntisity[2]=nStartPoint[2];
	nStartPointWithIntisity[3]=gAortaMean[0]; 
	nStartPointWithIntisity[4]=gAortaStdDev[0];
// get end point
	nEndPoint[0] = vModelPoints[vModelPoints.size()-1].x;
	nEndPoint[1] = vModelPoints[vModelPoints.size()-1].y;
	nEndPoint[2] = vModelPoints[vModelPoints.size()-1].z;
	
	// create an instance for Model-based minimal path with the correction
	miiMinPathModel *iMinPathModelCrct = new miiMinPathModel(nImWX, nImWY, nImWZ,nImgSpacing,\
		imgReadRaw.GetImageInfo(),chResultPath, true);//creat a class point
	
	
	// get current time(unit: s)
	int nET = GetTime();
	int nMaxItrNum = 2000000;
	int nSaveItrNum=200;
	int nItrNum=0;
	int nItrNumSum=0;

	/************** for Model-based method with the corrction *******************/
	cout << "\nModel-based method with the correction ..." << endl;

	// set the model
	iMinPathModelCrct->SetModelPoints(vModelPoints, nStartPoint);
	//set the model with intensity Add by JDQ
	//iMinPathModelCrct->SetModelPointsWithIntensity(vModelPointsWithIntensity,nStartPointWithIntisity);//transfor the intensity into public model intensity vector
	iMinPathModelCrct->SetModelPointsWithIntensity(vModelPointsWithIntensity,nStartPointWithIntisity);
	//set the mean and stddev of the start point for atlas image
	iMinPathModelCrct->SetMeanStdDevAtlasImg(gAortaMean[0], gAortaStdDev[0]);

	//set the mean and stddev of the start point for unseen image
	iMinPathModelCrct->SetMeanStdDevUnseenImg(gAortaMean[1], gAortaStdDev[1]);

	// set the first mean and std-dev of Unseen Image
	iMinPathModelCrct->SetMeanStdDev(fUnseenCAMean, fUnseenCAStdDev);

	
	// FMM initialization
	//iMinPathModelCrct->NewFastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), \
		sNormVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);
	//iMinPathModelCrct->FastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), \
		sNormVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);//sNormimg here 3 means ( Intensity + Vesselness + Direction)
	iMinPathModelCrct->ModifiedFastMarchingInit(sWinImg, imgReadRaw.GetImageInfo(), \
		sWinVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);
	// save starting point to "txt"
	WriteFileMC << "{" << nStartPoint[0] << ", " << nStartPoint[1] << ", " << nStartPoint[2] << "}(ModelCorrection)" << "->" << endl;
	
	for (int i = 0; i < nMaxSgmtNum; i++)
	{			
		cout << "Segment: " << i+1 << endl;
		// get current time(unit: s)
		nET = GetTime();
		
		// fast-marching-method evolution
		//nItrNum = iMinPathModelCrct->FastMarchingEvolution(nMaxItrNum,nItrNumSum);
		//nItrNum = iMinPathModelCrct->NewFastMarchingEvolution(nMaxItrNum,nItrNumSum);
		//nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolution(nMaxItrNum,nItrNumSum,sWinImg);
		nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModel(nMaxItrNum,nItrNumSum,sWinImg);

	    nItrNumSum=nItrNumSum+nItrNum;
		// calculate the elapsed time for FMM evolution
		int nUseTime = GetTime() - nET;

		if (nItrNum < nMaxItrNum && nItrNum > 0)
		{
			// find the minimal path by using back-propagation
			cout << "Finding minimal path..." << endl;
			if (iMinPathModelCrct->FindMinPath(nStartPoint, 2000))
			{
				bRunResult = true;
			}
			else
			{
				bRunResult = false;
				// output the elapsed time (uint: s)
				cout << "Elapsed time: " << nUseTime << "s" << endl;
				cout << "Finding min-path failed!" << endl;
				break;
			}			
		}
		else
		{
			bRunResult = false;
			// output the elapsed time (uint: s)
			cout << "Elapsed time: " << nUseTime << "s" << endl;
			cout << "FMM failed!" << endl;
			break;
		}

		// output the elapsed time (uint: s)
		cout << "Elapsed time: " << nUseTime << "s" << endl;

		// save the text results
		cout << "Saving results..." << endl;
		float fMCEndPointWorld[3];
		//iMinPathModelCrct->GetEndPoint(nMCEndPoint);
		iMinPathModelCrct->GetEndPointWorld(fMCEndPointWorld);
		WriteFileMC << i+1 << "/" << nMaxSgmtNum << ": ";
		WriteFileMC <<right<<fixed<<setfill('0')<<setprecision(1) \
			<< "Pos{" << fMCEndPointWorld[0] << ", " << fMCEndPointWorld[1] << ", " << fMCEndPointWorld[2] << "}, ";
		WriteFileMC << "Elapsed time" << "(" << nUseTime << "s), ";
		if (bRunResult)
			WriteFileMC << "Successful!" << endl;
		else
			WriteFileMC << "Failed!" << endl;

		// save the image results
		
/*		stringstream ss; // c++中用于格式化字符串的流
ss << value;
string strTemp;
ss >> strTemp; // 此时strTemp变为"100"
		string strFilenameMC = string(chResultPath);
		strFilenameMC = strFilenameMC + "/MCImageLine" + string(chTemp) + ".nii.gz";
		iMinPathModelCrct->GenMinPathGraph(imgReadVsls.GetImageInfo(), sVslsData, strFilenameMC, 0);

		strFilenameMC = string(chResultPath);
		strFilenameMC = strFilenameMC + "/MCLine" + string(chTemp) + ".nii.gz";
		iMinPathModelCrct->GenMinPathGraph(imgReadVsls.GetImageInfo(), sVslsData, strFilenameMC, 1);
*/
		// save U in minimal path
		/*
		string strFilenameU2 = string(chResultPath);
		strFilenameU2 = strFilenameU2 + "/ModelCrctU" + string(chTemp) + ".nii.gz";
		iMinPathModelCrct->SaveU2Image(imgReadRaw.GetImageInfo(), strFilenameU2);
		*/
		// save CA to vtk file
		char chTemp[25];
		_itoa_s(i, chTemp, 10);
		int nLen = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".vtk") + 1;
		char *chFileName2 = (char *)malloc(nLen);
		strcpy(chFileName2, chResultPath);
		strcat(chFileName2, "/MCLine");
		strcat(chFileName2, chTemp);
		strcat(chFileName2, ".vtk");
		iMinPathModelCrct->WriteCA2Vtk(chFileName2);
		free(chFileName2);

		
		int nLen1 = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".txt") + 1;
		char *chFileName3 = (char *)malloc(nLen1);
		strcpy(chFileName3, chResultPath);
		strcat(chFileName3, "/MCLine");
		strcat(chFileName3, chTemp);
		strcat(chFileName3, ".txt");
		iMinPathModelCrct->WriteCA2Txt(chFileName3);
		free(chFileName3);
		
		
/*
		nLen = strlen(chResultPath) + strlen("/CrctModel") + strlen(chTemp) + strlen(".txt") + 1;
		chFileName3 = (char *)malloc(nLen);
		strcpy(chFileName3, chResultPath);
		strcat(chFileName3, "/CrctModel");
		strcat(chFileName3, chTemp);
		strcat(chFileName3, ".txt");
		iMinPathModelCrct->GetModelPoints(vModelPoints);
		WriteVtk(vModelPoints, chFileName3);		
		free(chFileName3);
*/
		cout << "Results saved!\n" << endl;
		/*if (!iMinPathModelCrct->ReFastMarchingInit())
		{
			break;
		}*/
		/*if (!iMinPathModelCrct->ModifiedReFastMarchingInit())
		{
			break;
		}*/
		if (!iMinPathModelCrct->ModifiedReFastMarchingInitdontMoveModel())
		{
			break;
		}
		/*if (!iMinPathModelCrct->NewReFastMarchingInit()||!iMinPathModelCrct->UpdateUnseenImgIntensity(sNormImg))
		{
			break;
		}*/
		// update the whole intensity parameterAdd by JDQ
	}

	delete iMinPathModelCrct;
	WriteFileMC.close();
	
	// release memory
	delete[] sWinImg;
	delete[] sWinVsls;
	
	imgReadRaw.ReleaseMem();
	imgReadVsls.ReleaseMem();

	cout << "Run over!\n" << endl;
	return 1;	
}


int GetTime()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}


 