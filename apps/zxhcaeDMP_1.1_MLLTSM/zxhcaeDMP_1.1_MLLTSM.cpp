

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
	// read 1st number
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	strctTempPoint.x =atof(strNum.c_str());
	//fMean=1024*fMean/1896;//Add by JDQ
	// read 2nd number
	nStart = nEnd + 1;
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	strctTempPoint.y =atof(strNum.c_str());
	// read 3rd number
	nStart = nEnd + 1;
	nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	strNum.assign(strLine, nStart, nEnd - nStart);
	strctTempPoint.z =atof(strNum.c_str());
	cl.push_back(strctTempPoint);
    return true;
}
bool ReadandSetMolMeanFromTxt(char *chModelMapCAFileName,char *chFileName,vector< miiCNode<double, float> > &vModelPointsWithIntensity)
{
	if(vModelPointsWithIntensity.size()>0)vModelPointsWithIntensity.clear ();
	vector<double>vModelIntensity;
	ifstream inMolIntdata(chFileName);
	double data;string strLine="";
	bool bsizeeque=true;
	if(inMolIntdata)
	{

		while(getline(inMolIntdata,strLine,'\n'))
		{ 

			data = atof(strLine.c_str());
			vModelIntensity.push_back(data);

		}
	}
	else
	{
		cout << "Unable to open txt-file: " << endl;
		cout << chFileName << endl;
		return false; //exit(1); // terminate with error
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
	{
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

	void WriteModCA2Txt(char *chFileName,const vector<double>m_sMolPontIntsFinl)
	{
	ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = m_sMolPontIntsFinl.size();
	float fmodints;	
	for (int i = 0; i < nPointNum; i++)
	{
	 fmodints= m_sMolPontIntsFinl[i];
     WriteFileTxt <<right<<fixed<<setfill('0')<<setprecision(4)<< fmodints<< '\n';

	}	
	
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
	fUnseenCAMean=gAortaMean[1];//fAtlsCAMean-gAortaMean[0]+gAortaMean[1];
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
	//fAtlsCAMean=vModelPointsWithIntensity[0].val;
	float inten=0;
	for(int i=0;i<vModelPointsWithIntensity.size();i++)
	{
		inten=inten+vModelPointsWithIntensity[i].val;

	}

	fAtlsCAMean=inten/(vModelPointsWithIntensity.size());
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
void RandomDefPoint(float RefStart[3])
{
	int RAD_2=2;
		float fFi=0,fTheta=0,R=0;
		srand((unsigned)time(NULL));
		fFi=(float)M_PI* (rand()%360);
		float fRandVector[3]={0,0,0};
			fTheta=(float)M_PI*(rand()%180);
	R=(float)RAD_2*(rand() % 20)/19.0;
	fRandVector[0]=R*sin(fTheta)*cos(fFi);
	fRandVector[1]=R*sin(fTheta)*sin(fFi);
	fRandVector[2]=R*cos(fTheta);
	RefStart[0]=RefStart[0]+fRandVector[0];
	RefStart[1]=RefStart[1]+fRandVector[1];
	RefStart[2]=RefStart[2]+fRandVector[2];
}
bool ReadRefPointTxt(char *filename,vector< miiCNode<double, float> > &cl)
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

int main(int argc, char *argv[])//Use Resolution Mask to extraction
{
	
	
	//if( argc < 12 )
	//{
	//	std::cerr << "Usage: " << std::endl;
	//	std::cerr << "miiMinPathModelApp	vesselness	image	atlas_vesselx	results	atlas_vessel0	\
	//		atlas_vessel1	atlas_mean	unseen_mean	atlas_image	atlas_CA	segment_num" << std::endl;
	//	return -1;
	//}
	//set file path and name
	/*
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
	//
 //   string strFileNameRaw = string(argv[1]);
	//string strFileNameVsls = string(argv[2]);
	//
	//
	//char *chModelMapCAFileName = argv[3];
	//char *chResultPath = argv[4];

	//char *chRCAFileName = argv[5];
	//char *chLADFileName = argv[6];	

	//char *chAtlsMeanStdDevFileName = argv[7];
	//char *chUnseenMeanStdDevFileName = argv[8];
	//char *chMolPointsIntensityFileName=argv[9];
	//char *chFileNameAtlsCA = argv[10];

	//

	//int nMaxSgmtNum = atoi(argv[11]);
	//



	///////////////////////
//	//****SimuModelExtraction***//
//	string strFileNameRaw =  "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/Resmol_image02.nii";
//	string strFileNameVsls = "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/Remoeb_image02_vesselness.nii";
//	char *chModelMapCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/SetintensiyByOnePoint/Simulation19.vtk";//simu-model points(JDQ)
//	char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu19";
//	char *chRCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/vessel0/reference.vtk";
//	char *chLADFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/vessel1/reference.vtk";	
//	char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_mean_std.txt";
//	char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_to_atlas02_mean_std.txt";
//	char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/SetintensiyByOnePoint/Simuintensity19.txt";
//	int nMaxSgmtNum =30;
//	//string strFileNameRaw =  "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/00to01ResoImage/Resmol_image01.nii";
//	//string strFileNameVsls = "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/00to01ResoImage/Remoeb_image01_vesselness.nii";
//	//char *chModelMapCAFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image01_vessel0.vtk";//simu-model points(JDQ)
//	//char *chResultPath = "F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod00_to_unseen01_results/vessel0";
//	//char *chRCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset01/vessel0/reference.vtk";
//	//char *chLADFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset01/vessel1/reference.vtk";	
//	//char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_mean_std.txt";
//	//char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_to_atlas01_mean_std.txt";
//	//char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/SetintensiyByOnePoint/Simuintensity19.txt";
//	//int nMaxSgmtNum =30;
//
//
//	// define variables for Mean and Std-Deviation
//	double fAtlsCAMean, fAtlsCAStdDev;
//	double gAortaMean[2], gAortaStdDev[2];
//	double fUnseenCAMean, fUnseenCAStdDev; //
//	double dDiffUnseenToAtlasAortaMean;
//	double dRatioUnseenToAtlasAortaStdDev;
//	
//	ReadMeanStdDevFromTxt(chAtlsMeanStdDevFileName, gAortaMean[0], gAortaStdDev[0]);
//	ReadMeanStdDevFromTxt(chUnseenMeanStdDevFileName, gAortaMean[1], gAortaStdDev[1]);
////	gAortaMean[0] --- atlas
////	gAortaMean[1] --- unseen image
////	gAortaStdDev[0] --- atlas
////	gAortaStdDev[1] --- unseen image
//
////read model points intensity
//	
//	vector< miiCNode<double, float> > vModelPointsWithIntensity;
//
//	if(!ReadandSetMolMeanFromTxt(chModelMapCAFileName,chMolPointsIntensityFileName,vModelPointsWithIntensity))//modified method to get the model-points intensity based on input intensity file.
//	{
//		std::cerr << "size of model points and model points intensity are not the same !"; 
//		return -1;
//	}
//	//CalcAtlsCAMeanStdDev(strFileNameAtls, chFileNameAtlsCA, fAtlsCAMean, fAtlsCAStdDev);//calculate the mean and stddev of the model image
//	//CalcUnseenMeanStdDev(gAortaMean, gAortaStdDev, fAtlsCAMean, fAtlsCAStdDev, fUnseenCAMean, fUnseenCAStdDev);//calculate the mean and stddev of the useen image
//	ModifiedCalcAtlsCAMeanStdDev(vModelPointsWithIntensity, fAtlsCAMean);//set the first model mean intensity.
//	ModifiedCalcUnseenMeanStdDev(gAortaMean,gAortaStdDev,fAtlsCAMean,fUnseenCAMean,fUnseenCAStdDev);//calculate the first mean std of unseen image.
//	// for Model with the correction
//	string strMCFilename = string(chResultPath);
//	strMCFilename = strMCFilename + "/" + string("MethodMC.txt");
//	ofstream WriteFileMC(strMCFilename);
//
//	bool bRunResult;
//
//	//read a "nifti" file
//	zxhImageDataT<short> imgReadVsls, imgReadRaw,imgReadMolRaw,imgReadLineMask; //Change by JDQ
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
//	/*if( zxh::OpenImage( &imgReadLineMask, strFileNameLineMask ) == false )
//	{
//		std::cerr << "LineMask  image(nifti-file) is not found!"; 
//		return -1;
//	}*/
//
//	// get the image data
//	const short *sImData = imgReadRaw.GetImageData();
//	// get the image spacing
//	float nImgSpacing[4];//Add by JDQ
//	imgReadRaw.GetImageSpacing(nImgSpacing[0],nImgSpacing[1],nImgSpacing[2],nImgSpacing[3] );//Add by JDQ
//	// get the vesselness data 
//	const short *sVslsData = imgReadVsls.GetImageData();	
//	
//	//get the model image data
//	const short *sMolData = imgReadMolRaw.GetImageData();
//
//	//get the line mask data
//	//const short *sLineMaskData = imgReadLineMask.GetImageData();
//	//get the image size
//	int nImWX, nImWY, nImWZ, nImWT;
//	imgReadRaw.GetImageSize(nImWX, nImWY, nImWZ, nImWT);
//	
//	// normalize
//	short *sWinImg = new short[nImWX * nImWY * nImWZ];
//	short *sWinVsls = new short[nImWX * nImWY * nImWZ];
//
//	short sRawMin = 0, sRawMax = 0;
//	short sMolRawMin = 0, sMolRawMax = 0;
//	short sVslsMin = 0, sVslsMax = 200;
//	short sVslsMax1=0;short sVslsMin1=10000;
//	short sImgMax1=0;short sImgMin1=1000;
//	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
//	{
//		if (sVslsData[i] > sVslsMax1)
//			sVslsMax1=sVslsData[i] ;
//		if (sVslsData[i] < sVslsMin1)
//			sVslsMin1=sVslsData[i] ;
//	
//		if (sImData[i]> sImgMax1)
//			sImgMax1=sImData[i]; //sRawMin = sNormImg[i];
//		if (sImData[i] < sImgMin1)
//			sImgMin1=sImData[i] ;
//	}	
//	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
//	{
//		if (sVslsData[i] > sVslsMax)
//			sWinVsls[i] = sVslsMax;
//		else
//			sWinVsls[i] = sVslsData[i];
//	
//		if (sImData[i] < sRawMin)
//			sWinImg[i] = sRawMin; //sRawMin = sNormImg[i];
//		else
//		sWinImg[i] = sImData[i];
//	}	
//
//	//short sVslsMax1=0;short sVslsMin1=10000;
//	//short sImgMax1=0;short sImgMin1=1000;
//	//for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
//	//{
//	//	if (sWinVsls[i] > sVslsMax1)
//	//		sVslsMax1=sWinVsls[i] ;
//	//	if (sWinVsls[i] < sVslsMin1)
//	//		sVslsMin1=sWinVsls[i] ;
//	//
//	//	if (sWinImg[i]> sImgMax1)
//	//		sImgMax1=sWinImg[i]; //sRawMin = sNormImg[i];
//	//	if (sWinImg[i] < sImgMin1)
//	//		sImgMin1=sWinImg[i] ;
//	//}	
//
//	
//	
//	// change the orientation for test!!!
///*
//	short *sTempData = new short[nImWX * nImWY * nImWZ];
////	zxhImageData imgTempData ;
////	imgTempData.NewImage(imgReadRaw.GetDimension(),imgReadRaw.GetImageSize(), imgReadRaw.GetImageSpacing(), imgReadRaw.GetImageInfo() ) ; 
//	float fCord[3] = {0};
//	int nCord[3] = {0} ;
//	for (int iz = 0; iz < nImWZ; iz++)
//	{
//		for (int iy = 0; iy < nImWY; iy++)
//		{
//			for (int ix = 0; ix < nImWX; ix++)
//			{
//				fCord[0] = ix; fCord[1] = iy; fCord[2] = iz;
//				imgReadRaw.GetImageInfo()->ImageToWorld(fCord);
//				imgReadVsls.GetImageInfo()->WorldToImage(fCord);
//				nCord[0] = (int)(fCord[0] + 0.5);
//				nCord[1] = (int)(fCord[1] + 0.5);
//				nCord[2] = (int)(fCord[2] + 0.5);
//
////				imgReadRaw.GetImageInfo()->GridToIndex()
////				imgTempData.SetPixelByGreyscale( ix,iy,iz,0, imgReadVsls.GetPixelGreyscaleClosest(nCord[0],nCord[1],nCord[2],0) ) ; 
//	
//				sTempData[iz * nImWY * nImWX + iy * nImWX + ix] = \
//					sNormVsls[nCord[2] * nImWY * nImWX + nCord[1] * nImWX + nCord[0]];
//			}
//		}
//	}
//	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
//	{
////		sNormVsls[i] = sTempData[i];
//	}
//	delete[] sTempData;
//*/
//	// read model points from "*.vtk" file
//	vector< miiCNode<double, float> > vModelPoints;
//	ReadVtk(chModelMapCAFileName, vModelPoints);
//	//SetMolPointsIntensity(chModelMapCAFileName,strFileNameAtls,vModelPointsWithIntensity);//set the model-points intensity based on neighbour points intensiy.
//	// read the mapped RCA and LAD for calculating the starting point at the model
//	vector< miiCNode<double, float> > vRCAPoints;
//	ReadVtk(chRCAFileName, vRCAPoints);
//	vector< miiCNode<double, float> > vLADPoints;
//	ReadVtk(chLADFileName, vLADPoints);
//
//	float nStartPoint[3];
//	float nStartPointWithIntisity[5];
//	float nEndPoint[3];
//	// get starting point
//	nStartPoint[0] = (vRCAPoints[0].x + vLADPoints[0].x) / 2;
//	nStartPoint[1] = (vRCAPoints[0].y + vLADPoints[0].y) / 2;
//	nStartPoint[2] = (vRCAPoints[0].z + vLADPoints[0].z) / 2;
//	nStartPointWithIntisity[0]=nStartPoint[0];
//	nStartPointWithIntisity[1]=nStartPoint[1];
//	nStartPointWithIntisity[2]=nStartPoint[2];
//	nStartPointWithIntisity[3]=gAortaMean[0]; 
//	nStartPointWithIntisity[4]=gAortaStdDev[0];
//// get end point
//	nEndPoint[0] = vModelPoints[vModelPoints.size()-1].x;
//	nEndPoint[1] = vModelPoints[vModelPoints.size()-1].y;
//	nEndPoint[2] = vModelPoints[vModelPoints.size()-1].z;
//	
//	// create an instance for Model-based minimal path with the correction
//	miiMinPathModel *iMinPathModelCrct = new miiMinPathModel(nImWX, nImWY, nImWZ,nImgSpacing,\
//		imgReadRaw.GetImageInfo(),chResultPath, true);//creat a class point
//	
//	
//	// get current time(unit: s)
//	int nET = GetTime();
//	int nMaxItrNum = 2000000;
//	int nSaveItrNum=200;
//	int nItrNum=0;
//	int nItrNumSum=0;
//
//	/************** for Model-based method with the corrction *******************/
//	cout << "\nModel-based method with the correction ..." << endl;
//
//	// set the model
//	iMinPathModelCrct->SetModelPoints(vModelPoints, nStartPoint);
//	//set the model with intensity Add by JDQ
//	//iMinPathModelCrct->SetModelPointsWithIntensity(vModelPointsWithIntensity,nStartPointWithIntisity);//transfor the intensity into public model intensity vector
//	iMinPathModelCrct->SetModelPointsWithIntensity(vModelPointsWithIntensity,nStartPointWithIntisity);
//	//set the mean and stddev of the start point for atlas image
//	iMinPathModelCrct->SetMeanStdDevAtlasImg(gAortaMean[0], gAortaStdDev[0]);
//
//	//set the mean and stddev of the start point for unseen image
//	iMinPathModelCrct->SetMeanStdDevUnseenImg(gAortaMean[1], gAortaStdDev[1]);
//
//	// set the first mean and std-dev of Unseen Image
//	iMinPathModelCrct->SetMeanStdDev(fUnseenCAMean, fUnseenCAStdDev);
//
//	
//	// FMM initialization
//	//iMinPathModelCrct->NewFastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), \
//		sNormVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);
//	//iMinPathModelCrct->FastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), \
//		sNormVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);//sNormimg here 3 means ( Intensity + Vesselness + Direction)
//	iMinPathModelCrct->ModifiedFastMarchingInit(sWinImg, imgReadRaw.GetImageInfo(), \
//		sWinVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);
//	// save starting point to "txt"
//	WriteFileMC << "{" << nStartPoint[0] << ", " << nStartPoint[1] << ", " << nStartPoint[2] << "}(ModelCorrection)" << "->" << endl;
//	
//	for (int i = 0; i < nMaxSgmtNum; i++)
//	{			
//		cout << "Segment: " << i+1 << endl;
//		// get current time(unit: s)
//		nET = GetTime();
//		//Output Vector point position (Point C and Point D )
//		iMinPathModelCrct->CoutSandEndPosi(chResultPath);
//		
//		// fast-marching-method evolution
//		//nItrNum = iMinPathModelCrct->FastMarchingEvolution(nMaxItrNum,nItrNumSum);
//		//nItrNum = iMinPathModelCrct->NewFastMarchingEvolution(nMaxItrNum,nItrNumSum);
//		//nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolution(nMaxItrNum,nItrNumSum,sWinImg);
//		//nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModel(nMaxItrNum,nItrNumSum,sWinImg);
//		//nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelWithLineMask(nMaxItrNum,nItrNumSum,sWinImg,sLineMaskData);
//		nItrNum =iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNew(nMaxItrNum,nItrNumSum,sWinImg,imgReadRaw.GetImageInfo(),nStartPoint,2000);
//		//nItrNum =iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMask(nMaxItrNum,nItrNumSum,sWinImg,imgReadRaw.GetImageInfo(),nStartPoint,2000,sLineMaskData);
//		//nItrNum =iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelwithMask(nMaxItrNum,nItrNumSum,sWinImg,imgReadRaw.GetImageInfo(),nStartPoint,2000,sLineMaskData);
//		nItrNumSum=nItrNumSum+nItrNum;
//		
//		// calculate the elapsed time for FMM evolution
//		int nUseTime = GetTime() - nET;
//
//		if (nItrNum < nMaxItrNum && nItrNum > 0)
//		{
//			// find the minimal path by using back-propagation
//			cout << "Finding minimal path..." << endl;
//			if (iMinPathModelCrct->FindMinPathWithintensiy(imgReadRaw.GetImageInfo(),nStartPoint,sWinImg ,2000))//if (iMinPathModelCrct->FindMinPath(nStartPoint, 2000))
//			{
//				bRunResult = true;
//			}
//			else
//			{
//				bRunResult = false;
//				// output the elapsed time (uint: s)
//				cout << "Elapsed time: " << nUseTime << "s" << endl;
//				cout << "Finding min-path failed!" << endl;
//				break;
//			}			
//		}
//		else
//		{
//			bRunResult = false;
//			// output the elapsed time (uint: s)
//			cout << "Elapsed time: " << nUseTime << "s" << endl;
//			cout << "FMM failed!" << endl;
//			break;
//		}
//
//		// output the elapsed time (uint: s)
//		cout << "Elapsed time: " << nUseTime << "s" << endl;
//
//		// save the text results
//		cout << "Saving results..." << endl;
//		float fMCEndPointWorld[3];
//		//iMinPathModelCrct->GetEndPoint(nMCEndPoint);
//		iMinPathModelCrct->GetEndPointWorld(fMCEndPointWorld);
//		WriteFileMC << i+1 << "/" << nMaxSgmtNum << ": ";
//		WriteFileMC <<right<<fixed<<setfill('0')<<setprecision(1) \
//			<< "Pos{" << fMCEndPointWorld[0] << ", " << fMCEndPointWorld[1] << ", " << fMCEndPointWorld[2] << "}, ";
//		WriteFileMC << "Elapsed time" << "(" << nUseTime << "s), ";
//		if (bRunResult)
//			WriteFileMC << "Successful!" << endl;
//		else
//			WriteFileMC << "Failed!" << endl;
//
//		// save the image results
//
//		/*		stringstream ss; // c++中用于格式化字符串的流
//		ss << value;
//		string strTemp;
//		ss >> strTemp; // 此时strTemp变为"100"
//		string strFilenameMC = string(chResultPath);
//		strFilenameMC = strFilenameMC + "/MCImageLine" + string(chTemp) + ".nii.gz";
//		iMinPathModelCrct->GenMinPathGraph(imgReadVsls.GetImageInfo(), sVslsData, strFilenameMC, 0);
//
//		strFilenameMC = string(chResultPath);
//		strFilenameMC = strFilenameMC + "/MCLine" + string(chTemp) + ".nii.gz";
//		iMinPathModelCrct->GenMinPathGraph(imgReadVsls.GetImageInfo(), sVslsData, strFilenameMC, 1);
//		*/
//		// save U in minimal path
//		/*
//		string strFilenameU2 = string(chResultPath);
//		strFilenameU2 = strFilenameU2 + "/ModelCrctU" + string(chTemp) + ".nii.gz";
//		iMinPathModelCrct->SaveU2Image(imgReadRaw.GetImageInfo(), strFilenameU2);
//		*/
//		// save CA to vtk file
//		char chTemp[25];
//		_itoa_s(i, chTemp, 10);
//		int nLen0 = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".vtk") + 1;
//		char *chFileName0 = (char *)malloc(nLen0);
//		strcpy(chFileName0, chResultPath);
//		strcat(chFileName0, "/MCLine");
//		strcat(chFileName0, chTemp);
//		strcat(chFileName0, ".vtk");
//		iMinPathModelCrct->WriteCA2Vtk(chFileName0);
//		free(chFileName0);
//
//
//		int nLen1 = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".txt") + 1;
//		char *chFileName1 = (char *)malloc(nLen1);
//		strcpy(chFileName1, chResultPath);
//		strcat(chFileName1, "/MCLine");
//		strcat(chFileName1, chTemp);
//		strcat(chFileName1, ".txt");
//		iMinPathModelCrct->WriteCA2Txt(chFileName1);
//		free(chFileName1);
//
//		int nLen2 = strlen(chResultPath) + strlen("/MCLine_intensity") + strlen(chTemp) + strlen(".txt") + 1;
//		char *chFileName2 = (char *)malloc(nLen2);
//		strcpy(chFileName2, chResultPath);
//		strcat(chFileName2, "/MCLine_intensity");
//		strcat(chFileName2, chTemp);
//		strcat(chFileName2, ".txt");
//		iMinPathModelCrct->WriteCAIntTxt(chFileName2);
//		free(chFileName2);
//
//
//
//
//
//		/*
//		nLen = strlen(chResultPath) + strlen("/CrctModel") + strlen(chTemp) + strlen(".txt") + 1;
//		chFileName3 = (char *)malloc(nLen);
//		strcpy(chFileName3, chResultPath);
//		strcat(chFileName3, "/CrctModel");
//		strcat(chFileName3, chTemp);
//		strcat(chFileName3, ".txt");
//		iMinPathModelCrct->GetModelPoints(vModelPoints);
//		WriteVtk(vModelPoints, chFileName3);		
//		free(chFileName3);
//		*/
//		cout << "Results saved!\n" << endl;
//		/*if (!iMinPathModelCrct->ReFastMarchingInit())
//		{
//		break;
//		}*/
//		/*if (!iMinPathModelCrct->ModifiedReFastMarchingInit())
//		{
//		break;
//		}*/
//		if (!iMinPathModelCrct->ModifiedReFastMarchingInitdontMoveModel())
//		{
//			break;
//		}
//			
//		
//		/*if (!iMinPathModelCrct->NewReFastMarchingInit()||!iMinPathModelCrct->UpdateUnseenImgIntensity(sNormImg))
//		{
//		break;
//		}*/
//		// update the whole intensity parameterAdd by JDQ
//	}
//
//	delete iMinPathModelCrct;
//	WriteFileMC.close();
//
//	// release memory
//	delete[] sWinImg;
//	delete[] sWinVsls;
//
//	imgReadRaw.ReleaseMem();
//	imgReadVsls.ReleaseMem();
//
//	cout << "Run over!\n" << endl;
//	return 1;	
//
////***SimumodelExtraction**///


//***MapModelfrom training data***//
	
	//if( argc < 12 )
	//{
	//	std::cerr << "Usage: " << std::endl;
	//	std::cerr << "miiMinPathModelApp	vesselness	image	atlas_vesselx	results	atlas_vessel0	\
	//		atlas_vessel1	atlas_mean	unseen_mean	atlas_image	atlas_CA	segment_num" << std::endl;
	//	return -1;
	//}
	//set file path and name
	/*
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
	//
 //   string strFileNameRaw = string(argv[1]);
	//string strFileNameVsls = string(argv[2]);
	//
	//
	//char *chModelMapCAFileName = argv[3];
	//char *chResultPath = argv[4];

	//char *chRCAFileName = argv[5];
	//char *chLADFileName = argv[6];	

	//char *chAtlsMeanStdDevFileName = argv[7];
	//char *chUnseenMeanStdDevFileName = argv[8];
	//char *chMolPointsIntensityFileName=argv[9];
	//char *chFileNameAtlsCA = argv[10];

	//

	//int nMaxSgmtNum = atoi(argv[11]);
	//

	//string strFileNameRaw =  "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/Resmol_image02.nii";
	//string strFileNameVsls = "F:/Coronary_0/MultiResolution_results/DataMaskOutLung/WithWindoingVesselness200/SetOriIntensity/Remoeb_image02_vesselness.nii";
	//char *chModelMapCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/SetintensiyByOnePoint/Simulation19.vtk";//simu-model points(JDQ)
	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu19";
	//char *chRCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/vessel0/reference.vtk";
	//char *chLADFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset02/vessel1/reference.vtk";	
	//char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_mean_std.txt";
	//char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas02_to_atlas02_mean_std.txt";
	//char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset02/vessel0/Simu_reference/SetintensiyByOnePoint/Simuintensity19.txt";
	//int nMaxSgmtNum =30;






//for script file***
int o=0;
if( argc < 11 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << "zxhcaeDMP_1.1_MLLT	imageRaw (.nii)  vesselness	(.nii)	atlas_vesselx	results	atlas_vessel	\
				detected-ostial atlas_mean	unseen_mean	ModelPointsIntensity	segment_num  vecmaxdist" << std::endl;
		return -1;
	}
string strFileNameRaw = string(argv[1]);
string strFileNameVsls= string(argv[2]);
char *chModelMapCAFileName=argv[3];
char *chResultPath=argv[4];
char *chAtlasCAFileName =argv[5];
char *chREFFileName =argv[6];
char *chAtlsMeanStdDevFileName = argv[7];
char *chUnseenMeanStdDevFileName =argv[8];
char *chMolPointsIntensityFileName=argv[9];
int nMaxSgmtNum =atoi(argv[10]);
char *chREFEndFileName =argv[11];//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset07/vessel1/pointE.txt";


//vessel 0
//string strFileNameRaw ="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset01/Resmol_image00.nii";;
//string strFileNameVsls="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset01/Remoeb_image00_vesselness.nii";
//char *chModelMapCAFileName= "F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image01_meancenterline_v0.vtk";
//char *chResultPath= "F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResultsSM/mod01_to_unseen00_results/meanimg01v0model";
//char *chAtlasCAFileName ="F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image01_meancenterline_v0.vtk";
//char *chREFFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/rightostium_world.txt";
//char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas01_mean_std.txt";
//char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_mean_std.txt";  
//char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/mean_centerline_intensity/image01_01mline_v0_inten.txt";
//int nMaxSgmtNum =40;
//char *chREFEndFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel0/pointE.txt";



//vessel0,vessel1,vessel2
//string strFileNameRaw ="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset01/Resmol_image00.nii";;
//string strFileNameVsls="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset01/Remoeb_image00_vesselness.nii";
//char *chModelMapCAFileName= "F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image01_meancenterline_v3.vtk";
//char *chResultPath= "F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResultsSM/mod01_to_unseen00_results/meanimg01v3model";
//char *chAtlasCAFileName ="F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image01_meancenterline_v3.vtk";
//char *chREFFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/leftostium_world.txt";
//char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas01_mean_std.txt";
//char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_mean_std.txt";  
//char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/mean_centerline_intensity/image01_01mline_v3_inten.txt";
//int nMaxSgmtNum =30;
//char *chREFEndFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel3/pointE.txt";

//***for script file


//for fast test briefly***

//if( argc < 3 )
//	{
//	cerr << "usage: " << endl;
//	cerr << "zxhcaeDMP_1.1_MLL.exe	modelimage-number unseen-image-number vessel-number" << endl;
//	return -1;
//	} 
//
//char *chmodelNUM=argv[1];
//string strmodelNUM;
//int imodelnum=atoi(chmodelNUM);
//if (imodelnum<10)
//{
//	strmodelNUM="0"+string(chmodelNUM);
//}
//else
//	strmodelNUM=string("10");//more than 10;
//
//char *chunseenNUM=argv[2];
//string strunseenNUM;
//int iunseennum=atoi(chunseenNUM);
//if (iunseennum<10)
//{
//	strunseenNUM="0"+string(chunseenNUM);
//}
//else
//	strunseenNUM=string("10");//more than 10;
//
//char *chvessellNUM=argv[3];
//string strvesselNUM=string(chvessellNUM);
//int ivesselNUM=atoi(chvessellNUM);
//string strFileNameRaw ="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset"+strmodelNUM+"/vessel"+strvesselNUM+"/Resmol_image"+strunseenNUM+"_new1.nii";
//string strFileNameVsls ="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset"+strmodelNUM+"/vessel"+strvesselNUM+"/Remoeb_image"+strunseenNUM+"_vesselness_new1.nii";
//string strModelMapCAFileName ="F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image"+strmodelNUM+"_meancenterline_v"+strvesselNUM+".vtk";
//string strResultPath = "F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod"+strmodelNUM+"_to_unseen"+strunseenNUM+"_new_results1/meanimg"+strmodelNUM+"v"+strvesselNUM+"model";
//string strAtlasCAFileName,strREFFileName;
//if(ivesselNUM<1)
//{
//	strAtlasCAFileName ="F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image"+strmodelNUM+"_meancenterline_v0.vtk";
//	strREFFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset"+strunseenNUM+"/rightostium_world.txt";
//}
//else
//{
//	strAtlasCAFileName ="F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image"+strmodelNUM+"_meancenterline_v1.vtk";
//	strREFFileName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset"+strunseenNUM+"/leftostium_world.txt";
//}
//
//string strAtlsMeanStdDevFileName ="F:/Coronary_0/meanstd/atlas"+strmodelNUM+"_mean_std.txt";
//string strUnseenMeanStdDevFileName ="F:/Coronary_0/meanstd/atlas"+strmodelNUM+"_to_atlas"+strunseenNUM+"_mean_std.txt";  
//string strMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/mean_centerline_intensity/image"+strmodelNUM+"_"+strmodelNUM+"mline_v"+strvesselNUM+"_inten.txt";
//char *chModelMapCAFileName=new char[1024];
//char *chResultPath=new char[1024];
//char *chAtlasCAFileName=new char[1024];
//char *chREFFileName=new char[1024];
//char *chAtlsMeanStdDevFileName=new char[1024];
//char *chUnseenMeanStdDevFileName=new char[1024];
//char *chMolPointsIntensityFileName=new char[1024];
//stringstream sstr;
//sstr.clear();
//sstr<<strModelMapCAFileName;
//sstr>>chModelMapCAFileName;
//sstr.clear();
//sstr<<strResultPath;
//sstr>>chResultPath;
//sstr.clear();
//sstr<<strAtlasCAFileName;
//sstr>>chAtlasCAFileName;
//sstr.clear();
//sstr<<strREFFileName;
//sstr>>chREFFileName;
//sstr.clear();
//sstr<<strAtlsMeanStdDevFileName;
//sstr>>chAtlsMeanStdDevFileName;
//sstr.clear();
//sstr<<strUnseenMeanStdDevFileName;
//sstr>>chUnseenMeanStdDevFileName;
//sstr.clear();
//sstr<<strMolPointsIntensityFileName;
//sstr>>chMolPointsIntensityFileName;





int t=0;		

//***for fast test briefly

//string strFileNameRaw =  "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset01/vessel0/Resmol_image00.nii";
//string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/training/dataset01/vessel0/Remoeb_image00_vesselness.nii";
//char *chModelMapCAFileName = "F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image01_meancenterline_v0.vtk";
//char *chResultPath = "F:/Coronary_0/trainningdataZXHCAEDMP/LowResoResults/mod01_to_unseen00_resultsD/meanimg01v0model";
//char *chAtlasCAFileName = "F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image01_meancenterline_v0.vtk";//rightostium(0) while vessel0;others leftostium
////char *chLADFileName = "F:/Coronary_0/Coronary_Niessen/mean_centerline_2014_01_17/image03_meancenterline_v1.vtk";
//char *chREFFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/rightostium_world.txt";//rightostium while vessel0;others leftostium
////char *chREFLFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset01/leftostium_world.txt";
//char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas01_mean_std.txt";
//char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas01_to_atlas00_mean_std.txt";  
//char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/mean_centerline_intensity/image01_01mline_v0_inten.txt";
////string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image03.nii";
//int nMaxSgmtNum =30;	
	//***Comment
	//char *chModelMapCAFileName = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/SimuModel/dataset03/vessel2/Simu_reference/Simulation0.vtk";//simu-model points(JDQ)
	//char *chResultPath = "F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/SimuModelResult/dataset03/vessel2/Simu0";
	/*char *chRCAFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image03_vessel0.vtk";
	char *chLADFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image03_vessel1.vtk";	*/
	//char *chMolPointsIntensityFileName="F:/Coronary_0/Coronary_Niessen/ProcessByJDQ/SimuModel/dataset03/vessel2/Simu_reference/Simuintensity0.txt";
	//char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel3/reference.vtk";

	//Comment**//
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
	zxhImageDataT<short> imgReadVsls, imgReadRaw,imgReadMolRaw,imgReadLineMask; //Change by JDQ
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

	/*if( zxh::OpenImage( &imgReadLineMask, strFileNameLineMask ) == false )
	{
		std::cerr << "LineMask  image(nifti-file) is not found!"; 
		return -1;
	}*/

	// get the image data
	const short *sImData = imgReadRaw.GetImageData();
	// get the image spacing
	float nImgSpacing[4];//Add by JDQ
	imgReadRaw.GetImageSpacing(nImgSpacing[0],nImgSpacing[1],nImgSpacing[2],nImgSpacing[3] );//Add by JDQ
	// get the vesselness data 
	const short *sVslsData = imgReadVsls.GetImageData();	
	
	//get the model image data
	const short *sMolData = imgReadMolRaw.GetImageData();

	//get the line mask data
	//const short *sLineMaskData = imgReadLineMask.GetImageData();
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

		//if (sImData[i] < sRawMin)
		//	sWinImg[i] = sRawMin; //sRawMin = sNormImg[i];
		//else
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
	vector< miiCNode<double, float> > vAtlasPoints;
	ReadVtk(chAtlasCAFileName, vAtlasPoints);
	/*vector< miiCNode<double, float> > vLADPoints;
	ReadVtk(chLADFileName, vLADPoints);*/
	//read the estimated reference
	vector< miiCNode<double, float> > vREFPoints;
	ReadPointTxt(chREFFileName, vREFPoints);
	//read ref endpoints
	vector< miiCNode<double, float> > vREFEndPoints;
	ReadRefPointTxt(chREFEndFileName, vREFEndPoints);
	//ReadPointTxt(chREFLFileName, vREFLPoints);
	float nStartPoint[3];
	float nStartPointWithIntisity[5];
	float nEndPoint[3];
	float nREFStartPoint[3];
	// get starting point
	//start point is the mean value of leftostium_world and rightostium_world
	/*nStartPoint[0] = (vRCAPoints[0].x + vLADPoints[0].x) / 2;
	nStartPoint[1] = (vRCAPoints[0].y + vLADPoints[0].y) / 2;
	nStartPoint[2] = (vRCAPoints[0].z + vLADPoints[0].z) / 2;*/
	//when extracting vessel0
	/*nStartPoint[0] = vRCAPoints[0].x;
	nStartPoint[1] = vRCAPoints[0].y;
	nStartPoint[2] = vRCAPoints[0].z;*/
	// get the estimated start point of reference,vessel0
	/*nREFStartPoint[0] =vREFRPoints[0].x;
	nREFStartPoint[1] =vREFRPoints[0].y;
	nREFStartPoint[2] =vREFRPoints[0].z;*/
	//when extracting vessel1,vessel2,vessel3
	nStartPoint[0] = vAtlasPoints[0].x;
	nStartPoint[1] = vAtlasPoints[0].y;
	nStartPoint[2] = vAtlasPoints[0].z;
   // get the estimated start point of reference,vessel1,2,3
	nREFStartPoint[0] =vREFPoints[0].x;
	nREFStartPoint[1] =vREFPoints[0].y;
	nREFStartPoint[2] =vREFPoints[0].z;
	//RandomDefPoint(nREFStartPoint);//deform the refstart point randomly
	//set startpoint with intensity
	nStartPointWithIntisity[0]=nStartPoint[0];
	nStartPointWithIntisity[1]=nStartPoint[1];
	nStartPointWithIntisity[2]=nStartPoint[2];
	nStartPointWithIntisity[3]=gAortaMean[0]; 
	nStartPointWithIntisity[4]=gAortaStdDev[0];

// get end point
	nEndPoint[0] =vREFEndPoints[0].x;
	nEndPoint[1] =vREFEndPoints[0].y;
	nEndPoint[2] =vREFEndPoints[0].z;
	
	// create an instance for Model-based minimal path with the correction
	miiMinPathModel *iMinPathModelCrct = new miiMinPathModel(nImWX, nImWY, nImWZ,nImgSpacing,\
		imgReadRaw.GetImageInfo(),chResultPath,true);//creat a class point
	
	
	// get current time(unit: s)
	int nET = GetTime();
	int nMaxItrNum = 500000;
	int nMaxSumItrNum=1000000;
	int nSaveItrNum=200;
	int nItrNum=0;
	int nItrNumSum=0;

	/************** for Model-based method with the corrction *******************/
	cout << "\nModel-based method with the correction ..." << endl;

	// set the model
	//iMinPathModelCrct->SetModelPoints(vModelPoints, nStartPoint);
	iMinPathModelCrct->SetModelPointsStartPointCorrect(vModelPoints, nStartPoint,nREFStartPoint);//move start model points 
	//set the model with intensity Add by JDQ
	//iMinPathModelCrct->SetModelPointsWithIntensity(vModelPointsWithIntensity,nStartPointWithIntisity);//add the intensity into public model intensity vector
	iMinPathModelCrct->SetModelPointsWithIntensityCorrect(vModelPointsWithIntensity,nStartPointWithIntisity,nREFStartPoint);//add the intensity into public model intensity vector;and move the whole model
	//iMinPathModelCrct->SetModelPointsWithIntensity(vModelPointsWithIntensity,nStartPointWithIntisity);
	//set the mean and stddev of the start point for atlas image
	iMinPathModelCrct->SetMeanStdDevAtlasImg(gAortaMean[0], gAortaStdDev[0]);

	//set the mean and stddev of the start point for unseen image
	iMinPathModelCrct->SetMeanStdDevUnseenImg(gAortaMean[1], gAortaStdDev[1]);

	// set the first mean and std-dev of Unseen Image
	iMinPathModelCrct->SetMeanStdDev(fUnseenCAMean, fUnseenCAStdDev);

	nStartPoint[0] =nREFStartPoint[0];
	nStartPoint[1] =nREFStartPoint[1];
	nStartPoint[2] =nREFStartPoint[2];

	// FMM initialization
	//iMinPathModelCrct->NewFastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), \
		sNormVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);
	//iMinPathModelCrct->FastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), \
		sNormVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);//sNormimg here 3 means ( Intensity + Vesselness + Direction)
	iMinPathModelCrct->ModifiedFastMarchingInitLLSM(sWinImg, imgReadRaw.GetImageInfo(), \
		sWinVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);
	// save starting point to "txt"
	WriteFileMC << "{" << nStartPoint[0] << ", " << nStartPoint[1] << ", " << nStartPoint[2] << "}(ModelCorrection)" << "->" << endl;
	
	for (int i = 0; i < nMaxSgmtNum; i++)
	{			
		cout << "Segment: " << i+1 << endl;
		// get current time(unit: s)
		nET = GetTime();
		//Output Vector point position (Point C and Point D )
		//iMinPathModelCrct->CoutSandEndPosi(chResultPath);
		
		// fast-marching-method evolution
		//nItrNum = iMinPathModelCrct->FastMarchingEvolution(nMaxItrNum,nItrNumSum);
		//nItrNum = iMinPathModelCrct->NewFastMarchingEvolution(nMaxItrNum,nItrNumSum);
		//nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolution(nMaxItrNum,nItrNumSum,sWinImg);
		//nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModel(nMaxItrNum,nItrNumSum,sWinImg);
		//nItrNum = iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelWithLineMask(nMaxItrNum,nItrNumSum,sWinImg,sLineMaskData);
		nItrNum =iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLLSM(nMaxItrNum,nItrNumSum,sWinImg,imgReadRaw.GetImageInfo(),nStartPoint,2000);
		//nItrNum =iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMask(nMaxItrNum,nItrNumSum,sWinImg,imgReadRaw.GetImageInfo(),nStartPoint,2000,sLineMaskData);
		//nItrNum =iMinPathModelCrct->ModifiedFastMarchingEvolutiondontMoveModelwithMask(nMaxItrNum,nItrNumSum,sWinImg,imgReadRaw.GetImageInfo(),nStartPoint,2000,sLineMaskData);
		nItrNumSum=nItrNumSum+nItrNum;
		
		// calculate the elapsed time for FMM evolution
		int nUseTime = GetTime() - nET;

		if (nItrNum < nMaxItrNum && nItrNum > 0&&nItrNumSum<nMaxSumItrNum)
		{
			// find the minimal path by using back-propagation
			cout << "Finding minimal path..." << endl;
			if (iMinPathModelCrct->FindMinPathWithintensiyLL(imgReadRaw.GetImageInfo(),nStartPoint,sWinImg ,2000))//if (iMinPathModelCrct->FindMinPath(nStartPoint, 2000))
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
		int nLen0 = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".vtk") + 1;
		char *chFileName0 = (char *)malloc(nLen0);
		strcpy(chFileName0, chResultPath);
		strcat(chFileName0, "/MCLine");
		strcat(chFileName0, chTemp);
		strcat(chFileName0, ".vtk");
		iMinPathModelCrct->WriteCA2Vtk(chFileName0);
		free(chFileName0);


		int nLen1 = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".txt") + 1;
		char *chFileName1 = (char *)malloc(nLen1);
		strcpy(chFileName1, chResultPath);
		strcat(chFileName1, "/MCLine");
		strcat(chFileName1, chTemp);
		strcat(chFileName1, ".txt");
		iMinPathModelCrct->WriteCA2Txt(chFileName1);
		free(chFileName1);

		int nLen2 = strlen(chResultPath) + strlen("/MCLine_intensity") + strlen(chTemp) + strlen(".txt") + 1;
		char *chFileName2 = (char *)malloc(nLen2);
		strcpy(chFileName2, chResultPath);
		strcat(chFileName2, "/MCLine_intensity");
		strcat(chFileName2, chTemp);
		strcat(chFileName2, ".txt");
		iMinPathModelCrct->WriteCAIntTxt(chFileName2);
		free(chFileName2);





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
		if (!iMinPathModelCrct->ModifiedReFastMarchingInitdontMoveModelSM())
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
//***MapModelfrom training data***//
}
//int main(int argc, char* argv[])
//{
//	zxh::echo_zxh(argc,argv);
//	if( argc == 1 )
//	{
//		std::cout<<"  zxhcaeDMP [options]        \n";
//		Help();
//		return 0 ;
//	}
//	if( argc == 2 && strcmp( argv[1],"-H" )==0 )
//	{
//		std::cout<<"  zxhcaeDMP [options]        \n";
//		HELP();
//		return -1 ;
//	}
//	if( glbVerboseOutput>0 )
//	{
//		std::cout<<"\n * \n * zxhcaeDMP, version of 1.0  \n * \n";
//		zxh::echo_arguments( argc, argv ) ;
//	} 
//	return zxhcaeDMP_main(argc,argv);
//} 

int GetTime()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}


 