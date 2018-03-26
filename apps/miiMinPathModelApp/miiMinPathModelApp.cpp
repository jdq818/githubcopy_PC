
#include "miiMinPathModel.h"

#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkPolyLine.h"

#include <time.h>

#define TAB_CHAR	9

int GetTime()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}

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

int main(int argc, char *argv[])
{/**/
	/*
	if( argc < 10 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << "miiMinPathModelApp	vesselness	image	atlas_vesselx	results	atlas_vessel0	\
			atlas_vessel1	atlas_mean	unseen_mean	atlas_image	atlas_CA	max_segment	segment_num" << std::endl;
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

	int nMaxSgmtNum = 5;
	
*/
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

	int nMaxSgmtDist = atoi(argv[11]);

	int nMaxSgmtNum = atoi(argv[12]);
	*/
	string strFileNameVsls = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset01/image01_vesselness.nii"; 
	string strFileNameRaw = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset01/image01.nii";

	char *chModelMapCAFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image01_vessel0.vtk";//model points(JDQ)
	//char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen01_results";
	char *chResultPath = "F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen01_results_phys";
	char *chRCAFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image01_vessel0.vtk";
	char *chLADFileName = "F:/Coronary_0/Exp4_WhLabelReg/atlas00_to_image01_vessel1.vtk";	

	char *chAtlsMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_mean_std.txt";
	char *chUnseenMeanStdDevFileName = "F:/Coronary_0/meanstd/atlas00_to_atlas01_mean_std.txt";

	string strFileNameAtls = "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	char *chFileNameAtlsCA = "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel0/reference.vtk";

	int nMaxSgmtDist = 100;
	
	int nMaxSgmtNum = 10;

	// define variables for Mean and Std-Deviation
	

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
	CalcUnseenMeanStdDev(gAortaMean, gAortaStdDev, fAtlsCAMean, fAtlsCAStdDev, fUnseenCAMean, fUnseenCAStdDev);//calculate the mean and stddev of the useen image

	// for Model with the correction
	string strMCFilename = string(chResultPath);
	strMCFilename = strMCFilename + "/" + string("MethodMC.txt");
	ofstream WriteFileMC(strMCFilename);

	bool bRunResult;

	//read a "nifti" file
	zxhImageDataT<short> imgReadVsls, imgReadRaw; 
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
	
	//get the image size
	int nImWX, nImWY, nImWZ, nImWT;
	imgReadRaw.GetImageSize(nImWX, nImWY, nImWZ, nImWT);

	// normalize
	short *sNormImg = new short[nImWX * nImWY * nImWZ];
	short *sNormVsls = new short[nImWX * nImWY * nImWZ];
	short sRawMin = 0, sRawMax = 0;
	short sVslsMin = 0, sVslsMax = 50;

	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	{
		if (sVslsData[i] > sVslsMax)
			sNormVsls[i] = sVslsMax;
		else
			sNormVsls[i] = sVslsData[i];

		sNormImg[i] = sImData[i];

		if (sNormImg[i] < sRawMin)
			sNormImg[i] = sRawMin; //sRawMin = sNormImg[i];
		if (sNormImg[i] > sRawMax)
			sRawMax = sNormImg[i];

		if (sNormVsls[i] < sVslsMin)
			sVslsMin = sNormVsls[i];
		if (sNormVsls[i] > sVslsMax)
			sVslsMax = sNormVsls[i];
	}	

	for (int i = 0; i < nImWX * nImWY * nImWZ; i++)
	{
		sNormImg[i] = 1024 * ((double)(sNormImg[i] - sRawMin) / (sRawMax - sRawMin));
		sNormVsls[i] = 1024 * ((double)(sNormVsls[i] - sVslsMin) / (sVslsMax - sVslsMin));		
	}

	fUnseenCAMean = 1024 * (fUnseenCAMean / (sRawMax - sRawMin));
	fUnseenCAStdDev = 1024 * (fUnseenCAStdDev / (sRawMax - sRawMin));
	
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

	// read the mapped RCA and LAD for calculating the starting point at the model
	vector< miiCNode<double, float> > vRCAPoints;
	ReadVtk(chRCAFileName, vRCAPoints);
	vector< miiCNode<double, float> > vLADPoints;
	ReadVtk(chLADFileName, vLADPoints);

	float nStartPoint[3];
	float nEndPoint[3];
	// get starting point
	nStartPoint[0] = (vRCAPoints[0].x + vLADPoints[0].x) / 2;
	nStartPoint[1] = (vRCAPoints[0].y + vLADPoints[0].y) / 2;
	nStartPoint[2] = (vRCAPoints[0].z + vLADPoints[0].z) / 2;

	// create an instance for Model-based minimal path with the correction
	miiMinPathModel *iMinPathModelCrct = new miiMinPathModel(nImWX, nImWY, nImWZ,nImgSpacing,\
		imgReadRaw.GetImageInfo(),chResultPath, true);//creat a class point

	// get end point
	nEndPoint[0] = vModelPoints[vModelPoints.size()-1].x;
	nEndPoint[1] = vModelPoints[vModelPoints.size()-1].y;
	nEndPoint[2] = vModelPoints[vModelPoints.size()-1].z;

	// get current time(unit: s)
	int nET = GetTime();
	int nMaxItrNum = 2000000;
	int nItrNum;

	/************** for Model-based method with the corrction *******************/
	cout << "\nModel-based method with the correction ..." << endl;

	// set the model
	iMinPathModelCrct->SetModelPoints(vModelPoints, nStartPoint);

	// set the mean and std-dev
	iMinPathModelCrct->SetMeanStdDev(fUnseenCAMean, fUnseenCAStdDev);

	// FMM initialization
	iMinPathModelCrct->FastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), \
		sNormVsls, imgReadVsls.GetImageInfo(), nStartPoint, nEndPoint, 3);//sNormimg here 3 means ( Intensity + Vesselness + Direction)

	// save starting point to "txt"
	WriteFileMC << "{" << nStartPoint[0] << ", " << nStartPoint[1] << ", " << nStartPoint[2] << "}(ModelCorrection)" << "->" << endl;
	
	for (int i = 0; i < nMaxSgmtNum; i++)
	{			
		cout << "Segment: " << i+1 << endl;
		// get current time(unit: s)
		nET = GetTime();
		// fast-marching-method evolution
		nItrNum = iMinPathModelCrct->FastMarchingEvolution(nMaxItrNum);

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
		float nMCEndPoint[3];
		iMinPathModelCrct->GetEndPoint(nMCEndPoint);
		WriteFileMC << i+1 << "/" << nMaxSgmtNum << ": ";
		WriteFileMC <<right<<fixed<<setfill('0')<<setprecision(1) \
			<< "Pos{" << nMCEndPoint[0] << ", " << nMCEndPoint[1] << ", " << nMCEndPoint[2] << "}, ";
		WriteFileMC << "Elapsed time" << "(" << nUseTime << "s), ";
		if (bRunResult)
			WriteFileMC << "Successful!" << endl;
		else
			WriteFileMC << "Failed!" << endl;

		// save the image results
		char chTemp[2];
		_itoa_s(i, chTemp, 10);
/*		
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
		int nLen = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".vtk") + 1;
		char *chFileName2 = (char *)malloc(nLen);
		strcpy(chFileName2, chResultPath);
		strcat(chFileName2, "/MCLine");
		strcat(chFileName2, chTemp);
		strcat(chFileName2, ".vtk");
		iMinPathModelCrct->WriteCA2Vtk(chFileName2);
		free(chFileName2);
/*
		nLen = strlen(chResultPath) + strlen("/CrctModel") + strlen(chTemp) + strlen(".vtk") + 1;
		chFileName2 = (char *)malloc(nLen);
		strcpy(chFileName2, chResultPath);
		strcat(chFileName2, "/CrctModel");
		strcat(chFileName2, chTemp);
		strcat(chFileName2, ".vtk");
		iMinPathModelCrct->GetModelPoints(vModelPoints);
		WriteVtk(vModelPoints, chFileName2);		
		free(chFileName2);
*/
		cout << "Results saved!\n" << endl;

		if (!iMinPathModelCrct->ReFastMarchingInit())
		{
			break;
		}
	}

	delete iMinPathModelCrct;
	WriteFileMC.close();
	
	// release memory
	delete[] sNormImg;
	delete[] sNormVsls;
	imgReadRaw.ReleaseMem();
	imgReadVsls.ReleaseMem();

	cout << "Run over!\n" << endl;
	return 1;	
}
