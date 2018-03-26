
#include "miiMinPathOrg.h"

#include <time.h>

int GetTime()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}

int main(int argc, char *argv[])
{
	if( argc < 9 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << "miiMinPathOrgApp	image	vesselness start_x	start_y	start_z	\
					 end_x	end_y	end_z	result_dir" << std::endl;
		return -1;
	}

	string strFileNameRaw = string(argv[1]);
	string strFileNameVsls = string(argv[2]);
	
	float nStartPoint[3];
	nStartPoint[0] = atof(argv[3]);
	nStartPoint[1] = atof(argv[4]);
	nStartPoint[2] = atof(argv[5]);

	float nEndPoint[3];
	nEndPoint[0] = atof(argv[6]);
	nEndPoint[1] = atof(argv[7]);
	nEndPoint[2] = atof(argv[8]);

	char *chResultPath = argv[9];
	//set file path and name
// 	string strFileNameVsls = "E:/UI_Projects/lib/coronary/training/dataset00/image00_vesselness_ort.nii";
// 	string strFileNameRaw = "E:/UI_Projects/lib/coronary/training/dataset00/image00.nii";

	//read a "nifti" file
	zxhImageDataT<short> imgReadVsls, imgReadRaw; 
	if( zxh::OpenImage( &imgReadRaw, strFileNameRaw ) == false )
	{
		std::cerr << " Raw image(nifti-file) is not found!"; 
		return -1;
	}

	if( zxh::OpenImage( &imgReadVsls, strFileNameVsls ) == false )
	{
		std::cerr << "Vesselness image(nifti-file) is not found!"; 
		return -1;
	}

	// get the image data
	const short *sImData = imgReadRaw.GetImageData();

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
			sRawMin = sNormImg[i];
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

	// create an instance for minimal path 
	miiMinPathOrg *iMinPath = new miiMinPathOrg(nImWX, nImWY, nImWZ, imgReadRaw.GetImageInfo());

	// fast-marching-method initiation
	iMinPath->FastMarchingInit(sNormImg, imgReadRaw.GetImageInfo(), nStartPoint, nEndPoint, 1);

	// get current time (unit: s)
	int nET = GetTime();

	int nMaxItrNum = 20000000;

	cout << "Finding the path ..." << endl;

	// fast-marching-method evolution
	int nItrNum = iMinPath->FastMarchingEvolution(nMaxItrNum);

	// calculate the elapsed time for FMM evolution
	int nUseTime = GetTime() - nET;

	// output the elapsed time (uint: s)
	cout << "elapsed time:" << nUseTime << "s" << endl;

	if (nItrNum < nMaxItrNum)
	{
		// find the minimal path by using back-propagation
		if (iMinPath->FindMinPath(5000))
		{
			int nLen = strlen(chResultPath) + strlen("/Line") + strlen(".vtk") + 1;
			char *chVtkFileName = (char *)malloc(nLen);
			strcpy(chVtkFileName, chResultPath);
			strcat(chVtkFileName, "/Line");
			strcat(chVtkFileName, ".vtk");
			iMinPath->WriteCA2Vtk(chVtkFileName);

			// save the results
//			iMinPath->GenMinPathGraph(imgReadVsls.GetImageInfo(), sVslsData, "ImageLine.nii.gz", 0);

			nLen = strlen(chResultPath) + strlen("/Line") + strlen(".nii.gz") + 1;
			chVtkFileName = (char *)malloc(nLen);
			strcpy(chVtkFileName, chResultPath);
			strcat(chVtkFileName, "/Line");
			strcat(chVtkFileName, ".nii.gz");
			iMinPath->WriteCA2Vtk(chVtkFileName);
			iMinPath->GenMinPathGraph(imgReadVsls.GetImageInfo(), sVslsData, chVtkFileName, 1);

			cout << "Successful!" << endl;
		}
		else
		{
			cout << "Failed!" << endl;
		}
	}
	else
	{
		cout << "Failed!" << endl;
	}

	// release memory	
	imgReadRaw.ReleaseMem();
	imgReadVsls.ReleaseMem();
	delete iMinPath;
	delete[] sNormImg;
	delete[] sNormVsls;

	return 1;	
}
