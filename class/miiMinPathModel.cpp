/*******************************************************************
*
*    DESCRIPTION:  
*
*    AUTHOR: L. L.
*
*    HISTORY: 1.0
*
*    DATE: 2012-12-11
*
*******************************************************************/

/************************************
* upwind   
-------------------%%%%%++++++
-----------------%%%%%++++++++
---------------%%%%%++++++++++
-------------%%%%%++++++++++++
-----------%%%%%++++++++++++++
---------%%%%%++++++++++++++++
-------%%%%%++++++++++++++++++
-----%%%%%++++++++++++++++++++
* notes: 
* '-' --> 'alive', set 0
* '%' --> 'trial', set 1, Narrow Band
* '+' --> 'far', set 2
*
************************************/

/** include files **/
#include "miiMinPathModel.h"

//
miiMinPathModel::miiMinPathModel(int nImgWX, int nImgWY, int nImgWZ,float nImgSpacing[4],const zxhImageInfo *pBaseImgInfo, char *chResultPath, \
	 bool bReadMeanFlg):m_bReadMeanFlg(bReadMeanFlg), m_bFistEvlFlg(false),iMinDisBetweenVec(3), iMaxDisBetweenVec(12),miiMinPath(nImgWX,nImgWY,nImgWZ,nImgSpacing[0])
	//miiMinPathModel::miiMinPathModel(int nImgWX, int nImgWY, int nImgWZ, const zxhImageInfo *pBaseImgInfo, int nMaxSgmtDist, bool bReadMeanFlg):m_bReadMeanFlg(bReadMeanFlg), m_nModelPos(0), m_nFMMEvlNum(0), m_bFistEvlFlg(false)//Change by JDQ
{
	m_nImgWX = nImgWX;
	m_nImgWY = nImgWY;
	m_nImgWZ = nImgWZ;
	m_nImgSpacing[0]=nImgSpacing[0];//add by JDQ
	m_nImgSpacing[1]=nImgSpacing[1];//add by JDQ
	m_nImgSpacing[2]=nImgSpacing[2];//add by JDQ
	m_nImgSpacing[3]=nImgSpacing[3];//add by JDQ
	m_pBaseImgInfo = pBaseImgInfo;
	m_zxhBadPointImg.NewImage( m_pBaseImgInfo) ; 
	m_zxhMinNodeImg.NewImage( m_pBaseImgInfo);
	m_zxhBadPointMaskImg.NewImage( m_pBaseImgInfo);
	m_zxhNBUvalueImg.NewImage( m_pBaseImgInfo);
	m_zxhNBSpeedvalueImg.NewImage( m_pBaseImgInfo);
	m_zxhNBPvalueImg.NewImage( m_pBaseImgInfo);
	m_zxhNBSvalueImg.NewImage( m_pBaseImgInfo);
	m_zxhNBDvalueImg.NewImage( m_pBaseImgInfo);
	m_zxhNBD_SegvalueImg.NewImage( m_pBaseImgInfo);
	m_dU = new double[nImgWX * nImgWY * nImgWZ];
	m_nFmMap = new int[nImgWX * nImgWY * nImgWZ];
	m_dP2 = new double[nImgWX * nImgWY * nImgWZ];


	// create a min-heap for Narrow Band
	m_iMinHeap = new miiMinHeap<double>();	

	//m_nMaxDist2 = nMaxSgmtDist * nMaxSgmtDist;
	m_nMinDist2=20*20;//add by JDQ
	m_nMinDistmm2=1*1;
	m_nEndDistmm=1;//add by JDQ
	m_nKPCosTheta=0.8660;//add by JDQ
	//
	m_nBadpointSum=0;
	m_fModelVecWorld[3]=(0,0,0);
	m_fVesselVecWorld[3]=(0,0,0);
	m_iFarPValue=0;
	m_iAlivePValue=50;
	m_iActivePValue=20;
	m_iBadPointValue=200;
	m_iMinimalUPValue=100;
	m_nSaveFmItr=100000;
	m_cResultPath=chResultPath;
	IsEndPoint=false;
	IsCurEndPoint=false;
	bIsCurSMEndPoint=false;
	///

}

//
miiMinPathModel::~miiMinPathModel()
{
}

// Function Name: FastMarchingInit()
//
// Parameters: *sImgData: the pointer for the data of 3D image(normalized)
//			   *sVslsData: the pointer for the vesselness results				
//				nStartPoint[]: the start point 
//				nEndPoint[]: the end point
//				nMethod: the option for the potential function  
//						'0' - Gradient
//						'1' - Intensity + Vesselness
//						'2' - Intensity + Vesselness + Direction
//				*pImageInfo: the image info from nifti-file
//
// Description: initiate the FMM's data, including Narrow Band, distance function, and FMM map.  
//
// Returns: 
//
void miiMinPathModel::FastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, \
	float nStartPoint[], float nEndPoint[], int nMethod)
{
	// FMM initialization
	FastMarchingInitBase(sImgData, nStartPoint, nEndPoint);

	m_iSgmtPoint.x = m_iStartPoint.x;
	m_iSgmtPoint.y = m_iStartPoint.y;
	m_iSgmtPoint.z = m_iStartPoint.z;

	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iModelEndPoint.x = m_iEndPoint.x;
	m_iModelEndPoint.y = m_iEndPoint.y;
	m_iModelEndPoint.z = m_iEndPoint.z;	

	m_iSgmtPointWorld.x = m_iStartPointWorld.x;
	m_iSgmtPointWorld.y = m_iStartPointWorld.y;
	m_iSgmtPointWorld.z = m_iStartPointWorld.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iModelEndPointWorld.x = m_iEndPointWorld.x;
	m_iModelEndPointWorld.y = m_iEndPointWorld.y;
	m_iModelEndPointWorld.z = m_iEndPointWorld.z;	
	// potential function initialization
	PotentialFunction(sImgData, pImgInfo, sVslsData, pVslsInfo, nMethod);

	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return;
	}
}
void miiMinPathModel::NewFastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, \
	float nStartPoint[], float nEndPoint[], int nMethod)
{
	// FMM initialization
	FastMarchingInitBase(sImgData, nStartPoint, nEndPoint);

	m_iSgmtPoint.x = m_iStartPoint.x;
	m_iSgmtPoint.y = m_iStartPoint.y;
	m_iSgmtPoint.z = m_iStartPoint.z;

	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iModelEndPoint.x = m_iEndPoint.x;
	m_iModelEndPoint.y = m_iEndPoint.y;
	m_iModelEndPoint.z = m_iEndPoint.z;	

	m_iSgmtPointWorld.x = m_iStartPointWorld.x;
	m_iSgmtPointWorld.y = m_iStartPointWorld.y;
	m_iSgmtPointWorld.z = m_iStartPointWorld.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iModelEndPointWorld.x = m_iEndPointWorld.x;
	m_iModelEndPointWorld.y = m_iEndPointWorld.y;
	m_iModelEndPointWorld.z = m_iEndPointWorld.z;	
	

	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return;
	}
	VesselnessInit(pImgInfo, sVslsData, pVslsInfo, m_sVes,m_sNormVes);//Add by JDQ
	UnseenImgIntensitySet(sImgData,m_sUnseenImgIntensity);//Add by JDQ
}
void miiMinPathModel::ModifiedFastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, \
	float nStartPoint[], float nEndPoint[], int nMethod)
{
	// FMM initialization
	FastMarchingInitBase(sImgData, nStartPoint, nEndPoint);

	m_iSgmtPoint.x = m_iStartPoint.x;
	m_iSgmtPoint.y = m_iStartPoint.y;
	m_iSgmtPoint.z = m_iStartPoint.z;

	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iModelEndPoint.x = m_iEndPoint.x;
	m_iModelEndPoint.y = m_iEndPoint.y;
	m_iModelEndPoint.z = m_iEndPoint.z;	

	m_iSgmtPointWorld.x = m_iStartPointWorld.x;
	m_iSgmtPointWorld.y = m_iStartPointWorld.y;
	m_iSgmtPointWorld.z = m_iStartPointWorld.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iModelEndPointWorld.x = m_iEndPointWorld.x;
	m_iModelEndPointWorld.y = m_iEndPointWorld.y;
	m_iModelEndPointWorld.z = m_iEndPointWorld.z;	


	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return;
	}
	VesselnessInit(pImgInfo, sVslsData, pVslsInfo, m_sVes,m_sNormVes);//Add by JDQ
}
void miiMinPathModel::ModifiedFastMarchingInitLL(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, \
	float nStartPoint[], float nEndPoint[], int nMethod)
{
	// FMM initialization
	FastMarchingInitBase(sImgData, nStartPoint, nEndPoint);

	m_iSgmtPoint.x = m_iStartPoint.x;
	m_iSgmtPoint.y = m_iStartPoint.y;
	m_iSgmtPoint.z = m_iStartPoint.z;

	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iModelEndPoint.x = m_iEndPoint.x;
	m_iModelEndPoint.y = m_iEndPoint.y;
	m_iModelEndPoint.z = m_iEndPoint.z;	

	m_iSgmtPointWorld.x = m_iStartPointWorld.x;
	m_iSgmtPointWorld.y = m_iStartPointWorld.y;
	m_iSgmtPointWorld.z = m_iStartPointWorld.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iModelEndPointWorld.x = m_iEndPointWorld.x;
	m_iModelEndPointWorld.y = m_iEndPointWorld.y;
	m_iModelEndPointWorld.z = m_iEndPointWorld.z;	


	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return;
	}
	m_fCurrModelCArchLengthD = 0 ; 
	for (int i = 1; i <=m_nOneSegModelVectorEndPos ; i++)
	{ 
		m_fCurrModelCArchLengthD += sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1])); 
	}
	VesselnessInit(pImgInfo, sVslsData, pVslsInfo, m_sVes,m_sNormVes);//Add by JDQ

}
void miiMinPathModel::ModifiedFastMarchingInitLLSM(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, \
	float nStartPoint[], float nEndPoint[], int nMethod)
{
	// FMM initialization
	FastMarchingInitBaseSM(sImgData, nStartPoint, nEndPoint);

	m_iSgmtPoint.x = m_iStartPoint.x;
	m_iSgmtPoint.y = m_iStartPoint.y;
	m_iSgmtPoint.z = m_iStartPoint.z;

	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iModelEndPoint.x = m_iEndPoint.x;
	m_iModelEndPoint.y = m_iEndPoint.y;
	m_iModelEndPoint.z = m_iEndPoint.z;	

	m_iSgmtPointWorld.x = m_iStartPointWorld.x;
	m_iSgmtPointWorld.y = m_iStartPointWorld.y;
	m_iSgmtPointWorld.z = m_iStartPointWorld.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iModelEndPointWorld.x = m_iEndPointWorld.x;
	m_iModelEndPointWorld.y = m_iEndPointWorld.y;
	m_iModelEndPointWorld.z = m_iEndPointWorld.z;	


	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return;
	}
	VesselnessInit(pImgInfo, sVslsData, pVslsInfo, m_sVes,m_sNormVes);//Add by JDQ

}
// Function Name: ReFastMarchingInit()
//
// Parameters: bCrctFlag: the flag for the model correction
//						'true' - with the correction
//						'false' - without the correction
//
// Description: get the segmented point  
//
// Returns: 
//
bool miiMinPathModel::ReFastMarchingInit(bool bCrctFlag)
{
	// correct the CA model
	/*	
	if (bCrctFlag)
	{
		CorrectModelPoint();
	}
	*/
	if (bCrctFlag)
	{
		CorrectModelPointBasedSegLength();
	}
	
	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return false;
	}


	return true;
}
bool miiMinPathModel::ModifiedReFastMarchingInit(bool bCrctFlag)
{
	// correct the CA model
	/*	
	if (bCrctFlag)
	{
		CorrectModelPoint();
	}
	*/
	if (bCrctFlag)
	{
		CorrectModelPointBasedSegLength();
	}
	int ModelVectorStartPosi=m_nModelPos;
	if(m_nFMMEvlNum == 1)
	{
	m_nRealOriModelPointsPosi=m_nModelPos;
	m_nRealOriSegModelPointsWorld=m_vModelPointsWorld[m_nModelPos];
	
	}//m_nRealOriSegPointsWorld is fixed
	m_nRealOriModelPointsWorld=m_vModelPointsWorld[m_nRealOriModelPointsPosi];//Set the real start points of model.m_nRealOriModelPointsWorld is changing after one model correction.
	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return false;
	}
	int ModelVectorEndPosi=m_nModelPos;

	UpdateMeanStdOfUnseanImg( ModelVectorStartPosi, ModelVectorEndPosi);//put the model points into model image to calculate the meanstd of next segmentation add by JDQ 
	//ModifiedUpdateMeanStdOfUnseanImg( ModelVectorStartPosi, ModelVectorEndPosi);//read the mean intensity from the intensity.txt and std of unseen image is the gAortaStdDev[1] add by JDQ
	//ModifiedUpdateMeanStdOfUnseanImgFB( ModelVectorStartPosi, ModelVectorEndPosi);
	return true;
}
bool miiMinPathModel::ModifiedReFastMarchingInitdontMoveModel(bool bCrctFlag)
{
	
	
	if( IsEndPoint== true )//vector CD is the last model vector
	{
		cout << "Model vector has reached the model end-point!\n" << endl;
		ofstream PosiFile;
		int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt");
		char *chFileName = (char *)malloc(nLen);
		strcpy(chFileName, m_cResultPath);
		strcat(chFileName, "/Posi");
		strcat(chFileName, ".txt");
		PosiFile.open(chFileName,ios::app);
		if(!PosiFile)
		{cerr<<"open error!"<<endl;
		exit(1);
		}
		PosiFile<<"Model vector has reached the model end-point!"<<"\n";
		PosiFile.close();
		return false ; 
	}
	//if( IsCurEndPoint== true )//current sgment point is near from model end-point
	//{
	//	cout << "Current point has reached the model end-point!\n" << endl;
	//	ofstream PosiFile;
	//	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt");
	//	char *chFileName = (char *)malloc(nLen);
	//	strcpy(chFileName, m_cResultPath);
	//	strcat(chFileName, "/Posi");
	//	strcat(chFileName, ".txt");
	//	PosiFile.open(chFileName,ios::app);
	//	if(!PosiFile)
	//	{cerr<<"open error!"<<endl;
	//	exit(1);
	//	}
	//	PosiFile<<"Current point has reached the model end-point!"<<"\n";
	//	PosiFile.close();
	//	return false ; 
	//}
	

	if (bCrctFlag)
	{
		CorrectModelPointBasedSegLengthdontMoveModel();//set m_nOneSegModelVectorStartPos
	}


	if(m_nFMMEvlNum == 1)
	{
		m_nRealOriModelPointsPosi=m_nOneSegModelVectorStartPos;
		m_nRealOriSegModelPointsWorld=m_vModelPointsWorld[m_nOneSegModelVectorStartPos];
		m_nRealOriSegPointWorld=m_iSgmtPointWorld;
	}//m_nRealOriSegPointsWorld is fixed
	m_nRealOriModelPointsWorld=m_vModelPointsWorld[m_nRealOriModelPointsPosi];//Set the real start points of model.m_nRealOriModelPointsWorld is changing after one model correction.


	// update 'm_nModelPos'
	if (!UpdateModelSgmtPointNew())//set m_nOneSegModelVectorEndPos
	{
		return false;
	}
	
	//	if (!UpdateModelSgmtPointCD())//set m_nOneSegModelVectorEndPos
	//{
	//	return false;
	//}

	if(m_nOneSegModelVectorStartPos==m_vModelPointsWorld.size()-1)//if point C is the end point in model points vector,exchange the start and end.
	{
		int temp=m_nOneSegModelVectorEndPos;
		m_nOneSegModelVectorEndPos=m_nOneSegModelVectorStartPos;
		m_nOneSegModelVectorStartPos=temp;
		IsEndPoint=true;
	}

	m_fCurrModelCArchLength = 0 ; 
	for (int i = 1; i <=m_nOneSegModelVectorStartPos ; i++)
	{ 
		m_fCurrModelCArchLength += sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1])); 
	}

	m_fCurrModelCArchLengthD = 0 ; 
	for (int i = 1; i <=m_nOneSegModelVectorEndPos ; i++)
	{ 
		m_fCurrModelCArchLengthD += sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1])); 
	}
	if( !IsEndPoint== true )
	{

		//***Simu-Model***///
		//UpdateMeanStdOfUnseanImg( m_nOneSegModelVectorStartPos, m_nOneSegModelVectorEndPos);//put the model points into model image to calculate the meanstd of next segmentation add by JDQ 
				//***Simu-Model***///
		//ModifiedUpdateMeanStdOfUnseanImg( ModelVectorStartPosi, ModelVectorEndPosi);//read the mean intensity from the intensity.txt and std of unseen image is the gAortaStdDev[1] add by JDQ
		//ModifiedUpdateMeanStdOfUnseanImgFB( ModelVectorStartPosi, ModelVectorEndPosi);
	}

	return true;
}

	bool miiMinPathModel::DFMReFastMarchingWOBP(bool bCrctFlag)
{
	
	if( IsEndPoint== true )//vector CD is the last model vector
	{
		cout << "Model vector has reached the model end-point!\n" << endl;
		ofstream PosiFile;
		int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt");
		char *chFileName = (char *)malloc(nLen);
		strcpy(chFileName, m_cResultPath);
		strcat(chFileName, "/Posi");
		strcat(chFileName, ".txt");
		PosiFile.open(chFileName,ios::app);
		if(!PosiFile)
		{cerr<<"open error!"<<endl;
		exit(1);
		}
		PosiFile<<"Model vector has reached the model end-point!"<<"\n";
		PosiFile.close();
		return false ; 
	}
	//if( IsCurEndPoint== true )//current sgment point is near from model end-point
	//{
	//	cout << "Current point has reached the model end-point!\n" << endl;
	//	ofstream PosiFile;
	//	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt");
	//	char *chFileName = (char *)malloc(nLen);
	//	strcpy(chFileName, m_cResultPath);
	//	strcat(chFileName, "/Posi");
	//	strcat(chFileName, ".txt");
	//	PosiFile.open(chFileName,ios::app);
	//	if(!PosiFile)
	//	{cerr<<"open error!"<<endl;
	//	exit(1);
	//	}
	//	PosiFile<<"Current point has reached the model end-point!"<<"\n";
	//	PosiFile.close();
	//	return false ; 
	//}
	

	if (bCrctFlag)
	{
		DFM_CorrModel_WOBP();//set m_nOneSegModelVectorStartPos
	}


	if(m_nFMMEvlNum == 1)
	{
		m_nRealOriModelPointsPosi=m_nOneSegModelVectorStartPos;
		m_nRealOriSegModelPointsWorld=m_vModelPointsWorld[m_nOneSegModelVectorStartPos];
		m_nRealOriSegPointWorld=m_iSgmtPointWorld;
	}//m_nRealOriSegPointsWorld is fixed
	m_nRealOriModelPointsWorld=m_vModelPointsWorld[m_nRealOriModelPointsPosi];//Set the real start points of model.m_nRealOriModelPointsWorld is changing after one model correction.


	// update 'm_nModelPos'
	if (!UpdateModelSgmtPointNew())//set m_nOneSegModelVectorEndPos
	{
		return false;
	}
	
	//	if (!UpdateModelSgmtPointCD())//set m_nOneSegModelVectorEndPos
	//{
	//	return false;
	//}

	if(m_nOneSegModelVectorStartPos==m_vModelPointsWorld.size()-1)//if point C is the end point in model points vector,exchange the start and end.
	{
		int temp=m_nOneSegModelVectorEndPos;
		m_nOneSegModelVectorEndPos=m_nOneSegModelVectorStartPos;
		m_nOneSegModelVectorStartPos=temp;
		IsEndPoint=true;
	}

	m_fCurrModelCArchLength = 0 ; 
	for (int i = 1; i <=m_nOneSegModelVectorStartPos ; i++)
	{ 
		m_fCurrModelCArchLength += sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1])); 
	}

	m_fCurrModelCArchLengthD = 0 ; 
	for (int i = 1; i <=m_nOneSegModelVectorEndPos ; i++)
	{ 
		m_fCurrModelCArchLengthD += sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1])); 
	}
	if( !IsEndPoint== true )
	{

		//***Simu-Model***///
		//UpdateMeanStdOfUnseanImg( m_nOneSegModelVectorStartPos, m_nOneSegModelVectorEndPos);//put the model points into model image to calculate the meanstd of next segmentation add by JDQ 
				//***Simu-Model***///
		//ModifiedUpdateMeanStdOfUnseanImg( ModelVectorStartPosi, ModelVectorEndPosi);//read the mean intensity from the intensity.txt and std of unseen image is the gAortaStdDev[1] add by JDQ
		//ModifiedUpdateMeanStdOfUnseanImgFB( ModelVectorStartPosi, ModelVectorEndPosi);
	}

	return true;
}

bool miiMinPathModel::ModifiedReFastMarchingInitdontMoveModelSM(bool bCrctFlag)
{
	// correct the CA model
	/*	
	if (bCrctFlag)
	{
		CorrectModelPoint();
	}
	*/
	
	// QJD todo
	
	if( bIsCurSMEndPoint== true )//current sgment point is near from model end-point
	{
		cout << "Current point has reached the manual selected end-point!\n" << endl;
		return false ; 
	}
	if (bCrctFlag)
	{
		CorrectModelPointBasedSegLengthdontMoveModelSM();//set m_nOneSegModelVectorStartPos
	}


	if(m_nFMMEvlNum == 1)
	{
		m_nRealOriModelPointsPosi=m_nOneSegModelVectorStartPos;
		m_nRealOriSegModelPointsWorld=m_vModelPointsWorld[m_nOneSegModelVectorStartPos];
		m_nRealOriSegPointWorld=m_iSgmtPointWorld;
	}//m_nRealOriSegPointsWorld is fixed
	m_nRealOriModelPointsWorld=m_vModelPointsWorld[m_nRealOriModelPointsPosi];//Set the real start points of model.m_nRealOriModelPointsWorld is changing after one model correction.


	// update 'm_nModelPos'
	if (!UpdateModelSgmtPointNewSM())//set m_nOneSegModelVectorEndPos
	{
		return false;
	}
	
	//	if (!UpdateModelSgmtPointCD())//set m_nOneSegModelVectorEndPos
	//{
	//	return false;
	//}

	if(m_nOneSegModelVectorStartPos==m_vModelPointsWorld.size()-1)//if point C is the end point in model points vector,exchange the start and end.
	{
		int temp=m_nOneSegModelVectorEndPos;
		m_nOneSegModelVectorEndPos=m_nOneSegModelVectorStartPos;
		m_nOneSegModelVectorStartPos=temp;
	}
	if( !IsEndPoint== true )
	{

		//***Simu-Model***///
		//UpdateMeanStdOfUnseanImg( m_nOneSegModelVectorStartPos, m_nOneSegModelVectorEndPos);//put the model points into model image to calculate the meanstd of next segmentation add by JDQ 
				//***Simu-Model***///
		//ModifiedUpdateMeanStdOfUnseanImg( ModelVectorStartPosi, ModelVectorEndPosi);//read the mean intensity from the intensity.txt and std of unseen image is the gAortaStdDev[1] add by JDQ
		//ModifiedUpdateMeanStdOfUnseanImgFB( ModelVectorStartPosi, ModelVectorEndPosi);
	}

	return true;
}
bool miiMinPathModel::NewReFastMarchingInit(bool bCrctFlag)
{
	// correct the CA model
	/*	
	if (bCrctFlag)
	{
		CorrectModelPoint();
	}
	*/
	if (bCrctFlag)
	{
		CorrectModelPointBasedSegLength();
	}
	int ModelVectorStartPosi=m_nModelPos;
	// update 'm_nModelPos'
	if (!UpdateModelSgmtPoint())
	{
		return false;
	}
	int ModelVectorEndPosi=m_nModelPos;

	UpdateMeanStdOfUnseanImg( ModelVectorStartPosi, ModelVectorEndPosi);
	
	return true;
}
bool miiMinPathModel::UpdateMeanStdOfUnseanImg( int ModelVectorStartPosi, int ModelVectorEndPosi)
{
	double SumIntensity=0,fAtlsCAMean=0,fAtlsCAStdDev=0;
	int nCount=0;
	for(int i=ModelVectorStartPosi+1;i<=ModelVectorEndPosi;i++)
	{
		SumIntensity=SumIntensity+m_vModelPointsWorldWithIntensity[i].val;
		nCount++;
	}
	double fAtlsAortaMean=m_dAtlasStartMean;
	double fUnseenAortaMean=m_dUnseenStartMean;
	fAtlsCAMean=SumIntensity/(nCount+ 0.00000001);
	SumIntensity=0;
	double fMeanDiff = fAtlsAortaMean - fAtlsCAMean;
	m_fMean=fUnseenAortaMean - fMeanDiff;
	/*for(int i=ModelVectorStartPosi+1;i<=ModelVectorEndPosi;i++)
	{
		if (m_vModelPointsWorldWithIntensity[i].val> 10)
				{
					double temp =m_vModelPointsWorldWithIntensity[i].val- fAtlsCAMean;
					temp *= temp;
					SumIntensity += temp;
				} 
	}
	fAtlsCAStdDev= sqrt(SumIntensity/ (nCount + 1));*/
	double fUnseenAortaStdDev=m_dUnseenStartStdDev;
	double fAtlsAortaStdDev=m_dAtlasStartStdDev;
	m_fStdDev=fUnseenAortaStdDev;
  /* m_fStdDev=1024*fAtlsCAStdDev * fUnseenAortaStdDev / (fAtlsAortaStdDev + 0.00000001)/(m_dRawMax- m_dRawMin);	*/
		

	return true;
}
bool miiMinPathModel::ModifiedUpdateMeanStdOfUnseanImg( int ModelVectorStartPosi, int ModelVectorEndPosi)
{
	double fAtlsCAMean=m_vModelPointsWorldWithIntensity[ModelVectorEndPosi].val;
	
	
	double fAtlsAortaMean=m_dAtlasStartMean;
	double fUnseenAortaMean=m_dUnseenStartMean;
	double fMeanDiff = fAtlsAortaMean - fAtlsCAMean;
	double fUnseenCAMean=fUnseenAortaMean-fMeanDiff;
    double fUnseenCAStdDev=m_dUnseenStartStdDev;
	/*double fUnseenCAMean=(fUnseenAortaMean-fMeanDiff);
    double fUnseenCAStdDev=m_dUnseenStartStdDev;*/
     SetMeanStdDev(fUnseenCAMean, fUnseenCAStdDev);
	return true;
}
bool miiMinPathModel::ModifiedUpdateMeanStdOfUnseanImgFB( int ModelVectorStartPosi, int ModelVectorEndPosi)
{
	float fsearchrange=5;
	int isearchrange=(int)(fsearchrange/m_meandis_of_coronarymodel+0.5);
	int i=ModelVectorEndPosi;
		int iF=i-isearchrange;
		int iB=i+isearchrange;
		if (iF<0)iF=0;
		if (iB>m_vModelPointsWorldWithIntensity.size()-1)iB=m_vModelPointsWorldWithIntensity.size()-1;
		double temp=0;
	    int iCount=0;
		for (int j=iF;j<=iB;j++)
		{
			temp=temp+m_vModelPointsWorldWithIntensity[j].val;
		   iCount++;
		}
	double	fAtlsCAMean =temp/iCount;
	double fAtlsAortaMean=m_dAtlasStartMean;
	double fUnseenAortaMean=m_dUnseenStartMean;
	double fMeanDiff = fAtlsAortaMean - fAtlsCAMean;
	double fUnseenCAMean=fUnseenAortaMean-fMeanDiff;
    double fUnseenCAStdDev=m_dUnseenStartStdDev;
	/*double fUnseenCAMean=(fUnseenAortaMean-fMeanDiff);
    double fUnseenCAStdDev=m_dUnseenStartStdDev;*/
     SetMeanStdDev(fUnseenCAMean, fUnseenCAStdDev);
	return true;
}


// Function Name: UpdateUnseenImgIntensity()
//
// Parameters: intensity data
//
// Description: update the whole value of intensity data which is used for calculating potential function
//
// Returns: 
//
bool miiMinPathModel::UpdateUnseenImgIntensity(const short *sImgData)
{
	UnseenImgIntensitySet(sImgData,m_sUnseenImgIntensity);
	return true;
}
// Function Name: FastMarchingEvolution()
//
// Parameters: 
//
// Description: fast-marching-method evolution
//
// Returns: iteration number
//
int miiMinPathModel::FastMarchingEvolution(int nMaxItrNum,int nItrNumSum)
{
	// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
if((nItrNumSum+m_nFmItr)%m_nSaveFmItr==0)
		{
		char chTemp[25];
		_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		char *chFileName4 = (char *)malloc(nLen2);
		strcpy(chFileName4, m_cResultPath);
		strcat(chFileName4, "/MinNodeNBMinPath");
		strcat(chFileName4, chTemp);
		strcat(chFileName4, ".nii.gz");
		zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
		free(chFileName4);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
	
		/*float n=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		float n1=CalcDistmm2(m_vModelPointsWorld[1],m_vModelPointsWorld[2]);*/

		ofstream WriteFileTxt;
	  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/original/potentialnotchangeinupwind/m_dWP.txt",ios::app);//Add by JDQ
	  WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," <<iMinNode.y << "," << iMinNode.z<<"," <<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<iMinNode.val<<")"<<"\n ";//add by JDQ
	 WriteFileTxt.close();
		// search 6-neighbors of the minimum
		
		for (int i = 0; i < 6; i++)
		{
			
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float theta;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					theta = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					theta = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				UpWind1(nx, ny, nz, theta);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
				
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nModelPos].x-m_iStartPointWorld.x;
					 CurMVector[1]=m_vModelPointsWorld[m_nModelPos].y-m_iStartPointWorld.y;
					 CurMVector[2]=m_vModelPointsWorld[m_nModelPos].z-m_iStartPointWorld.z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					
					if (blSegPointFind) 
					{

						return m_nFmItr;
					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);
				}
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		{
						
			return m_nFmItr;
		}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::FastMarchingEvolution(int nMaxItrNum)
{
	// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;

	// iterate for FMM
	while (true)
	{
		m_nFmItr++;

		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		/*if(m_nFmItr%m_nSaveFmItr==0)
		{
		char chTemp[25];
		_itoa_s(m_nFmItr, chTemp, 10);
		int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".vtk") + 1;
		char *chFileName4 = (char *)malloc(nLen2);
		strcpy(chFileName4, m_cResultPath);
		strcat(chFileName4, "/MinNodeNBMinPath");
		strcat(chFileName4, chTemp);
		
		SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		free(chFileName4);
		}*/
	
		/*float n=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		float n1=CalcDistmm2(m_vModelPointsWorld[1],m_vModelPointsWorld[2]);*/
		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float theata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					theata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					theata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWind(nx, ny, nz, theata);
				UpWind1(nx, ny, nz, theata);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
				
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nModelPos].x-m_iStartPointWorld.x;
					 CurMVector[1]=m_vModelPointsWorld[m_nModelPos].y-m_iStartPointWorld.y;
					 CurMVector[2]=m_vModelPointsWorld[m_nModelPos].z-m_iStartPointWorld.z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					
					if (blSegPointFind) 
					{

						return m_nFmItr;
					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);
				}
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		{
		
						
			return m_nFmItr;
		}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::NewFastMarchingEvolution(int nMaxItrNum,int nItrNumSum)
{
	// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;


		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		/*if((nItrNumSum+m_nFmItr)%m_nSaveFmItr==0)
		{
		char chTemp[25];
		_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		char *chFileName4 = (char *)malloc(nLen2);
		strcpy(chFileName4, m_cResultPath);
		strcat(chFileName4, "/MinNodeNBMinPath");
		strcat(chFileName4, chTemp);
		strcat(chFileName4, ".nii.gz");
		SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		free(chFileName4);
		}*/
	
		/*float n=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		float n1=CalcDistmm2(m_vModelPointsWorld[1],m_vModelPointsWorld[2]);*/
		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				//NewUpWindPlus(nx, ny, nz, costheata);
				NewUpWindTimePlusMixed(nx, ny, nz, costheata);
				
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
				
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nModelPos].x-m_iStartPointWorld.x;
					 CurMVector[1]=m_vModelPointsWorld[m_nModelPos].y-m_iStartPointWorld.y;
					 CurMVector[2]=m_vModelPointsWorld[m_nModelPos].z-m_iStartPointWorld.z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					
					if (blSegPointFind) 
					{

						return m_nFmItr;
					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);
				}
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		{
		
						
			return m_nFmItr;
		}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolution(int nMaxItrNum,int nItrNumSum,const short *sImgData)
{
	// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;

		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		m_nSum=nItrNumSum+m_nFmItr;
	/*	if((m_nSum)%m_nSaveFmItr==0&&(nItrNumSum+m_nFmItr)<=200000)
		{
		char chTemp[25];
		_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		char *chFileName4 = (char *)malloc(nLen2);
		strcpy(chFileName4, m_cResultPath);
		strcat(chFileName4, "/MinNodeNBMinPath");
		strcat(chFileName4, chTemp);
		strcat(chFileName4, ".nii.gz");
		SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		free(chFileName4);
		}*/
	
		/*float n=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		float n1=CalcDistmm2(m_vModelPointsWorld[1],m_vModelPointsWorld[2]);*/
	
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/withintensityupdated/m_dWP.txt",ios::app);//Add by JDQ
	 // WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," <<iMinNode.y << "," << iMinNode.z<<"," <<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<iMinNode.val<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				//NewUpWindPlus(nx, ny, nz, costheata);
				//NewUpWindTime(nx, ny, nz, costheata);
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				//ModifiedUpWindPlusTimeMixed(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
				
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nModelPos].x-m_iStartPointWorld.x;
					 CurMVector[1]=m_vModelPointsWorld[m_nModelPos].y-m_iStartPointWorld.y;
					 CurMVector[2]=m_vModelPointsWorld[m_nModelPos].z-m_iStartPointWorld.z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						/*CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						/*CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
						  CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);
						
					}
					
					if (blSegPointFind) 
					{

						return m_nFmItr;
					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);
				}
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		{
		
						
			return m_nFmItr;
		}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}

	int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModel(int nMaxItrNum,int nItrNumSum,const short *sImgData)
{
	// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
      /* if((m_nSum)%m_nSaveFmItr==0&&(nItrNumSum+m_nFmItr)<=200000)
		{
			char chTemp[25];
			_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);
		}*/
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		
		
	
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/withintensityupdated/m_dWP.txt",ios::app);//Add by JDQ
	 // WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," <<iMinNode.y << "," << iMinNode.z<<"," <<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<iMinNode.val<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		//if(m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]!=0)
		//	double m=0;
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	 //WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," << iMinNode.y << "," << iMinNode.z<<","<<m_dSims[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dSimd[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				//NewUpWindPlus(nx, ny, nz, costheata);
				//NewUpWindTime(nx, ny, nz, costheata);
				/*if(!IsEndPoint)
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				else
				ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);*/
				//ModifiedUpWindPlusTimeMixed(nx, ny, nz, costheata,sImgData);
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				//ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				 m_zxhMinNodeImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
			
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
					 CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
					 CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						/*CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);
						
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						/*CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/

						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
						  m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);
						
					}
					
					if (blSegPointFind) 
					{
						char chTemp[25];
						_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
						int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen2);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						return m_nFmItr;
					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				}
				
				//add by jdq
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointC(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
      /* if((m_nSum)%m_nSaveFmItr==0&&(nItrNumSum+m_nFmItr)<=200000)
		{
			char chTemp[25];
			_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);
		}*/
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		m_nSum=nItrNumSum+m_nFmItr;
		
	
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/withintensityupdated/m_dWP.txt",ios::app);//Add by JDQ
	 // WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," <<iMinNode.y << "," << iMinNode.z<<"," <<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<iMinNode.val<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		//if(m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]!=0)
		//	double m=0;
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	 //WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," << iMinNode.y << "," << iMinNode.z<<","<<m_dSims[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dSimd[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				//NewUpWindPlus(nx, ny, nz, costheata);
				//NewUpWindTime(nx, ny, nz, costheata);
				/*if(!IsEndPoint)
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				else
				ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);*/
				//ModifiedUpWindPlusTimeMixed(nx, ny, nz, costheata,sImgData);
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				//ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				  m_zxhMinNodeImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
			
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
					 CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
					 CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						/*CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
						
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						/*CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/

						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
						  m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
						
					}
					
					if (blSegPointFind) 
					{
						char chTemp[25];
						_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
						int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen2);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						return m_nFmItr;
					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				}
				
				//add by jdq
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}

int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithMask(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
       /*if((m_nSum)%m_nSaveFmItr==0&&(nItrNumSum+m_nFmItr)<=200000)
		{
			char chTemp[25];
			_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);

		}*/
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		 m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		m_nSum=nItrNumSum+m_nFmItr;
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/withintensityupdated/m_dWP.txt",ios::app);//Add by JDQ
	 // WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," <<iMinNode.y << "," << iMinNode.z<<"," <<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<iMinNode.val<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		//if(m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]!=0)
		//	double m=0;
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	 //WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," << iMinNode.y << "," << iMinNode.z<<","<<m_dSims[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dSimd[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
			{
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				//NewUpWindPlus(nx, ny, nz, costheata);
				//NewUpWindTime(nx, ny, nz, costheata);
				/*if(!IsEndPoint)
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				else
				ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);*/
				//ModifiedUpWindPlusTimeMixed(nx, ny, nz, costheata,sImgData);
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				//ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				 //sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iAlivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
			
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
					 CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
					 CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						/*CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                         m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
						 //m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						/*CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/

						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
						  m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
						//m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
					}
					
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen2);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						//output badpoints into an image
						char chTemp1[25];
						_itoa_s(m_nBadpointSum, chTemp1, 10);
						int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						char *chFileName5 = (char *)malloc(nLen3);
						strcpy(chFileName5, m_cResultPath);
						strcat(chFileName5, "/Badpoints");
						strcat(chFileName5, chTemp1);
						strcat(chFileName5, ".nii.gz");
						zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						free(chFileName5);
						return m_nFmItr;


					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				}
				
				//add by jdq
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}

int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNew(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		m_nSum=nItrNumSum+m_nFmItr;
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];

			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);

				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
			
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						//char chTemp[25];
						//_itoa_s(nItrNumSum, chTemp, 10);
						//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName4 = (char *)malloc(nLen2);
						//strcpy(chFileName4, m_cResultPath);
						//strcat(chFileName4, "/MinNodeNBMinPath");
						//strcat(chFileName4, chTemp);
						//strcat(chFileName4, ".nii.gz");
						//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						//free(chFileName4);
						////output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for


		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}

	
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=60000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
			
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
				float costheta1,costheta2;	
			
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
				    costheta1 = CalcTheata(nx, ny, nz, true);
					ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
				}

				// upwind
			
			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen4);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						//output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						//output dimD
						int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName8 = (char *)malloc(nLen8);
						strcpy(chFileName8, m_cResultPath);
						strcat(chFileName8, "/NBDValue");
						strcat(chFileName8, chTemp);
						strcat(chFileName8, ".nii.gz");
						zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						free(chFileName8);
						
						//output the model vector as an nii.gz image
						zxhImageData zxhModelpointsImg;
						zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName2 = (char *)malloc(nLen2);
						strcpy(chFileName2, m_cResultPath);
						strcat(chFileName2, "/ModelVector");
						strcat(chFileName2, chTemp);
						strcat(chFileName2, ".nii.gz");
						zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						free(chFileName2);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for

		// reach around the end point of the model 
		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);


		if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		{
			IsCurEndPoint=true;

			return m_nFmItr;
		}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum||m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaF(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=60000000;
	//int iMaxSumSelNum=60000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
			
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
				float costheta1,costheta2;	
			   int DeaFNUM=6;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
				    costheta1 = CalcTheata(nx, ny, nz, true);
					ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
				}
				else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
				{
					//m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
				}
				// upwind
			
			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						/*char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen4);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);*/
						//output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						//output dimD
						//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName8 = (char *)malloc(nLen8);
						//strcpy(chFileName8, m_cResultPath);
						//strcat(chFileName8, "/NBDValue");
						//strcat(chFileName8, chTemp);
						//strcat(chFileName8, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName8);
						
						//output the model vector as an nii.gz image
						//zxhImageData zxhModelpointsImg;
						//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName2 = (char *)malloc(nLen2);
						//strcpy(chFileName2, m_cResultPath);
						//strcat(chFileName2, "/ModelVector");
						//strcat(chFileName2, chTemp);
						//strcat(chFileName2, ".nii.gz");
						//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						//free(chFileName2);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for

		// reach around the end point of the model 
		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);


		if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		{
			IsCurEndPoint=true;

			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum||m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	//int iMaxSumSelNum=60000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ

			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];

				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2;	
					int DeaFNUM=6;
					// the candidate point in aorta
					if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
					{
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz, true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
					{
						//m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					else
					{
						m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);


					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						// check whether segment point 
						bool blSegPointFind = false ;  
						CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
						if (blSegPointFind) 
						{//output alive points and trial points into an image
							/*char chTemp[25];
							_itoa_s(nItrNumSum, chTemp, 10);
							int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
							char *chFileName4 = (char *)malloc(nLen4);
							strcpy(chFileName4, m_cResultPath);
							strcat(chFileName4, "/MinNodeNBMinPath");
							strcat(chFileName4, chTemp);
							strcat(chFileName4, ".nii.gz");
							zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
							free(chFileName4);*/
							//output badpoints into an image
							//char chTemp1[25];
							//_itoa_s(m_nBadpointSum, chTemp1, 10);
							//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
							//char *chFileName5 = (char *)malloc(nLen3);
							//strcpy(chFileName5, m_cResultPath);
							//strcat(chFileName5, "/Badpoints");
							//strcat(chFileName5, chTemp1);
							//strcat(chFileName5, ".nii.gz");
							//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
							//free(chFileName5);
							//output dimD
							//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
							//char *chFileName8 = (char *)malloc(nLen8);
							//strcpy(chFileName8, m_cResultPath);
							//strcat(chFileName8, "/NBDValue");
							//strcat(chFileName8, chTemp);
							//strcat(chFileName8, ".nii.gz");
							//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
							//free(chFileName8);

							//output the model vector as an nii.gz image
							//zxhImageData zxhModelpointsImg;
							//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
							//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
							//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
							//char *chFileName2 = (char *)malloc(nLen2);
							//strcpy(chFileName2, m_cResultPath);
							//strcat(chFileName2, "/ModelVector");
							//strcat(chFileName2, chTemp);
							//strcat(chFileName2, ".nii.gz");
							//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
							//free(chFileName2);
							return m_nFmItr; 
						}
					}
					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
					//add by jdq
				}
			}//for

			// reach around the end point of the model 
			float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);


			if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
			{
				IsCurEndPoint=true;

				return m_nFmItr;
			}

			if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
			{
				return m_nFmItr;
			}

			if (m_nFmItr > nMaxItrNum||m_nSum>iMaxSumSelNum)
			{
				return m_nFmItr;
			}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_D(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		//{
		//	char chTemp[25];
		//	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%50000==0&&m_nFMVedPoiNUM>500000)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{

						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2;	
						int DeaFNUM=10;
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz, true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							//m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						//set the prameter's values as image intensity
						double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,1000*iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);


						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							// check whether segment point 
							bool blSegPointFind = false ;  
							CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2;	
						int DeaFNUM=6;
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz, true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							//m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						//set the prameter's values as image intensity
						double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,1000*iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);


						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							// check whether segment point 
							bool blSegPointFind = false ;  
							CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}

int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==178&&iMinNode.y==142&&iMinNode.z==190)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						
						float currentvec[3]={0,0,0};
						currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							/*costheta1 = CalcTheataN(nx, ny, nz,currentvec);
							costheta2= CalcSegThetaN(nx, ny, nz,currentvec);*/
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float currentvec[3]={0,0,0};
						currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							//m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}

int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==178&&iMinNode.y==142&&iMinNode.z==190)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==178&&iMinNode.y==142&&iMinNode.z==190)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}

int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC_DMP(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==178&&iMinNode.y==142&&iMinNode.z==190)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV_DMP(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV_DMP(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_DMP(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV_DMP(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV_DMP(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_DMP(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}

int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC_BPD(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		if(m_zxhBadPointMaskImg.GetPixelGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0)==ZXH_Foreground) 
			continue;
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==161&&iMinNode.y==121&&iMinNode.z==183)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
			        	double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>1000000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
							//{
							//	  CoutSandEndPosi();
							//}
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
							string strSpeedname = string(m_cResultPath);
							strSpeedname = strSpeedname + "/" + string("Speed.txt");
							ofstream WriteFileSP(strSpeedname);
							for(int i=0;i<m_vSpNorValue.size();i++)
							{
							WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
							} 
							WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);

							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image

								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);

								//if(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast>100)
							//{
									return 1; 
								//}
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}
int miiMinPathModel::DFM_WOBP(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		if(m_zxhBadPointMaskImg.GetPixelGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0)==ZXH_Foreground) 
		{
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			continue;

		}
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==161&&iMinNode.y==121&&iMinNode.z==183)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
						if(nx==253&&ny==257&&nz==268)
						{
							int n=0;
						short x=m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0);
						int mmm=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						}
				if (m_nFMMEvlNum<13)
				{

						if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
								//ModifiedUpWind_WOBP(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
									ModifiedUpWind_WOBP(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
									
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
			        	double p=sqrt(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]);
						double d=100*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPoint_DFM_WOP(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=30;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						
										ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						double p=sqrt(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]);
						double d=100*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>1000000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
							//{
							//	  CoutSandEndPosi();
							//}
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
							string strSpeedname = string(m_cResultPath);
							strSpeedname = strSpeedname + "/" + string("Speed.txt");
							ofstream WriteFileSP(strSpeedname);
							for(int i=0;i<m_vSpNorValue.size();i++)
							{
							WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
							} 
							WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);

							CheckIsSgmtPoint_DFM_WOP(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image

								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);

								//if(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast>100)
							//{
									return 1; 
								//}
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}
int miiMinPathModel::DFM_KT(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>0)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
			CoutSandEndPosi();
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		if(m_zxhBadPointMaskImg.GetPixelGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0)==ZXH_Foreground) 
			continue;
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==161&&iMinNode.y==121&&iMinNode.z==183)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModUpWind_MultiST(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						if(iNewNode.x==180&&iNewNode.y==145&&iNewNode.z==189)
				             int m=0;
						
						//set the prameter's values as image intensity
						double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModUpWind_MultiST(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}
int miiMinPathModel::DFM_KT_NewNeig(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };
		/*int gNbr[26][3] = {
			 {-1, 1, 1}, \
			 {-1, 1, 0}, \
			 {-1, 1, -1}, \
			 {-1, 0, 1}, \
             {-1, 0, 0}, \
			 {-1, 0, -1}, \
			 {-1, -1, 1}, \
			 {-1, -1, 0}, \
			 {-1, -1, -1}, \
					   { 0, 1, 1}, \
					   { 0,1, 0}, \
					   { 0, 1, -1}, \
					   { 0, 0,1}, \
					   { 0, 0, -1},\
					   { 0, -1, 1}, \
					   { 0, -1, 0}, \
					   { 0, -1, -1}, \
					    
			 {1, 1, 1}, \
			 {1, 1, 0}, \
			 {1, 1, -1}, \
			 {1, 0, 1}, \
             {1, 0, 0}, \
			 {1, 0, -1}, \
			 {1, -1, 1}, \
			 {1, -1, 0}, \
			 {1, -1, -1}, \
		};*/

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//cout<<"Visiting "<<m_nFMVedPoiNUM<<" Points"<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>0)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
			//CoutSandEndPosi();
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		if(m_zxhBadPointMaskImg.GetPixelGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0)==ZXH_Foreground) 
			continue;
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%10000==0&&m_nFMVedPoiNUM>0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==161&&iMinNode.y==121&&iMinNode.z==183)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if(nx < 0)nx=0;
				if(nx >= m_nImgWX)nx=m_nImgWX;
				if(ny < 0)ny=0;
				if(ny >= m_nImgWY)ny=m_nImgWY;
				if(nz < 0)nz=0;
				if(nz >= m_nImgWZ)nz=m_nImgWZ;
				if(nx==94&&ny==118&&nz==69)
				int m=0;
				if(nx==94&&ny==116&&nz==68)
				int m2=0;
				int is=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				if (m_nFMMEvlNum<3)
				{

					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 0)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>=1)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModUpWind_MultiST_test(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				
						
						//set the prameter's values as image intensity
						float p=1/(sqrt(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x])+0.000001);
						float d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						float d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						float s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						float s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetPixelByGreyscale(iNewNode.x,iNewNode.y,iNewNode.z,0,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
						//cout<<"Calcu "<<i<<" neighbour Points before if"<<nx<<","<<ny<<","<<nz<<","<<endl;//add by jdq
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=25;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModUpWind_MultiST_test(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
							//cout<<"Calcu u value of "<<nx<<","<<ny<<","<<nz<<","<<endl;//add by jdq
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModUpWind_MultiST_VslsOnly(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						

						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						float p=1/(sqrt(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x])+0.000001);
						float d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						float d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						float s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						float s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
					   m_zxhNBUvalueImg.SetPixelByGreyscale(iNewNode.x,iNewNode.y,iNewNode.z,0,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
						
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
						
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						
					}//if

				}
			}//for

	}//while

	return 1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC_BPD_O(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		if(m_zxhBadPointMaskImg.GetPixelGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0)==ZXH_Foreground) 
			continue;
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==178&&iMinNode.y==142&&iMinNode.z==190)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							//select the right vessel vector
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD_O(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							SelectVV(nx,ny,nz,u1,currentvec);
							costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
							//costheta3=CalcTheataN(nx, ny, nz,currentvec);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta1,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodeu1vector[3]={0,0,0};
							fCurMinNodeu1vector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodeu1vector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodeu1vector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							float fCurMinNodevector[3]={0,0,0};
							SelectVV(iMinNode.x,iMinNode.y,iMinNode.z,fCurMinNodeu1vector,fCurMinNodevector);
							CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD_O(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_3S(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
		//char chTemp[25];
		//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName5 = (char *)malloc(nLen5);
		//strcpy(chFileName5, m_cResultPath);
		//strcat(chFileName5, "/NBUValue");
		//strcat(chFileName5, chTemp);
		//strcat(chFileName5, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName9 = (char *)malloc(nLen9);
		//strcpy(chFileName9, m_cResultPath);
		//strcat(chFileName9, "/NBSPeedValue");
		//strcat(chFileName9, chTemp);
		//strcat(chFileName9, ".nii.gz");
		//zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName9);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%20000000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
 	//	  CoutSandEndPosi();
		//}

			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
			if(iMinNode.x==178&&iMinNode.y==142&&iMinNode.z==190)
				int m=0;
			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&!m_bFistEvlFlg)
						{
							if(SelectVV_3S(nx,ny,nz,u1,currentvec))//the angle between u1 and current vector are seperated into 3 sections
							{
								//if abs(angle)
								costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							    costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
								ModifiedUpWindPlusTime2CosthetaN(nx, ny,nz,costheta1,costheta2,costheta1,sImgData);
							}
							else
							{
								costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
							    costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
								ModifiedUpWindPlusTime2CosthetaN_3S(nx, ny,nz,costheta1,costheta2,costheta1,sImgData);
							}
						}
						// upwind
						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
					{
						
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2,costheta3;	
						int DeaFNUM=6;
						float u1[3]={0,0,0};
						u1[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
						u1[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
						u1[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
						float currentvec[3]={0,0,0};
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz,true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							if(SelectVV_3S(nx,ny,nz,u1,currentvec))//the angle between u1 and current vector are seperated into 3 sections
							{
								//if abs(angle)
								costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
								costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
								ModifiedUpWindPlusTime2CosthetaN(nx, ny,nz,costheta1,costheta2,costheta1,sImgData);
							}
							else
							{
								costheta1 = CalcThetaMV(nx, ny, nz,currentvec);
								costheta2= CalcThetaMLV(nx, ny, nz,currentvec);
								ModifiedUpWindPlusTime2CosthetaN_3S(nx, ny,nz,costheta1,costheta2,costheta1,sImgData);
							}
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						
						
						//set the prameter's values as image intensity
						/*double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);*/
						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							/*ConvS(iNewNode,pImageInfo,nStartPoint);
							if(m_nLn>0&&m_nLn%8000==0)
							{
								string strSpeedname = string(m_cResultPath);
								strSpeedname = strSpeedname + "/" + string("Speed.txt");
								ofstream WriteFileSP(strSpeedname);
								for(int i=0;i<m_vSpNorValue.size();i++)
								{
									WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
								} 
								WriteFileSP.close();
							}*/
							// check whether segment point 
							bool blSegPointFind = false ;  
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
							
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}


int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DNB(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	m_fMaxSPlength=0;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	m_nFMMTimeS= GetTimeNow();
	bool bSBP=false;
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
			char chTemp[25];
			_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);
			//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			//char *chFileName5 = (char *)malloc(nLen5);
			//strcpy(chFileName5, m_cResultPath);
			//strcat(chFileName5, "/NBUValue");
			//strcat(chFileName5, chTemp);
			//strcat(chFileName5, ".nii.gz");
			//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			//free(chFileName5);
			int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName9 = (char *)malloc(nLen9);
			strcpy(chFileName9, m_cResultPath);
			strcat(chFileName9, "/NBSPeedValue");
			strcat(chFileName9, chTemp);
			strcat(chFileName9, ".nii.gz");
			zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		    free(chFileName9);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		// if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%200000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//	CoutSandEndPosi();
		//}

		int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			if (m_nFMMEvlNum<3)
			{
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
				{


					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta
					if (m_nFMMEvlNum == 1)
					{
						// calculate the theata between model and vessel
						//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
						costheta1 = CalcTheata(nx, ny, nz,true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else if (m_nFMMEvlNum>1)
					{

						// calculate the theata between model and vessel
						/*costheta1 = CalcTheataN(nx, ny, nz,currentvec);
						costheta2= CalcSegThetaN(nx, ny, nz,currentvec);*/
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						//output convegence
						//ConvS(iNewNode,pImageInfo,nStartPoint);
						/*if(m_nLn==8000)
						{
						string strSpeedname = string(m_cResultPath);
						strSpeedname = strSpeedname + "/" + string("Speed.txt");
						ofstream WriteFileSP(strSpeedname);
						for(int i=0;i<m_vSpNorValue.size();i++)
						{
						WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
						} 
						WriteFileSP.close();
						}*/
						// check whether segment point 
						/*ConvS(iNewNode,pImageInfo,nStartPoint);
						if(m_nLn==8000)
						{
						string strSpeedname = string(m_cResultPath);
						strSpeedname = strSpeedname + "/" + string("Speed.txt");
						ofstream WriteFileSP(strSpeedname);
						for(int i=0;i<m_vSpNorValue.size();i++)
						{
						WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
						} 
						WriteFileSP.close();
						}*/
						int nFMMTime = GetTimeNow()-m_nFMMTimeS;
						bool blSegPointFind = false ; 
						if(nFMMTime>=300&&!bSBP)//never SBP ,after 300seconds SBP
						{
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else if(bSBP&&nFMMTime>60)//SBP before, 60s SBP once
						{  
							bSBP=false;
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else //
						{
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
								return 1; 
						}	

					}//if



					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 

				}//if

				//add by jdq
			}//if

			else//m_nFMMEvlNum>=3
			{

				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					int DeaFNUM=6;
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta

					if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM)
					{
						//m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					else
					{

						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						//output convegence
						/*ConvS(iNewNode,pImageInfo,nStartPoint);
						if(m_nLn==8000)
						{
						string strSpeedname = string(m_cResultPath);
						strSpeedname = strSpeedname + "/" + string("Speed.txt");
						ofstream WriteFileSP(strSpeedname);
						for(int i=0;i<m_vSpNorValue.size();i++)
						{
						WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
						} 
						WriteFileSP.close();
						}*/
						// check whether segment point 
						int nFMMTime = GetTime()-m_nFMMTimeS;bool blSegPointFind = false ; 
						if(nFMMTime>=300&&!bSBP)
						{
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else if(bSBP&&nFMMTime>60)
						{  
							bSBP=false;
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else 
						{
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
								return 1; 
						}	

					}//if

					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
				}//if going into mask
			}//else m_nFMMEvlNum>=3


		}//for

	}//while

	return -1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DNB_WORec(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	m_fMaxSPlength=0;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	m_nFMMTimeS= GetTimeNow();
	bool bSBP=false;
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
			char chTemp[25];
			_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);
			//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			//char *chFileName5 = (char *)malloc(nLen5);
			//strcpy(chFileName5, m_cResultPath);
			//strcat(chFileName5, "/NBUValue");
			//strcat(chFileName5, chTemp);
			//strcat(chFileName5, ".nii.gz");
			//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			//free(chFileName5);
			int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName9 = (char *)malloc(nLen9);
			strcpy(chFileName9, m_cResultPath);
			strcat(chFileName9, "/NBSPeedValue");
			strcat(chFileName9, chTemp);
			strcat(chFileName9, ".nii.gz");
			zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		    free(chFileName9);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		// if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%200000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//	CoutSandEndPosi();
		//}

		int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			if (m_nFMMEvlNum<3)
			{
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
				{


					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta
					if (m_nFMMEvlNum == 1)
					{
						// calculate the theata between model and vessel
						//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
						costheta1 = CalcTheata(nx, ny, nz,true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else if (m_nFMMEvlNum>1)
					{

						// calculate the theata between model and vessel
						/*costheta1 = CalcTheataN(nx, ny, nz,currentvec);
						costheta2= CalcSegThetaN(nx, ny, nz,currentvec);*/
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						//output convegence
						ConvS(iNewNode,pImageInfo,nStartPoint);
						/*if(m_nLn==8000)
						{
						string strSpeedname = string(m_cResultPath);
						strSpeedname = strSpeedname + "/" + string("Speed.txt");
						ofstream WriteFileSP(strSpeedname);
						for(int i=0;i<m_vSpNorValue.size();i++)
						{
						WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
						} 
						WriteFileSP.close();
						}*/
						// check whether segment point 
						/*ConvS(iNewNode,pImageInfo,nStartPoint);
						if(m_nLn==8000)
						{
						string strSpeedname = string(m_cResultPath);
						strSpeedname = strSpeedname + "/" + string("Speed.txt");
						ofstream WriteFileSP(strSpeedname);
						for(int i=0;i<m_vSpNorValue.size();i++)
						{
						WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
						} 
						WriteFileSP.close();
						}*/
						int nFMMTime = GetTimeNow()-m_nFMMTimeS;
						bool blSegPointFind = false ; 
						if(nFMMTime>=60&&!bSBP)//never SBP ,after 300seconds SBP
						{
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP_WORec(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else //
						{
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
								return 1; 
						}	

					}//if



					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 

				}//if

				//add by jdq
			}//if

			else//m_nFMMEvlNum>=3
			{

				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					int DeaFNUM=3;
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta

					if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM)
					{
						//m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					else
					{

						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						//output convegence
						ConvS(iNewNode,pImageInfo,nStartPoint);
						/*if(m_nLn==8000)
						{
						string strSpeedname = string(m_cResultPath);
						strSpeedname = strSpeedname + "/" + string("Speed.txt");
						ofstream WriteFileSP(strSpeedname);
						for(int i=0;i<m_vSpNorValue.size();i++)
						{
						WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
						} 
						WriteFileSP.close();
						}*/
						// check whether segment point 
						int nFMMTime = GetTime()-m_nFMMTimeS;bool blSegPointFind = false ; 
						if(nFMMTime>=60&&!bSBP)
						{
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP_WORec(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else 
						{
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
								return 1; 
						}	

					}//if

					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
				}//if going into mask
			}//else m_nFMMEvlNum>=3


		}//for

	}//while

	return -1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DNBCV(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	m_fMaxSPlength=0;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	m_nFMMTimeS= GetTimeNow();
	bool bSBP=false;
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		{
			char chTemp[25];
			_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);
			//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			//char *chFileName5 = (char *)malloc(nLen5);
			//strcpy(chFileName5, m_cResultPath);
			//strcat(chFileName5, "/NBUValue");
			//strcat(chFileName5, chTemp);
			//strcat(chFileName5, ".nii.gz");
			//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			//free(chFileName5);
			int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName9 = (char *)malloc(nLen9);
			strcpy(chFileName9, m_cResultPath);
			strcat(chFileName9, "/NBSPeedValue");
			strcat(chFileName9, chTemp);
			strcat(chFileName9, ".nii.gz");
			zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName9);
			int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName5 = (char *)malloc(nLen5);
			strcpy(chFileName5, m_cResultPath);
			strcat(chFileName5, "/NBUValue");
			strcat(chFileName5, chTemp);
			strcat(chFileName5, ".nii.gz");
			zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		// if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		if(m_nFMVedPoiNUM%200000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//	CoutSandEndPosi();
		//}

		int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			if (m_nFMMEvlNum<3)
			{
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
				{


					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta
					if (m_nFMMEvlNum == 1)
					{
						// calculate the theata between model and vessel
						//costheta1 = CalcTheataN(nx, ny, nz,currentvec,true);
						costheta1 = CalcTheata(nx, ny, nz,true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else if (m_nFMMEvlNum>1)
					{

						// calculate the theata between model and vessel
						/*costheta1 = CalcTheataN(nx, ny, nz,currentvec);
						costheta2= CalcSegThetaN(nx, ny, nz,currentvec);*/
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						//output convegence
						//ConvS(iNewNode,pImageInfo,nStartPoint);
						/*if(m_nLn==8000)
						{
						string strSpeedname = string(m_cResultPath);
						strSpeedname = strSpeedname + "/" + string("Speed.txt");
						ofstream WriteFileSP(strSpeedname);
						for(int i=0;i<m_vSpNorValue.size();i++)
						{
						WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
						} 
						WriteFileSP.close();
						}*/
						// check whether segment point 
						ConvS(iNewNode,pImageInfo,nStartPoint);
						bool blSegPointFind = false ; 
						int nFMMTime = GetTimeNow()-m_nFMMTimeS;
						if(m_nLn>20000)
						{
							string strSpeedname = string(m_cResultPath);
							strSpeedname = strSpeedname + "/" + string("Speed.txt");
							ofstream WriteFileSP(strSpeedname);
							for(int i=0;i<m_vSpNorValue.size();i++)
							{
								WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
							} 
							WriteFileSP.close();


							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else if(bSBP&&nFMMTime>60)//SBP before, 60s SBP once
						{  
							bSBP=false;
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
					
						//if(nFMMTime>=300&&!bSBP)//never SBP ,after 300seconds SBP
						//{
						//	int nNumVP=0;
						//	float fP1=0.2;
						//	std::vector<int> input = *new std::vector<int>();
						//	SBPInit(input);
						//	if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
						//	{
						//		delete[] m_nSBMap;
						//		return 1; 
						//	}
						//	else
						//	{
						//		delete[] m_nSBMap;
						//		bSBP=true;
						//	}
						//}
						//else if(bSBP&&nFMMTime>60)//SBP before, 60s SBP once
						//{  
						//	bSBP=false;
						//	int nNumVP=0;
						//	float fP1=0.2;
						//	std::vector<int> input = *new std::vector<int>();
						//	SBPInit(input);
						//	if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
						//	{
						//		delete[] m_nSBMap;
						//		return 1; 
						//	}
						//	else
						//	{
						//		delete[] m_nSBMap;
						//		bSBP=true;
						//	}
						//}
						else //
						{
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
								return 1; 
						}	

					}//if



					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 

				}//if

				//add by jdq
			}//if

			else//m_nFMMEvlNum>=3
			{

				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					int DeaFNUM=6;
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta

					if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM)
					{
						//m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					else
					{

						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						//output convegence
						ConvS(iNewNode,pImageInfo,nStartPoint);
						bool blSegPointFind = false ; 
						int nFMMTime = GetTime()-m_nFMMTimeS;
						if(m_nLn>20000)
						{
							string strSpeedname = string(m_cResultPath);
							strSpeedname = strSpeedname + "/" + string("Speed.txt");
							ofstream WriteFileSP(strSpeedname);
							for(int i=0;i<m_vSpNorValue.size();i++)
							{
								WriteFileSP<< i<< " " << m_vSpNorValue[i].val <<endl;
							} 
							WriteFileSP.close();
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						else if(bSBP&&nFMMTime>60)
						{  
							bSBP=false;
							int nNumVP=0;
							float fP1=0.2;
							std::vector<int> input = *new std::vector<int>();
							SBPInit(input);
							if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
							{
								delete[] m_nSBMap;
								return 1; 
							}
							else
							{
								delete[] m_nSBMap;
								bSBP=true;
							}
						}
						
						// check whether segment point 
						//int nFMMTime = GetTime()-m_nFMMTimeS;
						//if(nFMMTime>=300&&!bSBP)
						//{
						//	int nNumVP=0;
						//	float fP1=0.2;
						//	std::vector<int> input = *new std::vector<int>();
						//	SBPInit(input);
						//	if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
						//	{
						//		delete[] m_nSBMap;
						//		return 1; 
						//	}
						//	else
						//	{
						//		delete[] m_nSBMap;
						//		bSBP=true;
						//	}
						//}
						//else if(bSBP&&nFMMTime>60)
						//{  
						//	bSBP=false;
						//	int nNumVP=0;
						//	float fP1=0.2;
						//	std::vector<int> input = *new std::vector<int>();
						//	SBPInit(input);
						//	if(SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000))//if(SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000))
						//	{
						//		delete[] m_nSBMap;
						//		return 1; 
						//	}
						//	else
						//	{
						//		delete[] m_nSBMap;
						//		bSBP=true;
						//	}
						//}
						else 
						{
							//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							float fCurMinNodevector[3]={0,0,0};
							fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
							fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
							fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
							CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
							if (blSegPointFind) 
								return 1; 
						}	

					}//if

					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
				}//if going into mask
			}//else m_nFMMEvlNum>=3


		}//for

	}//while

	return -1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DFMP(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM,const zxhImageDataT<short>&imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	m_fMaxSPlength=0;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	int nFMMTimeS = GetTime();
	bool bSBP=false;
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		//{
		//	char chTemp[25];
		//	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
			//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			//char *chFileName5 = (char *)malloc(nLen5);
			//strcpy(chFileName5, m_cResultPath);
			//strcat(chFileName5, "/NBUValue");
			//strcat(chFileName5, chTemp);
			//strcat(chFileName5, ".nii.gz");
			//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			//free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		// if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}
		int nFMMTime = GetTime()-nFMMTimeS;bool blSegPointFind = false ; 
		if(m_nFMVedPoiNUM%200000==0)
		{
			cout<<"Visited Points Number:"<<m_nFMVedPoiNUM<<";"<<"The maxiumum Number:"<<MasknoneZeroSumNUM<<endl;
		}
		

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
		//if(m_nFMVedPoiNUM%m_nSaveFmItr==0&&m_nFMVedPoiNUM>100000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//	CoutSandEndPosi();
		//}

		int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetPixelByGreyscale(iMinNode.x,iMinNode.y,iMinNode.z,0,m_iAlivePValue);
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			if (m_nFMMEvlNum<3)
			{
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
				{
					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta
					if (m_nFMMEvlNum == 1)
					{
						costheta1 = CalcTheata(nx, ny, nz,true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else if (m_nFMMEvlNum>1)
					{
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					// upwind
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						// check whether segment point 
						bool blSegPointFind = false ; 
						//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
						float fCurMinNodevector[3]={0,0,0};
						fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
						fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
						fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
						CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
						if (blSegPointFind) 
							return 1; 
					}//if

					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
				}//if
			}//if

			else//m_nFMMEvlNum>=3
			{
			
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetPixelGreyscale(nx,ny,nz,0)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&imgReadLineMask.GetPixelGreyscale(nx,ny,nz,0)!=0)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float costheta1,costheta2,costheta3;	
					int DeaFNUM=6;
					float currentvec[3]={0,0,0};
					currentvec[0]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].x;
					currentvec[1]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].y;
					currentvec[2]=vUVectorPoint[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx].z;
					// the candidate point in aorta
					if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM)
					{
						//m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						costheta3=CalcTheataN(nx, ny, nz,currentvec);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2CosthetaN(nx, ny, nz, costheta1,costheta2,costheta3,sImgData);
					}
					else
					{
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					// upwind
					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					//set the prameter's values as image intensity
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching not including revisted points.
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
						// check whether segment point 
						bool blSegPointFind = false ; 
						//CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,2000,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
						float fCurMinNodevector[3]={0,0,0};
						fCurMinNodevector[0]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].x;
						fCurMinNodevector[1]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].y;
						fCurMinNodevector[2]=vUVectorPoint[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x].z;
						CheckIsSgmtPointAndBadPointMLLN(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,fCurMinNodevector,blSegPointFind) ;
						if (blSegPointFind) 
							return 1; 
					}	//if 

					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
				}//if
			}//else


		}//for

	}//while

	return -1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_D_LR(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int MasknoneZeroSumNUM)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=6000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMVedPoiNUM>80000)
		//{
		//	char chTemp[25];
		//	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
        // if m_vNarrowBand is empty error
		if(m_vNarrowBand.size()==1&&m_nFmItr>1)
		{
			cout<<"The NarrawBand is empty!";
			return -1;
		}
		if(m_nFMVedPoiNUM>6*MasknoneZeroSumNUM)
		{
			cout<<"all points in mask have been visited";
			return -1;
		}

		//int w=50;
		//if(m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
	 //    int x=0;
		//}
		//if((m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)%5000==0&&(m_nFMVedPoiNUM-m_nFMVedPoiNUMLast)>20000)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiL();
		//}
		//int w=100;

		//if(m_vNarrowBand.size()%w==0&&m_nFMMEvlNum>2&&m_vNarrowBand.size()<500)//((m_nFMVedPoiNUM)%m_nSaveFmItr==0&&m_nFMMEvlNum>10&&m_vNarrowBand.size()<500)
		//{
		//CoutSandEndPosiE();
		//}
			int m=1;
			miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
			// update FMM map by setting the minimum's location in FMM map as 'alive' 
			m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
			m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ

			// search 6-neighbors of the minimum 
			for (int i = 0; i < 6; i++)
			{
				int nx = iMinNode.x + gNbr[i][0];
				int ny = iMinNode.y + gNbr[i][1];
				int nz = iMinNode.z + gNbr[i][2];
				if (m_nFMMEvlNum<3)
				{
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{

						m_nFMVedPoiNUM++;//The number of points visited by fastmarching including revisted points.
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2;	
						int DeaFNUM=6;
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz, true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							//m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						//set the prameter's values as image intensity
						double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);


						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							
							m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							// check whether segment point 
							bool blSegPointFind = false ;  
							CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 


						//add by jdq
					}//if
				}//if
				else//m_nFMMEvlNum>=3
				{
			/*		short b=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					short m=m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					short lm=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					if(b!=ZXH_Foreground)
					{
						int p=0;
					}
					if(m!= FMM_ALIVE)
					{
						int p=0;
					}
					if(lm!=0)
					{
						int p=0;
					}*/
					if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
						&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
					{
						m_nFMVedPoiNUM++;//The number of points visited by fastmarching
						int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
						// for fist direction in aorta
						float costheta1,costheta2;	
						int DeaFNUM=6;
						// the candidate point in aorta
						if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
						{
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz, true);
							ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						}
						else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
						{
							//m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						else
						{
							m_bFistEvlFlg = true;
							// calculate the theata between model and vessel
							costheta1 = CalcTheata(nx, ny, nz);
							costheta2= CalcSegTheta(nx, ny, nz);
							//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
							ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
						}
						// upwind


						//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
						iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
						iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
						//set the prameter's values as image intensity
						double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
						double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
						double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
						if(p>33000)
							p=30000;
						if(s>33000)
							s=30000;
						if(d>33000)
							d=30000;
						m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
						m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
						m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
						m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);


						//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
						if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
						{
							m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
							m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
							m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
							// check whether segment point 
							bool blSegPointFind = false ;  
							CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,MasknoneZeroSumNUM,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
							if (blSegPointFind) 
							{//output alive points and trial points into an image
								/*char chTemp[25];
								_itoa_s(nItrNumSum, chTemp, 10);
								int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
								char *chFileName4 = (char *)malloc(nLen4);
								strcpy(chFileName4, m_cResultPath);
								strcat(chFileName4, "/MinNodeNBMinPath");
								strcat(chFileName4, chTemp);
								strcat(chFileName4, ".nii.gz");
								zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
								free(chFileName4);*/
								//output badpoints into an image
								//char chTemp1[25];
								//_itoa_s(m_nBadpointSum, chTemp1, 10);
								//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
								//char *chFileName5 = (char *)malloc(nLen3);
								//strcpy(chFileName5, m_cResultPath);
								//strcat(chFileName5, "/Badpoints");
								//strcat(chFileName5, chTemp1);
								//strcat(chFileName5, ".nii.gz");
								//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
								//free(chFileName5);
								//output dimD
								//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName8 = (char *)malloc(nLen8);
								//strcpy(chFileName8, m_cResultPath);
								//strcat(chFileName8, "/NBDValue");
								//strcat(chFileName8, chTemp);
								//strcat(chFileName8, ".nii.gz");
								//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
								//free(chFileName8);

								//output the model vector as an nii.gz image
								
								//zxhImageData zxhModelpointsImg;
								//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
								//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
								//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
								//char *chFileName2 = (char *)malloc(nLen2);
								//strcpy(chFileName2, m_cResultPath);
								//strcat(chFileName2, "/ModelVector");
								//strcat(chFileName2, chTemp);
								//strcat(chFileName2, ".nii.gz");
								//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
								//free(chFileName2);
								return 1; 
							}
						}
						else //'trial'
						{
							UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
						} 
						//add by jdq
					}//if

				}
			}//for

	}//while

	return 1;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLLDirH(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=500000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
			
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
				float costheta1,costheta2;	
			
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
				    costheta1 = CalcTheata(nx, ny, nz, true);
					ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
				}

				// upwind
			
			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen4);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						//output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						//output dimD
						int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName8 = (char *)malloc(nLen8);
						strcpy(chFileName8, m_cResultPath);
						strcat(chFileName8, "/NBDValue");
						strcat(chFileName8, chTemp);
						strcat(chFileName8, ".nii.gz");
						zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						free(chFileName8);
						
						//output the model vector as an nii.gz image
						zxhImageData zxhModelpointsImg;
						zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName2 = (char *)malloc(nLen2);
						strcpy(chFileName2, m_cResultPath);
						strcat(chFileName2, "/ModelVector");
						strcat(chFileName2, chTemp);
						strcat(chFileName2, ".nii.gz");
						zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						free(chFileName2);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for


		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum&&m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}

int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLLSM(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=1000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nSum)%m_nSaveFmItr==0&&m_nSum>10000)
		{
			char chTemp[25];
			_itoa_s(m_nSum, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);
			int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName5 = (char *)malloc(nLen5);
			strcpy(chFileName5, m_cResultPath);
			strcat(chFileName5, "/NBUValue");
			strcat(chFileName5, chTemp);
			strcat(chFileName5, ".nii.gz");
			zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName5);
			int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName6 = (char *)malloc(nLen6);
			strcpy(chFileName6, m_cResultPath);
			strcat(chFileName6, "/NBPValue");
			strcat(chFileName6, chTemp);
			strcat(chFileName6, ".nii.gz");
			zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName6);
			int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName7 = (char *)malloc(nLen7);
			strcpy(chFileName7, m_cResultPath);
			strcat(chFileName7, "/NBSValue");
			strcat(chFileName7, chTemp);
			strcat(chFileName7, ".nii.gz");
			zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName7);
			int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName8 = (char *)malloc(nLen8);
			strcpy(chFileName8, m_cResultPath);
			strcat(chFileName8, "/NBDValue");
			strcat(chFileName8, chTemp);
			strcat(chFileName8, ".nii.gz");
			zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName8);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
			
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
				float costheta1,costheta2;	
			
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
				    costheta1 = CalcTheata(nx, ny, nz, true);
					ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
				}
				else if (m_nFMMEvlNum < 10 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2CosthetaSM(nx, ny, nz, costheta1,costheta2,sImgData);
				}
				// upwind
			
			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,10*iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLLSM(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						//char chTemp[25];
						//_itoa_s(m_nSum, chTemp, 10);
						//int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName4 = (char *)malloc(nLen4);
						//strcpy(chFileName4, m_cResultPath);
						//strcat(chFileName4, "/MinNodeNBMinPath");
						//strcat(chFileName4, chTemp);
						//strcat(chFileName4, ".nii.gz");
						//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						//free(chFileName4);
						////output badpoints into an image
						////char chTemp1[25];
						////_itoa_s(m_nBadpointSum, chTemp1, 10);
						////int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						////char *chFileName5 = (char *)malloc(nLen3);
						////strcpy(chFileName5, m_cResultPath);
						////strcat(chFileName5, "/Badpoints");
						////strcat(chFileName5, chTemp1);
						////strcat(chFileName5, ".nii.gz");
						////zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						////free(chFileName5);
						////output dimD
						//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName8 = (char *)malloc(nLen8);
						//strcpy(chFileName8, m_cResultPath);
						//strcat(chFileName8, "/NBDValue");
						//strcat(chFileName8, chTemp);
						//strcat(chFileName8, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName8);
						//
						////output the model vector as an nii.gz image
						//zxhImageData zxhModelpointsImg;
						//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName2 = (char *)malloc(nLen2);
						//strcpy(chFileName2, m_cResultPath);
						//strcat(chFileName2, "/ModelVector");
						//strcat(chFileName2, chTemp);
						//strcat(chFileName2, ".nii.gz");
						//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						//free(chFileName2);
						//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen5);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/NBUValue");
						//strcat(chFileName5, chTemp);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName5);
						//int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName6 = (char *)malloc(nLen6);
						//strcpy(chFileName6, m_cResultPath);
						//strcat(chFileName6, "/NBPValue");
						//strcat(chFileName6, chTemp);
						//strcat(chFileName6, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName6);
						//int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName7 = (char *)malloc(nLen7);
						//strcpy(chFileName7, m_cResultPath);
						//strcat(chFileName7, "/NBSValue");
						//strcat(chFileName7, chTemp);
						//strcat(chFileName7, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName7);


						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for


		float nEndDistmm2 = CalcDistmm2(m_fSMEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//	bIsCurSMEndPoint=true;
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum&&m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMask(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		m_nSum=nItrNumSum+m_nFmItr;
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];

			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
			{
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);

				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				//NewUpWindPlus(nx, ny, nz, costheata);
				//NewUpWindTime(nx, ny, nz, costheata);
				/*if(!IsEndPoint)
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				else
				ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);*/
				//ModifiedUpWindPlusTimeMixed(nx, ny, nz, costheata,sImgData);
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				//ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d=1/(0.5*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						//char chTemp[25];
						//_itoa_s(nItrNumSum, chTemp, 10);
						//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName4 = (char *)malloc(nLen2);
						//strcpy(chFileName4, m_cResultPath);
						//strcat(chFileName4, "/MinNodeNBMinPath");
						//strcat(chFileName4, chTemp);
						//strcat(chFileName4, ".nii.gz");
						//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						//free(chFileName4);
						////output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for


		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}
		
		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=60000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
			{
			
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
				float costheta1,costheta2;	
			
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
				    costheta1 = CalcTheata(nx, ny, nz, true);
					ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
				}

				// upwind
			
			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=10*m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=10*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						/*char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen4);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);*/
						//output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						//output dimD
						//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName8 = (char *)malloc(nLen8);
						//strcpy(chFileName8, m_cResultPath);
						//strcat(chFileName8, "/NBDValue");
						//strcat(chFileName8, chTemp);
						//strcat(chFileName8, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName8);
						
						//output the model vector as an nii.gz image
						//zxhImageData zxhModelpointsImg;
						//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName2 = (char *)malloc(nLen2);
						//strcpy(chFileName2, m_cResultPath);
						//strcat(chFileName2, "/ModelVector");
						//strcat(chFileName2, chTemp);
						//strcat(chFileName2, ".nii.gz");
						//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						//free(chFileName2);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for


		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum&&m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
	}
	int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_D(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=2000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			if(m_nFMMEvlNum<3)
			{
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
					float costheta1,costheta2;	

					// the candidate point in aorta
					if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
					{
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz, true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else
					{
						m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
					}

					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

						/*if (m_nFMMEvlNum == 1)
						{
						m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
						}
						else
						{
						m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
						}

						if (blSegPointFind) 
						{
						return m_nFmItr;
						}*/ 
						// check whether segment point 
						bool blSegPointFind = false ;  
						CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
						if (blSegPointFind) 
						{//output alive points and trial points into an image
							/*char chTemp[25];
							_itoa_s(nItrNumSum, chTemp, 10);
							int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
							char *chFileName4 = (char *)malloc(nLen4);
							strcpy(chFileName4, m_cResultPath);
							strcat(chFileName4, "/MinNodeNBMinPath");
							strcat(chFileName4, chTemp);
							strcat(chFileName4, ".nii.gz");
							zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
							free(chFileName4);*/
							//output badpoints into an image
							//char chTemp1[25];
							//_itoa_s(m_nBadpointSum, chTemp1, 10);
							//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
							//char *chFileName5 = (char *)malloc(nLen3);
							//strcpy(chFileName5, m_cResultPath);
							//strcat(chFileName5, "/Badpoints");
							//strcat(chFileName5, chTemp1);
							//strcat(chFileName5, ".nii.gz");
							//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
							//free(chFileName5);
							//output dimD
							//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
							//char *chFileName8 = (char *)malloc(nLen8);
							//strcpy(chFileName8, m_cResultPath);
							//strcat(chFileName8, "/NBDValue");
							//strcat(chFileName8, chTemp);
							//strcat(chFileName8, ".nii.gz");
							//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
							//free(chFileName8);

							//output the model vector as an nii.gz image
							//zxhImageData zxhModelpointsImg;
							//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
							//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
							//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
							//char *chFileName2 = (char *)malloc(nLen2);
							//strcpy(chFileName2, m_cResultPath);
							//strcat(chFileName2, "/ModelVector");
							//strcat(chFileName2, chTemp);
							//strcat(chFileName2, ".nii.gz");
							//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
							//free(chFileName2);
							return m_nFmItr; 
						}
					}
					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
					//add by jdq
				}//if()
			}//if m_nFMMEvlNum
			else if(m_nFMMEvlNum>2)
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
					float costheta1,costheta2;	

					// the candidate point in aorta
					if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
					{
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz, true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else
					{
						m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
					}

					// upwind

			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen4);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						//output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						//output dimD
						//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName8 = (char *)malloc(nLen8);
						//strcpy(chFileName8, m_cResultPath);
						//strcat(chFileName8, "/NBDValue");
						//strcat(chFileName8, chTemp);
						//strcat(chFileName8, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName8);
						
						//output the model vector as an nii.gz image
						//zxhImageData zxhModelpointsImg;
						//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName2 = (char *)malloc(nLen2);
						//strcpy(chFileName2, m_cResultPath);
						//strcat(chFileName2, "/ModelVector");
						//strcat(chFileName2, chTemp);
						//strcat(chFileName2, ".nii.gz");
						//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						//free(chFileName2);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}//else m_nFMMEvlNum
			
		}//for


		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);


		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum&&m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
	}
	int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_DeaF(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=60000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
			{
			
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
				float costheta1,costheta2;	
			      int DeaFNUM=5;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz, true);
					ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
				}
				else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
				{
					//m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
				}

				// upwind
			
			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen4);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						//output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						//output dimD
						//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName8 = (char *)malloc(nLen8);
						//strcpy(chFileName8, m_cResultPath);
						//strcat(chFileName8, "/NBDValue");
						//strcat(chFileName8, chTemp);
						//strcat(chFileName8, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName8);
						
						//output the model vector as an nii.gz image
						//zxhImageData zxhModelpointsImg;
						//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName2 = (char *)malloc(nLen2);
						//strcpy(chFileName2, m_cResultPath);
						//strcat(chFileName2, "/ModelVector");
						//strcat(chFileName2, chTemp);
						//strcat(chFileName2, ".nii.gz");
						//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						//free(chFileName2);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for


		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum&&m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
	}
	int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_DeaF_D(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=2000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		//if((m_nSum)%m_nSaveFmItr==0&&m_nSum>50000)
		//{
		//	char chTemp[25];
		//	_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//	int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName4 = (char *)malloc(nLen2);
		//	strcpy(chFileName4, m_cResultPath);
		//	strcat(chFileName4, "/MinNodeNBMinPath");
		//	strcat(chFileName4, chTemp);
		//	strcat(chFileName4, ".nii.gz");
		//	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName4);
		//	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName5 = (char *)malloc(nLen5);
		//	strcpy(chFileName5, m_cResultPath);
		//	strcat(chFileName5, "/NBUValue");
		//	strcat(chFileName5, chTemp);
		//	strcat(chFileName5, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName5);
		//	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName6 = (char *)malloc(nLen6);
		//	strcpy(chFileName6, m_cResultPath);
		//	strcat(chFileName6, "/NBPValue");
		//	strcat(chFileName6, chTemp);
		//	strcat(chFileName6, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName6);
		//	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName7 = (char *)malloc(nLen7);
		//	strcpy(chFileName7, m_cResultPath);
		//	strcat(chFileName7, "/NBSValue");
		//	strcat(chFileName7, chTemp);
		//	strcat(chFileName7, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName7);
		//	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//	char *chFileName8 = (char *)malloc(nLen8);
		//	strcpy(chFileName8, m_cResultPath);
		//	strcat(chFileName8, "/NBDValue");
		//	strcat(chFileName8, chTemp);
		//	strcat(chFileName8, ".nii.gz");
		//	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//	free(chFileName8);
		//}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			if(m_nFMMEvlNum<3)
			{
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
					float costheta1,costheta2;	

					int DeaFNUM=5;
					// the candidate point in aorta
					if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
					{
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz, true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
					{
						//m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					else
					{
						m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
					}

					// upwind


					//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
					iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
					iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
					double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
					double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
					if(p>33000)
						p=30000;
					if(s>33000)
						s=30000;
					if(d>33000)
						d=30000;
					m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
					m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
					m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
					m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
					//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
					if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
					{
						m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
						m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

						/*if (m_nFMMEvlNum == 1)
						{
						m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
						}
						else
						{
						m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
						}

						if (blSegPointFind) 
						{
						return m_nFmItr;
						}*/ 
						// check whether segment point 
						bool blSegPointFind = false ;  
						CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
						if (blSegPointFind) 
						{//output alive points and trial points into an image
							char chTemp[25];
							_itoa_s(nItrNumSum, chTemp, 10);
							int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
							char *chFileName4 = (char *)malloc(nLen4);
							strcpy(chFileName4, m_cResultPath);
							strcat(chFileName4, "/MinNodeNBMinPath");
							strcat(chFileName4, chTemp);
							strcat(chFileName4, ".nii.gz");
							zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
							free(chFileName4);
							//output badpoints into an image
							//char chTemp1[25];
							//_itoa_s(m_nBadpointSum, chTemp1, 10);
							//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
							//char *chFileName5 = (char *)malloc(nLen3);
							//strcpy(chFileName5, m_cResultPath);
							//strcat(chFileName5, "/Badpoints");
							//strcat(chFileName5, chTemp1);
							//strcat(chFileName5, ".nii.gz");
							//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
							//free(chFileName5);
							//output dimD
							//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
							//char *chFileName8 = (char *)malloc(nLen8);
							//strcpy(chFileName8, m_cResultPath);
							//strcat(chFileName8, "/NBDValue");
							//strcat(chFileName8, chTemp);
							//strcat(chFileName8, ".nii.gz");
							//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
							//free(chFileName8);

							//output the model vector as an nii.gz image
							//zxhImageData zxhModelpointsImg;
							//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
							//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
							//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
							//char *chFileName2 = (char *)malloc(nLen2);
							//strcpy(chFileName2, m_cResultPath);
							//strcat(chFileName2, "/ModelVector");
							//strcat(chFileName2, chTemp);
							//strcat(chFileName2, ".nii.gz");
							//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
							//free(chFileName2);
							return m_nFmItr; 
						}
					}
					else //'trial'
					{
						UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
					} 
					//add by jdq
				}//if()
			}//if m_nFMMEvlNum
			else if(m_nFMMEvlNum>2)
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
					&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
				{

					int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
					// for fist direction in aorta
					float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
					float costheta1,costheta2;	
					int DeaFNUM=5;
					// the candidate point in aorta
					if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
					{
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz, true);
						ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					}
					else if (m_nFMMEvlNum>1&&m_nFMMEvlNum < DeaFNUM &&!m_bFistEvlFlg)
					{
						//m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					else
					{
						m_bFistEvlFlg = true;
						// calculate the theata between model and vessel
						costheta1 = CalcTheata(nx, ny, nz);
						costheta2= CalcSegTheta(nx, ny, nz);
						//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
						ModifiedUpWindPlusTime2Costheta_DeaF(nx, ny, nz, costheta1,costheta2,sImgData);
					}
					// upwind
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLL(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						char chTemp[25];
						_itoa_s(nItrNumSum, chTemp, 10);
						int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						char *chFileName4 = (char *)malloc(nLen4);
						strcpy(chFileName4, m_cResultPath);
						strcat(chFileName4, "/MinNodeNBMinPath");
						strcat(chFileName4, chTemp);
						strcat(chFileName4, ".nii.gz");
						zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						free(chFileName4);
						//output badpoints into an image
						//char chTemp1[25];
						//_itoa_s(m_nBadpointSum, chTemp1, 10);
						//int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen3);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/Badpoints");
						//strcat(chFileName5, chTemp1);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						//free(chFileName5);
						//output dimD
						//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName8 = (char *)malloc(nLen8);
						//strcpy(chFileName8, m_cResultPath);
						//strcat(chFileName8, "/NBDValue");
						//strcat(chFileName8, chTemp);
						//strcat(chFileName8, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName8);
						
						//output the model vector as an nii.gz image
						//zxhImageData zxhModelpointsImg;
						//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName2 = (char *)malloc(nLen2);
						//strcpy(chFileName2, m_cResultPath);
						//strcat(chFileName2, "/ModelVector");
						//strcat(chFileName2, chTemp);
						//strcat(chFileName2, ".nii.gz");
						//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						//free(chFileName2);
						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}//else m_nFMMEvlNum
			
		}//for


		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum&&m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
	}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLLSM(int nMaxItrNum,int nItrNumSum,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData)
{miiCNode<double> iMinNode, iNewNode;

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	int iMaxSumSelNum=60000000;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;
		m_nSum++;
		//output narraw band as image
		if((m_nSum)%m_nSaveFmItr==0&&m_nSum>10000)
		{
			char chTemp[25];
			_itoa_s(m_nSum, chTemp, 10);
			int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName4 = (char *)malloc(nLen2);
			strcpy(chFileName4, m_cResultPath);
			strcat(chFileName4, "/MinNodeNBMinPath");
			strcat(chFileName4, chTemp);
			strcat(chFileName4, ".nii.gz");
			zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName4);
			int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName5 = (char *)malloc(nLen5);
			strcpy(chFileName5, m_cResultPath);
			strcat(chFileName5, "/NBUValue");
			strcat(chFileName5, chTemp);
			strcat(chFileName5, ".nii.gz");
			zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName5);
			int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName6 = (char *)malloc(nLen6);
			strcpy(chFileName6, m_cResultPath);
			strcat(chFileName6, "/NBPValue");
			strcat(chFileName6, chTemp);
			strcat(chFileName6, ".nii.gz");
			zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName6);
			int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName7 = (char *)malloc(nLen7);
			strcpy(chFileName7, m_cResultPath);
			strcat(chFileName7, "/NBSValue");
			strcat(chFileName7, chTemp);
			strcat(chFileName7, ".nii.gz");
			zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName7);
			int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
			char *chFileName8 = (char *)malloc(nLen8);
			strcpy(chFileName8, m_cResultPath);
			strcat(chFileName8, "/NBDValue");
			strcat(chFileName8, chTemp);
			strcat(chFileName8, ".nii.gz");
			zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
			free(chFileName8);
		}
		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		//if(m_nFMMEvlNum>10&&m_vNarrowBand.size()<100)
		//{
		//char chTemp[25];
		//_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		//int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		//char *chFileName4 = (char *)malloc(nLen2);
		//strcpy(chFileName4, m_cResultPath);
		//strcat(chFileName4, "/MinNodeNBMinPath");
		//strcat(chFileName4, chTemp);
		//strcat(chFileName4, ".nii.gz");
		//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		//free(chFileName4);
		//}
			int m=1;
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		m_zxhMinNodeImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iAlivePValue);//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		
		// search 6-neighbors of the minimum 
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&&m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx)!=ZXH_Foreground&&m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
			{
			
				int m=m_zxhBadPointMaskImg.GetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx);
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
				float costheta1,costheta2;	
			
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
				    costheta1 = CalcTheata(nx, ny, nz, true);
					ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
				}
				else if (m_nFMMEvlNum < 10 &&!m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2Costheta(nx, ny, nz, costheta1,costheta2,sImgData);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheta1 = CalcTheata(nx, ny, nz);
					costheta2= CalcSegTheta(nx, ny, nz);
					//ModifiedUpWindPlusTime(nx, ny, nz, costheta1,sImgData);
					ModifiedUpWindPlusTime2CosthetaSM(nx, ny, nz, costheta1,costheta2,sImgData);
				}
				// upwind
			
			
				//ModifiedUpWindTime(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
                double p=1/(m_dP2[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.000001);
				double d=10*m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double d1=m_dSimd[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s1=m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x];
				double s=1/(0.5*m_dSims[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]+0.0001);
				if(p>33000)
					p=30000;
				if(s>33000)
					s=30000;
				if(d>33000)
					d=30000;
				m_zxhNBUvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,10*iNewNode.val);
				m_zxhNBPvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,p);
				m_zxhNBSvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,s);
				m_zxhNBDvalueImg.SetImageData(iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x,d);
				//sAMinNodeImg[iNewNode.z * m_nImgWY * m_nImgWX + iNewNode.y * m_nImgWX + iNewNode.x]=m_iActivePValue;//iNewNode.val;
				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					m_zxhMinNodeImg.SetImageData(nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx,m_iActivePValue);//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ

					/*if (m_nFMMEvlNum == 1)
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
					m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}

					if (blSegPointFind) 
					{
					return m_nFmItr;
					}*/ 
					// check whether segment point 
					bool blSegPointFind = false ;  
					CheckIsSgmtPointAndBadPointMLLSM(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath,blSegPointFind) ;// todo CheckIsSgmtPointAndBadPoint(iMinNode,pImageInfo,nStartPoint,sImgData,nMaxPath) ; 
					if (blSegPointFind) 
					{//output alive points and trial points into an image
						//char chTemp[25];
						//_itoa_s(m_nSum, chTemp, 10);
						//int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName4 = (char *)malloc(nLen4);
						//strcpy(chFileName4, m_cResultPath);
						//strcat(chFileName4, "/MinNodeNBMinPath");
						//strcat(chFileName4, chTemp);
						//strcat(chFileName4, ".nii.gz");
						//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4 );
						//free(chFileName4);
						////output badpoints into an image
						////char chTemp1[25];
						////_itoa_s(m_nBadpointSum, chTemp1, 10);
						////int nLen3 = strlen(m_cResultPath) + strlen("/Badpoints") + strlen(chTemp1) + strlen(".nii.gz") + 1;
						////char *chFileName5 = (char *)malloc(nLen3);
						////strcpy(chFileName5, m_cResultPath);
						////strcat(chFileName5, "/Badpoints");
						////strcat(chFileName5, chTemp1);
						////strcat(chFileName5, ".nii.gz");
						////zxh::SaveImage( &m_zxhBadPointImg, chFileName5 ) ;//SaveImage(m_pBaseImgInfo, sABadPointImg, chFileName5);
						////free(chFileName5);
						////output dimD
						//int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName8 = (char *)malloc(nLen8);
						//strcpy(chFileName8, m_cResultPath);
						//strcat(chFileName8, "/NBDValue");
						//strcat(chFileName8, chTemp);
						//strcat(chFileName8, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName8);
						//
						////output the model vector as an nii.gz image
						//zxhImageData zxhModelpointsImg;
						//zxhModelpointsImg.NewImage(m_pBaseImgInfo);
						//MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
						//int nLen2 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName2 = (char *)malloc(nLen2);
						//strcpy(chFileName2, m_cResultPath);
						//strcat(chFileName2, "/ModelVector");
						//strcat(chFileName2, chTemp);
						//strcat(chFileName2, ".nii.gz");
						//zxh::SaveImage( &zxhModelpointsImg, chFileName2 );
						//free(chFileName2);
						//int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName5 = (char *)malloc(nLen5);
						//strcpy(chFileName5, m_cResultPath);
						//strcat(chFileName5, "/NBUValue");
						//strcat(chFileName5, chTemp);
						//strcat(chFileName5, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName5);
						//int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName6 = (char *)malloc(nLen6);
						//strcpy(chFileName6, m_cResultPath);
						//strcat(chFileName6, "/NBPValue");
						//strcat(chFileName6, chTemp);
						//strcat(chFileName6, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName6);
						//int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
						//char *chFileName7 = (char *)malloc(nLen7);
						//strcpy(chFileName7, m_cResultPath);
						//strcat(chFileName7, "/NBSValue");
						//strcat(chFileName7, chTemp);
						//strcat(chFileName7, ".nii.gz");
						//zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
						//free(chFileName7);


						return m_nFmItr; 
					}
				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);//keeping in order in narrowband.
				} 
				//add by jdq
			}
		}//for


		float nEndDistmm2 = CalcDistmm2(m_fSMEndPointWorld,dfMinNodeWorld);

		// reach around the end point of the model 
		//if (nEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//	bIsCurSMEndPoint=true;
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum&&m_nSum>iMaxSumSelNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
int miiMinPathModel::ModifiedFastMarchingEvolutiondontMoveModelWithLineMask(int nMaxItrNum,int nItrNumSum,const short *sImgData,const short *sLineMaskData)
{
	// define a variable for the minimum in Narrow Band
	miiCNode<double> iMinNode, iNewNode;
	
	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };

	m_nFmItr = 0;
	m_nFMMEvlNum++;
	cout<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";//Add by JDQ
	// iterate for FMM
	while (true)
	{
		m_nFmItr++;

		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);
		miiCNode<double,float> dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;
		//sAMinNodeImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iAlivePValue;//Add by JDQ
		m_nSum=nItrNumSum+m_nFmItr;
		//out put the intermediate results
		/*if((m_nSum)%m_nSaveFmItr==0&&(nItrNumSum+m_nFmItr)<=200000)
		{
		char chTemp[25];
		_itoa_s(nItrNumSum+m_nFmItr, chTemp, 10);
		int nLen2 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath") + strlen(chTemp) + strlen(".nii.gz") + 1;
		char *chFileName4 = (char *)malloc(nLen2);
		strcpy(chFileName4, m_cResultPath);
		strcat(chFileName4, "/MinNodeNBMinPath");
		strcat(chFileName4, chTemp);
		strcat(chFileName4, ".nii.gz");
		SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
		free(chFileName4);
		}*/
	
		/*float n=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		float n1=CalcDistmm2(m_vModelPointsWorld[1],m_vModelPointsWorld[2]);*/
	
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/withintensityupdated/m_dWP.txt",ios::app);//Add by JDQ
	 // WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," <<iMinNode.y << "," << iMinNode.z<<"," <<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<iMinNode.val<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		//if(m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]!=0)
		//	double m=0;
		//ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	 //WriteFileTxt<< "iMinNode"<<m_nFmItr<<"("<<iMinNode.x << "," << iMinNode.y << "," << iMinNode.z<<","<<m_dSims[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_sVes[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dSimd[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_dP2[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<< ","<<m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<","<<m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];
			//short m=sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE&&sLineMaskData[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]!=0)
			{
				// for fist direction in aorta
				float nInitDsit2 = CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

				float nMinNodeDsit2=CalcDistmm2(m_vModelPointsWorld[0],dfMinNodeWorld);
	
				float costheata;
				// the candidate point in aorta
				if (m_nFMMEvlNum == 1 && nMinNodeDsit2 < (nInitDsit2 * 4) && !m_bFistEvlFlg)
				{
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz, true);
				}
				else
				{
					m_bFistEvlFlg = true;
					// calculate the theata between model and vessel
					costheata = CalcTheata(nx, ny, nz);
				}

				// upwind
				//UpWindPlusTime(nx, ny, nz, theata);
				//UpWindPlus(nx, ny, nz, theata);
				//NewUpWindPlus(nx, ny, nz, costheata);
				//NewUpWindTime(nx, ny, nz, costheata);
				/*if(!IsEndPoint)
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				else
				ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);*/
				//ModifiedUpWindPlusTimeMixed(nx, ny, nz, costheata,sImgData);
				ModifiedUpWindPlusTime(nx, ny, nz, costheata,sImgData);
				//ModifiedUpWindPlusTimeLast(nx, ny, nz, costheata,sImgData);
				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					
						m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
					//sAMinNodeImg[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]=m_iActivePValue;//Add by JDQ
					
					bool blSegPointFind=false;

	              
	         		float  CurMVector[4]={0};
					float  CurVSLVector[4]={0};
	
					 CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
					 CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
					 CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
					 CurVSLVector[0]=m_iLastStartPointWorld.x-m_iStartPointWorld.x;
					 CurVSLVector[1]=m_iLastStartPointWorld.y-m_iStartPointWorld.y;
					 CurVSLVector[2]=m_iLastStartPointWorld.z-m_iStartPointWorld.z;
					 CurVSLVector[3]=CalcDistmm2(m_iLastStartPointWorld,m_iStartPointWorld);
					if (m_nFMMEvlNum == 1)
					{ /* 
					    CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
					   	m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
						*/
						/*CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);

                          m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);
						
					}
					else
					{  /*
					    CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						 m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
						*/ 
						/*CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/

						  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
						  m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);
						
					}
					
					if (blSegPointFind) 
					{

						return m_nFmItr;
					}
					/*if (m_nFMMEvlNum == 1)
					{
                          m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVec1Distmm2,nAmm2);
					}
					else
					{
						  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,m_nModelVecODistmm2,nAmm2);
					}
					
					if (blSegPointFind) 
					{
						return m_nFmItr;
					}*/

				}
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);
				}
			}
		}//for

		float nModelEndDistmm2 = CalcDistmm2(m_iModelEndPointWorld,dfMinNodeWorld);
		
		// reach around the end point of the model 
		//if (nModelEndDistmm2 < m_nMinDistmm2)//if (nModelEndDist2 < (20 * 20))//Change by JDQ
		//{
		//
		//				
		//	return m_nFmItr;
		//}

		if (m_nFmItr > m_nImgWX * m_nImgWY * (m_nImgWZ - 1))
		{
			return m_nFmItr;
		}

		if (m_nFmItr > nMaxItrNum)
		{
			return m_nFmItr;
		}
	}//while

	return m_nFmItr;
}
// Function Name: UpWind()
//
// Parameters: x, y, z: coordinate corespond left(right),up(down),front(back)
//
// Description: add the direction to the cost function
//
// Returns: 
//
void miiMinPathModel::UpWind(int x, int y, int z, double theata)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// add the direction to the cost function
	double dTheataMean = 1, dTheataStdDev = 0.5;

	if (theata <= 0)
	{
		double dSim = - (((0.1 - dTheataMean)*(0.1 - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] / (dSim * dSim);
	}
	else if (theata > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		double dSim = - (((0.866 - dTheataMean)*(0.866 - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] / (dSim * dSim);
	}
	else
	{
		double dSim = - (((theata - dTheataMean)*(theata - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] / (dSim * dSim);//(dSim)^2?
	}

	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	ofstream WriteFileTxt;
	  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/original/potentialnotchangeinupwind/m_dWP.txt",ios::app);//Add by JDQ
	  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ")"<<"\n ";//add by JDQ
	 WriteFileTxt.close();
	if (m_nSum==500000)
		 double m2=1;
}
// Function Name: UpWind()
//
// Parameters: x, y, z: coordinate corespond left(right),up(down),front(back)
//
// Description: add the direction to the cost function
//
// Returns: 
//
void miiMinPathModel::UpWind1(int x, int y, int z, double theata)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// add the direction to the cost function
	double dTheataMean = 1, dTheataStdDev = 0.5;

	if (theata <= 0)
	{
		double dSim = - (((0.1 - dTheataMean)*(0.1 - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] / (dSim * dSim);
		
      
	
		
		}
	}
	else if (theata > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		double dSim = - (((0.866 - dTheataMean)*(0.866 - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] / (dSim * dSim);
			
		}
	}
	else
	{
		double dSim = - (((theata - dTheataMean)*(theata - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] / (dSim * dSim);
			
		}
	}
	
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	ofstream WriteFileTxt;
	  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/original/potentialnotchangeinupwind/m_dWP.txt",ios::app);//Add by JDQ
	  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ")"<<"\n ";//add by JDQ
	 WriteFileTxt.close();
	if (m_nSum==500000)
		 double m2=1;
}
// Function Name: UpWind()
//
// Parameters: x, y, z: coordinate corespond left(right),up(down),front(back)
//
// Description: add the direction to the cost function
//
// Returns: 
//
void miiMinPathModel::UpWindPlus(int x, int y, int z, double costheta)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// add the direction to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5;
	float lambdad=1;
	if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSim = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSim = exp(dSim);
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] + lambdad * dSim );
	 
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
}
void miiMinPathModel::UpWindPlusTime(int x, int y, int z, double theata)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// add the direction to the cost function
	double dTheataMean = 1, dTheataStdDev = 0.5;

	if (theata <= 0)
	{
		double dSim = - (((0.1 - dTheataMean)*(0.1 - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			1/(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+2*dSim);
	}
	else if (theata > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		double dSim = - (((0.866 - dTheataMean)*(0.866 - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			1/(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+2*dSim);
	}
	else
	{
		double dSim = - (((theata - dTheataMean)*(theata - dTheataMean)) / (dTheataStdDev * dTheataStdDev)) / 2;
		dSim = exp(dSim);
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			1/(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+2*dSim);
	}

	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
}

	void miiMinPathModel::NewUpWindPlus(int x, int y, int z, double costheta)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// add the direction to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5;
	float lambdad=1;
	if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSim = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSim = exp(dSim);

	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(m_sUnseenImgIntensity[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+ m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+dSim);
	double m=1/(m_sUnseenImgIntensity[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+ m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+dSim);
	
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
}

	void miiMinPathModel::ModifiedUpWindTime(int x, int y, int z, double costheta,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.00000001,omiga=0.0001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;double dSimd=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(dImgVal==-104)
		double s=1;
	if ( dImgVal < m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta <= 0)
	{
		dSimd = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
		}	
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		(1/((dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+epsilion)))*(1/((dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+epsilion)));
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}
	else if (costheta > 0.9397)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd = - (((0.9397 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);	
		}
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	(1/((dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+epsilion)))*(1/((dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+epsilion)));
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}
	else
	{
		dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
		}
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		(1/((dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+epsilion)))*(1/((dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+epsilion)));
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}//(1-lambdas)*(dSimd*dSimd)
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	if(m<0)
		int m=0;
	if(dSimd<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
// Function Name: NewUpWindPlus1()
//
// Parameters: x, y, z: coordinate corespond left(right),up(down),front(back)
//
// Description: calculate the intensity and direction parameter to the potential calculation
//
// Returns: 
//
void miiMinPathModel::ModifiedUpWindPlusTime(int x, int y, int z, double costheta,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;double dSimd=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	
	if ( dImgVal < m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta <= 0)
	{
		dSimd = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
		}	
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}
	else if (costheta > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);	
		}
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}
	else
	{
		dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
		}
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}//(1-lambdas)*(dSimd*dSimd)
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	if(m<0)
		int m=0;
	if(dSimd<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}

	// Function Name: NewUpWindPlus1()
//
// Parameters: x, y, z: coordinate corespond left(right),up(down),front(back)
//
// Description: calculate the intensity and direction parameter to the potential calculation
//
// Returns: 
//
void miiMinPathModel::ModifiedUpWindPlusTime2Costheta(int x, int y, int z, float costheta1,float costheta2,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;double dSimd=0;double dSimd1=0;double dSimd2=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}//(1-lambdas)*(dSimd*dSimd)
	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta2 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}//(1-lambdas)*(dSimd*dSimd)
	//dSimd=(dSimd1+dSimd2)/2;
	dSimd=dSimd1*dSimd1+dSimd2*dSimd2;

	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	if(m<0)
		int m=0;
	if(dSimd1<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}

	void miiMinPathModel::ModifiedUpWind_WOBP(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=1;
	double dSims=0;
	double dSimd=0;
	double dSimd1=0;
	double dSimd2=0;
	double dSimd3=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	//else if (costheta1 > 0.5)//cos(theta)=0.866  --->theta=60 degree.
	{
		//dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}


	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta2 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}

	if (costheta3 == 0)//current vessel director and last start vector
	{
		dSimd3 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else if (abs(costheta3) > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd3 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else
	{
		dSimd3 = - (((abs(costheta3) - dCosTheataMean)*(abs(costheta3) - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
		
	}
	//calculate the total direction parameter
	dSimd=dSimd1*dSimd1;
	//dSimd=dSimd3*dSimd3;
	//dSimd=dSimd1*dSimd1;
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	if(m<0)
		int m=0;
	if(dSimd1<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
void miiMinPathModel::ModifiedUpWindPlusTime2CosthetaN(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;
	double dSimd=0;
	double dSimd1=0;
	double dSimd2=0;
	double dSimd3=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	//else if (costheta1 > 0.5)//cos(theta)=0.866  --->theta=60 degree.
	{
		//dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}


	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta2 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}

	if (costheta3 == 0)//current vessel director and last start vector
	{
		dSimd3 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else if (abs(costheta3) > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd3 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else
	{
		dSimd3 = - (((abs(costheta3) - dCosTheataMean)*(abs(costheta3) - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
		
	}
	//calculate the total direction parameter
	dSimd=dSimd1*dSimd1;
	//dSimd=dSimd3*dSimd3;
	//dSimd=dSimd1*dSimd1;	
	int v=m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;

	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	double m=sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]);

	if(m<0)
		int m=0;
	if(dSimd1<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
void miiMinPathModel::ModifiedUpWindPlusTime2CosthetaN_SP(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	double ub = min(u_u, u_d);//y
	double uc = min(u_b, u_f);//z
	double ddela=0;
	double ddelb=0;
	double ddelc=0;
	double ddelmax=0;
	double ddelmed=0;
	double ddelmin=0;

	double ddeltx=m_nImgSpacing[0];
	double ddelty=m_nImgSpacing[1];
	double ddeltz=m_nImgSpacing[2];

	double dMax, dMin, dMid;
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dMin = ub;
		ddelmin=ddelty;

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dMin = ua;
		ddelmin=ddeltx;
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dMax = uc;
		ddelmax=ddeltz;
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dMin = uc;
		ddelmin=ddeltz;
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	ddela=ddelmax; ddelb=ddelmed;ddelc=ddelmin;

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;
	double dSimd=0;
	double dSimd1=0;
	double dSimd2=0;
	double dSimd3=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}


	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta2 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}

	if (costheta3 == 0)//current vessel director and last start vector
	{
		dSimd3 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else if (abs(costheta3) > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd3 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else
	{
		dSimd3 = - (((abs(costheta3) - dCosTheataMean)*(abs(costheta3) - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
		
	}
	//calculate the total direction parameter
	dSimd=dSimd1*dSimd1;
	//dSimd=dSimd3*dSimd3;
	//dSimd=dSimd1*dSimd1;
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	//solve 3 x
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(da3, db3, dc3, dRoots_1);
	//solve 2 x
	double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	bool bSolv_2 = QuadraticRoots(da2, db2, dc2, dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > ua && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > ub && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + uc;
	}
  //cout<< "("<<ddela << "," << ddelb << "," << ddelc<< ")"<<"\n ";
}
	void miiMinPathModel::ModUpWind_MultiST_test(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	// get the intensity to the cost function
	float dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=600;
	float dSims=0;
	float dSimd=0;
	float dSimd1=0;
	float dImgVal = (float)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}
	//calculate the total direction parameter
	dSimd=dSimd1*dSimd1;
	//DFM potential
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	float vsls=vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	//VBMP potential
	/*float dSimv1=vslweight*m_sNormVes[y * m_nImgWX + x]+0.01;
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(dSimv1+epsilion))*(1/(dSimv1+epsilion));
	}	
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimv1;*/
	float dP2=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	float dP=sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]);
	if(dP2>2500)
		int xx=7;
	//Multi Stencils
	MultiStencils(x,y,z,dP2);
	
  //cout<< "("<<ddela << "," << ddelb << "," << ddelc<< ")"<<"\n ";
}
	void miiMinPathModel::ModUpWind_MultiST(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;
	double dSimd=0;
	double dSimd1=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}
	//calculate the total direction parameter
	dSimd=dSimd1*dSimd1;
	//DFM potential
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	//VBMP potential
	/*float dSimv1=vslweight*m_sNormVes[y * m_nImgWX + x]+0.01;
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(dSimv1+epsilion))*(1/(dSimv1+epsilion));
	}	
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimv1;*/
	double dP2=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(dP2>2500)
		int xx=7;
	//Multi Stencils
	MultiStencils(x,y,z,dP2);
	
  //cout<< "("<<ddela << "," << ddelb << "," << ddelc<< ")"<<"\n ";
}
	void miiMinPathModel::ModUpWind_MultiST_VslsOnly(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;
	double dSimd=0;
	double dSimd1=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}
	//calculate the total direction parameter
	dSims = 1;
   dSimd=1;
	//DFM potential
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	//VBMP potential
	/*float dSimv1=vslweight*m_sNormVes[y * m_nImgWX + x]+0.01;
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(dSimv1+epsilion))*(1/(dSimv1+epsilion));
	}	
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimv1;*/
	double dP2=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(dP2>2500)
		int xx=7;
	//Multi Stencils
	MultiStencils(x,y,z,dP2);
	
  //cout<< "("<<ddela << "," << ddelb << "," << ddelc<< ")"<<"\n ";
}
void miiMinPathModel::FormTria(miiCNode<float,int> mLN,miiCNode<float,int> mRN,miiCNode<float,int> mUN,miiCNode<float,int> mDN,miiCNode<float,int> mFN,miiCNode<float,int> mBN,miiCNode<float,int>&mMaxN, miiCNode<float,int>&mMedN,miiCNode<float,int> &mMinN)
{
	// find the points with minmal value on each axix
	miiCNode<float,int> mLRN,mUDN,mFBN;
	if(mLN.val>mRN.val)
	{
		mLRN=mRN;
	}
	else
	{
		mLRN=mLN;
	}
	if(mUN.val>mDN.val)
	{
		mUDN=mDN;
	}
	else
	{
		mUDN=mUN;
	}	
	if(mFN.val>mBN.val)
	{
		mFBN=mBN;
	}
	else
	{
		mFBN=mFN;
	}
	//order the three minimal points
	if(mLRN.val>mUDN.val)
	{
		mMaxN=mLRN;
		mMinN=mUDN;
	}
	else
	{
		mMaxN=mUDN;
		mMinN=mLRN;
	}
	if(mFBN.val>mMaxN.val)
	{
		mMedN=mMaxN;
		mMaxN=mFBN;
	}
	else if(mFBN.val>mMinN.val)
	{
		mMedN=mFBN;
	}
	else
	{
		mMedN=mMinN;
		mMinN=mFBN;
	}

}
void miiMinPathModel::CalVecDelta(miiCNode<float,int> mCN,miiCNode<float,int>mMaxN, miiCNode<float,int>mMedN,miiCNode<float,int> mMinN,float fDel[3],float fCos[3])
{
	float daP[3]={0,0,0};
	float dbP[3]={0,0,0};
	float dcP[3]={0,0,0};
	daP[0]=mMaxN.x;daP[1]=mMaxN.y;daP[2]=mMaxN.z;
	dbP[0]=mMedN.x;dbP[1]=mMedN.y;dbP[2]=mMedN.z;
	dcP[0]=mMinN.x;dcP[1]=mMinN.y;dcP[2]=mMinN.z;
	//
	float daPworld[3]={daP[0],daP[1],daP[2]};
	float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	float dxyzPworld[3]={mCN.x,mCN.y,mCN.z};
	
	m_pBaseImgInfo->ImageToWorld(daPworld);
	m_pBaseImgInfo->ImageToWorld(dbPworld);
	m_pBaseImgInfo->ImageToWorld(dcPworld);
	m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	// a b c as the vector
	float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	float ddela=zxh::VectorOP_Magnitude(dva,3);
	float ddelb=zxh::VectorOP_Magnitude(dvb,3);
	float ddelc=zxh::VectorOP_Magnitude(dvc,3);
	float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	float cc=zxh::VectorOP_Cosine(dva,dvc,3);

	//store delta into arrays
	fDel[0]=ddela;
	fDel[1]=ddelb;
	fDel[2]=ddelc;
	//store cosine values into arrays
	fCos[0]=ca;
	fCos[1]=cb;
	fCos[2]=cc;
}
void miiMinPathModel::MultiStencils(int x, int y, int z,float dP2)
{
	// 6 stencils including all the neibour points
	float u[6]={0,0,0,0,0,0};
	u[0]=Stencil1_2(x,y,z,dP2);
	//cout<<"Stencil 1 down!"<<endl;
	//u[1]=Stencil2_2(x,y,z,dP2);
	//cout<<"Stencil 2 down!"<<endl;
	//u[2]=Stencil3_2(x,y,z,dP2);
	//cout<<"Stencil 3 down!"<<endl;
	//u[3]=Stencil4_2(x,y,z,dP2);
	//cout<<"Stencil 4 down!"<<endl;
	//u[4]=Stencil5_2(x,y,z,dP2);
	//cout<<"Stencil 5 down!"<<endl;
	//u[5]=Stencil6_2(x,y,z,dP2);
	//cout<<"Stencil 6 down!"<<endl;
	float fu=10000000000;
	for (int i=0;i<1;i++)
	{
		if(u[i]<fu)
		fu=u[i];
	}
	if(fu>100)
	{
		fu=100;
	}
	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=fu;
	int m=0;
}

double miiMinPathModel::Stencil1_2(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y + 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y - 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;
	if (u_l_x < 0 || u_l_x >= m_nImgWX)
	{    
		if(u_l_x < 0)u_l_x=x;
		if(u_l_x >= m_nImgWX)u_l_x=x;
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
	{
		if(u_r_x < 0)u_r_x=x;
		if(u_r_x >= m_nImgWX)u_r_x=x;
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		{
		if(u_u_y < 0)u_u_y=y;
		if(u_u_y >= m_nImgWY)u_u_y=y;
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		{
		if(u_d_y < 0)u_d_y=y;
		if(u_d_y >= m_nImgWY)u_d_y=y;
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		{
		if(u_b_z < 0)u_b_z=z;
		if(u_b_z >= m_nImgWZ)u_b_z=z;
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		{
		if(u_f_z < 0)u_f_z=z;
		if(u_f_z >= m_nImgWZ)u_f_z=z;
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	miiCNode<float,int> mCN;
	mCN.x=x;
	mCN.y=y;
	mCN.z=z;
	//store the six points with u
	miiCNode<float,int> mLN,mRN,mUN,mDN,mFN,mBN;
	//left point
	mLN.x=u_l_x;
	mLN.y=u_l_y;
	mLN.z=u_l_z;
	mLN.val=u_l;
		//right point
	mRN.x=u_r_x;
	mRN.y=u_r_y;
	mRN.z=u_r_z;
	mRN.val=u_r;
		//up point
	mUN.x=u_u_x;
	mUN.y=u_u_y;
	mUN.z=u_u_z;
	mUN.val=u_u;
		//down point
	mDN.x=u_d_x;
	mDN.y=u_d_y;
	mDN.z=u_d_z;
	mDN.val=u_d;
		//front point
	mFN.x=u_f_x;
	mFN.y=u_f_y;
	mFN.z=u_f_z;
	mFN.val=u_f;
		//back point
	mBN.x=u_b_x;
	mBN.y=u_b_y;
	mBN.z=u_b_z;
	mBN.val=u_b;
	//order the max med and min; generate a triagle
	miiCNode<float,int>mMaxN,mMedN,mMinN;
	FormTria(mLN,mRN,mUN,mDN,mFN,mBN,mMaxN,mMedN,mMinN);
	//calculate the three vectors and deltas;
	float fDelta[3]={0,0,0};
	float fCos[3]={0,0,0};
	CalVecDelta(mCN,mMaxN,mMedN,mMinN,fDelta,fCos);
	////select the minimal value of each axis;
	//double ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	//double dualr[3]={0,0,0};//ua,dualr,dua are in one group
	//if(u_l>u_r)
	//{
	//	ua=u_r;
	//	dualr[0]=u_r_x;
	//	dualr[1]=u_r_y;
	//	dualr[2]=u_r_z;
	//}
	//else
	//{
	//	ua=u_l;
	//	dualr[0]=u_l_x;
	//	dualr[1]=u_l_y;
	//	dualr[2]=u_l_z;
	//}
	//double ub = min(u_u, u_d);//y
	//double dubud[3]={0,0,0};//ub,dubud,dub are in one group
	//if(u_u>u_d)
	//{
	//	ub=u_d;
	//	dubud[0]=u_d_x;
	//	dubud[1]=u_d_y;
	//	dubud[2]=u_d_z;
	//}
	//else
	//{
	//	ub=u_u;
	//	dubud[0]=u_u_x;
	//	dubud[1]=u_u_y;
	//	dubud[2]=u_u_z;
	//}
	//double uc = min(u_b, u_f);//z
	//double ducbf[3]={0,0,0};//ub,dubud,dub are in one group
	//if(u_b>u_f)
	//{
	//	uc=u_f;
	//	ducbf[0]=u_f_x;
	//	ducbf[1]=u_f_y;
	//	ducbf[2]=u_f_z;
	//}
	//else
	//{
	//	ub=u_u;
	//	ducbf[0]=u_b_x;
	//	ducbf[1]=u_b_y;
	//	ducbf[2]=u_b_z;
	//}
	////order the three minimal values and store their positions
	//double ddela1=0;
	//double ddelb1=0;
	//double ddelc1=0;
	//double ddelmax=0;
	//double ddelmed=0;
	//double ddelmin=0;

	//double ddeltx=m_nImgSpacing[0];
	//double ddelty=m_nImgSpacing[1];
	//double ddeltz=m_nImgSpacing[2];

	//double dMax, dMin, dMid;//the max u values
	//double dmaxP[3]={0,0,0};
	//double dmedP[3]={0,0,0};
	//double dminP[3]={0,0,0};
	//double daP[3]={0,0,0};
	//double dbP[3]={0,0,0};
	//double dcP[3]={0,0,0};
	//if (ua > ub) 
	//{
	//	dMax = ua;
	//	ddelmax=ddeltx;
	//	dmaxP[0]=dualr[0];
	//	dmaxP[1]=dualr[1];
	//	dmaxP[2]=dualr[2];
	//	dMin = ub;
	//	ddelmin=ddelty;
	//	dminP[0]=dubud[0];
	//	dminP[1]=dubud[1];
	//	dminP[2]=dubud[2];

	//}
	//else
	//{
	//	dMax = ub;
	//	ddelmax=ddelty;
	//	dmaxP[0]=dubud[0];
	//	dmaxP[1]=dubud[1];
	//	dmaxP[2]=dubud[2];
	//	dMin = ua;
	//	ddelmin=ddeltx;
	//	dminP[0]=dualr[0];
	//	dminP[1]=dualr[1];
	//	dminP[2]=dualr[2];
	//}
	//if (uc > dMax)
	//{

	//	dMid = dMax;
	//	ddelmed=ddelmax;
	//	dmedP[0]=dmaxP[0];
	//	dmedP[1]=dmaxP[1];
	//	dmedP[2]=dmaxP[2];
	//	dMax = uc;
	//	ddelmax=ddeltz;
	//	dmaxP[0]=ducbf[0];
	//	dmaxP[1]=ducbf[1];
	//	dmaxP[2]=ducbf[2];
	//}
	//else if (uc > dMin)
	//{
	//	dMid = uc;
	//	ddelmed=ddeltz;
	//	dmedP[0]=ducbf[0];
	//	dmedP[1]=ducbf[1];
	//	dmedP[2]=ducbf[2];
	//}
	//else
	//{
	//	dMid = dMin;
	//	ddelmed=ddelmin;
	//	dmedP[0]=dminP[0];
	//	dmedP[1]=dminP[1];
	//	dmedP[2]=dminP[2];
	//	dMin = uc;
	//	ddelmin=ddeltz;
	//	dminP[0]=ducbf[0];
	//	dminP[1]=ducbf[1];
	//	dminP[2]=ducbf[2];
	//}

	//ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	//ddela1=ddelmax; 
	// ddelb1=ddelmed;
	// ddelc1=ddelmin;
	//daP[0]=dmaxP[0];daP[1]=dmaxP[1];daP[2]=dmaxP[2];
	//dbP[0]=dmedP[0];dbP[1]=dmedP[1];dbP[2]=dmedP[2];
	//dcP[0]=dminP[0];dcP[1]=dminP[1];dcP[2]=dminP[2];
	////
	//float daPworld[3]={daP[0],daP[1],daP[2]};
	//float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	//float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	//float dxyzPworld[3]={x,y,z};
	//
	//m_pBaseImgInfo->ImageToWorld(daPworld);
	//m_pBaseImgInfo->ImageToWorld(dbPworld);
	//m_pBaseImgInfo->ImageToWorld(dcPworld);
	//m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	//// a b c as the vector
	//float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	//float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	//float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	//double ddela=zxh::VectorOP_Magnitude(dva,3);
	//double ddelb=zxh::VectorOP_Magnitude(dvb,3);
	//double ddelc=zxh::VectorOP_Magnitude(dvc,3);
	//float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	//float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	//float cc=zxh::VectorOP_Cosine(dva,dvc,3);
	////define the three unit vectors
	//float dva1[3]={0,0,0};
	//dva1[0]=(float)daP[0]-x;
	//dva1[1]=(float)daP[1]-y;
	//dva1[2]=(float)daP[2]-z;
	//float dvb1[3]={0,0,0};
	//dvb1[0]=(float)dbP[0]-x;
	//dvb1[1]=(float)dbP[1]-y;
	//dvb1[2]=(float)dbP[2]-z;
	//float dvc1[3]={0,0,0};
	//dvc1[0]=(float)dcP[0]-x;
	//dvc1[1]=(float)dcP[1]-y;
	//dvc1[2]=(float)dcP[2]-z;
	////normalize
	//zxh::VectorOP_Normalise(dva1,3);
	//zxh::VectorOP_Normalise(dvb1,3);
	//zxh::VectorOP_Normalise(dvc1,3);
	////
	//float ca1=zxh::VectorOP_Cosine(dva1,dvb1,3);
	//float cb1=zxh::VectorOP_Cosine(dvb1,dvc1,3);
	//float cc1=zxh::VectorOP_Cosine(dva1,dvc1,3);
	////
     float ua=mMaxN.val;
	 float ub=mMedN.val;
	 float uc=mMinN.val;
	 float ca=fCos[0];
	 float cb=fCos[1];
	 float cc=fCos[2];
	 float ddela=fDelta[0];
	 float ddelb=fDelta[1];
	 float ddelc=fDelta[2];
	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);

	//double dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	//double dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	//double dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	////solve the equation
	//double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	//bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	////solve 2 x
	//double dA2=(1-pow(cc,2))/pow(ddelb,2)+1/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc);
	//double dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+1/pow(ddelc,2)*pow(uc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	//double dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	//double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	//if (dRtMax_3 > ua && bSolv_3 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_3;
	//}
	//else if (dRtMax_2 > ub && bSolv_2 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
	////solve 3 x
	//double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//double dRoots_0[2] = {0}, dRoots_1[2] = {0};
	//
	//bool bSolv_0 = QuadraticRoots(da3, db3, dc3, dRoots_0);
	////solve 2 x
	//double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//bool bSolv_1 = QuadraticRoots(da2, db2, dc2, dRoots_1);

	//double dRtMax_0 = max(dRoots_0[0], dRoots_0[1]);
	//double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);

	//if (dRtMax_0 > ua && bSolv_0 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_0;
	//}
	//else if (dRtMax_1 > ub && bSolv_1 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
}
double miiMinPathModel::Stencil2(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y + 1, u_u_z = z+1;
	int u_d_x = x,     u_d_y = y - 1, u_d_z = z-1;
	int u_f_x = x,     u_f_y = y-1,     u_f_z = z+1;
	int u_b_x = x,     u_b_y = y+1,     u_b_z = z-1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	//select the minimal value of each axis;
	//double ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	float dualr[3]={0,0,0};//ua,dualr,dua are in one group
	float ua=0;
	float ub=0;
	float uc=0;
	if(u_l>u_r)
	{
		ua=u_r;
		dualr[0]=u_r_x;
		dualr[1]=u_r_y;
		dualr[2]=u_r_z;
	}
	else
	{
		ua=u_l;
		dualr[0]=u_l_x;
		dualr[1]=u_l_y;
		dualr[2]=u_l_z;
	}
	//double ub = min(u_u, u_d);//y
	float dubud[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_u>u_d)
	{
		ub=u_d;
		dubud[0]=u_d_x;
		dubud[1]=u_d_y;
		dubud[2]=u_d_z;
	}
	else
	{
		ub=u_u;
		dubud[0]=u_u_x;
		dubud[1]=u_u_y;
		dubud[2]=u_u_z;
	}
	
	float ducbf[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_b>u_f)
	{
		uc=u_f;
		ducbf[0]=u_f_x;
		ducbf[1]=u_f_y;
		ducbf[2]=u_f_z;
	}
	else
	{
		ub=u_u;
		ducbf[0]=u_b_x;
		ducbf[1]=u_b_y;
		ducbf[2]=u_b_z;
	}
	//order the three minimal values and store their positions
	float ddelmax=0;
	float ddelmed=0;
	float ddelmin=0;
	float ddeltx=m_nImgSpacing[0];
	float ddelty=sqrt(pow(m_nImgSpacing[1],2)+pow(m_nImgSpacing[2],2));
	float ddeltz=sqrt(pow(m_nImgSpacing[1],2)+pow(m_nImgSpacing[2],2));
	float dMax, dMin, dMid;//the max u values
	float dmaxP[3]={0,0,0};
	float dmedP[3]={0,0,0};
	float dminP[3]={0,0,0};
	float daP[3]={0,0,0};
	float dbP[3]={0,0,0};
	float dcP[3]={0,0,0};
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dmaxP[0]=dualr[0];
		dmaxP[1]=dualr[1];
		dmaxP[2]=dualr[2];
		dMin = ub;
		ddelmin=ddelty;
		dminP[0]=dubud[0];
		dminP[1]=dubud[1];
		dminP[2]=dubud[2];

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dmaxP[0]=dubud[0];
		dmaxP[1]=dubud[1];
		dmaxP[2]=dubud[2];
		dMin = ua;
		ddelmin=ddeltx;
		dminP[0]=dualr[0];
		dminP[1]=dualr[1];
		dminP[2]=dualr[2];
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dmedP[0]=dmaxP[0];
		dmedP[1]=dmaxP[1];
		dmedP[2]=dmaxP[2];
		dMax = uc;
		ddelmax=ddeltz;
		dmaxP[0]=ducbf[0];
		dmaxP[1]=ducbf[1];
		dmaxP[2]=ducbf[2];
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
		dmedP[0]=ducbf[0];
		dmedP[1]=ducbf[1];
		dmedP[2]=ducbf[2];
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dmedP[0]=dminP[0];
		dmedP[1]=dminP[1];
		dmedP[2]=dminP[2];
		dMin = uc;
		ddelmin=ddeltz;
		dminP[0]=ducbf[0];
		dminP[1]=ducbf[1];
		dminP[2]=ducbf[2];
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	daP[0]=dmaxP[0];daP[1]=dmaxP[1];daP[2]=dmaxP[2];
	dbP[0]=dmedP[0];dbP[1]=dmedP[1];dbP[2]=dmedP[2];
	dcP[0]=dminP[0];dcP[1]=dminP[1];dcP[2]=dminP[2];
	//
	float daPworld[3]={daP[0],daP[1],daP[2]};
	float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	float dxyzPworld[3]={x,y,z};
	
	m_pBaseImgInfo->ImageToWorld(daPworld);
	m_pBaseImgInfo->ImageToWorld(dbPworld);
	m_pBaseImgInfo->ImageToWorld(dcPworld);
	m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	// a b c as the vector
	float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	float ddela=zxh::VectorOP_Magnitude(dva,3);
	float ddelb=zxh::VectorOP_Magnitude(dvb,3);
	float ddelc=zxh::VectorOP_Magnitude(dvc,3);
	float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	float cc=zxh::VectorOP_Cosine(dva,dvc,3);
	//

	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);

	//double dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	//double dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	//double dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	////solve the equation
	//double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	//bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	////solve 2 x
	//double dA2=(1-pow(cc,2))/pow(ddelb,2)+1/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc);
	//double dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+1/pow(ddelc,2)*pow(uc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	//double dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	//double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	//if (dRtMax_3 > ua && bSolv_3 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_3;
	//}
	//else if (dRtMax_2 > ub && bSolv_2 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
	////solve 3 x
	//double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//double dRoots_0[2] = {0}, dRoots_1[2] = {0};
	//
	//bool bSolv_0 = QuadraticRoots(da3, db3, dc3, dRoots_0);
	////solve 2 x
	//double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//bool bSolv_1 = QuadraticRoots(da2, db2, dc2, dRoots_1);

	//double dRtMax_0 = max(dRoots_0[0], dRoots_0[1]);
	//double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);

	//if (dRtMax_0 > ua && bSolv_0 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_0;
	//}
	//else if (dRtMax_1 > ub && bSolv_1 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
}
double miiMinPathModel::Stencil2_2(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y + 1, u_u_z = z+1;
	int u_d_x = x,     u_d_y = y - 1, u_d_z = z-1;
	int u_f_x = x,     u_f_y = y-1,     u_f_z = z+1;
	int u_b_x = x,     u_b_y = y+1,     u_b_z = z-1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	/*if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];*/
		if (u_l_x < 0 || u_l_x >= m_nImgWX)
	{    
		if(u_l_x < 0)u_l_x=x;
		if(u_l_x >= m_nImgWX)u_l_x=x;
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
	{
		if(u_r_x < 0)u_r_x=x;
		if(u_r_x >= m_nImgWX)u_r_x=x;
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY||u_u_z < 0 || u_u_z >= m_nImgWZ)
		{
		if(u_u_y < 0)u_u_y=y;
		if(u_u_y >= m_nImgWY)u_u_y=y;
		if(u_u_z < 0)u_u_z=z;
		if(u_u_z >= m_nImgWZ)u_u_z=z;
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY||u_d_z < 0 || u_d_z >= m_nImgWZ)
		{
		if(u_d_y < 0)u_d_y=y;
		if(u_d_y >= m_nImgWY)u_d_y=y;
		if(u_d_z < 0)u_d_z=z;
		if(u_d_z >= m_nImgWZ)u_d_z=z;
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ||u_b_y < 0 || u_b_y >= m_nImgWY)
		{
		if(u_b_z < 0)u_b_z=z;
		if(u_b_z >= m_nImgWZ)u_b_z=z;
		if(u_b_y < 0)u_b_y=y;
		if(u_b_y >= m_nImgWY)u_b_y=y;
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ||u_f_y < 0 || u_f_y >= m_nImgWY)
		{
		if(u_f_z < 0)u_f_z=z;
		if(u_f_z >= m_nImgWZ)u_f_z=z;
		if(u_f_y < 0)u_f_y=y;
		if(u_f_y >= m_nImgWY)u_f_y=y;
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];
	//store the center point
	miiCNode<float,int> mCN;
	mCN.x=x;
	mCN.y=y;
	mCN.z=z;
	//store the six points with u
	miiCNode<float,int> mLN,mRN,mUN,mDN,mFN,mBN;
	//left point
	mLN.x=u_l_x;
	mLN.y=u_l_y;
	mLN.z=u_l_z;
	mLN.val=u_l;
		//right point
	mRN.x=u_r_x;
	mRN.y=u_r_y;
	mRN.z=u_r_z;
	mRN.val=u_r;
		//up point
	mUN.x=u_u_x;
	mUN.y=u_u_y;
	mUN.z=u_u_z;
	mUN.val=u_u;
		//down point
	mDN.x=u_d_x;
	mDN.y=u_d_y;
	mDN.z=u_d_z;
	mDN.val=u_d;
		//front point
	mFN.x=u_f_x;
	mFN.y=u_f_y;
	mFN.z=u_f_z;
	mFN.val=u_f;
		//back point
	mBN.x=u_b_x;
	mBN.y=u_b_y;
	mBN.z=u_b_z;
	mBN.val=u_b;
	//order the max med and min; generate a triagle
	miiCNode<float,int>mMaxN,mMedN,mMinN;
	FormTria(mLN,mRN,mUN,mDN,mFN,mBN,mMaxN,mMedN,mMinN);
	//calculate the three vectors and deltas;
	float fDelta[3]={0,0,0};
	float fCos[3]={0,0,0};
	CalVecDelta(mCN,mMaxN,mMedN,mMinN,fDelta,fCos);
     float ua=mMaxN.val;
	 float ub=mMedN.val;
	 float uc=mMinN.val;
	 float ca=fCos[0];
	 float cb=fCos[1];
	 float cc=fCos[2];
	 float ddela=fDelta[0];
	 float ddelb=fDelta[1];
	 float ddelc=fDelta[2];
	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);
	
}
double miiMinPathModel::Stencil3(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z-1;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z+1;
	int u_u_x = x,     u_u_y = y + 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y - 1, u_d_z = z;
	int u_f_x = x-1,     u_f_y = y,     u_f_z = z+1;
	int u_b_x = x+1,     u_b_y = y,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	//select the minimal value of each axis;
	float ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	float dualr[3]={0,0,0};//ua,dualr,dua are in one group
	if(u_l>u_r)
	{
		ua=u_r;
		dualr[0]=u_r_x;
		dualr[1]=u_r_y;
		dualr[2]=u_r_z;
	}
	else
	{
		ua=u_l;
		dualr[0]=u_l_x;
		dualr[1]=u_l_y;
		dualr[2]=u_l_z;
	}
	float ub = min(u_u, u_d);//y
	float dubud[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_u>u_d)
	{
		ub=u_d;
		dubud[0]=u_d_x;
		dubud[1]=u_d_y;
		dubud[2]=u_d_z;
	}
	else
	{
		ub=u_u;
		dubud[0]=u_u_x;
		dubud[1]=u_u_y;
		dubud[2]=u_u_z;
	}
	float uc = min(u_b, u_f);//z
	float ducbf[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_b>u_f)
	{
		uc=u_f;
		ducbf[0]=u_f_x;
		ducbf[1]=u_f_y;
		ducbf[2]=u_f_z;
	}
	else
	{
		ub=u_u;
		ducbf[0]=u_b_x;
		ducbf[1]=u_b_y;
		ducbf[2]=u_b_z;
	}
	//order the three minimal values and store their positions
	float ddela1=0;
	float ddelb1=0;
	float ddelc1=0;
	float ddelmax=0;
	float ddelmed=0;
	float ddelmin=0;
	float ddeltx=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[2],2));
	float ddelty=m_nImgSpacing[1];
	float ddeltz=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[2],2));
	float dMax, dMin, dMid;//the max u values
	float dmaxP[3]={0,0,0};
	float dmedP[3]={0,0,0};
	float dminP[3]={0,0,0};
	float daP[3]={0,0,0};
	float dbP[3]={0,0,0};
	float dcP[3]={0,0,0};
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dmaxP[0]=dualr[0];
		dmaxP[1]=dualr[1];
		dmaxP[2]=dualr[2];
		dMin = ub;
		ddelmin=ddelty;
		dminP[0]=dubud[0];
		dminP[1]=dubud[1];
		dminP[2]=dubud[2];

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dmaxP[0]=dubud[0];
		dmaxP[1]=dubud[1];
		dmaxP[2]=dubud[2];
		dMin = ua;
		ddelmin=ddeltx;
		dminP[0]=dualr[0];
		dminP[1]=dualr[1];
		dminP[2]=dualr[2];
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dmedP[0]=dmaxP[0];
		dmedP[1]=dmaxP[1];
		dmedP[2]=dmaxP[2];
		dMax = uc;
		ddelmax=ddeltz;
		dmaxP[0]=ducbf[0];
		dmaxP[1]=ducbf[1];
		dmaxP[2]=ducbf[2];
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
		dmedP[0]=ducbf[0];
		dmedP[1]=ducbf[1];
		dmedP[2]=ducbf[2];
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dmedP[0]=dminP[0];
		dmedP[1]=dminP[1];
		dmedP[2]=dminP[2];
		dMin = uc;
		ddelmin=ddeltz;
		dminP[0]=ducbf[0];
		dminP[1]=ducbf[1];
		dminP[2]=ducbf[2];
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	daP[0]=dmaxP[0];daP[1]=dmaxP[1];daP[2]=dmaxP[2];
	dbP[0]=dmedP[0];dbP[1]=dmedP[1];dbP[2]=dmedP[2];
	dcP[0]=dminP[0];dcP[1]=dminP[1];dcP[2]=dminP[2];
	ddela1=ddelmax;
	 ddelb1=ddelmed;
	 ddelc1=ddelmin;
	//
	float daPworld[3]={daP[0],daP[1],daP[2]};
	float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	float dxyzPworld[3]={x,y,z};
	
	m_pBaseImgInfo->ImageToWorld(daPworld);
	m_pBaseImgInfo->ImageToWorld(dbPworld);
	m_pBaseImgInfo->ImageToWorld(dcPworld);
	m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	// a b c as the vector
	float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	float ddela=zxh::VectorOP_Magnitude(dva,3);
	float ddelb=zxh::VectorOP_Magnitude(dvb,3);
	float ddelc=zxh::VectorOP_Magnitude(dvc,3);
	float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	float cc=zxh::VectorOP_Cosine(dva,dvc,3);
	//define the three unit vectors
	float dva1[3]={0,0,0};
	dva1[0]=(float)daP[0]-x;
	dva1[1]=(float)daP[1]-y;
	dva1[2]=(float)daP[2]-z;
	float dvb1[3]={0,0,0};
	dvb1[0]=(float)dbP[0]-x;
	dvb1[1]=(float)dbP[1]-y;
	dvb1[2]=(float)dbP[2]-z;
	float dvc1[3]={0,0,0};
	dvc1[0]=(float)dcP[0]-x;
	dvc1[1]=(float)dcP[1]-y;
	dvc1[2]=(float)dcP[2]-z;
	//normalize
	zxh::VectorOP_Normalise(dva1,3);
	zxh::VectorOP_Normalise(dvb1,3);
	zxh::VectorOP_Normalise(dvc1,3);
	//
	float ca1=zxh::VectorOP_Cosine(dva1,dvb1,3);
	float cb1=zxh::VectorOP_Cosine(dvb1,dvc1,3);
	float cc1=zxh::VectorOP_Cosine(dva1,dvc1,3);
	//

	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);

	//double dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	//double dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	//double dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	////solve the equation
	//double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	//bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	////solve 2 x
	//double dA2=(1-pow(cc,2))/pow(ddelb,2)+1/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc);
	//double dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+1/pow(ddelc,2)*pow(uc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	//double dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	//double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	//if (dRtMax_3 > ua && bSolv_3 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_3;
	//}
	//else if (dRtMax_2 > ub && bSolv_2 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
	////solve 3 x
	//double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//double dRoots_0[2] = {0}, dRoots_1[2] = {0};
	//
	//bool bSolv_0 = QuadraticRoots(da3, db3, dc3, dRoots_0);
	////solve 2 x
	//double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//bool bSolv_1 = QuadraticRoots(da2, db2, dc2, dRoots_1);

	//double dRtMax_0 = max(dRoots_0[0], dRoots_0[1]);
	//double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);

	//if (dRtMax_0 > ua && bSolv_0 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_0;
	//}
	//else if (dRtMax_1 > ub && bSolv_1 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
}
double miiMinPathModel::Stencil3_2(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	////select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z-1;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z+1;
	int u_u_x = x,     u_u_y = y + 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y - 1, u_d_z = z;
	int u_f_x = x-1,     u_f_y = y,     u_f_z = z+1;
	int u_b_x = x+1,     u_b_y = y,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	//if (u_l_x < 0 || u_l_x >= m_nImgWX)
	//	u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	//else
	//	u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	//if (u_r_x < 0 || u_r_x >= m_nImgWX)
	//	u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	//else
	//	u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	//if (u_u_y < 0 || u_u_y >= m_nImgWY)
	//	u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	//else
	//	u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	//if (u_d_y < 0 || u_d_y >= m_nImgWY)
	//	u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	//else
	//	u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	//if (u_b_z < 0 || u_b_z >= m_nImgWZ)
	//	u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	//else
	//	u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	//if (u_f_z < 0 || u_f_z >= m_nImgWZ)
	//	u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	//else
	//	u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];
	if (u_l_x < 0 || u_l_x >= m_nImgWX||u_l_z < 0 || u_l_z >= m_nImgWZ)
	{    
		if(u_l_x < 0)u_l_x=x;
		if(u_l_x >= m_nImgWX)u_l_x=x;
		if(u_l_z < 0)u_l_z=z;
		if(u_l_z >= m_nImgWZ)u_l_z=z;
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX||u_r_z < 0 || u_r_z >= m_nImgWZ)
	{
		if(u_r_x < 0)u_r_x=x;
		if(u_r_x >= m_nImgWX)u_r_x=x;
		if(u_r_z < 0)u_r_z=z;
		if(u_r_z >= m_nImgWZ)u_r_z=z;
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
	{
		if(u_u_y < 0)u_u_y=y;
		if(u_u_y >= m_nImgWY)u_u_y=y;
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];
	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		{
		if(u_d_y < 0)u_d_y=y;
		if(u_d_y >= m_nImgWY)u_d_y=y;
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ||u_b_x < 0 || u_b_x >= m_nImgWX)
	{
		if(u_b_z < 0)u_b_z=z;
		if(u_b_z >= m_nImgWZ)u_b_z=z;
		if(u_b_x < 0)u_b_x=x;
		if(u_b_x >= m_nImgWX)u_b_x=x;
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ||u_f_x < 0 || u_f_x >= m_nImgWX)
	{
		if(u_f_z < 0)u_f_z=z;
		if(u_f_z >= m_nImgWZ)u_f_z=z;
		if(u_f_x < 0)u_f_x=x;
		if(u_f_x >= m_nImgWX)u_f_x=x;
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];
	//store the center point
	miiCNode<float,int> mCN;
	mCN.x=x;
	mCN.y=y;
	mCN.z=z;
	if(x==188&&y==157&&z==188)
		int ms=0;
	//store the six points with u
	miiCNode<float,int> mLN,mRN,mUN,mDN,mFN,mBN;
	//left point
	mLN.x=u_l_x;
	mLN.y=u_l_y;
	mLN.z=u_l_z;
	mLN.val=u_l;
		//right point
	mRN.x=u_r_x;
	mRN.y=u_r_y;
	mRN.z=u_r_z;
	mRN.val=u_r;
		//up point
	mUN.x=u_u_x;
	mUN.y=u_u_y;
	mUN.z=u_u_z;
	mUN.val=u_u;
		//down point
	mDN.x=u_d_x;
	mDN.y=u_d_y;
	mDN.z=u_d_z;
	mDN.val=u_d;
		//front point
	mFN.x=u_f_x;
	mFN.y=u_f_y;
	mFN.z=u_f_z;
	mFN.val=u_f;
		//back point
	mBN.x=u_b_x;
	mBN.y=u_b_y;
	mBN.z=u_b_z;
	mBN.val=u_b;
	//order the max med and min; generate a triagle
	miiCNode<float,int>mMaxN,mMedN,mMinN;
	FormTria(mLN,mRN,mUN,mDN,mFN,mBN,mMaxN,mMedN,mMinN);
	//calculate the three vectors and deltas;
	float fDelta[3]={0,0,0};
	float fCos[3]={0,0,0};
	CalVecDelta(mCN,mMaxN,mMedN,mMinN,fDelta,fCos);
     float ua=mMaxN.val;
	 float ub=mMedN.val;
	 float uc=mMinN.val;
	 float ca=fCos[0];
	 float cb=fCos[1];
	 float cc=fCos[2];
	 float ddela=fDelta[0];
	 float ddelb=fDelta[1];
	 float ddelc=fDelta[2];
	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);
	
}
double miiMinPathModel::Stencil4(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y+1,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y-1,	  u_r_z = z;
	int u_u_x = x+1,     u_u_y = y + 1, u_u_z = z;
	int u_d_x = x-1,     u_d_y = y - 1, u_d_z = z;
	int u_f_x = x,     u_f_y = y,     u_f_z = z+1;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	//select the minimal value of each axis;
	float ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	float dualr[3]={0,0,0};//ua,dualr,dua are in one group
	if(u_l>u_r)
	{
		ua=u_r;
		dualr[0]=u_r_x;
		dualr[1]=u_r_y;
		dualr[2]=u_r_z;
	}
	else
	{
		ua=u_l;
		dualr[0]=u_l_x;
		dualr[1]=u_l_y;
		dualr[2]=u_l_z;
	}
	float ub = min(u_u, u_d);//y
	float dubud[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_u>u_d)
	{
		ub=u_d;
		dubud[0]=u_d_x;
		dubud[1]=u_d_y;
		dubud[2]=u_d_z;
	}
	else
	{
		ub=u_u;
		dubud[0]=u_u_x;
		dubud[1]=u_u_y;
		dubud[2]=u_u_z;
	}
	float uc = min(u_b, u_f);//z
	float ducbf[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_b>u_f)
	{
		uc=u_f;
		ducbf[0]=u_f_x;
		ducbf[1]=u_f_y;
		ducbf[2]=u_f_z;
	}
	else
	{
		ub=u_u;
		ducbf[0]=u_b_x;
		ducbf[1]=u_b_y;
		ducbf[2]=u_b_z;
	}
	//order the three minimal values and store their positions
	float ddela1=0;
	float ddelb1=0;
	float ddelc1=0;
	float ddelmax=0;
	float ddelmed=0;
	float ddelmin=0;
	float ddeltx=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[1],2));
	float ddelty=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[1],2));
	float ddeltz=m_nImgSpacing[2];
	float dMax, dMin, dMid;//the max u values
	float dmaxP[3]={0,0,0};
	float dmedP[3]={0,0,0};
	float dminP[3]={0,0,0};
	float daP[3]={0,0,0};
	float dbP[3]={0,0,0};
	float dcP[3]={0,0,0};
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dmaxP[0]=dualr[0];
		dmaxP[1]=dualr[1];
		dmaxP[2]=dualr[2];
		dMin = ub;
		ddelmin=ddelty;
		dminP[0]=dubud[0];
		dminP[1]=dubud[1];
		dminP[2]=dubud[2];

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dmaxP[0]=dubud[0];
		dmaxP[1]=dubud[1];
		dmaxP[2]=dubud[2];
		dMin = ua;
		ddelmin=ddeltx;
		dminP[0]=dualr[0];
		dminP[1]=dualr[1];
		dminP[2]=dualr[2];
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dmedP[0]=dmaxP[0];
		dmedP[1]=dmaxP[1];
		dmedP[2]=dmaxP[2];
		dMax = uc;
		ddelmax=ddeltz;
		dmaxP[0]=ducbf[0];
		dmaxP[1]=ducbf[1];
		dmaxP[2]=ducbf[2];
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
		dmedP[0]=ducbf[0];
		dmedP[1]=ducbf[1];
		dmedP[2]=ducbf[2];
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dmedP[0]=dminP[0];
		dmedP[1]=dminP[1];
		dmedP[2]=dminP[2];
		dMin = uc;
		ddelmin=ddeltz;
		dminP[0]=ducbf[0];
		dminP[1]=ducbf[1];
		dminP[2]=ducbf[2];
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	daP[0]=dmaxP[0];daP[1]=dmaxP[1];daP[2]=dmaxP[2];
	dbP[0]=dmedP[0];dbP[1]=dmedP[1];dbP[2]=dmedP[2];
	dcP[0]=dminP[0];dcP[1]=dminP[1];dcP[2]=dminP[2];
	ddela1=ddelmax;
	 ddelb1=ddelmed;
	 ddelc1=ddelmin;
	//
	float daPworld[3]={daP[0],daP[1],daP[2]};
	float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	float dxyzPworld[3]={x,y,z};
	
	m_pBaseImgInfo->ImageToWorld(daPworld);
	m_pBaseImgInfo->ImageToWorld(dbPworld);
	m_pBaseImgInfo->ImageToWorld(dcPworld);
	m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	// a b c as the vector
	float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	float ddela=zxh::VectorOP_Magnitude(dva,3);
	float ddelb=zxh::VectorOP_Magnitude(dvb,3);
	float ddelc=zxh::VectorOP_Magnitude(dvc,3);
	float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	float cc=zxh::VectorOP_Cosine(dva,dvc,3);
	//define the three unit vectors
	float dva1[3]={0,0,0};
	dva1[0]=(float)daP[0]-x;
	dva1[1]=(float)daP[1]-y;
	dva1[2]=(float)daP[2]-z;
	float dvb1[3]={0,0,0};
	dvb1[0]=(float)dbP[0]-x;
	dvb1[1]=(float)dbP[1]-y;
	dvb1[2]=(float)dbP[2]-z;
	float dvc1[3]={0,0,0};
	dvc1[0]=(float)dcP[0]-x;
	dvc1[1]=(float)dcP[1]-y;
	dvc1[2]=(float)dcP[2]-z;
	//normalize
	zxh::VectorOP_Normalise(dva1,3);
	zxh::VectorOP_Normalise(dvb1,3);
	zxh::VectorOP_Normalise(dvc1,3);
	//
	float ca1=zxh::VectorOP_Cosine(dva1,dvb1,3);
	float cb1=zxh::VectorOP_Cosine(dvb1,dvc1,3);
	float cc1=zxh::VectorOP_Cosine(dva1,dvc1,3);
	//

	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);

	//double dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	//double dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	//double dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	////solve the equation
	//double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	//bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	////solve 2 x
	//double dA2=(1-pow(cc,2))/pow(ddelb,2)+1/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc);
	//double dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+1/pow(ddelc,2)*pow(uc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	//double dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	//double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	//if (dRtMax_3 > ua && bSolv_3 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_3;
	//}
	//else if (dRtMax_2 > ub && bSolv_2 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
	////solve 3 x
	//double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//double dRoots_0[2] = {0}, dRoots_1[2] = {0};
	//
	//bool bSolv_0 = QuadraticRoots(da3, db3, dc3, dRoots_0);
	////solve 2 x
	//double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//bool bSolv_1 = QuadraticRoots(da2, db2, dc2, dRoots_1);

	//double dRtMax_0 = max(dRoots_0[0], dRoots_0[1]);
	//double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);

	//if (dRtMax_0 > ua && bSolv_0 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_0;
	//}
	//else if (dRtMax_1 > ub && bSolv_1 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
}
double miiMinPathModel::Stencil4_2(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y+1,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y-1,	  u_r_z = z;
	int u_u_x = x+1,     u_u_y = y + 1, u_u_z = z;
	int u_d_x = x-1,     u_d_y = y - 1, u_d_z = z;
	int u_f_x = x,     u_f_y = y,     u_f_z = z+1;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;
	/*
	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];*/
		if (u_l_x < 0 || u_l_x >= m_nImgWX||u_l_y < 0 || u_l_y >= m_nImgWY)
	{    
		if(u_l_x < 0)u_l_x=x;
		if(u_l_x >= m_nImgWX)u_l_x=x;
		if(u_l_y < 0)u_l_y=y;
		if(u_l_y >= m_nImgWY)u_l_y=y;
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX||u_r_y < 0 || u_r_y >= m_nImgWY)
	{
		if(u_r_x < 0)u_r_x=x;
		if(u_r_x >= m_nImgWX)u_r_x=x;
		if(u_r_y < 0)u_r_y=y;
		if(u_r_y >= m_nImgWY)u_r_y=y;
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY||u_u_x < 0 || u_u_x >= m_nImgWX)
		{
		if(u_u_y < 0)u_u_y=y;
		if(u_u_y >= m_nImgWY)u_u_y=y;
		if(u_u_x < 0)u_u_x=x;
		if(u_u_x >= m_nImgWX)u_u_x=x;
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY||u_d_x < 0 || u_d_x >= m_nImgWX)
		{
		if(u_d_y < 0)u_d_y=y;
		if(u_d_y >= m_nImgWY)u_d_y=y;
		if(u_d_x < 0)u_d_x=x;
		if(u_d_x >= m_nImgWX)u_d_x=x;
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		{
		if(u_b_z < 0)u_b_z=z;
		if(u_b_z >= m_nImgWZ)u_b_z=z;
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		{
		if(u_f_z < 0)u_f_z=z;
		if(u_f_z >= m_nImgWZ)u_f_z=z;
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];
	//store the center point
	miiCNode<float,int> mCN;
	mCN.x=x;
	mCN.y=y;
	mCN.z=z;
	if(x==188&&y==157&&z==188)
		int ms=0;
	//store the six points with u
	miiCNode<float,int> mLN,mRN,mUN,mDN,mFN,mBN;
	//left point
	mLN.x=u_l_x;
	mLN.y=u_l_y;
	mLN.z=u_l_z;
	mLN.val=u_l;
		//right point
	mRN.x=u_r_x;
	mRN.y=u_r_y;
	mRN.z=u_r_z;
	mRN.val=u_r;
		//up point
	mUN.x=u_u_x;
	mUN.y=u_u_y;
	mUN.z=u_u_z;
	mUN.val=u_u;
		//down point
	mDN.x=u_d_x;
	mDN.y=u_d_y;
	mDN.z=u_d_z;
	mDN.val=u_d;
		//front point
	mFN.x=u_f_x;
	mFN.y=u_f_y;
	mFN.z=u_f_z;
	mFN.val=u_f;
		//back point
	mBN.x=u_b_x;
	mBN.y=u_b_y;
	mBN.z=u_b_z;
	mBN.val=u_b;
	//order the max med and min; generate a triagle
	miiCNode<float,int>mMaxN,mMedN,mMinN;
	FormTria(mLN,mRN,mUN,mDN,mFN,mBN,mMaxN,mMedN,mMinN);
	//calculate the three vectors and deltas;
	float fDelta[3]={0,0,0};
	float fCos[3]={0,0,0};
	CalVecDelta(mCN,mMaxN,mMedN,mMinN,fDelta,fCos);
     float ua=mMaxN.val;
	 float ub=mMedN.val;
	 float uc=mMinN.val;
	 float ca=fCos[0];
	 float cb=fCos[1];
	 float cc=fCos[2];
	 float ddela=fDelta[0];
	 float ddelb=fDelta[1];
	 float ddelc=fDelta[2];
	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);
	
}
double miiMinPathModel::Stencil5(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y+1,	  u_l_z = z+1;  
	int u_r_x = x + 1, u_r_y = y-1,	  u_r_z = z-1;
	int u_u_x = x+1,     u_u_y = y + 1, u_u_z = z-1;
	int u_d_x = x-1,     u_d_y = y - 1, u_d_z = z+1;
	int u_f_x = x+1,     u_f_y = y,     u_f_z = z+1;
	int u_b_x = x-1,     u_b_y = y,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	//select the minimal value of each axis;
	float ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	float dualr[3]={0,0,0};//ua,dualr,dua are in one group
	if(u_l>u_r)
	{
		ua=u_r;
		dualr[0]=u_r_x;
		dualr[1]=u_r_y;
		dualr[2]=u_r_z;
	}
	else
	{
		ua=u_l;
		dualr[0]=u_l_x;
		dualr[1]=u_l_y;
		dualr[2]=u_l_z;
	}
	float ub = min(u_u, u_d);//y
	float dubud[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_u>u_d)
	{
		ub=u_d;
		dubud[0]=u_d_x;
		dubud[1]=u_d_y;
		dubud[2]=u_d_z;
	}
	else
	{
		ub=u_u;
		dubud[0]=u_u_x;
		dubud[1]=u_u_y;
		dubud[2]=u_u_z;
	}
	float uc = min(u_b, u_f);//z
	float ducbf[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_b>u_f)
	{
		uc=u_f;
		ducbf[0]=u_f_x;
		ducbf[1]=u_f_y;
		ducbf[2]=u_f_z;
	}
	else
	{
		ub=u_u;
		ducbf[0]=u_b_x;
		ducbf[1]=u_b_y;
		ducbf[2]=u_b_z;
	}
	//order the three minimal values and store their positions
	float ddela1=0;
	float ddelb1=0;
	float ddelc1=0;
	float ddelmax=0;
	float ddelmed=0;
	float ddelmin=0;
	float ddeltx=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[1],2));
	float ddelty=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[1],2));
	float ddeltz=m_nImgSpacing[2];
	float dMax, dMin, dMid;//the max u values
	float dmaxP[3]={0,0,0};
	float dmedP[3]={0,0,0};
	float dminP[3]={0,0,0};
	float daP[3]={0,0,0};
	float dbP[3]={0,0,0};
	float dcP[3]={0,0,0};
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dmaxP[0]=dualr[0];
		dmaxP[1]=dualr[1];
		dmaxP[2]=dualr[2];
		dMin = ub;
		ddelmin=ddelty;
		dminP[0]=dubud[0];
		dminP[1]=dubud[1];
		dminP[2]=dubud[2];

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dmaxP[0]=dubud[0];
		dmaxP[1]=dubud[1];
		dmaxP[2]=dubud[2];
		dMin = ua;
		ddelmin=ddeltx;
		dminP[0]=dualr[0];
		dminP[1]=dualr[1];
		dminP[2]=dualr[2];
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dmedP[0]=dmaxP[0];
		dmedP[1]=dmaxP[1];
		dmedP[2]=dmaxP[2];
		dMax = uc;
		ddelmax=ddeltz;
		dmaxP[0]=ducbf[0];
		dmaxP[1]=ducbf[1];
		dmaxP[2]=ducbf[2];
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
		dmedP[0]=ducbf[0];
		dmedP[1]=ducbf[1];
		dmedP[2]=ducbf[2];
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dmedP[0]=dminP[0];
		dmedP[1]=dminP[1];
		dmedP[2]=dminP[2];
		dMin = uc;
		ddelmin=ddeltz;
		dminP[0]=ducbf[0];
		dminP[1]=ducbf[1];
		dminP[2]=ducbf[2];
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	daP[0]=dmaxP[0];daP[1]=dmaxP[1];daP[2]=dmaxP[2];
	dbP[0]=dmedP[0];dbP[1]=dmedP[1];dbP[2]=dmedP[2];
	dcP[0]=dminP[0];dcP[1]=dminP[1];dcP[2]=dminP[2];
	ddela1=ddelmax;
	 ddelb1=ddelmed;
	 ddelc1=ddelmin;
	//
	float daPworld[3]={daP[0],daP[1],daP[2]};
	float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	float dxyzPworld[3]={x,y,z};
	
	m_pBaseImgInfo->ImageToWorld(daPworld);
	m_pBaseImgInfo->ImageToWorld(dbPworld);
	m_pBaseImgInfo->ImageToWorld(dcPworld);
	m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	// a b c as the vector
	float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	float ddela=zxh::VectorOP_Magnitude(dva,3);
	float ddelb=zxh::VectorOP_Magnitude(dvb,3);
	float ddelc=zxh::VectorOP_Magnitude(dvc,3);
	float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	float cc=zxh::VectorOP_Cosine(dva,dvc,3);
	//define the three unit vectors
	float dva1[3]={0,0,0};
	dva1[0]=(float)daP[0]-x;
	dva1[1]=(float)daP[1]-y;
	dva1[2]=(float)daP[2]-z;
	float dvb1[3]={0,0,0};
	dvb1[0]=(float)dbP[0]-x;
	dvb1[1]=(float)dbP[1]-y;
	dvb1[2]=(float)dbP[2]-z;
	float dvc1[3]={0,0,0};
	dvc1[0]=(float)dcP[0]-x;
	dvc1[1]=(float)dcP[1]-y;
	dvc1[2]=(float)dcP[2]-z;
	//normalize
	zxh::VectorOP_Normalise(dva1,3);
	zxh::VectorOP_Normalise(dvb1,3);
	zxh::VectorOP_Normalise(dvc1,3);
	//
	float ca1=zxh::VectorOP_Cosine(dva1,dvb1,3);
	float cb1=zxh::VectorOP_Cosine(dvb1,dvc1,3);
	float cc1=zxh::VectorOP_Cosine(dva1,dvc1,3);
	//

	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);

	//double dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	//double dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	//double dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	////solve the equation
	//double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	//bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	////solve 2 x
	//double dA2=(1-pow(cc,2))/pow(ddelb,2)+1/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc);
	//double dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+1/pow(ddelc,2)*pow(uc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	//double dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	//double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	//if (dRtMax_3 > ua && bSolv_3 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_3;
	//}
	//else if (dRtMax_2 > ub && bSolv_2 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
	////solve 3 x
	//double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//double dRoots_0[2] = {0}, dRoots_1[2] = {0};
	//
	//bool bSolv_0 = QuadraticRoots(da3, db3, dc3, dRoots_0);
	////solve 2 x
	//double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//bool bSolv_1 = QuadraticRoots(da2, db2, dc2, dRoots_1);

	//double dRtMax_0 = max(dRoots_0[0], dRoots_0[1]);
	//double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);

	//if (dRtMax_0 > ua && bSolv_0 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_0;
	//}
	//else if (dRtMax_1 > ub && bSolv_1 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
}
double miiMinPathModel::Stencil5_2(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y+1,	  u_l_z = z+1;  
	int u_r_x = x + 1, u_r_y = y-1,	  u_r_z = z-1;
	int u_u_x = x+1,     u_u_y = y + 1, u_u_z = z-1;
	int u_d_x = x-1,     u_d_y = y - 1, u_d_z = z+1;
	int u_f_x = x+1,     u_f_y = y,     u_f_z = z+1;
	int u_b_x = x-1,     u_b_y = y,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;
	/*
	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];*/
	if (u_l_x < 0 || u_l_x >= m_nImgWX)
	{    
		if(u_l_x < 0)u_l_x=x;
		if(u_l_x >= m_nImgWX)u_l_x=x;
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
	{
		if(u_r_x < 0)u_r_x=x;
		if(u_r_x >= m_nImgWX)u_r_x=x;
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
	{
		if(u_u_y < 0)u_u_y=y;
		if(u_u_y >= m_nImgWY)u_u_y=y;
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
	{
		if(u_d_y < 0)u_d_y=y;
		if(u_d_y >= m_nImgWY)u_d_y=y;
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
	{
		if(u_b_z < 0)u_b_z=z;
		if(u_b_z >= m_nImgWZ)u_b_z=z;
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
	{
		if(u_f_z < 0)u_f_z=z;
		if(u_f_z >= m_nImgWZ)u_f_z=z;
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];
	//store the center point
	miiCNode<float,int> mCN;
	mCN.x=x;
	mCN.y=y;
	mCN.z=z;
	if(x==188&&y==157&&z==188)
		int ms=0;
	//store the six points with u
	miiCNode<float,int> mLN,mRN,mUN,mDN,mFN,mBN;
	//left point
	mLN.x=u_l_x;
	mLN.y=u_l_y;
	mLN.z=u_l_z;
	mLN.val=u_l;
		//right point
	mRN.x=u_r_x;
	mRN.y=u_r_y;
	mRN.z=u_r_z;
	mRN.val=u_r;
		//up point
	mUN.x=u_u_x;
	mUN.y=u_u_y;
	mUN.z=u_u_z;
	mUN.val=u_u;
		//down point
	mDN.x=u_d_x;
	mDN.y=u_d_y;
	mDN.z=u_d_z;
	mDN.val=u_d;
		//front point
	mFN.x=u_f_x;
	mFN.y=u_f_y;
	mFN.z=u_f_z;
	mFN.val=u_f;
		//back point
	mBN.x=u_b_x;
	mBN.y=u_b_y;
	mBN.z=u_b_z;
	mBN.val=u_b;
	//order the max med and min; generate a triagle
	miiCNode<float,int>mMaxN,mMedN,mMinN;
	FormTria(mLN,mRN,mUN,mDN,mFN,mBN,mMaxN,mMedN,mMinN);
	//calculate the three vectors and deltas;
	float fDelta[3]={0,0,0};
	float fCos[3]={0,0,0};
	CalVecDelta(mCN,mMaxN,mMedN,mMinN,fDelta,fCos);
     float ua=mMaxN.val;
	 float ub=mMedN.val;
	 float uc=mMinN.val;
	 float ca=fCos[0];
	 float cb=fCos[1];
	 float cc=fCos[2];
	 float ddela=fDelta[0];
	 float ddelb=fDelta[1];
	 float ddelc=fDelta[2];
	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);
	
}
double miiMinPathModel::Stencil6(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z+1;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z-1;
	int u_u_x = x-1,     u_u_y = y + 1, u_u_z = z-1;
	int u_d_x = x+1,     u_d_y = y - 1, u_d_z = z+1;
	int u_f_x = x+1,     u_f_y = y+1,     u_f_z = z+1;
	int u_b_x = x-1,     u_b_y = y-1,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	//select the minimal value of each axis;
	float ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	float dualr[3]={0,0,0};//ua,dualr,dua are in one group
	if(u_l>u_r)
	{
		ua=u_r;
		dualr[0]=u_r_x;
		dualr[1]=u_r_y;
		dualr[2]=u_r_z;
	}
	else
	{
		ua=u_l;
		dualr[0]=u_l_x;
		dualr[1]=u_l_y;
		dualr[2]=u_l_z;
	}
	float ub = min(u_u, u_d);//y
	float dubud[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_u>u_d)
	{
		ub=u_d;
		dubud[0]=u_d_x;
		dubud[1]=u_d_y;
		dubud[2]=u_d_z;
	}
	else
	{
		ub=u_u;
		dubud[0]=u_u_x;
		dubud[1]=u_u_y;
		dubud[2]=u_u_z;
	}
	float uc = min(u_b, u_f);//z
	float ducbf[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_b>u_f)
	{
		uc=u_f;
		ducbf[0]=u_f_x;
		ducbf[1]=u_f_y;
		ducbf[2]=u_f_z;
	}
	else
	{
		ub=u_u;
		ducbf[0]=u_b_x;
		ducbf[1]=u_b_y;
		ducbf[2]=u_b_z;
	}
	//order the three minimal values and store their positions
	float ddela1=0;
	float ddelb1=0;
	float ddelc1=0;
	float ddelmax=0;
	float ddelmed=0;
	float ddelmin=0;
	float ddeltx=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[1],2));
	float ddelty=sqrt(pow(m_nImgSpacing[0],2)+pow(m_nImgSpacing[1],2));
	float ddeltz=m_nImgSpacing[2];
	float dMax, dMin, dMid;//the max u values
	float dmaxP[3]={0,0,0};
	float dmedP[3]={0,0,0};
	float dminP[3]={0,0,0};
	float daP[3]={0,0,0};
	float dbP[3]={0,0,0};
	float dcP[3]={0,0,0};
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dmaxP[0]=dualr[0];
		dmaxP[1]=dualr[1];
		dmaxP[2]=dualr[2];
		dMin = ub;
		ddelmin=ddelty;
		dminP[0]=dubud[0];
		dminP[1]=dubud[1];
		dminP[2]=dubud[2];

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dmaxP[0]=dubud[0];
		dmaxP[1]=dubud[1];
		dmaxP[2]=dubud[2];
		dMin = ua;
		ddelmin=ddeltx;
		dminP[0]=dualr[0];
		dminP[1]=dualr[1];
		dminP[2]=dualr[2];
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dmedP[0]=dmaxP[0];
		dmedP[1]=dmaxP[1];
		dmedP[2]=dmaxP[2];
		dMax = uc;
		ddelmax=ddeltz;
		dmaxP[0]=ducbf[0];
		dmaxP[1]=ducbf[1];
		dmaxP[2]=ducbf[2];
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
		dmedP[0]=ducbf[0];
		dmedP[1]=ducbf[1];
		dmedP[2]=ducbf[2];
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dmedP[0]=dminP[0];
		dmedP[1]=dminP[1];
		dmedP[2]=dminP[2];
		dMin = uc;
		ddelmin=ddeltz;
		dminP[0]=ducbf[0];
		dminP[1]=ducbf[1];
		dminP[2]=ducbf[2];
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	daP[0]=dmaxP[0];daP[1]=dmaxP[1];daP[2]=dmaxP[2];
	dbP[0]=dmedP[0];dbP[1]=dmedP[1];dbP[2]=dmedP[2];
	dcP[0]=dminP[0];dcP[1]=dminP[1];dcP[2]=dminP[2];
	ddela1=ddelmax;
	 ddelb1=ddelmed;
	 ddelc1=ddelmin;
	//
	float daPworld[3]={daP[0],daP[1],daP[2]};
	float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	float dxyzPworld[3]={x,y,z};
	
	m_pBaseImgInfo->ImageToWorld(daPworld);
	m_pBaseImgInfo->ImageToWorld(dbPworld);
	m_pBaseImgInfo->ImageToWorld(dcPworld);
	m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	// a b c as the vector
	float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	float ddela=zxh::VectorOP_Magnitude(dva,3);
	float ddelb=zxh::VectorOP_Magnitude(dvb,3);
	float ddelc=zxh::VectorOP_Magnitude(dvc,3);
	float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	float cc=zxh::VectorOP_Cosine(dva,dvc,3);
	//define the three unit vectors
	float dva1[3]={0,0,0};
	dva1[0]=(float)daP[0]-x;
	dva1[1]=(float)daP[1]-y;
	dva1[2]=(float)daP[2]-z;
	float dvb1[3]={0,0,0};
	dvb1[0]=(float)dbP[0]-x;
	dvb1[1]=(float)dbP[1]-y;
	dvb1[2]=(float)dbP[2]-z;
	float dvc1[3]={0,0,0};
	dvc1[0]=(float)dcP[0]-x;
	dvc1[1]=(float)dcP[1]-y;
	dvc1[2]=(float)dcP[2]-z;
	//normalize
	zxh::VectorOP_Normalise(dva1,3);
	zxh::VectorOP_Normalise(dvb1,3);
	zxh::VectorOP_Normalise(dvc1,3);
	//
	float ca1=zxh::VectorOP_Cosine(dva1,dvb1,3);
	float cb1=zxh::VectorOP_Cosine(dvb1,dvc1,3);
	float cc1=zxh::VectorOP_Cosine(dva1,dvc1,3);
	//

	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);

	//double dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	//double dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	//double dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
	//	2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	////solve the equation
	//double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	//bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	////solve 2 x
	//double dA2=(1-pow(cc,2))/pow(ddelb,2)+1/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc);
	//double dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc/pow(ddelc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*(-ub-uc);
	//double dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+1/pow(ddelc,2)*pow(uc,2)+\
	//	2*(-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	//double dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	//double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	//if (dRtMax_3 > ua && bSolv_3 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_3;
	//}
	//else if (dRtMax_2 > ub && bSolv_2 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
	////solve 3 x
	//double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//double dRoots_0[2] = {0}, dRoots_1[2] = {0};
	//
	//bool bSolv_0 = QuadraticRoots(da3, db3, dc3, dRoots_0);
	////solve 2 x
	//double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//bool bSolv_1 = QuadraticRoots(da2, db2, dc2, dRoots_1);

	//double dRtMax_0 = max(dRoots_0[0], dRoots_0[1]);
	//double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);

	//if (dRtMax_0 > ua && bSolv_0 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_0;
	//}
	//else if (dRtMax_1 > ub && bSolv_1 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
}
double miiMinPathModel::Stencil6_2(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

		//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z+1;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z-1;
	int u_u_x = x-1,     u_u_y = y + 1, u_u_z = z-1;
	int u_d_x = x+1,     u_d_y = y - 1, u_d_z = z+1;
	int u_f_x = x+1,     u_f_y = y+1,     u_f_z = z+1;
	int u_b_x = x-1,     u_b_y = y-1,     u_b_z = z - 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;
	/*
	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];*/
		if (u_l_x < 0 || u_l_x >= m_nImgWX)
	{    
		if(u_l_x < 0)u_l_x=x;
		if(u_l_x >= m_nImgWX)u_l_x=x;
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
	{
		if(u_r_x < 0)u_r_x=x;
		if(u_r_x >= m_nImgWX)u_r_x=x;
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		{
		if(u_u_y < 0)u_u_y=y;
		if(u_u_y >= m_nImgWY)u_u_y=y;
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		{
		if(u_d_y < 0)u_d_y=y;
		if(u_d_y >= m_nImgWY)u_d_y=y;
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		{
		if(u_b_z < 0)u_b_z=z;
		if(u_b_z >= m_nImgWZ)u_b_z=z;
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		{
		if(u_f_z < 0)u_f_z=z;
		if(u_f_z >= m_nImgWZ)u_f_z=z;
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	}
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];
	//store the center point
	miiCNode<float,int> mCN;
	mCN.x=x;
	mCN.y=y;
	mCN.z=z;
	if(x==188&&y==157&&z==188)
		int ms=0;
	//store the six points with u
	miiCNode<float,int> mLN,mRN,mUN,mDN,mFN,mBN;
	//left point
	mLN.x=u_l_x;
	mLN.y=u_l_y;
	mLN.z=u_l_z;
	mLN.val=u_l;
		//right point
	mRN.x=u_r_x;
	mRN.y=u_r_y;
	mRN.z=u_r_z;
	mRN.val=u_r;
		//up point
	mUN.x=u_u_x;
	mUN.y=u_u_y;
	mUN.z=u_u_z;
	mUN.val=u_u;
		//down point
	mDN.x=u_d_x;
	mDN.y=u_d_y;
	mDN.z=u_d_z;
	mDN.val=u_d;
		//front point
	mFN.x=u_f_x;
	mFN.y=u_f_y;
	mFN.z=u_f_z;
	mFN.val=u_f;
		//back point
	mBN.x=u_b_x;
	mBN.y=u_b_y;
	mBN.z=u_b_z;
	mBN.val=u_b;
	//order the max med and min; generate a triagle
	miiCNode<float,int>mMaxN,mMedN,mMinN;
	FormTria(mLN,mRN,mUN,mDN,mFN,mBN,mMaxN,mMedN,mMinN);
	//calculate the three vectors and deltas;
	float fDelta[3]={0,0,0};
	float fCos[3]={0,0,0};
	CalVecDelta(mCN,mMaxN,mMedN,mMinN,fDelta,fCos);
     float ua=mMaxN.val;
	 float ub=mMedN.val;
	 float uc=mMinN.val;
	 float ca=fCos[0];
	 float cb=fCos[1];
	 float cc=fCos[2];
	 float ddela=fDelta[0];
	 float ddelb=fDelta[1];
	 float ddelc=fDelta[2];
	return SolveEq(ua,ub,uc,ca,cb,cc,ddela,ddelb,ddelc,dP2);
	
}
	double miiMinPathModel::SolveEq(float ua,float ub,float uc,float ca,float cb,float cc,float ddela,float ddelb,float ddelc,float dP2)
	{
		
	float dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	float dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
		2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	float dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
		2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	float dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
		2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//solve the equation
	double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	//solve 2 x
	float dA2=(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
		2*(ca*cc-cb)/(ddelb*ddelc);
	float dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
		2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	float dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
		2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	float dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	float dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	float droot=0;
	if (dRtMax_3 > ua && bSolv_3 == true)
	{
		droot= dRtMax_3;
	}
	else if (dRtMax_2 > ub && bSolv_2 == true)
	{
		droot= dRtMax_2;
	} 
	else
	{
		droot = sqrt(dP2*dD1/(1-pow(ca,2)))*ddelc + uc;
	}
	return droot;
	}
void miiMinPathModel::Stencil1(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points

	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	float ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	float ub = min(u_u, u_d);//y
	float uc = min(u_b, u_f);//z
	float ddela=0;
	float ddelb=0;
	float ddelc=0;
	float ddelmax=0;
	float ddelmed=0;
	float ddelmin=0;

	float ddeltx=m_nImgSpacing[0];
	float ddelty=m_nImgSpacing[1];
	float ddeltz=m_nImgSpacing[2];

	float dMax, dMin, dMid;
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dMin = ub;
		ddelmin=ddelty;

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dMin = ua;
		ddelmin=ddeltx;
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dMax = uc;
		ddelmax=ddeltz;
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dMin = uc;
		ddelmin=ddeltz;
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	ddela=ddelmax; ddelb=ddelmed;ddelc=ddelmin;
	//solve 3 x
	double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(da3, db3, dc3, dRoots_1);
	//solve 2 x
	double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	bool bSolv_2 = QuadraticRoots(da2, db2, dc2, dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > ua && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > ub && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(dP2) + uc;
	}
}
void miiMinPathModel::Stencil1_1(int x, int y, int z,float dP2)
{
	// stencil 1 including four neighbour points
	double hx=m_nImgSpacing[0];
	double hy=m_nImgSpacing[1];
	double hz=m_nImgSpacing[2];

	//select the u and order the u values decreasing
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y + 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y - 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	float u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (float)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (float)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (float)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (float)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (float)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (float)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (float)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	//select the minimal value of each axis;
	float ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	float dualr[3]={0,0,0};//ua,dualr,dua are in one group
	if(u_l>u_r)
	{
		ua=u_r;
		dualr[0]=u_r_x;
		dualr[1]=u_r_y;
		dualr[2]=u_r_z;
	}
	else
	{
		ua=u_l;
		dualr[0]=u_l_x;
		dualr[1]=u_l_y;
		dualr[2]=u_l_z;
	}
	float ub = min(u_u, u_d);//y
	float dubud[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_u>u_d)
	{
		ub=u_d;
		dubud[0]=u_d_x;
		dubud[1]=u_d_y;
		dubud[2]=u_d_z;
	}
	else
	{
		ub=u_u;
		dubud[0]=u_u_x;
		dubud[1]=u_u_y;
		dubud[2]=u_u_z;
	}
	float uc = min(u_b, u_f);//z
	float ducbf[3]={0,0,0};//ub,dubud,dub are in one group
	if(u_b>u_f)
	{
		uc=u_f;
		ducbf[0]=u_f_x;
		ducbf[1]=u_f_y;
		ducbf[2]=u_f_z;
	}
	else
	{
		ub=u_u;
		ducbf[0]=u_b_x;
		ducbf[1]=u_b_y;
		ducbf[2]=u_b_z;
	}
	//order the three minimal values and store their positions
	float ddela1=0;
	float ddelb1=0;
	float ddelc1=0;
	float ddelmax=0;
	float ddelmed=0;
	float ddelmin=0;

	float ddeltx=m_nImgSpacing[0];
	float ddelty=m_nImgSpacing[1];
	float ddeltz=m_nImgSpacing[2];

	float dMax, dMin, dMid;//the max u values
	float dmaxP[3]={0,0,0};
	float dmedP[3]={0,0,0};
	float dminP[3]={0,0,0};
	float daP[3]={0,0,0};
	float dbP[3]={0,0,0};
	float dcP[3]={0,0,0};
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dmaxP[0]=dualr[0];
		dmaxP[1]=dualr[1];
		dmaxP[2]=dualr[2];
		dMin = ub;
		ddelmin=ddelty;
		dminP[0]=dubud[0];
		dminP[1]=dubud[1];
		dminP[2]=dubud[2];

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dmaxP[0]=dubud[0];
		dmaxP[1]=dubud[1];
		dmaxP[2]=dubud[2];
		dMin = ua;
		ddelmin=ddeltx;
		dminP[0]=dualr[0];
		dminP[1]=dualr[1];
		dminP[2]=dualr[2];
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dmedP[0]=dmaxP[0];
		dmedP[1]=dmaxP[1];
		dmedP[2]=dmaxP[2];
		dMax = uc;
		ddelmax=ddeltz;
		dmaxP[0]=ducbf[0];
		dmaxP[1]=ducbf[1];
		dmaxP[2]=ducbf[2];
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
		dmedP[0]=ducbf[0];
		dmedP[1]=ducbf[1];
		dmedP[2]=ducbf[2];
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dmedP[0]=dminP[0];
		dmedP[1]=dminP[1];
		dmedP[2]=dminP[2];
		dMin = uc;
		ddelmin=ddeltz;
		dminP[0]=ducbf[0];
		dminP[1]=ducbf[1];
		dminP[2]=ducbf[2];
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	ddela1=ddelmax; 
	 ddelb1=ddelmed;
	 ddelc1=ddelmin;
	daP[0]=dmaxP[0];daP[1]=dmaxP[1];daP[2]=dmaxP[2];
	dbP[0]=dmedP[0];dbP[1]=dmedP[1];dbP[2]=dmedP[2];
	dcP[0]=dminP[0];dcP[1]=dminP[1];dcP[2]=dminP[2];
	//
	float daPworld[3]={daP[0],daP[1],daP[2]};
	float dbPworld[3]={dbP[0],dbP[1],dbP[2]};
	float dcPworld[3]={dcP[0],dcP[1],dcP[2]};
	float dxyzPworld[3]={x,y,z};
	
	m_pBaseImgInfo->ImageToWorld(daPworld);
	m_pBaseImgInfo->ImageToWorld(dbPworld);
	m_pBaseImgInfo->ImageToWorld(dcPworld);
	m_pBaseImgInfo->ImageToWorld(dxyzPworld);
	// a b c as the vector
	float dva[3]={daPworld[0]-dxyzPworld[0],daPworld[1]-dxyzPworld[1],daPworld[2]-dxyzPworld[2]};
	float dvb[3]={dbPworld[0]-dxyzPworld[0],dbPworld[1]-dxyzPworld[1],dbPworld[2]-dxyzPworld[2]};
	float dvc[3]={dcPworld[0]-dxyzPworld[0],dcPworld[1]-dxyzPworld[1],dcPworld[2]-dxyzPworld[2]};
	float ddela=zxh::VectorOP_Magnitude(dva,3);
	float ddelb=zxh::VectorOP_Magnitude(dvb,3);
	float ddelc=zxh::VectorOP_Magnitude(dvc,3);
	float ca=zxh::VectorOP_Cosine(dva,dvb,3);
	float cb=zxh::VectorOP_Cosine(dvb,dvc,3);
	float cc=zxh::VectorOP_Cosine(dva,dvc,3);
	//define the three unit vectors
	float dva1[3]={0,0,0};
	dva1[0]=(float)daP[0]-x;
	dva1[1]=(float)daP[1]-y;
	dva1[2]=(float)daP[2]-z;
	float dvb1[3]={0,0,0};
	dvb1[0]=(float)dbP[0]-x;
	dvb1[1]=(float)dbP[1]-y;
	dvb1[2]=(float)dbP[2]-z;
	float dvc1[3]={0,0,0};
	dvc1[0]=(float)dcP[0]-x;
	dvc1[1]=(float)dcP[1]-y;
	dvc1[2]=(float)dcP[2]-z;
	//normalize
	zxh::VectorOP_Normalise(dva1,3);
	zxh::VectorOP_Normalise(dvb1,3);
	zxh::VectorOP_Normalise(dvc1,3);
	//
	float ca1=zxh::VectorOP_Cosine(dva1,dvb1,3);
	float cb1=zxh::VectorOP_Cosine(dvb1,dvc1,3);
	float cc1=zxh::VectorOP_Cosine(dva1,dvc1,3);
	//


	double dD1=1-pow(ca,2)-pow(cb,2)-pow(cc,2)+2*ca*cb*cc;
	double dA3=(1-pow(cb,2))/pow(ddela,2)+(1-pow(cc,2))/pow(ddelb,2)+(1-pow(ca,2))/pow(ddelc,2)+\
		2*(cb*cc-ca)/(ddela*ddelb)+2*(ca*cb-cc)/(ddela*ddelc)+2*(ca*cc-cb)/(ddelb*ddelc);
	double dB3=-2*ua*(1-pow(cb,2))/pow(ddela,2)-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc*(1-pow(ca,2))/pow(ddelc,2)+\
		2*(cb*cc-ca)/(ddela*ddelb)*(-ub-ua)+2*(ca*cb-cc)/(ddela*ddelc)*(-ua-uc)+2*(ca*cc-cb)/(ddelb*ddelc)*(-ub-uc);
	double dC3=(1-pow(cb,2))/pow(ddela,2)*pow(ua,2)+(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+(1-pow(ca,2))/pow(ddelc,2)*pow(uc,2)+\
		2*(cb*cc-ca)/(ddela*ddelb)*ua*ub+2*(ca*cb-cc)/(ddela*ddelc)*ua*uc+2*(ca*cc-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	//solve the equation
	double dRoots_3[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_3 = QuadraticRoots(dA3, dB3, dC3, dRoots_3);
	//solve 2 x
	double dA2=(1-pow(cc,2))/pow(ddelb,2)+1/pow(ddelc,2)+\
		2*(-cb)/(ddelb*ddelc);
	double dB2=-2*ub*(1-pow(cc,2))/pow(ddelb,2)-2*uc/pow(ddelc,2)+\
		2*(-cb)/(ddelb*ddelc)*(-ub-uc);
	double dC2=(1-pow(cc,2))/pow(ddelb,2)*pow(ub,2)+1/pow(ddelc,2)*pow(uc,2)+\
		2*(-cb)/(ddelb*ddelc)*ub*uc-dP2*dD1;
	bool bSolv_2 = QuadraticRoots(dA2, dB2, dC2, dRoots_2);

	double dRtMax_3 = max(dRoots_3[0], dRoots_3[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);
	if (dRtMax_3 > ua && bSolv_3 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_3;
	}
	else if (dRtMax_2 > ub && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(dP2) + uc;
	}
	////solve 3 x
	//double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//double dRoots_0[2] = {0}, dRoots_1[2] = {0};
	//
	//bool bSolv_0 = QuadraticRoots(da3, db3, dc3, dRoots_0);
	////solve 2 x
	//double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	//double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	//double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-dP2;
	//bool bSolv_1 = QuadraticRoots(da2, db2, dc2, dRoots_1);

	//double dRtMax_0 = max(dRoots_0[0], dRoots_0[1]);
	//double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);

	//if (dRtMax_0 > ua && bSolv_0 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_0;
	//}
	//else if (dRtMax_1 > ub && bSolv_1 == true)
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	//} 
	//else
	//{
	//	m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	//		sqrt(dP2) + uc;
	//}
}
void miiMinPathModel::ModUpWind_MultSp(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double ua = min(u_l, u_r);//x;find the triangle with the minimal m_dU£¨u£© value
	double ub = min(u_u, u_d);//y
	double uc = min(u_b, u_f);//z
	double ddela=0;
	double ddelb=0;
	double ddelc=0;
	double ddelmax=0;
	double ddelmed=0;
	double ddelmin=0;

	double ddeltx=m_nImgSpacing[0];
	double ddelty=m_nImgSpacing[1];
	double ddeltz=m_nImgSpacing[2];

	double dMax, dMin, dMid;
	if (ua > ub) 
	{
		dMax = ua;
		ddelmax=ddeltx;
		dMin = ub;
		ddelmin=ddelty;

	}
	else
	{
		dMax = ub;
		ddelmax=ddelty;
		dMin = ua;
		ddelmin=ddeltx;
	}
	if (uc > dMax)
	{

		dMid = dMax;
		ddelmed=ddelmax;
		dMax = uc;
		ddelmax=ddeltz;
	}
	else if (uc > dMin)
	{
		dMid = uc;
		ddelmed=ddeltz;
	}
	else
	{
		dMid = dMin;
		ddelmed=ddelmin;
		dMin = uc;
		ddelmin=ddeltz;
	}

	ua = dMax; ub = dMid; uc = dMin;//order the trangle points  by their m_dU(u) 
	ddela=ddelmax; ddelb=ddelmed;ddelc=ddelmin;

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;
	double dSimd=0;
	double dSimd1=0;
	double dSimd2=0;
	double dSimd3=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}


	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta2 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}

	if (costheta3 == 0)//current vessel director and last start vector
	{
		dSimd3 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else if (abs(costheta3) > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd3 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else
	{
		dSimd3 = - (((abs(costheta3) - dCosTheataMean)*(abs(costheta3) - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
		
	}
	//calculate the total direction parameter
	dSimd=dSimd1*dSimd1;
	//dSimd=dSimd3*dSimd3;
	//dSimd=dSimd1*dSimd1;
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	//solve 3 x
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	double da3=1/pow(ddela,2)+1/pow(ddelb,2)+1/pow(ddelc,2);
	double db3=-2*(1/pow(ddela,2)*ua+1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	double dc3=1/pow(ddela,2)*ua*ua+1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(da3, db3, dc3, dRoots_1);
	//solve 2 x
	double da2=1/pow(ddelb,2)+1/pow(ddelc,2);
	double db2=-2*(1/pow(ddelb,2)*ub+1/pow(ddelc,2)*uc);
	double dc2=1/pow(ddelb,2)*ub*ub+1/pow(ddelc,2)*uc*uc-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	bool bSolv_2 = QuadraticRoots(da2, db2, dc2, dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > ua && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > ub && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + uc;
	}
  //cout<< "("<<ddela << "," << ddelb << "," << ddelc<< ")"<<"\n ";
}
void miiMinPathModel::ModifiedUpWindPlusTime2CosthetaN_3S(int x, int y, int z, float costheta1,float costheta2,float costheta3,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;
	double dSimd=0;
	double dSimd1=0;
	double dSimd2=0;
	double dSimd3=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}


	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta2 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}

	if (costheta3 == 0)//current vessel director and last start vector
	{
		dSimd3 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else if (abs(costheta3) > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd3 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
	}
	else
	{
		dSimd3 = - (((abs(costheta3) - dCosTheataMean)*(abs(costheta3) - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd3 = exp(dSimd3);
		
	}
	//calculate the total direction parameter
	//dSimd=dSimd1*dSimd1+dSimd2*dSimd2;
	dSimd=dSimd3*dSimd3;
	//dSimd=dSimd1*dSimd1;
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	if(m<0)
		int m=0;
	if(dSimd1<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}

void miiMinPathModel::ModifiedUpWindPlusTime2Costheta_DeaF(int x, int y, int z, float costheta1,float costheta2,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;double dSimd=0;double dSimd1=0;double dSimd2=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;

	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}//(1-lambdas)*(dSimd*dSimd)
	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta2 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}//(1-lambdas)*(dSimd*dSimd)
	//dSimd=(dSimd1+dSimd2)/2;
	
	//dSimd=dSimd1*dSimd1+dSimd2*dSimd2;
	
   dSims = 1;
   dSimd=1;
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	if(m<0)
		int m=0;
	if(dSimd1<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
void miiMinPathModel::ModifiedUpWindPlusTime2CosthetaSM(int x, int y, int z, float costheta1,float costheta2,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.000001,omiga=0.000001;
	float lambdas=0.75,lambdav=1,lambdad=1,vslweight=200;
	double dSims=0;double dSimd=0;double dSimd1=0;double dSimd2=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if ( (dImgVal < m_fMean) && dImgVal ==m_fMean)
	{
		dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
		dSims = exp(dSims);
	}
	else
		dSims = 1;
	// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
	float mindSimd=- (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	mindSimd=exp(dSimd);
	float maxdSimd=- (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	maxdSimd=exp(maxdSimd);
	if (costheta1 <= 0)//current vessel director and current model vector
	{
		dSimd1 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else if (costheta1 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd1 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
	}
	else
	{
		dSimd1 = - (((costheta1 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd1 = exp(dSimd1);
		
	}//(1-lambdas)*(dSimd*dSimd)
	if (costheta2 <= 0)//current vessel director and last start vector
	{
		dSimd2 = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else if (costheta2 > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		dSimd2 = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
	}
	else
	{
		dSimd2 = - (((costheta2 - dCosTheataMean)*(costheta1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd2 = exp(dSimd2);
		
	}//(1-lambdas)*(dSimd*dSimd)
	//dSimd=(dSimd1+dSimd2)/2;
	//dSimd=dSimd1*dSimd1+dSimd2*dSimd2;
	dSimd=maxdSimd;
	if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
	{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*vslweight*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd+(1-lambdas)*dSimd+epsilion))+omiga);
	}	
	/*else
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
	else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/dSimd;
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd1;
	double m=m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	if(m<0)
		int m=0;
	if(dSimd1<0.001)
		int m=1;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
	//     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}



	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	if(m1!=0)
		double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	//	ofstream WriteFileTxt;
	// WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	//WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	//WriteFileTxt.close();

	//cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
void miiMinPathModel::ModifiedUpWindPlusTimeLast(int x, int y, int z, double costheta,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.001,omiga=0.01;
	float lambdas=0.25,lambdav=1,lambdad=1;
	double dSims=0;double dSimd=0;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
				if ( dImgVal < m_fMean)
				{
					dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSims = exp(dSims);
				}
				else
					dSims = 1;
// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
if (costheta <= 0)
	{
		 dSimd = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
		}	
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		(1/(((lambdas)*dSims*200*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*0+(1-lambdas)*dSimd*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*200*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*0+(1-lambdas)*dSimd*dSimd+epsilion))+omiga);
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}
	else if (costheta > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		 dSimd = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);	
		}
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		(1/(((lambdas)*dSims*200*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*0+(1-lambdas)*dSimd*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*200*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*0+(1-lambdas)*dSimd*dSimd+epsilion))+omiga);
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}
	else
	{
	    dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
		}
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);*/
		if(m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] == FMM_FAR)
		{m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		(1/(((lambdas)*dSims*200*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*0+(1-lambdas)*dSimd*dSimd+epsilion))+omiga)*(1/(((lambdas)*dSims*200*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*0+(1-lambdas)*dSimd*dSimd+epsilion))+omiga);
		}	
		/*else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/((1-lambdas)*(dSimd*dSimd));*/
		else
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]/(dSimd*dSimd);
	}//(1-lambdas)*(dSimd*dSimd)
	m_dSims[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSims;
	m_dSimd[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]=dSimd;
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
 //     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}

    
	
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
		if(m1!=0)
			double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	 //	ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/testing/m_dWP.txt",ios::app);//Add by JDQ
	 //WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<dSims<<","<<m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<dSimd<<","<<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_nFmMap[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
	
		 //cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
// Function Name: ModifiedUpWindPlusTimeMixed()
//
// Parameters: x, y, z: coordinate corespond left(right),up(down),front(back),costheta,image intensity
//
// Description: add the intensity and direction parameter to the potential calculation; potential function doesn't change in upwind ;potential function =1/( v*s+d+e)
//
// Returns: 
//
void miiMinPathModel::ModifiedUpWindPlusTimeMixed(int x, int y, int z, double costheta,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.001;
	float lambdas=0.5,lambdav=1,lambdad=1;
	double dSims;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
				if ( dImgVal < m_fMean)
				{
					dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSims = exp(dSims);
				}
				else
					dSims = 1;
// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
if (costheta <= 0)
	{
		
		double dSimd = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/((lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd*dSimd+epsilion));*/
	  m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	    /*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims+m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+dSimd+0.01);*/
	
	}
	else if (costheta > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		double dSimd = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		
		/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/((lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd*dSimd+epsilion));*/
		/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims+m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+dSimd+0.01);*/
		 m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	}
	else
	{
		double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims+m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+dSimd+0.01);*/
		/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/((lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd*dSimd+epsilion));*/
		 m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	}
	
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
 //     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}

    
	
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
		if(m1!=0)
			double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	 //	ofstream WriteFileTxt;
	 // WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/svplusdplusepsilonwithoutintensityupdating/m_dWP.txt",ios::app);//Add by JDQ
	 // WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ")"<<"\n ";//add by JDQ
	 //WriteFileTxt.close();
	
		 //cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
void miiMinPathModel::ModifiedUpWindPlusTime1(int x, int y, int z, double costheta,const short *sImgData)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// get the intensity to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.001;
	float lambdas=0.5,lambdav=1,lambdad=1;
	double dSims;
	double dImgVal = (double)sImgData[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
				if ( dImgVal < m_fMean)
				{
					dSims = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSims = exp(dSims);
				}
				else
					dSims = 1;
// add the direction to the cost function
	/*if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSimd = exp(dSimd);*/
if (costheta <= 0)
	{
		double dSimd = - (((0.1 - dCosTheataMean)*(0.1 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
	
	}
	else if (costheta > 0.866)//cos(theta)=0.866  --->theta=30 degree.
	{
		double dSimd = - (((0.866 - dCosTheataMean)*(0.866 - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);	
		
	}
	else
	{
		double dSimd = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
		dSimd = exp(dSimd);
		
		m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
		
	}
	
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*dSimd*dSimd+dSimd*dSimd);
double m2=1/(dSims*m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+1);*/
	//if(x>100&&x<300&&y>70&&y<270&&z>120&&z<320&&m_sVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]!=0)
	//{
 //     ofstream WriteFileTxt;
	//  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotential/m_dWP.txt",ios::app);//Add by JDQ
	//  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ")"<<"\n ";//add by JDQ
	// WriteFileTxt.close();
	//}

    
	
	/*m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
		1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	double m=1/(lambdas*dSims*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSimd+epsilion);
	*/
	double m1=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
		if(m1!=0)
			double m2=m1;
	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
	 	ofstream WriteFileTxt;
	  WriteFileTxt.open ("F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/ModifiedPotentialTmePlusMixed/originalpotentialnotchangedinupwind/withintensityupdated/m_dWP.txt",ios::app);//Add by JDQ
	  WriteFileTxt<< "("<<x << "," << y << "," << z<<"," <<m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<<","<<m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]<< ")"<<"\n ";//add by JDQ
	 WriteFileTxt.close();
	
		 //cout<< "("<<x << "," << y << "," << z<<"," <<m2<< ")"<<"\n ";//add by JDQ
  
}
void miiMinPathModel::NewUpWindTimePlusMixed(int x, int y, int z, double costheta)
{
	int u_l_x = x - 1, u_l_y = y,	  u_l_z = z;  
	int u_r_x = x + 1, u_r_y = y,	  u_r_z = z;
	int u_u_x = x,     u_u_y = y - 1, u_u_z = z;
	int u_d_x = x,     u_d_y = y + 1, u_d_z = z;
	int u_b_x = x,     u_b_y = y,     u_b_z = z - 1;
	int u_f_x = x,     u_f_y = y,     u_f_z = z + 1;
	double u_l, u_r, u_u, u_d, u_b, u_f;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_l = (double)m_dU[u_l_z * m_nImgWY * m_nImgWX + u_l_y * m_nImgWX + u_l_x];

	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_r = (double)m_dU[u_r_z * m_nImgWY * m_nImgWX + u_r_y * m_nImgWX + u_r_x];

	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_u = (double)m_dU[u_u_z * m_nImgWY * m_nImgWX + u_u_y * m_nImgWX + u_u_x];

	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_d = (double)m_dU[u_d_z * m_nImgWY * m_nImgWX + u_d_y * m_nImgWX + u_d_x];

	if (u_b_z < 0 || u_b_z >= m_nImgWZ)
		u_b = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_b = (double)m_dU[u_b_z * m_nImgWY * m_nImgWX + u_b_y * m_nImgWX + u_b_x];

	if (u_f_z < 0 || u_f_z >= m_nImgWZ)
		u_f = (double)m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];
	else
		u_f = (double)m_dU[u_f_z * m_nImgWY * m_nImgWX + u_f_y * m_nImgWX + u_f_x];

	double a = min(u_l, u_r);//find the triangle with the minimal m_dU£¨u£© value
	double b = min(u_u, u_d);
	double c = min(u_b, u_f);

	double dMax, dMin, dMid;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}
	if (c > dMax)
	{
		dMid = dMax;
		dMax = c;
	}
	else if (c > dMin)
	{
		dMid = c;
	}
	else
	{
		dMid = dMin;
		dMin = c;
	}

	a = dMax; b = dMid; c = dMin;//order the trangle points  by their m_dU(u) 

	// add the direction to the cost function
	double dCosTheataMean = 1, dCosTheataStdDev = 0.5,epsilion=0.001;
	float lambdas=0.5,lambdav=1,lambdad=1;
	if( costheta < 0.1 ) costheta = 0.1 ; 
	if( costheta > 0.866 ) costheta = 0.866 ; 
	double dSim = - (((costheta - dCosTheataMean)*(costheta - dCosTheataMean)) / (dCosTheataStdDev * dCosTheataStdDev)) / 2;
	dSim = exp(dSim);
	m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
	1/(lambdas*m_sUnseenImgIntensity[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSim*dSim+epsilion);
	double m1=	1/(lambdas*m_sUnseenImgIntensity[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]*m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]+(1-lambdas)*dSim*dSim+epsilion);
	double m2=m_sNormVes[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x];

	double dRoots_1[2] = {0}, dRoots_2[2] = {0};
	bool bSolv_1 = QuadraticRoots(3, -2*(a+b+c), \
		a*a+b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_1);

	bool bSolv_2 = QuadraticRoots(2, -2*(b+c), \
		b*b+c*c-m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x], dRoots_2);

	double dRtMax_1 = max(dRoots_1[0], dRoots_1[1]);
	double dRtMax_2 = max(dRoots_2[0], dRoots_2[1]);

	if (dRtMax_1 > a && bSolv_1 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_1;
	}
	else if (dRtMax_2 > b && bSolv_2 == true)
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = dRtMax_2;
	} 
	else
	{
		m_dU[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x] = \
			sqrt(m_dP2[z * m_nImgWY * m_nImgWX + y * m_nImgWX + x]) + c;
	}
}


bool miiMinPathModel::ConvS(miiCNode<double,int> iNewNode,const zxhImageInfo *pImageInfo,float nStartPoint[3])
{
	Speed_BackTrackCV(iNewNode,pImageInfo,nStartPoint,15);

	NormSpeed();
	return true;
}
bool miiMinPathModel::NormSpeed()
{
	float maxspeedvalue=0;
	float fsum=0;
	int nsumNUM=0;
	int nstart=m_vSpValue.size()-1000;
	if(nstart<0)nstart=0;
	int nend=m_vSpValue.size()-1;
	miiCNode<double,int> miTempNode=m_vSpValue[nend];
	for(int j=nstart;j<=nend;j++)//calculate the pre-1000 value
	{
		if(m_vSpValue[j].val>maxspeedvalue)
		{
			maxspeedvalue=m_vSpValue[j].val;
			
		}
		fsum=fsum+m_vSpValue[j].val;
		nsumNUM++;
	}
	float fMeanSpeed=fsum/nsumNUM;
	m_fMinSpeed=zxh::minf(m_fMinSpeed,fMeanSpeed);
	miTempNode.val=fMeanSpeed/(maxspeedvalue+0.001);	
	if(miTempNode.val<m_fMinSpeed)
	{
		m_nLn++;
	}
	else
	{
		m_nLn=0;
	}
	m_vSpNorValue.push_back(miTempNode);
	return true;
}
// Function Name: PotentialFunction()
//
// Parameters: *sImgData: the pointer for the data of 3D image
//				nMethod: the option for the potential function 
//						'1' - Gradient
//						'2' - Vesselness + Intensity
//						'3' - Vesselness + Intensity + Direction(*)
//                      '4' - Vesselness + Intensity + Direction(+)                     
// Description: the potential function
//
// Returns: 
//
bool miiMinPathModel::PotentialFunction(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, int nMethod)
{
	if (sImgData == NULL)
	{
		std::cout << "The pointer cannot be null!" <<endl;
		return false;
	}

	switch (nMethod)
	{
	case 1:
		break;
	case 2:
		Vesselness2Invs(sImgData, pImgInfo, sVslsData, pVslsInfo, m_dP2);
		break;
	case 3:
		VesselnessSim2Invs(sImgData, pImgInfo, sVslsData, pVslsInfo, m_dP2);
		break;
	case 4:
		VesselnessSim2InvsPlus(sImgData, pImgInfo, sVslsData, pVslsInfo, m_dP2);
	case 5:
		VesselnessSim2InvsPlusTimeMixed(sImgData, pImgInfo, sVslsData, pVslsInfo, m_dP2);
	default:
		break;
	}

	return true;	
}

// Function Name: Vesselness2Invs()
//
// Parameters: *sImgData: the pointer for the data of 3D image
//			   *dImgG: the pointer for the inverse of the square of Gradient
//
// Description: calculate the inverse of the square of Vesselness 
//
// Returns: 
//
bool miiMinPathModel::Vesselness2Invs(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2)
{
	if (sImgData == NULL || sVslsData == NULL || dP2 == NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}
  //char *chResultPath ="F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/original/m_dP.txt";//Add byJDQ
  //    ofstream WriteFileTxt(chResultPath);//Add by JDQ
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = \
					1 / ((double)sVslsData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] * \
					sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] + 0.0000000001);
       
  /*   WriteFileTxt << "("<<ix << "," << iy << "," << iz<<"," <<dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]<< ")"<<"\n ";
	*/
			}
		}
	}

	return true;
}

// Function Name: VesselnessSim2Invs()
//
// Parameters: *sImgData: the pointer for the data of 3D image
//			   *dImgG: the pointer for the inverse of the square of Gradient
//
// Description: calculate the inverse of the square of Vesselness and the intensity similarty 
//
// Returns: 
//
bool miiMinPathModel::VesselnessSim2Invs(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2)
{
	if (sImgData == NULL || sVslsData == NULL || dP2 == NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}

	if (!m_bReadMeanFlg)
	{
		GetMeanStdDev(sImgData, m_fMean, m_fStdDev);
	}	

	double dSim;
	double dImgVal;
	//char *chResultPath ="F:/Coronary_0/miiMinPathModelApp_results/mod00_to_unseen02_results_physnoimgres/ori_vessel0/Thetalower30/NearestTangentPlane/WithMask/original/m_dP.txt";//Add byJDQ
 //     ofstream WriteFileTxt(chResultPath);//Add by JDQ
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				dImgVal = (double)sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				if ( dImgVal < m_fMean)
				{
					dSim = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSim = exp(dSim);
				}
				else
					dSim = 1;

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
				if (nCord[1] >= m_nImgWY)
					nCord[1] = m_nImgWY - 1;
				if (nCord[2] >= m_nImgWZ)
					nCord[2] = m_nImgWZ - 1;
				//dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = \
					1 / (((double)sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]]) * \
					dSim+((double)sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]])+ 1);
				//dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = \
					1 / (((double)sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]]) * \
					dSim+dSim+ 1);
				//dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = \
					1 / (((double)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]]) + 1) * \
					dSim);
				dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = \
					1 / ((double)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]]) * \
					dSim + 1);
				//if(ix>100&&ix<300&&iy>70&&iy<270&&iz>120&&iz<320&&sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]]!=0)
				// WriteFileTxt << "("<<ix << "," << iy << "," << iz<<"," <<dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]<< ")"<<"\n ";//add by JDQ
			}
		}
	}

 	return true;
}
bool miiMinPathModel::VesselnessSim2InvsPlus(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2)
{
	if (sImgData == NULL || sVslsData == NULL || dP2 == NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}

	if (!m_bReadMeanFlg)
	{
		GetMeanStdDev(sImgData, m_fMean, m_fStdDev);
	}	

	float MaxVsls=0;
	double dMaxSim=0;
	double dSim;
	double dImgVal;
	float lambdas=1, lambdav=1, epislon=ZXH_FloatInfinitesimal ;
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				
				dImgVal = (double)sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				if ( dImgVal < m_fMean)
				{
					dSim = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSim = exp(dSim);
				}
				else
					dSim = 1;
				dMaxSim=zxh::maxf(dSim,dMaxSim);
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
				dImgVal = (double)sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				if ( dImgVal < m_fMean)
				{
					dSim = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSim = exp(dSim)/dMaxSim;
				}
				else
					dSim = 1;

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
				if (nCord[1] >= m_nImgWY)
					nCord[1] = m_nImgWY - 1;
				if (nCord[2] >= m_nImgWZ)
					nCord[2] = m_nImgWZ - 1;
				dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = ( lambdav * (double)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]])/MaxVsls + \
					lambdas * dSim+ epislon );
			}
		}
	}

	return true;
}
bool miiMinPathModel::VesselnessSim2InvsPlusTimeMixed(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2)
{
	if (sImgData == NULL || sVslsData == NULL || dP2 == NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}

	if (!m_bReadMeanFlg)
	{
		GetMeanStdDev(sImgData, m_fMean, m_fStdDev);
	}	

	short MaxVsls=0;
	double dMaxSim=0;
	double dSim;
	double dImgVal;
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				
				dImgVal = (double)sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				if ( dImgVal < m_fMean)
				{
					dSim = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSim = exp(dSim);
				}
				else
					dSim = 1;
				dMaxSim=zxh::maxf(dSim,dMaxSim);
				MaxVsls=zxh::maxf(sVslsData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix],MaxVsls);
			}
		}
	}
	
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				dImgVal = (double)sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				if ( dImgVal < m_fMean)
				{
					dSim = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSim = exp(dSim)/dMaxSim;
				}
				else
					dSim = 1;

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
				if (nCord[1] >= m_nImgWY)
					nCord[1] = m_nImgWY - 1;
				if (nCord[2] >= m_nImgWZ)
					nCord[2] = m_nImgWZ - 1;
				dP2[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = (((double)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]])/MaxVsls)*\
					dSim+1);
			}
		}
	}

	return true;
}
bool miiMinPathModel::NewPotentialFunction(const short *sImgData, const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, int nMethod)
{
	if (sImgData == NULL)
	{
		std::cout << "The pointer cannot be null!" <<endl;
		return false;
	}

	switch (nMethod)
	{
	case 1:
		break;

	case 2:
		VesselnessInit(pImgInfo, sVslsData, pVslsInfo, m_sVes,m_sNormVes);
	    UnseenImgIntensitySet(sImgData,m_sUnseenImgIntensity);
	default:
		break;
	}

	return true;	
}
bool miiMinPathModel::VesselnessInit(const zxhImageInfo *pImgInfo, \
	const short *sVslsData, const zxhImageInfo *pVslsInfo, double *m_sVes,double *m_sNormVes)
{
	if (sVslsData == NULL || m_sVes == NULL ||m_sNormVes==NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}
		float MaxVsls=0;
	
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
				if (nCord[0] <0)
					nCord[0] = 0;
				if (nCord[1] >= m_nImgWY)
					nCord[1] = m_nImgWY - 1;
				if (nCord[1] < 0)
					nCord[1] = 0;
				if (nCord[2] >= m_nImgWZ)
					nCord[2] = m_nImgWZ - 1;
				if (nCord[2] < 0)
					nCord[2] = 0;
				m_sVes[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]=(double)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]]);
				m_sNormVes[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]=(double)(sVslsData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]])/MaxVsls;
			}
		}
	}

	return true;
}
bool miiMinPathModel::UnseenImgIntensitySet(const short *sImgData, double *m_sUnseenImgIntensity)
{
	if (sImgData == NULL ||m_sUnseenImgIntensity== NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}

	if (!m_bReadMeanFlg)
	{
		GetMeanStdDev(sImgData, m_fMean, m_fStdDev);
	}	
	double dSim;
	double dImgVal;
	for (int iz = 0; iz < m_nImgWZ; iz++)
	for (int iy = 0; iy < m_nImgWY; iy++)
	for (int ix = 0; ix < m_nImgWX; ix++)
		{
				
				dImgVal = (double)sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				if ( dImgVal < m_fMean)
				{
					dSim = - (((dImgVal - m_fMean)*(dImgVal - m_fMean)) / (m_fStdDev * m_fStdDev + 0.00000001)) / 2;
					dSim = exp(dSim);
				}
				else
					dSim = 1;
			m_sUnseenImgIntensity[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]=dSim;
	}

	return true;
}

// Function Name: GetMeanStdDev()
//
// Parameters: *sImgData: the pointer for the data of 3D image
//			   fMean: the mean of the CA's intensity
//			   fStdDev: the std-deviation of the CA's intensity
//
// Description: get the mean and std-deviation of the CA's intensity from the neighbor of the start point 
//
// Returns: 
//
bool miiMinPathModel::GetMeanStdDev(const short *sImgData, double &fMean, double &fStdDev)
{
	if (sImgData == NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}

	int nSearchRange = 10;
	int nCount = 0;
	double nSum = 0;
	int nXmin = m_iStartPoint.x - nSearchRange;
	int nXmax = m_iStartPoint.x + nSearchRange;

	int nYmin = m_iStartPoint.y - nSearchRange;
	int nYmax = m_iStartPoint.y + nSearchRange;

	int nZmin = m_iStartPoint.z - nSearchRange;
	int nZmax = m_iStartPoint.z + nSearchRange;

	// calculate mean
	for (int iz = nZmin; iz < nZmax; iz++)
	{
		for (int iy = nYmin; iy < nYmax; iy++)
		{
			for (int ix = nXmin; ix < nXmax; ix++)
			{
				if (sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] > 0)
				{
					nCount++;
					nSum += (double)sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				}
			}
		}
	}

	fMean = nSum / (nCount + 0.00000000001);
	nSum = 0;

	// calculate deviation
	for (int iz = nZmin; iz < nZmax; iz++)
	{
		for (int iy = nYmin; iy < nYmax; iy++)
		{
			for (int ix = nXmin; ix < nXmax; ix++)
			{
				if (sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] > 0)
				{
					double temp = sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] - fMean;
					temp *= temp;
					nSum += temp;
				}
			}
		}
	}

	fStdDev = sqrt(nSum / (nCount + 1));

	return true;
}

// Function Name: SetModelPoints()
//
// Parameters: vModelPoints: the original CA's model
//			   nModelStartPoint: the start point from the aorta center which is automatic detected
//
// Description: set the model for FMM from the external file
//
// Returns: 
//
bool miiMinPathModel::SetModelPoints(vector< miiCNode<double, float> > vModelPoints, float nModelStartPoint[3])
{
	if (vModelPoints.size() == 0)
	{
		cout << "The coronary model is null!" << endl;
		return false;
	}

	if (m_vModelPointsWorld.size() > 0)
	{
		m_vModelPoints.clear();
	}

	

	// convert the coordinate from world to image
	float fCord[3];
	float DistModelCur=0;
	fCord[0] = nModelStartPoint[0];
	fCord[1] = nModelStartPoint[1];
	fCord[2] = nModelStartPoint[2];
	miiCNode<double, float> fTempPointWorld;
	miiCNode<double, int> iTempPoint;
	fTempPointWorld.x =fCord[0];
	fTempPointWorld.y =fCord[1];
	fTempPointWorld.z =fCord[2];
	m_iOrgStartPointWorld.x=fCord[0];
	m_iOrgStartPointWorld.y=fCord[1];
	m_iOrgStartPointWorld.z=fCord[2];
	m_vModelPointsWorld.push_back(fTempPointWorld);
	m_pBaseImgInfo->WorldToImage(fCord);
	iTempPoint.x = (int)(fCord[0]+0.5);
	iTempPoint.y = (int)(fCord[1]+0.5);
	iTempPoint.z = (int)(fCord[2]+0.5);
	if (iTempPoint.x >= (m_nImgWX))
		iTempPoint.x = m_nImgWX - 1;
	if (iTempPoint.y >= m_nImgWY)
		iTempPoint.y = m_nImgWY - 1;
	if (iTempPoint.z >= m_nImgWZ)
		iTempPoint.z = m_nImgWZ - 1;
	m_vModelPoints.push_back(iTempPoint);
	float vPrevPointWorld[3] = {0};
    float vCurrPointWorld[3] = {0} ;
	for (int i =1; i < vModelPoints.size(); i++)
	{
		// convert the coordinate from world to image
		fCord[0] = vModelPoints[i].x;
		fCord[1] = vModelPoints[i].y;
		fCord[2] = vModelPoints[i].z;
		fTempPointWorld.x =fCord[0];
	    fTempPointWorld.y =fCord[1];
	    fTempPointWorld.z =fCord[2];
		
		m_vModelPointsWorld.push_back(fTempPointWorld);
		m_pBaseImgInfo->WorldToImage(fCord);
		iTempPoint.x = (int)(fCord[0]+0.5);
		iTempPoint.y = (int)(fCord[1]+0.5);
		iTempPoint.z = (int)(fCord[2]+0.5);
		if (iTempPoint.x >= m_nImgWX)
			iTempPoint.x = m_nImgWX - 1;
		if (iTempPoint.y >= m_nImgWY)
			iTempPoint.y = m_nImgWY - 1;
		if (iTempPoint.z >= m_nImgWZ)
			iTempPoint.z = m_nImgWZ - 1;
		if (iTempPoint.x <0)
			iTempPoint.x = 0;
		if (iTempPoint.y <0)
			iTempPoint.y = 0;
		if (iTempPoint.z <0)
			iTempPoint.z = 0;
		//std::cout<<iTempPoint.x<<","<<iTempPoint.y<<","<<iTempPoint.z<<"\n";
		m_vModelPoints.push_back(iTempPoint);

		vCurrPointWorld[0] = m_vModelPointsWorld[i].x ;
		vCurrPointWorld[1] = m_vModelPointsWorld[i].y ;
		vCurrPointWorld[2] = m_vModelPointsWorld[i].z ;
		vPrevPointWorld[0] = m_vModelPointsWorld[i-1].x ;
		vPrevPointWorld[1] = m_vModelPointsWorld[i-1].y ;
		vPrevPointWorld[2] = m_vModelPointsWorld[i-1].z ;

        if (i>1)
		{
			DistModelCur=DistModelCur+zxh::VectorOP_Distance(vPrevPointWorld, vCurrPointWorld, 3 ) ;  
		}

	}
	//calculate the mean distance of the model points
	m_meandis_of_coronarymodel=DistModelCur/(m_vModelPointsWorld.size()-2);
	return true;
}
bool miiMinPathModel::SetModelPointsStartPointCorrect(vector< miiCNode<double, float> > vModelPoints, float fModelStartPoint[3],float fREFStartPoint[3])
{
	if (vModelPoints.size() == 0)
	{
		cout<< "The coronary model is null!" << endl;
		return false;
	}

	if (m_vModelPointsWorld.size() > 0)
	{
		m_vModelPoints.clear();
	}

	

	// convert the coordinate from world to image
	float fCord[3],fMoveCord[3];
	float DistModelCur=0;
	fCord[0] = fModelStartPoint[0];
	fCord[1] = fModelStartPoint[1];
	fCord[2] = fModelStartPoint[2];
	fMoveCord[0] = fREFStartPoint[0]-fModelStartPoint[0];
	fMoveCord[1] = fREFStartPoint[1]-fModelStartPoint[1];
	fMoveCord[2] = fREFStartPoint[2]-fModelStartPoint[2];
	miiCNode<double, float> fTempPointWorld;
	miiCNode<double, int> iTempPoint;
	fTempPointWorld.x =fCord[0];
	fTempPointWorld.y =fCord[1];
	fTempPointWorld.z =fCord[2];
	m_iOrgStartPointWorld.x=fCord[0];
	m_iOrgStartPointWorld.y=fCord[1];
	m_iOrgStartPointWorld.z=fCord[2];
	m_vModelPointsWorld.push_back(fTempPointWorld);
	m_pBaseImgInfo->WorldToImage(fCord);
	iTempPoint.x = (int)(fCord[0]+0.5);
	iTempPoint.y = (int)(fCord[1]+0.5);
	iTempPoint.z = (int)(fCord[2]+0.5);
	if (iTempPoint.x >= (m_nImgWX))
		iTempPoint.x = m_nImgWX - 1;
	if (iTempPoint.y >= m_nImgWY)
		iTempPoint.y = m_nImgWY - 1;
	if (iTempPoint.z >= m_nImgWZ)
		iTempPoint.z = m_nImgWZ - 1;
	m_vModelPoints.push_back(iTempPoint);
	float vPrevPointWorld[3] = {0};
    float vCurrPointWorld[3] = {0} ;
	for (int i =1; i < vModelPoints.size(); i++)
	{
		// convert the coordinate from world to image
		fCord[0] = vModelPoints[i].x;
		fCord[1] = vModelPoints[i].y;
		fCord[2] = vModelPoints[i].z;
		fTempPointWorld.x =fCord[0];
	    fTempPointWorld.y =fCord[1];
	    fTempPointWorld.z =fCord[2];
		
		m_vModelPointsWorld.push_back(fTempPointWorld);
		m_pBaseImgInfo->WorldToImage(fCord);
		iTempPoint.x = (int)(fCord[0]+0.5);
		iTempPoint.y = (int)(fCord[1]+0.5);
		iTempPoint.z = (int)(fCord[2]+0.5);
		if (iTempPoint.x >= m_nImgWX)
			iTempPoint.x = m_nImgWX - 1;
		if (iTempPoint.y >= m_nImgWY)
			iTempPoint.y = m_nImgWY - 1;
		if (iTempPoint.z >= m_nImgWZ)
			iTempPoint.z = m_nImgWZ - 1;
		if (iTempPoint.x <0)
			iTempPoint.x = 0;
		if (iTempPoint.y <0)
			iTempPoint.y = 0;
		if (iTempPoint.z <0)
			iTempPoint.z = 0;
		//std::cout<<iTempPoint.x<<","<<iTempPoint.y<<","<<iTempPoint.z<<"\n";
		m_vModelPoints.push_back(iTempPoint);

		

		if (i>1)
		{
			vCurrPointWorld[0] = m_vModelPointsWorld[i].x ;
			vCurrPointWorld[1] = m_vModelPointsWorld[i].y ;
			vCurrPointWorld[2] = m_vModelPointsWorld[i].z ;
			vPrevPointWorld[0] = m_vModelPointsWorld[i-1].x ;
			vPrevPointWorld[1] = m_vModelPointsWorld[i-1].y ;
			vPrevPointWorld[2] = m_vModelPointsWorld[i-1].z ;

			DistModelCur=DistModelCur+zxh::VectorOP_Distance(vPrevPointWorld, vCurrPointWorld, 3 ) ;  
		}

	}
	for (int i =1; i < m_vModelPointsWorld.size(); i++)
	{
		vCurrPointWorld[0] = m_vModelPointsWorld[i].x ;
		vCurrPointWorld[1] = m_vModelPointsWorld[i].y ;
		vCurrPointWorld[2] = m_vModelPointsWorld[i].z ;
		vPrevPointWorld[0] = m_vModelPointsWorld[i-1].x ;
		vPrevPointWorld[1] = m_vModelPointsWorld[i-1].y ;
		vPrevPointWorld[2] = m_vModelPointsWorld[i-1].z ;

		m_ftoallength_model=m_ftoallength_model+zxh::VectorOP_Distance(vPrevPointWorld, vCurrPointWorld, 3 ) ;  
	}
	for (int i =0; i < m_vModelPointsWorld.size(); i++)
	{
		m_vModelPointsWorld[i].x=m_vModelPointsWorld[i].x+fMoveCord[0];
		m_vModelPointsWorld[i].y=m_vModelPointsWorld[i].y+fMoveCord[1];
		m_vModelPointsWorld[i].z=m_vModelPointsWorld[i].z+fMoveCord[2];
	}
	//calculate the mean distance of the model points
	m_meandis_of_coronarymodel=DistModelCur/(m_vModelPointsWorld.size()-2);

	
	return true;
}
bool miiMinPathModel::SetModelPointsWithIntensity(vector< miiCNode<double, float> > vModelPointsWithIntensity, float nModelStartPointWithIntensity[5])
{
	if (vModelPointsWithIntensity.size() == 0)
	{
		cout << "The coronary model is null!" << endl;
		return false;
	}

	if (m_vModelPointsWorldWithIntensity.size() > 0)
	{
		m_vModelPointsWorldWithIntensity.clear();
	}

	// convert the coordinate from world to image
	float fCord[5];
	fCord[0] = nModelStartPointWithIntensity[0];
	fCord[1] = nModelStartPointWithIntensity[1];
	fCord[2] = nModelStartPointWithIntensity[2];
	fCord[3]=nModelStartPointWithIntensity[3];
	fCord[4]=nModelStartPointWithIntensity[4];
	miiCNode<double, float> fTempPointWorld;
	miiCNode<double, int> iTempPoint;
	fTempPointWorld.x =fCord[0];
	fTempPointWorld.y =fCord[1];
	fTempPointWorld.z =fCord[2];
	fTempPointWorld.val=fCord[3];
	m_vModelPointsWorldWithIntensity.push_back(fTempPointWorld);
	for (int i =1; i < vModelPointsWithIntensity.size(); i++)
	{
		

		// copy the model points and its intensity
		fCord[0] = vModelPointsWithIntensity[i].x;
		fCord[1] = vModelPointsWithIntensity[i].y;
		fCord[2] = vModelPointsWithIntensity[i].z;
		fCord[3]=vModelPointsWithIntensity[i].val;
		fTempPointWorld.x =fCord[0];
	    fTempPointWorld.y =fCord[1];
	    fTempPointWorld.z =fCord[2];
		fTempPointWorld.val =fCord[3];
		m_vModelPointsWorldWithIntensity.push_back(fTempPointWorld);

      	}
	
	
	return true;
}

	bool miiMinPathModel::SetModelPointsWithIntensityCorrect(vector< miiCNode<double, float> > vModelPointsWithIntensity, float nModelStartPointWithIntensity[5],float fREFStartPoint[3])
	{
		if (vModelPointsWithIntensity.size() == 0)
		{
			cout << "The coronary model is null!" << endl;
			return false;
		}

		if (m_vModelPointsWorldWithIntensity.size() > 0)
		{
			m_vModelPointsWorldWithIntensity.clear();
		}
		// convert the coordinate from world to image
		float fCord[5],fMoveCord[3];;
		fCord[0] = nModelStartPointWithIntensity[0];
		fCord[1] = nModelStartPointWithIntensity[1];
		fCord[2] = nModelStartPointWithIntensity[2];
		fCord[3]=nModelStartPointWithIntensity[3];
		fCord[4]=nModelStartPointWithIntensity[4];
		//calculate the move vector
		fMoveCord[0] = fREFStartPoint[0]-nModelStartPointWithIntensity[0];
		fMoveCord[1] = fREFStartPoint[1]-nModelStartPointWithIntensity[1];
		fMoveCord[2] = fREFStartPoint[2]-nModelStartPointWithIntensity[2];
		miiCNode<double, float> fTempPointWorld;
		miiCNode<double, int> iTempPoint;
		fTempPointWorld.x =fCord[0];
		fTempPointWorld.y =fCord[1];
		fTempPointWorld.z =fCord[2];
		fTempPointWorld.val=fCord[3];
		m_vModelPointsWorldWithIntensity.push_back(fTempPointWorld);
		for (int i =1; i < vModelPointsWithIntensity.size(); i++)
		{


			// copy the model points and its intensity
			fCord[0] = vModelPointsWithIntensity[i].x;
			fCord[1] = vModelPointsWithIntensity[i].y;
			fCord[2] = vModelPointsWithIntensity[i].z;
			fCord[3]=vModelPointsWithIntensity[i].val;
			fTempPointWorld.x =fCord[0];
			fTempPointWorld.y =fCord[1];
			fTempPointWorld.z =fCord[2];
			fTempPointWorld.val =fCord[3];
			m_vModelPointsWorldWithIntensity.push_back(fTempPointWorld);

		}
		for (int i =0; i < m_vModelPointsWorldWithIntensity.size(); i++)
		{
			m_vModelPointsWorldWithIntensity[i].x=m_vModelPointsWorldWithIntensity[i].x+fMoveCord[0];
			m_vModelPointsWorldWithIntensity[i].y=m_vModelPointsWorldWithIntensity[i].y+fMoveCord[1];
			m_vModelPointsWorldWithIntensity[i].z=m_vModelPointsWorldWithIntensity[i].z+fMoveCord[2];
		}

		return true;
	}

// Function Name: GetModelPoints()
//
// Parameters: 
//
// Description: get the corrected model 
//
// Returns: 
//
bool miiMinPathModel::GetModelPoints(vector< miiCNode<double, float> > &vModelPoints)
{
	if (vModelPoints.size() > 0)
	{
		vModelPoints.clear();
	}

	for (int i = 0; i < m_vModelPoints.size(); i++)
	{
		// convert the coordinate from world to image
		float fCord[3];
		fCord[0] = (float)(m_vModelPoints[i].x);
		fCord[1] = (float)(m_vModelPoints[i].y);
		fCord[2] = (float)(m_vModelPoints[i].z);
		m_pBaseImgInfo->ImageToWorld(fCord);
		miiCNode<double, float> iTempPoint;
		iTempPoint.x = fCord[0];
		iTempPoint.y = fCord[1];
		iTempPoint.z = fCord[2];

		vModelPoints.push_back(iTempPoint);
	}

	return true;
}

bool miiMinPathModel::CorrectModelPoint()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model of the coronary is null!" << endl;
		return false;
	}

	// correct the model
	float fDist[3];

	// find the nearest point in the model for the candidate point of FMM
	float fMinDistmm2 = 100000000;//int nMinDist2 = 100000000;Change by JDQ
	float nDistPPmm2=0;int nMinPos=0;//int nDistPP2, nMinPos;Change by JDQ
	for (int i = 0; i< m_vModelPointsWorld.size(); i++)
	{
		nDistPPmm2=CalcDistmm2(m_iSgmtPointWorld,m_vModelPointsWorld[i]);
		if (nDistPPmm2 <fMinDistmm2)
		{
			fMinDistmm2 = nDistPPmm2;
			nMinPos = i;
			m_nModelPos = i;//Add by JDQ

		}
	}
	//CoutNPosi(m_nModelPos);
	// correct the model£¬NnDist is the vector from the modelpoint to unseen imagepiont
	fDist[0] = m_iSgmtPointWorld.x - m_vModelPointsWorld[nMinPos].x;
	fDist[1] = m_iSgmtPointWorld.y - m_vModelPointsWorld[nMinPos].y;
	fDist[2] = m_iSgmtPointWorld.z - m_vModelPointsWorld[nMinPos].z;

	for (int i = 0; i < m_vModelPointsWorld.size(); i++)
	{
		double temp = (double)m_vModelPointsWorld[i].x;
		temp = temp + (double)fDist[0]; // / ((i - m_nModelPos) / 1000 + 1);
		m_vModelPointsWorld[i].x = (float)temp;

		temp = (double)m_vModelPointsWorld[i].y;
		temp = temp + (double)fDist[1]; // / ((i - m_nModelPos) / 1000 + 1);
		m_vModelPointsWorld[i].y =(float)temp;

		temp = (double)m_vModelPointsWorld[i].z;
		temp = temp + (double)fDist[2]; // / ((i - m_nModelPos) / 1000 + 1);
		m_vModelPointsWorld[i].z =(float)temp;
	}
	
	return true;
}
// Function Name: CorrectModelPointBasedSegLength()
//
// Parameters: 
//
// Description: correct the model after Getting Point C;
//
// Returns: 
//
bool miiMinPathModel::CorrectModelPointBasedSegLength()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model of the coronary is null!" << endl;
		return false;
	}
	int nMinPos=0;float fDist[3]={0};
	FindModelPointCbasedOnGeodesicDistance(&nMinPos);
	
	m_nModelPos =nMinPos;
	//CoutNPosi(m_nModelPos);
	fDist[0] = m_iSgmtPointWorld.x - m_vModelPointsWorld[nMinPos].x;
	fDist[1] = m_iSgmtPointWorld.y - m_vModelPointsWorld[nMinPos].y;
	fDist[2] = m_iSgmtPointWorld.z - m_vModelPointsWorld[nMinPos].z;

	for (int i = 0; i < m_vModelPointsWorld.size(); i++)
	{
		double temp = (double)m_vModelPointsWorld[i].x;
		temp = temp + (double)fDist[0]; // / ((i - m_nModelPos) / 1000 + 1);
		m_vModelPointsWorld[i].x = (float)temp;

		temp = (double)m_vModelPointsWorld[i].y;
		temp = temp + (double)fDist[1]; // / ((i - m_nModelPos) / 1000 + 1);
		m_vModelPointsWorld[i].y =(float)temp;

		temp = (double)m_vModelPointsWorld[i].z;
		temp = temp + (double)fDist[2]; // / ((i - m_nModelPos) / 1000 + 1);
		m_vModelPointsWorld[i].z =(float)temp;
	}
	
	return true;

}
bool miiMinPathModel::DFM_CorrModel_WOBP()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model of the coronary is null!" << endl;
		return false;
	}
	int nMinPos=0;float fDist[3]={0};
	//FindModelPointCbasedOnGeodesicDistanceSkipCertainPonts(&nMinPos);
	//FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpoint(&nMinPos);
	DFM_FindC(&nMinPos);
	m_nOneSegModelVectorLastStartPos=m_nOneSegModelVectorStartPos;
	m_nOneSegModelVectorStartPos =nMinPos;

	//fDist[0] = m_iSgmtPointWorld.x - m_vModelPointsWorld[nMinPos].x;
	//fDist[1] = m_iSgmtPointWorld.y - m_vModelPointsWorld[nMinPos].y;
	//fDist[2] = m_iSgmtPointWorld.z - m_vModelPointsWorld[nMinPos].z;

	//for (int i = 0; i < m_vModelPointsWorld.size(); i++)
	//{
	//	double temp = (double)m_vModelPointsWorld[i].x;
	//	temp = temp + (double)fDist[0]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].x = (float)temp;

	//	temp = (double)m_vModelPointsWorld[i].y;
	//	temp = temp + (double)fDist[1]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].y =(float)temp;

	//	temp = (double)m_vModelPointsWorld[i].z;
	//	temp = temp + (double)fDist[2]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].z =(float)temp;
	//}

	return true;


}

// Function Name: FindModelPointCbasedOnGeodesicDistance
//
// Parameters: Position
//
// Description: Find the Point C which has the minimal distance to the tangent plane of segment point.
//
// Returns: 
//
bool miiMinPathModel:: FindModelPointCbasedOnGeodesicDistance(int *Posi)//Find the nearest model point position to the current segment point.
{
	float MDistmm=0;int Nposi=0;
	float fMinDistmm=1000000;

	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{

		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=m_fBackTrackDistmm)
		{ 
			Nposi=i;
			break ;
		} 
	}

	if( Nposi==0 )
	{
		*Posi = m_vModelPointsWorld.size()-1;
	}
	else
	{
		float TanVec[3]={m_fVesselVecWorld[0],m_fVesselVecWorld[1],m_fVesselVecWorld[2]};
		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*Posi = Nposi ;
			return false ; // 
		}
		else
			for (int i =Nposi-DistRangeAroundCPoint/m_meandis_of_coronarymodel; i <=Nposi+DistRangeAroundCPoint/m_meandis_of_coronarymodel; i++) // ¼ÇµÃÐÞ¸Ä5Îª±äÁ¿£¬Í¬Ê±Nposi-5/m_meandis_of_coronarymodel ºÍ Nposi+5/m_meandis_of_coronarymodelÐèÒªÅÐ¶ÏÔÚ0µ½size-1Ö®¼ä£¬Òò´Ë½¨ÒéÄãÐ´¸öº¯ÊýÀ´ÐÞ¸ÄÕâ¸öindex
			{
				if (i<0)i=0;
				if (i>m_vModelPointsWorld.size()-1)i=m_vModelPointsWorld.size()-1;
				float fSgmPointWorld[3]={m_iSgmtPointWorld.x,m_iSgmtPointWorld.y,m_iSgmtPointWorld.z};
				 
				float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);

				if (MPTTDistmm<fMinDistmm)
				{
					fMinDistmm=MPTTDistmm;
					 *Posi = i; 
				} 
			}
	}
	return true ; 
}
bool miiMinPathModel:: DFM_FindC(int *Posi)//Find the nearest model point position to the current segment point.
{
	
	int Pposi=0;
	float fMinNDistmm=1000000;
	//find the minimal position of point C based on the length of the backtrack length.
	float MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=m_fBackTrackDistmmSkipSomePonts)
		{ 
			Pposi=i;
			break ;
		}
	} 
	
	float d=CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	float fSgmPointWorld[3]={m_fSmoothSgmtPointWorld.x,m_fSmoothSgmtPointWorld.y,m_fSmoothSgmtPointWorld.z};
	float TanVec[3]={m_fSmoothVesselVecWorld[0],m_fSmoothVesselVecWorld[1],m_fSmoothVesselVecWorld[2]};
	//float TanVec[3]={m_fSmoothSgmtPointWorld.x-m_iStartPointWorld.x,m_fSmoothSgmtPointWorld.y-m_iStartPointWorld.y,m_fSmoothSgmtPointWorld.z-m_iStartPointWorld.z};
	float fFSgmPointWorld[3]={m_fForwardSmoothSgmtPointWorld.x,m_fForwardSmoothSgmtPointWorld.y,m_fForwardSmoothSgmtPointWorld.z};
	int TPosi=0,NPosi=0;
	if( Pposi==0||CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
	{
 		NPosi = m_vModelPointsWorld.size()-1; 

	}
	else
	{

		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*Posi = Pposi ;
			return false ; // 
		}
		else
		{
			int iStartOfNposi = Pposi-2.0/(m_meandis_of_coronarymodel+0.001); 
			int iEndOfNposi = Pposi+2.0/(m_meandis_of_coronarymodel+0.001); 
			if (iStartOfNposi<1) 
				iStartOfNposi = 1;
			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
				iEndOfNposi = m_vModelPointsWorld.size()-1; 

			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
			{
				float Fcord[3];
				Fcord[0]=m_vModelPointsWorld[i].x;
				Fcord[1]=m_vModelPointsWorld[i].y;
				Fcord[2]=m_vModelPointsWorld[i].z;
				float MPTNDistmm=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
				if (MPTNDistmm<fMinNDistmm)
				{
					fMinNDistmm=MPTNDistmm;
					NPosi = i; 
				} 
			}
		}
	}
	*Posi=NPosi;
	float fMEndDistmm2=CalcDistmm2(m_vModelPointsWorld[*Posi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	if (fMEndDistmm2<MINDISTTOPOINT2)
	{
		cout << "Reached the model end-point!\n" << endl;
		*Posi = m_vModelPointsWorld.size()-1; 
	}
	return true ; 
}

//void miiMinPathModel:: FindModelPointCbasedOnGeodesicDistance(int *Posi)//Find the nearest point (C Point) to the current segment point.
//{
//	float MDistmm=0;int Nposi=0;
//	 float fMinDistmm=1000000;
//	
//	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
// {
// 
// MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
//   if(MDistmm>=m_fBackTrackDistmm)
//   {
//    
//    Nposi=i;
//    break ;
//   } 
//}
// 
//	for (int i =Nposi-12/m_meandis_of_coronarymodel; i <=Nposi+12/m_meandis_of_coronarymodel; i++)
//  {
//	  float TanVec[3]={m_fVesselVecWorld[0],m_fVesselVecWorld[1],m_fVesselVecWorld[2]};
//	  float fSgmPointWorld[3]={m_iSgmtPointWorld.x,m_iSgmtPointWorld.y,m_iSgmtPointWorld.z};
//	  if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
//		  *Posi=m_nLModelPos;
//	  else
//	  {
//	  float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);
//	
//	  if (MPTTDistmm<fMinDistmm)
//		{
//			fMinDistmm=MPTTDistmm;
//			m_nLModelPos=i;
//			 *Posi = i;
//
//		}
//	  }
//  }
//    if( Nposi==0 ) 
//   {
//	   *Posi = m_vModelPointsWorld.size()-1;
//       m_nLModelPos=m_vModelPointsWorld.size()-1;
//	}
//}
// Function Name: ModelPointToTanPlaneDistmm
//
// Parameters: Position ,Tangent direction of the segment point,Segment point
//
// Description: Calculate the distance from the model point to the tangent plane.
//
// Returns:
//
float miiMinPathModel::ModelPointToTanPlaneDistmm(int i,float TanVec[3],float fSgmPointWorld[3])
{
	 if(i<0)i=0;
	 if(i>m_vModelPointsWorld.size()-1)i=m_vModelPointsWorld.size()-1;
	 float NDistmm;
	 float UMToSegPointVec[3]={0,0,0};//UMToSegPointVec means the undermined vector from model point to segment point;
	 UMToSegPointVec[0]=fSgmPointWorld[0]-m_vModelPointsWorld[i].x;
	 UMToSegPointVec[1]=fSgmPointWorld[1]-m_vModelPointsWorld[i].y;
	 UMToSegPointVec[2]=fSgmPointWorld[2]-m_vModelPointsWorld[i].z;
	 float fUNNorm=sqrt(UMToSegPointVec[0]*UMToSegPointVec[0]+UMToSegPointVec[1]*UMToSegPointVec[1]+UMToSegPointVec[2]*UMToSegPointVec[2]);
	 float fm_vVVWNorm=sqrt(TanVec[0]*TanVec[0]+TanVec[1]*TanVec[1]+TanVec[2]*TanVec[2]);
	 NDistmm=abs((UMToSegPointVec[0]*TanVec[0]+UMToSegPointVec[1]*TanVec[1]+UMToSegPointVec[2]*TanVec[2]))/( fm_vVVWNorm);

	 return NDistmm;
}

float miiMinPathModel::ModelPointToTanPlaneDistmm1(int i,float TanVec[3],float fSgmPointWorld[3])//the the point before the segment point is used to calculate the distance
{
	 if(i<0)i=0;
	 if(i>m_vModelPointsWorld.size()-1)i=m_vModelPointsWorld.size()-1;
	 float NDistmm;
	 float UMToSegPointVec[3]={0,0,0};//UMToSegPointVec means the undermined vector from model point to segment point;
	 UMToSegPointVec[0]=m_vModelPointsWorld[i].x-fSgmPointWorld[0];
	 UMToSegPointVec[1]=m_vModelPointsWorld[i].y-fSgmPointWorld[1];
	 UMToSegPointVec[2]=m_vModelPointsWorld[i].z-fSgmPointWorld[2];
	 float fUNNorm=sqrt(UMToSegPointVec[0]*UMToSegPointVec[0]+UMToSegPointVec[1]*UMToSegPointVec[1]+UMToSegPointVec[2]*UMToSegPointVec[2]);
	 float fm_vVVWNorm=sqrt(TanVec[0]*TanVec[0]+TanVec[1]*TanVec[1]+TanVec[2]*TanVec[2]);
	 NDistmm=abs((UMToSegPointVec[0]*TanVec[0]+UMToSegPointVec[1]*TanVec[1]+UMToSegPointVec[2]*TanVec[2]))/( fm_vVVWNorm);
	 return NDistmm;
}

bool miiMinPathModel::CorrectModelPointBasedSegLengthdontMoveModel()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model of the coronary is null!" << endl;
		return false;
	}
	int nMinPos=0;float fDist[3]={0};
	//FindModelPointCbasedOnGeodesicDistanceSkipCertainPonts(&nMinPos);
	//FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpoint(&nMinPos);
	FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpointDc(&nMinPos);
	m_nOneSegModelVectorLastStartPos=m_nOneSegModelVectorStartPos;
	m_nOneSegModelVectorStartPos =nMinPos;

	//fDist[0] = m_iSgmtPointWorld.x - m_vModelPointsWorld[nMinPos].x;
	//fDist[1] = m_iSgmtPointWorld.y - m_vModelPointsWorld[nMinPos].y;
	//fDist[2] = m_iSgmtPointWorld.z - m_vModelPointsWorld[nMinPos].z;

	//for (int i = 0; i < m_vModelPointsWorld.size(); i++)
	//{
	//	double temp = (double)m_vModelPointsWorld[i].x;
	//	temp = temp + (double)fDist[0]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].x = (float)temp;

	//	temp = (double)m_vModelPointsWorld[i].y;
	//	temp = temp + (double)fDist[1]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].y =(float)temp;

	//	temp = (double)m_vModelPointsWorld[i].z;
	//	temp = temp + (double)fDist[2]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].z =(float)temp;
	//}

	return true;

}
bool miiMinPathModel::CorrectModelPointBasedSegLengthdontMoveModelSM()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model of the coronary is null!" << endl;
		return false;
	}
	int nMinPos=0;float fDist[3]={0};
	//FindModelPointCbasedOnGeodesicDistanceSkipCertainPonts(&nMinPos);
	FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpointSM(&nMinPos);
	m_nOneSegModelVectorLastStartPos=m_nOneSegModelVectorStartPos;
	m_nOneSegModelVectorStartPos =nMinPos;
	m_fCurrModelCArchLength = 0 ; 
	for (int i = 1; i <=m_nOneSegModelVectorStartPos ; i++)
	{ 
		m_fCurrModelCArchLength += sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1])); 
	}
	//fDist[0] = m_iSgmtPointWorld.x - m_vModelPointsWorld[nMinPos].x;
	//fDist[1] = m_iSgmtPointWorld.y - m_vModelPointsWorld[nMinPos].y;
	//fDist[2] = m_iSgmtPointWorld.z - m_vModelPointsWorld[nMinPos].z;

	//for (int i = 0; i < m_vModelPointsWorld.size(); i++)
	//{
	//	double temp = (double)m_vModelPointsWorld[i].x;
	//	temp = temp + (double)fDist[0]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].x = (float)temp;

	//	temp = (double)m_vModelPointsWorld[i].y;
	//	temp = temp + (double)fDist[1]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].y =(float)temp;

	//	temp = (double)m_vModelPointsWorld[i].z;
	//	temp = temp + (double)fDist[2]; // / ((i - m_nModelPos) / 1000 + 1);
	//	m_vModelPointsWorld[i].z =(float)temp;
	//}

	return true;

}
// Function Name: FindModelPointCbasedOnGeodesicDistance
//
// Parameters: Position
//
// Description: Find the Point C which has the minimal distance to the tangent plane of segment point.
//
// Returns: 
//
bool miiMinPathModel:: FindModelPointCbasedOnGeodesicDistanceSkipCertainPonts(int *Posi)//Find the nearest model point position to the current segment point.
{
	float MDistmm=0;int Pposi=0;
	float fMinDistmm=1000000;
	//find the point C based on the length of the backtrack length.
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{

		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=m_fBackTrackDistmmSkipSomePonts)
		{ 
			Pposi=i;
			break ;
		} 
	} 
	float d=CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	if( Posi==0||CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
	{
		*Posi = m_vModelPointsWorld.size()-1; 

	}
	else
	{
		float TanVec[3]={m_fVesselVecWorld[0],m_fVesselVecWorld[1],m_fVesselVecWorld[2]};
		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*Posi = Pposi ;
			return false ; // 
		}
		else
		{
			int iStartOfNposi = Pposi-DistRangeAroundCPoint/m_meandis_of_coronarymodel ; 
			int iEndOfNposi = Pposi + DistRangeAroundCPoint/m_meandis_of_coronarymodel; 
			if (iStartOfNposi<0) 
				iStartOfNposi = 0;
			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
				iEndOfNposi = m_vModelPointsWorld.size()-1; 

			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
			{

				float fSgmPointWorld[3]={m_iSgmtPointWorld.x,m_iSgmtPointWorld.y,m_iSgmtPointWorld.z};
				float fFSgmPointWorld[3]={m_fForwardSgmtPointWorld.x,m_fForwardSgmtPointWorld.y,m_fForwardSgmtPointWorld.z};
				float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);

				if (MPTTDistmm<fMinDistmm)
				{
					fMinDistmm=MPTTDistmm;
					*Posi = i; 
				} 
			}
		}
	}

	return true ; 
}

	bool miiMinPathModel:: FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpoint(int *Posi)//Find the nearest model point position to the current segment point.
{
	
	int PMinPosi=m_vModelPointsWorld.size()-1,PMaxPos=m_vModelPointsWorld.size()-1,Pposi=0;
	float fMinTDistmm=1000000,fMinNDistmm=1000000;
	//find the minimal position of point C based on the length of the backtrack length.
	float MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=mZXH_MinPencentageForCPointLength*m_fBackTrackDistmmSkipSomePonts)
		{ 
			PMinPosi=i;
			break ;
		}
	} 
	MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=mZXH_MaxPencentageForCPointLength*m_fBackTrackDistmmSkipSomePonts)
		{ 
			PMaxPos=i;
			break ;
		}
	} 
	Pposi=(int)(PMinPosi+PMaxPos)/2;
	float d=CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	float fSgmPointWorld[3]={m_fSmoothSgmtPointWorld.x,m_fSmoothSgmtPointWorld.y,m_fSmoothSgmtPointWorld.z};
	float TanVec[3]={m_fSmoothVesselVecWorld[0],m_fSmoothVesselVecWorld[1],m_fSmoothVesselVecWorld[2]};
	//float TanVec[3]={m_fSmoothSgmtPointWorld.x-m_iStartPointWorld.x,m_fSmoothSgmtPointWorld.y-m_iStartPointWorld.y,m_fSmoothSgmtPointWorld.z-m_iStartPointWorld.z};
	float fFSgmPointWorld[3]={m_fForwardSmoothSgmtPointWorld.x,m_fForwardSmoothSgmtPointWorld.y,m_fForwardSmoothSgmtPointWorld.z};
	int TPosi=0,NPosi=0;
	if( Pposi==0||CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
	{
		*Posi = m_vModelPointsWorld.size()-1; 

	}
	else
	{

		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*Posi = Pposi ;
			return false ; // 
		}
		else
		{
			int iStartOfNposi = PMinPosi-DistRangeAroundCPoint/m_meandis_of_coronarymodel ; 
			int iEndOfNposi = PMaxPos + DistRangeAroundCPoint/m_meandis_of_coronarymodel; 
			if (iStartOfNposi<1) 
				iStartOfNposi = 1;
			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
				iEndOfNposi = m_vModelPointsWorld.size()-1; 

			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
			{


				float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);
				float Fcord[3];
				Fcord[0]=m_vModelPointsWorld[i].x;
				Fcord[1]=m_vModelPointsWorld[i].y;
				Fcord[2]=m_vModelPointsWorld[i].z;
				float MPTNDistmm=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
				if (MPTTDistmm<fMinTDistmm)
				{
					fMinTDistmm=MPTTDistmm;
					TPosi = i; 
				} 
				if (MPTNDistmm<fMinNDistmm)
				{
					fMinNDistmm=MPTNDistmm;
					NPosi = i; 
				} 
			}
		}
	}
	//float Fcord[3];
	//Fcord[0]=m_vModelPointsWorld[TPosi].x;
	//Fcord[1]=m_vModelPointsWorld[TPosi].y;
	//Fcord[2]=m_vModelPointsWorld[TPosi].z;
	//float fTDist=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
	//if (fMinTDistmm<1&&fTDist<fMinNDistmm)
	//{
	//		*Posi=TPosi;
	//	
	//}
	//else
	*Posi=NPosi;
	float fMEndDistmm2=CalcDistmm2(m_vModelPointsWorld[*Posi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	if (fMEndDistmm2<MINDISTTOPOINT2)
	{
		*Posi = m_vModelPointsWorld.size()-1; 
	}
	/*if(fMinTDistmm<1&&fMinNDistmm<fTDist)
		*Posi=NPosi;
	else
		*Posi=TPosi;*/
	/*float vesseltangent[3]={0,0,0};
	float modelvector[3]={0,0,0};
	float modeltosegvector[3]={0,0,0};
	float modelctolasdVec[3]={0,0,0};
	float modelctolascVec[3]={0,0,0};
	modelvector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	modelvector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	modelvector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	modeltosegvector[0]=fSgmPointWorld[0]-m_vModelPointsWorld[*Posi].x;
	modeltosegvector[1]=fSgmPointWorld[1]-m_vModelPointsWorld[*Posi].y;
	modeltosegvector[2]=fSgmPointWorld[2]-m_vModelPointsWorld[*Posi].z;
	modelctolasdVec[0]=m_vModelPointsWorld[*Posi].x-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x;
	modelctolasdVec[1]=m_vModelPointsWorld[*Posi].y-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y;
	modelctolasdVec[2]=m_vModelPointsWorld[*Posi].z-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z;
	modelctolascVec[0]=m_vModelPointsWorld[*Posi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	modelctolascVec[1]=m_vModelPointsWorld[*Posi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	modelctolascVec[2]=m_vModelPointsWorld[*Posi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	float modelvectordist=sqrt(zxh::VectorOP_DotProduct(modelvector,modelvector,3));
	float modeltosegvectordist=sqrt(zxh::VectorOP_DotProduct(modeltosegvector,modeltosegvector,3));
	float modelctolasdVeddist=sqrt(zxh::VectorOP_DotProduct(modelctolasdVec,modelctolasdVec,3));
	float modelctolasdVecdist=sqrt(zxh::VectorOP_DotProduct(modelctolascVec,modelctolascVec,3));*/
	/*float costheta=zxh::VectorOP_Cosine( modelvector,TanVec, 3);
	float dist1=sqrt(zxh::VectorOP_DotProduct(modelvector,modelvector,3));
	float dist2=costheta*dist1;*/
	return true ; 
}
bool miiMinPathModel:: FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpointDc(int *Posi)//Find the nearest model point position to the current segment point.
{
	
	int Pposi=0;
	float fMinNDistmm=1000000;
	//find the minimal position of point C based on the length of the backtrack length.
	float MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=m_fBackTrackDistmmSkipSomePonts)
		{ 
			Pposi=i;
			break ;
		}
	} 
	
	float d=CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	float fSgmPointWorld[3]={m_fSmoothSgmtPointWorld.x,m_fSmoothSgmtPointWorld.y,m_fSmoothSgmtPointWorld.z};
	float TanVec[3]={m_fSmoothVesselVecWorld[0],m_fSmoothVesselVecWorld[1],m_fSmoothVesselVecWorld[2]};
	//float TanVec[3]={m_fSmoothSgmtPointWorld.x-m_iStartPointWorld.x,m_fSmoothSgmtPointWorld.y-m_iStartPointWorld.y,m_fSmoothSgmtPointWorld.z-m_iStartPointWorld.z};
	float fFSgmPointWorld[3]={m_fForwardSmoothSgmtPointWorld.x,m_fForwardSmoothSgmtPointWorld.y,m_fForwardSmoothSgmtPointWorld.z};
	int TPosi=0,NPosi=0;
	if( Pposi==0||CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
	{
 		NPosi = m_vModelPointsWorld.size()-1; 

	}
	else
	{

		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*Posi = Pposi ;
			return false ; // 
		}
		else
		{
			int iStartOfNposi = Pposi-2.0/(m_meandis_of_coronarymodel+0.001); 
			int iEndOfNposi = Pposi+2.0/(m_meandis_of_coronarymodel+0.001); 
			if (iStartOfNposi<1) 
				iStartOfNposi = 1;
			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
				iEndOfNposi = m_vModelPointsWorld.size()-1; 

			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
			{
				float Fcord[3];
				Fcord[0]=m_vModelPointsWorld[i].x;
				Fcord[1]=m_vModelPointsWorld[i].y;
				Fcord[2]=m_vModelPointsWorld[i].z;
				float MPTNDistmm=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
				if (MPTNDistmm<fMinNDistmm)
				{
					fMinNDistmm=MPTNDistmm;
					NPosi = i; 
				} 
			}
		}
	}
	*Posi=NPosi;
	float fMEndDistmm2=CalcDistmm2(m_vModelPointsWorld[*Posi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	if (fMEndDistmm2<MINDISTTOPOINT2)
	{
		cout << "Reached the model end-point!\n" << endl;
		*Posi = m_vModelPointsWorld.size()-1; 
	}
	return true ; 
}
bool miiMinPathModel:: FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpointSM(int *Posi)//Find the nearest model point position to the current segment point.
{
	
	int PMinPosi=m_vModelPointsWorld.size()-1,PMaxPos=m_vModelPointsWorld.size()-1,Pposi=0;
	float fMinTDistmm=1000000,fMinNDistmm=1000000;
	//find the minimal position of point C based on the length of the backtrack length.
	float MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=mZXH_MinPencentageForCPointLength*m_fBackTrackDistmmSkipSomePonts)
		{ 
			PMinPosi=i;
			break ;
		}
	} 
	MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=mZXH_MaxPencentageForCPointLength*m_fBackTrackDistmmSkipSomePonts)
		{ 
			PMaxPos=i;
			break ;
		}
	} 
	Pposi=(int)(PMinPosi+PMaxPos)/2;
	float d=CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	float fSgmPointWorld[3]={m_fSmoothSgmtPointWorld.x,m_fSmoothSgmtPointWorld.y,m_fSmoothSgmtPointWorld.z};
	float TanVec[3]={m_fSmoothVesselVecWorld[0],m_fSmoothVesselVecWorld[1],m_fSmoothVesselVecWorld[2]};
	//float TanVec[3]={m_fSmoothSgmtPointWorld.x-m_iStartPointWorld.x,m_fSmoothSgmtPointWorld.y-m_iStartPointWorld.y,m_fSmoothSgmtPointWorld.z-m_iStartPointWorld.z};
	float fFSgmPointWorld[3]={m_fForwardSmoothSgmtPointWorld.x,m_fForwardSmoothSgmtPointWorld.y,m_fForwardSmoothSgmtPointWorld.z};
	int TPosi=0,NPosi=0;
	if( Pposi==0||CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
	{
		*Posi = m_vModelPointsWorld.size()-1; 

	}
	else
	{

		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*Posi = Pposi ;
			return false ; // 
		}
		else
		{
			int iStartOfNposi = PMinPosi-DistRangeAroundCPoint/m_meandis_of_coronarymodel ; 
			int iEndOfNposi = PMaxPos + DistRangeAroundCPoint/m_meandis_of_coronarymodel; 
			if (iStartOfNposi<1) 
				iStartOfNposi = 1;
			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
				iEndOfNposi = m_vModelPointsWorld.size()-1; 

			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
			{


				float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);
				float Fcord[3];
				Fcord[0]=m_vModelPointsWorld[i].x;
				Fcord[1]=m_vModelPointsWorld[i].y;
				Fcord[2]=m_vModelPointsWorld[i].z;
				float MPTNDistmm=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
				if (MPTTDistmm<fMinTDistmm)
				{
					fMinTDistmm=MPTTDistmm;
					TPosi = i; 
				} 
				if (MPTNDistmm<fMinNDistmm)
				{
					fMinNDistmm=MPTNDistmm;
					NPosi = i; 
				} 
			}
		}
	}
	//float Fcord[3];
	//Fcord[0]=m_vModelPointsWorld[TPosi].x;
	//Fcord[1]=m_vModelPointsWorld[TPosi].y;
	//Fcord[2]=m_vModelPointsWorld[TPosi].z;
	//float fTDist=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
	//if (fMinTDistmm<1&&fTDist<fMinNDistmm)
	//{
	//		*Posi=TPosi;
	//	
	//}
	//else
			*Posi=NPosi;
float fMEndDistmm2=CalcDistmm2(m_vModelPointsWorld[*Posi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
if (fMEndDistmm2<MINDISTTOPOINT2)
{
	*Posi = m_vModelPointsWorld.size()-1; 
}

	/*if(fMinTDistmm<1&&fMinNDistmm<fTDist)
		*Posi=NPosi;
	else
		*Posi=TPosi;*/
	/*float vesseltangent[3]={0,0,0};
	float modelvector[3]={0,0,0};
	float modeltosegvector[3]={0,0,0};
	float modelctolasdVec[3]={0,0,0};
	float modelctolascVec[3]={0,0,0};
	modelvector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	modelvector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	modelvector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	modeltosegvector[0]=fSgmPointWorld[0]-m_vModelPointsWorld[*Posi].x;
	modeltosegvector[1]=fSgmPointWorld[1]-m_vModelPointsWorld[*Posi].y;
	modeltosegvector[2]=fSgmPointWorld[2]-m_vModelPointsWorld[*Posi].z;
	modelctolasdVec[0]=m_vModelPointsWorld[*Posi].x-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x;
	modelctolasdVec[1]=m_vModelPointsWorld[*Posi].y-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y;
	modelctolasdVec[2]=m_vModelPointsWorld[*Posi].z-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z;
	modelctolascVec[0]=m_vModelPointsWorld[*Posi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	modelctolascVec[1]=m_vModelPointsWorld[*Posi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	modelctolascVec[2]=m_vModelPointsWorld[*Posi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	float modelvectordist=sqrt(zxh::VectorOP_DotProduct(modelvector,modelvector,3));
	float modeltosegvectordist=sqrt(zxh::VectorOP_DotProduct(modeltosegvector,modeltosegvector,3));
	float modelctolasdVeddist=sqrt(zxh::VectorOP_DotProduct(modelctolasdVec,modelctolasdVec,3));
	float modelctolasdVecdist=sqrt(zxh::VectorOP_DotProduct(modelctolascVec,modelctolascVec,3));*/
	/*float costheta=zxh::VectorOP_Cosine( modelvector,TanVec, 3);
	float dist1=sqrt(zxh::VectorOP_DotProduct(modelvector,modelvector,3));
	float dist2=costheta*dist1;*/
	return true ; 
}
// Function Name: UpdateModelSgmtPoint()
//
// Parameters: 
//
// Description: get the position of the segmented point in the model, 
//		update the value of "m_nModelPos"
//
// Returns: 
//
bool miiMinPathModel::UpdateModelSgmtPoint()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model is null!" << endl;
		return false;
	}

	/*if (m_nModelPos >= m_vModelPointsWorld.size() - 1)
	{
	return false;
	}*/

	// update and save start point	all points are seen on the unseen image
	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iStartPoint.x = m_iSgmtPoint.x;
	m_iStartPoint.y = m_iSgmtPoint.y;
	m_iStartPoint.z = m_iSgmtPoint.z;

	m_iStartPointWorld.x = m_iSgmtPointWorld.x;
	m_iStartPointWorld.y = m_iSgmtPointWorld.y;
	m_iStartPointWorld.z = m_iSgmtPointWorld.z;
	// find the model position corresponding to the segmented point of FMM
	int i=0;
	int Pixdist=(int)(iMinDisBetweenVec/m_meandis_of_coronarymodel+0.5);
	
	bool bPosFound = false ; 
	if (m_nOneSegModelVectorStartPos==0)
	{
		m_nOneSegModelVectorEndPos=1;
		bPosFound =  true;
	}
	else if(m_nOneSegModelVectorStartPos==m_vModelPointsWorld.size()-1)//Point C is the end point
	{

		for (i = m_nOneSegModelVectorStartPos - Pixdist; i >0 && bPosFound==false ; i--)//Add by JDQ
		{
			if( GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{
				m_nOneSegModelVectorEndPos = i;
				//Add by JDQ
				//float n=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_iStartPointWorld);
				bPosFound =  true;

			}
		}
	}
	else
	{


		for (i = m_nOneSegModelVectorStartPos + Pixdist; i < m_vModelPointsWorld.size() && bPosFound==false ; i++)//Add by JDQ
		{
			float m=GetModPointCos(i);
			float l=GetModVectorCos(i);
			if( GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{

				//Add by JDQ
				float n=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[i]);
				bPosFound =  true;
				float m=GetModPointCos(i);
				float l=GetModVectorCos(i);
				break;

			}
		}
		m_nOneSegModelVectorEndPos = i;

	}
	if (i == m_vModelPointsWorld.size())
	{
		m_nOneSegModelVectorEndPos = i-1;
		return true;
	}


	return bPosFound ;
}
// Function Name: UpdateModelSgmtPointNew()
//
// Parameters: 
//
// Description: get the position of the segmented point in the model, 
//		update the value of "m_nModelPos"
//
// Returns: 
//
bool miiMinPathModel::UpdateModelSgmtPointNew()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model is null!" << endl;
		return false;
	}

	/*if (m_nModelPos >= m_vModelPointsWorld.size() - 1)
	{
	return false;
	}*/

	// update and save start point	all points are seen on the unseen image
	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iStartPoint.x = m_iSgmtPoint.x;
	m_iStartPoint.y = m_iSgmtPoint.y;
	m_iStartPoint.z = m_iSgmtPoint.z;

	m_iStartPointWorld.x = m_iSgmtPointWorld.x;
	m_iStartPointWorld.y = m_iSgmtPointWorld.y;
	m_iStartPointWorld.z = m_iSgmtPointWorld.z;
	// find the model position corresponding to the segmented point of FMM
	int i=0;
	int Pixdist=(int)(iMinDisBetweenVec/m_meandis_of_coronarymodel+0.5);
	int PixdistMax=(int)(iMaxDisBetweenVec/m_meandis_of_coronarymodel+0.5);
	bool bPosFound = false ; 
	if (m_nOneSegModelVectorStartPos==0)
	{
		m_nOneSegModelVectorEndPos=1;
		bPosFound =  true;
	}
	else if(m_nOneSegModelVectorStartPos==m_vModelPointsWorld.size()-1)//Point C is the end point
	{

		for (i = m_nOneSegModelVectorStartPos - Pixdist; i >0 && bPosFound==false ; i--)//Add by JDQ
		{

			float fdist=sqrt(CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[i]));
			if( fdist>iMaxDisBetweenVec||GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{
				m_nOneSegModelVectorEndPos = i;
				//Add by JDQ
				//float n=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_iStartPointWorld);
				bPosFound =  true;
				return bPosFound ;

			}
		}
	}
	else
	{


		for (i = m_nOneSegModelVectorStartPos + Pixdist; i < m_vModelPointsWorld.size() && bPosFound==false ; i++)//Add by JDQ
		{
			float m=GetModPointCos(i);
			float l=GetModVectorCos(i);
			float fdist=sqrt(CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[i]));
			if( fdist>iMaxDisBetweenVec||GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{

				//Add by JDQ
				m_nOneSegModelVectorEndPos=i;
				bPosFound =  true;
				float m=GetModPointCos(i);
				float l=GetModVectorCos(i);
				return bPosFound ;

			}


		}

		m_nOneSegModelVectorEndPos =m_vModelPointsWorld.size()-1;//cannot find Point D ,set point D is the end point of the model

		bPosFound =  true;
		return bPosFound ;

	}
	

	return bPosFound ;
}
bool miiMinPathModel::UpdateModelSgmtPointNewSM()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model is null!" << endl;
		return false;
	}

	/*if (m_nModelPos >= m_vModelPointsWorld.size() - 1)
	{
	return false;
	}*/

	// update and save start point	all points are seen on the unseen image
	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iStartPoint.x = m_iSgmtPoint.x;
	m_iStartPoint.y = m_iSgmtPoint.y;
	m_iStartPoint.z = m_iSgmtPoint.z;

	m_iStartPointWorld.x = m_iSgmtPointWorld.x;
	m_iStartPointWorld.y = m_iSgmtPointWorld.y;
	m_iStartPointWorld.z = m_iSgmtPointWorld.z;
	// find the model position corresponding to the segmented point of FMM
	int i=0;
	int Pixdist=(int)(iMinDisBetweenVec/m_meandis_of_coronarymodel+0.5);
	int PixdistMax=(int)(iMaxDisBetweenVec/m_meandis_of_coronarymodel+0.5);
	bool bPosFound = false ; 
	if (m_nOneSegModelVectorStartPos==0)
	{
		m_nOneSegModelVectorEndPos=1;
		bPosFound =  true;
	}
	else if(m_nOneSegModelVectorStartPos==m_vModelPointsWorld.size()-1)//Point C is the end point
	{

		for (i = m_nOneSegModelVectorStartPos - Pixdist; i >0 && bPosFound==false ; i--)//Add by JDQ
		{

			float fdist=sqrt(CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[i]));
			if( fdist>iMaxDisBetweenVec||GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{
				m_nOneSegModelVectorEndPos = i;
				//Add by JDQ
				//float n=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_iStartPointWorld);
				bPosFound =  true;

			}
		}
	}
	else
	{


		for (i = m_nOneSegModelVectorStartPos + Pixdist; i < m_vModelPointsWorld.size() && bPosFound==false ; i++)//Add by JDQ
		{
			float m=GetModPointCos(i);
			float l=GetModVectorCos(i);
			float fdist=sqrt(CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[i]));
			if( fdist>iMaxDisBetweenVec||GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{

				//Add by JDQ
				m_nOneSegModelVectorEndPos=i;
				bPosFound =  true;
				float m=GetModPointCos(i);
				float l=GetModVectorCos(i);

				return bPosFound ;

			}


		}

		m_nOneSegModelVectorEndPos =m_vModelPointsWorld.size()-1;
		bPosFound =  true;
		return bPosFound ;


	}

	 return bPosFound ;
}
bool miiMinPathModel::UpdateModelSgmtPointCD()
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model is null!" << endl;
		return false;
	}

	/*if (m_nModelPos >= m_vModelPointsWorld.size() - 1)
	{
	return false;
	}*/

	// update and save start point	all points are seen on the unseen image
	m_iLastStartPoint.x = m_iStartPoint.x;
	m_iLastStartPoint.y = m_iStartPoint.y;
	m_iLastStartPoint.z = m_iStartPoint.z;

	m_iLastStartPointWorld.x = m_iStartPointWorld.x;
	m_iLastStartPointWorld.y = m_iStartPointWorld.y;
	m_iLastStartPointWorld.z = m_iStartPointWorld.z;

	m_iStartPoint.x = m_iSgmtPoint.x;
	m_iStartPoint.y = m_iSgmtPoint.y;
	m_iStartPoint.z = m_iSgmtPoint.z;

	m_iStartPointWorld.x = m_iSgmtPointWorld.x;
	m_iStartPointWorld.y = m_iSgmtPointWorld.y;
	m_iStartPointWorld.z = m_iSgmtPointWorld.z;
	// find the model position corresponding to the segmented point of FMM
	int i=0;
	int Pixdist=(int)(iMinDisBetweenVec/m_meandis_of_coronarymodel+0.5);
	int PixdistMax=(int)(iMaxDisBetweenVec/m_meandis_of_coronarymodel+0.5);
	bool bPosFound = false ; 
	if (m_nOneSegModelVectorStartPos==0)
	{
		m_nOneSegModelVectorEndPos=1;
		bPosFound =  true;
	}
	else if(m_nOneSegModelVectorStartPos==m_vModelPointsWorld.size()-1)//Point C is the end point
	{

		for (i = m_nOneSegModelVectorStartPos - Pixdist; i >0 && bPosFound==false ; i--)//Add by JDQ
		{

			float fdist=sqrt(CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[i]));
			if( fdist>iMaxDisBetweenVec||GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{
				m_nOneSegModelVectorEndPos = i;
				//Add by JDQ
				//float n=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_iStartPointWorld);
				bPosFound =  true;

			}
		}
	}

	else
	{


		for (i = m_nOneSegModelVectorStartPos + Pixdist; i < m_vModelPointsWorld.size() && bPosFound==false ; i++)//Add by JDQ
		{
			float m=GetModPointCos(i);
			float l=GetModVectorCos(i);
			float fdist=sqrt(CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[i]));
			if( fdist>iMaxDisBetweenVec||GetModPointCos(i)<m_nKPCosTheta||GetModVectorCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
				// if( GetModPointCos(i)<m_nKPCosTheta)//theta at current point <30¡ãthen seg when theta >30¡ã
			{

				//Add by JDQ
				
				bPosFound =  true;
				float m=GetModPointCos(i);
				float l=GetModVectorCos(i);
				break;

			}
             
		}
		m_nOneSegModelVectorEndPos = m_nOneSegModelVectorStartPos+Pixdist;
		bPosFound =  true;

	}
	if (i == m_vModelPointsWorld.size())
	{
		m_nOneSegModelVectorEndPos = i-1;
		return true;
	}


	return bPosFound ;
}


// Function Name: CalcTheata()
//
// Parameters: cx: the coordinate x
//			   cy: the coordinate y
//			   cz: the coordinate z
//			   bInitFlg: the flag for FMM evolution in aorta
//
// Description: 
//
// Returns: 
//
float miiMinPathModel::CalcTheata(int cx, int cy, int cz, bool bInitFlg)
{
	float vModel[3], vVessel[3];

	if (bInitFlg)//caculate the vector first
	{
		vModel[0] = m_vModelPointsWorld[1].x - m_vModelPointsWorld[0].x;
		vModel[1] = m_vModelPointsWorld[1].y - m_vModelPointsWorld[0].y;
		vModel[2] = m_vModelPointsWorld[1].z - m_vModelPointsWorld[0].z;
	}
	else//caculate the vector from the second segment.
	{
		/*vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x - m_iStartPointWorld.x;
		vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y - m_iStartPointWorld.y;
		vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z - m_iStartPointWorld.z;*/
		vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
		vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
		vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	}
	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	vVessel[0] = fCord[0] - m_iStartPointWorld.x;
	vVessel[1] = fCord[1] - m_iStartPointWorld.y;
	vVessel[2] = fCord[2] - m_iStartPointWorld.z;
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fModelNorm = sqrt((double)(vModel[0]*vModel[0] + vModel[1]*vModel[1] + vModel[2]*vModel[2]));
	double fVesselNorm = sqrt((double)(vVessel[0]*vVessel[0] + vVessel[1]*vVessel[1] + vVessel[2]*vVessel[2]));

	float theata = (double)(vModel[0] * vVessel[0] + vModel[1] * vVessel[1] + vModel[2] * vVessel[2]) \
		/ ( fModelNorm * fVesselNorm);//"theata" here means  cos(theata)

	return theata;
}
// Function Name: CalcThetaMV()
//
// Parameters: cx: the coordinate x
//			   cy: the coordinate y
//			   cz: the coordinate z
//			   currentvec: the current vector
//
// Description: 
//
// Returns: 
//
float miiMinPathModel::CalcThetaMV(int cx, int cy, int cz,float currentvec[3])
{
	float vModel[3], vVessel[3];

	//caculate the vector from the second segment.
	
		vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
		vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
		vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	miiCNode<double, int> nEndPos,nStartPos;
	nEndPos=WorldTransToImagePoint(m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
	nStartPos=WorldTransToImagePoint(m_vModelPointsWorld[m_nOneSegModelVectorStartPos]);
	//float fCord[3];
	//fCord[0] = (float)cx;
	//fCord[1] = (float)cy;
	//fCord[2] = (float)cz;
	//m_pBaseImgInfo->ImageToWorld(fCord);
	//vVessel[0] = fCord[0] - m_iStartPointWorld.x;
	//vVessel[1] = fCord[1] - m_iStartPointWorld.y;
	//vVessel[2] = fCord[2] - m_iStartPointWorld.z;
	//set vector u1 as vessel vetor
	vVessel[0] =currentvec[0];
	vVessel[1] =currentvec[1];
	vVessel[2] =currentvec[2];
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fModelNorm = sqrt((double)(vModel[0]*vModel[0] + vModel[1]*vModel[1] + vModel[2]*vModel[2]));
	double fVesselNorm = sqrt((double)(vVessel[0]*vVessel[0] + vVessel[1]*vVessel[1] + vVessel[2]*vVessel[2]));

	float theata = (double)(vModel[0] * vVessel[0] + vModel[1] * vVessel[1] + vModel[2] * vVessel[2]) \
		/ ( fModelNorm * fVesselNorm+0.0000001);//"theata" here means  cos(theata)
	return theata;
}

// Function Name: CalcThetaMV_DMP()
//
// Parameters: cx: the coordinate x
//			   cy: the coordinate y
//			   cz: the coordinate z
//			   currentvec: the current vector
//
// Description: _DMP 
//
// Returns: 
//
float miiMinPathModel::CalcThetaMV_DMP(int cx, int cy, int cz,float currentvec[3])
{
	float vModel[3], vVessel[3];

	//caculate the vector from the second segment.
	
		vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
		vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
		vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	miiCNode<double, int> nEndPos,nStartPos;
	nEndPos=WorldTransToImagePoint(m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
	nStartPos=WorldTransToImagePoint(m_vModelPointsWorld[m_nOneSegModelVectorStartPos]);
	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	vVessel[0] = fCord[0] - m_iStartPointWorld.x;
	vVessel[1] = fCord[1] - m_iStartPointWorld.y;
	vVessel[2] = fCord[2] - m_iStartPointWorld.z;
	//set vector u1 as vessel vetor
	/*vVessel[0] =currentvec[0];
	vVessel[1] =currentvec[1];
	vVessel[2] =currentvec[2];*/
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fModelNorm = sqrt((double)(vModel[0]*vModel[0] + vModel[1]*vModel[1] + vModel[2]*vModel[2]));
	double fVesselNorm = sqrt((double)(vVessel[0]*vVessel[0] + vVessel[1]*vVessel[1] + vVessel[2]*vVessel[2]));

	float theata = (double)(vModel[0] * vVessel[0] + vModel[1] * vVessel[1] + vModel[2] * vVessel[2]) \
		/ ( fModelNorm * fVesselNorm+0.0000001);//"theata" here means  cos(theata)
	return theata;
}
// Function Name: CalcThetaMLV()
//
// Parameters: cx: the coordinate x
//			   cy: the coordinate y
//			   cz: the coordinate z
//			   currentvec: the current vector
//
// Description: calculate the cosine value of theta between current vector(from u1) and last vessel vector (estimated value)
//
// Returns: 
//
float miiMinPathModel::CalcThetaMLV(int cx, int cy, int cz,float currentvec[3])
{
	float gLastVesselDirect[3], gCurrentVessel[3];
	gLastVesselDirect[0] = m_iStartPointWorld.x - m_iLastStartPointWorld.x;
	gLastVesselDirect[1] = m_iStartPointWorld.y - m_iLastStartPointWorld.y;
	gLastVesselDirect[2] = m_iStartPointWorld.z - m_iLastStartPointWorld.z;

	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	/*gCurrentVessel[0] = fCord[0] - m_iStartPointWorld.x;
	gCurrentVessel[1] = fCord[1] - m_iStartPointWorld.y;
	gCurrentVessel[2] = fCord[2] - m_iStartPointWorld.z;*/
	//set u1 as vessel vector
	gCurrentVessel[0] =currentvec[0];
	gCurrentVessel[1] =currentvec[1];
	gCurrentVessel[2] =currentvec[2];
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fLastVesDirectNorm = sqrt((double)(gLastVesselDirect[0]*gLastVesselDirect[0] + gLastVesselDirect[1]*gLastVesselDirect[1] + gLastVesselDirect[2]*gLastVesselDirect[2]));
	double fCurrVesDirectNorm = sqrt((double)(gCurrentVessel[0]*gCurrentVessel[0] + gCurrentVessel[1]*gCurrentVessel[1] + gCurrentVessel[2]*gCurrentVessel[2]));

	double theata = (double)(gLastVesselDirect[0] * gCurrentVessel[0] + gLastVesselDirect[1] * gCurrentVessel[1] + gLastVesselDirect[2] * gCurrentVessel[2])\
		/ ( fLastVesDirectNorm * fCurrVesDirectNorm+0.0000001);//"theata" here means  cos(theata)

	return (float)theata;
}
float miiMinPathModel::CalcThetaMLV_DMP(int cx, int cy, int cz,float currentvec[3])
{
	float gLastVesselDirect[3], gCurrentVessel[3];
	gLastVesselDirect[0] = m_iStartPointWorld.x - m_iLastStartPointWorld.x;
	gLastVesselDirect[1] = m_iStartPointWorld.y - m_iLastStartPointWorld.y;
	gLastVesselDirect[2] = m_iStartPointWorld.z - m_iLastStartPointWorld.z;

	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	gCurrentVessel[0] = fCord[0] - m_iStartPointWorld.x;
	gCurrentVessel[1] = fCord[1] - m_iStartPointWorld.y;
	gCurrentVessel[2] = fCord[2] - m_iStartPointWorld.z;
	//set u1 as vessel vector
	/*gCurrentVessel[0] =currentvec[0];
	gCurrentVessel[1] =currentvec[1];
	gCurrentVessel[2] =currentvec[2];*/
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fLastVesDirectNorm = sqrt((double)(gLastVesselDirect[0]*gLastVesselDirect[0] + gLastVesselDirect[1]*gLastVesselDirect[1] + gLastVesselDirect[2]*gLastVesselDirect[2]));
	double fCurrVesDirectNorm = sqrt((double)(gCurrentVessel[0]*gCurrentVessel[0] + gCurrentVessel[1]*gCurrentVessel[1] + gCurrentVessel[2]*gCurrentVessel[2]));

	double theata = (double)(gLastVesselDirect[0] * gCurrentVessel[0] + gLastVesselDirect[1] * gCurrentVessel[1] + gLastVesselDirect[2] * gCurrentVessel[2])\
		/ ( fLastVesDirectNorm * fCurrVesDirectNorm+0.0000001);//"theata" here means  cos(theata)

	return (float)theata;
}
float miiMinPathModel::CalcTheataN(int cx, int cy, int cz,float currentvec[3],bool bInitFlg)
{
	float vModel[3], vVessel[3];

	if (bInitFlg)//caculate the vector first
	{
		vModel[0] = m_vModelPointsWorld[1].x - m_vModelPointsWorld[0].x;
		vModel[1] = m_vModelPointsWorld[1].y - m_vModelPointsWorld[0].y;
		vModel[2] = m_vModelPointsWorld[1].z - m_vModelPointsWorld[0].z;
	}
	else//caculate the vector from the second segment.
	{
		vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
		vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
		vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	}
	//float fCord[3];
	//fCord[0] = (float)cx;
	//fCord[1] = (float)cy;
	//fCord[2] = (float)cz;
	//m_pBaseImgInfo->ImageToWorld(fCord);
	//vVessel[0] = fCord[0] - m_iStartPointWorld.x;
	//vVessel[1] = fCord[1] - m_iStartPointWorld.y;
	//vVessel[2] = fCord[2] - m_iStartPointWorld.z;
	//set vector u1 as vessel vetor
	vVessel[0] =currentvec[0];
	vVessel[1] =currentvec[1];
	vVessel[2] =currentvec[2];
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fModelNorm = sqrt((double)(vModel[0]*vModel[0] + vModel[1]*vModel[1] + vModel[2]*vModel[2]));
	double fVesselNorm = sqrt((double)(vVessel[0]*vVessel[0] + vVessel[1]*vVessel[1] + vVessel[2]*vVessel[2]));

	float theata = (double)(vModel[0] * vVessel[0] + vModel[1] * vVessel[1] + vModel[2] * vVessel[2]) \
		/ ( fModelNorm * fVesselNorm+0.0000001);//"theata" here means  cos(theata)

	return theata;
}
float miiMinPathModel::CalcSegTheta(int cx, int cy, int cz)
{
	float gLastVesselDirect[3], gCurrentVessel[3];
	gLastVesselDirect[0] = m_iStartPointWorld.x - m_iLastStartPointWorld.x;
	gLastVesselDirect[1] = m_iStartPointWorld.y - m_iLastStartPointWorld.y;
	gLastVesselDirect[2] = m_iStartPointWorld.z - m_iLastStartPointWorld.z;

	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	gCurrentVessel[0] = fCord[0] - m_iStartPointWorld.x;
	gCurrentVessel[1] = fCord[1] - m_iStartPointWorld.y;
	gCurrentVessel[2] = fCord[2] - m_iStartPointWorld.z;
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fLastVesDirectNorm = sqrt((double)(gLastVesselDirect[0]*gLastVesselDirect[0] + gLastVesselDirect[1]*gLastVesselDirect[1] + gLastVesselDirect[2]*gLastVesselDirect[2]));
	double fCurrVesDirectNorm = sqrt((double)(gCurrentVessel[0]*gCurrentVessel[0] + gCurrentVessel[1]*gCurrentVessel[1] + gCurrentVessel[2]*gCurrentVessel[2]));

	double theata = (double)(gLastVesselDirect[0] * gCurrentVessel[0] + gLastVesselDirect[1] * gCurrentVessel[1] + gLastVesselDirect[2] * gCurrentVessel[2]) \
		/ ( fLastVesDirectNorm * fCurrVesDirectNorm);//"theata" here means  cos(theata)

	return (float)theata;
}
float miiMinPathModel::CalcSegThetaN(int cx, int cy, int cz,float currentvec[3])
{
	float gLastVesselDirect[3], gCurrentVessel[3];
	gLastVesselDirect[0] = m_iStartPointWorld.x - m_iLastStartPointWorld.x;
	gLastVesselDirect[1] = m_iStartPointWorld.y - m_iLastStartPointWorld.y;
	gLastVesselDirect[2] = m_iStartPointWorld.z - m_iLastStartPointWorld.z;

	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	/*gCurrentVessel[0] = fCord[0] - m_iStartPointWorld.x;
	gCurrentVessel[1] = fCord[1] - m_iStartPointWorld.y;
	gCurrentVessel[2] = fCord[2] - m_iStartPointWorld.z;*/
	//set u1 as vessel vector
	gCurrentVessel[0] =currentvec[0];
	gCurrentVessel[1] =currentvec[1];
	gCurrentVessel[2] =currentvec[2];
	//record the direction vector
	//return zxh::VectorOP_Cosine( vModel, vVessel, 3 ) ; 
	double fLastVesDirectNorm = sqrt((double)(gLastVesselDirect[0]*gLastVesselDirect[0] + gLastVesselDirect[1]*gLastVesselDirect[1] + gLastVesselDirect[2]*gLastVesselDirect[2]));
	double fCurrVesDirectNorm = sqrt((double)(gCurrentVessel[0]*gCurrentVessel[0] + gCurrentVessel[1]*gCurrentVessel[1] + gCurrentVessel[2]*gCurrentVessel[2]));

	double theata = (double)(gLastVesselDirect[0] * gCurrentVessel[0] + gLastVesselDirect[1] * gCurrentVessel[1] + gLastVesselDirect[2] * gCurrentVessel[2])\
		/ ( fLastVesDirectNorm * fCurrVesDirectNorm+0.0000001);//"theata" here means  cos(theata)

	return (float)theata;
}

bool miiMinPathModel::SelectVV(int cx, int cy,int cz,float u1[3],float currentvec[3])
{
	float fEstVV[3];//estimated vessel vector
	float u1Cont[3];//contrary of u1
	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	fEstVV[0] = fCord[0] - m_iStartPointWorld.x;
	fEstVV[1] = fCord[1] - m_iStartPointWorld.y;
	fEstVV[2] = fCord[2] - m_iStartPointWorld.z;
	//calclulate the costheta of u1 and estimated vessel vector
	double du1Norm = sqrt((double)(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2]));
	double dEstVVNorm = sqrt((double)(fEstVV[0]*fEstVV[0] + fEstVV[1]*fEstVV[1] + fEstVV[2]*fEstVV[2]));
	double costheta = (double)(u1[0] * fEstVV[0] + u1[1] * fEstVV[1] + u1[2] * fEstVV[2])\
		/ ( du1Norm * dEstVVNorm+0.0000001);//"theata" here means  cos(theata)
	if (costheta>=0)
	{
		currentvec[0]=u1[0];
		currentvec[1]=u1[1];
		currentvec[2]=u1[2];
	}
	else
	{
		currentvec[0]=-u1[0];
		currentvec[1]=-u1[1];
		currentvec[2]=-u1[2];
	}

	return true;
}

bool miiMinPathModel::SelectVV_3S(int cx, int cy,int cz,float u1[3],float currentvec[3])
{
	float fEstVV[3];//estimated vessel vector
	float u1Cont[3];//contrary of u1
	float fCord[3];
	fCord[0] = (float)cx;
	fCord[1] = (float)cy;
	fCord[2] = (float)cz;
	m_pBaseImgInfo->ImageToWorld(fCord);
	fEstVV[0] = fCord[0] - m_iStartPointWorld.x;
	fEstVV[1] = fCord[1] - m_iStartPointWorld.y;
	fEstVV[2] = fCord[2] - m_iStartPointWorld.z;
	//calclulate the costheta of u1 and estimated vessel vector
	double du1Norm = sqrt((double)(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2]));
	double dEstVVNorm = sqrt((double)(fEstVV[0]*fEstVV[0] + fEstVV[1]*fEstVV[1] + fEstVV[2]*fEstVV[2]));
	double costheta = (double)(u1[0] * fEstVV[0] + u1[1] * fEstVV[1] + u1[2] * fEstVV[2])\
		/ ( du1Norm * dEstVVNorm+0.0000001);//"theata" here means  cos(theata)
	if (costheta>=0.866)
	{
		currentvec[0]=u1[0];
		currentvec[1]=u1[1];
		currentvec[2]=u1[2];
		return true;
	}
	else if(costheta<=-0.866)
	{
		currentvec[0]=-u1[0];
		currentvec[1]=-u1[1];
		currentvec[2]=-u1[2];
		return true;
	}
	else
	{
		currentvec[0]=u1[0];
		currentvec[1]=u1[1];
		currentvec[2]=u1[2];
		return false;
	}


	return true;
}

int miiMinPathModel::FindnxyPoiC(miiCNode<double,float>dnxyzWorld)
{
	int NPosi=0;
	float fMinNDistmm=1000000;
	float fnxyzWorld[3]={dnxyzWorld.x,dnxyzWorld.y,dnxyzWorld.z};
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		float Fcord[3];
		Fcord[0]=m_vModelPointsWorld[i].x;
		Fcord[1]=m_vModelPointsWorld[i].y;
		Fcord[2]=m_vModelPointsWorld[i].z;
		float MPTNDistmm=zxh::VectorOP_Distance(Fcord,fnxyzWorld,3);
		if (MPTNDistmm<fMinNDistmm)
		{
			fMinNDistmm=MPTNDistmm;
			NPosi = i; 
		} 
	}
	return NPosi;
}
// Function Name: SearchStartPoint()
//
// Parameters: 
//
// Description: by using the threshold segmentation to find the start point
//
// Returns: 
//
bool miiMinPathModel::SearchStartPoint(const short *sRawImg, int nOldStartPoint[3], int nNewStartPoint[3], int nSearchRange)
{
/*	int nSearchRange = 100;*/

	int nXL = nOldStartPoint[0] - nSearchRange;
	int nXR = nOldStartPoint[0] + nSearchRange;
	int nYL = nOldStartPoint[1] - nSearchRange;
	int nYR = nOldStartPoint[1] + nSearchRange;
	int nZL = nOldStartPoint[2] - nSearchRange;
	int nZR = nOldStartPoint[2] + nSearchRange;

	// avoid out of range
	if (nXL < 0) nXL = 0;
	if (nXR > m_nImgWX) nXR = m_nImgWX;
	if (nYL < 0) nYL = 0;
	if (nYR > m_nImgWY) nYR = m_nImgWY;
	if (nZL < 0) nZL = 0;
	if (nZR > m_nImgWZ) nZR = m_nImgWZ;

	unsigned long nCX = 0, nCY = 0, nCZ = 0;
	unsigned int nCount = 0;

	for (int iz = nZL; iz < nZR; iz++)
	{
		for (int iy = nYL; iy < nYR; iy++)
		{
			for (int ix = nXL; ix < nXR; ix++)
			{
				if (sRawImg[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] > 1200)
				{
					nCX += ix;
					nCY += iy;
					nCZ += iz;

					nCount++;
				}
			}
		}
	}

	// 
	nNewStartPoint[0] = nCX / nCount;
	nNewStartPoint[1] = nCY / nCount;
	nNewStartPoint[2] = nCZ / nCount;

// 	m_vModelPoints[0].x = m_iStartPoint.x;
// 	m_vModelPoints[1].y = m_iStartPoint.y;
// 	m_vModelPoints[2].z = m_iStartPoint.z;

	m_iSgmtPoint.x = m_iStartPoint.x;
	m_iSgmtPoint.y = m_iStartPoint.y;
	m_iSgmtPoint.z = m_iStartPoint.z;

	return true;
}
// Function Name: GetModPointCos()
//
// Parameters: position
//
// Description: get the cosine value of angel between Fvector(forward point to i )  and Bvector(i to backforward point)
//
// Returns: cosine value
//
float miiMinPathModel::GetModPointCos(int i)//Whole Add by JDQ
{
	miiCNode<double, float> fSecondModelPointWorld=m_vModelPointsWorld[1];
	miiCNode<double, float> fLastModelPointWorld=m_vModelPointsWorld[m_vModelPointsWorld.size()-1];
	/*if((CalcDistmm2(fSecondModelPointWorld,m_vModelPointsWorld[i]))<ZXHJDQCAE_FDSETPLENGTH*ZXHJDQCAE_FDSETPLENGTH||(CalcDistmm2(fLastModelPointWorld,m_vModelPointsWorld[i])<ZXHJDQCAE_FDSETPLENGTH*ZXHJDQCAE_FDSETPLENGTH))
	{
	return 1;
	}*/
	//Fi means in front of i
	int Fi=i-MODELVECTanLENGTH/m_meandis_of_coronarymodel;
	int Bi=i+MODELVECTanLENGTH/m_meandis_of_coronarymodel;
	Fi=CorrectTruePositionOfCoronaryModel(Fi);
	Bi=CorrectTruePositionOfCoronaryModel(Bi);
	float vModelF[3], vModelB[3];
	vModelF[0] = m_vModelPointsWorld[i].x - m_vModelPointsWorld[Fi].x;
	vModelF[1] = m_vModelPointsWorld[i].y - m_vModelPointsWorld[Fi].y;
	vModelF[2] = m_vModelPointsWorld[i].z - m_vModelPointsWorld[Fi].z;


	vModelB[0] =m_vModelPointsWorld[Bi].x-m_vModelPointsWorld[i].x;
	vModelB[1] =m_vModelPointsWorld[Bi].y-m_vModelPointsWorld[i].y;
	vModelB[2] =m_vModelPointsWorld[Bi].z-m_vModelPointsWorld[i].z;

	double fModelFNorm = sqrt((double)(vModelF[0]*vModelF[0] + vModelF[1]*vModelF[1] + vModelF[2]*vModelF[2]));
	double fModelBNorm = sqrt((double)(vModelB[0]*vModelB[0] + vModelB[1]*vModelB[1] + vModelB[2]*vModelB[2]));

	float costheta = (float)(vModelF[0] * vModelB[0] + vModelF[1]*vModelB[1] + vModelF[2] * vModelB[2])/( fModelFNorm * fModelBNorm);//
	//std::cout<<"PointCos"<<costheta<<"\n";
	return costheta;

}
// Function Name: GetModVectorCos()
//
// Parameters: position
//
// Description: get the cosine value of angel between tangent vector of point i and tangent vetor of m_nPosi
//
// Returns: cosine value
//
float miiMinPathModel::GetModVectorCos(int i)
{
	int tanDistPix=MODELVECTanLENGTH/m_meandis_of_coronarymodel;
	float NPointTangentVector[3];
	int forwardi = i+tanDistPix ; 
	int backwardi = i-tanDistPix ; 
	int forwardP = m_nOneSegModelVectorStartPos+tanDistPix ; 
	int backwardP =m_nOneSegModelVectorStartPos-tanDistPix ; 
	forwardi = CorrectTruePositionOfCoronaryModel( forwardi ) ; 
	backwardi=CorrectTruePositionOfCoronaryModel( backwardi ) ;
	forwardP = CorrectTruePositionOfCoronaryModel( forwardP) ; 
	backwardP=CorrectTruePositionOfCoronaryModel( backwardP) ;
	if (m_nOneSegModelVectorStartPos==0)//m_vModelPointsWorld[0]is the start segment point;m_fModelVecWorld is the current vector from start point to current model point
	{
		m_fModelVecWorld[0] = m_vModelPointsWorld[1].x - m_vModelPointsWorld[0].x;
		m_fModelVecWorld[1] = m_vModelPointsWorld[1].y - m_vModelPointsWorld[0].y;
		m_fModelVecWorld[2] = m_vModelPointsWorld[1].z - m_vModelPointsWorld[0].z;
	}
	else
	{
		m_fModelVecWorld[0]=m_vModelPointsWorld[forwardP].x-m_vModelPointsWorld[backwardP].x;
		m_fModelVecWorld[1]=m_vModelPointsWorld[forwardP].y-m_vModelPointsWorld[backwardP].y;
		m_fModelVecWorld[2]=m_vModelPointsWorld[forwardP].z-m_vModelPointsWorld[backwardP].z;
	}

	NPointTangentVector[0] =m_vModelPointsWorld[forwardi].x-m_vModelPointsWorld[backwardi].x; 
	NPointTangentVector[1] =m_vModelPointsWorld[forwardi].y-m_vModelPointsWorld[backwardi].y; 
	NPointTangentVector[2] =m_vModelPointsWorld[forwardi].z-m_vModelPointsWorld[backwardi].z; 

	double fModelSVecNorm = sqrt((double)(m_fModelVecWorld[0]*m_fModelVecWorld[0] + m_fModelVecWorld[1]*m_fModelVecWorld[1] +m_fModelVecWorld[2]*m_fModelVecWorld[2]));
	double fModelTVecNorm = sqrt((double)(NPointTangentVector[0]*NPointTangentVector[0] +NPointTangentVector[1]*NPointTangentVector[1] + NPointTangentVector[2]*NPointTangentVector[2]));

	float costheta = (float)(m_fModelVecWorld[0]*NPointTangentVector[0] + m_fModelVecWorld[1]*NPointTangentVector[1] +m_fModelVecWorld[2]*NPointTangentVector[2])/( fModelSVecNorm* fModelTVecNorm);//
	//std::cout<<"VectorCos"<<costheta<<"\n";
	return costheta;
}



/**********************************---End Line---**************************************/
