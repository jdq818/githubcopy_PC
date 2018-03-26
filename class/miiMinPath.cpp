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
#include "miiMinPath.h"
int GetTimeNow()
{
	time_t now_time;
	now_time = time(NULL);
	return now_time;
}
//
miiMinPath::miiMinPath(int nImgWX, int nImgWY, int nImgWZ,float nImgSpacing):\
					mZXH_MinLengthForBadPointCheck(50), mZXH_PencentageForBadPointLength(0.5), m_nImgWX(nImgWX), m_nImgWY(nImgWY), m_nImgWZ(nImgWZ),m_nImgSpacing(),m_nModelPos(0),m_nOneSegModelVectorStartPos(0),m_nOneSegModelVectorEndPos(0), m_nFMMEvlNum(0),ZXHJDQCAE_FDSETPLENGTH(1),VESELLSmgTanLENGTH(3),MODELVECTanLENGTH(3),DistRangeAroundCPoint(10),m_meandis_of_coronarymodel(0),MINDISTTOPOINT2(1),mZXH_MinPencentageForCPointLength(0.5),mZXH_MaxPencentageForCPointLength(1.5),m_ftoallength_model(0)
						//miiMinPath::miiMinPath(int nImgWX, int nImgWY, int nImgWZ):m_nImgWX(nImgWX), m_nImgWY(nImgWY), m_nImgWZ(nImgWZ)//Change by JDQ
{
	m_fCurrModelCArchLength = 0 ; 
	m_fMaxSPlength=0;
	m_dU = new double[nImgWX * nImgWY * nImgWZ];
	m_nFmMap = new int[nImgWX * nImgWY * nImgWZ];
	m_dP2 = new double[nImgWX * nImgWY * nImgWZ];
	
	//sAMinNodeImg=new short[nImgWX * nImgWY * nImgWZ];//add by JDQ
	//sABadPointImg=new short[nImgWX * nImgWY * nImgWZ];
	//sBadPointMask=new short[nImgWX * nImgWY * nImgWZ];//add by JDQ
	m_sUnseenImgIntensity=new double[nImgWX * nImgWY * nImgWZ];//add by JDQ
	m_sVes=new double[nImgWX * nImgWY * nImgWZ] ;//add by JDQ
	m_dSims=new double[nImgWX * nImgWY * nImgWZ];
	m_dSimd=new double[nImgWX * nImgWY * nImgWZ];
	m_sNormVes=new double[nImgWX * nImgWY * nImgWZ] ;//add by JDQ

	// create a min-heap for Narrow Band
	m_nOneSegModelVectorStartPos=0;
	m_nOneSegModelVectorEndPos=1;
	m_iMinHeap = new miiMinHeap<>();
	m_nSum=0;
	m_nSumFmItr=0;
	m_nFMVedPoiNUM=1;
	m_nFMVedPoiNUMLast=0;
	m_nLn=0;
	m_nMAXLn=0;
	m_fMinSpeed=0.05;
	m_nFMMTimeS=0;
}

//
miiMinPath::~miiMinPath()
{
	if (m_dU != NULL)
		delete[] m_dU;

	if (m_nFmMap != NULL)	
		delete[] m_nFmMap;
	

	if (m_dP2 != NULL)
		delete[] m_dP2;

	if (m_iMinHeap != NULL)
		delete m_iMinHeap;
	
	if (m_sUnseenImgIntensity != NULL)
		delete[] m_sUnseenImgIntensity;
	if (m_sVes != NULL)
		delete[] m_sVes;
	if (m_dSims != NULL)
		delete[] m_dSims;
	if (m_dSimd != NULL)
		delete[] m_dSimd;
	if (m_sNormVes != NULL)
		delete[] m_sNormVes;
}

// Function Name: FastMarchingInit()
//
// Parameters: *sImgData: the pointer for the data of 3D image
//				nMethod: the option for the potential function 
//						'0' - Gradient
//						'1' -
//
// Description: initiate the FMM's data, including Narrow Band, distance function, and FMM map.  
//
// Returns: 
//
void miiMinPath::FastMarchingInitBase(const short *sImgData, float nStartPoint[], float nEndPoint[])
{
	float nCord[3];

	// set starting point
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	m_iStartPointWorld.x =nStartPoint[0]; 
	m_iStartPointWorld.y =nStartPoint[1]; 
	m_iStartPointWorld.z =nStartPoint[2]; 
	m_iStartPointWorld.val = 0;
	m_pBaseImgInfo->WorldToImage(nCord);	
	m_iStartPoint.x = (int)(nCord[0]+0.5);
	m_iStartPoint.y = (int)(nCord[1]+0.5);
	m_iStartPoint.z = (int)(nCord[2]+0.5);
	m_iStartPoint.val = 0;
	if (m_iStartPoint.x >= m_nImgWX)
		m_iStartPoint.x = m_nImgWX - 1;
	if (m_iStartPoint.x < 0)
		m_iStartPoint.x = 0;
	if (m_iStartPoint.y >= m_nImgWY)
		m_iStartPoint.y = m_nImgWY - 1;
	if (m_iStartPoint.y < 0)
		m_iStartPoint.y = 0;
	if (m_iStartPoint.z >= m_nImgWZ)
		m_iStartPoint.z = m_nImgWZ - 1;
	if (m_iStartPoint.z < 0)
		m_iStartPoint.z = 0;

	// set end point
	nCord[0] = nEndPoint[0];
	nCord[1] = nEndPoint[1];
	nCord[2] = nEndPoint[2];
	m_iEndPointWorld.x =nEndPoint[0]; 
	m_iEndPointWorld.y =nEndPoint[1]; 
	m_iEndPointWorld.z =nEndPoint[2]; 
	m_iEndPointWorld.val = 0;
	m_pBaseImgInfo->WorldToImage(nCord);	
	m_iEndPoint.x = (int)(nCord[0]+0.5);
	m_iEndPoint.y = (int)(nCord[1]+0.5);
	m_iEndPoint.z = (int)(nCord[2]+0.5);
	m_iEndPoint.val = 0;
	if (m_iEndPoint.x >= m_nImgWX)
		m_iEndPoint.x = m_nImgWX - 1;
	if (m_iEndPoint.x < 0)
		m_iEndPoint.x = 0;
	if (m_iEndPoint.y >= m_nImgWY)
		m_iEndPoint.y = m_nImgWY - 1;
	if (m_iEndPoint.y < 0)
		m_iEndPoint.y = 0;
	if (m_iEndPoint.z >= m_nImgWZ)
		m_iEndPoint.z = m_nImgWZ - 1;
	if (m_iEndPoint.z < 0)
		m_iEndPoint.z = 0;
	// initiate the number of the iteration for FMM
	m_nFmItr = 0;

	// check narrow band
	if (m_vNarrowBand.size() > 0)
	{
		m_vNarrowBand.clear();
	}
	if (m_vSpNorValue.size() > 0)
	{
		m_vSpNorValue.clear();
	}
	// add the starting point to Narrow Band
	m_iStartPoint.val=0;
	m_vNarrowBand.push_back(m_iStartPoint);
	m_vNarrowBand.push_back(m_iStartPoint);

	// build the min-heap
	m_iMinHeap->BuildMinHeap(m_vNarrowBand);

	// initiate the distance function U and FMM map
	for (int i = 0; i <  m_nImgWX * m_nImgWY * m_nImgWZ; i++)
	{
		m_dU[i] = FMM_INF;
		m_dP2[i]=0;
		m_sVes[i]=0;
        m_sNormVes[i]=0;;
		m_nFmMap[i] = FMM_FAR;
		m_dSimd[i]=0;
	
	}

	for(int nx=0;nx<m_nImgWX;nx++)
		for(int ny=0;ny<m_nImgWY;ny++)
			for(int nz=0;nz<m_nImgWZ;nz++)
			{

				m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iFarPValue);//sAMinNodeImg[i]=m_iFarPValue;
				m_zxhBadPointImg.SetPixelByGreyscale(nx,ny,nz,0,m_iFarPValue);//sABadPointImg[i]=m_iFarPValue
				m_zxhBadPointMaskImg.SetPixelByGreyscale(nx,ny,nz,0,0); //sBadPointMask[i]=0;
				m_zxhNBUvalueImg.SetPixelByGreyscale(nx,ny,nz,0,-10);
				m_zxhNBSpeedvalueImg.SetPixelByGreyscale(nx,ny,nz,0,-10.0);
				m_zxhNBPvalueImg.SetPixelByGreyscale(nx,ny,nz,0,-10);
				m_zxhNBSvalueImg.SetPixelByGreyscale(nx,ny,nz,0,-10);
				m_zxhNBDvalueImg.SetPixelByGreyscale(nx,ny,nz,0,-10);
				m_zxhNBD_SegvalueImg.SetPixelByGreyscale(nx,ny,nz,0,-10);

			}


	// set the distance from starting point to starting point as zero 
	m_dU[m_iStartPoint.z * m_nImgWY * m_nImgWX + m_iStartPoint.y * m_nImgWX + m_iStartPoint.x] = 0;

	// set starting point as 'trail' in FMM map
	m_nFmMap[m_iStartPoint.z * m_nImgWY * m_nImgWX + m_iStartPoint.y * m_nImgWX + m_iStartPoint.x] = FMM_TRIAL;

//	PotentialFunction(sImgData, nMethod);
}
void miiMinPath::FastMarchingInitBaseSM(const short *sImgData, float nStartPoint[], float nEndPoint[])
{
	float nCord[3];

	// set starting point
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	m_iStartPointWorld.x =nStartPoint[0]; 
	m_iStartPointWorld.y =nStartPoint[1]; 
	m_iStartPointWorld.z =nStartPoint[2]; 
	m_iStartPointWorld.val = 0;
	m_pBaseImgInfo->WorldToImage(nCord);	
	m_iStartPoint.x = (int)(nCord[0]+0.5);
	m_iStartPoint.y = (int)(nCord[1]+0.5);
	m_iStartPoint.z = (int)(nCord[2]+0.5);
	m_iStartPoint.val = 0;
	if (m_iStartPoint.x >= m_nImgWX)
		m_iStartPoint.x = m_nImgWX - 1;
	if (m_iStartPoint.x < 0)
		m_iStartPoint.x = 0;
	if (m_iStartPoint.y >= m_nImgWY)
		m_iStartPoint.y = m_nImgWY - 1;
	if (m_iStartPoint.y < 0)
		m_iStartPoint.y = 0;
	if (m_iStartPoint.z >= m_nImgWZ)
		m_iStartPoint.z = m_nImgWZ - 1;
	if (m_iStartPoint.z < 0)
		m_iStartPoint.z = 0;

	// set end point
	nCord[0] = nEndPoint[0];
	nCord[1] = nEndPoint[1];
	nCord[2] = nEndPoint[2];
	m_iEndPointWorld.x =nEndPoint[0]; 
	m_iEndPointWorld.y =nEndPoint[1]; 
	m_iEndPointWorld.z =nEndPoint[2]; 
	m_iEndPointWorld.val = 0;
	//set the unique endpoint from output file
	m_fSMEndPointWorld.x =nEndPoint[0]; 
	m_fSMEndPointWorld.y =nEndPoint[1]; 
	m_fSMEndPointWorld.z =nEndPoint[2]; 
	m_fSMEndPointWorld.val = 0;
	m_pBaseImgInfo->WorldToImage(nCord);	
	m_iEndPoint.x = (int)(nCord[0]+0.5);
	m_iEndPoint.y = (int)(nCord[1]+0.5);
	m_iEndPoint.z = (int)(nCord[2]+0.5);
	m_iEndPoint.val = 0;
	if (m_iEndPoint.x >= m_nImgWX)
		m_iEndPoint.x = m_nImgWX - 1;
	if (m_iEndPoint.x < 0)
		m_iEndPoint.x = 0;
	if (m_iEndPoint.y >= m_nImgWY)
		m_iEndPoint.y = m_nImgWY - 1;
	if (m_iEndPoint.y < 0)
		m_iEndPoint.y = 0;
	if (m_iEndPoint.z >= m_nImgWZ)
		m_iEndPoint.z = m_nImgWZ - 1;
	if (m_iEndPoint.z < 0)
		m_iEndPoint.z = 0;
	// initiate the number of the iteration for FMM
	m_nFmItr = 0;

	// check narrow band
	if (m_vNarrowBand.size() > 0)
	{
		m_vNarrowBand.clear();
	}

	// add the starting point to Narrow Band
	m_vNarrowBand.push_back(m_iStartPoint);
	m_vNarrowBand.push_back(m_iStartPoint);

	// build the min-heap
	m_iMinHeap->BuildMinHeap(m_vNarrowBand);

	// initiate the distance function U and FMM map
	for (int i = 0; i < (m_nImgWX * m_nImgWY * m_nImgWZ); i++)
	{
		m_dU[i] = FMM_INF;
		m_dP2[i]=0;
		m_sVes[i]=0;
        m_sNormVes[i]=0;;
		m_nFmMap[i] = FMM_FAR;
		m_dSimd[i]=0;
		m_zxhMinNodeImg.SetImageData(i,m_iFarPValue);//sAMinNodeImg[i]=m_iFarPValue;
		m_zxhBadPointImg.SetImageData(i,m_iFarPValue);//sABadPointImg[i]=m_iFarPValue
		m_zxhBadPointMaskImg.SetImageData(i,0); //sBadPointMask[i]=0;
		m_zxhNBUvalueImg.SetImageData(i,-10);
		m_zxhNBPvalueImg.SetImageData(i,-10);
		m_zxhNBSvalueImg.SetImageData(i,-10);
		m_zxhNBDvalueImg.SetImageData(i,-10);
	}

	// set the distance from starting point to starting point as zero 
	m_dU[m_iStartPoint.z * m_nImgWY * m_nImgWX + m_iStartPoint.y * m_nImgWX + m_iStartPoint.x] = 0;

	// set starting point as 'trail' in FMM map
	m_nFmMap[m_iStartPoint.z * m_nImgWY * m_nImgWX + m_iStartPoint.y * m_nImgWX + m_iStartPoint.x] = FMM_TRIAL;

//	PotentialFunction(sImgData, nMethod);
}


// Function Name: UpdateNarrowBandVal()
//
// Parameters: iNode: the pixel needed to update in narrow band
//
// Description: update the distance 
//
// Returns: 
//
bool miiMinPath::UpdateNarrowBandVal(miiCNode<> iNode)
{
	for (int i = 1; i < m_vNarrowBand.size(); i++)
	{
		if (m_vNarrowBand[i].x == iNode.x && m_vNarrowBand[i].y == iNode.y \
			&& m_vNarrowBand[i].z == iNode.z)
		{
			if (iNode.val > m_vNarrowBand[i].val)
			{
				m_vNarrowBand[i].val = iNode.val;
				m_iMinHeap->MinHeapify(m_vNarrowBand, i);
			}
			else
			{
				m_iMinHeap->HeapDecreaseKey(m_vNarrowBand, i, iNode);
			}			

			return true;
		}		
	}

	return false;
}

// Function Name: UpWind()
//
// Parameters: x, y, z: coordinate
//
// Description:   
//
// Returns: 
//
void miiMinPath::UpWind(int x, int y, int z)
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

	double a = min(u_l, u_r);
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

	a = dMax; b = dMid; c = dMin;

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

// Function Name: QuadraticRoots()
//
// Parameters: a*(x^2) + b*x + c =0
//			   dRoots[2]: two rooters
//
// Description: calculate Quadratic function 
//
// Returns: 
//
bool miiMinPath::QuadraticRoots(double a, double b, double c, double dRoots[2])
{
	double delta = b * b - 4 * a * c;

	if (delta < 0)
	{
		return false;
	} 
	else
	{
		delta = sqrt(delta);
		dRoots[0] = (-b + delta) / (2 * a + 0.0000001);
		dRoots[1] = (-b - delta) / (2 * a + 0.0000001);
	}

	return true;
}

// Function Name: FindMinPath()
//
// Parameters: nStartPoint[3]: the start point
//			   nMaxPath: the max length of the path for preventing overtime
//
// Description: find the minimal path by using back-propagation in the distant map
//
// Returns: 
//
bool miiMinPath::FindMinPath(float nStartPoint[3], int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
				{ 1, 0, 0}, \
				{ 0,-1, 0}, \
				{ 0, 1, 0}, \
				{ 0, 0,-1}, \
				{ 0, 0, 1} };

	// define the end point as the source point
	m_iEndPoint=WorldTransToImagePoint(m_iEndPointWorld);
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;

	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}
	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
    miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = nCord[0];
	iOrgStartPoint.y = nCord[1];
	iOrgStartPoint.z = nCord[2];
	CorrectImagePos(iOrgStartPoint);
	/*if (iOrgStartPoint.x >= m_nImgWX)
		iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
		iOrgStartPoint.y = m_nImgWY - 1;
	if (iOrgStartPoint.z >= m_nImgWZ)
		iOrgStartPoint.z = m_nImgWZ - 1;*/


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				// save last point (starting point)
				vSgmtMinPath.push_back(iOrgStartPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					/*if (j>0)
				   {
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
				
				    }*/
             
					m_zxhMinNodeImg.SetImageData(vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x,m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=m_iMinimalUPValue;//Add by JDQ for output purpose
					m_vMinPath.push_back(vSgmtMinPath[j]);
				}  
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
            m_fBackTrackDistmm=BackTrackDistmm;
			m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
			//GetVesselVec();
			return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;

		// save the minimal value
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld=ImageTransToWorldPoint(iTempPoint);
		vSgmtMinPath.push_back(iTempPoint);
		
        float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
		/* 
		//calculate the distance in pixel level
		int nDist = ( nMinX- iOrgStartPoint.x) * (nMinX - iOrgStartPoint.x) \
			+ (nMinY - iOrgStartPoint.y) * (nMinY - iOrgStartPoint.y) \
			+ (nMinZ - iOrgStartPoint.z) * (nMinZ - iOrgStartPoint.z);
		
		//calculate the dist in mm level
		 (nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
			+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
			+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
		*/
		nDistmm = sqrt((double)nDistmm);
        
		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{
			
		
		 /*
			float nDistmm = (nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
			+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
			+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
		 */
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
				m_zxhMinNodeImg.SetImageData(vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x,m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=100;//Add by JDQ
				/*if (j>0)
				{
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
				}*/
				
				 
			} 
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
			m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
			m_fBackTrackDistmm=BackTrackDistmm;//does not skip any points.
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
	       {

		     MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		
	       }
			m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
			//GetVesselVec();
			return true;
		}		
	}//while
}

bool miiMinPath::FindMinPathWithintensiy(const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	fEndopointVec[0]=m_iEndPointWorld.x;
	fEndopointVec[1]=m_iEndPointWorld.y;
	fEndopointVec[2]=m_iEndPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	m_iEndPoint.x=int(fEndopointVec[0]+0.5);
	m_iEndPoint.y=int(fEndopointVec[1]+0.5);
	m_iEndPoint.z=int(fEndopointVec[2]+0.5);
	m_iEndPoint.val=m_iEndPointWorld.val;
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;
	zxhImageData SmoothLinePointsImg;
	SmoothLinePointsImg.NewImage(pImageInfo);
	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}
	if (m_sMolPontInts.size() > 0)
	{
		m_sMolPontInts.clear();
	}

	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);
	/*if (iOrgStartPoint.x >= m_nImgWX)
	iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
	iOrgStartPoint.y = m_nImgWY - 1;
	if (iOrgStartPoint.z >= m_nImgWZ)
	iOrgStartPoint.z = m_nImgWZ - 1;*/


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{


					m_zxhMinNodeImg.SetImageData(vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x,m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=m_iMinimalUPValue;//Add by JDQ for output purpose
					m_vMinPath.push_back(vSgmtMinPath[j]);
					m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				m_fBackTrackDistmm=BackTrackDistmm;
				m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
				m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(m_fBackTrackDistmmSkipSomePonts>m_fCurrModelCArchLength)
				{
					//SmoothPoints(SmoothLinePointsImg);//calculate the mean coordinate every "m_zxhSmoothBasedArcLengthmm"

				}
				//GetNextSegMeanInteStd(vSgmtMinPath,sImgData);
				GetVesselVec(pImageInfo);
				GetVesselVecSmooth(pImageInfo);
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);

		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
		/* 
		//calculate the distance in pixel level
		int nDist = ( nMinX- iOrgStartPoint.x) * (nMinX - iOrgStartPoint.x) \
		+ (nMinY - iOrgStartPoint.y) * (nMinY - iOrgStartPoint.y) \
		+ (nMinZ - iOrgStartPoint.z) * (nMinZ - iOrgStartPoint.z);

		//calculate the dist in mm level
		(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
		+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
		+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
		*/
		nDistmm = sqrt((double)nDistmm);

		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{


			/*
			float nDistmm = (nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
			+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
			+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
			*/
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
				m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				m_zxhMinNodeImg.SetImageData(vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x,m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=100;//Add by JDQ
				/*if (j>0)
				{
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
				}*/


			} 
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			m_fBackTrackDistmm=BackTrackDistmm;
			m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
			m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			//GetNextSegMeanInteStd(vSgmtMinPath,sImgData);
			GetVesselVec(pImageInfo);
			GetVesselVecSmooth(pImageInfo);
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
			{

				MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));

			}


			return true;
		}		
	}//while
}
bool miiMinPath::FindMinPathWithintensiyLL_NewNeig(const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	//int gNbr[26][3] = {
	//		 {-1, 1, 1}, \
	//		 {-1, 1, 0}, \
	//		 {-1, 1, -1}, \
	//		 {-1, 0, 1}, \
 //            {-1, 0, 0}, \
	//		 {-1, 0, -1}, \
	//		 {-1, -1, 1}, \
	//		 {-1, -1, 0}, \
	//		 {-1, -1, -1}, \
	//				   { 0, 1, 1}, \
	//				   { 0,1, 0}, \
	//				   { 0, 1, -1}, \
	//				   { 0, 0,1}, \
	//				   { 0, 0, -1},\
	//				   { 0, -1, 1}, \
	//				   { 0, -1, 0}, \
	//				   { 0, -1, -1}, \
	//				    
	//		 {1, 1, 1}, \
	//		 {1, 1, 0}, \
	//		 {1, 1, -1}, \
	//		 {1, 0, 1}, \
 //            {1, 0, 0}, \
	//		 {1, 0, -1}, \
	//		 {1, -1, 1}, \
	//		 {1, -1, 0}, \
	//		 {1, -1, -1}, \
	//	};
	//int gNbr[18][3] = {
	//		 {-1, 0, 1}, \
 //            {-1, 0, 0}, \
	//		 {-1, 0, -1}, \
	//		 {-1, 1, 0}, \
	//		 {-1, -1, 0}, \
	//				   { 0, 1, 1}, \
	//				   { 0,1, 0}, \
	//				   { 0, 1, -1}, \
	//				   { 0, 0,1}, \
	//				   { 0, 0, -1},\
	//				   { 0, -1, 1}, \
	//				   { 0, -1, 0}, \
	//				   { 0, -1, -1},\
	//		 {1, 1, 0}, \
	//		 {1, 0, 1}, \
 //            {1, 0, 0}, \
	//		 {1, 0, -1}, \
	//		 {1, -1, 0}, \
	//	};
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };
	// define the end point as the source point

	float fEndopointVec[3];
	fEndopointVec[0]=m_iEndPointWorld.x;
	fEndopointVec[1]=m_iEndPointWorld.y;
	fEndopointVec[2]=m_iEndPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	m_iEndPoint.x=int(fEndopointVec[0]+0.5);
	m_iEndPoint.y=int(fEndopointVec[1]+0.5);
	m_iEndPoint.z=int(fEndopointVec[2]+0.5);
	m_iEndPoint.val=m_iEndPointWorld.val;
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;
	zxhImageData SmoothLinePointsImg;
	SmoothLinePointsImg.NewImage(pImageInfo);
	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}//m_vMinPathWorld
	if (m_vMinPathWorld.size() > 0)
	{
		m_vMinPathWorld.clear();
	}
	if (m_sMolPontInts.size() > 0)
	{
		m_sMolPontInts.clear();
	}

	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];
	fOrgStartPointWorld.val=0;

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);
	iOrgStartPoint.val=0;
	CorrectImagePos(iOrgStartPoint);
	/*if (iOrgStartPoint.x >= m_nImgWX)
	iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
	iOrgStartPoint.y = m_nImgWY - 1;
	if (iOrgStartPoint.z >= m_nImgWZ)
	iOrgStartPoint.z = m_nImgWZ - 1;*/


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			cout<<nPathNum<<" Points have been tracked"<<endl;
				CoutSandEndPosi();//output parameter images
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{


					m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(m_nFMMEvlNum+1)*m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=m_iMinimalUPValue;//Add by JDQ for output purpose
					m_vMinPath.push_back(vSgmtMinPath[j]);
					m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
					miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
					m_vMinPathWorld.push_back(fSgmPointWorld);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				m_fBackTrackDistmm=BackTrackDistmm;
				m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
				m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(m_fBackTrackDistmmSkipSomePonts>1.25*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
				if(m_fBackTrackDistmmSkipSomePonts>m_fCurrModelCArchLength)
				{
					//SmoothPoints(SmoothLinePointsImg);//calculate the mean coordinate every "m_zxhSmoothBasedArcLengthmm"

				}
				GetNextSegMeanInteStd(m_vMinPath,sImgData);
				GetVesselVec(pImageInfo);
				GetVesselVecSmooth(pImageInfo);
				//CoutSandEndPosi();//output parameter images
				m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
				//SmoothPath();
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);

		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
	
		nDistmm = sqrt((double)nDistmm);

		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{

			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
				m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(m_nFMMEvlNum+1)*m_iMinimalUPValue);
				miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				m_vMinPathWorld.push_back(fSgmPointWorld);
			} 

			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			m_fBackTrackDistmm=BackTrackDistmm;
			m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
			m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			if(m_fBackTrackDistmmSkipSomePonts>1.25*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
			GetNextSegMeanInteStd(m_vMinPath,sImgData);
			GetVesselVec(pImageInfo);
			GetVesselVecSmooth(pImageInfo);
			//CoutSandEndPosi();//output parameter images
			m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
			//SmoothPath();
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
			{

				MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));

			}


			return true;
		}		
	}//while
}
bool miiMinPath::FindMinPathWithintensiyLL(const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	fEndopointVec[0]=m_iEndPointWorld.x;
	fEndopointVec[1]=m_iEndPointWorld.y;
	fEndopointVec[2]=m_iEndPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	m_iEndPoint.x=int(fEndopointVec[0]+0.5);
	m_iEndPoint.y=int(fEndopointVec[1]+0.5);
	m_iEndPoint.z=int(fEndopointVec[2]+0.5);
	m_iEndPoint.val=m_iEndPointWorld.val;
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;
	zxhImageData SmoothLinePointsImg;
	SmoothLinePointsImg.NewImage(pImageInfo);
	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}//m_vMinPathWorld
	if (m_vMinPathWorld.size() > 0)
	{
		m_vMinPathWorld.clear();
	}
	if (m_sMolPontInts.size() > 0)
	{
		m_sMolPontInts.clear();
	}

	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];
	fOrgStartPointWorld.val=0;

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);
	iOrgStartPoint.val=0;
	CorrectImagePos(iOrgStartPoint);
	/*if (iOrgStartPoint.x >= m_nImgWX)
	iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
	iOrgStartPoint.y = m_nImgWY - 1;
	if (iOrgStartPoint.z >= m_nImgWZ)
	iOrgStartPoint.z = m_nImgWZ - 1;*/


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			cout<<nPathNum<<" Points have been tracked"<<endl;
				CoutSandEndPosi();//output parameter images
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{


					m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(m_nFMMEvlNum+1)*m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=m_iMinimalUPValue;//Add by JDQ for output purpose
					m_vMinPath.push_back(vSgmtMinPath[j]);
					m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
					miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
					m_vMinPathWorld.push_back(fSgmPointWorld);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				m_fBackTrackDistmm=BackTrackDistmm;
				m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
				m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(m_fBackTrackDistmmSkipSomePonts>1.25*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
				if(m_fBackTrackDistmmSkipSomePonts>m_fCurrModelCArchLength)
				{
					//SmoothPoints(SmoothLinePointsImg);//calculate the mean coordinate every "m_zxhSmoothBasedArcLengthmm"

				}
				GetNextSegMeanInteStd(m_vMinPath,sImgData);
				GetVesselVec(pImageInfo);
				GetVesselVecSmooth(pImageInfo);
				CoutSandEndPosi();//output parameter images
				m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
				//SmoothPath();
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);

		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
	
		nDistmm = sqrt((double)nDistmm);

		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{

			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
				m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(m_nFMMEvlNum+1)*m_iMinimalUPValue);
				miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				m_vMinPathWorld.push_back(fSgmPointWorld);
			} 

			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			m_fBackTrackDistmm=BackTrackDistmm;
			m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
			m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			if(m_fBackTrackDistmmSkipSomePonts>1.25*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
			GetNextSegMeanInteStd(m_vMinPath,sImgData);
			GetVesselVec(pImageInfo);
			GetVesselVecSmooth(pImageInfo);
			CoutSandEndPosi();//output parameter images
			m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
			//SmoothPath();
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
			{

				MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));

			}


			return true;
		}		
	}//while
}
bool miiMinPath::FindMinPathWithintensiyLLV(const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	fEndopointVec[0]=m_iEndPointWorld.x;
	fEndopointVec[1]=m_iEndPointWorld.y;
	fEndopointVec[2]=m_iEndPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	m_iEndPoint.x=int(fEndopointVec[0]+0.5);
	m_iEndPoint.y=int(fEndopointVec[1]+0.5);
	m_iEndPoint.z=int(fEndopointVec[2]+0.5);
	m_iEndPoint.val=m_iEndPointWorld.val;
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;
	zxhImageData SmoothLinePointsImg;
	SmoothLinePointsImg.NewImage(pImageInfo);
	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}//m_vMinPathWorld
	if (m_vMinPathWorld.size() > 0)
	{
		m_vMinPathWorld.clear();
	}
	if (m_sMolPontInts.size() > 0)
	{
		m_sMolPontInts.clear();
	}

	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];
	fOrgStartPointWorld.val=0;

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);
	iOrgStartPoint.val=0;
	CorrectImagePos(iOrgStartPoint);
	/*if (iOrgStartPoint.x >= m_nImgWX)
	iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
	iOrgStartPoint.y = m_nImgWY - 1;
	if (iOrgStartPoint.z >= m_nImgWZ)
	iOrgStartPoint.z = m_nImgWZ - 1;*/


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			cout<<nPathNum<<" Points have been tracked"<<endl;
				CoutSandEndPosi();//output parameter images
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{


					m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(m_nFMMEvlNum+1)*m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=m_iMinimalUPValue;//Add by JDQ for output purpose
					m_vMinPath.push_back(vSgmtMinPath[j]);
					m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
					miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
					m_vMinPathWorld.push_back(fSgmPointWorld);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				m_fBackTrackDistmm=BackTrackDistmm;
				m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
				m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(m_fBackTrackDistmmSkipSomePonts>1.25*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
				if(m_fBackTrackDistmmSkipSomePonts>m_fCurrModelCArchLength)
				{
					//SmoothPoints(SmoothLinePointsImg);//calculate the mean coordinate every "m_zxhSmoothBasedArcLengthmm"

				}
				GetNextSegMeanInteStd(m_vMinPath,sImgData);
				GetVesselVec(pImageInfo);
				GetVesselVecSmooth(pImageInfo);
				CoutSandEndPosiV();//output parameter images
				m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
				//SmoothPath();
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);

		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
	
		nDistmm = sqrt((double)nDistmm);

		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{

			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
				m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(m_nFMMEvlNum+1)*m_iMinimalUPValue);
				miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				m_vMinPathWorld.push_back(fSgmPointWorld);
			} 

			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			m_fBackTrackDistmm=BackTrackDistmm;
			m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
			m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			if(m_fBackTrackDistmmSkipSomePonts>1.25*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
			GetNextSegMeanInteStd(m_vMinPath,sImgData);
			GetVesselVec(pImageInfo);
			GetVesselVecSmooth(pImageInfo);
			CoutSandEndPosiV();//output parameter images
			m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
			//SmoothPath();
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
			{

				MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));

			}


			return true;
		}		
	}//while
}
bool miiMinPath::FindMinPathWithintensiyLL_O(const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	miiCNode<>iEndPoint;
	fEndopointVec[0]=m_iEndPointWorld.x;
	fEndopointVec[1]=m_iEndPointWorld.y;
	fEndopointVec[2]=m_iEndPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	iEndPoint.x=int(fEndopointVec[0]+0.5);
	iEndPoint.y=int(fEndopointVec[1]+0.5);
	iEndPoint.z=int(fEndopointVec[2]+0.5);
	iEndPoint.val=m_iEndPointWorld.val;
	int nMinX = iEndPoint.x;
	int nMinY = iEndPoint.y;
	int nMinZ = iEndPoint.z;
	int nCentX = iEndPoint.x;
	int nCentY = iEndPoint.y;
	int nCentZ = iEndPoint.z;
	int nx, ny, nz;
	vector<miiCNode<>> vSgmtMinPath;
	vector<miiCNode<>> vMinPath;
	vector<miiCNode<double, float>>vMinPathWorld;
	if (vSgmtMinPath.size() > 0)
	{
		vSgmtMinPath.clear();
	}//m_vMinPathWorld
	if (vMinPathWorld.size() > 0)
	{
		vMinPathWorld.clear();
	}
	// save the source point
	vSgmtMinPath.push_back(iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;
	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			cout<<nPathNum<<" Points have been tracked"<<endl;
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				vSgmtMinPath.push_back(iTempPoint);
				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					vMinPath.push_back(vSgmtMinPath[j]);
					miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
					vMinPathWorld.push_back(fSgmPointWorld);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				m_fBackTrackDistmm=BackTrackDistmm;
				m_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
				m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(m_fBackTrackDistmmSkipSomePonts>1.5*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
				GetNextSegMeanInteStd(vMinPath,sImgData);
				GetVesselVecLocal(vMinPath,pImageInfo);
				GetVesselVecSmoothLocal(vMinPath,pImageInfo);
				CoutSandEndPosi();//output parameter images
				m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
				//SmoothPath();
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);
		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
		nDistmm = sqrt((double)nDistmm);
		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);
			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				vMinPath.push_back(vSgmtMinPath[j]);
				miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				vMinPathWorld.push_back(fSgmPointWorld);
			} 

			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			m_fBackTrackDistmm=BackTrackDistmm;
			m_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
			m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			if(m_fBackTrackDistmmSkipSomePonts>1.5*m_ftoallength_model)
				{
					cout<<"The result centerline has reached its maximum lenghth"<<endl;
					return false;
				}
			GetNextSegMeanInteStd(vMinPath,sImgData);
			GetVesselVecLocal(vMinPath,pImageInfo);
			GetVesselVecSmoothLocal(vMinPath,pImageInfo);
			CoutSandEndPosi();//output parameter images
			m_nFMVedPoiNUMLast=m_nFMVedPoiNUM;
			return true;
		}		
	}//while
}
bool miiMinPath::FindMinPathWithintensiyLLSM(const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	fEndopointVec[0]=m_iEndPointWorld.x;
	fEndopointVec[1]=m_iEndPointWorld.y;
	fEndopointVec[2]=m_iEndPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	m_iEndPoint.x=int(fEndopointVec[0]+0.5);
	m_iEndPoint.y=int(fEndopointVec[1]+0.5);
	m_iEndPoint.z=int(fEndopointVec[2]+0.5);
	m_iEndPoint.val=m_iEndPointWorld.val;
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;
	zxhImageData SmoothLinePointsImg;
	SmoothLinePointsImg.NewImage(pImageInfo);
	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}
	if (m_sMolPontInts.size() > 0)
	{
		m_sMolPontInts.clear();
	}

	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);
	/*if (iOrgStartPoint.x >= m_nImgWX)
	iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
	iOrgStartPoint.y = m_nImgWY - 1;
	if (iOrgStartPoint.z >= m_nImgWZ)
	iOrgStartPoint.z = m_nImgWZ - 1;*/


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{


					m_zxhMinNodeImg.SetImageData(vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x,m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=m_iMinimalUPValue;//Add by JDQ for output purpose
					m_vMinPath.push_back(vSgmtMinPath[j]);
					m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				m_fBackTrackDistmm=BackTrackDistmm;
				m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
				m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(m_fBackTrackDistmmSkipSomePonts>m_fCurrModelCArchLength)
				{
					//SmoothPoints(SmoothLinePointsImg);//calculate the mean coordinate every "m_zxhSmoothBasedArcLengthmm"

				}
				GetNextSegMeanInteStd(m_vMinPath,sImgData);
				GetVesselVecLocal(m_vMinPath,pImageInfo);
				GetVesselVecSmoothLocal(m_vMinPath,pImageInfo);
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);

		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
		/* 
		//calculate the distance in pixel level
		int nDist = ( nMinX- iOrgStartPoint.x) * (nMinX - iOrgStartPoint.x) \
		+ (nMinY - iOrgStartPoint.y) * (nMinY - iOrgStartPoint.y) \
		+ (nMinZ - iOrgStartPoint.z) * (nMinZ - iOrgStartPoint.z);

		//calculate the dist in mm level
		(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
		+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
		+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
		*/
		nDistmm = sqrt((double)nDistmm);

		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{


			/*
			float nDistmm = (nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
			+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
			+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
			*/
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
				m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				m_zxhMinNodeImg.SetImageData(vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x,m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=100;//Add by JDQ
				/*if (j>0)
				{
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
				}*/


			} 
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			m_fBackTrackDistmm=BackTrackDistmm;
			m_meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
			m_fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			GetNextSegMeanInteStd(m_vMinPath,sImgData);
			GetVesselVec(pImageInfo);
			GetVesselVecSmooth(pImageInfo);
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
			{

				MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));

			}


			return true;
		}		
	}//while
}
bool miiMinPath::FindMinPathWorld(float nStartPoint[3], int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
				{ 1, 0, 0}, \
				{ 0,-1, 0}, \
				{ 0, 1, 0}, \
				{ 0, 0,-1}, \
				{ 0, 0, 1} };

	// define the end point as the source point
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;

	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}
	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = (int)(nCord[0]+0.5);
	iOrgStartPoint.y = (int)(nCord[1]+0.5);
	iOrgStartPoint.z = (int)(nCord[2]+0.5);
	if (iOrgStartPoint.x >= m_nImgWX)
		iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
		iOrgStartPoint.y = m_nImgWY - 1;
	if (iOrgStartPoint.z >= m_nImgWZ)
		iOrgStartPoint.z = m_nImgWZ - 1;


	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				// save last point (starting point)
				vSgmtMinPath.push_back(iOrgStartPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					m_vMinPath.push_back(vSgmtMinPath[j]);
				if (j>0)
				{
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld);
				}
				}
				m_meandis_of_coronaryvessel=BackTrackDistmm/(vSgmtMinPath.size() - 1);
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;

		// save the minimal value
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld=ImageTransToWorldPoint(iTempPoint);
		vSgmtMinPath.push_back(iTempPoint);
		
        float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
		/* 
		//calculate the distance in pixel level
		int nDist = ( nMinX- iOrgStartPoint.x) * (nMinX - iOrgStartPoint.x) \
			+ (nMinY - iOrgStartPoint.y) * (nMinY - iOrgStartPoint.y) \
			+ (nMinZ - iOrgStartPoint.z) * (nMinZ - iOrgStartPoint.z);
		
		//calculate the dist in mm level
		 (nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
			+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
			+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
		*/
		nDistmm = sqrt((double)nDistmm);
        
		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{
			
		
		 /*
			float nDistmm = (nMinX - iOrgStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - iOrgStartPoint.x)*m_nImgSpacing[0] \
			+ (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]* (nMinY -iOrgStartPoint.y)*m_nImgSpacing[1]\
			+ (nMinZ - iOrgStartPoint.z)*m_nImgSpacing[2]* (nMinZ -iOrgStartPoint.z)*m_nImgSpacing[2];
		 */
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
					if (j>0)
				{
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld);
				}
			}
			m_meandis_of_coronaryvessel=BackTrackDistmm/(vSgmtMinPath.size() - 1);
			return true;
		}		
	}
}

// Function Name: FindMinPath()
//
// Parameters: nMaxPath: the max length of the path for preventing overtime
//
// Description: find the minimal path by using back-propagation in the distant map
//
// Returns: 
//
bool miiMinPath::FindMinPath(int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nMinZ = m_iEndPoint.z;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;
	int nCentZ = m_iEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<double>> vSgmtMinPath;

	// 	if (m_vMinPath.size() > 0)
	// 	{
	// 		m_vMinPath.clear();
	// 	}
	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<double> iTempPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == m_iStartPoint.x && ny == m_iStartPoint.y && nz == m_iStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				// save last point (starting point)
				vSgmtMinPath.push_back(m_iStartPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					m_vMinPath.push_back(vSgmtMinPath[j]);
				if (j>0)
				{
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld);
				}
				}
				m_meandis_of_coronaryvessel=BackTrackDistmm/(vSgmtMinPath.size() - 1);
				 m_fBackTrackDistmm=BackTrackDistmm;
				int FDist=MODELVECTanLENGTH/m_meandis_of_coronaryvessel;
			int vSgmtNum=vSgmtMinPath.size() - 1;
			
			m_fVesselVecWorld[0]=vSgmtMinPath[vSgmtNum].x-vSgmtMinPath[vSgmtNum-FDist].x;
			m_fVesselVecWorld[1]=vSgmtMinPath[vSgmtNum].y-vSgmtMinPath[vSgmtNum-FDist].y;
			m_fVesselVecWorld[2]=vSgmtMinPath[vSgmtNum].z-vSgmtMinPath[vSgmtNum-FDist].z;
				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;

		// save the minimal value
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];

		vSgmtMinPath.push_back(iTempPoint);
		float nDistmm=CalcDistmm2(fTempPointWorld,m_iStartPointWorld);
		/*int nDist = (nMinX - m_iStartPoint.x) * (nMinX - m_iStartPoint.x) \
			+ (nMinY - m_iStartPoint.y) * (nMinY - m_iStartPoint.y) \
			+ (nMinZ - m_iStartPoint.z) * (nMinZ - m_iStartPoint.z);

		*/   //Add by JDQ
          nDistmm = sqrt((double)nDistmm);
		if (nDistmm < m_nEndDistmm)//if (nDist < 0.7) Add by JDQ
		{
		 
			//calculate the dist in mm level  Add by JDQ
		 /*float nDistmm = (nMinX - m_iStartPoint.x)*m_nImgSpacing[0]*m_nImgSpacing[0]*(nMinX - m_iStartPoint.x)*m_nImgSpacing[0] \
			+ (nMinY - m_iStartPoint.y)*m_nImgSpacing[1]* (nMinY - m_iStartPoint.y)*m_nImgSpacing[1]\
			+ (nMinZ - m_iStartPoint.z)*m_nImgSpacing[2]* (nMinZ - m_iStartPoint.z)*m_nImgSpacing[2];*/
		 
			// save last point (starting point)
			vSgmtMinPath.push_back(m_iStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
				if (j>0)
				{
				fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
				BackTrackDistmm=BackTrackDistmm+CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld);
				}
			}
			m_meandis_of_coronaryvessel=BackTrackDistmm/(vSgmtMinPath.size() - 1);
			m_fBackTrackDistmm=BackTrackDistmm;
			int FDist=MODELVECTanLENGTH/m_meandis_of_coronaryvessel;
			int vSgmtNum=vSgmtMinPath.size() - 1;
			m_fVesselVecWorld[0]=vSgmtMinPath[vSgmtNum].x-vSgmtMinPath[vSgmtNum-FDist].x;
			m_fVesselVecWorld[1]=vSgmtMinPath[vSgmtNum].y-vSgmtMinPath[vSgmtNum-FDist].y;
			m_fVesselVecWorld[2]=vSgmtMinPath[vSgmtNum].z-vSgmtMinPath[vSgmtNum-FDist].z;
			return true;
		}		
	}
}
float miiMinPath::BackTrackInNB(miiCNode<double> iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return -1;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	miiCNode<double,int> miEndPoint;
	
	miEndPoint.x=iMinNode.x;
	miEndPoint.y=iMinNode.y;
	miEndPoint.z=iMinNode.z;
	miEndPoint.val=iMinNode.val;
	int nMinX = miEndPoint.x;
	int nMinY = miEndPoint.y;
	int nMinZ = miEndPoint.z;
	int nCentX = miEndPoint.x;
	int nCentY = miEndPoint.y;
	int nCentZ = miEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;
	vector<miiCNode<double,int>> mvMinPath;
	vector<miiCNode<double,float>> mvMinPathWorld;
	zxhImageData SmoothLinePointsImg;
	SmoothLinePointsImg.NewImage(pImageInfo);
	if (mvMinPath.size() > 0)
	{
		mvMinPath.clear();
	}//m_vMinPathWorld
	if (mvMinPathWorld.size() > 0)
	{
		mvMinPathWorld.clear();
	}
	// save the source point
	vSgmtMinPath.push_back(miEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return -1;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				vSgmtMinPath.push_back(iTempPoint);
				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					mvMinPath.push_back(vSgmtMinPath[j]);
					miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
					mvMinPathWorld.push_back(fSgmPointWorld);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				float mfBackTrackDistmm=BackTrackDistmm;
				float mmeandis_of_coronaryvessel=BackTrackDistmm/(mvMinPath.size() - 1);
				float mfBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(mfBackTrackDistmmSkipSomePonts>m_fCurrModelCArchLengthD*1.5)
				{
					return -1;
				}
				else
				{
				return mfBackTrackDistmmSkipSomePonts;
				}
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);

		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
	
		nDistmm = sqrt((double)nDistmm);

		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{

			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				mvMinPath.push_back(vSgmtMinPath[j]);
				miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				mvMinPathWorld.push_back(fSgmPointWorld);
			} 

			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			float mfBackTrackDistmm=BackTrackDistmm;
			float mmeandis_of_coronaryvessel=BackTrackDistmm/(mvMinPath.size() - 1);
			float mfBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.

			if(mfBackTrackDistmmSkipSomePonts>m_fCurrModelCArchLengthD*1.5)
			{
				return -1;
			}
			else
			{
				return mfBackTrackDistmmSkipSomePonts;
			}
		}		
	}//while

}
float miiMinPath::BackTrackInNBWOL(miiCNode<double> iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return -1;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	miiCNode<double,int> miEndPoint;
	
	miEndPoint.x=iMinNode.x;
	miEndPoint.y=iMinNode.y;
	miEndPoint.z=iMinNode.z;
	miEndPoint.val=iMinNode.val;
	int nMinX = miEndPoint.x;
	int nMinY = miEndPoint.y;
	int nMinZ = miEndPoint.z;
	int nCentX = miEndPoint.x;
	int nCentY = miEndPoint.y;
	int nCentZ = miEndPoint.z;
	int nx, ny, nz;

	vector<miiCNode<>> vSgmtMinPath;
	vector<miiCNode<double,int>> mvMinPath;
	vector<miiCNode<double,float>> mvMinPathWorld;
	zxhImageData SmoothLinePointsImg;
	SmoothLinePointsImg.NewImage(pImageInfo);
	if (mvMinPath.size() > 0)
	{
		mvMinPath.clear();
	}//m_vMinPathWorld
	if (mvMinPathWorld.size() > 0)
	{
		mvMinPathWorld.clear();
	}
	// save the source point
	vSgmtMinPath.push_back(miEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
	miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return -1;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					mvMinPath.push_back(vSgmtMinPath[j]);
					miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
					mvMinPathWorld.push_back(fSgmPointWorld);
				}  
				BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				float mfBackTrackDistmm=BackTrackDistmm;
				float mmeandis_of_coronaryvessel=BackTrackDistmm/(mvMinPath.size() - 1);
				float mfBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				if(mfBackTrackDistmmSkipSomePonts>1.5*m_fCurrModelCArchLengthD)
				{
					return -1;
				}
				else
				{
				return mfBackTrackDistmmSkipSomePonts;
				}
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);

		float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
	
		nDistmm = sqrt((double)nDistmm);

		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{

			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				mvMinPath.push_back(vSgmtMinPath[j]);
				miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
				mvMinPathWorld.push_back(fSgmPointWorld);
			} 

			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			float mfBackTrackDistmm=BackTrackDistmm;
			float mmeandis_of_coronaryvessel=BackTrackDistmm/(mvMinPath.size() - 1);
			float mfBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.

			if(mfBackTrackDistmmSkipSomePonts>1.5*m_fCurrModelCArchLengthD)
			{
				return -1;
			}
			else
			{
				return mfBackTrackDistmmSkipSomePonts;
			}
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
			{

				MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));

			}

		}		
	}//while

}


bool miiMinPath::Speed_BackTrack(miiCNode<double> &iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return -1;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	miiCNode<double,int> miEndPoint,iOrgStartPoint,iTempPoint;
	miiCNode<double,float> fOrgStartPointWorld,fTempPointWorld;
	miEndPoint.x=iMinNode.x;
	miEndPoint.y=iMinNode.y;
	miEndPoint.z=iMinNode.z;
	miEndPoint.val=iMinNode.val;
	int nMinX = miEndPoint.x;
	int nMinY = miEndPoint.y;
	int nMinZ = miEndPoint.z;
	int nCentX = miEndPoint.x;
	int nCentY = miEndPoint.y;
	int nCentZ = miEndPoint.z;
	int nx, ny, nz;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);

	for(int i=0;i<nMaxPath;i++)
	{

		for (int j = 0; j < 6; j++)
		{
			nx = nCentX + gNbr[j][0];
			ny = nCentY + gNbr[j][1];
			nz = nCentZ + gNbr[j][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = iMinNode.x;
				iTempPoint.y = iMinNode.y;
				iTempPoint.z = iMinNode.z;
				iTempPoint.val = (float)m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				m_zxhNBSpeedvalueImg.SetImageData(iTempPoint.z * m_nImgWY * m_nImgWX + iTempPoint.y * m_nImgWX + iTempPoint.x,iTempPoint.val);
                m_vSpValue.push_back(iTempPoint);
				return true;
				
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for i
		
		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
	}//for j
	iTempPoint.x = iMinNode.x;
	iTempPoint.y = iMinNode.y;
	iTempPoint.z = iMinNode.z;
	iTempPoint.val = (float)m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]-m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
	m_zxhNBSpeedvalueImg.SetImageData(iTempPoint.z * m_nImgWY * m_nImgWX + iTempPoint.y * m_nImgWX + iTempPoint.x,iTempPoint.val);
	m_vSpValue.push_back(iTempPoint);
	return true;
}
bool miiMinPath::Speed_BackTrackCV(miiCNode<double> &iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return -1;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	miiCNode<double,int> miEndPoint,iOrgStartPoint,iTempPoint;
	miiCNode<double,float> fOrgStartPointWorld,fTempPointWorld;
	miEndPoint.x=iMinNode.x;
	miEndPoint.y=iMinNode.y;
	miEndPoint.z=iMinNode.z;
	miEndPoint.val=iMinNode.val;
	int nMinX = miEndPoint.x;
	int nMinY = miEndPoint.y;
	int nMinZ = miEndPoint.z;
	int nCentX = miEndPoint.x;
	int nCentY = miEndPoint.y;
	int nCentZ = miEndPoint.z;
	int nx, ny, nz;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);

	for(int i=0;i<nMaxPath;i++)
	{

		for (int j = 0; j < 6; j++)
		{
			nx = nCentX + gNbr[j][0];
			ny = nCentY + gNbr[j][1];
			nz = nCentZ + gNbr[j][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = iMinNode.x;
				iTempPoint.y = iMinNode.y;
				iTempPoint.z = iMinNode.z;
				iTempPoint.val = (float)(i+1)/(m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]+0.000001);
				float fDeltaU=(float)m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				m_zxhNBSpeedvalueImg.SetImageData(iTempPoint.z * m_nImgWY * m_nImgWX + iTempPoint.y * m_nImgWX + iTempPoint.x,fDeltaU);
				m_vSpValue.push_back(iTempPoint);
				return true;

			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for i
		
		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
	}//for j
	iTempPoint.x = iMinNode.x;
	iTempPoint.y = iMinNode.y;
	iTempPoint.z = iMinNode.z;
	iTempPoint.val =(float)nMaxPath/(m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]-m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX]+0.000001);
	float fDeltaU=(float)(m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]-m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX]);
	m_zxhNBSpeedvalueImg.SetImageData(iTempPoint.z * m_nImgWY * m_nImgWX + iTempPoint.y * m_nImgWX + iTempPoint.x,fDeltaU);
	m_vSpValue.push_back(iTempPoint);
	return true;
}
bool miiMinPath::BackTrack15(miiCNode<double> &iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],miiCNode<double,int> &iEOUPoint,int nLPath)
{
	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return -1;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// define the end point as the source point

	float fEndopointVec[3];
	miiCNode<double,int> miEndPoint,iOrgStartPoint,iTempPoint;
	miiCNode<double,float> fOrgStartPointWorld,fTempPointWorld;
	miEndPoint.x=iMinNode.x;
	miEndPoint.y=iMinNode.y;
	miEndPoint.z=iMinNode.z;
	miEndPoint.val=iMinNode.val;
	int nMinX = miEndPoint.x;
	int nMinY = miEndPoint.y;
	int nMinZ = miEndPoint.z;
	int nCentX = miEndPoint.x;
	int nCentY = miEndPoint.y;
	int nCentZ = miEndPoint.z;
	int nx, ny, nz;

	// set starting point
	float nCord[3];
	float fCord[3];
	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);

	CorrectImagePos(iOrgStartPoint);

	for(int i=0;i<nLPath;i++)
	{

		for (int j = 0; j < 6; j++)
		{
			nx = nCentX + gNbr[j][0];
			ny = nCentY + gNbr[j][1];
			nz = nCentZ + gNbr[j][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				
				
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = (float)(i+1)/(m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]+0.000001);
				float fDeltaU=(float)m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				//m_zxhNBSpeedvalueImg.SetImageData(iTempPoint.z * m_nImgWY * m_nImgWX + iTempPoint.y * m_nImgWX + iTempPoint.x,fDeltaU);
                //m_vSpValue.push_back(iTempPoint);
				iEOUPoint=iTempPoint;
				return true;
				
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
				double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for i
		
		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
	}//for j
	iTempPoint.val = (float)(nLPath+1)/(m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]-m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX]+0.000001);
	iEOUPoint=iTempPoint;
	float fDeltaU=(float)(m_dU[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]-m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX]);
	//m_zxhNBSpeedvalueImg.SetImageData(iTempPoint.z * m_nImgWY * m_nImgWX + iTempPoint.y * m_nImgWX + iTempPoint.x,fDeltaU);
	return true;
}
// Function Name: GenMinPathGraph()
//
// Parameters: *pImageInfo: the image info from NIFTI
//			   *sRawData: the image volume
//				strFileName: the saved filename
//				nMethod: The type of the saved file
//						'0' - "ImageLine"
//						'1' - "Line"
// Description: Save the CA path or the image including the CA path to NIFTI file 
//
// Returns: 
//
bool miiMinPath::GenMinPathGraph(const zxhImageInfo *pImageInfo, const short *sImgData, \
	string strFileName, int nMethod)
{
	if (m_vMinPath.size() == 0)
	{
		cout << "The path of CA is null in the stage of the saving nifti!" << endl;
		return false;
	}

	short *sTrsmData = new short[m_nImgWX * m_nImgWY * m_nImgWZ];
	float fCord[3] = {0};
	int nCord[3] = {0} ;
	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				if (nMethod == 0)
				{
					fCord[0] = ix; fCord[1] = iy; fCord[2] = iz;
					m_pBaseImgInfo->ImageToWorld(fCord);
					pImageInfo->WorldToImage(fCord);
					nCord[0] = (int)(fCord[0] + 0.5);
					nCord[1] = (int)(fCord[1] + 0.5);
					nCord[2] = (int)(fCord[2] + 0.5);
					if (nCord[0] >= m_nImgWX)
						nCord[0] = m_nImgWX - 1;
					if (nCord[1] >= m_nImgWY)
						nCord[1] = m_nImgWY - 1;
					if (nCord[2] >= m_nImgWZ)
						nCord[2] = m_nImgWZ - 1;

					sTrsmData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = \
						sImgData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]];
				}
				else
				{
					sTrsmData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = 0; 
				}	

			}
		}
	}

	for (int i = 0; i < m_vMinPath.size(); i++)
	{
/*		float fCord[3];
		fCord[0] = (float)m_vMinPath[i].x;
		fCord[1] = (float)m_vMinPath[i].y;
		fCord[2] = (float)m_vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fCord);
		pImageInfo->WorldToImage(fCord);
		int nCord[3];
		nCord[0] = (int)(fCord[0]+0.5);
		nCord[1] = (int)(fCord[1]+0.5);
		nCord[2] = (int)(fCord[2]+0.5);
		if (nCord[0] >= m_nImgWX)
			nCord[0] = m_nImgWX - 1;
		if (nCord[1] >= m_nImgWY)
			nCord[1] = m_nImgWY - 1;
		if (nCord[2] >= m_nImgWZ)
			nCord[2] = m_nImgWZ - 1; */
		sTrsmData[m_vMinPath[i].z * m_nImgWY * m_nImgWX + \
			m_vMinPath[i].y * m_nImgWX + m_vMinPath[i].x] = 1000;
	}

	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };

	// label the starting and end points
	for (int i = 0; i < 6; i++)
	{
		int nx_1 = m_vMinPath[0].x + gNbr[i][0];
		int ny_1 = m_vMinPath[0].y + gNbr[i][1];
		int nz_1 = m_vMinPath[0].z + gNbr[i][2];

		int nx_2 = m_vMinPath[m_vMinPath.size()-1].x + gNbr[i][0];
		int ny_2 = m_vMinPath[m_vMinPath.size()-1].y + gNbr[i][1];
		int nz_2 = m_vMinPath[m_vMinPath.size()-1].z + gNbr[i][2];
/*
		float fCord[3];
		fCord[0] = (float)nx_1;
		fCord[1] = (float)ny_1;
		fCord[2] = (float)nz_1;
		m_pBaseImgInfo->ImageToWorld(fCord);
		pImageInfo->WorldToImage(fCord);
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

		sTrsmData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]] = 2000;
*/
		sTrsmData[nz_1 * m_nImgWY * m_nImgWX + ny_1 * m_nImgWX + nx_1] = 2000;
/*
		fCord[0] = (float)nx_2;
		fCord[1] = (float)ny_2;
		fCord[2] = (float)nz_2;
		m_pBaseImgInfo->ImageToWorld(fCord);
		pImageInfo->WorldToImage(fCord);
		nCord[0] = (int)(fCord[0]+0.5);
		nCord[1] = (int)(fCord[1]+0.5);
		nCord[2] = (int)(fCord[2]+0.5);
		if (nCord[0] >= m_nImgWX)
			nCord[0] = m_nImgWX - 1;
		if (nCord[1] >= m_nImgWY)
			nCord[1] = m_nImgWY - 1;
		if (nCord[2] >= m_nImgWZ)
			nCord[2] = m_nImgWZ - 1;
		sTrsmData[nCord[2] * m_nImgWY * m_nImgWX + nCord[1] * m_nImgWX + nCord[0]] = 2000;
		*/
		sTrsmData[nz_2 * m_nImgWY * m_nImgWX + ny_2 * m_nImgWX + nx_2] = 2000;
	}

	// save the image data
	SaveImage(m_pBaseImgInfo, sTrsmData, strFileName);

	delete[] sTrsmData;

	return true;
}


// Function Name: SaveImage()
//
// Parameters: *pImageInfo: the image info from NIFTI
//			   *sRawData: the image volume
//				strFileName: the saved filename
//
// Description: Save the image including the CA path to NIFTI file 
//
// Returns: 
//
bool miiMinPath::SaveImage(const zxhImageInfo * pImageInfo, short *sImgData, string strFileName)
{
	if (pImageInfo == NULL || sImgData == NULL)
	{
		std::cerr << "The pointer of image is null!\n";
		return false;
	}

	//define the variables
	short m_pixel;

	// 1 new image 
	zxhImageDataT<short> imgTest1; 

	imgTest1.NewImage( pImageInfo->Dimension, pImageInfo->Size, pImageInfo->Spacing, pImageInfo ) ; 
	const int *image_size = pImageInfo->Size;
	for( int ip=0; ip<imgTest1.GetNumberOfPixels(); ++ip )
		imgTest1.SetImageData( ip, sImgData[ip] ) ; 
	/*for( int i = 0; i < image_size[2]; i++ )
		for( int j = 0; j < image_size[1]; j++ )
			for( int k = 0; k < image_size[0]; k++ )
			{
				//get a pixel's value
				m_pixel = sImgData[i * image_size[1] * image_size[0] + j * image_size[0] + k];

				imgTest1.SetPixelByGreyscale(k, j, i, 0, m_pixel);
			}*/

	zxh::SaveImage(&imgTest1, strFileName); 

	return true;
}

bool miiMinPath::SaveImageD(const zxhImageInfo * pImageInfo, double *sImgData, string strFileName)
{
	if (pImageInfo == NULL || sImgData == NULL)
	{
		std::cerr << "The pointer of image is null!\n";
		return false;
	}

	//define the variables
	double m_pixel;

	// 1 new image 
	zxhImageDataT<double> imgTest1; 

	imgTest1.NewImage( pImageInfo->Dimension, pImageInfo->Size, pImageInfo->Spacing, pImageInfo ) ; 
	const int *image_size = pImageInfo->Size;
	//for( int ip=0; ip<imgTest1.GetNumberOfPixels(); ++ip )
	//	imgTest1.SetImageData( ip, sImgData[ip] ) ; 
	for( int i = 0; i < image_size[2]; i++ )
		for( int j = 0; j < image_size[1]; j++ )
			for( int k = 0; k < image_size[0]; k++ )
			{
				if(k==50&&j==49&&i==79)
				int m=0;
				//get a pixel's value
				m_pixel = sImgData[i * image_size[1] * image_size[0] + j * image_size[0] + k];

				imgTest1.SetPixelByGreyscale(k, j, i, 0, m_pixel);
			}

	zxh::SaveImage(&imgTest1, strFileName); 

	return true;
}
// Function Name: WriteCA2Vtk()
//
// Parameters: *pImageInfo: the image info from the header of Nifti
//			   *chFileName: the file name for saving vtk file 
//
// Description: 
//
// Returns: 
//
void miiMinPath::WriteCA2Vtk(char *chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	int nPointNum = m_vMinPath.size();

	float fImgPixel[3];	

	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = m_vMinPath[i].x;
		fImgPixel[1] = m_vMinPath[i].y;
		fImgPixel[2] = m_vMinPath[i].z;

		m_pBaseImgInfo->ImageToWorld(fImgPixel);
		iPoints->InsertNextPoint(fImgPixel[0], fImgPixel[1], fImgPixel[2]);
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


	// Function Name: WriteCA2Vtk()
//
// Parameters: *pImageInfo: the image info from the header of Nifti
//			   *chFileName: the file name for saving vtk file 
//
// Description: 
//
// Returns: 
//
void miiMinPath::WriteCA2VtkMML(char *chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	int nPointNum = m_vModelPointsWorld.size();

	float fImgPixel[3];	

	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = m_vModelPointsWorld[i].x;
		fImgPixel[1] = m_vModelPointsWorld[i].y;
		fImgPixel[2] = m_vModelPointsWorld[i].z;

		
		iPoints->InsertNextPoint(fImgPixel[0], fImgPixel[1], fImgPixel[2]);
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

void miiMinPath::WriteCA2Txt(char *chFileName)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = m_vMinPath.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = m_vMinPath[i].x;
		fImgPixel[1] = m_vMinPath[i].y;
		fImgPixel[2] = m_vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fImgPixel);	
		WriteFileTxt << right<<fixed<<setprecision(4) <<-fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] << '\n';

	}	

}


void miiMinPath::WriteCAIntTxt(char *chFileName)
	{
	ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = m_sMolPontInts.size();
	for (int i = 0; i < nPointNum; i++)
	{
     WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4) <<m_sMolPontInts[i] <<'\n';

	}	
	
	}
void miiMinPath::WriteCA2Vtk_O(const char *chFileName,const vector<miiCNode<double,float> >vMinPathWorld)
{

	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	int nPointNum = vMinPathWorld.size();

	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = vMinPathWorld[i].x;
		fImgPixel[1] = vMinPathWorld[i].y;
		fImgPixel[2] = vMinPathWorld[i].z;
		iPoints->InsertNextPoint(fImgPixel[0], fImgPixel[1], fImgPixel[2]);
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
void miiMinPath::WriteCA2Txt_O(const char *chFileName,const vector<miiCNode<double,float> >vMinPathWorld)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = vMinPathWorld.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
		fImgPixel[0] = vMinPathWorld[i].x;
		fImgPixel[1] = vMinPathWorld[i].y;
		fImgPixel[2] = vMinPathWorld[i].z;
		WriteFileTxt << right<<fixed<<setprecision(4) <<-fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] << '\n';

	}	

}
// Function Name: SaveU2Image()
//
// Parameters: *pImageInfo: the image info from NIFTI
//				strFileName: the saved filename
//
// Description: Save U value of FMM to NIFTI file 
//
// Returns: 
//
bool miiMinPath::SaveU2Image(const zxhImageInfo * pImageInfo, string strFileName)
{
	if (pImageInfo == NULL || m_dU == NULL)
	{
		std::cerr << "The pointer of image is null!\n";
		return false;
	}

	//define the variables
	short m_pixel;

	// 1 new image 
	zxhImageDataT<short> imgTest1; 

	imgTest1.NewImage( pImageInfo->Dimension, pImageInfo->Size, pImageInfo->Spacing, pImageInfo ) ; 
	const int *image_size = pImageInfo->Size;
	for( int iz = 0; iz < image_size[2]; iz++ )
	{
		for( int iy = 0; iy < image_size[1]; iy++ )
		{
			for( int ix = 0; ix < image_size[0]; ix++ )
			{
				//get a pixel's value
				m_pixel = m_dU[iz * image_size[1] * image_size[0] + iy * image_size[0] + ix];

				imgTest1.SetPixelByGreyscale(ix, iy, iz, 0, m_pixel);
			}
		}
	}

	zxh::SaveImage(&imgTest1, strFileName); 
	imgTest1.ReleaseMem();
	return true;
}
// Function Name: CalcDistmm2()
//                 
// Parameters:  nPixPoint1:The first Point
//              nPixPoint2:The second Point
//
// Description: calculate the mm distance between nPixPoint1 and nPixPoint2
//
// Returns: 
//
float miiMinPath::CalcDistmm2(miiCNode<double, float> nPixPoint1,miiCNode<double, float> nPixPoint2)
{  //wrong, delete this function zxh noted
	float Distmm2 = (nPixPoint1.x - nPixPoint2.x)*(nPixPoint1.x -nPixPoint2.x) + 
			(nPixPoint1.y -nPixPoint2.y)*(nPixPoint1.y - nPixPoint2.y) + 
			(nPixPoint1.z - nPixPoint2.z)*(nPixPoint1.z - nPixPoint2.z);
	return Distmm2;
}
int miiMinPath::FindImgMidPoint_DistanceThresholdPlusVectorAngle(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[])
{
	// calculate the distance(in mm level)Add by JDQ
     miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_vModelPointsWorld[0]);
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SMVector[3]={0};//current model vector
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SMVector, 3 );
	//Get the tangential direction of the m_iStartPoint
	
	if (m_nFMMEvlNum == 1)
	{
	//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
	//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
	//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2 )//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	} 
	else 
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	//else if (m_nFMMEvlNum<3)
 //  {
	//if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
	//	{

	//		m_iEndPointWorld.x = dfMinNodeWorld.x;
	//		m_iEndPointWorld.y = dfMinNodeWorld.y;
	//		m_iEndPointWorld.z = dfMinNodeWorld.z;

	//		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
	//		m_iSgmtPointWorld.y =dfMinNodeWorld.y;
	//		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
	//		*yFind=true;
	//	}
	//}
	//else 
	//{
	//	if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660) //theta<30бу
	//	{

	//		m_iEndPointWorld.x = dfMinNodeWorld.x;
	//		m_iEndPointWorld.y = dfMinNodeWorld.y;
	//		m_iEndPointWorld.z = dfMinNodeWorld.z;

	//		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
	//		m_iSgmtPointWorld.y =dfMinNodeWorld.y;
	//		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
	//		*yFind=true;
	//	}
	//}
		
	return m_nFmItr;
}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngle(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[])
{
	// calculate the distance(in mm level)Add by JDQ
     miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgMDistmm2=CalcDistmm2(m_nRealOriModelPointsWorld,dfMinNodeWorld);
	float COrgVDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriModelPointsWorld);//The m_nRealOriModelPointsWorld is changing after one model correction.
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SMVector[3]={0};//current model vector
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SMVector, 3 );
	//Get the tangential direction of the m_iStartPoint
	
	if (m_nFMMEvlNum == 1)
	{
	//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
	//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
	//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2 )//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<3)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgMDistmm2>=COrgVDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=3)
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgMDistmm2>=COrgVDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	//else if (m_nFMMEvlNum<3)
 //  {
	//if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
	//	{

	//		m_iEndPointWorld.x = dfMinNodeWorld.x;
	//		m_iEndPointWorld.y = dfMinNodeWorld.y;
	//		m_iEndPointWorld.z = dfMinNodeWorld.z;

	//		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
	//		m_iSgmtPointWorld.y =dfMinNodeWorld.y;
	//		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
	//		*yFind=true;
	//	}
	//}
	//else 
	//{
	//	if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660) //theta<30бу
	//	{

	//		m_iEndPointWorld.x = dfMinNodeWorld.x;
	//		m_iEndPointWorld.y = dfMinNodeWorld.y;
	//		m_iEndPointWorld.z = dfMinNodeWorld.z;

	//		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
	//		m_iSgmtPointWorld.y =dfMinNodeWorld.y;
	//		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
	//		*yFind=true;
	//	}
	//}
		
	return m_nFmItr;
}

	int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[])
{
	// calculate the distance(in mm level)Add by JDQ
     miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SMVector[3]={0};//current model vector
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SMVector, 3 );
	//Get the tangential direction of the m_iStartPoint
	
	if (m_nFMMEvlNum == 1&&(!IsEndPoint))
	{
	//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
	//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
	//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<3&&(!IsEndPoint))
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=3&&(!IsEndPoint))
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		} 
	}
	else if (m_nFMMEvlNum>=3&&IsEndPoint)

	{
			if (CSVDistmm2 >9&&costheta>0.8660) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}  
	}

	//else if (m_nFMMEvlNum<3)
 //  {
	//if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
	//	{

	//		m_iEndPointWorld.x = dfMinNodeWorld.x;
	//		m_iEndPointWorld.y = dfMinNodeWorld.y;
	//		m_iEndPointWorld.z = dfMinNodeWorld.z;

	//		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
	//		m_iSgmtPointWorld.y =dfMinNodeWorld.y;
	//		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
	//		*yFind=true;
	//	}
	//}
	//else 
	//{
	//	if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660) //theta<30бу
	//	{

	//		m_iEndPointWorld.x = dfMinNodeWorld.x;
	//		m_iEndPointWorld.y = dfMinNodeWorld.y;
	//		m_iEndPointWorld.z = dfMinNodeWorld.z;

	//		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
	//		m_iSgmtPointWorld.y =dfMinNodeWorld.y;
	//		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
	//		*yFind=true;
	//	}
	//}
		
	return m_nFmItr;
}
float miiMinPath::ModelPointToTanPlaneDistmm(int i,float TanVec[3],float fSgmPointWorld[3])
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


bool miiMinPath::CheckBadPoint_FindPointC(miiCNode<double, float> fSgmtPointWorld,float fmeandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3],int *iPosi )
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model of the coronary is null!" << endl;
		return false;
	}
	float fDist[3]={0};
	float MDistmm=0;int Nposi=0;
	float fMinDistmm=1000000;
	//find the point C based on the length of the backtrack length.
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{

		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=fBackTrackDistmmSkipSomePonts)
		{ 
			Nposi=i;
			break ;
		} 
	}  
	float d=CalcDistmm2(m_vModelPointsWorld[Nposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	if( Nposi==0||CalcDistmm2(m_vModelPointsWorld[Nposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
	{
		*iPosi = m_vModelPointsWorld.size()-1; 

	}
	else
	{
		float TanVec[3]={fVesselVecWorld[0],fVesselVecWorld[1],fVesselVecWorld[2]};
		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*iPosi = Nposi ;
			return false ; // 
		}
		else
		{
			int iStartOfNposi = Nposi-DistRangeAroundCPoint/m_meandis_of_coronarymodel ; 
			int iEndOfNposi = Nposi + DistRangeAroundCPoint/m_meandis_of_coronarymodel; 
			if (iStartOfNposi<0) 
				iStartOfNposi = 0;
			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
				iEndOfNposi = m_vModelPointsWorld.size()-1; 

			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
			{

				float fSgmPointWorld[3]={fSgmtPointWorld.x,fSgmtPointWorld.y,fSgmtPointWorld.z};
				//float fFSgmPointWorld[3]={m_fForwardSgmtPointWorld.x,m_fForwardSgmtPointWorld.y,m_fForwardSgmtPointWorld.z};
				float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);

				if (MPTTDistmm<fMinDistmm)
				{
					fMinDistmm=MPTTDistmm;
					*iPosi = i; 
				} 
			}
		}
	}
	return true;
}
bool miiMinPath::CheckBadPoint_FindPointCSmoothSeg(miiCNode<double, float> fSgmtPointWorld,float fm_meandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,miiCNode<double, float> fSmoothSgmtPointWorld,miiCNode<double, float> fForwardSmoothSgmtPointWorld,float fSmoothVesselVecWorld[3],int *Posi )
{
	
	int PMinPosi=m_vModelPointsWorld.size()-1,PMaxPos=m_vModelPointsWorld.size()-1,Pposi=0;
	float fMinTDistmm=1000000,fMinNDistmm=1000000;
	//find the minimal position of point C based on the length of the backtrack length.
	float MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=mZXH_MinPencentageForCPointLength*fBackTrackDistmmSkipSomePonts)
		{ 
			PMinPosi=i;
			break ;
		}
	} 
	//find the maximal position of point C based on the length of the backtrack length.
	MDistmm=0;
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=mZXH_MaxPencentageForCPointLength*fBackTrackDistmmSkipSomePonts)
		{ 
			PMaxPos=i;
			break ;
		}
	} 
	Pposi=(int)(PMinPosi+PMaxPos)/2;
	float d=CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	float fSgmPointWorld[3]={fSmoothSgmtPointWorld.x,fSmoothSgmtPointWorld.y,fSmoothSgmtPointWorld.z};
	float TanVec[3]={fSmoothVesselVecWorld[0],fSmoothVesselVecWorld[1],fSmoothVesselVecWorld[2]};
	//float TanVec[3]={m_fSmoothSgmtPointWorld.x-m_iStartPointWorld.x,m_fSmoothSgmtPointWorld.y-m_iStartPointWorld.y,m_fSmoothSgmtPointWorld.z-m_iStartPointWorld.z};
	float fFSgmPointWorld[3]={fForwardSmoothSgmtPointWorld.x,fForwardSmoothSgmtPointWorld.y,fForwardSmoothSgmtPointWorld.z};
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
			if (iStartOfNposi<0) 
				iStartOfNposi = 0;
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
	float Fcord[3];
	Fcord[0]=m_vModelPointsWorld[TPosi].x;
	Fcord[1]=m_vModelPointsWorld[TPosi].y;
	Fcord[2]=m_vModelPointsWorld[TPosi].z;
	float fTDist=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
	/*if (fMinTDistmm<1)
	{
		if(fTDist<fMinNDistmm)
			*Posi=TPosi;
		else
			*Posi=NPosi;
	}
	else
		return false;*/
	if(fMinTDistmm<1&&fMinNDistmm<fTDist)
		*Posi=NPosi;
	else
		*Posi=TPosi;
	float vesseltangent[3]={0,0,0};
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
	float modelctolasdVecdist=sqrt(zxh::VectorOP_DotProduct(modelctolascVec,modelctolascVec,3));
	/*float costheta=zxh::VectorOP_Cosine( modelvector,TanVec, 3);
	float dist1=sqrt(zxh::VectorOP_DotProduct(modelvector,modelvector,3));
	float dist2=costheta*dist1;*/
	return true ; 
}

	bool miiMinPath::CheckBadPoint_FindPointCNPosiSmoothSeg(miiCNode<double, float> fSgmtPointWorld,float fm_meandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,miiCNode<double, float> fSmoothSgmtPointWorld,miiCNode<double, float> fForwardSmoothSgmtPointWorld,float fSmoothVesselVecWorld[3],int *Posi )
{
	
	int PMinPosi=m_vModelPointsWorld.size()-1,PMaxPos=m_vModelPointsWorld.size()-1,Pposi=0,NPosi=0;
	float fMinTDistmm=1000000,fMinNDistmm=1000000;
	float fSgmPointWorld[3]={fSmoothSgmtPointWorld.x,fSmoothSgmtPointWorld.y,fSmoothSgmtPointWorld.z};
	float TanVec[3]={fSmoothVesselVecWorld[0],fSmoothVesselVecWorld[1],fSmoothVesselVecWorld[2]};
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
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
	*Posi=NPosi;
	
//	//**find the nearest distance position for point C to tangent plane**//
//	//find the minimal position of point C based on the length of the backtrack length.
//	float MDistmm=0;
//	
//	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
//	{
//		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
//		if(MDistmm>=mZXH_MinPencentageForCPointLength*fBackTrackDistmmSkipSomePonts)
//		{ 
//			PMinPosi=i;
//			break ;
//		}
//
//	} 
//	//find the maximal position of point C based on the length of the backtrack length.
//	MDistmm=0;
//	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
//	{
//		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
//		if(MDistmm>=mZXH_MaxPencentageForCPointLength*fBackTrackDistmmSkipSomePonts)
//		{ 
//			PMaxPos=i;
//			break ;
//		}
//	} 
//	Pposi=(int)(PMinPosi+PMaxPos)/2;
//	float d=CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
//	
//	//float TanVec[3]={m_fSmoothSgmtPointWorld.x-m_iStartPointWorld.x,m_fSmoothSgmtPointWorld.y-m_iStartPointWorld.y,m_fSmoothSgmtPointWorld.z-m_iStartPointWorld.z};
//	float fFSgmPointWorld[3]={fForwardSmoothSgmtPointWorld.x,fForwardSmoothSgmtPointWorld.y,fForwardSmoothSgmtPointWorld.z};
//	int TPosi=0;
//	if( Pposi==0||CalcDistmm2(m_vModelPointsWorld[Pposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
//	{
//		*Posi = m_vModelPointsWorld.size()-1; 
//
//	}
//	else
//	{
//
//		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
//		{
//			*Posi = Pposi ;
//			return false ; // 
//		}
//		else
//		{
//			int iStartOfNposi = PMinPosi-DistRangeAroundCPoint/m_meandis_of_coronarymodel ; 
//			int iEndOfNposi = PMaxPos + DistRangeAroundCPoint/m_meandis_of_coronarymodel; 
//			if (iStartOfNposi<0) 
//				iStartOfNposi = 0;
//			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
//				iEndOfNposi = m_vModelPointsWorld.size()-1; 
//
//			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
//			{
//
//
//				float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);
//			
//				if (MPTTDistmm<fMinTDistmm)
//				{
//					fMinTDistmm=MPTTDistmm;
//					TPosi = i; 
//				} 
//				
//			}
//		}
//	}
////**find the nearest distance position for point C to tangent plane**// 
	


	//float Fcord[3];
	//Fcord[0]=m_vModelPointsWorld[TPosi].x;
	//Fcord[1]=m_vModelPointsWorld[TPosi].y;
	//Fcord[2]=m_vModelPointsWorld[TPosi].z;
	//float fTDist=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
	//if (fMinTDistmm<1)
	//{
	//	if(fTDist<fMinNDistmm)
	//		*Posi=TPosi;
	//	else
	//		*Posi=NPosi;
	//}
	//else
	//	return false;
	/*if(fMinTDistmm<1&&fMinNDistmm<fTDist)
		*Posi=NPosi;
	else
		*Posi=TPosi;*/
	//float vesseltangent[3]={0,0,0};
	//float modelvector[3]={0,0,0};
	//float modeltosegvector[3]={0,0,0};
	//float modelctolasdVec[3]={0,0,0};
	//float modelctolascVec[3]={0,0,0};
	//modelvector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	//modelvector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	//modelvector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//modeltosegvector[0]=fSgmPointWorld[0]-m_vModelPointsWorld[*Posi].x;
	//modeltosegvector[1]=fSgmPointWorld[1]-m_vModelPointsWorld[*Posi].y;
	//modeltosegvector[2]=fSgmPointWorld[2]-m_vModelPointsWorld[*Posi].z;
	//modelctolasdVec[0]=m_vModelPointsWorld[*Posi].x-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x;
	//modelctolasdVec[1]=m_vModelPointsWorld[*Posi].y-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y;
	//modelctolasdVec[2]=m_vModelPointsWorld[*Posi].z-m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z;
	//modelctolascVec[0]=m_vModelPointsWorld[*Posi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	//modelctolascVec[1]=m_vModelPointsWorld[*Posi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	//modelctolascVec[2]=m_vModelPointsWorld[*Posi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//float modelvectordist=sqrt(zxh::VectorOP_DotProduct(modelvector,modelvector,3));
	//float modeltosegvectordist=sqrt(zxh::VectorOP_DotProduct(modeltosegvector,modeltosegvector,3));
	//float modelctolasdVeddist=sqrt(zxh::VectorOP_DotProduct(modelctolasdVec,modelctolasdVec,3));
	//float modelctolasdVecdist=sqrt(zxh::VectorOP_DotProduct(modelctolascVec,modelctolascVec,3));
	/*float costheta=zxh::VectorOP_Cosine( modelvector,TanVec, 3);
	float dist1=sqrt(zxh::VectorOP_DotProduct(modelvector,modelvector,3));
	float dist2=costheta*dist1;*/
	return true ; 
}
bool miiMinPath::FindMinPathForBpoints(miiCNode<double, float> fSgmtPointWorld,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath,float &fm_meandis_of_coronaryvessel,float &fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3])
{


	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}
	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
				{ 1, 0, 0}, \
				{ 0,-1, 0}, \
				{ 0, 1, 0}, \
				{ 0, 0,-1}, \
				{ 0, 0, 1} };

	// define the end point as the source point
	float fEndopointVec[3];
	miiCNode<double, int> diEndPoint;
	miiCNode<double,float> fForwardSgmtPointWorld;
	fEndopointVec[0]=fSgmtPointWorld.x;
	fEndopointVec[1]=fSgmtPointWorld.y;
	fEndopointVec[2]=fSgmtPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	diEndPoint.x=int(fEndopointVec[0]+0.5);
	diEndPoint.y=int(fEndopointVec[1]+0.5);
	diEndPoint.z=int(fEndopointVec[2]+0.5);
	diEndPoint.val=fSgmtPointWorld.val;
	int nMinX = diEndPoint.x;
	int nMinY = diEndPoint.y;
	int nMinZ = diEndPoint.z;
	int nCentX = diEndPoint.x;
	int nCentY = diEndPoint.y;
	int nCentZ = diEndPoint.z;
	int nx, ny, nz;
	vector<miiCNode<>> vSgmtMinPath;
	vector<miiCNode<>> vMinPath;
	vector<short> sMolPontInts;
	// save the source point
	vSgmtMinPath.push_back(diEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
    miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	float fBackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	

	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);
	
	CorrectImagePos(iOrgStartPoint);

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					vMinPath.push_back(vSgmtMinPath[j]);
					sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				}  
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
            fBackTrackDistmm=BackTrackDistmm;
			fm_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
			fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
			GetVesselVecfromVessel(pImageInfo,vMinPath,&fForwardSgmtPointWorld,fVesselVecWorld);
			return false;
			 }

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);
		
        float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
		
		nDistmm = sqrt((double)nDistmm);
        
		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				vMinPath.push_back(vSgmtMinPath[j]);
				sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
			} 
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			fBackTrackDistmm=BackTrackDistmm;
			fm_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
			fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			GetVesselVecfromVessel(pImageInfo,vMinPath,&fForwardSgmtPointWorld,fVesselVecWorld);
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
	       {
		     MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
	       }
			return false;
		}		
	}//while

}
bool miiMinPath::FindMinPathForBpoints_test(miiCNode<double, float> fSgmtPointWorld,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath,float *fm_meandis_of_coronaryvessel,float *fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3])
{


	if (m_dU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}
	// define a array for neighbor searching
	int gNbr[6][3] = { {-1, 0, 0}, \
				{ 1, 0, 0}, \
				{ 0,-1, 0}, \
				{ 0, 1, 0}, \
				{ 0, 0,-1}, \
				{ 0, 0, 1} };

	// define the end point as the source point
	float fEndopointVec[3];
	miiCNode<double, int> diEndPoint;
	miiCNode<double,float> fForwardSgmtPointWorld;
	fEndopointVec[0]=fSgmtPointWorld.x;
	fEndopointVec[1]=fSgmtPointWorld.y;
	fEndopointVec[2]=fSgmtPointWorld.z;
	pImageInfo->WorldToImage(fEndopointVec);
	diEndPoint.x=int(fEndopointVec[0]+0.5);
	diEndPoint.y=int(fEndopointVec[1]+0.5);
	diEndPoint.z=int(fEndopointVec[2]+0.5);
	diEndPoint.val=fSgmtPointWorld.val;
	int nMinX = diEndPoint.x;
	int nMinY = diEndPoint.y;
	int nMinZ = diEndPoint.z;
	int nCentX = diEndPoint.x;
	int nCentY = diEndPoint.y;
	int nCentZ = diEndPoint.z;
	int nx, ny, nz;
	vector<miiCNode<>> vSgmtMinPath;
	vector<miiCNode<>> vMinPath;
	vector<short> sMolPontInts;
	// save the source point
	vSgmtMinPath.push_back(diEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;
    miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	float BackTrackDistmm=0;
	float fBackTrackDistmm=0;
	int nPathNum = 0;

	// set starting point
	float nCord[3];
	float fCord[3];
	

	nCord[0] = nStartPoint[0]; 
	nCord[1] = nStartPoint[1]; 
	nCord[2] = nStartPoint[2]; 
	fOrgStartPointWorld.x=nCord[0];
	fOrgStartPointWorld.y=nCord[1];
	fOrgStartPointWorld.z=nCord[2];

	m_pBaseImgInfo->WorldToImage(nCord);	
	iOrgStartPoint.x = int(nCord[0]+0.5);
	iOrgStartPoint.y = int(nCord[1]+0.5);
	iOrgStartPoint.z =int(nCord[2]+0.5);
	
	CorrectImagePos(iOrgStartPoint);

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 6; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];
			nz = nCentZ + gNbr[i][2];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.z = nz;
				iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				vSgmtMinPath.push_back(iTempPoint);

				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					vMinPath.push_back(vSgmtMinPath[j]);
					sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				}  
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
            fBackTrackDistmm=BackTrackDistmm;
			*fm_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
			*fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
			GetVesselVecfromVessel(pImageInfo,vMinPath,&fForwardSgmtPointWorld,fVesselVecWorld);
			return false;
			 }

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			{
				// find the minimal value around the center point 
				if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				{
					nMinX = nx;
					nMinY = ny;
					nMinZ = nz;
				}
			}
		}//for

		nCentX = nMinX;
		nCentY = nMinY;
		nCentZ = nMinZ;
		fCord[0]=(float)nMinX;
		fCord[1]=(float)nMinY;
		fCord[2]=(float)nMinZ;
		// save the minimal value
		m_pBaseImgInfo->ImageToWorld(fCord);
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.z = nMinZ;
		iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		fTempPointWorld.x=fCord[0];
		fTempPointWorld.y=fCord[1];
		fTempPointWorld.z=fCord[2];
		vSgmtMinPath.push_back(iTempPoint);
		
        float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);
		
		nDistmm = sqrt((double)nDistmm);
        
		if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		{
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				vMinPath.push_back(vSgmtMinPath[j]);
				sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
			} 
			BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			fBackTrackDistmm=BackTrackDistmm;
			*fm_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
			*fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			GetVesselVecfromVessel(pImageInfo,vMinPath,&fForwardSgmtPointWorld,fVesselVecWorld);
			float MDistmm=0;
			float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
	       {
		     MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
	       }
			return false;
		}		
	}//while

}
 bool miiMinPath::FindMinPathForBpointsWithSmoothSeg(miiCNode<double, float> fSgmtPointWorld,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath,float &fm_meandis_of_coronaryvessel,float &fBackTrackDistmmSkipSomePonts,miiCNode<double, float> &fSmoothSgmtPointWorld,miiCNode<double, float> &fForwardSmoothSgmtPointWorld,float fSmoothVesselVecWorld[3])
 {

	 if (m_dU == NULL || m_nFmItr == 0)
	 {
		 cout << "Cannot find minimal path!" << endl;
		 return false;
	 }
	 // define a array for neighbor searching
	 int gNbr[6][3] = { {-1, 0, 0}, \
	 { 1, 0, 0}, \
	 { 0,-1, 0}, \
	 { 0, 1, 0}, \
	 { 0, 0,-1}, \
	 { 0, 0, 1} };

	 // define the end point as the source point
	 float fEndopointVec[3];
	 miiCNode<double, int> diEndPoint;
	 miiCNode<double,float> fForwardSgmtPointWorld;
	 fEndopointVec[0]=fSgmtPointWorld.x;
	 fEndopointVec[1]=fSgmtPointWorld.y;
	 fEndopointVec[2]=fSgmtPointWorld.z;
	 pImageInfo->WorldToImage(fEndopointVec);
	 diEndPoint.x=int(fEndopointVec[0]+0.5);
	 diEndPoint.y=int(fEndopointVec[1]+0.5);
	 diEndPoint.z=int(fEndopointVec[2]+0.5);
	 diEndPoint.val=fSgmtPointWorld.val;
	 int nMinX = diEndPoint.x;
	 int nMinY = diEndPoint.y;
	 int nMinZ = diEndPoint.z;
	 int nCentX = diEndPoint.x;
	 int nCentY = diEndPoint.y;
	 int nCentZ = diEndPoint.z;
	 int nx, ny, nz;
	 vector<miiCNode<>> vSgmtMinPath;
	 vector<miiCNode<>> vMinPath;
	 vector<short> sMolPontInts;
	 // save the source point
	 vSgmtMinPath.push_back(diEndPoint);

	 miiCNode<> iTempPoint, iOrgStartPoint;
	 miiCNode<double,float>fTempPointWorld,fSgmFPointWorld,fSgmBPointWorld,fOrgStartPointWorld;
	 float BackTrackDistmm=0;
	 float fBackTrackDistmm=0;
	 int nPathNum = 0;

	 // set starting point
	 float nCord[3];
	 float fCord[3];


	 nCord[0] = nStartPoint[0]; 
	 nCord[1] = nStartPoint[1]; 
	 nCord[2] = nStartPoint[2]; 
	 fOrgStartPointWorld.x=nCord[0];
	 fOrgStartPointWorld.y=nCord[1];
	 fOrgStartPointWorld.z=nCord[2];

	 m_pBaseImgInfo->WorldToImage(nCord);	
	 iOrgStartPoint.x = int(nCord[0]+0.5);
	 iOrgStartPoint.y = int(nCord[1]+0.5);
	 iOrgStartPoint.z =int(nCord[2]+0.5);

	 CorrectImagePos(iOrgStartPoint);

	 while (1)
	 {
		 nPathNum++;
		 if (nPathNum > nMaxPath)
		 {
			 return false;
		 }

		 for (int i = 0; i < 6; i++)
		 {
			 nx = nCentX + gNbr[i][0];
			 ny = nCentY + gNbr[i][1];
			 nz = nCentZ + gNbr[i][2];

			 // reach the starting point
			 if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
			 {
				 iTempPoint.x = nx;
				 iTempPoint.y = ny;
				 iTempPoint.z = nz;
				 iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				 vSgmtMinPath.push_back(iTempPoint);

				 for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				 {
					 vMinPath.push_back(vSgmtMinPath[j]);
					 sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
				 }  
				 BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);
				 fBackTrackDistmm=BackTrackDistmm;
				 fm_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
				 fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points
				GetSmoothVesselVecfromVessel(pImageInfo,vMinPath,fSmoothSgmtPointWorld,fForwardSmoothSgmtPointWorld,fSmoothVesselVecWorld);
				 return false;
			 }

			 // prevent out of boundary
			 if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
			 {
				 // find the minimal value around the center point 
				 if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
					 m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
				 {
					 nMinX = nx;
					 nMinY = ny;
					 nMinZ = nz;
				 }
			 }
		 }//for

		 nCentX = nMinX;
		 nCentY = nMinY;
		 nCentZ = nMinZ;
		 fCord[0]=(float)nMinX;
		 fCord[1]=(float)nMinY;
		 fCord[2]=(float)nMinZ;
		 // save the minimal value
		 m_pBaseImgInfo->ImageToWorld(fCord);
		 iTempPoint.x = nMinX;
		 iTempPoint.y = nMinY;
		 iTempPoint.z = nMinZ;
		 iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
		 fTempPointWorld.x=fCord[0];
		 fTempPointWorld.y=fCord[1];
		 fTempPointWorld.z=fCord[2];
		 vSgmtMinPath.push_back(iTempPoint);

		 float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);

		 nDistmm = sqrt((double)nDistmm);

		 if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
		 {
			 // save last point (starting point)
			 vSgmtMinPath.push_back(iOrgStartPoint);

			 for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
			 {
				 vMinPath.push_back(vSgmtMinPath[j]);
				 sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
			 } 
			 BackTrackDistmm=fCalcMinPathLength(vSgmtMinPath);//does not skip any points
			 fBackTrackDistmm=BackTrackDistmm;
			 fm_meandis_of_coronaryvessel=BackTrackDistmm/(vMinPath.size() - 1);
			 fBackTrackDistmmSkipSomePonts=fCalcMinPathLengthSkipSomePonts(vSgmtMinPath);//whole minimal path length skipping some points.
			 //GetVesselVecfromVessel(pImageInfo,vMinPath,&fForwardSgmtPointWorld,fVesselVecWorld);
			 GetSmoothVesselVecfromVessel(pImageInfo,vMinPath,fSmoothSgmtPointWorld,fForwardSmoothSgmtPointWorld,fSmoothVesselVecWorld);
			 float MDistmm=0;
			 float ndist=sqrt(CalcDistmm2(m_iOrgStartPointWorld,m_vModelPointsWorld[m_nOneSegModelVectorEndPos]));//model straight length from start(0) to current position.
			 for (int i = 1; i <=m_nOneSegModelVectorEndPos; i++)//model arc length from start(0) to current position.
			 {
				 MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
			 }
			 return false;
		 }		
	 }//while

 }
 
bool miiMinPath::IsSegmentBPoint(miiCNode<double, float> fSgmtPointWorld,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath)
{
	
	float fBackTrackDistmmSkipSomePonts=0;
	float fmeandis_of_coronaryvessel=0;
	float fVesselVecWorld[3];
	float fSmoothVesselVecWorld[3];
	miiCNode<double, float> fSmoothSgmtPointWorld;
	miiCNode<double, float> fForwardSmoothSgmtPointWorld;
	int iPosi=0;
	//FindMinPathForBpoints(fSgmtPointWorld,pImageInfo, nStartPoint, sImgData,nMaxPath,fmeandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fVesselVecWorld);
	FindMinPathForBpointsWithSmoothSeg(fSgmtPointWorld,pImageInfo,nStartPoint,sImgData,nMaxPath,fmeandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fSmoothSgmtPointWorld,fForwardSmoothSgmtPointWorld,fSmoothVesselVecWorld);
	float fCurrModelCArchLength = 0 ; 
	//CheckBadPoint_FindPointC(fSgmtPointWorld,fmeandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fVesselVecWorld,&iPosi );
	//CheckBadPoint_FindPointCSmoothSeg(fSgmtPointWorld,fmeandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fSmoothSgmtPointWorld,fForwardSmoothSgmtPointWorld,fSmoothVesselVecWorld,&iPosi );
	CheckBadPoint_FindPointCNPosiSmoothSeg(fSgmtPointWorld,fmeandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fSmoothSgmtPointWorld,fForwardSmoothSgmtPointWorld,fSmoothVesselVecWorld,&iPosi );
	if( m_fCurrModelCArchLength < mZXH_MinLengthForBadPointCheck )
		return false ;

	if(iPosi< m_nOneSegModelVectorStartPos * mZXH_PencentageForBadPointLength ) //m_nOneSegModelVectorStartPos)
		return true;

	return false;

}
bool miiMinPath::IsSegmentBPoint_test(miiCNode<double, float> fSgmtPointWorld,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath)
{
	float fBackTrackDistmmSkipSomePonts=0;
	float fm_meandis_of_coronaryvessel=0;
	float fVesselVecWorld[3];
	int iPosi=0;
	FindMinPathForBpoints(fSgmtPointWorld,pImageInfo, nStartPoint, sImgData,nMaxPath,fm_meandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fVesselVecWorld);
	CheckBadPoint_FindPointC(fSgmtPointWorld,fm_meandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fVesselVecWorld,&iPosi);//FindPointC(fSgmtPointWorld,fm_meandis_of_coronaryvessel,fBackTrackDistmmSkipSomePonts,fVesselVecWorld,&iPosi);
	if( m_fCurrModelCArchLength < mZXH_MinLengthForBadPointCheck )
		return false ;
	if(iPosi<m_nOneSegModelVectorStartPos*0.5)
		return true;
	return false;

}
	int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
		// calculate the distance(in mm level)Add by JDQ
     miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SMVector[3]={0};//current model vector
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SMVector, 3 );
	//Get the tangential direction of the m_iStartPoint
	
	if (m_nFMMEvlNum == 1&&(!IsEndPoint))
	{
	//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
	//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
	//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<3&&(!IsEndPoint))
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=3&&(!IsEndPoint))
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
		if ( CSVDistmm2 > nVectorDistmm2&&costheta>0.8660 ) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			      *yFind=true;
			//else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			//{
			//	UpdateBadPointsNeighbour(iMinNode);
			//	m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
			//	m_nBadpointSum++;
			//	m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
			//	m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			//}
		}
		//else if(m_iSgmtPointWorld.x!=m_iOrgStartPointWorld.x&&m_iSgmtPointWorld.y!=m_iOrgStartPointWorld.y&&m_iSgmtPointWorld.z!=m_iOrgStartPointWorld.z&& CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		//{  
		//	// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
		//	if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
		//	{
		//	m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
		//	m_nBadpointSum++;
		//	UpdateBadPointsNeighbour(iMinNode);
		//	m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
		//	m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
		//	}
		//}
		 
	}
	else if (m_nFMMEvlNum>=3&&IsEndPoint)

	{
			if (CSVDistmm2 >9&&costheta>0.8660) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}  
	}
		return m_nFmItr;
}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	SVector[0]=(SLVVector[0]+SMVector[0])/2;
	SVector[1]=(SLVVector[1]+SMVector[1])/2;
	SVector[2]=(SLVVector[2]+SMVector[2])/2;
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	//Get the tangential direction of the m_iStartPoint

	if (m_nFMMEvlNum == 1&&(!IsEndPoint))
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<3&&(!IsEndPoint))
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=3&&(!IsEndPoint))
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costheta>0.75 ) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign  orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
	
	else if (m_nFMMEvlNum>=3&&IsEndPoint)

	{
		if (CSVDistmm2 >9&&costheta>0.7) //theta<30бу
		{ 
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}  
	}
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector

	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	SVector[0]=(SLVVector[0]+SMVector[0])/2;
	SVector[1]=(SLVVector[1]+SMVector[1])/2;
	SVector[2]=(SLVVector[2]+SMVector[2])/2;


	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	//Get the C position for the current iMinNode
	int iPosi;
	float fPointCTan[3]={0,0,0};
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	//Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	GetPointCTangentVector(iPosi,fPointCTan);
	//Calclulate the distance from current iMinNode to the PointC tangent plane;
	float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

	if (m_nFMMEvlNum == 1&&fPctoTanPlanDist<1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<3&&fPctoTanPlanDist<1)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=3)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&fPctoTanPlanDist<1&&costheta>0.5 ) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
	
	
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQ(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	//vetor from last startpoint to current startpoint on the vessel
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	//average the two vetor
	SVector[0]=(SLVVector[0]+SMVector[0])/2;
	SVector[1]=(SLVVector[1]+SMVector[1])/2;
	SVector[2]=(SLVVector[2]+SMVector[2])/2;
	//current startpoint to current selecting point on the vessel
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	//current model vector lengthmm2
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );

	m_zxhNBD_SegvalueImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,10*costheta);
	//Get the C position for the current iMinNode
	int iPosi;
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	float fPointDTanStart[3]={0,0,0};
	GetPointCTangentVector(m_nOneSegModelVectorStartPos,fPointDTanStart);
	//calculate vector the iposi to the start point on the model
	float fMStartToiPosi[3]={0,0,0};
	fMStartToiPosi[0]=m_vModelPointsWorld[iPosi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	fMStartToiPosi[1]=m_vModelPointsWorld[iPosi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	fMStartToiPosi[2]=m_vModelPointsWorld[iPosi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//caclulate the cosine value between fMStartToiPosi and SNVVector
	float costheta1=zxh::VectorOP_Cosine( fMStartToiPosi,SNVVector, 3 );
	////Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	//float fPointCTan[3]={0,0,0};
	//GetPointCTangentVector(iPosi,fPointCTan);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

 //     //Get the tangential direction of corresponding PointD on the model of the current model start point on the vessel;
	//float fPointDTanEnd[3]={0,0,0};
	//GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPdtoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointDTanEnd);

	float fconstcosine=0.8660;
	if (m_nFMMEvlNum == 1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<3&&costheta1>0)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta>fconstcosine&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=3)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costheta1>0&&costheta>fconstcosine ) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if( CSVDistmm2 > nVectorDistmm2&&costheta1>0)

		{
			fconstcosine=zxh::maxf(0.75*fconstcosine,0);
			if (costheta>fconstcosine ) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
	
	
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel(approximate)
	float SNVVectorU1[3]={0};//Current Segpoint  to startpoint on the vessel from u1
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	
	//vetor from last startpoint to current startpoint on the vessel
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	//average the two vetor
	SVector[0]=(SLVVector[0]+SMVector[0])/2;
	SVector[1]=(SLVVector[1]+SMVector[1])/2;
	SVector[2]=(SLVVector[2]+SMVector[2])/2;
	//current startpoint to current selecting point on the vessel
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;

	//current model vector lengthmm2 and angel
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	m_zxhNBD_SegvalueImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,10*costheta);
	//set the u1 as current vessel vetor
	SNVVectorU1[0]=fCurMinNodevector[0];
	SNVVectorU1[1]=fCurMinNodevector[1];
	SNVVectorU1[2]=fCurMinNodevector[2];
	float costhetaU1=abs(zxh::VectorOP_Cosine( SNVVectorU1,SVector, 3 ));
	//Get the tangent of point D and calculate the distance coming from current point
	
	float fPointDTanEnd[3]={0,0,0};
	GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	float fSegPoint[3]={0};
	fSegPoint[0]=dfMinNodeWorld.x;
	fSegPoint[1]=dfMinNodeWorld.y;
	fSegPoint[2]=dfMinNodeWorld.z;
	float TPTanDistmm=ModelPointToTanPlaneDistmm(m_nOneSegModelVectorEndPos,fPointDTanEnd,fSegPoint);

	//combine two thetas
	float costheta2=abs(zxh::VectorOP_Cosine( fPointDTanEnd,SNVVectorU1, 3 ));
	costheta2=costhetaU1;
	//calculate vector the iposi to the start point on the model
	int iPosi;
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	float fMStartToiPosi[3]={0,0,0};
	fMStartToiPosi[0]=m_vModelPointsWorld[iPosi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	fMStartToiPosi[1]=m_vModelPointsWorld[iPosi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	fMStartToiPosi[2]=m_vModelPointsWorld[iPosi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//caclulate the cosine value between fMStartToiPosi and SNVVector
	float costheta1=zxh::VectorOP_Cosine( fMStartToiPosi,SMVector, 3 );
	////Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	//float fPointCTan[3]={0,0,0};
	//GetPointCTangentVector(iPosi,fPointCTan);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

 //     //Get the tangential direction of corresponding PointD on the model of the current model start point on the vessel;
	//float fPointDTanEnd[3]={0,0,0};
	//GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPdtoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointDTanEnd);

	float fconstcosine=0.8660;
	float fconstcosineU1=0.8660;
	if (m_nFMMEvlNum == 1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<3&&costheta1>0)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costheta2>fconstcosineU1&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=3)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costheta1>0&&costheta2>fconstcosineU1) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if( CSVDistmm2 > nVectorDistmm2&&costheta1>0)

		{
			//fconstcosine=zxh::maxf(0.75*fconstcosine,0);
			//fconstcosineU1=zxh::maxf(0.75*fconstcosineU1,0.5);
			if (costheta2>fconstcosineU1) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel(approximate)
	float SNVVectorU1[3]={0};//Current Segpoint  to startpoint on the vessel from u1
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	
	//vetor from last startpoint to current startpoint on the vessel
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	//average the two vetor
	SVector[0]=SMVector[0];
	SVector[1]=SMVector[1];
	SVector[2]=SMVector[2];
	//current startpoint to current selecting point on the vessel
	SNVVector[0]=fCurMinNodevector[0];
	SNVVector[1]=fCurMinNodevector[1];
	SNVVector[2]=fCurMinNodevector[2];

	//current model vector lengthmm2 and angel
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	m_zxhNBD_SegvalueImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,10*costheta);
	//set the u1 as current vessel vetor
	SNVVectorU1[0]=fCurMinNodevector[0];
	SNVVectorU1[1]=fCurMinNodevector[1];
	SNVVectorU1[2]=fCurMinNodevector[2];
	float costhetaU1=zxh::VectorOP_Cosine( SNVVectorU1,SVector, 3 );
	//Get the tangent of point D and calculate the distance coming from current point
	
	float fPointDTanEnd[3]={0,0,0};
	GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	float fSegPoint[3]={0};
	fSegPoint[0]=dfMinNodeWorld.x;
	fSegPoint[1]=dfMinNodeWorld.y;
	fSegPoint[2]=dfMinNodeWorld.z;
	float TPTanDistmm=ModelPointToTanPlaneDistmm(m_nOneSegModelVectorEndPos,fPointDTanEnd,fSegPoint);

	//combine two thetas
	float costheta2=abs(zxh::VectorOP_Cosine( fPointDTanEnd,SNVVectorU1, 3 ));
	//calculate vector the iposi to the start point on the model
	int iPosi;
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	float fMStartToiPosi[3]={0,0,0};
	fMStartToiPosi[0]=m_vModelPointsWorld[iPosi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	fMStartToiPosi[1]=m_vModelPointsWorld[iPosi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	fMStartToiPosi[2]=m_vModelPointsWorld[iPosi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//caclulate the cosine value between fMStartToiPosi and SNVVector
	float costheta1=zxh::VectorOP_Cosine( fMStartToiPosi,SMVector, 3 );
	////Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	//float fPointCTan[3]={0,0,0};
	//GetPointCTangentVector(iPosi,fPointCTan);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

 //     //Get the tangential direction of corresponding PointD on the model of the current model start point on the vessel;
	//float fPointDTanEnd[3]={0,0,0};
	//GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPdtoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointDTanEnd);

	float fconstcosine=0.8660;
	float fconstcosineU1=0.8660;
	if (m_nFMMEvlNum == 1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<6)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=6)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if( CSVDistmm2 > nVectorDistmm2&&costheta1>0)

		{
			//fconstcosine=zxh::maxf(0.75*fconstcosine,0);
			//fconstcosineU1=zxh::maxf(0.75*fconstcosineU1,0.5);
			if (costhetaU1>fconstcosineU1) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_DMP(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel(approximate)
	float SNVVectorU1[3]={0};//Current Segpoint  to startpoint on the vessel from u1
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	
	//vetor from last startpoint to current startpoint on the vessel
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	//average the two vetor
	SVector[0]=SMVector[0];
	SVector[1]=SMVector[1];
	SVector[2]=SMVector[2];
	//current startpoint to current selecting point on the vessel
	//SNVVector[0]=fCurMinNodevector[0];
	//SNVVector[1]=fCurMinNodevector[1];
	//SNVVector[2]=fCurMinNodevector[2];
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	//current model vector lengthmm2 and angel
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	m_zxhNBD_SegvalueImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,10*costheta);
	//set the u1 as current vessel vetor
	//SNVVectorU1[0]=fCurMinNodevector[0];
	//SNVVectorU1[1]=fCurMinNodevector[1];
	//SNVVectorU1[2]=fCurMinNodevector[2];
	SNVVectorU1[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVectorU1[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVectorU1[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	float costhetaU1=zxh::VectorOP_Cosine( SNVVectorU1,SVector, 3 );
	//Get the tangent of point D and calculate the distance coming from current point
	
	float fPointDTanEnd[3]={0,0,0};
	GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	float fSegPoint[3]={0};
	fSegPoint[0]=dfMinNodeWorld.x;
	fSegPoint[1]=dfMinNodeWorld.y;
	fSegPoint[2]=dfMinNodeWorld.z;
	float TPTanDistmm=ModelPointToTanPlaneDistmm(m_nOneSegModelVectorEndPos,fPointDTanEnd,fSegPoint);

	//combine two thetas
	float costheta2=abs(zxh::VectorOP_Cosine( fPointDTanEnd,SNVVectorU1, 3 ));
	//calculate vector the iposi to the start point on the model
	int iPosi;
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	float fMStartToiPosi[3]={0,0,0};
	fMStartToiPosi[0]=m_vModelPointsWorld[iPosi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	fMStartToiPosi[1]=m_vModelPointsWorld[iPosi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	fMStartToiPosi[2]=m_vModelPointsWorld[iPosi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//caclulate the cosine value between fMStartToiPosi and SNVVector
	float costheta1=zxh::VectorOP_Cosine( fMStartToiPosi,SMVector, 3 );
	////Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	//float fPointCTan[3]={0,0,0};
	//GetPointCTangentVector(iPosi,fPointCTan);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

 //     //Get the tangential direction of corresponding PointD on the model of the current model start point on the vessel;
	//float fPointDTanEnd[3]={0,0,0};
	//GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPdtoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointDTanEnd);

	float fconstcosine=0.8660;
	float fconstcosineU1=0.8660;
	if (m_nFMMEvlNum == 1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<6)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=6)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if( CSVDistmm2 > nVectorDistmm2&&costheta1>0)

		{
			//fconstcosine=zxh::maxf(0.75*fconstcosine,0);
			//fconstcosineU1=zxh::maxf(0.75*fconstcosineU1,0.5);
			if (costhetaU1>fconstcosineU1) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
		return m_nFmItr;

}

int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel(approximate)
	float SNVVectorU1[3]={0};//Current Segpoint  to startpoint on the vessel from u1
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	
	//vetor from last startpoint to current startpoint on the vessel
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	//average the two vetor
	SVector[0]=SMVector[0];
	SVector[1]=SMVector[1];
	SVector[2]=SMVector[2];
	//current startpoint to current selecting point on the vessel
	SNVVector[0]=fCurMinNodevector[0];
	SNVVector[1]=fCurMinNodevector[1];
	SNVVector[2]=fCurMinNodevector[2];

	//current model vector lengthmm2 and angel
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	m_zxhNBD_SegvalueImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,10*costheta);
	//set the u1 as current vessel vetor
	SNVVectorU1[0]=fCurMinNodevector[0];
	SNVVectorU1[1]=fCurMinNodevector[1];
	SNVVectorU1[2]=fCurMinNodevector[2];
	float costhetaU1=zxh::VectorOP_Cosine( SNVVectorU1,SVector, 3 );
	//Get the tangent of point D and calculate the distance coming from current point
	
	float fPointDTanEnd[3]={0,0,0};
	GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	float fSegPoint[3]={0};
	fSegPoint[0]=dfMinNodeWorld.x;
	fSegPoint[1]=dfMinNodeWorld.y;
	fSegPoint[2]=dfMinNodeWorld.z;
	float TPTanDistmm=ModelPointToTanPlaneDistmm(m_nOneSegModelVectorEndPos,fPointDTanEnd,fSegPoint);

	//combine two thetas
	float costheta2=abs(zxh::VectorOP_Cosine( fPointDTanEnd,SNVVectorU1, 3 ));
	//calculate vector the iposi to the start point on the model
	int iPosi;
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	float fMStartToiPosi[3]={0,0,0};
	fMStartToiPosi[0]=m_vModelPointsWorld[iPosi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	fMStartToiPosi[1]=m_vModelPointsWorld[iPosi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	fMStartToiPosi[2]=m_vModelPointsWorld[iPosi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//caclulate the cosine value between fMStartToiPosi and SNVVector
	float costheta1=zxh::VectorOP_Cosine( fMStartToiPosi,SMVector, 3 );
	////Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	//float fPointCTan[3]={0,0,0};
	//GetPointCTangentVector(iPosi,fPointCTan);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

 //     //Get the tangential direction of corresponding PointD on the model of the current model start point on the vessel;
	//float fPointDTanEnd[3]={0,0,0};
	//GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPdtoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointDTanEnd);

	float fconstcosine=0.8660;
	float fconstcosineU1=0.8660;
	if (m_nFMMEvlNum == 1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<6)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=6)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour_BPD(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if( CSVDistmm2 > nVectorDistmm2&&costheta1>0)

		{
			//fconstcosine=zxh::maxf(0.75*fconstcosine,0);
			//fconstcosineU1=zxh::maxf(0.75*fconstcosineU1,0.5);
			if (costhetaU1>fconstcosineU1) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour_BPD(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour_BPD(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DFM_WOBP(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel(approximate)
	float SNVVectorU1[3]={0};//Current Segpoint  to startpoint on the vessel from u1
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	
	//vetor from last startpoint to current startpoint on the vessel
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	//average the two vetor
	SVector[0]=SMVector[0];
	SVector[1]=SMVector[1];
	SVector[2]=SMVector[2];
	//current startpoint to current selecting point on the vessel
	SNVVector[0]=fCurMinNodevector[0];
	SNVVector[1]=fCurMinNodevector[1];
	SNVVector[2]=fCurMinNodevector[2];

	//current model vector lengthmm2 and angel
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	m_zxhNBD_SegvalueImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,10*costheta);
	//set the u1 as current vessel vetor
	SNVVectorU1[0]=fCurMinNodevector[0];
	SNVVectorU1[1]=fCurMinNodevector[1];
	SNVVectorU1[2]=fCurMinNodevector[2];
	float costhetaU1=zxh::VectorOP_Cosine( SNVVectorU1,SVector, 3 );
	//Get the tangent of point D and calculate the distance coming from current point
	
	float fPointDTanEnd[3]={0,0,0};
	GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	float fSegPoint[3]={0};
	fSegPoint[0]=dfMinNodeWorld.x;
	fSegPoint[1]=dfMinNodeWorld.y;
	fSegPoint[2]=dfMinNodeWorld.z;
	float TPTanDistmm=ModelPointToTanPlaneDistmm(m_nOneSegModelVectorEndPos,fPointDTanEnd,fSegPoint);

	//combine two thetas
	float costheta2=abs(zxh::VectorOP_Cosine( fPointDTanEnd,SNVVectorU1, 3 ));
	//calculate vector the iposi to the start point on the model
	int iPosi;
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	float fMStartToiPosi[3]={0,0,0};
	fMStartToiPosi[0]=m_vModelPointsWorld[iPosi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	fMStartToiPosi[1]=m_vModelPointsWorld[iPosi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	fMStartToiPosi[2]=m_vModelPointsWorld[iPosi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//caclulate the cosine value between fMStartToiPosi and SNVVector
	float costheta1=zxh::VectorOP_Cosine( fMStartToiPosi,SMVector, 3 );
	////Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	//float fPointCTan[3]={0,0,0};
	//GetPointCTangentVector(iPosi,fPointCTan);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

 //     //Get the tangential direction of corresponding PointD on the model of the current model start point on the vessel;
	//float fPointDTanEnd[3]={0,0,0};
	//GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPdtoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointDTanEnd);

	float fconstcosine=0.8660;
	float fconstcosineU1=0.8660;
	if (m_nFMMEvlNum == 1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<4)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=4)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour_WBPD(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if( CSVDistmm2 > nVectorDistmm2&&costheta1>0)

		{
			//fconstcosine=zxh::maxf(0.75*fconstcosine,0);
			//fconstcosineU1=zxh::maxf(0.75*fconstcosineU1,0.5);
			if (costhetaU1>fconstcosineU1) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour_WBPD(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour_WBPD(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD_O(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel(approximate)
	float SNVVectorU1[3]={0};//Current Segpoint  to startpoint on the vessel from u1
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector
	
	//vetor from last startpoint to current startpoint on the vessel
	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	//average the two vetor
	SVector[0]=SMVector[0];
	SVector[1]=SMVector[1];
	SVector[2]=SMVector[2];
	//current startpoint to current selecting point on the vessel
	SNVVector[0]=fCurMinNodevector[0];
	SNVVector[1]=fCurMinNodevector[1];
	SNVVector[2]=fCurMinNodevector[2];

	//current model vector lengthmm2 and angel
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	m_zxhNBD_SegvalueImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,10*costheta);
	//set the u1 as current vessel vetor
	SNVVectorU1[0]=fCurMinNodevector[0];
	SNVVectorU1[1]=fCurMinNodevector[1];
	SNVVectorU1[2]=fCurMinNodevector[2];
	float costhetaU1=zxh::VectorOP_Cosine( SNVVectorU1,SVector, 3 );
	//Get the tangent of point D and calculate the distance coming from current point
	
	float fPointDTanEnd[3]={0,0,0};
	GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	float fSegPoint[3]={0};
	fSegPoint[0]=dfMinNodeWorld.x;
	fSegPoint[1]=dfMinNodeWorld.y;
	fSegPoint[2]=dfMinNodeWorld.z;
	float TPTanDistmm=ModelPointToTanPlaneDistmm(m_nOneSegModelVectorEndPos,fPointDTanEnd,fSegPoint);

	//combine two thetas
	float costheta2=abs(zxh::VectorOP_Cosine( fPointDTanEnd,SNVVectorU1, 3 ));
	//calculate vector the iposi to the start point on the model
	int iPosi;
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	float fMStartToiPosi[3]={0,0,0};
	fMStartToiPosi[0]=m_vModelPointsWorld[iPosi].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	fMStartToiPosi[1]=m_vModelPointsWorld[iPosi].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	fMStartToiPosi[2]=m_vModelPointsWorld[iPosi].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	//caclulate the cosine value between fMStartToiPosi and SNVVector
	float costheta1=zxh::VectorOP_Cosine( fMStartToiPosi,SMVector, 3 );
	////Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	//float fPointCTan[3]={0,0,0};
	//GetPointCTangentVector(iPosi,fPointCTan);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);

 //     //Get the tangential direction of corresponding PointD on the model of the current model start point on the vessel;
	//float fPointDTanEnd[3]={0,0,0};
	//GetPointCTangentVector(m_nOneSegModelVectorEndPos,fPointDTanEnd);
	////Calclulate the distance from current iMinNode to the PointC tangent plane;
	//float fPdtoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointDTanEnd);

	float fconstcosine=0.8660;
	float fconstcosineU1=0.8660;
	if (m_nFMMEvlNum == 1)
	{
		//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
		//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
		//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2)//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;

			*yFind=true;
		}
	} 
	else  if (m_nFMMEvlNum<6)
	{
		if (CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	else if (m_nFMMEvlNum>=6)
	{ 
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;

		
		if ( CSVDistmm2 > nVectorDistmm2&&costhetaU1>fconstcosineU1&&costheta1>0) //theta<30бу to determine whether is a segment point
		{  
			// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				*yFind=true;
			else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			{
				UpdateBadPointsNeighbour_BPD(iMinNode);
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}
		else if( CSVDistmm2 > nVectorDistmm2&&costheta1>0)

		{
			//fconstcosine=zxh::maxf(0.75*fconstcosine,0);
			//fconstcosineU1=zxh::maxf(0.75*fconstcosineU1,0.5);
			if (costhetaU1>fconstcosineU1) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour_BPD(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
		}
		else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
		{  
			// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
			if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
			{
				m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
				m_nBadpointSum++;
				UpdateBadPointsNeighbour_BPD(iMinNode);
				m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
				m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
			}
		}

	}
		return m_nFmItr;

}
int miiMinPath::ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLLSM(miiCNode<double, int> iMinNode,bool *yFind,float CurMVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath)
{
	// calculate the distance(in mm level)Add by JDQ
	miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float COrgVDistmm2=CalcDistmm2(m_iOrgStartPointWorld,dfMinNodeWorld);//Original Vessel Point.
	float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorEndPos],m_iOrgStartPointWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// float COrgVDistmm2=CalcDistmm2(m_nRealOriSegPointWorld,dfMinNodeWorld);//first  Vessel segment Point.
	//float COrgMDistmm2=CalcDistmm2(m_vModelPointsWorld[m_nModelPos],m_nRealOriSegModelPointsWorld);//Original Model Point;The m_nRealOriModelPointsWorld is changing after one model correction.
	// 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);
	float SLVVector[3]={0};//from last startpoint to current startpoint 
	float SNVVector[3]={0};//Current Segpoint  to startpoint on the vessel
	float SVector[3]={0};
	float SMVector[3]={0};//current model vector

	SLVVector[0]=m_iStartPointWorld.x-m_iLastStartPointWorld.x;
	SLVVector[1]=m_iStartPointWorld.y-m_iLastStartPointWorld.y;
	SLVVector[2]=m_iStartPointWorld.z-m_iLastStartPointWorld.z;
	
	SMVector[0]=CurMVector[0];
	SMVector[1]=CurMVector[1];
	SMVector[2]=CurMVector[2];
	SVector[0]=(SLVVector[0]+SMVector[0])/2;
	SVector[1]=(SLVVector[1]+SMVector[1])/2;
	SVector[2]=(SLVVector[2]+SMVector[2])/2;


	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	float nVectorDistmm2=CurMVector[3];
	float costheta=zxh::VectorOP_Cosine( SNVVector,SVector, 3 );
	//Get the C position for the current iMinNode
	int iPosi;
	float fPointCTan[3]={0,0,0};
	FindPointCNearPosi(dfMinNodeWorld,iPosi);
	//Get the tangential direction of corresponding PointC on the model with the iMinNode on the vessel;
	GetPointCTangentVector(iPosi,fPointCTan);
	//Calclulate the distance from current iMinNode to the PointC tangent plane;
	float fPctoTanPlanDist=CalcDistFromPointCTanPlan(dfMinNodeWorld,iPosi,fPointCTan);
	//calculate the distance from current point to unique endpoint
	float CEVDistmm2=CalcDistmm2(dfMinNodeWorld,m_fSMEndPointWorld);
	
	if (CEVDistmm2<3)
	{
		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y =dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
		bIsCurSMEndPoint=true;
		*yFind=true;

	}
	else//CEVDistmm2>=1
	{

		if (m_nFMMEvlNum == 1&&fPctoTanPlanDist<1)
		{
			//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
			//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
			//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
			if (CSVDistmm2 > nVectorDistmm2)//
			{
				//first tangential direction have already initialized
				m_iEndPointWorld.x = dfMinNodeWorld.x;
				m_iEndPointWorld.y =dfMinNodeWorld.y;
				m_iEndPointWorld.z = dfMinNodeWorld.z;

				m_iSgmtPointWorld.x =dfMinNodeWorld.x;
				m_iSgmtPointWorld.y =dfMinNodeWorld.y;
				m_iSgmtPointWorld.z = dfMinNodeWorld.z;

				*yFind=true;
			}
		} 
		else  if (m_nFMMEvlNum<3&&fPctoTanPlanDist<1)
		{
			if (CSVDistmm2 > nVectorDistmm2&&costheta>0.8660&& COrgVDistmm2>=COrgMDistmm2) //theta<30бу
			{

				m_iEndPointWorld.x = dfMinNodeWorld.x;
				m_iEndPointWorld.y = dfMinNodeWorld.y;
				m_iEndPointWorld.z = dfMinNodeWorld.z;

				m_iSgmtPointWorld.x = dfMinNodeWorld.x;
				m_iSgmtPointWorld.y =dfMinNodeWorld.y;
				m_iSgmtPointWorld.z = dfMinNodeWorld.z;

				*yFind=true;
			}
		}
		else if (m_nFMMEvlNum>=3)
		{ 
			miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
			oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
			oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y = dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;


			if ( CSVDistmm2 > nVectorDistmm2&&fPctoTanPlanDist<1&&costheta>0.8660 ) //theta<30бу to determine whether is a segment point
			{  
				// is sgm point candidate, if yes return, else assign orignal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if(!IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
					*yFind=true;
				else  // segment point not Found, assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				{
					UpdateBadPointsNeighbour(iMinNode);
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}
			else if(CSVDistmm2 > nVectorDistmm2 ) // to determine whether is a bad point 
			{  
				// NOT segment point, (1) check whether bad point; (2) assign orginal values back to m_iEndPointWorld and m_iSgmtPointWorld
				if( IsSegmentBPoint(m_iSgmtPointWorld,pImageInfo,nStartPoint,sImgData ,nMaxPath))
				{
					m_zxhBadPointImg.SetImageData(iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x,m_iBadPointValue);//sABadPointImg[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x]=m_iBadPointValue;
					m_nBadpointSum++;
					UpdateBadPointsNeighbour(iMinNode);
					m_iEndPointWorld.SetValueFrom( oriEndPointWorld  ) ;  // return to original value
					m_iSgmtPointWorld.SetValueFrom( oriSgmtPointWorld  ) ;
				}
			}

		}

	}//CEVDistmm2>=1
		return m_nFmItr;

}
bool miiMinPath::UpdateBadPointsNeighbour(miiCNode<double, int> iSgmtPoint)//update the neighbourpoints of Bad points and  make a mask used to kill some points in narraw band
{
	// define a array for neighbor searching

	int gNbr[6][3] = { {-1, 0, 0}, \
				{ 1, 0, 0}, \
				{ 0,-1, 0}, \
				{ 0, 1, 0}, \
				{ 0, 0,-1}, \
				{ 0, 0, 1} };
	int nCentX = iSgmtPoint.x;
	int nCentY = iSgmtPoint.y;
	int nCentZ = iSgmtPoint.z;
	int nx, ny, nz;
	for (int i = 0; i < 6; i++)
	{
		nx = nCentX + gNbr[i][0];
		ny = nCentY + gNbr[i][1];
		nz = nCentZ + gNbr[i][2];
		if(nx<0)nx=0;
		if(nx>=m_nImgWX)nx=m_nImgWX-1;
		if(ny<0)ny=0;
		if(ny>=m_nImgWY)ny=m_nImgWY-1;
		if(nz<0)nz=0;
		if(nz>=m_nImgWZ)nz=m_nImgWZ-1;
		if(m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) 
			m_zxhBadPointMaskImg.SetPixelByGreyscale(nx,ny,nz,0,ZXH_Foreground);//sBadPointMask[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = ZXH_Foreground; 
	}
	return true;

}
bool miiMinPath::UpdateBadPointsNeighbour_WBPD(miiCNode<double, int> iSgmtPoint)//update the neighbourpoints of Bad points and  make a mask used to kill some points in narraw band
{
	// define a array for neighbor searching

	int gNbr[6][3] = { {-1, 0, 0}, \
				{ 1, 0, 0}, \
				{ 0,-1, 0}, \
				{ 0, 1, 0}, \
				{ 0, 0,-1}, \
				{ 0, 0, 1} };
	int nCentX = iSgmtPoint.x;
	int nCentY = iSgmtPoint.y;
	int nCentZ = iSgmtPoint.z;
	int nx, ny, nz;
	for (int i = 0; i < 6; i++)
	{
		nx = nCentX + gNbr[i][0];
		ny = nCentY + gNbr[i][1];
		nz = nCentZ + gNbr[i][2];
		if(nx<0)nx=0;
		if(nx>=m_nImgWX)nx=m_nImgWX-1;
		if(ny<0)ny=0;
		if(ny>=m_nImgWY)ny=m_nImgWY-1;
		if(nz<0)nz=0;
		if(nz>=m_nImgWZ)nz=m_nImgWZ-1;
		if(m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] ==FMM_FAR) 
		{
			m_zxhBadPointMaskImg.SetPixelByGreyscale(nx,ny,nz,0,ZXH_Foreground);//sBadPointMask[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = ZXH_Foreground; 
			m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] =FMM_ALIVE;
		}
		if(m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] ==FMM_TRIAL) 
		{
			for (int i = 1; i < m_vNarrowBand.size(); i++)
			{
				if (m_vNarrowBand[i].x == nx && m_vNarrowBand[i].y == ny \
					&& m_vNarrowBand[i].z == nz)
				{
                  m_iMinHeap->MinHeapDelete(m_vNarrowBand,i);
				}

			}
			
			m_zxhMinNodeImg.SetPixelByGreyscale(nx,ny,nz,0,m_iAlivePValue);
			m_zxhBadPointMaskImg.SetPixelByGreyscale(nx,ny,nz,0,ZXH_Foreground);//sBadPointMask[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = ZXH_Foreground; 
			m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_ALIVE;
		}
	}
	return true;

}
bool miiMinPath::UpdateBadPointsNeighbour_BPD(miiCNode<double, int> iSgmtPoint)//update the neighbourpoints of Bad points and  make a mask used to kill some points in narraw band
{
	// define a array for neighbor searching

	int gNbr[6][3] = { {-1, 0, 0}, \
				{ 1, 0, 0}, \
				{ 0,-1, 0}, \
				{ 0, 1, 0}, \
				{ 0, 0,-1}, \
				{ 0, 0, 1} };
	int nCentX = iSgmtPoint.x;
	int nCentY = iSgmtPoint.y;
	int nCentZ = iSgmtPoint.z;
	int nx, ny, nz;
	for (int i = 0; i < 6; i++)
	{
		nx = nCentX + gNbr[i][0];
		ny = nCentY + gNbr[i][1];
		nz = nCentZ + gNbr[i][2];
		if(nx<0)nx=0;
		if(nx>=m_nImgWX)nx=m_nImgWX-1;
		if(ny<0)ny=0;
		if(ny>=m_nImgWY)ny=m_nImgWY-1;
		if(nz<0)nz=0;
		if(nz>=m_nImgWZ)nz=m_nImgWZ-1;
		if(m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] !=FMM_ALIVE) 
			m_zxhBadPointMaskImg.SetPixelByGreyscale(nx,ny,nz,0,ZXH_Foreground);//sBadPointMask[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = ZXH_Foreground; 
	}
	return true;

}
int miiMinPath::FindImgMidPoint_DistanceThresholdPlusTriangle(miiCNode<double, int> iMinNode,bool *yFind, float CurMVector[],float CurVSLVector[])
{ //fMinNode, fMinNodelWorld
	// calculate the distance(in mm level)Add by JDQ
     miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	 
	float CSVDistmm2=CalcDistmm2(dfMinNodeWorld,m_iStartPointWorld);//distance from Current point to Startpoint on the vessel
	float LSVDistmm2=CurVSLVector[3];//distance from Last startpoint to Startpoint on the vessel
	float CLVDistmm2=CalcDistmm2(m_iLastStartPointWorld,dfMinNodeWorld);//distance from Current point to Last Startpoint on the vessel
	float nVectorDistmm2=CurMVector[3];
	float SLVector[3];//Current Segpoint  to startpoint on the vessel
	float SNVVector[3];//current model vector
	//Get the tangential direction of the m_iStartPoint
	SLVector[0]=CurVSLVector[0];
	SLVector[1]=CurVSLVector[1];
	SLVector[2]=CurVSLVector[2];
	SNVVector[0]=dfMinNodeWorld.x-m_iStartPointWorld.x;
	SNVVector[1]=dfMinNodeWorld.y-m_iStartPointWorld.y;
	SNVVector[2]=dfMinNodeWorld.z-m_iStartPointWorld.z;
	float costheta=zxh::VectorOP_Cosine( SLVector,SNVVector, 3 );
	if (m_nFMMEvlNum == 1)
	{
	//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
	//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
	//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (CSVDistmm2 > nVectorDistmm2 )//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	} 
	else
	{
		if (CSVDistmm2 > nVectorDistmm2 &&costheta<-0.8660)//theta>120 (CSVDistmm2 +LSVDistmm2-CLVDistmm2 )<=sqrt(CSVDistmm2)*sqrt(LSVDistmm2)*(-0.5)
			                                      // &&(CLVDistmm2 -(CSVDistmm2 + LSVDistmm2))/(sqrt(CSVDistmm2))*(sqrt(LSVDistmm2))>=0
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	return m_nFmItr;

}

int miiMinPath::FindImgMidPoint(miiCNode<double, int> iMinNode,bool *yFind,float nVectorDistmm2)
{ //fMinNode, fMinNodelWorld
	// calculate the distance(in mm level)Add by JDQ
    miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iMinNode);
	float nDistmm2=CalcDistmm2(m_iStartPointWorld,dfMinNodeWorld);
	//Get the tangential direction of the m_iStartPoint
	
	if (m_nFMMEvlNum == 1)
	{
	//m_fModelVecWorld[0]=m_vModelPoints[1].x-m_vModelPoints[0].x;
	//m_fModelVecWorld[1]=m_vModelPoints[1].y-m_vModelPoints[0].x;
	//m_fModelVecWorld[2]=m_vModelPoints[1].z-m_vModelPoints[0].x;
		if (nDistmm2 > nVectorDistmm2 )//
		{
			//first tangential direction have already initialized
			m_iEndPointWorld.x =dfMinNodeWorld.x;
			m_iEndPointWorld.y =dfMinNodeWorld.y;
			m_iEndPointWorld.z =dfMinNodeWorld.z;

			m_iSgmtPointWorld.x =dfMinNodeWorld.x;
			m_iSgmtPointWorld.y =dfMinNodeWorld.y;
			m_iSgmtPointWorld.z =dfMinNodeWorld.z;
			*yFind=true;
		}
	} 
	else
	{
		if (nDistmm2 > nVectorDistmm2 )
		{

			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y = dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			*yFind=true;
		}
	}
	return m_nFmItr;

}
void miiMinPath::CheckIsSgmtPointAndBadPoint(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	if (m_nFMMEvlNum == 1)
	{ /* 
	  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
	  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
	  */
		/*CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		 // m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
		//m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);

	}
	else
	{  /*
	   CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
	   m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
	   */ 
		/*CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
		m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		  // m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
		   m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
			//m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);
	}
}
void miiMinPath::CheckIsSgmtPointAndBadPointMLL(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	if (m_nFMMEvlNum == 1)
	{ 
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQ(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);

	}
	else
	{  

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQ(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
	}
}
void miiMinPath::CheckIsSgmtPointAndBadPointMLLN(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	if (m_nFMMEvlNum == 1)
	{ 
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,fCurMinNodevector,nMaxPath);

	}
	else
	{  

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData,fCurMinNodevector,nMaxPath);
	}
}

void miiMinPath::CheckIsSgmtPointAndBadPointMLLN_VSPSC(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	if (m_nFMMEvlNum == 1)
	{ 
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,fCurMinNodevector,nMaxPath);

	}
	else
	{  

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData,fCurMinNodevector,nMaxPath);
	}
}
void miiMinPath::CheckIsSgmtPointAndBadPointMLLN_VSPSC_DMP(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	if (m_nFMMEvlNum == 1)
	{ 
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_DMP(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,fCurMinNodevector,nMaxPath);

	}
	else
	{  

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_DMP(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData,fCurMinNodevector,nMaxPath);
	}
}

void miiMinPath::CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	if (m_nFMMEvlNum == 1)
	{ 
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,fCurMinNodevector,nMaxPath);

	}
	else
	{  

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData,fCurMinNodevector,nMaxPath);
	}
}
void miiMinPath::CheckIsSgmtPoint_DFM_WOP(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	if (m_nFMMEvlNum == 1)
	{ 
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=ModifiedFindImgMidPoint_DFM_WOBP(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,fCurMinNodevector,nMaxPath);

	}
	else
	{  

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		m_nFmItr=ModifiedFindImgMidPoint_DFM_WOBP(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData,fCurMinNodevector,nMaxPath);
	}
}


void miiMinPath::CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD_O(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	
	if (m_nFMMEvlNum == 1)
	{ 
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,fCurMinNodevector,nMaxPath);

	}
	else
	{  

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData,fCurMinNodevector,nMaxPath);
	}
}
void miiMinPath::CheckIsSgmtPointAndBadPointMLLSM(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,bool &blSegPointFind)
{

	float  CurMVector[4]={0};
	CurMVector[0]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	CurMVector[1]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	CurMVector[2]=m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z-m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	if (m_nFMMEvlNum == 1)
	{ /* 
	  CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
	  m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind, CurMVector,CurVSLVector);
	  */
		/*CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/
		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[0],m_vModelPointsWorld[1]);
		 // m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
		m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLLSM(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
		//m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);

	}
	else
	{  /*
	   CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
	   m_nFmItr=FindImgMidPoint_DistanceThresholdPlusTriangle(iMinNode,&blSegPointFind,CurMVector,CurVSLVector);
	   */ 
		/*CurMVector[3]=CalcDistmm2(m_iStartPointWorld,m_vModelPointsWorld[m_nModelPos]);
		m_nFmItr=FindImgMidPoint_DistanceThresholdPlusVectorAngle(iMinNode,&blSegPointFind,CurMVector);*/

		CurMVector[3]=CalcDistmm2(m_vModelPointsWorld[m_nOneSegModelVectorStartPos],m_vModelPointsWorld[m_nOneSegModelVectorEndPos]);
		  // m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
		   m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLLSM(iMinNode,&blSegPointFind,CurMVector,pImageInfo,nStartPoint,sImgData ,nMaxPath);
			//m_nFmItr=ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(iMinNode,&blSegPointFind,CurMVector);
	}
}
void miiMinPath::SmoothPoints(zxhImageData &SmoothLinePointsImg)
{
	float BackTrackDistmm=fCalcMinPathLength(m_vMinPath);
	float ZXHSegLenghmm=3;
	float meandis_of_coronaryvessel=BackTrackDistmm/(m_vMinPath.size() - 1);
	miiCNode<double,float>fSgmFPointWorld,fSgmSPointWorld;
	float ZXHSegLengthmm=2;
	
	float dist=0;
	int num=0;
	for(int k=0;k<m_vMinPath.size();k+num)
	{
		for(int i=k;i<m_vMinPath.size();i+k)
		{
			fSgmSPointWorld=ImageTransToWorldPoint(m_vMinPath[k]);
			fSgmFPointWorld=ImageTransToWorldPoint(m_vMinPath[i]);
			dist=dist+sqrt(CalcDistmm2(fSgmSPointWorld,fSgmFPointWorld));
			if(dist>ZXHSegLenghmm)
			{
				num=i;
				break;
			}
		} 
	}

}
bool miiMinPath::SelectBPasFinalSegPoint(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,int nMaxPath)//find the position of the point with maximum length
{
	cout<<"Selecting best point by backtracking"<<endl;
	if(m_nFMMEvlNum==16)
	{
		int x=0;
	}
	bool bLVP=false;//if it is the last SBP;
	int nNBSumNUM=m_vNarrowBand.size();
	if(nNumVP==nNBSumNUM)//visit all the points in nb
	{
		return false;
	}
	float fMaxDist=0;
	int nBPNum=0;
	if(nNBSumNUM-nNumVP>1000)//if the number >100 random selceting
	{
		int NUM2=fP1*nNBSumNUM;
		srand((unsigned)time(NULL));
		for(int i=0;i<=NUM2;i++)
		{
			int RNUM=rand()%nNBSumNUM;
			int s=m_nSBMap[RNUM];
			if(s==1)
			{
				RNUM=rand()%nNBSumNUM;
			}
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(m_vNarrowBand[RNUM]);
			if(s!=0)
			{
				int n=0;
			}
			if(m_nSBMap[RNUM]!=1)
			{
				float fcvlength=BackTrackInNB(m_vNarrowBand[RNUM],pImageInfo,nStartPoint,nMaxPath);
				nNumVP++;
				m_nSBMap[RNUM]=1;
				if(fcvlength<0)
				{
					continue;
				}
			
				if (fcvlength>fMaxDist)
				{
					fMaxDist=fcvlength;
					nBPNum=RNUM;
				}
			}
			
			
		}//for
		
	}
	else//number<1000 visited every point
	{
		for(int i=0;i<nNBSumNUM;i++)
		{
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(m_vNarrowBand[i]);
			
				if(m_nSBMap[i]!=1)
				{
					float fcvlength=BackTrackInNB(m_vNarrowBand[i],pImageInfo,nStartPoint,nMaxPath);
					if(fcvlength<0)
					{
						
						continue;
					}
					if (fcvlength>fMaxDist)
					{
						fMaxDist=fcvlength;
						nBPNum=i;
					}
				}

		}//for i
		bLVP=true;
	}
	int nNonezeroNUM=0;
	for(int i = 0; i <nNBSumNUM; i++)
	{
		if(m_nSBMap[i]==1)
			nNonezeroNUM++;
	}
	float fD=5;
	float fP=0.2;
	miiCNode<double, int> iNode=m_vNarrowBand[nBPNum];
	bool bP=false;
	
		
	if(SelectSP(iNode,pImageInfo,nStartPoint,2000,fD,fP,bP))//after getting the position, continue searching around the position
	{
		float fcvlength=BackTrackInNB(iNode,pImageInfo,nStartPoint,nMaxPath);
		if(fcvlength>=m_fCurrModelCArchLengthD)//if > LOD , then get the target point
		{
			m_fCurrModelCArchLengthD=fcvlength;
			miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
			oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
			oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iNode);
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y = dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			return true;
		}
		else// record the current maximum lengh and point
		{
			if(fcvlength>m_fMaxSPlength)
			{
				m_fMaxSPlength=fcvlength;
				m_miiMaxSP=iNode;
			}
			fP1=zxh::minf(fP1*1.2,1);
			if(bLVP)
			{
				return false;
			}
			else
			{
				SelectBPasFinalSegPoint(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,2000);
			}
		}
	}

	else
	{
		return false;
	}
}

bool miiMinPath::SelectBPasFinalSegPoint_RandomNoOv(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int>&input,int nMaxPath)//find the position of the point with maximum length
{
	cout<<"Selecting best point by backtracking"<<endl;
	if(m_nFMMEvlNum==16)
	{
		int x=0;
	}
	bool bLVP=false;//if it is the last SBP;
	int nNBSumNUM=m_vNarrowBand.size();
	if(input.empty())//visit all the points in nb
	{
		return false;
	}
	float fMaxDist=0;
	int nBPNum=0;
	int NUM2=fP1*nNBSumNUM;
	vector<int> vRandomNum;
	if(!GetRandom(NUM2,nNBSumNUM-nNumVP,input,vRandomNum))
		bLVP=true;
	for(int i=0;i<vRandomNum.size();i++)//get the position
	{
		int RNUM=vRandomNum[i];
		miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(m_vNarrowBand[RNUM]);
		float fcvlength=BackTrackInNB(m_vNarrowBand[RNUM],pImageInfo,nStartPoint,nMaxPath);
		nNumVP++;
		m_nSBMap[RNUM]=1;
		if(fcvlength<0)
		{
			continue;
		}

		if (fcvlength>fMaxDist)
		{
			fMaxDist=fcvlength;
			nBPNum=RNUM;
		}
	}//for
	int nNonezeroNUM=0;
	float fD=5;
	float fP=0.2;
	miiCNode<double, int> iNode=m_vNarrowBand[nBPNum];
	bool bP=false;
	
		
	if(SelectSP(iNode,pImageInfo,nStartPoint,2000,fD,fP,bP))//after getting the position, continue searching around the position
	{
		float fcvlength=BackTrackInNB(iNode,pImageInfo,nStartPoint,nMaxPath);
		if(fcvlength>=m_fCurrModelCArchLengthD)//if > LOD , then get the target point
		{
			m_fCurrModelCArchLengthD=fcvlength;
			miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
			oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
			oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iNode);
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y = dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			return true;
		}
		else// record the current maximum lengh and point
		{
			if(fcvlength>m_fMaxSPlength)
			{
				m_fMaxSPlength=fcvlength;
				m_miiMaxSP=iNode;
			}
			fP1=zxh::minf(fP1*1.2,1);
			if(bLVP)
			{
				return false;
			}
			else
			{
				return SelectBPasFinalSegPoint_RandomNoOv(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000);
			}
		}
	}

	else
	{
		return false;
	}
}
bool miiMinPath::SelectBPasFinalSegPoint_RandomNoOv15(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int>&input,int nMaxPath)//find the position of the point with maximum length
{
	cout<<"Selecting best point by backtracking"<<endl;
	if(m_nFMMEvlNum==16)
	{
		int x=0;
	}
	bool bLVP=false;//if it is the last SBP;
	int nNBSumNUM=m_vNarrowBand.size();
	if(input.empty())//visit all the points in nb
	{
		return false;
	}
	float fMaxDist=0;
	int nBPNum=0;
	int NUM2=fP1*nNBSumNUM;
	vector<int> vRandomNum;
	if(!GetRandom(NUM2,nNBSumNUM-nNumVP,input,vRandomNum))
		bLVP=true;
	for(int i=0;i<vRandomNum.size();i++)//get the position
	{
		int RNUM=vRandomNum[i];
		miiCNode<double, int>iEOUP;
		BackTrack15(m_vNarrowBand[RNUM],pImageInfo,nStartPoint,iEOUP,15);
		nNumVP++;
		m_nSBMap[RNUM]=1;
		if(iEOUP.val>0.05)
		{
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iEOUP);
			float fcvlength=BackTrackInNB(iEOUP,pImageInfo,nStartPoint,nMaxPath);

			if(fcvlength<0)
			{
				continue;
			}

			if (fcvlength>fMaxDist)
			{
				fMaxDist=fcvlength;
				nBPNum=RNUM;
			}
		}
	}//for
	int nNonezeroNUM=0;
	float fD=5;
	float fP=0.2;
	miiCNode<double, int> iNode=m_vNarrowBand[nBPNum];
	bool bP=false;
	
		
	if(SelectSP(iNode,pImageInfo,nStartPoint,2000,fD,fP,bP))//after getting the position, continue searching around the position
	{
		float fcvlength=BackTrackInNBWOL(iNode,pImageInfo,nStartPoint,nMaxPath);
		if(fcvlength>=m_fCurrModelCArchLengthD)//if > LOD , then get the target point
		{
			m_fCurrModelCArchLengthD=fcvlength;
			miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
			oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
			oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iNode);
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y = dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			return true;
		}
		else// record the current maximum lengh and point
		{
			if(fcvlength>m_fMaxSPlength)
			{
				m_fMaxSPlength=fcvlength;
				m_miiMaxSP=iNode;
			}
			fP1=zxh::minf(fP1*1.2,1);
			if(bLVP)
			{
				m_nFMMTimeS=GetTimeNow();
				return false;
			}
			else
			{
				return SelectBPasFinalSegPoint_RandomNoOv15(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000);
			}
		}
	}

	else
	{
		return false;
	}
}

bool miiMinPath::SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int>&input,int nMaxPath)//find the position of the point with maximum length
{
	cout<<"Selecting best point by backtracking"<<endl;
	if(m_nFMMEvlNum==16)
	{
		int x=0;
	}
	bool bLVP=false;//if it is the last SBP;
	int nNBSumNUM=m_vNarrowBand.size();
	if(input.empty())//visit all the points in nb
	{
		return false;
	}
	float fMaxDist=0;
	int nBPNum=0;
	int NUM2=fP1*nNBSumNUM;
	vector<int> vRandomNum;
	if(!GetRandom(NUM2,nNBSumNUM-nNumVP,input,vRandomNum))
		bLVP=true;
	for(int i=0;i<vRandomNum.size();i++)//get the position
	{
		int RNUM=vRandomNum[i];
		miiCNode<double, int>iEOUP;
		BackTrack15(m_vNarrowBand[RNUM],pImageInfo,nStartPoint,iEOUP,15);
		nNumVP++;
		m_nSBMap[RNUM]=1;
		if(iEOUP.val>0.05)
		{
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iEOUP);
			float fcvlength=BackTrackInNB(iEOUP,pImageInfo,nStartPoint,nMaxPath);

			if(fcvlength<0)
			{
				continue;
			}

			if (fcvlength>fMaxDist)
			{
				fMaxDist=fcvlength;
				nBPNum=RNUM;
			}
		}
	}//for
	int nNonezeroNUM=0;
	float fD=5;
	float fP=0.2;
	miiCNode<double, int> iNode=m_vNarrowBand[nBPNum];
	bool bP=false;
	float fcvlength=BackTrackInNBWOL(iNode,pImageInfo,nStartPoint,nMaxPath);
	cout<<"Current Length:"<<fcvlength<<";"<<"Maximum Length:"<<m_fMaxSPlength<<";"<<"LOD:"<<m_fCurrModelCArchLengthD<<endl;
 	if(fcvlength>=m_fCurrModelCArchLengthD)//if > LOD , then get the target point
	{
		m_fCurrModelCArchLengthD=fcvlength;
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 
		miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iNode);
		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
		return true;
	}
	else// record the current maximum lengh and point
	{
		if(fcvlength>m_fMaxSPlength)
		{
			m_fMaxSPlength=fcvlength;
			m_miiMaxSP=iNode;
		}
		fP1=zxh::minf(fP1*1.2,1);
		if(bLVP)
		{
			cout<<"5 minutes more FM"<<endl;
			m_nFMMTimeS=GetTimeNow();
			return false;
		}
		else
		{
			return SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000);
		}
	}
}
bool miiMinPath::SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP_WORec(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int>&input,int nMaxPath)//find the position of the point with maximum length
{
	cout<<"Selecting best point by backtracking"<<endl;
	if(m_nFMMEvlNum==16)
	{
		int x=0;
	}
	bool bLVP=false;//if it is the last SBP;
	int nNBSumNUM=m_vNarrowBand.size();
	if(input.empty())//visit all the points in nb
	{
		return false;
	}
	float fMaxDist=0;
	int nBPNum=0;
	int NUM2=fP1*nNBSumNUM;
	vector<int> vRandomNum;
	if(!GetRandom(NUM2,nNBSumNUM-nNumVP,input,vRandomNum))
		bLVP=true;
	for(int i=0;i<vRandomNum.size();i++)//get the position
	{
		int RNUM=vRandomNum[i];
		miiCNode<double, int>iEOUP;
		BackTrack15(m_vNarrowBand[RNUM],pImageInfo,nStartPoint,iEOUP,15);
		nNumVP++;
		m_nSBMap[RNUM]=1;
		if(iEOUP.val>0.5)
		{
			float fcvlength=BackTrackInNB(m_vNarrowBand[RNUM],pImageInfo,nStartPoint,nMaxPath);

			if(fcvlength<0)
			{
				continue;
			}

			if (fcvlength>fMaxDist)
			{
				fMaxDist=fcvlength;
				nBPNum=RNUM;
			}
		}
	}//for
	int nNonezeroNUM=0;
	float fD=5;
	float fP=0.2;
	miiCNode<double, int> iNode=m_vNarrowBand[nBPNum];
	bool bP=false;
	float fcvlength=BackTrackInNBWOL(iNode,pImageInfo,nStartPoint,nMaxPath);
	cout<<"Current Length:"<<fcvlength<<";"<<"Maximum Length:"<<m_fMaxSPlength<<";"<<"LOD:"<<m_fCurrModelCArchLengthD<<";"<<"WL:"<<m_ftoallength_model<<endl;
 	if(fcvlength>=m_fCurrModelCArchLengthD)//if > LOD , then get the target point
	{
		m_fCurrModelCArchLengthD=fcvlength;
		miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
		oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
		oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 
		miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iNode);
		m_iEndPointWorld.x = dfMinNodeWorld.x;
		m_iEndPointWorld.y = dfMinNodeWorld.y;
		m_iEndPointWorld.z = dfMinNodeWorld.z;

		m_iSgmtPointWorld.x = dfMinNodeWorld.x;
		m_iSgmtPointWorld.y = dfMinNodeWorld.y;
		m_iSgmtPointWorld.z = dfMinNodeWorld.z;
		return true;
	}
	else// record the current maximum lengh and point
	{

		
		if(bLVP)
		{
			cout<<"The best segmentpoint can't be find!Return the Maximal Length segmentpoint!"<<endl;

			m_nFMMTimeS=GetTimeNow();
			m_fCurrModelCArchLengthD=fcvlength;
			miiCNode<double,float> oriEndPointWorld, oriSgmtPointWorld ; 
			oriEndPointWorld.SetValueFrom( m_iEndPointWorld ) ; 
			oriSgmtPointWorld.SetValueFrom( m_iSgmtPointWorld ) ;  // store the original value 
			miiCNode<double, float>dfMinNodeWorld=ImageTransToWorldPoint(iNode);
			m_iEndPointWorld.x = dfMinNodeWorld.x;
			m_iEndPointWorld.y = dfMinNodeWorld.y;
			m_iEndPointWorld.z = dfMinNodeWorld.z;

			m_iSgmtPointWorld.x = dfMinNodeWorld.x;
			m_iSgmtPointWorld.y = dfMinNodeWorld.y;
			m_iSgmtPointWorld.z = dfMinNodeWorld.z;
			return true;
		}
		else
		{

			if(fcvlength>m_fMaxSPlength)
			{
				m_fMaxSPlength=fcvlength;
				m_miiMaxSP=iNode;
			}
			fP1=zxh::minf(fP1*1.2,1);
			return SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP_WORec(iMinNode,sImgData,pImageInfo,nStartPoint,nNumVP,fP1,input,2000);
		}
	}
}
bool miiMinPath::GetRandom(int NUM2,int nNBSumNUM,vector<int>&input,vector<int>&vRandomNum)//generate random index number without overlaping
{
	if(input.empty())
	{
		return false;
	}
	srand((int)time(NULL));

	vector<int> output = *new vector<int>();

	int end = input.size();

	do{   
		if(input.empty())
		{
			vRandomNum=output;
			return false;
		}
		vector<int>::iterator iter = input.begin();
		int num = rand()%end;
		iter = iter+num;
		output.push_back(*iter);
		input.erase(iter);
		end--;
	}while(output.size()<NUM2);
	vRandomNum=output;
	return true;
}
bool miiMinPath::SelectSP(miiCNode<double, int> &iNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,float fD, float fP,bool bP)
{
	//generate ball-like kernal for limiting the scope later on 
	
	int Bx=(fD+1)/m_nImgSpacing[0];
	int By=(fD+1)/m_nImgSpacing[1];
	int Bz=(fD+1)/m_nImgSpacing[2];
	vector<miiCNode<double,int>>vNbr;
	if(vNbr.size()>0)
	{
	vNbr.clear();
	}
	for(int ix=-Bx;ix<=Bx;ix++)
		for(int iy=-Bx;iy<=By;iy++)
			for(int iz=-Bz;iz<=Bz;iz++)
			{
				if(ix*ix+iy*iy+iz*iz<=(fD+1)*(fD+1))

				{
					miiCNode<double,int> temppoins;
					temppoins.x=ix;
					temppoins.y=iy;
					temppoins.z=iz;
					vNbr.push_back(temppoins);
				}
			}
	vector<miiCNode<double, int>> vLocNBPoiSet;
	if(vLocNBPoiSet.size()>0)
	{
		vLocNBPoiSet.clear();
	}
	return NeiBour(iNode,pImageInfo,fD,fP,bP,nStartPoint,nMaxPath,vLocNBPoiSet,vNbr);
}

bool miiMinPath::NeiBour(miiCNode<double, int> &iNode,const zxhImageInfo *pImageInfo,float fD,float fP,bool bP,float nStartPoint[3],int nMaxPath,vector<miiCNode<double, int>> &vLocNBPoiSet,vector<miiCNode<double, int>> &vNbr)
{

	vLocNBPoiSet.push_back(iNode);
	float fCenNode[3]={0};
	fCenNode[0]=(float)iNode.x;
	fCenNode[1]=(float)iNode.y;
	fCenNode[2]=(float)iNode.z;
	m_pBaseImgInfo->ImageToWorld(fCenNode);
	float lxy=sqrt(m_nImgSpacing[0]*m_nImgSpacing[0]+m_nImgSpacing[1]*m_nImgSpacing[1]);
	float lxz=sqrt(m_nImgSpacing[0]*m_nImgSpacing[0]+m_nImgSpacing[2]*m_nImgSpacing[2]);
	float lyz=sqrt(m_nImgSpacing[1]*m_nImgSpacing[1]+m_nImgSpacing[2]*m_nImgSpacing[2]);
	float minlxyxz=zxh::minf(lxy,lxz);
	float minlxyz=zxh::minf(lyz,minlxyxz);

	for (int i = 0; i < vNbr.size(); i++)
	{
		int nx = iNode.x+vNbr[i].x;
		int ny = iNode.y+vNbr[i].y;
		int nz = iNode.z+vNbr[i].z;
		if(nx<=m_nImgWX&&ny<=m_nImgWY&&nz<=m_nImgWZ)
		{
			if(nx>m_nImgWX-1)nx=m_nImgWX-1;
			if(ny>m_nImgWY-1)ny=m_nImgWY-1;
			if(nz>m_nImgWZ-1)nz=m_nImgWZ-1;
			if(nx<0)nx=0;
			if(ny<0)ny=0;
			if(nz<0)nz=0;

			float fCord[3]={0};
			fCord[0]=(float)nx;
			fCord[1]=(float)ny;
			fCord[2]=(float)nz;
			m_pBaseImgInfo->ImageToWorld(fCord);
			float dist=zxh::VectorOP_Distance(fCord,fCenNode,3);

			if(dist<fD)
			{
				if((m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx]== FMM_TRIAL)||nx==m_nImgWX-1||ny==m_nImgWY-1||nz==m_nImgWZ-1||nx==0||ny==0||nz==0)
				{
					miiCNode<double, int> itempnode;
					itempnode.x=nx;
					itempnode.y=ny;
					itempnode.z=nz;
					vLocNBPoiSet.push_back(itempnode);
				}
			}
		}
	}
	int nNBNUM=vLocNBPoiSet.size();
	float fMaxDist=0;
	int nBPNum=0;
	float fcvlength=0;
	int nSN=fP*nNBNUM;
	//find nBPNum in the "trial points"set.
	if(nNBNUM==1)
	{
		return true;
	}
	else if(nNBNUM<=500)
	{
		
		for(int i=0;i<nNBNUM;i++)
		{

			fcvlength=BackTrackInNBWOL(vLocNBPoiSet[i],pImageInfo,nStartPoint,nMaxPath);
			if(fcvlength<0)continue;
			if (fcvlength>fMaxDist)
			{
				fMaxDist=fcvlength;
				nBPNum=i;
			}
		}

		iNode=vLocNBPoiSet[nBPNum];
		return true;
	}
	else
	{
		srand((unsigned)time(NULL));
		fcvlength=BackTrackInNBWOL(vLocNBPoiSet[0],pImageInfo,nStartPoint,nMaxPath);

		if (fcvlength>fMaxDist)
		{
			fMaxDist=fcvlength;
			nBPNum=0;
		}
		for(int i=0;i<nSN;i++)//find the first position of the node
		{
			int I=rand()%nNBNUM;
			fcvlength=BackTrackInNBWOL(vLocNBPoiSet[I],pImageInfo,nStartPoint,nMaxPath);
			if(fcvlength<0)continue;
			if (fcvlength>fMaxDist)
			{
				fMaxDist=fcvlength;
				nBPNum=I;
			}
		}
	}
	miiCNode<double, int> iOldNode;
	iOldNode=iNode;
	if(fMaxDist>m_fCurrModelCArchLengthD)
	{
		iNode=vLocNBPoiSet[nBPNum];
		m_fCurrModelCArchLengthD=fMaxDist;
	}
	//judge the nBPNum point
	if(iOldNode.x==vLocNBPoiSet[nBPNum].x&&iOldNode.y==vLocNBPoiSet[nBPNum].y&&iOldNode.z==vLocNBPoiSet[nBPNum].z)
	{
		if(fP==1)
		{
			
			for(int i=0;i<nNBNUM;i++)
			{

				fcvlength=BackTrackInNBWOL(vLocNBPoiSet[i],pImageInfo,nStartPoint,nMaxPath);
				if(fcvlength<0)continue;
				if (fcvlength>fMaxDist)
				{
					fMaxDist=fcvlength;
					nBPNum=i;
				}
			}
			
				iNode=vLocNBPoiSet[nBPNum];
			
			return true;
		}
		else
		{
			fP=zxh::minf(1,fP*1.2);
			SelectSP(iNode,pImageInfo,nStartPoint,2000,fD,fP,bP);
		}
	}

	else
	{
		float lxy=sqrt(m_nImgSpacing[0]*m_nImgSpacing[0]+m_nImgSpacing[1]*m_nImgSpacing[1]);
		float lxz=sqrt(m_nImgSpacing[0]*m_nImgSpacing[0]+m_nImgSpacing[2]*m_nImgSpacing[2]);
		float lyz=sqrt(m_nImgSpacing[1]*m_nImgSpacing[1]+m_nImgSpacing[2]*m_nImgSpacing[2]);
		float minlxyxz=zxh::minf(lxy,lxz);
		float minlxyz=zxh::minf(lyz,minlxyxz);
		fD=zxh::maxf(minlxyz,fD*0.8);
		fP=zxh::minf(1,fP*1.2);
		iNode=vLocNBPoiSet[nBPNum];
		SelectSP(iNode,pImageInfo,nStartPoint,2000,fD,fP,bP);

	}
	
	
}

bool miiMinPath::FindPointC(miiCNode<double, float> fSgmtPointWorld,float fmeandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3],int *iPosi)
{
	if (m_vModelPointsWorld.size() == 0)
	{
		cerr << "The model of the coronary is null!" << endl;
		return false;
	}
	float fDist[3]={0};
	float MDistmm=0;int Nposi=0;
	float fMinDistmm=1000000;
	//find the point C based on the length of the backtrack length.
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{

		MDistmm=MDistmm+sqrt(CalcDistmm2(m_vModelPointsWorld[i],m_vModelPointsWorld[i-1]));
		if(MDistmm>=fBackTrackDistmmSkipSomePonts)
		{ 
			Nposi=i;
			break ;
		} 
	} 
	float d=CalcDistmm2(m_vModelPointsWorld[Nposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1]);
	if( Nposi==0||CalcDistmm2(m_vModelPointsWorld[Nposi],m_vModelPointsWorld[m_vModelPointsWorld.size()-1])<MINDISTTOPOINT2)//cant find the right position or very near to endpoint
	{
		*iPosi = m_vModelPointsWorld.size()-1; 
		
	}
	else
	{
		float TanVec[3]={fVesselVecWorld[0],fVesselVecWorld[1],fVesselVecWorld[2]};
		if (TanVec[0]==0&TanVec[1]==0&TanVec[2]==0)
		{
			*iPosi = Nposi ;
			return false ; // 
		}
		else
		{
			int iStartOfNposi = Nposi-DistRangeAroundCPoint/m_meandis_of_coronarymodel ; 
			int iEndOfNposi = Nposi + DistRangeAroundCPoint/m_meandis_of_coronarymodel; 
			if (iStartOfNposi<0) 
				iStartOfNposi = 0;
			if (iEndOfNposi>m_vModelPointsWorld.size()-1) 
				iEndOfNposi = m_vModelPointsWorld.size()-1; 

			for( int i=iStartOfNposi; i<=iEndOfNposi; ++i )//calculate the distance of position i to the tangent plane in a range of distance.
			{
				
				float fSgmPointWorld[3]={fSgmtPointWorld.x,fSgmtPointWorld.y,fSgmtPointWorld.z};
				//float fFSgmPointWorld[3]={m_fForwardSgmtPointWorld.x,m_fForwardSgmtPointWorld.y,m_fForwardSgmtPointWorld.z};
				float MPTTDistmm=ModelPointToTanPlaneDistmm(i,TanVec,fSgmPointWorld);

				if (MPTTDistmm<fMinDistmm)
				{
					fMinDistmm=MPTTDistmm;
					 *iPosi = i; 
				} 
			}
		}
	}
	return true;
}
bool miiMinPath::FindPointCNearPosi(miiCNode<double, float>dfMinNodeWorld,int &iPosi)
{
	int NPosi=0;
	float fMinNDistmm=1000000;
	float fSgmPointWorld[3]={dfMinNodeWorld.x,dfMinNodeWorld.y,dfMinNodeWorld.z};
	for (int i = 1; i < m_vModelPointsWorld.size(); i++)
	{
		float Fcord[3];
		Fcord[0]=m_vModelPointsWorld[i].x;
		Fcord[1]=m_vModelPointsWorld[i].y;
		Fcord[2]=m_vModelPointsWorld[i].z;
		float MPNDistmm=zxh::VectorOP_Distance(Fcord,fSgmPointWorld,3);
		if (MPNDistmm<fMinNDistmm)
		{
			fMinNDistmm=MPNDistmm;
			NPosi = i; 
		} 
	}
	iPosi=NPosi;
	return true;
}
bool miiMinPath::GetPointCTangentVector(int i,float fPointCTan[3])
{
	//Fi means in front of i
	int Fi=i-MODELVECTanLENGTH/m_meandis_of_coronarymodel;
	int Bi=i+MODELVECTanLENGTH/m_meandis_of_coronarymodel;
	Fi=CorrectTruePositionOfCoronaryModel(Fi);
	Bi=CorrectTruePositionOfCoronaryModel(Bi);
	float vModelF[3], vModelB[3];
	//the forward vector
	vModelF[0] = m_vModelPointsWorld[i].x - m_vModelPointsWorld[Fi].x;
	vModelF[1] = m_vModelPointsWorld[i].y - m_vModelPointsWorld[Fi].y;
	vModelF[2] = m_vModelPointsWorld[i].z - m_vModelPointsWorld[Fi].z;

	// the backward vector
	vModelB[0] =m_vModelPointsWorld[Bi].x-m_vModelPointsWorld[i].x;
	vModelB[1] =m_vModelPointsWorld[Bi].y-m_vModelPointsWorld[i].y;
	vModelB[2] =m_vModelPointsWorld[Bi].z-m_vModelPointsWorld[i].z;

	//the tangent vector
	fPointCTan[0]=(vModelF[0]+vModelB[0])/2;
	fPointCTan[1]=(vModelF[1]+vModelB[1])/2;
	fPointCTan[2]=(vModelF[2]+vModelB[2])/2;
	return true;
}
float miiMinPath::CalcDistFromPointCTanPlan(miiCNode<double, float>fSgmPointWorld,int i,float TanVec[3])
{
	 if(i<0)i=0;
	 if(i>m_vModelPointsWorld.size()-1)i=m_vModelPointsWorld.size()-1;
	 float NDistmm;
	 float UMToSegPointVec[3]={0,0,0};//UMToSegPointVec means the undermined vector from model point to segment point;
	 UMToSegPointVec[0]=fSgmPointWorld.x-m_vModelPointsWorld[i].x;
	 UMToSegPointVec[1]=fSgmPointWorld.y-m_vModelPointsWorld[i].y;
	 UMToSegPointVec[2]=fSgmPointWorld.z-m_vModelPointsWorld[i].z;
	 float fUNNorm=sqrt(UMToSegPointVec[0]*UMToSegPointVec[0]+UMToSegPointVec[1]*UMToSegPointVec[1]+UMToSegPointVec[2]*UMToSegPointVec[2]);
	 float fm_vVVWNorm=sqrt(TanVec[0]*TanVec[0]+TanVec[1]*TanVec[1]+TanVec[2]*TanVec[2]);
	 NDistmm=abs((UMToSegPointVec[0]*TanVec[0]+UMToSegPointVec[1]*TanVec[1]+UMToSegPointVec[2]*TanVec[2]))/( fm_vVVWNorm);

	 return NDistmm;
}
// Function Name: GetNextSegMeanInteStd()
//                 
// Parameters:  vSgmtMinPath: minimal path vector
//              sImgData: Input image intensity array
//
// Description: calculate the mean intensity ,std based on the segment point
//
// Returns: 
//
void miiMinPath::GetNextSegMeanInteStd(vector<miiCNode<double,int>> vSgmtMinPath,const short *sImgData)
{
	miiCNode<double,int>fSgmPoint,TempPoint;

	float ForwNUM=0;
	float ForwDistmm_5=5;
	ForwNUM=ForwDistmm_5/(m_meandis_of_coronaryvessel+0.1);
	float BackTrackDistmm=0;
	double IntenSum=0;
	int STNUM=vSgmtMinPath.size()-ForwNUM;
	if (STNUM<=0) STNUM=0;
	int n=0;
	for (int j=STNUM;j<vSgmtMinPath.size();j++)
	{
		TempPoint=vSgmtMinPath[j];
		short m=sImgData[TempPoint.z * m_nImgWY * m_nImgWX + TempPoint.y * m_nImgWX + TempPoint.x];
		if(m<0) m=0.0;
		IntenSum=IntenSum+m;
		n++;

	}
	m_fMean=IntenSum/n;
	double fUnseenAortaStdDev=m_dUnseenStartStdDev;
	m_fStdDev=fUnseenAortaStdDev;
}

bool miiMinPath::BoundaryCorrect(int *PointPos,int ImgNewVslsSize[4])
{
	for (int i=0;i<3;i++)
	{
	PointPos[i]=zxh::maxf(0,PointPos[i]);
	PointPos[i]=zxh::minf(ImgNewVslsSize[i]-1,PointPos[i]);
	}
	return true;
}
void miiMinPath::MapCurrentVectortoImgandOutput(zxhImageDataT<short> &zxhModelpointsImg)
{
	int zxhModelpointsImgSize[4]={1};
	zxhModelpointsImg.GetImageSize(zxhModelpointsImgSize[0],zxhModelpointsImgSize[1],zxhModelpointsImgSize[2],zxhModelpointsImgSize[3]);
	for(int it=0;it<zxhModelpointsImgSize[3];++it)
		for(int iz=0;iz<zxhModelpointsImgSize[2];++iz)
			for(int iy=0;iy<zxhModelpointsImgSize[1];++iy)
				for(int ix=0;ix<zxhModelpointsImgSize[0];++ix)
				{
					zxhModelpointsImg.SetPixelByGreyscale(ix,iy,iz,it,0);
				}

				for (int i=0;i<m_vModelPointsWorld.size();i++)
				{
					float PointPosWorld[ImageDimensionMax]={0};
					int PointPos[4]={0};
					PointPosWorld[0]=m_vModelPointsWorld[i].x;
					PointPosWorld[1]=m_vModelPointsWorld[i].y;
					PointPosWorld[2]=m_vModelPointsWorld[i].z;
					zxhModelpointsImg.GetImageInfo()->WorldToImage(PointPosWorld);
					PointPos[0]=zxh::round(PointPosWorld[0]);
					PointPos[1]=(int)(PointPosWorld[1]+0.5);
					PointPos[2]=(int)(PointPosWorld[2]+0.5);
					BoundaryCorrect(PointPos,zxhModelpointsImgSize); 
					int ZXHModelPointsInte=15;
					if(i>=m_nOneSegModelVectorStartPos&&i<=m_nOneSegModelVectorEndPos)

						zxhModelpointsImg.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXHModelPointsInte);
					else
						zxhModelpointsImg.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);
				}
				
}
void miiMinPath::CoutSandEndPosi()//cout the position of the start point and end point
{
	ofstream PosiFile;
	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt")+1;//error occur in matlab script if there is not "+1"
	char *chFileName = new char[nLen];
	strcpy(chFileName, m_cResultPath);
	strcat(chFileName, "/Posi");
	strcat(chFileName, ".txt");
	PosiFile.open(chFileName,ios::app);
	if(!PosiFile)
	{cerr<<"open error!"<<endl;
	exit(1);
	}
	PosiFile<<"VisitedPointNUM:"<<m_nFMVedPoiNUM<<"CPosi:"<<m_nOneSegModelVectorStartPos<<";"<<"DPosi:"<<m_nOneSegModelVectorEndPos<<";"<<"Diff:"<<m_nOneSegModelVectorEndPos-m_nOneSegModelVectorStartPos<<";"<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";
	PosiFile.close();
	delete[] chFileName;

	char chTemp[25];
	char chTempNum[25];
	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
	_itoa_s(m_nFMMEvlNum-1, chTempNum, 10);
	//output the narrawband points
   int nLen4 = strlen(m_cResultPath) + strlen("/VistPons") + strlen(chTemp) +strlen("_")+strlen(chTempNum)+ strlen(".nii.gz") + 1;
	char *chFileName4 = (char *)malloc(nLen4);
	strcpy(chFileName4, m_cResultPath);
	strcat(chFileName4, "/VistPons");
	strcat(chFileName4, chTemp);
	strcat(chFileName4, "_");
	strcat(chFileName4, chTempNum);
	strcat(chFileName4, ".nii.gz");
	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);
	////SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////
	//////int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName9 = (char *)malloc(nLen9);
	//////strcpy(chFileName9, m_cResultPath);
	//////strcat(chFileName9, "/NBSPeedValue");
	//////strcat(chFileName9, chTemp);
	//////strcat(chFileName9, ".nii.gz");
	//////zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//////free(chFileName9);
	//////
	//////int noneZeroNUM=0;
	//////for(int nx=0;nx<m_nImgWX;nx++)
	//////	for(int ny=0;ny<m_nImgWY;ny++)
	//////		for(int nz=0;nz<m_nImgWZ;nz++)
	//////		{
	//////			short tempint=m_zxhMinNodeImg.GetPixelGreyscale(nx,ny,nz,0);
	//////			if (tempint!=0)
	//////				noneZeroNUM++;
	//////		}
	//////string strFileName4=chFileName4;
	//////cout<<"None-Zero Voxel Number in NB:"<<noneZeroNUM<<endl;
	//////free(chFileName4);
	////////
	////////output the model vector as an nii.gz image
	//////zxhImageDataT<short> zxhModelpointsImg;
	//////zxhModelpointsImg.NewImage(m_pBaseImgInfo);
	//////MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
	////////OUTPUT MODEL VECTOR AS IMAGE
	//////
	//////int nLen3 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName3 = (char *)malloc(nLen3);
	//////strcpy(chFileName3, m_cResultPath);
	//////strcat(chFileName3, "/ModelVector");
	//////strcat(chFileName3, chTemp);
	//////strcat(chFileName3, ".nii.gz");
	//////zxh::SaveImage( &zxhModelpointsImg, chFileName3 );
	//////free(chFileName3);
	//////
	//
	//	double fm=m_zxhNBUvalueImg.GetPixelGreyscale(50,49,79,0);
	//	m_zxhNBUvalueImg.SetPixelByGreyscale(50,49,79,0,fm);
	//	double um=m_dU[79 * m_nImgWY * m_nImgWX + 49 * m_nImgWX + 50];
	//int nLen5 = strlen(m_cResultPath) + strlen("/UV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//char *chFileName5 = (char *)malloc(nLen5);
	//strcpy(chFileName5, m_cResultPath);
	//strcat(chFileName5, "/UV");
	//strcat(chFileName5, chTemp);
	//strcat(chFileName5, ".nii.gz");
	//string strchFileName=chFileName5;
	//zxh::SaveImage<float>(&m_zxhNBUvalueImg,strchFileName);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//free(chFileName5);
	/*
	int nLen6 = strlen(m_cResultPath) + strlen("/NBPV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName6 = (char *)malloc(nLen6);
	strcpy(chFileName6, m_cResultPath);
	strcat(chFileName6, "/NBPV");
	strcat(chFileName6, chTemp);
	strcat(chFileName6, ".nii.gz");
	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName6);*/
	////
	////int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	////char *chFileName7 = (char *)malloc(nLen7);
	////strcpy(chFileName7, m_cResultPath);
	////strcat(chFileName7, "/NBSValue");
	////strcat(chFileName7, chTemp);
	////strcat(chFileName7, ".nii.gz");
	////zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////free(chFileName7);
	////int nLen8 = strlen(m_cResultPath) + strlen("/NBDV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	////char *chFileName8 = (char *)malloc(nLen8);
	////strcpy(chFileName8, m_cResultPath);
	////strcat(chFileName8, "/NBDV");
	////strcat(chFileName8, chTemp);
	////strcat(chFileName8, ".nii.gz");
	////zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////free(chFileName8);
	//
 //   int nLen9 = strlen(m_cResultPath) + strlen("/NBD_Segvalue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//char *chFileName9 = (char *)malloc(nLen9);
	//strcpy(chFileName9, m_cResultPath);
	//strcat(chFileName9, "/NBD_Segvalue");
	//strcat(chFileName9, chTemp);
	//strcat(chFileName9, ".nii.gz");
	//zxh::SaveImage( &m_zxhNBD_SegvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//free(chFileName9);

	//zxhModelpointsImg.ReleaseMem();



}
void miiMinPath::CoutSandEndPosi2()//cout the position of the start point and end point
{
	ofstream PosiFile;
	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt")+1;//error occur in matlab script if there is not "+1"
	char *chFileName = new char[nLen];
	strcpy(chFileName, m_cResultPath);
	strcat(chFileName, "/Posi");
	strcat(chFileName, ".txt");
	PosiFile.open(chFileName,ios::app);
	if(!PosiFile)
	{cerr<<"open error!"<<endl;
	exit(1);
	}
	PosiFile<<"VisitedPointNUM:"<<m_nFMVedPoiNUM<<"CPosi:"<<m_nOneSegModelVectorStartPos<<";"<<"DPosi:"<<m_nOneSegModelVectorEndPos<<";"<<"Diff:"<<m_nOneSegModelVectorEndPos-m_nOneSegModelVectorStartPos<<";"<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";
	PosiFile.close();
	delete[] chFileName;

	//char chTemp[25];
	//char chTempNum[25];
	//_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
	//_itoa_s(m_nFMMEvlNum-1, chTempNum, 10);
	//output the narrawband points
 //  int nLen4 = strlen(m_cResultPath) + strlen("/VistPons") + strlen(chTemp) +strlen("_")+strlen(chTempNum)+ strlen(".nii.gz") + 1;
	//char *chFileName4 = (char *)malloc(nLen4);
	//strcpy(chFileName4, m_cResultPath);
	//strcat(chFileName4, "/VistPons");
	//strcat(chFileName4, chTemp);
	//strcat(chFileName4, "_");
	//strcat(chFileName4, chTempNum);
	//strcat(chFileName4, ".nii.gz");
	//zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);
	////SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////
	//////int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName9 = (char *)malloc(nLen9);
	//////strcpy(chFileName9, m_cResultPath);
	//////strcat(chFileName9, "/NBSPeedValue");
	//////strcat(chFileName9, chTemp);
	//////strcat(chFileName9, ".nii.gz");
	//////zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//////free(chFileName9);
	//////
	//////int noneZeroNUM=0;
	//////for(int nx=0;nx<m_nImgWX;nx++)
	//////	for(int ny=0;ny<m_nImgWY;ny++)
	//////		for(int nz=0;nz<m_nImgWZ;nz++)
	//////		{
	//////			short tempint=m_zxhMinNodeImg.GetPixelGreyscale(nx,ny,nz,0);
	//////			if (tempint!=0)
	//////				noneZeroNUM++;
	//////		}
	//////string strFileName4=chFileName4;
	//////cout<<"None-Zero Voxel Number in NB:"<<noneZeroNUM<<endl;
	//////free(chFileName4);
	////////
	////////output the model vector as an nii.gz image
	//////zxhImageDataT<short> zxhModelpointsImg;
	//////zxhModelpointsImg.NewImage(m_pBaseImgInfo);
	//////MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
	////////OUTPUT MODEL VECTOR AS IMAGE
	//////
	//////int nLen3 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName3 = (char *)malloc(nLen3);
	//////strcpy(chFileName3, m_cResultPath);
	//////strcat(chFileName3, "/ModelVector");
	//////strcat(chFileName3, chTemp);
	//////strcat(chFileName3, ".nii.gz");
	//////zxh::SaveImage( &zxhModelpointsImg, chFileName3 );
	//////free(chFileName3);
	//////
	//
	//	double fm=m_zxhNBUvalueImg.GetPixelGreyscale(50,49,79,0);
	//	m_zxhNBUvalueImg.SetPixelByGreyscale(50,49,79,0,fm);
	//	double um=m_dU[79 * m_nImgWY * m_nImgWX + 49 * m_nImgWX + 50];
	//int nLen5 = strlen(m_cResultPath) + strlen("/UV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//char *chFileName5 = (char *)malloc(nLen5);
	//strcpy(chFileName5, m_cResultPath);
	//strcat(chFileName5, "/UV");
	//strcat(chFileName5, chTemp);
	//strcat(chFileName5, ".nii.gz");
	//string strchFileName=chFileName5;
	//zxh::SaveImage<float>(&m_zxhNBUvalueImg,strchFileName);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//free(chFileName5);
	/*
	int nLen6 = strlen(m_cResultPath) + strlen("/NBPV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName6 = (char *)malloc(nLen6);
	strcpy(chFileName6, m_cResultPath);
	strcat(chFileName6, "/NBPV");
	strcat(chFileName6, chTemp);
	strcat(chFileName6, ".nii.gz");
	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName6);*/
	////
	////int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	////char *chFileName7 = (char *)malloc(nLen7);
	////strcpy(chFileName7, m_cResultPath);
	////strcat(chFileName7, "/NBSValue");
	////strcat(chFileName7, chTemp);
	////strcat(chFileName7, ".nii.gz");
	////zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////free(chFileName7);
	////int nLen8 = strlen(m_cResultPath) + strlen("/NBDV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	////char *chFileName8 = (char *)malloc(nLen8);
	////strcpy(chFileName8, m_cResultPath);
	////strcat(chFileName8, "/NBDV");
	////strcat(chFileName8, chTemp);
	////strcat(chFileName8, ".nii.gz");
	////zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////free(chFileName8);
	//
 //   int nLen9 = strlen(m_cResultPath) + strlen("/NBD_Segvalue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//char *chFileName9 = (char *)malloc(nLen9);
	//strcpy(chFileName9, m_cResultPath);
	//strcat(chFileName9, "/NBD_Segvalue");
	//strcat(chFileName9, chTemp);
	//strcat(chFileName9, ".nii.gz");
	//zxh::SaveImage( &m_zxhNBD_SegvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//free(chFileName9);

	//zxhModelpointsImg.ReleaseMem();



}
void miiMinPath::CoutSandEndPosiV()//cout the position of the start point and end point
{
	ofstream PosiFile;
	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt")+1;//error occur in matlab script if there is not "+1"
	char *chFileName = new char[nLen];
	strcpy(chFileName, m_cResultPath);
	strcat(chFileName, "/Posi");
	strcat(chFileName, ".txt");
	PosiFile.open(chFileName,ios::app);
	if(!PosiFile)
	{cerr<<"open error!"<<endl;
	exit(1);
	}
	PosiFile<<"VisitedPointNUM:"<<m_nFMVedPoiNUM<<"CPosi:"<<m_nOneSegModelVectorStartPos<<";"<<"DPosi:"<<m_nOneSegModelVectorEndPos<<";"<<"Diff:"<<m_nOneSegModelVectorEndPos-m_nOneSegModelVectorStartPos<<";"<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";
	PosiFile.close();
	delete[] chFileName;

	/*char chTemp[25];
	char chTempNum[25];
	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
	_itoa_s(m_nFMMEvlNum-1, chTempNum, 10);
	//output the narrawband points
   int nLen4 = strlen(m_cResultPath) + strlen("/VistPons") + strlen(chTemp) +strlen("_")+strlen(chTempNum)+ strlen(".nii.gz") + 1;
	char *chFileName4 = (char *)malloc(nLen4);
	strcpy(chFileName4, m_cResultPath);
	strcat(chFileName4, "/VistPons");
	strcat(chFileName4, chTemp);
	strcat(chFileName4, "_");
	strcat(chFileName4, chTempNum);
	strcat(chFileName4, ".nii.gz");
	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);*/
	////SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////
	//////int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName9 = (char *)malloc(nLen9);
	//////strcpy(chFileName9, m_cResultPath);
	//////strcat(chFileName9, "/NBSPeedValue");
	//////strcat(chFileName9, chTemp);
	//////strcat(chFileName9, ".nii.gz");
	//////zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//////free(chFileName9);
	//////
	//////int noneZeroNUM=0;
	//////for(int nx=0;nx<m_nImgWX;nx++)
	//////	for(int ny=0;ny<m_nImgWY;ny++)
	//////		for(int nz=0;nz<m_nImgWZ;nz++)
	//////		{
	//////			short tempint=m_zxhMinNodeImg.GetPixelGreyscale(nx,ny,nz,0);
	//////			if (tempint!=0)
	//////				noneZeroNUM++;
	//////		}
	//////string strFileName4=chFileName4;
	//////cout<<"None-Zero Voxel Number in NB:"<<noneZeroNUM<<endl;
	//////free(chFileName4);
	////////
	////////output the model vector as an nii.gz image
	//////zxhImageDataT<short> zxhModelpointsImg;
	//////zxhModelpointsImg.NewImage(m_pBaseImgInfo);
	//////MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
	////////OUTPUT MODEL VECTOR AS IMAGE
	//////
	//////int nLen3 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName3 = (char *)malloc(nLen3);
	//////strcpy(chFileName3, m_cResultPath);
	//////strcat(chFileName3, "/ModelVector");
	//////strcat(chFileName3, chTemp);
	//////strcat(chFileName3, ".nii.gz");
	//////zxh::SaveImage( &zxhModelpointsImg, chFileName3 );
	//////free(chFileName3);
	//////
	//
	//	double fm=m_zxhNBUvalueImg.GetPixelGreyscale(50,49,79,0);
	//	m_zxhNBUvalueImg.SetPixelByGreyscale(50,49,79,0,fm);
	//	double um=m_dU[79 * m_nImgWY * m_nImgWX + 49 * m_nImgWX + 50];
	/*int nLen5 = strlen(m_cResultPath) + strlen("/UV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName5 = (char *)malloc(nLen5);
	strcpy(chFileName5, m_cResultPath);
	strcat(chFileName5, "/UV");
	strcat(chFileName5, chTemp);
	strcat(chFileName5, ".nii.gz");
	string strchFileName=chFileName5;
	zxh::SaveImage<double>(&m_zxhNBUvalueImg,strchFileName);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName5);
	int nLen6 = strlen(m_cResultPath) + strlen("/NBPV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName6 = (char *)malloc(nLen6);
	strcpy(chFileName6, m_cResultPath);
	strcat(chFileName6, "/NBPV");
	strcat(chFileName6, chTemp);
	strcat(chFileName6, ".nii.gz");
	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName6);*/
	////
	////int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	////char *chFileName7 = (char *)malloc(nLen7);
	////strcpy(chFileName7, m_cResultPath);
	////strcat(chFileName7, "/NBSValue");
	////strcat(chFileName7, chTemp);
	////strcat(chFileName7, ".nii.gz");
	////zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////free(chFileName7);
	////int nLen8 = strlen(m_cResultPath) + strlen("/NBDV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	////char *chFileName8 = (char *)malloc(nLen8);
	////strcpy(chFileName8, m_cResultPath);
	////strcat(chFileName8, "/NBDV");
	////strcat(chFileName8, chTemp);
	////strcat(chFileName8, ".nii.gz");
	////zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////free(chFileName8);
	//
 //   int nLen9 = strlen(m_cResultPath) + strlen("/NBD_Segvalue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//char *chFileName9 = (char *)malloc(nLen9);
	//strcpy(chFileName9, m_cResultPath);
	//strcat(chFileName9, "/NBD_Segvalue");
	//strcat(chFileName9, chTemp);
	//strcat(chFileName9, ".nii.gz");
	//zxh::SaveImage( &m_zxhNBD_SegvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//free(chFileName9);

	//zxhModelpointsImg.ReleaseMem();



}
void miiMinPath::CoutSandEndPosi1()//cout the position of the start point and end point
{
	ofstream PosiFile;
	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt")+1;//error occur in matlab script if there is not "+1"
	char *chFileName = new char[nLen];
	strcpy(chFileName, m_cResultPath);
	strcat(chFileName, "/Posi");
	strcat(chFileName, ".txt");
	PosiFile.open(chFileName,ios::app);
	if(!PosiFile)
	{cerr<<"open error!"<<endl;
	exit(1);
	}
	PosiFile<<"VisitedPointNUM:"<<m_nFMVedPoiNUM<<"CPosi:"<<m_nOneSegModelVectorStartPos<<";"<<"DPosi:"<<m_nOneSegModelVectorEndPos<<";"<<"Diff:"<<m_nOneSegModelVectorEndPos-m_nOneSegModelVectorStartPos<<";"<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";
	PosiFile.close();
	delete[] chFileName;

	char chTemp[25];
	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
	//output the narrawband points
   int nLen4 = strlen(m_cResultPath) + strlen("/VistPons") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName4 = (char *)malloc(nLen4);
	strcpy(chFileName4, m_cResultPath);
	strcat(chFileName4, "/VistPons");
	strcat(chFileName4, chTemp);
	strcat(chFileName4, ".nii.gz");
	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);////SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	////
	//////int nLen9 = strlen(m_cResultPath) + strlen("/NBSPeedValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName9 = (char *)malloc(nLen9);
	//////strcpy(chFileName9, m_cResultPath);
	//////strcat(chFileName9, "/NBSPeedValue");
	//////strcat(chFileName9, chTemp);
	//////strcat(chFileName9, ".nii.gz");
	//////zxh::SaveImage( &m_zxhNBSpeedvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//////free(chFileName9);
	//////
	//////int noneZeroNUM=0;
	//////for(int nx=0;nx<m_nImgWX;nx++)
	//////	for(int ny=0;ny<m_nImgWY;ny++)
	//////		for(int nz=0;nz<m_nImgWZ;nz++)
	//////		{
	//////			short tempint=m_zxhMinNodeImg.GetPixelGreyscale(nx,ny,nz,0);
	//////			if (tempint!=0)
	//////				noneZeroNUM++;
	//////		}
	//////string strFileName4=chFileName4;
	//////cout<<"None-Zero Voxel Number in NB:"<<noneZeroNUM<<endl;
	//////free(chFileName4);
	////////
	////////output the model vector as an nii.gz image
	//////zxhImageDataT<short> zxhModelpointsImg;
	//////zxhModelpointsImg.NewImage(m_pBaseImgInfo);
	//////MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
	////////OUTPUT MODEL VECTOR AS IMAGE
	//////
	//////int nLen3 = strlen(m_cResultPath) + strlen("/ModelVector") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//////char *chFileName3 = (char *)malloc(nLen3);
	//////strcpy(chFileName3, m_cResultPath);
	//////strcat(chFileName3, "/ModelVector");
	//////strcat(chFileName3, chTemp);
	//////strcat(chFileName3, ".nii.gz");
	//////zxh::SaveImage( &zxhModelpointsImg, chFileName3 );
	//////free(chFileName3);
	//////
	//
	//	double fm=m_zxhNBUvalueImg.GetPixelGreyscale(50,49,79,0);
	//	m_zxhNBUvalueImg.SetPixelByGreyscale(50,49,79,0,fm);
	//	double um=m_dU[79 * m_nImgWY * m_nImgWX + 49 * m_nImgWX + 50];
	int nLen5 = strlen(m_cResultPath) + strlen("/UV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName5 = (char *)malloc(nLen5);
	strcpy(chFileName5, m_cResultPath);
	strcat(chFileName5, "/UV");
	strcat(chFileName5, chTemp);
	strcat(chFileName5, ".nii.gz");
	string strchFileName=chFileName5;
	zxh::SaveImage(&m_zxhNBUvalueImg,strchFileName);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName5);
	int nLen6 = strlen(m_cResultPath) + strlen("/NBPV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName6 = (char *)malloc(nLen6);
	strcpy(chFileName6, m_cResultPath);
	strcat(chFileName6, "/NBPV");
	strcat(chFileName6, chTemp);
	strcat(chFileName6, ".nii.gz");
	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName6);
	////
	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName7 = (char *)malloc(nLen7);
	strcpy(chFileName7, m_cResultPath);
	strcat(chFileName7, "/NBSValue");
	strcat(chFileName7, chTemp);
	strcat(chFileName7, ".nii.gz");
	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName7);
	int nLen8 = strlen(m_cResultPath) + strlen("/NBDV") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName8 = (char *)malloc(nLen8);
	strcpy(chFileName8, m_cResultPath);
	strcat(chFileName8, "/NBDV");
	strcat(chFileName8, chTemp);
	strcat(chFileName8, ".nii.gz");
	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName8);
	//
 //   int nLen9 = strlen(m_cResultPath) + strlen("/NBD_Segvalue") + strlen(chTemp) + strlen(".nii.gz") + 1;
	//char *chFileName9 = (char *)malloc(nLen9);
	//strcpy(chFileName9, m_cResultPath);
	//strcat(chFileName9, "/NBD_Segvalue");
	//strcat(chFileName9, chTemp);
	//strcat(chFileName9, ".nii.gz");
	//zxh::SaveImage( &m_zxhNBD_SegvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	//free(chFileName9);

	//zxhModelpointsImg.ReleaseMem();



}
void miiMinPath::CoutSandEndPosiE()//cout the position of the start point and end point
{
	ofstream PosiFile;
	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt")+1;//error occur in matlab script if there is not "+1"
	char *chFileName = new char[nLen];
	strcpy(chFileName, m_cResultPath);
	strcat(chFileName, "/Posi");
	strcat(chFileName, ".txt");
	PosiFile.open(chFileName,ios::app);
	if(!PosiFile)
	{cerr<<"open error!"<<endl;
	exit(1);
	}
	PosiFile<<"VisitedPointNUM:"<<m_nFMVedPoiNUM<<"CPosi:"<<m_nOneSegModelVectorStartPos<<";"<<"DPosi:"<<m_nOneSegModelVectorEndPos<<";"<<"Diff:"<<m_nOneSegModelVectorEndPos-m_nOneSegModelVectorStartPos<<";"<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";
	PosiFile.close();
	delete[] chFileName;
	char chTemp[25];
	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
	//output the narrawband points
	int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath_E") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName4 = (char *)malloc(nLen4);
	strcpy(chFileName4, m_cResultPath);
	strcat(chFileName4, "/MinNodeNBMinPath_E");
	strcat(chFileName4, chTemp);
	strcat(chFileName4, ".nii.gz");
	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName4);

	//
	//output the model vector as an nii.gz image
	zxhImageDataT<short> zxhModelpointsImg;
	zxhModelpointsImg.NewImage(m_pBaseImgInfo);
	MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
	//OUTPUT MODEL VECTOR AS IMAGE
	
	int nLen3 = strlen(m_cResultPath) + strlen("/ModelVector_E") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName3 = (char *)malloc(nLen3);
	strcpy(chFileName3, m_cResultPath);
	strcat(chFileName3, "/ModelVector_E");
	strcat(chFileName3, chTemp);
	strcat(chFileName3, ".nii.gz");
	zxh::SaveImage( &zxhModelpointsImg, chFileName3 );
	free(chFileName3);
	
	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue_E") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName5 = (char *)malloc(nLen5);
	strcpy(chFileName5, m_cResultPath);
	strcat(chFileName5, "/NBUValue_E");
	strcat(chFileName5, chTemp);
	strcat(chFileName5, ".nii.gz");
	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName5);
	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue_E") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName6 = (char *)malloc(nLen6);
	strcpy(chFileName6, m_cResultPath);
	strcat(chFileName6, "/NBPValue_E");
	strcat(chFileName6, chTemp);
	strcat(chFileName6, ".nii.gz");
	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName6);
	
	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue_E") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName7 = (char *)malloc(nLen7);
	strcpy(chFileName7, m_cResultPath);
	strcat(chFileName7, "/NBSValue_E");
	strcat(chFileName7, chTemp);
	strcat(chFileName7, ".nii.gz");
	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName7);
	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue_E") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName8 = (char *)malloc(nLen8);
	strcpy(chFileName8, m_cResultPath);
	strcat(chFileName8, "/NBDValue_E");
	strcat(chFileName8, chTemp);
	strcat(chFileName8, ".nii.gz");
	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName8);
	int nLen9 = strlen(m_cResultPath) + strlen("/NBD_Segvalue_E") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName9 = (char *)malloc(nLen9);
	strcpy(chFileName9, m_cResultPath);
	strcat(chFileName9, "/NBD_Segvalue_E");
	strcat(chFileName9, chTemp);
	strcat(chFileName9, ".nii.gz");
	zxh::SaveImage( &m_zxhNBD_SegvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName9);
	zxhModelpointsImg.ReleaseMem();
}
void miiMinPath::CoutSandEndPosiL()//cout the position of the start point and end point
{
	ofstream PosiFile;
	int nLen = strlen(m_cResultPath) + strlen("/Posi") +strlen(".txt")+1;//error occur in matlab script if there is not "+1"
	char *chFileName = new char[nLen];
	strcpy(chFileName, m_cResultPath);
	strcat(chFileName, "/Posi");
	strcat(chFileName, ".txt");
	PosiFile.open(chFileName,ios::app);
	if(!PosiFile)
	{cerr<<"open error!"<<endl;
	exit(1);
	}
	PosiFile<<"VisitedPointNUM:"<<m_nFMVedPoiNUM<<"CPosi:"<<m_nOneSegModelVectorStartPos<<";"<<"DPosi:"<<m_nOneSegModelVectorEndPos<<";"<<"Diff:"<<m_nOneSegModelVectorEndPos-m_nOneSegModelVectorStartPos<<";"<<"m_fMean:"<<m_fMean<<";"<<"m_fStdDev:"<<m_fStdDev<<"\n";
	PosiFile.close();
	delete[] chFileName;
	char chTemp[25];
	_itoa_s(m_nFMVedPoiNUM, chTemp, 10);
	//output the narrawband points
	int nLen4 = strlen(m_cResultPath) + strlen("/MinNodeNBMinPath_L") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName4 = (char *)malloc(nLen4);
	strcpy(chFileName4, m_cResultPath);
	strcat(chFileName4, "/MinNodeNBMinPath_L");
	strcat(chFileName4, chTemp);
	strcat(chFileName4, ".nii.gz");
	zxh::SaveImage( &m_zxhMinNodeImg, chFileName4);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName4);

	//
	//output the model vector as an nii.gz image
	zxhImageDataT<short> zxhModelpointsImg;
	zxhModelpointsImg.NewImage(m_pBaseImgInfo);
	MapCurrentVectortoImgandOutput(zxhModelpointsImg);//map modelpoints to an new image as rawimage,and mark current vetor and then output.
	//OUTPUT MODEL VECTOR AS IMAGE
	
	int nLen3 = strlen(m_cResultPath) + strlen("/ModelVector_L") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName3 = (char *)malloc(nLen3);
	strcpy(chFileName3, m_cResultPath);
	strcat(chFileName3, "/ModelVector_L");
	strcat(chFileName3, chTemp);
	strcat(chFileName3, ".nii.gz");
	zxh::SaveImage( &zxhModelpointsImg, chFileName3 );
	free(chFileName3);
	
	int nLen5 = strlen(m_cResultPath) + strlen("/NBUValue_L") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName5 = (char *)malloc(nLen5);
	strcpy(chFileName5, m_cResultPath);
	strcat(chFileName5, "/NBUValue_L");
	strcat(chFileName5, chTemp);
	strcat(chFileName5, ".nii.gz");
	zxh::SaveImage( &m_zxhNBUvalueImg, chFileName5);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName5);
	int nLen6 = strlen(m_cResultPath) + strlen("/NBPValue_L") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName6 = (char *)malloc(nLen6);
	strcpy(chFileName6, m_cResultPath);
	strcat(chFileName6, "/NBPValue_L");
	strcat(chFileName6, chTemp);
	strcat(chFileName6, ".nii.gz");
	zxh::SaveImage( &m_zxhNBPvalueImg, chFileName6);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName6);
	
	int nLen7 = strlen(m_cResultPath) + strlen("/NBSValue_L") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName7 = (char *)malloc(nLen7);
	strcpy(chFileName7, m_cResultPath);
	strcat(chFileName7, "/NBSValue_L");
	strcat(chFileName7, chTemp);
	strcat(chFileName7, ".nii.gz");
	zxh::SaveImage( &m_zxhNBSvalueImg, chFileName7);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName7);
	int nLen8 = strlen(m_cResultPath) + strlen("/NBDValue_L") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName8 = (char *)malloc(nLen8);
	strcpy(chFileName8, m_cResultPath);
	strcat(chFileName8, "/NBDValue_L");
	strcat(chFileName8, chTemp);
	strcat(chFileName8, ".nii.gz");
	zxh::SaveImage( &m_zxhNBDvalueImg, chFileName8);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName8);
	int nLen9 = strlen(m_cResultPath) + strlen("/NBD_Segvalue_L") + strlen(chTemp) + strlen(".nii.gz") + 1;
	char *chFileName9 = (char *)malloc(nLen9);
	strcpy(chFileName9, m_cResultPath);
	strcat(chFileName9, "/NBD_Segvalue_L");
	strcat(chFileName9, chTemp);
	strcat(chFileName9, ".nii.gz");
	zxh::SaveImage( &m_zxhNBD_SegvalueImg, chFileName9);//SaveImage(m_pBaseImgInfo, sAMinNodeImg, chFileName4);
	free(chFileName9);
	zxhModelpointsImg.ReleaseMem();
}
void miiMinPath::SmoothPath()
{
	float fStepDist=m_fBackTrackDistmm/100;


	if (m_vMinPathSegWorld.size() > 0)
	{
		m_vMinPathSegWorld.clear();
	}
	if (m_vMinPathSeg.size() > 0)
	{
		m_vMinPathSeg.clear();
	}
	m_vMinPathSegWorld.push_back(m_vMinPathWorld[0]);
	m_vMinPathSeg.push_back(m_vMinPath[0]);
	int idistNum_2mm=2/m_meandis_of_coronaryvessel;
	for(int i=1;i<101;i++)//100 seg point
	{
		int SegNUM=0;
		float dist=0;
		for(int j=0;j<m_vMinPathWorld.size();j++)
		{
			dist=dist+sqrt(CalcDistmm2(m_vMinPathWorld[j],m_vMinPathWorld[j+1]));
			if(dist>=i*fStepDist)
			{
				SegNUM=j;
				break;
			}
		}

		int SegStart=SegNUM-idistNum_2mm;
		int SegEnd=SegNUM+idistNum_2mm;
		if (SegStart<=0)SegStart=0;
		if (SegEnd>=0)SegEnd=m_vMinPathWorld.size()-1;
		float fSegSum[3]={0};
		int n=0;
		for (int k=SegStart;k<=SegEnd;k++)
		{
			fSegSum[0]=fSegSum[0]+m_vMinPathWorld[k].x;
			fSegSum[1]=fSegSum[1]+m_vMinPathWorld[k].y;
			fSegSum[2]=fSegSum[2]+m_vMinPathWorld[k].z;
			n++;
		}
		miiCNode<double,float> mtempPointWorld;
		miiCNode<double,int> mtempPoint;
		mtempPointWorld.x=fSegSum[0]/n;
		mtempPointWorld.y=fSegSum[1]/n;
		mtempPointWorld.z=fSegSum[2]/n;
		m_vMinPathSegWorld.push_back(mtempPointWorld);
		float nCord[3];
		float fCord[3];
		nCord[0] =mtempPointWorld.x; 
		nCord[1] =mtempPointWorld.y; 
		nCord[2] =mtempPointWorld.z; 
		m_pBaseImgInfo->WorldToImage(nCord);	
		mtempPoint.x = int(nCord[0]+0.5);
		mtempPoint.y = int(nCord[1]+0.5);
		mtempPoint.z =int(nCord[2]+0.5);
		CorrectImagePos(mtempPoint);
		m_vMinPathSeg.push_back(m_vMinPath[SegNUM]);
	}
}
bool miiMinPath::UpdateVSPSet()
{
	if (m_vVSPSet.size() == 0)
	{
		m_vVSPSet.clear();
	}
	m_vVSPSet.push_back(m_iEndPointWorld);
	return true;
}
bool miiMinPath::UpdateCurveSet(const zxhImageInfo *pImageInfo,float nStartPoint[3],char *chResultPath)
{
	if(m_vVSPSet.size()==0)
	{
		cout << "The Vessel segment point set is empty!" << endl;
		return false;
	}
	v_vCurveSet.clear();
	miiCNode<double, int> iEndPoint;
		vector<miiCNode<double,int>> vSgmtMinPath;
        vector<miiCNode<double,float> >vMinPathWorld;
		miiCNode<double,int> iTempPoint, iOrgStartPoint;
		miiCNode<double,float>fTempPointWorld,fOrgStartPointWorld;
	cout << "There are "<<m_vVSPSet.size() <<" Vessel Segment point(s) in..." << endl;
	for(int l=0;l<m_vVSPSet.size();l++)//back tracking from vessel segment point set as curves
	{
			vSgmtMinPath.clear();
			vMinPathWorld.clear();
		// define a array for neighbor searching
		int gNbr[6][3] = { {-1, 0, 0}, \
		{ 1, 0, 0}, \
		{ 0,-1, 0}, \
		{ 0, 1, 0}, \
		{ 0, 0,-1}, \
		{ 0, 0, 1} };

		// define the end point as the source point

		float fEndopointVec[3]={0,0,0};
		fEndopointVec[0]=m_vVSPSet[l].x;
		fEndopointVec[1]=m_vVSPSet[l].y;
		fEndopointVec[2]=m_vVSPSet[l].z;
		pImageInfo->WorldToImage(fEndopointVec);
		iEndPoint.x=int(fEndopointVec[0]+0.5);
		iEndPoint.y=int(fEndopointVec[1]+0.5);
		iEndPoint.z=int(fEndopointVec[2]+0.5);
		iEndPoint.val=m_vVSPSet[l].val;
		int nMinX = iEndPoint.x;
		int nMinY = iEndPoint.y;
		int nMinZ = iEndPoint.z;
		int nCentX = iEndPoint.x;
		int nCentY = iEndPoint.y;
		int nCentZ = iEndPoint.z;
		int nx, ny, nz;
		// save the source point
		vSgmtMinPath.push_back(iEndPoint);

		// set starting point
		float nCord[3]={0,0,0};
		float fCord[3]={0,0,0};
		nCord[0] = nStartPoint[0]; 
		nCord[1] = nStartPoint[1]; 
		nCord[2] = nStartPoint[2]; 
		fOrgStartPointWorld.x=nCord[0];
		fOrgStartPointWorld.y=nCord[1];
		fOrgStartPointWorld.z=nCord[2];

		pImageInfo->WorldToImage(nCord);	
		iOrgStartPoint.x = int(nCord[0]+0.5);
		iOrgStartPoint.y = int(nCord[1]+0.5);
		iOrgStartPoint.z =int(nCord[2]+0.5);
		CorrectImagePos(iOrgStartPoint);
		cout<<"befor while"<<endl;
		while (1)
		{
			for (int i = 0; i < 6; i++)
			{
				nx = nCentX + gNbr[i][0];
				ny = nCentY + gNbr[i][1];
				nz = nCentZ + gNbr[i][2];

				// reach the starting point
				if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y && nz == iOrgStartPoint.z)
				{
					iTempPoint.x = nx;
					iTempPoint.y = ny;
					iTempPoint.z = nz;
					iTempPoint.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					vSgmtMinPath.push_back(iTempPoint);
					for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
					{
						m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(l+1)*m_iMinimalUPValue);//sAMinNodeImg[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]=m_iMinimalUPValue;//Add by JDQ for output purpose
						//m_sMolPontInts.push_back(sImgData[vSgmtMinPath[j].z * m_nImgWY * m_nImgWX + vSgmtMinPath[j].y * m_nImgWX + vSgmtMinPath[j].x]);
						miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
						vMinPathWorld.push_back(fSgmPointWorld);
					}  
					cout<<"in while1"<<endl;
					//CoutSandEndPosi();//output parameter images
					break;
				}

				// prevent out of boundary
				if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ)
				{
					// find the minimal value around the center point 
					double m=m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];
					double n=m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
					if (m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] <= \
						m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX])
					{
						nMinX = nx;
						nMinY = ny;
						nMinZ = nz;
					}
				}
			}//for
				cout<<"after for"<<endl;
			nCentX = nMinX;
			nCentY = nMinY;
			nCentZ = nMinZ;
			fCord[0]=(float)nMinX;
			fCord[1]=(float)nMinY;
			fCord[2]=(float)nMinZ;
			// save the minimal value
			pImageInfo->ImageToWorld(fCord);
			iTempPoint.x = nMinX;
			iTempPoint.y = nMinY;
			iTempPoint.z = nMinZ;
			iTempPoint.val = m_dU[nMinZ * m_nImgWY * m_nImgWX + nMinY * m_nImgWX + nMinX];
			fTempPointWorld.x=fCord[0];
			fTempPointWorld.y=fCord[1];
			fTempPointWorld.z=fCord[2];
			vSgmtMinPath.push_back(iTempPoint);

			float nDistmm =CalcDistmm2(fTempPointWorld,fOrgStartPointWorld);

			nDistmm = sqrt((double)nDistmm);
			cout<<"after for1"<<endl;
			if (nDistmm < m_nEndDistmm) //if (nDist < 3) Add by JDq
			{
				cout<<"in while2"<<endl;
				// save last point (starting point)
				vSgmtMinPath.push_back(iOrgStartPoint);
					cout<<"in while3"<<endl;
				for (int j = vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					cout<<"in while4"<<endl;
					m_zxhMinNodeImg.SetPixelByGreyscale(vSgmtMinPath[j].x,vSgmtMinPath[j].y,vSgmtMinPath[j].z,0,(l+1)*m_iMinimalUPValue);
					cout<<"in while5"<<endl;
					miiCNode<double,float> fSgmPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
					cout<<"in while6"<<endl;
					vMinPathWorld.push_back(fSgmPointWorld);
				} 
				
				//CoutSandEndPosi();//output parameter image
				break;
			}		
		}//while
		cout<<"after while"<<endl;
		v_vCurveSet.push_back(vMinPathWorld);
	}//for l	
	
	return true;
}
bool miiMinPath::Outputvtktxt(const char *chResultPath)
{
	const char *chLResultPath=chResultPath;
	for(int l=0;l<v_vCurveSet.size();l++)
	{
		//write the result curve into files;	
		//char chTemp[25];
		//_itoa_s(l, chTemp, 10);
		//int nLen0 = strlen(chResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".vtk") + 1;
		//char *chFileName0=NULL;
		//chFileName0=new char[nLen0];
		//assert(chFileName0!= NULL);
		//strcpy(chFileName0, chResultPath);
		//strcat(chFileName0, "/MCLine");
		//strcat(chFileName0, chTemp);
		//strcat(chFileName0, ".vtk");
		//cout<<"befor while2"<<endl;
		//WriteCA2Vtk_O(chFileName0,v_vCurveSet[l]);
		//delete[]chFileName0;
		//chFileName0=NULL;

		//char *chFileName1=NULL;
		//chFileName1=new char[nLen0];
		//assert(chFileName1!= NULL);
		//strcpy(chFileName1, chResultPath);
		//strcat(chFileName1, "/MCLine");
		//strcat(chFileName1, chTemp);
		//strcat(chFileName1, ".txt");
		//WriteCA2Txt_O(chFileName1,v_vCurveSet[l]);
		//delete[]chFileName1;
		//chFileName1=NULL;
		/////////////
		char chTemp[25];
		_itoa_s(l, chTemp, 10);
		int nLen0 = strlen(chLResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".vtk") + 1;
		char *chFileName0 = (char *)calloc(1,nLen0);
		strcpy(chFileName0, chLResultPath);
		strcat(chFileName0, "/MCLine");
		strcat(chFileName0, chTemp);
		strcat(chFileName0, ".vtk");
		WriteCA2Vtk_O(chFileName0,v_vCurveSet[l]);
		free(chFileName0);

		int nLen1 = strlen(chLResultPath) + strlen("/MCLine") + strlen(chTemp) + strlen(".txt") + 1;
		char *chFileName1 = (char *)calloc(1,nLen1);
		strcpy(chFileName1, chLResultPath);
		strcat(chFileName1, "/MCLine");
		strcat(chFileName1, chTemp);
		strcat(chFileName1, ".txt");
		WriteCA2Txt_O(chFileName1,v_vCurveSet[l]);
		free(chFileName1);

		cout << "The "<<l+1<<"-th of "<<v_vCurveSet.size()<<" curve is saved in "<<chLResultPath<<"successfully."<< endl;

	}//for
	return true;
}
