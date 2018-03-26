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
#include "miiMinPathOrg.h"

//
miiMinPathOrg::miiMinPathOrg(int nImgWX, int nImgWY, int nImgWZ, const zxhImageInfo *pBaseImgInfo)
{
	m_nImgWX = nImgWX;
	m_nImgWY = nImgWY;
	m_nImgWZ = nImgWZ;
	m_pBaseImgInfo = pBaseImgInfo;
	m_dU = new double[nImgWX * nImgWY * nImgWZ];
	m_nFmMap = new int[nImgWX * nImgWY * nImgWZ];
	m_dP2 = new double[nImgWX * nImgWY * nImgWZ];

	// create a min-heap for Narrow Band
	m_iMinHeap = new miiMinHeap<double>();
}

//
miiMinPathOrg::~miiMinPathOrg()
{
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
void miiMinPathOrg::FastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
			float nStartPoint[], float nEndPoint[], int nMethod)
{
	// FMM initialization
	FastMarchingInitBase(sImgData, nStartPoint, nEndPoint);

	// potential function initialization
	PotentialFunction(sImgData, nMethod);
}

// Function Name: FastMarchingEvolution()
//
// Parameters: 
//
// Description: fast-marching-method evolution
//
// Returns: iteration number
//
int miiMinPathOrg::FastMarchingEvolution(int nMaxItrNum)
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

	// iterate for FMM
	while (true)
	{
		m_nFmItr++;

		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);

		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;

		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];

			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				UpWind(nx, ny, nz);

				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
				} 
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);
				}
			}
		}

		if (iMinNode.x == m_iEndPoint.x && iMinNode.y == m_iEndPoint.y && iMinNode.z == m_iEndPoint.z)
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
	}

	return m_nFmItr;
}
// Function Name: FastMarchingEvolution()
//
// Parameters: 
//
// Description: fast-marching-method evolution
//
// Returns: iteration number
//
int miiMinPathOrg::FastMarchingEvolution(int nMaxItrNum,int i)
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

	// iterate for FMM
	while (true)
	{
		m_nFmItr++;

		// extract and remove the minimum from Narrow Band 
		iMinNode = m_iMinHeap->HeapExtractMin(m_vNarrowBand);

		// update FMM map by setting the minimum's location in FMM map as 'alive' 
		m_nFmMap[iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x] = FMM_ALIVE;

		// search 6-neighbors of the minimum
		for (int i = 0; i < 6; i++)
		{
			int nx = iMinNode.x + gNbr[i][0];
			int ny = iMinNode.y + gNbr[i][1];
			int nz = iMinNode.z + gNbr[i][2];

			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY && nz >= 0 && nz < m_nImgWZ \
				&& m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] != FMM_ALIVE)
			{
				UpWind(nx, ny, nz);

				iNewNode.x = nx; iNewNode.y = ny; iNewNode.z = nz; 
				iNewNode.val = m_dU[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx];

				if (m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] == FMM_FAR) // 'far'
				{
					m_nFmMap[nz * m_nImgWY * m_nImgWX + ny * m_nImgWX + nx] = FMM_TRIAL;
					m_iMinHeap->MinHeapInsert(m_vNarrowBand, iNewNode);
				} 
				else //'trial'
				{
					UpdateNarrowBandVal(iNewNode);
				}
			}
		}

		if (iMinNode.x == m_iEndPoint.x && iMinNode.y == m_iEndPoint.y && iMinNode.z == m_iEndPoint.z)
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
	}

	return m_nFmItr;
}
bool miiMinPathOrg::PotentialFunction(const short *sImgData, int nMethod)
{
	if (sImgData == NULL)
	{
		std::cout << "The pointer cannot be null!" <<endl;
		return false;
	}

	switch (nMethod)
	{
	case 1:
		ImgGradient2Invs(sImgData, m_dP2);
		break;
	default:
		break;
	}

	return true;	
}

// Function Name: ImgGradient2Invs()
//
// Parameters: *sImgData: the pointer for the data of 3D image
//			   *dImgG: the pointer for the inverse of the square of Gradient
//
// Description: calculate the inverse of the square of Gradient 
//
// Returns: 
//
bool miiMinPathOrg::ImgGradient2Invs(const short *sImgData, double *dImgG)
{
	if (sImgData == NULL || dImgG == NULL)
	{
		std::cout << "The pointer is null!" << endl;
		return false;
	}

	double Fx, Fy, Fz;

	int nx, ny, nz;

	for (int iz = 0; iz < m_nImgWZ; iz++)
	{
		for (int iy = 0; iy < m_nImgWY; iy++)
		{
			for (int ix = 0; ix < m_nImgWX; ix++)
			{
				if (iz == m_nImgWZ - 1)	nz = iz - 1;
				else nz = iz + 1;

				if (iy == m_nImgWY - 1)	ny = iy - 1;
				else ny = iy + 1;

				if (ix == m_nImgWX - 1)	nx = ix - 1;
				else nx = ix + 1;

				if (ix == 275 && iy == 258 && iz == 161)
				{
					int temp = 1;
				}

				Fx = (double)(sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + nx] \
					- sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]);
				Fy = (double)(sImgData[iz * m_nImgWY * m_nImgWX + ny * m_nImgWX + ix] \
					- sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]);
				Fz = (double)(sImgData[nz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] \
					- sImgData[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix]);

				dImgG[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] = \
					(1 / (Fx * Fx + Fy * Fy + Fz * Fz + 0.0000000000001));

				if (dImgG[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix] > 0)
				{
					double temp = dImgG[iz * m_nImgWY * m_nImgWX + iy * m_nImgWX + ix];
				}
			}
		}
	}

	return true;
}
