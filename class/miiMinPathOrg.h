/*#pragma once*/
#ifndef miiMinPathOrg_h
#define miiMinPathOrg_h

/** include files **/
#include "miiMinPath.h"

class miiMinPathOrg: public miiMinPath
{
public:
	///
	miiMinPathOrg(int nImgWX, int nImgWY, int nImgWZ, const zxhImageInfo *pBaseImgInfo);

	///
	~miiMinPathOrg();

	/// the initiation of fast-marching-method(FMM)
	void FastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
		float nStartPoint[], float nEndPoint[], int nMethod);


	/// 
	
	int FastMarchingEvolution(int nMaxItrNum);
	int FastMarchingEvolution(int nMaxItrNum,int i);//Add by JDQ
private:
	bool PotentialFunction(const short *sImgData, int nMethod = 1);

	bool ImgGradient2Invs(const short *sImgData, double *dImgG);

};

#endif
