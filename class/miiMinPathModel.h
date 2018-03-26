/*#pragma once*/
#ifndef miiMinPathModel_h
#define miiMinPathModel_h

/** include files **/
#include "miiMinPath.h"
#define ZXHJDQCAE_FDSETPLENGTHMM 5 
int GetTime() ;
class miiMinPathModel: public miiMinPath
{
public:
	///
	miiMinPathModel(int nImgWX, int nImgWY, int nImgWZ,float nImgSpacing[4],const zxhImageInfo *pBaseImgInfo, char *chResultPath,\
		 bool bReadMeanFlg = false);
	///miiMinPathModel(int nImgWX, int nImgWY, int nImgWZ, const zxhImageInfo *pBaseImgInfo, int nMaxSgmtDist, bool bReadMeanFlg = false);
	///Change by JDQ

	///
	~miiMinPathModel();	 


	/// the initiation of fast-marching-method for model
	void FastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, \
		float nStartPoint[], float nEndPoint[], int nMethod = 4);
	/// the new initiation of fast-marching-method for the model
	void NewFastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, \
		float nStartPoint[], float nEndPoint[], int nMethod = 4);
	void ModifiedFastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, \
		float nStartPoint[], float nEndPoint[], int nMethod = 4);
	void ModifiedFastMarchingInitLL(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, \
		float nStartPoint[], float nEndPoint[], int nMethod = 4);
	void ModifiedFastMarchingInitLLSM(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, \
		float nStartPoint[], float nEndPoint[], int nMethod = 4);
	/// the re-initiation of fast-marching-method for model
	bool ReFastMarchingInit(bool bCrctFlag = true);
	bool ModifiedReFastMarchingInit(bool bCrctFlag = true);
	bool ModifiedReFastMarchingInitdontMoveModel(bool bCrctFlag = true);
	bool DFMReFastMarchingWOBP(bool bCrctFlag = true);
	
	
	bool ModifiedReFastMarchingInitdontMoveModelSM(bool bCrctFlag = true);
	bool NewReFastMarchingInit(bool bCrctFlag = true);

	
	/// 
	int FastMarchingEvolution(int nMaxItrNum,int i);
	int FastMarchingEvolution(int nMaxItrNum);
	int NewFastMarchingEvolution(int nMaxItrNum,int i);
	int ModifiedFastMarchingEvolution(int nMaxItrNum,int i,const short *sImgData);
	int ModifiedFastMarchingEvolutiondontMoveModel(int nMaxItrNum,int i,const short *sImgData);//model line doesn't move no matter correction.
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointC(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath = 2000);//model line doesn't move no matter correction.with find mininal path;point C;
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNew(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaF(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_D(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask

	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_D_LR(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL zxhcaeDMP_1.1_MLLM
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_3S(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;by 3 conditions

	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;vessel segmentpoint selection conditions
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC_DMP(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;vessel segmentpoint selection conditions;change into DMP

	
	
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC_BPD(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;vessel segmentpoint selection conditions;Bad points detection
		int DFM_WOBP(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;vessel segmentpoint selection conditions;Bad points detection
	
	
	
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel_VSPSC_BPD_O(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;vessel segmentpoint selection conditions;Bad points detection

	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DNB(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;SBP until find the best;Recursion
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DNB_WORec(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;SBP find the best once;without Recursion
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DNBCV(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;SBP;Covergence
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DFMP(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;copye from SBP above,but except it.
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLLDirH(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLLSM(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath);//Semi mehthod;The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL
	
	
	///with mask below
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMask(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//model line doesn't move no matter correction.with find mininal path;point C;with line mask filter
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//model line doesn't move no matter correction.with find mininal path;point C;with line mask filter
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_D(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//model line doesn't move no matter correction.with find mininal path;point C;with line mask filter;with the new intersection mask
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_DM(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//model line doesn't move no matter correction.with find mininal path;point C;with line mask filter;with the new intersection mask
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_DeaF(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//model line doesn't move no matter correction.with find mininal path;point C;with line mask filter;Deacrease the weight of direction and intensity
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_DeaF_D(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//model line doesn't move no matter correction.with find mininal path;point C;with line mask filter;Deacrease the weight of direction and intensity;with the new intersection mask
	int ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLLSM(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//Semi method;model line doesn't move no matter correction.with find mininal path;point C;with line mask filter


	int ModifiedFastMarchingEvolutiondontMoveModelwithMask(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const short *sLineMaskData);//model line doesn't move no matter correction.with find mininal path;point C;with line mask filte
	int ModifiedFastMarchingEvolutiondontMoveModelWithLineMask(int nMaxItrNum,int i,const short *sImgData,const short *sLineMaskData);//model line doesn't move no matter correction.With line Mask filter
///
	int DFM_KT(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;vessel segmentpoint selection conditions;Bad points detection
	///
	int DFM_KT_NewNeig(int nMaxItrNum,int i,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,const zxhImageDataT<short> &imgReadLineMask,vector< miiCNode<double, float> > &vUVectorPoint);//The New one comes from with mask;model line doesn't move no matter correction.with find mininal path;point C;Get developing information from LL,Decrease the weight of direction and intensity;copy from MLL but with mask;proceed the intersection mask;select suitable vessel vector;vessel segmentpoint selection conditions;Bad points detection
	
	///
	bool SetModelPoints(vector< miiCNode<double, float> > vModelPoints, float nModelStartPoint[3]);
	///
	bool SetModelPointsStartPointCorrect(vector< miiCNode<double, float> > vModelPoints, float fModelStartPoint[3],float fREFStartPoint[3]);

	bool SetModelPointsWithIntensity(vector< miiCNode<double, float> > vModelPointsWorldWithIntensity, float nModelStartPointWithIntensity[5]);
	bool SetModelPointsWithIntensityCorrect(vector< miiCNode<double, float> > vModelPointsWorldWithIntensity, float nModelStartPointWithIntensity[5],float fREFStartPoint[3]);

	
	///set the intensity of startpoint of unseean image
	inline bool SetMeanStdDevUnseenImg(double  UnseenMean, double UnseenStdDev){

		m_dUnseenStartMean=UnseenMean;
		m_dUnseenStartStdDev=UnseenStdDev;
		return true;
	};
	inline bool SetMeanStdDevAtlasImg(double  AtlasMean, double AtlasStdDev){

		m_dAtlasStartMean=AtlasMean;
		m_dAtlasStartStdDev=AtlasStdDev;
		return true;
	};
	inline void SBPInit(vector<int>&input){
		m_nSBMap = new int[m_vNarrowBand.size()]; //remember deltete the array

		for (int i = 0; i < m_vNarrowBand.size(); i++)//initiallize the Selecting Best point
		{
			m_nSBMap[i]=0;
			input.push_back(i);
		}
	};
	///
	bool GetModelPoints(vector< miiCNode<double, float> > &vModelPoints);

	///
	bool CorrectModelPoint();
	///
	bool CorrectModelPointBasedSegLength();
	///
	bool CorrectModelPointBasedSegLengthdontMoveModel();
	///
	bool DFM_CorrModel_WOBP();

	///
	bool CorrectModelPointBasedSegLengthdontMoveModelSM();
	///
	bool UpdateMeanStdOfUnseanImg(int StartPosi,int EndPosi);//put the model points into model image to calculate the meanstd of next segmentation add by JDQ 
	///read the mean intensity from the intensity.txt and std of unseen image is the gAortaStdDev[1] ,and update directly add by JDQ
	bool ModifiedUpdateMeanStdOfUnseanImg(int StartPosi,int EndPosi);
	///read the mean intensity from the intensity.txt and std of unseen image is the gAortaStdDev[1],but update between an rangeadd by JDQ
	bool ModifiedUpdateMeanStdOfUnseanImgFB(int StartPosi,int EndPosi);
	///
	bool FindModelPointCbasedOnGeodesicDistance(int *Posi);
	///
	
	
	///
	bool DFM_FindC(int *Posi);
	///
	bool FindModelPointCbasedOnGeodesicDistanceSkipCertainPonts(int *Posi);
	///
	bool FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpoint(int *Posi);//the segmentpoint tangent vetor is smoothed one
	///
	bool FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpointDc(int *Posi);//the segmentpoint tangent vetor is smoothed one;decrease the range of point C
	///
	bool FindModelPointCbasedOnGeodesicDistanceSkipCertainPontsSmoothSegpointSM(int *Posi);//the segmentpoint tangent vetor is smoothed one
	///
	float ModelPointToTanPlaneDistmm(int i,float TanVec[3],float fSmgPoint[3]);
	///
	float ModelPointToTanPlaneDistmm1(int i,float TanVec[3],float fForwardSmgPoint[3]);
	///
	bool SearchStartPoint(const short *sRawImg, int nOldStartPoint[3], int nNewStartPoint[3], int nSearchRange = 100);

	///
	bool UpdateUnseenImgIntensity(const short *sImgData);
	float GetModPointCos(int nModelPos);//Add by JDQ
	float GetModVectorCos(int nModelPos);//Add by JDQ
	inline void GetEndPoint(float nEndPoint[3]) {
		float fCord[3];
		fCord[0] = (float)(m_iSgmtPoint.x);
		fCord[1] = (float)(m_iSgmtPoint.y); 
		fCord[2] = (float)(m_iSgmtPoint.z);
		m_pBaseImgInfo->ImageToWorld(fCord);
		nEndPoint[0] = fCord[0];
		nEndPoint[1] = fCord[1];
		nEndPoint[2] = fCord[2];
	};
	
	inline void GetEndPointWorld(float nEndPointWorld[3]) {
		float fCord[3];
	     nEndPointWorld[0]= (float)(m_iSgmtPointWorld.x);
		 nEndPointWorld[1]= (float)(m_iSgmtPointWorld.y);
		 nEndPointWorld[2]= (float)(m_iSgmtPointWorld.z);
	};
	inline void SetMeanStdDev(double fMean, double fStdDev) {
		m_fMean = fMean;
		m_fStdDev = fStdDev;
		
	};
	inline void SetMaxVeclenth(float maxveclenth) {
		iMaxDisBetweenVec=maxveclenth;
		
	};

	
	
private:
	//
	bool m_bReadMeanFlg;
	
	


	// define the variables for detecting automatically the segmented points   
	// define the max distance between two segmented points
	int m_nMaxDist2;
	//// define the max distance between two segmented points on mm level
	int m_nMinDist2;//Add by JDQ
	float m_nMinDistmm2;
	int m_nSaveFmItr;
	
	


	

	//
	

	//define max/min mean intensity and stdDeve of coronary artery of unseen image.
	double m_dRawMax,m_dRawMin;

	bool m_bFistEvlFlg;
	float iMinDisBetweenVec;
	float iMaxDisBetweenVec;
	
	//
	void UpWind(int x, int y, int z, double theata);//Upwind means potential changes in upwind
	void UpWind1(int x, int y, int z, double theata);//Upwind1 means potential does not change in upwind
	void UpWindPlus(int x, int y, int z, double theata);
	void UpWindPlusTime(int x, int y, int z, double theata);

	void NewUpWindPlus(int x, int y, int z, double theata);
	void NewUpWindTimePlusMixed(int x, int y, int z, double theata);
	void ModifiedUpWindTime(int x, int y, int z, double theata,const short *sImgData);//ModifiedUpWindTime means potential changes in upwind

	void ModifiedUpWindPlusTime(int x, int y, int z, double theata,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind
	void ModifiedUpWindPlusTime2Costheta(int x, int y, int z, float theta1,float theta2,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind,using two costheta
	void ModifiedUpWindPlusTime2CosthetaN(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind,using two costheta
	void ModifiedUpWind_WOBP(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind,using two costheta
	
	void ModUpWind_MultSp(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//
	void ModifiedUpWindPlusTime2CosthetaN_SP(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind,using two costheta
	void ModifiedUpWindPlusTime2CosthetaN_3S(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind,using two costheta

	void ModifiedUpWindPlusTime2CosthetaSM(int x, int y, int z, float theta1,float theta2,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind,using two costheta;weaken the function of direction
	void ModifiedUpWindPlusTime2Costheta_DeaF(int x, int y, int z, float theta1,float theta2,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind,using two costheta

	void ModifiedUpWindPlusTimeLast(int x, int y, int z, double theata,const short *sImgData);//ModifiedUpWindPlusTime means potential changes in upwind
	void ModifiedUpWindPlusTime1(int x, int y, int z, double theata,const short *sImgData);// ModifiedUpWindPlusTime1 means potential does not change in upwind
	void ModifiedUpWindPlusTimeMixed(int x, int y, int z, double theata,const short *sImgData);//ModifiedUpWindPlusTimeMixed means potential doesn't change in upwind ;potential function =1/( v*s+d+e)
	
	void ModUpWind_MultiST(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//multi stencils fast marching
	void ModUpWind_MultiST_test(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//multi stencils fast marching
	void ModUpWind_MultiST_VslsOnly(int x, int y, int z, float theta1,float theta2,float theta3,const short *sImgData);//multi stencils fast marching
	void MultiStencils(int x, int y, int z,float P2);

	void FormTria(miiCNode<float,int> mLN,miiCNode<float,int> mRN,miiCNode<float,int> mUN,miiCNode<float,int> mDN,miiCNode<float,int> mFN,miiCNode<float,int> mBN,miiCNode<float,int>&mMaxN, miiCNode<float,int>&mMedN,miiCNode<float,int> &mMinN);
	void CalVecDelta(miiCNode<float,int>mCN,miiCNode<float,int>mMaxN, miiCNode<float,int>mMedN,miiCNode<float,int> mMinN,float fDel[3],float fCos[3]);
	void Stencil1(int x, int y, int z,float P2);//original
	void Stencil1_1(int x, int y, int z,float P2);//stencil1
	double Stencil1_2(int x, int y, int z,float P2);//stencil1 for 2type
	double Stencil2(int x, int y, int z,float P2);//stencil 2
	double Stencil2_2(int x, int y, int z,float P2);//stencil 2
	double Stencil3(int x, int y, int z,float P2);//stencil 3
	double Stencil3_2(int x, int y, int z,float P2);//stencil 3
	double Stencil4(int x, int y, int z,float P2);//stencil 4
	double Stencil4_2(int x, int y, int z,float P2);//stencil 3
	double Stencil5(int x, int y, int z,float P2);//stencil 5
	double Stencil5_2(int x, int y, int z,float P2);//stencil 3
	double Stencil6(int x, int y, int z,float P2);//stencil 6
	double Stencil6_2(int x, int y, int z,float P2);//stencil 3


	double SolveEq(float ua,float ub,float uc,float ca,float cb,float cc,float ddela,float ddelb,float ddelc,float dP2);
	bool ConvS(miiCNode<double,int> iNewNode,const zxhImageInfo *pImageInfo,float nStartPoint[3]);
	bool NormSpeed();
	bool PotentialFunction(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, int nMethod);
	
	bool ImgGradient2Invs(const short *sImgData, double *dImgG);

	bool Vesselness2Invs(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2);

	bool VesselnessSim2Invs(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2);
	bool VesselnessSim2InvsPlus(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2);
	bool VesselnessSim2InvsPlusTimeMixed(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, double *dP2);
	bool NewPotentialFunction(const short *sImgData, const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, int nMethod);
	bool VesselnessInit(const zxhImageInfo *pImgInfo, \
		const short *sVslsData, const zxhImageInfo *pVslsInfo, double *m_sVes,double *m_sNormVes);
	bool UnseenImgIntensitySet(const short *sImgData,double *m_sUnseenImgIntensity);
	
	bool GetMeanStdDev(const short *sImgData, double &dMean, double &StdDev);

	bool UpdateModelSgmtPoint();
	bool UpdateModelSgmtPointNew();//limit the max length between point C and point D
	bool UpdateModelSgmtPointCD();//limit the max length between point C and point D;add the point D condition
	bool UpdateModelSgmtPointNewSM();//limit the max length between point C and point D
	
	//
	
	float CalcTheata(int cx, int cy, int cz, bool bInitFlg = false);
	float CalcSegTheta(int cx, int cy, int cz);// Calculate the cosine value of the theta between laststartpointtocurrent and currentstartpointtodfNode.
	float CalcTheataN(int cx, int cy, int cz,float currentvec[3] ,bool bInitFlg = false);//using u1
	float CalcSegThetaN(int cx, int cy, int cz,float currentvec[3]);//using u1
	int FindnxyPoiC(miiCNode<double,float>dnxyzWorld);

	float CalcThetaMV(int cx, int cy, int cz,float currentvec[3]);//
	float CalcThetaMV_DMP(int cx, int cy, int cz,float currentvec[3]);//
	float CalcThetaMLV(int cx, int cy, int cz,float currentvec[3]);//using u1 as current vessel vector; last start to current start point as last vessel vector
	float CalcThetaMLV_DMP(int cx, int cy, int cz,float currentvec[3]);//using u1 as current vessel vector; last start to current start point as last vessel vector
	bool SelectVV(int cx, int cy,int cz,float u1[3],float currentvec[3]);//Select right direction of vessel vector at (cx,cy)

	bool SelectVV_3S(int cx, int cy,int cz,float u1[3],float currentvec[3]);//Select right direction of vessel vector at (cx, cy, cz)

};

#endif
