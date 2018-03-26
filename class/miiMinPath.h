/*#pragma once*/
#ifndef miiMinPath_h
#define miiMinPath_h

// the definition of fast-marching-method
#define FMM_ALIVE 0
#define FMM_TRIAL 1
#define FMM_FAR 2

// define a value for infinity
#define FMM_INF 800000000

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>

#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkPolyLine.h"

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"

#include "miiMinHeap.h"
using namespace std;
int GetTimeNow() ;
class miiMinPath
{
public:
	///
	miiMinPath(){};
	miiMinPath(int nImgWX, int nImgWY, int nImgWZ,float nImgSpacing);
	///miiMinPath(int nImgWX, int nImgWY, int nImgWZ);
	///Change by JDQ

	///
	virtual ~miiMinPath();	 

// 	virtual void FastMarchingInit(const short *sImgData, const zxhImageInfo *pImgInfo, \
// 		float nStartPoint[], float nEndPoint[], int nMethod) = 0;

	/// 
	virtual int FastMarchingEvolution(int nMaxItrNum,int i) = 0;
	virtual int FastMarchingEvolution(int nMaxItrNum) = 0;
	///
	bool FindMinPath(float nStartPoint[3], int nMaxPath = 2000);
	///
	bool FindMinPathWithintensiy(const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);
	///
	bool FindMinPathWithintensiyLL(const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);
	///
	bool FindMinPathWithintensiyLLV(const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);
	///
	bool FindMinPathWithintensiyLL_NewNeig(const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);
	//
	bool FindMinPathWithintensiyLL_O(const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);
	//
	bool FindMinPathWithintensiyLLSM(const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);

	///
	float BackTrackInNBWOL(miiCNode<double> iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath = 2000);
	///
	float BackTrackInNB(miiCNode<double> iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath = 2000);
	///
	bool Speed_BackTrack(miiCNode<double> &iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nLPath = 15);
	///
	bool Speed_BackTrackCV(miiCNode<double> &iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nLPath = 15);
	///
	bool BackTrack15(miiCNode<double> &iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],miiCNode<double,int> &iEOUPoint,int nLPath = 15);
	///
	bool FindMinPathWorld(float nStartPoint[3], int nMaxPath);
	///
	bool FindMinPath(int nMaxPath = 2000);
	///
	float CalcDistmm2(miiCNode<double, float> nPixPoint1,miiCNode<double, float> nPixPoint2);//Add by JDQ
	///
	void SmoothPoints(zxhImageData &SmoothLinePointsImg);//calculate the mean coordinate every "m_zxhSmoothBasedArcLengthmm"
	///
	bool GenMinPathGraph(const zxhImageInfo *pImageInfo, const short *sImgData, \
		string strFileName = "abc.nii.gz", int nMethod = 0);

	///
	void WriteCA2Vtk(char *chFileName);

	void WriteCA2VtkMML(char *chFileName);
	void WriteCA2Txt(char *chFileName);

	void WriteCAIntTxt(char *chFileName);
	
	///
	bool SaveU2Image(const zxhImageInfo * pImageInfo, string strFileName);
	///
	int FindImgMidPoint(miiCNode<double, int> iMinNode,bool *blSegPointFind=false,float nVectorDistmm2=0);//Add by JDQ
	int FindImgMidPoint_DistanceThresholdPlusTriangle(miiCNode<double, int> iMinNode,bool *blSegPointFind,float CurMVector[4],float CurVSLVector[4]);
	int FindImgMidPoint_DistanceThresholdPlusVectorAngle(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4]);
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngle(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4]);//2conditions after 2rd segmentation
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMove(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4]);//2conditions after 2rd segmentation
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointB(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath);//bad points are skipped here.
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath);//bad points are skipped here.Modified from LL
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQ(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath);//bad points are skipped here.Modified from LL;Modiified byJDQ
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath);//bad points are skipped here.Modified from LL;Modiified byJDQ
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath);//bad points are skipped here.Modified from LL;Modiified byJDQ
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_DMP(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath);//bad points are skipped here.Modified from LL;Modiified byJDQ


	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath);//bad points are skipped here.Modified from LL;Modiified byJDQ
		int ModifiedFindImgMidPoint_DFM_WOBP(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath);//bad points are skipped here.Modified from LL;Modiified byJDQ
	
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLLSM(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath);//bad points are skipped here.Modified from LL
	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC_BPD_O(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,float fCurMinNodevector[3],int nMaxPath);//bad points are skipped here.Modified from LL;Modiified byJDQ


	int ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBAsSeg(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4],const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData ,int nMaxPath);//bad points are skipped here.
	void CheckIsSgmtPointAndBadPoint(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,bool &blSegPointFind);
	void CheckIsSgmtPointAndBadPointMLL(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,bool &blSegPointFind);
	void CheckIsSgmtPointAndBadPointMLLN(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind);
	void CheckIsSgmtPointAndBadPointMLLN_VSPSC(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind);
	void CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind);
void CheckIsSgmtPoint_DFM_WOP(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind);
	
	void CheckIsSgmtPointAndBadPointMLLN_VSPSC_DMP(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind);
	void CheckIsSgmtPointAndBadPointMLLN_VSPSC_BPD_O(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,float fCurMinNodevector[3],bool &blSegPointFind);


	void CheckIsSgmtPointAndBadPointMLLSM(miiCNode<double, int>iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],const short *sImgData,int nMaxPath,bool &blSegPointFind);
	//int FindImgMidPoint_DistanceThresholdPlusVectorAngle2Factors(miiCNode<double, int> iMinNode,bool *blSegPointFind,float SLNVector[4]);
	bool IsSegmentBPoint(miiCNode<double, float> iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);//judge the segmentpoint a bad point
	bool IsSegmentBPoint_test(miiCNode<double, float> iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath = 2000);
	///
	bool SelectBPasFinalSegPoint(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,int nMaxPath = 2000);
	///without overlap
	bool SelectBPasFinalSegPoint_RandomNoOv(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int> &input,int nMaxPath = 2000);
	///
	bool SelectBPasFinalSegPoint_RandomNoOv15(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int> &input,int nMaxPath = 2000);
	///
	bool SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int> &input,int nMaxPath = 2000);
	///
	bool SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP_WORec(miiCNode<double, int> &iMinNode,const short *sImgData,const zxhImageInfo *pImageInfo,float nStartPoint[3],int &nNumVP,float &fP1,vector<int> &input,int nMaxPath = 2000);
	///
	bool GetRandom(int NUM2,int nNBSumNUM,vector<int>&input,vector<int>&vRandomNum);
	///
	bool SelectSP(miiCNode<double, int> &iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3],int nMaxPath,float fD, float fP,bool bP);

	bool NeiBour(miiCNode<double, int> &iMinNode,const zxhImageInfo *pImageInfo,float fD,float fP,bool bP,float nStartPoint[3],int nMaxPath,vector<miiCNode<double, int>> &vLocNBPoiSet,vector<miiCNode<double, int>> &vNbr);
	bool FindMinPathForBpoints(miiCNode<double, float> fSgmtPointWorld,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath,float &fm_meandis_of_coronaryvessel,float &fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3]);
    ///
    bool FindMinPathForBpointsWithSmoothSeg(miiCNode<double, float> fSgmtPointWorld,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath,float &fm_meandis_of_coronaryvessel,float &fBackTrackDistmmSkipSomePonts,miiCNode<double, float> &fSmoothSgmtPointWorld,miiCNode<double, float> &fForwardSmoothSgmtPointWorld,float fSmoothVesselVecWorld[3]);
	///
	bool FindMinPathForBpoints_test(miiCNode<double, float> iMinNode,const zxhImageInfo *pImageInfo,float nStartPoint[3], const short *sImgData,int nMaxPath,float *fmeandis_of_coronaryvessel,float *fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3]);
	///
	bool FindPointC(miiCNode<double, float> fSgmtPointWorld,float fmeandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3],int *iPosi);
	///
	bool FindPointCNearPosi(miiCNode<double, float>dfMinNodeWorld,int &iPosi);//get the point C for the dfMinNodeWorld;
	///
	bool GetPointCTangentVector(int iPosi,float fPointCTan[3]);
	///
    float CalcDistFromPointCTanPlan(miiCNode<double, float>dfMinNodeWorld,int fPointC,float TanfPointCTan[3]);
	///
	bool CheckBadPoint_FindPointC(miiCNode<double, float> fSgmtPointWorld,float fm_meandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,float fVesselVecWorld[3],int *iPosi );
	///
	bool CheckBadPoint_FindPointCSmoothSeg(miiCNode<double, float> fSgmtPointWorld,float fm_meandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,miiCNode<double, float> fSmoothSgmtPointWorld,miiCNode<double, float> fForwardSmoothSgmtPointWorld,float fSmoothVesselVecWorld[3],int *iPosi );
	///
	bool CheckBadPoint_FindPointCNPosiSmoothSeg(miiCNode<double, float> fSgmtPointWorld,float fm_meandis_of_coronaryvessel,float fBackTrackDistmmSkipSomePonts,miiCNode<double, float> fSmoothSgmtPointWorld,miiCNode<double, float> fForwardSmoothSgmtPointWorld,float fSmoothVesselVecWorld[3],int *iPosi );
	///
	bool UpdateBadPointsNeighbour(miiCNode<double, int> iSgmtPoint);
	//
	bool UpdateBadPointsNeighbour_BPD(miiCNode<double, int> iSgmtPoint);//with the **_BPD
	///
		bool UpdateBadPointsNeighbour_WBPD(miiCNode<double, int> iSgmtPoint);//with the **_BPD
	///
	float ModelPointToTanPlaneDistmm(int i,float TanVec[3],float fSmgPoint[3]);
	///
	void CoutSandEndPosi();
	///
	void CoutSandEndPosiV();
	///
	void CoutSandEndPosi1();
	///
		void CoutSandEndPosi2();
	///
	void CoutSandEndPosiE();//if narrawband's points number is small coutE
	///
	void CoutSandEndPosiL();//if narrawband's points number is small coutE
	///
	
	//map modelpoints to an new image as rawimage,and mark current vetor and then output.
	void MapCurrentVectortoImgandOutput(zxhImageData &zxhModelpointsImg);
	///
	bool BoundaryCorrect(int *PointPos,int ImgNewVslsSize[4]);
	///
	void SmoothPath();
	///
	bool UpdateVSPSet();
	///
	bool UpdateCurveSet(const zxhImageInfo *pImgInfo,float nStartPoint[3],char *chResultPath);
	///
	bool Outputvtktxt(const char *chResultPath);
	///
	void WriteCA2Vtk_O(const char *chFileName,const vector<miiCNode<double,float> > vMinPathWorld);
		///
    void WriteCA2Txt_O(const char *chFileName,const vector<miiCNode<double,float> > vMinPathWorld);
	inline void GetStartPoint(float nStartPoint[3]) {
		float fCord[3];
		fCord[0] = (float)(m_iStartPoint.x); 
		fCord[1] = (float)(m_iStartPoint.y); 
		fCord[2] = (float)(m_iStartPoint.z);
		m_pBaseImgInfo->ImageToWorld(fCord);
		nStartPoint[0] = fCord[0];
		nStartPoint[1] = fCord[1];
		nStartPoint[2] = fCord[2];
		if (nStartPoint[0] >= m_nImgWX)
			nStartPoint[0] = m_nImgWX - 1;
		if (nStartPoint[1] >= m_nImgWY)
			nStartPoint[1] = m_nImgWY - 1;
		if (nStartPoint[2] >= m_nImgWZ)
			nStartPoint[2] = m_nImgWZ - 1;
	};

	inline void SetStartPoint(float nStartPoint[3]) {
		float fCord[3];
		fCord[0] = nStartPoint[0]; fCord[1] = nStartPoint[1]; fCord[2] = nStartPoint[2];
		m_pBaseImgInfo->WorldToImage(fCord);
		m_iStartPoint.x = (int)(fCord[0]+0.5);
		m_iStartPoint.y = (int)(fCord[1]+0.5);
		m_iStartPoint.z = (int)(fCord[2]+0.5);
		if (m_iStartPoint.x >= m_nImgWX)
			m_iStartPoint.x = m_nImgWX - 1;
		if (m_iStartPoint.y >= m_nImgWY)
			m_iStartPoint.y = m_nImgWY - 1;
		if (m_iStartPoint.z >= m_nImgWZ)
			m_iStartPoint.z = m_nImgWZ - 1;
	};

	inline void GetEndPoint(float nEndPoint[3]) {
		float fCord[3];
		fCord[0] = (float)(m_iEndPoint.x);
		fCord[1] = (float)(m_iEndPoint.y); 
		fCord[2] = (float)(m_iEndPoint.z);
		m_pBaseImgInfo->ImageToWorld(fCord);
		nEndPoint[0] = fCord[0];
		nEndPoint[1] = fCord[1];
		nEndPoint[2] = fCord[2];
		if (nEndPoint[0] >= m_nImgWX)
			nEndPoint[0] = m_nImgWX - 1;
		if (nEndPoint[1] >= m_nImgWY)
			nEndPoint[1] = m_nImgWY - 1;
		if (nEndPoint[2] >= m_nImgWZ)
			nEndPoint[2] = m_nImgWZ - 1;
	};
	inline void CorrectImagePos(miiCNode<double, int> &diPoint ){
		if (diPoint.x >= m_nImgWX)
			diPoint.x = m_nImgWX - 1;
	    if (diPoint.y >= m_nImgWY)
			diPoint.y = m_nImgWY - 1;
     	if (diPoint.z >= m_nImgWZ)
			diPoint.z = m_nImgWZ - 1;
		if (diPoint.x <= 0)
			diPoint.x = 0;
		if (diPoint.y <= 0)
			diPoint.y = 0;
		if (diPoint.z <= 0)
			diPoint.z = 0;
	};
	inline miiCNode<double, float> ImageTransToWorldPoint(miiCNode<double, int>&diPoint ){
        float fCord[3];
		CorrectImagePos(diPoint);
		fCord[0] = (float)(diPoint.x);
		fCord[1] = (float)(diPoint.y); 
		fCord[2] = (float)(diPoint.z);
		m_pBaseImgInfo->ImageToWorld(fCord);
		miiCNode<double, float> dfPoint;
		dfPoint.x= fCord[0];
		dfPoint.y = fCord[1];
		dfPoint.z= fCord[2];
		dfPoint.val=diPoint.val;
		return dfPoint;

	};
	inline miiCNode<double, int> WorldTransToImagePoint(miiCNode<double, float>dfPointWorld ){
        float fCord[3];

		fCord[0] = dfPointWorld.x;
		fCord[1] = dfPointWorld.y; 
		fCord[2] = dfPointWorld.z;
		m_pBaseImgInfo->WorldToImage(fCord);

		miiCNode<double, int> diPoint;
		diPoint.x= int(fCord[0]);
		diPoint.y =int(fCord[1]);
		diPoint.z= int(fCord[2]);
		diPoint.val=dfPointWorld.val;
	     CorrectImagePos(diPoint);

		return diPoint;

	};
	inline int CorrectTruePositionOfCoronaryModel(int i)
	{
		if(i<1) return 1;
		//Fi means in front of i
		if(i>m_vModelPointsWorld.size()-1) return m_vModelPointsWorld.size()-1;
		else 
			return i;
	}
//inline void  CoutNPosi(int i){
////ofstream NPosiFile("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu7/Posi.txt",ios::app);
//ofstream NPosiFile("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu0/Posi.txt",ios::app);
//
//if(!NPosiFile)
//   {cerr<<"open error!"<<endl;
//   exit(1);
//   }
//NPosiFile<<"NMPosi:"<<i<<";"<<"";
//NPosiFile.close();
//}
inline void GetVesselVecLocal(vector<miiCNode<double,int>> vMinPath,const zxhImageInfo *pImageInfo){

	int FDist=VESELLSmgTanLENGTH/m_meandis_of_coronaryvessel;
	int vSgmtNum=vMinPath.size() - 1;
	int FPoint=vSgmtNum-FDist;
	if (FPoint<0)FPoint=0;
	float fCord[3];
	float fFCord[3];
	fCord[0]=vMinPath[vSgmtNum].x;
	fCord[1]=vMinPath[vSgmtNum].y;
	fCord[2]=vMinPath[vSgmtNum].z;
	m_pBaseImgInfo->ImageToWorld(fCord);
	fFCord[0]=vMinPath[FPoint].x;
	fFCord[1]=vMinPath[FPoint].y;
	fFCord[2]=vMinPath[FPoint].z;
	m_pBaseImgInfo->ImageToWorld(fFCord);
	//set the point before segmentpoint
	m_fForwardSgmtPointWorld.x=fFCord[0];
	m_fForwardSgmtPointWorld.y=fFCord[1];
	m_fForwardSgmtPointWorld.z=fFCord[2];
	m_fVesselVecWorld[0]=fCord[0]-fFCord[0];
	m_fVesselVecWorld[1]=fCord[1]-fFCord[1];
	m_fVesselVecWorld[2]=fCord[2]-fFCord[2];
}
inline void GetVesselVec(const zxhImageInfo *pImageInfo){

	int FDist=VESELLSmgTanLENGTH/m_meandis_of_coronaryvessel;
	int vSgmtNum=m_vMinPath.size() - 1;
	int FPoint=vSgmtNum-FDist;
	if (FPoint<0)FPoint=0;
	float fCord[3];
	float fFCord[3];
	fCord[0]=m_vMinPath[vSgmtNum].x;
	fCord[1]=m_vMinPath[vSgmtNum].y;
	fCord[2]=m_vMinPath[vSgmtNum].z;
	m_pBaseImgInfo->ImageToWorld(fCord);
	fFCord[0]=m_vMinPath[FPoint].x;
	fFCord[1]=m_vMinPath[FPoint].y;
	fFCord[2]=m_vMinPath[FPoint].z;
	m_pBaseImgInfo->ImageToWorld(fFCord);
	//set the point before segmentpoint
	m_fForwardSgmtPointWorld.x=fFCord[0];
	m_fForwardSgmtPointWorld.y=fFCord[1];
	m_fForwardSgmtPointWorld.z=fFCord[2];
	m_fVesselVecWorld[0]=fCord[0]-fFCord[0];
	m_fVesselVecWorld[1]=fCord[1]-fFCord[1];
	m_fVesselVecWorld[2]=fCord[2]-fFCord[2];
}
inline void GetVesselVecSmooth(const zxhImageInfo *pImageInfo){

	//calculate the mean status around segmentpoint as the status of sgementpoint
	float Segsmoothdist1mm=1;
	float FSegsmoothdist2mm=2;
	int SegsmoothNum=Segsmoothdist1mm/m_meandis_of_coronaryvessel;
	int FSegsmoothNum=FSegsmoothdist2mm/m_meandis_of_coronaryvessel;
	int vSgmtNum=m_vMinPath.size() - 1;
	int forwardi=vSgmtNum-SegsmoothNum;
	if(forwardi<0)forwardi=0;
	float fTempworldSum[3]={0,0,0};
	for(int i=forwardi;i<=vSgmtNum;i++)
	{  
		float fTemp[3];
		fTemp[0]=m_vMinPath[i].x;
		fTemp[1]=m_vMinPath[i].y;
		fTemp[2]=m_vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fTemp);
		fTempworldSum[0]=fTempworldSum[0]+fTemp[0];
		fTempworldSum[1]=fTempworldSum[1]+fTemp[1];
		fTempworldSum[2]=fTempworldSum[2]+fTemp[2];
	}
	m_fSmoothSgmtPointWorld.x=fTempworldSum[0]/(vSgmtNum-forwardi+1);
	m_fSmoothSgmtPointWorld.y=fTempworldSum[1]/(vSgmtNum-forwardi+1);
	m_fSmoothSgmtPointWorld.z=fTempworldSum[2]/(vSgmtNum-forwardi+1);
	//calculate the mean status around front segmentpoint as the status of front sgementpoint
	float fTempworldSum1[3]={0,0,0};
	int forwardFSegsmoothPosi=vSgmtNum-FSegsmoothNum-SegsmoothNum;
	if(forwardFSegsmoothPosi<0)forwardFSegsmoothPosi=0;
	int forwardBSegsmoothPosi=vSgmtNum-FSegsmoothNum+SegsmoothNum;
	if(forwardBSegsmoothPosi>m_vMinPath.size() - 1)forwardBSegsmoothPosi=m_vMinPath.size() - 1;
	for(int i=forwardFSegsmoothPosi;i<=forwardBSegsmoothPosi;i++)
	{
		float fTemp[3];
		fTemp[0]=m_vMinPath[i].x;
		fTemp[1]=m_vMinPath[i].y;
		fTemp[2]=m_vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fTemp);
		fTempworldSum1[0]=fTempworldSum1[0]+fTemp[0];
		fTempworldSum1[1]=fTempworldSum1[1]+fTemp[1];
		fTempworldSum1[2]=fTempworldSum1[2]+fTemp[2];
	}
	m_fForwardSmoothSgmtPointWorld.x=fTempworldSum1[0]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	m_fForwardSmoothSgmtPointWorld.y=fTempworldSum1[1]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	m_fForwardSmoothSgmtPointWorld.z=fTempworldSum1[2]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	//set the tangent vector at segmentpoint
	m_fSmoothVesselVecWorld[0]=m_fSmoothSgmtPointWorld.x-m_fForwardSmoothSgmtPointWorld.x;
	m_fSmoothVesselVecWorld[1]=m_fSmoothSgmtPointWorld.y-m_fForwardSmoothSgmtPointWorld.y;
	m_fSmoothVesselVecWorld[2]=m_fSmoothSgmtPointWorld.z-m_fForwardSmoothSgmtPointWorld.z;
	//calculate the angel of current modelvector and smoothvessel vector
	float vModel[3]={0,0,0};
	vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	float costheta=zxh::VectorOP_Cosine( m_fSmoothVesselVecWorld,vModel, 3 );
}
inline void GetVesselVecSmoothLocal(vector<miiCNode<double,int>> vMinPath,const zxhImageInfo *pImageInfo){

	//calculate the mean status around segmentpoint as the status of sgementpoint
	float Segsmoothdist1mm=1;
	float FSegsmoothdist2mm=2;
	int SegsmoothNum=Segsmoothdist1mm/m_meandis_of_coronaryvessel;
	int FSegsmoothNum=FSegsmoothdist2mm/m_meandis_of_coronaryvessel;
	int vSgmtNum=vMinPath.size() - 1;
	int forwardi=vSgmtNum-SegsmoothNum;
	if(forwardi<0)forwardi=0;
	float fTempworldSum[3]={0,0,0};
	for(int i=forwardi;i<=vSgmtNum;i++)
	{  
		float fTemp[3];
		fTemp[0]=vMinPath[i].x;
		fTemp[1]=vMinPath[i].y;
		fTemp[2]=vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fTemp);
		fTempworldSum[0]=fTempworldSum[0]+fTemp[0];
		fTempworldSum[1]=fTempworldSum[1]+fTemp[1];
		fTempworldSum[2]=fTempworldSum[2]+fTemp[2];
	}
	m_fSmoothSgmtPointWorld.x=fTempworldSum[0]/(vSgmtNum-forwardi+1);
	m_fSmoothSgmtPointWorld.y=fTempworldSum[1]/(vSgmtNum-forwardi+1);
	m_fSmoothSgmtPointWorld.z=fTempworldSum[2]/(vSgmtNum-forwardi+1);
	//calculate the mean status around front segmentpoint as the status of front sgementpoint
	float fTempworldSum1[3]={0,0,0};
	int forwardFSegsmoothPosi=vSgmtNum-FSegsmoothNum-SegsmoothNum;
	if(forwardFSegsmoothPosi<0)forwardFSegsmoothPosi=0;
	int forwardBSegsmoothPosi=vSgmtNum-FSegsmoothNum+SegsmoothNum;
	if(forwardBSegsmoothPosi>vMinPath.size() - 1)forwardBSegsmoothPosi=vMinPath.size() - 1;
	for(int i=forwardFSegsmoothPosi;i<=forwardBSegsmoothPosi;i++)
	{
		float fTemp[3];
		fTemp[0]=vMinPath[i].x;
		fTemp[1]=vMinPath[i].y;
		fTemp[2]=vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fTemp);
		fTempworldSum1[0]=fTempworldSum1[0]+fTemp[0];
		fTempworldSum1[1]=fTempworldSum1[1]+fTemp[1];
		fTempworldSum1[2]=fTempworldSum1[2]+fTemp[2];
	}
	m_fForwardSmoothSgmtPointWorld.x=fTempworldSum1[0]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	m_fForwardSmoothSgmtPointWorld.y=fTempworldSum1[1]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	m_fForwardSmoothSgmtPointWorld.z=fTempworldSum1[2]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	//set the tangent vector at segmentpoint
	m_fSmoothVesselVecWorld[0]=m_fSmoothSgmtPointWorld.x-m_fForwardSmoothSgmtPointWorld.x;
	m_fSmoothVesselVecWorld[1]=m_fSmoothSgmtPointWorld.y-m_fForwardSmoothSgmtPointWorld.y;
	m_fSmoothVesselVecWorld[2]=m_fSmoothSgmtPointWorld.z-m_fForwardSmoothSgmtPointWorld.z;
	//calculate the angel of current modelvector and smoothvessel vector
	float vModel[3]={0,0,0};
	vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	float costheta=zxh::VectorOP_Cosine( m_fSmoothVesselVecWorld,vModel, 3 );
}
inline void GetVesselVecfromVessel(const zxhImageInfo *pImageInfo,vector<miiCNode<>> vMinPath,miiCNode<double,float> *fForwardSgmtPointWorld,float fVesselVecWorld[3]){

	int FDist=VESELLSmgTanLENGTH/m_meandis_of_coronaryvessel;
	int vSgmtNum=vMinPath.size() - 1;
	int FSgmtNum=0;
	float fCord[3];
	float fFCord[3];
	fCord[0]=vMinPath[vSgmtNum].x;
	fCord[1]=vMinPath[vSgmtNum].y;
	fCord[2]=vMinPath[vSgmtNum].z;
	FSgmtNum=vSgmtNum-FDist;
	if(FSgmtNum<0)FSgmtNum=0;
	m_pBaseImgInfo->ImageToWorld(fCord);
	fFCord[0]=vMinPath[FSgmtNum].x;
	fFCord[1]=vMinPath[FSgmtNum].y;
	fFCord[2]=vMinPath[FSgmtNum].z;
	m_pBaseImgInfo->ImageToWorld(fFCord);
	//set the point before segmentpoint
	fForwardSgmtPointWorld->x=fFCord[0];
	fForwardSgmtPointWorld->y=fFCord[1];
	fForwardSgmtPointWorld->z=fFCord[2];
	fVesselVecWorld[0]=fCord[0]-fFCord[0];
	fVesselVecWorld[1]=fCord[1]-fFCord[1];
	fVesselVecWorld[2]=fCord[2]-fFCord[2];
}
inline void GetSmoothVesselVecfromVessel(const zxhImageInfo *pImageInfo,vector<miiCNode<>> vMinPath,miiCNode<double, float> &fSmoothSgmtPointWorld,miiCNode<double,float> &fForwardSmoothSgmtPointWorld,float fSmoothVesselVecWorld[3]){

	//calculate the mean status around segmentpoint as the status of sgementpoint
	float Segsmoothdist2mm=2;
	float FSegsmoothdist3mm=3;
	int SegsmoothNum=Segsmoothdist2mm/m_meandis_of_coronaryvessel;
	int FSegsmoothNum=FSegsmoothdist3mm/m_meandis_of_coronaryvessel;
	int vSgmtNum=vMinPath.size() - 1;
	int forwardi=vSgmtNum-SegsmoothNum;
	if(forwardi<0)forwardi=0;
	float fTempworldSum[3]={0,0,0};
	for(int i=forwardi;i<=vSgmtNum;i++)
	{  
		float fTemp[3];
		fTemp[0]=vMinPath[i].x;
		fTemp[1]=vMinPath[i].y;
		fTemp[2]=vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fTemp);
		fTempworldSum[0]=fTempworldSum[0]+fTemp[0];
		fTempworldSum[1]=fTempworldSum[1]+fTemp[1];
		fTempworldSum[2]=fTempworldSum[2]+fTemp[2];
	}
	fSmoothSgmtPointWorld.x=fTempworldSum[0]/(vSgmtNum-forwardi+1);
	fSmoothSgmtPointWorld.y=fTempworldSum[1]/(vSgmtNum-forwardi+1);
	fSmoothSgmtPointWorld.z=fTempworldSum[2]/(vSgmtNum-forwardi+1);
	//calculate the mean status around front segmentpoint as the status of front sgementpoint
	float fTempworldSum1[3]={0,0,0};
	int forwardFSegsmoothPosi=vSgmtNum-FSegsmoothNum-SegsmoothNum;
	if(forwardFSegsmoothPosi<0)forwardFSegsmoothPosi=0;
	int forwardBSegsmoothPosi=vSgmtNum-FSegsmoothNum+SegsmoothNum;
	if(forwardBSegsmoothPosi>vMinPath.size() - 1)forwardBSegsmoothPosi=vMinPath.size() - 1;
	for(int i=forwardFSegsmoothPosi;i<=forwardBSegsmoothPosi;i++)
	{
		float fTemp[3];
		fTemp[0]=vMinPath[i].x;
		fTemp[1]=vMinPath[i].y;
		fTemp[2]=vMinPath[i].z;
		m_pBaseImgInfo->ImageToWorld(fTemp);
		fTempworldSum1[0]=fTempworldSum1[0]+fTemp[0];
		fTempworldSum1[1]=fTempworldSum1[1]+fTemp[1];
		fTempworldSum1[2]=fTempworldSum1[2]+fTemp[2];
	}
	fForwardSmoothSgmtPointWorld.x=fTempworldSum1[0]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	fForwardSmoothSgmtPointWorld.y=fTempworldSum1[1]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	fForwardSmoothSgmtPointWorld.z=fTempworldSum1[2]/(forwardBSegsmoothPosi-forwardFSegsmoothPosi+1);
	//set the tangent vector at segmentpoint
	fSmoothVesselVecWorld[0]=fSmoothSgmtPointWorld.x-fForwardSmoothSgmtPointWorld.x;
	fSmoothVesselVecWorld[1]=fSmoothSgmtPointWorld.y-fForwardSmoothSgmtPointWorld.y;
	fSmoothVesselVecWorld[2]=fSmoothSgmtPointWorld.z-fForwardSmoothSgmtPointWorld.z;
	//calculate the angel of current modelvector and smoothvessel vector
	float vModel[3]={0,0,0};
	vModel[0] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].x -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].x;
	vModel[1] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].y -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].y;
	vModel[2] = m_vModelPointsWorld[m_nOneSegModelVectorEndPos].z -m_vModelPointsWorld[m_nOneSegModelVectorStartPos].z;
	float costheta=zxh::VectorOP_Cosine( m_fSmoothVesselVecWorld,vModel, 3 );
}
inline float fCalcMinPathLength(vector<miiCNode<double,int>> vSgmtMinPath){
	miiCNode<double,float>fSgmFPointWorld,fSgmBPointWorld;
	float BackTrackDistmm=0;

	for (int j=1;j<vSgmtMinPath.size();j++)
	{
	fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
	fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
	BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
	}

return BackTrackDistmm;
}
inline float fCalcMinPathLengthSkipSomePonts(vector<miiCNode<double,int>> vSgmtMinPath){//skip some points to calculate the total length of the minimal path.
	miiCNode<double,float>fSgmFPointWorld,fSgmBPointWorld;
	float BackTrackDistmm=0.0;
	float fdensityofminmalpath=0.0;
	float CERTAINLENGTH=2.0;
	int BackPintPosi=0;
	for (int j=1;j<vSgmtMinPath.size();j++)
	{
	fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
	fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j-1]);
	BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
	}
	fdensityofminmalpath=BackTrackDistmm/(vSgmtMinPath.size()-1+0.001);
	BackTrackDistmm=0;
	BackPintPosi=(int)(CERTAINLENGTH/fdensityofminmalpath+0.5);
	int i=0;
	for (int j=0;(j+BackPintPosi)<vSgmtMinPath.size();)
	{
	
	fSgmBPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j+BackPintPosi]);
	fSgmFPointWorld=ImageTransToWorldPoint(vSgmtMinPath[j]);
	BackTrackDistmm=BackTrackDistmm+sqrt(CalcDistmm2(fSgmBPointWorld,fSgmFPointWorld));
	j=j+BackPintPosi;
	}
	

	
return BackTrackDistmm;
}
void GetNextSegMeanInteStd(vector<miiCNode<double,int>> vSgmtMinPath,const short *sImgData);

//void CoutVPosi(int i)
//{
//ofstream VPosiFile;
////VPosiFile.open("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu7/Posi.txt",ios::app);
//VPosiFile.open("F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu0/Posi.txt",ios::app);
//if(!VPosiFile)
//   {cerr<<"open error!"<<endl;
//   exit(1);
//   }
//VPosiFile<<"VPosi:"<<i<<";"<<"\n";
//VPosiFile.close();
//}

	void CoutsAMinNodeImg(short *sAMinNodeImg);
	bool SaveImage(const zxhImageInfo * pImageInfo, short *sImgData, string strFileName);
	bool SaveImageD(const zxhImageInfo * pImageInfo, double *sImgData, string strFileName);
public:
	zxhImageData m_zxhMinNodeImg; //short *sAMinNodeImg;	//add by JDQ
	zxhImageData m_zxhBadPointImg;//short *sABadPointImg;	// zxhImageData m_imgABadPointImg
	zxhImageData m_zxhBadPointMaskImg;//short *sBadPointMask;	//to set badpoints'neighbour points which is far as mask;
	//zxhImageData m_zxhModelpointsImg;//to map the modelpoints to raw image and mark current vector
	zxhImageDataF m_zxhNBUvalueImg;//short *sBadPointMask;	//to set u value as an image;
	zxhImageDataF m_zxhNBPvalueImg;//short *sBadPointMask;	//to set badpoints'neighbour points which is far as mask;
	zxhImageDataF m_zxhNBSvalueImg;//dimS value'sreciprocal
	zxhImageDataF m_zxhNBDvalueImg;//dimD value.reciprocal
	zxhImageDataF m_zxhNBD_SegvalueImg;//dimD in selecting segmentpoint value.reciprocal
	zxhImageDataF m_zxhNBSpeedvalueImg;//
	char *m_cResultPath;
protected:
	// the image size
	int m_nImgWX, m_nImgWY, m_nImgWZ;
	// the model image size
	int m_nMolImgWX, m_nMolImgWY, m_nMolImgWZ;//Add by JDQ
	float m_nImgSpacing[4];//Add by JDQ
	// image info
	const zxhImageInfo *m_pBaseImgInfo;

	// define the starting and end point of the coronary 
	miiCNode<> m_iStartPoint, m_iEndPoint;
	miiCNode<> m_iLastStartPoint;
	miiCNode<> m_iModelEndPoint;
	// define an instance for saving the position of the segmented point
	miiCNode<double> m_iSgmtPoint;
	// define the distance function U in Upwind
	miiCNode<double,float> m_iOrgStartPointWorld,m_iStartPointWorld, m_iEndPointWorld,m_fSMEndPointWorld;
	miiCNode<double,float> m_iLastStartPointWorld;
	miiCNode<double,float> m_iModelEndPointWorld;
	miiCNode<double,float> m_iSgmtPointWorld;
	miiCNode<double,float> m_fForwardSgmtPointWorld;//the  points in front of sgmtpointworld

	miiCNode<double,float> m_fSmoothSgmtPointWorld;//the smooth point sgmtpointworld
	miiCNode<double,float> m_fForwardSmoothSgmtPointWorld;//the smooth point in front of sgmtpointworld
	double *m_dU;
	//define the minimallengthof bad point for checking
	float mZXH_MinLengthForBadPointCheck; 
	//define the percentage for badpoint length
	float mZXH_PencentageForBadPointLength; 
	//define the minimal percentage of length while selecting C Points;
	float mZXH_MinPencentageForCPointLength;
	//define the maximal percentage of length while selecting C Points;
	float mZXH_MaxPencentageForCPointLength; 
	// define the map of the fast-marching-method for 
	// saving 'alive', 'trial(Narrow Band)', and 'far'
	int *m_nFmMap;
	//
	int *m_nSBMap;
	// potential of minimal path
	double *m_dP2;
	//define the array for MinNode NB and MinPath
	
	//define the intensity of unseen image(s)
	double *m_sUnseenImgIntensity;
	double *m_dSims;
	//define the vesselness
	double *m_sVes;
	double *m_sNormVes;
	//define the direction
	double *m_dSimd;
	//define the number of visited points
	int m_nLn;
	//define the maximum speed value of visited points
	float m_nMAXLn;
	//define the minimum speed value of visited points
	float m_fMinSpeed;
	// define the min-heap for NarrowBand
	miiMinHeap<> *m_iMinHeap;

	// define the vector of Narrow Band for 'trail' in fast-marching-method
	vector<miiCNode<>> m_vNarrowBand;
	//define the vector of speed about convergence
	vector<miiCNode<>> m_vSpValue;
	//define the vector of normalized speed about convergence
	vector<miiCNode<>> m_vSpNorValue;
	//
	// define a vector for saving the points of minimal path
	vector<miiCNode<>> m_vMinPath;//minimal path points vector
	// define a vector for saving the points of minimal path
	vector<miiCNode<double,float> > m_vMinPathWorld;//minimal path points world vector
	// define a vector for saving the points of minimal path
	vector<miiCNode<>> m_vMinPathSeg;//minimal path points segpoit
	// define a vector for saving the seg points of minimal path
	vector<miiCNode<double,float> > m_vMinPathSegWorld;
	//
	vector<miiCNode<double,float>> m_vVSPSet;//vessel segment point set
	//minimal path points intensity vector
	vector<short> m_sMolPontInts;//backtrack path intensity
	// model points
	vector< miiCNode<double, int> > m_vModelPoints;
	vector< miiCNode<double, float> > m_vModelPointsWorld;//Add by JDQ
	vector< miiCNode<double, float> >m_vModelPointsWorldWithIntensity;//Add by JDQ

	
	// the position of the model points
	int m_nModelPos;
	//the positon of the start model points in one segmentation
	int m_nOneSegStartModelPos;
	//the positon of the last start model points in one segmentation
	int m_nOneSegModelVectorLastStartPos;
	
	///the positon of the start model points in one segmentation
	int m_nOneSegModelVectorStartPos;
	/// arch length of the model from the start point to Point C
	float m_fCurrModelCArchLength ;
	/// arch length of the model from the start point to Point D
	float m_fCurrModelCArchLengthD ;
	///record the maximum length of Selected point in NB
	float m_fMaxSPlength;
	///
	miiCNode<double, int> m_miiMaxSP;
	// the positon of the end model points in one segmentation
	int m_nOneSegModelVectorEndPos;
	// define the variable for recording the number of FMM iteration
	int m_nFmItr;
	int m_nSumFmItr;
	float m_nEndDistmm;//Add by JDQ
	int m_nFMMEvlNum;
	int m_nSum;
	int m_nFMVedPoiNUM;
	int m_nFMVedPoiNUMLast;
	int m_nBadpointSum;
	int m_nRealOriModelPointsPosi;
	miiCNode<double, float> m_nRealOriModelPointsWorld;//this point is always change after one model correction
	miiCNode<double, float> m_nRealOriSegModelPointsWorld;//this model point doesnt change after model correction
	miiCNode<double, float> m_nRealOriSegPointWorld;//the first segment point;
	float m_nKPLambda;//Add by JDQ
	float m_nKPEpsilon;//Add by JDQ
	float m_nKPCosTheta;//Add by JDQ
	float ZXHJDQCAE_FDSETPLENGTH;//Add by JDQ
	float MODELVECTanLENGTH;
	float VESELLSmgTanLENGTH;
	float DistRangeAroundCPoint;
	float MINDISTTOPOINT2;
	float m_ftoallength_model;
	float m_meandis_of_coronarymodel;//Add by JDQ
	float m_meandis_of_coronaryvessel;//Add by JDQ
	float m_fModelVecWorld[3];//model direction vector
	float m_fVesselVecWorld[3];//vessel direction vector
	float m_fSmoothVesselVecWorld[3];//vessel direction vector
	float m_fBackTrackDistmm;
	float m_fBackTrackDistmmSkipSomePonts;
	int m_iFarPValue;
	int m_iAlivePValue;
	int m_iActivePValue;
	int m_iBadPointValue;
	int m_iMinimalUPValue;
	bool IsEndPoint;//whether or not Point C is the Endpoint
	bool IsCurEndPoint;//whether or not Current point is the Endpoint
	bool bIsCurSMEndPoint;//wheter or not current point is near the semi-endpoint
	/// the initiation of fast-marching-method(FMM)
	void FastMarchingInitBase(const short *sImgData, float nStartPoint[], float nEndPoint[]);
	void FastMarchingInitBaseSM(const short *sImgData, float nStartPoint[], float nEndPoint[]);
    // define mean intensity and stdDeve of coronary artery of unseen image.
	double m_fMean, m_fStdDev;
    //define mean intensity and stdDev of start point for the unseen image
	double m_dUnseenStartMean, m_dUnseenStartStdDev;
	//
	int m_nFMMTimeS;
	//
	double m_dAtlasStartMean,m_dAtlasStartStdDev;
	
	//
	void UpWind(int x, int y, int z);
	//
	
	bool UpdateNarrowBandVal(miiCNode<>);

	bool QuadraticRoots(double a, double b, double c, double dRoots[2]);	
	///calculate the mm level distance from the pixel level
	vector<vector<miiCNode<double,float> >> v_vCurveSet;
	
};

#endif
