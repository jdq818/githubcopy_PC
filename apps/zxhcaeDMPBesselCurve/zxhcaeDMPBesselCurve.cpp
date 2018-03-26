
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"


#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"


#include "vtkUnstructuredGrid.h"

// for read
#include "vtkUnstructuredGridReader.h"

// for write
#include "vtkUnstructuredGridWriter.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"
#include "miiMinHeap.h"

#define TAB_CHAR	9
using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
	void cvPoint(float px,float py,float pz);
}PointCordTypeDef;
void PointCordTypeDef::cvPoint(float px,float py,float pz)
{
	x=px;
	y=py;
	z=pz;
};
typedef struct
{
	float fMean;
	float fStdVar;
	float fMaxErr;
}ErrTypeDef;

void ReadVtk(char *chFileName, vector<PointCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find VTK-file!" << endl;
		return;
	}
	if (!PointCord.empty())
	{
		PointCord.clear();
	}

	vtkSmartPointer<vtkUnstructuredGridReader> iVtkReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	iVtkReader->SetFileName( chFileName );
	iVtkReader->Update();

	vtkSmartPointer<vtkUnstructuredGrid> iGridRead = iVtkReader->GetOutput();

	int nPointNum = iGridRead->GetMaxCellSize();

	double dCord[3];
	PointCordTypeDef strctTempPoint;

	for (int i = 0; i < nPointNum; i++)
	{
		iGridRead->GetPoint(i, dCord);
		strctTempPoint.x = dCord[0];
		strctTempPoint.y = dCord[1];
		strctTempPoint.z = dCord[2];
		PointCord.push_back(strctTempPoint);
	}
}
bool ReadTxt(char *chFileName,vector<short>&sMolPontInts)
{
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
     sMolPontInts.push_back(data);
         }
	return true;

   }
}

void WriteVtk(vector< PointCordTypeDef > PointCord, char* chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	/*	int nPointNum = PointCord.size();*/

	int nPointNum = PointCord.size();

	for (int i = 0; i < nPointNum; i++)
	{
		iPoints->InsertNextPoint(PointCord[i].x, PointCord[i].y, PointCord[i].z);
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
void WriteTxt(vector< PointCordTypeDef > PointCord, char* chFileName)
{
    ofstream WriteFileTxt(chFileName);
	int nPointNum = PointCord.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
	 fImgPixel[0] = PointCord[i].x;
	 fImgPixel[1] = PointCord[i].y;
	 fImgPixel[2] = PointCord[i].z;
     WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4)<< -fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] << "\n ";

	}	
}
void WriteCAIntTxt(vector<short> sMolPontInts,char* chFileName)
{
ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = sMolPontInts.size();
	for (int i = 0; i < nPointNum; i++)
	{
     WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4) <<sMolPontInts[i] <<"\n";

	}	
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
float bezier2funcX(float uu,PointCordTypeDef *controlP){  
   float part0 = controlP[0].x * uu * uu;  
   float part1 = 2 * controlP[1].x * uu *(1 - uu);  
   float part2 = controlP[2].x * (1 - uu) * (1 - uu);     
   return part0 + part1 + part2 ;   
} 
float bezier2funcY(float uu,PointCordTypeDef *controlP){  
   float part0 = controlP[0].y * uu * uu;  
   float part1 = 2 * controlP[1].y * uu *(1 - uu);  
   float part2 = controlP[2].y * (1 - uu) * (1 - uu);     
   return part0 + part1 + part2 ;   
}  
float bezier2funcZ(float uu,PointCordTypeDef *controlP){  
   float part0 = controlP[0].z * uu * uu;  
   float part1 = 2 * controlP[1].z * uu *(1 - uu);  
   float part2 = controlP[2].z * (1 - uu) * (1 - uu);     
   return part0 + part1 + part2 ;   
}  
float bezier3funcX(float uu,PointCordTypeDef *controlP){  
   float part0 = controlP[0].x * uu * uu * uu;  
   float part1 = 3 * controlP[1].x * uu * uu * (1 - uu);  
   float part2 = 3 * controlP[2].x * uu * (1 - uu) * (1 - uu);  
   float part3 = controlP[3].x * (1 - uu) * (1 - uu) * (1 - uu);     
   return part0 + part1 + part2 + part3;   
}      
float bezier3funcY(float uu,PointCordTypeDef *controlP){  
   float part0 = controlP[0].y * uu * uu * uu;  
   float part1 = 3 * controlP[1].y * uu * uu * (1 - uu);  
   float part2 = 3 * controlP[2].y * uu * (1 - uu) * (1 - uu);  
   float part3 = controlP[3].y * (1 - uu) * (1 - uu) * (1 - uu);     
   return part0 + part1 + part2 + part3;   
}
float bezier3funcZ(float uu,PointCordTypeDef *controlP){  
   float part0 = controlP[0].z * uu * uu * uu;  
   float part1 = 3 * controlP[1].z * uu * uu * (1 - uu);  
   float part2 = 3 * controlP[2].z * uu * (1 - uu) * (1 - uu);  
   float part3 = controlP[3].z * (1 - uu) * (1 - uu) * (1 - uu);     
   return part0 + part1 + part2 + part3;   
}  
void createCurve0(PointCordTypeDef originPoint[4],int originCount,vector<PointCordTypeDef> &curvePoint)
{  
    //生成4控制点，产生贝塞尔曲线  
	PointCordTypeDef controlPoint2[4];
	int i=0;
	controlPoint2[0] = originPoint[i];  
	controlPoint2[1] = originPoint[i+1];
	controlPoint2[2] = originPoint[i+2];
	controlPoint2[3] =originPoint[i+3];
	float u = 1;  
	while(u >= 0){  
		float px = bezier3funcX(u,controlPoint2);  
		float py = bezier3funcY(u,controlPoint2);
		float pz= bezier3funcZ(u,controlPoint2);
		//u的步长决定曲线的疏密  
		u -= 0.005;  
		PointCordTypeDef tempP;
		tempP.cvPoint(px,py,pz);
		//存入曲线点   
		curvePoint.push_back(tempP);  
	}     
}  
void createCurve1(PointCordTypeDef originPoint[4],int originCount,vector<PointCordTypeDef> &curvePoint)
{  
    //控制点收缩系数 ，经调试0.6较好，CvPoint是opencv的，可自行定义结构体(x,y)  
    float scale = 0.6;  
    PointCordTypeDef *midpoints=new PointCordTypeDef[originCount];  
    //生成中点       
    for(int i = 0 ;i < originCount-1 ; i++){      
        int nexti = i + 1;  
        midpoints[i].x = (originPoint[i].x + originPoint[nexti].x)/2.0;  
        midpoints[i].y = (originPoint[i].y + originPoint[nexti].y)/2.0;  
		midpoints[i].z = (originPoint[i].z + originPoint[nexti].z)/2.0;
    }      
    //平移中点  
    PointCordTypeDef *extrapoints=new PointCordTypeDef[2 * originCount];   
    for(int i = 1 ;i < originCount-1 ; i++)
	{  
         int nexti = i + 1;  
         int backi = i - 1;  
         PointCordTypeDef midinmid;  
         midinmid.x = (midpoints[i].x + midpoints[backi].x)/2.0;  
         midinmid.y = (midpoints[i].y + midpoints[backi].y)/2.0;  
		 midinmid.z = (midpoints[i].z + midpoints[backi].z)/2.0;
         float offsetx = originPoint[i].x - midinmid.x;  
         float offsety = originPoint[i].y - midinmid.y;
		 float offsetz = originPoint[i].z - midinmid.z;
         int extraindex = 2 * i;  
         extrapoints[extraindex].x = midpoints[backi].x + offsetx;  
         extrapoints[extraindex].y = midpoints[backi].y + offsety;
		 extrapoints[extraindex].z = midpoints[backi].z + offsetz;  
         //朝 originPoint[i]方向收缩   
         float addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;  
         float addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;
		 float addz = (extrapoints[extraindex].z - originPoint[i].z) * scale;  
         extrapoints[extraindex].x = originPoint[i].x + addx;  
         extrapoints[extraindex].y = originPoint[i].y + addy;  
         extrapoints[extraindex].z = originPoint[i].z + addz;  
         int extranexti = extraindex + 1;  
         extrapoints[extranexti].x = midpoints[i].x + offsetx;  
         extrapoints[extranexti].y = midpoints[i].y + offsety;
		 extrapoints[extranexti].z = midpoints[i].z + offsetz;
         //朝 originPoint[i]方向收缩   
         addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;  
         addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;
		 addz = (extrapoints[extranexti].z - originPoint[i].z) * scale;  
         extrapoints[extranexti].x = originPoint[i].x + addx;  
         extrapoints[extranexti].y = originPoint[i].y + addy;  
         extrapoints[extranexti].z = originPoint[i].z + addz;  
    }      
     
	PointCordTypeDef controlPoint1[3];  
	//生成3控制点，产生第一段贝塞尔曲线  
	    int i1=0;
		controlPoint1[0] = originPoint[i1];  
		int extraindex1 = i1+2;  
		controlPoint1[1] = extrapoints[extraindex1];  
		int nexti1 = i1 + 1;  
		controlPoint1[2] = originPoint[nexti1];      
		float u1 = 1;  
		while(u1 >= 0){  
			float px = bezier2funcX(u1,controlPoint1);  
			float py = bezier2funcY(u1,controlPoint1);
			float pz= bezier2funcZ(u1,controlPoint1);
			//u的步长决定曲线的疏密  
			u1 -= 0.005;  
			PointCordTypeDef tempP;
			tempP.cvPoint(px,py,pz);
			//存入曲线点   
			curvePoint.push_back(tempP);  
		}      
	
    PointCordTypeDef controlPoint2[4];  
    //生成4控制点，产生贝塞尔曲线  
    for(int i = 1 ;i < originCount-2 ; i++){  
           controlPoint2[0] = originPoint[i];  
           int extraindex = 2 * i;  
           controlPoint2[1] = extrapoints[extraindex + 1];  
           int extranexti = extraindex + 2;  
           controlPoint2[2] = extrapoints[extranexti];  
           int nexti = i + 1;  
           controlPoint2[3] = originPoint[nexti];      
           float u = 1;  
           while(u >= 0){  
               float px = bezier3funcX(u,controlPoint2);  
               float py = bezier3funcY(u,controlPoint2);
			   float pz= bezier3funcZ(u,controlPoint2);
               //u的步长决定曲线的疏密  
               u -= 0.005;  
               PointCordTypeDef tempP;
			   tempP.cvPoint(px,py,pz);
               //存入曲线点   
               curvePoint.push_back(tempP);  
           }      
    }  
		PointCordTypeDef controlPoint3[3];  
	//生成3控制点，产生最后一段贝塞尔曲线  
	    int i3=originCount-2;
		controlPoint3[0] = originPoint[i3];  
		int extraindex3 = i3+3;  
		controlPoint3[1] = extrapoints[extraindex3];  
		int nexti3 = i3 + 1;  
		controlPoint3[2] = originPoint[nexti3];      
		float u3 = 1;  
		while(u3 >= 0){  
			float px = bezier2funcX(u3,controlPoint3);  
			float py = bezier2funcY(u3,controlPoint3);
			float pz= bezier2funcZ(u3,controlPoint3);
			//u的步长决定曲线的疏密  
			u3 -= 0.005;  
			PointCordTypeDef tempP;
			tempP.cvPoint(px,py,pz);
			//存入曲线点   
			curvePoint.push_back(tempP);  
		}    
		delete[]midpoints;
		delete[]extrapoints;
}  
//void createCurve(PointCordTypeDef originPoint[4],int originCount,vector<PointCordTypeDef> &curvePoint)
//{  
//    //控制点收缩系数 ，经调试0.6较好，CvPoint是opencv的，可自行定义结构体(x,y)  
//    float scale = 0.6;  
//    PointCordTypeDef *midpoints=new PointCordTypeDef[originCount];  
//    //生成中点       
//    for(int i = 0 ;i < originCount ; i++){      
//        int nexti = (i + 1) % originCount;  
//        midpoints[i].x = (originPoint[i].x + originPoint[nexti].x)/2.0;  
//        midpoints[i].y = (originPoint[i].y + originPoint[nexti].y)/2.0;  
//    }      
//      
//    //平移中点  
//    PointCordTypeDef *extrapoints=new PointCordTypeDef[2 * originCount];   
//    for(int i = 0 ;i < originCount ; i++)
//	{  
//         int nexti = (i + 1) % originCount;  
//         int backi = (i + originCount - 1) % originCount;  
//         PointCordTypeDef midinmid;  
//         midinmid.x = (midpoints[i].x + midpoints[backi].x)/2.0;  
//         midinmid.y = (midpoints[i].y + midpoints[backi].y)/2.0;  
//         int offsetx = originPoint[i].x - midinmid.x;  
//         int offsety = originPoint[i].y - midinmid.y;  
//         int extraindex = 2 * i;  
//         extrapoints[extraindex].x = midpoints[backi].x + offsetx;  
//         extrapoints[extraindex].y = midpoints[backi].y + offsety;  
//         //朝 originPoint[i]方向收缩   
//         int addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;  
//         int addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;  
//         extrapoints[extraindex].x = originPoint[i].x + addx;  
//         extrapoints[extraindex].y = originPoint[i].y + addy;  
//           
//         int extranexti = (extraindex + 1)%(2 * originCount);  
//         extrapoints[extranexti].x = midpoints[i].x + offsetx;  
//         extrapoints[extranexti].y = midpoints[i].y + offsety;  
//         //朝 originPoint[i]方向收缩   
//         addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;  
//         addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;  
//         extrapoints[extranexti].x = originPoint[i].x + addx;  
//         extrapoints[extranexti].y = originPoint[i].y + addy;  
//           
//    }      
//      
//    PointCordTypeDef controlPoint[4];  
//    //生成4控制点，产生贝塞尔曲线  
//    for(int i = 0 ;i < originCount ; i++){  
//           controlPoint[0] = originPoint[i];  
//           int extraindex = 2 * i;  
//           controlPoint[1] = extrapoints[extraindex + 1];  
//           int extranexti = (extraindex + 2) % (2 * originCount);  
//           controlPoint[2] = extrapoints[extranexti];  
//           int nexti = (i + 1) % originCount;  
//           controlPoint[3] = originPoint[nexti];      
//           float u = 1;  
//           while(u >= 0){  
//               int px = bezier3funcX(u,controlPoint);  
//               int py = bezier3funcY(u,controlPoint);  
//               //u的步长决定曲线的疏密  
//               u -= 0.005;  
//               PointCordTypeDef tempP = cvPoint(px,py,);  
//               //存入曲线点   
//               curvePoint.push_back(tempP);  
//           }      
//    }  
//}  
//三次贝塞尔曲线  

int main(int argc, char *argv[])

{
	//**for release
	if( argc < 2 )
	{
		cerr << "Usage: " << endl;
		cerr << "zxhcaeDMPBesselCurve	refpath resultpath" << endl;
		return -1;
	}
	char *REFRFileNameHrd =argv[1];
	char *chResultPathName =argv[2];
	//**for release
	//**for debug
	/*char *REFRFileNameHrd ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel0/point";
	char *chResultPathName ="F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel0/RefFromPointBessel1";*/
	//**for debug
	int iREFRFileLen = strlen(REFRFileNameHrd) + strlen("S.txt");
	// define variables 
	char *chLetter[4]={"S.txt","A.txt","B.txt","E.txt"};
	//read the four points selected from reference (pointS,pointA,pointB,pointE)
	char *REFRFileName;
	vector< miiCNode<double, float> > vREFfourPoints;
	PointCordTypeDef originPoint[4];
	for (int i=0;i<4;i++)
	{	
		REFRFileName = (char*)malloc(iREFRFileLen);
		strcpy(REFRFileName, REFRFileNameHrd);
		strcat(REFRFileName, chLetter[i]);
		ifstream str(REFRFileName);
		ReadPointTxt(REFRFileName, vREFfourPoints);
		PointCordTypeDef PTTemp;
		PTTemp.x=vREFfourPoints[i].x;
		PTTemp.y=vREFfourPoints[i].y;
		PTTemp.z=vREFfourPoints[i].z;
		originPoint[i]=PTTemp;
	}
	//order the A B
	float fPointS[3]={originPoint[0].x,originPoint[0].y,originPoint[0].z};
	float fPointA[3]={originPoint[1].x,originPoint[1].y,originPoint[1].z};
	float fPointB[3]={originPoint[2].x,originPoint[2].y,originPoint[2].z};
	float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
	float fSBDist=zxh::VectorOP_Distance(fPointS,fPointB,3);
	if(fSADist>fSBDist)
	{
		PointCordTypeDef PTemp=originPoint[1];
		originPoint[1]=originPoint[2];
		originPoint[2]=PTemp;
		cout<<"Current Order is S B A E;"<<endl;
	}
	vector<PointCordTypeDef>curvePoint;
	int originCount=4;
	createCurve1(originPoint,originCount,curvePoint);
	//createCurve0(originPoint,originCount,curvePoint);
	//vector<miiCNode<double,float>>vCurvePointWorld;
	//miiCNode<double,float> ptemp;
	//for (int i=0;i<curvePoint.size();i++)
	//{
	//	ptemp.x=curvePoint[i].x;
	//	ptemp.y=curvePoint[i].y;
	//	ptemp.z=curvePoint[i].z;
	//	vCurvePointWorld.push_back(ptemp);
	//}
	// save lenthen-path
	int iCurveNameLen = strlen(chResultPathName) + strlen(".txt");
	char *CurveNamevtk = (char*)malloc(iCurveNameLen);
	strcpy(CurveNamevtk, chResultPathName);
	strcat(CurveNamevtk, ".vtk");
	WriteVtk(curvePoint, CurveNamevtk);
	//
	char *CurveNametxt = (char*)malloc(iCurveNameLen);
	strcpy(CurveNametxt, chResultPathName);
	strcat(CurveNametxt, ".txt");
	WriteTxt(curvePoint, CurveNametxt);
	cout << "Bessel curve fitting is down!" << endl;
}