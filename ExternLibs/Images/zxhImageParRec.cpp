
#include "zxhImageParRec.h"
#include <string>
#include <sstream>
#include <iostream>
//
////namespace zxh
////{
//
////constructor
////template<class ZXH_PixelTypeDefault>
//zxhImageParRec::zxhImageParRec(void)
//{
//}
//
////template<class ZXH_PixelTypeDefault>
//zxhImageParRec::~zxhImageParRec(void)
//{
//}
//
//// copy from nminReadParRecHamburg.c
////template<class ZXH_PixelTypeDefault>
//int zxhImageParRec::GrepLineFromFilePlus(int iplus,FILE *fpascii,char *matchstring,char *matchline)
//{
//  char dummystring[max_buffer_length];
//  int match = 1; // match turns to zero if matchstring is found in ascii file
//  char *retValue;
//
//  if (!fpascii) {
//    fprintf(stderr, "__grepLineFromFile : Invalid file pointer\n");
//    return -1;
//  }
//  fseek(fpascii, 0, 0);	// start at first line of file
//  retValue = dummystring; // any dummy pointer
//  while (!feof(fpascii) && retValue && match) {
//    retValue = fgets(dummystring, max_buffer_length, fpascii);
//    match    = strncmp( dummystring, matchstring, strlen(matchstring));
//  }
//  if (feof(fpascii) || !retValue) { // returns non-zero is eof is reached
//    return -1;
//  }
//  if (!match) {
//    for(int i=0;i<iplus;++i)
//		fgets(dummystring, max_buffer_length, fpascii);
//    strcpy(matchline, dummystring);
//    return 1;
//  }
//
//  return 1;
//}
//
//// copy from nminReadParRecHamburg.c
////template<class ZXH_PixelTypeDefault>
//int zxhImageParRec::GrepLineFromFile(FILE *fpascii,char *matchstring,char *matchline)
//{
//  char dummystring[max_buffer_length];
//  int match = 1; // match turns to zero if matchstring is found in ascii file
//  char *retValue;
//
//  if (!fpascii) {
//    fprintf(stderr, "__grepLineFromFile : Invalid file pointer\n");
//    return -1;
//  }
//  fseek(fpascii, 0, 0);	// start at first line of file
//  retValue = dummystring; // any dummy pointer
//  while (!feof(fpascii) && retValue && match) {
//    retValue = fgets(dummystring, max_buffer_length, fpascii);
//    match    = strncmp( dummystring, matchstring, strlen(matchstring));
//  }
//  if (feof(fpascii) || !retValue) { // returns non-zero is eof is reached
//    return -1;
//  }
//  if (!match) {
//    strcpy(matchline, dummystring);
//    return 1;
//  }
//
//  return 1;
//}
//
//// copy form nminReadParRecHamburg
////template<class ZXH_PixelTypeDefault>
//bool zxhImageParRec::OpenImageParRec(const char *pFilePar,const char* pFileRec)
//{
//	// begin extract the information from the Par file
//	float ap,	//for coronal(slice_orientation=3)
//		fh,		//for axail(slice_orientation=1)
//		rl;		//for sagital(slice_orientation=2)
//	int slice_orientation;//  slice orientation ( TRA=1/SAG=2/COR=3 )        (integer)
//	float trig_time0,trig_time1,spacing_x,spacing_y;
//	int no_phases,no_slices,x_reso,y_reso;
//	int slice_first,phase_first;
//	// Image file
//	FILE *fp = 0, *fp_par = 0;
//	int bitsPerPixel = 16;
//	char  line[max_buffer_length];
//
//	// Open file
//	fp_par = fopen(pFilePar, "r");
//
//	if (fp_par == 0)
//	{
//		fprintf(stderr, "loadParRecImage: Unable to open file %s\n", pFilePar);
//		return false;;
//	}
//
//	//get image series information
//	std::istringstream istrs;
//	std::string sline;
//	int ipos_data;
//	if(GrepLineFromFile( fp_par, ".    Max. number of cardiac phases", line)==1)
//        sscanf(line, "%*s %*s %*s %*s %*s %*s %*s  %d\n", &no_phases);
//	if(GrepLineFromFile( fp_par, ".    Max. number of slices/locations", line)==1)
//		sscanf(line, "%*s %*s %*s %*s %*s %*s  %d\n", &no_slices);
//	/*GrepLineFromFile( fp_par,".    Image pixel size [8 or 16 bits]",line);
//	sline	= line;
//	ipos_data	= sline.find(":",0);
//	istrs.str(sline.substr(ipos_data+1,sline.length()-ipos_data-1));
//	istrs>>bitsPerPixel;
//	GrepLineFromFile( fp_par, ".    Recon resolution (x, y)",line);
//	sline	= line;
//	ipos_data	= sline.find(":",0);
//	istrs.str(sline.substr(ipos_data+1,sline.length()-ipos_data-1));
//	istrs>>x_reso;istrs>>y_reso;  this is wrong when image recontructed*/
//	GrepLineFromFile( fp_par, ".    FOV (ap,fh,rl) [mm]",line);
//	sline	= line;
//	ipos_data	= (int)sline.find(":",0);
//	istrs.str(sline.substr(ipos_data+1,sline.length()-ipos_data-1));
//	istrs>>ap;istrs>>fh;istrs>>rl;
//
//	GrepLineFromFilePlus(2,fp_par, "#sl ec dyn ph ty", line);
//	///           1   2   3   4  5   6   7   8    9  0 1   2   3    4   5  6   7   8    9   0  1   2  3    4   5   6  7   8  9   0   1  2   3
//	sscanf(line, "%d %*s %*s %d %*s %*s %*s %*s %*s %d %d %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d %*s %*s %*s %*s %*s %*s %f",
//		&slice_first,&phase_first,&x_reso,&y_reso,&slice_orientation,&trig_time0);
//	if(no_phases*no_slices>1)
//	{
//		GrepLineFromFilePlus(3,fp_par, "#sl ec dyn ph ty", line);
//		///           1   2   3   4  5   6   7   8    9  0 1   2   3    4   5  6   7   8    9   0  1   2  3    4   5   6  7   8  9  0   1  2   3
//		sscanf(line, "%d %*s %*s %d %*s %*s %*s %*s %*s %d %d %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d %*s %*s %f %f %*s %*s %f",
//			&slice_first,&phase_first,&x_reso,&y_reso,&slice_orientation,&spacing_x,&spacing_y,&trig_time1);
//		if(phase_first>1)
//		{
//			m_afSpacing[3]=trig_time1-trig_time0;
//		}//else zxhtodo
//	}
//	fclose(fp_par);
//
//	m_afSpacing[0]=spacing_x;
//	m_afSpacing[1]=spacing_y;
//	switch(slice_orientation)
//	{
//	case 1://axail
//		m_afSpacing[2]=fh/(float)no_slices;
//		break;
//	case 2://sagital
//		m_afSpacing[2]=rl/(float)no_slices;
//		break;
//	default://3 coronal
//		m_afSpacing[2]=ap/(float)no_slices;
//	}
//
//	// Initialize number of in-plane pixels and total Pixels
//	m_aiSize[0]=x_reso;
//	m_aiSize[1]=y_reso;
//	m_aiSize[2]=no_slices;
//	m_aiSize[3]=no_phases;
//	m_iResolution=x_reso*y_reso;
//	m_iVolume=m_iResolution*no_slices;
//
//	//m_afSpacing[3]=(float)miniseconds_cardiac_phases/(float)m_aiSize[3];
//	if(m_aiSize[3]>1)//&&m_aiSize[2]>1)
//		m_iDimension=4;
//	else if(m_aiSize[3]==1&&m_aiSize[2]>1)
//		m_iDimension=3;
//	else if(m_aiSize[3]==1&&m_aiSize[2]==1)
//		m_iDimension=2;
//	// new data
//	m_pImageData	= new ParRecType[m_iVolume*m_aiSize[3]];
//	if(m_pImageData==0)
//	{
//		std::cerr<<"Open Image:"<<pFilePar<<"Error:can not new image data\n";
//		return false;
//	}
//	char * cpBuffer;
//	ParRecType  * ipBuffer;
//	char *buffer	= new char[m_aiSize[0]*m_aiSize[1]+1];//only for when bitsPerPixel==8
//
//	//get raw data
//	fp = fopen(pFileRec, "r");
//	if (fp == 0)
//	{
//		fprintf(stderr, "loadParRecImage: Unable to open file %s\n", pFileRec);
//		return false;;
//	}
//	// store data
//	if(slice_first>1)
//	{
//		for (int iph=0;iph<no_phases;++iph)
//		{
//			for (int islc=0;islc<no_slices;++islc)
//			{
//				unsigned long iOffset	= iph*m_iVolume+islc*m_iResolution;
//				if(bitsPerPixel==16)
//				{
//					cpBuffer	= (char*)(m_pImageData+iOffset);
//					fseek(fp,iOffset*2,SEEK_SET);
//					fread(cpBuffer,2,m_iResolution,fp);
//				}
//				else //bitsPerPixel==8
//				{
//					fseek(fp,iOffset,SEEK_SET);
//					ipBuffer	= m_pImageData+iOffset;
//					fread(buffer,1,m_iResolution,fp);
//					for(unsigned int i=0;i<m_iResolution;++i)
//						ipBuffer[i]	= (ParRecType)buffer[i];
//				}
//			}
//		}
//	}else //phase_first, normally, phase_first
//	{
//		for(int islc=0;islc<no_slices;++islc)
//		{
//			for(int iph=0;iph<no_phases;++iph)
//			{
//				unsigned long int iOffsetOfFile	= (iph+islc*no_phases)*m_iResolution;
//				unsigned long int iOffsetOfImage	= (iph*m_iVolume+islc*m_iResolution);
//				if(bitsPerPixel==16)
//				{
//					cpBuffer	= (char*)(m_pImageData+iOffsetOfImage);
//					fseek(fp,iOffsetOfFile*2,SEEK_SET);
//					fread(cpBuffer,2,m_iResolution,fp);
//				}
//				else //bitsPerPixel==8
//				{
//					fseek(fp,iOffsetOfFile,SEEK_SET);
//					ipBuffer	= m_pImageData+iOffsetOfImage;
//					fread(buffer,1,m_aiSize[0]*m_aiSize[1],fp);
//					for(int i=0;i<m_aiSize[0]*m_aiSize[1];++i)
//						ipBuffer[i]	= (ParRecType)buffer[i];
//				}
//			}
//		}
//	}
//	delete []buffer;
//	// Close file
//	fclose(fp);
//	//fprintf(stdout,"Open Par_Rec:%s Image Sucess!\n",pFilePar);
//	std::cout<<"Open Par_Rec:"<<pFilePar<<" Image Success!\n" ;
//	return true;
//}
////template<class ZXH_PixelTypeDefault>
//bool zxhImageParRec::OpenImage(const char*pFilename)
//{
//	std::string sFile	= pFilename;
//	unsigned int ipos=sFile.length()-4;
//	std::string sPar	= sFile.substr(0,ipos)+".PAR";
//	std::string sRec	= sFile.substr(0,ipos)+".REC";
//	return OpenImageParRec(sPar.c_str(),sRec.c_str());
//}

//}//end of namespace






