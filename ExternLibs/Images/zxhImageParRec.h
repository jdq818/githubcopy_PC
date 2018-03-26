
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhImageParRec.h    $
  Language:  C++
  Date:      $Date: From  2006-11 $
  Version:   $Revision: 2.0 $

=========================================================================*/

#ifndef zxhImageParRec_h
#define zxhImageParRec_h

#include "zxhImageData.h"
#include <stdio.h>

//namespace zxh
//{

///
/// \class zxhImageParRec
/// \brief REC_PAR image data operation
///
/// \ingroup zxhImageData
///
typedef ZXH_PixelTypeDefault ParRecType;
class zxhImageParRec:
	public zxhImageDataT<ParRecType>
{
public:
	///constract
	zxhImageParRec(void){};

	///\return
	virtual ~zxhImageParRec(void){};
	 
//
//	/// open image
//	virtual bool OpenImage(const char*pFilename);
//
//	/// read REC_PAR image, if Rec=0 then uses pFilePar-"PAR"+"REC"
//	bool OpenImageParRec(const char *pFilePar,const char* pFileRec);
//
//
//protected:
//	/// for read REC_PAR
//	int GrepLineFromFilePlus(int iplus,FILE *fpascii,char *matchstring,char *matchline);
//	/// for read REC_PAR
//	int GrepLineFromFile(FILE *fpascii,char *matchstring,char *matchline);
//
//	// image info 0:x,1:y,2:z,3:phase 	int m_aiResolution4D[4];
//	// for spacing[3] is temporal dimension space, unit is minisecond	float m_afSpacing4D[4];
//
//	/// constant value
//	const static int max_buffer_length	= 2048;
//	///constant value
//	const static int miniseconds_cardiac_phases = 800;
};

//};//end of namespace
#endif //zxhImageParRec_h


