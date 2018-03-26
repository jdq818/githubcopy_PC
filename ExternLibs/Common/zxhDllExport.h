
/*=========================================================================

  Program:   ZXH Registration Software
  Author:	 Xiahai Zhuang
  Module:    $RCSfle: zxhDllExport.h    $
  Language:  C++
  Date:      $Date: From  2004-01 $
  Version:   $Revision: 1.0, 2.0 $

=========================================================================*/
// .NAME zxhDllExport - maitain marcro for dllexport
// .SECTION Description
//

#ifndef _zxhDllExport_h
#define _zxhDllExport_h

//#define ZXH_FOR_DLL_EXPORT//if compile dll then define
//#define ZXH_FOR_DLL_INPORT //header for using dll

#ifdef ZXH_FOR_DLL_EXPORT//if compile dll then define
#define ZXH_DLL_EXPORT __declspec( dllexport )// compile dll
#else

#ifdef ZXH_FOR_DLL_INPORT //header for using dll
#define ZXH_DLL_EXPORT __declspec( dllimport )//use dll
#else
#define ZXH_DLL_EXPORT
#endif


 // error message id define, e.g: std::cerr<<"error: system error "<<ZXH_SYSTEMERROR_zxhMetricBase_ResampleTestImagesSpacingPhysIntervel <<"\n"; 
#define ZXH_SYSTEMERROR_zxhMetricBase_ResampleTestImagesSpacingPhysIntervel 1 

#endif //ZXH_FOR_DLL_EXPORT

#endif//_zxhDllExport_h
