
/*=========================================================================

  Program:   class jdqPoint
  Author:	 Jia, Dengqiang
  Module:    $RCSfle: jdqPoint.h    $
  Language:  C++
  Date:      $Date: From  2017-10$
  Version:   $Revision: 1.0$
			 $Revision: 1.1$

  Update log: 2017-10-22 some new member functions were included

=========================================================================*/
//dijkstra.h
//Classes to determine connections between two paths
//Jia, Dengqiang
//cat2008util.h
//Utilities for the Coronary Artery Tracking Challenge 2008 software
//See http://cat08.bigr.nl/ for details
//Theo van Walsum
//Michiel Schaap
//Coert Metz

#ifndef __Util_H__
#define __Util_H__

#include <string>
#include <sstream>
#include <iterator>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

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

namespace jdq2017 {

//readScores
//**************
//Function to read interobserver overlap variabilities
//  for score computation from a text file
//IN
//  filename: filename to read scores from
//OUT:
//  overlap_io: interobserver overlap variabilities for all vessels
bool readScores(const char *filename,
                std::map< int, std::map< int, std::vector<double> > > &overlap_io)
{
  std::ifstream  input(filename);

  if (!input) return false;

  int linenr = 0;
  int ds_nr = -1;
  int vs_nr = 3;

  // Read interobserver overlap variabilities for all vessels
  while (input)
  {
    ++linenr;
    std::string line;
    if (!std::getline(input, line)) break;
    if (line[0] == '#') continue;
    std::istringstream linestream(line.c_str());
    std::vector<double> nrs;
    int ds, vs;
    linestream >> ds >> vs;
    //Split line on spaces/tabs/newlines
    std::copy(std::istream_iterator<double>(linestream), 
      std::istream_iterator<double>(), std::back_inserter(nrs));
    //Check input
    if (nrs.size() != 3)
    {
      std::cerr << "Error reading line " << linenr << " from file " 
        << filename  << std::endl;
      return false;
    }
    if (overlap_io[ds][vs] != std::vector<double>() )
    {
      std::cerr << "Error reading line " << linenr << " from file " 
        << filename << ": duplicate entry" << std::endl;
      return false;
    }
   overlap_io[ds][vs] = nrs;
  }
  return true;
}

//pathLength
//**************
//Function to compute Euclidian length of a path/centerline
//IN
//  path: vector containing path points
//OUT:
//  length of the path (sum of all line segment lengths)
template <typename Point>
double pathLength(const std::vector<Point> &path)
{
  if (path.size() < 2) return 0.0;
  double cumlen = 0.0;
  for (typename std::vector<Point>::const_iterator i = path.begin()+1;
       i != path.end(); 
       ++i)
       cumlen += (*(i-1) - *i).length();
  return cumlen;
}

//readReference
//**************
//Read reference standard centerline from textfile
//IN
//  filename: filename to read centerline from
//OUT:
//  ref: vector with points of the path
//  rad: vector with radii attached to path points
//  io: vector with interobserver accuracy variability attached to path points
template <typename Point>
bool readReference(const char *filename,
                   std::vector<Point> &ref,
                   std::vector<double> &rad,
                   std::vector<double> &io)
{
    std::ifstream  input(filename);
    
    int linenr = 0;
    while (input)
    {
      ++linenr;
      std::string line;
      //Read line from textfile
      if (!std::getline(input, line)) break;
      if (line[0] == '#') continue;
      std::istringstream linestream(line.c_str());
      std::vector<double> nrs;
      //Split line on spaces/tabs/newlines
      std::copy(std::istream_iterator<double>(linestream), 
                std::istream_iterator<double>(), std::back_inserter(nrs));
      //Check input
      if (nrs.size() != 5)
      {
        std::cerr << "Error reading line " << linenr << " from file " 
                  << filename << std::endl;
        return false;
      }
      //Add point, radius and io to corresponding output variable
      ref.push_back(Point(nrs[0], nrs[1], nrs[2]));
      rad.push_back(nrs[3]);
      io.push_back(nrs[4]); 
    }
    return true;
}

//readCenterline
//**************
//Read centerline from textfile
//IN
//  filename: filename to read centerline from
//  nrElementsPerLine: number of elements per line (if > 4, every element on the
//    line after the third one will be skipped
//OUT:
//  cl: vector with points of the path
//  io: vector with interobserver accuracy variability attached to path points
template <typename Point>
bool readCenterline(const char *filename,
                    std::vector<Point> &cl,
                    int nrElementsPerLine = 3)
{
    std::ifstream  input(filename);
    
    int linenr = 0;
    while (input)
    {
      ++linenr;
      std::string line;
      //Read line from file
      if (!std::getline(input, line)) break;
      if (line[0] == '#') continue;
      std::istringstream linestream(line.c_str());
      std::vector<double> nrs;
      //Split line on spaces/tabs/newlines
      std::copy(std::istream_iterator<double>(linestream), 
                std::istream_iterator<double>(), std::back_inserter(nrs));
      //Check input: skip empty lines ...
      if (nrs.size() == 0) continue;
      // Check input: count number of elements on line ...
      if (nrs.size() != nrElementsPerLine)
      {
        std::cerr << "Error reading line " << linenr << " from file " 
                  << filename << std::endl;
        return false;
      }
      //Add point to centerline
      cl.push_back(Point(nrs[0], nrs[1], nrs[2]));
    }
    return true;
}

//readCenterlinevtk
//**************
//Read centerline from textfile
//IN
//  filename: filename to read centerline from
//  nrElementsPerLine: number of elements per line (if > 4, every element on the
//    line after the third one will be skipped
//OUT:
//  cl: vector with points of the path
//  io: vector with interobserver accuracy variability attached to path points
template <typename Point>
bool readCenterlinevtk(const char *filename,
                    std::vector<Point> &cl,
                    int nrElementsPerLine = 3)
{
	if (filename == NULL)
	{
		cout << "Cannot find VTK-file!" << endl;
		return false;
	}
	if (!cl.empty())
	{
		cl.clear();
	}

	vtkSmartPointer<vtkUnstructuredGridReader> iVtkReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	iVtkReader->SetFileName( filename );
	iVtkReader->Update();

	vtkSmartPointer<vtkUnstructuredGrid> iGridRead = iVtkReader->GetOutput();

	int nPointNum = iGridRead->GetMaxCellSize();

	double dCord[3];
	Point strctTempPoint;

	for (int i = 0; i < nPointNum; i++)
	{
		iGridRead->GetPoint(i, dCord);
		strctTempPoint._x = dCord[0];
		strctTempPoint._y =dCord[1];
		strctTempPoint._z = dCord[2];
		cl.push_back(strctTempPoint);
	}
		//copy of origin
	/*
    std::ifstream  input(filename);
    
    int linenr = 0;
    while (input)
    {
      ++linenr;
      std::string line;
      //Read line from file
      if (!std::getline(input, line)) break;
      if (line[0] == '#') continue;
      std::istringstream linestream(line.c_str());
      std::vector<double> nrs;
      //Split line on spaces/tabs/newlines
      std::copy(std::istream_iterator<double>(linestream), 
                std::istream_iterator<double>(), std::back_inserter(nrs));
      //Check input: skip empty lines ...
      if (nrs.size() == 0) continue;
      // Check input: count number of elements on line ...
      if (nrs.size() != nrElementsPerLine)
      {
        std::cerr << "Error reading line " << linenr << " from file " 
                  << filename << std::endl;
        return false;
      }
      //Add point to centerline
      cl.push_back(Point(nrs[0], nrs[1], nrs[2]));
    }
	*/ 
	//copy of origin
    return true;
}

//startCLFrom
//**************
//Calculate clipping point based on disk at start of reference
//IN
//  s: first point reference standard
//  n: normal of clipping plane at s
//  rad: radius first point reference standard
//  cl: centerline for which clipping needs to be determined
//OUT:
//  index where centerline starts, clipping should be performed up to this point
template <typename Point>
int startCLFrom(const Point &s, const Point &n, double rad,
                const std::vector<Point> &cl)
{
  int startAt = -1;
  typename std::vector<Point>::const_iterator i(cl.begin());
  typename std::vector<Point>::const_iterator j(i+1);
  while ( startAt < 0 && j != cl.end() )
  {
    // test whether line segment intersects plane of disc
    if ( (*i-s).dot(n) * (*j-s).dot(n) < 0.0 )
    {
      // possibly found intersection, linearly interpolate position at disk
      double a = (*i-s).dot(n);
      double b = -1.0*(*j-s).dot(n);
      Point pos = b/(a+b)* *i + a/(a+b)* *j;
      
      // test distance of intersection to start pos: should be less than 2 x 
      // radius
      if ( (pos-s).length() < 2.0*rad )
        startAt = int(j - cl.begin());
    }
    ++i, ++j;
  }
  return (startAt >= 0) ? startAt : 0;
}

//ResamplePaths
//**************
//Structs to resample paths equidistantly
template <typename Point >
struct ResamplePaths
{
  std::vector< Point > &resultPath()  { return _resultPath; }
  std::vector< double > &resultRadius()  { return _resultRadius; }
  //Resampling operator
  void operator() (const std::vector<Point> &path, double samplingDistance);
  
  void operator() (const std::vector<Point> &path, const std::vector<double> &radius, double samplingDistance);
  
private:
  //Resulting resampled centerline
  std::vector< Point > _resultPath;
  std::vector< double > _resultRadius;
};

//operator()
//*************
//Path resample operator
//IN
//  path: path to perform resampling on
//  samplingDistance: samplingDistance output path in mm
template <typename Point >
void ResamplePaths<Point>::
operator() ( const std::vector<Point> &path,
             double samplingDistance)
{
  //Check for valid samplingdistance
  if (samplingDistance == 0.0)
  {
    _resultPath = path;
    return;
  }

  //Resampled result path is stored in _resultPath
  _resultPath.clear();

  double curLen = 0.0;
  double sampleNr = 0;
  int curIdx = 0;

  while (curIdx < int(path.size())-1)
  {
    // resample up to nextIdx
    double nextLen = curLen + (path[curIdx+1]-path[curIdx]).length();
	double cnLen=nextLen-curLen;
    while (sampleNr*samplingDistance < nextLen)
    {
      // linearly interpolate between curIdx at curLen, and curIdx at nextLen
      double dist = sampleNr * samplingDistance;
      double a = (nextLen-dist) / (nextLen - curLen); // weight for curIdx
      double b = (dist - curLen) / (nextLen - curLen);// weight for curIdx+1

      _resultPath.push_back(a*path[curIdx] + b*path[curIdx+1]);
      ++sampleNr;
    }
    curLen = nextLen;
    ++curIdx;
  }
  //std::cout<<_resultPath.size()<<std::endl;
  // add last sample
  _resultPath.push_back(path.back());
}

//operator()
//*************
//Path resample operator
//IN
//  path: path to perform resampling on
// radius: radius at every point in path
//  samplingDistance: samplingDistance output path in mm
template <typename Point >
void ResamplePaths<Point>::
operator() ( const std::vector<Point> &path,
             const std::vector<double> &radius,
             double samplingDistance)
{
  //Check for valid samplingdistance
  if (samplingDistance == 0.0)
  {
    _resultPath = path;
    _resultRadius = radius;
    return;
  }

  //Resampled result path is stored in _resultPath
  _resultPath.clear();
  _resultRadius.clear();

  double curLen = 0.0;
  double sampleNr = 0;
  int curIdx = 0;

  while (curIdx < int(path.size())-1)
  {
    // resample up to nextIdx
    double nextLen = curLen + (path[curIdx+1]-path[curIdx]).length();
    while (sampleNr*samplingDistance < nextLen)
    {
      // linearly interpolate between curIdx at curLen, and curIdx at nextLen
      double dist = sampleNr * samplingDistance;
      double a = (nextLen-dist) / (nextLen - curLen); // weight for curIdx
      double b = (dist - curLen) / (nextLen - curLen);// weight for curIdx+1
      _resultPath.push_back(a*path[curIdx] + b*path[curIdx+1]);
      _resultRadius.push_back(a*radius[curIdx] + b*radius[curIdx+1]); 
      ++sampleNr;
    }
    curLen = nextLen;
    ++curIdx;
  }

  // add last sample
  _resultPath.push_back(path.back());
  _resultRadius.push_back(radius.back());
}

//writeReference
//**************
//Write reference standard to file
//IN
//  filename: file to write to
//  ref: reference centerline to write
//  rad: radii of centerline to write
//  io: io of centerline to write
//OUT:
//  returns true if writing succeeded
template <typename Point>
bool writeReference(const char *filename,
                    std::vector<Point> &ref,
                    std::vector<double> &rad,
                    std::vector<double> &io) {
  // Open file
  std::ofstream file_op(filename);

  // Check if file is open
  if (file_op.is_open()) {
    // Write path to file
    int i=0;
    for (typename std::vector<Point>::const_iterator it = ref.begin(); it != ref.end(); ++it, ++i) {
      Point v = *it;
      file_op << v[0] << " " << v[1] << " " << v[2] << " " << rad[i] << " " << io[i] << std::endl;
    }
    file_op.close();
    return true;
  } else {
    return false;
  }  
}

}

#endif



