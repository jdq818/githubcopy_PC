
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

#ifndef __point3D_H__
#define __point3D_H__

#include <cmath>

namespace jdq2017 {

struct point3D
{
  double _x, _y, _z;
  point3D() : _x(0.0), _y(0.0), _z(0.0) {}
  point3D(double x, double y, double z) : _x(x), _y(y), _z(z) {}
  point3D(const point3D &p) : _x(p._x), _y(p._y), _z(p._z) {}
  point3D &operator=(const point3D &p) { 
    if (this != &p) _x=p._x, _y=p._y, _z=p._z;
	return *this;
  }
  double operator[](int i) const { 
    if (i==0) return _x;
    else if (i==1) return _y; 
    else if (i==2) return _z; 
    else return 0;
  }
  point3D operator+(const point3D &p) const { return point3D(_x+p._x, _y+p._y, _z+p._z); }
  point3D operator-(const point3D &p) const { return point3D(_x-p._x, _y-p._y, _z-p._z); }
  point3D operator*(double f) const { return point3D(_x*f, _y*f, _z*f); }
  point3D operator/(double f) const { return point3D(_x/f, _y/f, _z/f); }
  point3D &operator+=(const point3D &p) { _x+=p._x, _y+=p._y, _z+=p._z; return *this; }
  double dot(const point3D &p) const { return _x*p._x+_y*p._y+_z*p._z; }
  double length() const { return std::sqrt(_x*_x+_y*_y+_z*_z); }
  void normalize() {
	  double f = 1.0/length();
	  _x /= f;
	  _y /= f;
	  _z /= f;
  }
};

point3D operator*(double f, const point3D &p)
{
	return p*f;
}

}
#endif

