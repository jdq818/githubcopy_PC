
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
#ifndef __twoPaths_h__
#define __twoPaths_h__

#include <algorithm>
#include <utility>
#include <vector>
#include <limits>
#include <set>

namespace jdq2017 {

// Typedef for pairing
typedef std::pair< int, int > Connection;

// Functor that determines the length of a connection
// Maintains a pointer to both input paths, and requires
// (binary) subtraction operator and length() member function
// for Point datatype.
template <typename Point>
struct ConnectionLength
{
  ConnectionLength(const std::vector<Point> &A,
    const std::vector<Point> &B) 
    : _a(&A), _b(&B) {}

  double operator() (const Connection &c)
  {
    return ((*_a)[c.first]-(*_b)[c.second]).length();
  }

  const std::vector<Point> *_a;
  const std::vector<Point> *_b;

  ConnectionLength & operator= (const ConnectionLength &rhs)
  {
    _a = rhs._a;
    _b = rhs._b;
    return *this;
  }
};


//
// Class to determine the "connection" between two paths.
template <typename Point>
class Dijkstra {

public:

  // Operator, determines the optimal connection between two sets of
  // nD points, based on either minimal surface or minimal total
  // connection length. Defaults to minimal connection length.
  double operator() (const std::vector< Point > &a,
                     const std::vector< Point > &b,
                     bool minimalSurface = false);


  // Returns the last pairing found.
  const std::vector< Connection > &pairing() const { return _shortestPath; }

  // Determines the distance given a pairing: connection length
  double edgeLength (const std::vector< Point > &a,
                         const std::vector< Point > &b, 
                         const std::vector< Connection > &pairing)
  { return distance(a,b,pairing,false); }

  // Determines the distance given a pairing: total triangle surface area.
  double surfaceArea(const std::vector< Point > &a,
                         const std::vector< Point > &b, 
                         const std::vector< Connection > &pairing)
  { return distance(a,b,pairing,true); }


private:

  // Helper enum for Connection Status
  typedef enum { Done, OnFront, NotVisited } NodeStatus;

  // Found path
  std::vector< Connection >  _shortestPath;
  
  // Returns maximum value for double
  static double distMax() { 
    return std::numeric_limits<double>::max(); 
  }

  // Determines the distance for a given pairing.
  double distance(const std::vector< Point > &a,
    const std::vector< Point > &b, 
    const std::vector< Connection > &pairing, 
    bool minimalSurface);

  // Switches between two distance measures (for edges): triangle area
  // or edge length. 
  double distance(const Point &p1, const Point &p2, const Point & p3,
                  bool minimalSurface)
  {
    return minimalSurface ? triangleArea(p1,p2,p3) : edgeLength(p2,p3);
  }

  // Determines length of an edge. Three arguments are used because of
  // the triangleArea function.
  double edgeLength(const Point & p2,
                    const Point & p3)
  {
    return (p2-p3).length();
  }

  //////////
  // Determines surface area of a triangle (for 2D, 3D, and
  // probably more dimensions ...
  double triangleArea(const Point &p1, 
                      const Point &p2,
                      const Point &p3)
  {
    // Base of the triangle
    Point base21(p2 - p1);
    // Normalized vector along the base
    Point base21norm(base21);
    base21norm.normalize();

    // One of the other sides
    Point side23(p2 - p3);
    // Height vector: part of side23 that is orthogonal to the base
    Point height(side23 - side23.dot(base21norm)*base21norm);

    // Surface area: height*base/2
    return 0.5*height.length()*base21.length();
  }

};

// Function implementing Dijkstra 
template <typename Point>
double Dijkstra<Point>::
operator()( const std::vector< Point > &a,
            const std::vector< Point > &b, bool minimalSurface )
{

  // Initialize local vars
  _shortestPath.clear();

  // Sanity check
  if (a.empty() || b.empty()) return -1.0;

  // --------------------------------------------------
  // Temporary storage: 

  // 2D Array of all nodes, contains node status
  std::vector< std::vector< NodeStatus > > 
    nodeStatus( a.size(), 
                std::vector<NodeStatus>( b.size(), 
                                            NotVisited ) );

  // 2D array of current distances (for Done and OnFront nodes)
  std::vector< std::vector< double > > 
    distanceMap( a.size(), 
                 std::vector<double>( b.size(), distMax() ) );

  // 2D array of pointers to previous node, for tracing back the path
  std::vector< std::vector< unsigned char > > 
    prevPointer( a.size(), 
                 std::vector<unsigned char>( b.size(), 0 ) );

  typedef std::pair< double, Connection > QueueElement;

  std::set< QueueElement > priorityQueue;

  // Initialize queue with first element
  double dist = distance(a[0], b[0], a[0], minimalSurface);
  priorityQueue.insert( QueueElement( dist, Connection(0,0) ) );
  nodeStatus[0][0] = OnFront;
  distanceMap[0][0] = dist;

  // Keep popping nodes while there is still a node to pop, and while
  // the distance to the next node is less then the distance to the
  // end point.
  while ( !priorityQueue.empty() && 
          priorityQueue.begin()->first < distanceMap.back().back() )
    {
      // pop front, and iterate until node is found that is not Done
      // yet, and that has a better distance than the current stored
      // distance.
      Connection queueElem;
      double dist;
      do {
        queueElem = priorityQueue.begin()->second;
        dist = priorityQueue.begin()->first;
        priorityQueue.erase(priorityQueue.begin());
      } while ( !priorityQueue.empty() &&
                nodeStatus[queueElem.first][queueElem.second] == Done );

      // Finished if the queue is empty, and no valid node is found.
      if ( nodeStatus[queueElem.first][queueElem.second] == Done ) break;

      // Finished if we are guaranteed that the path to the end point
      // can not be improved anymore
      if ( dist > distanceMap.back().back() ) break;

      // Process the node: update nodeStatus and distanceMap datastructures.
      nodeStatus[queueElem.first][queueElem.second] = Done;
      distanceMap[queueElem.first][queueElem.second] = dist;

      // Process neighbours:
      if ( queueElem.second < int(b.size())-1 )
        {
          double newDist = dist + distance(b[queueElem.second],
			                               a[queueElem.first],
                                           b[queueElem.second+1], 
                                           minimalSurface);
          
          switch ( nodeStatus[queueElem.first][queueElem.second+1] ) {
          case Done: continue;
          case OnFront:
            if (newDist >= distanceMap[queueElem.first][queueElem.second+1]) continue;
          case NotVisited: // and if better then current on front
            nodeStatus[queueElem.first][queueElem.second+1] = OnFront;
            distanceMap[queueElem.first][queueElem.second+1] = newDist;
            prevPointer[queueElem.first][queueElem.second+1] = 2;
            priorityQueue.insert( QueueElement( newDist, 
                                                Connection( queueElem.first, 
                                                            queueElem.second+1) ) );
          }
        }
      if ( queueElem.first < int(a.size())-1 )
        {
          double newDist = dist + distance(a[queueElem.first],
                                             b[queueElem.second],
                                             a[queueElem.first+1],
                                             minimalSurface);
          
          switch ( nodeStatus[queueElem.first+1][queueElem.second] ) {
          case Done: continue;
          case OnFront:
            if (newDist >= distanceMap[queueElem.first+1][queueElem.second]) continue;
          case NotVisited: // and if better then current on front
            nodeStatus[queueElem.first+1][queueElem.second] = OnFront;
            distanceMap[queueElem.first+1][queueElem.second] = newDist;
            prevPointer[queueElem.first+1][queueElem.second] = 1;
            priorityQueue.insert( QueueElement( newDist, 
                                                Connection( queueElem.first+1,
                                                            queueElem.second ) ) );
          }
        }
    }

  // Trace back path from end.
  std::vector< Connection >  revPath; // reversed path
  revPath.push_back( Connection(int(a.size()-1), int(b.size()-1)) );
  while ( revPath.back().first ||  revPath.back().second)
    {
      switch ( prevPointer[revPath.back().first][revPath.back().second] )
        {
        case 1:
          revPath.push_back( Connection ( revPath.back().first-1, 
                                          revPath.back().second ) );
          break;
        case 2:
          revPath.push_back( Connection ( revPath.back().first, 
                                    revPath.back().second-1 ) );
          break;
        default:
        case 0:
          exit(1);
          // can not happen!
        }
    }

  // Reverse the traced back path.
  std::reverse_copy( revPath.begin(), revPath.end(),
                     back_inserter(_shortestPath) );

  // Return the minimal distance.
  return distanceMap.back().back();
}


template <typename Point>
double Dijkstra<Point>::
distance   (const std::vector< Point > &a,
            const std::vector< Point > &b, 
            const std::vector< Connection > &pairing, 
            bool minimalSurface)
{
  // Determine distance given  a pairing...
  if (pairing.empty() || a.empty() || b.empty()) return -1.0;

  double dist = distance(a[pairing[0].first], b[pairing[0].second], a[pairing[0].first], minimalSurface);

  for (unsigned int i = 1; i < pairing.size(); ++i)
  {
    if (pairing[i-1].first == pairing[i].first)
      dist += distance(b[pairing[i-1].second], a[pairing[i].first],  b[pairing[i].second], minimalSurface);
    else
      dist += distance(a[pairing[i-1].first], b[pairing[i].second], a[pairing[i].first], minimalSurface);
  }
  return dist;
}

}
#endif


