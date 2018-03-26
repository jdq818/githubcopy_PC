#ifndef miiMinHeap_h
#define miiMinHeap_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#define INF	8000000000000 

using namespace std;

template <class T1 = double, class T2 = int>
// T1-the coordinate value, T2-the coordinate
class miiCNode
{
public:
	miiCNode(){};
	~miiCNode(){};
	T1 val;
	T2 x, y, z;
	void SetValueFrom( miiCNode<T1,T2> & src ) ;
};
template <class T1, class T2>
void miiCNode<T1,T2>::SetValueFrom( miiCNode<T1,T2> & src ) 
{
	val = src.val ; 
	x = src.x ; 
	y = src.y ; 
	z = src.z ;
};
 

template <class T1 = double, class T2 = int>
class miiMinHeap
{
public:
	miiMinHeap();
	miiMinHeap(vector<miiCNode<T1, T2>> &vArray);
	~miiMinHeap();

	bool BuildMinHeap(vector< miiCNode<T1, T2> > &vArray);
 	bool HeapSort(vector< miiCNode<T1, T2> > &vArray);
	miiCNode<T1, T2> HeapExtractMin(vector< miiCNode<T1, T2> > &vArray);
// 	bool HeapMinimum(vector< miiCNode<T> > &vArray);
 	bool HeapDecreaseKey(vector< miiCNode<T1, T2> > &vArray, int nIdx, miiCNode<T1, T2> iNewNode);
 	bool MinHeapInsert(vector< miiCNode<T1, T2> > &vArray, miiCNode<T1, T2> iNewNode);
	bool MinHeapify(vector< miiCNode<T1, T2> > &vArray, int nNode);
		bool MinHeapDelete(vector< miiCNode<T1, T2> > &vArray, int nIdx);
private:
	int m_nHeapSize;

	_inline int LeftNode(int nNode) { return (nNode * 2);}
	_inline int RightNode(int nNode) { return (nNode * 2 + 1);}
	_inline int ParentNode(int nNode) { return (nNode / 2);}
};

template <class T1, class T2>
miiMinHeap<T1, T2>::miiMinHeap()
{
}

template <class T1, class T2>
miiMinHeap<T1, T2>::miiMinHeap(vector<miiCNode<T1, T2>> &vArray)
{
	if (vArray.size() <= 1)
	{
		cerr << "The heap size cannot be less than 1!" << endl;
		exit(1);
	}

	// build the min-heap
	BuildMinHeap(vArray);
}

template <class T1, class T2>
miiMinHeap<T1, T2>::~miiMinHeap()
{

}
template <class T1, class T2>
bool miiMinHeap<T1, T2>::BuildMinHeap(vector< miiCNode<T1, T2> > &vArray)
{
	if (vArray.size() <= 1)
	{
		cerr << "The heap size cannot be less than 1!" << endl;
		return false;
	}

	// get the heap size
	m_nHeapSize = (int)vArray.size() - 1;

	int i = m_nHeapSize / 2;

	while (i > 0)
	{
		MinHeapify(vArray, i);
		i--;
	}

	return true;
}

template <class T1, class T2>
bool miiMinHeap<T1, T2>::MinHeapify(vector< miiCNode<T1, T2> > &vArray, int nNode)
{
	if (vArray.size() <= 1)
	{
		cerr << "The heap size cannot be less than 1!" << endl;
		return false;
	}

	// get the index of left node
	int nLeft = LeftNode(nNode);

	// get the index of right node
	int nRight = RightNode(nNode);

	// define a variable for saving the index of the minimal node
	int nMinNode = 0;

	if (nLeft <= m_nHeapSize)
	{
		if (vArray[nLeft].val < vArray[nNode].val)
		{
			nMinNode = nLeft;
		}
		else
		{
			nMinNode = nNode;
		}
	} 


	if (nRight <= m_nHeapSize)
	{
		if (vArray[nRight].val < vArray[nMinNode].val)
		{
			nMinNode = nRight;
		}
	} 

	if (nMinNode != 0 && nMinNode != nNode)
	{
		miiCNode<T1, T2> vTemp = vArray[nNode];
		vArray[nNode] = vArray[nMinNode];
		vArray[nMinNode] = vTemp;

		MinHeapify(vArray, nMinNode);
	}

	return true;
}

template <class T1, class T2>
bool miiMinHeap<T1, T2>::HeapSort(vector< miiCNode<T1, T2> > &vArray)
{
	if (vArray.size() <= 0)
	{
		std::cerr << "The heap size cannot be less than 1!" << endl;
		return false;
	}

	BuildMinHeap(vArray);

	miiCNode<T> iTempNode;

	for (int i = m_nHeapSize; i >= 2; i--)
	{
		iTempNode = vArray[1];
		vArray[1] = vArray[i];
		vArray[i] = iTempNode;

		m_nHeapSize--;

		MinHeapify(vArray, 1);
	}

	m_nHeapSize = vArray.size() - 1;

	return true;
}

template <class T1, class T2>
miiCNode<T1, T2> miiMinHeap<T1, T2>::HeapExtractMin(vector< miiCNode<T1, T2> > &vArray)
{
	if (vArray.size() == 1)
	{
		return vArray[0];
	}
	miiCNode<T1, T2> iMinNode = vArray[1];

	vArray[1] = vArray[vArray.size() - 1];

	m_nHeapSize--;

	MinHeapify(vArray, 1);

	vArray.erase(vArray.end() - 1);

	return iMinNode;
}

template <class T1, class T2>
bool miiMinHeap<T1, T2>::HeapDecreaseKey(vector< miiCNode<T1, T2> > &vArray, int nIdx, miiCNode<T1, T2> iNewNode)
{
	if (vArray.size() <= 0)
	{
		cerr << "The heap size cannot be less than 1!" << endl;
		return false;
	}

	if (iNewNode.val > vArray[nIdx].val)
	{
		cerr << "New key value is larger than current key value" << endl;
		return false;
	}

	vArray[nIdx] = iNewNode;

	miiCNode<T1, T2> iTempNode;
	while (nIdx > 1 && vArray[ParentNode(nIdx)].val > vArray[nIdx].val)
	{
		iTempNode = vArray[ParentNode(nIdx)];
		vArray[ParentNode(nIdx)] = vArray[nIdx];
		vArray[nIdx] = iTempNode;

		nIdx = ParentNode(nIdx);
	}

	return true;
}

	template <class T1, class T2>
bool miiMinHeap<T1, T2>::MinHeapDelete(vector< miiCNode<T1, T2> > &vArray, int nIdx)
{
	if (vArray.size() <= 0)
	{
		cerr << "The heap size cannot be less than 1!" << endl;
		return false;
	}

	miiCNode<T1, T2> iMinNode = vArray[1];
	HeapDecreaseKey(vArray, nIdx, iMinNode);
	HeapExtractMin(vArray);
	return true;
}

template <class T1, class T2>
bool miiMinHeap<T1, T2>::MinHeapInsert(vector< miiCNode<T1, T2> > &vArray, miiCNode<T1, T2> iNewNode)
{
	if (vArray.size() <= 0)
	{
		std::cerr << "The heap size cannot be less than 1!" << endl;
		return false;
	}

	m_nHeapSize++;

	vArray.push_back(iNewNode);
	vArray[m_nHeapSize].val = INF;

	HeapDecreaseKey(vArray, m_nHeapSize, iNewNode);

	return true;
}

#endif