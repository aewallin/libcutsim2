/*

  Main class and structures for DC

  Copyright (C) 2011  Tao Ju

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef OCTREE_H
#define OCTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "GeoCommon.hpp"
#include "eigen.hpp"
#include "dualcontouring_tables.hpp"

// these required for intersection-free algorithm?
// #include "HashMap.hpp"
// #include "intersection.hpp"

// Clamp all minimizers to be inside the cell
//#define CLAMP

// If SOG vertices are relative to each cell
//#define SOG_RELATIVE

//#define EDGE_TEST_CONVEXITY
//#define EDGE_TEST_FLIPDIAGONAL
//#define EDGE_TEST_NEW

//#define TESS_UNIFORM
//#define TESS_NONE

// old definiton was: 0== InternalNode, 1== LeafNode, 2==PseudoLeafNode
enum NodeType { INTERNAL, LEAF, PSEUDOLEAF };

/* Tree nodes */
class OctreeNode {
public:
    OctreeNode(){};
    virtual ~OctreeNode(){};
    virtual NodeType getType() = 0; // 0== InternalNode, 1== LeafNode, 2==PseudoLeafNode
};

class InternalNode : public OctreeNode {
public: // no signs, height, len, or QEF stored for internal node
	OctreeNode * child[8] ;
	InternalNode ()  {
		for ( int i = 0 ; i < 8 ; i ++ )
			child[i] = NULL ;
	};
	virtual ~InternalNode() {
		for ( int i = 0 ; i < 8 ; i ++ ) {
			if (child[i] != NULL) {
				delete child[i];
			}
		}
	}
	NodeType getType ( ) { return INTERNAL; };
};

class QEFMixin {
protected:
	unsigned char signs;
public:
	char height; // depth
	float mp[3]; // this is the minimizer point of the QEF
	int index; // vertex index in PLY file
	float ata[6], atb[3], btb; // QEF data
	
	void clearQEF() {
		for ( int i = 0 ; i < 6 ; i ++ )
			ata[i] = 0 ;

		for ( int i = 0 ; i < 3 ; i ++ ) {
			mp[i] = 0;
			atb[i] = 0 ;
		}

		btb = 0 ;
	}
	void setQEF(float ata1[6], float atb1[3], float btb1, float mp1[3] ) {
		for ( int i = 0 ; i < 6 ; i ++ )
			ata[i] = ata1[i] ;

		for ( int i = 0 ; i < 3 ; i ++ ) {
			mp[i] = mp1[i] ;
			atb[i] = atb1[i] ;
		}
		btb = btb1 ;
	}
	int getSign ( int index ) { return (( signs >> index ) & 1 ); };
};

class LeafNode : public OctreeNode, public QEFMixin {

public:
	virtual ~LeafNode() {};
	
	// Construction
	LeafNode( int ht, unsigned char sg, float coord[3] )  {
		height = ht ;
		signs = sg ;
		clearQEF();
		index = -1 ;
	};

	// Construction by QEF
	// each edge of the cube can have an intersection point
	// so we can give up to 12 intersection points with 12 normal-vectors
	// specify number of intersections in numint
	//
	// st is the minimum bounding-box point
	// st + (1,1,1)*len is the maximum bounding-box point
	LeafNode( int ht, unsigned char sg, int st[3], int len, int numint, float inters[12][3], float norms[12][3] ) {
		height = ht;
		signs = sg;
		index = -1;
		clearQEF();
		
		float pt[3] ={0,0,0} ;
		if ( numint > 0 ) {
			for ( int i = 0 ; i < numint ; i ++ ) {
				float* norm = norms[i] ;
				float* p = inters[i] ;
				// printf("Norm: %f, %f, %f Pts: %f, %f, %f\n", norm[0], norm[1], norm[2], p[0], p[1], p[2] ) ;

				// QEF
				ata[ 0 ] += (float) ( norm[ 0 ] * norm[ 0 ] );
				ata[ 1 ] += (float) ( norm[ 0 ] * norm[ 1 ] );
				ata[ 2 ] += (float) ( norm[ 0 ] * norm[ 2 ] );
				ata[ 3 ] += (float) ( norm[ 1 ] * norm[ 1 ] );
				ata[ 4 ] += (float) ( norm[ 1 ] * norm[ 2 ] );
				ata[ 5 ] += (float) ( norm[ 2 ] * norm[ 2 ] );
				double pn = p[0] * norm[0] + p[1] * norm[1] + p[2] * norm[2] ;
				atb[ 0 ] += (float) ( norm[ 0 ] * pn ) ;
				atb[ 1 ] += (float) ( norm[ 1 ] * pn ) ;
				atb[ 2 ] += (float) ( norm[ 2 ] * pn ) ;
				btb += (float) pn * (float) pn ;
				// Minimizer
				pt[0] += p[0] ;
				pt[1] += p[1] ;
				pt[2] += p[2] ;
			}
			// we minimize towards the average of all intersection points
			pt[0] /= numint ;
			pt[1] /= numint ;
			pt[2] /= numint ;
			// Solve
			float mat[10] ;
			BoundingBoxf * box = new BoundingBoxf();
			box->begin.x = (float) st[0] ;
			box->begin.y = (float) st[1] ;
			box->begin.z = (float) st[2] ;
			box->end.x = (float) st[0] + len ;
			box->end.y = (float) st[1] + len ;
			box->end.z = (float) st[2] + len ;
			
			// eigen.hpp
			// calculate minimizer point, and return error
			// QEF: ata, atb, btb
			// pt is the average of the intersection points
			// mp is the result
			// box is a bounding-box for this node
			// mat is storage for calcPoint() ?
			// float error =  // never used
			calcPoint( ata, atb, btb, pt, mp, box, mat ) ;

#ifdef CLAMP // Clamp all minimizers to be inside the cell
			if ( mp[0] < st[0] || mp[1] < st[1] || mp[2] < st[2] || // mp is outside bounding-box min-pt
				mp[0] > st[0] + len || mp[1] > st[1] + len || mp[2] > st[2] + len ) // mp is outside bounding-box max-pt
			{
				mp[0] = pt[0] ; // reject mp by calcPoint, instead clamp solution to the mass-center
				mp[1] = pt[1] ;
				mp[2] = pt[2] ;
			}
#endif
		}
		else {
			printf("Number of edge intersections in this leaf cell is zero!\n") ;
			mp[0] = st[0] + len / 2;
			mp[1] = st[1] + len / 2;
			mp[2] = st[2] + len / 2;
		}
	};

	NodeType getType ( ) { return LEAF ; };
};

// leaf, but not at max depth
// created by merging child-nodes
class PseudoLeafNode : public OctreeNode, public QEFMixin {
public:
	OctreeNode * child[8] ; // Children 
	
	virtual ~PseudoLeafNode() {
		for ( int i = 0 ; i < 8 ; i ++ ) {
			if (child[i] != NULL)
				delete child[i];
		}
	}
	
	// Construction, without QEF
	PseudoLeafNode ( int ht, unsigned char sg, float coord[3] )  {
		height = ht;
		signs = sg;
		clearQEF();
		for ( int i = 0 ; i < 3 ; i ++ ) 
			mp[i] = coord[i] ;

		for ( int i = 0 ; i < 8 ; i ++ ) 
			child[i] = NULL ;

		index = -1 ;
	};

	// construction with QEF
	PseudoLeafNode ( int ht, unsigned char sg, float ata1[6], float atb1[3], float btb1, float mp1[3] )  {
		height = ht ;
		signs = sg ;
		setQEF(ata1, atb1, btb1, mp1);
		for ( int i = 0 ; i < 8 ; i ++ )
			child[i] = NULL ;

		index = -1 ;
	};
	NodeType getType ( ) { return PSEUDOLEAF ; };
};





/**
 * Class for building and processing an octree
 */
class Octree {
public:
	OctreeNode* root ;
	int dimen; 	   // Length of grid
	int maxDepth;
	//int hasQEF;    // used in simplify(), can be removed
	int faceVerts, edgeVerts;
	int actualTris ; // number of triangles produced by cellProcContour()
	int founds, news ;
public:
	Octree();
	~Octree();
	void simplify( float thresh );
	void genContour( char* fname );
	//void genContourNoInter ( char* fname ) ; // not called from main() ?
	//void genContourNoInter2 ( char* fname ) ;

	void countNodes(int result[3]);
	void countNodes( OctreeNode* node, int result[3] );
private:
	float simplify_threshold;
	OctreeNode* simplify( OctreeNode* node, int st[3], int len, float thresh ) ;

	//void readSOG ( char* fname ) ; // read SOG file
	//OctreeNode* readSOG ( FILE* fin, int st[3], int len, int ht, float origin[3], float range ) ;
	//void readDCF ( char* fname ) ; // read DCF file
	//OctreeNode* readDCF ( FILE* fin, int st[3], int len, int ht ) ;

// Contouring
	void generateVertexIndex( OctreeNode* node, int& offset, FILE* fout ) ; // not used by NoInter2-functions?

	void cellProcContour ( OctreeNode* node, FILE* fout ) ;
	void faceProcContour ( OctreeNode* node[2], int dir, FILE* fout ) ;
	void edgeProcContour ( OctreeNode* node[4], int dir, FILE* fout ) ;
	void processEdgeWrite ( OctreeNode* node[4], int dir, FILE* fout ) ;
	
	// for cuonting number of vertices/faces, prior to actual contouring
	void cellProcCount ( OctreeNode* node, int& nverts, int& nfaces ) ;
	void faceProcCount ( OctreeNode* node[2], int dir, int& nverts, int& nfaces ) ;
	void edgeProcCount ( OctreeNode* node[4], int dir, int& nverts, int& nfaces ) ;
	void processEdgeCount ( OctreeNode* node[4], int dir, int& nverts, int& nfaces ) ;

/* not used !?
	void cellProcContourNoInter( OctreeNode* node, int st[3], int len, HashMap2* hash, TriangleList* list, int& numTris ) ;
	void faceProcContourNoInter( OctreeNode* node[2], int st[3], int len, int dir, HashMap2* hash, TriangleList* list, int& numTris ) ;
	void edgeProcContourNoInter( OctreeNode* node[4], int st[3], int len, int dir, HashMap2* hash, TriangleList* list, int& numTris ) ;
	void processEdgeNoInter( OctreeNode* node[4], int st[3], int len, int dir, HashMap2* hash, TriangleList* list, int& numTris ) ;
*/

/*
	void cellProcContourNoInter2( OctreeNode* node, int st[3], int len, HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts ) ;
	void faceProcContourNoInter2( OctreeNode* node[2], int st[3], int len, int dir, HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts ) ;
	void edgeProcContourNoInter2( OctreeNode* node[4], int st[3], int len, int dir, HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts ) ;
	void processEdgeNoInter2( OctreeNode* node[4], int st[3], int len, int dir, HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts ) ;
*/

	/**
	 *  Non-intersecting test and tesselation
	 */
	
	//int testFace( int st[3], int len, int dir, float v1[3], float v2[3] ) ;
	//int testEdge( int st[3], int len, int dir, OctreeNode* node[4], float v[4][3] ) ;
	
	// not called?
	// void makeFaceVertex( int st[3], int len, int dir, OctreeNode* node1, OctreeNode* node2, float v[3] ) ;
	
	//void makeEdgeVertex( int st[3], int len, int dir, OctreeNode* node[4], float mp[4][3], float v[3] ) ;
};


#endif
