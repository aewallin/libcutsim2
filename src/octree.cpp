/*

  Implementations of Octree member functions.

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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <iostream>
#include <cassert>

#include "octree.hpp"
#include "PLYWriter.hpp"

Octree::Octree() {

}

void Octree::simplify( float thresh ) {
	int st[3] = {0,0,0};
	this->root = simplify( this->root, st, this->dimen, thresh ) ;
}

// count number of INTERNAL, LEAF, PSEUDOLEAF nodes
void Octree::countNodes(int result[3]) {
	for (int i=0;i<3;i++)
		result[i]=0;
	countNodes(this->root, result);
}
// recursive node counter
void Octree::countNodes( OctreeNode* node, int result[3] ) {
	switch( node->getType() ) {
		case INTERNAL:
			result[0]++;
			InternalNode* inode;
			inode = (InternalNode*)node;
			for (int i=0;i<8;i++) {
				if (inode->child[i]!=NULL)
					countNodes( inode->child[i], result);
			}
			break;
		case PSEUDOLEAF:
			result[1]++;
			PseudoLeafNode* pnode;
			pnode = (PseudoLeafNode*)node;
			for (int i=0;i<8;i++) {
				if (pnode->child[i]!=NULL)
					countNodes( pnode->child[i], result);
			}
			break;
		case LEAF:
			result[2]++;
			break;
	}
}


// simplify by collapsing nodes where the parent node QEF solution is good enough
OctreeNode* Octree::simplify( OctreeNode* node, int st[3], int len, float thresh ) {
	if ( node == NULL )
		return NULL ;

	if ( node->getType() == INTERNAL ) {
		InternalNode* inode = (InternalNode*)node ;
		int simple = 1;

		// QEF data
		float ata[6] = { 0, 0, 0, 0, 0, 0 };
		float atb[3] = { 0, 0, 0 } ;
		float pt[3] = { 0, 0, 0 } ;
		//float mp[3] = { 0, 0, 0 } ;
		float btb = 0 ;
		int signs[8] = {-1,-1,-1,-1,-1,-1,-1,-1} ;
		int midsign = -1 ;

		// child data
		int nlen = len / 2 ;
		int nst[3] ;
		int ec = 0 ;
		int ht ;

		for ( int i = 0 ; i < 8 ; i ++ ) { // recurse into tree
			nst[0] = st[0] + vertMap[i][0] * nlen ;
			nst[1] = st[1] + vertMap[i][1] * nlen ;
			nst[2] = st[2] + vertMap[i][2] * nlen ;

			inode->child[i] = simplify( inode->child[i], nst, nlen, thresh ) ;
			
			if ( inode->child[i] != NULL ) {
				if ( inode->child[i]->getType() == INTERNAL ) {
					simple = 0 ;
				}
				else if ( inode->child[i]->getType() == LEAF ) { // sum child leaf QEFs
					LeafNode* lnode = (LeafNode *) inode->child[i] ;
					ht = lnode->height ;

					for ( int j = 0 ; j < 6 ; j ++ )
						ata[j] += lnode->ata[j] ; 

					for ( int j = 0 ; j < 3 ; j ++ ) {
						atb[j] += lnode->atb[j] ;
						pt[j] += lnode->mp[j] ;
					}
					if ( lnode->mp[0] == 0 )
						printf("%f %f %f, Height: %d\n", lnode->mp[0], lnode->mp[1], lnode->mp[2], ht) ;

					btb += lnode->btb ;
					ec++ ; // QEF count (?)

					midsign = lnode->getSign( 7 - i ) ;
					signs[i] = lnode->getSign( i ) ;
				}
				else { // pseudoleaf
					assert( inode->child[i]->getType() == PSEUDOLEAF );
					PseudoLeafNode* pnode = (PseudoLeafNode *) inode->child[i];
					ht = pnode->height ;

					for ( int j = 0 ; j < 6 ; j ++ )
						ata[j] += pnode->ata[j] ;

					for ( int j = 0 ; j < 3 ; j ++ ) {
						atb[j] += pnode->atb[j] ;
						pt[j] += pnode->mp[j] ;
					}
					btb += pnode->btb ;
					ec ++ ;

					midsign = pnode->getSign( 7 - i ) ;
					signs[i] = pnode->getSign( i ) ;
				}
			}
		} // all QEFs summed 

		if ( simple ) { // some non-INTERNAL children (?)
			if ( ec == 0 ) { // no QEFs found/summed above ( all childs INTERNAL ?)
				//printf("deleting INTERNAL node because all children INTERNAL\n");
				delete node;
				return NULL;
			}
			else {
				pt[0] = pt[0] / ec; // average of summed points
				pt[1] = pt[1] / ec;
				pt[2] = pt[2] / ec;
				//if ( pt[0] < st[0] || pt[1] < st[1] || pt[2] < st[2] ||
				//	pt[0] > st[0] + len || pt[1] > st[1] + len || pt[2] > st[2] + len )
				//{ // pt outside node cube
					//printf("Out! %f %f %f, Box: (%d %d %d) Len: %d ec: %d\n", pt[0], pt[1], pt[2], st[0],st[1],st[2],len,ec) ;
				//}

				unsigned char sg = 0 ;
				for ( int i = 0 ; i < 8 ; i ++ ) {
					if ( signs[i] == 1 )
						sg |= ( 1 << i ) ;
					else if ( signs[i] == -1 ) {  // Undetermined, use center sign instead
						if ( midsign == 1 )
							sg |= ( 1 << i ) ;
						else if ( midsign == -1 )
							printf("Wrong!");
					}
				}

				// Solve QEF for parent node
				float mat[10];
				BoundingBoxf* box = new BoundingBoxf();
				box->begin.x = (float) st[0] ;
				box->begin.y = (float) st[1] ;
				box->begin.z = (float) st[2] ;
				box->end.x = (float) st[0] + len ;
				box->end.y = (float) st[1] + len ;
				box->end.z = (float) st[2] + len ;
				// pt is the average of child-nodes
				// mp is the new solution point
				float mp[3] = { 0, 0, 0 } ;
				float error = calcPoint( ata, atb, btb, pt, mp, box, mat ) ;
#ifdef CLAMP
				if ( mp[0] < st[0] || mp[1] < st[1] || mp[2] < st[2] || // mp is outside boudning-box
					mp[0] > st[0] + len || mp[1] > st[1] + len || mp[2] > st[2] + len ) {
					mp[0] = pt[0] ;
					mp[1] = pt[1] ;
					mp[2] = pt[2] ;
				}
#endif
				if ( error <= thresh ) { // if parent QEF solution is good enough
					PseudoLeafNode* pnode = new PseudoLeafNode( ht+1, sg, ata, atb, btb, mp ) ;
					delete inode ;
					return pnode ;
				}
				else { // QEF solution not good enough
					return node ;
				}
			}
			
		} else { // simple == 0
			return node ;
		}
	} else { // type != INTERNAL, we can't collapse LEAF or PSEUDOLEAF, so quit
		return node ;
	}
}


/*
void Octree::readDCF( char* fname ) {
	FILE* fin = fopen( fname, "rb" ) ;
	if ( fin == NULL )
		printf("Can not open file %s.\n", fname) ;
	
	// Process header
	char version[10] ;
	fread( version, sizeof( char ), 10, fin ) ;
	if ( strcmp( version, "multisign" ) != 0 ) {
		printf("Wrong DCF version.\n") ;
		exit(0) ;
	}
	
	fread( &(this->dimen), sizeof( int ), 1, fin ) ;
	fread( &(this->dimen), sizeof( int ), 1, fin ) ;
	fread( &(this->dimen), sizeof( int ), 1, fin ) ;
	this->maxDepth = 0 ;
	int temp = 1 ;
	while ( temp < this->dimen ) {
		maxDepth ++ ;
		temp <<= 1 ;
	}
	printf("Dimensions: %d Depth: %d\n", this->dimen, maxDepth ) ;

	// Recursive reader
	int st[3] = {0, 0, 0} ;
	this->root = readDCF( fin, st, dimen, maxDepth ) ;

	// optional octree simplification
	if (simplify_threshold > 0 ) {
		std::cout << "Simplifying with threshold " << simplify_threshold << "\n";
		int nodecount[3];
		countNodes( nodecount );
		std::cout << " Before simplify: Internal " << nodecount[0] << "\tPseudo " << nodecount[1] << "\tLeaf " << nodecount[2] << "\n";
		simplify( simplify_threshold );
		countNodes( nodecount );
		std::cout << "  After simplify: Internal " << nodecount[0] << "\tPseudo " << nodecount[1] << "\tLeaf " << nodecount[2] << "\n";

	}
	printf("Done reading.\n") ;	
	fclose( fin ) ;
}*/

/*
// only InternalNode and LeafNode returned by this function
// st cube corner (x,y,z)
// len cube side length
// ht node depth (root has 0, leaf has maxDepth)
OctreeNode* Octree::readDCF( FILE* fin, int st[3], int len, int height ) {
	OctreeNode* rvalue = NULL ;

	int type ;
	fread( &type, sizeof( int ), 1, fin ) ;  // Get type
	// printf("%d %d (%02d, %02d, %02d) NodeType: %d\n", ht, len, st[0], st[1], st[2], type);

	if ( type == 0 ) { // Internal node
		rvalue = new InternalNode() ;
		int child_len = len / 2 ; // len of child node is half that of parent
		int child_st[3] ;

		for ( int i = 0 ; i < 8 ; i ++ ) {  // create eight child nodes
			child_st[0] = st[0] + vertMap[i][0] * child_len; // child st position
			child_st[1] = st[1] + vertMap[i][1] * child_len;
			child_st[2] = st[2] + vertMap[i][2] * child_len;
			((InternalNode *)rvalue)->child[i] = readDCF( fin, child_st, child_len, height - 1 ) ; // height is one less than parent
		}
		return rvalue ;
	}
	
	else if ( type == 1 ) { // Empty node, 
		short sg ;
		fread( &sg, sizeof( short ), 1, fin ) ; // signs not used??
		assert( rvalue == NULL );
		return rvalue ;
	}
	
	else if ( type == 2 ) { // Leaf node
		short rsg[8] ;
		fread( rsg, sizeof( short ), 8, fin ) ;
		unsigned char sg = 0 ;
		for ( int i = 0 ; i < 8 ; i ++ ) { // set signs
			if ( rsg[i] != 0 )
				sg |= ( 1 << i ) ;
		}

		// intersections and normals
		float inters[12][3], norms[12][3] ;
		int numinters = 0;
		for ( int i = 0 ; i < 12 ; i ++ ) { // potentially there are 12 intersections
			int num ;
			fread( &num, sizeof( int ), 1, fin ) ;
			if ( num > 0 ) {
				for ( int j = 0 ; j < num ; j ++ ) {
					float off ;
					fread( &off, sizeof( float ), 1, fin ) ;
					fread( norms[numinters], sizeof( float ), 3, fin ) ;
					int dir = i / 4 ;
					int base = edgevmap[ i ][ 0 ] ;
					inters[numinters][0] = st[0] + vertMap[base][0] * len ;
					inters[numinters][1] = st[1] + vertMap[base][1] * len ;
					inters[numinters][2] = st[2] + vertMap[base][2] * len ;
					inters[numinters][dir] += off ;
					numinters ++ ;
				}
			}
		}
		
		if ( numinters > 0 )
			rvalue = new LeafNode( height, sg, st, len, numinters, inters, norms ) ;
		else
			rvalue = NULL ;
		
		return rvalue ;
	}
	else {
		printf("Wrong! Type: %d\n", type);
		exit(-1);
	}
}
*/




// original algorithm, may produce intersecting polygons?
void Octree::genContour( char* fname ) {
	int numTris = 0 ;
	int numVertices = 0 ;

	FILE* fout = fopen ( fname, "wb" ) ;
	cellProcCount ( root, numVertices, numTris ) ;
	printf("Vertices counted: %d Triangles counted: %d \n", numVertices, numTris ) ;
	PLYWriter::writeHeader( fout, numVertices, numTris ) ;
	int offset = 0; // start of vertex index

	clock_t start = clock();
	generateVertexIndex( root, offset, fout );  // write vertices to file, populate node->index
	printf("Wrote %d vertices to file\n", offset ) ;

	actualTris = 0 ;
	cellProcContour( this->root, fout ) ; // a single call to root runs algorithm on entire tree
	clock_t finish = clock();
	printf("Time used: %f seconds.\n", (float) (finish - start) / (float) CLOCKS_PER_SEC ) ;
	printf("Actual triangles written: %d\n", actualTris ) ;
	fclose( fout ) ;
}

// this writes out octree vertices to the PLY file
// each vertex gets an index, which is stored in node->index
void Octree::generateVertexIndex( OctreeNode* node, int& offset, FILE* fout ) {
	NodeType type = node->getType() ;

	if ( type == INTERNAL ) { // Internal node, recurse into tree
		InternalNode* inode = ( (InternalNode* ) node ) ;
		for ( int i = 0 ; i < 8 ; i ++ ) {
			if ( inode->child[i] != NULL )
				generateVertexIndex( inode->child[i], offset, fout ) ;
		}
	}
	else if ( type == LEAF ) { // Leaf node
		LeafNode* lnode = ((LeafNode *) node) ;
		PLYWriter::writeVertex( fout, lnode->mp ) ; // write out mp
		lnode->index = offset;
		offset++;
	}
	else if ( type == PSEUDOLEAF ) { // Pseudo leaf node
		PseudoLeafNode* pnode = ((PseudoLeafNode *) node) ;
		PLYWriter::writeVertex( fout, pnode->mp ) ; // write out mp
		pnode->index = offset;
		offset++;
	}
}

// cellProcContour( this->root ) is the entry-point to the entire algorithm
void Octree::cellProcContour( OctreeNode* node, FILE* fout )  {
	if ( node == NULL )
		return ;

	int type = node->getType() ;

	if ( type == INTERNAL ) { // internal node
		InternalNode* inode = (( InternalNode * ) node );
		for ( int i = 0 ; i < 8 ; i ++ ) // 8 Cell calls on children
			cellProcContour( inode->child[ i ], fout );

		for ( int i = 0 ; i < 12 ; i ++ ) {  // 12 face calls, faces between each child node
			int c[ 2 ] = { cellProcFaceMask[ i ][ 0 ], cellProcFaceMask[ i ][ 1 ] };
			OctreeNode* fcd[2];
			fcd[0] = inode->child[ c[0] ] ;
			fcd[1] = inode->child[ c[1] ] ;
			faceProcContour( fcd, cellProcFaceMask[ i ][ 2 ], fout ) ;
		}

		for ( int i = 0 ; i < 6 ; i ++ ) {  // 6 edge calls
			int c[ 4 ] = { cellProcEdgeMask[ i ][ 0 ], cellProcEdgeMask[ i ][ 1 ], cellProcEdgeMask[ i ][ 2 ], cellProcEdgeMask[ i ][ 3 ] };
			OctreeNode* ecd[4] ;
			for ( int j = 0 ; j < 4 ; j ++ )
				ecd[j] = inode->child[ c[j] ] ;

			edgeProcContour( ecd, cellProcEdgeMask[ i ][ 4 ], fout ) ;
		}
	}
}

// node[2] are the two nodes that share a face
// dir comes from cellProcFaceMask[i][2]  where i=0..11
void Octree::faceProcContour ( OctreeNode* node[2], int dir, FILE* fout )  {
	// printf("I am at a face! %d\n", dir ) ;
	if ( ! ( node[0] && node[1] ) ) {
		// printf("I am none.\n") ;
		return ;
	}

	NodeType type[2] = { node[0]->getType(), node[1]->getType() } ;

	if ( type[0] == INTERNAL || type[1] == INTERNAL ) { // both nodes internal
		// 4 face calls
		OctreeNode* fcd[2] ;
		for ( int i = 0 ; i < 4 ; i ++ ) {
			int c[2] = { faceProcFaceMask[ dir ][ i ][ 0 ], faceProcFaceMask[ dir ][ i ][ 1 ] };
			for ( int j = 0 ; j < 2 ; j ++ ) {
				if ( type[j] > 0 )
					fcd[j] = node[j];
				else 
					fcd[j] = ((InternalNode *) node[ j ] )->child[ c[j] ];
			}
			faceProcContour( fcd, faceProcFaceMask[ dir ][ i ][ 2 ], fout ) ;
		}

		// 4 edge calls
		int orders[2][4] = {{ 0, 0, 1, 1 }, { 0, 1, 0, 1 }} ;
		OctreeNode* ecd[4] ;
			
		for ( int i = 0 ; i < 4 ; i ++ ) {
			int c[4] = { faceProcEdgeMask[ dir ][ i ][ 1 ], faceProcEdgeMask[ dir ][ i ][ 2 ],
						 faceProcEdgeMask[ dir ][ i ][ 3 ], faceProcEdgeMask[ dir ][ i ][ 4 ] };
			int* order = orders[ faceProcEdgeMask[ dir ][ i ][ 0 ] ] ;

			for ( int j = 0 ; j < 4 ; j ++ ) {
				if ( type[order[j]] > 0 )
					ecd[j] = node[order[j]] ;
				else
					ecd[j] = ( (InternalNode *) node[ order[ j ] ] )->child[ c[j] ] ;
			}
			edgeProcContour( ecd, faceProcEdgeMask[ dir ][ i ][ 5 ], fout ) ;
		}
//		printf("I am done.\n") ;
	}
	else {
//		printf("I don't have any children.\n") ;
	}
}

// a common edge between four nodes in node[4]
// "dir" comes from cellProcEdgeMask
void Octree::edgeProcContour ( OctreeNode* node[4], int dir, FILE* fout ) {
	if ( ! ( node[0] && node[1] && node[2] && node[3] ) )
		return;

	NodeType type[4] = { node[0]->getType(), node[1]->getType(), node[2]->getType(), node[3]->getType() } ;

	if ( type[0] != INTERNAL && type[1] != INTERNAL  && type[2] != INTERNAL && type[3] != INTERNAL ) {
		processEdgeWrite( node, dir, fout ) ; // a face (quad?) is output
	} else {
		// 2 edge calls
		OctreeNode* ecd[4] ;
		for ( int i = 0 ; i < 2 ; i ++ ) {
			int c[ 4 ] = { edgeProcEdgeMask[ dir ][ i ][ 0 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 1 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 2 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 3 ] } ;

			for ( int j = 0 ; j < 4 ; j ++ ) {
				if ( type[j] > 0 )
					ecd[j] = node[j] ;
				else
					ecd[j] = ((InternalNode *) node[j])->child[ c[j] ] ;
			}

			edgeProcContour( ecd, edgeProcEdgeMask[ dir ][ i ][ 4 ], fout ) ;
		}

	}
}

// this writes out a face to the PLY file
// vertices already exist in the file
// so here we write out topology only, i.e. sets of indices that form a face
void Octree::processEdgeWrite ( OctreeNode* node[4], int dir, FILE* fout )  {
	// Get minimal cell
	int  minht = this->maxDepth+1, mini = -1 ;
	int ind[4], sc[4];
	//int flip[4] = {0,0,0,0} ;
	int flip2;
	for ( int i = 0 ; i < 4 ; i ++ ) {
		if ( node[i]->getType() == LEAF ) {
			LeafNode* lnode = ((LeafNode *) node[i]) ;
			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( lnode->height < minht ) {
				minht = lnode->height ;
				mini = i ;
				if ( lnode->getSign(c1) > 0 )
					flip2 = 1 ;
				else
					flip2 = 0 ;
			}
			ind[i] = lnode->index ;

			if ( lnode->getSign( c1 ) == lnode->getSign( c2 ) )
				sc[ i ] = 0 ;
			else
				sc[ i ] = 1 ;

//				if ( lnode->getSign(c1) > 0 )
//				{
//					flip[ i ] = 1 ;
//				}
			// }
		} else {
			assert( node[i]->getType() == PSEUDOLEAF );
			PseudoLeafNode* pnode = ((PseudoLeafNode *) node[i]) ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( pnode->height < minht ) {
				minht = pnode->height;
				mini = i;
				if ( pnode->getSign(c1) > 0 )
					flip2 = 1 ;
				else
					flip2 = 0 ;
			}
			ind[i] = pnode->index ;

			if ( pnode->getSign( c1 ) == pnode->getSign( c2 ) )
				sc[ i ] = 0 ;
			else
				sc[ i ] = 1 ;

//				if ( pnode->getSign(c1) > 0 )
//				{
//					flip[ i ] = 1 ;
//					flip2 = 1 ;
//				}
//			}
		}

	}

	if ( sc[ mini ] == 1 ) { // condition for any triangle output?
		if ( flip2 == 0 ) {
			actualTris ++ ;
			if ( ind[0] == ind[1] ) { // two indices same, so output triangle
				int tind[] = { ind[0], ind[3], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else if ( ind[1] == ind[3] ) {
				int tind[] = { ind[0], ind[1], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else if ( ind[3] == ind[2] ) {
				int tind[] = { ind[0], ind[1], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else if ( ind[2] == ind[0] ) {
				int tind[] = { ind[1], ind[3], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else { // all indices unique, so output a quad by outputting two triangles
				int tind1[] = { ind[0], ind[1], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind1 ) ;
				int tind2[] = { ind[0], ind[3], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind2 ) ;
				actualTris ++ ; // two triangles, so add one here also
			}
		} else {
			actualTris ++ ;
			if ( ind[0] == ind[1] ) {
				int tind[] = { ind[0], ind[2], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else if ( ind[1] == ind[3] ) {
				int tind[] = { ind[0], ind[2], ind[1] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else if ( ind[3] == ind[2] ) {
				int tind[] = { ind[0], ind[3], ind[1] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else if ( ind[2] == ind[0] ) {
				int tind[] = { ind[1], ind[2], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			} else {
				int tind1[] = { ind[0], ind[3], ind[1] } ;
				PLYWriter::writeFace( fout, 3, tind1 ) ;
				int tind2[] = { ind[0], ind[2], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind2 ) ;
				actualTris++; // two triangles, so add one here also
			}
		}
	}
}


// used initially for counting number of vertices
// genContour calls cellProcCount(root)
void Octree::cellProcCount( OctreeNode* node, int& nverts, int& nfaces )  {
	if ( node == NULL )
		return ;

	int type = node->getType() ;

	if (type != INTERNAL)
		nverts ++ ; // leaf or pseudoleaf produce a vertex
	else { // recurse into tree
		InternalNode* inode = (( InternalNode * ) node ) ;
		
		for ( int i = 0 ; i < 8 ; i ++ ) // 8 Cell calls
			cellProcCount( inode->child[ i ], nverts, nfaces ) ;

		OctreeNode* fcd[2] ;
		for ( int i = 0 ; i < 12 ; i ++ ) {  // 12 face calls
			int c[ 2 ] = { cellProcFaceMask[ i ][ 0 ], cellProcFaceMask[ i ][ 1 ] };
			fcd[0] = inode->child[ c[0] ] ;
			fcd[1] = inode->child[ c[1] ] ;
			faceProcCount( fcd, cellProcFaceMask[ i ][ 2 ], nverts, nfaces ) ;
		}
		
		OctreeNode* ecd[4] ;
		for ( int i = 0 ; i < 6 ; i ++ ) { // 6 edge calls
			int c[ 4 ] = { cellProcEdgeMask[ i ][ 0 ], cellProcEdgeMask[ i ][ 1 ], cellProcEdgeMask[ i ][ 2 ], cellProcEdgeMask[ i ][ 3 ] };
			for ( int j = 0 ; j < 4 ; j ++ )
				ecd[j] = inode->child[ c[j] ] ;
			edgeProcCount( ecd, cellProcEdgeMask[ i ][ 4 ], nverts, nfaces ) ;
		}
	}
}

void Octree::faceProcCount ( OctreeNode* node[2], int dir, int& nverts, int& nfaces ) {
	if ( ! ( node[0] && node[1] ) )
		return ;

	int type[2] = { node[0]->getType(), node[1]->getType() } ;

	if ( type[0] == INTERNAL || type[1] == INTERNAL ) {
		// 4 face calls
		OctreeNode* fcd[2] ;
		for ( int i = 0 ; i < 4 ; i ++ )
		{
			int c[2] = { faceProcFaceMask[ dir ][ i ][ 0 ], faceProcFaceMask[ dir ][ i ][ 1 ] };
			for ( int j = 0 ; j < 2 ; j ++ ) {
				if ( type[j] != INTERNAL )
					fcd[j] = node[j] ;
				else
					fcd[j] = ((InternalNode *) node[ j ] )->child[ c[j] ] ;
			}
			faceProcCount( fcd, faceProcFaceMask[ dir ][ i ][ 2 ], nverts, nfaces ) ;
		}

		int orders[2][4] = {{ 0, 0, 1, 1 }, { 0, 1, 0, 1 }} ;
		OctreeNode* ecd[4] ;
			
		for ( int i = 0 ; i < 4 ; i ++ ) {  // 4 edge calls
			int c[4] = { faceProcEdgeMask[ dir ][ i ][ 1 ], faceProcEdgeMask[ dir ][ i ][ 2 ],
						 faceProcEdgeMask[ dir ][ i ][ 3 ], faceProcEdgeMask[ dir ][ i ][ 4 ] };
			int* order = orders[ faceProcEdgeMask[ dir ][ i ][ 0 ] ] ;

			for ( int j = 0 ; j < 4 ; j ++ ) {
				if ( type[order[j]] != INTERNAL )
					ecd[j] = node[order[j]] ;
				else
					ecd[j] = ( (InternalNode *) node[ order[ j ] ] )->child[ c[j] ] ;
			}

			edgeProcCount( ecd, faceProcEdgeMask[ dir ][ i ][ 5 ], nverts, nfaces ) ;
		}
	}
}

void Octree::edgeProcCount ( OctreeNode* node[4], int dir, int& nverts, int& nfaces ) 
{
	if ( ! ( node[0] && node[1] && node[2] && node[3] ) )
		return ;

	int type[4] = { node[0]->getType(), node[1]->getType(), node[2]->getType(), node[3]->getType() } ;
	if ( type[0] != INTERNAL && type[1] != INTERNAL && type[2] != INTERNAL && type[3] != INTERNAL )
		processEdgeCount( node, dir, nverts, nfaces ) ;
	else {
		// 2 edge calls
		OctreeNode* ecd[4] ;
		for ( int i = 0 ; i < 2 ; i ++ ) {
			int c[ 4 ] = { edgeProcEdgeMask[ dir ][ i ][ 0 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 1 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 2 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 3 ] } ;

			for ( int j = 0 ; j < 4 ; j ++ ) {
				if ( type[j] != INTERNAL )
					ecd[j] = node[j] ;
				else
					ecd[j] = ((InternalNode *) node[j])->child[ c[j] ] ;
			}
			edgeProcCount( ecd, edgeProcEdgeMask[ dir ][ i ][ 4 ], nverts, nfaces ) ;
		}
	}
}

void Octree::processEdgeCount ( OctreeNode* node[4], int dir, int& nverts, int& nfaces )  {
	// Get minimal cell
	int i,   minht = maxDepth+1, mini = -1 ;
	//int ind[4]; // set but not used
	int sc[4];
	//flip[4] = {0,0,0,0} ; // set byt not used
	for ( i = 0 ; i < 4 ; i ++ ) {
		if ( node[i]->getType() == 1 ) {
			LeafNode* lnode = ((LeafNode *) node[i]) ;

			if ( lnode->height < minht ) {
				minht = lnode->height ;
				mini = i ;
			}
			//ind[i] = lnode->index ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( lnode->getSign( c1 ) == lnode->getSign( c2 ) ) {
				sc[ i ] = 0 ;
			}
			else {
				sc[ i ] = 1 ;
				//if ( lnode->getSign(c1) > 0 )
				//	flip[ i ] = 1 ;
			}
		}
		else {
			PseudoLeafNode* pnode = ((PseudoLeafNode *) node[i]) ;
			if ( pnode->height < minht ) {
				minht = pnode->height ;
				mini = i ;
			}
			//ind[i] = pnode->index ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( pnode->getSign( c1 ) == pnode->getSign( c2 ) ) {
				sc[ i ] = 0 ;
			}
			else {
				sc[ i ] = 1 ;
				//if ( pnode->getSign(c1) > 0 )
				//	flip[ i ] = 1 ;
			}
		}
	}

	if ( sc[ mini ] == 1 ) {
		nfaces ++ ;
		if ( node[0] != node[1] && node[1] != node[3] && node[3] != node[2] && node[2] != node[0] )
			nfaces ++ ; // quad, so two triangles
	}

}

