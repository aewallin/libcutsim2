/*
  Copyright (C) 2011 Tao Ju

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
#include "octree.hpp"

#include <math.h>
#include <iostream>

//#include <boost/program_options.hpp>
//namespace po = boost::program_options;

//#define ALLOW_INTERSECTION

/*	Parameters
 *	argv[1]:	name of input file (.dcf format)
 *	argv[2]:	name of output file (.ply format)
 *	argv[3]:	(OPTIONAL) name of secondary output file (.ply format) 
 *              when using dual contouring, storing self-intersecting triangles.
*/

int main( int args, char* argv[] ) {
	Octree* tree = new Octree();
	int counts[3];
	tree->countNodes(counts);
	std::cout << " Internal " << counts[0] << "\tPseudo " << counts[1] << "\tLeaf " << counts[2] << "\n";

	delete tree;
}

