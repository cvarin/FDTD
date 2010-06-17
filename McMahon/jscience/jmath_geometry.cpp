/*
  Copyright (C) 2006 - 2008  Jeffrey M. McMahon

  This file is part of JSCIENCE.

  JSCIENCE is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  JSCIENCE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with JSCIENCE.  If not, see <http://www.gnu.org/licenses/>.
*/


//************************************************************************
//------------------------------------------------------------------------
//
// NAME: JMATH_GEOMETRY
//
// VERSION: Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_geometry.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 11/20/08
//
//	DESC:
//
//	NOTES: 	i. 
//
//
//	LIST OF SUBROUTINES:
//
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************


//************************************************************************
//			INCLUDE FILES
//************************************************************************
#include "jmath_geometry.h"


//************************************************************************
//			CUSTOM STRUCTURES
//************************************************************************

//************************************************************************
//			GLOBAL VARIABLES
//************************************************************************

int gauss_npts_2d[4], gauss_npts_3d[9];
double gauss_coord_1d[20][20], gauss_weight_1d[20][20], gauss_coord_2d[4][7][3], gauss_weight_2d[4][7], gauss_coord_3d[9][45][4], gauss_weight_3d[9][45];
int igauss_generated = 0;

//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
double get_area(TRIANGLE tri)
{

  double area;
  double lside1, lside2, lside3;

  // 1 - 2: 1
  lside1 = sqrt((tri.v[1].x - tri.v[0].x)*(tri.v[1].x - tri.v[0].x) + (tri.v[1].y - tri.v[0].y)*(tri.v[1].y - tri.v[0].y) + (tri.v[1].z - tri.v[0].z)*(tri.v[1].z - tri.v[0].z));

  // 2 - 3: 2
  lside2 = sqrt((tri.v[2].x - tri.v[1].x)*(tri.v[2].x - tri.v[1].x) + (tri.v[2].y - tri.v[1].y)*(tri.v[2].y - tri.v[1].y) + (tri.v[2].z - tri.v[1].z)*(tri.v[2].z - tri.v[1].z));

  // 3 - 1: 3
  lside3 = sqrt((tri.v[0].x - tri.v[2].x)*(tri.v[0].x - tri.v[2].x) + (tri.v[0].y - tri.v[2].y)*(tri.v[0].y - tri.v[2].y) + (tri.v[0].z - tri.v[2].z)*(tri.v[0].z - tri.v[2].z));

  area = 0.25*sqrt(4*lside1*lside1*lside2*lside2 - (lside1*lside1 + lside2*lside2 - lside3*lside3)*(lside1*lside1 + lside2*lside2 - lside3*lside3));

  return area;
}



//========================================================================
//========================================================================
//
//	NAME:	VECTOR_DOUBLE get_tri_nhat(TRIANGLE tri_info)
//	DESC:	Get the outward pointing unit normal vector for a triangle
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. the normal vector points outwards from (0, 0, 0)
//
//
//========================================================================
//========================================================================
VECTOR_DOUBLE get_tri_nhat(TRIANGLE tri_info)
{
  VECTOR_DOUBLE nhat;

  //=========================================================
  // GET THE TWO VECTORS DEFINING THE TRIANGLE
  //=========================================================
  // vectors defining the triangle
  // i.) two vectors define the triangle surface, one
  // pointing from node 1 to node 2 and another pointing
  // from node 1 to node 3
  VECTOR_DOUBLE vec_node12, vec_node13, tri_mid, test_pt;

  // let us now get the vector that points from node 1 -> node 2
  vec_node12.x = tri_info.v[1].x - tri_info.v[0].x;
  vec_node12.y = tri_info.v[1].y - tri_info.v[0].y;
  vec_node12.z = tri_info.v[1].z - tri_info.v[0].z;  
  
  // let us now get the vector that points from node 1 -> node 3
  vec_node13.x = tri_info.v[2].x - tri_info.v[0].x;
  vec_node13.y = tri_info.v[2].y - tri_info.v[0].y;
  vec_node13.z = tri_info.v[2].z - tri_info.v[0].z;  

  //---------------------------------------------------------
  // GET THE NON-UNITY NORMAL VECTOR
  //---------------------------------------------------------
  // i.) get the normal vector by taking the cross product of
  // vec_node12 and vec_node13
  nhat = cross_product(vec_node12, vec_node13);


  //---------------------------------------------------------
  // NORMALIZE THE NORMAL VECTOR
  //---------------------------------------------------------
  double length = sqrt(nhat.x*nhat.x + nhat.y*nhat.y + nhat.z*nhat.z);

  nhat.x /= length;
  nhat.y /= length;
  nhat.z /= length;

  //---------------------------------------------------------
  // DIRECTIONALIZE THE NORMAL VECTOR OUTWARDS
  //---------------------------------------------------------
  // now make sure that the normal vector points outwards
  // (assuming the surfaces form a sphere around the origin (0, 0, 0))

  //double tri_xmid, tri_ymid, tri_zmid, dist_tri, testx, testy, testz, dist_point;

  // get the centroid of the triangle
  tri_mid = (tri_info.v[0] + tri_info.v[1] + tri_info.v[2])/3.0;
//  tri_mid.x = (tri_info.v[0].x + tri_info.v[1].x + tri_info.v[2].x)/3.0;
//  tri_mid.y = (tri_info.v[0].y + tri_info.v[1].y + tri_info.v[2].y)/3.0;
//  tri_mid.z = (tri_info.v[0].z + tri_info.v[1].z + tri_info.v[2].z)/3.0;

  // get the distance of the centroid from the origin
  //dist_tri = sqrt((tri_xmid)*(tri_xmid) + (tri_ymid)*(tri_ymid) + (tri_zmid)*(tri_zmid));

  // get a point 1 unit off the surface of the triangle 
  //testx = tri_xmid + nhat.x;
  //testy = tri_ymid + nhat.y;  
  //testz = tri_zmid + nhat.z; 

  //dist_point = sqrt((testx)*(testx) + (testy)*(testy) + (testz)*(testz));

  test_pt = tri_mid + nhat;

  // if the normal vector points inwards reverse it
  if(get_length(test_pt) < get_length(tri_mid))
  {
    nhat *= -1.0;
    //nhat.x = 0.0 - nhat.x;
    //nhat.y = 0.0 - nhat.y;
    //nhat.z = 0.0 - nhat.z;
  }

  return nhat;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
double get_tri_circumcircle_radius(TRIANGLE tri)
{
  double r;
  VECTOR_DOUBLE vec12, vec23, vec31;

  vec12 = tri.v[1] - tri.v[2];
  vec23 = tri.v[2] - tri.v[3];
  vec31 = tri.v[3] - tri.v[1];

  r = get_length(vec12)*get_length(vec23)*get_length(vec31)/(2.0*get_length(cross_product(vec12, vec23)));

  return r;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
void get_tri_circumcircle_center(TRIANGLE tri, VECTOR_DOUBLE &center)
{
  VECTOR_DOUBLE vec12, vec13, vec21, vec23, vec31, vec32;
  double denom, talpha, tbeta, tgamma;

  vec12 = tri.v[1] - tri.v[2];
  vec13 = tri.v[1] - tri.v[3];
  vec21 = tri.v[2] - tri.v[1];
  vec23 = tri.v[2] - tri.v[3];
  vec31 = tri.v[3] - tri.v[1];
  vec32 = tri.v[3] - tri.v[2];

  // GET THE DENOMINATOR OF ALL PARAMETERS
  denom = 2.0*get_length(cross_product(vec12, vec23))*get_length(cross_product(vec12, vec23));

  talpha = get_length(vec23)*get_length(vec23)*dot_product(vec12, vec13)/denom;

  tbeta = get_length(vec13)*get_length(vec13)*dot_product(vec21, vec23)/denom;

  tgamma = get_length(vec12)*get_length(vec12)*dot_product(vec31, vec32)/denom;

  center = tri.v[1]*talpha + tri.v[2]*tbeta + tri.v[3]*tgamma;

/*
//&&&&&&&&&&&&&&&&

  VECTOR_DOUBLE vec_12, vec_13, vec_23;

  vec_12 = tri.position[1] - tri.position[2];
  vec_13 = tri.position[1] - tri.position[3];
  vec_23 = tri.position[2] - tri.position[3];

  double a = get_length(vec_12);
  double b = get_length(vec_13);
  double c = get_length(vec_23);

  double denom2 = 2.0*(a*a*b*b + a*a*c*c + b*b*c*c) - (pow(a,4)+pow(b,4)+pow(c,4));

  double acoeff, bcoeff, ccoeff;
  acoeff = a*a*(b*b + c*c - a*a);
  bcoeff = b*b*(a*a + c*c - b*b);
  ccoeff = c*c*(a*a + b*b - c*c);

  VECTOR_DOUBLE center2 = tri.position[3]*acoeff + tri.position[2]*bcoeff + tri.position[1]*ccoeff;
  center2 /= denom2;

  cout << "x1: " << center.y << endl;
  cout << "x2: " << center2.y << endl;
  int rgh;
  cin >> rgh;

//&&&&&&&&&&&&&&&&&&&&&&&&
*/

  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//		iii. as long as the nodes are properly labeled (right-hand volume
//		rule) the volume should be positive.  It is up to the user to 
//		check the volume themself where appropriate
//
//
//========================================================================
//========================================================================
double get_volume(TETRAHEDRON tet)
{
  return ( (tet.v[2].x*tet.v[1].y*tet.v[0].z - tet.v[0].z*tet.v[3].x*tet.v[1].y - tet.v[1].x*tet.v[2].y*tet.v[0].z + tet.v[3].x*tet.v[2].y*tet.v[0].z + tet.v[1].x*tet.v[3].y*tet.v[0].z - tet.v[2].x*tet.v[3].y*tet.v[0].z - tet.v[2].x*tet.v[0].y*tet.v[1].z + tet.v[3].x*tet.v[0].y*tet.v[1].z + tet.v[0].x*tet.v[2].y*tet.v[1].z - tet.v[3].x*tet.v[2].y*tet.v[1].z - tet.v[0].x*tet.v[3].y*tet.v[1].z + tet.v[2].x*tet.v[3].y*tet.v[1].z + tet.v[1].x*tet.v[0].y*tet.v[2].z - tet.v[3].x*tet.v[0].y*tet.v[2].z - tet.v[0].x*tet.v[1].y*tet.v[2].z + tet.v[3].x*tet.v[1].y*tet.v[2].z + tet.v[0].x*tet.v[3].y*tet.v[2].z - tet.v[1].x*tet.v[3].y*tet.v[2].z - tet.v[1].x*tet.v[0].y*tet.v[3].z + tet.v[2].x*tet.v[0].y*tet.v[3].z + tet.v[0].x*tet.v[1].y*tet.v[3].z - tet.v[2].x*tet.v[1].y*tet.v[3].z - tet.v[0].x*tet.v[2].y*tet.v[3].z + tet.v[1].x*tet.v[2].y*tet.v[3].z)/6.0 );
}


/*
!! ALTERNATIVE FORMULATION !!
//========================================================================
//========================================================================
//
//	NAME:	double tetrahedron_volume_calc(double x[], double y[], double z[])
//	DESC:	
//
//	NOTES: 	i. we are going to use the formula from wikipedia:
//
//			if the four vertices are defined by a, b, c, d then we can use
//			(we will assume a, b, c, and d correlate to 1, 2, 3, and 4):
//		
//			V = (1/6)*|(a-d) . ((b-d) X (c-d))|
//
//		ii. the volume we obtain here will never be negative because of the absolute 
//		value sign (in FEM we have to be careful though because we want to define
//		things such that the tetrahedra have positive volume without the abs sign)
//		
//
//
//
//
//========================================================================
//========================================================================
double tetrahedron_volume_calc(tet_t tet)
{

  double volume = 0.0;

  VECTOR_DOUBLE amd, bmd, cmd, cp;
  double dp;

  VECTOR_DOUBLE *tet_nodes;
  tet_nodes=new VECTOR_DOUBLE[5];
  int i;
  for(i=0;i<4;++i)
  {
    tet_nodes[i+1].x=tet.v[i].x;
    tet_nodes[i+1].y=tet.v[i].y;
    tet_nodes[i+1].z=tet.v[i].z;
  }

  // ENTER THE LOOP TO CALCULATE THE VOLUME
  // i. this should only be run a maximum of twice (i.e. twice if we obtain
  // a negative volume on the first pass through)
  do
  {

    // RESET THE VOLUME
    // i. this is in case we calculated a negative volume on the first pass through
    volume = 0.0;
 
    // GET IMPORTANT VECTORS
    // i. these are all the vectors in the notes above (d-a, d-b, d-c)
    amd = tet_nodes[1] - tet_nodes[4];
    bmd = tet_nodes[2] - tet_nodes[4];
    cmd = tet_nodes[3] - tet_nodes[4];

    // GET THE CROSS PRODUCT BETWEEN d-b & d-c
    cp = cross_product(bmd, cmd);

    // GET THE DOT PRODUCT BETWEEN d-a & (d-b X d-c)
    dp = dot_product(amd, cp);

    
    // GET THE VOLUME
    // i. we essentially have this, we just have to multiply by the appropriate factor
    // ii. !!! according to the wikipedia formula we should take the absolute value of the dot
    // product, see the notes above for why I am not doing this
    // iii. !!! without the fabs I think we consistently get stuck in an infinite loop here 
    // (at least by only swapping nodes 1 & 2 below)
    volume = fabs(dp)/6.0;
    //volume = dp/6.0;


    // SWAP NODES 1 & 2
    // i. this is in case we have calculated a negative volume and we take a second pass through.  
    // ii. if we swap two nodes then we will change the sign (!!! I think, if not we could get 
    // stuck in an infinite loop here) of the volume calculated
    vector_swap(tet_nodes[1], tet_nodes[2]);


  // if we calculate a negative volume then restart with the swapped nodes
  } while (volume < 0.0);


  delete [] tet_nodes;

  return volume;

}
*/


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//
//========================================================================
//========================================================================
void get_tet_info(TETRAHEDRON tet, double ae[4], double be[4], double ce[4], double de[4], double &volume_e)
{

  ae[0] = tet.v[1].x*tet.v[2].y*tet.v[3].z - tet.v[1].x*tet.v[3].y*tet.v[2].z + tet.v[2].x*tet.v[3].y*tet.v[1].z - tet.v[2].x*tet.v[1].y*tet.v[3].z + tet.v[3].x*tet.v[1].y*tet.v[2].z - tet.v[3].x*tet.v[2].y*tet.v[1].z;
  ae[1] = 0.0 - tet.v[2].x*tet.v[3].y*tet.v[0].z + tet.v[2].x*tet.v[0].y*tet.v[3].z - tet.v[3].x*tet.v[0].y*tet.v[2].z + tet.v[3].x*tet.v[2].y*tet.v[0].z - tet.v[0].x*tet.v[2].y*tet.v[3].z + tet.v[0].x*tet.v[3].y*tet.v[2].z;
  ae[2] = tet.v[3].x*tet.v[0].y*tet.v[1].z - tet.v[3].x*tet.v[1].y*tet.v[0].z + tet.v[0].x*tet.v[1].y*tet.v[3].z - tet.v[0].x*tet.v[3].y*tet.v[1].z + tet.v[1].x*tet.v[3].y*tet.v[0].z - tet.v[1].x*tet.v[0].y*tet.v[3].z;
  ae[3] = 0.0 - tet.v[0].x*tet.v[1].y*tet.v[2].z + tet.v[0].x*tet.v[2].y*tet.v[1].z - tet.v[1].x*tet.v[2].y*tet.v[0].z + tet.v[1].x*tet.v[0].y*tet.v[2].z - tet.v[2].x*tet.v[0].y*tet.v[1].z + tet.v[2].x*tet.v[1].y*tet.v[0].z;

  be[0] = tet.v[2].y*tet.v[1].z - tet.v[3].y*tet.v[1].z - tet.v[1].y*tet.v[2].z + tet.v[3].y*tet.v[2].z + tet.v[1].y*tet.v[3].z - tet.v[2].y*tet.v[3].z;
  be[1] = tet.v[3].y*tet.v[0].z - tet.v[2].y*tet.v[0].z + tet.v[0].y*tet.v[2].z - tet.v[3].y*tet.v[2].z - tet.v[0].y*tet.v[3].z + tet.v[2].y*tet.v[3].z;
  be[2] = tet.v[1].y*tet.v[0].z - tet.v[3].y*tet.v[0].z - tet.v[0].y*tet.v[1].z + tet.v[3].y*tet.v[1].z + tet.v[0].y*tet.v[3].z - tet.v[1].y*tet.v[3].z;
  be[3] = tet.v[2].y*tet.v[0].z - tet.v[1].y*tet.v[0].z + tet.v[0].y*tet.v[1].z - tet.v[2].y*tet.v[1].z - tet.v[0].y*tet.v[2].z + tet.v[1].y*tet.v[2].z;

  ce[0] = tet.v[3].x*tet.v[1].z - tet.v[2].x*tet.v[1].z + tet.v[1].x*tet.v[2].z - tet.v[3].x*tet.v[2].z - tet.v[1].x*tet.v[3].z + tet.v[2].x*tet.v[3].z;
  ce[1] = tet.v[2].x*tet.v[0].z - tet.v[3].x*tet.v[0].z - tet.v[0].x*tet.v[2].z + tet.v[3].x*tet.v[2].z + tet.v[0].x*tet.v[3].z - tet.v[2].x*tet.v[3].z;
  ce[2] = tet.v[3].x*tet.v[0].z - tet.v[1].x*tet.v[0].z + tet.v[0].x*tet.v[1].z - tet.v[3].x*tet.v[1].z - tet.v[0].x*tet.v[3].z + tet.v[1].x*tet.v[3].z;
  ce[3] = tet.v[1].x*tet.v[0].z - tet.v[2].x*tet.v[0].z - tet.v[0].x*tet.v[1].z + tet.v[2].x*tet.v[1].z + tet.v[0].x*tet.v[2].z - tet.v[1].x*tet.v[2].z;

  de[0] = tet.v[2].x*tet.v[1].y - tet.v[3].x*tet.v[1].y - tet.v[1].x*tet.v[2].y + tet.v[3].x*tet.v[2].y + tet.v[1].x*tet.v[3].y - tet.v[2].x*tet.v[3].y; 
  de[1] = 0.0 - tet.v[2].x*tet.v[0].y + tet.v[3].x*tet.v[0].y + tet.v[0].x*tet.v[2].y - tet.v[3].x*tet.v[2].y - tet.v[0].x*tet.v[3].y + tet.v[2].x*tet.v[3].y;
  de[2] = tet.v[1].x*tet.v[0].y - tet.v[3].x*tet.v[0].y - tet.v[0].x*tet.v[1].y + tet.v[3].x*tet.v[1].y + tet.v[0].x*tet.v[3].y - tet.v[1].x*tet.v[3].y;
  de[3] = 0.0 - tet.v[1].x*tet.v[0].y + tet.v[2].x*tet.v[0].y + tet.v[0].x*tet.v[1].y - tet.v[2].x*tet.v[1].y - tet.v[0].x*tet.v[2].y + tet.v[1].x*tet.v[2].y;

  volume_e = get_volume(tet);


  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//
//========================================================================
//========================================================================
void get_tri_info(TRIANGLE tri, double ae[3], double be[3], double ce[3], double &area, VECTOR_DOUBLE &tri_nhat, VECTOR_DOUBLE &vec_u, VECTOR_DOUBLE &vec_v)
{
  int i;
  VECTOR_DOUBLE vec_node12, vec_node13;
  double u_coord[3], v_coord[3], areatemp;

  // GET VECTORS DEFINING TRIANGLE PLANE
  vec_node12 = tri.v[1] - tri.v[0];
  vec_node13 = tri.v[2] - tri.v[0];

  // GET NORMAL VECTOR
  tri_nhat = get_tri_nhat(tri);
  area = get_area(tri);

  // DEFINE NODE 1 AS (0, 0) IN THE NEW COORD SYSTEM
  u_coord[0] = 0.0;
  v_coord[0] = 0.0;


  for(i = 0; i < 2; ++i)
  {
    // THE u COORDINATE POINTS FROM node 1 --> node 2
    vec_u = vec_node12;
    normalize(vec_u);

    // the v coordinate vector will be perpendicular to the u vector and the normal vector
    vec_v = cross_product(tri_nhat, vec_u);
    normalize(vec_v); 

    // now let us get the coordinates of the three nodes defining this face in the 
    // (u, v) coordinate system

    if(i == 1)
    {
      // node 2 is at (dot product of vec_node12 with vec_u, dot product of vec_node12 with vec_v (should be 0))
      u_coord[1] = dot_product(vec_node12, vec_u);
      v_coord[1] = dot_product(vec_node12, vec_v);

      // node 3 is at (dot product of vec_node13 with vec_u, dot product of vec_node13 with vec_v)
      u_coord[2] = dot_product(vec_node13, vec_u);
      v_coord[2] = dot_product(vec_node13, vec_v);    
    }
    else
    {
      // REVERSE THE NODE ORDER
      // node 2 is at (dot product of vec_node12 with vec_u, dot product of vec_node12 with vec_v (should be 0))
      u_coord[2] = dot_product(vec_node12, vec_u);
      v_coord[2] = dot_product(vec_node12, vec_v);

      // node 3 is at (dot product of vec_node13 with vec_u, dot product of vec_node13 with vec_v)
      u_coord[1] = dot_product(vec_node13, vec_u);
      v_coord[1] = dot_product(vec_node13, vec_v);  
    }

    // remembering u is like our x coordinate and v is like our y coordinate we
    // can get the "be" and "ce" values from Jin's book, pg. 96
    be[0] = v_coord[1] - v_coord[2];
    be[1] = v_coord[2] - v_coord[0];
    be[2] = v_coord[0] - v_coord[1];

    ce[0] = u_coord[2] - u_coord[1];
    ce[1] = u_coord[0] - u_coord[2];
    ce[2] = u_coord[1] - u_coord[0];

    areatemp = 0.5*(be[0]*ce[1] - be[1]*ce[0]);
    if(areatemp > 0.0)
    {
      break;
    }
  }

    // GET ae
    ae[0] = u_coord[1]*v_coord[2] - v_coord[1]*u_coord[2];
    ae[1] = u_coord[2]*v_coord[0] - v_coord[2]*u_coord[0];
    ae[2] = u_coord[0]*v_coord[1] - v_coord[0]*u_coord[1];

  // DO A QUICK ERROR CHECK FOR LABELING IN THIS COORDINATE SYSTEM
  //double areatest = 0.5*(be[0]*ce[1] - be[1]*ce[0]);
  if( fabs( (area - areatemp)/area ) > 1.01)
  {
    cout << "Error in (u, v) coordinate system for triangle!" << endl;
    cout << "area: " << area << "   areatest: " << areatemp << endl;
    exit(1);
  }
       
  return;
}

//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//
//========================================================================
//========================================================================
void get_simplex_coord_tri(TRIANGLE tri, VECTOR_DOUBLE r, double &zeta1, double &zeta2, double &zeta3)
{
  // MAKE SURE DESIRED POINT IS IN PLANE OF TRIANGLE
  VECTOR_DOUBLE rp;
  get_projpt_tri(tri, r, rp);

  double aetri[3], betri[3], cetri[3], tri_area, ru, rv;
  VECTOR_DOUBLE tri_nhat, tri_vec_u, tri_vec_v;

  // GET THE INFORMATION ABOUT THE TRIANGLE
  get_tri_info(tri, aetri, betri, cetri, tri_area, tri_nhat, tri_vec_u, tri_vec_v);

  // GET THE COORDINATES OF THE PROJECTED r COORDINATE IN (u, v)
  // i. !!! should this be r + tri.v[0] ???
  ru = dot_product(rp - tri.v[0], tri_vec_u);
  rv = dot_product(rp - tri.v[0], tri_vec_v);

  // GET THE AREA (SIMPLEX) COORDINATES OF r0
  zeta1 = (1.0/(2.0*tri_area))*( aetri[0] + betri[0]*ru + cetri[0]*rv );
  zeta2 = (1.0/(2.0*tri_area))*( aetri[1] + betri[1]*ru + cetri[1]*rv );
  zeta3 = (1.0/(2.0*tri_area))*( aetri[2] + betri[2]*ru + cetri[2]*rv );

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//
//========================================================================
//========================================================================
void get_cart_coord_tri(TRIANGLE tri, double L[3], VECTOR_DOUBLE &r)
{

  double aetri[3], betri[3], cetri[3], tri_area, r1u, r1v, r2u, r2v, r3u, r3v, u, v;
  VECTOR_DOUBLE tri_nhat, tri_vec_u, tri_vec_v, r1, r2, r3;

  // GET THE INFORMATION ABOUT THE TRIANGLE
  get_tri_info(tri, aetri, betri, cetri, tri_area, tri_nhat, tri_vec_u, tri_vec_v);

  // GET THE COORDINATES OF THE VERTICES IN (u, v)
  r1u = dot_product(tri.v[0] - tri.v[0], tri_vec_u);
  r1v = dot_product(tri.v[0] - tri.v[0], tri_vec_v);

  r2u = dot_product(tri.v[1] - tri.v[0], tri_vec_u);
  r2v = dot_product(tri.v[1] - tri.v[0], tri_vec_v);

  r3u = dot_product(tri.v[2] - tri.v[0], tri_vec_u);
  r3v = dot_product(tri.v[2] - tri.v[0], tri_vec_v);

  // GET THE (u, v) COORDINATE OF r
  u = L[0]*r1u + L[1]*r2u + L[2]*r3u;
  v = L[0]*r1v + L[1]*r2v + L[2]*r3v;

  // GET THE (x, y, z) COORDINATE OF r FROM (u, v)
  // i. !!! should this be tri.v[0] ???
  r = tri_vec_u*u + tri_vec_v*v + tri.v[0];

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//
//========================================================================
//========================================================================
void get_projpt_tri(TRIANGLE tri, VECTOR_DOUBLE r, VECTOR_DOUBLE &r0)
{
  VECTOR_DOUBLE tri_nhat = get_tri_nhat(tri);
  double u = dot_product(tri_nhat, (tri.v[2] - r))/dot_product(tri_nhat, tri_nhat);
  r0 = r + tri_nhat*u;

  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//
//========================================================================
//========================================================================
void get_simplex_coord_tet(TETRAHEDRON tet, VECTOR_DOUBLE r, double zeta[4])
{

  double aetet[4], betet[4], cetet[4], detet[4], tet_volume;

  // GET THE INFORMATION ABOUT THE TRIANGLE
  get_tet_info(tet, aetet, betet, cetet, detet, tet_volume);

  // GET THE AREA (SIMPLEX) COORDINATES OF r0
  zeta[0] = (1.0/(6.0*tet_volume))*( aetet[0] + betet[0]*r.x + cetet[0]*r.y + detet[0]*r.z);
  zeta[1] = (1.0/(6.0*tet_volume))*( aetet[1] + betet[1]*r.x + cetet[1]*r.y + detet[1]*r.z);
  zeta[2] = (1.0/(6.0*tet_volume))*( aetet[2] + betet[2]*r.x + cetet[2]*r.y + detet[2]*r.z);
  zeta[3] = (1.0/(6.0*tet_volume))*( aetet[3] + betet[3]*r.x + cetet[3]*r.y + detet[3]*r.z);

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//
//========================================================================
//========================================================================
void get_cart_coord_tet(TETRAHEDRON tet, double zeta[4], VECTOR_DOUBLE &r)
{
  double aetet[4], betet[4], cetet[4], detet[4], tet_volume;

  // GET THE INFORMATION ABOUT THE TRIANGLE
  get_tet_info(tet, aetet, betet, cetet, detet, tet_volume);

  r = tet.v[0]*zeta[0] + tet.v[1]*zeta[1] + tet.v[2]*zeta[2] + tet.v[3]*zeta[3];

  return;
}

//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//		ii. ! 0 based notation
//		iii. !!! change these to a set of external files to read these in
//
//========================================================================
//========================================================================
void generate_gauss()
{

  //=========================================================
  // 1D
  //=========================================================
  gauss_coord_1d[0][0] = 0.0;
  gauss_weight_1d[0][0] = 2.0;

  gauss_coord_1d[1][0] = -0.5773502692;
  gauss_coord_1d[1][1] = 0.5773502692;
  gauss_weight_1d[1][0] = 1.0;
  gauss_weight_1d[1][1] = 1.0;

  gauss_coord_1d[2][0] =  -0.7745966692;
  gauss_coord_1d[2][1] = 0.0;
  gauss_coord_1d[2][2] = 0.7745966692;
  gauss_weight_1d[2][0] = 0.5555555556;
  gauss_weight_1d[2][1] = 0.8888888889;
  gauss_weight_1d[2][2] = 0.5555555556;

  gauss_coord_1d[3][0] = -0.8611363116;
  gauss_coord_1d[3][1] = -0.3399810436;
  gauss_coord_1d[3][2] = 0.3399810436;
  gauss_coord_1d[3][3] = 0.8611363116;
  gauss_weight_1d[3][0] = 0.3478548451;
  gauss_weight_1d[3][1] = 0.6521451549;
  gauss_weight_1d[3][2] = 0.6521451549;
  gauss_weight_1d[3][3] = 0.3478548451;

  gauss_coord_1d[4][0] = -0.9061798459;
  gauss_coord_1d[4][1] = -0.5384693101;
  gauss_coord_1d[4][2] = 0.0;
  gauss_coord_1d[4][3] = 0.5384693101;
  gauss_coord_1d[4][4] = 0.9061798459;
  gauss_weight_1d[4][0] = 0.2369268851;
  gauss_weight_1d[4][1] = 0.4786286705;
  gauss_weight_1d[4][2] = 0.5688888889;
  gauss_weight_1d[4][3] = 0.4786286705;
  gauss_weight_1d[4][4] = 0.2369268851;

  gauss_coord_1d[5][0] = -0.9324695142;
  gauss_coord_1d[5][1] = -0.6612093865;
  gauss_coord_1d[5][2] = -0.2386191861;
  gauss_coord_1d[5][3] = 0.2386191861;
  gauss_coord_1d[5][4] = 0.6612093865;
  gauss_coord_1d[5][5] = 0.9324695142;
  gauss_weight_1d[5][0] = 0.1713244924;
  gauss_weight_1d[5][1] = 0.3607615730;
  gauss_weight_1d[5][2] = 0.4679139346;
  gauss_weight_1d[5][3] = 0.4679139346;
  gauss_weight_1d[5][4] = 0.3607615730;
  gauss_weight_1d[5][5] = 0.1713244924;

  gauss_coord_1d[6][0] = -0.9491079123;
  gauss_coord_1d[6][1] = -0.7415311856;
  gauss_coord_1d[6][2] = -0.4058451514;
  gauss_coord_1d[6][3] = 0.0;
  gauss_coord_1d[6][4] = 0.4058451514;
  gauss_coord_1d[6][5] = 0.7415311856;
  gauss_coord_1d[6][6] = 0.9491079123;
  gauss_weight_1d[6][0] = 0.1294849662;
  gauss_weight_1d[6][1] = 0.2797053915;
  gauss_weight_1d[6][2] = 0.3818300505;
  gauss_weight_1d[6][3] = 0.4179591837;
  gauss_weight_1d[6][4] = 0.3818300505;
  gauss_weight_1d[6][5] = 0.2797053915;
  gauss_weight_1d[6][6] = 0.1294849662;


  gauss_coord_1d[9][0] = -0.973906528517;
  gauss_coord_1d[9][1] = -0.865063366689;
  gauss_coord_1d[9][2] = -0.679409568299;
  gauss_coord_1d[9][3] = -0.433395394129;
  gauss_coord_1d[9][4] = -0.148874338982;
  gauss_coord_1d[9][5] = 0.148874338982;
  gauss_coord_1d[9][6] = 0.433395394129;
  gauss_coord_1d[9][7] = 0.679409568299;
  gauss_coord_1d[9][8] = 0.865063366689;
  gauss_coord_1d[9][9] = 0.973906528517;
  gauss_weight_1d[9][0] = 0.0666713443087;
  gauss_weight_1d[9][1] = 0.149451349151;
  gauss_weight_1d[9][2] = 0.219086362516;
  gauss_weight_1d[9][3] = 0.26926671931;
  gauss_weight_1d[9][4] = 0.295524224715;
  gauss_weight_1d[9][5] = 0.295524224715;
  gauss_weight_1d[9][6] = 0.26926671931;
  gauss_weight_1d[9][7] = 0.219086362516;
  gauss_weight_1d[9][8] = 0.149451349151;
  gauss_weight_1d[9][9] = 0.0666713443087;



  gauss_coord_1d[19][0] = -0.993128599185;
  gauss_coord_1d[19][1] = -0.963971927278;
  gauss_coord_1d[19][2] = -0.912234428251;
  gauss_coord_1d[19][3] = -0.839116971822;
  gauss_coord_1d[19][4] = -0.74633190646;
  gauss_coord_1d[19][5] = -0.636053680727;
  gauss_coord_1d[19][6] = -0.510867001951;
  gauss_coord_1d[19][7] = -0.373706088715;
  gauss_coord_1d[19][8] = -0.227785851142;
  gauss_coord_1d[19][9] = -0.0765265211335;
  gauss_coord_1d[19][10] = 0.0765265211335;
  gauss_coord_1d[19][11] = 0.227785851142;
  gauss_coord_1d[19][12] = 0.373706088715;
  gauss_coord_1d[19][13] = 0.510867001951;
  gauss_coord_1d[19][14] = 0.636053680727;
  gauss_coord_1d[19][15] = 0.74633190646;
  gauss_coord_1d[19][16] = 0.839116971822;
  gauss_coord_1d[19][17] = 0.912234428251;
  gauss_coord_1d[19][18] = 0.963971927278;
  gauss_coord_1d[19][19] = 0.993128599185;
  gauss_weight_1d[19][0] = 0.0176140070678;
  gauss_weight_1d[19][1] = 0.0406014298819;
  gauss_weight_1d[19][2] = 0.0626720482976;
  gauss_weight_1d[19][3] = 0.0832767415506;
  gauss_weight_1d[19][4] = 0.101930119826;
  gauss_weight_1d[19][5] = 0.118194531969;
  gauss_weight_1d[19][6] = 0.131688638458;
  gauss_weight_1d[19][7] = 0.142096109327;
  gauss_weight_1d[19][8] = 0.149172986482;
  gauss_weight_1d[19][9] = 0.15275338714;
  gauss_weight_1d[19][10] = 0.15275338714;
  gauss_weight_1d[19][11] = 0.149172986482;
  gauss_weight_1d[19][12] = 0.142096109327;
  gauss_weight_1d[19][13] = 0.131688638458;
  gauss_weight_1d[19][14] = 0.118194531969;
  gauss_weight_1d[19][15] = 0.101930119826;
  gauss_weight_1d[19][16] = 0.0832767415506;
  gauss_weight_1d[19][17] = 0.0626720482976;
  gauss_weight_1d[19][18] = 0.0406014298819;
  gauss_weight_1d[19][19] = 0.0176140070678;


  //=========================================================
  // 2D TRIANGLES
  //=========================================================

  gauss_npts_2d[0] = 1;
  gauss_coord_2d[0][0][0] = 1.0/3.0;
  gauss_coord_2d[0][0][1] = 1.0/3.0;
  gauss_coord_2d[0][0][2] = 1.0/3.0;
  gauss_weight_2d[0][0] = 1.0;

  gauss_npts_2d[1] = 3;
  gauss_coord_2d[1][0][0] = 0.5;
  gauss_coord_2d[1][1][0] = 0.0;
  gauss_coord_2d[1][2][0] = 0.5;
  gauss_coord_2d[1][0][1] = 0.5;
  gauss_coord_2d[1][1][1] = 0.5;
  gauss_coord_2d[1][2][1] = 0.0;
  gauss_coord_2d[1][0][2] = 0.0;
  gauss_coord_2d[1][1][2] = 0.5;
  gauss_coord_2d[1][2][2] = 0.5;
  gauss_weight_2d[1][0] = 1.0/3.0;
  gauss_weight_2d[1][1] = 1.0/3.0;
  gauss_weight_2d[1][2] = 1.0/3.0;

  gauss_npts_2d[2] = 4;
  gauss_coord_2d[2][0][0] = 0.33333333;
  gauss_coord_2d[2][1][0] = 0.73333333;
  gauss_coord_2d[2][2][0] = 0.13333333;
  gauss_coord_2d[2][3][0] = 0.13333333;
  gauss_coord_2d[2][0][1] = 0.33333333;
  gauss_coord_2d[2][1][1] = 0.13333333;
  gauss_coord_2d[2][2][1] = 0.73333333;
  gauss_coord_2d[2][3][1] = 0.13333333;
  gauss_coord_2d[2][0][2] = 0.33333333;
  gauss_coord_2d[2][1][2] = 0.13333333;
  gauss_coord_2d[2][2][2] = 0.13333333;
  gauss_coord_2d[2][3][2] = 0.73333333;
  gauss_weight_2d[2][0] = -0.5625;
  gauss_weight_2d[2][1] = 0.52083333;
  gauss_weight_2d[2][2] = 0.52083333;
  gauss_weight_2d[2][3] = 0.52083333;

  gauss_npts_2d[3] = 7;
  gauss_coord_2d[3][0][0] = 0.33333333;
  gauss_coord_2d[3][1][0] = 0.5;
  gauss_coord_2d[3][2][0] = 0.0;
  gauss_coord_2d[3][3][0] = 0.5;
  gauss_coord_2d[3][4][0] = 1.0;
  gauss_coord_2d[3][5][0] = 0.0;
  gauss_coord_2d[3][6][0] = 0.0;
  gauss_coord_2d[3][0][1] = 0.33333333;
  gauss_coord_2d[3][1][1] = 0.5;
  gauss_coord_2d[3][2][1] = 0.5;
  gauss_coord_2d[3][3][1] = 0.0;
  gauss_coord_2d[3][4][1] = 0.0;
  gauss_coord_2d[3][5][1] = 1.0;
  gauss_coord_2d[3][6][1] = 0.0;
  gauss_coord_2d[3][0][2] = 0.33333333;
  gauss_coord_2d[3][1][2] = 0.0;
  gauss_coord_2d[3][2][2] = 0.5;
  gauss_coord_2d[3][3][2] = 0.5;
  gauss_coord_2d[3][4][2] = 0.0;
  gauss_coord_2d[3][5][2] = 0.0;
  gauss_coord_2d[3][6][2] = 1.0;
  gauss_weight_2d[3][0] = 0.45;
  gauss_weight_2d[3][1] = 0.13333333;
  gauss_weight_2d[3][2] = 0.13333333;
  gauss_weight_2d[3][3] = 0.13333333;
  gauss_weight_2d[3][4] = 0.05;
  gauss_weight_2d[3][5] = 0.05;
  gauss_weight_2d[3][6] = 0.05;

  //=========================================================
  // 3D TETRAHEDRA
  //=========================================================

  gauss_npts_3d[0] = 1;
  gauss_coord_3d[0][0][0] = 0.25;
  gauss_coord_3d[0][0][1] = 0.25;
  gauss_coord_3d[0][0][2] = 0.25;
  gauss_coord_3d[0][0][3] = 0.25;
  gauss_weight_3d[0][0] = 1.0;

  gauss_npts_3d[1] = 4;
  gauss_coord_3d[1][0][0] = 0.5585410197;
  gauss_coord_3d[1][1][0] = 0.1381966011;
  gauss_coord_3d[1][2][0] = 0.1381966011;
  gauss_coord_3d[1][3][0] = 0.1381966011;
  gauss_coord_3d[1][0][1] = 0.1381966011;
  gauss_coord_3d[1][1][1] = 0.5585410197;
  gauss_coord_3d[1][2][1] = 0.1381966011;
  gauss_coord_3d[1][3][1] = 0.1381966011;
  gauss_coord_3d[1][0][2] = 0.1381966011;
  gauss_coord_3d[1][1][2] = 0.1381966011;
  gauss_coord_3d[1][2][2] = 0.5585410197;
  gauss_coord_3d[1][3][2] = 0.1381966011;
  gauss_coord_3d[1][0][3] = 0.1381966011;
  gauss_coord_3d[1][1][3] = 0.1381966011;
  gauss_coord_3d[1][2][3] = 0.1381966011;
  gauss_coord_3d[1][3][3] = 0.5585410197;
  gauss_weight_3d[1][0] = 0.25;
  gauss_weight_3d[1][1] = 0.25;
  gauss_weight_3d[1][2] = 0.25;
  gauss_weight_3d[1][3] = 0.25;

  gauss_npts_3d[2] = 5;
  gauss_coord_3d[2][0][0] = 0.25;
  gauss_coord_3d[2][1][0] = 0.5;
  gauss_coord_3d[2][2][0] = 0.1666666667;
  gauss_coord_3d[2][3][0] = 0.1666666667;
  gauss_coord_3d[2][4][0] = 0.1666666667;
  gauss_coord_3d[2][0][1] = 0.25;
  gauss_coord_3d[2][1][1] = 0.1666666667;
  gauss_coord_3d[2][2][1] = 0.5;
  gauss_coord_3d[2][3][1] = 0.1666666667;
  gauss_coord_3d[2][4][1] = 0.1666666667;
  gauss_coord_3d[2][0][2] = 0.25;
  gauss_coord_3d[2][1][2] = 0.1666666667;
  gauss_coord_3d[2][2][2] = 0.1666666667;
  gauss_coord_3d[2][3][2] = 0.5;
  gauss_coord_3d[2][4][2] = 0.1666666667;
  gauss_coord_3d[2][0][3] = 0.25;
  gauss_coord_3d[2][1][3] = 0.1666666667;
  gauss_coord_3d[2][2][3] = 0.1666666667;
  gauss_coord_3d[2][3][3] = 0.1666666667;
  gauss_coord_3d[2][4][3] = 0.5;
  gauss_weight_3d[2][0] = -0.8;
  gauss_weight_3d[2][1] = 0.45;
  gauss_weight_3d[2][2] = 0.45;
  gauss_weight_3d[2][3] = 0.45;
  gauss_weight_3d[2][4] = 0.45;



  gauss_npts_3d[4] = 14;
  gauss_coord_3d[4][0][0] = 0.0673422422;
  gauss_coord_3d[4][1][0] = 0.3108859193;
  gauss_coord_3d[4][2][0] = 0.3108859193;
  gauss_coord_3d[4][3][0] = 0.3108859193;
  gauss_coord_3d[4][4][0] = 0.7217942491;
  gauss_coord_3d[4][5][0] = 0.0927352503;
  gauss_coord_3d[4][6][0] = 0.0927352503;
  gauss_coord_3d[4][7][0] = 0.0927352503;
  gauss_coord_3d[4][8][0] = 0.4544962959;
  gauss_coord_3d[4][9][0] = 0.4544962959;
  gauss_coord_3d[4][10][0] = 0.4544962959;
  gauss_coord_3d[4][11][0] = 0.0455037041;
  gauss_coord_3d[4][12][0] = 0.0455037041;
  gauss_coord_3d[4][13][0] = 0.0455037041;

  gauss_coord_3d[4][0][1] = 0.3108859193;
  gauss_coord_3d[4][1][1] = 0.0673422422;
  gauss_coord_3d[4][2][1] = 0.3108859193;
  gauss_coord_3d[4][3][1] = 0.3108859193;
  gauss_coord_3d[4][4][1] = 0.0927352503;
  gauss_coord_3d[4][5][1] = 0.7217942491;
  gauss_coord_3d[4][6][1] = 0.0927352503;
  gauss_coord_3d[4][7][1] = 0.0927352503;
  gauss_coord_3d[4][8][1] = 0.4544962959;
  gauss_coord_3d[4][9][1] = 0.0455037041;
  gauss_coord_3d[4][10][1] = 0.0455037041;
  gauss_coord_3d[4][11][1] = 0.4544962959;
  gauss_coord_3d[4][12][1] = 0.4544962959;
  gauss_coord_3d[4][13][1] = 0.0455037041;

  gauss_coord_3d[4][0][2] = 0.3108859193;
  gauss_coord_3d[4][1][2] = 0.3108859193;
  gauss_coord_3d[4][2][2] = 0.0673422422;
  gauss_coord_3d[4][3][2] = 0.3108859193;
  gauss_coord_3d[4][4][2] = 0.0927352503;
  gauss_coord_3d[4][5][2] = 0.0927352503;
  gauss_coord_3d[4][6][2] = 0.7217942491;
  gauss_coord_3d[4][7][2] = 0.0927352503;
  gauss_coord_3d[4][8][2] = 0.0455037041;
  gauss_coord_3d[4][9][2] =  0.4544962959;
  gauss_coord_3d[4][10][2] = 0.0455037041;
  gauss_coord_3d[4][11][2] = 0.4544962959;
  gauss_coord_3d[4][12][2] = 0.0455037041;
  gauss_coord_3d[4][13][2] = 0.4544962959;

  gauss_coord_3d[4][0][3] = 0.3108859193;
  gauss_coord_3d[4][1][3] = 0.3108859193;
  gauss_coord_3d[4][2][3] = 0.3108859193;
  gauss_coord_3d[4][3][3] = 0.0673422422;
  gauss_coord_3d[4][4][3] = 0.0927352503;
  gauss_coord_3d[4][5][3] = 0.0927352503;
  gauss_coord_3d[4][6][3] = 0.0927352503;
  gauss_coord_3d[4][7][3] = 0.7217942491;
  gauss_coord_3d[4][8][3] = 0.0455037041;
  gauss_coord_3d[4][9][3] = 0.0455037041;
  gauss_coord_3d[4][10][3] = 0.4544962959;
  gauss_coord_3d[4][11][3] = 0.0455037041;
  gauss_coord_3d[4][12][3] = 0.4544962959;
  gauss_coord_3d[4][13][3] = 0.4544962959;

  gauss_weight_3d[4][0] = 0.1126879257;
  gauss_weight_3d[4][1] = 0.1126879257;
  gauss_weight_3d[4][2] = 0.1126879257;
  gauss_weight_3d[4][3] = 0.1126879257;
  gauss_weight_3d[4][4] = 0.0734930431;
  gauss_weight_3d[4][5] = 0.0734930431;
  gauss_weight_3d[4][6] = 0.0734930431;
  gauss_weight_3d[4][7] = 0.0734930431;
  gauss_weight_3d[4][8] = 0.0425460207;
  gauss_weight_3d[4][9] = 0.0425460207;
  gauss_weight_3d[4][10] = 0.0425460207;
  gauss_weight_3d[4][11] = 0.0425460207;
  gauss_weight_3d[4][12] = 0.0425460207;
  gauss_weight_3d[4][13] = 0.0425460207;

  gauss_npts_3d[8] = 45;

  int j, m;

  // LOOP OVER ALL INTEGRATION ORDERS
  //for(j = 0; j <= max_p; ++j)
  //{
    j = 9;

    // FINISH SETTING UP FILENAME
    char filename1_3d_absc[] = "./gauss_legendre_rules/abscissas_3d_p";
    char filename1_3d_weights[] = "./gauss_legendre_rules/weights_3d_p";
    char filename_dat[] = ".dat";
    char filename_order[5];
    char filename_absc[50];
    char filename_weights[50];

    // COPY ALL FILENAMES TOGETHER
    sprintf(filename_order,"%i",j); 
    strcpy(filename_absc, filename1_3d_absc);
    strcat(filename_absc, filename_order);
    strcat(filename_absc, filename_dat);
    strcpy(filename_weights, filename1_3d_weights);
    strcat(filename_weights, filename_order);
    strcat(filename_weights, filename_dat);

    // OPEN FILES
    ifstream file_absc, file_weights;
    //file_absc.open(filename_absc, ios::in);
    //file_weights.open(filename_weights, ios::in);
    file_absc.open("./abscissas_3d_p9.dat", ios::in);
    file_weights.open("./weights_3d_p9.dat", ios::in);

    // READ IN ABSCISSAS AND WEIGHTS
    for(m = 0; m < gauss_npts_3d[j-1]; ++m)
    {
      // READ IN ABSCICCAS
      file_absc >> gauss_coord_3d[j-1][m][0];
      file_absc >> gauss_coord_3d[j-1][m][1];
      file_absc >> gauss_coord_3d[j-1][m][2];
      gauss_coord_3d[j-1][m][3] = 1.0 - gauss_coord_3d[j-1][m][0] - gauss_coord_3d[j-1][m][1] - gauss_coord_3d[j-1][m][2];
      // READ IN WEIGHTS
      file_weights >> gauss_weight_3d[j-1][m];
    } // ++m

    // CLOSE FILE
    file_absc.close();
    file_weights.close();
  //}



  igauss_generated = 1;

  return;
}

//**********************************************
    // FIND THE CENTROID OF ELEMENT
  //  for(i=1;i<=4;++i)
  //  {
   //   vert = element[e].vertex[i];

  //    x[i]=node_positions[(vert-1)*3 + 0];
  //    y[i]=node_positions[(vert-1)*3 + 1];
  //    z[i]=node_positions[(vert-1)*3 + 2];
  //  }

  //  dpos.x=x[4];
  //  dpos.y=y[4];
  //  dpos.z=z[4];
   
  //  ra.x=x[1]-dpos.x;
  //  rb.x=x[2]-dpos.x;
   // rc.x=x[3]-dpos.x;

  //  ra.y=y[1]-dpos.y;
  //  rb.y=y[2]-dpos.y;
  //  rc.y=y[3]-dpos.y;

   // ra.z=z[1]-dpos.z;
   // rb.z=z[2]-dpos.z;
   // rc.z=z[3]-dpos.z;

  //  centroid.x=(ra.x+rb.x+rc.x)/4.0+dpos.x;
  //  centroid.y=(ra.y+rb.y+rc.y)/4.0+dpos.y;
  //  centroid.z=(ra.z+rb.z+rc.z)/4.0+dpos.z;
//**********************************************

bool point_in_tet(tet_t tet, VECTOR_DOUBLE r)
{
  double tol = 1.0e-10;

  int s, t;
  tet_t tettemp;
  double vol_sum = 0.0;

  //=========================================================
  // DIVIDE TET INTO 4 SUBTETS USING r AND SUM VOLUMES
  //=========================================================
  for(s = 0; s < 4; ++s)
  {
    // ADD IN THE TET POINTS
    for(t = 0; t < 4; ++t)
    {
      if(t != s)
      {
        tettemp.v[t] = tet.v[t];
      }
    }

    // ADD IN r
    tettemp.v[s] = r;
 
    // CALCULATE AND ADD IN THE SUBTETRAHEDRON VOLUME
    // i. we don't care about node ordering here
    vol_sum += fabs(get_volume(tettemp));
  }


  //=========================================================
  // COMPARE VOLUMES & RETURN
  //=========================================================
  if(fabs(vol_sum/get_volume(tet) - 1.0) < tol)
  {
    return true;
  }

  return false;
}



//========================================================================
//========================================================================
//
//	NAME:	void assemble_global_system(int mpi_nprocessors)
//	DESC:	forms the preconditioner matrix
//
//	NOTES:	i. this subroutine works on kmatrix_proc
//
//
//========================================================================
//========================================================================
void get_centroid_tet(TETRAHEDRON tet, VECTOR_DOUBLE &r)
{

  double L[4];
  L[0] = L[1] = L[2] = L[3] = 0.25;

  get_cart_coord_tet(tet, L, r);

  return;
}


