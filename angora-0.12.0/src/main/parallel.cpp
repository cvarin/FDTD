/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//Defines functions that handle information exchange with other nodes in the network

#include "headers.h"

#include "parallel.h"

extern Array<double,3> Hx,Hy,Hz;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_SubComm,MPI_CartSubComm;
#endif
extern int rank_behind,rank_front,rank_left,rank_right,rank_below,rank_above;
extern int rank, nodes;
extern int rank_x, rank_y, rank_z;
extern int nodes_x, nodes_y, nodes_z;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;
extern Array<double,2> SendBuf_Hx_y,SendBuf_Hz_y,RecvBuf_Hx_y,RecvBuf_Hz_y,
 SendBuf_Hy_x,SendBuf_Hz_x,RecvBuf_Hy_x,RecvBuf_Hz_x,
 SendBuf_Hx_z,SendBuf_Hy_z,RecvBuf_Hx_z,RecvBuf_Hy_z;
#ifndef MPI_DISABLE
extern MPI_Status Status;
#endif

///** OVERRIDDEN OSTREAM **/
//class MPI_streambuf : public streambuf
//{
//private:
//  int xsputn (char_type* s, streamsize n)
//  {
//  	if (_rank==0) cout.write(s, n);
//    return n;
//  }
//
//  const int _rank;
//};
//
//ostream MPI_cout(&MPIbuf);
///** OVERRIDDEN OSTREAM **/


void init_parallel()
{//assigns the group number, rank number and group size
	//A node only has to know these things:
	//1) Its subgroup index (=0...NumOfGroups-1)
	//2) Its rank in the subgroup, and the total number of nodes in the subgroup

	//Before anything, determine the global number of nodes, and rank in MPI_COMM_WORLD
	//these will not be used globally in the code, so define local
	int rank_global,nodes_global;
#ifndef MPI_DISABLE
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_global);
	MPI_Comm_size(MPI_COMM_WORLD, &nodes_global);
#else
	// no MPI, only one node
	rank_global = 0;
	nodes_global = 1;
#endif

// 	//divide the global communicator into subgroups running independent grids
// 	//Create the subcommunicators
// 	div_t partition = div(rank_global,number_of_runs);
// 	GridIndex = partition.rem;
// 	//split the global communicator according to the subgroup index
// 	MPI_Comm_split(MPI_COMM_WORLD,GridIndex,rank,&MPI_SubComm);
//
// 	//finally, assign the rank number and group size for a single grid (these will be used throughout the code)
// 	MPI_Comm_rank(MPI_SubComm, &rank);
// 	MPI_Comm_size(MPI_SubComm, &nodes);

	// in the following, we abandon the idea of multiple grids
	//we will consider multiple runs, which will happen inside a loop in FDTD.cpp
	//GridIndex will now indicate run index (still called GridIndex, though)

	//define subcommunicator
	//use the common value 0 instead of GridIndex to include all nodes in the same subcommunicator
#ifndef MPI_DISABLE
	MPI_Comm_split(MPI_COMM_WORLD,0,rank,&MPI_SubComm);
#endif

	//finally, assign the rank number and group size for a single grid (these will be used throughout the code)
	//since the common value 0 has been used above, rank and nodes must be equal to rank_global and nodes_global
#ifndef MPI_DISABLE
	MPI_Comm_rank(MPI_SubComm, &rank);
	MPI_Comm_size(MPI_SubComm, &nodes);
#else
	// no MPI, only one node
	rank = 0;
	nodes = 1;
#endif

}

void exchangeH()
//	Handles the H-field exchange before update
{
#ifndef MPI_DISABLE
	enum FieldComponentTag{HxyTag=1,HzyTag=2,HyxTag=3,HzxTag=4, HxzTag=5, HyzTag=6};
	enum LimitTag{LowLimitTag=1,HighLimitTag=2};

	//Tangential H-field values are exchanged in every direction.
	//The tangential H-values that are 1/2 cell away from the grid the boundary in the INWARD direction
	//are sent to the neighboring node. The corresponding tangential H-field values received from
	//the neighboring node are placed 1/2 cell away from the grid boundary in the OUTWARD direction.
	//These values are used to update the E-field components on the grid boundary.

	// exchange Hy with the node behind this one
	SendBuf_Hy_x = Hy(iback,Range(jleft,jright+1),Range(klower,kupper));
	MPI_Sendrecv(SendBuf_Hy_x.data(),SendBuf_Hy_x.size(),MPI_DOUBLE,rank_behind,HyxTag,
		RecvBuf_Hy_x.data(),RecvBuf_Hy_x.size(),MPI_DOUBLE,rank_behind,HyxTag,MPI_CartSubComm,&Status);
	Hy(iback-1,Range(jleft,jright+1),Range(klower,kupper))=RecvBuf_Hy_x;
	// exchange Hz with the node behind this one
	SendBuf_Hz_x = Hz(iback,Range(jleft,jright),Range(klower,kupper+1));
	MPI_Sendrecv(SendBuf_Hz_x.data(),SendBuf_Hz_x.size(),MPI_DOUBLE,rank_behind,HzxTag,
		RecvBuf_Hz_x.data(),RecvBuf_Hz_x.size(),MPI_DOUBLE,rank_behind,HzxTag,MPI_CartSubComm,&Status);
	Hz(iback-1,Range(jleft,jright),Range(klower,kupper+1))=RecvBuf_Hz_x;

	// exchange Hy with the node in front of this one
	SendBuf_Hy_x = Hy(ifront,Range(jleft,jright+1),Range(klower,kupper));
	MPI_Sendrecv(SendBuf_Hy_x.data(),SendBuf_Hy_x.size(),MPI_DOUBLE,rank_front,HyxTag,
		RecvBuf_Hy_x.data(),RecvBuf_Hy_x.size(),MPI_DOUBLE,rank_front,HyxTag,MPI_CartSubComm,&Status);
	Hy(ifront+1,Range(jleft,jright+1),Range(klower,kupper))=RecvBuf_Hy_x;
	// exchange Hz with the node in front of this one
	SendBuf_Hz_x = Hz(ifront,Range(jleft,jright),Range(klower,kupper+1));
	MPI_Sendrecv(SendBuf_Hz_x.data(),SendBuf_Hz_x.size(),MPI_DOUBLE,rank_front,HzxTag,
		RecvBuf_Hz_x.data(),RecvBuf_Hz_x.size(),MPI_DOUBLE,rank_front,HzxTag,MPI_CartSubComm,&Status);
	Hz(ifront+1,Range(jleft,jright),Range(klower,kupper+1))=RecvBuf_Hz_x;

	// exchange Hx with the node on the left of this one
	SendBuf_Hx_y = Hx(Range(iback,ifront+1),jleft,Range(klower,kupper));
	MPI_Sendrecv(SendBuf_Hx_y.data(),SendBuf_Hx_y.size(),MPI_DOUBLE,rank_left,HxyTag,
		RecvBuf_Hx_y.data(),RecvBuf_Hx_y.size(),MPI_DOUBLE,rank_left,HxyTag,MPI_CartSubComm,&Status);
	Hx(Range(iback,ifront+1),jleft-1,Range(klower,kupper))=RecvBuf_Hx_y;
	// exchange Hz with the node on the left of this one
	SendBuf_Hz_y = Hz(Range(iback,ifront),jleft,Range(klower,kupper+1));
	MPI_Sendrecv(SendBuf_Hz_y.data(),SendBuf_Hz_y.size(),MPI_DOUBLE,rank_left,HzyTag,
		RecvBuf_Hz_y.data(),RecvBuf_Hz_y.size(),MPI_DOUBLE,rank_left,HzyTag,MPI_CartSubComm,&Status);
	Hz(Range(iback,ifront),jleft-1,Range(klower,kupper+1))=RecvBuf_Hz_y;

	// exchange Hx with the node on the right of this one
	SendBuf_Hx_y = Hx(Range(iback,ifront+1),jright,Range(klower,kupper));
	MPI_Sendrecv(SendBuf_Hx_y.data(),SendBuf_Hx_y.size(),MPI_DOUBLE,rank_right,HxyTag,
		RecvBuf_Hx_y.data(),RecvBuf_Hx_y.size(),MPI_DOUBLE,rank_right,HxyTag,MPI_CartSubComm,&Status);
	Hx(Range(iback,ifront+1),jright+1,Range(klower,kupper))=RecvBuf_Hx_y;
	// exchange Hz with the node on the right of this one
	SendBuf_Hz_y = Hz(Range(iback,ifront),jright,Range(klower,kupper+1));
	MPI_Sendrecv(SendBuf_Hz_y.data(),SendBuf_Hz_y.size(),MPI_DOUBLE,rank_right,HzyTag,
		RecvBuf_Hz_y.data(),RecvBuf_Hz_y.size(),MPI_DOUBLE,rank_right,HzyTag,MPI_CartSubComm,&Status);
	Hz(Range(iback,ifront),jright+1,Range(klower,kupper+1))=RecvBuf_Hz_y;

	// exchange Hx with the node below this one
	SendBuf_Hx_z = Hx(Range(iback,ifront+1),Range(jleft,jright),klower);
	MPI_Sendrecv(SendBuf_Hx_z.data(),SendBuf_Hx_z.size(),MPI_DOUBLE,rank_below,HxzTag,
		RecvBuf_Hx_z.data(),RecvBuf_Hx_z.size(),MPI_DOUBLE,rank_below,HxzTag,MPI_CartSubComm,&Status);
	Hx(Range(iback,ifront+1),Range(jleft,jright),klower-1)=RecvBuf_Hx_z;
	// exchange Hy with the node below this one
	SendBuf_Hy_z = Hy(Range(iback,ifront),Range(jleft,jright+1),klower);
	MPI_Sendrecv(SendBuf_Hy_z.data(),SendBuf_Hy_z.size(),MPI_DOUBLE,rank_below,HyzTag,
		RecvBuf_Hy_z.data(),RecvBuf_Hy_z.size(),MPI_DOUBLE,rank_below,HyzTag,MPI_CartSubComm,&Status);
	Hy(Range(iback,ifront),Range(jleft,jright+1),klower-1)=RecvBuf_Hy_z;

	// exchange Hx with the node above this one
	SendBuf_Hx_z = Hx(Range(iback,ifront+1),Range(jleft,jright),kupper);
	MPI_Sendrecv(SendBuf_Hx_z.data(),SendBuf_Hx_z.size(),MPI_DOUBLE,rank_above,HxzTag,
		RecvBuf_Hx_z.data(),RecvBuf_Hx_z.size(),MPI_DOUBLE,rank_above,HxzTag,MPI_CartSubComm,&Status);
	Hx(Range(iback,ifront+1),Range(jleft,jright),kupper+1)=RecvBuf_Hx_z;
	// exchange Hy with the node above this one
	SendBuf_Hy_z = Hy(Range(iback,ifront),Range(jleft,jright+1),kupper);
	MPI_Sendrecv(SendBuf_Hy_z.data(),SendBuf_Hy_z.size(),MPI_DOUBLE,rank_above,HyzTag,
		RecvBuf_Hy_z.data(),RecvBuf_Hy_z.size(),MPI_DOUBLE,rank_above,HyzTag,MPI_CartSubComm,&Status);
	Hy(Range(iback,ifront),Range(jleft,jright+1),kupper+1)=RecvBuf_Hy_z;
#endif
	// Obviously, if no MPI (only one node), no field exchange is done.
	// Tangential H-field values 1/2 cell away from the grid boundary are always left zero.
}

void MPI_exit(const int& exitcode)
{//parallel version of "exit(#exitcode#)" that calls MPI_Finalize, if MPI is enabled
//#ifndef MPI_DISABLE
//	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Finalize();
//#endif
	exit(exitcode);
}
