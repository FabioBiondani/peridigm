/*! \file Peridigm_NeighborhoodData.hpp */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#ifndef PERIDIGM_NEIGHBORHOODDATA_HPP
#define PERIDIGM_NEIGHBORHOODDATA_HPP

// #include "mpi.h"

namespace PeridigmNS {

class NeighborhoodData {

public:

  NeighborhoodData() 
    : numOwnedPoints(0), ownedIDs(0), neighborhoodListSize(0), neighborhoodList(0), neighborhoodPtr(0), overlapNeighborhoodListSize(0), overlapNeighborhoodList(0), specularBondPositions(0){}

  NeighborhoodData(const NeighborhoodData& other)
    : numOwnedPoints(0), ownedIDs(0), neighborhoodListSize(0), neighborhoodList(0), neighborhoodPtr(0), overlapNeighborhoodListSize(0), overlapNeighborhoodList(0), specularBondPositions(0)
  {
    SetNumOwned(other.NumOwnedPoints());
    SetNeighborhoodListSize(other.NeighborhoodListSize());
    SetOverlapNeighborhoodListSize(other.OverlapNeighborhoodListSize());
    memcpy(ownedIDs, other.ownedIDs, numOwnedPoints*sizeof(int));
    memcpy(neighborhoodPtr, other.neighborhoodPtr, numOwnedPoints*sizeof(int));
    memcpy(neighborhoodList, other.neighborhoodList, neighborhoodListSize*sizeof(int));
    memcpy(overlapNeighborhoodList, other.overlapNeighborhoodList, overlapNeighborhoodListSize*sizeof(int));
  }

  ~NeighborhoodData(){
	if(ownedIDs != 0)
	  delete[] ownedIDs;
	if(neighborhoodList != 0)
	  delete[] neighborhoodList;
    if(neighborhoodPtr != 0)
      delete[] neighborhoodPtr;
	if(overlapNeighborhoodList != 0)
	  delete[] overlapNeighborhoodList;
  }

  void SetNumOwned(int numOwned){
	numOwnedPoints = numOwned;
	if(ownedIDs != 0)
	  delete[] ownedIDs;
	ownedIDs = new int[numOwned];
    if(neighborhoodPtr != 0)
      delete[] neighborhoodPtr;
    neighborhoodPtr = new int[numOwned];
	if(neighborhoodListSize != 0)
        specularBondPositions = new int[neighborhoodListSize-numOwnedPoints];
  }

  void SetNeighborhoodListSize(int neighborhoodSize){
	neighborhoodListSize = neighborhoodSize;
	if(neighborhoodList != 0)
	  delete[] neighborhoodList;
	neighborhoodList = new int[neighborhoodListSize];
	if(numOwnedPoints != 0)
        specularBondPositions = new int[neighborhoodListSize-numOwnedPoints];
  }

  int NumOwnedPoints() const{
	return numOwnedPoints;
  }

  int* OwnedIDs() const{
	return ownedIDs;
  }

  int* NeighborhoodPtr() const{
	return neighborhoodPtr;
  }

  int NeighborhoodListSize() const{
	return neighborhoodListSize;
  }

  int* NeighborhoodList() const{
	return neighborhoodList;
  }

  void SetOverlapNeighborhoodListSize(int overlapneighborhoodSize){
	overlapNeighborhoodListSize = overlapneighborhoodSize;
	if(overlapNeighborhoodList != 0)
	  delete[] overlapNeighborhoodList ;
	overlapNeighborhoodList  = new int[overlapNeighborhoodListSize];
  }

  int OverlapNeighborhoodListSize() const{
	return overlapNeighborhoodListSize;
  }

  int* OverlapNeighborhoodList() const{
	return overlapNeighborhoodList;
  }

  void SetSpecularBondPositions(Teuchos::RCP<const Epetra_BlockMap> overlapScalarBondMap){
    int GID1;
    for (int i=0;i<(neighborhoodListSize-numOwnedPoints);++i){
        for(int j=0;j<overlapScalarBondMap->NumMyElements();++j){
            if ((i-overlapScalarBondMap->FirstPointInElement(j))<overlapScalarBondMap->ElementSize(j)) {
                GID1 = overlapScalarBondMap->GID(j);
                break;
            }
        }
        int GID2 = *(overlapNeighborhoodList+i);
        int LID2 = overlapScalarBondMap->LID(GID2);
        int numNeighborsSpecularPoint = overlapScalarBondMap->ElementSize(LID2);
        int firstpoint = overlapScalarBondMap->FirstPointInElement(LID2);
        int j=0;
        for (j=0;j<numNeighborsSpecularPoint;++j){
            if (*(overlapNeighborhoodList+firstpoint+j)==GID1){
                *(specularBondPositions+i) = firstpoint+j;
                break;
            }
        }
//         int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//         std::cout << "Proc: " << rank << "  Pos: " << i << "  GID1: " << GID1 << "  GID2: " << GID2 << "  Specu: " << firstpoint+j << std::endl;
    }
  }

  int SpecularBondPositionsSize() const{
	return (neighborhoodListSize-numOwnedPoints);
  }

  int* SpecularBondPositions() const{
	return specularBondPositions;
  }

  double memorySize() const{
    int sizeInBytes =
      (numOwnedPoints + 2*neighborhoodListSize + overlapNeighborhoodListSize + 2 )*sizeof(int) + 4*sizeof(int*);
    double sizeInMegabytes = sizeInBytes/1048576.0;
    return sizeInMegabytes;
  }

protected:
  int numOwnedPoints;
  int numOverlapPoints;
  int* ownedIDs;
  int neighborhoodListSize;
  int* neighborhoodList;
  int* neighborhoodPtr;
  int overlapNeighborhoodListSize;
  int* overlapNeighborhoodList;
  int* specularBondPositions;
};

}

#endif // PERIDIGM_NEIGHBORHOODDATA_HPP
