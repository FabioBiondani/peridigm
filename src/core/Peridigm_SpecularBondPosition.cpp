/*! \file Peridigm_Discretization.cpp */

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

#include "Peridigm_SpecularBondPosition.hpp"
#include <vector>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Ionit_Initializer.h>
#include <iostream>

using namespace std;

PeridigmNS::SpecularBondPosition::SpecularBondPosition(const Teuchos::RCP<const Epetra_Comm>& comm, Teuchos::RCP<const Epetra_BlockMap> map1D, Teuchos::RCP<const Epetra_BlockMap> map1Doverlap, Teuchos::RCP<const Epetra_BlockMap> bondMap, int neighborListSize,int* neighborList){

    // vector of number of neighbors
    Epetra_Vector numNeigh(*map1D);
    double* numNeighPtr;
    numNeigh.ExtractView(&numNeighPtr);
    int i=0;
    while (  i<neighborListSize ){
        *(numNeighPtr++)+= *(neighborList+i);
        i+= *(neighborList+i)+1;
    }

    // vector of number of neighbors with ghosted points
    Epetra_Import importerMap1D(*map1Doverlap,*map1D);
    Epetra_Vector numNeighOverlap(*map1Doverlap);
    numNeighOverlap.Import(numNeigh,importerMap1D,Insert);

    // create bondMap with ghosted points
    int* elementSizeList  = new int[map1Doverlap->NumMyElements()];
    int* myGlobalElements = new int[map1Doverlap->NumMyElements()];
    for(int i=0;i<map1Doverlap->NumMyElements();++i){
        *(elementSizeList+i)  = numNeighOverlap[i];
        *(myGlobalElements+i) = map1Doverlap->GID(i);
    }

    bondMapOverlap = Teuchos::rcp(new Epetra_BlockMap(-1, map1Doverlap->NumMyElements(), myGlobalElements, elementSizeList, 0, *comm));
    delete[] myGlobalElements;
    delete[] elementSizeList;

    // 
    Epetra_Vector NeighborsGID(*bondMap);
    int j=0;
    for (int i=0;i<map1D->NumMyElements();++i,++j){
        int numNeigh_i = numNeigh[i];
        for(int k=0;k<numNeigh_i;++k,++j){
            NeighborsGID[j-i]=map1Doverlap->GID( *(neighborList+j+1) );
        }
    }

    //
    Epetra_Import importerBondMap(*bondMapOverlap,*bondMap);
    NeighborsGIDoverlap = Teuchos::rcp(new Epetra_Vector(*bondMapOverlap));
    NeighborsGIDoverlap->Import(NeighborsGID,importerBondMap,Insert);
//     cout << NeighborsGIDoverlap << endl;
}
// // // /*
// // //     //
// // //     specularBondPos = Teuchos::rcp(new Epetra_Vector(*bondMap));
// // //     double* specularBondPosPtr;
// // //     specularBondPos->ExtractView(&specularBondPosPtr);
// // //     double* NeighborsGIDoverlapPtr;
// // //     NeighborsGIDoverlap->ExtractView(&NeighborsGIDoverlapPtr);
// // //     int GID1;
// // //     for (int i=0;i<(neighborListSize-map1D->NumMyElements());++i){
// // //         for(int j=0;j<map1D->NumMyElements();++j){
// // //             if (i>=bondMap->FirstPointInElement(j)) {
// // //                 GID1 = map1D->GID(j);
// // //                 continue;
// // //             }
// // //         }
// // //         int GID2 = NeighborsGID[i];
// // //         int LID2 = map1Doverlap->LID(GID2);
// // //         int numNeighborsSpecularPoint = numNeighOverlap[LID2];
// // //         int firstpoint = bondMapOverlap->FirstPointInElement(LID2);
// // //         for (int j=0;j<numNeighborsSpecularPoint;++j){
// // //             if (*(NeighborsGIDoverlapPtr+j)==GID1){
// // //                 *(specularBondPosPtr+i) = firstpoint+j;
// // //                 continue;
// // //             }
// // //         }
// // //     }
// // // //     cout << *specularBondPos << endl;
// // // 
// // // }*/

