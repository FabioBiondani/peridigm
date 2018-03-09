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
#include <iostream>

PeridigmNS::SpecularBondPosition::SpecularBondPosition(int numOwnedPoints,int* neighborhoodPtr ,int neighborListSize,int* neighborList){
  
  specularPosNeighList = new int[neighborListSize-numOwnedPoints];
  int ID=-1;
  int k=0;
  int i=0;
  while( i<neighborListSize ){
//       std::cout << i << "  " << *(neighborList+i) << std::endl;
//       std::cout << (ID<numOwnedPoints) << " " << (i==neighborhoodPtr[ID+1]) << std::endl;
      if ((ID<numOwnedPoints) && (i==*(neighborhoodPtr+ID+1))) {
          ++ID;
//           std::cout << "ID:  " << ID << std::endl;
      } else {
          int NeighPtr = *(neighborhoodPtr+ *(neighborList+i));
//           std::cout << "NP:    " << NeighPtr << std::endl;
          int j=1; int m=0;
          while (m==0) {
              if (j==(*(neighborhoodPtr+ *(neighborList+i)+1)-NeighPtr)){
                  break;
                  std::cout << "Pointer to specular node not found" << std::endl;
              } else if (*(neighborList+NeighPtr+j)==ID){
                  m=1;
                  *(specularPosNeighList+k)=NeighPtr+j-*(neighborList+i)-1;
                  ++k;
//                   std::cout << "SPNL:    " << NeighPtr+j-*(neighborList+i)-1 << std::endl;
              } else ++j;
          }
      }
      ++i;
  }

}

