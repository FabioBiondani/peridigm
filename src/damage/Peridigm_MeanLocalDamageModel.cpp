/*! \file Peridigm_MeanLocalDamageModel.cpp */

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

#include "Peridigm_MeanLocalDamageModel.hpp"
#include "Peridigm_Field.hpp"

#include <boost/math/constants/constants.hpp>
const double pi = boost::math::constants::pi<double>();

using namespace std;

PeridigmNS::MeanLocalDamageModel::MeanLocalDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params), m_volumeFieldId(-1), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_localdamageFieldId(-1), m_bondDamageFieldId(-1)/*, m_bondsLeftFieldId(-1)*/, m_damageFieldId(-1)
{
  m_criticalLocalDamage = params.get<double>("Critical Local Damage");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_localdamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Local_Damage");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");
//   m_bondsLeftFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bonds_Left");
  m_damageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
  
  
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_localdamageFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
//   m_fieldIds.push_back(m_bondsLeftFieldId);
  m_fieldIds.push_back(m_damageFieldId);
}

PeridigmNS::MeanLocalDamageModel::~MeanLocalDamageModel()
{
}

void
PeridigmNS::MeanLocalDamageModel::initialize(const double dt,
                                                   const int numOwnedPoints,
                                                   const int* ownedIDs,
                                                   const int* neighborhoodList,
                                                   PeridigmNS::DataManager& dataManager) const
{
  double *bondDamage, /**BondsLeft,*/ *damage;
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
//   dataManager.getData(m_bondsLeftFieldId, PeridigmField::STEP_NP1)->ExtractView(&BondsLeft);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);

  // Initialize damage to zero
  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID /*, ++BondsLeft*/){
	int nodeID = ownedIDs[iID];
    damage[nodeID] = 0.0;
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
//     *BondsLeft = numNeighbors;
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      bondDamage[bondIndex++] = 0.0;
	}
  }
}

void
PeridigmNS::MeanLocalDamageModel::computeDamage(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{
  double *volume, *x, *y;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  
  double *LocalDamageN, *LocalDamageNP1;
  dataManager.getData(m_localdamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&LocalDamageNP1);
  dataManager.getData(m_localdamageFieldId, PeridigmField::STEP_N)->ExtractView(&LocalDamageN);
  
  double *DamageNP1;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&DamageNP1);
  
  double *bondDamageN, *bondDamageNP1;
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamageN);

//   double *BondsLeftNP1;
//   dataManager.getData(m_bondsLeftFieldId, PeridigmField::STEP_NP1)->ExtractView(&BondsLeftNP1);


  const int *neighPtr = neighborhoodList;
  double* localDamageOverlap= LocalDamageNP1;
  double* volumeOverlap= volume;
  double* xOverlap = x;
  double* bondDamageOverlapN= bondDamageN;
  double* bondDamageOverlapNP1= bondDamageNP1;
  double meanDa1;
  double meanDa2;
  double Da;
  double DaP;
  double radius;  // mean radius of the element
  double radiusP;
  double distance;
  double mcoordX;
  double mcoordY;
  double mcoordZ;
  double mcoordXP;
  double mcoordYP;
  double mcoordZP;
  
//   int BondsLeft;
 
  for(int p=0;p<numOwnedPoints;p++,localDamageOverlap++,volumeOverlap++,xOverlap+=3/*,BondsLeftNP1++*/){
    int numNeigh = *neighPtr; neighPtr++;
    Da=*localDamageOverlap;
    radius=std::pow(*volumeOverlap*3/(4*pi),1.0/3.0);
    mcoordX=*xOverlap;
    mcoordY=*(xOverlap+1);
    mcoordZ=*(xOverlap+2);
    
//     BondsLeft=numNeigh;
    for(int n=0;n<numNeigh;n++,neighPtr++,bondDamageOverlapNP1++,bondDamageOverlapN++){
      int localId = *neighPtr;
      DaP = LocalDamageNP1[localId];
      radiusP = std::pow(volume[localId]*3/(4*pi),1.0/3.0);
      mcoordXP = x[3*localId];
      mcoordYP = x[3*localId+1];
      mcoordZP = x[3*localId+2];
      distance     = std::pow(std::pow(mcoordX-mcoordXP,2.0)+std::pow(mcoordY-mcoordYP,2.0)+std::pow(mcoordZ-mcoordZP,2.0),1.0/2.0);
      
      meanDa1 = Da + radius/distance*(DaP-Da);
      meanDa2 = Da + (distance-radiusP)/distance*(DaP-Da);
      if ((meanDa1>m_criticalLocalDamage) || (meanDa2>m_criticalLocalDamage)) {
          *bondDamageOverlapNP1=1.;
          //cout << Da << "  " << DaP << "  " << meanDa1 << "  " << meanDa2 << "  " << radius << "  " << radiusP << "  " << distance << "  " << endl;
      }
      else {*bondDamageOverlapNP1=*bondDamageOverlapN;}
//       if (*bondDamageOverlapNP1==1.)
//           BondsLeft-=1;
    }
//     *BondsLeftNP1 = BondsLeft;
  }

  //  Update the element damage (percent of bonds broken)
  double totalDamage;
  int neighborhoodListIndex(0), bondIndex(0);
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int nodeId = ownedIDs[iID];
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
	totalDamage = 0.0;
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
	  totalDamage += bondDamageNP1[bondIndex++];
	}
	if(numNeighbors > 0)
	  totalDamage /= numNeighbors;
	else
	  totalDamage = 0.0;
 	DamageNP1[nodeId] = totalDamage;
  }
}
