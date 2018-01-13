/*! \file Peridigm_TupekCorrespondenceDamageModel.cpp */

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

#include "Peridigm_TupekCorrespondenceDamageModel.hpp"
#include "Peridigm_Field.hpp"

using namespace std;

PeridigmNS::TupekCorrespondenceDamageModel::TupekCorrespondenceDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params), m_damageFieldId(-1), m_bondDamageFieldId(-1)
{

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_damageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
}

PeridigmNS::TupekCorrespondenceDamageModel::~TupekCorrespondenceDamageModel()
{
}

void
PeridigmNS::TupekCorrespondenceDamageModel::initialize(const double dt,
                                                   const int numOwnedPoints,
                                                   const int* ownedIDs,
                                                   const int* neighborhoodList,
                                                   PeridigmNS::DataManager& dataManager) const
{
  double *damage, *bondDamage;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  // Initialize damage to zero
  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int nodeID = ownedIDs[iID];
    damage[nodeID] = 0.0;
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      bondDamage[bondIndex++] = 0.0;
	}
  }
}

void
PeridigmNS::TupekCorrespondenceDamageModel::computeDamage(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{
  double *damage, *bondDamageN, *bondDamageNP1;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamageN);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);

  int neighborhoodListIndex(0), bondIndex(0);
  int nodeId, numNeighbors, neighborID, iID, iNID;

  // Set the bond damage to the previous value
  *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N));

  // Update the bond damage
  // Break bonds if the extension is greater than the critical extension

  for(iID=0 ; iID<numOwnedPoints ; ++iID){
	nodeId = ownedIDs[iID];
	numNeighbors = neighborhoodList[neighborhoodListIndex++];
	for(iNID=0 ; iNID<numNeighbors ; ++iNID){
	  neighborID = neighborhoodList[neighborhoodListIndex++];
      if((damage[nodeId] > 0.2) || (damage[neighborID] > 0.2) || bondDamageN[bondIndex] == 1.0){
        bondDamageNP1[bondIndex] = 1.0;
      }
      bondIndex += 1;
    }
  }
  
}
