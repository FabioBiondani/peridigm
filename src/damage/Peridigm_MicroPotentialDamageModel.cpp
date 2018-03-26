/*! \file Peridigm_MicroPotentialDamageModel.cpp */

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

#include <cmath>
#include "Peridigm_MicroPotentialDamageModel.hpp"
#include "Peridigm_Field.hpp"

using namespace std;

PeridigmNS::MicropotentialDamageModel::MicropotentialDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params), m_Jintegral(0.0), m_modelCoordinatesFieldId(-1), m_horizonFieldId(-1), m_damageFieldId(-1), m_bondDamageFieldId(-1), m_deltaTemperatureFieldId(-1),m_microPotentialFieldId(-1), m_specularBondPositionFieldId(-1)
{
  obj_Jintegral.set(params,"J_integral");
  m_Jintegral= obj_Jintegral.compute(0.0);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_horizonFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_damageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");
  m_deltaTemperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");
  m_microPotentialFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Micro-Potential");
  m_specularBondPositionFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Specular_Bond_Position");


  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_deltaTemperatureFieldId);
  m_fieldIds.push_back(m_microPotentialFieldId);
  m_fieldIds.push_back(m_specularBondPositionFieldId);
}

PeridigmNS::MicropotentialDamageModel::~MicropotentialDamageModel()
{
}

void
PeridigmNS::MicropotentialDamageModel::initialize(const double dt,
                                                   const int numOwnedPoints,
                                                   const int* ownedIDs,
                                                   const int* neighborhoodList,
                                                   PeridigmNS::DataManager& dataManager) const
{
  double *damage, *bondDamage/*, *BondsLeft*/;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  // Initialize damage to zero
  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID /*, ++BondsLeft*/){
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
PeridigmNS::MicropotentialDamageModel::computeDamage(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{
  double *x, *horizon, *damage, *bondDamageN, *bondDamageNP1, *deltaTemperature, *miPot, *specu;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamageN);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
  dataManager.getData(m_microPotentialFieldId, PeridigmField::STEP_N)->ExtractView(&miPot);
  dataManager.getData(m_specularBondPositionFieldId, PeridigmField::STEP_NONE)->ExtractView(&specu);

  double trialDamage(0.0), totalDamage;
  int neighborhoodListIndex(0), bondIndex(0);
  int nodeId, numNeighbors, neighborID, iID, iNID;
  double nodeInitialX[3], initialDistance, bond_Jintegral;

  // Set the bond damage to the previous value
  *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N));

  // Update the bond damage
  // Break bonds if the extension is greater than the critical extension

  for(iID=0 ; iID<numOwnedPoints ; ++iID){
	nodeId = ownedIDs[iID];
	nodeInitialX[0] = x[nodeId*3];
	nodeInitialX[1] = x[nodeId*3+1];
	nodeInitialX[2] = x[nodeId*3+2];
    double localT = *(deltaTemperature+nodeId);
    numNeighbors = neighborhoodList[neighborhoodListIndex++];
//     *BondsLeftNP1 = numNeighbors;
	for(iNID=0 ; iNID<numNeighbors ; ++iNID){
	  neighborID = neighborhoodList[neighborhoodListIndex++];
      initialDistance = 
        distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
                 x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);
      double neighT = *(deltaTemperature+neighborID);
      
      bond_Jintegral = obj_Jintegral.compute((localT+neighT)/2.0);

      double m_criticalMicroPotential = 4.0/(m_pi*pow(*(horizon+nodeId),4.0))*bond_Jintegral;
//       double m_criticalMicroPotential = 12.0/(11.0*m_pi*pow(*(horizon+nodeId),4.0))*bond_Jintegral;
//       double m_criticalMicroPotential = 5.0/(m_pi*pow(*(horizon+nodeId),5.0))*initialDistance*bond_Jintegral;

      double bondMicroPotential = miPot[bondIndex] ;

      if (bondMicroPotential != miPot[int(specu[bondIndex])] ) cout << "MALISSSIMOOOOOOOOOOOOO" << endl;

      trialDamage = 0.0;
      if(bondMicroPotential > m_criticalMicroPotential)
        trialDamage = 1.0;
      if(trialDamage > bondDamageNP1[bondIndex])
        bondDamageNP1[bondIndex] = trialDamage;
      bondIndex += 1;
    }
  }

  //  Update the element damage (percent of bonds broken)

  neighborhoodListIndex = 0;
  bondIndex = 0;
  for(iID=0 ; iID<numOwnedPoints ; ++iID){
	nodeId = ownedIDs[iID];
	numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
	totalDamage = 0.0;
	for(iNID=0 ; iNID<numNeighbors ; ++iNID){
	  totalDamage += bondDamageNP1[bondIndex++];
	}
	if(numNeighbors > 0)
	  totalDamage /= numNeighbors;
	else
	  totalDamage = 0.0;
 	damage[nodeId] = totalDamage;
  }
}
