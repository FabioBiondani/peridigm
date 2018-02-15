/*! \file Peridigm_ElasticMaterial.cpp */

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

#include "Peridigm_JohnsonCookOrdinaryMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "jc_ordinary.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::JohnsonCookOrdinaryMaterial::JohnsonCookOrdinaryMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_bulkModulus(0.0), m_shearModulus(0.0), m_density(0.0), m_alpha(0.0), m_horizon(0.0),
    m_applyThermalStrains(false),
    m_computePartialStress(false),
    m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_volumeFieldId(-1), m_damageFieldId(-1), m_weightedVolumeFieldId(-1), m_dilatationFieldId(-1), m_modelCoordinatesFieldId(-1),
    m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_partialStressFieldId(-1), m_bondDamageFieldId(-1),
    m_deltaTemperatureFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");

  if(params.isParameter("Thermal Expansion Coefficient")){
    m_alpha = params.get<double>("Thermal Expansion Coefficient");
    m_applyThermalStrains = true;
  }

  if(params.isParameter("Compute Partial Stress"))
    m_computePartialStress = params.get<bool>("Compute Partial Stress");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Volume");
  m_damageFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Damage");
  m_weightedVolumeFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Weighted_Volume");
  m_dilatationFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Dilatation");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Bond_Damage");
  if(m_applyThermalStrains)
    m_deltaTemperatureFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Temperature_Change");
  if(m_computePartialStress)
    m_partialStressFieldId         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_weightedVolumeFieldId);
  m_fieldIds.push_back(m_dilatationFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  if(m_applyThermalStrains)
    m_fieldIds.push_back(m_deltaTemperatureFieldId);
  if(m_computePartialStress)
    m_fieldIds.push_back(m_partialStressFieldId);
}

PeridigmNS::JohnsonCookOrdinaryMaterial::~JohnsonCookOrdinaryMaterial()
{
}

void
PeridigmNS::JohnsonCookOrdinaryMaterial::initialize(const double dt,
                                        const int numOwnedPoints,
                                        const int* ownedIDs,
                                        const int* neighborhoodList,
                                        PeridigmNS::DataManager& dataManager)
{
  // Extract pointers to the underlying data
  double *xOverlap,  *cellVolumeOverlap, *weightedVolume;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolumeOverlap);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);

  MATERIAL_EVALUATION::computeWeightedVolume(xOverlap,cellVolumeOverlap,weightedVolume,numOwnedPoints,neighborhoodList,m_horizon);

}

void
PeridigmNS::JohnsonCookOrdinaryMaterial::computeForce(const double dt,
                                          const int numOwnedPoints,
                                          const int* ownedIDs,
                                          const int* neighborhoodList,
                                          PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  if(m_computePartialStress)
    dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *weightedVolume, *dilatation, *bondDamage, *force, *deltaTemperature, *partialStress;

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);
  deltaTemperature = NULL;
  if(m_applyThermalStrains)
    dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
  partialStress = NULL;
  if(m_computePartialStress)
    dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&partialStress);

  MATERIAL_EVALUATION::computeDilatation(x,y,weightedVolume,cellVolume,bondDamage,dilatation,neighborhoodList,numOwnedPoints,m_horizon,m_OMEGA,m_alpha,deltaTemperature);
  MATERIAL_EVALUATION::computeInternalForceJohnsonCookOrdinary(x,y,weightedVolume,cellVolume,dilatation,bondDamage,force,partialStress,neighborhoodList,numOwnedPoints,m_bulkModulus,m_shearModulus,m_horizon,m_alpha,deltaTemperature);
}

// void
// PeridigmNS::JohnsonCookOrdinaryMaterial::computeStoredElasticEnergyDensity(const double dt,
//                                                                const int numOwnedPoints,
//                                                                const int* ownedIDs,
//                                                                const int* neighborhoodList,
//                                                                PeridigmNS::DataManager& dataManager) const
// {
//   // This function is intended to be called from a compute class.
//   // The compute class should have already created the Stored_Elastic_Energy_Density field id.
//   int storedElasticEnergyDensityFieldId = PeridigmNS::FieldManager::self().getFieldId("Stored_Elastic_Energy_Density");
// 
//   double *x, *y, *cellVolume, *weightedVolume, *dilatation, *storedElasticEnergyDensity, *bondDamage;
//   dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
//   dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
//   dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
//   dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
//   dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
//   dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
//   dataManager.getData(storedElasticEnergyDensityFieldId, PeridigmField::STEP_NONE)->ExtractView(&storedElasticEnergyDensity);
// 
//   double *deltaTemperature = NULL;
//   if(m_applyThermalStrains)
//     dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
// 
//   int iID, iNID, numNeighbors, nodeId, neighborId;
//   double omega, nodeInitialX[3], nodeCurrentX[3];
//   double initialDistance, currentDistance, deviatoricExtension, neighborBondDamage;
//   double nodeDilatation, alpha, temp;
// 
//   int neighborhoodListIndex(0), bondIndex(0);
//   for(iID=0 ; iID<numOwnedPoints ; ++iID){
// 
//     nodeId = ownedIDs[iID];
//     nodeInitialX[0] = x[nodeId*3];
//     nodeInitialX[1] = x[nodeId*3+1];
//     nodeInitialX[2] = x[nodeId*3+2];
//     nodeCurrentX[0] = y[nodeId*3];
//     nodeCurrentX[1] = y[nodeId*3+1];
//     nodeCurrentX[2] = y[nodeId*3+2];
//     nodeDilatation = dilatation[nodeId];
//     alpha = 15.0*m_shearModulus/weightedVolume[nodeId];
// 
//     temp = 0.0;
// 
//     numNeighbors = neighborhoodList[neighborhoodListIndex++];
//     for(iNID=0 ; iNID<numNeighbors ; ++iNID){
//       neighborId = neighborhoodList[neighborhoodListIndex++];
//       neighborBondDamage = bondDamage[bondIndex++];
//       initialDistance = 
//         distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
//                  x[neighborId*3], x[neighborId*3+1], x[neighborId*3+2]);
//       currentDistance = 
//         distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
//                  y[neighborId*3], y[neighborId*3+1], y[neighborId*3+2]);
//       if(m_applyThermalStrains)
// 	currentDistance -= m_alpha*deltaTemperature[nodeId]*initialDistance;
//       deviatoricExtension = (currentDistance - initialDistance) - nodeDilatation*initialDistance/3.0;
//       omega=m_OMEGA(initialDistance,m_horizon);
//       temp += (1.0-neighborBondDamage)*omega*deviatoricExtension*deviatoricExtension*cellVolume[neighborId];
//     }
//     storedElasticEnergyDensity[nodeId] = 0.5*m_bulkModulus*nodeDilatation*nodeDilatation + 0.5*alpha*temp;
//   }
// }
