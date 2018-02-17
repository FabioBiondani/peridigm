/*! \file Peridigm_JohnsonCookOrdinaryMaterial.cpp */

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
    m_bulkModulus(0.0), m_shearModulus(0.0), m_alpha(0.0),  m_density(0.0), m_horizon(0.0),
    m_applySurfaceCorrectionFactor(false), m_applyThermalStrains(false),
    m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_volumeFieldId(-1), m_damageFieldId(-1), m_weightedVolumeFieldId(-1), m_dilatationFieldId(-1), m_modelCoordinatesFieldId(-1),
    m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1), m_surfaceCorrectionFactorFieldId(-1),
    m_deltaTemperatureFieldId(-1),
    m_ElasticEnergyDensityFieldId(-1),m_VonMisesStressFieldId(-1),
    m_deviatoricPlasticExtensionFieldId(-1),m_equivalentPlasticStrainFieldId(-1),m_accumulatedPlasticStrainFieldId(-1),m_LocalDamageFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  obj_bulkModulus.set(params);
  obj_shearModulus.set(params);
  m_bulkModulus = obj_bulkModulus.compute(0.0);
  m_shearModulus = obj_shearModulus.compute(0.0);
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");
  
  m_MeltingTemperature = params.get<double>("Melting Temperature");
  m_ReferenceTemperature = params.get<double>("Reference Temperature");
  m_A  = params.get<double>("Constant A");
  m_N  = params.get<double>("Constant N");
  m_B  = params.get<double>("Constant B");
  m_C  = params.get<double>("Constant C");
  m_M  = params.get<double>("Constant M");
  m_D1 = params.get<double>("Constant D1");
  m_D2 = params.get<double>("Constant D2");
  m_D3 = params.get<double>("Constant D3");
  m_D4 = params.get<double>("Constant D4");
  m_D5 = params.get<double>("Constant D5");
  m_DC = params.get<double>("Constant DC");

  

  if(params.isParameter("Thermal Expansion Coefficient")){
    obj_alphaVol.set(params,"Thermal Expansion Coefficient");
    m_alpha= obj_alphaVol.compute(0.0);
    m_applyThermalStrains = true;
  }
  
  if(params.isParameter("Apply Shear Correction Factor"))
    m_applySurfaceCorrectionFactor = params.get<bool>("Apply Shear Correction Factor");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Volume");
  m_damageFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Damage");
  m_weightedVolumeFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Weighted_Volume");
  m_dilatationFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Dilatation");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Bond_Damage");
  m_surfaceCorrectionFactorFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Surface_Correction_Factor");
  if(m_applyThermalStrains)
    m_deltaTemperatureFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Temperature_Change");

  m_ElasticEnergyDensityFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Elastic_Energy_Density");
  m_VonMisesStressFieldId                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Von_Mises_Stress");

  m_deviatoricPlasticExtensionFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Deviatoric_Plastic_Extension");
  m_equivalentPlasticStrainFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");
  m_accumulatedPlasticStrainFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Accumulated_Plastic_Strain");
  m_LocalDamageFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Local_Damage");


  
  
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_weightedVolumeFieldId);
  m_fieldIds.push_back(m_dilatationFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_surfaceCorrectionFactorFieldId);
  if(m_applyThermalStrains)
    m_fieldIds.push_back(m_deltaTemperatureFieldId);
  
  m_fieldIds.push_back(m_ElasticEnergyDensityFieldId);
  m_fieldIds.push_back(m_VonMisesStressFieldId);
  
  m_fieldIds.push_back(m_deviatoricPlasticExtensionFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_accumulatedPlasticStrainFieldId);
  m_fieldIds.push_back(m_LocalDamageFieldId);
  
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

  dataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  if(m_applySurfaceCorrectionFactor){
    Epetra_Vector temp(*dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1));
	int lengthYOverlap = temp.MyLength();
	double  *yOverlap,  *surfaceCorrectionFactor;
	temp.ExtractView(&yOverlap);
	dataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->ExtractView(&surfaceCorrectionFactor);
    MATERIAL_EVALUATION::computeShearCorrectionFactor(numOwnedPoints,lengthYOverlap,xOverlap,yOverlap,cellVolumeOverlap,weightedVolume,neighborhoodList,m_horizon,surfaceCorrectionFactor);
  }

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
  dataManager.getData(m_ElasticEnergyDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_VonMisesStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *weightedVolume, *dilatation, *bondDamage, *scf, *force, *deltaTemperature;

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->ExtractView(&scf);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);
  deltaTemperature = NULL;
  if(m_applyThermalStrains)
    dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
  
  double *W, *sigmaVM;
  
  dataManager.getData(m_ElasticEnergyDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&W);
  dataManager.getData(m_VonMisesStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&sigmaVM);
  
  double *edN, *eqpsN, *dapsN, *DaN, *edNP1, *eqpsNP1, *dapsNP1, *DaNP1;
  dataManager.getData(m_deviatoricPlasticExtensionFieldId, PeridigmField::STEP_N)->ExtractView(&edN);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&eqpsN);
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&dapsN);
  dataManager.getData(m_LocalDamageFieldId, PeridigmField::STEP_N)->ExtractView(&DaN);
  dataManager.getData(m_deviatoricPlasticExtensionFieldId, PeridigmField::STEP_NP1)->ExtractView(&edNP1);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&eqpsNP1);
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&dapsNP1);
  dataManager.getData(m_LocalDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&DaNP1);


  MATERIAL_EVALUATION::computeDilatation(x,y,weightedVolume,cellVolume,bondDamage,dilatation,neighborhoodList,numOwnedPoints,m_horizon,m_OMEGA,m_alpha,deltaTemperature);
  MATERIAL_EVALUATION::computeInternalForceJohnsonCookOrdinary(
      x,
      y,
      weightedVolume,
      cellVolume,
      dilatation,
      bondDamage,
      scf,
      force,
      neighborhoodList,
      numOwnedPoints,
      obj_bulkModulus,
      obj_shearModulus,
      obj_alphaVol,
      m_horizon,
      W,
      sigmaVM,
      edN,
      edNP1,
      eqpsN,
      eqpsNP1,
      dapsN,
      dapsNP1,
      DaN,
      DaNP1,
      deltaTemperature,
      dt,
      m_MeltingTemperature,
      m_ReferenceTemperature,
      m_A,
      m_N,
      m_B,
      m_C,
      m_M,
      m_D1,
      m_D2,
      m_D3,
      m_D4,
      m_D5,
      m_DC
);
}
