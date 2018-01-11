/*! \file Peridigm_ThermalBB_JCCorrMaterial.cpp */

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

#include "Peridigm_ThermalBB_JCCorrMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "JC_correspondence.h"
#include "nonlocal_thermal_diffusion.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ThermalBB_JCCorrMaterial::ThermalBB_JCCorrMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_MeltingTemperature(0.0),m_ReferenceTemperature(0.0),m_A(0.0),m_N(0.0),m_B(0.0),m_C(0.0),m_M(0.0),
    m_D1(0.0),m_D2(0.0),m_D3(0.0),m_D4(0.0),m_D5(0.0),m_DC(0.0),
    m_unrotatedRateOfDeformationFieldId(-1), m_unrotatedCauchyStressFieldId(-1), m_vonMisesStressFieldId(-1),
    m_equivalentPlasticStrainFieldId(-1),    m_accumulatedPlasticStrainFieldId(-1),    m_DamageFieldId(-1),
    m_deltaTemperatureFieldId(-1),  m_heatFlowFieldId(-1)

{
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
  
  obj_termCond.set(params,"Thermal Conductivity");
  m_horizon = params.get<double>("Horizon");
  
  
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_unrotatedCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_vonMisesStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Von_Mises_Stress");
  m_equivalentPlasticStrainFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");
  m_accumulatedPlasticStrainFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Accumulated_Plastic_Strain");
  m_DamageFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_deltaTemperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");
  m_heatFlowFieldId     = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Heat_Flow");

  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_vonMisesStressFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_accumulatedPlasticStrainFieldId);
  m_fieldIds.push_back(m_DamageFieldId);
  m_fieldIds.push_back(m_deltaTemperatureFieldId);
}

PeridigmNS::ThermalBB_JCCorrMaterial::~ThermalBB_JCCorrMaterial()
{
}

void
PeridigmNS::ThermalBB_JCCorrMaterial::initialize(const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             PeridigmNS::DataManager& dataManager)
{

  PeridigmNS::CorrespondenceMaterial::initialize(dt,
                                                  numOwnedPoints,
                                                  ownedIDs,
                                                  neighborhoodList,
                                                  dataManager);

  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
}








void
PeridigmNS::ThermalBB_JCCorrMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               PeridigmNS::DataManager& dataManager) const
{
  double *unrotatedCauchyStressN, *unrotatedCauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressN);

  double *unrotatedRateOfDeformation;
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);
  
  double *vonMisesStressN, *vonMisesStressNP1;
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&vonMisesStressNP1);
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_N)->ExtractView(&vonMisesStressN);
  
  double *equivalentPlasticStrainN, *equivalentPlasticStrainNP1;
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&equivalentPlasticStrainNP1);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&equivalentPlasticStrainN);

  double *accumulatedPlasticStrainN, *accumulatedPlasticStrainNP1;
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&accumulatedPlasticStrainNP1);
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&accumulatedPlasticStrainN);

  double *DamageN, *DamageNP1;
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&DamageNP1);
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_N)->ExtractView(&DamageN);

  double *deltaTemperatureN, *deltaTemperatureNP1;
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperatureNP1);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_N)->ExtractView(&deltaTemperatureN);
  
//   *deltaTemperatureNP1=m_ReferenceTemperature;
  CORRESPONDENCE::updateJohnsonCookCauchyStress(unrotatedRateOfDeformation, 
                                                unrotatedCauchyStressN, 
                                                unrotatedCauchyStressNP1, 
                                                vonMisesStressNP1,
                                                equivalentPlasticStrainN, 
                                                accumulatedPlasticStrainN, 
                                                DamageN, 
                                                equivalentPlasticStrainNP1, 
                                                accumulatedPlasticStrainNP1, 
                                                DamageNP1, 
                                                numOwnedPoints, 
                                                obj_bulkModulus, 
                                                obj_shearModulus,
                                                obj_alphaVol,
                                                deltaTemperatureN,
                                                deltaTemperatureNP1,
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

void
PeridigmNS::ThermalBB_JCCorrMaterial::computeHeatFlow(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{

//             Zero out the heat flow
    dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

//     Extract pointers to the underlying data
    double *x, *y, *cellVolume, *bondDamage, *heatFlow, *deltaTemperature;

    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
    dataManager.getData(m_DamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage); ////////////////////
    dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->ExtractView(&heatFlow);
    dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
    
    
    MATERIAL_EVALUATION::computeHeatFlow(x,
                                         y,
                                         cellVolume,
                                         bondDamage,
                                         heatFlow,
                                         neighborhoodList,
                                         numOwnedPoints,
                                         obj_termCond,
                                         m_horizon,
                                         deltaTemperature);
}


