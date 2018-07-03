/*! \file Peridigm_JC_CorrespondenceMaterial.cpp */

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

#include "Peridigm_JohnsonCookCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "jc_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include "correspondence.h"
#include <math.h>

using namespace std;

PeridigmNS::JohnsonCookCorrespondenceMaterial::JohnsonCookCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_MeltingTemperature(0.0),m_ReferenceTemperature(0.0),m_A(0.0),m_N(0.0),m_B(0.0),m_C(0.0),m_M(0.0),m_doteqps0(1.0),
    m_Beta(0.0),
    m_unrotatedRateOfDeformationFieldId(-1), m_unrotatedCauchyStressFieldId(-1), m_vonMisesStressFieldId(-1),
    m_equivalentPlasticStrainFieldId(-1),m_bondDamageFieldId(-1),
    m_deltaTemperatureFieldId(-1),
    m_cumulativeHeatFieldId(-1)
{
  if (params.isParameter("Constant A")){
      m_A  = params.get<double>("Constant A");
      m_N  = params.get<double>("Constant N");
      m_B  = params.get<double>("Constant B");
      m_C  = params.get<double>("Constant C");
      m_M  = params.get<double>("Constant M");
      m_MeltingTemperature = params.get<double>("Melting Temperature");
      m_ReferenceTemperature = params.get<double>("Reference Temperature");
      if (params.isParameter("Reference Strain Rate"))
          m_doteqps0 = params.get<double>("Reference Strain Rate");
  } else {
      m_A = 1e200;
      m_N = 0.0;
      m_B = 0.0;
      m_C = 0.0;
      m_M = 1e-100;
      m_MeltingTemperature = 1e100;
      m_ReferenceTemperature = 0.0;
  }
  
  if (params.isParameter("Beta"))
      m_Beta = params.get<double>("Beta");
  else if (params.isParameter("Beta Taylor-Quinney"))
      m_Beta = params.get<double>("Beta Taylor-Quinney");
  else
      m_Beta = 1.0;

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_vonMisesStressFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Von_Mises_Stress");
  m_equivalentPlasticStrainFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");
  m_bondDamageFieldId                 = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_deltaTemperatureFieldId           = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");
  
  m_cumulativeHeatFieldId       = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cumulative_Adiabatic_Heat");

  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_vonMisesStressFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_deltaTemperatureFieldId);
  m_fieldIds.push_back(m_cumulativeHeatFieldId);
}

PeridigmNS::JohnsonCookCorrespondenceMaterial::~JohnsonCookCorrespondenceMaterial()
{
}

void
PeridigmNS::JohnsonCookCorrespondenceMaterial::initialize(const double dt,
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
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
}

void
PeridigmNS::JohnsonCookCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               const int* neighborhoodList,
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

  double *deltaTemperatureN, *deltaTemperatureNP1;
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperatureNP1);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_N)->ExtractView(&deltaTemperatureN);

  CORRESPONDENCE::updateJohnsonCookCauchyStress(
      unrotatedRateOfDeformation, 
      unrotatedCauchyStressN, 
      unrotatedCauchyStressNP1, 
      vonMisesStressNP1,
      equivalentPlasticStrainN, 
      equivalentPlasticStrainNP1, 
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
      m_doteqps0
  );
  
  double *cumulativeHeatN, *cumulativeHeatNP1;
  dataManager.getData(m_cumulativeHeatFieldId, PeridigmField::STEP_N)->ExtractView(&cumulativeHeatN);
  dataManager.getData(m_cumulativeHeatFieldId, PeridigmField::STEP_NP1)->ExtractView(&cumulativeHeatNP1);

  for(int iID=0 ; iID<numOwnedPoints ; ++iID, ++cumulativeHeatN, ++cumulativeHeatNP1, ++vonMisesStressNP1, ++equivalentPlasticStrainN, ++equivalentPlasticStrainNP1){
      *cumulativeHeatNP1 = *cumulativeHeatN + *vonMisesStressNP1 * (*equivalentPlasticStrainNP1 - *equivalentPlasticStrainN);
  }
}
