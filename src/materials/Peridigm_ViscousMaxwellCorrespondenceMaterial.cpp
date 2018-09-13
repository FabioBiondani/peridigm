/*! \file Peridigm_ViscousMaxwellCorrespondenceMaterial.cpp */

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

#include "Peridigm_ViscousMaxwellCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "viscousmaxwell_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ViscousMaxwellCorrespondenceMaterial::ViscousMaxwellCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : ViscousCorrespondenceMaterial(params), m_analysisHasThermal(false), m_internalVariablesFieldId(-1), m_unrotatedCauchyStressFieldId(-1), m_deltaTemperatureFieldId(-1)
{
  obj_bulkModulus.set(params);
  obj_shearModulus.set(params);
  obj_lambda.set(params,"Maxwell Model Stiffness Ratio");
  obj_tau.set(params,"Relaxation Time");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  if (fieldManager.hasField("Temperature_Change")||params.isParameter("Thermal Expansion Coefficient"))
    m_analysisHasThermal = true;
  m_unrotatedCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_internalVariablesFieldId     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Internal_Variables");
  if (m_analysisHasThermal)
    m_deltaTemperatureFieldId      = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");

  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);  
  m_fieldIds.push_back(m_internalVariablesFieldId);
  if(m_analysisHasThermal)
    m_fieldIds.push_back(m_deltaTemperatureFieldId);  
}

PeridigmNS::ViscousMaxwellCorrespondenceMaterial::~ViscousMaxwellCorrespondenceMaterial()
{
}

void
PeridigmNS::ViscousMaxwellCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               const int* neighborhoodList,
                                                               PeridigmNS::DataManager& dataManager) const
{
  double *unrotatedCauchyStressN, *unrotatedCauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressN);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
  double *unrotatedViscousCauchyStress;
  dataManager.getData(m_unrotatedViscousCauchyStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedViscousCauchyStress);
  double *internalVariablesN, *internalVariablesNP1;
  dataManager.getData(m_internalVariablesFieldId, PeridigmField::STEP_N)->ExtractView(&internalVariablesN);
  dataManager.getData(m_internalVariablesFieldId, PeridigmField::STEP_NP1)->ExtractView(&internalVariablesNP1);
  double *damage;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  double *deltaTemperature(NULL);
  if (m_analysisHasThermal)
    dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);

  CORRESPONDENCE::computeViscousMaxwellCauchyStress(unrotatedCauchyStressN,
                                                    unrotatedCauchyStressNP1,
                                                    internalVariablesN,
                                                    internalVariablesNP1,
                                                    unrotatedViscousCauchyStress,
                                                    damage,
                                                    deltaTemperature,
                                                    numOwnedPoints,
                                                    obj_lambda,
                                                    obj_tau,
                                                    dt);
}
