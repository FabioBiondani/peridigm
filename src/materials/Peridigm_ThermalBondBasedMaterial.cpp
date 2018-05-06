/*! \file Peridigm_ThermalBondBasedMaterial.cpp */

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

#include "Peridigm_ThermalBondBasedMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "thermal_bondbased.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <iostream>

using namespace std;

PeridigmNS::ThermalBondBasedMaterial::ThermalBondBasedMaterial(const Teuchos::ParameterList& params)
  : ThermalMaterial(params),
	m_horizon(0.0),
	m_thermalShock(false),
	m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
	m_volumeFieldId(-1), m_damageFieldId(-1), m_weightedVolumeFieldId(-1), m_dilatationFieldId(-1), m_modelCoordinatesFieldId(-1),
	m_coordinatesFieldId(-1), m_bondDamageFieldId(-1),
	m_heatFlowFieldId(-1),m_deltaTemperatureFieldId(-1)
{
  //! TODO Add meaningful asserts on material properties.
	m_horizon = params.get<double>("Horizon");
    
    obj_specificHeat.set(params,"Specific Heat");
    obj_termCond.set(params,"Thermal Conductivity");
	if(m_thermalShock) obj_convectionConstant.set(matparams,"Convection Constant");

	PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
	m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
	m_damageFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
	m_weightedVolumeFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
	m_dilatationFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
	m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
	m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
	m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
	m_heatFlowFieldId        		 = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Heat_Flow");
	m_deltaTemperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");

	m_fieldIds.push_back(m_volumeFieldId);
	m_fieldIds.push_back(m_damageFieldId);
	m_fieldIds.push_back(m_weightedVolumeFieldId);
	m_fieldIds.push_back(m_dilatationFieldId);
	m_fieldIds.push_back(m_modelCoordinatesFieldId);
	m_fieldIds.push_back(m_coordinatesFieldId);
	m_fieldIds.push_back(m_bondDamageFieldId);
	m_fieldIds.push_back(m_heatFlowFieldId);
	m_fieldIds.push_back(m_deltaTemperatureFieldId);
}

PeridigmNS::ThermalBondBasedMaterial::~ThermalBondBasedMaterial()
{
}

void
PeridigmNS::ThermalBondBasedMaterial::initialize(const double dt,
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
PeridigmNS::ThermalBondBasedMaterial::computeHeatFlow(	const double dt,
															const int numOwnedPoints,
															const int* ownedIDs,
															const int* neighborhoodList,
															PeridigmNS::DataManager& dataManager) const
{

// 			Zero out the heat flow
	dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

// 	Extract pointers to the underlying data
	double *x, *y, *cellVolume, *bondDamage, *heatFlow, *deltaTemperature;

	dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
	dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
	dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
	dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
	dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->ExtractView(&heatFlow);
	dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
    
    
	MATERIAL_EVALUATION::computeBondBasedHeatFlow(	x,
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
