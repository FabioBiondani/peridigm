/*! \file Peridigm_ViscousCorrespondenceMaterial.cpp */

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

#include "Peridigm_BaseViscousCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>

using namespace std;

PeridigmNS::ViscousCorrespondenceMaterial::ViscousCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : ViscousMaterial(params),
    m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_horizonFieldId(-1), m_volumeFieldId(-1), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1), m_damageFieldId(-1),
    m_deformationGradientFieldId(-1), m_shapeTensorInverseFieldId(-1), m_rotationTensorFieldId(-1), m_partialStressFieldId(-1),
    m_unrotatedViscousCauchyStressFieldId(-1), m_viscousCauchyStressFieldId(-1),
    m_singularityFieldId(-1), m_singularityDetachment(true)
{
    
	PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    m_horizonFieldId                    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
    m_volumeFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
    m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
    m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
    m_bondDamageFieldId                 = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
    m_damageFieldId                     = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");

    m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");
    m_shapeTensorInverseFieldId         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Shape_Tensor_Inverse");
    m_rotationTensorFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
    m_partialStressFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");

    m_unrotatedViscousCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Viscous_Cauchy_Stress");
    m_viscousCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Viscous_Cauchy_Stress");

	m_fieldIds.push_back(m_horizonFieldId);
	m_fieldIds.push_back(m_volumeFieldId);
	m_fieldIds.push_back(m_modelCoordinatesFieldId);
	m_fieldIds.push_back(m_coordinatesFieldId);
	m_fieldIds.push_back(m_forceDensityFieldId);
	m_fieldIds.push_back(m_bondDamageFieldId);
	m_fieldIds.push_back(m_damageFieldId);

	m_fieldIds.push_back(m_deformationGradientFieldId);
	m_fieldIds.push_back(m_shapeTensorInverseFieldId);
	m_fieldIds.push_back(m_rotationTensorFieldId);
	m_fieldIds.push_back(m_partialStressFieldId);

    m_fieldIds.push_back(m_unrotatedViscousCauchyStressFieldId);
    m_fieldIds.push_back(m_viscousCauchyStressFieldId);
    
    if (params.isParameter("Singularity Detachment")){
      m_singularityDetachment  = params.get<bool>("Singularity Detachment");
    }
    if (m_singularityDetachment){
      m_singularityFieldId                = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::CONSTANT, "Correspondence_Singularity");
    }

}

PeridigmNS::ViscousCorrespondenceMaterial::~ViscousCorrespondenceMaterial()
{
}

void
PeridigmNS::ViscousCorrespondenceMaterial::initialize(const double dt,
											const int numOwnedPoints,
											const int* ownedIDs,
											const int* neighborhoodList,
											PeridigmNS::DataManager& dataManager)
{
  dataManager.getData(m_unrotatedViscousCauchyStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_viscousCauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_viscousCauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
}


void
PeridigmNS::ViscousCorrespondenceMaterial::computeForce(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager) const
{
  // Zero out the viscous forces
  dataManager.getData(m_unrotatedViscousCauchyStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  double *horizon, *volume, *modelCoordinates;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);

  double *damage, *bondDamage, *singu;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  if(m_singularityDetachment){
      dataManager.getData(m_singularityFieldId, PeridigmField::STEP_NONE)->ExtractView(&singu);
  }else singu = nullptr;
  
  // Evaluate the Cauchy stress using the routine implemented in the derived class (specific correspondence material model)
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  // computeCauchyStress() typically uses the following fields which are accessed via the DataManager:
  //   Input:  unrotated rate-of-deformation tensor
  //   Input:  unrotated Cauchy stress at step N
  //   Input:  internal state data (managed in the derived class)
  //   Output: unrotated Cauchy stress at step N+1
  computeCauchyStress(dt, numOwnedPoints, neighborhoodList, dataManager);

  // rotate back to the Eulerian frame
  double *rotationTensorNP1, *unrotatedViscousCauchyStressNP1, *viscousCauchyStressNP1;
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);
  dataManager.getData(m_unrotatedViscousCauchyStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedViscousCauchyStressNP1);
  dataManager.getData(m_viscousCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&viscousCauchyStressNP1);

  CORRESPONDENCE::rotateCauchyStress(rotationTensorNP1,
                                     unrotatedViscousCauchyStressNP1,
                                     viscousCauchyStressNP1,
                                     numOwnedPoints);

  // Cauchy stress is now updated and in the rotated state.  Proceed with
  // conversion to Piola-Kirchoff and force-vector states.

  double *forceDensity;
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  double *partialStress;
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&partialStress);

  double *delta = horizon;
  double* stress = viscousCauchyStressNP1;

  double *defGrad, *shapeTensorInv;
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&defGrad);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInv);

  double *modelCoordinatesPtr, *neighborModelCoordinatesPtr,
  *forceDensityPtr, *neighborForceDensityPtr, *partialStressPtr;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double TX, TY, TZ, omega, vol, neighborVol, jacobianDeterminant;
  int numNeighbors, neighborIndex;
  int bondIndex(0);

  vector<double> defGradInvVector(9), piolaStressVector(9), tempVector(9);
  double* defGradInv = &defGradInvVector[0];
  double* piolaStress = &piolaStressVector[0];
  double* temp = &tempVector[0];

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID, 
          ++delta, defGrad+=9, stress+=9, shapeTensorInv+=9, ++damage, ++singu){

    numNeighbors = *neighborListPtr; neighborListPtr++;
    if ((m_singularityDetachment)&&(*singu==1.0)){
        neighborListPtr+=numNeighbors;
        bondIndex+=numNeighbors;
        continue;
    }

    // first Piola-Kirchhoff stress = J * cauchyStress * defGrad^-T

    // Invert the deformation gradient and store the determinant
    CORRESPONDENCE::Invert3by3Matrix(defGrad, jacobianDeterminant, defGradInv);

    //P = J * \sigma * F^(-T)
    CORRESPONDENCE::MatrixMultiply(false, true, jacobianDeterminant, stress, defGradInv, piolaStress);

    // Inner product of Piola stress and the inverse of the shape tensor
    CORRESPONDENCE::MatrixMultiply(false, false, (1.0-*damage), piolaStress, shapeTensorInv, temp);

    // Loop over the neighbors and compute contribution to force densities
    modelCoordinatesPtr = modelCoordinates + 3*iID;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondIndex++){

      neighborIndex = *neighborListPtr;
      neighborModelCoordinatesPtr = modelCoordinates + 3*neighborIndex;

      undeformedBondX = *(neighborModelCoordinatesPtr)   - *(modelCoordinatesPtr);
      undeformedBondY = *(neighborModelCoordinatesPtr+1) - *(modelCoordinatesPtr+1);
      undeformedBondZ = *(neighborModelCoordinatesPtr+2) - *(modelCoordinatesPtr+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      omega = m_OMEGA(undeformedBondLength, *delta)*(1.0-bondDamage[bondIndex]);
      TX = omega * ( *(temp)   * undeformedBondX + *(temp+1) * undeformedBondY + *(temp+2) * undeformedBondZ );
      TY = omega * ( *(temp+3) * undeformedBondX + *(temp+4) * undeformedBondY + *(temp+5) * undeformedBondZ );
      TZ = omega * ( *(temp+6) * undeformedBondX + *(temp+7) * undeformedBondY + *(temp+8) * undeformedBondZ );

      vol = volume[iID];
      neighborVol = volume[neighborIndex];

      forceDensityPtr = forceDensity + 3*iID;
      neighborForceDensityPtr = forceDensity + 3*neighborIndex;

      // Sum the viscous force densities into the force densities
      *(forceDensityPtr)            += TX * neighborVol;
      *(forceDensityPtr+1)          += TY * neighborVol;
      *(forceDensityPtr+2)          += TZ * neighborVol;
      *(neighborForceDensityPtr)    -= TX * vol;
      *(neighborForceDensityPtr+1)  -= TY * vol;
      *(neighborForceDensityPtr+2)  -= TZ * vol;

      partialStressPtr = partialStress + 9*iID;
      *(partialStressPtr)   += TX*undeformedBondX*neighborVol;
      *(partialStressPtr+1) += TX*undeformedBondY*neighborVol;
      *(partialStressPtr+2) += TX*undeformedBondZ*neighborVol;
      *(partialStressPtr+3) += TY*undeformedBondX*neighborVol;
      *(partialStressPtr+4) += TY*undeformedBondY*neighborVol;
      *(partialStressPtr+5) += TY*undeformedBondZ*neighborVol;
      *(partialStressPtr+6) += TZ*undeformedBondX*neighborVol;
      *(partialStressPtr+7) += TZ*undeformedBondY*neighborVol;
      *(partialStressPtr+8) += TZ*undeformedBondZ*neighborVol;
      
    }
  }
}




