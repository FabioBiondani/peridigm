/*! \file Peridigm_CorrespondenceMaterial.cpp */

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

#include "Peridigm_CorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "correspondence.h"
#include <Teuchos_Assert.hpp>
#include <iostream>

using namespace std;

PeridigmNS::CorrespondenceMaterial::CorrespondenceMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0), m_hourglassCoefficient(0.0),
    m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_horizonFieldId(-1), m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_velocitiesFieldId(-1), 
    m_hourglassForceDensityFieldId(-1), m_forceDensityFieldId(-1),
    m_bondDamageFieldId(-1), m_damageFieldId(-1),
    m_deformationGradientFieldId(-1),
    m_shapeTensorInverseFieldId(-1),
    m_leftStretchTensorFieldId(-1),
    m_rotationTensorFieldId(-1), 
    m_unrotatedCauchyStressFieldId(-1),
    m_cauchyStressFieldId(-1), 
    m_unrotatedRateOfDeformationFieldId(-1),
    m_partialStressFieldId(-1),
    m_specularBondPositionFieldId(-1),
    m_microPotentialFieldId(-1),
    m_useSpecularBondPositions(false)
{
  //! \todo Add meaningful asserts on material properties.
  obj_bulkModulus.set(params);
  obj_shearModulus.set(params);
  m_bulkModulus = obj_bulkModulus.compute(0.0);
  m_shearModulus = obj_shearModulus.compute(0.0);
  //m_bulkModulus = calculateBulkModulus(params);
  //m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");
  m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");

  obj_alphaVol.set(params,"Thermal Expansion Coefficient");
  m_alphaVol= obj_alphaVol.compute(0.0);

  if (params.isParameter("Use Specular Bond Position")){
      m_useSpecularBondPositions  = params.get<bool>("Use Specular Bond Position");
  }

  
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticCorrespondence material model.\n");
//   TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticCorrespondence material model.\n");
//   TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticCorrespondence material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_horizonFieldId                    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_velocitiesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_hourglassForceDensityFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Hourglass_Force_Density");
  m_bondDamageFieldId                 = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_damageFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");
  m_leftStretchTensorFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor");
  m_rotationTensorFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  m_shapeTensorInverseFieldId         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Shape_Tensor_Inverse");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_cauchyStressFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_partialStressFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");
  if (m_useSpecularBondPositions){
      m_specularBondPositionFieldId       = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Specular_Bond_Position");
      m_microPotentialFieldId             = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Micro-Potential");
  }

  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocitiesFieldId);
  m_fieldIds.push_back(m_hourglassForceDensityFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_leftStretchTensorFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_partialStressFieldId);
  if (m_useSpecularBondPositions){
      m_fieldIds.push_back(m_specularBondPositionFieldId);
      m_fieldIds.push_back(m_microPotentialFieldId);
  }
}

PeridigmNS::CorrespondenceMaterial::~CorrespondenceMaterial()
{
}

void
PeridigmNS::CorrespondenceMaterial::initialize(const double dt,
                                               const int numOwnedPoints,
                                               const int* ownedIDs,
                                               const int* neighborhoodList,
                                               PeridigmNS::DataManager& dataManager)
{
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *leftStretchTensorN;
  double *leftStretchTensorNP1;
  double *rotationTensorN;
  double *rotationTensorNP1;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);

  //Initialize the left stretch and rotation tenor to the identity matrix
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorNP1, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorNP1, numOwnedPoints);

  //Initialize the inverse of the shape tensor and the deformation gradient
  double *volume;
  double *horizon;
  double *modelCoordinates;
  double *coordinates;
  double *shapeTensorInverse;
  double *deformationGradient;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);

  double *bondDamage;
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamage);


  int shapeTensorReturnCode = 
    CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                               horizon,
                                                                               modelCoordinates,
                                                                               coordinates,
                                                                               shapeTensorInverse,
                                                                               deformationGradient,
                                                                               bondDamage,
                                                                               neighborhoodList,
                                                                               numOwnedPoints);

  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::initialize() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
//   TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);
  if (shapeTensorReturnCode != 0) cout << shapeTensorErrorMessage;


}

void
PeridigmNS::CorrespondenceMaterial::computeForce(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces and partial stress
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *horizon, *volume, *modelCoordinates, *coordinatesN, *coordinatesNP1, *velocities, *shapeTensorInverse, *deformationGradient;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&coordinatesN);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);

  double *damage, *bondDamage;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  // Compute the inverse of the shape tensor and the approximate deformation gradient
  // The approximate deformation gradient will be used by the derived class (specific correspondence material model)
  // to compute the Cauchy stress.
  // The inverse of the shape tensor is stored for later use after the Cauchy stress calculation
  int shapeTensorReturnCode = 
    CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                               horizon,
                                                                               modelCoordinates,
                                                                               coordinatesNP1,
                                                                               shapeTensorInverse,
                                                                               deformationGradient,
                                                                               bondDamage,
                                                                               neighborhoodList,
                                                                               numOwnedPoints);
  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
//   TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);
  if (shapeTensorReturnCode != 0) cout << shapeTensorErrorMessage;

  double *leftStretchTensorN, *leftStretchTensorNP1, *rotationTensorN, *rotationTensorNP1, *unrotatedRateOfDeformation;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);

  // Compute left stretch tensor, rotation tensor, and unrotated rate-of-deformation.
  // Performs a polar decomposition via Flanagan & Taylor (1987) algorithm.
  int rotationTensorReturnCode = CORRESPONDENCE::computeUnrotatedRateOfDeformationAndRotationTensor(volume,
                                                                                                    horizon,
                                                                                                    modelCoordinates, 
                                                                                                    velocities, 
                                                                                                    deformationGradient,
                                                                                                    shapeTensorInverse,
                                                                                                    leftStretchTensorN,
                                                                                                    rotationTensorN,
                                                                                                    leftStretchTensorNP1,
                                                                                                    rotationTensorNP1,
                                                                                                    unrotatedRateOfDeformation,
                                                                                                    bondDamage,
                                                                                                    neighborhoodList, 
                                                                                                    numOwnedPoints, 
                                                                                                    dt);
  string rotationTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
  rotationTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
//   TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode != 0, rotationTensorErrorMessage);
  if (rotationTensorReturnCode != 0) cout << rotationTensorErrorMessage;
  
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
  double *unrotatedCauchyStressNP1, *cauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);

  CORRESPONDENCE::rotateCauchyStress(rotationTensorNP1,
                                     unrotatedCauchyStressNP1,
                                     cauchyStressNP1,
                                     numOwnedPoints);

  // Cauchy stress is now updated and in the rotated state.  Proceed with
  // conversion to Piola-Kirchoff and force-vector states.

  double *forceDensity;
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  double *partialStress;
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&partialStress);

  double *delta = horizon;
  double* stress = cauchyStressNP1;
  double* shapeTensorInv = shapeTensorInverse;
  double* defGrad = deformationGradient;

  double *modelCoordinatesPtr, *coordinatesPtrN, *coordinatesPtrNP1, *neighborModelCoordinatesPtr, *neighborCoordinatesPtrN, *neighborCoordinatesPtrNP1, *forceDensityPtr, *neighborForceDensityPtr, *partialStressPtr;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double TX, TY, TZ, omega, vol, neighborVol, jacobianDeterminant, deltaDeformedBondX, deltaDeformedBondY, deltaDeformedBondZ;
  int numNeighbors, neighborIndex;
  int bondIndex(0);
  
  string matrixInversionErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to invert deformation gradient.\n";
  matrixInversionErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";

  vector<double> defGradInvVector(9), piolaStressVector(9), tempVector(9);
  double* defGradInv = &defGradInvVector[0];
  double* piolaStress = &piolaStressVector[0];
  double* temp = &tempVector[0];

  double *specu, *miPotNP1, *miPotNP1overlap;
  if (m_useSpecularBondPositions){
      dataManager.getData(m_specularBondPositionFieldId, PeridigmField::STEP_NONE)->ExtractView(&specu);
      dataManager.getData(m_microPotentialFieldId, PeridigmField::STEP_NP1)->ExtractView(&miPotNP1);
  }
  miPotNP1overlap = miPotNP1;

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID, damage++,
          ++delta, defGrad+=9, stress+=9, shapeTensorInv+=9){

    // first Piola-Kirchhoff stress = J * cauchyStress * defGrad^-T

    // Invert the deformation gradient and store the determinant
    int matrixInversionReturnCode =
      CORRESPONDENCE::Invert3by3Matrix(defGrad, jacobianDeterminant, defGradInv);
//     TEUCHOS_TEST_FOR_EXCEPT_MSG(matrixInversionReturnCode != 0, matrixInversionErrorMessage);
    if (matrixInversionReturnCode != 0) cout << matrixInversionErrorMessage;
    
    //P = J * \sigma * F^(-T)
    CORRESPONDENCE::MatrixMultiply(false, true, jacobianDeterminant, stress, defGradInv, piolaStress);

    // Inner product of Piola stress and the inverse of the shape tensor
    CORRESPONDENCE::MatrixMultiply(false, false, (1.0-*damage), piolaStress, shapeTensorInv, temp);

    // Loop over the neighbors and compute contribution to force densities
    modelCoordinatesPtr = modelCoordinates + 3*iID;
    coordinatesPtrN     = coordinatesN        + 3*iID;
    coordinatesPtrNP1   = coordinatesNP1      + 3*iID;
    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondIndex++, specu++, miPotNP1++){

      neighborIndex = *neighborListPtr;
      neighborModelCoordinatesPtr = modelCoordinates + 3*neighborIndex;
      neighborCoordinatesPtrN     = coordinatesN     + 3*neighborIndex;
      neighborCoordinatesPtrNP1   = coordinatesNP1   + 3*neighborIndex;
      
      undeformedBondX = *(neighborModelCoordinatesPtr)   - *(modelCoordinatesPtr);
      undeformedBondY = *(neighborModelCoordinatesPtr+1) - *(modelCoordinatesPtr+1);
      undeformedBondZ = *(neighborModelCoordinatesPtr+2) - *(modelCoordinatesPtr+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      omega = m_OMEGA(undeformedBondLength, *delta)*(1-bondDamage[bondIndex]);
      TX = omega * ( *(temp)   * undeformedBondX + *(temp+1) * undeformedBondY + *(temp+2) * undeformedBondZ );
      TY = omega * ( *(temp+3) * undeformedBondX + *(temp+4) * undeformedBondY + *(temp+5) * undeformedBondZ );
      TZ = omega * ( *(temp+6) * undeformedBondX + *(temp+7) * undeformedBondY + *(temp+8) * undeformedBondZ );

      vol = volume[iID];
      neighborVol = volume[neighborIndex];

      forceDensityPtr = forceDensity + 3*iID;
      neighborForceDensityPtr = forceDensity + 3*neighborIndex;

      *(forceDensityPtr)   += TX * neighborVol;
      *(forceDensityPtr+1) += TY * neighborVol;
      *(forceDensityPtr+2) += TZ * neighborVol;
      *(neighborForceDensityPtr)   -= TX * vol;
      *(neighborForceDensityPtr+1) -= TY * vol;
      *(neighborForceDensityPtr+2) -= TZ * vol;

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
      
      deltaDeformedBondX = (*(neighborCoordinatesPtrNP1)   - *(coordinatesPtrNP1))   - (*(neighborCoordinatesPtrN)   - *(coordinatesPtrN));
      deltaDeformedBondY = (*(neighborCoordinatesPtrNP1+1) - *(coordinatesPtrNP1+1)) - (*(neighborCoordinatesPtrN+1) - *(coordinatesPtrN+1));
      deltaDeformedBondZ = (*(neighborCoordinatesPtrNP1+2) - *(coordinatesPtrNP1+2)) - (*(neighborCoordinatesPtrN+2) - *(coordinatesPtrN+2));

      if (m_useSpecularBondPositions){
          int specuId = int(*specu);
          *miPotNP1+=                TX*deltaDeformedBondX + TY*deltaDeformedBondY + TZ*deltaDeformedBondZ;
          miPotNP1overlap[specuId]+= TX*deltaDeformedBondX + TY*deltaDeformedBondY + TZ*deltaDeformedBondZ;
      }
    }
  }

  // Compute hourglass forces for stabilization of low-energy and/or zero-energy modes
  dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *hourglassForceDensity;
  dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&hourglassForceDensity);

  // \todo HOURGLASS FORCES ARE NOT OUTPUT TO EXODUS CORRECTLY BECAUSE THEY ARE NOT ASSEMBLED ACROSS PROCESSORS.
  //       They are summed into the force vector below, and the force vector is assembled across processors,
  //       so the calculation runs correctly, but the hourglass output is off.

  CORRESPONDENCE::computeHourglassForce(volume,
                                        horizon,
                                        modelCoordinates,
                                        coordinatesNP1,
                                        deformationGradient,
                                        hourglassForceDensity,
                                        bondDamage,
                                        neighborhoodList,
                                        numOwnedPoints,
                                        m_bulkModulus,
                                        m_hourglassCoefficient);

  // Sum the hourglass force densities into the force densities
  Teuchos::RCP<Epetra_Vector> forceDensityVector = dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1);
  Teuchos::RCP<Epetra_Vector> hourglassForceDensityVector = dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1);
  forceDensityVector->Update(1.0, *hourglassForceDensityVector, 1.0);
}
