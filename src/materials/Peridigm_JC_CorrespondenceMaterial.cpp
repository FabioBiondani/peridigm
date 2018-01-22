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

#include "Peridigm_JC_CorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "JC_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include "correspondence.h"

#include <boost/math/constants/constants.hpp>
const double pi = boost::math::constants::pi<double>();

using namespace std;

PeridigmNS::JC_CorrespondenceMaterial::JC_CorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_MeltingTemperature(0.0),m_ReferenceTemperature(0.0),m_A(0.0),m_N(0.0),m_B(0.0),m_C(0.0),m_M(0.0),
    m_D1(0.0),m_D2(0.0),m_D3(0.0),m_D4(0.0),m_D5(0.0),m_DC(0.0),
    m_unrotatedRateOfDeformationFieldId(-1), m_unrotatedCauchyStressFieldId(-1), m_vonMisesStressFieldId(-1),
    m_equivalentPlasticStrainFieldId(-1),    m_accumulatedPlasticStrainFieldId(-1),    m_DamageFieldId(-1), 
    m_deltaTemperatureFieldId(-1), m_BondsLeftFieldId(-1), m_DissipationFieldId(-1)


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
  
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_unrotatedCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_vonMisesStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Von_Mises_Stress");
  m_equivalentPlasticStrainFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");
  m_accumulatedPlasticStrainFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Accumulated_Plastic_Strain");
  m_DamageFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_BondsLeftFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bonds Left");
  m_deltaTemperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");
  m_DissipationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dissipation");
  
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_vonMisesStressFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_accumulatedPlasticStrainFieldId);
  m_fieldIds.push_back(m_DamageFieldId);
  m_fieldIds.push_back(m_BondsLeftFieldId);
  m_fieldIds.push_back(m_deltaTemperatureFieldId);
  m_fieldIds.push_back(m_DissipationFieldId);
}

PeridigmNS::JC_CorrespondenceMaterial::~JC_CorrespondenceMaterial()
{
}

void
PeridigmNS::JC_CorrespondenceMaterial::initialize(const double dt,
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
  dataManager.getData(m_BondsLeftFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_BondsLeftFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_DissipationFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_DissipationFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
}

void
PeridigmNS::JC_CorrespondenceMaterial::computeCauchyStress(const double dt,
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

  double *accumulatedPlasticStrainN, *accumulatedPlasticStrainNP1;
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&accumulatedPlasticStrainNP1);
  dataManager.getData(m_accumulatedPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&accumulatedPlasticStrainN);

  double *DamageN, *DamageNP1;
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&DamageNP1);
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_N)->ExtractView(&DamageN);

  double *deltaTemperatureN, *deltaTemperatureNP1;
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperatureNP1);
  dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_N)->ExtractView(&deltaTemperatureN);

  double *DissipationNP1;
  dataManager.getData(m_DissipationFieldId, PeridigmField::STEP_NP1)->ExtractView(&DissipationNP1);
  
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
                                                m_DC,
                                                DissipationNP1
                                               );
  
  
}


void
PeridigmNS::JC_CorrespondenceMaterial::computeForce(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces and partial stress
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *horizon, *volume, *modelCoordinates, *coordinates, *velocities, *shapeTensorInverse, *deformationGradient;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);

  double *DamageN, *DamageNP1;
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&DamageNP1);
  dataManager.getData(m_DamageFieldId, PeridigmField::STEP_N)->ExtractView(&DamageN);

  double *bondDamageN, *bondDamageNP1;
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamageN);

  double *BondsLeftNP1;
  dataManager.getData(m_BondsLeftFieldId, PeridigmField::STEP_NP1)->ExtractView(&BondsLeftNP1);

  
  const int *neighPtr = neighborhoodList;

  // 
  double* bondDamageOverlapN= bondDamageN;
  double* bondDamageOverlapNP1= bondDamageNP1;
  for(int p=0;p<numOwnedPoints;p++){
    int numNeigh = *neighPtr; neighPtr++;
    for(int n=0;n<numNeigh;n++,neighPtr++,bondDamageOverlapNP1++,bondDamageOverlapN++){
      *bondDamageOverlapNP1=*bondDamageOverlapN;
    }
  }

  
  // Compute the inverse of the shape tensor and the approximate deformation gradient
  // The approximate deformation gradient will be used by the derived class (specific correspondence material model)
  // to compute the Cauchy stress.
  // The inverse of the shape tensor is stored for later use after the Cauchy stress calculation
  int shapeTensorReturnCode = 
    CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                               horizon,
                                                                               modelCoordinates,
                                                                               coordinates,
                                                                               shapeTensorInverse,
                                                                               deformationGradient,
                                                                               bondDamageNP1,
                                                                               neighborhoodList,
                                                                               numOwnedPoints);
  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);

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
                                                                                                    bondDamageNP1,
                                                                                                    neighborhoodList, 
                                                                                                    numOwnedPoints, 
                                                                                                    dt);
  string rotationTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
  rotationTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode != 0, rotationTensorErrorMessage);

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

  
  
  // FIND BROKEN BONDS
  const int *neighPtr2 = neighborhoodList;
  double* DamageOverlap= DamageNP1;
  double* volumeOverlap= volume;
  double* modelcoordinatesOverlap = coordinates;
  double meanDa1;
  double meanDa2;
  double Da;
  double DaP;
  double radius;  // mean radius of the element
  double radiusP;
  double distance;
  double mcoordX;
  double mcoordY;
  double mcoordZ;
  double mcoordXP;
  double mcoordYP;
  double mcoordZP;
  
  int BondsLeft;
  bondDamageOverlapN= bondDamageN;
  bondDamageOverlapNP1= bondDamageNP1;

  for(int p=0;p<numOwnedPoints;p++,DamageOverlap++,volumeOverlap++,modelcoordinatesOverlap+=3,BondsLeftNP1++){
    int numNeigh = *neighPtr2; neighPtr2++;
    Da=*DamageOverlap;
    radius=pow(*volumeOverlap*3/(4*pi),1/3);
    mcoordX=*modelcoordinatesOverlap;
    mcoordY=*(modelcoordinatesOverlap+1);
    mcoordZ=*(modelcoordinatesOverlap+2);
    
    BondsLeft=numNeigh;
    for(int n=0;n<numNeigh;n++,neighPtr2++,bondDamageOverlapNP1++,bondDamageOverlapN++){
      int localId = *neighPtr2;
      DaP = DamageNP1[localId];
      radiusP = pow(volume[localId]*3/(4*pi),1/3);
      mcoordXP = modelCoordinates[3*localId];
      mcoordYP = modelCoordinates[3*localId+1];
      mcoordZP = modelCoordinates[3*localId+2];
      distance = pow(pow(mcoordX-mcoordXP,2)+pow(mcoordY-mcoordYP,2)+pow(mcoordZ-mcoordZP,2),1/2);
      
      meanDa1 = Da + radius/distance*(DaP-Da);
      meanDa2 = Da + (distance-radiusP)/distance*(DaP-Da);
      if ((meanDa1>m_DC) || (meanDa2>m_DC)) {
          *bondDamageOverlapNP1=1.;
//           cout << *DamageNP1 << "  " << meanDa << "  " << Da << "  " << DaP << "  " << "\n";
          BondsLeft-=1;
      }
      else {*bondDamageOverlapNP1=*bondDamageOverlapN;}
    }
    *BondsLeftNP1 = BondsLeft;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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

  double *modelCoordinatesPtr, *neighborModelCoordinatesPtr, *forceDensityPtr, *neighborForceDensityPtr, *partialStressPtr;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double TX, TY, TZ, omega, vol, neighborVol, jacobianDeterminant;
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

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID, 
          ++delta, defGrad+=9, stress+=9, shapeTensorInv+=9){

    // first Piola-Kirchhoff stress = J * cauchyStress * defGrad^-T

    // Invert the deformation gradient and store the determinant
    int matrixInversionReturnCode =
      CORRESPONDENCE::Invert3by3Matrix(defGrad, jacobianDeterminant, defGradInv);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(matrixInversionReturnCode != 0, matrixInversionErrorMessage);
    
    //P = J * \sigma * F^(-T)
    CORRESPONDENCE::MatrixMultiply(false, true, jacobianDeterminant, stress, defGradInv, piolaStress);

    // Inner product of Piola stress and the inverse of the shape tensor
    CORRESPONDENCE::MatrixMultiply(false, false, 1.0, piolaStress, shapeTensorInv, temp);

    // Loop over the neighbors and compute contribution to force densities
    modelCoordinatesPtr = modelCoordinates + 3*iID;
    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondIndex++){

      neighborIndex = *neighborListPtr;
      neighborModelCoordinatesPtr = modelCoordinates + 3*neighborIndex;
      
      undeformedBondX = *(neighborModelCoordinatesPtr)   - *(modelCoordinatesPtr);
      undeformedBondY = *(neighborModelCoordinatesPtr+1) - *(modelCoordinatesPtr+1);
      undeformedBondZ = *(neighborModelCoordinatesPtr+2) - *(modelCoordinatesPtr+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      omega = m_OMEGA(undeformedBondLength, *delta)*(1-bondDamageNP1[bondIndex]);
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
                                        coordinates,
                                        deformationGradient,
                                        hourglassForceDensity,
                                        bondDamageNP1,
                                        neighborhoodList,
                                        numOwnedPoints,
                                        m_bulkModulus,
                                        m_hourglassCoefficient);

  // Sum the hourglass force densities into the force densities
  Teuchos::RCP<Epetra_Vector> forceDensityVector = dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1);
  Teuchos::RCP<Epetra_Vector> hourglassForceDensityVector = dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1);
  forceDensityVector->Update(1.0, *hourglassForceDensityVector, 1.0);
}
