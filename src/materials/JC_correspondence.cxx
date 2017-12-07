//! \file JC_correspondence.cxx

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

#include "JC_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <math.h>

namespace CORRESPONDENCE {

template<typename ScalarT>
void updateJohnsonCookCauchyStress
(
const ScalarT* unrotatedRateOfDeformation, 
const ScalarT* cauchyStressN, 
ScalarT* cauchyStressNP1, 
ScalarT* vonMisesStress,
const ScalarT* equivalentPlasticStrainN,
const ScalarT* accumulatedPlasticStrainN,
const ScalarT* DamageN,
ScalarT* equivalentPlasticStrainNP1,
ScalarT* accumulatedPlasticStrainNP1,
ScalarT* DamageNP1,
const int numPoints, 
const double bulkMod,
const double shearMod,
const double dt,
const ScalarT* TemperatureN,
const double MeltingTemperature,
const double ReferenceTemperature,
const double constA,
const double constN,
const double constB,
const double constC,
const double constM,
const double constD1,
const double constD2,
const double constD3,
const double constD4,
const double constD5
)
{
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* stressN = cauchyStressN;
  ScalarT* stressNP1 = cauchyStressNP1;

  ScalarT* vmStress = vonMisesStress;

  const ScalarT* eqpsN = equivalentPlasticStrainN;
  const ScalarT* dapsN = accumulatedPlasticStrainN;
  const ScalarT* DaN = DamageN;
  
  ScalarT hmlgT = (*TemperatureN - ReferenceTemperature) / (MeltingTemperature - ReferenceTemperature) ; // Homologous Temperature
  
  ScalarT* eqpsNP1 = equivalentPlasticStrainNP1;
  ScalarT* dapsNP1 = accumulatedPlasticStrainNP1;
  ScalarT* DaNP1 = DamageNP1;
  
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];
  ScalarT deviatoricStressN[9];
  ScalarT deviatoricStressMagnitudeN;

  ScalarT deviatoricStressNP1[9];
  ScalarT deviatoricStressMagnitudeNP1;

  ScalarT tempA[9];
  ScalarT tempB[9];

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldStressN;
  
  
  for(int iID=0 ; iID<numPoints ; ++iID, rateOfDef+=9, stressN+=9,
        stressNP1+=9, ++vmStress,++eqpsN,++eqpsNP1,
      ++dapsN,++dapsNP1,++DaN,++DaNP1){

      //strainInc = dt * rateOfDef
      for (int i = 0; i < 9; i++) {
          strainInc[i] = *(rateOfDef+i)*dt;
          deviatoricStrainInc[i] = strainInc[i];
      }

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //Compute an elastic ``trial stress''
      for (int i = 0; i < 9; i++) {
          *(stressNP1+i) = *(stressN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(stressNP1) += bulkMod*dilatationInc;
      *(stressNP1+4) += bulkMod*dilatationInc;
      *(stressNP1+8) += bulkMod*dilatationInc;

      sphericalStressNP1 = (*(stressNP1) + *(stressNP1+4) + *(stressNP1+8))/3.0;

      // Compute the ``trial'' von Mises stress
      for (int i = 0; i < 9; i++) {
          deviatoricStressNP1[i] = *(stressNP1+i);
      }
      deviatoricStressNP1[0] -= sphericalStressNP1;
      deviatoricStressNP1[4] -= sphericalStressNP1;
      deviatoricStressNP1[8] -= sphericalStressNP1;

      // Compute \sigma_ij * \sigma_ij
      tempScalar = 0.0;
      for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
              tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
      }

      *vmStress = sqrt(3.0/2.0*tempScalar);
      
      // update strains and damage
      *eqpsNP1 = *eqpsN;
      *dapsNP1 = *dapsN;
      *DaNP1 = *DaN;

      
      // calculate actual yield stress with previous equivalent plastic strain and damage (actual limit of elastic behaviour)
      // consider daps_dot=0
      yieldStressN = (1-DaN)*(constA+constB*dapsN^constN)*(1+0)^constC*(1-hmlgT^constM);
      
      //If true, the step is plastic and we need to return to the yield
      //surface.  
      if(*vmStress > yieldStressN) {
          
          // FIND PLASTIC STRAIN AND DAMAGE USING NEWTON'S METHOD
          
          
          // eqps : equivalent plastic strain
          // feps : failure equivalent plastic strain
          // daps : damage accumulated plastic strain
          // daps_dot : daps time derived
          // Da   : damage
          // sigmaA : hydrostatic stress ratio
          // hmlgT : homologous temperature
          
          // daps_dotNP1 = (dapsNP1-dapsN)/dt = daps_dotNP1[eqpsNP1]
          // eqps_dotNP1 = (eqpsNP1-eqpsN)/dt = eqps_dotNP1[eqpsNP1]
          // daps_dotNP1 = eqps_dotNP1*(1-DaNP1)
          // dapsNP1 = dapsN + (eqpsNP1-eqpsN)*(1-DaNP1) = dapsNP1[eqpsNP1,DaNP1]
          
          // yieldStressNP1 = (1-DaNP1)*(constA+constB*dapsNP1^constN)*(1+daps_dotNP1)^constC*(1-hmlgT^constM) = yieldStressNP1[DaNP1,eqpsNP1]
          // scalarDeviatoricStrainInc = scalarDeviatoricEnergyDensity / deviatoricStressMagnitudeNP1
          // fun1 = scalarDeviatoricStrainInc - eqpsNP1 + eqpsN -1/(2*shearMod)*(sqrt(2/3)*yieldStressNP1-yieldStressN) = fun1[DaNP1,eqpsNP1] = 0
          
          // sigmaA_NP1 = sphericalStressNP1/yieldStressNP1 = sigmaA_NP1[DaNP1,eqpsNP1] // how to obtain sphericalStress during iterations?
          
          // fepsNP1 = (constD1+constD2*exp(constD3*sigmaA_NP1))*(1+eqps_dotNP1)^constD4*(1+constD5hmlgT)
          
          // DaNP1 = DaNP1 + (eqpsNP1-eqpsN)/fepsNP1
          // fun2 = DaNP1 - DaN - (eqpsNP1-eqpsN)/fepsNP1 = 0
          
          // params: dt hmlgT shearMod constA constN constB constC constM constD1 constD2 constD3 constD4 constD5 / yieldStressN DaN eqpsN dapsN / scalarDeviatoricStrainInc sphericalStressNP1
          
          // EXPRESS fun1 and fun2 function of UNKNOWNS DaNP1 and eqpsNP1
          
          // fun1 = scalarDeviatoricStrainInc + eqpsN - eqpsNP1 + (yieldStressN - (2^(1/2)*3^(1/2)*(hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)/3)/(2*shearMod)
          // fun2 = DaNP1 - DaN + (eqpsN - eqpsNP1)/((1 - (eqpsN - eqpsNP1)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)))*(constD5hmlgT + 1))
          
          // DERIVE fun1 and fun2 wrt DaNP1 and eqpsNP1 to obtain the Jacobian for the Newton's method
          
          // fun1_eqpsNP1 = ((2^(1/2)*3^(1/2)*constB*constN*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^(constN - 1)*(hmlgT^constM - 1)*(DaNP1 - 1)^2*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)/3 + (2^(1/2)*3^(1/2)*constC*(hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)^2*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^(constC - 1))/(3*dt))/(2*shearMod) - 1
          // fun1_DaNP1 = -((2^(1/2)*3^(1/2)*(hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)/3 + (2^(1/2)*3^(1/2)*constC*(hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(eqpsN - eqpsNP1)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^(constC - 1))/(3*dt) + (2^(1/2)*3^(1/2)*constB*constN*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^(constN - 1)*(hmlgT^constM - 1)*(eqpsN - eqpsNP1)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)/3)/(2*shearMod)
          
          // fun2_eqpsNP1 = - 1/((1 - (eqpsN - eqpsNP1)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)))*(constD5hmlgT + 1)) - (constD4*(eqpsN - eqpsNP1))/(dt*(1 - (eqpsN - eqpsNP1)/dt)^(constD4 + 1)*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)))*(constD5hmlgT + 1)) - (constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC))*((constC*constD3*sphericalStressNP1)/(dt*(hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^(constC + 1)) + (constB*constD3*constN*sphericalStressNP1*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^(constN - 1))/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)^2*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC))*(eqpsN - eqpsNP1))/((1 - (eqpsN - eqpsNP1)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)))^2*(constD5hmlgT + 1))
          // fun2_DaNP1 = (constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC))*(eqpsN - eqpsNP1)*((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)^2*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC) + (constC*constD3*sphericalStressNP1*(eqpsN - eqpsNP1))/(dt*(hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^(constC + 1)) + (constB*constD3*constN*sphericalStressNP1*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^(constN - 1)*(eqpsN - eqpsNP1))/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)^2*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)))/((1 - (eqpsN - eqpsNP1)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(dapsN + (eqpsN - eqpsNP1)*(DaNP1 - 1))^constN)*(DaNP1 - 1)*(((eqpsN - eqpsNP1)*(DaNP1 - 1))/dt + 1)^constC)))^2*(constD5hmlgT + 1)) + 1
          
          
          // Avoid divide-by-zero
          deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));
          
          // Compute \S_ij * \epsilon_inc_ij
          tempScalar = 0.0;
          for (int j = 0; j < 3; j++) {
              for (int i = 0; i < 3; i++) {
                  tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStrainInc[i+3*j];
              }
          }

          ScalarT scalarDeviatoricStrainInc = tempScalar / deviatoricStressMagnitudeNP1;
          
          ScalarT fun1;
          ScalarT fun2;
          ScalarT fun1_eqps;
          ScalarT fun2_eqps;
          ScalarT fun1_Da;
          ScalarT fun2_Da;
          ScalarT daps_dotNP1;
          ScalarT yieldStressNP1;
          
          ScalarT eqps_it = *eqpsN;
          ScalarT Da_it   = *DaN;
          
          
          ScalarT detM, invM11, invM12, invM21, invM22;
  
          
          fun1 = scalarDeviatoricStrainInc + *eqpsN - eqps_it + (yieldStressN - (2^(1/2)*3^(1/2)*(hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)/3)/(2*shearMod);
          
          fun2 = Da_it - *DaN + (*eqpsN - eqps_it)/((1 - (*eqpsN - eqps_it)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)))*(constD5*(hmlgT) + 1)) ;
          
          while ( ((fun1^2+fun2^2)^0.5) > 1e-7 ) {
              
              fun1_eqps = ((2^(1/2)*3^(1/2)*constB*constN*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^(constN - 1)*(hmlgT^constM - 1)*(Da_it - 1)^2*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)/3 + (2^(1/2)*3^(1/2)*constC*(hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)^2*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^(constC - 1))/(3*dt))/(2*shearMod) - 1 ;
              
              fun1_Da = -((2^(1/2)*3^(1/2)*(hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)/3 + (2^(1/2)*3^(1/2)*constC*(hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(*eqpsN - eqps_it)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^(constC - 1))/(3*dt) + (2^(1/2)*3^(1/2)*constB*constN*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^(constN - 1)*(hmlgT^constM - 1)*(*eqpsN - eqps_it)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)/3)/(2*shearMod) ;
              
              
              fun2_eqps = - 1/((1 - (*eqpsN - eqps_it)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)))*(constD5*(hmlgT) + 1)) - (constD4*(*eqpsN - eqps_it))/(dt*(1 - (*eqpsN - eqps_it)/dt)^(constD4 + 1)*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)))*(constD5*(hmlgT) + 1)) - (constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC))*((constC*constD3*sphericalStressNP1)/(dt*(hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^(constC + 1)) + (constB*constD3*constN*sphericalStressNP1*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^(constN - 1))/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)^2*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC))*(*eqpsN - eqps_it))/((1 - (*eqpsN - eqps_it)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)))^2*(constD5*(hmlgT) + 1)) ;
              
              fun2_Da = (constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC))*(*eqpsN - eqps_it)*((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)^2*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC) + (constC*constD3*sphericalStressNP1*(*eqpsN - eqps_it))/(dt*(hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^(constC + 1)) + (constB*constD3*constN*sphericalStressNP1*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^(constN - 1)*(*eqpsN - eqps_it))/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)^2*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)))/((1 - (*eqpsN - eqps_it)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)))^2*(constD5*(hmlgT) + 1)) + 1 ;
              
              
              /* 
               *  /      \   /                    \   / /         \   /          \ \   /     \
               *  | fun1 |   | fun1_Da  fun1_eqps |   | |  Da_it  |   |  Da_old  | |   |  0  |
               *  |      | + |                    | * | |         | - |          | | = |     |
               *  | fun2 |   | fun2_Da  fun2_eqps |   | | eqps_it |   | eqps_old | |   |  0  |
               *  \      /   \                    /   \ \         /   \          / /   \     /
               *     -f     =          M             *   (   x_it    -    x_old   )
               *    x_it   =  -inv(M) * f + x_old
               */
              
              detM   = fun1_Da*fun2_eqps - fun2_Da*fun1_eqps ;
              invM11 =  fun2_eqps / detM ;
              invM22 =  fun1_Da   / detM ;
              invM12 = -fun2_Da   / detM ;
              invM21 = -fun1_eqps / detM ;
              
              Da_it   = -invM11*fun1 - invM12*fun2 + Da_it;
              eqps_it = -invM21*fun1 - invM22*fun2 + eqps_it;
              
              
              // update fun1 and fun2 with the new value of eqps_it and Da_it
              fun1 = scalarDeviatoricStrainInc + *eqpsN - eqps_it + (yieldStressN - (2^(1/2)*3^(1/2)*(hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)/3)/(2*shearMod);
              
              fun2 = Da_it - *DaN + (*eqpsN - eqps_it)/((1 - (*eqpsN - eqps_it)/dt)^constD4*(constD1 + constD2*exp((constD3*sphericalStressNP1)/((hmlgT^constM - 1)*(constA + constB*(*dapsN + (*eqpsN - eqps_it)*(Da_it - 1))^constN)*(Da_it - 1)*(((*eqpsN - eqps_it)*(Da_it - 1))/dt + 1)^constC)))*(constD5*(hmlgT) + 1));
              
          };
          if ((Da_it<=*DaN) && (eqps_it<=*eqpsN)) {
              // not damaged more
          }
          else if ((Da_it>=*DaN) && (eqps_it>=*eqpsN)) {
              *eqpsNP1=eqps_it;
              *DaNP1=Da_it;
              *dapsNP1 = *dapsN + (*eqpsNP1-*eqpsN)*(1-*DaNP1);
              daps_dotNP1 = (*dapsNP1-*dapsN)/dt;
              yieldStressNP1 = (1-*DaNP1)*(constA+constB*(*dapsNP1)^constN)*(1+daps_dotNP1)^constC*(1-hmlgT^constM);
              
              
              //For perfectly plastic deformations we just have a constant factor to 
              //multiply the trial deviatoric stress
              tempScalar = yieldStressN/yieldStressNP1;
              
              // Return the deviatoric stress to the yield surface
              for (int i = 0; i < 9; i++) {
                  deviatoricStressNP1[i] *= tempScalar; 
                  *(stressNP1+i) = deviatoricStressNP1[i];
              }
              
              // Update the Cauchy Stress
              *stressNP1 += sphericalStressNP1;
              *(stressNP1+4) += sphericalStressNP1;
              *(stressNP1+8) += sphericalStressNP1;
              
              // Update the von Mises stress now that the state of stress is on the
              // yield surface
              tempScalar = 0.0;
              for (int j = 0; j < 3; j++) {
                  for (int i = 0; i < 3; i++) {
                      tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
                  }
              }
              
              *vmStress = yieldStressNP1;
          }
          else {
          // error
          };
          
          
          
          
      } else {
          // The step is elastic
      };

  }
}

// Explicit template instantiation for double
template void updateJohnsonCookCauchyStress<double>
(
const double* unrotatedRateOfDeformation, 
const double* cauchyStressN, 
double* cauchyStressNP1, 
double* vonMisesStress,
const double* equivalentPlasticStrainN,
const double* accumulatedPlasticStrainN,
const double* DamageN,
double* equivalentPlasticStrainNP1,
double* accumulatedPlasticStrainNP1,
double* DamageNP1,
const int numPoints, 
const double bulkMod,
const double shearMod,
const double dt,
const double* TemperatureN,
const double MeltingTemperature,
const double ReferenceTemperature,
const double constA,
const double constN,
const double constB,
const double constC,
const double constM,
const double constD1,
const double constD2,
const double constD3,
const double constD4,
const double constD5
)

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void updateJohnsonCookCauchyStress<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformation, 
const Sacado::Fad::DFad<double>* cauchyStressN, 
Sacado::Fad::DFad<double>* cauchyStressNP1, 
Sacado::Fad::DFad<double>* vonMisesStress,
const Sacado::Fad::DFad<double>* equivalentPlasticStrainN,
const Sacado::Fad::DFad<double>* accumulatedPlasticStrainN,
const Sacado::Fad::DFad<double>* DamageN,
Sacado::Fad::DFad<double>* equivalentPlasticStrainNP1,
Sacado::Fad::DFad<double>* accumulatedPlasticStrainNP1,
Sacado::Fad::DFad<double>* DamageNP1,
const int numPoints, 
const double bulkMod,
const double shearMod,
const double dt,
Sacado::Fad::DFad<double>* TemperatureN,
const double MeltingTemperature,
const double ReferenceTemperature,
const double constA,
const double constN,
const double constB,
const double constC,
const double constM,
const double constD1,
const double constD2,
const double constD3,
const double constD4,
const double constD5
);

}
