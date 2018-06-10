//! \file viscousmaxwell_ordinary.cxx

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

#include <cmath>
#include <iostream>
#include "viscousmaxwell_correspondence.h"
#include "material_utilities.h"

using std::cout;
using std::endl;
namespace CORRESPONDENCE {

void computeViscousMaxwellCauchyStress
  (
   const double *unrotatedCauchyStressN,
   const double *unrotatedCauchyStressNP1,
   const double *internalVariablesN,
   double *internalVariablesNP1,
   double *unrotatedViscousCauchyStress,
   const double* damage,
   const double* deltaTemperature,
   int numOwnedPoints,
   PeridigmNS::Material::TempDepConst obj_lambda,
   PeridigmNS::Material::TempDepConst obj_tau,
   const double dt
)
{
    double m_tau(1.0), m_lambda(1.0);
    double c1(0.0), decay(0.0), beta_i(0.0);
    if (deltaTemperature==nullptr){
        m_tau = obj_tau.compute(0.0);
        m_lambda = obj_lambda.compute(0.0);
        if ((m_lambda!=0.0)&&(m_tau!=0.0)){
            c1 = m_tau / dt;
            decay = exp(-1.0/c1);
            beta_i=1.-c1*(1.-decay);
        }
    }
    
    const double *deltaT    = deltaTemperature;
    const double* stressN   = unrotatedCauchyStressN;
    const double* stressNP1 = unrotatedCauchyStressNP1;
    const double* internalN = internalVariablesN;
    double* internalNP1     = internalVariablesNP1;
    double* vstress         = unrotatedViscousCauchyStress;

    double devstressN[9];
    double devstressNP1[9];
    double isostress;

    if ((m_lambda!=0.0)&&(m_tau!=0.0))
	for(int p=0;p<numOwnedPoints;p++, deltaT++, damage++, stressN+=9, stressNP1+=9, internalN+=9, internalNP1+=9, vstress+=9){

        if(deltaTemperature){
            m_tau = obj_tau.compute(*deltaT);
            m_lambda = obj_lambda.compute(*deltaT);
            if ((m_lambda==0.0)||(m_tau==0.0)) continue;
            c1 = m_tau / dt;
            decay = exp(-1.0/c1);
            beta_i=1.-c1*(1.-decay);
        }

        if (*damage==1.0) continue;

        isostress = ( stressN[0]   + stressN[4]   + stressN[8]   )/3.0;
		for(int n=0;n<9;n++){
            devstressN[n] = stressN[n];
        }
        devstressN[0] -= isostress;
        devstressN[4] -= isostress;
        devstressN[8] -= isostress;

        isostress = ( stressNP1[0] + stressNP1[4] + stressNP1[8] )/3.0;
		for(int n=0;n<9;n++){
            devstressNP1[n] = stressNP1[n];
        }
        devstressNP1[0] -= isostress;
        devstressNP1[4] -= isostress;
        devstressNP1[8] -= isostress;

        
        for(int n=0;n<9;n++){
            *(internalNP1+n) = m_lambda * devstressN[n] * (1-decay) +
                               *(internalN+n) * decay  +
                               m_lambda * beta_i * (devstressNP1[n]-devstressN[n]);
            
            *(vstress+n) = m_lambda * devstressNP1[n] - *(internalNP1+n);
        }
	}
}

}

