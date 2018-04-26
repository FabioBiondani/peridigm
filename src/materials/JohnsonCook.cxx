//! \file JohnsonCook.cxx

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



#include <Sacado.hpp>
#include <math.h>
#include <JohnsonCook.h>


namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void JohnsonCookSolve
(
const ScalarT vmStressTrial,
const ScalarT* eqpsN,
ScalarT* eqpsNP1,
ScalarT* yieldStress,
const double shearModulus,
const double constA,
const double constN,
const double constB,
const double constC,
const double pow_hmlgT_M,
const double dt
)
{
    // FIND PLASTIC STRAIN INCREMENT USING NEWTON'S METHOD

    // eqps : equivalent plastic strain
    // Deqps : increment of equivalent plastic strain
    // hmlgT : homologous temperature
    
    // EXPRESS yield function function of UNKNOWN eqps
    // DERIVE yield function wrt deqps to obtain the Jacobian for the Newton's method

    ScalarT lambda = 0.0;
    ScalarT yieldFunction = 1.0;
    ScalarT yieldFunction_lambda;
    ScalarT teqps;
    ScalarT teqps_lambda;
    
    ScalarT pow_eqps_nM1;
    ScalarT pow_eqps_n;
    ScalarT pow_1teqps_C;
    ScalarT pow_1teqps_CM1;
    
    ScalarT yieldStress_lambda;
    
    
    int it=0;
    int max_it = 20;
    while (std::abs(yieldFunction) > 1e-7){
        if (it>0)
        {
            // x_it   =  -inv(M) * f + x_old
            lambda   = -yieldFunction/yieldFunction_lambda + lambda;
        }

        ++it;
        *eqpsNP1= *eqpsN + lambda; // equivalent plastic strain
        teqps = lambda/dt; // time derived eqps
        teqps_lambda= 1./dt;

        pow_eqps_n = 0.;
        if (*eqpsNP1>0.){pow_eqps_n = +(pow(*eqpsNP1,constN));}
        else if (*eqpsNP1<0.){pow_eqps_n = +(pow(-*eqpsNP1,constN));}
        
        pow_eqps_nM1 = 0.;
        if (*eqpsNP1>0.){pow_eqps_nM1 = +(pow(*eqpsNP1,constN-1.0));}
        else if(*eqpsNP1<0.){pow_eqps_nM1 = +(pow(-*eqpsNP1,constN-1.0));}
        
        pow_1teqps_C = 0.;
        if (1.0+teqps>0.){pow_1teqps_C = +(pow(1.0+teqps,constC));}
        else if (1.0+teqps<0.){pow_1teqps_C = +(pow(-1.0-teqps,constC));}
        
        pow_1teqps_CM1 = 0.;
        if (1.0+teqps>0.){pow_1teqps_CM1 = +(pow(1.0+teqps,constC-1.0));}
        else if (1.0+teqps<0.){pow_1teqps_CM1 = +(pow(-1.0-teqps,constC-1.0));}
        
        
        *yieldStress =
        (constA+constB*pow_eqps_n)*
        pow_1teqps_C*
        (1-pow_hmlgT_M);
        
        yieldStress_lambda = 
        pow_1teqps_CM1*
        (1-pow_hmlgT_M)*
        ( +(constB*constN*pow_eqps_nM1)*(1+teqps)
        +(constA+constB*pow_eqps_n)*constC*teqps_lambda   );
        
        yieldFunction = (vmStressTrial - 3.*shearModulus*lambda - (*yieldStress))/constA; // adimensional
        yieldFunction_lambda = (-3*shearModulus-yieldStress_lambda)/constA;
        //std::cout << "it=" << it <<"   yieldFunction=" << yieldFunction << "   lambda=" << lambda << "\n";
        if (it==max_it){
            std::cout << "WARNING: NOT-CONVERGED PLASTIC STRAIN LOOP:" <<  "   yieldFunction=" << yieldFunction << "   lambda=" << lambda << "   trialstress=" << vmStressTrial << std::endl;
            yieldFunction=0.0;
            *eqpsNP1= *eqpsN; // equivalent plastic strain
//             std::cout << "  vmStressTrial=" << vmStressTrial ;
//             std::cout << "\n";
        }
    }
    if (lambda<0.0) {
        std::cout << "WARNING: PLASTIC STRAIN LOOP CONVERGED TO A NEGATIVE PLASTIC STRAIN" << std::endl << "         setting lambda to zero." << std::endl ;
        *eqpsNP1= *eqpsN; // equivalent plastic strain
    }
}

// Explicit template instantiation for double
template void JohnsonCookSolve<double>
(
const double vmStressTrial,
const double* eqpsN,
double* eqpsNP1,
double* yieldStress,
const double shearModulus,
const double constA,
const double constN,
const double constB,
const double constC,
const double pow_hmlgT_M,
const double dt
);

/* Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void JohnsonCookSolve<Sacado::Fad::DFad<double>>
(
const Sacado::Fad::DFad<double> vmStressTrial,
const Sacado::Fad::DFad<double>* eqpsN,
Sacado::Fad::DFad<double>* eqpsNP1,
Sacado::Fad::DFad<double>* yieldStress,
const double shearModulus,
const double constA,
const double constN,
const double constB,
const double constC,
const double pow_hmlgT_M,
const double dt
);

}
