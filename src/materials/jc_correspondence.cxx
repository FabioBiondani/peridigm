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

#include <Teuchos_Assert.hpp>
#include "jc_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include "Peridigm_Material.hpp"
#include <Sacado.hpp>
#include <math.h>
#include <iostream>
#include "JohnsonCook.h"


namespace CORRESPONDENCE {


    
template<typename ScalarT>
void updateJohnsonCookCauchyStress
(
const ScalarT* unrotatedRateOfDeformation, 
const ScalarT* cauchyStressN, 
ScalarT* cauchyStressNP1, 
ScalarT* vonMisesStressNP1,
const ScalarT* equivalentPlasticStrainN,
ScalarT* equivalentPlasticStrainNP1,
const int numPoints, 
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
const double* deltaTemperatureN,
const double* deltaTemperatureNP1,
const double dt,
const double MeltingTemperature,
const double ReferenceTemperature,
const double constA,
const double constN,
const double constB,
const double constC,
const double constM
)
{
    const ScalarT* rateOfDef = unrotatedRateOfDeformation;
    const ScalarT* stressN = cauchyStressN;
    ScalarT* stressNP1 = cauchyStressNP1;
    
    ScalarT* vmStress = vonMisesStressNP1;
    
    const ScalarT* eqpsN = equivalentPlasticStrainN;
    
    ScalarT hydroStressN;
    ScalarT deviatoricStressN[9];
    
    ScalarT* eqpsNP1 = equivalentPlasticStrainNP1;
    
    
    ScalarT strainInc[9];
    ScalarT deviatoricStrainInc[9];
    
    ScalarT deviatoricStressNP1[9];
    
    ScalarT dilatationInc;
    ScalarT hydroStressNP1;
    ScalarT tempScalar;
    
    ScalarT vmStressTrial;
    
    ScalarT lambda;
    ScalarT yieldStress;
    
    ScalarT pow_eqps_n;
    double hmlgT;
    double pow_hmlgT_M;
    
    double bulkModN;
    double shearModN;
    double alphaN;
    double bulkModNP1;
    double shearModNP1;
    double alphaNP1;
//     double YoungModNP1;
//     double PoissRatNP1;
    
    double ThermalExpansionStrain;
    
    for(int iID=0 ; iID<numPoints ; ++iID, rateOfDef+=9, stressN+=9, stressNP1+=9, ++vmStress,
        ++eqpsN,   ++eqpsNP1,
        ++deltaTemperatureN,    ++deltaTemperatureNP1
        ){
        
        // temperatures
        hmlgT = (*deltaTemperatureNP1 - ReferenceTemperature) / (MeltingTemperature - ReferenceTemperature) ; // Homologous Temperature

        //Tempdouble = *Temperature;
        bulkModN    =obj_bulkModulus.compute(*deltaTemperatureN);
        shearModN   =obj_shearModulus.compute(*deltaTemperatureN);
        alphaN      =obj_alphaVol.compute(*deltaTemperatureN);
        bulkModNP1  =obj_bulkModulus.compute(*deltaTemperatureNP1);
        shearModNP1 =obj_shearModulus.compute(*deltaTemperatureNP1);
        alphaNP1    =obj_alphaVol.compute(*deltaTemperatureNP1);

        //strainInc = dt * rateOfDef
        for (int i = 0; i < 9; i++) {
            strainInc[i] = *(rateOfDef+i)*dt;
        }
        
        // Thermal isovolumetric expansion
        ThermalExpansionStrain = (alphaNP1+alphaN)/2* (*deltaTemperatureNP1-(*deltaTemperatureN));
        strainInc[0] -= ThermalExpansionStrain;
        strainInc[4] -= ThermalExpansionStrain;
        strainInc[8] -= ThermalExpansionStrain;

        for (int i = 0; i < 9; i++) {
            deviatoricStrainInc[i] = strainInc[i];
        }
        
        //dilatation
        dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

        //deviatoric strain
        deviatoricStrainInc[0] -= dilatationInc/3.0;
        deviatoricStrainInc[4] -= dilatationInc/3.0;
        deviatoricStrainInc[8] -= dilatationInc/3.0;

        
        
        hydroStressN=(*(stressN) + *(stressN+4) + *(stressN+8))/3.0;
        for (int i = 0; i < 9; i++) {
            *(deviatoricStressN+i) = *(stressN+i);
        }
        *(deviatoricStressN) -= hydroStressN;
        *(deviatoricStressN+4) -= hydroStressN;
        *(deviatoricStressN+8) -= hydroStressN;
        
        
        //Compute an elastic ``trial stress''
        for (int i = 0; i < 9; i++) {
            *(deviatoricStressNP1+i) = *(deviatoricStressN+i)*shearModNP1/shearModN + deviatoricStrainInc[i]*2.0*shearModNP1;
        }
        hydroStressNP1=hydroStressN*bulkModNP1/bulkModN+bulkModNP1*dilatationInc;
        
        
        for (int i = 0; i < 9; i++) {
            *(stressNP1+i) = *(deviatoricStressNP1+i);
        }
        *(stressNP1) += hydroStressNP1;
        *(stressNP1+4) += hydroStressNP1;
        *(stressNP1+8) += hydroStressNP1;


        // Compute \sigma_ij * \sigma_ij
        tempScalar = 0.0;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
                tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
            }
        }

        vmStressTrial = sqrt(3.0/2.0*tempScalar);
        
        if ((hmlgT<0.) && (constM<1.)) pow_hmlgT_M=hmlgT;
        else pow_hmlgT_M=pow(hmlgT,constM);
        if (*eqpsN>0.) pow_eqps_n = +(pow(*eqpsN,constN));
        else pow_eqps_n = 0.;
        
        yieldStress = // actual yield stress if step is elastic
                (constA+constB*pow_eqps_n)*
                (1-pow_hmlgT_M);

        //If true, the step is plastic and we need to return to the yield
        //surface.
        if( vmStressTrial - yieldStress >= 0 ) {
            
            MATERIAL_EVALUATION::JohnsonCookSolve(
                vmStressTrial,
                eqpsN,
                eqpsNP1,
                &yieldStress,
                shearModNP1,
                constA,
                constN,
                constB,
                constC,
                pow_hmlgT_M,
                dt
            );
            
            lambda=*eqpsNP1-*eqpsN;
            
            if (lambda>=0.) {
                // RADIAL RETURN
                // update deviatoric stress
                // vmStress = yieldStress in the new condition
                *vmStress = yieldStress;
                tempScalar = *vmStress/vmStressTrial;
                for (int i = 0; i < 9; i++) {
                    deviatoricStressNP1[i] *= tempScalar; 
                    *(stressNP1+i) = deviatoricStressNP1[i];
                }
                *stressNP1 += hydroStressNP1;
                *(stressNP1+4) += hydroStressNP1;
                *(stressNP1+8) += hydroStressNP1;

//                 tempScalar = 0.0;  for (int j = 0; j < 3; j++) { for (int i = 0; i < 3; i++) { tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j]; } }
//                 std::cout << "Von Mises Stress: " << pow(1.5*tempScalar,.5) << "   " << yieldStress << std::endl;

            }
            else{
                std::cout << "ERROR: Negative delta plastic epsilon after loop\n" <<  "  Delta plastic strain:  " << lambda << "\n";
                TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "")
            } // end if lambda
        } else { // The step is elastic
            *eqpsNP1= *eqpsN;
            *vmStress = vmStressTrial;
        }; // end if yield
    } // end for points
}//end updateJohnsonCookCauchyStress

// Explicit template instantiation for double
template void updateJohnsonCookCauchyStress<double>
(
const double* unrotatedRateOfDeformation, 
const double* cauchyStressN, 
double* cauchyStressNP1,
double* vonMisesStressNP1,
const double* equivalentPlasticStrainN,
double* equivalentPlasticStrainNP1,
const int numPoints, 
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
const double* deltaTemperatureN,
const double* deltaTemperatureNP1,
const double dt,
const double MeltingTemperature,
const double ReferenceTemperature,
const double constA,
const double constN,
const double constB,
const double constC,
const double constM
);

/* Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void updateJohnsonCookCauchyStress<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformation, 
const Sacado::Fad::DFad<double>* cauchyStressN, 
Sacado::Fad::DFad<double>* cauchyStressNP1, 
Sacado::Fad::DFad<double>* vonMisesStressNP1,
const Sacado::Fad::DFad<double>* equivalentPlasticStrainN,
Sacado::Fad::DFad<double>* equivalentPlasticStrainNP1,
const int numPoints, 
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
const double* deltaTemperatureN,
const double* deltaTemperatureNP1,
const double dt,
const double MeltingTemperature,
const double ReferenceTemperature,
const double constA,
const double constN,
const double constB,
const double constC,
const double constM
);

}




