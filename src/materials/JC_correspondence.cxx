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
#include "Peridigm_Material.hpp"
#include <Sacado.hpp>
#include <math.h>
#include <iostream>


namespace CORRESPONDENCE {

    
template<typename ScalarT>
void updateJohnsonCookCauchyStress
(
const ScalarT* unrotatedRateOfDeformation, 
const ScalarT* cauchyStressN, 
ScalarT* cauchyStressNP1, 
ScalarT* vonMisesStressNP1,
const ScalarT* equivalentPlasticStrainN,
const ScalarT* accumulatedPlasticStrainN,
const ScalarT* DamageN,
ScalarT* equivalentPlasticStrainNP1,
ScalarT* accumulatedPlasticStrainNP1,
ScalarT* DamageNP1,
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
const double constM,
const double constD1,
const double constD2,
const double constD3,
const double constD4,
const double constD5,
const double constDC,
ScalarT* DissipationNP1
)
{
    const ScalarT* rateOfDef = unrotatedRateOfDeformation;
    const ScalarT* stressN = cauchyStressN;
    ScalarT* stressNP1 = cauchyStressNP1;
    
    ScalarT* vmStress = vonMisesStressNP1;
    
    const ScalarT* eqpsN = equivalentPlasticStrainN;
    const ScalarT* dapsN = accumulatedPlasticStrainN;
    const ScalarT* DaN = DamageN;
    
    
    ScalarT hydroStressN;
    ScalarT deviatoricStressN[9];
    
    ScalarT* eqpsNP1 = equivalentPlasticStrainNP1;
    ScalarT* dapsNP1 = accumulatedPlasticStrainNP1;
    ScalarT* DaNP1 = DamageNP1;
    
    
    ScalarT strainInc[9];
    ScalarT deviatoricStrainInc[9];
    
    ScalarT deviatoricStressNP1[9];
    
    ScalarT dilatationInc;
    ScalarT hydroStressNP1;
    ScalarT tempScalar;
    
    ScalarT yieldStressHat0;
    
    
    int it;
    
    ScalarT Deqps;
    ScalarT teqps=0;
    ScalarT teqps_Deqps;
    ScalarT pow_eqps_n;
    ScalarT pow_eqps_nM1;
    ScalarT pow_1teqps_C;
    ScalarT pow_1teqps_Cm1;
    ScalarT yieldStressHat;
    ScalarT yieldStressHat_Deqps;
    ScalarT fun1;
    ScalarT fun1_Deqps;
    
    ScalarT hmlgT;
    ScalarT pow_hmlgT_M;
    ScalarT tdaps;
    ScalarT pow_1tdaps_D4;
    ScalarT pow_1tdaps_D4M1;
    ScalarT tdaps_Da;
    ScalarT hydroStress;
    ScalarT yieldStress;
    ScalarT frs;
    ScalarT frs_Da;
    ScalarT fun2;
    ScalarT fun2_Da;
    
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
        ++eqpsN,   ++eqpsNP1,   ++dapsN,   ++dapsNP1,   ++DaN,   ++DaNP1,
        ++deltaTemperatureN,    ++deltaTemperatureNP1, ++DissipationNP1
        ){
        
        // temperatures
        hmlgT = (*deltaTemperatureNP1 - ReferenceTemperature) / (MeltingTemperature - ReferenceTemperature) ; // Homologous Temperature

        //Tempdouble = *Temperature;
        bulkModN    =obj_bulkModulus.compute(*deltaTemperatureN);
        shearModN   =obj_shearModulus.compute(*deltaTemperatureN);
        alphaN      =obj_alphaVol.compute(*deltaTemperatureN);
        shearModNP1 =obj_shearModulus.compute(*deltaTemperatureNP1);
        bulkModNP1  =obj_bulkModulus.compute(*deltaTemperatureNP1);
        alphaNP1    =obj_alphaVol.compute(*deltaTemperatureNP1);

        
        if (*DaN==1.){
//             std::cout << iID << " Entirely damaged\n";
            *DaNP1=*DaN;
            *eqpsNP1=*eqpsN;
            *dapsNP1=*dapsN;
            *vmStress=0.;
            for (int i = 0; i < 9; i++) {
                stressNP1[i] = 0.;
            }
            continue;
        }

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
            *(deviatoricStressNP1+i) = *(deviatoricStressN+i)*shearModNP1/shearModN + (1.-*DaN)*deviatoricStrainInc[i]*2.0*shearModNP1;
        }
        hydroStressNP1=hydroStressN*bulkModNP1/bulkModN+(1.-*DaN)*bulkModNP1*dilatationInc;
        
        
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

        *vmStress = sqrt(3.0/2.0*tempScalar);
        
        // update strains and damage
        Deqps=0.0;
        *eqpsNP1 = *eqpsN;
        *dapsNP1 = *dapsN;
        *DaNP1 = *DaN;
        
        
        if ((hmlgT<0.) && (constM<1.)){
            pow_hmlgT_M=hmlgT;
            
        }
        else{
            pow_hmlgT_M=pow(hmlgT,constM);
        }
        
        yieldStressHat0 = // without considering damage
                (constA+constB*pow(*eqpsN,constN))*
                (1-pow_hmlgT_M);
        //std::cout << "\n\n\n" << *vmStress/(1.-*DaN) << ">" << yieldStressHat0 << "\n";

        //If true, the step is plastic and we need to return to the yield
        //surface.
        if( *vmStress/(1.-*DaN) >= yieldStressHat0 ) {
            
            //std::cout << "\n" << *vmStress/(1.-*DaN) << ">" << yieldStressHat0 << "\n";
            //std::cout << "DaN: " << *DaN << "   dapsN:" << *dapsN << "\n";
            // FIND PLASTIC STRAIN AND DAMAGE USING NEWTON'S METHOD
            
            // eqps : equivalent plastic strain
            // Deqps : increment of equivalent plastic strain
            // daps : damage accumulated plastic strain
            // Da   : damage
            // hmlgT : homologous temperature
                        
            
            // EXPRESS fun1 function of UNKNOWN eqps
            // DERIVE fun1 wrt deqps to obtain the Jacobian for the Newton's method
            
            fun1=1;
            it=0;
            while (std::abs(fun1) > 1e-7) {
                if (it>0)
                {
                    // x_it   =  -inv(M) * f + x_old
                    Deqps   = -fun1/fun1_Deqps + Deqps;
                }
                
                ++it;
                *eqpsNP1= *eqpsN + Deqps; // equivalent plastic strain
                teqps = Deqps/dt; // time derived eqps
                teqps_Deqps= 1./dt;
            
                pow_eqps_nM1 = 0.;
                if (*eqpsNP1>0.){pow_eqps_nM1 = +(pow(*eqpsNP1,constN-1));}
                
                pow_eqps_n = 0.;
                if (*eqpsNP1>0.){pow_eqps_n = +(pow(*eqpsNP1,constN));}
                
                pow_1teqps_C = 0.;
                if (1+teqps>0.){pow_1teqps_C = +(pow(1+teqps,constC));}
                pow_1teqps_Cm1 = 0.;
                if (1+teqps>0.){pow_1teqps_Cm1 = +(pow(1+teqps,constC-1));}
                
                
                yieldStressHat =
                (constA+constB*pow_eqps_n)*
                pow_1teqps_C*
                (1-pow_hmlgT_M);
                
                yieldStressHat_Deqps = 
                pow_1teqps_Cm1*
                (1-pow_hmlgT_M)*
                ( +(constB*constN*pow_eqps_nM1)*(1+teqps)
                +(constA+constB*pow_eqps_n)*constC*teqps_Deqps   );
                
                fun1 = ((*vmStress)/(1-*DaN) - 3.*shearModNP1*Deqps - yieldStressHat)/constA; // adimensional
                fun1_Deqps = (-3*shearModNP1-yieldStressHat_Deqps)/constA;
                //std::cout << "it=" << it <<"   fun1=" << fun1 << "   Deqps=" << Deqps << "\n";
                if (it==20){fun1=0;
                    std::cout << "WARNING: NOT-CONVERGED PLASTIC STRAIN LOOP:" <<  "   fun1=" << fun1 << "   Deqps=" << Deqps <<  "   iID=" << iID+1;
                    std::cout << "  vmStress=" << *vmStress ;
                    std::cout << "\n";
                    for (int i = 0; i < 9; i++) {
                        std::cout << *(stressN+i) << "  "<< *(stressNP1+i) << "  " << *(rateOfDef+i) << "  ";
                    }
                    std::cout << "\n";
                }
            }
            
            if (Deqps>=0.) {
                // EXPRESS fun2 function of UNKNOWN Da
                // DERIVE fun2 wrt Da to obtain the Jacobian for the Newton's method
                fun2=1;
                it=0;
                while (std::abs(fun2) > 1e-7) {
                    if (it>0)
                    {
                        // x_it   =  -inv(M) * f + x_old
                        *DaNP1   = -fun2/fun2_Da + *DaNP1;
                    }
                    ++it;
                    tdaps = teqps/(1-*DaNP1); // time derived daps
                    tdaps_Da = 1/pow(1-*DaNP1,2)*teqps;
                    *dapsNP1 = tdaps*dt+*dapsN; // damage accumulated plastic strain
                    
                    yieldStress = yieldStressHat*(1.-*DaNP1);
                    
                    hydroStress = hydroStressNP1*
                    (1.-*DaNP1)/(1.-*DaN);
                    
                    if (tdaps>0.){
                        pow_1tdaps_D4=pow(1+tdaps,constD4);
                        pow_1tdaps_D4M1=pow(1+tdaps,constD4-1);
                    }
                    else{
                        pow_1tdaps_D4=1;
                        pow_1tdaps_D4M1=1;
                    }
                            
                    frs = (constD1+constD2*exp(constD3*hydroStress/yieldStress))*pow_1tdaps_D4*(1+constD5*hmlgT); // fracture strain
                    frs_Da = (constD1+constD2*exp(constD3*hydroStress/yieldStress))*constD4*pow_1tdaps_D4M1*tdaps_Da*(1+constD5*hmlgT);
                    
//                     std::cout << frs << "  " << (constD1+constD2*exp(constD3*hydroStress/yieldStress)) << "  " << hydroStress/yieldStress << "  " << "\n";
                    
                    fun2 = *DaNP1-*DaN -constDC/(frs*(1-*DaNP1))*Deqps ;
                    fun2_Da = 1 + constDC/pow((1-*DaNP1)*frs,2)*(-frs+(1-*DaNP1)*frs_Da)*Deqps;
                    if(*vmStress/(1.-*DaN) > 5*yieldStressHat0){
                    //std::cout << "it=" << it <<"   fun2=" << fun2 << "   Da=" << *DaNP1 << "\n";
                    }
                    if (it==20){*DaNP1=1;fun2=0.;
                        std::cout << "WARNING: NOT-CONVERGED DAMAGE LOOP:" <<  "   Imposed damage=" << *DaNP1 << "   iID=" << iID+1;
                        std::cout << "   fun2=" << fun2 <<"  vmStress=" << *vmStress ;
                        std::cout << "\n";
                    for (int i = 0; i < 9; i++) {
                        std::cout << *(stressN+i) << "  "<< *(stressNP1+i) << "  " << *(rateOfDef+i) << "  ";
                    }
                    std::cout << "\n";
                        
                    }
            
                }
                if (*DaNP1>1.){*DaNP1=1.;}
                if (*DaNP1>=*DaN){
                    // RADIAL RETURN
                    // update deviatoric stress
                    // vmStress = yieldStress in the new condition
                    tempScalar = yieldStressHat*(1.-*DaNP1)/(*vmStress);
                    for (int i = 0; i < 9; i++) {
                        deviatoricStressNP1[i] *= tempScalar; 
                        *(stressNP1+i) = deviatoricStressNP1[i];
                    }
                    // update hydrostatic stress
                    hydroStressNP1*=(1-*DaNP1)/(1-*DaN);
                    *stressNP1 += hydroStressNP1;
                    *(stressNP1+4) += hydroStressNP1;
                    *(stressNP1+8) += hydroStressNP1;
                     
                    // Check if same Von Mises stress obtained
                    // Update the von Mises stress now that the state of stress is on the
                    // yield surface
                    tempScalar = 0.0;
                    for (int j = 0; j < 3; j++) {
                        for (int i = 0; i < 3; i++) {
                            tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
                        }
                    }
                    //std::cout << "Von Mises Stress: " << pow(1.5*tempScalar,.5) << "   " << yieldStress;
                    
                    *vmStress = yieldStress;
                }
                else
                {
                    std::cout << "WARNING: Negative Delta Damage:" << "   DaNP1:  " << *DaNP1 << "   DaN:  " << *DaN << "\n";
                    *DaNP1=*DaN;
                    *dapsNP1=*dapsN;
                } // end if damage
            }
            else{
                std::cout << "WARNING: Negative delta plastic epsilon\n" <<  "  Delta plastic strain:  " << Deqps << "\n";
                Deqps=0;
                *eqpsNP1= *eqpsN + Deqps;
            } // end if Deqps
            
        } else {
            yieldStress = (1.-*DaNP1)*yieldStressHat0;
            // The step is elastic
        }; // end if yield
        
        
        
        
//         // compute DissipationNP1 (Damage conjugate)
//         YoungModNP1 = 9*bulkModNP1*shearModNP1/(3*bulkModNP1+shearModNP1);
//         PoissRatNP1 = (3*bulkModNP1-2*shearModNP1)/(2*(3*bulkModNP1+shearModNP1));
//         if (*vmStress==0) *DissipationNP1=0.0;
//         else *DissipationNP1 = pow(*vmStress,2) / (2*YoungModNP1*pow(1-*DaNP1,2)) *( 2/3*(1.+PoissRatNP1) +3*(1-2*PoissRatNP1)*pow(hydroStressNP1/(*vmStress),2));
//         //cout << hydroStressNP1 << "  " << *vmStress << endl;
        
        
        
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
const double* accumulatedPlasticStrainN,
const double* DamageN,
double* equivalentPlasticStrainNP1,
double* accumulatedPlasticStrainNP1,
double* DamageNP1,
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
const double constM,
const double constD1,
const double constD2,
const double constD3,
const double constD4,
const double constD5,
const double constDC,
double* DissipationNP1
);

/* Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void updateJohnsonCookCauchyStress<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformation, 
const Sacado::Fad::DFad<double>* cauchyStressN, 
Sacado::Fad::DFad<double>* cauchyStressNP1, 
Sacado::Fad::DFad<double>* vonMisesStressNP1,
const Sacado::Fad::DFad<double>* equivalentPlasticStrainN,
const Sacado::Fad::DFad<double>* accumulatedPlasticStrainN,
const Sacado::Fad::DFad<double>* DamageN,
Sacado::Fad::DFad<double>* equivalentPlasticStrainNP1,
Sacado::Fad::DFad<double>* accumulatedPlasticStrainNP1,
Sacado::Fad::DFad<double>* DamageNP1,
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
const double constM,
const double constD1,
const double constD2,
const double constD3,
const double constD4,
const double constD5,
const double constDC,
Sacado::Fad::DFad<double>* DissipationNP1
);

}
