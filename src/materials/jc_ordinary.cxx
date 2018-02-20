//! \file jc_ordinary.cxx

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

#include "jc_ordinary.h"
#include "Peridigm_Material.hpp"
#include <cmath>
#include <Sacado.hpp>
#include "jc_ordinary.h"
#include "material_utilities.h"

#include <math.h>
#include <iostream>


namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeInternalForceJohnsonCookOrdinary
(
const double* xOverlap,
const ScalarT* yOverlap,
const double* mOwned,
const double* volumeOverlap,
const ScalarT* dilatationOwned,
const double* bondDamage,
const double* scfOwned,
ScalarT* fInternalOverlap,
const int*  localNeighborList,
int numOwnedPoints,
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
double horizon,
ScalarT* ElasticEnergyDensity,
const ScalarT* microPotentialN,
ScalarT* microPotentialNP1,
ScalarT* VonMisesStress,
const ScalarT* deviatoricPlasticExtensionN,
ScalarT* deviatoricPlasticExtensionNP1,
const ScalarT* EquivalentPlasticStrainN,
ScalarT* EquivalentPlasticStrainNP1,
const ScalarT* AccumulatedPlasticStrainN,
ScalarT* AccumulatedPlasticStrainNP1,
const ScalarT* LocalDamageN,
ScalarT* LocalDamageNP1,
const double* deltaTemperature,
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
const double constDC
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	double K = obj_bulkModulus.compute(0.0);
	double MU = obj_shearModulus.compute(0.0);
    double thermalExpansionCoefficient;
    if(deltaTemperature) thermalExpansionCoefficient = obj_alphaVol.compute(0.0);

	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
    const double *deltaT = deltaTemperature;
	const double *m = mOwned;
	const double *v = volumeOverlap;
    const double *scf = scfOwned;
	const ScalarT *theta = dilatationOwned;
	ScalarT *fOwned = fInternalOverlap;
    
    const ScalarT *muW_N = microPotentialN;
    ScalarT *muW_NP1 = microPotentialNP1;
    ScalarT *W   = ElasticEnergyDensity;
    ScalarT Wd_hat;
    ScalarT Wi_hat;
    ScalarT W0_hat;
    ScalarT F_hat;

    ScalarT *VMstress = VonMisesStress;
    const ScalarT *edpN = deviatoricPlasticExtensionN;
    ScalarT *edpNP1 = deviatoricPlasticExtensionNP1;
    const ScalarT *eqpsN = EquivalentPlasticStrainN;
    ScalarT *eqpsNP1 = EquivalentPlasticStrainNP1;
    const ScalarT *dapsN = AccumulatedPlasticStrainN;
    ScalarT *dapsNP1 = AccumulatedPlasticStrainNP1;
    const ScalarT *DaN = LocalDamageN;
    ScalarT *DaNP1 = LocalDamageNP1;
    
    double hmlgT;
    
    ScalarT pow_eqps_n;
    ScalarT pow_hmlgT_M;
    ScalarT yieldStress0;


//     ScalarT Deqps;
//     ScalarT teqps=0;
//     ScalarT teqps_Deqps;
// 
//     ScalarT pow_eqps_nM1;
//     ScalarT pow_1teqps_C;
//     ScalarT pow_1teqps_Cm1;
//     ScalarT yieldStressHat_Deqps;
//     ScalarT fun1;
//     ScalarT fun1_Deqps;
// 
//     ScalarT tdaps;
//     ScalarT pow_1tdaps_D4;
//     ScalarT pow_1tdaps_D4M1;
//     ScalarT tdaps_Da;
//     ScalarT hydroStress;


	const int *neighPtr = localNeighborList;
    const int *neighPtrOverLap = localNeighborList;
    const double *bondDamageOverLap = bondDamage;
    const ScalarT *edpN_OverLap = edpN;
    ScalarT *edpNP1_OverLap = edpNP1;
	double cellVolume, alpha, X_dx, X_dy, X_dz, zeta, omega;
	ScalarT Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz, e/*, c1*/;
    ScalarT ed;
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, W++, VMstress++, eqpsN++, eqpsNP1++, dapsN++, dapsNP1++, DaN++, DaNP1++, deltaT++, m++, theta++, scf++){
        
        if(deltaTemperature){
            K    = obj_bulkModulus.compute(*deltaT);
            MU   = obj_shearModulus.compute(*deltaT);
            thermalExpansionCoefficient = obj_alphaVol.compute(*deltaT);
        }

        
		int numNeigh = *neighPtr; neighPtr++; neighPtrOverLap++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		alpha = 15.0*MU/(*m);
		alpha *= (*scf);
		double selfCellVolume = v[p];
                
        Wd_hat=0.0;
		for(int n=0;n<numNeigh;n++,neighPtrOverLap++,bondDamageOverLap++,edpN_OverLap++,edpNP1_OverLap++){
            *edpNP1_OverLap=*edpN_OverLap;

			int localId = *neighPtrOverLap;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yOverlap[3*localId];
			X_dx = XP[0]-X[0];
			X_dy = XP[1]-X[1];
			X_dz = XP[2]-X[2];
			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
            e = dY - zeta;
            if(deltaTemperature){
//                 std::cout << *deltaT << std::endl;
              e -= thermalExpansionCoefficient*(*deltaT)*zeta;
            }
			omega = scalarInfluenceFunction(zeta,horizon);
            ed = e - *theta/3*zeta - *edpNP1_OverLap; // deviatoric Extension
            
            // compute deviatoric energy density
            Wd_hat += (1.0-*bondDamageOverLap)* alpha/2 * ed * omega * ed * cellVolume;
		}
		
		*VMstress = sqrt(6*MU*Wd_hat);
        
        if(deltaTemperature) hmlgT = (*deltaT - ReferenceTemperature) / (MeltingTemperature - ReferenceTemperature) ; // Homologous Temperature
        else hmlgT=0.0;
        
        // update strains and damage
        *eqpsNP1 = *eqpsN;
        *dapsNP1 = *dapsN;
        *DaNP1 = *DaN;
        
        if ((hmlgT<0.) && (constM<1.)) pow_hmlgT_M=hmlgT;
        else pow_hmlgT_M=pow(hmlgT,constM);
        
        pow_eqps_n = 0.;
        if (*eqpsN>0.) pow_eqps_n = +(pow(*eqpsN,constN));

        
        yieldStress0 = // without considering damage, and teqps=0
          (constA+constB*pow_eqps_n)*
          (1-pow_hmlgT_M);
        W0_hat = // yielding energy density
          pow(yieldStress0,2.0)/(6*MU);
          
        F_hat = Wd_hat - W0_hat ;

        if (F_hat > 0.0){
            // std::cout << "YIELDING: " << yieldStress0 << " " << ed << " " << *VMstress <<  std::endl;
        }
        else{ // step is elastic
        }


        Wd_hat=0.0;
        for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,edpN++,edpNP1++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yOverlap[3*localId];
			X_dx = XP[0]-X[0];
			X_dy = XP[1]-X[1];
			X_dz = XP[2]-X[2];
			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
            e = dY - zeta;
            if(deltaTemperature)
              e -= thermalExpansionCoefficient*(*deltaT)*zeta;
			omega = scalarInfluenceFunction(zeta,horizon);
            ed = e - *theta/3*zeta - *edpNP1; // deviatoric Extension
            t = (1.0-*bondDamage)*(3.0*K*(*theta)/(*m)* omega * zeta + omega * alpha * ed);
			// c1 = omega*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
// 			c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0);
// 			t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);
            
			fx = t * Y_dx / dY;
			fy = t * Y_dy / dY;
			fz = t * Y_dz / dY;

			*(fOwned+0) += fx*cellVolume;
			*(fOwned+1) += fy*cellVolume;
			*(fOwned+2) += fz*cellVolume;
			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;

            // compute deviatoric energy density
            Wd_hat += (1.0-*bondDamage)* alpha/2 * ed * omega * ed * cellVolume;
        }
        
        Wi_hat = K/2 * (*theta) * (*theta) ;
		*W = (1.0-*DaNP1)*(Wi_hat + Wd_hat) ;

	}
}

/** Explicit template instantiation for double. */
template void computeInternalForceJohnsonCookOrdinary<double>
(
const double* xOverlap,
const double* yOverlap,
const double* mOwned,
const double* volumeOverlap,
const double* dilatationOwned,
const double* bondDamage,
const double* scfOwned,
double* fInternalOverlap,
const int*  localNeighborList,
int numOwnedPoints,
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
double horizon,
double* ElasticEnergyDensity,
const double* microPotentialN,
double* microPotentialNP1,
double* VonMisesStress,
const double* deviatoricPlasticExtensionN,
double* deviatoricPlasticExtensionNP1,
const double* EquivalentPlasticStrainN,
double* EquivalentPlasticStrainNP1,
const double* AccumulatedPlasticStrainN,
double* AccumulatedPlasticStrainNP1,
const double* LocalDamageN,
double* LocalDamageNP1,
const double* deltaTemperature,
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
const double constDC
 );

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceJohnsonCookOrdinary<Sacado::Fad::DFad<double> >
(
const double* xOverlap,
const Sacado::Fad::DFad<double>* yOverlap,
const double* mOwned,
const double* volumeOverlap,
const Sacado::Fad::DFad<double>* dilatationOwned,
const double* bondDamage,
const double* scfOwned,
Sacado::Fad::DFad<double>* fInternalOverlap,
const int*  localNeighborList,
int numOwnedPoints,
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
double horizon,
Sacado::Fad::DFad<double>* ElasticEnergyDensity,
const Sacado::Fad::DFad<double>* microPotentialN,
Sacado::Fad::DFad<double>* microPotentialNP1,
Sacado::Fad::DFad<double>* VonMisesStress,
const Sacado::Fad::DFad<double>* deviatoricPlasticExtensionN,
Sacado::Fad::DFad<double>* deviatoricPlasticExtensionNP1,
const Sacado::Fad::DFad<double>* EquivalentPlasticStrainN,
Sacado::Fad::DFad<double>* EquivalentPlasticStrainNP1,
const Sacado::Fad::DFad<double>* AccumulatedPlasticStrainN,
Sacado::Fad::DFad<double>* AccumulatedPlasticStrainNP1,
const Sacado::Fad::DFad<double>* LocalDamageN,
Sacado::Fad::DFad<double>* LocalDamageNP1,
const double* deltaTemperature,
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
const double constDC
);

}
