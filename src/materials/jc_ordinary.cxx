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

#include <Teuchos_Assert.hpp>
#include "jc_ordinary.h"
#include "Peridigm_Material.hpp"
#include <cmath>
#include <Sacado.hpp>
#include "jc_ordinary.h"
#include "material_utilities.h"
#include <math.h>
#include <iostream>
#include "JohnsonCook.h"


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
ScalarT* VonMisesStress,
const ScalarT* deviatoricPlasticExtensionN,
ScalarT* deviatoricPlasticExtensionNP1,
const ScalarT* EquivalentPlasticStrainN,
ScalarT* EquivalentPlasticStrainNP1,
ScalarT* deviatoricForceDensity,
const double* deltaTemperature,
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
const double horizon,
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
    
    ScalarT Wd;
    
    ScalarT *vmStress = VonMisesStress;
    const ScalarT *edpN = deviatoricPlasticExtensionN;
    ScalarT *edpNP1 = deviatoricPlasticExtensionNP1;
    const ScalarT *eqpsN = EquivalentPlasticStrainN;
    ScalarT *eqpsNP1 = EquivalentPlasticStrainNP1;
    ScalarT *td = deviatoricForceDensity;
    
    double hmlgT;
    
    ScalarT vmStressTrial;
    
    ScalarT pow_eqps_n;
    double pow_hmlgT_M;
    
    ScalarT yieldStress;
    
    ScalarT lambda;
    ScalarT tempScalar;

    const int *neighPtr = localNeighborList;
    double cellVolume, alpha, X_dx, X_dy, X_dz, zeta, omega;
    ScalarT Y_dx, Y_dy, Y_dz, dY, t, ti, fx, fy, fz, e;
    ScalarT ed;

    const int *neighPtrOverLap = localNeighborList;
    const ScalarT *edpN_OverLap = deviatoricPlasticExtensionN;
    ScalarT *edpNP1_OverLap = deviatoricPlasticExtensionNP1;
    ScalarT *td_OverLap = deviatoricForceDensity;
    const double *bondDamageOverLap = bondDamage;


    for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, vmStress++, eqpsN++, eqpsNP1++, deltaT++, m++, theta++, scf++){
        
        if(deltaTemperature){
            hmlgT = (*deltaT - ReferenceTemperature) / (MeltingTemperature - ReferenceTemperature) ; // Homologous Temperature
            K    = obj_bulkModulus.compute(*deltaT);
            MU   = obj_shearModulus.compute(*deltaT);
            thermalExpansionCoefficient = obj_alphaVol.compute(*deltaT);
        }
        else hmlgT=0.0;

		int numNeigh = *neighPtr; neighPtr++; neighPtrOverLap++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		alpha = 15.0*MU/(*m);
		alpha *= (*scf);
		double selfCellVolume = v[p];

        Wd=0.0;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,edpN++,edpNP1++,td++){
            *edpNP1=*edpN;
            

			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yOverlap[3*localId];
            X_dx = XP[0]-X[0];  X_dy = XP[1]-X[1];  X_dz = XP[2]-X[2];
            zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
            Y_dx = YP[0]-Y[0];  Y_dy = YP[1]-Y[1];  Y_dz = YP[2]-Y[2];
            dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
            e = dY - zeta;
            if(deltaTemperature) e -= thermalExpansionCoefficient*(*deltaT)*zeta;
			omega = scalarInfluenceFunction(zeta,horizon);
            ed = e - *theta/3*zeta - *edpNP1; // deviatoric Extension

            *td = (1.0-*bondDamage)*(omega * alpha * ed);

            if (*edpN/zeta>0.2) std::cout << *edpN << "  " << *edpNP1 << std::endl;

            // compute deviatoric energy density
            Wd += (1.0-*bondDamage)* alpha/2 * ed * omega * ed * cellVolume;
		}

		vmStressTrial = sqrt(6*MU*Wd);
        
        if ((hmlgT<0.) && (constM<1.)) pow_hmlgT_M=hmlgT;
        else pow_hmlgT_M=pow(hmlgT,constM);
        if (*eqpsN>0.) pow_eqps_n = +(pow(*eqpsN,constN));
        else pow_eqps_n = 0.;

        yieldStress = // actual yield stress if step is elastic
                (constA+constB*pow_eqps_n)*
                (1-pow_hmlgT_M);

        if( vmStressTrial - yieldStress >= 0 ) {
            MATERIAL_EVALUATION::JohnsonCookSolve(
                vmStressTrial,eqpsN,eqpsNP1,&yieldStress,
                MU,constA,constN,constB,constC,pow_hmlgT_M,dt);

            lambda = *eqpsNP1-*eqpsN ;

            if (lambda>=0.) {
                // RADIAL RETURN
                // update deviatoric stress
                // vmStress = yieldStress in the new condition
                *vmStress = yieldStress;
                
                tempScalar = *vmStress / vmStressTrial;
//                 Wd=0.0;
                for(int n=0;n<numNeigh;n++,neighPtrOverLap++,bondDamageOverLap++,edpN_OverLap++,edpNP1_OverLap++,td_OverLap++){

                    int localId = *neighPtrOverLap;
                    cellVolume = v[localId];
                    const double *XP = &xOverlap[3*localId];
                    const ScalarT *YP = &yOverlap[3*localId];
                    X_dx = XP[0]-X[0];  X_dy = XP[1]-X[1];  X_dz = XP[2]-X[2];
                    Y_dx = YP[0]-Y[0];  Y_dy = YP[1]-Y[1];  Y_dz = YP[2]-Y[2];
                    dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
                    zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
                    omega = scalarInfluenceFunction(zeta,horizon);

                    *td_OverLap *= tempScalar;
                    ti= (1.0-*bondDamage)*(3.0*K*(*theta)/(*m)* omega * zeta);
                    t = ti + *td_OverLap;

                    fx = t * Y_dx / dY; fy = t * Y_dy / dY; fz = t * Y_dz / dY;

                    *(fOwned+0) += fx*cellVolume;
                    *(fOwned+1) += fy*cellVolume;
                    *(fOwned+2) += fz*cellVolume;
                    fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
                    fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
                    fInternalOverlap[3*localId+2] -= fz*selfCellVolume;

                    *edpNP1_OverLap = *edpN_OverLap + *m/5.0 *(*td_OverLap)/omega / (*vmStress) * lambda ;
                    
//                     // compute deviatoric energy density
//                     e = dY - zeta;
//                     if(deltaTemperature) e -= thermalExpansionCoefficient*(*deltaT)*zeta;
//                     ed = e - *theta/3*zeta - *edpNP1_OverLap; // deviatoric Extension
//                     Wd += (1.0-*bondDamageOverLap)* alpha/2 * ed * omega * ed * cellVolume;
                }
//                 cout << "yielding " << sqrt(6*MU*Wd) << " " << *vmStress << endl;
            } else{
                std::cout << "ERROR: Negative delta plastic epsilon after loop\n" <<  "  Delta plastic strain:  " << lambda << std::endl;
                TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "")
            } // end if lambda
        } else { // The step is elastic
            *eqpsNP1= *eqpsN;
            *vmStress = vmStressTrial;

            for(int n=0; n<numNeigh; n++,neighPtrOverLap++,bondDamageOverLap++,edpN_OverLap++,edpNP1_OverLap++,td_OverLap++) {

                int localId = *neighPtrOverLap;
                cellVolume = v[localId];
                const double *XP = &xOverlap[3*localId];
                const ScalarT *YP = &yOverlap[3*localId];
                X_dx = XP[0]-X[0];  X_dy = XP[1]-X[1];  X_dz = XP[2]-X[2];
                Y_dx = YP[0]-Y[0];  Y_dy = YP[1]-Y[1];  Y_dz = YP[2]-Y[2];
                dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
                zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
                omega = scalarInfluenceFunction(zeta,horizon);

                ti= (1.0-*bondDamage)*(3.0*K*(*theta)/(*m)* omega * zeta);
                t = ti + *td_OverLap;

                fx = t * Y_dx / dY;  fy = t * Y_dy / dY;  fz = t * Y_dz / dY;

                *(fOwned+0) += fx*cellVolume;  
                *(fOwned+1) += fy*cellVolume;
                *(fOwned+2) += fz*cellVolume;
                fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
                fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
                fInternalOverlap[3*localId+2] -= fz*selfCellVolume;

            }
//             cout << "elastic " << sqrt(6*MU*Wd) << " " << *vmStress << endl;
            
        }; // end if yield

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
double* VonMisesStress,
const double* deviatoricPlasticExtensionN,
double* deviatoricPlasticExtensionNP1,
const double* EquivalentPlasticStrainN,
double* EquivalentPlasticStrainNP1,
double* deviatoricForceDensity,
const double* deltaTemperature,
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
const double horizon,
const double dt,
const double MeltingTemperature,
const double ReferenceTemperature,
const double constA,
const double constN,
const double constB,
const double constC,
const double constM
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
Sacado::Fad::DFad<double>* VonMisesStress,
const Sacado::Fad::DFad<double>* deviatoricPlasticExtensionN,
Sacado::Fad::DFad<double>* deviatoricPlasticExtensionNP1,
const Sacado::Fad::DFad<double>* EquivalentPlasticStrainN,
Sacado::Fad::DFad<double>* EquivalentPlasticStrainNP1,
Sacado::Fad::DFad<double>* deviatoricForceDensity,
const double* deltaTemperature,
PeridigmNS::Material::BulkMod obj_bulkModulus,
PeridigmNS::Material::ShearMod obj_shearModulus,
PeridigmNS::Material::TempDepConst obj_alphaVol,
const double horizon,
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
