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

#include <cmath>
#include <Sacado.hpp>
#include "jc_ordinary.h"
#include "material_utilities.h"

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
		double BULK_MODULUS,
		double SHEAR_MODULUS,
        double horizon,
        ScalarT* ElasticEnergyDensity,
        ScalarT* DeviatoricElasticEnergyDensity,
        ScalarT* VonMisesStress,
        ScalarT* EquivalentStrain,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;

	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
    const double *deltaT = deltaTemperature;
	const double *m = mOwned;
	const double *v = volumeOverlap;
    const double *scf = scfOwned;
	const ScalarT *theta = dilatationOwned;
	ScalarT *fOwned = fInternalOverlap;
    ScalarT *W  = ElasticEnergyDensity;
    ScalarT *Wd = DeviatoricElasticEnergyDensity;
    ScalarT *sigmaVM = VonMisesStress;
    ScalarT *eps_eq = EquivalentStrain;

	const int *neighPtr = localNeighborList;
	double cellVolume, alpha, X_dx, X_dy, X_dz, zeta, omega;
	ScalarT Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz, e, c1;
    ScalarT ed;
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, W++, Wd++, sigmaVM++, eps_eq++, deltaT++, m++, theta++, scf++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		alpha = 15.0*MU/(*m);
		alpha *= (*scf);
		double selfCellVolume = v[p];
        *Wd=0;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
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
			// c1 = omega*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
			c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0);
			t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);
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
            ed = e - *theta/3;
            *Wd += (1.0-*bondDamage)* alpha/2 * ed * omega * ed * cellVolume;

		}
		*W = K/2 * (*theta) * (*theta) + *Wd;
        *sigmaVM = sqrt(6*MU*(*Wd));
        *eps_eq = *sigmaVM/(3*MU);

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
		double BULK_MODULUS,
		double SHEAR_MODULUS,
        double horizon,
        double* ElasticEnergyDensity,
        double* DeviatoricElasticEnergyDensity,
        double* VonMisesStress,
        double* EquivalentStrain,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
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
		double BULK_MODULUS,
		double SHEAR_MODULUS,
        double horizon,
        Sacado::Fad::DFad<double>* ElasticEnergyDensity,
        Sacado::Fad::DFad<double>* DeviatoricElasticEnergyDensity,
        Sacado::Fad::DFad<double>* VonMisesStress,
        Sacado::Fad::DFad<double>* EquivalentStrain,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
);

}
