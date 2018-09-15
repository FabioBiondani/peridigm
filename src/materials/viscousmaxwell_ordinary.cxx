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
#include "viscousmaxwell_ordinary.h"
#include "material_utilities.h"

using std::cout;
using std::endl;
namespace MATERIAL_EVALUATION {

void computeInternalForceViscousMaxwell
  (
   const double dt,
   const double horizon,
   const double *xOverlap,
   const double *yNOverlap,
   const double *yNP1Overlap,
   const double *mOwned,
   const double* volumeOverlap,
   const double* dilatationOwnedN,
   const double* dilatationOwnedNP1,
   const double* bondDamage,
   const double *edbN,
   double *edbNP1,
   double *fInternalOverlap,
   const int*  localNeighborList,
   int numOwnedPoints,
   PeridigmNS::Material::BulkMod obj_bulkModulus,
   PeridigmNS::Material::ShearMod obj_shearModulus,
   PeridigmNS::Material::TempDepConst obj_lambda,
   PeridigmNS::Material::TempDepConst obj_tau,
   PeridigmNS::Material::TempDepConst obj_alphaVol,
   bool applyThermalStrains,
   bool temperatureDependence,
   const double* deltaTemperatureN,
   const double* deltaTemperatureNP1
)
{

    double MU(0.0), tau(1.0), lambda(1.0), alphaVolN(0.0) , alphaVolNP1(0.0);
    double c1(0.0), decay(0.0), beta_i(0.0);
    MU = obj_shearModulus.compute(0.0);
    tau = obj_tau.compute(0.0);
    lambda = obj_lambda.compute(0.0);
    if ((lambda!=0.0)&&(tau!=0.0)){
        c1 = tau / dt;
        decay = exp(-1.0/c1);
        beta_i=1.-c1*(1.-decay);
    }
    if(applyThermalStrains){
        alphaVolN   = obj_alphaVol.compute(0.0);
        alphaVolNP1 = alphaVolN;
    }
    

	/*
	 * Compute processor local contribution to internal force
	 */
	double omega;

	const double *xOwned = xOverlap;
	const double *yNOwned = yNOverlap;
	const double *yNP1Owned = yNP1Overlap;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const double *thetaN = dilatationOwnedN;
	const double *thetaNP1 = dilatationOwnedNP1;
	double *fOwned = fInternalOverlap;
    const double *deltaT_N   = deltaTemperatureN;
    const double *deltaT_NP1 = deltaTemperatureNP1;

	const int *neighPtr = localNeighborList;
	double cellVolume, dx, dy, dz, zeta, dYN, dYNP1, t, td, edN, edNP1, delta_ed;
    if ((lambda!=0.0)&&(tau!=0.0))
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yNOwned +=3, yNP1Owned +=3, fOwned+=3, m++, thetaN++, thetaNP1++, deltaT_N++, deltaT_NP1++){

        if(temperatureDependence){
            MU   = obj_shearModulus.compute(*deltaT_NP1);
            tau = obj_tau.compute(*deltaT_NP1);
            lambda = obj_lambda.compute(*deltaT_NP1);
            if ((lambda==0.0)||(tau==0.0)) continue;
            
            c1 = tau / dt;
            decay = exp(-1.0/c1);
            beta_i=1.-c1*(1.-decay);
            if(applyThermalStrains){
                alphaVolN   = obj_alphaVol.compute(*deltaT_N);
                alphaVolNP1 = obj_alphaVol.compute(*deltaT_NP1);
            }
        }

        int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const double *YN = yNOwned;
		const double *YNP1 = yNP1Owned;
		double weightedVolume = *m;
		double dilatationN   = *thetaN;
		double dilatationNP1 = *thetaNP1;
		double alpha = 15.0*MU/weightedVolume;
		double selfCellVolume = v[p];
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,edbN++,edbNP1++){
			
			double damageNP1 = (1.0-*bondDamage);
            if (damageNP1==0.0) continue;
            
            int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP    = &xOverlap[3*localId];
			const double *YPN   = &yNOverlap[3*localId];
			const double *YPNP1 = &yNP1Overlap[3*localId];
			dx = XP[0]-X[0];
			dy = XP[1]-X[1];
			dz = XP[2]-X[2];
			zeta = sqrt(dx*dx+dy*dy+dz*dz);

			/*
			 * volumetric scalar state
			 */
			double eiN   = dilatationN * zeta / 3.0;
			double eiNP1 = dilatationNP1 * zeta / 3.0;

			/*
			 * COMPUTE edN
			 */
			dx = YPN[0]-YN[0];
			dy = YPN[1]-YN[1];
			dz = YPN[2]-YN[2];
			dYN = sqrt(dx*dx+dy*dy+dz*dz);
			/*
			 */
			edN = (dYN - zeta) - eiN;
            if(applyThermalStrains) edN -= alphaVolN*(*deltaT_N)*zeta;

            
			/*
			 * COMPUTE edNP1
			 */
			dx = YPNP1[0]-YNP1[0];
			dy = YPNP1[1]-YNP1[1];
			dz = YPNP1[2]-YNP1[2];
			dYNP1 = sqrt(dx*dx+dy*dy+dz*dz);
			edNP1 = (dYNP1 - zeta) - eiNP1;
            if(applyThermalStrains) edNP1 -= alphaVolNP1*(*deltaT_NP1)*zeta;

			/*
			 * Increment to deviatoric extension state
			 */
			delta_ed = edNP1-edN;
			/*
			 * Integrate back extension state forward in time
			 */
			*edbNP1 = edN * (1-decay) + (*edbN)*decay  + beta_i * delta_ed;

			/*
			 * Compute deviatoric force state
			 */
            omega = scalarInfluenceFunction(zeta,horizon);
			td = lambda * alpha * omega * ( edNP1 - *edbNP1 );

            /*
			 * Note that damage has already been applied once to 'td' (through ed) above.
			 */
			t = damageNP1 * td;
			double fx = t * dx / dYNP1;
			double fy = t * dy / dYNP1;
			double fz = t * dz / dYNP1;

			*(fOwned+0) += fx*cellVolume;
			*(fOwned+1) += fy*cellVolume;
			*(fOwned+2) += fz*cellVolume;
			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
		}
	}
}

}

