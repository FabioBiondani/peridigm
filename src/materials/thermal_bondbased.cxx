//! \file bondbased_thermal_diffusion.cxx

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
#include <algorithm>
#include <vector>
#include <Sacado.hpp>
#include "thermal_bondbased.h"
#include "material_utilities.h"
#include <boost/math/constants/constants.hpp>
#include <iostream>
// ----------------------------------HEAT FLOW---------------------------------
namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeBondBasedHeatFlow
(
	const double*  xOverlap,
	const ScalarT* yOverlap,
	const double* volumeOverlap,
	const double* bondDamage,
	ScalarT* heatFlowOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
    PeridigmNS::Material::TempDepConst obj_thermalConductivity,
	double horizon,
	ScalarT* deltaTemperatureOverlap,
    bool temperatureDependence,
	ScalarT* TCF
)
{
	/*
	 * Compute processor local contribution to internal heat flux
	 */
	double K_T = obj_thermalConductivity.compute(0.0);
    double microConductivity;
	const double PI_G = boost::math::constants::pi<double>();
	const double *xOwned = xOverlap;
// 	const ScalarT *yOwned = yOverlap;
	const ScalarT *deltaTemperatureOwned = deltaTemperatureOverlap;
	const double *v = volumeOverlap;
	ScalarT *heatFlowOwned = heatFlowOverlap;
	const int *neighPtr = localNeighborList;
	double cellVolume;
// 	double omega;
	double X_dx, X_dy, X_dz, zeta;
// 	ScalarT Y_dx, Y_dy, Y_dz, dY;
    ScalarT dT, q1;
    double deltaTdouble=0;
    
    bool boolTCF=false; // Thermal Correction Factor
    if (TCF!=NULL)
        boolTCF=true;
    
// 	loop over all the nodes
	for(int p=0;p<numOwnedPoints;p++, deltaTemperatureOwned++, xOwned +=3,/* yOwned +=3,*/ heatFlowOwned++, TCF++){
		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
// 		const ScalarT *Y = yOwned;
		const ScalarT *deltaT = deltaTemperatureOwned;
        
//      local thermal conductivity
        deltaTdouble = convT2double(*deltaT);
        if(temperatureDependence)
            K_T = obj_thermalConductivity.compute(deltaTdouble);
        
        microConductivity = 6 * K_T /( PI_G * horizon*horizon*horizon*horizon);
        
// 		loop over the horizon region
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const ScalarT *deltaTP = &deltaTemperatureOverlap[localId];
			const double *XP = &xOverlap[3*localId];
// 			const ScalarT *YP = &yOverlap[3*localId];
			X_dx = XP[0]-X[0];
			X_dy = XP[1]-X[1];
			X_dz = XP[2]-X[2];
			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
// 			Y_dx = YP[0]-Y[0];
// 			Y_dy = YP[1]-Y[1];
// 			Y_dz = YP[2]-Y[2];
// 			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
			dT = *deltaTP - *deltaT;
			q1 = (1-*bondDamage)*microConductivity*dT/zeta ;
            if (boolTCF) q1 *= *TCF;
			*heatFlowOwned += q1 *cellVolume;
		}
	}
}

/** Explicit template instantiation for double. */
template void computeBondBasedHeatFlow<double>
(
	const double*  xOverlap,
	const double* yOverlap,
	const double* volumeOverlap,
	const double* bondDamage,
	double* heatFlowOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
    PeridigmNS::Material::TempDepConst obj_thermalConductivity,
    double horizon,
	double* deltaTemperatureOverlap,
    bool temperatureDependence,
	double* TCF
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeBondBasedHeatFlow<Sacado::Fad::DFad<double> >
(
	const double*  xOverlap,
	const Sacado::Fad::DFad<double>* yOverlap,
	const double* volumeOverlap,
	const double* bondDamage,
	Sacado::Fad::DFad<double>* heatFlowOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
    PeridigmNS::Material::TempDepConst obj_thermalConductivity,
	double horizon,
	Sacado::Fad::DFad<double>* deltaTemperatureOverlap,
    bool temperatureDependence,
	Sacado::Fad::DFad<double>* TCF
);

}
