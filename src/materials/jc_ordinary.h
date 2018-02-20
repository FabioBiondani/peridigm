//! \file jc_ordinary.h

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
#ifndef JCORDINARY_H
#define JCORDINARY_H

#include "Peridigm_Material.hpp"

namespace MATERIAL_EVALUATION {

//! Computes contributions to the internal force resulting from owned points.
template<typename ScalarT>
void computeInternalForceJohnsonCookOrdinary
(
const double* xOverlapPtr,
const ScalarT* yOverlapPtr,
const double* mOwned,
const double* volumeOverlapPtr,
const ScalarT* dilatationOwned,
const double* bondDamage,
const double* scfOwned,
ScalarT* fInternalOverlapPtr,
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
);

}

#endif // JCORDINARY_H
