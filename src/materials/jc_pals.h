#include <math.h>
#include <vector>
#include "damagepals.h"

#include "Peridigm_InfluenceFunction.hpp"



namespace MATERIAL_EVALUATION {

namespace PALS {

void computeInternalForceJohnsonCookPals
(
	const double *xOverlap,
	const double *yOverlapN,
	const double *yOverlapNP1,
	const double *volumeOverlap,
	const std::vector<const double *>& _omega_multipliers,
	const double *omega_constant,
	const std::vector<const double *>& _sigma_multipliers,
	const double *sigma_constant,
	const double *dilatation,
	const double *pals_pressure,
	double *fInternalOverlap,
	const int *localNeighborList,
	int numOwnedPoints,
    PeridigmNS::Material::BulkMod obj_bulkModulus,
    PeridigmNS::Material::ShearMod obj_shearModulus,
	double horizon,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
    const double* bondDamage,
    const bool useSpecularBondPosition,
    const double* specularBondPosition,
    double* microPotentialN,
    double* microPotentialNP1,
 
    double* VonMisesStress,
    const double* deviatoricPlasticExtensionN,
    double* deviatoricPlasticExtensionNP1,
    const double* EquivalentPlasticStrainN,
    double* EquivalentPlasticStrainNP1,
    double* deviatoricForceDensity,

    const double dt,
    const double MeltingTemperature,
    const double ReferenceTemperature,
    const double constA,
    const double constN,
    const double constB,
    const double constC,
    const double constM,
    const double doteps0
);

}

}


