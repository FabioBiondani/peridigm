#include <math.h>
#include <vector>
#include "pals.h"

#include "Peridigm_InfluenceFunction.hpp"



namespace MATERIAL_EVALUATION {

namespace PALS {


void
compute_lagrange_multipliers
(
	const double *xOverlap,
	const double *volumeOverlap,
	int num_owned_points,
	const int *localNeighborList,
	double horizon,
	std::vector<double *>& sig_owned,
	double *beta_sig_owned,
	std::vector<double *>& tau_owned,
	double *beta_tau_owned,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
    const double* bondDamage
);

void
compute_lagrange_multipliers
(
	const double *xOverlap,
	const double *volumeOverlap,
	int num_owned_points,
	const int *localNeighborList,
	double horizon,
	std::vector<double *>& sig_owned,
	double *beta_sig_owned,
	std::vector<double *>& tau_owned,
	double *beta_tau_owned,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
    const double* damageN,
    const double* damageNP1,
    const double* bondDamage
);

void
compute_lagrange_multipliers_point
(
	const double *X,
	const double *xOverlap,
	const double *volumeOverlap,
	const int* neigh,
	double h,
	double *omega_multipliers,
	double *omega_constant,
	double *sigma_multipliers,
	double *sigma_constant,
	const FunctionPointer dilation_influence_function,
	const FunctionPointer deviatoric_influence_function,
    const double* bondDamage
);

void computeWeightedVolume
(
	const double *xOverlap,
	const double *volumeOverlap,
	const std::vector<const double *>& _sigma_multipliers,
	const double *sigma_constant,
	double *weighted_volume,
	int myNumPoints,
	const int* localNeighborList,
	double horizon,
    const double* damageN,
    const double* damageNP1,
    const double* bondDamage,
	const FunctionPointer SIGMA_0=PeridigmNS::InfluenceFunction::self().getInfluenceFunction()
);

void computeDilatationAndPalsPressure
(
	const double *xOverlap,
	const double *yOverlap,
	const double *volumeOverlap,
	const std::vector<const double *>& _omega_multipliers,
	const double *omega_constant,
	const std::vector<const double *>& _sigma_multipliers,
	const double *sigma_constant,
	const double *weighted_volume,
	double *dilatation,
	double *pals_pressure,
	const int *localNeighborList,
	int numOwnedPoints,
	double BULK_MODULUS,
	double SHEAR_MODULUS,
	double horizon,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
    const double* bondDamage
);

void computeInternalForceDamagePals
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
	double BULK_MODULUS,
	double SHEAR_MODULUS,
	double horizon,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
    const double* bondDamage,
    const bool useSpecularBondPosition,
    const double* specularBondPosition,
    double* microPotentialN,
    double* microPotentialNP1
);

}

}


