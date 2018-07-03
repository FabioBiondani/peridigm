#include "Peridigm_Material.hpp"
#include "jc_pals.h"
#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include "material_utilities.h"
#include "JohnsonCook.h"

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
)
{
	/*
	 * Compute processor local contribution to internal force
	 */
//    	double BULK_MODULUS = obj_bulkModulus.compute(0.0);
	double SHEAR_MODULUS = obj_shearModulus.compute(0.0);

	//double K = BULK_MODULUS;
	double TWO_MU = 2.0 * SHEAR_MODULUS;

	const double *xOwned = xOverlap;
	const double *yOwnedN = yOverlapN;
	const double *yOwnedNP1 = yOverlapNP1;
	double lambda_X[NUM_LAGRANGE_MULTIPLIERS];
	const double *oc=omega_constant;
	double tau_X[NUM_LAGRANGE_MULTIPLIERS];
	const double *sc=sigma_constant;
	double *fOwned = fInternalOverlap;

	double bond[3];
	double a, b, c;
	double xi;
	double YNP1_dx, YNP1_dy, YNP1_dz;
	double dYNP1;
	double YN_dx, YN_dy, YN_dz;
	double dY_dx, dY_dy, dY_dz;
	double e;
	double eps;
	double omega, sigma, omegaOrdinary;
	double t;
	double fx, fy, fz;
	double cell_volume;
	const double *theta = dilatation;
	const double *p = pals_pressure;
    
    double Wd;
    
    double *vmStress = VonMisesStress;
    const double *edpN = deviatoricPlasticExtensionN;
    double *edpNP1 = deviatoricPlasticExtensionNP1;
    const double *eqpsN = EquivalentPlasticStrainN;
    double *eqpsNP1 = EquivalentPlasticStrainNP1;
    double *td = deviatoricForceDensity;
    
    const double *specu = specularBondPosition;
    double *miPotNP1 = microPotentialNP1;
    double *miPotNP1_Overlap = microPotentialNP1;

//     double hmlgT(0.0);
    double pow_hmlgT_M(0.0);

    double vmStressTrial;
    double pow_eqps_n;    
    double yieldStress, lambda, tempScalar, ti;

	const int *neighPtr = localNeighborList;
    
    const int *neighPtr_ = localNeighborList;
    const double *edpN_ = deviatoricPlasticExtensionN;
    double *edpNP1_ = deviatoricPlasticExtensionNP1;
    double *td_ = deviatoricForceDensity;
    const double *bondDamage_ = bondDamage;

	for(int q=0; q<numOwnedPoints;q++, xOwned+=3, yOwnedN+=3, yOwnedNP1+=3, fOwned+=3, vmStress++, eqpsN++, eqpsNP1++, oc++, sc++, theta++, p++){

		int numNeigh = *neighPtr; neighPtr++; neighPtr_++;
		const double *X = xOwned;
		const double *YN = yOwnedN;
		const double *YNP1 = yOwnedNP1;
		// Collect computed Lagrange multipliers for this point
		for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
			lambda_X[i]=_omega_multipliers[i][q];
			tau_X[i]=_sigma_multipliers[i][q];
		}

		double self_cell_volume = volumeOverlap[q];
        
        Wd=0.0;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,edpN++,edpNP1++,td++){
            *edpNP1=*edpN;

            int localId = *neighPtr;
			cell_volume= volumeOverlap[localId];
			const double *XP = &xOverlap[3*localId];
			const double *YNP1P = &yOverlapNP1[3*localId];
			bond[0] = XP[0]-X[0];  bond[1] = XP[1]-X[1];  bond[2] = XP[2]-X[2];
            a = bond[0];           b = bond[1];           c = bond[2];
			xi = sqrt(a*a+b*b+c*c);
            YNP1_dx = YNP1P[0]-YNP1[0];  YNP1_dy = YNP1P[1]-YNP1[1];  YNP1_dz = YNP1P[2]-YNP1[2];
            dYNP1 = sqrt(YNP1_dx*YNP1_dx+YNP1_dy*YNP1_dy+YNP1_dz*YNP1_dz);
            e = dYNP1 - xi;

			pals_influence<deviatoric_influence> SIGMA(SIGMA_0,*sc,tau_X);
            sigma = SIGMA(bond,horizon);

            eps = e - *theta/3*xi - *edpNP1; // deviatoric Extension
            *td = (1.0-*bondDamage)*(sigma * TWO_MU * eps);

            // compute deviatoric energy density
            Wd += (1.0-*bondDamage)* SHEAR_MODULUS * eps * sigma * eps * cell_volume;
        }
        
        vmStressTrial = sqrt(6*SHEAR_MODULUS*Wd);

        if (*eqpsN>0.) pow_eqps_n = +(pow(*eqpsN,constN));
        else pow_eqps_n = 0.;

        yieldStress = // actual yield stress if step is elastic
                (constA+constB*pow_eqps_n)*
                (1-pow_hmlgT_M);

        if( vmStressTrial - yieldStress >= 0 ) {
            MATERIAL_EVALUATION::JohnsonCookSolve(
                vmStressTrial,eqpsN,eqpsNP1,&yieldStress,
                SHEAR_MODULUS,constA,constN,constB,constC,pow_hmlgT_M,doteps0,dt);

            lambda = *eqpsNP1-*eqpsN ;

            if (lambda>=0.) {
                // RADIAL RETURN
                // update deviatoric stress
                // vmStress = yieldStress in the new condition
                *vmStress = yieldStress;
                
                tempScalar = *vmStress / vmStressTrial;
                
                for(int n=0;n<numNeigh;n++,neighPtr_++,bondDamage_++,edpN_++,edpNP1_++,td_++,specu++,miPotNP1++){

                    int localId = *neighPtr_;
                    cell_volume = volumeOverlap[localId];
                    const double *XP = &xOverlap[3*localId];
                    bond[0] = XP[0]-X[0];
                    bond[1] = XP[1]-X[1];
                    bond[2] = XP[2]-X[2];
                    a = bond[0]; b = bond[1]; c = bond[2];
                    xi = sqrt(a*a+b*b+c*c);
                    const double *YNP1P = &yOverlapNP1[3*localId];
                    YNP1_dx = YNP1P[0]-YNP1[0];  YNP1_dy = YNP1P[1]-YNP1[1];  YNP1_dz = YNP1P[2]-YNP1[2];
                    dYNP1 = sqrt(YNP1_dx*YNP1_dx+YNP1_dy*YNP1_dy+YNP1_dz*YNP1_dz);

                    pals_influence<dilatation_influence> OMEGA(OMEGA_0,*oc,lambda_X);
                    omega = OMEGA(bond,horizon);

                    *td_ *= tempScalar;
                    ti= (1.0-*bondDamage_)*(*p)*omega*xi;
                    t = ti + *td_;

                    fx = t * YNP1_dx / dYNP1; fy = t * YNP1_dy / dYNP1; fz = t * YNP1_dz / dYNP1;

                    *(fOwned+0) += fx*cell_volume;
                    *(fOwned+1) += fy*cell_volume;
                    *(fOwned+2) += fz*cell_volume;
                    fInternalOverlap[3*localId+0] -= fx*self_cell_volume;
                    fInternalOverlap[3*localId+1] -= fy*self_cell_volume;
                    fInternalOverlap[3*localId+2] -= fz*self_cell_volume;

                    *edpNP1_ = *edpN_ + 3.0/2.0 *(*td_)/omega / (*vmStress) * lambda ;


                    if(useSpecularBondPosition){
                        int specuId = int(*specu);
                        const double *YNP = &yOverlapN[3*localId];
                        YN_dx = YNP[0]-YN[0];  YN_dy = YNP[1]-YN[1];  YN_dz = YNP[2]-YN[2];
                        dY_dx = YNP1_dx - YN_dx; dY_dy = YNP1_dy - YN_dy;  dY_dz = YNP1_dz - YN_dz;

                        *miPotNP1+=                 fx*dY_dx + fy*dY_dy + fz*dY_dz ;
                        miPotNP1_Overlap[specuId]+= fx*dY_dx + fy*dY_dy + fz*dY_dz ;
                    }
                }
            }else{
                std::cout << "ERROR: Negative delta plastic epsilon after loop\n" <<  "  Delta plastic strain:  " << lambda << std::endl;
                TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "")
            } // end if lambda
        }else 
        { // The step is elastic
            *eqpsNP1= *eqpsN;
            *vmStress = vmStressTrial;

            for(int n=0; n<numNeigh; n++,neighPtr_++,bondDamage_++,edpN_++,edpNP1_++,td_++,specu++,miPotNP1++) {

                int localId = *neighPtr_;
                cell_volume = volumeOverlap[localId];
                const double *XP = &xOverlap[3*localId];
                bond[0] = XP[0]-X[0];
                bond[1] = XP[1]-X[1];
                bond[2] = XP[2]-X[2];
                a = bond[0]; b = bond[1]; c = bond[2];
                xi = sqrt(a*a+b*b+c*c);
                const double *YNP1P = &yOverlapNP1[3*localId];
                YNP1_dx = YNP1P[0]-YNP1[0];  YNP1_dy = YNP1P[1]-YNP1[1];  YNP1_dz = YNP1P[2]-YNP1[2];
                dYNP1 = sqrt(YNP1_dx*YNP1_dx+YNP1_dy*YNP1_dy+YNP1_dz*YNP1_dz);

                pals_influence<dilatation_influence> OMEGA(OMEGA_0,*oc,lambda_X);
                omega = OMEGA(bond,horizon);

                ti= (1.0-*bondDamage_)*(*p)*omega*xi;
                t = ti + *td_;

                fx = t * YNP1_dx / dYNP1;  fy = t * YNP1_dy / dYNP1;  fz = t * YNP1_dz / dYNP1;

                *(fOwned+0) += fx*cell_volume;
                *(fOwned+1) += fy*cell_volume;
                *(fOwned+2) += fz*cell_volume;
                fInternalOverlap[3*localId+0] -= fx*self_cell_volume;
                fInternalOverlap[3*localId+1] -= fy*self_cell_volume;
                fInternalOverlap[3*localId+2] -= fz*self_cell_volume;

                if(useSpecularBondPosition){
                    int specuId = int(*specu);
                    const double *YNP = &yOverlapN[3*localId];
                    YN_dx = YNP[0]-YN[0];  YN_dy = YNP[1]-YN[1];  YN_dz = YNP[2]-YN[2];
                    dY_dx = YNP1_dx - YN_dx; dY_dy = YNP1_dy - YN_dy;  dY_dz = YNP1_dz - YN_dz;

                    *miPotNP1+=                 fx*dY_dx + fy*dY_dy + fz*dY_dz ;
                    miPotNP1_Overlap[specuId]+= fx*dY_dx + fy*dY_dy + fz*dY_dz ;
                }

            }
            
         }; // end if yield

	}

}

}


}
