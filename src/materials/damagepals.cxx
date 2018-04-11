#include "damagepals.h"
#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

namespace PALS {

//Lapack linear equations
//http://www.netlib.org/lapack/lug/node38.html

// Factorize general matrix
// http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f_source.html

// Solve general matrix
// http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f_source.html

//Cholesky factorize (dpotrf)
//http://www.netlib.org/lapack/explore-html/df/da8/VARIANTS_2cholesky_2RL_2dpotrf_8f_source.html

//Solve using factorize (dpotrs)
//http://www.math.utah.edu/software/lapack/lapack-d/dpotrs.html

extern "C" {
	/*
	 * Lapack
	 * Cholesky factorize dense matrix
	 * Cholesky backsolve dense matrix
	 */
	void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
	void dpotrs_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);

	/*
	 * Lapack
	 * Factorize general 'dense' matrix
	 */
	void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
	void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
}

using std::setw;
using std::vector;

void
compute_lagrange_multipliers
(
	const double *xOverlap,
	const double *volumeOverlap,
	int num_owned_points,
	const int *localNeighborList,
	double horizon,
	vector<double *>& omega_multipliers,
	double *omega_constants,
	vector<double *>& sigma_multipliers,
	double *sigma_constants,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
    const double *bondDamage
)
{

	const double *X=xOverlap;
	const int *neigh_X=localNeighborList;
	int num_neigh_X=*neigh_X;
	double *oc_X=omega_constants;
	double *sc_X=sigma_constants;
	double omega_X[NUM_LAGRANGE_MULTIPLIERS], sigma_X[NUM_LAGRANGE_MULTIPLIERS];
	for(int p=0;p<num_owned_points;p++,X+=3,oc_X++,sc_X++){

		// evaluate lagrange multipliers and normalizing constants for point X
		//print_point_3d(std::cout,p,X);
		compute_lagrange_multipliers_point(X,xOverlap,volumeOverlap,neigh_X,horizon,omega_X,oc_X,sigma_X,sc_X,OMEGA_0,SIGMA_0,bondDamage);

		// Collect computed Lagrange multipliers for this point
		for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
			omega_multipliers[i][p]=omega_X[i];
			sigma_multipliers[i][p]=sigma_X[i];
		}
		//print_N_vector(std::cout,"lagrange_multipliers dilatation",6,sig_X);
		//print_N_vector(std::cout,"lagrange_multipliers deviatoric",6,tau_X);

		// move to next point
		num_neigh_X=*neigh_X;
		neigh_X+=(num_neigh_X+1);
        bondDamage+=num_neigh_X;
	}

}

void
compute_lagrange_multipliers
(
	const double *xOverlap,
	const double *volumeOverlap,
	int num_owned_points,
	const int *localNeighborList,
	double horizon,
	vector<double *>& omega_multipliers,
	double *omega_constants,
	vector<double *>& sigma_multipliers,
	double *sigma_constants,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
    const double* damageN,
    const double* damageNP1,
    const double *bondDamage
)
{

	const double *X=xOverlap;
	const int *neigh_X=localNeighborList;
	int num_neigh_X=*neigh_X;
	double *oc_X=omega_constants;
	double *sc_X=sigma_constants;
	double omega_X[NUM_LAGRANGE_MULTIPLIERS], sigma_X[NUM_LAGRANGE_MULTIPLIERS];
	for(int p=0;p<num_owned_points;p++,X+=3,oc_X++,sc_X++, damageN++, damageNP1++){

        if (*damageNP1 > *damageN){
            if (*damageNP1<=DAMAGE_TOLERANCE){
                // evaluate lagrange multipliers and normalizing constants for point X
                //print_point_3d(std::cout,p,X);
                compute_lagrange_multipliers_point(X,xOverlap,volumeOverlap,neigh_X,horizon,omega_X,oc_X,sigma_X,sc_X,OMEGA_0,SIGMA_0,bondDamage);

                // Collect computed Lagrange multipliers for this point
                for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
                    omega_multipliers[i][p]=omega_X[i];
                    sigma_multipliers[i][p]=sigma_X[i];
                }
                //print_N_vector(std::cout,"lagrange_multipliers dilatation",6,sig_X);
                //print_N_vector(std::cout,"lagrange_multipliers deviatoric",6,tau_X);
            }
            else{
                for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
                    omega_multipliers[i][p]=0.0;
                    sigma_multipliers[i][p]=0.0;
                }
            }
        }

		// move to next point
		num_neigh_X=*neigh_X;
		neigh_X+=(num_neigh_X+1);
        bondDamage+=num_neigh_X;
	}

}


void
compute_lagrange_multipliers_point
(
	const double *X,
	const double *xOverlap,
	const double *volumeOverlap,
	const int *neighborhood,
	double h,
	double *omega_multipliers,
	double *omega_constant,
	double *sigma_multipliers,
	double *sigma_constant,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0,
	const double* bondDamage
)
{
	/*
	 * Compute Lagrange multipliers for point X: 6 each (omega_multipliers, sigma_multipliers)
	 * INPUT
	 * X: 3D point
	 * q: neighbors of X
	 * vol_q: volume of each neighbor q
	 * num_neigh: number of neighbors @ X
	 * h: horizon
	 * OUTPUT
	 * omega_multipliers[6]: Lagrange multipliers for dilatation
	 * sigma_multipliers[6]: Lagrange multipliers for deviatoric
	 * omega_constant[1]: normalizing constant for trial dilatation influence function
	 * sigma_constant[1]: normalizing constant for trial deviatoric influence function
	 */

	/*
	 * RHS side vectors of linear problems
	 */
	double rhs_dil[NUM_LAGRANGE_MULTIPLIERS];
	double rhs_dev[NUM_LAGRANGE_MULTIPLIERS];

	// Matrix at point
	double k[NUM_LAGRANGE_MULTIPLIERS*NUM_LAGRANGE_MULTIPLIERS];
	double k_dev[NUM_LAGRANGE_MULTIPLIERS*NUM_LAGRANGE_MULTIPLIERS];
	// Initialize RHS vectors and Matrices
	for(int i=0,p=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
		rhs_dil[i]=0.0; rhs_dev[i]=0.0;
		for(int j=0;j<NUM_LAGRANGE_MULTIPLIERS;j++,p++){
			k[p]=0.0;
			k_dev[p]=0.0;
		}
	}

	/*
	 * Normalize Gaussian influence function at point
	 */
    const int *neigh=neighborhood;
	int num_neigh=*neigh; neigh++;
// 	double m(0.0),sum_dev(0.0);
	const double x=*X, y=*(X+1), z=*(X+2);
	for(int iq=0;iq<num_neigh;iq++){

		int neigh_local_id=*neigh; neigh++;
		const double *q=xOverlap+3*neigh_local_id;
		double vol=volumeOverlap[neigh_local_id];
        double damage=bondDamage[neigh_local_id];
		double bond[3]={*(q+0)-x,*(q+1)-y,*(q+2)-z};
		double a=bond[0], b=bond[1], c=bond[2];
		double a2=a*a, b2=b*b, c2=c*c;
		double ab=a*b, ac=a*c, bc=b*c;
		/*
		 * Components of extension state
		 * NOTE
		 * Each of these extension state has a bond in the denominator;
		 * But the dilatation matrix is formed by the product of the bond
		 * with the extension state and hence the bond cancels.
		 */
		double e_k[]={a2,b2,c2,2*ab,2*ac,2*bc};

		/*
		 * Symmetric dilatation matrix;
		 */
		for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
			double e_i=e_k[i];
			int col=i*NUM_LAGRANGE_MULTIPLIERS;
			for(int j=0;j<NUM_LAGRANGE_MULTIPLIERS;j++){
				double e_j=e_k[j];
				double kij=e_i*e_j*vol;
				k[col+j]+=(1.0-damage)*kij;
			}
		}

		/*
		 * Components of deviatoric extension state
		 */
		double xi2=a2+b2+c2;
		double r=std::sqrt(xi2);
		double avg=r/3.0;
		double epsilon_k[]={a2/r-avg,b2/r-avg,c2/r-avg,2*ab/r,2*ac/r,2*bc/r};

		/*
		 * Symmetric deviatoric matrix;
		 */
		for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
			double epsilon_i=epsilon_k[i];
			double epsilon_i2=epsilon_i*epsilon_i;
			int col=i*NUM_LAGRANGE_MULTIPLIERS;
			for(int j=0;j<NUM_LAGRANGE_MULTIPLIERS;j++){
				double epsilon_j=epsilon_k[j];
				double epsilon_j2=epsilon_j*epsilon_j;
				double kij=epsilon_i2*epsilon_j2*vol;
				k_dev[col+j]+=(1.0-damage)*kij;
			}
		}

		/*
		 * Reference influence function evaluation
		 */
		double omega_0=OMEGA_0(r,h);
		double sigma_0=SIGMA_0(r,h);

		// sums for normalization of OMEGA_0 and SIGMA_0
		// m+=omega_0*xi2*vol;
		// double epsilon_0=2*(ab+bc+ac)/r;
		// sum_dev+=sigma_0*epsilon_0*epsilon_0*vol;

		/*
		 * RHS vectors
		 */
		for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
			rhs_dil[i]+=(1.0-damage)*omega_0*e_k[i]*vol;
			rhs_dev[i]+=(1.0-damage)*sigma_0*epsilon_k[i]*epsilon_k[i]*vol;
		}

	}

   // normalizing scalars
   // const double D=3.0;
   // const double norm_omega_0=D/m;
   // const double norm_sigma_0=6.0/sum_dev;

	// finalize RHS vectors
	const double c_1=1.0, c_2=1.0;
	rhs_dil[0]=1.0-c_1*rhs_dil[0];
	rhs_dil[1]=1.0-c_1*rhs_dil[1];
	rhs_dil[2]=1.0-c_1*rhs_dil[2];
	rhs_dil[3]=0.0-c_1*rhs_dil[3];
	rhs_dil[4]=0.0-c_1*rhs_dil[4];
	rhs_dil[5]=0.0-c_1*rhs_dil[5];

	const double c23=2.0/3.0;
	rhs_dev[0]=c23-c_2*rhs_dev[0];
	rhs_dev[1]=c23-c_2*rhs_dev[1];
	rhs_dev[2]=c23-c_2*rhs_dev[2];
	rhs_dev[3]=2.0-c_2*rhs_dev[3];
	rhs_dev[4]=2.0-c_2*rhs_dev[4];
	rhs_dev[5]=2.0-c_2*rhs_dev[5];

	/*
	 * ISSUE is that dx, dy, or dz may be zero;
	 * For 2D meshes in x-y plane, this occurs because all of the
	 * z-coordinate values have the same constant value and hence dz=0.
	 * The linear problem in these cases is ill-defined; Following code fixes.
	 */
	int i_diag[]={0,1*6+2-1,2*6+3-1,3*6+4-1,4*6+5-1,5*6+6-1};
	double trace_k=0;
	for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++)
		trace_k+=k[i_diag[i]];
	double small=1.0e-15;
	for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
		int d=i_diag[i];
		if((trace_k==0.0) || (k[d]/trace_k < small)) {
			/*
			 * Zero out row/col for column i
			 * Remember that k is column major
			 * 'r' tracks across  row 'i'
			 * 'c' tracks down column 'i'
			 */
			int N=NUM_LAGRANGE_MULTIPLIERS;
			for(int j=0,r=i,c=i*N;j<N;j++,r+=N,c++){
				k[r]=0.0;
				k[c]=0.0;
				k_dev[r]=0.0;
				k_dev[c]=0.0;
			}
			// Set diagonal to 1.0
			k[d]=1.0;
			k_dev[d]=1.0;
			// Set rhs vectors to zero
			rhs_dil[i]=0.0;
			rhs_dev[i]=0.0;
		}
	}

    int solveError;
	//print_symmetrix_6x6(std::cout,"K_DIL",k);
	//print_N_vector(std::cout,"RHS_DIL",NUM_LAGRANGE_MULTIPLIERS,rhs_dil);
	solveError = solve_linear_problem(k,rhs_dil);
	//print_N_vector(std::cout,"RHS_DIL",NUM_LAGRANGE_MULTIPLIERS,rhs_dil);

//	print_symmetrix_6x6(std::cout,"K_DEV",k_dev);
//	print_N_vector(std::cout,"RHS_DEV",NUM_LAGRANGE_MULTIPLIERS,rhs_dev);
    solveError+= solve_linear_problem(k_dev,rhs_dev);
//	print_N_vector(std::cout,"RHS_DEV",NUM_LAGRANGE_MULTIPLIERS,rhs_dev);

	// Copy solution into output vectors
    if (solveError ==0)
        for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
            omega_multipliers[i]=rhs_dil[i];
            sigma_multipliers[i]=rhs_dev[i];
        }
    else
        for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
            omega_multipliers[i]=0.0;
            sigma_multipliers[i]=0.0;
        }

	/*
	 * This code computes scaling factors for the composite influence
	 * functions.  Generally, it turns out the scaling for 'OMEGA' is 1.0
	 * but scaling for SIGMA is a non-trivial number.
	 * ALSO -- note that the 'pals_influence' function definition can
	 * be adjusted to scale on OMEGA_0 and SIGMA_0 only.  Its a bit
	 * subtle here -- the functions below are computed on the assumption
	 * that both parts of the influence function are scaled, ie, these
	 * are computed using the 'composite/whole' function.  Assuming that
	 * there are no bugs, it turns out that scaling 'SIGMA' does
	 * not always work;  Turning scaling off has no effect on 'OMEGA' since 
    * the scaling is implicitly correct with the choice of 'matching 
    * deformations' used here.  However, scaling on the deviatoric 
    * piece does not generally work but it does work when scaling is 
    * not used.
	 */
   //const double one(1.0);
   //pals_influence<dilatation_influence> OMEGA(OMEGA_0,one,rhs_dil);
   //pals_influence<deviatoric_influence> SIGMA(SIGMA_0,one,rhs_dev);
   //const double norm_omega=compute_normalizing_constant_point(OMEGA,X,xOverlap,volumeOverlap,neighborhood,h);
   //const double norm_sigma=compute_normalizing_constant_point(SIGMA,X,xOverlap,volumeOverlap,neighborhood,h);
	//*omega_constant=norm_omega;
	//*sigma_constant=norm_sigma;

   /*
    * Turn scaling off;
    */
	*omega_constant=1.0;
	*sigma_constant=1.0;

}

void computeWeightedVolume
(
	const double *xOverlap,
	const double *volumeOverlap,
	const vector<const double *>& sigma_multipliers,
	const double *sigma_constant,
	double *weighted_volume,
	int numOwnedPoints,
	const int* localNeighborList,
	double horizon,
    const double* damageN,
    const double* damageNP1,
    const double* bondDamage,
	const FunctionPointer SIGMA_0
)
{
	double bond[3];
	const double *xOwned = xOverlap;
	double tau_X[NUM_LAGRANGE_MULTIPLIERS];
	const double *sc=sigma_constant;
	double *m=weighted_volume;
	const int *neighPtr = localNeighborList;
	double cellVolume;
	for(int q=0; q<numOwnedPoints;q++, xOwned+=3, sc++, m++, damageN++, damageNP1++){
		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
        if((*damageNP1 == *damageN)&&(*damageN>=DAMAGE_TOLERANCE)){
            neighPtr+=numNeigh;
            bondDamage+=numNeigh;
            continue;
        }
        else{
            *m=0.0;

            // Collect computed Lagrange multipliers for this point
            for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
                tau_X[i]=sigma_multipliers[i][q];
            }

            for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
                int localId = *neighPtr;
                cellVolume = volumeOverlap[localId];
                const double *XP = &xOverlap[3*localId];
                bond[0]=XP[0]-X[0];
                bond[1]=XP[1]-X[1];
                bond[2]=XP[2]-X[2];
                double a = bond[0];
                double b = bond[1];
                double c = bond[2];
                double xi2 = a*a+b*b+c*c;
                pals_influence<deviatoric_influence> SIGMA(SIGMA_0,*sc,tau_X);
            double sigma = SIGMA(bond,horizon);
                *m += (1.0-*bondDamage)*sigma*xi2*cellVolume;
            }
        }
	}
}

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
){

	double K = BULK_MODULUS;
	double TWO_MU = 2.0 * SHEAR_MODULUS;
	double bond[3];
	const double *xOwned = xOverlap;
	const double *yOwned = yOverlap;
	double lambda_X[NUM_LAGRANGE_MULTIPLIERS];
	const double *oc=omega_constant;
	double tau_X[NUM_LAGRANGE_MULTIPLIERS];
	const double *sc=sigma_constant;
	const double *m=weighted_volume;
	double *theta = dilatation;
	double *p = pals_pressure;
	double cellVolume;
	const int *neighPtr = localNeighborList;
	for(int q=0; q<numOwnedPoints;q++, xOwned+=3, yOwned+=3, oc++, sc++, m++, theta++, p++){
		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const double *Y = yOwned;
		// Collect computed Lagrange multipliers for this point
		for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
			lambda_X[i]=_omega_multipliers[i][q];
			tau_X[i]=_sigma_multipliers[i][q];
		}
		*theta = double(0.0);
		*p = double(0.0);
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = volumeOverlap[localId];
			const double *XP = &xOverlap[3*localId];
			const double *YP = &yOverlap[3*localId];
			bond[0]=XP[0]-X[0];
			bond[1]=XP[1]-X[1];
			bond[2]=XP[2]-X[2];
			double a = bond[0];
			double b = bond[1];
			double c = bond[2];
			double xi2 = a*a+b*b+c*c;
			a = YP[0]-Y[0];
			b = YP[1]-Y[1];
			c = YP[2]-Y[2];
			double dY = a*a+b*b+c*c;
			double x = sqrt(xi2);
			double e = sqrt(dY)-x;
			pals_influence<dilatation_influence> OMEGA(OMEGA_0,*oc,lambda_X);
			pals_influence<deviatoric_influence> SIGMA(SIGMA_0,*sc,tau_X);
         double omega = OMEGA(bond,horizon);
         double sigma = SIGMA(bond,horizon);
         double omega_x=omega*x;
         double sigma_x=sigma*x;
         *theta+=(1.0-*bondDamage)*omega_x*e*cellVolume;
         *p+=-(1.0-*bondDamage)*(TWO_MU*sigma_x/3.0)*e*cellVolume;
		}
		// Final piece of pals_pressure requires dilatation that is only ready here
		*p+=(K+TWO_MU*(*m)/9.0)*(*theta);
	}
}

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
    double* microPotentialNP1,
    double* weightedVolume
)
{
	/*
	 * Compute processor local contribution to internal force
	 */
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
    double *m=weightedVolume;

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
    
    const double *specu = specularBondPosition;
    double *miPotN   = microPotentialN;
    double *miPotNP1 = microPotentialNP1;
    double *miPotNP1_Overlap = microPotentialNP1;

    
	const int *neighPtr = localNeighborList;
	for(int q=0; q<numOwnedPoints;q++, xOwned+=3, yOwnedN+=3, yOwnedNP1+=3, fOwned+=3, oc++, sc++, theta++, p++,m++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const double *YN = yOwnedN;
		const double *YNP1 = yOwnedNP1;
		// Collect computed Lagrange multipliers for this point
		for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
			lambda_X[i]=_omega_multipliers[i][q];
			tau_X[i]=_sigma_multipliers[i][q];
		}

		double self_cell_volume = volumeOverlap[q];
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,specu++,miPotN++,miPotNP1++){
			int localId = *neighPtr;
			const double *XP = &xOverlap[3*localId];
			const double *YNP1P = &yOverlapNP1[3*localId];
			bond[0] = XP[0]-X[0];
			bond[1] = XP[1]-X[1];
			bond[2] = XP[2]-X[2];
			a = bond[0];
			b = bond[1];
			c = bond[2];
			xi = sqrt(a*a+b*b+c*c);
			YNP1_dx = YNP1P[0]-YNP1[0];
			YNP1_dy = YNP1P[1]-YNP1[1];
			YNP1_dz = YNP1P[2]-YNP1[2];
			dYNP1 = sqrt(YNP1_dx*YNP1_dx+YNP1_dy*YNP1_dy+YNP1_dz*YNP1_dz);
            e = dYNP1 - xi;
            eps=e-(*theta)*xi/3.0;
			pals_influence<dilatation_influence> OMEGA(OMEGA_0,*oc,lambda_X);
			pals_influence<deviatoric_influence> SIGMA(SIGMA_0,*sc,tau_X);
            omega = OMEGA(bond,horizon);
            sigma = SIGMA(bond,horizon);


            t=(1.0-*bondDamage)*((*p)*omega*xi+TWO_MU*sigma*eps);
			fx = t * YNP1_dx / dYNP1;
			fy = t * YNP1_dy / dYNP1;
			fz = t * YNP1_dz / dYNP1;
			cell_volume= volumeOverlap[localId];
			*(fOwned+0) += fx*cell_volume;
			*(fOwned+1) += fy*cell_volume;
			*(fOwned+2) += fz*cell_volume;
			fInternalOverlap[3*localId+0] -= fx*self_cell_volume;
			fInternalOverlap[3*localId+1] -= fy*self_cell_volume;
			fInternalOverlap[3*localId+2] -= fz*self_cell_volume;

            if(useSpecularBondPosition){
//                 omegaOrdinary = scalarInfluenceFunction(xi,horizon);
//                 double alpha = 15*SHEAR_MODULUS/(*m);
//                 t=(1.0-*bondDamage)*omegaOrdinary*((*p)*xi+alpha*eps);
                fx = t * YNP1_dx / dYNP1;
                fy = t * YNP1_dy / dYNP1;
                fz = t * YNP1_dz / dYNP1;
                int specuId = int(*specu);
                const double *YNP = &yOverlapN[3*localId];
                YN_dx = YNP[0]-YN[0];  YN_dy = YNP[1]-YN[1];  YN_dz = YNP[2]-YN[2];
                dY_dx = YNP1_dx - YN_dx; dY_dy = YNP1_dy - YN_dy;  dY_dz = YNP1_dz - YN_dz;
                
                *miPotNP1                 += fx*dY_dx + fy*dY_dy + fz*dY_dz ;
                miPotNP1_Overlap[specuId] += fx*dY_dx + fy*dY_dy + fz*dY_dz ;
            }

		}

	}

}

}


}
