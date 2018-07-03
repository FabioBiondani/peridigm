/*
 * Peridigm_PALS_Model.hpp
 *
 *  Created on: Aug 21, 2013
 *      Author: jamitch
 */

#ifndef PERIDIGM_JOHNSONCOOKPALSMATERIAL_HPP_
#define PERIDIGM_JOHNSONCOOKPALSMATERIAL_HPP_

#include "Peridigm_Material.hpp"
#include "Peridigm_InfluenceFunction.hpp"

namespace PeridigmNS {

class JohnsonCookPalsMaterial : public Material {
public:

	typedef PeridigmNS::InfluenceFunction::functionPointer FunctionPointer;

	//! Constructor.
	JohnsonCookPalsMaterial(const Teuchos::ParameterList & params);

	//! Destructor.
	virtual ~JohnsonCookPalsMaterial();

	//! Return name of material type
	virtual std::string Name() const { return("Johnson-Cook Pals"); }

	//! Returns the density of the material.
	virtual double Density() const { return m_density; }

	//! Returns the bulk modulus of the material.
	virtual double BulkModulus() const { return m_bulkModulus; }

	//! Returns the shear modulus of the material.
	virtual double ShearModulus() const { return m_shearModulus; }

	//! Returns the horizon.
	virtual double Horizon() const { return m_horizon; }

	//! Returns a vector of field IDs corresponding to the variables associated with the material.
	virtual std::vector<int> FieldIds() const { return m_fieldIds; }

	//! Initialized data containers and computes weighted volume.
	virtual void
	initialize(const double dt,
			  const int numOwnedPoints,
			  const int* ownedIDs,
			  const int* neighborhoodList,
			  PeridigmNS::DataManager& dataManager);

	//! Evaluate the internal force.
	virtual void
	computeForce(const double dt,
		 const int numOwnedPoints,
		 const int* ownedIDs,
		 const int* neighborhoodList,
				PeridigmNS::DataManager& dataManager) const;

    //! Compute stored elastic density energy.
    virtual void
    computeStoredElasticEnergyDensity(const double dt,
                                      const int numOwnedPoints,
                                      const int* ownedIDs,
                                      const int* neighborhoodList,
                                      PeridigmNS::DataManager& dataManager) const;

 protected:


   // material parameters
   double m_bulkModulus;
   double m_shearModulus;
   BulkMod obj_bulkModulus;
   ShearMod obj_shearModulus;
   double m_density;
   double m_horizon;

   // Influence functions
   FunctionPointer m_OMEGA_0;
   FunctionPointer m_SIGMA_0;

    double m_MeltingTemperature;
    double m_ReferenceTemperature;
    double m_A;
    double m_N;
    double m_B;
    double m_C;
    double m_M;
    double m_doteqps0;

   // field spec ids for all relevant data
   std::vector<int> m_fieldIds;
   int m_volumeFieldId;
   int m_weightedVolumeFieldId;
   int m_normalizedWeightedVolumeFieldId;
   int m_dilatationFieldId;
   int m_palsPressureFieldId;
   int m_modelCoordinatesFieldId;
   int m_coordinatesFieldId;
   int m_forceDensityFieldId;
   int m_damageFieldId;
   int m_bondDamageFieldId;

   const int num_lagrange_multipliers;
   int m_dilatationNormalizationFieldId;
   int m_deviatoricNormalizationFieldId;
   std::vector<int> m_dilatationLagrangeMultiplersFieldIds;
   std::vector<int> m_deviatoricLagrangeMultiplersFieldIds;

    int m_VonMisesStressFieldId;
    
    int m_deviatoricPlasticExtensionFieldId;
    int m_equivalentPlasticStrainFieldId;
    int m_deviatoricForceDensityFieldId;

   int m_specularBondPositionFieldId;
   int m_microPotentialFieldId;

   bool m_useSpecularBondPosition;


 };


}


#endif /* PERIDIGM_PALS_MODEL_HPP_ */
