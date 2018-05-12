#ifndef PERIDIGM_VISCOUSMATERIAL_HPP
#define PERIDIGM_VISCOUSMATERIAL_HPP

#include "Peridigm_Material.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Peridigm_DataManager.hpp"

namespace PeridigmNS{
	
	class ViscousMaterial {
	public:
        Teuchos::ParameterList matparams;
		
// 		! Constructor
		ViscousMaterial() {};
		ViscousMaterial(const Teuchos::ParameterList & params) {
            matparams = params;
        };
		
// 		! Destructor
		virtual ~ViscousMaterial(){};
		
        //! Returns a vector of field IDs corresponding to the variables associated with the material.
        virtual std::vector<int> FieldIds() const = 0;

        //! Initialize the material model.
        virtual void
        initialize(const double dt,
                   const int numOwnedPoints,
                   const int* ownedIDs,
                   const int* neighborhoodList,
                   PeridigmNS::DataManager& dataManager) {}
        //! Evaluate the internal force.
        virtual void
        computeForce(const double dt,
                     const int numOwnedPoints,
                     const int* ownedIDs,
                     const int* neighborhoodList,
                     PeridigmNS::DataManager& dataManager) const = 0;
	};

}
#endif // PERIDIGM_VISCOUSMATERIAL_HPP
