#ifndef PERIDIGM_BASETHERMALMATERIAL_HPP
#define PERIDIGM_BASETHERMALMATERIAL_HPP

#include "Peridigm_Material.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Peridigm_DataManager.hpp"

namespace PeridigmNS{
	
	class ThermalMaterial {
	public:
        Teuchos::ParameterList matparams;
		
// 		! Constructor
		ThermalMaterial() {};
		ThermalMaterial(const Teuchos::ParameterList & params) {
            matparams = params;
        };
		
// 		! Destructor
		virtual ~ThermalMaterial(){};
		
        //! Returns a vector of field IDs corresponding to the variables associated with the material.
        virtual std::vector<int> FieldIds() const = 0;

        //! Initialize the material model.
        virtual void
        initialize(const double dt,
                   const int numOwnedPoints,
                   const int* ownedIDs,
                   const int* neighborhoodList,
                   PeridigmNS::DataManager& dataManager) {}

        virtual void
		computeHeatFlow(const double dt,
						const int numOwnedPoints,
						const int* ownedIDs,
						const int* neighborhoodList, 
						PeridigmNS::DataManager& dataManager) const = 0;
	};

}
#endif // PERIDIGM_BASETHERMALMATERIAL_HPP
