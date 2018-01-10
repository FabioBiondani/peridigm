#ifndef PERIDIGM_PURETHERMALMATERIAL_HPP
#define PERIDIGM_PURETHERMALMATERIAL_HPP

#include <Teuchos_ParameterList.hpp>
#include "Peridigm_DataManager.hpp"


namespace PeridigmNS{
	
	class PureThermalMaterial{
	public: 
		
// 		! Constructor
		PureThermalMaterial(const Teuchos::ParameterList & params) {};
		
// 		! Destructor
		virtual ~PureThermalMaterial(){};
		
		virtual void 
		computeHeatFlow(const double dt,
						const int numOwnedPoints,
						const int* ownedIDs,
						const int* neighborhoodList, 
						PeridigmNS::DataManager& dataManager) const = 0;
	};

}
#endif // PERIDIGM_PURETHERMALMATERIAL_HPP
