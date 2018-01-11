#ifndef PERIDIGM_THERMALMATERIAL_HPP
#define PERIDIGM_THERMALMATERIAL_HPP

#include "Peridigm_Material.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Peridigm_DataManager.hpp"

namespace PeridigmNS{
	
	class ThermalMaterial {
	public: 
		
// 		! Constructor
		ThermalMaterial() {};
		
// 		! Destructor
		virtual ~ThermalMaterial(){};
		
		virtual void 
		computeHeatFlow(const double dt,
						const int numOwnedPoints,
						const int* ownedIDs,
						const int* neighborhoodList, 
						PeridigmNS::DataManager& dataManager) const = 0;
	};

}
#endif // PERIDIGM_THERMALMATERIAL_HPP
