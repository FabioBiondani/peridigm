//! \file Peridigm_JohnsonCookCorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_JOHNSONCOOKCORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_JOHNSONCOOKCORRESPONDENCEMATERIAL_HPP

#include "Peridigm_CorrespondenceMaterial.hpp"
#include <map>

namespace PeridigmNS {

  class JC_CorrespondenceMaterial : public CorrespondenceMaterial{
  public:

	//! Constructor.
    JC_CorrespondenceMaterial(const Teuchos::ParameterList& params);

    //! Destructor.
    virtual ~JC_CorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Johnson-Cook Correspondence"); }

    //! Initialize the derived class
    virtual void initialize(const double dt, 
                            const int numOwnedPoints, 
                            const int* ownedIDs,
                            const int* neighborhoodList,
                            PeridigmNS::DataManager& dataManager);

    //! Evaluate the Cauchy stress.
    virtual void computeCauchyStress(const double dt,
                                     const int numOwnedPoints,
                                     const int* neighborhoodList,
                                     PeridigmNS::DataManager& dataManager) const;
    
//     //! Evaluate the internal force.
//     void computeForce(const double dt,
//                       const int numOwnedPoints,
//                       const int* ownedIDs,
//                       const int* neighborhoodList,
//                       PeridigmNS::DataManager& dataManager) const;



    //! Returns the requested material property
    //! A dummy method here.
    virtual double lookupMaterialProperty(const std::string keyname) const {return 0.0;}

  protected:
    
    // material properties
    double m_MeltingTemperature;
    double m_ReferenceTemperature;
    double m_A;
    double m_N;
    double m_B;
    double m_C;
    double m_M;
    double m_D1;
    double m_D2;
    double m_D3;
    double m_D4;
    double m_D5;
    double m_DC;

    // field spec ids for all relevant data
    int m_unrotatedRateOfDeformationFieldId;
    int m_unrotatedCauchyStressFieldId;
    int m_vonMisesStressFieldId;
    int m_equivalentPlasticStrainFieldId;
    int m_accumulatedPlasticStrainFieldId;
    int m_LocalDamageFieldId;
    int m_bondDamageFieldId;
    int m_deltaTemperatureFieldId;
    int m_DissipationFieldId;
  };
}

#endif // PERIDIGM_JOHNSONCOOKCORRESPONDENCEMATERIAL_HPP
