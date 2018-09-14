//! \file Peridigm_ViscousMaxwellCorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_VISCOUSMAXWELLCORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_VISCOUSMAXWELLCORRESPONDENCEMATERIAL_HPP

#include "Peridigm_BaseViscousCorrespondenceMaterial.hpp"
#include "Peridigm_Material.hpp"

namespace PeridigmNS {

  class ViscousMaxwellCorrespondenceMaterial : public ViscousCorrespondenceMaterial{
  public:

	//! Constructor.
    ViscousMaxwellCorrespondenceMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~ViscousMaxwellCorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Viscous Maxwell Correspondence"); }

    //! Returns a vector of field IDs corresponding to the variables associated with the material.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Evaluate the Cauchy stress.
    virtual void computeCauchyStress(const double dt,
                                     const int numOwnedPoints,
                                     const int* neighborhoodList,
                                     PeridigmNS::DataManager& dataManager) const;

  protected:

    bool m_analysisHasThermal;
    int m_internalVariablesFieldId;
    int m_unrotatedCauchyStressFieldId;
    int m_deltaTemperatureFieldId;

    // material parameters
    Material::BulkMod obj_bulkModulus;
    Material::ShearMod obj_shearModulus;
    Material::TempDepConst obj_lambda;
    Material::TempDepConst obj_tau;
    bool m_temperatureDependence;
  };
}

#endif // PERIDIGM_VISCOUSMAXWELLCORRESPONDENCEMATERIAL_HPP
