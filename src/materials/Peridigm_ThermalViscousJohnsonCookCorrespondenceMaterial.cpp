/*! \file Peridigm_ThermalViscousJohnsonCookCorrespondenceMaterial.cpp */

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

#include "Peridigm_ThermalBondBasedMaterial.hpp"
#include "Peridigm_ViscousJohnsonCookCorrespondenceMaterial.hpp"
#include "Peridigm_ThermalViscousJohnsonCookCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::ThermalViscousJohnsonCookCorrespondenceMaterial::ThermalViscousJohnsonCookCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : ViscousJohnsonCookCorrespondenceMaterial(params), ThermalBondBasedMaterial(params)
{
    m_fieldIds = ViscousJohnsonCookCorrespondenceMaterial::m_fieldIds;
    for(int m=0;m<int(ThermalBondBasedMaterial::m_fieldIds.size());m++){
        int thermalFieldId = ThermalBondBasedMaterial::m_fieldIds[m];
        bool found = false;
        for(int n=0;n<int(ViscousJohnsonCookCorrespondenceMaterial::m_fieldIds.size());n++){
            int mechaFieldId = ViscousJohnsonCookCorrespondenceMaterial::m_fieldIds[n];
            if (thermalFieldId==mechaFieldId){
                found = true;
                break;
            }
        }
        if (!found) m_fieldIds.push_back(thermalFieldId);
    }
}

PeridigmNS::ThermalViscousJohnsonCookCorrespondenceMaterial::~ThermalViscousJohnsonCookCorrespondenceMaterial()
{}

void
PeridigmNS::ThermalViscousJohnsonCookCorrespondenceMaterial::initialize(const double dt,
                                        const int numOwnedPoints,
                                        const int* ownedIDs,
                                        const int* neighborhoodList,
                                        PeridigmNS::DataManager& dataManager)
{
    ViscousJohnsonCookCorrespondenceMaterial::initialize(dt,
                                            numOwnedPoints,
                                            ownedIDs,
                                            neighborhoodList,
                                            dataManager);
    ThermalBondBasedMaterial::initialize(dt,
                                         numOwnedPoints,
                                         ownedIDs,
                                         neighborhoodList,
                                         dataManager);    
}
