/*! \file Peridigm_ViscousJohnsonCookCorrespondenceMaterial.cpp */

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

#include "Peridigm_JohnsonCookCorrespondenceMaterial.hpp"
#include "Peridigm_ViscousMaxwellCorrespondenceMaterial.hpp"
#include "Peridigm_ViscousJohnsonCookCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::ViscousJohnsonCookCorrespondenceMaterial::ViscousJohnsonCookCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : JohnsonCookCorrespondenceMaterial(params), ViscousMaxwellCorrespondenceMaterial(params)
{
    m_fieldIds = JohnsonCookCorrespondenceMaterial::m_fieldIds;
    for(int m=0;m<int(ViscousMaxwellCorrespondenceMaterial::m_fieldIds.size());m++){
        int viscousFieldId = ViscousMaxwellCorrespondenceMaterial::m_fieldIds[m];
        bool found = false;
        for(int n=0;n<int(JohnsonCookCorrespondenceMaterial::m_fieldIds.size());n++){
            int mechaFieldId = JohnsonCookCorrespondenceMaterial::m_fieldIds[n];
            if (viscousFieldId==mechaFieldId){
                found = true;
                break;
            }
        }
        if (!found) m_fieldIds.push_back(viscousFieldId);
    }
}

PeridigmNS::ViscousJohnsonCookCorrespondenceMaterial::~ViscousJohnsonCookCorrespondenceMaterial()
{}

void
PeridigmNS::ViscousJohnsonCookCorrespondenceMaterial::initialize(const double dt,
                                        const int numOwnedPoints,
                                        const int* ownedIDs,
                                        const int* neighborhoodList,
                                        PeridigmNS::DataManager& dataManager)
{
    JohnsonCookCorrespondenceMaterial::initialize(dt,
                                              numOwnedPoints,
                                              ownedIDs,
                                              neighborhoodList,
                                              dataManager);
    ViscousMaxwellCorrespondenceMaterial::initialize(dt,
                                                   numOwnedPoints,
                                                   ownedIDs,
                                                   neighborhoodList,
                                                   dataManager);    
}

void
PeridigmNS::ViscousJohnsonCookCorrespondenceMaterial::computeForce(const double dt,
                                                                  const int numOwnedPoints,
                                                                  const int* ownedIDs,
                                                                  const int* neighborhoodList,
                                                                  PeridigmNS::DataManager& dataManager) const
{
  JohnsonCookCorrespondenceMaterial::computeForce(dt,
                                              numOwnedPoints,
                                              ownedIDs,
                                              neighborhoodList,
                                              dataManager);
  ViscousMaxwellCorrespondenceMaterial::computeForce(dt,
                                                   numOwnedPoints,
                                                   ownedIDs,
                                                   neighborhoodList,
                                                   dataManager);
}
