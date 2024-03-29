add_subdirectory (unit_test)

# Collection of libs that we need to link with
set(PD_MATERIAL_SOURCES
    ../core/Peridigm_InfluenceFunction.cpp
    elastic.cxx
    elastic_bond_based.cxx
    elastic_plastic.cxx
    elastic_plastic_hardening.cxx
    viscoelastic.cxx
    correspondence.cxx
    elastic_correspondence.cxx
    elastic_plastic_correspondence.cxx
    jc_correspondence.cxx
    isotropic_hardening_correspondence.cxx
    viscoplastic_needleman_correspondence.cxx
    material_utilities.cxx
    nonlocal_diffusion.cxx
    nonlocal_thermal_diffusion.cxx
    thermal_bondbased.cxx
    pals.cxx
    damagepals.cxx
    jc_pals.cxx
    jc_ordinary.cxx
    JohnsonCook.cxx
    viscousmaxwell_correspondence.cxx
    viscousmaxwell_ordinary.cxx
    viscous_bond_based.cxx
)

# Optional source files for Sandia internal development
IF(PERIDIGM_SANDIA_INTERNAL)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} partial_stress.cxx)
ENDIF()
IF(PERIDIGM_PV)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} elastic_pv.cxx correspondence_pv.cxx linear_lps_pv.cxx)
ENDIF()
IF(PERIDIGM_KOKKOS)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} elastic_kokkos.cxx)
ENDIF()

# Optional source files for CJL development
IF(PERIDIGM_CJL)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} lammi_concrete_model.cxx)
ENDIF()

# PdMaterialUtilities library
add_library(PdMaterialUtilities ${PD_MATERIAL_SOURCES})

IF (INSTALL_PERIDIGM)
   install(TARGETS PdMaterialUtilities EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()
