# CMake generated Testfile for 
# Source directory: /home/fabio/local/peridigm/src/materials/unit_test
# Build directory: /home/fabio/local/peridigm/build/src/materials/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(utPeridigm_ElasticMaterial "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_ElasticMaterial")
add_test(utPeridigm_MultiphysicsElasticMaterial "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_MultiphysicsElasticMaterial")
subdirs(twoPoint_SLS_Relaxation)
subdirs(utPeridigm_ElasticPlastic)
