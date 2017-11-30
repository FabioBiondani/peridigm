# CMake generated Testfile for 
# Source directory: /home/fabio/local/peridigm/src/io/unit_test
# Build directory: /home/fabio/local/peridigm/build/src/io/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(utPeridigm_ProximitySearch_np1 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_ProximitySearch")
add_test(utPeridigm_ProximitySearch_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_ProximitySearch")
add_test(utPeridigm_ProximitySearch_np3 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "3" "./utPeridigm_ProximitySearch")
add_test(utPeridigm_ProximitySearch_np4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "4" "./utPeridigm_ProximitySearch")
add_test(utPeridigm_ProximitySearch_np5 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "5" "./utPeridigm_ProximitySearch")
add_test(utPeridigm_SearchTree "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_SearchTree")
add_test(ut_kdtree_nn_search "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_kdtree_nn_search")
