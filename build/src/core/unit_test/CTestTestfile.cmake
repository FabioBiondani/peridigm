# CMake generated Testfile for 
# Source directory: /home/fabio/local/peridigm/src/core/unit_test
# Build directory: /home/fabio/local/peridigm/build/src/core/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(utPeridigm_State "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_State")
add_test(utPeridigm_State_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_State")
