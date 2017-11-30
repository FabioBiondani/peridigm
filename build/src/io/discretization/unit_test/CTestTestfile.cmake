# CMake generated Testfile for 
# Source directory: /home/fabio/local/peridigm/src/io/discretization/unit_test
# Build directory: /home/fabio/local/peridigm/build/src/io/discretization/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(utPeridigm_PdQuickGridDiscretization "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_PdQuickGridDiscretization")
add_test(utPeridigm_PdQuickGridDiscretization_MPI_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_PdQuickGridDiscretization_MPI_np2")
add_test(utPeridigm_ExodusDiscretization "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_ExodusDiscretization")
add_test(utPeridigm_ExodusDiscretization_MPI_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_ExodusDiscretization")
add_test(utPeridigm_GeometryUtils "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_GeometryUtils")
