# CMake generated Testfile for 
# Source directory: /home/fabio/local/peridigm/src/compute/unit_test
# Build directory: /home/fabio/local/peridigm/build/src/compute/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(utPeridigm_Compute_Force "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_Compute_Force")
add_test(utPeridigm_Compute_Force_MPI_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_Compute_Force")
add_test(utPeridigm_Compute_Angular_Momentum "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_Compute_Angular_Momentum")
add_test(utPeridigm_Compute_Angular_Momentum_MPI_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_Compute_Angular_Momentum")
add_test(utPeridigm_Compute_Linear_Momentum "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_Compute_Linear_Momentum")
add_test(utPeridigm_Compute_Linear_Momentum_MPI_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_Compute_Linear_Momentum")
add_test(utPeridigm_Compute_Kinetic_Energy "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utPeridigm_Compute_Kinetic_Energy")
add_test(utPeridigm_Compute_Kinetic_Energy_MPI_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./utPeridigm_Compute_Kinetic_Energy")
