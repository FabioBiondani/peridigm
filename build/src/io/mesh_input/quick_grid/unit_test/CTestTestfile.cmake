# CMake generated Testfile for 
# Source directory: /home/fabio/local/peridigm/src/io/mesh_input/quick_grid/unit_test
# Build directory: /home/fabio/local/peridigm/build/src/io/mesh_input/quick_grid/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ut_Q2CylinderRingHorizon "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_Q2CylinderRingHorizon")
add_test(ut_QuickGridHorizon "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_QuickGridHorizon")
add_test(ut_QuickGrid_solidCylinder_np1 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_QuickGrid_solidCylinder_np1")
add_test(ut_QuickGrid_solidCylinder_np2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./ut_QuickGrid_solidCylinder_np2")
add_test(ut_QuickGridMPI_np2_3x3x2 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./ut_QuickGridMPI_np2_3x3x2")
add_test(ut_QuickGridMPI_np3 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "3" "./ut_QuickGridMPI_np3")
add_test(ut_SmallMeshCylinder_np4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "4" "./ut_SmallMeshCylinder_np4")
add_test(ut_SmallMeshCylinder "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_SmallMeshCylinder")
