# CMake generated Testfile for 
# Source directory: /home/fabio/local/peridigm/src/io/pdneigh/unit_test
# Build directory: /home/fabio/local/peridigm/build/src/io/pdneigh/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ut_frameset_2x2x1_np4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "4" "./ut_frameset_2x2x1_np4")
add_test(ut_neighborhood_list "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_neighborhood_list")
add_test(ut_Y-Z_crack "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_Y-Z_crack")
add_test(ut_Y-Z_crack_np4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "4" "./ut_Y-Z_crack_np4")
add_test(ut_Y-Z_crack_jacobian_np4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "4" "./ut_Y-Z_crack_jacobian_np4")
add_test(ut_axisAlignedBoundaryPoints "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_axisAlignedBoundaryPoints")
add_test(ut_NeighborhoodProcessorFrame "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_NeighborhoodProcessorFrame")
add_test(ut_QuickGrid_loadBal_np2_3x1x1 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./ut_QuickGrid_loadBal_np2_3x1x1")
add_test(ut_QuickGrid_loadBal_np2_4x4x4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "2" "./ut_QuickGrid_loadBal_np2_4x4x4")
add_test(ut_QuickGrid_loadBal_np8_4x4x4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "8" "./ut_QuickGrid_loadBal_np8_4x4x4")
add_test(utFinitePlane "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utFinitePlane")
add_test(utReloadBalance_np4 "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "mpiexec" "-np" "4" "./utReloadBalance_np4")
add_test(ut_reLoadBalance "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_reLoadBalance")
add_test(utZoltanQuery_pointSizeInBytes "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utZoltanQuery_pointSizeInBytes")
add_test(utZoltanQuery_packUnpackPointsMultiFunction "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./utZoltanQuery_packUnpackPointsMultiFunction")
add_test(ut_point "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_point")
add_test(ut_kdtree_rectangular_search "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_kdtree_rectangular_search")
add_test(ut_kdtree_spherical_search "python" "/home/fabio/local/peridigm/build/scripts/run_unit_test.py" "./ut_kdtree_spherical_search")
