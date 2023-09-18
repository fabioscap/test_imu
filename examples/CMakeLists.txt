add_executable(integrate_trajectory integrate_trajectory.cpp)
target_link_libraries(integrate_trajectory 
                      srrg2_solver_core_library
                      # srrg_geometry is header-only
                      # eigen3 is header-only
                      synthetic
                      common
                      variables_and_factors
)

add_executable(test_ test.cpp)
target_link_libraries(test_ 
                      srrg2_solver_core_library
                      # srrg_geometry is header-only
                      # eigen3 is header-only
                      synthetic
                      common
                      variables_and_factors
)