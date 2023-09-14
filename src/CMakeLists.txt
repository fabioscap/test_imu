add_executable(nodo_prova main.cpp)
target_link_libraries(nodo_prova 
                      srrg2_solver_core_library
                      # srrg_geometry is header-only
                      # eigen3 is header-only
)