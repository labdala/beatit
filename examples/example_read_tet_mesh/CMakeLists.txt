add_executable(example_read_tet_mesh main.cpp)
message("=== Adding Example : example_read_tet_mesh")

set_target_properties(example_read_tet_mesh PROPERTIES  OUTPUT "example_read_tet_mesh.x")


target_link_libraries(example_read_tet_mesh beatit)
target_link_libraries(example_read_tet_mesh ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")
#CMAKE_SOURCE_DIR
#target_include_directories(example_read_tet_mesh ${CMAKE_SOURCE_DIR})

SET_TARGET_PROPERTIES(example_read_tet_mesh PROPERTIES LINKER_LANGUAGE CXX)
