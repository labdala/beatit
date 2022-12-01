SET(TESTNAMEV example_laryssa_calc_mesh_length)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_laryssa_calc_mesh_length")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} /nas/longleaf/home/laryssa/dogwood/sfw/linux/libmesh/1.7.0/1.7.0-opt/lib/libtimpi_opt.so.3)

#link_directories(${LIBMESH_DIR}/lib)


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)


SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )
