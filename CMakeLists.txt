cmake_minimum_required(VERSION 3.12)
Project("Mesh2SDF")

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries

# SET(CMAKE_BUILD_TYPE Coverage)
SET(CMAKE_BUILD_TYPE Release)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/../bin)

#These flags might not work on every system, especially the release flags, comment out as needed
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

# load Eigen
find_package(Eigen3 3.3 REQUIRED)

# load assimp
message("\n-- Locating assimp...")
find_package(assimp QUIET)
if(assimp_FOUND)
  message("--  assimp located!\n")
else()
  message("--  assimp NOT located!\n")
endif()

add_library(mesh2sdf MeshSDF.cpp util.cpp)
target_link_libraries(mesh2sdf Eigen3::Eigen)
set_target_properties(mesh2sdf PROPERTIES PUBLIC_HEADER "MeshSDF.hpp;array3.hpp;types.hpp;util.hpp")

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} mesh2sdf)

add_executable(QuerySDF mainQuery.cpp)
target_link_libraries(QuerySDF mesh2sdf)

install(TARGETS mesh2sdf EXPORT mesh2sdf-targets
  RUNTIME DESTINATION bin/${PROJECT_NAME}
  LIBRARY DESTINATION lib/${PROJECT_NAME}
  PUBLIC_HEADER DESTINATION include/${PROJECT_NAME})

if (assimp_FOUND)
  target_compile_definitions(${PROJECT_NAME} PRIVATE HAVE_ASSIMP)
  target_compile_definitions(QuerySDF PRIVATE HAVE_ASSIMP)
  target_link_libraries(${PROJECT_NAME} assimp)
  target_link_libraries(QuerySDF assimp)
endif()

if (USE_DOUBLE_PRECISION)
message("\n-- Compiling with double precision...")
  target_compile_definitions(mesh2sdf PRIVATE MESH2SDF_DOUBLE_PRECISION)
  target_compile_definitions(${PROJECT_NAME} PRIVATE MESH2SDF_DOUBLE_PRECISION)
  target_compile_definitions(QuerySDF PRIVATE MESH2SDF_DOUBLE_PRECISION)
endif()