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

message(CHECK_START "Fetching Eigen3")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
# note: To disable eigen tests,
# you should put this code in a add_subdirectory to avoid to change
# BUILD_TESTING for your own project too since variables are directory
# scoped
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
FetchContent_MakeAvailable(Eigen)

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")

# load assimp
message("\n-- Locating assimp...")
find_package(assimp QUIET)
if(assimp_FOUND)
  message("--  assimp located!\n")
else()
  message("--  assimp NOT located!\n")
endif()


add_executable(${PROJECT_NAME} main.cpp  MeshSDF.cpp util.cpp)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

add_executable(QuerySDF mainQuery.cpp MeshSDF.cpp util.cpp)
target_link_libraries(QuerySDF Eigen3::Eigen)

if (assimp_FOUND)
  target_compile_definitions(${PROJECT_NAME} PRIVATE HAVE_ASSIMP)
  target_compile_definitions(QuerySDF PRIVATE HAVE_ASSIMP)
  target_link_libraries(${PROJECT_NAME} assimp)
  target_link_libraries(QuerySDF assimp)
endif()