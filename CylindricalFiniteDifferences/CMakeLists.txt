cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_COMPILER_FORCED TRUE)
set(CMAKE_CXX_COMPILER_FORCED TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
project(kevin CXX)

INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/include
  )

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/unsupported/)
include_directories(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/supported/)
INCLUDE_DIRECTORIES(
  "${CMAKE_SOURCE_DIR}/hazel/src/"
)
find_package(LAPACK)
find_package(BLAS)

if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()


find_package(OpenMP)

ADD_LIBRARY(node ${CMAKE_SOURCE_DIR}/src/Node.cpp)
ADD_LIBRARY(element ${CMAKE_SOURCE_DIR}/src/Element.cpp)
ADD_LIBRARY(reference ${CMAKE_SOURCE_DIR}/src/reference_configuration.cpp)
ADD_LIBRARY(material ${CMAKE_SOURCE_DIR}/src/material_class.cpp)
ADD_LIBRARY(lc ${CMAKE_SOURCE_DIR}/src/LocalConstructor.cpp)
ADD_LIBRARY(elasticproblem ${CMAKE_SOURCE_DIR}/src/ElasticProblem.cpp)
ADD_LIBRARY(dofhandler ${CMAKE_SOURCE_DIR}/src/DOFHandler.cpp)
add_library(hvector ${CMAKE_SOURCE_DIR}/hazel/src/Vector.cpp)
add_library(hmatrix ${CMAKE_SOURCE_DIR}/hazel/src/Matrix.cpp)



ADD_EXECUTABLE(run ${CMAKE_SOURCE_DIR}/src/main.cpp)
if(OpenMP_CXX_FOUND)
    
    
endif()

TARGET_LINK_LIBRARIES(hvector ${lapackblas_libraries})
TARGET_LINK_LIBRARIES(hmatrix ${lapackblas_libraries})
TARGET_LINK_LIBRARIES(hvector hmatrix)
TARGET_LINK_LIBRARIES(hmatrix hvector)

TARGET_LINK_LIBRARIES(element node)
TARGET_LINK_LIBRARIES(element reference)
TARGET_LINK_LIBRARIES(lc element)
TARGET_LINK_LIBRARIES(lc material)
TARGET_LINK_LIBRARIES(elasticproblem node)
TARGET_LINK_LIBRARIES(elasticproblem element)
TARGET_LINK_LIBRARIES(elasticproblem reference)
TARGET_LINK_LIBRARIES(elasticproblem material)
TARGET_LINK_LIBRARIES(elasticproblem lc)

TARGET_LINK_LIBRARIES(dofhandler node)
TARGET_LINK_LIBRARIES(dofhandler reference)
TARGET_LINK_LIBRARIES(dofhandler element)
TARGET_LINK_LIBRARIES(dofhandler material)
TARGET_LINK_LIBRARIES(elasticproblem dofhandler)

TARGET_LINK_LIBRARIES(elasticproblem hvector)
TARGET_LINK_LIBRARIES(elasticproblem hmatrix)

TARGET_LINK_LIBRARIES(run node)
TARGET_LINK_LIBRARIES(run element)
TARGET_LINK_LIBRARIES(run reference)
TARGET_LINK_LIBRARIES(run material)
TARGET_LINK_LIBRARIES(run elasticproblem)





# Define a distclean target to remove every generated file:
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
  COMMENT "distclean invoked"
  )