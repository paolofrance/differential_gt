cmake_minimum_required(VERSION 3.1)

project(differential_gt)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
roscpp
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)


catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ${PROJECT_NAME}
 CATKIN_DEPENDS 
 DEPENDS EIGEN3
)

include_directories( include ${catkin_INCLUDE_DIRS} )



add_library(
cgt src/differential_gt/cgt.cpp
)
target_link_libraries(cgt  ${catkin_LIBRARIES} Eigen3::Eigen)

add_library(
ncgt src/differential_gt/ncgt.cpp
)
target_link_libraries(ncgt  ${catkin_LIBRARIES} Eigen3::Eigen)

install(
  TARGETS cgt
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)
install(
  DIRECTORY include/cgt/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

install(
  TARGETS ncgt
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)
install(
  DIRECTORY include/ncgt/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

add_executable(dgt_node src/differential_gt/dgt_node.cpp)
target_link_libraries(dgt_node cgt ncgt ${catkin_LIBRARIES})


