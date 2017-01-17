cmake_minimum_required(VERSION 2.8)
 
PROJECT(main)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(main MACOSX_BUNDLE main.h main.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(main ${VTK_LIBRARIES})
else()
  target_link_libraries(main vtkHybrid vtkWidgets)
endif()