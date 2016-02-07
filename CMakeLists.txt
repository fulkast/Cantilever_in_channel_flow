# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( LB2D_Framework_ )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

FIND_PACKAGE( OpenMP)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(${GLEW_INCLUDE_DIRS})
link_libraries(${GLEW_LIBRARIES})
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

SET(BMP thirdparty/EasyBMP/EasyBMP.cpp
             thirdparty/EasyBMP/EasyBMP.h
             thirdparty/EasyBMP/EasyBMP_BMP.h
             thirdparty/EasyBMP/EasyBMP_DataStructures.h
             thirdparty/EasyBMP/EasyBMP_VariousBMPutilities.h)


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11 -DUSE_OPENGL_VISUALIZATION -lgomp -lglut -lGLU -lGL -lGLEW -lX11 ")

#VISFLAG = -DUSE_OPENGL_VISUALIZATION
#LDFLAGS = -lgomp -lglut -lGLU -lGL -lGLEW -lX11 

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

	file(GLOB Dir1_Sources RELATIVE "third_party/EasyBMP/" "*.cpp")

	add_library(BMP third_party/EasyBMP/EasyBMP.cpp)

	target_include_directories (BMP PUBLIC /third_party/EasyBMP)

	add_executable(LB2D main.cpp ${Dir1_Sources} geometry_2D.cpp)

  include( CGAL_CreateSingleSourceCGALProgram )
	
  #create_single_source_cgal_program( "main.cpp" )

	target_link_libraries(LB2D LINK_PUBLIC ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} BMP)


else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

