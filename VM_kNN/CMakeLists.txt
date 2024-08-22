# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

project( VM_kNN )


set (CGAL_DIR "/home/yduan/local/CGAL-5.6.1")
set (GMP_LIBRARIES "/home/yduan/local/gmp-6.3.0/lib")
set (GMP_INCLUDE_DIR "/home/yduan/local/gmp-6.3.0/include")
set (MPFR_LIBRARIES "/home/yduan/local/mpfr-4.2.1/lib")
set (MPFR_INCLUDE_DIR "/home/yduan/local/mpfr-4.2.1/include")
set (Boost_INCLUDE_DIR "/home/yduan/local/boost_1_85_0")
# set (Boost_INCLUDE_DIR "/usr/lmp/boost-1.69.0/include")


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()



# include for local directory
include_directories(${GMP_INCLUDE_DIR})
include_directories(${MPFR_INCLUDE_DIR})
# include for local package
link_directories(${GMP_LIBRARIES})
link_directories(${MPFR_LIBRARIES})


option(LINK_STATIC_LIBS
       "Link the executable to static libraries"
       ON)


if (LINK_STATIC_LIBS)
  # link static libraries on unix.
  if (UNIX)
    set(my_flag "-static-libgcc -static-libstdc++ -Wl,-Bstatic -lmpfr -lgmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${my_flag}")
    message(${CMAKE_CXX_FLAGS})
  endif(UNIX)
endif(LINK_STATIC_LIBS)


# Creating entries for target: topoVM
# ############################

add_executable( VM_kNN  gsd.c io2D.cpp comn.cpp main.cpp)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS VM_kNN )

# Link the executable to CGAL and third-party libraries
target_link_libraries(VM_kNN PRIVATE CGAL::CGAL gmp)
