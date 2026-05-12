# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

project( kNN_VM )


# set (LOCAL_DIR "/home/sli/yduan/local")
# set (CGAL_DIR "${LOCAL_DIR}/CGAL-5.6.2")
#set (GMP_LIBRARIES "${LOCAL_DIR}/gmp/lib")
#set (GMP_INCLUDE_DIR "${LOCAL_DIR}/gmp/include")
#set (MPFR_LIBRARIES "${LOCAL_DIR}/mpfr/lib")
#set (MPFR_INCLUDE_DIR "${LOCAL_DIR}/mpfr/include")
#set (Boost_INCLUDE_DIR "${LOCAL_DIR}/boost_1_87_0")


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
# include_directories(${GMP_INCLUDE_DIR})
# include_directories(${MPFR_INCLUDE_DIR})
# include for local package
# link_directories(${GMP_LIBRARIES})
# link_directories(${MPFR_LIBRARIES})


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


# Creating entries for target: kNN_VM
# ############################

add_executable( kNN_VM  gsd.c io2D.cpp comn.cpp main_conti.cpp vm_voro_2.cpp vm_knn.cpp)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS kNN_VM )

# Link the executable to CGAL and third-party libraries
target_link_libraries(kNN_VM PRIVATE CGAL::CGAL libgmpxx.a libgmp.a)
