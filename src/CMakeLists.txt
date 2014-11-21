set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

set(lrvTrack_SRCS_CPP lrvTrackOL.cpp
blobUtils.cpp 
lrvTrackBase.cpp
lrvTrackDebug.cpp
larvaDistanceMap.cpp
larvaObject.cpp
lrvTrackFit.cpp)

add_executable(lrvTrack ${lrvTrack_SRCS_CPP})

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_BINARY_DIR}")

if(Boost_USE_MULTITHREADED)
  find_package(Boost REQUIRED COMPONENTS log log_setup program_options filesystem system timer thread)
else(Boost_USE_MULTITHREADED)
  find_package(Boost REQUIRED COMPONENTS log log_setup program_options filesystem system timer)
endif(Boost_USE_MULTITHREADED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(lrvTrack ${Boost_LIBRARIES})
endif(Boost_FOUND)


# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/lrvTrackConfig.h.in"
  "${PROJECT_BINARY_DIR}/lrvTrackConfig.h"
  )

find_package( OpenCV REQUIRED COMPONENTS core highgui videoio video imgcodecs imgproc ml)
if(OpenCV_FOUND)
  include_directories(${OpenCV_INCLUDE_DIRS})
  link_directories(${OpenCV_LIBRARY_DIRS})
  target_link_libraries(lrvTrack ${OpenCV_LIBRARIES})
endif()

set(CvBlob_LIBRARIES cvblob)
find_package (CvBlob REQUIRED)
if (CvBlob_FOUND)
  MESSAGE("-- Found CVBlob: INCLUDE: ${CvBlob_INCLUDE_DIRS}, LIB: ${CvBlob_LIBRARIES}")
  include_directories(${CvBlob_INCLUDE_DIRS})
  target_link_libraries (lrvTrack ${CvBlob_LIBRARIES})
endif (CvBlob_FOUND)

find_package( TBB REQUIRED )
# If you want FindTBB to fail if TBB is not found, use this form:
# find_package(TBB REQUIRED)
if(NOT TBB_FOUND)
  MESSAGE (STATUS "TBB not found. Project will not be built.")
else(NOT TBB_FOUND)
  # Here goes the actual project
  include_directories(${TBB_INCLUDE_DIRS})
  link_directories(${TBB_LIBRARY_DIRS})
  target_link_libraries (lrvTrack ${TBB_LIBRARIES})
endif(NOT TBB_FOUND)

add_library(alglib alglib/alglibinternal.cpp
alglib/alglibmisc.cpp
alglib/ap.cpp
alglib/integration.cpp
alglib/interpolation.cpp
alglib/linalg.cpp
alglib/optimization.cpp
alglib/solvers.cpp
alglib/specialfunctions.cpp
alglib/statistics.cpp
)
target_link_libraries (lrvTrack alglib)
