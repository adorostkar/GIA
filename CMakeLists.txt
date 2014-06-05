# top level CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
FIND_PACKAGE(deal.II 8.0 REQUIRED
      HINTS
        ${DEAL_II_DIR} $ENV{DEAL_II_DIR} )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(GIA)

SET(CLEAN_UP_FILES
  *gmv *gnuplot *gpl *eps *pov *vtk *vtu *ucd *.d2 *dat *.log *.m *.1
)


INCLUDE_DIRECTORIES(include)

FILE(GLOB headers ./include/*.h)
FILE(GLOB sources ./source/*.cc)
ADD_EXECUTABLE(elastic ${sources} ${headers})

option(LOG_RUN ON)
IF(LOG_RUN)
    ADD_DEFINITIONS(-DLOGRUN)
ENDIF(LOG_RUN)

ADD_CUSTOM_TARGET(clear
  COMMAND rm ${CLEAN_UP_FILES}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
TARGET_LINK_LIBRARIES(elastic /opt/local/lib/libboost_program_options-mt.a)
DEAL_II_SETUP_TARGET(elastic)
