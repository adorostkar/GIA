# top level CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
FIND_PACKAGE(deal.II 8.0 REQUIRED
      HINTS
        ${DEAL_II_DIR} $ENV{DEAL_II_DIR} )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(DimitarElastic)

SET(CLEAN_UP_FILES
  *gmv *gnuplot *gpl *eps *pov *vtk *ucd *.d2 data* *.dat *.log *.m
)


INCLUDE_DIRECTORIES(include)

FILE(GLOB headers ./include/*.h)
FILE(GLOB sources ./source/*.cc)
ADD_EXECUTABLE(elastic ${sources} ${headers})
#ADD_EXECUTABLE(elastic ${sources})

DEAL_II_SETUP_TARGET(elastic)
