cmake_minimum_required (VERSION 2.8)
project (HESSIO)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# Allow the user to choose a build type (using ccmake or cmake-gui you'll get a drop-down box with the following options)
set( BuildConfig "CTA_MAX_SC" CACHE STRING
     "Configuration to use for building hessio")

set_property( CACHE BuildConfig PROPERTY STRINGS
              HESS_PHASE_1
              HESS_PHASE_2
              HESS_PHASE_3
              CTA
              CTA_ULTRA
              CTA_ULTRA=5
              CTA_MAX
              CTA_PROD1
              CTA_PROD2
              CTA_PROD3
              CTA_SC
              CTA_SC=1
              CTA_SC=2
              CTA_SC=3
              CTA_SC=4
              CTA_PROD2_SC
              CTA_PROD3_SC
              CTA_MINI2
              CTA_MINI3
              CTA_MAX_SC
              MEGAPIX)

message(STATUS "BuildConfig='${BuildConfig}'")
add_definitions(-D${BuildConfig})

# these should be detected using cmake's config stuff instead of hard-coded:
# add_definitions(-DHAVE_STD_VECTOR -DHAVE_STD_VALARRAY -DHAVE_STD_STRING -DHAVE_64BIT_INT -DSIXTY_FOUR_BITS)

# build the libraries and executables
add_subdirectory("hessioxxx/src")
install(DIRECTORY hessioxxx/include/ DESTINATION include/pyhessio)
