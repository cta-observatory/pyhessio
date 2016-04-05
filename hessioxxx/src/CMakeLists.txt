set( HESSIO_SOURCES
                    atmprof.c
                    current.c
                    dhsort.c
                    eventio.c
                    eventio_registry.c
                    fileopen.c
                    histogram.c
                    hconfig.c
                    moments.c
                    io_histogram.c 
                    io_history.c 
                    io_simtel.c
                    io_trgmask.c
                    straux.c 
                    warning.c 
                    io_hess.c)
include_directories( ${PROJECT_SOURCE_DIR}/include )

set( HESSIO_HEADERS
       ${PROJECT_SOURCE_DIR}/include/history.h
       ${PROJECT_SOURCE_DIR}/include/initial.h
       ${PROJECT_SOURCE_DIR}/include/io_basic.h 
       ${PROJECT_SOURCE_DIR}/include/io_histogram.h 
       ${PROJECT_SOURCE_DIR}/include/io_history.h
       ${PROJECT_SOURCE_DIR}/include/mc_tel.h
       ${PROJECT_SOURCE_DIR}/include/straux.h
       ${PROJECT_SOURCE_DIR}/include/warning.h 
       ${PROJECT_SOURCE_DIR}/include/io_hess.h
       ${PROJECT_SOURCE_DIR}/include/io_trgmask.h
       )

# Libraries

add_library( hessio SHARED ${HESSIO_SOURCES} ${HESSIO_HEADERS})

add_library( hessio++ SHARED
             EventIO.cc 
             ${PROJECT_SOURCE_DIR}/include/EventIO.hh
             ${HESSIO_SOURCES} ${HESSIO_INCLUDES} )


# C Executables

add_executable( testio testio.c )
target_link_libraries( testio hessio m )

add_executable( read_hess read_hess.c rec_tools.c user_analysis.c reconstruct.c  camera_image.c basic_ntuple.c)
target_link_libraries( read_hess hessio m )

add_executable( read_hess_nr read_hess_nr.c rec_tools.c camera_image.c )
target_link_libraries( read_hess_nr hessio m )

add_executable( split_hessio split_hessio.c )
target_link_libraries( split_hessio hessio m )

add_executable( check_trgmask check_trgmask.c )
target_link_libraries( check_trgmask hessio m )

add_executable( gen_trgmask gen_trgmask.c )
target_link_libraries( gen_trgmask hessio m )

add_executable( extract_hess extract_hess.c )
target_link_libraries( extract_hess hessio m )

add_executable( merge_simtel merge_simtel.c )
target_link_libraries( merge_simtel hessio m )

add_executable( list_histograms list_histograms.c )
target_link_libraries( list_histograms hessio m )

add_executable( add_histograms add_histograms.c )
target_link_libraries( add_histograms hessio m )

add_executable( fcat fcat.c )
target_link_libraries( fcat hessio m )


# C++ executables

add_executable( testio++ TestIO.cc )
target_link_libraries( testio++ hessio++ )

add_executable( statio statio.cc )
target_link_libraries( statio hessio++ )

add_executable( filterio filterio.cc )
target_link_libraries( filterio hessio++ )

