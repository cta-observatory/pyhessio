cmake_minimum_required (VERSION 2.8)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
include_directories( ${PROJECT_SOURCE_DIR}/hessioxxx/include )

set( HESSIO_HEADERS
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/history.h
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/initial.h
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/io_basic.h 
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/io_histogram.h 
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/io_history.h
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/mc_tel.h
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/straux.h
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/warning.h 
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/io_hess.h
       ${PROJECT_SOURCE_DIR}/hessioxxx/include/io_trgmask.h
       )

# Libraries

add_library( hessio SHARED ${HESSIO_SOURCES} ${HESSIO_HEADERS})

install (TARGETS hessio LIBRARY DESTINATION lib)
