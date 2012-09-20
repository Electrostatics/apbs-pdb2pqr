set(LIBRARY_NAME "apbs_pmgc")

message(STATUS "Adding ${LIBRARY_NAME}")

set(SOURCES "")
list(APPEND SOURCES buildAd.c)
list(APPEND SOURCES buildBd.c)
list(APPEND SOURCES buildGd.c)
list(APPEND SOURCES buildPd.c)
list(APPEND SOURCES cgd.c)
list(APPEND SOURCES gsd.c)
list(APPEND SOURCES matvecd.c)
list(APPEND SOURCES mgcsd.c)
list(APPEND SOURCES mgdrvd.c)
list(APPEND SOURCES mgsubd.c)
list(APPEND SOURCES mikpckd.c)
list(APPEND SOURCES mlinpckd.c)
list(APPEND SOURCES mypdec.c)
list(APPEND SOURCES newtond.c)
list(APPEND SOURCES newdrvd.c)
list(APPEND SOURCES powerd.c)
list(APPEND SOURCES smoothd.c)
list(APPEND SOURCES mgfasd.c)

message(STATUS "With source files ${SOURCES}")

set(LIB_DEPS "")
#ist(APPEND LIB_DEPS _lib)
list(APPEND LIB_DEPS ${APBS_LIBS})

message(STATUS "With library dependencies ${LIB_DEPS}")

add_library(${LIBRARY_NAME} SHARED ${SOURCES})

target_link_libraries(${LIBRARY_NAME} ${LIB_DEPS})

message(STATUS "Added ${LIBRARY_NAME}")
