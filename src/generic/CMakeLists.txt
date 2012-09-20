set(LIBRARY_NAME "apbs_generic")

message(STATUS "Adding ${LIBRARY_NAME}")

set(SOURCES "")
list(APPEND SOURCES nosh.c)
list(APPEND SOURCES mgparm.c)
list(APPEND SOURCES femparm.c)
list(APPEND SOURCES pbeparm.c)
list(APPEND SOURCES apolparm.c)
list(APPEND SOURCES vacc.c)
list(APPEND SOURCES valist.c)
list(APPEND SOURCES vatom.c)
list(APPEND SOURCES vpbe.c)
list(APPEND SOURCES vcap.c)
list(APPEND SOURCES vclist.c)
list(APPEND SOURCES vstring.c)
list(APPEND SOURCES vparam.c)
list(APPEND SOURCES vgreen.c)

message(STATUS "With source files ${SOURCES}")

set(LIB_DEPS "")
#ist(APPEND LIB_DEPS _lib)
list(APPEND LIB_DEPS ${APBS_LIBS})

message(STATUS "With library dependencies ${LIB_DEPS}")

add_library(${LIBRARY_NAME} SHARED ${SOURCES})

target_link_libraries(${LIBRARY_NAME} ${LIB_DEPS})

message(STATUS "Added ${LIBRARY_NAME}")
