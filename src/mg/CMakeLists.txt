set(LIBRARY_NAME "apbs_mg")

message(STATUS "Adding ${LIBRARY_NAME}")

set(SOURCES "")
list(APPEND SOURCES vgrid.c)
list(APPEND SOURCES vmgrid.c)
list(APPEND SOURCES vopot.c)
list(APPEND SOURCES vpmg.c)
list(APPEND SOURCES vpmgp.c)

message(STATUS "With source files ${SOURCES}")

set(LIB_DEPS "")
list(APPEND LIB_DEPS apbs_generic)
list(APPEND LIB_DEPS apbs_pmgc)
list(APPEND LIB_DEPS ${APBS_LIBS})

message(STATUS "With library dependencies ${LIB_DEPS}")

add_library(${LIBRARY_NAME} SHARED ${SOURCES})

target_link_libraries(${LIBRARY_NAME} ${LIB_DEPS})

message(STATUS "Added ${LIBRARY_NAME}")
