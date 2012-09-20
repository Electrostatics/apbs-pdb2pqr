set(LIBRARY_NAME "apbs_fem")

message(STATUS "Adding ${LIBRARY_NAME}")

set(SOURCES "")
list(APPEND SOURCES vfetk.c)
list(APPEND SOURCES vcsm.c)
list(APPEND SOURCES vpee.c)

message(STATUS "With source files ${SOURCES}")

set(LIB_DEPS "")
#ist(APPEND LIB_DEPS _lib)
list(APPEND LIB_DEPS ${APBS_LIBS})

message(STATUS "With library dependencies ${LIB_DEPS}")

add_library(${LIBRARY_NAME} SHARED ${SOURCES})

target_link_libraries(${LIBRARY_NAME} ${LIB_DEPS})

message(STATUS "Added ${LIBRARY_NAME}")
