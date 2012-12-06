ENABLE_LANGUAGE(Fortran)

#SET(APBS_AUX_PATH ${APBS_SRC_PATH}/src/CMakeFiles/apbs.dir)

INCLUDE_DIRECTORIES(. ${APBS_ROOT}/src ${APBS_INST_PATH}/include)
LINK_DIRECTORIES(${APBS_ROOT_PATH}/lib)

#set(LIBS "")
#list(APPEND LIBS "apbs_generic")
#list(APPEND LIBS "apbs_mg")
#list(APPEND LIBS "apbs_pmgc")

#ADD_LIBRARY(iapbs SHARED apbs_driver.c apbs_driver.h)
ADD_LIBRARY(iapbs STATIC apbs_driver.c apbs_driver.h)

ADD_EXECUTABLE(wrapper wrapper.f)

TARGET_LINK_LIBRARIES(wrapper 
 ${CMAKE_BINARY_DIR}/src/CMakeFiles/apbs.dir/routines.c.o
 iapbs
# ${MALOC_LIBRARY}
 apbs_mg
 apbs_generic
 apbs_pmgc)

#${LIBS}

INSTALL(TARGETS iapbs DESTINATION lib)
INSTALL(TARGETS wrapper DESTINATION bin)

