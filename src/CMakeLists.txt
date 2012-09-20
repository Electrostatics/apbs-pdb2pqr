add_executable(apbs main.c routines.c)

add_subdirectory(generic)
target_link_libraries(apbs apbs_generic)

add_subdirectory(pmgc)
target_link_libraries(apbs apbs_pmgc)

target_link_libraries(apbs apbs_mg)
add_subdirectory(mg)

if(ENABLE_FETK)
   target_link_libraries(apbs apbs_fem)
   add_subdirectory(fem)
endif(ENABLE_FETK)

target_link_libraries(apbs ${APBS_LIBS})

configure_file(
    ${APBS_ROOT}/src/config/apbscfg.h.in
    ${APBS_ROOT}/src/apbscfg.h
    )
