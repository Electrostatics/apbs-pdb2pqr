add_items(
    SOURCES
    cpbconcz2.cpp
    modules.cpp
    pbsolvercz.cpp
    solvationsub.cpp
    subfunction.cpp
    surfconcz.cpp
)

add_items(
    EXTERNAL_HEADERS
    modules.h
    Mat.h
    cpbconcz2.h
)

add_sublibrary(geoflow)
