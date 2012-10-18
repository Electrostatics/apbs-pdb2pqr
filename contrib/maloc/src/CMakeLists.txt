macro(add_items TYPE)

    # Make sure the TYPE argument is valid
    set(ALLOWABLE_TYPES SOURCES EXTERNAL_HEADERS INTERNAL_HEADERS)
    list(FIND ALLOWABLE_TYPES ${TYPE} VALID_TYPE)
    if(VALID_TYPE EQUAL -1)
        message(
            FATAL_ERROR
            "Passed an invalid type to add_items macro.  "
            "Must be one of SOURCES, EXTERNAL_HEADERS, or INTERNAL_HEADERS.  "
            "Was ${TYPE}"
        )
    endif()
    
    # The items to add will be any additional parameters beyond TYPE
    set(ITEMS ${ARGN})
    
    # Get the relative path of the current folder from src
    file(
        RELATIVE_PATH
        PATH
        "${CMAKE_SOURCE_DIR}/src"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    
    # Add each item successively
    foreach(ITEM ${ITEMS})
    
        # If this folder is not src, append the relative path to the item
        if(PATH)
            set(ITEM "${PATH}/${ITEM}")
        endif()
        
        # Add the item to the appropriately typed list
        list(APPEND ${TYPE} "${ITEM}")
    endforeach()
    
    # If this folder is not src, push the updated item list to the parent scope
    if(PATH)
        set(${TYPE} ${${TYPE}} PARENT_SCOPE)
    endif()
    
endmacro()


include_directories(${MALOC_ROOT}/src)
include_directories(${MALOC_ROOT}/src/base)
include_directories(${MALOC_ROOT}/src/psh)
include_directories(${MALOC_ROOT}/src/vsh)
include_directories(${MALOC_ROOT}/src/vsys)

add_subdirectory(base)
add_subdirectory(psh)
add_subdirectory(vsh)
add_subdirectory(vsys)

message(STATUS "Whith sources: ${SOURCES}")
message(STATUS "Whith external headers: ${EXTERNAL_HEADERS}")
message(STATUS "Whith internal headers: ${INTERNAL_HEADERS}")

add_library(maloc SHARED ${SOURCES} ${EXTERNAL_HEADERS} ${INTERNAL_HEADERS})
target_link_libraries(maloc ${MALOC_LIBS})

configure_file(
    ${MALOC_ROOT}/src/config/maloccf.h.in
    ${MALOC_ROOT}/src/maloccf.h
    )
    
INSTALL(FILES ${EXTERNAL_HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR})
INSTALL(TARGETS maloc DESTINATION ${LIBRARY_INSTALL_DIR})

