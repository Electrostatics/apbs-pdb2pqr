cmake_minimum_required (VERSION 3.13)

project(abps_pdb2pqr)

###############################################################################
#  Get the Current Version Number so that it can be embedded in files later
###############################################################################
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION _VERSION_DATA)
string(REPLACE "_" ";" _VERSION_LIST ${_VERSION_DATA})
list(GET _VERSION_LIST 0 MAJOR_VERSION)
list(GET _VERSION_LIST 1 MINOR_VERSION)
list(GET _VERSION_LIST 2 MICRO_VERSION)

###############################################################################
#  From: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
#
#  Use git to add the submodules, initialize them, and update them for the user
###############################################################################
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

add_subdirectory(apbs)