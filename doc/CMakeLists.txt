################################################################################
# TODO:  Generate license.html
################################################################################

find_package(Doxygen REQUIRED)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/Doxyfile.in
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/Doxyfile
    )

add_custom_target(
    doc
    ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/programmer/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/programmer
    COMMENT "Generating documentation with Doxygen" VERBATIM
)
