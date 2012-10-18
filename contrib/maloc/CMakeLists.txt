cmake_minimum_required (VERSION 2.6)

################################################################################
# TODO:
################################################################################

################################################################################
# Future Enchancements
# --------------------
# [ ] Adding a profiling mode set by a flag
# [ ] Add functionality for creating rpm and deb packages
################################################################################


################################################################################
# Test platforms
# --------------
# [ ] Ubuntu x86_64
# [ ] Ubuntu i386
# [ ] Redhat 5
# [ ] ArchLinux x86_64
# [ ] Mac OSX
# [ ] Windows 7 x86_64
# [ ] Windows 7 i386
################################################################################

include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckSymbolExists)
include(CheckTypeSize)
include(ExternalProject)

set(MALOC_VERSION "1.0.1")
set(FETK_VERSION "1.0.0")
set(PACKAGE_NAME "maloc")
set(PACKAGE_TARNAME "maloc")
set(PACKAGE_VERSION "1.0")
set(PACKAGE_STRING "maloc 1.0")
set(PACKAGE_BUGREPORT "mholst@math.ucsd.edu")

project(maloc)



################################################################################
# Set project paths                                                            #
################################################################################

message(STATUS "Setting project paths")

set(MALOC_ROOT ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${MALOC_ROOT}/bin)
set(LIBRARY_OUTPUT_PATH ${MALOC_ROOT}/lib)
#set(MALOC_BINARY ${EXECUTABLE_OUTPUT_PATH}/maloc)


set(INCLUDE_INSTALL_DIR include/maloc)
set(LIBRARY_INSTALL_DIR lib)
set(EXECUTABLE_INSTALL_DIR bin)



################################################################################
# Set up temporary files and directories                                       #
################################################################################

file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/temp)



################################################################################
# Set the lookup paths for external library dependencies                       #
################################################################################

message(STATUS "Setting lookup paths for libraries")

set(SYS_LIBPATHS "")
list(APPEND SYS_LIBPATHS /usr/lib)
list(APPEND SYS_LIBPATHS /usr/local/lib)
list(APPEND SYS_LIBPATHS /lib/x86_64-linux-gnu)

set(MALOC_LIBS "")



################################################################################
# Enable ansi pedantic compiling                                               #
################################################################################

option(ENABLE_PEDANTIC "Enable the pedantic ansi compilation" OFF)

if(ENABLE_PEDANTIC)
    ADD_DEFINITIONS("-Wall -pedantic -ansi")
    message(STATUS "Pedantic compilation enabled")
endif(ENABLE_PEDANTIC)



################################################################################
# Handle the math library dependency                                           #
################################################################################

message(STATUS "Checking for math library")

find_library(
    MATH_LIBRARY
    NAMES m
    PATHS ${SYS_LIBPATHS}
    DOC   "The math library"
    )

if(MATH_LIBRARY)
    list(APPEND MALOC_LIBS ${MATH_LIBRARY})
    message(STATUS "The math library was found: ${MATH_LIBRARY}")
else(MATH_LIBRARY)
    message(FATAL_ERROR "The math library (libm) is required for apbs")
endif(MATH_LIBRARY)




################################################################################
# Handle the liberty library dependency                                        #
################################################################################

message(STATUS "Checking for the liberty library")

find_library(
    LIBERTY_LIBRARY
    NAMES iberty
    PATHS ${SYS_LIBPATHS}
    DOC   "The liberty library"
    )

if(LIBERTY_LIBRARY)
    list(APPEND MALOC_LIBS ${LIBERTY_LIBRARY})
    message(STATUS "The liberty library was found: ${LIBERTY_LIBRARY}")
else(LIBERTY_LIBRARY)
    message(WARNING "Liberty library not found")
endif(LIBERTY_LIBRARY)

# TODO: Determine if the liberty library is ever needed any more



################################################################################
# Handle the nsl library dependency                                            #
################################################################################

find_library(
    NSL_LIBRARY
    NAMES nsl
    PATHS ${SYS_LIBPATHS}
    DOC   "The nsl library"
    )

if(NSL_LIBRARY)
    list(APPEND MALOC_LIBS ${NSL_LIBRARY})
    message(STATUS "The nsl library was found: ${NSL_LIBRARY}")
else(NSL_LIBRARY)
    message(WARNING "The nsl library was not found")
endif(NSL_LIBRARY)

# TODO: Determine if nsl library is ever needed any more



################################################################################
# Handle the socket library dependency                                         #
################################################################################

find_library(
    SOCKET_LIBRARY
    NAMES socket
    PATHS ${SYS_LIBPATHS}
    DOC   "The socket library"
    )

if(SOCKET_LIBRARY)
    list(APPEND MALOC_LIBS ${SOCKET_LIBRARY})
    message(STATUS "The socket library was found: ${SOCKET_LIBRARY}")
else(SOCKET_LIBRARY)
    message(WARNING "The socket library was not found")
endif(SOCKET_LIBRARY)

# TODO: Determine if socket library is ever needed any more



################################################################################
# Handle the thread library dependency                                         #
################################################################################

find_library(
    THREAD_LIBRARY
    NAMES thread
    PATHS ${SYS_LIBPATHS}
    DOC   "The thread library"
    )

if(THREAD_LIBRARY)
    list(APPEND MALOC_LIBS ${THREAD_LIBRARY})
    message(STATUS "The thread library was found: ${THREAD_LIBRARY}")
else(THREAD_LIBRARY)
    message(WARNING "The thread library was not found")
endif(THREAD_LIBRARY)

# TODO: Determine if thread library is ever needed any more



################################################################################
# Handle conditional inclusion of the readline library                         #
################################################################################

option(ENABLE_READLINE "Enable the readline library" OFF)

if(ENABLE_READLINE)
    message(STATUS "Checking for realine library")
    find_library(
        READLINE_LIBRARY
        NAMES readline
        PATHS ${SYS_LIBPATHS}
        )

    if(READLINE_LIBRARY)

        CHECK_INCLUDE_FILES("readline/readline.h" HAVE_READLINE_READLINE_H)
        CHECK_INCLUDE_FILES("readline/history.h" HAVE_READLINE_HISTORY_H)

        list(APPEND MALOC_LIBS ${READLINE_LIBRARY})
        set(HAVE_LIBREADLINE 1)
        message(STATUS "The readline library was found: ${READLINE_LIBRARY}")



        message(STATUS "Checking fore ncurses library")
        find_library(
            NCURSES_LIBRARY
            NAMES ncurses
            PATHS ${SYS_LIBPATHS}
            )

        if(NCURSES_LIBRARY)
            list(APPEND MALOC_LIBS ${NCURSES_LIBRARY})
            message(STATUS "The ncurses library was found: ${NCURSES_LIBRARY}")
        else(NCURSES_LIBRARY)
            message(WARNING "The ncurses library was not found")
        endif(NCURSES_LIBRARY)
    else(READLINE_LIBRARY)
        message(WARNING "The readline library was not found.  Did not enable readline")
    endif(READLINE_LIBRARY)
endif(ENABLE_READLINE)



################################################################################
# Handle conditional inclusion of the electric fence library                   #
################################################################################

option(ENABLE_EFENCE "Enable the electric fence library" OFF)

if(ENABLE_EFENCE)
    message(STATUS "Checking for efence library")
    find_library(EFENCE_LIBRARY NAMES efence PATHS ${SYS_LIBPATHS})
    if(EFENCE_LIBRARY)
        list(APPEND MALOC_LIBS ${EFENCE_LIBRARY})
        message(STATUS "The efence library was found: ${EFENCE_LIBRARY}")
    else(EFENCE_LIBRARY)
        message(WARNING "The efence library was not found. Did not enable efence")
    endif(EFENCE_LIBRARY)
endif(ENABLE_EFENCE)



################################################################################
# Check for required header files                                              #
################################################################################


set(CMAKE_EXTRA_INCLUDE_FILES stddef.h sys/stat.h sys/types.h)
CHECK_TYPE_SIZE("mode_t" MODE_T)
CHECK_TYPE_SIZE("pid_t" PID_T)
CHECK_TYPE_SIZE("size_t" SIZE_T)
set(CMAKE_EXTRA_INCLUDE_FILES)

if(NOT HAVE_MODE_T)
    set(mode_t "int")
endif(NOT HAVE_MODE_T)

if(NOT HAVE_PID_T)
    set(pid_t "int")
endif(NOT HAVE_PID_T)

if(NOT HAVE_SIZE_T)
    set(size_t "unsigned int")
endif(NOT HAVE_SIZE_T)

CHECK_INCLUDE_FILES("stdlib.h;stdarg.h;string.h;float.h" STDC_HEADERS)

if(NOT CYGWIN)

    CHECK_INCLUDE_FILES("unistd.h" HAVE_UNISTD_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required unistd.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("sys/stat.h" HAVE_SYS_STAT_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required sys/stat.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("sys/socket.h" HAVE_SYS_SOCKET_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required sys/socket.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("sys/time.h" HAVE_SYS_TIME_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required sys/time.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("sys/types.h" HAVE_SYS_TYPES_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required sys/types.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("sys/un.h" HAVE_SYS_UN_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required sys/un.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("sys/wait.h" HAVE_SYS_WAIT_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required sys/wait.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("fcntl.h" HAVE_FCNTL_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required rcntl.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("rpc/rpc.h" HAVE_RPC_RPC_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required rpc/rpc.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("netinet/in.h" HAVE_NETINET_IN_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required netinet/in.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("arpa/inet.h" HAVE_ARPA_INET_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required arpa/inet.h header")
    endif(NOT HAVE_UNISTD_H)

    CHECK_INCLUDE_FILES("netdb.h" HAVE_NETDB_H)
    if(NOT HAVE_UNISTD_H)
        message(FATAL_ERROR "Missing required netdb.h header")
    endif(NOT HAVE_UNISTD_H)
    list(APPEND MALOC_REQUIRED_HEADERS "netdb.h")

endif(NOT CYGWIN)



################################################################################
# Check for a few required functions and symbols                               #
################################################################################

CHECK_FUNCTION_EXISTS(xdrmem_create HAVE_XDR)

if(NOT HAVE_XDR)
    message(WARNING "The xdrmem_create function was not found")
endif(NOT HAVE_XDR)



CHECK_FUNCTION_EXISTS(getcwd HAVE_GETCWD)

if(NOT HAVE_GETCWD)
    message(WARNING "The getcwd function was not found")
endif(NOT HAVE_GETCWD)



CHECK_SYMBOL_EXISTS(O_NONBLOCK "fcntl.h" HAVE_O_NONBLOCK)

if(NOT HAVE_O_NONBLOCK)
    message(WARNING "The O_NONBLOCK symbol was not found")
endif(NOT HAVE_O_NONBLOCK)



################################################################################
# Handle conditional availability of macro embedding                           #
################################################################################

try_compile(
    HAVE_EMBED
    ${MALOC_ROOT}/build
    ${MALOC_ROOT}/src/config/embed_test.c
    )

# TODO: Determine if the EMBED macro is even used



################################################################################
# Handle conditional debug building                                            #
################################################################################

option(ENABLE_DEBUG "Enable debugging mode" OFF)

if(ENABLE_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    set(DEBUG 1)
    set(HAVE_DEBUG 1)
    message(STATUS "Debugging compilation enabled")
endif(ENABLE_DEBUG)



################################################################################
# Handle conditional building with MPI Support                                 #
################################################################################

option(ENABLE_MPI "Enable MPI parallelism" ON)

if(ENABLE_MPI)
    if(NOT ENABLE_DEBUG)
        message(STATUS "Checking for MPI")
        find_package(MPI)
        if(MPI_FOUND)
            message(STATUS "MPI support enabled")
            include_directories(${MPI_INCLUDE_PATH})
            list(APPEND MALOC_LIBS ${MPI_LIBRARIES})
        else(MPI_FOUND)
            message(WARNING "MPI was not found; disabling")
        endif(MPI_FOUND)
    else(NOT ENABLE_DEBUG)
        message(WARNING "MPI may not be enabled in debugging mode")
    endif(NOT ENABLE_DEBUG)
endif(ENABLE_MPI)



################################################################################
# Handle library checks for embedded unix environments in windows              #
################################################################################

if(MINGW)
    message(STATUS "Checking for wsock32 in MinGW environment")
    find_library(
        MINGW_WSOCK32
        NAMES wsock32
        PATHS ${SYS_LIBPATHS}
        DOC   "The wsock32 library"
        )

    if(MINGW_WSOCK32)
        message(STATUS "The wsock32 library was found: ${MINGW_HAS_WSOCK32}")
    else(MINGW_WSOCK32)
        message(FATAL_ERROR "The wsock32 library was not fond")
    endif(MINGW_WSOCK32)
endif(MINGW)



if(CYGWIN)
    message(STATUS "Checking fore rpclib in Cygwin environment")
    find_library(
        CYGWIN_RPCLIB
        NAMES rpclib
        PATHS ${SYS_LIBPATHS}
        DOC   "The rpclib library"
        )

    if(CYGWIN_RPCLIB)
        message(STATUS "The rpclib library was found: ${CYGWIN_HAS_RPCLIB}")
    else(CYGWIN_RPCLIB)
        message(FATAL_ERROR "The rpclib library was not found")
    endif(NOT CYGWIN_RPCLIB)



    message(STATUS "Checking for wsock32 in Cygwin environment")
    find_library(
        CYGWIN_WSOCK32
        NAMES wsock32
        PATHS ${SYS_LIBPATHS}
        DOC   "The wsock32 library"
        )

    if(CYGWIN_WSOCK32)
        message(STATUS "The wsock32 library was found: ${CYGWIN_WSOCK32}")
        list(APPPEND MALOC_LIBS ${CYGWIN_WSOCK32})
    else(CYGWIN_WSOCK32)
        message(FATAL_ERROR "The wsock32 library was not fond")
    endif(CYGWIN_WSOCK32)

    set(MSDOS 1)
    set(WIN32 1)
endif(CYGWIN)


################################################################################
# Build MALOCsources                                                           #
################################################################################

add_subdirectory(src)



################################################################################
# Clean up temporary files and directories                                     #
################################################################################

file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}/temp)

