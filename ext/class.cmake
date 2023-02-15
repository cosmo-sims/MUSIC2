cmake_minimum_required(VERSION 3.11)
include(FetchContent)
FetchContent_Declare(
    class
    GIT_REPOSITORY https://github.com/ohahn/class_public.git
    GIT_TAG master
    GIT_SHALLOW YES
    GIT_PROGRESS TRUE
    USES_TERMINAL_DOWNLOAD TRUE   # <---- this is needed only for Ninja
)

FetchContent_GetProperties(class)
if(NOT class_POPULATED)
    set(FETCHCONTENT_QUIET OFF)
    FetchContent_Populate(class)
    add_subdirectory(${class_SOURCE_DIR} ${class_BINARY_DIR})
endif()
