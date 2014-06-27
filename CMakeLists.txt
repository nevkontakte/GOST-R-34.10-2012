project(gost_ecc)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(src SRC_LIST)
include_directories(include)

add_library(${PROJECT_NAME} SHARED ${SRC_LIST})

include(ExternalProject)

ExternalProject_Add(signature_contest
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/signature_contest
    CONFIGURE_COMMAND ./configure
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)
