project(gost_ecc)
cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)
add_definitions(-Wall -Werror -Wpedantic)

aux_source_directory(src SRC_LIST)
aux_source_directory(test TEST_SRC_LIST)
aux_source_directory(production PRODUCTION_SRC_LIST)

include_directories(include)

find_library(CRYPTOPP_LIBRARY cryptopp)

add_library(${PROJECT_NAME} SHARED ${SRC_LIST} ${PRODUCTION_SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${CRYPTOPP_LIBRARY})

add_executable(${PROJECT_NAME}_test ${SRC_LIST} ${TEST_SRC_LIST})
target_link_libraries(${PROJECT_NAME}_test ${CRYPTOPP_LIBRARY})

include(ExternalProject)

ExternalProject_Add(signature_contest
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/signature_contest
    CONFIGURE_COMMAND ./configure
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)

file(GLOB_RECURSE HEADER_FILES ${PROJECT_SOURCE_DIR}/include/*.h) # Make QtCreator to show headers
file(GLOB_RECURSE VALIDATOR_FILES ${PROJECT_SOURCE_DIR}/ext/signature_contest/**) # Make QtCreator to show headers
add_custom_target(${PROJECT_NAME}_aux SOURCES ${HEADER_FILES} ${VALIDATOR_FILES})
