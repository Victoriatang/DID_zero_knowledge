cmake_minimum_required(VERSION 3.29)
project(BBSDID LANGUAGES CXX)

option(BBSAC_BUILD_TESTS "Build tests" ON)
option(BBSAC_WARNINGS_AS_ERRORS "Treat warnings as errors for project targets" OFF)
option(BBSAC_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(BBSAC_MCL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/mcl")
set(BBSAC_BICYCL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/bicycl")
if(NOT EXISTS "${BBSAC_MCL_DIR}/CMakeLists.txt" OR NOT EXISTS "${BBSAC_BICYCL_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "Could not find mcl/bicycl dependencies under include/.")
endif()

add_subdirectory("${BBSAC_MCL_DIR}" "${CMAKE_BINARY_DIR}/deps/mcl" EXCLUDE_FROM_ALL)
add_subdirectory("${BBSAC_BICYCL_DIR}" "${CMAKE_BINARY_DIR}/deps/bicycl" EXCLUDE_FROM_ALL)

get_target_property(MCL_INTERFACE_INCLUDES mcl INTERFACE_INCLUDE_DIRECTORIES)
if(MCL_INTERFACE_INCLUDES)
    set_property(TARGET mcl PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${MCL_INTERFACE_INCLUDES}")
endif()
get_target_property(BICYCL_INTERFACE_INCLUDES bicycl INTERFACE_INCLUDE_DIRECTORIES)
if(BICYCL_INTERFACE_INCLUDES)
    set_property(TARGET bicycl PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${BICYCL_INTERFACE_INCLUDES}")
endif()

find_library(GMP_LIBRARY gmp REQUIRED)
find_library(GMPXX_LIBRARY gmpxx REQUIRED)
find_package(OpenSSL REQUIRED)

add_library(bbsac_project_options INTERFACE)
target_compile_features(bbsac_project_options INTERFACE cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bbsac_project_options INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
    )
    if(BBSAC_WARNINGS_AS_ERRORS)
        target_compile_options(bbsac_project_options INTERFACE -Werror)
    endif()
    target_compile_options(bbsac_project_options INTERFACE
        $<$<CONFIG:Debug>:-O0 -g3 -fno-omit-frame-pointer>
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
        $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
    )
    if(BBSAC_ENABLE_SANITIZERS)
        target_compile_options(bbsac_project_options INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(bbsac_project_options INTERFACE -fsanitize=address,undefined)
    endif()
elseif(MSVC)
    target_compile_options(bbsac_project_options INTERFACE /W4 /permissive- /EHsc)
    if(BBSAC_WARNINGS_AS_ERRORS)
        target_compile_options(bbsac_project_options INTERFACE /WX)
    endif()
    target_compile_options(bbsac_project_options INTERFACE
        $<$<CONFIG:Debug>:/Od /Zi>
        $<$<CONFIG:Release>:/O2 /DNDEBUG>
        $<$<CONFIG:RelWithDebInfo>:/O2 /Zi /DNDEBUG>
    )
endif()
aux_source_directory(src/BBS_signature PROTOCOL_SOURCES)
#list(REMOVE_ITEM PROTOCOL_SOURCES "src/BBS_signature/prove_CL.cpp")
add_library(bbsac_core
#    src/protocol/protocol_unit.cpp
    ${PROTOCOL_SOURCES}
)

target_include_directories(bbsac_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)



target_compile_features(bbsac_core PUBLIC cxx_std_17)

target_link_libraries(bbsac_core
    PUBLIC
        bbsac_project_options
        mcl
        bicycl
        ${GMP_LIBRARY}
        ${GMPXX_LIBRARY}
        OpenSSL::Crypto
        OpenSSL::SSL
)

add_executable(bbsac_cli tests/BBS_test.cpp)
target_link_libraries(bbsac_cli PRIVATE bbsac_core)
set_target_properties(bbsac_cli PROPERTIES OUTPUT_NAME BBS-DID-Implementation)

if(BBSAC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()


