cmake_minimum_required(VERSION 3.20)

project(bigint-cpp
    VERSION 1.0.0
    DESCRIPTION "Arbitrary precision integer arithmetic library for C++"
    LANGUAGES CXX
)

# ============================================================================
# Project Options
# ============================================================================
# Only build tests/benchmarks by default when this is the top-level project
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(BIGINT_IS_TOP_LEVEL ON)
else()
    set(BIGINT_IS_TOP_LEVEL OFF)
endif()

option(BIGINT_BUILD_TESTS "Build bigint-cpp tests" ${BIGINT_IS_TOP_LEVEL})
option(BIGINT_BUILD_BENCHMARKS "Build bigint-cpp benchmarks" ${BIGINT_IS_TOP_LEVEL})

# ============================================================================
# C++ Standard
# ============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ============================================================================
# Compiler Warnings
# ============================================================================
add_library(bigint_warnings INTERFACE)
add_library(bigint::warnings ALIAS bigint_warnings)

if(MSVC)
    target_compile_options(bigint_warnings INTERFACE
        /W4
        /permissive-
    )
else()
    target_compile_options(bigint_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wconversion
        -Wsign-conversion
    )
endif()

# ============================================================================
# Main Library
# ============================================================================
add_library(bigint STATIC
    src/BigInt.cpp
    src/MontgomeryContext.cpp
)
add_library(bigint::bigint ALIAS bigint)

target_include_directories(bigint
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(bigint
    PRIVATE
        bigint::warnings
)

# ============================================================================
# Tests
# ============================================================================
if(BIGINT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ============================================================================
# Benchmarks
# ============================================================================
if(BIGINT_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# ============================================================================
# Install Configuration
# ============================================================================
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS bigint bigint_warnings
    EXPORT bigint-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT bigint-targets
    FILE bigint-targets.cmake
    NAMESPACE bigint::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bigint
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/bigint-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/bigint-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bigint
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/bigint-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/bigint-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/bigint-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bigint
)
