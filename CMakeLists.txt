cmake_minimum_required(VERSION 4.1.0)
project(task1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(BUILD_TESTS "Build GTest suite" OFF)

if(ENABLE_ASAN AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(MSVC)
        add_compile_options(/fsanitize=address)
        add_link_options(/fsanitize=address)
    else()

        add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address)
    endif()
endif()

include_directories(include)

add_executable(hello src/main.cpp)

if(BUILD_TESTS)
    find_package(GTest REQUIRED)
    enable_testing()
    set(TEST_EXECUTABLE ${PROJECT_NAME}-tests)
    add_executable(${TEST_EXECUTABLE} tests/main.cpp tests/tdma.cpp)
    target_compile_definitions(
        ${TEST_EXECUTABLE}
        PRIVATE $<$<CONFIG:Debug>:DEBUG>
    )
    target_link_libraries(
        ${TEST_EXECUTABLE}
        PRIVATE GTest::gtest GTest::gtest_main
    )
    include(GoogleTest)
    gtest_discover_tests(${TEST_EXECUTABLE})
endif()
