cmake_minimum_required(VERSION 3.20)

project(BP_LCGCD VERSION 0.1.0 LANGUAGES C)

add_executable(bp_lcgcd src/apps/main.c)
add_library(bplcgcd_gf2 STATIC src/core/gf2.c src/core/systematic.c)
target_include_directories(bplcgcd_gf2 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
add_library(bplcgcd_decoder STATIC src/decoders/decoder.c)
target_link_libraries(bplcgcd_decoder PUBLIC bplcgcd_gf2 bplcgcd_search)
target_link_libraries(bp_lcgcd PRIVATE bplcgcd_decoder)
add_library(bplcgcd_search STATIC src/search/tep.c src/search/local_constraints.c src/search/trellis.c src/search/viterbi.c src/search/slva.c)
target_link_libraries(bplcgcd_search PUBLIC bplcgcd_gf2)
add_library(bplcgcd_channel STATIC src/channel/channel.c)
target_link_libraries(bplcgcd_channel PUBLIC bplcgcd_gf2)
if(NOT MSVC)
    target_link_libraries(bplcgcd_channel PUBLIC m)
endif()
add_library(bplcgcd_simulation STATIC src/simulation/simulation.c)
target_link_libraries(bplcgcd_simulation PUBLIC bplcgcd_channel bplcgcd_decoder)
target_link_libraries(bp_lcgcd PRIVATE bplcgcd_simulation)

foreach(target IN ITEMS bp_lcgcd bplcgcd_gf2 bplcgcd_decoder bplcgcd_search bplcgcd_channel bplcgcd_simulation)
    set_target_properties(${target} PROPERTIES
        C_STANDARD 11
        C_STANDARD_REQUIRED YES
        C_EXTENSIONS NO
    )
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endforeach()

include(CTest)
if(BUILD_TESTING)
    add_test(NAME cli_smoke
        COMMAND "${CMAKE_COMMAND}"
            "-DCLI_EXECUTABLE=$<TARGET_FILE:bp_lcgcd>"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_smoke.cmake"
    )
    set_tests_properties(cli_smoke PROPERTIES TIMEOUT 10)
    add_test(NAME simulation_cli COMMAND "${CMAKE_COMMAND}"
        "-DCLI_EXECUTABLE=$<TARGET_FILE:bp_lcgcd>"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/simulation_cli.cmake")
    set_tests_properties(simulation_cli PROPERTIES TIMEOUT 15)

    foreach(test IN ITEMS gf2 systematic decoder tep gcd channel simulation local_constraints trellis viterbi slva)
        add_executable(${test}_test tests/${test}_test.c)
        if(test STREQUAL "simulation")
            target_link_libraries(${test}_test PRIVATE bplcgcd_simulation)
        elseif(test STREQUAL "channel")
            target_link_libraries(${test}_test PRIVATE bplcgcd_channel bplcgcd_decoder)
        elseif(test STREQUAL "decoder" OR test STREQUAL "gcd" OR test STREQUAL "local_constraints" OR test STREQUAL "trellis")
            target_link_libraries(${test}_test PRIVATE bplcgcd_decoder)
        elseif(test STREQUAL "tep" OR test STREQUAL "viterbi" OR test STREQUAL "slva")
            target_link_libraries(${test}_test PRIVATE bplcgcd_search)
        else()
            target_link_libraries(${test}_test PRIVATE bplcgcd_gf2)
        endif()
        set_target_properties(${test}_test PROPERTIES
            C_STANDARD 11 C_STANDARD_REQUIRED YES C_EXTENSIONS NO)
        if(MSVC)
            target_compile_options(${test}_test PRIVATE /W4)
        else()
            target_compile_options(${test}_test PRIVATE -Wall -Wextra -Wpedantic)
        endif()
        add_test(NAME ${test} COMMAND ${test}_test)
        set_tests_properties(${test} PROPERTIES TIMEOUT 10)
    endforeach()
endif()
