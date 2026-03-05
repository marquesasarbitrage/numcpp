enable_testing()

add_executable(testList
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/probability.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/objects.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/polysolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/interpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/gaussquad.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/stats.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/multstats.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/functions.cpp
)

target_link_libraries(testList numcpp)

add_test(NAME tests COMMAND testList)