enable_testing()

#add_executable(testList
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/probability.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/objects.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/polysolver.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/interpolation.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/gaussquad.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/stats.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/multstats.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/tests/functions.cpp
#)

foreach(test_file functions probability objects polysolver gaussquad interpolation stats multstats)
    add_executable(numcpptests-${test_file} tests/${test_file}.cpp)
    target_link_libraries(numcpptests-${test_file} numcpp)
    add_test(NAME ${test_file} COMMAND numcpptests-${test_file})
endforeach()

#target_link_libraries(testList numcpp)

#add_test(NAME tests COMMAND testList)