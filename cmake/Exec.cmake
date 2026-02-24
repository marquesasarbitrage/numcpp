add_executable(numcpp-probability ${CMAKE_CURRENT_SOURCE_DIR}/tests/probability.cpp)
target_link_libraries(numcpp-probability PUBLIC numcpp)

add_executable(numcpp-objects ${CMAKE_CURRENT_SOURCE_DIR}/tests/objects.cpp)
target_link_libraries(numcpp-objects PUBLIC numcpp)

add_executable(numcpp-polysolver ${CMAKE_CURRENT_SOURCE_DIR}/tests/polysolver.cpp)
target_link_libraries(numcpp-polysolver PUBLIC numcpp)

add_executable(numcpp-interpolation ${CMAKE_CURRENT_SOURCE_DIR}/tests/interpolation.cpp)
target_link_libraries(numcpp-interpolation PUBLIC numcpp)

add_executable(numcpp-gaussquad ${CMAKE_CURRENT_SOURCE_DIR}/tests/gaussquad.cpp)
target_link_libraries(numcpp-gaussquad PUBLIC numcpp)