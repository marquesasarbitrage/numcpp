add_library(numcpp INTERFACE)
target_link_libraries(numcpp INTERFACE Eigen3::Eigen)
target_include_directories(numcpp INTERFACE include)
