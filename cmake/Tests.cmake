include(CTest)

set(TESTS
    numcpp-probability      tests/probability.cpp
    numcpp-covmat      tests/covmat.cpp
    numcpp-polysolver      tests/polysolver.cpp
    numcpp-gaussquad      tests/gaussquad.cpp
    numcpp-interpolation      tests/interpolation.cpp
    numcpp-stats      tests/stats.cpp
    numcpp-multstats      tests/multstats.cpp
    numcpp-neldermead      tests/neldermead.cpp
    numcpp-bfgs      tests/bfgs.cpp
    numcpp-bisection      tests/bisection.cpp
    numcpp-nraphson      tests/nraphson.cpp
    numcpp-ols      tests/ols.cpp
    numcpp-stochprocess      tests/stochprocess.cpp
    numcpp-brent      tests/brent.cpp
)

while(TESTS)
    list(POP_FRONT TESTS test_name test_file)

    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} PRIVATE numcpp)

    add_test(NAME ${test_name} COMMAND ${test_name})
endwhile()