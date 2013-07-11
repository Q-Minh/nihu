SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -std=c++11" )

enable_testing()

# Placeholder expression test
add_executable(placeholderexp_test placeholderexp_test.cpp)

# Lambda metafunction test
add_executable(lambda_test lambda_test.cpp)

add_test(LambdaTest lambda_test "")
