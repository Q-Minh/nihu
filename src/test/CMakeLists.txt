SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -std=c++11" )

enable_testing()

# Domain test
add_executable(domain_test domain_test.cpp)

# Shape function test
add_executable(shape_test shape_test.cpp)

# Element test
add_executable(element_test element_test.cpp)

# Field test
add_executable(field_test field_test.cpp)

# Placeholder expression test
add_executable(placeholderexp_test placeholderexp_test.cpp)

# Lambda metafunction test
add_executable(lambda_test lambda_test.cpp)

add_test(LambdaTest lambda_test "")
