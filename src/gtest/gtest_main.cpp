#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
	int status = 0;
	
	testing::InitGoogleTest(&argc, argv);
	
	testing::GTEST_FLAG(filter) = "bounding_box.*";
	status |= RUN_ALL_TESTS();

	testing::GTEST_FLAG(filter) = "cluster.*";
	status |= RUN_ALL_TESTS();

	testing::GTEST_FLAG(filter) = "elem_center_iterator.*";
	status |= RUN_ALL_TESTS();

	testing::GTEST_FLAG(filter) = "field_center_iterator.*";
	status |= RUN_ALL_TESTS();

	testing::GTEST_FLAG(filter) = "laplace_2d_expansion.*";
	status |= RUN_ALL_TESTS();

	testing::GTEST_FLAG(filter) = "helmholtz_2d_wb_expansion.*";
	status |= RUN_ALL_TESTS();

	testing::GTEST_FLAG(filter) = "chebyshev.*";
	status |= RUN_ALL_TESTS();

#if 0

	testing::GTEST_FLAG(filter) = "laplace_3d_expansion.*";
	status |= RUN_ALL_TESTS();

#endif


	return status;
}
