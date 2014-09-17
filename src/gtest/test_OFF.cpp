#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include <gtest/gtest.h>

TEST(Off, Tag2ElemId)
{
	EXPECT_EQ(line_1_elem::id, nvert2elem_id(2, line_1_tag(), tria_1_tag()) );
	EXPECT_EQ(tria_1_elem::id, nvert2elem_id(3, tria_1_tag()) );
}

