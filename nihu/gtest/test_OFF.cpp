#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include <gtest/gtest.h>

TEST(Off, Tag2ElemId)
{
	EXPECT_EQ(NiHu::line_1_elem::id, nvert2elem_id(2, NiHu::line_1_tag(), NiHu::tria_1_tag()) );
	EXPECT_EQ(NiHu::tria_1_elem::id, nvert2elem_id(3, NiHu::tria_1_tag()) );
}

