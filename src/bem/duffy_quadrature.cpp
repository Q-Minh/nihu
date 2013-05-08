#include "duffy_quadrature.hpp"

unsigned const duffy_traits<tria_1_shape_set>::duffy_corner_indices[3][2+1] = {
	{1, /*|*/ 1, 2},
	{1, /*|*/ 2, 0},
	{1, /*|*/ 0, 1}
};

unsigned const duffy_traits<tria_1_shape_set>::duffy_face_indices[4+1] = {
	3, /*|*/ 0, 1, 2, 0
};


unsigned const duffy_traits<quad_1_shape_set>::duffy_corner_indices[4][3+1] = {
	{2, /*|*/ 1, 2, 3},
	{2, /*|*/ 2, 3, 0},
	{2, /*|*/ 3, 0, 1},
	{2, /*|*/ 0, 1, 2}
};

unsigned const duffy_traits<quad_1_shape_set>::duffy_face_indices[5+1] = {
	4, /*|*/ 0, 1, 2, 3, 0
};


unsigned const duffy_traits<tria_2_shape_set>::duffy_corner_indices[6][3+1] = {
	{1, /*|*/ 2, 4},
	{2, /*|*/ 2, 4, 0},
	{1, /*|*/ 4, 0},
	{2, /*|*/ 4, 0, 2},
	{1, /*|*/ 0, 2},
	{2, /*|*/ 0, 2, 4}
};

unsigned const duffy_traits<tria_2_shape_set>::duffy_face_indices[4+1] = {
	3, /*|*/ 0, 2, 4, 0
};


