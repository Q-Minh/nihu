#include "../bem/shapeset.hpp"

typedef constant_shape_set<tria_domain> tria_0_shape;
typedef isoparam_shape_set<quad_domain> quad_1_shape;
typedef isoparam_shape_set<tria_domain> tria_1_shape;
typedef parallelogram_shape_set para_shape;

#include <iostream>

int main(void)
{
	std::cout << "1st corner of tria_0: " << std::endl << *(tria_0_shape::corner_begin()) << std::endl << std::endl;
	std::cout << "1st corner of quad_1: " << std::endl << *(quad_1_shape::corner_begin()) << std::endl << std::endl;
	std::cout << "Shape values of quad_1 (0,0): " << std::endl <<
		quad_1_shape::eval_shape(quad_1_shape::xi_t(0.0,0.0)) << std::endl << std::endl;
	std::cout << "Shape derivatives of quad_1 (1.0,-0.75): " << std::endl <<
		quad_1_shape::eval_dshape(quad_1_shape::xi_t(1.0,-0.75)) << std::endl << std::endl;
	std::cout << "3rd corner of para_shape_set: " << std::endl <<
		*(para_shape::corner_begin()+2) << std::endl << std::endl;
	
	return 0;
}