#include <gtest/gtest.h>
#include "library/lib_shape.hpp"

typedef Eigen::Matrix<double, 1, 1> Vector1d;
typedef Eigen::Matrix<double, 2, 1> Vector2d;
typedef Eigen::Matrix<double, 3, 1> Vector3d;

Vector1d xi1;
Vector2d xi2;
Vector3d xi3;

typedef decltype(Vector1d::Zero()) z1_t;
typedef decltype(Vector2d::Zero()) z2_t;
typedef decltype(Vector3d::Zero()) z3_t;

template <class ShapeSet>
struct GeneralShapeTester
{
	static void eval(double corners[][ShapeSet::domain_t::dimension])
	{
		double const eps(1e-10);

		// test corner locations
		for (unsigned i = 0; i < ShapeSet::num_nodes; ++i)
			for (unsigned j = 0; j < ShapeSet::domain_t::dimension; ++j)
				EXPECT_NEAR(ShapeSet::corner_at(i)(j), corners[i][j], eps);

		// test shape functions in nodes
		for (unsigned i = 0; i < ShapeSet::num_nodes; ++i)
		{
			auto L = ShapeSet::template eval_shape<0>(ShapeSet::corner_at(i));
			EXPECT_LE((L - decltype(L)::Unit(i)).norm(), eps);
		}
		// test sum of shape functions and derivatives in random locations
		unsigned const N(1000);
		for (unsigned n = 0; n < N; ++n)
		{
			auto xi = ShapeSet::xi_t::Random();
			// sum of L should be one
			auto L = ShapeSet::template eval_shape<0>(xi);
			EXPECT_NEAR(L.sum(), 1., eps);
			// sum of dL should be 0
			auto dL = ShapeSet::template eval_shape<1>(xi);
			EXPECT_NEAR(dL.sum(), 0., eps);
			// sum of ddL should be 0
			auto ddL = ShapeSet::template eval_shape<2>(xi);
			EXPECT_NEAR(ddL.sum(), 0., eps);
		}
	}
};

// testing line_0_shape_set
static_assert(std::is_same<NiHu::line_0_shape_set::domain_t,
	NiHu::line_domain>::value, "Line 0 shape set domain failure");
static_assert(NiHu::line_0_shape_set::num_nodes == 1,
	"Line 0 shape set number of nodes failure");
static_assert(NiHu::line_0_shape_set::polynomial_order == 0,
	"Line 0 shape set polynomial order failure");
static_assert(NiHu::line_0_shape_set::jacobian_order == 0,
	"Line 0 shape set Jacobian order failure");
static_assert(NiHu::line_0_shape_set::id == 1201,
	"Line 0 shape set id failure");
static_assert(std::is_same<NiHu::line_0_shape_set::scalar_t, double>::value,
	"Line 0 shape set scalar failure");
static_assert(std::is_same<NiHu::line_0_shape_set::xi_t, Vector1d>::value,
	"Line 0 shape set xi_t failure");
static_assert(std::is_same<NiHu::line_0_shape_set::shape_t, Vector1d>::value,
	"Line 0 shape set shape value type failure");
static_assert(std::is_same<NiHu::line_0_shape_set::dshape_t, Vector1d>::value,
	"Line 0 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::line_0_shape_set::ddshape_t, Vector1d>::value,
	"Line 0 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::line_0_shape_set::eval_shape<0>(xi1)), Vector1d const &>::value,
	"Line 0 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::line_0_shape_set::eval_shape<1>(xi1)), z1_t>::value,
	"Line 0 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::line_0_shape_set::eval_shape<2>(xi1)), z1_t>::value,
	"Line 0 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::line_0_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof1> >::value,
		"Line 0 shape set poisition DOF vector failure");

TEST(ShapeSet, Line0){
	double corners[][1] = { {0.} };
	GeneralShapeTester<NiHu::line_0_shape_set>::eval(corners);
}



// testing line_1_shape_set
static_assert(std::is_same<NiHu::line_1_shape_set::domain_t,
	NiHu::line_domain>::value, "Line 1 shape set domain failure");
static_assert(NiHu::line_1_shape_set::num_nodes == 2,
	"Line 1 shape set number of nodes failure");
static_assert(NiHu::line_1_shape_set::polynomial_order == 1,
	"Line 1 shape set polynomial order failure");
static_assert(NiHu::line_1_shape_set::jacobian_order == 0,
	"Line 1 shape set Jacobian order failure");
static_assert(NiHu::line_1_shape_set::id == 1202,
	"Line 1 shape set id failure");
static_assert(std::is_same<NiHu::line_1_shape_set::scalar_t, double>::value,
	"Line 1 shape set scalar failure");
static_assert(std::is_same<NiHu::line_1_shape_set::xi_t, Vector1d>::value,
	"Line 1 shape set xi_t failure");
static_assert(std::is_same<NiHu::line_1_shape_set::shape_t, Vector2d>::value,
	"Line 1 shape set shape value type failure");
static_assert(std::is_same<NiHu::line_1_shape_set::dshape_t, Vector2d>::value,
	"Line 1 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::line_1_shape_set::ddshape_t,
	Vector2d>::value,
	"Line 1 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::line_1_shape_set::eval_shape<0>(xi1)), Vector2d>::value,
	"Line 1 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::line_1_shape_set::eval_shape<1>(xi1)), Vector2d const &>::value,
	"Line 1 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::line_1_shape_set::eval_shape<2>(xi1)), z2_t>::value,
	"Line 1 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::line_1_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof0, NiHu::dof0> >::value,
	"Line 1 shape set poisition DOF vector failure");

TEST(ShapeSet, Line1){
	double corners[][1] = { {-1.}, {1.} };
	GeneralShapeTester<NiHu::line_1_shape_set>::eval(corners);
}


// testing line_1_gauss_shape_set
TEST(ShapeSet, Line1Gauss) {
	// test nodal corners
	double corners[][1] = {
		{-1. / std::sqrt(3.0)}, {+1. / std::sqrt(3.0)}
	};
	GeneralShapeTester<NiHu::line_1_gauss_shape_set>::eval(corners);
}


// testing line_2_shape_set
static_assert(std::is_same<NiHu::line_2_shape_set::domain_t,
	NiHu::line_domain>::value, "Line 2 shape set domain failure");
static_assert(NiHu::line_2_shape_set::num_nodes == 3,
	"Line 2 shape set number of nodes failure");
static_assert(NiHu::line_2_shape_set::polynomial_order == 2,
	"Line 2 shape set polynomial order failure");
static_assert(NiHu::line_2_shape_set::jacobian_order == 1,
	"Line 2 shape set Jacobian order failure");
static_assert(NiHu::line_2_shape_set::id == 1203,
	"Line 2 shape set id failure");
static_assert(std::is_same<NiHu::line_2_shape_set::scalar_t, double>::value,
	"Line 2 shape set scalar failure");
static_assert(std::is_same<NiHu::line_2_shape_set::xi_t, Vector1d>::value,
	"Line 2 shape set xi_t failure");
static_assert(std::is_same<NiHu::line_2_shape_set::shape_t, Vector3d>::value,
	"Line 2 shape set shape value type failure");
static_assert(std::is_same<NiHu::line_2_shape_set::dshape_t, Vector3d>::value,
	"Line 2 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::line_2_shape_set::ddshape_t, Vector3d>::value,
	"Line 2 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::line_2_shape_set::eval_shape<0>(xi1)), Vector3d>::value,
	"Line 2 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::line_2_shape_set::eval_shape<1>(xi1)), Vector3d>::value,
	"Line 2 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::line_2_shape_set::eval_shape<2>(xi1)), Vector3d const &>::value,
	"Line 2 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::line_2_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof0, NiHu::dof1, NiHu::dof0> >::value,
	"Line 2 shape set poisition DOF vector failure");

TEST(ShapeSet, Line2){
	double corners[][1] = { {-1.}, {0.}, {1.} };
	GeneralShapeTester<NiHu::line_2_shape_set>::eval(corners);
}


// testing tria_0_shape_set
static_assert(std::is_same<NiHu::tria_0_shape_set::domain_t,
	NiHu::tria_domain>::value, "Tria 0 shape set domain failure");
static_assert(NiHu::tria_0_shape_set::num_nodes == 1,
	"Tria 0 shape set number of nodes failure");
static_assert(NiHu::tria_0_shape_set::polynomial_order == 0,
	"Tria 0 shape set polynomial order failure");
static_assert(NiHu::tria_0_shape_set::jacobian_order == 0,
	"Tria 0 shape set Jacobian order failure");
static_assert(NiHu::tria_0_shape_set::id == 2301,
	"Tria 0 shape set id failure");
static_assert(std::is_same<NiHu::tria_0_shape_set::scalar_t, double>::value,
	"Tria 0 shape set scalar failure");
static_assert(std::is_same<NiHu::tria_0_shape_set::xi_t, Vector2d>::value,
	"Tria 0 shape set xi_t failure");
static_assert(std::is_same<NiHu::tria_0_shape_set::shape_t, Vector1d>::value,
	"Tria 0 shape set shape value type failure");
static_assert(std::is_same<NiHu::tria_0_shape_set::dshape_t,
	Eigen::Matrix<double, 1, 2> >::value,
	"Tria 0 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::tria_0_shape_set::ddshape_t,
	Eigen::Matrix<double, 1, 3> >::value,
	"Tria 0 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::tria_0_shape_set::eval_shape<0>(xi2)), Vector1d const &>::value,
	"Tria 0 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::tria_0_shape_set::eval_shape<1>(xi2)), decltype(Eigen::Matrix<double, 1, 2>::Zero())>::value,
	"Tria 0 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::tria_0_shape_set::eval_shape<2>(xi2)), decltype(Eigen::Matrix<double, 1, 3>::Zero())>::value,
	"Tria 0 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::tria_0_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof2> >::value,
	"Tria 0 shape set poisition DOF vector failure");

TEST(ShapeSet, Tria0){
	double corners[][2] = { {1./3., 1./3.} };
	GeneralShapeTester<NiHu::tria_0_shape_set>::eval(corners);
}


// testing tria_1_shape_set
static_assert(std::is_same<NiHu::tria_1_shape_set::domain_t, NiHu::tria_domain>::value, "Tria 1 shape set domain failure");
static_assert(NiHu::tria_1_shape_set::num_nodes == 3,
	"Tria 1 shape set number of nodes failure");
static_assert(NiHu::tria_1_shape_set::polynomial_order == 1,
	"Tria 1 shape set polynomial order failure");
static_assert(NiHu::tria_1_shape_set::jacobian_order == 0,
	"Tria 1 shape set Jacobian order failure");
static_assert(NiHu::tria_1_shape_set::id == 2303,
	"Tria 1 shape set id failure");
static_assert(std::is_same<NiHu::tria_1_shape_set::scalar_t, double>::value,
	"Tria 1 shape set scalar failure");
static_assert(std::is_same<NiHu::tria_1_shape_set::xi_t, Vector2d>::value,
	"Tria 1 shape set xi_t failure");
static_assert(std::is_same<NiHu::tria_1_shape_set::shape_t, Vector3d>::value,
	"Tria 1 shape set shape value type failure");
static_assert(std::is_same<NiHu::tria_1_shape_set::dshape_t,
	Eigen::Matrix<double, 3, 2> >::value,
	"Tria 1 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::tria_1_shape_set::ddshape_t,
	Eigen::Matrix<double, 3, 3> >::value,
	"Tria 1 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::tria_1_shape_set::eval_shape<0>(xi2)), Vector3d>::value,
	"Tria 1 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::tria_1_shape_set::eval_shape<1>(xi2)), Eigen::Matrix<double, 3, 2> const &>::value,
	"Tria 1 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::tria_1_shape_set::eval_shape<2>(xi2)), decltype(Eigen::Matrix<double, 3, 3>::Zero())>::value,
	"Tria 1 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::tria_1_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof0, NiHu::dof0, NiHu::dof0> >::value,
	"Tria 1 shape set poisition DOF vector failure");

TEST(ShapeSet, Tria1){
	// test nodal corners
	double corners[][2] = { {0., 0.}, {1., 0.}, {0., 1.} };
	GeneralShapeTester<NiHu::tria_1_shape_set>::eval(corners);
}


// testing tria_2_shape_set
static_assert(std::is_same<NiHu::tria_2_shape_set::domain_t,
	NiHu::tria_domain>::value, "Tria 2 shape set domain failure");
static_assert(NiHu::tria_2_shape_set::num_nodes == 6,
	"Tria 2 shape set number of nodes failure");
static_assert(NiHu::tria_2_shape_set::polynomial_order == 2,
	"Tria 2 shape set polynomial order failure");
static_assert(NiHu::tria_2_shape_set::jacobian_order == 2,
	"Tria 2 shape set Jacobian order failure");
static_assert(NiHu::tria_2_shape_set::id == 2306, "Tria 2 shape set id failure");
static_assert(std::is_same<NiHu::tria_2_shape_set::scalar_t, double>::value,
	"Tria 2 shape set scalar failure");
static_assert(std::is_same<NiHu::tria_2_shape_set::xi_t, Vector2d>::value, "Tria 2 shape set xi_t failure");
static_assert(std::is_same<NiHu::tria_2_shape_set::shape_t,
	Eigen::Matrix<double, 6, 1> >::value,
	"Tria 2 shape set shape value type failure");
static_assert(std::is_same<NiHu::tria_2_shape_set::dshape_t,
	Eigen::Matrix<double, 6, 2> >::value,
	"Tria 2 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::tria_2_shape_set::ddshape_t,
	Eigen::Matrix<double, 6, 3> >::value,
	"Tria 2 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::tria_2_shape_set::eval_shape<0>(xi2)), Eigen::Matrix<double, 6, 1> >::value,
	"Tria 2 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::tria_2_shape_set::eval_shape<1>(xi2)), Eigen::Matrix<double, 6, 2> >::value,
	"Tria 2 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::tria_2_shape_set::eval_shape<2>(xi2)), Eigen::Matrix<double, 6, 3> >::value,
	"Tria 2 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::tria_2_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1> >::value,
	"Tria 2 shape set poisition DOF vector failure");

TEST(ShapeSet, Tria2){
	// test nodal corners
	double corners[][2] = { {0., 0.}, {.5, 0.}, {1., 0.}, {.5, .5}, {0., 1.}, {0., .5} };
	GeneralShapeTester<NiHu::tria_2_shape_set>::eval(corners);
}


// testing quad_0_shape_set
static_assert(std::is_same<NiHu::quad_0_shape_set::domain_t,
	NiHu::quad_domain>::value, "Quad 0 shape set domain failure");
static_assert(NiHu::quad_0_shape_set::num_nodes == 1,
	"Quad 0 shape set number of nodes failure");
static_assert(NiHu::quad_0_shape_set::polynomial_order == 0,
	"Quad 0 shape set polynomial order failure");
static_assert(NiHu::quad_0_shape_set::jacobian_order == 0,
	"Quad 0 shape set Jacobian order failure");
static_assert(NiHu::quad_0_shape_set::id == 2401,
	"Quad 0 shape set id failure");
static_assert(std::is_same<NiHu::quad_0_shape_set::scalar_t, double>::value,
	"Quad 0 shape set scalar failure");
static_assert(std::is_same<NiHu::quad_0_shape_set::xi_t, Vector2d>::value,
	"Quad 0 shape set xi_t failure");
static_assert(std::is_same<NiHu::quad_0_shape_set::shape_t, Vector1d>::value,
	"Quad 0 shape set shape value type failure");
static_assert(std::is_same<NiHu::quad_0_shape_set::dshape_t,
	Eigen::Matrix<double, 1, 2> >::value,
	"Quad 0 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::quad_0_shape_set::ddshape_t,
	Eigen::Matrix<double, 1, 3> >::value,
	"Quad 0 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::quad_0_shape_set::eval_shape<0>(xi2)), Vector1d const &>::value,
	"Quad 0 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::quad_0_shape_set::eval_shape<1>(xi2)), decltype(Eigen::Matrix<double, 1, 2>::Zero())>::value,
	"Quad 0 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::quad_0_shape_set::eval_shape<2>(xi2)), decltype(Eigen::Matrix<double, 1, 3>::Zero())>::value,
	"Quad 0 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::quad_0_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof2> >::value,
	"Quad 0 shape set poisition DOF vector failure");

TEST(ShapeSet, Quad0){
	// test nodal corners
	double corners[][2] = { {0., 0.} };
	GeneralShapeTester<NiHu::quad_0_shape_set>::eval(corners);
}


// testing quad_1_shape_set
static_assert(std::is_same<NiHu::quad_1_shape_set::domain_t,
	NiHu::quad_domain>::value, "Quad 1 shape set domain failure");
static_assert(NiHu::quad_1_shape_set::num_nodes == 4,
	"Quad 1 shape set number of nodes failure");
static_assert(NiHu::quad_1_shape_set::polynomial_order == 1,
	"Quad 1 shape set polynomial order failure");
static_assert(NiHu::quad_1_shape_set::jacobian_order == 1,
	"Quad 1 shape set Jacobian order failure");
static_assert(NiHu::quad_1_shape_set::id == 2404,
	"Quad 1 shape set id failure");
static_assert(std::is_same<NiHu::quad_1_shape_set::scalar_t, double>::value,
	"Quad 1 shape set scalar failure");
static_assert(std::is_same<NiHu::quad_1_shape_set::xi_t, Vector2d>::value,
	"Quad 1 shape set xi_t failure");
static_assert(std::is_same<NiHu::quad_1_shape_set::shape_t,
	Eigen::Matrix<double, 4, 1> >::value,
	"Quad 1 shape set shape value type failure");
static_assert(std::is_same<NiHu::quad_1_shape_set::dshape_t,
	Eigen::Matrix<double, 4, 2> >::value,
	"Quad 1 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::quad_1_shape_set::ddshape_t, Eigen::Matrix<double, 4, 3> >::value,
	"Quad 1 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::quad_1_shape_set::eval_shape<0>(xi2)), Eigen::Matrix<double, 4, 1> >::value,
	"Quad 1 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::quad_1_shape_set::eval_shape<1>(xi2)), Eigen::Matrix<double, 4, 2> >::value,
	"Quad 1 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::quad_1_shape_set::eval_shape<2>(xi2)), Eigen::Matrix<double, 4, 3> const &>::value,
	"Quad 1 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::quad_1_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof0, NiHu::dof0, NiHu::dof0, NiHu::dof0> >::value,
	"Quad 1 shape set poisition DOF vector failure");

TEST(ShapeSet, Quad1){
	// test nodal corners
	double corners[][2] = { {-1., -1.}, {1., -1.}, {1., 1.}, {-1., 1.} };
	GeneralShapeTester<NiHu::quad_1_shape_set>::eval(corners);
}


// testing quad_2_shape_set
static_assert(std::is_same<NiHu::quad_2_shape_set::domain_t,
	NiHu::quad_domain>::value, "Quad 2 shape set domain failure");
static_assert(NiHu::quad_2_shape_set::num_nodes == 9,
	"Quad 2 shape set number of nodes failure");
static_assert(NiHu::quad_2_shape_set::polynomial_order == 2,
	"Quad 2 shape set polynomial order failure");
static_assert(NiHu::quad_2_shape_set::jacobian_order == 3,
	"Quad 2 shape set Jacobian order failure");
static_assert(NiHu::quad_2_shape_set::id == 2409,
	"Quad 2 shape set id failure");
static_assert(std::is_same<NiHu::quad_2_shape_set::scalar_t, double>::value,
	"Quad 2 shape set scalar failure");
static_assert(std::is_same<NiHu::quad_2_shape_set::xi_t, Vector2d>::value,
	"Quad 2 shape set xi_t failure");
static_assert(std::is_same<NiHu::quad_2_shape_set::shape_t,
	Eigen::Matrix<double, 9, 1> >::value,
	"Quad 2 shape set shape value type failure");
static_assert(std::is_same<NiHu::quad_2_shape_set::dshape_t,
	Eigen::Matrix<double, 9, 2> >::value,
	"Quad 2 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::quad_2_shape_set::ddshape_t,
	Eigen::Matrix<double, 9, 3> >::value,
	"Quad 2 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::quad_2_shape_set::eval_shape<0>(xi2)), Eigen::Matrix<double, 9, 1>	>::value,
	"Quad 2 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::quad_2_shape_set::eval_shape<1>(xi2)), Eigen::Matrix<double, 9, 2> >::value,
	"Quad 2 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::quad_2_shape_set::eval_shape<2>(xi2)), Eigen::Matrix<double, 9, 3> >::value,
	"Quad 2 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::quad_2_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1, NiHu::dof2> >::value,
	"Quad 2 shape set poisition DOF vector failure");

TEST(ShapeSet, Quad2){
	// test nodal corners
	double corners[][2] = { {-1., -1.}, {0., -1.}, {1., -1.}, {1., 0.}, {1., 1.}, {0., 1.}, {-1., 1.}, {-1., 0.}, {0., 0.} };
	GeneralShapeTester<NiHu::quad_2_shape_set>::eval(corners);
}


// testing quad_28_shape_set
static_assert(std::is_same<NiHu::quad_28_shape_set::domain_t,
	NiHu::quad_domain>::value, "Quad 28 shape set domain failure");
static_assert(NiHu::quad_28_shape_set::num_nodes == 8,
	"Quad 28 shape set number of nodes failure");
static_assert(NiHu::quad_28_shape_set::polynomial_order == 2,
	"Quad 28 shape set polynomial order failure");
static_assert(NiHu::quad_28_shape_set::jacobian_order == 3,
	"Quad 28 shape set Jacobian order failure");
static_assert(NiHu::quad_28_shape_set::id == 2408,
	"Quad 28 shape set id failure");
static_assert(std::is_same<NiHu::quad_28_shape_set::scalar_t, double>::value,
	"Quad 28 shape set scalar failure");
static_assert(std::is_same<NiHu::quad_28_shape_set::xi_t, Vector2d>::value,
	"Quad 28 shape set xi_t failure");
static_assert(std::is_same<NiHu::quad_28_shape_set::shape_t,
	Eigen::Matrix<double, 8, 1> >::value,
	"Quad 28 shape set shape value type failure");
static_assert(std::is_same<NiHu::quad_28_shape_set::dshape_t,
	Eigen::Matrix<double, 8, 2> >::value,
	"Quad 28 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::quad_28_shape_set::ddshape_t,
	Eigen::Matrix<double, 8, 3> >::value,
	"Quad 28 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::quad_28_shape_set::eval_shape<0>(xi2)), Eigen::Matrix<double, 8, 1> >::value,
	"Quad 28 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::quad_28_shape_set::eval_shape<1>(xi2)), Eigen::Matrix<double, 8, 2> >::value,
	"Quad 28 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::quad_28_shape_set::eval_shape<2>(xi2)), Eigen::Matrix<double, 8, 3> >::value,
	"Quad 28 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::quad_28_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1, NiHu::dof0, NiHu::dof1> >::value,
	"Quad 28 shape set poisition DOF vector failure");

TEST(ShapeSet, Quad28){
	// test nodal corners
	double corners[][2] =
	{
		{-1., -1.}, {0., -1.}, {1., -1.}, {1., 0.},
		{1., 1.}, {0., 1.}, {-1., 1.}, {-1., 0.}
	};
	GeneralShapeTester<NiHu::quad_28_shape_set>::eval(corners);
}


// testing brick_0_shape_set
static_assert(std::is_same<NiHu::brick_0_shape_set::domain_t,
	NiHu::brick_domain>::value, "Brick 0 shape set domain failure");
static_assert(NiHu::brick_0_shape_set::num_nodes == 1,
	"Brick 0 shape set number of nodes failure");
static_assert(NiHu::brick_0_shape_set::polynomial_order == 0,
	"Brick 0 shape set polynomial order failure");
static_assert(NiHu::brick_0_shape_set::jacobian_order == 0,
	"Brick 0 shape set Jacobian order failure");
static_assert(NiHu::brick_0_shape_set::id == 3801,
	"Brick 0 shape set id failure");
static_assert(std::is_same<NiHu::brick_0_shape_set::scalar_t, double>::value,
	"Brick 0 shape set scalar failure");
static_assert(std::is_same<NiHu::brick_0_shape_set::xi_t, Vector3d>::value,
	"Brick 0 shape set xi_t failure");
static_assert(std::is_same<NiHu::brick_0_shape_set::shape_t, Vector1d>::value,
	"Brick 0 shape set shape value type failure");
static_assert(std::is_same<NiHu::brick_0_shape_set::dshape_t,
	Eigen::Matrix<double, 1, 3> >::value,
	"Brick 0 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::brick_0_shape_set::ddshape_t,
	Eigen::Matrix<double, 1, 6> >::value,
	"Brick 0 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::brick_0_shape_set::eval_shape<0>(xi3)), Vector1d const &>::value,
	"Brick 0 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::brick_0_shape_set::eval_shape<1>(xi3)), decltype(Eigen::Matrix<double, 1, 3>::Zero())>::value,
	"Brick 0 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::brick_0_shape_set::eval_shape<2>(xi3)), decltype(Eigen::Matrix<double, 1, 6>::Zero())>::value,
	"Brick 0 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::brick_0_shape_set::position_dof_vector,
	tmp::vector<NiHu::dof3> >::value,
	"Brick 0 shape set poisition DOF vector failure");

TEST(ShapeSet, Brick0){
	// test nodal corners
	double corners[][3] = { {0., 0., 0.} };
	GeneralShapeTester<NiHu::brick_0_shape_set>::eval(corners);
}


// testing brick_1_shape_set
static_assert(std::is_same<NiHu::brick_1_shape_set::domain_t,
	NiHu::brick_domain>::value, "Brick 1 shape set domain failure");
static_assert(NiHu::brick_1_shape_set::num_nodes == 8,
	"Brick 1 shape set number of nodes failure");
static_assert(NiHu::brick_1_shape_set::polynomial_order == 1,
	"Brick 1 shape set polynomial order failure");
static_assert(NiHu::brick_1_shape_set::jacobian_order == 1,
	"Brick 1 shape set Jacobian order failure");
static_assert(NiHu::brick_1_shape_set::id == 3808,
	"Brick 1 shape set id failure");
static_assert(std::is_same<NiHu::brick_1_shape_set::scalar_t, double>::value,
	"Brick 1 shape set scalar failure");
static_assert(std::is_same<NiHu::brick_1_shape_set::xi_t, Vector3d>::value,
	"Brick 1 shape set xi_t failure");
static_assert(std::is_same<NiHu::brick_1_shape_set::shape_t,
	Eigen::Matrix<double, 8, 1> >::value,
	"Brick 1 shape set shape value type failure");
static_assert(std::is_same<NiHu::brick_1_shape_set::dshape_t,
	Eigen::Matrix<double, 8, 3> >::value,
	"Brick 1 shape set shape derivative value type failure");
static_assert(std::is_same<NiHu::brick_1_shape_set::ddshape_t, Eigen::Matrix<double, 8, 6> >::value,
	"Brick 1 shape set shape 2nd derivative value type failure");
static_assert(std::is_same<decltype(NiHu::brick_1_shape_set::eval_shape<0>(xi3)), Eigen::Matrix<double, 8, 1> >::value,
	"Brick 1 shape set shape return type failure");
static_assert(std::is_same<decltype(NiHu::brick_1_shape_set::eval_shape<1>(xi3)), Eigen::Matrix<double, 8, 3>>::value,
	"Brick 1 shape set shape derivative return type failure");
static_assert(std::is_same<decltype(NiHu::brick_1_shape_set::eval_shape<2>(xi3)), Eigen::Matrix<double, 8, 6> >::value,
	"Brick 1 shape set shape second derivative return type failure");
static_assert(std::is_same<NiHu::brick_1_shape_set::position_dof_vector, tmp::vector<NiHu::dof0, NiHu::dof0, NiHu::dof0, NiHu::dof0, NiHu::dof0, NiHu::dof0, NiHu::dof0, NiHu::dof0> >::value,
	"Brick 1 shape set poisition DOF vector failure");

TEST(ShapeSet, Brick1){
	// test nodal corners
	double corners[][3] = {
		{-1., -1., -1.}, {1., -1., -1.}, {1., 1., -1.}, {-1., 1., -1.},
		{-1., -1., 1.}, {1., -1., 1.}, {1., 1., 1.}, {-1., 1., 1.}
	 };
	GeneralShapeTester<NiHu::brick_1_shape_set>::eval(corners);
}

