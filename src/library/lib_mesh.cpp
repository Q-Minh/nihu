#include "lib_mesh.hpp"

#include <boost/math/constants/constants.hpp>


namespace NiHu
{

NiHu::mesh<tmp::vector<NiHu::line_1_elem> > 
create_line_1_circle_mesh(double r, int N)
{
	using namespace boost::math::double_constants;
	
	typedef Eigen::Matrix<double, Eigen::Dynamic, 2> dMatrix;
	typedef Eigen::Matrix<unsigned, Eigen::Dynamic, 3> uMatrix;
	
	dMatrix nodes(N,2);
	uMatrix elements(N, 3);
	for (int i = 0; i < N; ++i)
	{
		double phi = i * two_pi / N;
		nodes(i,0) = r * std::cos(phi);
		nodes(i,1) = r * std::sin(phi);
		
		elements(i,0) = NiHu::line_1_elem::id;
		elements(i,1) = i % N;
		elements(i,2) = (i+1) % N;
	}
	
	return NiHu::create_mesh(nodes, elements, NiHu::line_1_tag());
}

} // namespace NiHu
