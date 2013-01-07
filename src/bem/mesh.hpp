#include <array>

template <unsigned nDim>
class Coord : public std::array<double, nDim> { };


class tria_tag;
class quad_tag;

template <class> class nNodes;
template<> struct nNodes<quad_tag> { static unsigned const value = 3; };
template<> struct nNodes<tria_tag> { static unsigned const value = 4; };

template <class ElemType>
class Elem : public std::array<unsigned, nNodes<ElemType>::value> { };

typedef Elem<quad_tag> QuadElem;
typedef Elem<tria_tag> TriaElem;

#include <vector>

template <class Elem>
class ElemContainer
{
protected:
	std::vector<Elem > elements;
};

template <unsigned nDim>
class Mesh : public ElemContainer<QuadElem>, ElemContainer<TriaElem>
{
private:
	std::vector<TriaElem > trias;

public:
	template <class Elem>
	void add_elem(Elem const &e) { ElemContainer<Elem>::elements.push_back(e); }

	template <class Elem>
	typename std::vector<Elem>::iterator begin(void) { return ElemContainer<Elem>::elements.begin(); }

	template <class Elem>
	typename std::vector<Elem>::iterator end(void) { return ElemContainer<Elem>::elements.end(); }
};

