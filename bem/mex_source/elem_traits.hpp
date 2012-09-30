template <class elemType>
struct elem_traits;

class TriaElem;
class QuadElem;

template<>
struct elem_traits<TriaElem>
{
    enum
    {
        nNodes = 3,
        isLinear = true
    };
};


template<>
struct elem_traits<QuadElem>
{
    enum
    {
        nNodes = 4,
        isLinear = false
    };
};
