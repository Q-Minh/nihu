#ifndef REAL_PART_TYPE_HPP_INCLUDED
#define REAL_PART_TYPE_HPP_INCLUDED

namespace NiHu
{

/// \brief metafunction returning the real scalar part type
/// \tparam T the input type
template <class T>
struct real_part_type
{
	typedef T type;
};

/// \brief specialisation of real_part_type for a complex scalar
/// \tparam T the real scalar type
template <class T>
struct real_part_type<std::complex<T> >
{
	typedef T type;
};

} // end of namespace NiHu

#endif // REAL_PART_TYPE_HPP_INCLUDED
