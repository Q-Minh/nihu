#ifndef BRICK_HPP_INCLUDED
#define BRICK_HPP_INCLUDED

#include "../tmp/bool.hpp"
#include "../tmp/vector.hpp"
#include "../tmp/algorithm.hpp"

/** \brief terminating type of brick types */
struct empty_wall
{
	/** \brief self-returning type */
	typedef empty_wall type;

	/** \brief constructor from everything */
	template <class...Args>
	empty_wall(Args...args) {};
};

/** \brief metafunction to explode a wall into a vector of bricks */
template <class wall>
struct wall_to_bricks : tmp::push_front<
	typename wall_to_bricks<typename wall::base_t>::type,
	typename wall::template wrap<empty_wall>::type
> {};

/** \brief specialisation of ::wall_to_bricks for the empty wall */
template <>
struct wall_to_bricks<empty_wall> : tmp::vector<> {};

/** \brief metafunction to build a wall from a vector of bricks */
template <class bricks>
struct bricks_to_wall
{
	typedef typename tmp::deref<typename tmp::begin<bricks>::type>::type::template wrap<
		typename bricks_to_wall<typename tmp::pop_front<bricks>::type>::type
	>::type type;
};

/** \brief specialisation of ::wall_to_bricks for the empty brick vector case */
template <>
struct bricks_to_wall<tmp::vector<> > : empty_wall {};


/** \brief merge walls A and B */
template <class wallA, class wallB>
struct merge : bricks_to_wall<
	typename tmp::unique<
		typename tmp::concatenate<
			typename wall_to_bricks<wallA>::type,
			typename wall_to_bricks<wallB>::type
		>::type
	>::type
> {};


template <template <class T> class brick, class wall = empty_wall>
struct build :
	public brick<wall>
{
	typedef build type;
	typedef wall base_t;
	template <class newWall> struct wrap : build<brick, newWall> {};

	template <class...Args>
	build(Args const &...args) :
		brick<wall>(args...)
	{
	}
};


#endif // BRICK_HPP_INCLUDED

