/**
 * \file brick.hpp
 * \brief definition of brick data classes
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef BRICK_HPP_INCLUDED
#define BRICK_HPP_INCLUDED

#include "../tmp/bool.hpp"
#include "../tmp/vector.hpp"
#include "../tmp/algorithm.hpp"

namespace internal
{
	/** \brief terminating type of wall types */
	struct empty_wall
	{
		/** \brief self-returning type */
		typedef empty_wall type;

		/** \brief constructor from everything */
		template <class...Args>
		empty_wall(Args...args) {}
	};

	/** \brief metafunction to explode a wall into a vector of bricks */
	template <class wall>
	struct wall_to_bricks : tmp::push_back<
		typename wall_to_bricks<typename wall::base_t>::type,
		typename wall::template wrap<empty_wall>::type
	> {};

	/** \brief specialisation of ::wall_to_bricks for the empty wall */
	template <>
	struct wall_to_bricks<empty_wall> : tmp::vector<> {};

	/** \brief helper metafunction to put a brick on the top of a wall */
	template <class brick, class wall>
	struct put_on : brick::template wrap<wall> {};

	/** \brief metafunction to build a wall from a vector of bricks */
	template <class bricks>
	struct bricks_to_wall : tmp::accumulate<
		bricks,
		empty_wall,
		put_on<tmp::_2, tmp::_1>
	> {};
}

/** \brief merge walls
 * \tparam wallA the first wall to merge
 * \tparam wallB the second wall to merge
 * \return a wall consisting of unique bricks in the same order
 */
template <class wallA, class wallB = internal::empty_wall>
struct merge : internal::bricks_to_wall<
	typename tmp::unique<
		typename tmp::concatenate<
			typename internal::wall_to_bricks<wallA>::type,
			typename internal::wall_to_bricks<wallB>::type
		>::type
	>::type
> {};

/** \brief convert a brick template and a wall into a wall
 * \tparam brick a brick template
 * \tparam wall a wall
 * \return a wall with the brick on the top
 */
template <template <class T> class brick, class wall = internal::empty_wall>
struct build :
	public brick<wall>
{
	/** \brief self-returning metafunction */
	typedef build type;
	/** \brief the supporting wall */
	typedef wall base_t;
	/** \brief metafunction to replace top brick to a new wall
	 * \tparam newWall the new wall where the top brick is to be placed
	 * \return the new wall with the brick on top
	 */
	template <class newWall>
	struct wrap :
		build<brick, newWall> {};

	/** constructor of wrapper class
	 * \tparam Args arbitrary paramerer types
	 * \param args the parameters
	 */
	template <class...Args>
	build(Args const &...args) :
		brick<wall>(args...)
	{
	}
};

#endif // BRICK_HPP_INCLUDED

