/**
 * \file  collection.hpp
 * \brief Collect a set of classes into a container
 */

#include <type_traits> // is_base_of

/**
 * \brief collect a set of classes by multiple variadic inheritance
 * \tparam Args the classes that are collected into a container
 */
template <class...Args>
class collect :
	public Args...
{
public:
	/** \brief constructor from individual arguments */
	collect(Args const &...args) : Args(args)...
	{
	}

	/** \brief default constructor */
	collect()
	{
	}
};

/**
 * \brief merge two collections by preserving unique member types only
 * \tparam C1 the first collection
 * \tparam C2 the second collection
 */
template <class C1, class C2>
struct merger;


namespace internal
{
	/** \brief merge a plain class into a collection if it is already contained */
	template <class C1, class D2, bool enable = std::is_base_of<C1, D2>::value>
	struct merger_impl
	{
		typedef D2 ret_type;

		static ret_type eval(C1 const &c1, D2 const &d2)
		{
			return d2;
		}
	};


	/** \brief merge a plain class into a collection if it is not contained */
	template <class C1, class...Args2>
	struct merger_impl<C1, collect<Args2...>, false>
	{
		typedef collect<Args2...> d2_type;
		typedef collect<C1, Args2...> ret_type;

		static ret_type eval(C1 const &c1, d2_type const &d2)
		{
			return ret_type(c1, static_cast<Args2 const &>(d2)...);
		}
	};
} // end of namespace internal



/**
 * \brief merge a plain class into a collection
 * \tparam C1 the class to merge
 * \tparam Args the collection to extend is collect<Args...>
 */
template <class C1, class...Args2>
struct merger<C1, collect<Args2...> > :
	internal::merger_impl<C1, collect<Args2...> >
{
};


/**
 * \brief merge a one-element collection into an other collection
 * \tparam A1 the one-element collection is collect<A1>
 * \tparam Args the collection to extend is collect<Args...>
 */
template <class A1, class...Args2>
struct merger<collect<A1>, collect<Args2...> >
{
	typedef collect<A1> d1_type;
	typedef collect<Args2...> d2_type;
	typedef typename merger<A1, d2_type>::ret_type ret_type;

	static ret_type eval(d1_type const &d1, d2_type const &d2)
	{
		return merger<A1, d2_type>::eval(d1, d2);
	}
};


/**
 * \brief merge a general collection into an other collection
 * \tparam A1 the one-element collection is collect<A1>
 * \tparam Args the collection to extend is collect<Args...>
 */
template <class A1, class...Args1, class...Args2>
struct merger<collect<A1, Args1...>, collect<Args2...> >
{
	typedef collect<A1, Args1...> d1_type;
	typedef collect<Args2...> d2_type;

	typedef typename merger<collect<Args1...>, d2_type>::ret_type inter_type;
	typedef typename merger<A1, inter_type>::ret_type ret_type;

	static ret_type eval(d1_type const &d1, d2_type const &d2)
	{
		return merger<A1, inter_type>::eval(d1,
			merger<collect<Args1...>, d2_type >::eval(d1, d2)
		);
	}
};


/** \brief factory function to merge two collections
 * \tparam Args1 the first collection members
 * \tparam Args2 the second collection members
 * \param [in] c1 the first collection
 * \param [in] c2 the second collection
 * \return the merged collection
 */
template <class...Args1, class...Args2>
typename merger<collect<Args1...>, collect<Args2...> >::ret_type
	merge_data(collect<Args1...> const &c1, collect<Args2...> const &c2)
{
	return merger<collect<Args1...>, collect<Args2>...>::eval(c1, c2);
}


