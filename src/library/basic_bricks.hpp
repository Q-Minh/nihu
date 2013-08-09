#ifndef BASIC_BRICKS_HPP_INCLUDED
#define BASIC_BRICKS_HPP_INCLUDED

/** \brief a brick representing a distance vector \f${\bf r} = {\bf y} - {\bf x}\f$
 * \tparam space the coordinate space the distance is defined over
 */
template <class space>
struct distance_vector_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_distance_vector(trial_input.get_x()-test_input.get_x())
		{
		}

		/** \brief return distance vector
		 * \return distance vector
		 */
		typename space::location_t const &get_distance_vector(void) const
		{
			return m_distance_vector;
		}

	private:
		typename space::location_t m_distance_vector;
	};
};


/** \brief a brick representing a scalar distance \f$r = |{\bf r}|\f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct distance_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_distance(wall::get_distance_vector().norm())
		{
		}

		/** \brief return distance
		 * \return scalar distance
		 */
		scalar const & get_distance(void) const
		{
			return m_distance;
		}

	private:
		scalar m_distance;
	};
};


/** \brief a brick representing a distance normal derivative \f$ r'_{n_y} \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct rdny_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_rdny(wall::get_distance_vector().dot(trial_input.get_unit_normal()) / wall::get_distance())
		{
		}

		/** \brief return distance
		 * \return scalar distance
		 */
		scalar const & get_rdny(void) const
		{
			return m_rdny;
		}

	private:
		scalar m_rdny;
	};
};


/** \brief a brick representing a distance normal derivative \f$ r'_{n_x} \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct rdnx_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_rdnx(-1.0 * wall::get_distance_vector().dot(test_input.get_unit_normal()) / wall::get_distance())
		{
		}

		/** \brief return distance
		 * \return scalar distance
		 */
		scalar const & get_rdnx(void) const
		{
			return m_rdnx;
		}

	private:
		scalar m_rdnx;
	};
};


#endif // BASIC_BRICKS_HPP_INCLUDED
