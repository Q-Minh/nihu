#ifndef MEX_MATRIX_HPP_INCLUDED
#define MEX_MATRIX_HPP_INCLUDED

#include <mex.h>

/** \brief container class of N complex matrices stored in Matlab format
 * \tparam N number of complex matrices
 */
template <size_t N>
class mex_complex_matrix;

/** \brief index proxy class of N complex matrices
 * \tparam N number of complex matrices
 */
template <size_t N>
class index_proxy
{
public:
	/** \brief constructor
	 * \param cont the parent container class
	 * \param row the row index
	 * \param col the column index
	 */
	index_proxy(mex_complex_matrix<N> cont, size_t row, size_t col) : m_parent(cont), m_row(row), m_col(col) {}

	/** \brief increase operator
	 * \param data the data to add to the conatiner
	 */
	template <class data_t>
	void operator +=(data_t const &data)
	{
		m_parent.increase(m_row, m_col, data);
	}

private:
	mex_complex_matrix<N> & m_parent; /**< \brief reference to the parent */
	size_t const m_row; /**< \brief the row index */
	size_t const m_col; /**< \brief the column index */
};


/** \brief container class of N complex matrices stored in Matlab format
 * \tparam N number of complex matrices
 */
template <size_t N>
class mex_complex_matrix
{
public:
	/** \brief constructor
	 * \param rows number of rows
	 * \param cols number of columns
	 */
	mex_complex_matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols)
	{
		for (size_t i = 0; i < N; ++i)
		{
			m_matrix[i] = mxCreateDoubleMatrix(m_rows, m_cols, mxCOMPLEX);
			m_real[i] = mxGetPr(m_matrix[i]);
			m_imag[i] = mxGetPi(m_matrix[i]);
		}
	}

	/** \brief increase a given element of all matrices with a vector data
	 * \param row row index
	 * \param col column index
	 * \param data the data to increase our matrix with
	 */
	template <class data_t>
	void increase(size_t row, size_t col, data_t data)
	{
		for (size_t i = 0; i <N; ++i)
		{
			m_real[i][row+m_rows*col] += data[i].real();
			m_imag[i][row+m_rows*col] += data[i].imag();
		}
	}

	/** \brief index operator that returns a proxy
	 * \param row row index
	 * \param col column index
	 * \return proxy object storing the parent container and the indices
	 */
	index_proxy<N> operator() (size_t row, size_t col)
	{
		return index_proxy<N>(*this, row, col);
	}

	/** \brief return the i-th Matlab matrix
	 * \param i matrix index
	 * \return the i-th matrix in Matlab format
	 */
	mxArray *get_matrix(size_t i = 0) const
	{
		return m_matrix[i];
	}

protected:
	size_t const m_rows; /**< \brief number of rows */
	size_t const m_cols; /**< \brief number of columns */
	mxArray *m_matrix[N]; /**< \brief array of Matlab matrices */
	double *m_real[N];	/**< \brief arrays of real data */
	double *m_imag[N];	/**< \brief arrays of imaginary data */
};

#endif // MEX_MATRIX_HPP_INCLUDED
