/**
 * \file mex_matrix.hpp
 * \brief A Matlab mex complex matrix interface
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */
#ifndef MEX_MATRIX_HPP_INCLUDED
#define MEX_MATRIX_HPP_INCLUDED

#include <mex.h>

/** \brief container class of a complex matrix stored in Matlab format */
class mex_complex_matrix;

/** \brief index proxy class of a complex matrix */
class index_proxy
{
public:
	/** \brief constructor
	 * \param matrix the parent container class
	 * \param row the row index
	 * \param col the column index
	 */
	index_proxy(mex_complex_matrix &matrix, int row, int col) : m_matrix(matrix), m_row(row), m_col(col) {}

	/** \brief increment operator
	 * \param data the data to add to the conatiner
	 */
	template <class data_t>
	void operator +=(data_t const &data) const
	{
		m_matrix.increment(m_row, m_col, data);
	}

private:
	mex_complex_matrix &m_matrix; /**< \brief reference to the parent */
	int const m_row; /**< \brief the row index */
	int const m_col; /**< \brief the column index */
};


/** \brief container class of a complex matrix stored in Matlab format */
class mex_complex_matrix
{
public:
	/** \brief constructor
	 * \param rows number of rows
	 * \param cols number of columns
	 */
	mex_complex_matrix(int rows, int cols) : m_rows(rows), m_cols(cols)
	{
		m_matrix = mxCreateDoubleMatrix(m_rows, m_cols, mxCOMPLEX);
		m_real = mxGetPr(m_matrix);
		m_imag = mxGetPi(m_matrix);
	}

	/** \brief increment a given element of the matrix with a complex data
	 * \param row row index
	 * \param col column index
	 * \param data the data to increase our matrix with
	 */
	template <class data_t>
	void increment(int row, int col, data_t data)
	{
		m_real[row+m_rows*col] += data.real();
		m_imag[row+m_rows*col] += data.imag();
	}

	/** \brief index operator that returns a proxy
	 * \param row row index
	 * \param col column index
	 * \return proxy object storing the parent container and the indices
	 */
	index_proxy operator() (int row, int col)
	{
		return index_proxy(*this, row, col);
	}

	/** \brief return the Matlab matrix
	 * \return the matrix in Matlab format
	 */
	mxArray *get_matrix() const
	{
		return m_matrix;
	}

protected:
	int const m_rows;	/**< \brief number of rows */
	int const m_cols;	/**< \brief number of columns */
	mxArray *m_matrix;	/**< \brief Matlab matrix */
	double *m_real;		/**< \brief array of real data */
	double *m_imag;		/**< \brief array of imaginary data */
};

#endif // MEX_MATRIX_HPP_INCLUDED
