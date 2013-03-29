/**
 * \file mex_matrix.hpp
 * \brief A Matlab mex matrix interface
 * \details The interface makes it possible to construct real and complex mex-type
 * matrices in native Matlab format, and use them in C++ as if they were
 * Eigen-type matrices.
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */
#ifndef MEX_MATRIX_HPP_INCLUDED
#define MEX_MATRIX_HPP_INCLUDED

#include <mex.h>

/** \brief container class of a real matrix stored in Matlab format */
class mex_real_matrix
{
public:
	/** \brief constructor
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 */
	mex_real_matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols)
	{
		m_matrix = mxCreateDoubleMatrix(m_rows, m_cols, mxREAL);
		m_real = mxGetPr(m_matrix);
	}


	/** \brief constructor
	 * \param [in] mex_matrix pointer to the native Matlab matrix format
	 */
	mex_real_matrix(mxArray const *mex_matrix) : m_matrix(mex_matrix)
	{
		m_rows = mxGetM(m_matrix);
		m_cols = mxGetN(m_matrix);
		m_real = mxGetPr(m_matrix);
	}

	double const &operator() (size_t row, size_t col) const
	{
		return m_real[row + m_rows * col];
	}


	double operator() (size_t row, size_t col)
	{
		return m_real[row + m_rows * col];
	}

	size_t rows(void) const { return m_rows; }
	size_t cols(void) const { return m_cols; }


	/** \brief return the Matlab matrix
	 * \return the matrix in Matlab format
	 */
	mxArray const *get_matrix() const
	{
		return m_matrix;
	}

protected:
	mxArray const *m_matrix;	/**< \brief Matlab matrix */
	size_t m_rows;			/**< \brief number of rows */
	size_t m_cols;			/**< \brief number of columns */
	double *m_real;		/**< \brief array of real data */
};



/** \brief container class of a complex matrix stored in Matlab format */
class mex_complex_matrix
{
public:
	/** \brief index proxy class of a complex matrix */
	class index_proxy
	{
	public:
		/** \brief constructor
		 * \param matrix the parent container class
		 * \param row the row index
		 * \param col the column index
		 */
		index_proxy(mex_complex_matrix &matrix, size_t row, size_t col)
			: m_matrix(matrix), m_row(row), m_col(col) {}

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
		size_t const m_row; /**< \brief the row index */
		size_t const m_col; /**< \brief the column index */
	};


	/** \brief constructor
	 * \param rows number of rows
	 * \param cols number of columns
	 */
	mex_complex_matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols)
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
	void increment(size_t row, size_t col, data_t const &data)
	{
		m_real[row+m_rows*col] += data.real();
		m_imag[row+m_rows*col] += data.imag();
	}

	/** \brief index operator that returns a proxy
	 * \param row row index
	 * \param col column index
	 * \return proxy object storing the parent container and the indices
	 */
	index_proxy operator() (size_t row, size_t col)
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
	size_t m_rows;	/**< \brief number of rows */
	size_t m_cols;	/**< \brief number of columns */
	mxArray *m_matrix;	/**< \brief Matlab matrix */
	double *m_real;		/**< \brief array of real data */
	double *m_imag;		/**< \brief array of imaginary data */
};


#endif // MEX_MATRIX_HPP_INCLUDED

