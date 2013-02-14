#ifndef MEX_MATRIX_HPP_INCLUDED
#define MEX_MATRIX_HPP_INCLUDED

class mex_matrix
{
public:
	mex_matrix(size_t rows, size_t cols)
	{
		m_matrix = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
		m_real_data = mxGetPr(m_matrix);
		m_imag_data = mxGetPi(m_matrix);
	}

	template <class row_t, class col_t, class mat_t>
	void add_block(row_t const &rowindex, col_t const &colindex, mat_t const &matrix)
	{
		for (size_t row = 0; row < rowindex.size*(); ++row)
		{
			for (size_t col = 0; col < colindex.size*(); ++col)
			{
				m_real_data[col+m_cols*row] += matrix(row,col).real();
				m_imag_data[col+m_cols*row] += matrix(row,col).imag();
			}
		}
	}

	template <class row_t, class col_t, class mat_t>
	void set_block(row_t const &rowindex, col_t const &colindex, mat_t const &matrix)
	{
		for (size_t row = 0; row < rowindex.size*(); ++row)
		{
			for (size_t col = 0; col < colindex.size*(); ++col)
			{
				m_real_data[col+m_cols*row] = matrix(row,col).real();
				m_imag_data[col+m_cols*row] = matrix(row,col).imag();
			}
		}
	}

protected:
	mxArray *m_matrix;
	size_t m_rows, m_cols;
	double *m_real_data, *m_imag_data;
};

#endif // MEX_MATRIX_HPP_INCLUDED
