/**
 * \file matrix_block.hpp
 * \author Peter Fiala and Peter Rucz, fiala@hit.bme.hu, rucz@hit.bme.hu
 * \brief implementation of a block matrix
 */


#include <cstddef> // for size_t

/**
 * \brief proxy class to represent a block of a matrix
 * \tparam Matrix the matrix type that is indexed
 * \tparam RowIndex the row index vector type
 * \tparam ColIndex the column index vector type
 */
template <class Matrix, class RowIndex, class ColIndex = RowIndex>
class matrix_block
{
public:
	matrix_block(Matrix &m, RowIndex const &rows, ColIndex const &cols) : m_matrix(m), m_rows(rows), m_cols(cols)
	{
	}

	template <class SubMatrix>
	void operator +=(SubMatrix const &rhs) const
	{
		for (int i = 0; i < m_rows.size(); ++i)
			for (int j = 0; j < m_cols.size(); ++j)
				m_matrix(m_rows(i), m_cols(j)) += rhs(i, j);
	}

protected:
	Matrix &m_matrix;
	RowIndex const &m_rows;
	ColIndex const &m_cols;
};


/**
 * \brief factory function of matrix_block
 * \details The function allows easy creation of block instances of an arbitrary matrix object
 * \tparam Matrix the matric type
 * \tparam RowIndex the row index vector type
 * \tparam ColIndex the column index vector type
 * \param matrix [in] the matrix object to index
 * \param rows [in] the rows to be selected
 * \param cols [in] the columns to be selected
 * \return the block proxy
 */
template <class Matrix, class RowIndex, class ColIndex>
matrix_block<Matrix, RowIndex, ColIndex> block(Matrix &matrix, RowIndex const &rows, ColIndex const &cols)
{
	return matrix_block<Matrix, RowIndex, ColIndex>(matrix, rows, cols);
}

/*

#include <iostream>

template <size_t R, size_t C>
class matrix
{
private:
	double data[R][C];
public:
	double &operator()(size_t r, size_t c)
	{
		return data[r][c];
	}

	double operator()(size_t r, size_t c) const
	{
		return data[r][c];
	}

};

template <size_t N>
class vector
{
private:
	size_t data[N];
public:
	size_t operator[] (size_t i) const
	{
		return data[i];
	}

	size_t &operator[] (size_t i)
	{
		return data[i];
	}

	size_t size(void) const { return N; }
};

int main(void)
{
	matrix<10,10> big;
	matrix<2,2> small;
	vector<3> rows;
	vector<3> cols;

	rows[0] = 0;
	rows[1] = 1;
	rows[2] = 2;

	cols[0] = 0;
	cols[1] = 1;
	cols[2] = 2;

	std::cout << rows[0];

	block(big, rows, cols) += small;

	return 0;
}

*/

