/**
 * \file matrix_block.hpp
 * \author Peter Fiala and Peter Rucz, fiala@hit.bme.hu, rucz@hit.bme.hu
 * \brief implementation of a block matrix
 */

template <class Matrix, class Index>
class matrix_block
{
public:
	matrix_block(Matrix &m, Index const &rows, Index const &cols) : m_matrix(m), m_rows(rows), m_cols(cols)
	{
	}
	
	template <class SubMatrix>
	void operator +=(SubMatrix const &rhs) const
	{
		for (size_t i = 0; i < m_rows.size(); ++i)
			for (size_t j = 0; j < m_cols.size(); ++j)
				m_matrix(m_rows[i], m_cols[j]) += rhs(i, j);
	}
	
protected:
	Matrix &m_matrix;
	Index const &m_rows;
	Index const &m_cols;
};


#include <iostream>

template <class Matrix, class Index>
matrix_block<Matrix, Index> block(Matrix &matrix, Index const &rows, Index const &cols)
{
	return matrix_block<Matrix, Index>(matrix, rows, cols);
}

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
