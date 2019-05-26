#ifndef MATRIX_H
#define MATRIX_H

////////////////////////////////////////////////////////////////////////////////
//
// FILE:        Matrix.h
// DESCRIPTION: contains template matrix class and implementation
// AUTHOR:      Dan Fabian
// DATE:        5/24/2019

#include <vector>
#include <thread>
#include <utility>
#include <cmath>
#include <iostream>
#include <cassert>

using std::vector;
using std::thread;
using std::pair; using std::make_pair;
using std::cout; using std::endl;

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX
template <typename T = int> // T must be some type of number
class Matrix {
public:
	// constructors
	Matrix();                                    // default
	Matrix(size_t n);                            // identity matrix of size n x n
	Matrix(size_t rows, size_t cols, T val = T(0)); // cells default to 0's
	Matrix(const Matrix<T>& orig);               // copy constructor
	Matrix(const vector<vector<T>>& orig);       // explicit constructor
	Matrix(const vector<T>& orig);               // vector constructor

	// binary operators
	Matrix<T>  operator*  (const T& rhs)          const;  // scalar multiplication
	Matrix<T>  operator*  (const Matrix<T>& rhs)  const;  // matrix multiplication, uses multi threading
	Matrix<T>  operator+  (const T& rhs)          const;  // adds rhs to every cell
	Matrix<T>  operator+  (const Matrix<T>& rhs)  const;  // matrix addition
	Matrix<T>  operator-  (const Matrix<T>& rhs)  const { return *this + (rhs * T(-1)); }
	Matrix<T>  operator-  (const T & rhs)         const { return *this + (rhs * T(-1)); }
	Matrix<T>  operator^  (const size_t & rhs)    const;  // matrix power
	vector<T>& operator[] (const size_t & val);           // used to adjust cells
	vector<T>  operator[] (const size_t & val)    const;  // const cell retrieval
	Matrix<T>& operator=  (const Matrix<T> & rhs);        // copies matrix
	bool       operator== (const Matrix<T> & rhs) const;
	bool       operator!= (const Matrix<T> & rhs) const;

	// unary operators
	Matrix<T>& operator+= (const Matrix<T> & rhs)     { *this = *this + rhs; return *this; }
	Matrix<T>& operator+= (const T & rhs)             { *this = *this + rhs; return *this; }
	Matrix<T>& operator*= (const Matrix<T> & rhs)     { *this = *this * rhs; return *this; }
	Matrix<T>& operator*= (const T & rhs)             { *this = *this * rhs; return *this; }
	Matrix<T>& operator-= (const Matrix<T> & rhs)     { *this = *this - rhs; return *this; }
	Matrix<T>& operator-= (const T & rhs)             { *this = *this - rhs; return *this; }
	Matrix<T>& operator^= (const T & rhs)             { *this = *this ^ rhs; return *this; }
	Matrix<T>  operator-  ()                    const { return *this* T(-1); }
	Matrix<T>  operator+  ()                    const { return *this; }

	// methods
	void                       print             ()                   const; // prints all cells
	void                       swap              (Matrix<T> & mat);           // constant time swap
	size_t                     getRows           ()                   const { return rows_; }
	size_t                     getCols           ()                   const { return cols_; }
	void                       multiplyRow       (size_t row, T val);        // multiplies each element in row by val
	void                       addRow            (size_t row, T val);        // adds a val to each element in row
	void                       rowOp             (size_t row, size_t rowToChange, T val);
	void                       rowSwap           (size_t row1, size_t row2); // swaps two rows
	void                       transpose         ();                         // transposes matrix
	Matrix<T>                  subMatrix         (size_t rowBegin, size_t rowEnd, size_t colBegin, size_t colEnd) const;
	Matrix<T>                  rowEchelon        ()                   const; // optimized for square matrices, returns matrix in row Echelon form
	Matrix<T>                  rowEchelonAnySize ()                   const; // returns matrix in row Echelon form
	pair<Matrix<T>, Matrix<T>> rowEchelon        (Matrix<T> ans)      const; // returns matrix in row Echelon form with an answers matrix
	T                          determinant       ()                   const; // calcs determinant with row echelon
	Matrix<T>                  inverse           ()                   const; // finds inverse of square matrix
	size_t                     rank              ()                   const; // finds rank of matrix

	// friends
	template <typename K>
	friend bool      checkDimSame(const Matrix<K> & left, const Matrix<K> & right);
	template <typename K>
	friend bool      checkDimMultiply(const Matrix<K> & left, const Matrix<K> & right);
	template <typename K>
	friend Matrix<K> operator*        (const K & left, const Matrix<K> & right);
	template <typename K>
	friend Matrix<K> operator+        (const K & left, const Matrix<K> & right);
	template <typename K>
	friend Matrix<K> operator-        (const K& left, const Matrix<K>& right);

private:
	vector<vector<T>> cells_;
	size_t            rows_;
	size_t            cols_;

};

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX constructors
////////////////////////////////////////
// completely empty matrix
template <typename T>
Matrix<T>::Matrix()
	: rows_(0), cols_(0) {}

////////////////////////////////////////
// creates an identity matrix of size n x n
template <typename T>
Matrix<T>::Matrix(size_t n)
	: rows_(n), cols_(n)
{
	for (size_t i = 0; i != n; ++i)
	{
		vector<T> tmp;
		for (size_t j = 0; j != n; ++j)
		{
			if (j == i) // diagnals
				tmp.push_back(T(1));
			else
				tmp.push_back(T(0));
		}
		cells_.push_back(tmp);
	}
}

////////////////////////////////////////
// creates a matrix with special size, initialized with a value
template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, T val)
	: rows_(rows), cols_(cols)
{
	for (size_t i = 0; i != rows_; ++i)
	{
		vector<T> tmp;
		for (size_t j = 0; j != cols_; ++j) tmp.push_back(val);
		cells_.push_back(tmp);
	}
}

////////////////////////////////////////
// copy constructor
template <typename T>
Matrix<T>::Matrix(const Matrix<T> & orig)
	: rows_(orig.rows_), cols_(orig.cols_), cells_(orig.cells_) {}

////////////////////////////////////////
// explicit constructor
template <typename T>
Matrix<T>::Matrix(const vector<vector<T>> & orig)
	: rows_(orig.size()), cols_(orig[0].size()), cells_(orig)
{
	// check that all columns are same size
	for (const auto& e : orig)
		assert(e.size() == cols_);
}

////////////////////////////////////////
// vector constructor
template <typename T>
Matrix<T>::Matrix(const vector<T> & orig)
	: rows_(orig.size()), cols_(1)
{
	for (size_t i = 0; i != rows_; ++i)
		cells_.push_back(vector<T>(1, orig[i]));
}

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX operators
////////////////////////////////////////
// used to adjust cells, example: matrix[3] returns fourth row
template <typename T>
vector<T>& Matrix<T>::operator[](const size_t & val)
{
	return cells_[val];
}

////////////////////////////////////////
// used to return cells, example: matrix[3] returns fourth row
template <typename T>
vector<T> Matrix<T>::operator[](const size_t & val) const
{
	return cells_[val];
}

////////////////////////////////////////
// overloaded assignment
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> & rhs)
{
	cells_ = rhs.cells_;
	rows_ = rhs.rows_;
	cols_ = rhs.cols_;
	return *this;
}

////////////////////////////////////////
// matrix addition
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> & rhs) const
{
	// make sure dims are the same for addition
	assert(checkDimSame(*this, rhs));

	Matrix<T> result(rows_, cols_);
	for (size_t i = 0; i != rows_; ++i)
		for (size_t j = 0; j != cols_; ++j)
			result[i][j] = (*this)[i][j] + rhs[i][j];

	return result;
}

////////////////////////////////////////
// addition with a single number
template <typename T>
Matrix<T> Matrix<T>::operator+(const T & rhs) const
{
	return *this + Matrix<T>(this->rows_, this->cols_, rhs);
}

////////////////////////////////////////
// matrix multiplication, uses multi threading
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> & rhs) const
{
	// make sure dims are good to multiply
	assert(checkDimMultiply(*this, rhs));

	Matrix<T> result(rows_, rhs.cols_);

	// lambda for calculating a range of rows in the result matrix, [rowBegin, rowEnd)
	auto calcRows = [&](size_t rowBegin, size_t rowEnd) {

		for (size_t row = rowBegin; row != rowEnd; ++row)
			for (size_t col = 0; col != result.cols_; ++col)
				for (size_t i = 0; i != cols_; ++i)
					result[row][col] += (*this)[row][i] * rhs[i][col];
	};

	// if both matrices have a dimension greater than num use multi threading
	size_t num = 10;
	if ((rows_ > num || cols_ > num) &&
		(rhs.getRows() > num || rhs.getCols() > num))
	{
		size_t threadsToUse = thread::hardware_concurrency();

		// divide workload
		vector<thread> threads;
		if (result.getRows() <= threadsToUse)
			for (size_t i = 0; i != result.getRows(); ++i)
				threads.push_back(thread(calcRows, i, i + 1));

		else // if more rows than threads to use divide up work
		{
			size_t rowsPerThread = static_cast<size_t>(floor(result.getRows() / threadsToUse));
			size_t leftOverRows = result.getRows() - (rowsPerThread * threadsToUse);

			// assign threads
			size_t i = 0;
			for (; i != threadsToUse; ++i)
				threads.push_back(thread(calcRows, i * rowsPerThread, (i + 1) * rowsPerThread));

			threads.push_back(thread(calcRows, i * rowsPerThread, i * rowsPerThread + leftOverRows));
		}

		// synchronize
		for (size_t i = 0; i != threads.size(); ++i)
			threads[i].join();
	}
	else // multi threading isn't necessary
		calcRows(0, result.getRows());

	return result;
}

////////////////////////////////////////
// scalar multiplication
template <typename T>
Matrix<T> Matrix<T>::operator*(const T & rhs) const
{
	Matrix<T> result(*this);
	for (size_t row = 0; row != result.rows_; ++row)
		for (size_t col = 0; col != result.cols_; ++col)
			result[row][col] *= rhs;

	return result;
}
////////////////////////////////////////
// comparing matrices
template <typename T>
bool Matrix<T>::operator==(const Matrix<T> & rhs) const
{
	if (!checkDimSame(*this, rhs))
		return false;

	for (size_t i = 0; i != rows_; ++i)
		for (size_t j = 0; j != cols_; ++j)
			if ((*this)[i][j] != rhs[i][j])
				return false;

	return true;
}

////////////////////////////////////////
// matrix power 
template <typename T>
Matrix<T> Matrix<T>::operator^(const size_t & rhs) const
{
	// must be square
	assert(cols_ == rows_);

	Matrix<T> result(rows_);
	for (size_t i = 0; i != rhs; ++i)
		result = result * (*this);

	return result;
}

////////////////////////////////////////
// comparing matrices
template <typename T>
bool Matrix<T>::operator!=(const Matrix<T> & rhs) const
{
	return !(*this == rhs);
}

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX methods
////////////////////////////////////////
// prints out matrix, used for testing
template <typename T>
void Matrix<T>::print() const
{
	for (size_t i = 0; i != rows_; ++i)
	{
		for (size_t j = 0; j != cols_; ++j) cout << cells_[i][j] << ' ';
		cout << endl;
	}
}

////////////////////////////////////////
// constant time swap
template <typename T>
void Matrix<T>::swap(Matrix<T> & mat)
{
	Matrix<T> tmp = mat;
	mat = *this;
	*this = tmp;
}

////////////////////////////////////////
// multiplies a whole row by a value
template <typename T>
void Matrix<T>::multiplyRow(size_t row, T val)
{
	for (auto& e : cells_[row])
		e *= val;
}

////////////////////////////////////////
// adds a value to a whole row
template <typename T>
void Matrix<T>::addRow(size_t row, T val)
{
	for (auto& e : cells_[row])
		e += val;
}

////////////////////////////////////////
// multiplies a row by a val then adds column-wise to the row to be changed
template <typename T>
void Matrix<T>::rowOp(size_t row, size_t rowToChange, T val)
{
	vector<T> tmp = cells_[row];
	for (auto& e : tmp)
		e *= val;

	for (size_t i = 0; i != cols_; ++i)
		cells_[rowToChange][i] += tmp[i];
}

////////////////////////////////////////
// swaps two rows in a matrix
template <typename T>
void Matrix<T>::rowSwap(size_t row1, size_t row2)
{
	if (row1 == row2) return;

	vector<T> tmp = (*this)[row1];
	(*this)[row1] = (*this)[row2];
	(*this)[row2] = tmp;
}

////////////////////////////////////////
// transposes matrix
template <typename T>
void Matrix<T>::transpose()
{
	Matrix<T> result(cols_, rows_);
	for (size_t row = 0; row != result.rows_; ++row)
		for (size_t col = 0; col != result.cols_; ++col)
			result[row][col] = (*this)[col][row];

	*this = result;
}

////////////////////////////////////////
// creates a sub-matrix, inclusive 
template <typename T>
Matrix<T> Matrix<T>::subMatrix(size_t rowBegin, size_t rowEnd, size_t colBegin, size_t colEnd) const
{
	assert(rowEnd > rowBegin);
	assert(colEnd > colBegin);

	Matrix<T> result(rowEnd - rowBegin + 1, colEnd - colBegin + 1);

	for (size_t i = 0; i != result.getRows(); ++i)
		for (size_t j = 0; j != result.getCols(); ++j)
			result[i][j] = (*this)[i + rowBegin][j + colBegin];

	return result;
}

////////////////////////////////////////
// returns matrix in row Echelon form, need to use a type that is able to work with fractional parts
// this version is optimized for square matrices
template <typename T>
Matrix<T> Matrix<T>::rowEchelon() const
{
	// must be square
	assert(rows_ == cols_);

	// make sure the matrix isn't just 0's
	assert((*this) != Matrix<T>(rows_, cols_, T(0)));

	Matrix<T> result(*this);
	for (size_t i = 0; i != rows_; ++i)
	{
		// check if rows need to be swapped
		// if leading term is 0, search through all rows to find a leading term not 0
		size_t j = i;
		for (; j != rows_ && result[j][i] == T(0); ++j) {}

		// check if j reached the end of all rows, then matrix can't be solved so throw error
		if (j == rows_) return result;

		// then swap
		result.rowSwap(i, j); // if i == j, nothing happens

		for (size_t row = i + 1; row != rows_; ++row)
		{
			T rowMultiplier = (result[row][i] / result[i][i]) * T(-1);
			result.rowOp(i, row, rowMultiplier);
		}
	}

	return result;
}

////////////////////////////////////////
// returns matrix in row Echelon form, need to use a type that is able to work with fractional parts
template <typename T>
Matrix<T> Matrix<T>::rowEchelonAnySize() const
{
	// make sure the matrix isn't just 0's
	assert((*this) != Matrix<T>(rows_, cols_, T(0)));

	Matrix<T> result(*this);
	for (size_t i = 0; i != rows_; ++i)
	{
		// check if rows need to be swapped
		// also search through all rows and columns to find a leading term 
		size_t leadingTerm = cols_ - 1;
		size_t rowToSwap = i;
		bool termFound = false;
		for (size_t j = i; j != cols_ && !termFound; ++j)
			for (size_t k = i; k != rows_ && !termFound; ++k)
				if (j <= leadingTerm && result[k][j] != T(0)) // if leading term found
				{
					termFound = true;
					leadingTerm = j;
					rowToSwap = k;
				}

		result.rowSwap(i, rowToSwap); // then swap rows if needed

		if (result[i][leadingTerm] == T(0)) return result;

		for (size_t row = i + 1; row != rows_; ++row)
		{
			T rowMultiplier = (result[row][leadingTerm] / result[i][leadingTerm]) * T(-1);
			result.rowOp(i, row, rowMultiplier);
		}
	}

	return result;
}

////////////////////////////////////////
// returns matrix in row Echelon form, need to use a type that is able to work with fractional parts
template <typename T>
pair<Matrix<T>, Matrix<T>> Matrix<T>::rowEchelon(Matrix<T> answers) const
{
	Matrix<T> equations(*this);

	// check if the matrix can actually be solved
	assert(cols_ == rows_);
	assert(answers.getRows() == equations.getRows());

	// make sure the matrix isn't just 0's
	assert(equations != Matrix<T>(equations.getRows(), equations.getCols(), T(0)));

	for (size_t i = 0; i != equations.getRows(); ++i)
	{
		// check if rows need to be swapped
		// if leading term is 0, search through all rows to find a leading term not 0
		size_t j = i;
		for (; j != equations.getRows() && equations[j][i] == T(0); ++j) {}

		// check if j reached the end of all rows, then matrix can't be solved so throw error
		assert(j != equations.getRows());

		// then swap
		equations.rowSwap(i, j); // if i == j, nothing happens
		answers.rowSwap(i, j);

		for (size_t row = i + 1; row != equations.getRows(); ++row)
		{
			T rowMultiplier = (equations[row][i] / equations[i][i]) * T(-1);
			equations.rowOp(i, row, rowMultiplier);
			answers.rowOp(i, row, rowMultiplier);
		}
	}

	return make_pair(equations, answers);
}

////////////////////////////////////////
// calculates determinant using Row Echelon
template <typename T>
T Matrix<T>::determinant() const
{
	// only for square matrices
	assert(rows_ == cols_);

	// check if its an identity matrix, then it's just 1
	if ((*this) == Matrix<T>(rows_)) return T(1);

	T determ = T(1);

	Matrix<T> result(*this);
	for (size_t i = 0; i != rows_; ++i)
	{
		// check if rows need to be swapped
		// if leading term is 0, search through all rows to find a leading term not 0
		size_t j = i;
		for (; j != rows_ && result[j][i] == T(0); ++j) {}

		// then we are done so exit loop
		if (j == rows_) break;

		// then swap
		if (i != j)
		{
			determ *= T(-1);
			result.rowSwap(i, j);
		}

		for (size_t row = i + 1; row != rows_; ++row)
		{
			T rowMultiplier = (result[row][i] / result[i][i]) * T(-1);
			result.rowOp(i, row, rowMultiplier);
		}
	}

	for (size_t i = 0; i != result.getRows(); ++i)
		determ *= result[i][i];

	return determ;
}

////////////////////////////////////////
// finds inverse of a square matrix
template <typename T>
Matrix<T> Matrix<T>::inverse() const
{
	// only for square matrices
	assert(rows_ == cols_);

	// check if inverse exists
	pair<Matrix<T>, Matrix<T>> echelonPair = this->rowEchelon(Matrix<T>(rows_));
	Matrix<T> origReduced = echelonPair.first;
	Matrix<T> result = echelonPair.second;
	for (size_t i = 0; i != origReduced.getRows(); ++i)
	{
		size_t j = 0;
		for (; j != origReduced.getRows() && origReduced[i][j] == T(0); ++j) {} // if a non zero value is found move to the next row
		assert(j != origReduced.getRows());
	}

	for (size_t i = origReduced.getRows() - 1; i >= 0 &&
		i < origReduced.getRows(); --i) // need to check that i < rows because size_t 
										// cycles back to positive if subtracting from 0
	{
		// normalize row
		result.multiplyRow(i, static_cast<T>(1) / origReduced[i][i]); // must normalize result first
		origReduced.multiplyRow(i, static_cast<T>(1) / origReduced[i][i]);

		for (size_t row = i - 1; row < origReduced.getRows(); --row)
		{
			T rowMultiplier = (origReduced[row][i] / origReduced[i][i]) * T(-1);
			origReduced.rowOp(i, row, rowMultiplier);
			result.rowOp(i, row, rowMultiplier);
		}
	}

	return result;
}

////////////////////////////////////////
// finds rank of a matrix
template <typename T>
size_t Matrix<T>::rank() const
{
	Matrix<T> origReduced = this->rowEchelonAnySize();
	size_t rank = 0;
	for (; rank != origReduced.getRows(); ++rank)
	{
		size_t j = 0;
		for (; j != origReduced.getCols() && origReduced[rank][j] == T(0); ++j) {} // if a non zero value is found move to the next row
		if (j == origReduced.getCols()) return rank;
	}

	return rank;
}

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX friends
////////////////////////////////////////
// checks if dimensions are the same
template <typename T>
bool checkDimSame(const Matrix<T> & left, const Matrix<T> & right)
{
	return left.rows_ == right.rows_ && left.cols_ == right.cols_;
}

////////////////////////////////////////
// checks if dimensions are good to multiply
template <typename T>
bool checkDimMultiply(const Matrix<T> & left, const Matrix<T> & right)
{
	return left.cols_ == right.rows_;
}

////////////////////////////////////////
// allows for multiplying by numbers in different orders
template <typename T>
Matrix<T> operator*(const T & left, const Matrix<T> & right)
{
	return right * left;
}

////////////////////////////////////////
// allows for adding by numbers in different orders
template <typename T>
Matrix<T> operator+(const T & left, const Matrix<T> & right)
{
	return right + left;
}

////////////////////////////////////////
// allows for subtracting by numbers in different orders
template <typename T>
Matrix<T> operator-(const T& left, const Matrix<T>& right)
{
	return (right - left) * T(-1);
}

#endif // MATRIX_H