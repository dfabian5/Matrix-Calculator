#ifndef SOLVER_H
#define SOLVER_H

////////////////////////////////////////////////////////////////////////////////
//
// FILE:        Calculator.h
// DESCRIPTION: contains functions for solving a system of linear equations, and UI for calculating
// AUTHOR:      Dan Fabian
// DATE:        5/25/2019

#include "Matrix.h"
#include "Fraction.h"
#include <chrono>
#include <utility>
#include <cassert>
#include <complex>
#include <iostream>

using namespace std::chrono;
using std::pair; using std::make_pair;
using std::complex;
using std::cin; using std::cout; using std::endl;

typedef high_resolution_clock Clock;

////////////////////////////////////////////////////////////////////////////////
//
// FUNCTIONS
////////////////////////////////////////
// creates a pair of square matrices, row echelon form and the corresponding answers matrix.
// input equations matrix must be square
template <typename T>
pair<Matrix<T>, Matrix<T>> rowEchelon(Matrix<T> equations, Matrix<T> answers)
{
	// check if the matrix can actually be solved
	assert(equations.getCols() == equations.getRows());
	assert(answers.getRows() == equations.getRows());
	assert(answers.getCols() == 1);
	
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
// returns a vector of all answers in order, uses rowEchelon function
// input equations matrix must be square
template <typename T>
vector<T> solve(Matrix<T> equations, Matrix<T> answers)
{
	pair<Matrix<T>, Matrix<T>> result = rowEchelon(equations, answers);
	Matrix<T> adjustedEqs = result.first;
	Matrix<T> adjustedAns = result.second;

	for (size_t i = equations.getRows() - 1; i >= 0 && 
		 i < equations.getRows(); --i) // need to check that i < rows because size_t 
									   // cycles back to positive if subtracting from 0
	{
		// normalize row
		adjustedAns.multiplyRow(i, static_cast<T>(1) / adjustedEqs[i][i]); // must normalize Ans first
		adjustedEqs.multiplyRow(i, static_cast<T>(1) / adjustedEqs[i][i]);

		for (size_t row = i - 1; row < equations.getRows(); --row)
		{
			T rowMultiplier = (adjustedEqs[row][i] / adjustedEqs[i][i]) * T(-1);
			adjustedEqs.rowOp(i, row, rowMultiplier);
			adjustedAns.rowOp(i, row, rowMultiplier);
		}
	}

	vector<T> finalAns;
	for (size_t i = 0; i != adjustedAns.getRows(); ++i)
		finalAns.push_back(adjustedAns[i][0]);

	return finalAns;
}

////////////////////////////////////////////////////////////////////////////////
//
// HELPER CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// get input from square matrix
template <typename T>
Matrix<T> input(bool isSquare = true)
{
	// get size of matrix
	size_t row;
	size_t col;
	if (!isSquare)
	{
		cout << "Enter number of rows in Matrix: ";
		cin >> row;
		cout << "Enter number of columns in Matrix: ";
		cin >> col;
		cout << endl;
	}
	else
	{
		cout << "Enter size of Matrix: ";
		cin >> row;
		cout << endl;

		col = row;
	}

	// enter matrix
	Matrix<T> mat(row, col);
	for (size_t i = 0; i != row; ++i)
	{
		cout << "ENTER ROW " << i << ':' << endl;
		for (size_t j = 0; j != col; ++j)
		{
			cout << "Enter value in column " << j << ": ";
			cin >> mat[i][j];
		}
		cout << endl;
	}

	return mat;
}

////////////////////////////////////////
// get input from square matrix with rational complex specialization
template <>
Matrix<complex<Fraction<int>>> input<complex<Fraction<int>>>(bool isSquare)
{
	// get size of matrix
	size_t row;
	size_t col;
	if (!isSquare)
	{
		cout << "Enter number of rows in Matrix: ";
		cin >> row;
		cout << "Enter number of columns in Matrix: ";
		cin >> col;
		cout << endl;
	}
	else
	{
		cout << "Enter size of Matrix: ";
		cin >> row;
		cout << endl;

		col = row;
	}

	// enter matrix
	Matrix<complex<Fraction<int>>> mat(row, col);
	for (size_t i = 0; i != row; ++i)
	{
		cout << "ENTER ROW " << i << ':' << endl;
		for (size_t j = 0; j != col; ++j)
		{
			cout << "Enter real part in column " << j << ": ";
			Fraction<int> real; cin >> real;
			cout << "Enter imaginary part in column " << j << ": ";
			Fraction<int> imag; cin >> imag;
			mat[i][j] = complex<Fraction<int>>(real, imag);
		}
		cout << endl;
	}

	return mat;
}

////////////////////////////////////////
// get input from square matrix with decimal complex specialization
template <>
Matrix<complex<double>> input<complex<double>>(bool isSquare)
{
	// get size of matrix
	size_t row;
	size_t col;
	if (!isSquare)
	{
		cout << "Enter number of rows in Matrix: ";
		cin >> row;
		cout << "Enter number of columns in Matrix: ";
		cin >> col;
		cout << endl;
	}
	else
	{
		cout << "Enter size of Matrix: ";
		cin >> row;
		cout << endl;

		col = row;
	}
	// enter matrix
	Matrix<complex<double>> mat(row, col);
	for (size_t i = 0; i != row; ++i)
	{
		cout << "ENTER ROW " << i << ':' << endl;
		for (size_t j = 0; j != col; ++j)
		{
			cout << "Enter real part in column " << j << ": ";
			double real; cin >> real;
			cout << "Enter imaginary part in column " << j << ": ";
			double imag; cin >> imag;
			mat[i][j] = complex<double>(real, imag);
		}
		cout << endl;
	}

	return mat;
}

////////////////////////////////////////
// print out computation time
void outputComputationTime(steady_clock::time_point t1, steady_clock::time_point t2)
{
	cout << "Computation took "
		<< duration_cast<std::chrono::microseconds>(t2 - t1).count()
		<< " microseconds" << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// LINEAR EQUATION CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// linear equation solver UI
template <typename T>
void solver()
{
	// get number of variables
	cout << "Enter number of variables to solve: ";
	size_t numOfVariables; cin >> numOfVariables;
	cout << endl;

	// enter equations
	vector<vector<T>> eqs;
	vector<vector<T>> ans;
	for (size_t i = 0; i != numOfVariables; ++i)
	{
		cout << "ENTER EQUATION " << i << ':' << endl;
		vector<T> tmpEq;
		vector<T> tmpAns;
		char varName = char(123 - numOfVariables);
		for (size_t j = 0; j != numOfVariables + 1; ++j)
		{
			T coef;
			if (j != numOfVariables)
			{
				cout << "Enter coefficient for ";
				cout << varName << ": ";
				varName += 1;

				cin >> coef;
				tmpEq.push_back(T(coef));
			}
			else
			{
				cout << "Enter answer to equation " << i << ": ";
				cin >> coef;
				tmpAns.push_back(T(coef));
			}
		}

		cout << endl;

		eqs.push_back(tmpEq);
		ans.push_back(tmpAns);
	}

	// find answers
	Matrix<T> eqsMatrix(eqs);
	Matrix<T> ansMatrix(ans);

	auto t1 = Clock::now();
	vector<T> finalAns = solve<T>(eqsMatrix, ansMatrix);
	auto t2 = Clock::now();

	// print answer
	char name = char(123 - numOfVariables);
	for (size_t i = 0; i != finalAns.size(); ++i)
	{
		cout << name << " = " << finalAns[i] << endl;
		name += 1;
	}
	cout << endl;

	cout << "Computation took "
		<< duration_cast<std::chrono::microseconds>(t2 - t1).count()
		<< " microseconds" << endl;
}

////////////////////////////////////////
// solver template specialization for complex numbers with rational coef
template <>
void solver<complex<Fraction<int>>>()
{
	typedef Fraction<int> T;

	// get number of variables
	cout << "Enter number of variables to solve: ";
	size_t numOfVariables; cin >> numOfVariables;
	cout << endl;

	// enter equations
	vector<vector<complex<T>>> eqs;
	vector<vector<complex<T>>> ans;
	for (size_t i = 0; i != numOfVariables; ++i)
	{
		cout << "ENTER EQUATION " << i << ':' << endl;
		vector<complex<T>> tmpEq;
		vector<complex<T>> tmpAns;
		char varName = char(123 - numOfVariables);
		for (size_t j = 0; j != numOfVariables + 1; ++j)
		{
			T realCoef;
			T imagCoef;
			if (j != numOfVariables)
			{
				cout << "Enter coeficient real part for ";
				cout << varName << ": ";
				cin >> realCoef;

				cout << "Enter coeficient imaginary part for ";
				cout << varName << ": ";
				cin >> imagCoef;

				varName += 1;
				tmpEq.push_back(complex<T>(realCoef, imagCoef));
			}
			else
			{
				cout << "Enter answer real part to equation " << i << ": ";
				cin >> realCoef;
				cout << "Enter answer imaginary part to equation " << i << ": ";
				cin >> imagCoef;
				tmpAns.push_back(complex<T>(realCoef, imagCoef));
			}
		}

		cout << endl;

		eqs.push_back(tmpEq);
		ans.push_back(tmpAns);
	}

	// find answers
	Matrix<complex<T>> eqsMatrix(eqs);
	Matrix<complex<T>> ansMatrix(ans);

	auto t1 = Clock::now();
	vector<complex<T>> finalAns = solve<complex<T>>(eqsMatrix, ansMatrix);
	auto t2 = Clock::now();

	// print answer
	char name = char(123 - numOfVariables);
	for (size_t i = 0; i != finalAns.size(); ++i)
	{
		cout << name << " = "
			<< finalAns[i].real() << " + " << finalAns[i].imag() << 'i'
			<< endl;
		name += 1;
	}
	cout << endl;

	cout << "Computation took "
		<< duration_cast<std::chrono::microseconds>(t2 - t1).count()
		<< " microseconds" << endl;
}

////////////////////////////////////////
// solver template specialization for complex numbers with decimal coef
template <>
void solver<complex<double>>()
{
	typedef double T;

	// get number of variables
	cout << "Enter number of variables to solve: ";
	size_t numOfVariables; cin >> numOfVariables;
	cout << endl;

	// enter equations
	vector<vector<complex<T>>> eqs;
	vector<vector<complex<T>>> ans;
	for (size_t i = 0; i != numOfVariables; ++i)
	{
		cout << "ENTER EQUATION " << i << ':' << endl;
		vector<complex<T>> tmpEq;
		vector<complex<T>> tmpAns;
		char varName = char(123 - numOfVariables);
		for (size_t j = 0; j != numOfVariables + 1; ++j)
		{
			T realCoef;
			T imagCoef;
			if (j != numOfVariables)
			{
				cout << "Enter coeficient real part for ";
				cout << varName << ": ";
				cin >> realCoef;

				cout << "Enter coeficient imaginary part for ";
				cout << varName << ": ";
				cin >> imagCoef;

				varName += 1;
				tmpEq.push_back(complex<T>(realCoef, imagCoef));
			}
			else
			{
				cout << "Enter answer real part to equation " << i << ": ";
				cin >> realCoef;
				cout << "Enter answer imaginary part to equation " << i << ": ";
				cin >> imagCoef;
				tmpAns.push_back(complex<T>(realCoef, imagCoef));
			}
		}

		cout << endl;

		eqs.push_back(tmpEq);
		ans.push_back(tmpAns);
	}

	// find answers
	Matrix<complex<T>> eqsMatrix(eqs);
	Matrix<complex<T>> ansMatrix(ans);

	auto t1 = Clock::now();
	vector<complex<T>> finalAns = solve<complex<T>>(eqsMatrix, ansMatrix);
	auto t2 = Clock::now();

	// print answer
	char name = char(123 - numOfVariables);
	for (size_t i = 0; i != finalAns.size(); ++i)
	{
		cout << name << " = "
			<< finalAns[i].real() << " + " << finalAns[i].imag() << 'i'
			<< endl;
		name += 1;
	}
	cout << endl;

	cout << "Computation took "
		<< duration_cast<std::chrono::microseconds>(t2 - t1).count()
		<< " microseconds" << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// RANK CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// rank solver UI
template <typename T>
void solveRank()
{
	Matrix<T> mat = input<T>(false);

	// find rank
	auto t1 = Clock::now();
	size_t rank = mat.rank();
	auto t2 = Clock::now();

	// print answer
	cout << "Rank: " << rank << endl;

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// rank solver UI with rational complex
template <>
void solveRank<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> mat = input<T>(false);

	// find rank
	auto t1 = Clock::now();
	size_t rank = mat.rank();
	auto t2 = Clock::now();

	// print answer
	cout << "Rank: " << rank << endl;

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// rank solver UI with decimal complex numbers
template <>
void solveRank<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> mat = input<T>(false);

	// find rank
	auto t1 = Clock::now();
	size_t rank = mat.rank();
	auto t2 = Clock::now();

	// print answer
	cout << "Rank: " << rank << endl;

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// DETERMINANT CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// determinant solver UI
template <typename T>
void solveDeterminant()
{
	Matrix<T> mat = input<T>();

	// find determinant
	auto t1 = Clock::now();
	T determinant = mat.determinant();
	auto t2 = Clock::now();

	// print answer
	cout << "Determinant: " << determinant << endl;

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// determinant solver UI with rational complex specialization
template <>
void solveDeterminant<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> mat = input<T>();

	// find determinant
	auto t1 = Clock::now();
	T determinant = mat.determinant();
	auto t2 = Clock::now();

	// print answer
	cout << "Determinant, output is in the form (real, imaginary): " << determinant << endl;

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// determinant solver UI decimal complex specialization
template <>
void solveDeterminant<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> mat = input<T>();

	// find determinant
	auto t1 = Clock::now();
	complex<double> determinant = mat.determinant();
	auto t2 = Clock::now();

	// print answer
	cout << "Determinant, output is in the form (real, imaginary): " 
		<< determinant.real() << " + " << determinant.imag() << 'i' 
		<< endl;

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// INVERSE CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// inverse solver UI
template <typename T>
void solveInverse()
{
	Matrix<T> mat = input<T>();

	// make sure inverse exists
	if (mat.determinant() == 0)
	{
		cout << "Inverse doesn't exist" << endl;
		return;
	}

	// find inverse
	auto t1 = Clock::now();
	Matrix<T> result = mat.inverse();
	auto t2 = Clock::now();

	// print answer
	cout << "Inverse: " <<  endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// inverse solver UI with complex rational number specialization
template <>
void solveInverse<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> mat = input<T>();

	// make sure inverse exists
	if (mat.determinant() == T(0))
	{
		cout << "Inverse doesn't exist" << endl;
		return;
	}

	// find inverse
	auto t1 = Clock::now();
	Matrix<T> result = mat.inverse();
	auto t2 = Clock::now();

	// print answer
	cout << "Inverse, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// inverse solver UI with complex rational number specialization
template <>
void solveInverse<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> mat = input<T>();

	// make sure inverse exists
	if (mat.determinant() == T(0))
	{
		cout << "Inverse doesn't exist" << endl;
		return;
	}

	// find inverse
	auto t1 = Clock::now();
	Matrix<T> result = mat.inverse();
	auto t2 = Clock::now();

	// print answer
	cout << "Inverse, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX POWER CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// matrix power solver UI
template <typename T>
void solvePower()
{
	Matrix<T> mat = input<T>();

	cout << "Enter power (whole numbers only): ";
	size_t power; cin >> power;

	// find power
	auto t1 = Clock::now();
	Matrix<T> result = mat ^ power;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer: " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix power solver UI with rational complex number specialization
template <>
void solvePower<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> mat = input<T>();

	cout << "Enter power (whole numbers only): ";
	size_t power; cin >> power;

	// find power
	auto t1 = Clock::now();
	Matrix<T> result = mat ^ power;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix power solver UI with real complex number specialization
template <>
void solvePower<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> mat = input<T>();

	cout << "Enter power (whole numbers only): ";
	size_t power; cin >> power;

	// find power
	auto t1 = Clock::now();
	Matrix<T> result = mat ^ power;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX TRANSPOSE CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// matrix transpose solver UI
template <typename T>
void solveTranspose()
{
	Matrix<T> mat = input<T>(false);

	// find transpose
	auto t1 = Clock::now();
	mat.transpose();
	auto t2 = Clock::now();

	// print answer
	cout << "Answer: " << endl;
	mat.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix transpose solver UI with rational complex specialization
template <>
void solveTranspose<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> mat = input<T>(false);

	// find transpose
	auto t1 = Clock::now();
	mat.transpose();
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	mat.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix transpose solver UI with decimal complex specialization
template <>
void solveTranspose<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> mat = input<T>(false);

	// find transpose
	auto t1 = Clock::now();
	mat.transpose();
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	mat.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX MULTIPLICATION CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// matrix multiplication solver UI
template <typename T>
void solveMultiply()
{
	Matrix<T> one = input<T>(false);
	Matrix<T> two = input<T>(false);

	// check if they can be multiplied
	if (!checkDimMultiply(one, two))
	{
		cout << "Cannot multiply these matrices together" << endl;
		return;
	}

	// find answer
	auto t1 = Clock::now();
	Matrix<T> result = one * two;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer: " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix multiplication solver UI rational complex specialization
template <>
void solveMultiply<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> one = input<T>(false);
	Matrix<T> two = input<T>(false);

	// check if they can be multiplied
	if (!checkDimMultiply(one, two))
	{
		cout << "Cannot multiply these matrices together" << endl;
		return;
	}

	// find answer
	auto t1 = Clock::now();
	Matrix<T> result = one * two;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix multiplication solver UI decimal complex specialization
template <>
void solveMultiply<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> one = input<T>(false);
	Matrix<T> two = input<T>(false);

	// check if they can be multiplied
	if (!checkDimMultiply(one, two))
	{
		cout << "Cannot multiply these matrices together" << endl;
		return;
	}

	// find answer
	auto t1 = Clock::now();
	Matrix<T> result = one * two;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// SCALAR MULTIPLICATION CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// scalar multiplication solver UI 
template <typename T>
void solveScalar()
{
	Matrix<T> one = input<T>(false);
	
	cout << "Enter scalar: ";
	T scalar; cin >> scalar;

	// find answer
	auto t1 = Clock::now();
	Matrix<T> result = one * scalar;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer: " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// scalar multiplication solver UI rational complex specialization
template <>
void solveScalar<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> one = input<T>(false);

	cout << "Enter real part of scalar: ";
	Fraction<int> scalarReal; cin >> scalarReal;
	cout << "Enter imaginary part of scalar: ";
	Fraction<int> scalarImag; cin >> scalarImag;
	T scalar(scalarReal, scalarImag);

	// find answer
	auto t1 = Clock::now();
	Matrix<T> result = one * scalar;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// scalar multiplication solver UI decimal complex specialization
template <>
void solveScalar<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> one = input<T>(false);

	cout << "Enter real part of scalar: ";
	double scalarReal; cin >> scalarReal;
	cout << "Enter imaginary part of scalar: ";
	double scalarImag; cin >> scalarImag;
	T scalar(scalarReal, scalarImag);

	// find answer
	auto t1 = Clock::now();
	Matrix<T> result = one * scalar;
	auto t2 = Clock::now();

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// ROW ECHELON CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// row echelon solver UI
template <typename T>
void solveRowEchelon()
{
	Matrix<T> mat = input<T>(false);

	// solve
	auto t1 = Clock::now();
	Matrix<T> result = mat.rowEchelon();
	auto t2 = Clock::now();

	// print answer
	cout << "Row Echelon Form: " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// row echelon solver UI rational complex specialization
template <>
void solveRowEchelon<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> mat = input<T>(false);

	// solve
	auto t1 = Clock::now();
	Matrix<T> result = mat.rowEchelon();
	auto t2 = Clock::now();

	// print answer
	cout << "Row Echelon Form, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// row echelon solver UI rational complex specialization
template <>
void solveRowEchelon<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> mat = input<T>(false);

	// solve
	auto t1 = Clock::now();
	Matrix<T> result = mat.rowEchelon();
	auto t2 = Clock::now();

	// print answer
	cout << "Row Echelon Form, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// MATRIX ADD/SUBTRACT CALCULATOR UI FUNCTIONS
////////////////////////////////////////
// matrix addition or subtraction solver UI
template <typename T>
void solveAdd()
{
	Matrix<T> one = input<T>(false);
	Matrix<T> two = input<T>(false);

	// check if they can be added
	if (!checkDimSame(one, two))
	{
		cout << "Cannot add or subtract these matrices together" << endl;
		return;
	}

	cout << "Add or Subtract? (a,s): ";
	char option; cin >> option;

	// calculate
	Matrix<T> result;
	auto t1 = Clock::now();
	auto t2 = Clock::now();
	if (option == 's')
	{
		t1 = Clock::now();
		result = one - two;
		t2 = Clock::now();
	}
	else
	{
		t1 = Clock::now();
		result = one + two;
		t2 = Clock::now();
	}

	// print answer
	cout << "Answer: " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix addition or subtraction solver UI with rational complex specialization
template <>
void solveAdd<complex<Fraction<int>>>()
{
	typedef complex<Fraction<int>> T;

	Matrix<T> one = input<T>(false);
	Matrix<T> two = input<T>(false);

	// check if they can be added
	if (!checkDimSame(one, two))
	{
		cout << "Cannot add or subtract these matrices together" << endl;
		return;
	}

	cout << "Add or Subtract? (a,s): ";
	char option; cin >> option;

	// calculate
	Matrix<T> result;
	auto t1 = Clock::now();
	auto t2 = Clock::now();
	if (option == 's')
	{
		t1 = Clock::now();
		result = one - two;
		t2 = Clock::now();
	}
	else
	{
		t1 = Clock::now();
		result = one + two;
		t2 = Clock::now();
	}

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////
// matrix addition or subtraction solver UI with decimal complex specialization
template <>
void solveAdd<complex<double>>()
{
	typedef complex<double> T;

	Matrix<T> one = input<T>(false);
	Matrix<T> two = input<T>(false);

	// check if they can be added
	if (!checkDimSame(one, two))
	{
		cout << "Cannot add or subtract these matrices together" << endl;
		return;
	}

	cout << "Add or Subtract? (a,s): ";
	char option; cin >> option;

	// calculate
	Matrix<T> result;
	auto t1 = Clock::now();
	auto t2 = Clock::now();
	if (option == 's')
	{
		t1 = Clock::now();
		result = one - two;
		t2 = Clock::now();
	}
	else
	{
		t1 = Clock::now();
		result = one + two;
		t2 = Clock::now();
	}

	// print answer
	cout << "Answer, output is in the form (real, imaginary): " << endl;
	result.print();

	outputComputationTime(t1, t2);
}

////////////////////////////////////////////////////////////////////////////////
//
// HELPER CALCULATOR SELECTOR
////////////////////////////////////////
// runs selected calculator
template <typename T>
void calcSelector(size_t calcSelection)
{
	if      (calcSelection == 0) solver<T>();
	else if (calcSelection == 1) solveRank<T>();
	else if (calcSelection == 2) solveDeterminant<T>();
	else if (calcSelection == 3) solveInverse<T>();
	else if (calcSelection == 4) solvePower<T>();
	else if (calcSelection == 5) solveTranspose<T>();
	else if (calcSelection == 6) solveMultiply<T>();
	else if (calcSelection == 7) solveScalar<T>();
	else if (calcSelection == 8) solveRowEchelon<T>();
	else                         solveAdd<T>();
}
#endif // SOLVER_H