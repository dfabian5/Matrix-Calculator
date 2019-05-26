////////////////////////////////////////////////////////////////////////////////
//
// FILE:        main.cpp
// DESCRIPTION: runs all tests and solves user entered systems
// AUTHOR:      Dan Fabian 
// DATE:        5/25/2019

#include "Calculator.h" // contains all needed headers

int main()
{
	cout << "ULTIMATE MATRIX CALCULATOR" << endl << endl;

	// main loop
	char another; // used to ask to enter another system
	do
	{
		cout << "Which field would you like to use? " << endl
			<< "0. Rationals" << endl
			<< "1. Reals" << endl
			<< "2. Complex with rational coeficients" << endl
			<< "3. Complex with decimal coeficients" << endl
			<< endl << "Selection: ";
		int fieldSelection; cin >> fieldSelection;
		cout << endl;

		cout << "Which calculator would you like to use? " << endl
			<< "0. Linear Equation Solver" << endl
			<< "1. Matrix Rank Solver" << endl
			<< "2. Matrix Determinant Solver" << endl
			<< "3. Matrix Inverse Solver" << endl
			<< "4. Matrix Power Solver" << endl
			<< "5. Matrix Transpose Solver" << endl
			<< "6. Matrix Multiplication Solver" << endl
			<< "7. Scalar Multiplication Solver" << endl
			<< "8. Matrix Addition/Subtracton Solver" << endl
			<< endl << "Selection: ";
		int calcSelection; cin >> calcSelection;

		// run calc
		if      (fieldSelection == 0) calcSelector<Fraction<int>>(calcSelection);
		else if (fieldSelection == 1) calcSelector<double>(calcSelection);
		else if (fieldSelection == 2) calcSelector<complex<Fraction<int>>>(calcSelection);
		else                          calcSelector<complex<double>>(calcSelection);

		cout << endl << "Enter another problem? (y/n): ";
		cin >> another;
		cout << endl;

	} while (another != 'n');
}
