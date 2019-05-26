# Matrix-Calculator Made By Dan Fabian
This project contains full classes for both Fractions and Matricies to use in the calculator.
The matrix class has optimizations such as multi threading for multiplication, a special Row
Echelon function for square matrices that heavily decreases the run time of finding inverses
and determinants.

## File Contents
  - "Calculator.h" contains all the calculator UI and linear equation solving functions
  - "Matrix.h" contains templated Matrix class and implementation
  - "Fraction" contains templated Fraction class and implementation
  - "main.cpp" contains main function...

## Calculator Can Work With:
  - Rational numbers
  - Real numbers
  - Complex numbers with rational coeficients
  - Complex numbers with real coeficients

## Calculator Can:
  - calculate answers to a system of linear equations
  - calculate a matrix's rank
  - calculate a matrix's determinant
  - calculate a matrix's inverse
  - calculate a matrix to a power
  - calculate a matrix's transpose
  - calculate a matrix's Row Echelon Form
  - multiply two matrices
  - scalar multiplication
  - add/subtract two matrices

## Notes
  - Matrix and Fraction classes are both templated for any type of number
  - There is no max size of matricies, the only limitation is your system architecture
