#ifndef FRACTION_H
#define FRACTION_H

////////////////////////////////////////////////////////////////////////////////
//
// FILE:        Fraction.h
// DESCRIPTION: contains template fraction class and implementation, also GCD function
// AUTHOR:      Dan Fabian
// DATE:        5/21/2019

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstddef> // used for size_t
#include <type_traits>

using std::ostream; using std::istream;
using std::is_same;

////////////////////////////////////////////////////////////////////////////////
//
// FRACTION
template <typename T = int> // T must be some type of number, DOUBLE'S NOT ALLOWED
class Fraction {
public:
	// prevents doubles and floats
	static_assert(!is_same<T, double>::value && !is_same<T, float>::value,
		"No doubles or floats allowed");

	// constructors
	Fraction() : numer_(T(0)), denom_(T(1)) {}
	Fraction(T num, T den);                                        // initalize values, makes sure denom != 0
	Fraction(T num) : numer_(num), denom_(T(1)) {}                 // initalize to a whole value
	Fraction(const Fraction<T>& orig);                             // copy constructor

	// methods
	void   simplify();                                             // simplify the fraction, uses Euclidean Alg
	double toDecimal() const { return double(numer_) / double(denom_); }
	void   inverse();                                              // flip fraction

	// comparisons
	bool operator== (Fraction<T> rhs) const;
	bool operator<  (Fraction<T> rhs) const;
	bool operator!= (Fraction<T> rhs) const { return !(*this == rhs); }
	bool operator<= (Fraction<T> rhs) const { return *this == rhs || *this < rhs; }
	bool operator>  (Fraction<T> rhs) const { return rhs < *this; }
	bool operator>= (Fraction<T> rhs) const { return *this == rhs || *this > rhs; }

	// binary operations
	Fraction<T>& operator= (const Fraction<T> & rhs);       // overloaded assignment
	Fraction<T>  operator+ (Fraction<T> rhs)         const; // addition
	Fraction<T>  operator* (const Fraction<T> & rhs) const; // multiplication
	Fraction<T>  operator- (const Fraction<T> & rhs) const; // subtraction (only works if T is allowed to be negative)
	Fraction<T>  operator/ (Fraction<T> rhs)         const; // division
	Fraction<T>  operator^ (const Fraction<T> & rhs) const; // powers

	// unary operators
	Fraction<T>& operator+= (const Fraction & rhs)      { *this = *this + rhs; return *this; }
	Fraction<T>& operator*= (const Fraction & rhs)      { *this = *this * rhs; return *this; }
	Fraction<T>& operator-= (const Fraction & rhs)      { *this = *this - rhs; return *this; }
	Fraction<T>& operator/= (const Fraction & rhs)      { *this = *this / rhs; return *this; }
	Fraction<T>& operator^= (const Fraction & rhs)      { *this = *this ^ rhs; return *this; }
	Fraction<T>  operator-  ()                    const { return *this* T(-1); }
	Fraction<T>  operator+  ()                    const { return *this; }

	// conversion operators
	explicit operator double() const { return this->toDecimal(); }

	// friend operators
	template <typename K>
	friend Fraction<K> operator+ (const K & lhs, const Fraction<K> & rhs);
	template <typename K>
	friend Fraction<K> operator* (const K & lhs, const Fraction<K> & rhs);
	template <typename K>
	friend Fraction<K> operator- (const K & lhs, const Fraction<K> & rhs);
	template <typename K>
	friend Fraction<K> operator/ (const K & lhs, const Fraction<K> & rhs);
	template <typename K>
	friend Fraction<K> operator^ (const K & lhs, const Fraction<K> & rhs);

	// input and output
	template <typename K>
	friend ostream& operator<< (ostream & out, const Fraction<K> & val);
	template <typename K>
	friend istream& operator>> (istream & in, Fraction<K> & val);

private:
	T numer_;
	T denom_;
};

////////////////////////////////////////////////////////////////////////////////
//
// GCD function
template <typename T>
T GCD(const T & a, const T & b)
{
	// use absolute values
	T aAbs = abs(a);
	T bAbs = abs(b);

	if (aAbs < bAbs)
		return aAbs == 0 ? bAbs : GCD<T>(aAbs, bAbs % aAbs);

	return bAbs == 0 ? aAbs : GCD<T>(bAbs, aAbs % bAbs);
}

////////////////////////////////////////////////////////////////////////////////
//
// FRACTION constructors
////////////////////////////////////////
// constructs fraction checking for a non zero denom
template <typename T>
Fraction<T>::Fraction(T num, T den)
{
	assert(den != T(0));
	numer_ = num;
	denom_ = den;
}

////////////////////////////////////////
// copy constructor
template <typename T>
Fraction<T>::Fraction(const Fraction<T> & orig)
	: numer_(orig.numer_), denom_(orig.denom_) {}

////////////////////////////////////////////////////////////////////////////////
//
// FRACTION methods
////////////////////////////////////////
// finds GCD then simplifies
template <typename T>
void Fraction<T>::simplify()
{
	T gcd = GCD<T>(numer_, denom_);

	// simplify
	numer_ /= gcd;
	denom_ /= gcd;

	// check if both are negative, or if denom is negative move it to numer
	if (numer_ < T(0) && denom_ < T(0) || denom_ < T(0))
	{
		numer_ *= T(-1);
		denom_ *= T(-1);
	}

	if (numer_ == T(0)) denom_ = T(1);
}

////////////////////////////////////////
// flips fraction
template <typename T>
void Fraction<T>::inverse()
{
	if (numer_ == T(0))
	{
		denom_ = T(1);
		return;
	}

	T tmp = numer_;
	numer_ = denom_;
	denom_ = tmp;

	this->simplify();
}

////////////////////////////////////////////////////////////////////////////////
//
// FRACTION comparisons
////////////////////////////////////////
// checks for equality
template <typename T>
bool Fraction<T>::operator==(Fraction<T> rhs) const
{
	Fraction<T> thisCopy(*this);
	thisCopy.simplify();
	rhs.simplify();

	return thisCopy.denom_ == rhs.denom_ && thisCopy.numer_ == rhs.numer_;
}

////////////////////////////////////////
// checks for less than
template <typename T>
bool Fraction<T>::operator<(Fraction<T> rhs) const
{
	Fraction<T> thisCopy(*this);
	thisCopy.simplify();
	rhs.simplify();

	return (thisCopy.numer_ * rhs.denom_) < (rhs.numer_ * thisCopy.denom_);
}

////////////////////////////////////////////////////////////////////////////////
//
// FRACTION operations
////////////////////////////////////////
// overloaded assignment
template <typename T>
Fraction<T>& Fraction<T>::operator=(const Fraction<T> & rhs)
{
	numer_ = rhs.numer_;
	denom_ = rhs.denom_;

	return *this;
}

////////////////////////////////////////
// fraction addition
template <typename T>
Fraction<T> Fraction<T>::operator+(Fraction<T> rhs) const
{
	Fraction<T> thisCopy(*this);
	thisCopy.simplify();
	rhs.simplify();

	Fraction<T> result((thisCopy.numer_ * rhs.denom_) + (rhs.numer_ * thisCopy.denom_),
		thisCopy.denom_ * rhs.denom_);

	result.simplify();
	return result;
}

////////////////////////////////////////
// fraction multiplication
template <typename T>
Fraction<T> Fraction<T>::operator*(const Fraction<T> & rhs) const
{
	Fraction<T> result(numer_ * rhs.numer_, denom_ * rhs.denom_);

	result.simplify();
	return result;
}

////////////////////////////////////////
// fraction subtraction
template <typename T>
Fraction<T> Fraction<T>::operator-(const Fraction<T> & rhs) const
{
	Fraction<T> result = *this + (rhs * T(-1));

	result.simplify();
	return result;
}

////////////////////////////////////////
// fraction division
template <typename T>
Fraction<T> Fraction<T>::operator/(Fraction<T> rhs) const
{
	assert(rhs != 0);

	rhs.inverse();
	Fraction<T> result = *this * rhs;

	result.simplify();
	return result;
}

////////////////////////////////////////
// fraction division
template <typename T>
Fraction<T> Fraction<T>::operator^(const Fraction<T> & rhs) const
{
	if (rhs == 0) return Fraction<T>(T(1));

	Fraction<T> result(static_cast<T>(pow(numer_, rhs.toDecimal())),
		static_cast<T>(pow(denom_, rhs.toDecimal())));

	result.simplify();
	return result;
}

////////////////////////////////////////////////////////////////////////////////
//
// FRACTION friend operations
////////////////////////////////////////
// addition with non-fraction lhs
template <typename T>
Fraction<T> operator+(const T & lhs, const Fraction<T> & rhs)
{
	return rhs + lhs;
}

////////////////////////////////////////
// multiplication with non-fraction lhs
template <typename T>
Fraction<T> operator*(const T & lhs, const Fraction<T> & rhs)
{
	return rhs * lhs;
}

////////////////////////////////////////
// subtraction with non-fraction lhs
template <typename T>
Fraction<T> operator-(const T & lhs, const Fraction<T> & rhs)
{
	return (rhs - lhs) * T(-1);
}

////////////////////////////////////////
// division with non-fraction lhs
template <typename T>
Fraction<T> operator/(const T & lhs, const Fraction<T> & rhs)
{
	return Fraction<T>(lhs) / rhs;
}

////////////////////////////////////////
// power with non-fraction lhs
template <typename T>
Fraction<T> operator^(const T & lhs, const Fraction<T> & rhs)
{
	return Fraction<T>(lhs) ^ rhs;
}

////////////////////////////////////////////////////////////////////////////////
//
// FRACTION output and input
////////////////////////////////////////
// outputs fraction
template <typename T>
ostream& operator<<(ostream & out, const Fraction<T> & val)
{
	if (val.denom_ != 1)
		out << '(' << val.numer_ << " / " << val.denom_ << ')';
	else
		out << val.numer_;

	return out;
}

// input fraction
template <typename T>
istream& operator>>(istream & in, Fraction<T> & val)
{
	T numer; in >> numer;

	val = Fraction<T>(numer);

	return in;
}

#endif // FRACTION_H

