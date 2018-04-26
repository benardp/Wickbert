/** @file Interval.h
*  @brief Interval class definition.
*  @author Matei Stroila, Bart Stander, John C. Hart
*
* Adapted by Hart from Stander's original 1997 ism code.
* This class Interval is a wrapper for the boost library's interval class.
* @todo Remove global functions like IExp().
*/


#ifndef INTERVAL_H
#define INTERVAL_H

#include <math.h>
#include "Infinite.h"
#include <boost/numeric/interval.hpp>
#include <boost/assign.hpp>

/** Tristate logic for interval conditions.
* \li I_TRUE - True for all possible members.
* \li I_FALSE - False for all possible members.
* \li I_MAYBE - True for some, false for others.
*/
enum I_TRISTATE { I_TRUE = 0, I_FALSE = -1, I_MAYBE = 1 };

/** Used to indicate result of division by zero
* \li INORMAL - Division resulted in a single interval - normally
* \li ISPLIT  - Division resulted in two separate intervals
* \li IDIV_BY_ZERO - Error since denominator was zero
*/
enum IDIF_TYPES { INORMAL = 0, ISPLIT, IDIV_BY_ZERO };

/** Complement of an I_TRISTATE variable.
* \li I_TRUE->I_FALSE
* \li I_FALSE->I_TRUE
* \li I_MAYBE->I_MAYBE
*/


using namespace boost::numeric;
using namespace interval_lib;
using namespace compare::certain;
template <class Type> struct my_arith_rounding: rounded_arith_opp<Type>
{
	double int_down(Type f) { return floor(f); }
	double int_up  (Type f) { return ceil(f); }
};

//forward declarations needed for the operator <<
template <class Type> class Interval;
template <class Type> std::ostream & operator<< (std::ostream& os, const Interval<Type>& m);	

/** Uses checking_no_nan (instead of checking_strict which doesn't allow
 *  zero-width intervals.
 */
template <class Type> class Interval : public interval<Type,policies<save_state<
	rounded_transc_exact<Type,my_arith_rounding<Type> > >,
	checking_no_nan<Type> > >
{

	typedef interval<Type,policies<save_state<
	rounded_transc_exact<Type,my_arith_rounding<Type> > >,
	checking_no_nan<Type> > > BoostInterval;
public:
	
	/// @name Constructors @{
	Interval(void){	};
    Interval(Type const &a) : BoostInterval(a)
	{	};
    Interval(Type const &l, Type const &h) : BoostInterval(l,h)
	{	};
	Interval(BoostInterval const &i) : BoostInterval(i)
	{	};
	Interval(Interval<Type> const &x) : BoostInterval()
	{
		this->assign(x.lower(),x.upper());
	};
    /// @}

	virtual ~Interval(void){	};
	
    /// @name Attribute Access and Assignment @{
	inline Type low(void) const{
		return this->lower();
	} 
    inline Type high(void) const{
		return this->upper();
	}
	
    operator Type() const { return (Type)(median(*this)); } ///< Type conversion
    
	inline void setInterval(Type const &low, Type const &high){
		this->assign(low,high);
	};

	Interval<Type>& operator=(const Interval<Type>& x){
		this->assign(x.lower(),x.upper());
		return *this;
	};


	/// @}
	
    /// @name Checks @{
    virtual void print() const;
    virtual void error() const;
	virtual void check_interval() const;
    bool thin() const;
    /// @}
	
    /// @name Unary Operators @{
    inline Type width() const{
#if 0
		return boost::numeric::width(*this);
#else
		return high() - low();
#endif
	}
	
	inline Type center() const{
		return median(*this);
	}
	
    Interval<Type> exp() const;
	
	inline Interval<Type> squared() const
	{
		return Interval<Type>(boost::numeric::square(*this));
	}
	
    Interval<Type> cubed() const;
    
	inline Interval<Type> log() const
	{
		return Interval<Type>(boost::numeric::log(*this));
	}
    
	Interval<Type> pow(Interval<Type> n) const;
	
	inline Interval<Type> pow(int n) const
	{
		return Interval<Type>(boost::numeric::pow(*this,n));
	}
	
	
	Interval<Type> sqrt() const;
    /// @}
	
    /// @name Binary Operators @{
	Interval<Type> intersection(Interval<Type> &b) const;
	Interval<Type> unionWith(Interval<Type> &b) const;
    /// @}
	
    /// @name Conditions @{
	
    bool contains(const Interval<Type> &a) const;
    bool containsProper(const Interval<Type> &a) const;
    bool overlaps(const Interval<Type> &a) const;
    bool disjoint(const Interval<Type> &a) const;
    bool isZero() const;
    bool isPositive() const;
    bool isNegative() const;
    bool isNonnegative() const;
    bool isNonpositive() const;
    bool isMixed() const;
    /// @}

/// Basic output
    friend std::ostream & operator<<<Type> (std::ostream& os, const Interval<Type>& m);	
		
};


/* arithmetic operators involving intervals */
template <class Type>
inline Interval<Type> operator+(const Interval<Type>& x){
	return Interval<Type>(boost::numeric::operator+(x));
};
template <class Type>
inline Interval<Type> operator-(const Interval<Type>& x){
	return Interval<Type>(boost::numeric::operator-(x));
};
template <class Type>
inline Interval<Type> operator+(const Interval<Type>& x, const Interval<Type>& y){
	return Interval<Type>(boost::numeric::operator+(x,y));
};
template <class Type>
inline Interval<Type> operator+(const Interval<Type>& x, const Type& y){
	return Interval<Type>(boost::numeric::operator+(x,y));
};
template <class Type>
inline Interval<Type> operator+(const Type& x, const Interval<Type>& y){
	return Interval<Type>(boost::numeric::operator+(x,y));
};

template <class Type>
inline Interval<Type> operator-(const Interval<Type>& x, const Interval<Type>& y){
	return Interval<Type>(boost::numeric::operator-(x,y));
};
template <class Type>
inline Interval<Type> operator-(const Interval<Type>& x, const Type& y){
	return Interval<Type>(boost::numeric::operator-(x,y));
};
template <class Type>
inline Interval<Type> operator-(const Type& x, const Interval<Type>& y){
	return Interval<Type>(boost::numeric::operator-(x,y));
};

template <class Type>
inline Interval<Type> operator*(const Interval<Type>& x, const Interval<Type>& y){
	return Interval<Type>(boost::numeric::operator*(x,y));
};
template <class Type>
inline Interval<Type> operator*(const Interval<Type>& x, const Type& y){
	return Interval<Type>(boost::numeric::operator*(x,y));
};
template <class Type>
inline Interval<Type> operator*(const Type& x, const Interval<Type>& y){
	return Interval<Type>(boost::numeric::operator*(x,y));
};

template <class Type>
inline Interval<Type> operator/(const Interval<Type>& x, const Interval<Type>& y){
	return Interval<Type>(boost::numeric::operator/(x,y));
};
template <class Type>
inline Interval<Type> operator/(const Interval<Type>& x, const Type& y){
	return Interval<Type>(boost::numeric::operator/(x,y));
};
template <class Type>
inline Interval<Type> operator/(const Type& r, const Interval<Type>& x){
	return Interval<Type>(boost::numeric::operator/(r,x));
};


typedef Interval<double> Intervald;

/** Full Interval division
* In this extended version, divide by zero is allowed. When this happens
* two semi-infinite intervals result.
* \param x Numerator.
* \param y Denominator.
* \param v Quotient.
* \param v2 Second part of quotient if a split occurs.
* \return
*   \li ISPLIT if a split occurs.
*   \li ()IDIV_BY_0, if that occurs.
*   \li INORMAL (nothing bad happened).
*
* IDinf, interval divide by infinity, was taken from John M. Snyder's
*       Thesis, California Institute of Technology, 1991
*/
template<class Type>
IDIF_TYPES IDinf(Interval<Type> &v, Interval<Type> &v2, 
                 const Interval<Type> &x, const Interval<Type> &y)
{
	Interval<Type> NearZero(1./(-INT_MAX),1./INT_MAX);
	
	if (y.isMixed()) 
    {
		if (x.isNegative()) 
        {
			if (NearZero.contains(Interval<Type>(y.high()))) 
            {
				v=Interval<Type>(x.high()/y.low(),Infinite<Type>::I_INF);
				return(IDIV_BY_ZERO);
            } 
			else if (NearZero.contains(Interval<Type>(y.low()))) 
            {
				v=Interval<Type>(Infinite<Type>::I_NEGINF,-x.high()/y.high());
				return(IDIV_BY_ZERO);
            } 
			else 
            {
				v=Interval<Type>(Infinite<Type>::I_NEGINF,x.high()/y.high());
				v2=Interval<Type>(x.high()/y.low(),Infinite<Type>::I_INF);
				return(ISPLIT);
            }
        } 
		else if (x.isPositive()) 
        {
			if (NearZero.contains(Interval<Type>(y.high()))) 
            {
				v = Interval<Type>(Infinite<Type>::I_NEGINF,x.low()/y.low());
				return(IDIV_BY_ZERO);
            } 
			else if (NearZero.contains(Interval<Type>(y.low()))) 
            {
				v = Interval<Type>(x.low()/y.high(),Infinite<Type>::I_INF);
				return(IDIV_BY_ZERO);
            } 
			else 
            {
				v = Interval<Type>(Infinite<Type>::I_NEGINF,x.low()/y.low());
				v2 = Interval<Type>(x.low()/y.high(),Infinite<Type>::I_INF);
				return(ISPLIT);
            }
        } 
		else 
        {
			/* numerator negative and positive */
			v=Interval<Type>(Infinite<Type>::I_NEGINF,
                             Infinite<Type>::I_INF);
			return(IDIV_BY_ZERO);
        }
    } 
	else 
    { /* else a normal divide */
		v = x / y;
		return(INORMAL); /* normal interval division */
    }

#ifdef _DEBUG
v.check_interval();
v2.check_interval();
#endif
}

/// Positive is True only if interval entirely positive.
template <class Type>
bool Interval<Type>::isPositive() const 
{
	return *this > (Type) 0;
}

/// isNegative is True only if interval entirely negative.
template <class Type>
bool Interval<Type>::isNegative() const 
{
	return *this < (Type) 0;
}

/// isNonnegative is True only if interval entirely positive.
template <class Type>
bool Interval<Type>::isNonnegative() const 
{
	return *this >= (Type) 0;
}

/// isNonpositive is True only if interval entirely negative.
template <class Type>
bool Interval<Type>::isNonpositive() const 
{
	return *this <= (Type) 0;
}

/*
template <class Type>
Interval<Type> Interval<Type>::pow(int n) const
{
	return Interval<Type>(boost::numeric::pow(*this,n));
}
*/

/// Return e to the interval power.
template <class Type>
Interval<Type> Interval<Type>::exp() const 
{
	return Interval<Type>(boost::numeric::exp(*this));
}

template <class Type>
Interval<Type> Interval<Type>::sqrt() const
{
	return Interval<Type>(boost::numeric::sqrt(*this));
}

/// Checks if interval is zero width
template <class Type>
bool Interval<Type>::thin() const 
{
	return singleton(*this);
}


/// Checks if the interval is pretty much [0,0]
template <class Type>
bool Interval<Type>::isZero() const 
{
	return (thin() && in_zero(*this));
}

/// Overlaps is True if a intersects me
template <class Type>
bool Interval<Type>::overlaps(const Interval<Type> &a) const 
{
	return overlap(*this,a);
}

/// Raises an interval to the power of another interval
template <class Type>
Interval<Type> Interval<Type>::pow(Interval<Type> n) const
{
	Interval<Type> retval;
	Interval<Type> t1,t2;
	
	if (this->lower() < 0.0)
    {
		std::cerr << "Error: Trying to pow() a partially negative interval: " <<
		this << std::endl;
    }
	else
    {
		bool _lowIsZero = (this->lower() == 0.0);
		if ((_lowIsZero) && (n.lower() <= 0.0))
        {
			std::cerr <<
            "Error: Trying to pow() with a partially negative exponent: " <<
			this << std::endl;
        }
		else
        {
			if (_lowIsZero)
            {
				t2 = Interval<Type>(this->upper());
				t1 = t2.log();
            }
			else
				t1 = (*this).log();
			
			t2 = t1 * n;
			if (_lowIsZero)
            {
				t1 = t2.exp();
				retval.setInterval(std::min(0.0,t1.lower()),
								   std::max(0.0,t1.upper()));
            }
			else
				retval = t2.exp();
        }
    }
	return retval;
}

/// Cubes the interval
template <class Type>
Interval<Type> Interval<Type>::cubed() const 
{
	return Interval<Type>(this->lower()*this->lower()*this->lower(),this->upper()*this->upper()*this->upper());
}

template <class Type>
Interval<Type> Interval<Type>::intersection(Interval<Type> &b) const
{
	return Interval<Type>(boost::numeric::intersect(*this, b));
}

template <class Type>
bool Interval<Type>::contains(const Interval<Type> &a) const
{
	return subset(a, *this); 
}

template <class Type>
bool Interval<Type>::containsProper(const Interval<Type> &a) const
{
	return proper_subset(a, *this); 
}

/// Disjoint is True if a does not intersect me
template <class Type>
bool Interval<Type>::disjoint(const Interval<Type> &a) const
{
	return this->lower() > a.upper() || this->upper() < a.lower();
}

/** True only if interval has positive and negative elements.
* \note [0,0] is not mixed.
*/
template <class Type>
bool Interval<Type>::isMixed() const 
{
	return this->lower() < 0.0 && this->upper() > 0.0;
}

/// Union of intervals = smallest interval containing both intervals.
template <class Type>
Interval<Type> Interval<Type>::unionWith(Interval &b) const
{
	return Interval<Type>(gmMin(this->lower(),b.lower()),gmMax(this->upper(),b.upper()));
}

/// Basic print out of the interval
template <class Type>
void Interval<Type>::print() const
{
	std::cout << this << std::endl;
}


/// Error function can be called when low > high
template <class Type>
void Interval<Type>::error() const 
{
	std::cerr << "Error: [" << this << std::endl;
}


/// Checks if lower < upper.
template <class Type>
void Interval<Type>::check_interval() const 
{
	if (this->lower() > this->upper()) 
		this->error();
}

/// Basic Output
template <class Type> 
std::ostream & operator<<(std::ostream& os, const Interval<Type>& m)
{
  os << "[ " << m.low() << "," << m.high() << "]";
  return os;
}	



#endif
