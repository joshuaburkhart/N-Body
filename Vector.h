/*						
  3-dimensional vectors				
						
  John Conery					
  CIS 455/555					
  University of Oregon				
*/

#ifndef _vector_h_
#define _vector_h_

#include <iostream>
using std::ostream;

// A Vector is a simple 3-tuple of doubles.  This class defines methods
// and operators for vector objects (e.g. adding a scalar to a vector,
// or calculating the norm of a vector), and it also defines several
// friend operators that create new vectors from existing vector objects
// (e.g. computing the sum of two vectors).

class Vector {
public:
  // constructors
  Vector();					// default constructor
  Vector(const Vector& v);			// copy constructor
  Vector(double x, double y, double z); 	// main constructor
  Vector(const double *msg);			// make vector from message

  // operators on existing vector objects
  Vector& operator=(const Vector& v); 		// vector = vector
  Vector& operator=(double a);			// vector = scalar
  Vector& operator+=(const Vector& v); 		// vector += vector
  Vector& operator+=(double a);			// vector += scalar
  Vector& operator-=(const Vector& v); 		// vector -= vector
  Vector& operator-=(double a);			// vector -= scalar
  Vector& operator*=(double a);			// vector += scalar
  double norm();				// Euclidean norm (wrt (0,0,0))
  void makeMessage(double *msg);	 	// makes an MPI message

  // friend operators that create new vector objects
  friend Vector operator+(const Vector &v1, const Vector& v2);
  friend Vector operator+(const Vector &v1, double a);
  friend Vector operator-(const Vector &v1, const Vector& v2);
  friend Vector operator-(const Vector &v1, double a);
  friend Vector operator*(const Vector &v1, double a);

  // Dot product (vector * vector = scalar)
  friend double operator*(const Vector &v1, const Vector& v2);

  // output operator
  friend ostream& operator<<(ostream& s, const Vector& v);

protected:
  double x;
  double y;
  double z;
};

#endif

