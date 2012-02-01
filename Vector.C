/*						
  Implementation of 3-dimensional vectors	
						
  John Conery					
  CIS 455/555					
  University of Oregon				
 */

#include "Vector.h"

#include <math.h>

// Default constructor -- set vector to zero (i.e. (0,0,0)):

Vector::Vector() {
  x = y = z = 0.0;
}

// Copy constructor -- the new vector has same (x,y,z) coordinates
// as an existing vector V:

Vector::Vector(const Vector& v) {
  x = v.x;
  y = v.y;
  z = v.z;
}

// Basic user constructor -- make a new vector from specified x, 
// y, and z coordinates:

Vector::Vector(double px, double py, double pz) {
  x = px;
  y = py;
  z = pz;
}

// Alternate constructor -- make a new vector from an array of 3
// doubles (intended for reconstituting vector arriving in an MPI message):

Vector::Vector(const double* v) {
  x = v[0];
  y = v[1];
  z = v[2];
}

// Assignment operator

Vector& Vector::operator=(const Vector& v) {
  x = v.x;
  y = v.y;
  z = v.z;
  return *this;
}

Vector& Vector::operator=(double a) {
  x = a;
  y = a;
  z = a;
  return *this;
}

// Assignment versions of arithemetic operations -- these operators modify
// existing vector objects.  They all return references to vectors so
// they can be used to define the corresponding stand-alone operators.

Vector& Vector::operator+=(const Vector& v) {
  x += v.x;
  y += v.y;
  z += v.z;
  return *this;
}

Vector& Vector::operator+=(double a) {
  x += a;
  y += a;
  z += a;
  return *this;
}

Vector& Vector::operator-=(const Vector& v) {
  x -= v.x;
  y -= v.y;
  z -= v.z;
  return *this;
}

Vector& Vector::operator-=(double a) {
  x -= a;
  y -= a;
  z -= a;
  return *this;
}

Vector& Vector::operator*=(double a) {
  x *= a;
  y *= a;
  z *= a;
  return *this;
}

// Compute the Euclidean norm of the vector:

double Vector::norm() {
  return sqrt(x*x + y*y + z*z);
}

// Fill an array of doubles with values from this vector; intended
// to be used for sending a vector object in an MPI message.

void Vector::makeMessage(double *msg) {
    msg[0] = x;
    msg[1] = y;
    msg[2] = z;
}

// Stand-alone versions of the basic operators -- they are defined
// as friends, but are implemented in terms of the corresponding
// assignment operators.  See "More Effective C++", Item 22.

Vector operator+(const Vector &v1, const Vector& v2) {
  return Vector(v1) += v2;
}

Vector operator+(const Vector &v1, double a) {
  return Vector(v1) += a;
}

Vector operator-(const Vector &v1, const Vector& v2) {
  return Vector(v1) -= v2;
}

Vector operator-(const Vector &v1, double a) {
  return Vector(v1) -= a;
}

Vector operator*(const Vector &v1, double a) {
  return Vector(v1) *= a;
}

// Vector multiplication (dot product), e.g. a = v0 * v1.

double operator*(const Vector& v1, const Vector& v2) {
  return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// Print a vector on a stream.  Print as three doubles
// separated by spaces (suitable for reading by Matlab):

ostream& operator<<(ostream& s, const Vector& v) {
  s << v.x << " " << v.y << " " << v.z;
  return s;
}

