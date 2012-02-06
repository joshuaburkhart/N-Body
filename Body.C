/*
  Implementation of Body class for N-body simulation

	*** Distributed version ***

	empty bodies for most methods -- add the code
	for the methods you need

*/

#include "Body.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

// G = gravitational constant, in newtons-kg^2/m^2

const double G = 6.67E-11;

// Constructors --

// The default constructor makes a body with no mass and 
// 3D zero vectors for position and velocity.

Body::Body() {
  _mass = 0.0;
  _position = 0.0;
  _velocity = 0.0;
  _force = 0.0;
}

// Create a body with specified mass, position, and velocity; the
// force is initalized to the zero vector.

Body::Body(double m, Vector x, Vector v) {
  _mass = m;
  _position = x;
  _velocity = v;
  _force = 0.0;
}

// Initialize a new body from the contents of an MPI message

Body::Body(const body_message *m) {

  /* your code here */

}

// "Getters" and "Setters"

double Body::mass() {
  return _mass;
}

void Body::setMass(double m) {
  _mass = m;
}

Vector Body::position() {
  return _position;
}

void Body::setPosition(const Vector &v) {
  _position = v;
}

Vector Body::velocity() {
  return _velocity;
}

void Body::setVelocity(const Vector &v) {
  _velocity = v;
}

Vector Body::force() {
  return _force;
}

void Body::setForce(const Vector &v) {
  _force = v;
}

// Instance Methods

void Body::clearForce() {

  /* your code here */

}

void Body::addForce(const Body &b) {

  /* your code here */

}

void Body::move(int dt) {

  /* your code here */

}

void Body::makeMessage(body_message *m) {

  /* your code here */

}

double* Body::message(double o[]){
  o[0]=_position.getX();
  o[1]=_position.getY();
  o[2]=_position.getZ();
  o[3]=_mass;
  return o;
} 

// Output operator 

ostream& operator<<(ostream& s, const Body& b) {
  s << b._position;
  return s;
}
