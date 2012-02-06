/*
  Header file for body object in N-body simulation
 */

#ifndef _body_h_
#define _body_h_

#include "Vector.h"
#include <iostream>
using std::ostream;

struct body_message {
  double id;
  double mass;
  double position[3];
  double velocity[3];
  double force[3];
};

const int body_message_size = sizeof(body_message)/sizeof(double);

class Body {
public:
  // constructors
  Body();
  Body(double m, Vector x, Vector v);
  Body(const body_message *m);

  // "getters" and "setters"
  double mass();
  void setMass(double m);
  Vector position();
  void setPosition(const Vector &v);
  Vector velocity();
  void setVelocity(const Vector &v);
  Vector force();
  void setForce(const Vector &v);

  // other instance methods
  void clearForce();
  void addForce(const Body &b);
  void move(int dt);
  void makeMessage(body_message *m);

  double* message(double o[]);

  // output operator
  friend ostream& operator<<(ostream& s, const Body& b);

protected:
  double _mass;
  Vector _position;
  Vector _velocity;
  Vector _force;
};

#endif
