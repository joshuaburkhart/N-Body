/*						
  Demo program for Vector class			
						
  John Conery					
  University of Oregon				
  March 1999					
*/

#include "Vector.h"

#include <iostream>
using std::cout;
using std::endl;

Vector foo(Vector a, Vector b) {
  return a + b*2;
}

int main(int argc, char *argv[]) {

  cout << "Constructors" << endl;
  Vector v0;
  Vector v1(1,2,3);
  cout << "Vector v0;        v0 = " << v0 << endl;
  cout << "Vector v1(1,2,3); v1 = " << v1 << endl;
  cout << endl;

  cout << "Assignments" << endl;
  Vector v2 = Vector(2,4,6);
  cout << "v2 = Vector(2,4,6) = " << v2 << endl;
  v0 = v1;
  cout << "v0 = v1; v0 = " << v0 << endl;
  v0 = 3;
  cout << "v0 = 3; v0 = " << v0 << endl;
  cout << endl;
  
  cout << "Assignment operators" << endl;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  v2 += v1;
  cout << "  v2 += v1; v2 = " << v2 << endl;
  v2 += 1;
  cout << "  v2 += 1; v2 = " << v2 << endl;
  v2 -= 1;
  cout << "  v2 -= 1; v2 = " << v2 << endl;
  v2 -= v1;
  cout << "  v2 -= v1; v2 = " << v2 << endl;
  v2 *= 2;
  cout << "  v2 *= 2; v2 = " << v2 << endl;
  cout << endl;
  
  cout << "Stand-alone operators" << endl;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  Vector v3;
  v3 = v1 + v2;
  cout << "  v3 = v1 + v2; v3 = " << v3 << endl;
  cout << "  v1 + v2 = " << v1 + v2 << endl;
  v3 = v1 + 10;
  cout << "  v3 = v1 + 10; v3 = " << v3 << endl;
  cout << "  v1 + 10 = " << v1 + 10 << endl;
  v3 = v1 - v2;
  cout << "  v3 = v1 - v2; v3 = " << v3 << endl;
  cout << "  v1 - v2 = " << v1 - v2 << endl;
  v3 = v1 - 10;
  cout << "  v3 = v1 - 10; v3 = " << v3 << endl;
  cout << "  v1 - 10 = " << v1 - 10 << endl;
  v3 = v1 * 10;
  cout << "  v3 = v1 * 10; v3 = " << v3 << endl;
  cout << "  v1 * 10 = " << v1 * 10 << endl;
  cout << endl;

  cout << "Inner product" << endl;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << "v1 * v2 = " << v1 * v2 << endl;
  cout << endl;

  cout << "Euclidean norm" << endl;
  cout << "v1.norm() = " << v1.norm() << endl;
  cout << endl;

  cout << "Passing vectors to functions (copy constructor)" << endl;
  cout << "foo(v1,v2) = v1 + 2*v2 = " << foo(v1,v2) << endl;
  cout << endl;

  cout << "MPI messages" << endl;
  double msg[3];
  cout << "v1 = " << v1 << endl;
  v1.makeMessage(msg);
  cout << "v1.makeMessage(msg); msg = ";
  cout << "[" << msg[0] << "," << msg[1] << "," << msg[2] << "]" << endl;
  msg[0] += 1000;
  msg[1] += 1000;
  msg[2] += 1000;
  Vector v4(msg);
  cout << "msg[i] += 1000; Vector v4(msg) = " << v4 << endl;
  cout << endl;
}
