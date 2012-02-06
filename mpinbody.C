/*
   N-body simulation, using particle-particle (PP) method.
   
   ** Distributed version -- add code to main() to compute
   ** interactions between bodies

   To verify the method, this program uses 10 bodies with masses
   of the sun and 9 planets.  The initial conditions were determined
   by downloading positions and velocities as of Jan 1, 1970 from the
   Solar System Database at JPL.

   John Conery
   CIS 455/555
   
   Updated Winter 2012 -- Use actual locations of planets (from JPL),
   print data in a format that can be read as an R frame object
*/

#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include "Vector.h"
#include "Body.h"
#include <mpi.h>

const int DT = 86459;     // time step = number of seconds in one day
const int T = 365;        // number of time steps to execute
const double G = 6.67E-11;

Body b[] = {
  /* sun */      Body( 1.9891E+30, Vector(0, 0, 0), Vector(0, 0, 0 ) ),
  /* mercury */  Body( 3.302E+23, Vector(3.83713E+10, 2.877025E+10, -1.175808E+09), Vector(-38787.67, 41093.05, 6918.461) ),
  /* venus */    Body( 4.8685E+24, Vector(-5.377313E+09, -1.085956E+11, -1.164748E+09), Vector(34741.48, -1865.747, -2031.506) ),
  /* earth */    Body( 5.9736E+24, Vector(-2.700743E+10, 1.446007E+11, 9686451), Vector(-29770.44, -5568.042, 0.3961261) ),
  /* mars */     Body( 6.4185E+23, Vector(1.983825E+11, 7.422924E+10, -3.334841E+09), Vector(-7557.626, 24761.27, 704.7457) ),
  /* jupiter */  Body( 1.89813E+27, Vector(-7.496502E+11, -3.201711E+11, 1.811155E+10), Vector(4982.522, -11417.83, -64.66531) ),
  /* saturn */   Body( 5.68319E+26, Vector(1.082806E+12, 8.510841E+11, -5.793461E+10), Vector(-6487.118, 7565.952, 125.4422) ),
  /* uranus */   Body( 8.68103E+25, Vector(-2.724616E+12, -2.894003E+11, 3.428801E+10), Vector(671.3469, -7099.093, -35.04028) ),
  /* neptune */  Body( 1.0241E+26, Vector(-2.328072E+12, -3.891086E+12, 1.337436E+11), Vector(4633.961, -2767.423, -49.57268) ),
  /* pluto */    Body( 1.314E+22, Vector(-4.551135E+12, 3.175277E+11, 1.282177E+12), Vector(635.998, -5762.115, 440.8821) )
};

string names[] = { "sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto" };

const int NP = sizeof(b)/sizeof(Body);

Vector calc_accel(Body i, double x, double y, double z, double jmass);
int output(int t, int h, int w, double buf[]);

int nprocs,myid;

int main(int argc, char *argv[]) {
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Status s;

  if(myid==0){
    int nn = ( sizeof(names) / sizeof(string *) );
    for (int i = 0; i < nn; i++) {
      cout << names[i] + "x "; 
      cout << names[i] + "y ";
      cout << names[i] + "z ";
    }
    cout << endl;
  }
  double rec_buffer[NP*4];
  for(int t = 0; t < T; t++){
    double send_buffer[4];
    b[myid].message(send_buffer);
    MPI_Allgather(send_buffer,4,MPI_DOUBLE,rec_buffer,4,MPI_DOUBLE,MPI_COMM_WORLD);
    if(myid==0){
      output(t,NP,3,rec_buffer); //put mass last in the message
    }
    Vector accel;   
    for(int i = 0; i < NP; i++){
      if(i != myid){
        accel += calc_accel(b[myid],rec_buffer[i*4+0],rec_buffer[i*4+1],rec_buffer[i*4+2],rec_buffer[i*4+3]);
      }
    }
    Vector newVeloc = b[myid].velocity() += accel * DT;
    Vector newPosit = b[myid].position() += newVeloc * DT;
    b[myid].setVelocity(newVeloc);
    b[myid].setPosition(newPosit);
  }
  return MPI_Finalize();
}

int output(int t, int h,int w,double buf[]){
  cout << t << " ";
  for(int i = 0; i < h; i++){
    for(int j = 0; j < w; j++){
      cout << buf[i*w+j] << " ";
    }
  }
  cout << endl;
}

Vector calc_accel(Body i, double x,double y, double z, double jmass){
  Vector j(x,y,z);
  Vector k = i.position() - j;
  double l = (k * (k * k)).norm();
  Vector m = (k * ((double) 1/l)) * -1 * G * jmass;
  return m;
}
