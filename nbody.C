#include "Body.h"

#define N 10
#define TIME_STEPS 10
#define DELTA_T 1

Vector calc_accel(Body i, Body j);

int main(int argc, int argv[]){

double coords[N][3];
coords[0] = {0.0,9.0,8.2};
coords[1] = {2.0,8.0,6.2};
coords[2] = {3.0,7.0,2.2};
coords[3] = {4.0,6.0,2.2};
coords[4] = {5.0,5.0,3.2};
coords[5] = {6.0,4.0,4.2};
coords[6] = {7.0,3.0,5.2};
coords[7] = {8.0,2.0,6.2};
coords[8] = {9.0,1.0,7.2};
coords[9] = {10.0,1.0,9.2};

double veloc[N][3];
veloc[0] = {0.0,9.0,3.2};
veloc[1] = {2.0,8.0,5.2};
veloc[2] = {3.0,7.0,1.2};
veloc[3] = {4.0,6.0,7.2};
veloc[4] = {5.0,5.0,8.2};
veloc[5] = {6.0,4.0,9.2};
veloc[6] = {7.0,3.0,0.2};
veloc[7] = {8.0,2.0,1.2};
veloc[8] = {9.0,1.0,2.2};
veloc[9] = {10.0,1.0,1.2};

Vector position[N];
Vector velocity[N];
Body system[N];

for(int i=0; i < N; i++){
  position[i]=new Vector(coords[i][0],coords[i][1],coords[i][2]);
  velocity[i]=new Vector(veloc[i][0],veloc[i][1],veloc[i][2]);
  system[i]=new Body((double) (3*i+i),position[i],velocity[i]);
}

for(int t=0; t < TIME_STEPS; t++){
  printf("time %i: ",t);
  for(int i=0; i < N; i++){
    printf("body %i ",i);
    Vector accel[N-1];
    for(int j=0; j < N; j++){
      if(i != j){
        accel[i] += calc_accel(system[i],system[j]);
      }
      //else continue
    }
    Vector newVeloc = system[i].velocity() += accel[i] * DELTA_T;
    Vector newPosit = system[i].position() += newVeloc * DELTA_T;
    system[i].setVelocity(newVeloc);
    system[i].setPosition(newPosit);
    cout << "(pos " << system[i].position() << ")";
    cout << "(vel " << system[i].velocity() << ")" << endl;
  }
}

return 0;
}

Vector calc_accel(Body i, Body j){

Vector k = i.position() - j.position();
Vector l = abs(k * k * k);
Vector m = -G * j.mass() * (k/l);

return m;
}
