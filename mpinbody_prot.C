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

#include <math.h>
using std:sqrt;

#include <stdlib.h>
#include <fstream>
using std::ifstream;

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

const int NP = sizeof(b)/sizeof(Body);

Vector calc_accel(Body i, double x, double y, double z, double jmass);
int output(int t, int h, int w, double buf[]);

int nprocs,myid;

int main(int argc, char *argv[]) {

MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Status s;

//each thread reads from file... only the lines it needs
/*
 * open file
 * allocate memory for bodies (first line contains num bodies)
 * float* bod_array = (float*) malloc numbodies / nprocs
 * int row = 0;
 * int bod_count = 0;
 * for(int i = 0; row < numbodies; i++){
 *  row = i * nprocs + myid
 *  bod_array[i] = read file line_number[row]
 *  bod_count++;
 * }
 *
 *put my data into an array of the form:
 *bod_count, mybod0.m, mybod0.x, mybod0.y, mybod0.z...
 *size of array will be size_of(bod_array) + 1 // for bod_count
 *make reception array nprocs * (size_of(bod_array) + 1)
 *mpi_allgather(send size_of(bod_array) + 1, receive size_of(bod_array) + 1)
 *
 *calculate new coords / velocities for my bodies from received array
 *for(each of my bodies){
 *  for(each of the received arrays){
 *    for(each of the values in each of the received arrays){
 *      mod_accel(my body, value)
 *    }
 *  }
 *}
 *
 */
  double *local_bod_array;
  int local_bod_idx=0;
  int line_number=0;
  int row=myid;
  int total_bod_count,bods_per_proc,local_bod_count;
  string line;
  ifstream input_data;
  input_data.open("galaxy.dat");
  if(input_data.is_open()){
    if(input_data.good()){
      input_data >> total_bod_count;
      bods_per_proc = total_bod_count/nprocs;
      printf("total_bod_count is %i!\n",total_bod_count);
      printf("bods_per_proc is %i!\n",bods_per_proc);
    }
    local_bod_array = (double*) malloc(sizeof(double) * 7 * bods_per_proc);
    getline(input_data,line);//moving to the first row of data
    while(input_data.good() && row < total_bod_count){
      getline(input_data,line);
      if(row==line_number){
        printf("I am thread %i at row %i, storing body %i\n",myid,row,local_bod_idx);
        size_t start = 0;
        size_t end = 0;
        for(int j=0;j<7;j++){
          end = line.find(" ",start+1);
          double dub = atof(line.substr(start,end).c_str());
          *(local_bod_array+(local_bod_idx*7)+j) = dub;
          printf("start: %lu, end: %lu, just added measurement %f\n",start,end,dub);
          start = end;
        }
        row+=nprocs;
        local_bod_idx++;
      }
      line_number++;
    }
  }
local_bod_count = local_bod_idx;//it's a better name for the var from here on down
printf("Thread %i holding %i bodies\n",myid,local_bod_count);
for(int i = 0; i < local_bod_count; i++){
printf("bod %i mass: %f\n",i,*(local_bod_array+(i*7)+0));
printf("bod %i x: %f\n",i,*(local_bod_array+(i*7)+1));
printf("bod %i y: %f\n",i,*(local_bod_array+(i*7)+2));
printf("bod %i z: %f\n",i,*(local_bod_array+(i*7)+3));
printf("bod %i vx: %f\n",i,*(local_bod_array+(i*7)+4));
printf("bod %i vy: %f\n",i,*(local_bod_array+(i*7)+5));
printf("bod %i vz: %f\n",i,*(local_bod_array+(i*7)+6));
}

int sendbuf_size=(4 * bods_per_proc) + 1;
printf("sendbuf_size is %i\n",sendbuf_size);

double *sendbuf,*recvbuf;
sendbuf = (double*) malloc(sizeof(double) * sendbuf_size);
recvbuf = (double*) malloc(sizeof(double) * sendbuf_size * nprocs); //different from total_bod_count!!

*(sendbuf)= (double) local_bod_count;

for(int t=0;t<T;t++){
  //sendbuf arrangement: 0->local_bod_count, 1->bod0.m, 2->bod0.x, 3->bod0.y, 4->bod0.z, 5->bod1.m, ... , (4n+1)->bodn.m, (4n+2)->bodn.x, (4n+3)->bodn.y, (4n+4)->bodn.z
  for(int i = 0; i < local_bod_count; i++){
    *(sendbuf+(i*4)+1)=*(local_bod_array+(i*7)+0);
    *(sendbuf+(i*4)+2)=*(local_bod_array+(i*7)+1);
    *(sendbuf+(i*4)+3)=*(local_bod_array+(i*7)+2);
    *(sendbuf+(i*4)+4)=*(local_bod_array+(i*7)+3);
  }

  printf("local_bod_count is %i\n",local_bod_count);

  for(int i = 0; i < sendbuf_size; i++){
    printf("sendbuf %i is %f\n",i,*(sendbuf+i));
  }

  MPI_Allgather(sendbuf,sendbuf_size,MPI_DOUBLE,recvbuf,sendbuf_size,MPI_DOUBLE,MPI_COMM_WORLD);
  
  printf("Listing all bod counts\n");
  for(int i=0; i<nprocs;i++){
    if(i != myid){//skip my own bodies
      int foreign_bod_count=(int)  *(recvbuf+(i * sendbuf_size));
      printf("sendbuf bod count %i is %f\n",i,foreign_bod_count);
      for(int j = 0; j < foreign_bod_count; j++){
        double foreign_m = *(recvbuf+(i * sendbuf_size)+(4 * j)+1);
        double foreign_x = *(recvbuf+(i * sendbuf_size)+(4 * j)+2);
        double foreign_y = *(recvbuf+(i * sendbuf_size)+(4 * j)+3);
        double foreign_z = *(recvbuf+(i * sendbuf_size)+(4 * j)+4);
        for(int k = 0; k < local_bod_count; k++){
          double local_m = *(local_bod_array+(k*7)+0);
          double local_x = *(local_bod_array+(k*7)+1);
          double local_y = *(local_bod_array+(k*7)+2);
          double local_z = *(local_bod_array+(k*7)+3);
          double local_vx = *(local_bod_array+(k*7)+4);
          double local_vy = *(local_bod_array+(k*7)+5);
          double local_vz = *(local_bod_array+(k*7)+6);
          // Vector k = i.position() - j;
          double dx = local_x - foreign_x;
          double dy = local_y - foreign_y;
          double dz = local_z - foreign_z;
          // double l = (k * (k * k)).norm();
          double dot_p = dx*dx + dy*dy + dz*dz;
          dx *= dot_p;
          dy *= dot_p;
          dz *= dot_p;
          double l = sqrt(dx*dx + dy*dy + dz*dz);
          // Vector m = (k * ((double) 1/l)) * -1 * G * jmass;
          double accel_x = dx * ((double) 1/l) * -1 * foreign_m;
          double accel_y = dy * ((double) 1/l) * -1 * foreign_m;
          double accel_z = dz * ((double) 1/l) * -1 * foreign_m;
          // Vector newVeloc = b[myid].velocity() += accel * DT;
          double new_vx = local_vx += accel_x * DT;
          double new_vy = local_vy += accel_y * DT;
          double new_vz = local_vz += accel_z * DT;
          // Vector newPosit = b[myid].position() += newVeloc * DT;
          double new_x = local_x += new_vx * DT;
          double new_y = local_y += new_vy * DT;         
          double new_z = local_z += new_vz * DT;
          // b[myid].setPosition(newPosit);
          *(local_bod_array+(k*7)+1) = new_x;
          *(local_bod_array+(k*7)+2) = new_y;
          *(local_bod_array+(k*7)+3) = new_z;
          // b[myid].setVelocity(newVeloc);
          *(local_bod_array+(k*7)+4) = new_vx;
          *(local_bod_array+(k*7)+5) = new_vy;
          *(local_bod_array+(k*7)+6) = new_vz;
        }   
      }
    }
  }
}
  
//for each timestep
//
//  do the allgather
//
//  iterate through and pull out values
//
//  recalculate coords for my bodies
//
//  rebuild send buffer 

/*
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
    //printf("myid: %i, x: %f y: %f z: %f mass: %f\n",myid,send_buffer[0],send_buffer[1],send_buffer[2],send_buffer[3]);
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
*/
  free(sendbuf);
  free(recvbuf);
  free(bod_array);
  return MPI_Finalize();
}

int output(int t, int h,int w,double buf[]){
  cout << t << " ";
  for(int i = 0; i < h; i++){
    for(int j = 0; j < w; j++){
      cout << buf[i*4+j] << " ";
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
