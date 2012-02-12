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
#include <stdlib.h>
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <mpi.h>

const int DT = 86459;     // time step = number of seconds in one day
const int T = 365;        // number of time steps to execute
double G = 6.67E-11;
int nprocs,myid;

int main(int argc, char *argv[]) {
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Status s;
  double *local_bod_array;
  int local_bod_idx=0;
  int line_number=0;
  int row=myid;
  int total_bod_count,bods_per_proc,local_bod_count;
  string line;
  ifstream input_data;
  input_data.open("/home13/jburkhar/N-Body/galaxy.dat");
  if(input_data.is_open()){
    if(input_data.good()){
      input_data >> total_bod_count;
      bods_per_proc = total_bod_count/nprocs;
    }
    local_bod_array = (double*) malloc(sizeof(double) * 7 * bods_per_proc);
    getline(input_data,line);//moving to the first row of data
    while(input_data.good() && row < total_bod_count){
      getline(input_data,line);
      if(row==line_number){
        size_t start = 0;
        size_t end = 0;
        for(int j=0;j<7;j++){
          end = line.find(" ",start+1);
          double dub = atof(line.substr(start,end).c_str());
          *(local_bod_array+(local_bod_idx*7)+j) = dub;
          start = end;
        }
        row+=nprocs;
        local_bod_idx++;
      }
      line_number++;
    }
  }
  local_bod_count = local_bod_idx;//it's a better name for the var from here on down
  int sendbuf_size=(4 * bods_per_proc) + 1;
  double *sendbuf,*recvbuf;
  sendbuf = (double*) malloc(sizeof(double) * sendbuf_size);
  recvbuf = (double*) malloc(sizeof(double) * sendbuf_size * nprocs); //different from total_bod_count!!
  *(sendbuf)= (double) local_bod_count;
  for(int t=0;t<T;t++){
    if(myid==0){
      cout << t << " ";
    }
    //sendbuf arrangement: 0->local_bod_count, 1->bod0.m, 2->bod0.x, 3->bod0.y, 4->bod0.z, 5->bod1.m, ... , (4n+1)->bodn.m, (4n+2)->bodn.x, (4n+3)->bodn.y, (4n+4)->bodn.z
    for(int i = 0; i < local_bod_count; i++){
      *(sendbuf+(i*4)+1)=*(local_bod_array+(i*7)+0);
      *(sendbuf+(i*4)+2)=*(local_bod_array+(i*7)+1);
      *(sendbuf+(i*4)+3)=*(local_bod_array+(i*7)+2);
      *(sendbuf+(i*4)+4)=*(local_bod_array+(i*7)+3);
    }
    MPI_Allgather(sendbuf,sendbuf_size,MPI_DOUBLE,recvbuf,sendbuf_size,MPI_DOUBLE,MPI_COMM_WORLD);
    for(int k = 0; k < local_bod_count; k++){
      double local_x;
      double local_y;
      double local_z;
      double local_vx;
      double local_vy;
      double local_vz;
      double accel_x=0;
      double accel_y=0;
      double accel_z=0;
      for(int i=0; i<nprocs;i++){
        int foreign_bod_count=(int)  *(recvbuf+(i * sendbuf_size));
        for(int j = 0; j < foreign_bod_count; j++){
          if(!(i == myid && j==k)){//skip bodies that are themselves
            double foreign_m = *(recvbuf+(i * sendbuf_size)+(4 * j)+1);
            double foreign_x = *(recvbuf+(i * sendbuf_size)+(4 * j)+2);
            double foreign_y = *(recvbuf+(i * sendbuf_size)+(4 * j)+3);
            double foreign_z = *(recvbuf+(i * sendbuf_size)+(4 * j)+4);
            double local_m = *(local_bod_array+(k*7)+0);
            local_x = *(local_bod_array+(k*7)+1);
            local_y = *(local_bod_array+(k*7)+2);
            local_z = *(local_bod_array+(k*7)+3);
            local_vx = *(local_bod_array+(k*7)+4);
            local_vy = *(local_bod_array+(k*7)+5);
            local_vz = *(local_bod_array+(k*7)+6);
            double dx = local_x - foreign_x;
            double dy = local_y - foreign_y;
            double dz = local_z - foreign_z;
            double dot_p = dx*dx + dy*dy + dz*dz;
            double tdx = dx * dot_p;
            double tdy = dy * dot_p;
            double tdz = dz * dot_p;
            double l = sqrt(tdx*tdx + tdy*tdy + tdz*tdz);
            double tax = dx * ((double) 1/l) * -1 * G * foreign_m;
            double tay = dy * ((double) 1/l) * -1 * G * foreign_m;
            double taz = dz * ((double) 1/l) * -1 * G * foreign_m;
            accel_x += tax;
            accel_y += tay;
            accel_z += taz;
          }   
        }
      }
      double new_vx = local_vx + accel_x * DT;
      double new_vy = local_vy + accel_y * DT;
      double new_vz = local_vz + accel_z * DT;
      double new_x = local_x + new_vx * DT;
      double new_y = local_y + new_vy * DT;         
      double new_z = local_z + new_vz * DT;
      *(local_bod_array+(k*7)+1) = new_x;
      *(local_bod_array+(k*7)+2) = new_y;
      *(local_bod_array+(k*7)+3) = new_z;
      *(local_bod_array+(k*7)+4) = new_vx;
      *(local_bod_array+(k*7)+5) = new_vy;
      *(local_bod_array+(k*7)+6) = new_vz;
    }
    //print coords for this timestep
    if(myid==0){
      for(int i = 0;i < bods_per_proc; i++){
        for(int j = 0; j < nprocs; j++){
          int tmp_bod_count=(int) *(recvbuf+(j * sendbuf_size)+0);
          if(i < tmp_bod_count){
            double tmp_x = *(recvbuf+(j * sendbuf_size)+(i * 4)+2);
            double tmp_y = *(recvbuf+(j * sendbuf_size)+(i * 4)+3);
            double tmp_z = *(recvbuf+(j * sendbuf_size)+(i * 4)+4);
            cout << tmp_x << " " << tmp_y << " " << tmp_z << " ";
          }
        }
      }
      cout << endl;
    }
  }  
  free(sendbuf);
  free(recvbuf);
  free(local_bod_array);
  return MPI_Finalize();
}
