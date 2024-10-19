#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
//////////// gcc main.c -o main -lm ////////////

/// Define vector type ///
typedef struct {
    double x; double y; double z;
} vector;
typedef struct {
    int neighbors[4];
} VectorArray;

/// Set Your simulation parameters ///
vector H = {0.0, 0.0, -10.0};
double initial_sx=1E-6, initial_sy=0.0, initial_sz=1.0;
double Jxx=0.3, Jyy=0.3, Jzz=0.3; 
double ku =-1; 

int time_step = 3000000;
int width=10, len=10;
double lumba = 0.05;
double gama = 1.76 * 1E11;
double dt = 1E-15;

/////////////////////////define vector caulation function////////////////////////////////////////
vector crossProduct(vector a, vector b) {
    vector result;
    result.x = a.y * b.z - a.z * b.y; result.y = a.z * b.x - a.x * b.z; result.z = a.x * b.y - a.y * b.x;
    return result;
}

vector scalarMultiply(vector a, double b) {
    vector result;
    result.x = a.x*b; result.y = a.y*b; result.z = a.z*b;
    return result;
}

vector scalarDivide(vector a, double b) {
    vector result;
    result.x = a.x/b; result.y = a.y/b; result.z = a.z/b;
    return result;
}


vector add(vector a, vector b) {
    vector result;
    result.x = a.x+b.x; result.y = a.y+b.y; result.z = a.z+b.z;
    return result;
}
vector minus(vector a, vector b) {
    vector result;
    result.x = a.x-b.x; result.y = a.y-b.y; result.z = a.z-b.z;
    return result;
}

vector normalized(vector a){
    a=scalarMultiply(a, 1.0/sqrt(a.x*a.x+a.y*a.y+a.z*a.z));
    return a;
}
double sech(double x) {
    return 1 / cosh(x);
}



vector f(double fspinx, double fspiny, double fspinz,double dtime, vector H_eff) {
    vector fspin = {fspinx, fspiny, fspinz};
    vector eq18_1 = crossProduct(fspin, H_eff);
    vector eq18_2 = crossProduct(scalarMultiply(fspin, lumba), eq18_1);
    vector d_spin = scalarMultiply(add(eq18_1, eq18_2), -1.0 * gama / (1.0 + lumba * lumba));
    return d_spin;
}

void set_initial_spin(vector* multi_spin,int num_spin){

  for(int i = 0; i < num_spin; i++) {
     multi_spin[i].x = initial_sx;
     multi_spin[i].y = initial_sy;
     multi_spin[i].z = initial_sz;
  }
}

void set_spin_position(vector* multi_spin_position,int num_spin, int width,int len){

  for(int i = 0; i < num_spin; i++) {
    multi_spin_position[i].x = i%width;
    multi_spin_position[i].y = i/width;
    multi_spin_position[i].z = 0.0;
  }
}


void find_spin_neighbor(VectorArray* neighbor_index,int num_spin,vector* multi_spin_position, int width,int len){

  for (int i = 0; i < num_spin; i++){
    //-x neighbor
    if (i%width==0){
       neighbor_index[i].neighbors[0] = i+width-1;
    } else{
       neighbor_index[i].neighbors[0] = i-1;
    }

    //+x neighbor
    if (i%width==4){
       neighbor_index[i].neighbors[1] = i-width+1;
    } else{
       neighbor_index[i].neighbors[1] = i+1;
    }
    //-y neighbor
    if (i/width==0){
       neighbor_index[i].neighbors[2] = i+width*(len-1);
    } else{
       neighbor_index[i].neighbors[2] = i-width;
    }
    //+y neighbor
    if (i/width==(len-1)){
       neighbor_index[i].neighbors[3] = i-(width*(len-1));
    } else{
       neighbor_index[i].neighbors[3] = i+width;
    }
  }
}
void calculate_H1_eff_app(vector* H_eff,VectorArray* neighbor_index,vector* multi_spin, int Nth_spin){
   int neighbor1=neighbor_index[Nth_spin].neighbors[0];
   int neighbor2=neighbor_index[Nth_spin].neighbors[1];
   int neighbor3=neighbor_index[Nth_spin].neighbors[2];
   int neighbor4=neighbor_index[Nth_spin].neighbors[3];
   H_eff[Nth_spin].x =H.x ; H_eff[Nth_spin].y =H.y ; H_eff[Nth_spin].z =H.z ;
}

void calculate_H2_eff_exc(vector* H_eff,VectorArray* neighbor_index,vector* multi_spin, int Nth_spin){
   int neighbor1=neighbor_index[Nth_spin].neighbors[0];
   int neighbor2=neighbor_index[Nth_spin].neighbors[1];
   int neighbor3=neighbor_index[Nth_spin].neighbors[2];
   int neighbor4=neighbor_index[Nth_spin].neighbors[3];
   H_eff[Nth_spin].x += Jxx* (multi_spin[neighbor1].x + 
	                       multi_spin[neighbor2].x +
			       multi_spin[neighbor3].x +
			       multi_spin[neighbor4].x ) ;

   H_eff[Nth_spin].y += Jyy* (multi_spin[neighbor1].y +
                               multi_spin[neighbor2].y +
                               multi_spin[neighbor3].y +
                               multi_spin[neighbor4].y ) ;
	                    
   H_eff[Nth_spin].z += Jzz* (multi_spin[neighbor1].z +
                               multi_spin[neighbor2].z +
                               multi_spin[neighbor3].z +
                               multi_spin[neighbor4].z ) ;
}

void calculate_H3_eff_ani(vector* H_eff,VectorArray* neighbor_index,vector* multi_spin, int Nth_spin){
   H_eff[Nth_spin].z += -ku*2.0*multi_spin[Nth_spin].z;
}


////////////////////////////////////////////////////////////////////////////////////////////////

int main(void)
{
  int num_spin = width*len;
  vector H_eff[num_spin];
  vector multi_spin[num_spin]; 
  vector multi_spin_position[num_spin];
  VectorArray neighbor_index[num_spin];
  vector totalspin={0.0,0.0,0.0};
  set_initial_spin(multi_spin,num_spin);
  set_spin_position(multi_spin_position,num_spin, width, len);
  find_spin_neighbor(neighbor_index,num_spin, multi_spin_position, width,len);
  vector analytic_spin={1.0 , 0.0, 0.0};
  vector rk1 [num_spin];
  vector rk2 [num_spin];
  vector rk3 [num_spin];
  vector rk4 [num_spin];

  ///open file///
  FILE *fptr;
  fptr = fopen("result.txt", "w");
  if (fptr == NULL) {
    printf("Error opening file.\n");
    return 1;
  }
 //////////////////////  Start Spin Iteration  //////////////////////

  for(int i = 0; i < time_step; i++) {
    if (i % 100 == 0) {
    fprintf(fptr, "%.16f %.16f %.16f %.16f \n", dt*i, totalspin.x, totalspin.y, totalspin.z);
    printf("%.16f %.16f %.16f %.16f result\n", dt*i, totalspin.x, totalspin.y, totalspin.z);
    }
    totalspin.x=0.0,totalspin.y=0.0,totalspin.z=0.0;
 ///  Run RK1  ///
    for (int j = 0; j < num_spin; j++) {
      calculate_H1_eff_app(H_eff, neighbor_index, multi_spin, j);
      calculate_H2_eff_exc(H_eff, neighbor_index, multi_spin, j);
      calculate_H3_eff_ani(H_eff, neighbor_index, multi_spin, j); 
      rk1[j] = f(multi_spin[j].x, multi_spin[j].y, multi_spin[j].z, 0.0,H_eff[j]);
      rk1[j] = scalarMultiply(rk1[j], dt);
    }
 ///  Run RK2  ///
    for (int j = 0; j < num_spin; j++) {
      calculate_H1_eff_app(H_eff, neighbor_index, rk1, j);
      calculate_H2_eff_exc(H_eff, neighbor_index, rk1, j);
      calculate_H3_eff_ani(H_eff, neighbor_index, multi_spin, j);
      rk2[j] = f(multi_spin[j].x + rk1[j].x/2.0, multi_spin[j].y + rk1[j].y/2.0, multi_spin[j].z + rk1[j].z/2.0, 0.5*dt,H_eff[j]);
      rk2[j] = scalarMultiply(rk2[j], dt);
    }
 ///  Run RK3  ///
    for (int j = 0; j < num_spin; j++) {
      calculate_H1_eff_app(H_eff, neighbor_index, rk2, j);
      calculate_H2_eff_exc(H_eff, neighbor_index, rk2, j);
      calculate_H3_eff_ani(H_eff, neighbor_index, multi_spin, j);
      rk3[j] = f(multi_spin[j].x + rk2[j].x/2.0, multi_spin[j].y + rk2[j].y/2.0, multi_spin[j].z + rk2[j].z/2.0, 0.5*dt,H_eff[j]);
      rk3[j] = scalarMultiply(rk3[j], dt);
    }
 ///  Run RK4  ///
    for (int j = 0; j < num_spin; j++) {
      calculate_H1_eff_app(H_eff, neighbor_index, rk3, j);
      calculate_H2_eff_exc(H_eff, neighbor_index, rk3, j);
      calculate_H3_eff_ani(H_eff, neighbor_index, multi_spin, j);
      rk4[j] = f(multi_spin[j].x + rk3[j].x, multi_spin[j].y + rk3[j].y, multi_spin[j].z + rk3[j].z, dt,H_eff[j]);
      rk4[j] = scalarMultiply(rk4[j], dt);
    }
 ///  Average RK1~RK4  and Do Spin Normalization of ///
    for (int j = 0; j < num_spin; j++) {
      multi_spin[j].x = multi_spin[j].x + (rk1[j].x + 2*rk2[j].x + 2*rk3[j].x + rk4[j].x)/6.0;
      multi_spin[j].y = multi_spin[j].y + (rk1[j].y + 2*rk2[j].y + 2*rk3[j].y + rk4[j].y)/6.0;
      multi_spin[j].z = multi_spin[j].z + (rk1[j].z + 2*rk2[j].z + 2*rk3[j].z + rk4[j].z)/6.0;
      multi_spin[j] = normalized(multi_spin[j]);
    }
 /// Record total Spin ///
    for (int j = 0; j < num_spin; j++) {
      totalspin.x+=multi_spin[j].x; totalspin.y+=multi_spin[j].y; totalspin.z+=multi_spin[j].z;
    }
  }
//////////////////////  End Spin Iteration  //////////////////////  
  fclose(fptr);
  return 0;
}



