#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
///////gcc spin.c -o spin -lm//// make sure the math modle can run
typedef struct {
    double x; double y; double z;
} vector;


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
vector H = {0.0, 0.0, 10.0};
double H_ext=10.0;
double lumba = 0.05;
double gama = 1.76 * 1E11; // Update gamma to the larger value
double dt = 1E-15;
double thirty=30.0;
//vector f(vector spin){
//
//}
vector f(double fspinx, double fspiny, double fspinz,double dtime) {
    vector fspin = {fspinx, fspiny, fspinz};
    vector eq18_1 = crossProduct(fspin, H);
    vector eq18_2 = crossProduct(scalarMultiply(fspin, lumba), eq18_1);
    vector d_spin = scalarMultiply(add(eq18_1, eq18_2), -1.0 * gama / (1.0 + lumba * lumba));
//    fspin = add(fspin, scalarMultiply(d_spin, dtime));
    return d_spin;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int main(void)
{
  double analytic_factor1 ,analytic_factor2;
  int time_step = 1000000;
  int width=2, len=2;
  int num_spin = width*len;
  vector multi_spin[num_spin]; 
  ///Use 2 array of spin, one for update t(x+1) spin, one for recording the t(x) spin information
  ///dynamic memory C, use pointer to record the information; sth like ptr = (int*) malloc(100 * sizeof(int));
  vector multi_spin_position[num_spin];
  vector totalspin={0.0,0.0,0.0};
  for(int i = 0; i < num_spin; i++) {
     multi_spin[i].x = 1.0;
     multi_spin[i].y = 0.0;
     multi_spin[i].z = 0.0;
//     printf("%.16f %.16f %.16f\n",multi_spin[i].x,multi_spin[i].y,multi_spin[i].z);
  }
  for(int i = 0; i < num_spin; i++) {
     multi_spin_position[i].x = i%width;
     multi_spin_position[i].y = i/width;
     multi_spin_position[i].z = 0.0;
     printf("%.16f %.16f %.16f\n",multi_spin_position[i].x,multi_spin_position[i].y,multi_spin[i].z);
  }

  vector analytic_spin={1.0 , 0.0, 0.0};
  vector rk1 ={0.0,0.0,0.0};
  vector rk2 ={0.0,0.0,0.0};
  vector rk3 ={0.0,0.0,0.0};
  vector rk4 ={0.0,0.0,0.0};

  /////////////////////open file//////////////////////////////
  FILE *fptr;
  fptr = fopen("result.txt", "w");
  if (fptr == NULL) {
    printf("Error opening file.\n");
    return 1;
  }
 //////////////////////iteration of spin/////////////////////

  for(int i = 0; i < time_step; i++) {
    analytic_factor1 = H_ext*lumba*gama*dt*(i)/(1+lumba*lumba);
    analytic_factor2 = H_ext*gama*dt*(i)/(1+lumba*lumba);

    analytic_spin.x = sech(analytic_factor1)*cos(analytic_factor2);
    analytic_spin.y = sech(analytic_factor1)*sin(analytic_factor2);
    analytic_spin.z = tanh(analytic_factor1);
    fprintf(fptr, "%.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", dt*i, totalspin.x, totalspin.y, totalspin.z, analytic_spin.x*num_spin,analytic_spin.y*num_spin,analytic_spin.z*num_spin);
    totalspin.x=0.0,totalspin.y=0.0,totalspin.z=0.0;

    for (int j = 0; j < num_spin; j++) {
      rk1 = f(multi_spin[j].x, multi_spin[j].y, multi_spin[j].z, 0.0);
      rk1 = scalarMultiply(rk1, dt);

      rk2 = f(multi_spin[j].x + rk1.x/2.0, multi_spin[j].y + rk1.y/2.0, multi_spin[j].z + rk1.z/2.0, 0.5*dt);
      rk2 = scalarMultiply(rk2, dt);

      rk3 = f(multi_spin[j].x + rk2.x/2.0, multi_spin[j].y + rk2.y/2.0, multi_spin[j].z + rk2.z/2.0, 0.5*dt);
      rk3 = scalarMultiply(rk3, dt);

      rk4 = f(multi_spin[j].x + rk3.x, multi_spin[j].y + rk3.y, multi_spin[j].z + rk3.z, dt);
      rk4 = scalarMultiply(rk4, dt);

      multi_spin[j].x = multi_spin[j].x + (rk1.x + 2*rk2.x + 2*rk3.x + rk4.x)/6.0;
      multi_spin[j].y = multi_spin[j].y + (rk1.y + 2*rk2.y + 2*rk3.y + rk4.y)/6.0;
      multi_spin[j].z = multi_spin[j].z + (rk1.z + 2*rk2.z + 2*rk3.z + rk4.z)/6.0;
      multi_spin[j] = normalized(multi_spin[j]);
    }
    for (int j = 0; j < num_spin; j++) {
      totalspin.x+=multi_spin[j].x; totalspin.y+=multi_spin[j].y; totalspin.z+=multi_spin[j].z;
    }
  }  
  
  fclose(fptr);
  return 0;
}



