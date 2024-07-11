#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
///////gcc spin.c -o spin -lm//// make sure the math modle can run
typedef struct {
    float x; float y; float z;
} vector;


/////////////////////////define vector caulation function////////////////////////////////////////
vector crossProduct(vector a, vector b) {
    vector result;
    result.x = a.y * b.z - a.z * b.y; result.y = a.z * b.x - a.x * b.z; result.z = a.x * b.y - a.y * b.x;
    return result;
}

vector scalarMultiply(vector a, float b) {
    vector result;
    result.x = a.x*b; result.y = a.y*b; result.z = a.z*b;
    return result;
}

vector add(vector a, vector b) {
    vector result;
    result.x = a.x+b.x; result.y = a.y+b.y; result.z = a.z+b.z;
    return result;
}
double sech(double x) {
    return 1 / cosh(x);
}
vector H = {0.0, 0.0, 10.0};
double H_ext=10.0;
double lumba = 0.1;
double gama = 1.76 * 1E11; // Update gamma to the larger value
double dt = 1E-15;

vector f(double fspinx, double fspiny, double fspinz,double dtime) {
    vector fspin = {fspinx, fspiny, fspinz};
    vector eq18_1 = crossProduct(fspin, H);
    vector eq18_2 = crossProduct(scalarMultiply(fspin, lumba), eq18_1);
    vector d_spin = scalarMultiply(add(eq18_1, eq18_2), -1 * gama / (1 + lumba * lumba));
    fspin = add(fspin, scalarMultiply(d_spin, dtime));
    return fspin;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int main(void)
{
  double analytic_factor1 ,analytic_factor2;
  int time_step = 400000;

  vector analytic_spin={1.0 , 0.0, 0.0};

  vector spin = {1.0, 0.0, 0.0};
  vector H = {0.0, 0.0, H_ext};
  vector d_spin = {0.0, 0.0, 0.0};
  vector cross_1 = {0.0, 0.0, 0.0};
  vector cross_2 = {0.0, 0.0, 0.0};

  /////////////////////open file//////////////////////////////
  FILE *fptr;
  fptr = fopen("result.txt", "w");
  if (fptr == NULL) {
    printf("Error opening file.\n");
    return 1;
  }
 //////////////////////iteration of spin/////////////////////

  for(int i = 0; i < time_step; i++) {
    spin = f(spin.x, spin.y, spin.z,dt);
    analytic_factor1 = H_ext*lumba*gama*1*dt*i/(1+lumba*lumba);
    analytic_factor2 = H_ext*gama*1*dt*i/(1+lumba*lumba);
    analytic_spin.x = sech(analytic_factor1)*cos(analytic_factor2);
    analytic_spin.y = sech(analytic_factor1)*sin(analytic_factor2);
    analytic_spin.z = tanh(analytic_factor1);
    
    fprintf(fptr, "%.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", dt*i, spin.x, spin.y, spin.z, analytic_spin.x,analytic_spin.y,analytic_spin.z);
  }  
  fclose(fptr);
  return 0;
}


