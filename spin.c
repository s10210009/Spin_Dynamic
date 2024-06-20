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
////////////////////////////////////////////////////////////////////////////////////////////////

int main(void)
{
  double dt = 1E-15;
  double gamma = 1.76*(1E11), lumba=0.1; 
  double analytic_factor1 ,analytic_factor2;
  double H_ext = 10.0;
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
    vector cross_1 = crossProduct(spin, H);
    vector cross_2 = crossProduct(scalarMultiply(spin,lumba), cross_1);
    d_spin=scalarMultiply( add(cross_1,cross_2), -1*gamma/(1+lumba*lumba));
    spin = add(spin, scalarMultiply(d_spin,dt) );

    analytic_factor1 = H_ext*lumba*gamma*1*dt*i/(1+lumba*lumba);
    analytic_factor2 = H_ext*gamma*1*dt*i/(1+lumba*lumba);
    analytic_spin.x = sech(analytic_factor1)*cos(analytic_factor2);
    analytic_spin.y = sech(analytic_factor1)*sin(analytic_factor2);
    analytic_spin.z = tanh(analytic_factor1);
    
    fprintf(fptr, "%.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", dt*i, spin.x, spin.y, spin.z, analytic_spin.x,analytic_spin.y,analytic_spin.z);
  }  
  fclose(fptr);
  return 0;
}


