#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
///////gcc spin.c -o spin -lm//// make sure the math modle can run
typedef struct {
    double x; double y; double z;
} vector;
typedef struct {
    int neighbors[4];
} VectorArray;


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

void set_initial_spin(vector* multi_spin,int num_spin){

  for(int i = 0; i < num_spin; i++) {
     multi_spin[i].x = 1.0;
     multi_spin[i].y = 0.0;
     multi_spin[i].z = 0.0;
     printf("%.16f %.16f %.16f\n",multi_spin[i].x,multi_spin[i].y,multi_spin[i].z);
  }
}

void set_spin_position(vector* multi_spin_position,int num_spin, int width,int len){

  for(int i = 0; i < num_spin; i++) {
    multi_spin_position[i].x = i%width;
    multi_spin_position[i].y = i/width;
    multi_spin_position[i].z = 0.0;
    printf("%.16f %.16f %.16f position\n",multi_spin_position[i].x,multi_spin_position[i].y,multi_spin_position[i].z);
  }
}


void find_spin_neighbor(VectorArray* neighbor_index,int num_spin,vector* multi_spin_position, int width,int len){

  for (int i = 0; i < num_spin; i++){
    //-x方向
    if (i%width==0){
       neighbor_index[i].neighbors[0] = i+width-1;
    } else{
       neighbor_index[i].neighbors[0] = i-1;
    }
    //+x方向
    if (i%width==4){
       neighbor_index[i].neighbors[1] = i-width+1;
    } else{
       neighbor_index[i].neighbors[1] = i+1;
    }
    //-y方向
    if (i/width==0){
       neighbor_index[i].neighbors[2] = i+width*(len-1);
    } else{
       neighbor_index[i].neighbors[2] = i-width;
    }
    //+y方向
    if (i/width==(len-1)){
       neighbor_index[i].neighbors[3] = i-(width*(len-1));
    } else{
       neighbor_index[i].neighbors[3] = i+width;
    }
//    printf("%d %d %d %d %d, neighbor index\n",i, neighbor_index[i].neighbors[0],neighbor_index[i].neighbors[1],neighbor_index[i].neighbors[2],neighbor_index[i].neighbors[3]);
  }
}

void calculate_H_eff(vector* H_eff,VectorArray* neighbor_index,vector* multi_spin, int Nth_spin){
   int neighbor1=neighbor_index[Nth_spin].neighbors[0];
   int neighbor2=neighbor_index[Nth_spin].neighbors[1];
   int neighbor3=neighbor_index[Nth_spin].neighbors[2];
   int neighbor4=neighbor_index[Nth_spin].neighbors[3];
   double Jxx=0.1;
   double Jyy=0.1;
   double Jzz=0.1;
   H_eff[Nth_spin].x =H.x - Jxx* (multi_spin[Nth_spin].x*multi_spin[neighbor1].x + 
	                          multi_spin[Nth_spin].x*multi_spin[neighbor2].x +
			          multi_spin[Nth_spin].x*multi_spin[neighbor3].x +
			          multi_spin[Nth_spin].x*multi_spin[neighbor4].x ) ;
   H_eff[Nth_spin].y =H.y - Jyy* (multi_spin[Nth_spin].y*multi_spin[neighbor1].y +
                                  multi_spin[Nth_spin].y*multi_spin[neighbor2].y +
                                  multi_spin[Nth_spin].y*multi_spin[neighbor3].y +
                                  multi_spin[Nth_spin].y*multi_spin[neighbor4].y ) ;
	                    
   H_eff[Nth_spin].z =H.z - Jzz* (multi_spin[Nth_spin].x*multi_spin[neighbor1].z +
                                  multi_spin[Nth_spin].x*multi_spin[neighbor2].z +
                                  multi_spin[Nth_spin].x*multi_spin[neighbor3].z +
                                  multi_spin[Nth_spin].x*multi_spin[neighbor4].z ) ;
//   printf("%d %d %d %d, neighbor index\n", neighbor1,neighbor2,neighbor3,neighbor4);
//   H_eff.x = multi_spin[Nth_spin].x*(multi_spin[].x)

}

////////////////////////////////////////////////////////////////////////////////////////////////

int main(void)
{
  double analytic_factor1 ,analytic_factor2;
  int time_step = 2;//1000000;
  int width=5, len=5;
  int num_spin = width*len;
  vector H_eff[num_spin];
  vector multi_spin[num_spin]; 
  ///Use 2 array of spin, one for update t(x+1) spin, one for recording the t(x) spin information
  ///dynamic memory C, use pointer to record the information; sth like ptr = (int*) malloc(100 * sizeof(int));
  vector multi_spin_position[num_spin];
  VectorArray neighbor_index[num_spin];
  vector totalspin={0.0,0.0,0.0};

  set_initial_spin(multi_spin,num_spin);
  for (int i = 0; i < num_spin; i++) {
     printf("Spin  %d: %.16f %.16f %.16f\n", i, multi_spin[i].x, multi_spin[i].y, multi_spin[i].z);
  }
  set_spin_position(multi_spin_position,num_spin, width, len);

  find_spin_neighbor(neighbor_index,num_spin, multi_spin_position, width,len);
  printf("%d  %d\n", neighbor_index[1].neighbors[0], neighbor_index[5].neighbors[3]);
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
    fprintf(fptr, "%.16f %.16f %.16f %.16f\n", dt*i, totalspin.x, totalspin.y, totalspin.z);
    totalspin.x=0.0,totalspin.y=0.0,totalspin.z=0.0;
    
    for (int j = 0; j < num_spin; j++) {
      calculate_H_eff(H_eff, neighbor_index, multi_spin, j);
      printf("%.16f %.16f %.16f \n",H_eff[j].x, H_eff[j].y,H_eff[j].z);


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



