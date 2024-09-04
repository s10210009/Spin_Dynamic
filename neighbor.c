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

int main(void)
{
  int width=5, len=3;
  int num_spin = width*len;
  vector multi_spin[num_spin];
  vector multi_spin_position[num_spin];
  vector totalspin={0.0,0.0,0.0};

  VectorArray neighbor_position[num_spin];
  for(int i = 0; i < num_spin; i++) {
    multi_spin_position[i].x = i%width;
    multi_spin_position[i].y = i/width;
    multi_spin_position[i].z = 0.0;
    printf("%.16f %.16f %.16f te\n",multi_spin_position[i].x,multi_spin_position[i].y,multi_spin_position[i].z);
  }
  VectorArray neighbor_index[num_spin];
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
    if (i/width==len){
       neighbor_index[i].neighbors[3] = i-width*(len-1);
    } else{
       neighbor_index[i].neighbors[3] = i+width;
    }
					     //index
					     //0   1   2   3   4
					     //5   6   7   8   9
					     //10  11  12  13  14
  }
//  for(int i = 0; i < num_spin; i++) {
//    for (int j = 0; j < 4; j++){

//    neighbor_position[i].neighbors[j].x = 0 ;
//    neighbor_position[i].neighbors[j].y = 0 ;
//    neighbor_position[i].neighbors[j].z = 0 ;
//    }
//  }

}

