## The program is to simulate spin dynamic under external magnetic field. The algorithm is followed by this paper, Atomistic spin model simulations of magnetic
nanomaterials, Condens. Matter 26 (2014) 103202. 
We use RK4 to simulate the change of spin with time. You could varies the parameter in spin.c

## Run the program by 

### gcc spin.c -o spin -lm

### ./spin
### python plot.py (or python3 plot.py)

Then, you could see the output figures.


The simulation results 

![spin](https://github.com/user-attachments/assets/c201f756-a954-4120-8fba-51581d7ab7ee)



The difference between simulation results and analytic soultion

![error](https://github.com/user-attachments/assets/1180458a-9fd8-453a-acad-43d2a3c73fb9)
