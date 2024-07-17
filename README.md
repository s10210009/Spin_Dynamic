The program is to simulate spin dynamic under external magnetic field. The algorithm is followed by this paper, Atomistic spin model simulations of magnetic
nanomaterials, Condens. Matter 26 (2014) 103202. We use RK4 to simulate the change of spin with time.
You could varies the parameter in spin.c
and run the program by 
gcc spin.c -o spin -lm
./spin
python plot.py (or python3 plot.py)

Then, you could see the output figures.
