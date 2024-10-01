import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('result.txt')
data = np.array(data)
data = np.transpose(data)
t1=2
t2=100000000
plt.figure()
plt.plot(data[0,t1:t2],data[1,t1:t2], label='sx')
#plt.plot(data[0,t1:t2],data[4,t1:t2], label='analtic_sx')
plt.legend()
plt.xlabel('t')
plt.savefig('spin_x.png')
plt.tight_layout()

plt.figure()
plt.plot(data[0,t1:t2],data[2,t1:t2], label='sy')
#plt.plot(data[0,t1:t2],data[5,t1:t2], label='analtic_sy')

plt.legend()
plt.xlabel('t')
plt.savefig('spin_y.png')
plt.tight_layout()

plt.figure()
plt.plot(data[0,t1:t2],data[3,t1:t2], label='sz')
#plt.plot(data[0,t1:t2],data[6,t1:t2], label='analtic_sz')

plt.legend()
plt.xlabel('t')
plt.savefig('spin_z.png')
plt.tight_layout()


plt.figure()
plt.plot(data[0,t1:t2],data[1,t1:t2], label='sx')
plt.plot(data[0,t1:t2],data[2,t1:t2], label='sy')
plt.plot(data[0,t1:t2],data[3,t1:t2], label='sz')
#plt.plot(data[0,t1:t2],data[6,t1:t2], label='analtic_sz')

plt.legend()
plt.xlabel('t')
plt.savefig('spin.png')
plt.tight_layout()


#plt.figure()
#plt.plot(data[0,t1:t2],data[1,t1:t2]-data[4,t1:t2], label='d_sx')
#plt.plot(data[0,t1:t2],data[2,t1:t2]-data[5,t1:t2], label='d_sy')
#plt.plot(data[0,t1:t2],data[3,t1:t2]-data[6,t1:t2], label='d_sz')
#plt.legend()
#plt.xlabel('t')
#plt.savefig('error.png')
#plt.tight_layout()
#
