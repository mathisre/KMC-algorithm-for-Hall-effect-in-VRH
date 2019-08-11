from __future__ import print_function
import os
import sys
import numpy as np

try:
	T1 = float(sys.argv[1])
	T2 = float(sys.argv[2])
	dT = float(sys.argv[3])
except IndexError, ValueError:
	print("Put T1 , T2 and dT")
	exit()

Np = int(round((T2-T1)/dT))
T_lin = np.linspace(T1, T2, Np+1)

print("Temperatures:")
print(T_lin)

if len(sys.argv) > 4:
	try: 
		n1 = int(sys.argv[4])
		n2 = int(sys.argv[5])
	except IndexError, ValueError:
		print("Can put n1, n2 in argv")
else:
	n = 10
	n1 = 0
	n2 = n
print("n1:", n1, "n2;", n2)
ctype = 7
triJumps = 1

screen = -1
Cint = 1

L = 100
fill = 0.5
loc = 2
maxJL = 5
directory = "../finalRandomPositions/DOS_ES/"

for i in range(n1,n2):
	for T in T_lin:
		Ex = T/5
		Nt = int(1*10**6)
		Hz = 0.0
		np.random.seed(int(round(Nt+(Ex+5)*(Hz+5)*(T+1)*(i+1)*loc*(screen+10)*(triJumps+1))))
		seed = np.random.randint(10**7)
		seedconfig = np.random.randint(10**7)
		seedstate = np.random.randint(10**7)
		statement = "./glatzprogram outpre=%s%d_ Nt=%d Ex=%.7f Hz=%.7f T=%.7f seed2=%d rsconfig=%d rsstate=%d screen=%.4f fill=%.4f L=%d ctype=%d loc=%f triJumps=%d maxJL=%d Cint=%d" % (directory,i,Nt,Ex,Hz,T,seed,seedconfig,seedstate,screen,fill,L,ctype, loc, triJumps, maxJL, Cint)
		print(statement)
		os.system(statement)
