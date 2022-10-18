import numpy as np
from sl import eliminacao_gaussiana
import time
start_time = time.time()


def int_pol(x, y):
    q = len(x)
    Q = np.zeros((q, q+1))
    for i in range(q):
        Q[i][q] = y[i]
      
    for i in range(q):
        for j in range(q):
            Q[i][j] = x[i]**j
    return Q
def pol_coef(x, y): 
    print(int_pol(x, y), "\n\n\n")
    P = eliminacao_gaussiana(int_pol(x, y))
    for i in range(len(P)):
        print(P[i],"X ^",i)
   
         


