                                                                             
import sys
import numpy as np
from int_pol import pol_coef



print("\n\nBem vindo ao programa de interpolação do Leão\n")

modo = int(input("gostaria de colocar os dados manualmente (0) ou upar os dados para o software? (1)"))
if(modo == 0):
    n = int(input("\n\n\nInsira o número de pontos\n\n"))
    x = np.zeros(n)
    y = np.zeros(n)
    for i in range(n):
        x[i] = float(input("insira o valor de x\n"))
        y[i] = float(input("nsira o valor de y\n"))
    print(x)
    print(y)

    pol_coef(x, y)
else:
    x = sys.argv[1]
    y = sys.argv[2]
    
    pol_coef(x, y)