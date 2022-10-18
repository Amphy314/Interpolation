import numpy as np
def solve_spline( x, y, derivatives=[], verbose=False ):
  """
  Resolve o sistema linear das splines para encontrar os coeficientes.

  Recebe:
    x,y         => Vetores com dados
    derivatives => Se recebe alguma coisas, calcula spline fixa, caso contrário, livre.
    verbose     => Mostra passo-a-passo
  Retorna:
    coefficients => Matriz com coeficientes da spline interpolada
                    Linha: intervalos
                    Coluna: potências
                    
                    [ [ a0, b0, c0, d0 ],
                      [ a1, b1, c1, d1 ],
                      ... ]
  """

  n = len(x)  # Número de pontos
  assert len(y) == n, "Vetores de dados com tamanhos diferentes!"



  a = y.copy()
  h = np.zeros(n-1)
  for i in range(n-1):
    h[i] = x[i+1] - x[i] 

  # Matrizes para resolver
  A = np.zeros( (n,n) )
  B = np.zeros( n )


  # Parte comum à Spline Fixa e Natural
  for i in range(1,n-1):
    A[i][i-1] = h[i-1]
    A[i][i]   = 2*( h[i-1] + h[i] )
    A[i][i+1] = h[i]
  
    B[i] = 3.0*(a[i+1] - a[i])/float(h[i]) - 3.0*(a[i] - a[i-1])/float(h[i-1])
  
  

  # Parte que depende do tipo da Spline
  # Tenta usar as derivadas, se não der, faz a spline natural
  try:
    df0 = derivatives[0]
    df1 = derivatives[1]

    A[0][0]   = 2*h[0]
    A[0][1]   = h[0]
    A[-1][-2] = h[-2]
    A[-1][-1] = 2*h[-2]
  
    B[0] = 3.0*(a[1] - a[0])/float(h[0]) - 3*df0
    B[-1] = 3*df1 - 3.0*(a[-1] - a[-2])/float(h[-2])
  
  # Spline natural
  except:
    A[0][0]   = 1.0
    A[-1][-1] = 1.0


  # Encontra os coeficientes
  c = np.linalg.solve( A, B )
  b = np.zeros(n-1)
  d = np.zeros(n-1)
  for i in range(n-1):
    b[i] = (a[i+1]-a[i])/float(h[i]) - h[i]*( 2*c[i] + c[i+1] )/3.0
    d[i] = (c[i+1]-c[i])/( 3.0 * h[i] )
  



  if verbose:
    
    print("Matriz de splines:")
    print()
    for i in range(len(A)):
      for j in range(len(A[i])):
        print(A[i][j],end=" ")
      print("| ",B[i])
    
    print("Vetores:")
    print(" i  |      x      |      y      |       h     |       a     |      b      |      c      |      d      |")
    fmt="%2d  " + "|  %9.6g  "*7
    for i in range(n-1):
      print(fmt%(i,x[i],y[i],h[i],a[i],b[i],c[i],d[i]))
    


  coef = np.zeros( (n-1,4) )
  for i in range(n-1):
    coef[i][0] = a[i]
    coef[i][1] = b[i]
    coef[i][2] = c[i]
    coef[i][3] = d[i]


  return coef





def calc_spline( x, X, coef ):
  """
  Calcula a spline no ponto X
  Os coeficientes estão na matriz coef
  
  Recebe:
    x => Vetor com dados
    X => Ponto (escalar ou vetor) onde calcular a spline
    coef => Coeficientes da spline, conforme calculados pela função solve_spline
  
  Retorna:
    valor ou vetor da spline calculada em X
  """
  
  y = 0
  
  try:
    n = len(X)
    
    y = np.zeros
    for i in range(n):
      y[i] = calc_spline( x, X[i], coef )
  
  except:
    
    # Encontra a posição de X dentro do vetor x
    k = x.searchsorted( X )
  
    if k>0      : k -= 1  # Use a função S do ponto da esquerda
    if k==len(x): k -= 1  # X[i] > x[n]
    
    H = X - x[k]
    ak = coef[k][0]
    bk = coef[k][1]
    ck = coef[k][2]
    dk = coef[k][3]
    y = ak + H*( bk + H*( ck + H*dk ) )
  
  
  return y