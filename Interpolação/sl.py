#! /usr/bin/python3
# -*- coding: utf-8 -*-
#
# Algoritmos para encontrar soluções de Sistemas de Equações Lineares
#_________________________________________________________
# Universidade Federal de Santa Catarina
# Departamento de Engenharias da Mobilidade
# Curso de Cálculo Numérico
# Prof. Alexandre Zabot
# https://zabot.paginas.ufsc.br
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import numpy as np
 
 

def eliminacao_gaussiana( A, verbose=False ):
  """
  Algoritmo da Eliminação Gaussiana
  
  Recebe:
    A       => Matriz aumentada de coeficientes
    verbose => Imprime passo a passo?
    
  Retorna:
    x  => Vetor solução
  """

  n = len(A)



  if verbose: print(A)
  
  
  # Faço o escalonamento para cada uma das colunas
  for c in range(n-1):
    if verbose:
      print ("\n\n----------------------------"       )
      print ("Eliminação Gaussiana na coluna %d." % c , end = '')
  
  
  
    # Procuro um Pivô
    p=c
    while A[p,c]==0 and p<n-1:
      p+=1
    if A[p,c]==0:
      exit("O Pivô é nulo, o algoritmo falha!")
  
  
    # Se for necessário, troca linhas      
    if p==c:
      if verbose:
        print("Pivô A[%d,%d] = %.6g, não preciso trocar as linhas:" % (p,c,A[p,c]))
    else:  # O pivô não está na linha atual, faço uma troca de linhas
      if verbose: 
        print("Pivô A[%d,%d] = %.6g, trocando as linhas %d <-> %d" % (p,c,A[p,c],p,c) )
      x = A[c].copy()
      A[c] = A[p]
      A[p] = x
      if verbose: print(A)
    
    
    # Faz o escalonamento para a coluna c
    if verbose: print( "" )
    for l in range(c+1,n):
      m = A[l,c]/A[c,c]
      if verbose: print("E_%d - %f E_%d ->  E_%d " % (l,m,c,l))
      A[l] = A[l] - m*A[c]
    if verbose: print(A)
  
  
  # Não posso iniciar a substituição regressiva
  if A[n-1,n-1]==0:
    exit("A[%d,%d] = 0, não existe solução única!" % (n-1,n-1) )
  
  
  # Substituição regressiva
  if verbose: 
    print("\n\n----------------------------")
    print("Substituição regressiva"         )
  x = [0]*n
  x[n-1] = A[n-1,n]/A[n-1,n-1]
  if verbose:
    print( "x_%d = %.6g / %.6g = %.6g" % (n-1, A[n-1,n], A[n-1,n-1], x[n-1]) )
  
  for i in range(n-2,-1,-1):
    s = 0
    strsoma = ""
    if verbose: print( "x_%d = " % (i), end = '')
    for j in range(i+1,n):
      s += A[i,j]*x[j]
      if verbose: strsoma += "%.6g x_%d + " % (A[i,j],j)
    x[i] = (A[i,n] - s)/A[i,i]
    if verbose:
      print ( "%.6g - (%s) )/ %.6g =  %.6g" % (A[i,n],strsoma,A[i,i],x[i]) )
  if verbose: print( x )
  


  return x





def print_tridiagonal( a,b,c,d ):
  """
  Imprimie uma matriz tridiagonal de modo formatado
  
  Recebe:
    Vetores da matriz, conforme definido na aula
    
  Retorna:
    Nada, apenas imprime saída
  """
  n = len(a)
  
  print (" %+11.6f %+11.6f" % (b[0],c[0]), end = '')
  print (" "*12*(n-2),"| %+11.6f" % d[0]  )
  
  for i in range(1,n-1):
    print (" "*12*(i-1), end = '')
    print ("%+11.6f %+11.6f %+11.6f" % (a[i],b[i],c[i]), end = '')
    print (" "*(12*(n-2-i)),"| %+11.6f" % d[i] )

  print (' '*12*(n-2), end = '')
  print (" %+11.6f %+11.6f | %+11.6f" % (a[-1],b[-1],d[-1]) )
  
  
  

def thomas( A, verbose=False ):
  """
  Algoritmo de Thomas
  https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  
  Recebe:
    A       => Matriz aumentada de coeficientes
    verbose => Imprime passo a passo?
    
  Retorna:
    x  => Vetor solução
  """  
  n = len(A)
  
  a = np.zeros(n)
  b = np.zeros(n)
  c = np.zeros(n)
  d = np.zeros(n)

  b[0] = A[0][0]
  c[0] = A[0][1]
  d[0] = A[0][-1]
  for i in range(1,n-1):
    a[i] = A[i][i-1]
    b[i] = A[i][i]
    c[i] = A[i][i+1]
    d[i] = A[i][-1]
  a[-1] = A[-1][-3]
  b[-1] = A[-1][-2]
  d[-1] = A[-1][-1]


  if verbose:
    print("Matriz original")
    print_tridiagonal( a,b,c,d )


  den = float( b[0] )
  c[0] /= den
  d[0] /= den
  for i in range(1,n):
    den = float( b[i] - a[i]*c[i-1] )
    c[i] /= den
    d[i]  = (d[i] - a[i]*d[i-1])/ den

  if verbose:
    print("Depois do passo progressivo")
    print_tridiagonal( a,b,c,d )
    
  x = np.copy(d)
  for i in range(n-2,-1,-1):
    x[i] -= c[i]*x[i+1]


  return x  








def decompoeLU( A, verbose=True ):
  """
  Decomposição LU
  Atenção: sobreescreve a matriz A
  
  Recebe:
    A       => Matriz de coeficientes
    verbose => Imprime passo a passo?
    
  Retorna:
    nada, a matriz de entrada é sobreescrita
  """
  n = len(A)
  
  
  if verbose: print(A)
  
  for j in range(n):  # Coluna
    for i in range(n): # linha
      
      # Parte L da matriz
      if i>=j:
        s = 0
        for k in range(0,j):
          s += A[i][k] * A[k][j]
        A[i][j] -= s

      
      # Parte U da matriz
      else:
        s = 0
        for k in range(0,i):
          s += A[i][k] * A[k][j]        
        
        A[i][j] = ( A[i][j] - s )/float(A[i][i])


  if verbose: print(A)
  



def resolveLU( LU, B, verbose=True ):
  """
  Resolve um sistema já decomposto em LU
  Atenção: sobreescreve a matriz A
  
  Recebe:
    LU      => Matriz LU decomposta
    B       => Vetor com termos independentes do Sistema Linear
    verbose => Imprime passo a passo?
    
  Retorna:
    x  => Vetor solução
  """
  n = len(LU)
  assert len(B)==n, "Matrizes de tamanhos diferentes"
  
  # Resolve Ly=b
  y = np.zeros(n,dtype=float)
  for i in range(n):
    s=0
    for j in range(0,i):
      s += LU[i][j]*y[j]
    y[i] = ( B[i] - s )/float(LU[i][i])
  
  if verbose: print("y = ",y)
  
  # Resolve Ux=y
  x = np.zeros(n,dtype=float)
  for i in range(n-1,-1,-1):
    s=0
    for j in range(i+1,n):
      s += LU[i][j]*x[j]
    x[i] = y[i] - s

  return x


def lu( A, B, verbose=True ):
  """
  Resolve um sistema pela decomposição LU
  Atenção: sobreescreve a matriz A
  
  Recebe:
    A       => Matriz de coeficientes
    B       => Vetor com termos independentes do Sistema Linear
    verbose => Imprime passo a passo?
    
  Retorna:
    x  => Vetor solução
  """
  decompoeLU( A, verbose )
  return resolveLU( A, B, verbose )
  









def gauss_seidel( A, X, TOL, Nmax, verbose=False ):
  """
  Resolve um sistema pelo Método de Gauss-Seidel
  
  Recebe:
    A       => Matriz aumentada de coeficientes
    X       => Chute inicial
    TOL     => Tolerância a se alcançar (b-a)/2 < tol
    Nmax    => Número máximo de iterações    
    verbose => Imprime passo a passo?
    
  Retorna:
    x  => Vetor solução
  """
  n = len(A)
  
  if verbose: print(A)
  
  
  # Cópias que guardam os valores dos passos atuais e anteriores
  x0 = X.copy()
  
  
  # Aplica o algoritmo
  convergiu=False
  k=0
  
  if verbose:
    print("  n |",end="")
    for i in range(n):
      print("        x_%d        |" % i,end="")
    print(" Erro absoluto |")

    print("  0 |",end="")
    for xi in X:
      print(" %+.10e |" % xi,end="")
    print("    -----      |")
    
    
  while k<Nmax and not convergiu:
      
    # Atualiza os termos
    for i in range(n):
      
      # Calcula as somas
      s1 = 0
      s2 = 0
      for j in range(n):
        if   j<i: s1 += A[i][j]*X[j]   # Usa os x já atualizados
        elif j>i: s2 += A[i][j]*x0[j]  # Usa os x antigos
      
      # Calcula o novo termo
      X[i] = ( -s1 - s2 + A[i][n] )/float(A[i][i])
        
    
    # Verifica se convergiu
    diff2 = (X-x0)*(X-x0)
    erro_abs = np.sqrt( diff2.sum() )
    if erro_abs < TOL:
      convergiu=True
    
    x0=X.copy()
    k+=1
    
    if verbose:
      print(" %2d |" % k,end="")
      for xi in X:
        print(" %+.10e |" % xi,end="")
      print("    %.1e    |" % erro_abs)
  
  if verbose:
    print("")
    if convergiu:
      print("Convergiu com %d iterações" % k )
    else:
      print("NÃO convergiu depois de %d iterações" % k)
  
  return X

