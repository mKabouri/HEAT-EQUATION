import numpy as np
import math

def somme(debut, fin, M):
    resultat=0
    for i in range (debut,fin):#fin<-i
        resultat+=M[fin-1,i]
    return resultat

def somme_produits(p,debut, fin, M):
    resultat2=0
    for i in range(debut,fin):
        resultat2+=M[p,i]*M[i,fin-1]
    return resultat2

def somme_produits2(debut, fin, M):
    resultat2=0
    for i in range(debut,fin):
        resultat2+=M[fin,i]*M[fin,i]
    return resultat2

def decomposition(A,n):
    T=np.zeros((n,n))
    T[0,0]=math.sqrt(A[0,0])
    for i in range (1,n):
        T[i,0]=A[0,i]/T[0,0]
    for k in range (1,n-1):
        T[k,k]=math.sqrt(A[k,k]-somme(1,k,T))
        for i in range (k,n):
            T[i,k]=(A[i,k]-somme_produits(i,1,k,T))/T[k,k]
    return T
                            
def test():
    A=np.zeros((3,3))
    for i in range(3):
        A[i,i]=2
    A[0,2]=0
    A[2,0]=0
    for i in range(3):
        for j in range(3):
            if(i!=j and (i>0 or j>0)):
               A[i,j]=-1
    B=decomposition(A,3)
    return B

def test2():
    A=np.zeros((4,4))
    for j in range(4):
        A[0,j]=1
    for j in range(4):
        A[j,0]=1
    for j in range(1,4):
        A[1,j]=5
    for j in  range(1,4):
        A[j,1]=5;
    k=2
    for h in range(2,4):
        A[k,h]=14
    for h in range(2,4):
        A[h,k]=14
    A[3,3]=15
    return decomposition(A,4)



print(test2())










