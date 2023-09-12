from string import ascii_letters
import numpy as np
from matplotlib import pyplot as plt
from sqlalchemy import desc
np.set_printoptions(linewidth=300, threshold=np.inf)

#####################################################################################################################################
# CHOLESKY
#####################################################################################################################################

def squared_sum(T, i):
	r = T[i][0]**2
	for k in range(1, i):
		r += T[i][k]**2
	return r

def prod_sum(T, i, j):
	r = T[i][0]*T[j][0]
	for k in range(1, i):
		r += T[i][k]*T[j][k]
	return r

# DENSE VERSION OF CHOLESKY FACTORISATION
def cholesky(A):
	n = A.shape[0]
	T = np.zeros((n, n))

	#* Calcule de la 1ere colonne
	T[0][0] = np.sqrt(A[0][0])
	for k in range(1, n):
		T[k][0] = A[0][k]/T[0][0]

	#* Calcule des colonnes restantes
	for i in range(1, n):
		T[i][i] = np.sqrt(A[i][i] - squared_sum(T, i))
		for j in range(i+1, n):
			T[j][i] = (A[i][j] - prod_sum(T, i, j))/T[i][i]
	return T

def cholesky_incomplete(A):
	n = A.shape[0]
	T = np.zeros((n, n))

	#* Calcule de la 1ere colonne
	T[0][0] = np.sqrt(A[0][0])
	for k in range(1, n):
		T[k][0] = A[0][k]/T[0][0]

	#* Calcule des colonnes restantes
	for i in range(1, n):
		T[i][i] = np.sqrt(A[i][i] - squared_sum(T, i))
		for j in range(i+1, n):
			if (A[i][j]!=0): #* pas de calcul pour les cellules nulles
				T[j][i] = (A[i][j] - prod_sum(T, i, j))/T[i][i]
	return T

#####################################################################################################################################
# HEAT EQUATIONS AND GRADIENT METHODS
#####################################################################################################################################

# Algorithm used to generate the tridiagonal matrix used in the heat equation
def heat(N):
	A = np.zeros((N*N, N*N), int)
	b = np.zeros((N*N, 1), int)
	for i in range(N*N):
		for j in range(N*N):
			dif = abs(i-j)
			#* Second condition is used to have the 0s at the edge of each block through the diagonal
			if ((abs(dif-N)==0 or dif == 1) and not((j%N == 0 and i == j - 1) or (i%N == 0 and j == i-1))):
				A[i][j] = 1

	np.fill_diagonal(A, -4)
	return [A, b]

# def solve_heat_eq(A, b, N):
# 	return  np.linalg.solve(((N+1)**2)*A, b)

def descent(L, b):
	n = b.shape[0]
	y = np.zeros((n, 1))
	y[0] = b[0]/L[0][0]
	for i in range(1, n):
		y[i] = b[i]
		for j in range(i):
			y[i] -= L[i][j]*y[j]
		y[i] /= L[i][i]
	return y

def ascent(U, y):
	n = y.shape[0]
	x = np.zeros((n, 1))
	x[n-1] = y[n-1]/U[n-1][n-1]
	for i in range(n-2, -1, -1):
		x[i] = y[i]
		for j in range(i+1, n):
			x[i] -= U[i][j]*x[j]
		x[i] /= U[i][i]
	return x

# Conjugate gradient method without the preconditioning
def conjgrad(A, b, x0):
	r = b-np.dot(A, x0)
	p = r
	rsold = np.dot(np.transpose(r), r)
	for i in range(1000000):
		Ap = np.dot(A, p)
		alpha = rsold/(np.dot(np.transpose(p), Ap))
		x0 = x0 + alpha*p;
		r = r - alpha*Ap;
		rsnew = np.dot(np.transpose(r), r)
		if np.sqrt(rsnew) < 10**(-10):
			break
		p = r + rsnew/rsold*p;
		rsold = rsnew
	return x0

# Conjugate gradient method using the preconditioning
def precond_conjgrad(A, b, x0):
	#* Initialisation
	L = cholesky_incomplete(A)
	r = b - np.dot(A, x0)
	#? resolution du systeme
	z = ascent(np.transpose(L), descent(L, r))
	p = z
	rzold = np.dot(np.transpose(r), z)

	#* boucle
	for i in range(1000000):
		Ap = np.dot(A, p)
		alpha = rzold/(np.dot(np.transpose(p), Ap))
		x0 = x0 + alpha*p
		r = r - alpha*Ap
		z = ascent(np.transpose(L), descent(L, r))
		rznew = np.dot(np.transpose(r), z)
		if np.sqrt(rznew) < 10**(-10):
			break
		p = z + rznew/rzold*p
		rzold = rznew
	return x0

#####################################################################################################################################
# TESTS OF THE FUNCTIONS
#####################################################################################################################################

def test_cholesky_dense():
	A = np.zeros((4, 4))
	A[-2,:] = 14
	A[:,-2] = 14
	A[:,1] = 5
	A[1,:] = 5
	A[-1, -1] = 15
	A[0,:] = 1
	A[:,0] = 1

	L = cholesky(A)
	print("Test de la fonction dense de cholesky")
	print("A :\n", A)
	print("L :\n", L)
	print("L.L' == A?:\n", np.dot(L, np.transpose(L))==A)

def test_conjgrad(N):
	print("Test of the conjugate gradient method")
	A,b = heat(N)
	b[:N] = -1
	x = conjgrad(A, b, np.zeros((A.shape[0], 1)))
	x_real = np.linalg.solve(A, b)
	print("Conjugate gradient method result is ok : ", x - x_real < 10**(-10))

def test_precond_conjgrad(N):
	print("Test of the preconditioned conjugate gradient method")
	A,b = heat(N)
	b[:N] = -1
	x = precond_conjgrad(-A, -b, np.zeros((A.shape[0], 1)))
	x_real = np.linalg.solve(A, b)
	print("Conjugate gradient method is ok : ", x - x_real < 10**(-10))

def radiator(N):
	A,b = heat(N)
	b.shape = (N, N)
	b[N//2][N//2] = -1
	b.shape = (N*N, 1)
	x = precond_conjgrad(-((N+1)**2)*A, -b, np.zeros((A.shape[0], 1)))
	plt.figure(1)
	plt.imshow(x.reshape(N, N), cmap="hot")

def heated_northern_wall(N):
	A,b = heat(N)
	b[:N] = -1
	x = precond_conjgrad(-((N+1)**2)*A, -b, np.zeros((A.shape[0], 1)))
	plt.figure(2)
	plt.imshow(x.reshape(N, N), cmap="hot")

#####################################################################################################################################
# MAIN FUNCTION
#####################################################################################################################################

if __name__=='__main__':
	test_cholesky_dense()
	test_conjgrad(4)
	test_precond_conjgrad(4)
	radiator(10)
	heated_northern_wall(10)
	plt.show()