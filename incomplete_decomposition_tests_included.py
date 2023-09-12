import numpy as np


#Implementation de la factorisation incomplete de Cholesky ci-dessus. La factori#sation est stockee sous la forme d'une matrice triangulaire inferieure, avec le#s elements du triangle superieur mis a zero.(comme convenu ci-contre)
#************************
def incomplete_decomposition(a):
	n = a.shape[0] #on recupere la dimension commune (matrice carre)
	for k in range (n):
		a[k,k] = np.sqrt(a[k,k]);#on commence par faire le meme
                #calcul que precedemment deduit des relations apres developpement du produit matriciel.
		for i in range (k,n):
		    if (a[i,k]!=0):#DIFFERENCE MAJEUR vis a vis de l'algorithme des 2 derniers questions ; la le calcul se fait que si la valeur en cette position precise par iteration (i,k parametrent la position courante) est pas nulle.
		        a[i,k] = a[i,k]/a[k,k]            
		for j in range(k,n):
		    for i in range(j-1,n):
		        if (a[i,j]!=0):#de meme apres avoir fini les premiers calculs comme explique precedemment pour ne pas avoir de dependances entre les calculs.
		            a[i,j] = a[i,j]-a[i,k]*a[j,k]
        for i in range(n):
            for j in range(i,n):
                a[i,j] = 0 #force par effet de bord le caractere garanti de la matrice en trigulaire superieur/inferieur selon a ce qu'on considere T ou bien la transposee ici ca sera t normalement comme avant; l'effet de bord permet de limiter la multiplication d'instances d'objets Array pour effet un control de la complexite.
        return a
#*****************************(opere par effet de bord cest pour cela
#les elements dans la partie haute sont remis a zero).//

#tests rassembles ici:>>       
if __name__=='__main__':
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
    print(incomplete_decomposition(A))

#explication de la difference majeure vis a vis de l'autre:
#En analyse numerique, une factorisation de Cholesky incomplete d'une matrice definie positive symetrique est une approximation de la factorisation de Cholesky.Une factorisation de Cholesky incomplete est souvent utilisee comme preconditionneur pour des algorithmes comme la methode du gradient conjugue.(La partie qui suit celle-ci en tiendra compte de ce lien et de son interet pour fournir une factorisation incomplete suivant cette derniere question de la partie1.)

#LE PROCEDE DE CALCUL PRECIS VA ETRE DEFINI DANS LE RAPPORT.

#Une facon de trouver une telle matrice K(telque KK* sera le preconditionneur dont on parlerait apres) UTILISEE PRINCIPALEMENT ICI est d'utiliser l'algorithme pour trouver la decomposition exacte de Cholesky, sauf que:


#toute entree est mise a zero si l'entree correspondante dans A est egalement nu
#lle. Cela donne une factorisation de Cholesky incomplete qui est aussi creuse q
#ue la matrice A.(et c'est le but pour avoir une approximation plus orientee ver
#s un calcul simple pour le preconditionneur en question qui va etre utilise pro
#chainement dans des algorithmes y compris dans ce projet (orientee resolution d
#'equation de la chaleur+animation..) l'algorithme traitant la methode du gradie
#nt conjugue.

#utilisation de la bibliotheque comme evoque avant:(direct mais ca permet pas de comprendre l'algo sur sa profondeur mais d'utiliser une operation pour une application future ce qui etait pas principalement notre but en se repartissant les taches sur la partie 1)
##____________________________________
#A = np.array([[1,-2j],[2j,5]])

#A->affichage dans la boucle d'interaction:(binaire interpreteur: le programme
#source python2/3)
#array([[ 1.+0.j, -0.-2.j],
 #      [ 0.+2.j,  5.+0.j]])

#L = np.linalg.cholesky(A)

#L->Log: donne la valeur en __stdout de la decomposition de Chlovesly.
##_____________________________________
