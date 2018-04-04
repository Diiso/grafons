import numpy as np
import matplotlib.pyplot as plt
import random
import math as math
from scipy import integrate


## Classes

class Graphon:
    def __init__(self, W):
        self.W = W
    def __call__(self, x, y):
        return self.W(x,y)

class GraphonExponentiel(Graphon):
    def __init__(self, g):
        def W(x,y):
            return (math.exp(g(x)+g(y)))/(1+math.exp(g(x)+g(y)))
        super().__init__(W)

class Graphe:
    def __init__(self, V, E):
        self.V = V
        self.E = E
    def add(self, e):
        self.E.append(e)
    def __call__(self, i, j):
        return (((i,j) in self.E) or ((j,i) in self.E))


## Fonctions diverses

def id(x):
    return x

def base_n(n, k, p):
    k1 = k
    T = []
    for i in range(p):
        T.append(k1//(n**(p-i-1)))
        k1 = k1%(n**(p-i-1))
    return T

def diffTableau(T):
    for i in range(len(T)):
        for j in range(len(T)):
            if i!=j and T[i]==T[j]:
                return False
    return True


## Méthodes d'intégration

#Méthode des rectangles
#Méthode des trapèzes
#Méthode de Simpson

#Utilisation de scipy : integrate.quad(f,a,b), integrate.romberg(f,a,b) + méthode des trapèzes toute faite

#Méthode de Monte-Carlo en dimension 1
def MC1(f,a,b,N):
    res = 0
    for i in range(N):
        res+=f(random.uniform(a,b))
    return res*(b-a)/N

#Méhode de Monte-Carlo bis (en dimension p)
def MCp(f,a,b,p,N=1000,var=False):  #integrale de f de a à b, en dimension p, par la méthode de Monte Carlo avec N pts
    X = []
    for i in range(p):
        X.append(a+(b-a)*np.random.random_sample(N))
    d = 1.0/((b-a)**p)
    somme = [0.0]*2

    for i in range(N):
        x = []
        for k in range(p):
            x.append(X[k][i])

        val = f(x)
        somme[0] += val
        somme[1] += val*val

    moyenne = somme[0]/(d*N)
    variance = somme[1]/(d*d*N)-moyenne*moyenne

    if var:
        return (moyenne,math.sqrt(variance/N)*1.96)
    return moyenne


## Fonctions pour la simulation

def generateurDeGraphesTests(k): #k est le nombre de noeuds
    T = [] #Notre tableau de graphes
    Etot = [] #Tableau d'aretes possibles
    V = [] #Contient les sommets, commun a tous les graphes

    for i in range(1,k+1):
        V.append(i)
    for i in range(1,k+1):
        for j in range(i+1,k+1):
            Etot.append((i,j))

    K = int(k**2-(k*(k+1)/2))
    for i in range(2**K):
        F = Graphe(V,[])
        b = bin(i+2**K)
        for j in range(K):
            if (b[len(b)-1-j] == "1"):
                F.add(Etot[j])
        T.append(F)

    return T

def nbInjectionsGraphes(F,G):
    assert(len(F.V)<=len(G.V))

    p = len(F.V)
    n = len(G.V)

    res = 0

    for i in range(n**p):
        Bi = base_n(n,i,p)
        if diffTableau(Bi):
            temp = 1
            for e in F.E:
                temp*=int(G(Bi[e[0]-1],Bi[e[1]-1]))
            res+=temp
    return res

def tauxInjectionsGraphes(F,G):
    assert (len(G.V)>=len(F.V))

    temp = math.factorial(len(G.V))/math.factorial(len(G.V)-len(F.V))
    return nbInjectionsGraphes(F,G)/temp

def tauxInjectionsGraphons(F,W):

    def f(x):
        res = 1
        for e in F.E:
            res *= W(x[e[0]-1],x[e[1]-1])
        return res

    return MCp(f,0,1,len(F.V))

def erreur(F,G,W,epsilon):
    res = True
    for i in range (len(F)):
        if (abs(tauxInjectionsGraphes(F[i],G)-tauxInjectionsGraphons(F[i],W))>epsilon):
            res = False
    return res

def simulation(W,n):    #renvoie G_n (de la suite qui tend vers W)
    V = []
    X = []

    for i in range(1,n+1):
        V.append(i)
        X.append(random.uniform(0,1))

    E = []
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            p = random.uniform(0,1)
            if p<=W(X[i-1],X[j-1]):
                E.append((i,j))

    return Graphe(V,E)

def CVps_temp(f,k,epsilon):

    W = GraphonExponentiel(f)
    F = generateurDeGraphesTests(k)

    res = False
    n = 10
    G = simulation(W,n)
    while ( n<1000 and not(res) and not(erreur(F,G,W,epsilon)) ):
        n += 1
        G = simulation(W,n)
    return erreur(F,G,W,epsilon)

def CVps(f,k,epsilon,N=1000):

    res = True

    for k in range(N):
        res = res and CVps_temp(f,k,epsilon)

    return res

def CVl_temp(f,k,n,N=1000):

    W = GraphonExponentiel(f)
    F = generateurDeGraphesTests(k)

    esp = []
    var = []
    X = []

    for i in range(len(F)):
        esp.append(0)
        var.append(0)
        X.append([])

    for j in range(N) :
        G = simulation(W,n)

        for i in range (len(F)):
            print(i,j)
            X[i].append(math.sqrt(n)*(tauxInjectionsGraphes(F[i],G)-tauxInjectionsGraphons(F[i],W)))
            esp[i] += X[i][j]

    for i in range(k):
        esp[i] *= 1/N

    for j in range(N):
        print(j)
        var[i] += abs(X[i][j]-esp[i])**2

    for i in range(k):
        var[i] *= 1/(N**2)


    plt.hist(X[1],500,normed=1)
    plt.show()

    #print(X)


    return esp,var




##############################################################################
############################    MAIN    #####################################
##############################################################################


print(CVps(id,2,5,1))
print(CVl_temp(id,3,50))



## Tests

"""
W = GraphonExponentiel(id)
V = [1,2]
E = [(1,2),(2,3)]
E2 = []

F = Graphe(V,E2)
print(F(2,3))
F.add((1,2))
print(F(1,3))
print (W(0,0))
generateurDeGraphesTests(4)
print(base_n(12,16,3))

V3 = [1,2,3]
E3 = [(1,2),(1,3)]
G = Graphe(V3,E3)
print(nbInjectionsGraphes(F,G))
print("Taux d'injections : ",tauxInjections(F,G)*100,"%")

U = simulation(W,10)
print(U.V)
print(U.E)

def f1(x):
    res = 0
    for i in range (len(x)):
        res += x[i]**2
    return res

def f2(x):
    res = 1
    for i in range(len(x)):
        res *= x[i]**2
    return res

print(MCp(f1,0,1,2))

"""
