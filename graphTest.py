import numpy as np
import matplotlib as plt
import random
import math as math

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
    for item in T:
        print(item.E)
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
                # print(res)
            res+=temp
    return res

def tauxInjections(F,G):
    return(nbInjectionsGraphes(F,G)/(math.factorial(len(G.V)/math.factorial(len(G.V)-len(F.V)))))

def simulation(W,n):
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




##############################################################################
############################    MAIN    #####################################
##############################################################################

W = GraphonExponentiel(id)
V = [1,2]
E = [(1,2),(2,3)]
E2 = []

F = Graphe(V,E2)
# print(F(2,3))
F.add((1,2))
# print(F(1,3))
# print (W(0,0))
# generateurDeGraphesTests(4)
# print(base_n(12,16,3))

V3 = [1,2,3]
E3 = [(1,2),(1,3)]
G = Graphe(V3,E3)
print(nbInjectionsGraphes(F,G))
print("Taux d'injections : ",tauxInjections(F,G)*100,"%")

U = simulation(W,10)
print(U.V)
print(U.E)


#Lorsque n tend vers l'infini, le teux d'injection de F dans Gn tend vers un bail 
#qui est le taux d'injection de F dans le graphon W
