import numpy as np
from sympy import symbols, lambdify
import matplotlib.pyplot as plt

def np_lambdify(varname, func):
    lamb = lambdify(varname, func, modules=['numpy'])
    if func.is_constant():
        return lambda t: np.full_like(t, lamb(t))
    else:
        return lambda t: lamb(np.array(t))

def aproksymacjaDyskretna():
    s = np.array([])
    t = np.array([])
    
    for k in range (0, stopien+1):
        for j in range (0, stopien+1):
            s = np.append(s, np.sum(np.power(xi, j+k))) 
            
    s = s.reshape(stopien+1, stopien+1) # macierz s
    
    for k in range (0, stopien+1):
        t = np.append(t, np.sum(np.multiply(np.power(xi, k), yi))) # wektor t
        
    s1 = np.linalg.inv(s) # macierz odwrotna macierzy s
    a = np.dot(s1,t) # wektor współczynników
    
    g = 0
    
    for i in range (0, len(a)):
        g += a[i]*x**i # funkcja aproksymacyjna g
    
    xx = np.linspace(xi[0]*0.9, xi[-1]*1.1, 1000) 
    yy = lambdify(x, g)(xx)
    
    plt.plot(xx, np.transpose(yy))
    plt.plot(xi, yi, '.')
    plt.show()

    return g

def aproksymacjaDyskretnaGS():
    v = np.array([])
    w = np.array([])
    v = np.append(v, 1)
    w = np.append(w, x**0)
    
    xii = np.array(xi) # przekształcam listę na np.array, żeby móc skorzystać z niej w operacji lambdify
    yii = np.array(yi)
        
    for k in range(2, stopien+2):
        v = np.append(v, x**(k-1))
        w = np.append(w, x**(k-1))
        for j in range(1, k):
            w[k-1] -= np.sum(np.multiply(lambdify(x, v[k-1])(xii), lambdify(x, w[j-1])(xii))) / np.sum(np.power(np_lambdify(x, w[j-1])(xii), 2)) * w[j-1]
    
    l = np.array([]) # iloczyny skalarne (w|g)
    c = np.array([]) # współczynniki przy w
    g = 0

    for i in range(0, stopien+1):
        l = np.append(l, np.sum(np.multiply(lambdify(x, w[i])(xii), yii)))
        c = np.append(c, l[i] / np.sum(np.power(np_lambdify(x, w[i])(xii), 2)))
        g += c[i] * w[i]
    
    xx = np.linspace(xi[0]*0.9, xi[-1]*1.1, 1000)
    yy = lambdify(x, g)(xx)
    
    plt.plot(xx, np.transpose(yy))
    plt.plot(xi, yi, '.')
    plt.show()
    
    return g

if __name__ == '__main__':

    metoda = input("Metoda aproksymacji (Dyskretna/DyskretnaGS):")
    stopien = int(input("Stopien wielomianu aproksymacyjnego: "))
    xi = [float(x) for x in input("Wartości xi: ").split()]
    yi = [float(x) for x in input("Wartości f(xi): ").split()]

    if len(xi) != len(yi):
        print("\nRóżna liczba argumentów i wartości funkcji! Wprowadź dane ponownie.")

    x = symbols('x') # definicja x jako zmiennej

    if metoda == "Dyskretna":
        aproksymacjaDyskretna()
    elif metoda == "DyskretnaGS":
        aproksymacjaDyskretnaGS()

