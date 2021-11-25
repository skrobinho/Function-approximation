import math
import numpy as np
from sympy import symbols, integrate, sympify
from sympy.plotting import plot

def aproksymacjaCiagla():
    s = np.array([])
    t = np.array([])
    
    for k in range (0, stopien+1):
        for j in range (0, stopien+1):
            s = np.append(s, integrate(x**(j+k), (x,a,b))) 
    
    s = s.reshape(stopien+1, stopien+1).astype(np.float32) # macierz s
        
    for k in range (0, stopien+1):
        t = np.append(t, integrate(f*x**k, (x,a,b))) # wektor t
    
    s1 = np.linalg.inv(s) # macierz odwrotna macierzy s
    c = np.dot(s1,t).astype(np.float32) # wektor współczynników
    
    g = 0
    
    for i in range (0, len(c)):
        g += c[i]*x**i # funkcja aproksymacyjna g
    
    # kreślenie wykresów funkcji f i g na podanym przedziale
    if ylim2>ylim1:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False)
        p2 = plot(g, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False, line_color='red')
    else:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False)
        p2 = plot(g, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False, line_color='red')
    
    p1.extend(p2)
    p1.show()

    return g

def aproksymacjaCiaglaGS():
    v = {}
    w = {}
    v["v1"] = 1
    w["w1"] = 1

    for k in range(2, stopien+2):
        v["v{0}".format(k)] = x**(k-1)
        w["w{0}".format(k)] = x**(k-1)
        for j in range(1, k):
            w["w{0}".format(k)] -= integrate((v['v{0}'.format(k)]*w['w{0}'.format(j)]),(x,a,b)) / integrate((w['w{0}'.format(j)]**2),(x,a,b)) * w['w{0}'.format(j)]

    l = {} # iloczyny skalarne (w|g)
    c = {} # współczynniki przy w
    g = 0

    for i in range(1, stopien+2):
        l["l{0}".format(i)] = integrate((f*w['w{0}'.format(i)]),(x,a,b))
        c["c{0}".format(i)] = l["l{0}".format(i)] / integrate((w['w{0}'.format(i)]**2),(x,a,b))
        g += c["c{0}".format(i)] * w["w{0}".format(i)]

    if ylim2>ylim1:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False)
        p2 = plot(g, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False, line_color='red')
    else:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False)
        p2 = plot(g, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False, line_color='red')
    
    p1.extend(p2)
    p1.show()
    
    return g

if __name__ == '__main__':
    
    metoda = input("Metoda aproksymacji (Ciagla/CiaglaGS):")
    stopien = int(input("Stopien wielomianu aproksymacyjnego: "))
    a = sympify(input("Granica dolna przedzialu: "))
    b = sympify(input("Granica gorna przedzialu: "))
    f = sympify(input("Funkcja pierwotna: "))

    x, e, pi = symbols('x e pi') # definicja x i e jako zmiennych
    f = f.subs(e, math.e).subs(pi, math.pi) # przypisanie zmiennej e wartości e
    a = float(a.subs(e, math.e).subs(pi, math.pi))
    b = float(b.subs(e, math.e).subs(pi, math.pi))

    ylim1 = f.subs(x,a) # wartość funkcji f w punkcie a
    ylim2 = f.subs(x,b)

    if ylim2>ylim1:
        ylim1 = 0.9*ylim1
        ylim2 = 1.1*ylim2
    else:
        ylim2 = 0.9*ylim2
        ylim1 = 1.1*ylim1

    if metoda == "Ciagla":
        aproksymacjaCiagla()
    elif metoda == "CiaglaGS":
        aproksymacjaCiaglaGS()
