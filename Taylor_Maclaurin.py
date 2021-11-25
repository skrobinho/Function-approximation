import math
import numpy as np
from sympy import symbols, lambdify, sympify, diff, evalf, expand
from sympy.plotting import plot

def szeregTaylora():
    w = f.subs(x,x0)
    fx = f
    
    for k in range(1, stopien+1):
        fx = diff(fx, x) 
        w += fx.subs(x, x0)/math.factorial(k)*(x-x0)**k
    
    w = expand(w.evalf())
    
    if ylim2>ylim1:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False)
        p2 = plot(w, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False, line_color='red')
    else:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False)
        p2 = plot(w, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False, line_color='red')
    
    p1.extend(p2)
    p1.show()
    
    return w

def szeregMaclaurina():
    x0 = 0
    w = f.subs(x,x0)
    fx = f
    
    
    for k in range(1, stopien+1):
        fx = diff(fx, x) 
        w += fx.subs(x, x0)/math.factorial(k)*(x-x0)**k
    
    w = expand(w.evalf())
    
    if ylim2>ylim1:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False)
        p2 = plot(w, xlim=(a,b), ylim=(float(ylim1), float(ylim2)), axis_center=(0,0), show=False, line_color='red')
    else:
        p1 = plot(f, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False)
        p2 = plot(w, xlim=(a,b), ylim=(float(ylim2), float(ylim1)), axis_center=(0,0), show=False, line_color='red')
    
    p1.extend(p2)
    p1.show()
        
    return w

if __name__ == '__main__':

    metoda = input("Metoda aproksymacji (Taylor/Maclaurin):")
    stopien = int(input("Stopien wielomianu aproksymacyjnego: "))
    a = sympify(input("Granica dolna przedzialu: "))
    b = sympify(input("Granica gorna przedzialu: "))
    f = sympify(input("Funkcja pierwotna: "))
    x0 = sympify(input("Punkt x0: "))

    x, e, pi = symbols('x e pi') # definicja x i e jako zmiennych
    f = f.subs(e, math.e).subs(pi, math.pi)
    a = float(a.subs(e, math.e).subs(pi, math.pi))
    b = float(b.subs(e, math.e).subs(pi, math.pi))
    x0 = x0.subs(e, math.e).subs(pi, math.pi)

    ylim1 = f.subs(x,a) # wartość funkcji f w punkcie a
    ylim2 = f.subs(x,b)

    if ylim2>ylim1:
        ylim1 = 0.9*ylim1
        ylim2 = 1.1*ylim2
    else:
        ylim2 = 0.9*ylim2
        ylim1 = 1.1*ylim1

    if metoda == "Taylor":
        szeregTaylora()
    elif metoda == "Maclaurin":
        szeregMaclaurina()