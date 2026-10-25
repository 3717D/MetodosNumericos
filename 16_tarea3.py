import sympy as sp 
import numpy as np 
import matplotlib.pyplot as plt

# --- Valores generales ---
x = sp.symbols('x')
f = 2*(x**3) - 6*(x**2) + 8*x - 5
fprima = sp.diff(f, x, simplify=True)
itmax = 100
eamin = 0.0001      # Porcentaje

# ------
xreal = 1.423853799 # VALOR OBTENIDO MEDIANTE CARDANO
# ------

fev = f.subs(x, xreal).evalf()
print(fev)
# >> -2.14784634522402e-10 
# = 0.000000000214784634522402 (cantidad muy pequeña, error dado por
# usar una calculadora muy básica)

# --- Valores bisección y regla falsa ---
xl = 1
xu = 2

# --- Valores Newton-Raphson ---
x0nr = 2

# --- Valores secante ---
x0s = 1
x1s = 2

def biseccion(f, xl, xu, iteracionesMax):
    listabi = []

    for i in range(iteracionesMax):
        xr = (xl + xu)/2
        fxr = f.subs(x,xr).evalf()
        fxl = f.subs(x,xl).evalf()

        if fxr * fxl < 0:
            xrant = xu
            xu = xr
        else: 
            xrant = xl
            xl = xr

        if i > 0: ea = abs((xr - xrant)/xr)*100
        else: ea = 100

        et = abs((xreal - xr)/xreal) * 100

        # No entiendo por qué la matriz sale en notación científica
        listabi.append((i, xr, ea, et))
        if ea < eamin: break

    return xr, listabi, i

def reglaFalsa(f, xl, xu, iteracionesMax):
    listarf = []

    for i in range(iteracionesMax):
        
        fxl = f.subs(x, xl).evalf()
        fxu = f.subs(x, xu).evalf()

        xr = xu + ((fxu * (xl - xu)) / (fxu - fxl))
        # Esta fórmula me salió por accidente, y es diferente a la vista en
        # clase. Lo tengo anotado en el cuaderno :P
        fxr = f.subs(x,xr).evalf()

        if fxr * fxu < 0: xl = xr 
        else: xu = xr

        if i > 0: ea = abs((xr - xrant)/xr)*100
        else: ea = 100

        et = abs((xreal - xr)/xreal) * 100

        listarf.append((i, xr, ea, et))
        xrant = xr

        if ea < eamin: break

    return xr, listarf, i

def newtonRaphson(f, fp, x0, iteracionesMax):
    listaNR = []

    for i in range(iteracionesMax):
        fx0 = f.subs(x,x0).evalf()
        fprimax0 = fprima.subs(x,x0).evalf()
        xi = x0 - (fx0 / fprimax0)

        if i > 0: ea = abs((xi - x0)/xi)*100
        else: ea = 100

        et = abs((xreal - xi)/xreal) * 100

        listaNR.append((i, xi, ea, et))

        x0 = xi

        if ea < eamin: break
    
    return x0, listaNR, i

def secante(f, x0, x1, iteracionesMax):
    listaSE = []
    
    for i in range(iteracionesMax):
        fx0 = f.subs(x, x0).evalf()
        fx1 = f.subs(x, x1).evalf()

        # xip1 = x_(i+1)
        x2 = x1 - ((fx1)*(x1 - x0) / (fx1 - fx0))

        if i > 0: ea = abs((x2 - x1)/x2)*100
        else: ea = 100

        et = abs((xreal - x2)/xreal) * 100

        listaSE.append((i, x2, ea, et))

        x0 = x1
        x1 = x2

        if ea < eamin: break
    
    return x2, listaSE, i

# --- Tabla bis ---
print("BISECCION")
xiBI, listaBI, iBI = biseccion(f, xl, xu, itmax)
res_BI = np.array(listaBI)
print(res_BI)

# --- Tabla ref ---
print(f"\nREGLA FALSA")
xiRF, listaRF, iRF = reglaFalsa(f, xl, xu, itmax)
res_RF = np.array(listaRF)
print(res_RF)

# --- Tabla NR ---
print(f"\nNEWTON-RAPHSON")
xiNR, listaNR, iNR = newtonRaphson(f, fprima, x0nr, itmax)
res_NR = np.array(listaNR)
print(res_NR)

# --- Tabla sec ---
print(f"\nSECANTE")
xiSE, listaSE, iSE = secante(f, x0s, x1s, itmax)
res_SE = np.array(listaSE)
print(res_SE)

fn = sp.lambdify(x, f, 'numpy')
xn = np.linspace(1.4, 1.5, 100)
# Tiene DEMASIADO zoom para que se aprecie mejor la separación

plt.figure()
plt.subplot(211)
plt.plot(xn, fn(xn), 'orange', label='f(x)', alpha=0.5)
plt.plot(xreal, fn(xreal), 'cx', label='Valor real', alpha=0.7)
plt.axhline(0, color='k')

plt.plot(xiBI, fn(xiBI), 'ro', label='Bisección',alpha=0.5)
plt.plot(xiRF, fn(xiRF), 'gd', label='Regla Falsa',alpha=0.5)
plt.plot(xiNR, fn(xiNR), 'bv', label='Newton-Raphson',alpha=0.5)
plt.plot(xiSE, fn(xiSE), 'ys', label='Secante',alpha=0.5)

plt.title('Raíces de la función')
#plt.xlabel('x')    # No se ven bien
#plt.ylabel('y')

plt.legend()
plt.grid(True)

plt.subplot(212)
# ----
iteBiseccion = [row[0] for row in listaBI]
errRBiseccion = [row[2] for row in listaBI]
# ----
iteRegFal = [row[0] for row in listaRF]
errRRegFal = [row[2] for row in listaRF]
# ----
iteNewRap = [row[0] for row in listaNR]
errRNewRap = [row[2] for row in listaNR]
# ----
iteSecante = [row[0] for row in listaSE]
errRSecante = [row[2] for row in listaSE]

plt.plot(iteBiseccion, errRBiseccion, 'o-', markersize=3,color='r', label='Bisección', lw=0.5)
plt.plot(iteRegFal, errRRegFal, 'd-', markersize=3,color='g', label='Reg Falsa', lw=0.5)
plt.plot(iteNewRap, errRNewRap, 'v-', markersize=3,color='b', label='Newt-Raph', lw=0.5)
plt.plot(iteSecante, errRSecante, 's-', markersize=3, color='y', label='Secante', lw=0.5)

plt.title("Error vs Iteraciones")
plt.xlabel('Iteración')
plt.ylabel('Porcentaje de error')


plt.legend()
plt.grid(True)
plt.show()

