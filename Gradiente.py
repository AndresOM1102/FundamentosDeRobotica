import sympy as sp
from f import *
from copy import copy


def cinematica_directa(q):
    T01 = dh(1,q[0],1,0)
    T12 = dh(1, q[1], 0, np.pi/2)
    T23 = dh(0, q[2], 1, 0)
    T34 = dh(0,q[3],1,0)
    Tf = T01@T12@T23@T34
    return Tf

def jacobian_position(q, delta=0.0001):
    J = np.zeros((3,len(q)))
    T = cinematica_directa(q)
    for i in range(len(q)):
        # Copiar la configuracion articular inicial
        dq = copy(q)
        dq[i] += delta
        # Transformacion homogenea luego del incremento (q+delta)
        T_inc = cinematica_directa(dq)
        # Aproximacion del Jacobiano de posicion usando diferencias finitas
        J[0:3,i]=(T_inc[0:3,3]-T[0:3,3])/delta
    return J


import numpy as np
cos = np.cos
sin = np.sin

q = np.array([0.2, 0.2,0.2,0.2])  # Valor inicial en el espacio articular
xd = np.array([0.52, 0.5,3.9])  # Valor deseado en el espacio cartesiano

epsilon = 1e-3  # Tolerancia del error
max_iter = 1000  # Máximo numero de iteraciones
alpha = 0.2
# Iteraciones: Método de Gradiente
norma_error = []
iteracion = []


for i in range(max_iter):

    Jacob = jacobian_position(q)
    pos = cinematica_directa(q)[0:3,3]
    e = xd - pos
    q = q + alpha*np.dot(Jacob.T, e)

    norma_error.append(np.linalg.norm(e))
    iteracion.append(i)

    # Condición de término
    if (np.linalg.norm(e) < epsilon):
        print("Convergió en la iteración ",i)
        break


print(q)


# Ploteo de figura :

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

plt.plot(iteracion, norma_error)
plt.title("Convergencia de la cinematica inversa"); plt.xlabel(" iteracion"); plt.ylabel("norma del error")
plt.grid()
plt.show()