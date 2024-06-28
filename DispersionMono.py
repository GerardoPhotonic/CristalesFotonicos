"""Este programa grafica la relación de dispersión de una monocapa de grafeno 
para distintos valores de potencial químico"""

#%% Definimos los parámetros físicos EDITABLES#################################
# Índice de refracción del medio 1
n1 = 1
# Índice de refracción del medio 2
n2 = 1
# Potencial químico
mu = 0.3
mu2=0.6
mu3=0.9
# Scattering
tau = 1e-12
###############################################################################

#%% Establecemos las librerías y scripts necesarias
import numpy as np
import matplotlib.pyplot as plt
import FuncionesSigmaWGrafeno as graf

#%% Definimos los parámetros computacionales necesarios
# Creamos un vector para la frecuencia y frecuencia angular
numDatos = 500
frecuencias = 15 * 1e12
nu = np.linspace(0, frecuencias, numDatos)
nu1 = nu * 1e-12  # Frecuencia normalizada a TeraHertz
w = 2 * np.pi * nu

#%% Definimos los parámetros físicos CONSTANTES
# La velocidad de la luz
c = 3e8  # [m]/[s]
# Impedancia del vacío
eta0 = 120 * np.pi  # [Ohm]
# Determinamos sigma
sigma1 = graf.SIGMATOTAL(w, mu, tau)
sigma2 = graf.SIGMATOTAL(w, mu2, tau)
sigma3 = graf.SIGMATOTAL(w, mu3, tau)
#Creamos q_x
qx=(w/c)*np.sqrt(1- ((2/(eta0*sigma1))**2) ) 
qx2=(w/c)*np.sqrt(1- ((2/(eta0*sigma2))**2) )
qx3=(w/c)*np.sqrt(1- ((2/(eta0*sigma3))**2) )

#%% Mostramos los resultados
plt.figure(figsize=(10, 7))
plt.plot( qx.real*1e-6 ,nu1 , 'r', linewidth=2.5, label='$\mu=0.3$')
plt.plot( qx2.real*1e-6 ,nu1, 'b', linewidth=2.5, label='$\mu=0.6$')
plt.plot( qx3.real*1e-6 ,nu1, 'g', linewidth=2.5, label='$\mu=0.9$')
plt.plot((w/c)*1e-6 ,nu1, 'k--', linewidth=1)
plt.plot((w/c)*1e-6*4 ,nu1, 'k--', linewidth=1)
plt.text(0.32, 14, 'aire',)
plt.text(0.64, 7, 'vidrio',)
plt.xlabel('$q_{x}$ real')
plt.ylabel('$f [THz]$')
plt.xlim(0,0.7)
plt.ylim(0,15)
plt.legend()
plt.grid()
plt.show()