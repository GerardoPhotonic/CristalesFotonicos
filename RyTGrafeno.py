#%%
"""
Este programa almacena la matriz de transferencia para una capa de grafeno con 
diferentes valores de los indices de refracción; una polarización tipo TM y un 
potenciál químico de 1ev a incidencia normal
"""

#%% Establecemos las librerías y scripts necesarias
import numpy as np
import matplotlib.pyplot as plt
import FuncionesSigmaWGrafeno as graf

#%% Definimos los parámetros computacionales necesarios
#Creamos un vector para la frecuencia y frecuencia angular
numDatos=2000
frecuencias=15*1e12
nu=np.linspace(0,frecuencias,numDatos)
nu1=nu*1e-12 #Frecuencia normalizada a TeraHertz
w=2*np.pi*nu

#%% Definimos las constantes físicas necesarias
#Indice de refracción del medio 1
n1=1
#Indice de refrcción del medio 2
n2=1
#La velocidad de la luz
c=3e8 #[m]/[s]
#Ángulo incidente
theta=0 #[rad] Incidencia normal
#Impedancia del vacío
eta0=120*np.pi #[Ohm]
#Aceptancia del medio como función del ángulo de incidencia y del índice de refracción
Z1=(n1**2)/(eta0*(n1**2-(n1*np.sin(theta))**2)**(1/2))
Z2=(n2**2)/(eta0*(n2**2-(n1*np.sin(theta))**2)**(1/2))
#Determinamos sigma
sigma=graf.SIGMAM(w)

#%% Determinamos la matriz M en función de Z1 Y Z2
m11=(1/2)*( 1+(sigma/Z1)+(Z2/Z1) )
m12=(1/2)*( 1+(sigma/Z1)-(Z2/Z1) )
m21=(1/2)*( 1-(sigma/Z1)-(Z2/Z1) )
m22=(1/2)*( 1-(sigma/Z1)+(Z2/Z1) )

#%% Determinamos los coeficientes de reflexión y transmisión.
r= (m21) / (m11)
t= ((1) / (m11))
R=abs(r)**2
T=abs(t)**2

#%% Mostramos los resultados
plt.figure(figsize=(10, 7))
plt.plot(nu1,R,'r',label='$|R|^{2}$') #Reflectancia
plt.plot(nu1,T,'g',label='$|T|^{2}$') #Transmitancia
plt.xlabel('Frecuencia (THz)')
plt.ylabel('$|R|^{2}, |T|^{2}$')
plt.xlim(0,15)
plt.ylim(0,1)
plt.xticks(range(0, 16, 3))
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
plt.text(3, 0.9, '$\mu_{c}=1 eV$',)
plt.legend()
plt.grid()
plt.show()