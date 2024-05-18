#%%
"""
Este programa almacena la matriz de transferencia para dos capa de grafeno con 
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
n2=3.193  #Pues es el medio entre capas y n=sqrt(eps_r)
#Indice de refrcción del medio 3
n3=1
#La velocidad de la luz
c=3e8 #[m]/[s]
#Ángulo incidente
theta=0 #[rad] Incidencia normal
#Impedancia del vacío
eta0=120*np.pi #[Ohm]
#Distancia entre capas
d=10e-6 #[m]
#Definimos la componente Kz del medio 2
K2z=(w/c)*(n2**2 - (n1*np.sin(theta))**2)**(1/2)
#Aceptancia del medio como función del ángulo de incidencia y del índice de refracción
Z1=(n1**2)/(eta0*(n1**2-(n1*np.sin(theta))**2)**(1/2))
Z2=(n2**2)/(eta0*(n2**2-(n1*np.sin(theta))**2)**(1/2))
Z3=(n3**2)/(eta0*(n3**2 - (n2*np.sin(theta)))**(1/2))
#Determinamos sigma
sigma=graf.SIGMAP(w)

#%% Determinamos la matriz M1 en función de Z1 Y Z2
m11=(1/2)*( 1+(sigma/Z1)+(Z2/Z1) )
m12=(1/2)*( 1+(sigma/Z1)-(Z2/Z1) )
m21=(1/2)*( 1-(sigma/Z1)-(Z2/Z1) )
m22=(1/2)*( 1-(sigma/Z1)+(Z2/Z1) )

#%% Determinamos la matriz M2 en función de Z2 Y Z3
m2_11=(1/2)*( 1+(sigma/Z2)+(Z3/Z2) )
m2_12=(1/2)*( 1+(sigma/Z2)-(Z3/Z2) )
m2_21=(1/2)*( 1-(sigma/Z2)-(Z3/Z2) )
m2_22=(1/2)*( 1-(sigma/Z2)+(Z3/Z2) )

#%% Determinamos la matríz exponencial
e11=np.exp(-1j*K2z*d)
e22=np.exp(1j*K2z*d)

#%% Determinamos la matriz total de transferencia
MT11=m11*e11*m2_11 + m12*e22*m2_21
MT12=m11*e11*m2_12 + m12*e22*m2_22
MT21=m21*e11*m2_11 + m22*e22*m2_21
MT22=m21*e11*m2_12 + m22*e22*m2_22

MTOTAL=np.array([[MT11, MT12], [MT21, MT22]])
#%% Determinamos los coeficientes de transmisión.
r= (MT21) / (MT11)
t= ((1) / (MT11))
T=abs(t)**2
R=abs(r)**2
#%% Mostramos los resultados
plt.figure(figsize=(10, 7))
plt.plot(nu1,T,'r',label='$|T|^{2}$') #Transmitancia
#plt.plot(nu1,R,'b',label='$|R|^{2}$') #Reflectancia
plt.xlabel('Frecuencia (THz)')
plt.ylabel('$|T|^{2}$')
plt.xlim(0,8)
plt.ylim(0,1)
plt.xticks(range(0, 9, 1))
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
#plt.text(4, 0.4, '$\mu_{c}=0.5 eV$',)
plt.legend()
plt.grid()
plt.show()