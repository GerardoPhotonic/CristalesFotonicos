"""
Este programa almacena las gráficas del espectro de reflexión de 8 capas de 
grafeno con la configuración de Otto. Un prisma, 8 capas de grafeno y un 
dieléctrico igual entre capas de grafeno 
"""

#%% Establecemos las libreria­as y scripts necesarias
import numpy as np
import FuncionesSigmaWGrafeno as graf
import matplotlib.pyplot as plt

#%% Definimos los parametros computacionales necesarios
#Numero de puntos de simulacion
numDatos=20000

#%% Definimos las constantes fi­sicas necesarias
#Indice de refraccion del prisma
n_pris=4
#Indice de refrccion del espacio de entrada
n_in=1
#Indice de refrccion del espacio entre capas de grafeno
n=1
#Indice de refraccion del espacio de salida
n_out=1
#Frecuencia 
nu=5e12
#Frecuencia angular
w=(2*np.pi)*nu
#La velocidad de la luz
c=2.99e8 #[m]/[s]
#Angulo incidente en grados
factor=np.pi/180
theta=np.linspace(0*factor, 89*factor, numDatos) #[rad]
#Impedancia del vaci­o
eta0=120 * np.pi #[Ohm]
#Distancia entre el prisma y la primer capa de grafeno
D=10e-6 #[m]
#Distancia entre capas de grafeno
d=1e-6 #[m]
#Llamamos a la funcion sigma
sigma=graf.SIGMAN(w)
#Numero de capas
capas=10 #Minimo 2
#Definimos la potencia
potencia=capas-1   # potencia = numero de capas - 1

#%% Calculamos las cantidades fisicas que dependen de las constantes
#Calculamos las componentes del vector de onda
Kpris=(w/c)*(((n_pris**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))
Kin=(w/c)*(((n_in**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2)) 
K=(w/c)*(((n**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))
Kout=(w/c)*(((n_out**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))

#Calculamos las impedancias
Yp=(n_pris)**2/(eta0*((n_pris**2 - (n_pris*np.sin(theta))**2)*(1+0j) )**(1/2))
Yin=(n_in)**2 / (eta0*((n_in**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2) )
Y=(n)**2 / (eta0*((n**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))
Yout=(n_out)/ (eta0*((n_out**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))

#%% Construimos las matrices 
#Matriz asociada al prisma
p11=(1/2)*np.ones(numDatos)
p12=1/(2*Yp)
p21=(-1/2)*np.ones(numDatos)
p22=1/(2*Yp) 
Ap=np.array([[p11,p12],[p21,p22]])

#Matriz asociada al medio 1
in11=np.cos(Kin*D)
in12=-(1j/Yin)*np.sin(Kin*D)
in21=-1j*Yin*np.sin(Kin*D)
in22=np.cos(Kin*D)
Ain=np.array([[in11,in12],[in21,in22]])

#Matriz asociada al medio 2
a11=np.cos(K*d)
a12=-(1j/Y)*np.sin(K*d)
a21=-1j*Y*np.sin(K*d) + (sigma*np.cos(K*d))
a22=np.cos(K*d) - (1j*(sigma/Y)*np.sin(K*d))
A=np.array([[a11,a12],[a21,a22]])

#Matriz asociada al medio 2 elevada a una potencia
Apot=np.zeros([2,2,numDatos], dtype=complex)
for i in range(len(theta)):
    Apot[:,:,i]=np.linalg.matrix_power(A[:,:,i], potencia)
    
#Matriz asociada al medio final
af11=np.ones(numDatos)   
af12=(-1)*np.ones(numDatos)
af21=Yout + sigma
af22=Yout - sigma
Af=np.array([[af11,af12],[af21,af22]])


#%% Realizamos el producto de matrices
F1=np.zeros([2,2,numDatos], dtype=complex)
for i in range(len(theta)):
    F1[:,:,i]=Ap[:,:,i]@ Ain[:,:,i]
    
F2=np.zeros([2,2,numDatos], dtype=complex)
for i in range(len(theta)):
    F2[:,:,i]=F1[:,:,i]@ Apot[:,:,i]

Ff=np.zeros([2,2,numDatos], dtype=complex)
for l in range(len(theta)):
    Ff[:,:,l]=F2[:,:,l]@ Af[:,:,l]
    
#%% Calculamos los coeficientes de reflexión
r=(Ff[1,0]) / (Ff[0,0])
R=abs(r)**2

#%% Mostramos resultados
plt.figure(figsize=(10, 7))
plt.plot(theta/factor,R,'b',linewidth=2.5,label='$|R|^{2}$')
plt.axvline(x=22.41, color='r', linestyle='--')
#plt.xlim([22,22.6])
#plt.ylim([0,0.4])
plt.xlabel('Angulo $ \phi $')
plt.ylabel('$|R|^{2}$')
plt.legend()
plt.grid()
plt.show()
