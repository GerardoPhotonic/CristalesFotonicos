"""
Este programa almacena los perfiles de campo de una bicapa de grafeno en 
en configuración ATR inmerso en aire
"""

#%% Establecemos las libreriÂ­as y scripts necesarias
import numpy as np
import FuncionesSigmaWGrafeno as graf
import matplotlib.pyplot as plt

#%% Definimos los parametros computacionales necesarios
#Numero de puntos de simulacion
numDatos=50000

#%% Definimos las constantes fi­sicas necesarias
#Indice de refraccion del prisma
n_pris=4
#Indice de refrccion del espacio de entrada
n_in=1
#Indice de refraccion del espacio de salida
n_out=1
#Frecuencia 
nu=5e12
#Frecuencia angular
w=2*np.pi*nu
#La velocidad de la luz
c=2.99e8 #[m]/[s]
#Angulo incidente en grados
factor=np.pi/180
theta=30.5*factor
#Impedancia del vaci­o
eta0=120 * np.pi #[Ohm]
#Distancia entre el prisma y la capa de grafeno
D=17e-6 #[m]
#Llamamos a la funcion sigma
sigma=graf.SIGMAN(w)

#%% Calculamos las cantidades fisicas que dependen de las constantes
#Calculamos las componentes del vector de onda
Kpris=(w/c)*(((n_pris**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))
Kin=(w/c)*(((n_in**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2)) 
Kout=(w/c)*(((n_out**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))

#Calculamos las impedancias
Yp=(n_pris)**2/(eta0*((n_pris**2 - (n_pris*np.sin(theta))**2)*(1+0j) )**(1/2))
Yin=(n_in)**2 / (eta0*((n_in**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2) )
Yout=(n_out)/ (eta0*((n_out**2 - (n_pris*np.sin(theta))**2)*(1+0j))**(1/2))
Y1=Yp/Yout

#%% Construimos la matriz M0
p11=(Y1/2)

p12=Y1/(2*Yp)
p21=(-Y1/2)
p22=Y1/(2*Yp)
Ap=np.array([[p11,p12],[p21,p22]])

in11=1
in12=-1
in21=Yin
in22=Yin
In=np.array([[in11,in12],[in21,in22]])

M0=Ap@In

#%% Construimos la matriz MD
exp11=np.exp(-1j*Kin*D)
exp12=0
exp21=0
exp22=np.exp(1j*Kin*D)
Exp=np.array([[exp11,exp12],[exp21,exp22]])

o11=(1/2)
o12=1/(2*Yin)
o21=(-1/2)
o22=1/(2*Yin)
Ao=np.array([[o11,o12],[o21,o22]])

sig11=1
sig12=-1
sig21=Yin + sigma
sig22=Yin - sigma
Sig=np.array([[sig11,sig12],[sig21,sig22]])

salida11=np.exp(1j*Kout*D)
salida12=0
salida21=0
salida22=np.exp(-1j*Kout*D)
Salida=np.array([[salida11,salida12],[salida21,salida22]])

MD=Exp@Ao@Sig@Salida

#%% Calculamos la matriz MF
M=M0@MD

#%% Determinamos r y t
t=1 / (M[0,0])
r= M[1,0]/M[0,0]

#%% Determinamos a2, b2
a1=MD[0,0]*t
b1=MD[1,0]*t

#%% Modelamos los campos
z1=np.linspace(-5e-6, 0, numDatos)
z2=np.linspace(0, D, numDatos)
z3=np.linspace(D, D+5e-6, numDatos)

perfil1=np.exp(1j*Kpris*z1) + r*np.exp(-1j*Kpris*z1)
perfil2=a1*np.exp(1j*Kin*z2) + b1*np.exp(-1j*Kin*z2)
perfil3=t*np.exp(1j*Kout*z3)

#%% Mostramos resultados
plt.figure(figsize=(10, 7))
#plt.plot(z1*1e6,(perfil1).real,'b',linewidth=2.5, label='Prisma')
plt.plot(z2*1e6,(perfil2).real,'r',linewidth=2.5, label='Aire')
plt.plot(z3*1e6, (perfil3).real,'g',linewidth=2.5, label='Aire')
#plt.axvline(x=0, color='k', linestyle='--')
plt.axvline(x=17, color='k', linestyle='-')
#plt.ylim([-0.3,0.3])
#plt.xlim([0.8e-5,1.2e-5])
plt.grid()
plt.xlabel('$z [\mu m]$')
plt.ylabel('Amplitudes del campo $H_{y}$')
plt.show()