#%%
"""
Este programa almacena la función del grafeno. La suma total de sigma
sigma total=sigma intra + sigma inter y todas las sigmas necesarias en los
proyectos futuros
"""
#%% Solicitamos las librerias necesarias
import numpy as np
import math as mt

#%% Definimos los parámetros computacionales necesarios
#Creamos un vector para la frecuencia y la frecuencia angular
numDatos=8000
frecuencias=15*1e12
f=np.linspace(0,frecuencias,numDatos)
f1=f*1e-12
w=2*np.pi*f

#%% Creamos la función que calcula sigma de w y mu=0.1 (GRAFICAS)
def SIGMAW1(w):
    mu1=0.1 #[eV]
    hbarra=6.582*1e-16
    factor= (1.0j*1.602*1e-19)/(np.pi*hbarra)
    intra=(mu1)/(hbarra*(w+1.0j*1e12))
    inter=(1/4)*np.log((2*mu1+hbarra*(w+1.0j*1e12))/(2*mu1-hbarra*(w+1.0j*1e12)))
    sigmaW1=(factor*intra) - (factor*inter)
    return sigmaW1

#%% Creamos la función que calcula sigma de w y mu=0.2 (GRAFICAS)
def SIGMAW2(w):
    mu2=0.2 #[eV]
    hbarra=6.582*1e-16
    factor= (1.0j*1.602*1e-19)/(np.pi*hbarra)
    intra=(mu2)/(hbarra*(w+1.0j*1e12))
    inter=(1/4)*np.log((2*mu2+hbarra*(w+1.0j*1e12))/(2*mu2-hbarra*(w+1.0j*1e12)))
    sigmaW2=(factor*intra) - (factor*inter)
    return sigmaW2

#%% Creamos la función que calcula sigma de w y mu=0.3 (GRAFICAS)
def SIGMAW3(w):
    mu3=0.3 #[eV]
    hbarra=6.582*1e-16
    factor= (1.0j*1.602*1e-19)/(np.pi*hbarra)
    intra=(mu3)/(hbarra*(w+1.0j*1e12))
    inter=(1/4)*np.log((2*mu3+hbarra*(w+1.0j*1e12))/(2*mu3-hbarra*(w+1.0j*1e12)))
    sigmaW3=(factor*intra) - (factor*inter)
    return sigmaW3

#%% Creamos la función grafeno (LA TRANSMISIÓN EN UNA BICAPA)
def SIGMAP(w):
    mu=1 #[eV]
    hbarra=6.582*1e-16
    factor= (1.0j*1.602*1e-19)/(np.pi*hbarra)
    intra=(mu)/(hbarra*(w+1.0j*2e12))
    inter=(1/4)*np.log((2*mu+hbarra*(w+1.0j*2e12))/(2*mu-hbarra*(w+1.0j*2e12)))
    sigmaP=(factor*intra) - (factor*inter)
    return sigmaP

#%% Creamos la función grafeno (RyT para una monocapa)
def SIGMAM(w):
    mu=1 #[eV]
    hbarra=6.582*1e-16
    factor= (1.0j*1.602*1e-19)/(np.pi*hbarra)
    intra=(mu)/(hbarra*(w+1.0j*2e12))
    inter=(1/4)*np.log((2*mu+hbarra*(w+1.0j*2e12))/(2*mu-hbarra*(w+1.0j*2e12)))
    sigmaW=(factor*intra) - (factor*inter)
    return sigmaW

#%% Creamos la función grafeno (T para el caso n-ésimo)
def SIGMAN(w):
    mu=0.8 #[eV]
    hbarra=6.582*1e-16
    factor= (1.0j*1.602*1e-19)/(np.pi*hbarra)
    intra=(mu)/(hbarra*(w+1.0j*1e12))
    inter=(1/4)*np.log((2*mu+hbarra*(w+1.0j*1e12))/(2*mu-hbarra*(w+1.0j*1e12)))
    sigmaW=(factor*intra) - (factor*inter)
    return sigmaW

#%% Creamos la función que calcula sigma de w con valores en T
def SIGMAT(w):
    T=300 #[Kelvin]
    muT=0.5 #[eV]
    hbarra=6.582e-16 #[eV*S]
    gamma=0.5e12 #[1/S] (1/1ps)
    ecuadrada= 2.56e-38 #[C]^2
    kb=8.617e-5 #[eV]/[K]
    
    factor1=(-1j*ecuadrada*kb*T)/(np.pi*hbarra**2*(w-1j*gamma))
    factor2=(muT)/(kb*T) + 2*np.log(mt.exp(-muT/(kb*T)) + 1)
    intra=factor1*factor2
    
    factor3=(-1j*ecuadrada)/(4*np.pi*hbarra)
    factor4=(2*muT-(w-1j*gamma)*hbarra)/(2*muT+(w-1j*gamma)*hbarra)
    inter=factor3*np.log(factor4)
    
    sigmaT=inter+intra
    return sigmaT
