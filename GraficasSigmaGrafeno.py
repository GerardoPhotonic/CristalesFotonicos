#%%
"""
Este programa grafica la función sigma(w) para valores a bajas temperaturas
empleando las funciones hechas en el código FuncionesSigmaWGrafeno.py e importándolas
como módulos. 
"""
#%% Llamamos a los scrips y librerías necesarias
import matplotlib.pyplot as plt
import numpy as np
import FuncionesSigmaWGrafeno as graf

#%% Definimos las constantes físicas necesarias
hbarra=6.582e-16 # [eV*s]
#Definimos el factor de norlmalización
sigma0=(1.602e-19)/(4*hbarra)

#%% Definimos los parametros computacionales necesarios 
#Creamos un vector para la frecuencia y frecuencia angular
numDatos=8000
frecuencias=1*1e15
f=np.linspace(0,frecuencias,numDatos)
f1=f*1e-12 #El mismo vector frecuencia pero normalizado para THz
w=2*np.pi*f

#%% Mostramos los resultados
plt.figure(figsize=(10, 7))

#Generamos el gráfico de la función normalizada con mu = 0.3
plt.plot(f1, graf.SIGMAW3(w).imag/sigma0, 'r--',)
plt.plot(f1, graf.SIGMAW3(w).real/sigma0, 'r')
plt.text(150, -1.5, '$\mu_{c}=0.3 eV$')

#Generamos el gráfico de la función normalizada con mu = 0.2
plt.plot(f1, graf.SIGMAW2(w).real/sigma0, 'b')
plt.plot(f1, graf.SIGMAW2(w).imag/sigma0, 'b--')
plt.text(100, -1.5, '$\mu_{c}=0.2 eV$')

#Generamos el gráfico de la función normalizada con mu = 0.1
plt.plot(f1, graf.SIGMAW1(w).real/sigma0, 'k', label='Parte real')
plt.legend()
plt.plot(f1, graf.SIGMAW1(w).imag/sigma0, 'k--', label='Parte imaginaria')
plt.legend()
plt.text(50, -1.5, '$\mu_{c}=0.1 eV$')

#Damos las especificaciónes del gráfico
plt.xlabel('Frecuencia (THz)')
plt.ylabel('$\sigma/ \sigma_{0}$')
plt.xlim(20,200)
plt.ylim(-2,2)
plt.grid()
plt.show()



