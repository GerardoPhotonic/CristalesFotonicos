import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

# Definir la función f(x, y) = (x - 3)^2 + (y - 2)^2
def f(x, y):
    return (x - 3) ** 2 + (y + 6) ** 2

# Definir las derivadas parciales usando diferencia central
def df_dx(x, y, h=1e-5):
    return (f(x + h, y) - f(x - h, y)) / (2 * h)

def df_dy(x, y, h=1e-5):
    return (f(x, y + h) - f(x, y - h)) / (2 * h)

# Implementar el algoritmo de descenso de gradiente para dos variables
def gradient_descent(starting_point, learning_rate, num_iterations):
    x, y = starting_point
    history = [(x, y)]
    
    for _ in range(num_iterations):
        grad_x = df_dx(x, y)
        grad_y = df_dy(x, y)
        x = x - learning_rate * grad_x
        y = y - learning_rate * grad_y
        history.append((x, y))
    
    return (x, y), history

# Parámetros del algoritmo
starting_point = (0, 0)  # Punto de inicio
learning_rate = 0.1  # Tasa de aprendizaje
num_iterations = 50  # Número de iteraciones

# Ejecutar el descenso de gradiente
min_point, history = gradient_descent(starting_point, learning_rate, num_iterations)

# Imprimir el resultado
print(f"El valor mínimo de (x, y) es: {min_point}")

# Graficar la función y el proceso de descenso de gradiente en 3D
x_values = np.linspace(-10, 10, 100)
y_values = np.linspace(-10, 10, 100)
X, Y = np.meshgrid(x_values, y_values)
Z = f(X, Y)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, abs(Z), cmap='viridis', alpha=0.8)

history = np.array(history)
Z_history = f(history[:, 0], history[:, 1])
ax.plot(history[:, 0], history[:, 1], Z_history, color='red', marker='o', label='Descenso de gradiente')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('f(x, y)') 
plt.show()
