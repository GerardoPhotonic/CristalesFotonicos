import numpy as np
import matplotlib.pyplot as plt

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

# Graficar la función y el proceso de descenso de gradiente
x_values = np.linspace(-10, 7, 100)
y_values = np.linspace(-10, 5, 100)
X, Y = np.meshgrid(x_values, y_values)
Z = f(X, Y)

plt.contourf(X, Y, Z, levels=50)
history = np.array(history)
plt.scatter(history[:, 0], history[:, 1], color='red', zorder=5)
plt.plot(history[:, 0], history[:, 1], color='red', label='Descenso de gradiente')
plt.xlabel('x')
plt.ylabel('y')
#plt.title('Descenso de Gradiente en f(x, y) = (x - 3)^2 + (y - 2)^2')
#plt.legend()
plt.show()
