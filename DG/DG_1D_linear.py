import numpy as np
import matplotlib.pyplot as pl
from scipy.integrate import quad
from scipy.misc import derivative

# 1D Scalar Conservation Equation solved using Discontinuous Galerkin method
# Only linear basis used 

c = 2 # flow speed
x_start = 0.0
x_end   = 1.0 
t_start = 0.0
t_end = 5.0
t_elements = 5000 # no of time steps
elements = 32 # no of space elements

x_array = np.linspace(x_start, x_end, elements)
t_array = np.linspace(t_start, t_end, t_elements)
delta_t = (t_end - t_start)/t_elements

def initial_value_function(x):
    return np.sin(2*np.pi*x)

def J (c, phi):
    return c*phi

def basis_0 (x, xl, xr):
    return (xr - x)/(xr - xl)

def basis_1 (x, xl, xr):
    return (x - xl)/(xr - xl)

def generate_mass_matrix (xl, xr):
    mass_matrix = [[0.0, 0.0], [0.0,0.0]]
    mass_matrix[0][0], err = quad(lambda x,a,b: basis_0(x, a, b)*basis_0(x, a, b), xl, xr, args=(xl,xr,))
    mass_matrix[0][1], err = quad(lambda x,a,b: basis_0(x, a, b)*basis_1(x, a, b), xl, xr, args=(xl,xr,))
    mass_matrix[1][0], err = quad(lambda x,a,b: basis_1(x, a, b)*basis_0(x, a, b), xl, xr, args=(xl,xr,))
    mass_matrix[1][1], err = quad(lambda x,a,b: basis_1(x, a, b)*basis_1(x, a, b), xl, xr, args=(xl,xr,))
    return np.array(mass_matrix)

def generate_stiffness_matrix (xl, xr):
    stiffness_matrix = [[0.0,0.0], [0.0,0.0]]
    stiffness_matrix[0][0], err = quad(lambda x,a,b: derivative(basis_0, x, args=(a,b,))*basis_0(x, a, b), xl, xr, args=(xl,xr))
    stiffness_matrix[0][1], err = quad(lambda x,a,b: derivative(basis_0, x, args=(a,b,))*basis_1(x, a, b), xl, xr, args=(xl,xr))
    stiffness_matrix[1][0], err = quad(lambda x,a,b: derivative(basis_1, x, args=(a,b,))*basis_0(x, a, b), xl, xr, args=(xl,xr))
    stiffness_matrix[1][1], err = quad(lambda x,a,b: derivative(basis_1, x, args=(a,b,))*basis_1(x, a, b), xl, xr, args=(xl,xr))
    return np.array(stiffness_matrix)

def plot_elements (x_elements, y_elements):
    for k in range(len(x_elements)):
        pl.plot(x_elements[k], y_elements[k])

initial_array = initial_value_function(x_array)
basis_weights = []
x_elements    = []
complete_basis_weights = []

for i in range(len(initial_array) - 1):
    basis_weights.append([initial_array[i], initial_array[i+1]])
    x_elements.append([x_array[i], x_array[i+1]])

complete_basis_weights.append(basis_weights)
plot_elements(x_elements, basis_weights)
pl.savefig('0.png')
pl.clf()

for i in range(1, len(t_array)):
    for j in range(len(x_array) - 1):

        mass_matrix      = generate_mass_matrix(x_array[j], x_array[j+1])
        stiffness_matrix = generate_stiffness_matrix(x_array[j], x_array[j+1])
        temp_matrix_1    = [0.0, 0.0]
        temp_matrix_1[0] = c*stiffness_matrix[0][0]*basis_weights[j][0] + c*stiffness_matrix[0][1]*basis_weights[j][1]
        temp_matrix_1[1] = c*stiffness_matrix[1][0]*basis_weights[j][0] + c*stiffness_matrix[1][1]*basis_weights[j][1]

        if (j == 0):
            temp_matrix_1[0] = temp_matrix_1[0] + J(c, basis_weights[-1][1])
        else:
            temp_matrix_1[0] = temp_matrix_1[0] + J(c, basis_weights[j - 1][1])

        temp_matrix_1[1] = temp_matrix_1[1] - J(c, basis_weights[j][1])
        basis_weights[j] = basis_weights[j] + delta_t*np.matmul(np.linalg.inv(mass_matrix), np.array(temp_matrix_1))

    if (i >= 100) and (i % 100 == 0):
        plot_elements(x_elements, basis_weights)
        pl.savefig(str(i) + '.png')
        print "Progress = " + str(((i * 1.0)/t_elements)*100.0) + "%"
        pl.clf()

    complete_basis_weights.append(basis_weights)