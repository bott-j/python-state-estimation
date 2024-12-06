#!/usr/bin/env python
"""simpleGaussNewton.py: Provides example in using the Gauss-Newton algorithm\
for linear and non-linear systems of equations. """

# Import thrid-party packages
import numpy as np

# Set the maximum number of iterations
MAX_ITER = 200

## Single non-linear equation
#  (one state variable, two measurements) 
x = np.transpose([[-100]])
z = np.transpose([[2]])

# For each iteration
for i in range(0, MAX_ITER):
     J = [[3*x[0][0]**2]]
     h_vk = np.transpose([[x[0][0]**3]])
     x_next = np.array(x) + np.linalg.inv(np.transpose(J)@J)@np.transpose(J)@(np.array(z) - np.array(h_vk)) 
     x = x_next

# Display the result
print("Example 1: x = " + str(x))

## Two independant equations, one non-linear one linear
#  (two state variables, two measurements)
x = np.transpose([[-100,-100]])
z = np.transpose([[1,1]])

# For each iteration
for i in range(0, MAX_ITER):
     J = [[3*x[0][0]**2,0],[0,3*x[1][0]**2]]
     h_vk = np.transpose([[x[0][0]**3,x[1][0]**3]])
     x_next = np.array(x) + np.linalg.inv(np.transpose(J)@J)@np.transpose(J)@(np.array(z) - np.array(h_vk)) 
     x = x_next

# Display the result
print("Example 2: x = " + str(x))

## Two non-linear simultaneuous equations
#  (two state variables, two measurements) 
x = np.transpose([[-100,-100]])
z = np.transpose([[1,2]])

# For each iteration
for i in range(0, MAX_ITER):
     J = [[3*x[0][0]**2,1],[2*x[0][0],1]]
     h_vk = np.transpose([[x[0][0]**3+x[1][0],x[0][0]**2+x[1][0]]])
     x_next = np.array(x) + np.linalg.inv(np.transpose(J)@J)@np.transpose(J)@(np.array(z) - np.array(h_vk)) 
     x = x_next

# Display the result
print("Example 3: x = " + str(x))
