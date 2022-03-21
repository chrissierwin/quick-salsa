'''Implement the steepest descent method for minimizing a convex 
quadratic function f(x) = .5*x^T Qx + b^T x in Python. Stop when 
f(x(k)) - f(x(k+1)) < error or the number of steps exceeds N. Use
Q, b, error, N and the initial guess x(0) as your inputs. Your 
output should be a table with the columns in the k-th row showing 
(i) the step number k; (ii) x(k); (iii) f(x(k)) - f(x(k-1)). Test 
your program with the following data:'''

'''Note: The method is called the method of steepest descent because 
for analytic g(z), constant phase contours are equivalent to steepest 
descent contours.'''


import numpy as np









# test program with the following data =============
error_tol = 1e-6

# test 1
Q1 = np.array([[2, -1,  0],
               [-1, 2, -1],
               [0, -1,  2]])
b1 = np.array([[1],
               [0],
               [1]])
x01 = np.array([[3],
                [5],
                [7]])
N1 = 30

# test 2
Q2 = np.array([[2, -1,  0],
               [-1, 2, -1],
               [0, -1,  2]])
b2 = np.array([[1],
               [0],
               [1]])
x02 = np.array([[-1],
                [ 2],
                [-3]])
N2 = 30

# test 3
Q3 = np.array([[3  -1, 0, 1],
               [-1, 4, 0, 2],
               [0,  0, 2, 1],
               [1,  2, 1, 5]])
b3 = np.array([[-1],
               [0],
               [1],
               [2]])
x03 = np.array([[1],
                [2],
                [3],
                [4]])
N3 = 40