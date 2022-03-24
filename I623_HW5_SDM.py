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
'''
# the function
def func(Q, b, x):
    #where Q is a matrix, b is a vector, and x is a vector
    return .5*(x.T @ Q @ x) + b.T @ x

# gradient of the function
def grad_func(Q, x):
    #where Q is a matrix and x is a vector
    return np.array(Q @ x).ravel()'''

# fix list error ( should be an array)

def SDM(Q, b, x, N, E):
    '''where Q is a matrix, b is a vector, x is a vector, and 
    E is the stopping condition that tells us when to stop searching. 
    Prints a table of results. '''

    print(' iteration \t x_k \t f(x_k) - f(x_(k+1))')
    print('-------------------------------------------------')
    xk = x
    for i in range(N):
        pk = -1*(Q @ x)                               # step direction at k
        #print(-1*(Q @ x) )
        ak = -1 * ((Q @ x).T @ pk) / (pk.T @ Q @ pk)    # step size at k
        xk1 = xk + (ak * pk)                          # update x_(k+1)

        fxk = (.5*(x.T @ Q @ x) + b.T @ x)    # original function at x_k
        fxk1 = (.5*(x.T @ Q @ x) + b.T @ x)   # original function at x_(k+1)

        fmin = fxk - fxk1                     # towards stopping condition

        print('  ', i, '\t', '[', round(xk[0,0],4), round(xk[1,0],4), round(xk[2,0],4), ']', '\t', round(fmin[0,0],4))  # print current results at k

        # check for stopping condition
        if (fmin < E):
            i = N
        else: 
            xk = xk1
'''
def SDM(Q, b, xk, N, E):
    #where Q is a matrix, b is a vector, x is a vector, and 
    #E is the stopping condition that tells us when to stop searching. 
    #Prints a table of results. 

    print(' iteration \t x_k \t f(x_k) - f(x_(k+1))')
    print('-------------------------------------------------')

    for i in range(N):
        pk = -1*grad_func(Q, xk)                               # step direction at k
        print(type(pk))
        print(type(grad_func(Q, xk)))
        print(grad_func(Q, xk).T @ pk)
        print((pk.T @ Q @ pk))
        ak = -1 * ((grad_func(Q, xk)).T @ pk) / (pk.T @ Q @ pk) # step size at k

        fxk = func(Q, b, xk)   # original function at x_k

        xk1 = xk + (ak * pk)   # update x_(k+1)
        fxk1 = func(Q, b, xk1) # original function at x_(k+1)

        fmin = fxk - fxk1      # towards stopping condition

        print('  ', i, '\t', '[', round(xk[0,0],4), round(xk[1,0],4), round(xk[2,0],4), ']', '\t', round(fmin[0,0],4))  # print current results at k

        # check for stopping condition
        if (fmin < E):
            i = N
        else: 
            xk = xk1
'''
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


# Test cases ============================
'''# test 1
print('\nTest 1: \n')
test1 = SDM(Q1, b1, x01, N1, error_tol)
#print('\nTest 1: \n', test1)
'''
# test 2
print('\nTest 2: \n')
test2 = SDM(Q2, b2, x02, N2, error_tol)

# test 3
print('\nTest 3: \n')
#test3 = SDM(Q3, b3, x03, N3, error_tol)
#print('\nTest 3: \n', test3)

# end