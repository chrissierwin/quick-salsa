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

def SDM(Q, b, x, N, E):
    '''where Q is a matrix, b is a vector, x is a vector, and 
    E is the stopping condition that tells us when to stop searching. 
    Prints a table of results. '''

    print(' iteration \t\t   x_k \t\t     f(x_k) - f(x_(k-1))')
    print('-------------------------------------------------------------------')
    xk = x; i=0        
    fxkm1 = (.5*(xk.T @ Q @ xk) + b.T @ xk)[0,0]                  # original function at x_k
    
    # iterate until stopping condition is met
    while i < (N+1):
        pk = -1*(Q @ xk + b)                                      # step direction at k
        ak = (-1*((Q @ xk + b).T @ pk) / (pk.T @ Q @ pk))[0,0]    # step size at k
        xk1 = xk + (ak * pk)                                      # update x_(k+1)

        fxk = (.5*(xk.T @ Q @ xk) + b.T @ xk)[0,0]                # original function at x_k
        fxkp1 = (.5*(xk1.T @ Q @ xk1) + b.T @ xk1)[0,0]           # original function at x_(k+1)

        fmin = fxk - fxkp1                                        # f(x_k) - f(x_(k+1))
        fminp = np.abs(fxk - fxkm1)                               # f(x_k) - f(x_(k-1))

        if len(b) == 3:
            print('  ', i, '\t', '[', '{:.4f}'.format(xk[0,0]), '{:.4f}'.format(xk[1,0]), '{:.4f}'.format(xk[2,0]), ']', '\t', '{:.4f}'.format(fminp))  # print current results at k
        elif len(b) == 4:
            print('  ', i, '\t', '[', '{:.4f}'.format(xk[0,0]), '{:.4f}'.format(xk[1,0]), '{:.4f}'.format(xk[2,0]), '{:.4f}'.format(xk[3,0]), ']', '\t', '{:.4f}'.format(fminp))  # print current results at k
        else:
            print('cannot print table for the given matrices')
        
        # check for stopping condition
        if (fmin < E):
            i = N+1
        else: 
            xk = xk1
            fxkm1 = fxk
            i += 1


# test program with the following data =============
error_tol = 0.000001 #1e-6

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
Q3 = np.array([3, -1, 0, 1, -1, 4, 0, 2, 0, 0, 2, 1, 1, 2, 1, 5]).reshape(4,4)
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
# test 1
print('\nTest 1: \n')
test1 = SDM(Q1, b1, x01, N1, error_tol)

# test 2
print('\nTest 2: \n')
test2 = SDM(Q2, b2, x02, N2, error_tol)

# test 3
print('\nTest 3: \n')
test3 = SDM(Q3, b3, x03, N3, error_tol)

# end