
def f(x):
    x = complex(x)
    return (((x+5)*x-9)*x-85)*x-136
def RootFinder (f, p0, p1, p2, tol, maxit=100):
    if (tol <= 0) or (maxit<1):
        print('Either tolerance=',tol,'or the max. number of iterations ', maxit, 'is illogical, try again')
        return 
    
    import numpy as np
    import cmath

    p0 = complex(p0)
    p1 = complex(p1)
    p2 = complex(p2)

    n=0 
    while n < maxit:
       f0= f(p0)
       f1= f(p1)
       f2= f(p2)
       h0= p0 - p2
       h1= p1 - p2
       denom = h0 * h1 * (h0 - h1)
       if denom ==0:
         print('Unsuitable starting values, pick distinct p0,p1,p2.')
         return 
       a= (h1 * (f0 - f2) - h0 * (f1  - f2)) / denom
       b= ((h0**2) * (f1-f2) - (h1**2) * (f0-f2)) / denom
       c=f2

       disc= b * b - 4*a*c
       sqrt_disc= cmath.sqrt(disc)

       s=np.sign(b.real)
       if s==0:
           s=1.0 
    
       denom_choice= b + s * sqrt_disc
       if denom_choice == 0: 
          denom_choice = b - s * sqrt_disc
          if denom_choice==0:                               
              print('iteration cannot be done')
              return 
       p3= p2 - 2 * c/denom_choice 
       f3= f(p3)

       n= n+1 
       print(n, 'p_n=', p3,'f(p_n)',f3)

       if abs(f3) < tol:
          print('method converged in', n, 'iterations')
          print('final approximation:',p3)
          print('|f(p_n)|=',abs(f3))
          return p3, n
       p0, p1, p2 = p1, p2, p3
    print('ERROR: method did not converge within the specified max of:',maxit, 'iterations')
    return n 
def fmt_complex_5(z):
    z= complex(z)
    if abs(z.imag) < 1e-14:
        return str(round(z.real, 5))
    
    if z.imag >= 0 :
        sign = '+' 
    else:
       sign= '-'
    real_part= str(round(z.real,5))
    imag_part= str(round(abs(z.imag), 5))
    return real_part+ sign + imag_part + 'j'    
tol = 1e-4

print('case A: p0=0, p1=1, p2=2, tol=1e-4')
rootA,itA= RootFinder(f, 0, 1, 2, tol)

print('final output A:')
if rootA is not None:
    valA = f(rootA)
    print('approx root 5dp:', fmt_complex_5(rootA))
    print('absolute |f(root)|:', f'{abs(valA):.5g}')
    print('iterations taken:', itA)
else:
    print('max iterations exceeded, not converged')
    print('iterations taken:', itA)
    
print('case B: p0=3, p1=4, p2=5, tol=1e-4 ')
rootB,itB = RootFinder(f, 3, 4, 5, tol)

print('final output B:')
if rootB is not None:
    valB = f(rootB)
    print('approx root 5dp:', fmt_complex_5(rootB))
    print('absolute |f(root)|:', f'{abs(valB):.5g}')
    print('iterations taken:', itB)
else:
    print('max iterations exceeded, not converged')
    print('iterations taken:', itB) 