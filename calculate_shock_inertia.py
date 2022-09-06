import numpy as np
def calculate_shock_inertia(d,x,y,z):

  #calculate inertia tensor
  dnorm = np.sum(d)

  flag_x2_tensor = 0

  print "flag_x2_tensor = ",flag_x2_tensor

  if(flag_x2_tensor==0):
    Ixx = np.sum( d*(y*y + z*z) )/dnorm
    Iyy = np.sum( d*(x*x + z*z) )/dnorm
    Izz = np.sum( d*(x*x+ y*y) )/dnorm
    Ixy = np.sum( d*(x*y) )/dnorm
    Ixz = np.sum( d*(x*z) )/dnorm
    Iyz = np.sum( d*(y*z) )/dnorm
  else:

    Ixx = np.sum( d*x**2 )/dnorm
    Iyy = np.sum( d*y**2)/dnorm
    Izz = np.sum( d*z**2 )/dnorm
    Ixy = np.sum( d*x*y )/dnorm
    Ixz = np.sum( d*(x*z) )/dnorm
    Iyz = np.sum( d*(y*z) )/dnorm

  #Ixx = total( (y^2 + z^2) ) 
  #Iyy = total( (x^2 + z^2) )
  #Izz = total( (x^2 + y^2) )
  #Ixy = total( (x*y) )
  #Ixz = total( (x*z) )
  #Iyz = total( (y*z) )


  Iyx = Ixy
  Izx = Ixz
  Izy = Iyz

  #print(Ixx)

  a = np.zeros((3,3), dtype=np.float32)

  
  a[0,0] =     Ixx
  a[0,1] = -1.*Ixy
  a[0,2] = -1.*Ixz

  
  a[1,0] = -1.*Iyx
  a[1,1] =     Iyy
  a[1,2] = -1.*Iyz

  a[2,0] = -1.*Izx
  a[2,1] = -1.*Izy
  a[2,2] =     Izz
  

  #print(A)

  flag_svdc = 1

  if(flag_svdc==0):

    '''
    #SVDC, A, W, U, V
    #print,"Eigendcomposition..."
    #w = la_eigenql(a, eigenvectors = u)
    w = eigenql(a, /ascending, eigenvectors = u)
    v = u
    s = "Eigendcomposition... %e %e " % (A(2,0)*u(0,2)+A(2,1)*u(1,2)+A(2,2)*u(2,2),w(2)*u(2,2),w)
    '''
  else:

    U,s,V = np.linalg.svd(a, full_matrices=True)

    #
    #The SVD is commonly written as a = U S V.H. 
    #The v returned by this function is V.H and u = U.
    #If U is a unitary matrix, it means that it satisfies U.H = inv(U).
    #The rows of v are the eigenvectors of a.H a. 
    #The columns of u are the eigenvectors of a a.H. 
    #For row i in v and column i in u, the corresponding eigenvalue is s[i]**2.



    #print(U.shape)
    #print(s.shape)
    #print(V.shape)


    #see http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html

    print("*****u******")
    print(U)
    print("*****s******")
    print(s)
    print("*****V******")
    print(V)
    print("*******")
    S = np.zeros((3,3), dtype=np.complex)
    S[:3, :3] = np.diag(s)
    print(S)
    #print(np.dot(U,np.dot(S,V)))
    #print(np.allclose(s, np.dot(U,np.dot(S,V))))
    print(np.allclose(a, np.dot(U, np.dot(S, V))))

    '''
    SVDC, A, W, U, V, /column

    #SVDC does not sort the values of w
    wmax = max(w,wi)  ;find the maximum value
    #sort w
    print,w(0),w(1),w(2)
    wi = sort(-1.0 * w)

    wtmp = w
    vtmp = v
    utmp = u
    for i=0,2 do begin
      wtmp(i)   = w(wi(i))
      utmp(*,i) = u(*,wi(i))
      vtmp(*,i) = v(*,wi(i))
    endfor

    w = wtmp
    v = vtmp
    u = utmp
    '''

  w = s
  return U, w, V
