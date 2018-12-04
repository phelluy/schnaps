import numpy as np
import scipy.sparse.linalg as spla
from numpy import linalg as LA

def M1toMin(x,R1):
	return 1/np.tanh(x) -1/x -R1

def M1toMinD(x):
	return -1/(np.sinh(x)*np.sinh(x)) + 1/(x*x)
		
def M1GetBa(R1):
	xk  = _NewtonInitGuess
	absR1 = np.absolute(R1)
	i=1
	error = 0
	if (absR1 > _NewtonInitGuess):
		xk = 1/(1-absR1)
		
	xkp1=xk
	while (i < _NewtonIter):
		xkp1 = xk - M1toMin(xk,absR1)/M1toMinD(xk)
		error = np.absolute(xkp1-xk)
		if error < _NewtonTol :
			if R1 > 0 :
				return xkp1
			else :
				return  -xkp1
		xk = xkp1
		i=i+1
	return 0

def DegToXYZ(nodePolar):
	newArray = np.zeros((_Nd,3))
	for i in range(0, _Nd):
		newArray[i][0] =  np.cos(np.pi*nodePolar[2*i]/180)*np.sin(np.pi*nodePolar[2*i+1]/180)
		newArray[i][1] =  np.sin(np.pi*nodePolar[2*i]/180)*np.sin(np.pi*nodePolar[2*i+1]/180)
		newArray[i][2] =  np.cos(np.pi*nodePolar[2*i+1]/180)
	return newArray

def Func(xk,w,I):
	jac = np.zeros((4,4));
	gk = np.zeros(4)
	F1 = 0; F2x = 0; F2y = 0; F2z = 0;
	DF2xbx=0; DF2xby = 0; DF2xbz = 0; DF2xa = 0;
	DF2ybx=0; DF2yby = 0; DF2ybz = 0; DF2ya = 0;
	DF2zbx=0; DF2zby = 0; DF2zbz = 0; DF2za = 0;
	
	for i in range(0,_Nd):
		fiwi =  _Weight*np.exp(xk[0]*_Node[i][0] + xk[1]*_Node[i][1] +xk[2]*_Node[i][2])
		fiwivx =  fiwi*_Node[i][0]
		fiwivy =  fiwi*_Node[i][1]
		fiwivz =  fiwi*_Node[i][2]
		#Gk elements
		F1 += fiwi
		F2x += fiwivx
		F2y += fiwivy
		F2z += fiwivz
	
		#Jacobian Elements
		DF2xbx += fiwivx*_Node[i][0]
		DF2ybx += fiwivy*_Node[i][0]
		DF2zbx += fiwivz*_Node[i][0]
		
		
		DF2xby += fiwivx*_Node[i][1]
		DF2yby += fiwivy*_Node[i][1]
		DF2zby += fiwivz*_Node[i][1]
		
		DF2xbz += fiwivx*_Node[i][2]
		DF2ybz += fiwivy*_Node[i][2]
		DF2zbz += fiwivz*_Node[i][2]
		
		DF2xa  +=  fiwivx
		DF2ya +=  fiwivy
		DF2za += fiwivz
		
		
	#Multiplication by a 	
	DF2xbx *= xk[3]
	DF2xby *= xk[3]
	DF2xbz *= xk[3]
	DF2ybx *= xk[3]
	DF2yby *= xk[3]
	DF2ybz *= xk[3]
	DF2zbx *= xk[3]
	DF2zby *= xk[3]
	DF2zbz *= xk[3]
	
	#Jacobian Assembly
	jac[0][0]=xk[3]*F2x ; jac[0][1]=xk[3]*F2y ; jac[0][2]=xk[3]*F2z ; jac[0][3]=F1;
	jac[1][0]=DF2xbx ; jac[1][1]=DF2xby ; jac[1][2]=DF2xbz ; jac[1][3]=DF2xa;
	jac[2][0]=DF2ybx ; jac[2][1]=DF2yby ; jac[2][2]=DF2ybz ; jac[2][3]=DF2ya;
	jac[3][0]=DF2zbx ; jac[3][1]=DF2zby ; jac[3][2]=DF2zbz ; jac[3][3]=DF2za;
	#Gk Assembly
	gk[0]=xk[3]*F1-w ; gk[1]=xk[3]*F2x-I[0]; gk[2]=xk[3]*F2y-I[1] ; gk[3]=xk[3]*F2z-I[2];
	# print "Cond(Jac) : ", LA.cond(jac, 2)
	# print jac
	return gk, jac

def M1GetB(x0,w,I):
  i=0
  xk=x0
  while i < _NewtonIter :
    gk,jac = Func(xk,w,I)
    # print "bk : ", xk
    sk = np.linalg.solve(jac, -gk)
    # sk = spla.gmres(jac, -gk, x0=None, tol=1e-05)[0]
    # if i > 100:
    xk = sk + xk
    # else:
    # xk = xk + 0.001*gk
    i=i+1
  return xk

		
_NewtonIter = 10
_NewtonTol = 1e-8
_NewtonInitGuess = 0.4

#Mesh loading 

# _NodePolar = np.array([ 0.000000000000000, 90.000000000000000, 180.000000000000000, 90.000000000000000, 90.000000000000000, 90.000000000000000, -90.000000000000000, 90.000000000000000, 90.000000000000000, 0.000000000000000, 90.000000000000000, 180.000000000000000])
# _Nd = len(_NodePolar)/2
# _Weight = 0.166666666666667
# _Node = DegToXYZ(_NodePolar)
_Node= np.loadtxt('./lebedev_005.txt')
_Nd = len(_Node)
_Weight =  4*np.pi/_Nd


if __name__ == '__main__':

  xk = np.zeros(4)
  I = np.zeros(3)

  #Input parameter
  w=1
  I[0] =0.10*w
  I[1] =0
  I[2] =0

  Inorm = np.sqrt(I[0]*I[0] + I[1]*I[1] + I[2]*I[2])
  print "Inorm :", Inorm

  print "R1 :", I[0]/w, I[1]/w, I[2]/w
  #Method coming from the continuous optimisation problem
  bx = M1GetBa(I[0]/w)
  by = M1GetBa(I[1]/w)
  bz = M1GetBa(I[2]/w)
  bnorm = np.sqrt(bx*bx+by*by+bz*bz)
  a=w*bnorm/(4*np.pi*np.sinh(bnorm))
  # print "b : ",bx,by,bz,a


  '''
  Solving the discrete problem 
  '''
  #Init Guess for MultiD NEwton Method
  xk[0] =  1/(1-np.absolute(I[0]/w)) #bx
  xk[1] = 1/(1-np.absolute(I[1]/w)) #by
  xk[2] = 1/(1-np.absolute(I[2]/w)) #bz
  bnorm = np.sqrt(xk[0]*xk[0] + xk[1]*xk[1] + xk[2]*xk[2])
  xk[3] = w*bnorm/(4*np.pi*np.sinh(bnorm)) #a

  csumm=0
  csumm1=0
  csumm2=0
  csumm3=0 
  dsumm=0
  dsumm1=0
  dsumm2=0
  dsumm3=0

  xkp1 = M1GetB(xk,w,I)
  print "bx by bz a continuous : ", bx,by,bz,a
  print "bx by bz a discrete : ", xkp1[0], xkp1[1], xkp1[2], xkp1[3]
  for i in range(0,_Nd):
    temp =  xkp1[3]*4*np.pi*_Node[i][3]*np.exp(xkp1[0]*_Node[i][0] + xkp1[1]*_Node[i][1] +xkp1[2]*_Node[i][2])
   
    # temp = 4*np.pi*_Node[i][3]
    # temp2 = 4*np.pi*_Node[i][3]
    temp2 =  a*4*np.pi*_Node[i][3]*np.exp(bx*_Node[i][0] + by*_Node[i][1] +bz*_Node[i][2])
    print "f pour i=",i," valeur : ",a*np.exp(bx*_Node[i][0] + by*_Node[i][1] +bz*_Node[i][2])
    csumm = csumm + temp2
    csumm1 = csumm1 + temp2*_Node[i][0]
    csumm2 = csumm2 + temp2*_Node[i][1]
    csumm3 = csumm3 + temp2*_Node[i][2]
    dsumm = dsumm + temp
    dsumm1 = dsumm1 + temp*_Node[i][0]
    dsumm2 = dsumm2 + temp*_Node[i][1]
    dsumm3 = dsumm3 + temp*_Node[i][2]
  print "Expected values", w,I[0],I[1],I[2]
  print "Continuous values", csumm, csumm1, csumm2,csumm3
  print "Continuous errors", abs(csumm-w), abs(csumm1-I[0]), abs(csumm2-I[1]),abs(csumm3-I[2])
  print "Discrete values", dsumm, dsumm1, dsumm2,dsumm3
  print "Discrete errors", abs(dsumm-w), abs(dsumm1-I[0]), abs(dsumm2-I[1]),abs(dsumm3-I[2])
  
	
