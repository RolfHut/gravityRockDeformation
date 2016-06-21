#python class which computes 4-pi Geodesy style normalized assocaiited Legendre functions and its derivative
#copyright R. Rietbroek
#Date: 20 June 2016
#License: pending opensource  license to be determined, all rights reserved
#todo: Performance (array accessing order) and stability checks
#literature:
#Holmes, Simon A., and Will E. Featherstone. "A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions." Journal of Geodesy 76.5 (2002): 279-299. 
#routines here are thought to be stable up to degree and order 2700 (after that under/overflow may occur)
#The performance of the routine can likely be improved by making better use of numpy array operations
import numpy as np
import math
#import imp
class Pnm:
	#Constructor precomputes some stuff woth square roots  based on the requested maximum degree (
	def  __init__(self,nmax):
		self.nmax=nmax
		self.Wnn=np.zeros(self.nmax+1)
		self.Wnm=np.zeros([self.nmax+1,self.nmax+1])
		self.Wnn[1]=math.sqrt(3)
		for n in range(2,self.nmax+1):
			self.Wnn[n]=math.sqrt((2*n+1)/(2.0*n))

		#also precompute  Wnm array
		for n in range(1,self.nmax+1):
			for m in range(0,n):
				self.Wnm[n][m]=math.sqrt((2*n+1.0)/(n+m)*(2*n-1.0)/(n-m))


    #compute the value of the associated legendre function at x (0..1) usually x is cos(colatitude)
	def at(self,cosTheta):
		assert cosTheta>=-1 and cosTheta<=1,"Error in Pnm.at(x): -1 <= x <= 1 not fulfilled"
		P=np.zeros([self.nmax+1,self.nmax+1])
		sinTheta=math.sqrt(1-cosTheta**2)
		
		#fill up diagonal with 1^-280 for stability concerns. This is the seed of the recursion
		for n in range(0,self.nmax+1):
			P[n][n]=1e-280
			
		#also compute the off-diagonal elements (apply recursion for fixed order, varying degree)
		sect=1e280
		for m in range(0,self.nmax):
			n=m+1
			P[n][m]=self.Wnm[n][m]*cosTheta*P[n-1][m]	
			##proceed recursion with 2 terms
			for n in range(m+2,self.nmax+1):
				P[n][m]=self.Wnm[n][m]*(cosTheta*P[n-1][m]-P[n-2][m]/self.Wnm[n-1][m])

			#now rescale the values back with 1e280*Pmm (this is an extra action is added for the sake of stability)
			P[m:self.nmax+1,m]*=sect
			#compute next scaled sectorial
			sect=sect*self.Wnn[m+1]*sinTheta
		return P

	def derivTheta(self,cosTheta):
		dP=np.zeros([self.nmax+1,self.nmax+1])
		P=self.at(cosTheta)
		sinTheta=math.sqrt(1-cosTheta**2)
		for n in range(0,self.nmax+1):
			for m in range(0,n):
				dP[n][m]=(cosTheta*n*P[n][m]-math.sqrt((n**2-m**2)*(2*n+1)/(2*n-1))*P[n-1][m])/sinTheta
		return dP
            
#Some stuff for testing the SphereHarm stuff

#analytical solution of Pnm for degree 5 and order 2
def P52(theta):
	sc=math.sqrt(242550/860160)
	return sc*(2*math.cos(theta)+math.cos(3*theta)-3*math.cos(5*theta))

#analytical solution of dPnm/dtheta for degree 5 and order 2
def dP52(theta):
	sc=math.sqrt(242550/860160)
	return -sc*(2*math.sin(theta)+3*math.sin(3*theta)-15*math.sin(5*theta))

def sanity_check():
	#do a runtime import of packages (which are not necessary for the other SphereHarm routines)
	import matplotlib.pyplot as plt
	
	#make a vector of theta values (but do not include poles, to avoid divs by zero in the derivative)
	step=0.01
	theta=np.arange(step/2,math.pi-step/2,step)
	#maximum degree high enough to have the module sweating a bit
	
	nmax=20
	Pleg=Pnm(nmax)
		#def at(self,cosTheta):
		#return self.__compute__(cosTheta)[0:self.nmax,0:self.nmax]
	Pvals=np.zeros([nmax+1,nmax+1,len(theta)])
	dPvals=np.zeros([nmax+1,nmax+1,len(theta)])
	ith=0
	for th in theta:
		Pvals[:,:,ith]=Pleg.at(math.cos(th))
		dPvals[:,:,ith]=Pleg.derivTheta(math.cos(th))
		ith+=1
	
	#get analytical solution of test case (degree 5 order 2):
	Pana=[P52(th) for th in theta]	
	dPana=[dP52(th) for th in theta]

#plot analytical solution together with numerical solution for degree 5 order 2
	
	plt.plot(theta,Pvals[5,2,:],'b--')
	plt.plot(theta,Pana,'ro')
	plt.plot(theta,dPvals[5,2,:],'g--')
	plt.plot(theta,dPana,'co')
	plt.show()
	
	
    
