#python class which computes 4-pi Geodesy style normalized assocaiited Legendre functions and its derivative
#copyright R. Rietbroek
#Date: 20 June 2016
#License: pending opensource  license to be determined, all rights reserved
#todo: Performance (array accessing order) and stability checks
import numpy as np
import math
class Pnm:
	#Constructor precomputes some stuff woth square roots  based on the requested maximum degree (
	def  __init__(self,nmax):
		self.nmax=nmax+1 #note we add an extra degree to also be able to solve for the derivative
		self.Wnn=np.zeros(self.nmax+1)
		self.Wnm=np.zeros([self.nmax+1,self.nmax+1])
		self.Wnn[1]=math.sqrt(3)
		for n in range(2,self.nmax+1):
			self.Wnn[n]=math.sqrt((2*n+1)/(2.0*n))

		#also precompute  Wnm array
		for n in range(1,self.nmax+1):
			for m in range(0,n):
				self.Wnm[n][m]=math.sqrt((2*n+1.0)/(n+m)*(2*n-1.0)/(n-m))


	def at(self,cosTheta):
		return self.__compute__(cosTheta)[0:self.nmax,0:self.nmax]
    #compute the value of the associated legendre function at x (0..1) usually x is cos(colatitude)
	def __compute__(self,cosTheta):
		assert cosTheta>=0 and cosTheta<=1,"Error in Pnm.at(x): 0 <= x <= 1 not fulfilled"
		P=np.zeros([self.nmax+1,self.nmax+1])
		sinTheta=math.sqrt(1-cosTheta**2)
		#fill up diagonal
		P[0][0]=1
		for n in range(1,self.nmax+1):
			#print(self.Wnn[n],n)
			P[n][n]=self.Wnn[n]*sinTheta*P[n-1][n-1]
		#also compute the off-diagonal elements (apply recursion for fixed order, varying degree)
		
		for m in range(0,self.nmax):
			#first value just off the diagonal
			n=m+1
			P[n][m]=self.Wnm[n][m]*cosTheta*P[n-1][m]
			##proceed recursion with 2 terms
			for n in range(m+2,self.nmax+1):
				P[n][m]=self.Wnm[n][m]*(cosTheta*P[n-1][m]-P[n-2][m]/self.Wnm[n-1][m])
		return P

	def derivTheta(self,cosTheta):
		dP=np.zeros([self.nmax,self.nmax])
		P=self.__compute__(cosTheta)
		sinTheta=math.sqrt(1-cosTheta**2)
		for n in range(1,self.nmax):
			#no need to do the sectorials though
			for m in range(0,n):
				dP[n][m]=(cosTheta*(n+1)*P[n][m]-math.sqrt((n+m+1)*(n-m+1))*math.sqrt((2*n+1)/(2*n+3))*P[n+1][m])/sinTheta
		return dP
