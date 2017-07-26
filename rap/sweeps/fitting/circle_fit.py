def circle_fit(self, Show_Plot = True):

	S21 = self.loop.z
	Freq = self.loop.freq

	LargeCircle = 10
	def pythag(m,n):
		'''compute pythagorean distance
		   sqrt(m*m + n*n)'''
		return np.sqrt(np.square(m) + np.square(n))

	def eigen2x2(a,b,c):
		'''a,b,c - matrix components 	[[a c]
										 [c d]]	
		   d1,d2 - eigen values where |d1| >= |d2|
		   (Vx,Vy) - unit eigen vector of d1,  Note: (-Vy,Vx) is eigen vector for d2
		'''
		disc = pythag(a-b,2*c) # discriminant
		d1 = max(a+b + disc, a+b - disc)/2
		d2 = (a*b-c*c)/d1

		if np.abs(a-d1) > np.abs(b-d1):
			f = pythag(c,d1-a)
			if f == 0.0:
				Vx = 1.
				Vy = 0.
			else:
				Vx = c/f
				Vy = (d1-a)/f
		else:
			f = pythag(c,d1-b)
			if f == 0.0:
				Vx = 1.
				Vy = 0.
			else:
				Vx = (d1-b)/f
				Vy = c/f					
		return d1,d2,Vx,Vy

	def F(x,y,a,b):
		''' computes and returns the value of the objective fuction.
		do this for the case of a large circle and a small circle  '''
		
		if (np.abs(a) < LargeCircle) and (np.abs(b) < LargeCircle): # Case of Small circle
			xx = x - a
			yy = y - b
			D = pythag(xx,yy)

			r = D.mean()

			return (np.square(D - r)).mean()
		else:	# Case of Large circle
			a0 = a - x.mean()
			b0 = b - y.mean()
			d = 1.0/pythag(a0,b0)
			dd = d*d
			s = b0*d
			c = a0*d

			xx = x - x.mean()
			yy = y - y.mean()
			z = np.square(xx) + np.square(yy)
			p = xx*c + yy*s
			t = d*z - 2.0*p
			g = t/(1+np.sqrt(1.+d*t))
			W = (z+p*g)/(2.0+d*g)
			Z = z

			return Z.mean() - W.mean()*(2.0+d*d*W.mean())

	def GradHessF(x,y,a,b):
		'''Compute gradient of F, GradF = [F1,F2] and Hessian of F, HessF = [[A11 A12]
																			  A12 A22]]
		at point p = [a,b].
		Note Hessian is symmetric. 
		'''
		if (np.abs(a) < LargeCircle) and (np.abs(b) < LargeCircle): # Case of Small circle
			xx = x - a
			yy = y - b
			r = pythag(xx,yy)
			u = xx/r
			v = yy/r

			Mr = r.mean()
			Mu = u.mean()
			Mv = v.mean()
			Muu = (u*u).mean()
			Mvv = (v*v).mean()
			Muv = (u*v).mean()
			Muur = (u*u/r).mean()
			Mvvr = (v*v/r).mean()
			Muvr = (u*v/r).mean()
		


			F1 = a + Mu * Mr - x.mean()
			F2 = b + Mv * Mr - y.mean()

			A11 = 1.0 - Mu * Mu - Mr * Mvvr
			A22 = 1.0 - Mv * Mv - Mr * Muur
			A12 = -1.0 * Mu * Mv + Mr * Muvr

		else:	# Case of Large circle
			a0 = a - x.mean()
			b0 = b - y.mean()
			d = 1.0/pythag(a0,b0)
			dd = d*d
			s = b0*d
			c = a0*d

			xx = x - x.mean()
			yy = y - y.mean()
			z = np.square(xx) + np.square(yy)
			p = xx*c + yy*s
			t = 2.0*p - d*z 
			w = np.sqrt(1.0-d*t)				
			g = -1.0*t/(1.0+w)
			g1 = 1.0/(1.0+d*g)
			gg1 = g*g1
			gg2 = g/(2.0 + d * g)
			aa = (xx+g*c)/w
			bb = (yy+g*s)/w	

			X = (xx*gg1).mean()
			Y = (yy*gg1).mean()
			R = (z+t*gg2).mean()
			T = (t*gg1).mean()
			W = (t*gg1*gg2).mean()	
			AA = (aa*aa*g1).mean()
			BB = (bb*bb*g1).mean()
			AB = (aa*bb*g1).mean()
			AG = (aa*gg1).mean()
			BG = (bb*gg1).mean()
			GG = (g*gg1).mean()	

			U = (T-b*W)*c*0.5 - X + R*c*0.5
			V = (T-b*W)*s*0.5 - Y + R*s*0.5

			F1 = d * ((dd*R*U - d*W*c + T*c)*0.5 - X)
			F2 = d * ((dd*R*V - d*W*s + T*s)*0.5 - Y)

			UUR = ((GG-R*0.5)*c + 2.0*(AG-U))*c + AA
			VVR = ((GG-R*0.5)*s + 2.0*(BG-V))*s + BB
			UVR = ((GG-R*0.5)*c + (AG-U))*s + (BG-V)*c + AB

			A11 = dd*(U*(2.0*c-dd*U) - R*s*s*0.5 - VVR*(1.0+dd*R*0.5))
			A22 = dd*(V*(2.0*s-dd*V) - R*c*c*0.5 - UUR*(1.0+dd*R*0.5))
			A12 = dd*(U*s + V*c + R*s*c*0.5 - dd*U*V + UVR*(1.0 + dd*R*0.5))
		return F1,F2,A11,A22,A12
	
	def sigma(x,y,loop):
		'''estimate of Sigma = square root of RSS divided by N
		gives the root-mean-square error of the geometric circle fit'''
		dx = x-loop.a
		dy = x-loop.b
		loop.sigma = (pythag(dx,dy)-loop.r).mean()
		return loop

	def CircleFitByChernovHoussam(x,y, init, lambda_init):
		import copy
		import sys


		REAL_EPSILON = sys.float_info.epsilon
		REAL_MAX = sys.float_info.max


		IterMAX=200
		check_line= True
		#dmin = 1.0

		ParLimit2 = 100.
		epsilon = 1.e+7*REAL_EPSILON
		factor1 = 32.
		factor2 = 32.
		ccc = 0.4
		factorUp = 10.
		factorDown = 0.1

		new = copy.copy(init)
		#new = sys.modules[__name__].loop() #This is how to access the loop class from inside this function
		#old = loop()

		new.s = F(x,y,init.a,init.b) # compute root mean square error
		F1,F2,A11,A22,A12 = GradHessF(x,y,init.a,init.b) # compute gradient vector and Hessian matrix
		new.Gx = F1
		new.Gy = F2
		new.g = pythag(F1,F2) # The gradient vector and its norm
		lambda_ = lambda_init
		sBest = gBest = progess = REAL_MAX

		enough = False
		i = 0
		ii = 0
		while not enough:
			if i > 0:
				# evaluate the progress made during the previous iteration
				progress = (np.abs(new.a - old.a)+np.abs(new.b - old.b))/(np.square(old.a) + np.square(old.b) + 1.0)
			old = copy.copy(new)

			i = i+1
			if i > IterMAX: #termination due to going over the limit
				enough = True
				break
			d1,d2,Vx,Vy = eigen2x2(A11,A22,A12) #eigendecomposition of the Hessian matrix
			dmin = min(d1,d2) #recording the smaller e-value
			AB = pythag(old.a,old.b) + 1.0 # approximation to the circle size
			# main stopping rule: terminate iterations if 
			# the gradient vector is small enough and the 
			# progress is not substantial 
			if (old.g < factor1*REAL_EPSILON) and (progress<epsilon):
				#print('primary stopping rule')
				enough = True
				break
			# secondary stopping rule (prevents some stupid cycling)
			if (old.s >= sBest) and (old.g >= gBest):
				print(old.s, sBest, old.g, gBest)
				#print('secondary stopping rule')
				enough = True
				break

			if (sBest > old.s):
				sBest = old.s  # updating the smallest value of the objective function found so far
			if (gBest > old.g): 
				gBest = old.g  # updating the smallest length of the gradient vector found so far

			G1 = Vx*F1 + Vy*F2  # rotating the gradient vector
			G2 = Vx*F2 - Vy*F1  # (expressing it in the eigensystem of the Hessian matrix)

			while not enough: # starting point of an "inner" iteration (adjusting lambda)
				# enforcing a lower bound on lambda that guarantees that
				# (i)  the augmented Hessian matrix is positive definite
				# (ii) the step is not too big (does not exceed a certain 
				# fraction of the circle size) the fraction is defined by 
				# the factor "ccc")
				if lambda_ < (np.abs(G1)/AB/ccc) - d1:
					lambda_ = np.abs(G1)/AB/ccc - d1
				if lambda_ < (np.abs(G2)/AB/ccc) - d2: 
					lambda_ = np.abs(G2)/AB/ccc - d2

				# compute the step (dX,dY) by using the current value of lambda
				dX = old.Gx*(Vx*Vx/(d1+lambda_)+Vy*Vy/(d2+lambda_)) + old.Gy*Vx*Vy*(1.0/(d1+lambda_)-1.0/(d2+lambda_))
				dY = old.Gx*Vx*Vy*(1.0/(d1+lambda_)-1.0/(d2+lambda_)) + old.Gy*(Vx*Vx/(d2+lambda_)+Vy*Vy/(d1+lambda_))

				# updating the loop parameter
				new.a = old.a - dX
				new.b = old.b - dY

				if (new.a==old.a) and (new.b==old.b): #if no change, terminate iterations
					enough  = True
					break

				#check if the circle is very large
				if np.abs(new.a)>ParLimit2 or np.abs(new.b)>ParLimit2:
					#when the circle is very large for the first time, check if 
					#the best fitting line gives the best fit

					if check_line:   # initially, check_line= True, then it is set to zero

						#compute scatter matrix
						dx = x - x.mean()
						dy = y - y.mean()
						Mxx = (dx*dx).sum()
						Myy = (dy*dy).sum()
						Mxy = (dy*dx).sum()
						dL1,dL2,VLx,VLy = eigen2x2(Mxx,Myy,Mxy)  # eigendecomposition of scatter matrix

						#compute the third mixed moment (after rotation of coordinates)
						dx = (x - x.mean())*VLx + (y - y.mean())*VLy
						dy = (y - y.mean())*VLx - (x - x.mean())*VLy
						Mxxy = (dx*dx*dy).sum()

						#rough estimate of the center to be used later to recoved from the wrong valley
						if Mxxy > 0.0:
							R = ParLimit2
						else:
							R = -ParLimit2

						aL = -VLy*R
						bL =  VLx*R                 
						check_line = False

					# check if the circle is in the wrong valley
					if (new.a*VLy - new.b*VLx)*R>0.0: 
						# switch to the rough circle estimate (precomupted earlier)
						new.a = aL;                 
						new.b = bL;                 
						new.s = F(x,y,new.a,new.b)    # compute the root-mean-square error
						
						# compute the gradient vector and Hessian matrix
						F1,F2,A11,A22,A12 = GradHessF(x,y,new.a,new.b)  

						# the gradient vector and its norm 
						new.Gx = F1;  
						new.Gy = F2;   
						new.g = pythag(F1,F2)  
						lambda_ = lambda_init     #reset lambda
						sBest = gBest = REAL_MAX  #reset best circle characteristics 
						break
				
				# compute the root-mean-square error
				new.s = F(x,y,new.a,new.b) 
				# compute the gradient vector and Hessian matrix
				F1,F2,A11,A22,A12 = GradHessF(x,y,new.a,new.b)

				# the gradient vector and its norm  
				new.Gx = F1  
				new.Gy = F2   
				new.g = pythag(F1,F2) 

				# check if improvement is gained
				if new.s < sBest*(1.0+factor2*REAL_EPSILON):  #yes, improvement
					lambda_ *= factorDown     # reduce lambda
					break 
				else:
					ii += 1
					if ii > IterMAX: #termination due to going over the limit
						enough = True
						break
					lambda_ *= factorUp #increace lambda
					continue
		

		old.r = pythag(x - old.a, y - old.b).mean() 
		old.outer_iterations = i
		old.inner_iterations = ii
		loop = old
		exit_code = 0
		if old.outer_iterations  > IterMAX:
			exit_code  = 1

		if old.inner_iterations  > IterMAX:
			exit_code = 2

		if (dmin <= 0.0) and (exit_code==0):
			exit_code  = 3

		loop.circle_fit_exit_code = exit_code
		loop = sigma(x,y,loop)

		return loop
	



	x = S21.real
	y = S21.imag


	self.loop.a =  0#guess.real#0
	self.loop.b =  0#guess.imag #0
	lambda_init = 0.001
	#self.loop = CircleFitByChernovHoussam(x,y, self.loop, lambda_init)
	if True: #self.loop.circle_fit_exit_code != 0:
		#print('Circle Fit Failed! Trying again...')
		#another initial guess
		norm = np.abs(S21[1:5].mean())
		S21 = S21/norm
		guess = np.mean(S21)
		self.loop.a =  guess.real#0
		self.loop.b =  guess.imag #0
		lambda_init = 0.001
		x = S21.real
		y = S21.imag
		self.loop = CircleFitByChernovHoussam(x,y, self.loop, lambda_init)
		self.loop.a = self.loop.a*norm
		self.loop.b = self.loop.b*norm
		self.loop.r = self.loop.r*norm
		self.loop.z = S21*norm

		if self.loop.circle_fit_exit_code != 0:
			print('!!!!!!!!!!!!!!    Circle Fit Failed Again! Giving Up...')

	if Show_Plot:
		fig, ax = self.plot_loop(show = False)[:2]		
		t = np.linspace(0, 2.0*np.pi, num=50, endpoint=True)
		j = np.complex(0,1); zc = self.loop.a + j*self.loop.b;  r = self.loop.r
		line = ax.plot(zc.real + r*np.cos(t),zc.imag + r*np.sin(t),'y-', label = 'Circle Fit')
		line = ax.plot([zc.real],[zc.imag],'yx', markersize = 10, markeredgewidth = 4, label = 'Center')
		ax.set_aspect('equal')
		plt.show()