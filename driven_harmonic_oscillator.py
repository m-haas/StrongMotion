#Returns displacement for a driven and damped harmonic oscillator - Linear System
#Where driving force is a seismogram
#Uses Newmark's method with linear acceleration method

class Newmark:
    '''
    Newmark solver for a SDOF system
    '''

    def __init__(self,m,xi,k):
        '''
        Given mass (m), damping ratio (xi) and stiffness
        sets up the SDOF system
        '''
        #mass
        self.m = m
        #damping ratio
        #self.xi = xi
        #stiffness
        self.k = k
        #SDOF parameters
        #c from xi
        self.c=xi*2.*(k*m)**0.5

    def set_init_conditions(self,u_0,ut_0,gamma,beta):
        '''
        Initial conditions for Newmark solver
        and assuming initial conditions for
        deformation (u0) and velocity (ut0)
        and solving approach (gamma) and (beta)
        '''
        self.u_0   = u_0
        self.ut_0  = ut_0
        self.gamma = gamma
        self.beta  = beta

    def groundmotion_response(self,p,dt):
        '''
        Linear system ground motion response
        for a force from an acceleration record (p)
        with samples with given time spacing (dt)
        NOTE: ONLY LINEAR SYSTEM i.e. fs = ku
        (Chopra p167)
        '''
        #initial acceleration
        utt_0 = (p[0] - self.c*self.ut_0 - self.k*self.u_0)/self.m

        #vectors for system (set to initial value)
        u   = [0  for i in range(len(p))]
        ut  = u
        utt = u
        u[0]   = self.u_0
        ut[0]  = self.ut_0
        utt[0] = utt_0

        #factors
        a0 = self.gamma/(self.beta*dt)
        a1 = 1./(self.beta*dt**2)
        a2 = 1./(self.beta*dt)
        a3 = self.gamma/self.beta
        a4 = 1./(2.*self.beta)
        a5 = dt*(self.gamma/(2.*self.beta) - 1)
        a6 = self.gamma*a2
        a7 = dt*(1-self.gamma/(2.*self.beta))

        #khat
        khat = self.k + a0*self.c + a1*self.m
        #a and b
        a = a2*self.m + a3*self.c
        b = a4*self.m + a5*self.c

        #go through time steps (except last one since we always update next one)
        for i in range(len(p)-1):
            #get incremental approximations
            #incremental gm forcing
            dphati = p[i+1] - p[i] + a*ut[i] + b*utt[i]
            #displacement increment
            du = dphati/khat
            #velocity increment
            dut = a6*du - a3*ut[i] + a7*utt[i]
            #acceleration increment
            dutt = a1*du - a2*ut[i] - a4*utt[i]

            #update displacement, velocity and acceleration
            u[i+1]   = u[i]   + du
            ut[i+1]  = ut[i]  + dut
            utt[i+1] = utt[i] + dutt

            #print table
            print "{:2.1f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f}".format(dt*i,p[i],utt[i],p[i+1]-p[i],dphati,du,dut,dutt,ut[i],u[i])

        return [u,ut,utt]

##mass [kg]
#m = 91153.0644091
##damping ratio [-]
#xi = 0.05
##stiffness [N/m]
#k = 3090534.95577
#
#solver = Newmark(m,xi,k)
#
##initial conditions
#u_0 = 0;
#ut_0 = 0;
#
##linear acceleration method
#gamma = 1./2
#beta = 1./6


#def dY_dt(Y, t=0):
#    """
#    Return the derivative
#    """
#    return numpy.array([ Y[1], -omega**2*numpy.sin(Y[0])-b*Y[1]+F0*numpy.sin(omega2*t)])
#
#from scipy import integrate
#
#t = numpy.linspace(0,50,20000)
#omega = 1
#omega2=0.6
#b=0.5 #damping
#
#F0=1.0 #forcing
#
#Y0=numpy.array([0.0, 0]) #starting conditions
#Y = integrate.odeint(dY_dt, Y0, t)
#
#plt.plot(t,Y[:,0],label="F0=0.0")
#plt.legend(loc="upper left")
#plt.xlabel("Time (sec)")
#plt.ylabel("Angle (rad)")

class CentralDifferenceMethod():
    """
    Central Difference Method to solve Equation of Motion
    """
    def __init__(self,m,xi,k):
        '''
        Given mass (m), damping ratio (xi) and stiffness
        sets up the SDOF system
        '''
        #mass
        self.m = m
        #damping ratio
        #self.xi = xi
        #stiffness
        self.k = k
        #SDOF parameters
        #c from xi
        self.c=xi*2.*(k*m)**0.5

    def groundmotion_response(self,p,dt,u_0,ut_0):
        '''
        Linear system ground motion response
        for a force from an acceleration record (p)
        with samples with given time spacing (dt)
        and assuming initial conditions for
        deformation (u0) and velocity (ut0)
        '''
        #initial acceleration
        utt_0 = (p[0] - self.c*ut_0 - self.k*u_0)/self.m
        uprev = u_0 - dt*ut_0 + (dt**2.)/2.*utt_0

        #vectors for system (set to initial value)
        u   = [uprev] + [0.   for i in range(len(p))]
        ut  = [0.]    + [0.  for i in range(len(p))]
        utt = [0.]    + [0. for i in range(len(p))]
        u[1] = u_0
        u[1] = ut_0
        utt[1] = utt_0

        #secant stiffness
        khat = self.m/(dt**2.) + self.c/(2.*dt)
        #a and b
        a = self.m/(dt**2.) - self.c/(2*dt)
        b = self.k - 2.*self.m/(dt**2.)

        #go through time steps
        #NOTE: u is shifted one up
        #      --> shift gm one up
        #          start at 1 stop at n-1)
        p = [0.]+list(p)
        for i in range(1,len(p)-1):
            #get incremental approximations
            #central gm forcing: pi - au[i-1] - bu[i]
            phati = p[i] - a*u[i-1] - b*u[i]
            #u[i+1]
            u[i+1] = phati/khat
            #ut[i]
            ut[i] = (u[i+1]-u[i-1])/(2.*dt)
            #utt[i]
            utt[i] = (u[i+1]-2.*u[i]+u[i-1])/(dt**2.)

            #print table
            print "{:2.1f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f}".format(dt*(i-1),p[i],u[i-1],u[i],phati,u[i+1])

        return [u[1:],ut[1:],utt[1:]]
