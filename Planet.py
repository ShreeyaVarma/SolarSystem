import matplotlib
matplotlib.use('TKAgg')

import math
import numpy as np
from numpy.linalg import norm


'''
Planet class

'''


class Planet(object):
 
    def __init__(self, name, mass, orbit, colour):

        self.name = name
        # mass in kg
        self.m = mass
        # orbital radius in m
        self.orbit = orbit
        # colour - need to strip trailing line return!
        self.c = colour.strip()
        # set year to zero
        self.year = 0

    def initialise(self, G, p):
        # inital position, initial coords = (orbit radius, 0)
        self.r = np.array([self.orbit, 0])
        # inital velocity, tangential to position
        # speed = sqrt(G*marsmass/r)
        if (self.orbit == 0.0):
            self.v = np.array([0, 0])
        else:
            vel = math.sqrt(G*p.m/self.orbit)
            self.v = np.array([0, vel])
        # intial accelatation, using gravitational force law
        if (self.orbit == 0.0):
            self.a = np.array([0, 0])
        else:
            self.a = self.updateAcc(G, p)
        # set acc_old = acc to start Beeman
        self.a_old = self.a

        #print(f"{self.name}:: pos: ({self.r[0]}, {self.r[1]}), vel = ({self.v[0]}, {self.v[1]} ")

    def updatePos(self, dt):
        # keep old position to check for year
        self.r_old = self.r
       
        # update position first: Beeman
        self.r = self.r + self.v*dt + (4*self.a - self.a_old)*dt*dt/6.0
       

    def updateVel(self, G, dt, p):
        # update velocity second: Beeman
        a_new = self.updateAcc(G, p)
        self.v = self.v + (2*a_new + 5*self.a - self.a_old)*dt/6.0
        # now update acc ready for next iteration
        self.a_old = self.a
        self.a = a_new
       
       
    def updateAcc(self, G, p):
        # update acc (gravitational force law)
        pos = self.r - p.r
        a = -G*p.m*pos/math.pow(norm(pos),3)
        #print(f"Acc: {a[0], a[1]}")

        return a


    def newYear(self):
        # update the year when the planet passes the +x axis
        if (self.r_old[1] < 0.0 and self.r[1] >= 0.0):
            self.year +=1
            return True
        else:
            return False


    # determine kinetic energy
    def kineticEnergy(self):
        # ke in J
        ke = (np.dot(self.v, self.v))*self.m/2
        return ke
