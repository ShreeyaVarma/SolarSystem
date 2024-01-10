from Planet import Planet

import matplotlib
matplotlib.use('TKAgg')

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy.linalg import norm


'''
Class to run the orbital simulation

'''


class Solar(object):

    def __init__(self):

        inputdata = []
        #filename = raw_input("Input file: ")
        filename = "parameters-solar.dat"
        filein = open(filename, "r")
        for line in filein.readlines():
            if (not line.startswith("#")):
                inputdata.append(line)
        filein.close()

        # simulation parameters
        self.niter = int(inputdata[0])
        self.dt = float(inputdata[1])
        self.G = float(inputdata[2])
        print(self.G)
        # list for mars and moons
        self.bodies = []
       
        # rest of input data is mars and moon data in four line 'chunks'
        # first entry must be mars
        for i in range(3, len(inputdata)-4, 4):
            name = inputdata[i]
            mass = float(inputdata[i+1])
            orbit = float(inputdata[i+2])
            colour = inputdata[i+3]
            self.bodies.append(Planet(name, mass, orbit, colour))
                               
        # set initial positions and velocities relative to sun
        # sum must be first element in bodies list!
        for i in range(0, len(self.bodies)):
            self.bodies[i].initialise(self.G, self.bodies[0])



    def init(self):
        # initialiser for animator
        return self.patches


    def animate(self, i):
        # keep track of time in earth years
        time = (i+1)*self.dt

        # update positions
        for j in range(0, len(self.bodies)):
            self.bodies[j].updatePos( self.dt)
            self.patches[j].center = self.bodies[j].r
           
        # then update velocities
        for j in range(0, len(self.bodies)):
            for k in range(0, len(self.bodies)):
                if (j != k):
                    self.bodies[j].updateVel(self.G, self.dt, self.bodies[k])

        # check year and print year if new year for any planet
        for j in range(0, len(self.bodies)):
            if (self.bodies[j].newYear()):
               print(self.bodies[j].name.strip() + " " + str(self.bodies[j].year) + " years = " + str(time) + " earth years")
               # in new year is earth year, also print total energy

               if (self.bodies[j].name.strip() == 'earth'):
                   # need to convert from earth masses AU^2 yr^-2 to kg m^2 s-2 (J)
                   # 1 earth mass = 5.97219e24 kg
                   # 1 AU = 1.496e+11 m
                   c =(5.97219e+24*1.496e+11*1.496e+11)/(3.154e+7*3.154e+7)
                   energy = self.energy()*c
                   self.time.append(time)
                   print('Time = ' + str(time) + ' earth years. Total energy = ' + '{:.3}'.format(energy) + ' J')
               

        return self.patches

    def energy(self):
        ke = 0.0
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if (k != j):
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double countin
        pe = pe / 2
        totEnergy = ke + pe
        print(f"Kinetic Energy: {ke}, Potential Energy: {pe}, Total Energy: {totEnergy}")

        self.ke.append(ke)
        self.pe.append(pe)
        self.tot_e.append(totEnergy)

        return totEnergy

    def calcTotalEnergy(self, i):
        ke = 0.0
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if (k != j):
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double countin
        pe = pe / 2
        totEnergy = ke + pe
        print('Time = ' + str(i) + ' iterations. Total energy = ' + '{:.3e}'.format(totEnergy))

    def plot_energies(self):
        figure, axes = plt.subplots(3)

        figure.tight_layout()

        print(self.tot_e)

        axes[0].plot(self.time, self.ke ,color="blue")
        axes[1].plot(self.time, self.pe, color="red")
        axes[2].plot(self.time, self.tot_e, color="black")

        plt.show()

    def run(self):

        # set up the plot components        
        fig = plt.figure()
        ax = plt.axes()

        # create an array for patches (planet and moons)
        self.patches = []
        self.ke = []
        self.pe = []
        self.tot_e = []
        self.time = []

        # get orbital radius of outermost moon to set size of orbiting bodies and of plot    
        # hacky - should really check to see which moon is outermost
        maxOrb = math.sqrt(np.dot(self.bodies[-1].r, self.bodies[-1].r))

        # add the planet and moons to the Axes and patches
        for i in range(0, len(self.bodies)):
            if (i == 0):
                self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.05*maxOrb, color = self.bodies[i].c, animated = True)))
            else:
                self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.02*maxOrb, color = self.bodies[i].c, animated = True)))
       
        # set up the axes
        # scale axes so circle looks like a circle and set limits with border b for prettier plot
        b = 1.2
        lim = maxOrb*b
        print(lim)
        ax.axis('scaled')
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
               
        anim = FuncAnimation(fig, self.animate, init_func = self.init, frames = self.niter, repeat = False, interval = 0.01, blit= True)

        plt.show()
