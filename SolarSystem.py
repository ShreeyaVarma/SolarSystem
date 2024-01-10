from ctypes import alignment
from numpy import array, dot
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sqrt, tan

class SolarSystem:
    gravitational_const = None
    furthest_revolving_body = None

    def __init__(self, timestep) -> None:
        self.bodies = []

        self.fig = plt.figure(0)
        self.axis = plt.axes()
        self.fig.tight_layout()

        self.timestep = timestep
        self.patches = []
        self.time_passed = 0

        self.ke_values = []
        self.pe_values = []
        self.total_e_values = []

        self.planetery_alignments = []

    def save_orbital_periods(self, orbital_periods):
        self.real_orbital_periods = orbital_periods

    def create_patch(self, body):
        patch_size = 0.15 if body.name.lower() == "sun" else min((body.mass/12 + 0.06) /2, 0.8)

        self.patches.append(self.axis.add_patch(plt.Circle(body.position, patch_size, color = body.colour, animated = True)))

    def add_body(self, body):
        self.bodies.append(body)

        if body.name.lower() == "sun":
            return

        if self.furthest_revolving_body == None:
            self.furthest_revolving_body = body

        elif body.orbital_radius > self.furthest_revolving_body.orbital_radius:
            self.furthest_revolving_body = body

    def calculate_energies(self):
        ke = 0.0
        pe = 0.0

        for idx, first_body in enumerate(self.bodies):
            ke += first_body.get_kinetic_energy()

            for second_body in self.bodies[idx+1:]:
                pe -= first_body.get_potential_energy(second_body, norm(first_body.position - second_body.position))

        return (ke, pe, ke + pe)

    def plot_energies(self):
        figure, axes = plt.subplots(3)

        time = [i * self.timestep for i in range(len(self.ke_values))]
        figure.tight_layout()

        [axes[i].set_xlabel("Time") for i in range(len(axes))]
        axes[0].set_ylabel("Kinetic Energy")
        axes[1].set_ylabel("Potential Energy")
        axes[2].set_ylabel("Total Energy")

        axes[0].plot(time, self.ke_values ,color="blue")
        axes[1].plot(time, self.pe_values, color="red")
        axes[2].plot(time, self.total_e_values, color="black")

        plt.show()

    def has_sat_reached_mars_earth(self):
        earth_idx = 0
        mars_idx = 0

        for idx, body in enumerate(self.bodies):
            if body.name.lower() == "mars":
                mars_idx = idx

            if body.name.lower() == "earth":
                earth_idx = idx

        distance = self.bodies[-1].position - self.bodies[mars_idx].position
        distance_mag = norm(distance)

        if distance_mag < self.sat_shortest_dis_to_mars and self.sat_shortest_dis_to_mars > 0.01:
            self.sat_shortest_dis_to_mars = distance_mag
            self.time_took_to_reach_mars = self.time_passed

        if self.sat_shortest_dis_to_mars <= 0.05:
            distance = self.bodies[-1].position - self.bodies[earth_idx].position
            distance_mag = norm(distance)

            if distance_mag < self.sat_shortest_dis_to_earth:
                self.sat_shortest_dis_to_earth = distance_mag

    def animate(self, i):
        self.time_passed += self.timestep

        for idx, body in enumerate(self.bodies):
            old_pos = body.position
            body.update_position(self.timestep, direct_euler = True)

            if(body.name.lower() != "sun" and str(type(body)) != "<class 'solar_system.SolarBody'>"):
                body.time_passed_in_year += self.timestep
                if body.has_year_passed(old_pos):
                    body.years_passed += 1

                    body.orbital_periods.append(body.time_passed_in_year)
                    body.time_passed_to_year = 0
                    print(f"{body.name}: Year Passed: {body.years_passed}\tEarth Time: {self.time_passed}")

            self.patches[idx].center = body.position

        for first_body in self.bodies:
            for second_body in self.bodies:
                if first_body != second_body:
                    first_body.update_velocity(self.timestep, second_body, direct_euler = True)

        self.has_sat_reached_mars_earth()

        if self.check_planetery_alignment():
            self.planetery_alignments.append(self.time_passed)

        energies = self.calculate_energies()
        self.ke_values.append(energies[0])
        self.pe_values.append(energies[1])
        self.total_e_values.append(energies[2])

        to_joules_const = (5.97219e+24 * 1.496e+11 * 1.496e+11) / (3.154e+7 * 3.154e+7)
        self.file.write(f"Total Energy: {round((energies[2] * to_joules_const), ndigits=3)} J\n")
       
        return self.patches

    def animate_only_sun_g(self, i):
        self.time_passed += self.timestep

        for idx, body in enumerate(self.bodies):
            old_pos = body.position
            body.update_position(self.timestep, direct_euler = False)

            if(body.name.lower() != "sun" and str(type(body)) != "<class 'solar_system.SolarBody'>"):
                body.time_passed_to_year += self.timestep
                if body.has_year_passed(old_pos):
                    body.years_passed += 1

                    body.orbital_periods.append(body.time_passed_to_year)
                    body.time_passed_to_year = 0
                    print(f"{body.name}: Year Passed: {body.years_passed}\tEarth Time: {self.time_passed}")

            self.patches[idx].center = body.position


        for first_body in self.bodies:
            for second_body in self.bodies:
                if (first_body != second_body) and (second_body.name.lower() == "sun"):
                    first_body.update_velocity(self.timestep, second_body, direct_euler = True)
                   
        self.has_sat_reached_mars_earth()

        energies = self.calculate_energies()
        self.ke_values.append(energies[0])
        self.pe_values.append(energies[1])
        self.total_e_values.append(energies[2])

        to_joules_const = (5.97219e+24 * 1.496e+11 * 1.496e+11) / (3.154e+7 * 3.154e+7)
        self.file.write(f"Total Energy: {round((energies[2] * to_joules_const), ndigits=3)} J\n")
       
        return self.patches

    def run_simulation(self, frames, only_sun_g =  False):
        sat_initial_vel = norm(self.bodies[-1].velocity)
        ax_lim = (-self.furthest_revolving_body.orbital_radius - 0.5, self.furthest_revolving_body.orbital_radius + 0.5)

        self.axis.set_xlim(*ax_lim)
        self.axis.set_ylim(*ax_lim)

        self.file = open("Total Energy.txt", 'w')

        self.sat_shortest_dis_to_mars = ax_lim[1]
        self.sat_shortest_dis_to_earth = ax_lim[1]
        self.time_took_to_reach_mars = 0

        if only_sun_g:
            animation = FuncAnimation(self.fig, self.animate_only_sun_g, frames = frames, repeat = False, interval = self.timestep, blit= True)

        else:
            animation = FuncAnimation(self.fig, self.animate, frames = frames, repeat = False, interval = self.timestep, blit= True)

        plt.show()

        self.file.close()

        file = open("Solar Output.txt", 'w')
        file.write(f"WITH INITIAL VELOCITY OF {-sat_initial_vel} units\n")
        file.write(f"Shortest Distance Between Satellite and Mars was: {self.sat_shortest_dis_to_mars}\n")

        time_to_mars = (int)(12*self.time_took_to_reach_mars)
        file.write(f"This took around {time_to_mars} months which is {time_to_mars - 7} more months then NASA's Perseverance\n")

        file.write(f"After going to marks Satelite manage to come close to earth with shortest distance of {self.sat_shortest_dis_to_earth}\n")

        file.writelines(self.compare_orbital_periods())

        for alignment_time in self.planetery_alignments:
            file.write(f"Planets Alinged after {alignment_time} earth year\n")

        file.close()

        print("\n OUTPUT IS WRITTEN OUT TO A TEXT FILE CALLED Solar Output.txt")

    def compare_orbital_periods(self):
        uncertainties = []
        for idx, body in enumerate(self.bodies[1:]):
            if len(body.orbital_periods) > 0:
                calculated_orbital_period = (sum(body.orbital_periods) / len(body.orbital_periods))
                uncertaintity = ((calculated_orbital_period - self.real_orbital_periods[idx]) / self.real_orbital_periods[idx]) * 100

                uncertainties.append(f"{body.name} Orbital Period Uncertaintity = {round(uncertaintity, 5)}%\n")

        return uncertainties

    def check_planetery_alignment(self):
        angles = []
        for body in self.bodies:
            if str(type(body)) != "<class 'solar_system.SolarBody'>":
                angles.append(tan( body.position[1] / body.position[0] ))

        for idx, angle in enumerate(angles[:-1]):
            if abs(angle - angles[idx + 1]) > 5:
                return False

        return True

class SolarBody:
    def __init__(self, solar_system,name, mass, orbital_radius, colour, position=array([0, 0]), velocity= array([0, 0]))-> None:
        self.solar_system = solar_system
        self.name = name
        self.mass = mass
        self.orbital_radius = orbital_radius

        self.colour = colour
        self.position = position
        self.velocity = velocity
        self.acceleration = array([0,0])
        self.old_acceleration = self.acceleration
        self.years_passed = 0
        self.orbital_periods = []

        self.solar_system.add_body(self)
        self.solar_system.create_patch(self)

    def update_position(self, timestep, direct_euler = False):
       
        if(not direct_euler):
            self.position = self. position + self.velocity * timestep + 1/6*(4*self.acceleration - self.old_acceleration)*(timestep**2)

        else:
            self.position = self.position + self.velocity * timestep

    def update_velocity(self, timestep, other_body, direct_euler = False):
        if(not direct_euler):
            new_acceleration = self.update_acceleartion(other_body)

            self.velocity = self.velocity + 1/6*(2*new_acceleration + 5*self.acceleration - self.old_acceleration)*timestep

            self.old_acceleration = self.acceleration
            self.acceleration = new_acceleration

        else:
            self.velocity = self.velocity + self.acceleration * timestep
            self.acceleration = self.update_acceleartion(other_body)

    def update_acceleartion(self, otherBody):        
        distance = otherBody.position - self.position
        distance_mag = norm(distance)

        acceleration = ((SolarSystem.gravitational_const * otherBody.mass ) / (distance_mag ** 3)) * distance

        return acceleration

    def get_kinetic_energy(self):
        return (1/2*self.mass) * (dot(self.velocity, self.velocity))

    def get_potential_energy(self, bigger_body, distance_btw_bodies):
        return (SolarSystem.gravitational_const * bigger_body.mass * self.mass) / distance_btw_bodies

    def has_year_passed(self, old_position):
        return (old_position[1] < 0 and self.position[1] > 0 and self.position[0] >= 0)
       
class Planet(SolarBody):
    def __init__(self, solar_system, name, mass, colour, orbital_radius) -> None:
        vel = self._calculate_initial_vel(solar_system.bodies[0].mass, orbital_radius)
        super().__init__(solar_system, name, mass,orbital_radius, colour, array([orbital_radius, 0]), velocity=vel)

        acc = self._calculate_initial_acc(solar_system.bodies[0], orbital_radius)
        self.acceleration = acc
        self.old_acceleration = acc

        self.time_passed_in_year = 0

    def _calculate_initial_vel(self, bigger_mass, distance_btw_bodies):
        if distance_btw_bodies == 0:
            return array([0,0])

        vel = sqrt((SolarSystem.gravitational_const * bigger_mass)/ distance_btw_bodies)
        return array([0,vel ])

    def _calculate_initial_acc(self, bigger_body, distance_btw_bodies):
        if distance_btw_bodies == 0:
            return array([0,0])

        return self.update_acceleartion(bigger_body)

class Sun(SolarBody):
    def __init__(self, solar_system, mass, colour= "yellow", position=array([0, 0]), velocity=array([0, 0])) -> None:
        super().__init__(solar_system, "Sun",mass,0, colour, position, velocity)
