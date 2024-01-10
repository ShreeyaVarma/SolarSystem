from numpy import array
from solar_system import Planet, Sun, SolarSystem, SolarBody
import matplotlib.pyplot as plt
from math import sqrt
 
solar_system = SolarSystem(0.001)

def add_bodies_to_solar_system():
    global sun

    planets_data = read_data("planet_details.txt")

    SolarSystem.gravitational_const = eval(planets_data[0])
    sun = Sun(solar_system, eval(planets_data[1]))
   
    orbital_periods = [eval(period) for period in planets_data[2].split(",")]
    solar_system.save_orbital_periods(orbital_periods)

    for planet_data in planets_data[3: -2]:
            name, colour, mass, orbital_radius = planet_data.split(",")
            Planet(solar_system, name.strip(), eval(mass), colour.strip(), eval(orbital_radius))

    name, colour, mass, orbital_radius = planets_data[-2].split(",")
    earth_radius = eval(planets_data[-1])
    v = sqrt((SolarSystem.gravitational_const * 2) / earth_radius)
    SolarBody(solar_system, name, eval(mass), eval(orbital_radius), colour.strip(), array([1,0]), array([0, -v * 2.93]))

def read_data(filename):
    planet_file = open(filename, "r")
    temp_planets_data = planet_file.readlines()
    planets_data = []

    for data in temp_planets_data:
        if (data[0] != '#' and data.strip() != ""):
            planets_data.append(data)

    return planets_data


    planet_file.close()
def main():
    add_bodies_to_solar_system()
   
    solar_system.run_simulation(100000)

    solar_system.compare_orbital_periods()

    solar_system.plot_energies()

if __name__ == "__main__":
    main()
