from orbit_decay_num_sim_v5 import orbit
import numpy as np
from decimal import Decimal as D
import matplotlib.pyplot as plt
import pandas as pd

init_alt = 409000		# Initial altitude, m

r_earth = D(np.average([6378137, 6356752]))  	# Radius of Earth, m
G = D(6.67*10**-11)  							# Gravitational constant
mass_earth = D(5.972*10**24)  					# Earth mass, kg
mu = G*mass_earth  								# Specific gravitational constant of Earth
circ_vel = np.sqrt(mu/(r_earth + init_alt))		# Orbital velocity for circular orbit at init_alt

init_vel = circ_vel		# Initial velocity, m/s

C_D = 2.2				# Drag coefficient (2.2 is assumed constant by Jacchia 70)
t_step = 2000			# Time step, s


masses = [5.0]						# A range of reasonable masses for a 3U cubesat
area_refs = [100, 300, 350]			# A range of cross-sectional areas for a 3U cubesat, cm^2

dict = {}

for i in range(len(masses)):
	for j in range(len(area_refs)):
		dict_kword = str(int(masses[i])) + "kg_" + str(int(area_refs[j])) + "cm2"
		run_name = "deorbit_" + dict_kword
		
		mass_sat = masses[i]
		A_ref = area_refs[j]/10000

		time, altitude = orbit(init_alt, init_vel, C_D, A_ref, mass_sat, t_step, drag_on=True, circular_orbit=False, error_testing=False, num_orbits=100000)
		
		df = pd.DataFrame([time, altitude]).transpose()
		#df.columns = ['time', 'altitude', 'radius', 'radial_vel', 'true_anomoly', 'ang_vel']
		df.columns = ['time (s)', 'altitude (m)']
		df.to_csv(run_name+".csv", index=False)

		dict[dict_kword] = (time, altitude, mass_sat, A_ref)


fig_no = 0
colors = ['blue', 'green', 'red', 'orange']

for i in masses:
	fig_no += 1
	plt.figure(fig_no)
	kwords = []

	for j in dict.keys():
		if j.find(str(int(i))+"kg") == 0:
			kwords.append(j)

	for j in kwords:
		plt.plot(dict.get(j)[0]/(2400*24), dict.get(j)[1]/1000)

	str_title = "Satellite Altitude Decay for Varying Reference Areas\n" + str(i) + " kg CubeSat, Time Step = " + str(t_step) + " sec."
	plt.title(str_title)
	plt.xlabel("Time (days)")
	plt.ylabel("Altitude (km)")
	plt.legend([str(int(i)) + " cm$^2$" for i in area_refs], loc='upper right')

plt.show()