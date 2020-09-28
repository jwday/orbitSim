from orbit_decay_num_sim_v5 import orbit
import numpy as np
from decimal import Decimal as D
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime

# This code assumes you start at 409 km and will not be going higher. Therefore, for vis-viva calculations, other orbital altitudes will be taken at perigee

start_time = datetime(2013, 3, 31, 12)		# Unless otherwise specified, start time is as specified 
alt_pe = [408250, 407000, 349000]			# Altitudes at perigee, m
r_earth = D(6378100)  						# Radius of Earth, m
G = D(6.67*10**-11)  						# Gravitational constant
mass_earth = D(5.972*10**24)  				# Earth mass, kg
mu = D(3.986*10**14)						# Specific gravitational constant of Earth, G*mass_earth
circ_vel = np.sqrt(mu/(r_earth + 409000))	# Orbital velocity for circular orbit at 409km (ISS)

# init_vel = circ_vel		# Initial velocity at perigee, m/s
# --- OR: Given an apogee and perigee, calculate the velocity at the init altitude using the Vis-viva equation
init_vel = []
for alt in alt_pe:
	ap = r_earth + 409000
	pe = r_earth + alt
	a = (ap+pe)/2
	r = ap
	vel_visviva = np.sqrt(mu*((2/r) - (1/a)))
	init_vel.append(vel_visviva)

delta_vs_reqd = [circ_vel - i for i in init_vel]

C_D = 2.2				# Drag coefficient (2.2 is assumed constant by Jacchia 70)
t_step = [2000]			# Time step, s


masses = [5.00]				# A range of reasonable masses for a 3U cubesat
area_refs = [100, 300, 350]			# A range of cross-sectional areas for a 3U cubesat, cm^2

# fig1, ax1 = plt.subplots()
# fig1.suptitle("Satellite Altitude Decay for Varying Reference Areas")
# plt.xlabel("Time (days)")
# plt.ylabel("Altitude (km)")

dict = {}

for i, vel in enumerate(init_vel):
	for j, area in enumerate(area_refs):
		for k, mass in enumerate(masses):
			for m, ts in enumerate(t_step):
				dict_kword = str(mass) + "kg_" + str(area) + "cm2_" + str(alt_pe[i]/1000) + "km_" + str(ts) + "sec_Cd" + str(C_D)
				run_name = "deorbit_" + dict_kword
				
				A_ref = area/10000	# Covert cm^2 to m^2

				# alt_pe and velocity must be taken at PERIGEE!!!
				time, altitude, states = orbit(alt_pe[i], vel, C_D, A_ref, mass, ts, start_time, drag_on=True, circular_orbit=False, error_testing=False, num_orbits=100000)
				
				df = pd.DataFrame([time, altitude, states[0], states[1], states[2], states[3]]).transpose()
				df.columns = ['Time (s)', 'Altitude (m)', 'Radius (m)', 'Radial Velocity (m/s)', 'True Anomoly (rad)', 'Angular Velocity (rad/s)']
				# df.columns = ['time (s)', 'altitude (m)']
				df.to_csv(run_name+".csv", index=False)

				dict[dict_kword] = (time, altitude, states)
				# ax1.plot(dict.get(dict_kword)[0]/(2400*24), dict.get(dict_kword)[1]/1000)

# plt.legend([key for key in dict.keys()], loc='upper right')

# fig_no = 0
# colors = ['blue', 'green', 'red', 'orange']

# for i, mass in enumerate(masses):
# 	for j, ts in enumerate(t_step):
		
# 		fig_no += 1
# 		plt.figure(fig_no)
# 		kwords = []

# 		# for j in dict.keys():
# 		# 	if j.find(str(int(i))+"kg") == 0:
# 		# 		kwords.append(j)

# 		for k in dict.keys():
# 			plt.plot(dict.get(k)[0]/(2400*24), dict.get(k)[1]/1000)

# 		str_title = "Satellite Altitude Decay for Varying Reference Areas"
# 		plt.title(str_title)
# 		plt.xlabel("Time (days)")
# 		plt.ylabel("Altitude (km)")
# 		plt.legend([str(int(i)) + " cm$^2$" for i in area_refs], loc='upper right')

# plt.show()