# Clohessy-Wiltshire
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm

z = 54.5										# Distance from center of station
mass = 5										# Mass of spacecraft, kg

r_earth = (np.average([6378137, 6356752]))  	# Radius of Earth, m
G = (6.67*10**-11)  							# Gravitational constant
mass_earth = (5.972*10**24)  					# Earth mass, kg
mu = G*mass_earth  								# Specific gravitational constant of Earth
mu = 3.986*10**14

at = 408000 + r_earth							# Semi-major axis of orbit (circular)
at = 6.781*10**6
n = math.sqrt(mu/at**3)

z_pp = (n**2)*z
thrust = z_pp*mass

print()
# print('Standard Gravitational Constant (\u03BC): {:.3e}'.format(mu))
print('Station-Keeping Distance: {:>5} meters'.format(z))
print('Station-Keeping Thrust: {:>6} \u03BCN'.format(round(thrust*10**6)))
print()


# ax2 = ax1.twinx()

# z_target = [15, 30, 45, 60]  # Target point distance from target spacecraft CoM
z_target = np.arange(0, 61, 0.1)		# Target position
d = np.arange(0, 5, 0.01)				# Stray distances
Z, D = np.meshgrid(z_target, d)
colors = ['r', 'y', 'g', 'b']

columns = ['{0} m'.format(i) for i in z_target]

dv_req = pd.DataFrame(index=d, columns=columns)
tof_between_burns = pd.DataFrame(index=d, columns=columns)


for i in z_target:
	# z0 = i + d								# Initial distance = observation target distance (z_target) + outward maximum allowable stray distance (d)
	t_to_min = np.arccos((i-d)/(i+d))/n		# Time to minimum distance (inward maximum stray distance from target, z_target - d)
	v_at_min = -n*np.sin(n*t_to_min)*(i+d)	# Velocity at minimum distance (velocity when at inward maximum stray distance from target)
	dv = -2*v_at_min						# delta-V required to exactly reverse velocity (hence factor of -2)
	# impulse = mass * dv						# Impulse required to exactly reverse velocity

	target_dist_str = '{0} m'.format(i)
	
	dv_req[target_dist_str] = dv
	tof_between_burns[target_dist_str] = 2*t_to_min

dv_tot = dv_req/(tof_between_burns/5520)


fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=True, dpi=300, figsize=(7.5, 5.5))
cs = ax1.contourf(Z, D, dv_req, levels=11, cmap=cm.Blues)
ax1.set_title('\u0394V Required to Stay Within Allowable Stray Distance (m/s)')
cbar = fig1.colorbar(cs, ax=ax1)

cs = ax2.contourf(Z, D, tof_between_burns, levels=11, cmap=cm.Oranges)
ax2.set_title('Time Between Required Burns (sec)')
cbar = fig1.colorbar(cs, ax=ax2)

cs = ax3.contourf(Z, D, dv_tot, levels=11, cmap=cm.Greens)
ax3.set_title('\u0394V Required per Orbit (m/s)')
cbar = fig1.colorbar(cs, ax=ax3)

# Create a big, empty subplot and label the axes on it instead of the individual subplots
# fig1.suptitle('Target Out-of-Plane Distance and Allowable Stray Distance')
fig1.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel('Inspection Target Out-of-Plane Distance (m)', x=0.4)
plt.ylabel('Allowable Stray Distance from Target (m)')

fig1.tight_layout(h_pad=1.0)
# plt.subplots_adjust(hspace=0.4)

plt.show()