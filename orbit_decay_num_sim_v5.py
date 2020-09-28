# Numerical Integration of Polar Eqns. of Motion for Orbital Motion

# =============================================================================
# Polar Equations of Motion
# =============================================================================
#	 r'' = r(th')^2 - mu/r^2			   	# Radial acceleration
# 	 th'' = -2(r')(th')/r + a_T/r			# Angular acceleration

# As a baseline, first-order method solely to test functionality, drag will be 
# constant on the object of interest to observe how the orbit decays
# Define variables to integrate:
#	Let:
# 		x[0] = r
#		x[1] = r'
#		x[2] = th
#		x[3] = th'
#
#	Therefore:
#		x[0]' = x[1]
#		x[1]' = x[0]*x[3]**2 - mu/(x[0]**2)
#		x[2]' = x[3]
#		x[3]' = -2*x[1]*x[3]/x[0] + a_T/x[0]


# =============================================================================
# Script begin
# =============================================================================
from __future__ import print_function
import numpy as np
import math
from scipy.integrate import ode
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
import decimal
from decimal import Decimal as D
from decimal_pi import pi as dec_pi
from atm_model_exp import nasa_eam_model, msise00_model
from mpmath import *
import pandas as pd
from datetime import datetime, timedelta

# =============================================================================
# Specify constants and simulation parameters
# =============================================================================

# -- Set Integration Parameters -----------------------------------------------

# 'dopri5' is 4th-order Dormand-Prince method with error calculated to 5th order,
# identical to the MATLAB solver 'ode45'. 'dopri5' has multiple options for step
# size adjustment, tolerance, max step size, etc...
# Defaults:
# 	rtol=1e-6, atol=1e-12,
#   nsteps=500,
#	safety=0.9,
#	ifactor=10.0

# 'dop853' is an 8th-order Dormand-Prince method with 5th order error and something
# something 3rd order. It's about 40% slower than dopri45 but gives better results.
# Options are identical to dopri45.
# Defaults:
# 	rtol=1e-6, atol=1e-12,
#   nsteps=500,
#	safety=0.9,
#	ifactor=6.0


def orbit(h_perigee, v_perigee, C_D, A_ref, mass_sat, t_step, start_time, **kwargs):
	drag_on = kwargs["drag_on"]
	circular_orbit = kwargs["circular_orbit"]
	error_testing = kwargs["error_testing"]
	sim_time = start_time
	
	if error_testing:
		circular_orbit=False
		drag_on=False
		num_orbits = kwargs["num_orbits"]
	else:
		pass


	print("")
	print("")
	if not error_testing:
		print("Deorbiting a " + str(int(mass_sat)) + "kg cubesat of ref. area " + str(int(A_ref*10000)) + " cm^2 (" + str(int(t_step)) + " sec time step)")
	else:
		print("Error testing over " + str(num_orbits) + " orbits")


	# =============================================================================
	# Define Numerical and Integration Parameters
	# =============================================================================	

	decimal.getcontext().prec=8

	integrator = 'dopri5'
	rtol=10**-3
	atol=10**-6
	nsteps=1000

	dt = t_step  # Time step (seconds)




	# =============================================================================
	# Define Constants
	# =============================================================================

	# ISS_alt = D(408)*1000  # Kilometers

	r_earth = D(6378100)  	# Radius of Earth, m
	G = D(6.67*10**-11)  							# Gravitational constant
	mass_earth = D(5.972*10**24)  					# Earth mass, kg
	# mu = G*mass_earth  								# Specific gravitational constant of Earth
	mu = D(3.986*10**14)



	# =============================================================================
	# Calculate remaining orbit characteristics
	# =============================================================================

	r_perigee = r_earth + D(h_perigee)						# Radius at perigee
	
	if circular_orbit:
		a = r_perigee										# Semi-major axis (circular)
		H = (a*mu)**D(0.5) 									# Angular momentum (circular)
		v_perigee = H/a  									# Orbital velocity (circular)
	else:
		a = 1/(2/r_perigee - (D(v_perigee)**D(2.0))/mu)  	# Semi-major axis	
	

	## DO I NEED THIS STUFF??
		# eccen = (v_perigee**2)/2 - mu/r_perigee  		# Eccentricity
		# orbit_h = r_perigee*v_perigee					# Specific angular momentum



		# # I'm defining the apoapsis and periapsis of the ISS orbit around Earth and
		# # using these two properties to determine the remaining orbit characteristics
		# # init_r_a = earth_rad + init_alt  # Apoapsis
		# # init_r_p = init_r_a  # Periapsis (for circular orbit, r_p = r_a)

		# init_a = (init_r_a + init_r_p)/2  # Semi-major axis
		# init_e = (init_r_a - init_r_p)/(init_r_a + init_r_p)  # Eccentricity
		# init_p = init_a*(1 - init_e**2)  # Semi-latus rectum

		# init_orbital_E = -mu/(2*init_a) # Energy
		# init_orbital_H = (init_p*mu)**D(0.5) # Angular momentum
		# init_orbital_vel = init_orbital_H/init_a  # Orbital velocity
		# init_orbital_period = 2*dec_pi()*((init_a**3)/mu)**D(0.5)  # Orbital period
	##



	# =============================================================================
	# Initial State Vector
	# =============================================================================

	x0 = np.array([float(r_perigee), 0, 0, float(v_perigee)/float(r_perigee)])
	t0 = 0



	# =============================================================================
	# State Function
	# =============================================================================

	def func(t, x, drag):
	    # Returns a 1x4 array of x_dots, aka the derivates of the EOMs above
	    x_dot = np.zeros(4)  # Pre-make the array

	    x_dot[0] = D(x[1])  																	# r'
	    x_dot[1] = D(x[0])*D(x[3])**2 - mu/D(x[0])**2  											# r''
	    x_dot[2] = D(x[3])  																	# th'
	    x_dot[3] = -2*D(x[1])*D(x[3])/D(x[0]) + drag(x[0], x[3], mass_sat, sim_time)/D(x[0])  	# th''
	  
	    x_dots = np.array([x_dot[0], x_dot[1], x_dot[2], x_dot[3]])
	    
	    return x_dots




	# =============================================================================
	# Drag Function
	# =============================================================================

	def drag(alt, ang_vel, mass_sat, sim_time):
		if drag_on and not error_testing:
			vel = D(alt)*(D(ang_vel) - D(7.2921159*10**-5))  	# Takes into account the angular velocity of the Earth
			# rho = nasa_eam_model(D(alt) - r_earth)[2]					# Density is calculated by passing altitude to this function
			rho = msise00_model((D(alt) - r_earth), sim_time)[2]		
			drag_force = D(0.5)*D(C_D)*D(A_ref)*rho*vel**2		# Drag force
			a_T = -drag_force/D(mass_sat)						# Tengential acceleration from F = ma

		else:
			a_T = 0

		return D(a_T)




	# =============================================================================
	# Integration Routine
	# =============================================================================
	r = ode(func)  											# Make the ODE object
	r.set_integrator(integrator, rtol=rtol, atol=atol, nsteps=nsteps)		# Set integrator and tolerance
	r.set_initial_value(x0, t0)  							# Set initial values from Initial State Vector
	r.set_f_params(drag)  									# Pass the 'drag' function as a parameter into the EOMs



	x = [[x0[0]], [x0[1]], [x0[2]], [x0[3]]] 				# Initialize state vector list for logging
	t = [t0]  												# Initialize time vector list for logging



	# This block sets up the "integration time remaining" estimator
	# -----------------------------------------------------------------------------

	comp_time_zero = time.time()  # Store the time at the beginning of integration
	comp_time = []
	comp_time.append(0)

	# -----------------------------------------------------------------------------




	# =============================================================================
	# Integrate until we crash (when altitude goes to zero)
	# =============================================================================

	def append_data(x, r, t):
		# Append to state vector list for logging
		x[0].append(float(r.y[0]))  	# r
		x[1].append(float(r.y[1]))  	# r'
		x[2].append(float(r.y[2]))  	# theta (aka true anomoly)
		x[3].append(float(r.y[3]))  	# theta'

		# Append to time vector list for logging
		t.append(r.t)  					# time

		return t, x


	N = 0		# Number of full orbits
	N_prev = 0
	orbit = [0]

	# These are only used for error testing
	radius = []
	true_anomoly = []
	estimated_periapsis = []


	while r.successful():
		r.integrate(D(r.t)+D(dt))
		th_home = 2*np.pi*N
		sim_time += timedelta(seconds=dt)

		# Log the number of full orbits made
		if (float(r.y[2]) - th_home) >= 2*np.pi:  
			N += 1
			orbit.append(N)

		else:
			pass

		# Do something different depending on the kwarg parameters passed to the function
		if drag_on and not error_testing and r.y[0] > r_earth:
			t, x = append_data(x, r, t)
			min_altitude = min(x[0]) - float(r_earth)
			print("\rTime: {:2.1f}      Min Altitude: {:2.1f}".format(r.t, min_altitude), end='')

		elif not drag_on and not error_testing and N < num_orbits:
			t, x = append_data(x, r, t)
			print("\rTime: {:2.1f}      Number of Orbits: {:d}".format(r.t, N), end='')


		elif error_testing and N <= num_orbits:
			t, x = append_data(x, r, t)
			print("\rTime: {:2.1f}      Number of Orbits: {:d}".format(r.t, N), end='')

			if N != N_prev:
				radius.append([float(x[0][-1]), float(x[0][-2]), float(x[0][-3]), float(x[0][-4])])
				true_anomoly.append([float(x[2][-1])-2*np.pi*N, float(x[2][-2])-2*np.pi*N, float(x[2][-3])-2*np.pi*N, float(x[2][-4])-2*np.pi*N])
				interpolated_orbit = interp1d(true_anomoly[-1], radius[-1], kind='cubic')

				estimated_periapsis.append(interpolated_orbit(0).tolist())

				N_prev = N
			else:
				pass
		else:
			break



	times = np.array(t)			# Converts time list to array
	altitude = np.array([i-float(r_earth) for i in x[0]])
	states = np.array(x)  		# Convert states list to array

	if not error_testing:
		return(times, altitude, states)
	else:
		plot_error(orbit, estimated_periapsis, num_orbits, r_perigee, t_step, integrator, atol, rtol)



def plot_error(orbit, estimated_periapsis, num_orbits, r_perigee, t_step, integrator, atol, rtol):
	error = [0] + [float(r_perigee) - i for i in estimated_periapsis]

	fig_err = plt.figure(figsize=(8,8))
	plt.plot(orbit[:-1], error)

	main_title = 'Variation (m) of Periapsis over ' + str(num_orbits) + ' Orbits'
	subtitle1 = 'Step Size: ' + str(t_step) + ' sec.     Precision: ' + str(10**-decimal.getcontext().prec)
	subtitle2 = 'Integrator: ' + integrator + '     atol: ' + str(atol) + '     rtol: ' + str(rtol)

	plt.title(main_title + '\n' + subtitle1 + '\n' + subtitle2)
	plt.xlabel('Orbit No.')
	plt.ylabel('Variation from Inital Value (m)')
	plt.show()
