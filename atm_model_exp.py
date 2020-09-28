# A simple, exponential model of upper atmosphere based on the Earth Atmosphere Model (NASA) based
# largely on the ideal gas model 
# Inputs: Altitude
# Outputs: Temperature, Pressure, Density

import sys
import math
from decimal import Decimal as D
from decimal_pi import exp as dec_exp
import msise00
from datetime import datetime


def nasa_eam_model(alt):
	if (alt < 11000): # meters
		T = D(15.04) - D(0.00649)*alt  						# Temperature (C)
		p = D(101.29)*((T+D(273.15))/D(288.08))**D(5.256) 	# Pressure (kPa)

	elif (alt >= 11000 and alt < 25000):
		T = -56.46											# Temperature (C)
		p = D(22.65)*dec_exp(D(1.73)-D(0.000157)*alt)		# Pressure (kPa)

	elif (alt >= 25000):
		T = D(-131.21) + D(0.00299)*alt						# Temperature (C)
		p = D(2.488)*((T+D(273.15))/D(216.6))**D(-11.388)	# Pressure (kPa)

	r = p/(D(0.2869)*(D(T)+D(273.15))) 						# Density (kg/m^3)

	return [T, p, r]


def msise00_model(alt, datetime_obj):
	atmos = msise00.run(time=datetime_obj, altkm=alt/1000, glat=65., glon=-148.)
	r = D(atmos['Total'].values[0,0,0,0])			# Denisty (kg/m^3)
	T = D(atmos['Tn'].values[0,0,0,0])				# Temperature (K)
	p = r*D(0.2869)*T								# Pressure (kPa)

	return [T-D(273.15), p, r]						# Return C, kPa, kg/m^3
	

if __name__ == '__main__':
	alt = float(sys.argv[1])
	print('')
	print('You want to know the atmospheric properties at %f meters?' % alt)
	print('')
	[T, p, r] = nasa_eam_model(alt)
	print('Temperature: %f C' % T)
	print('Pressure: %f kPa' % p)
	print('Density: %f kg/m^3' % r)