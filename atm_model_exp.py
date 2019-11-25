# A simple, exponential model of upper atmosphere based on the Earth Atmosphere Model (NASA) based
# largely on the ideal gas model 
# Inputs: Altitude
# Outputs: Temperature, Pressure, Density

import sys
import math
from decimal import Decimal as D
from decimal_pi import exp as dec_exp


def nasa_eam(alt):
	if (alt < 11000):
		T = D(15.04) - D(0.00649)*alt  # Temperature
		p = D(101.29)*((T+D(273.15))/D(288.08))**D(5.256)  # Pressure

	elif (alt >= 11000 and alt < 25000):
		T = -56.46
		p = D(22.65)*dec_exp(D(1.73)-D(0.000157)*alt)

	elif (alt >= 25000):
		T = D(-131.21) + D(0.00299)*alt
		p = D(2.488)*((T+D(273.15))/D(216.6))**D(-11.388)

	r = p/(D(0.2869)*(D(T)+D(273.15)))  # Density

	return [T, p, r]


if __name__ == '__main__':
	alt = float(sys.argv[1])
	print('')
	print('You want to know the atmospheric properties at %f meters?' % alt)
	print('')
	[T, p, r] = nasa_eam(alt)
	print('Temperature: %f C' % T)
	print('Pressure: %f kPa?' % p)
	print('Density: %f kg/m^3' % r)