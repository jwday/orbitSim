import numpy as np
from decimal import Decimal as D
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob

cwd = os.getcwd()
# data_files = [f for f in glob.glob('deorbit*') if os.path.isfile(f)]
data_files = [f for f in glob.glob('deorbit_5.0*') if os.path.isfile(f)]

fig1, ax1 = plt.subplots()
fig1.suptitle("Satellite Altitude Decay for Varying Reference Areas")
plt.xlabel("Time (days)")
plt.ylabel("Altitude (km)")

for f in data_files:
	df = pd.read_csv(f)
	ax1.plot(df['Time (s)']/(2400*24), df['Altitude (m)']/1000)
	
plt.legend([f for f in data_files], loc='upper right')
plt.show()