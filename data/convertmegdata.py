import scipy.io
import numpy as np
import pandas as pd


mat = scipy.io.loadmat('meg/3019/Resting_State/VEf_b_1_48.mat')

print(mat.keys())


signal = mat['VEf_b_1_48']

print("Shape:", signal.shape)
print("Type:", type(signal))




# Save to CSV
np.savetxt('meg_data_3.csv', signal, delimiter=',')
