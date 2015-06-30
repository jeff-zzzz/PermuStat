"""Calculate statistics on simulated example.
This program uses the functions in :file:`get_T.py` to calculate
the t-statistics of 2 arrays. The results are saved in bin.
the files used will be downloaded if they are not in the working directory
of the IPython engines.
"""

from __future__ import print_function

from IPython.parallel import Client 
import numpy as np
from timeit import default_timer as clock
from permustat import perm_t

# simulated data 
A = np.arange(120).reshape((2,3,20))  
B = np.arange(180).reshape((2,3,30))  

nSim = 150

# Connect to the IPython cluster
c = Client()
c[:].run('./permustat/perm_t.py')

# the number of engines
n = len(c)
nSimn = nSim/n
id0 = c.ids[0]
v = c[:]
v.block = True
#fetch the pi-files
#print("downloading %i files of pi"%n)
#v.map(fetch_pi_file, files[:n])
v.map(perm_t.get_Ts, A, B, nSimn)
print("done")

# Run 10m digits on 1 engine
t1 = clock()
Ts1 = c[id0].apply_sync(perm_t.get_Ts, A, B, nSim)
t2 = clock()
digits_per_second1 = 10.0e6/(t2-t1)
print("Digits per second (1 core, 10m digits):   ", digits_per_second1)


# Run n*10m digits on all engines
t1 = clock()
Ts_each = v.map(perm_t.get_Ts, A, B, nSimn)
Ts_all = reduce_freqs(Ts_each)
t2 = clock()
digits_per_second8 = n*10.0e6/(t2-t1)
print("Digits per second (%i engines, %i0m digits): "%(n,n), digits_per_second8)

print("Speedup: ", digits_per_second8/digits_per_second1)

# plot_two_digit_freqs(freqs150m)
# plt.title("2 digit sequences in %i0m digits of pi"%n)
# plt.show()

