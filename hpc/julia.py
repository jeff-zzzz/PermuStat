# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:18:45 2015
@author: jesong1126
julia.py
"""

z = -1.8-1.8j
print abs(z)

c = -0.62772-0.42193j
z = 0+0j
for n in range(9):
    z = z*z +c
    print "{}: z={:33}, abs(z)={:0.2f}, c={}".format(n, z, abs(z), c)

    
import time 
# area of complex space to investigate 
x1, x2, y1, y2 = -1.8, 1.8, -1.8, 1.8
c_real, c_imag = -0.62772, -.42193

def cal_pure_python(des_width, max_ier):
    x_step = (float(x2-x1)/float(des_width))
    y_step = (float(y2-y1)/float(des_width))
    x=[]
    y=[]
    ycoord = y2 
    while ycoord > y1:
        y.append(ycoord)
        ycoord += y_step 
    xcoord = x1
    while xcoord < x2:
        x.append(xcoord)
        xcoord += x_step 
        

zs = []
cs = []
for ycoord in y:
    for xcoord in x:
        zs.append(complex(xcoord,ycoord))
        cs.append(complex(c_real,c_imag))

print "Length of x: ", len(x)     
print "Total elements: ", len(zs)     
start_time = time.time()
output = calculate_z_serial_purepython(max_iter, zs,zc)
end_time = time.time()
secs = start_time - end_time
#for z in coordinates: 
#    for iteration in range(maxiter):
#        if abs(z) < 2.0 :
#            z = z*z +c
#        else: 
#            break
    

