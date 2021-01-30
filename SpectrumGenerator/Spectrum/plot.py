import numpy as np
import matplotlib.pyplot as plt
import os
with open('/home/rustin/Projets/JSC_PI/spectrum_generator/data/spectrum1.txt','r') as f:
    data_x,data_y = [],[]
    lines = f.readlines()
    count = 0
    for line in lines:
        count = count + 1
        if(count > 2):
        	atom = line.split(" ")
        	print(atom)
        	data_x.append(float(atom[1]))
        	data_y.append(float(atom[2]))
		#print(count)
plt.plot(data_x,data_y,'o',markersize=4)

plt.xlabel('x') 
plt.ylabel('y')

plt.savefig(os.path.join(os.path.expanduser('~'), 'Projets/JSC_PI/spectrum_generator/data/plots/', 'spectrum_ploted.png')) 
plt.show()
