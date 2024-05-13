# for analysis2 pi plot, binning genomic regions


import sys
import numpy as np


num_average = 5 # number of lines to average    


with open(sys.argv[1]) as infile:
    group = []
    counter = 0
    header = infile.readline().strip().split()
    header = header[:-1]
    print("\t".join(header))
    for line in infile:
        newline = line.strip().split()
        newline = newline[:-1]
        #print(newline) # get rid of seed number
        newline = list(map(float,newline))
        group.append(np.array(newline))

        if counter % num_average != 0:
            pass
        elif counter == 0:
            pass
        else:
            # average rows               
            avgs = np.mean(group, axis=0)
            print("\t".join(map(str,avgs)))
            group = [] # reset                  

        counter+=1
