import numpy as np

def scaling_law(a): #takes an asteroid diameter [km] and returns a basin diameter [km]
    return (359.55178234952507*a)**(0.6950181547012931)

def asteroid_diameter(filename): # takes a filename and outputs the asteroid diameter 
    file0    = open(filename,'r')
    time_0, time_f, save_step, R_p, M_p, R_imp, v_imp, coords, R_c = eval(file0.readline().replace('\n',''))
    file0.close()
    return 2*R_imp

def extract(file0): # quality of life function that is used to read the ejecta.data file
    return np.array(list(eval(file0.readline().replace('\n',''))))

def trim(line_values): # reads a line of Basin.txt
    
    # read the ':' separated list
    name = line_values[1]
    size = line_values[2]
    North= line_values[3]
    East = line_values[4]
    
    # remove the spaces on either side of the name without destroying spaces in the center of the name
    # i.e. South Pole-Aitkens
    lower,upper,length = 0,-1,len(name)
    for i in range(length):
        if name[i]==' ':
            continue
        else:
            lower = i
            break
    for i in range(length):
        if name[-i-1]==' ':
            continue
        else:
            upper = -i
            break
            
    # returns the name, diameter, and coordinates of the basin
    return name[lower:upper],float(size),float(North),float(East)
