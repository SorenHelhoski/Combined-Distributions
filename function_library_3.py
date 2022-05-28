import numpy as np

def Silverman(x):
    a= 4*np.std(x) / (3*len(x))
    return a**(1/5)


def transient_scaling_law(x,inverse=False): # x is in meters, returns with meters (diameters) (imp->crater by default)
    # in mks units
    velocity = 12_000.
    gravity = 1.62
    imp_angle = np.pi/2

    # in the same units as eachother
    imp_dens = 3.25
    targ_dens = 2.71
    
    if inverse:
        return  (x / ( 1.161 * (imp_dens/targ_dens)**(1/3) * (velocity)**0.44 * (gravity)**(-0.22) * (np.sin(imp_angle))**(1/3) ))**(1/0.78)
    
    return  1.161 * (imp_dens/targ_dens)**(1/3) * (x)**0.78 * (velocity)**0.44 * (gravity)**(-0.22) * (np.sin(imp_angle))**(1/3)

def scaling_law(a,inverse=False): #takes an asteroid diameter [km] and returns a basin diameter [km] (or the inverse)
    A = 3.484406473362318
    eta = -0.10879516054102951
    simple_complex = 20_000. #m
    
    if inverse:
        return 1/1000 * transient_scaling_law( (simple_complex**eta * a*1000 / A )**( 1/(1+eta)) , inverse=True)
        
    return 1/1000 * A * (simple_complex)**(-eta) * (transient_scaling_law(a*1000,inverse=False)) **(1+eta)

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
    size = line_values[2] # this will be useless when called for sample sites
    North= line_values[-2]
    East = line_values[-1]
    
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
            
    if line_values[0] == 'LB':   
        # returns the name, diameter, and coordinates of the basin
        return name[lower:upper],float(size),float(North),float(East)
    if line_values[0] == 'S':   
        # returns the name and coordinates of the sample
        return name[lower:upper],float(North),float(East)
