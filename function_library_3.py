import numpy as np
from random import randrange

def Silverman(x):
    a= 4/ (3*len(x))
    return a**(1/5)*np.std(x)

#===================================================================================

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
    #A = 1.35
    #eta = 0.086
    simple_complex = 19_000. #m
    
    if inverse:
        return 1/1000 * transient_scaling_law( (simple_complex**eta * a*1000 / A )**( 1/(1+eta)) , inverse=True)
        
    return 1/1000 * A * (simple_complex)**(-eta) * (transient_scaling_law(a*1000,inverse=False)) **(1+eta)

def asteroid_diameter(filename): # takes a filename and outputs the asteroid diameter 
    file0    = open(filename,'r')
    time_0, time_f, save_step, R_p, M_p, R_imp, v_imp = eval(file0.readline().replace('\n',''))
    file0.close()
    return 2*R_imp

#===================================================================================

import pandas as pd 

def extract(file0): # quality of life function that is used to read the ejecta.data file
    header, data = {},{}
    header['time_0'], header['time_f'], header['save_step'], header['R_p'], header['M_p'], header['R_imp'], header['v_imp'] = eval(file0.readline().replace('\n',''))
    
    data['used'] = np.array(list(eval(file0.readline().replace('\n',''))))
    data['material'] = np.array(list(eval(file0.readline().replace('\n',''))))
    data['volume']   = np.array(list(eval(file0.readline().replace('\n',''))))
    data['spread']   = np.array(list(eval(file0.readline().replace('\n',''))))
    data['mass']     = np.array(list(eval(file0.readline().replace('\n',''))))
    data['launch_time']     = np.array(list(eval(file0.readline().replace('\n',''))))
    data['launch_polar']    = np.array(list(eval(file0.readline().replace('\n',''))))
    data['launch_length']   = np.array(list(eval(file0.readline().replace('\n',''))))
    data['launch_height']   = np.array(list(eval(file0.readline().replace('\n',''))))
    data['launch_velocity'] = np.array(list(eval(file0.readline().replace('\n',''))))
    data['launch_angle']    = np.array(list(eval(file0.readline().replace('\n',''))))
    data['pressure']    = np.array(list(eval(file0.readline().replace('\n',''))))
    data['temperature'] = np.array(list(eval(file0.readline().replace('\n',''))))
    data['landing_time']   = np.array(list(eval(file0.readline().replace('\n',''))))
    data['landing_polar']  = np.array(list(eval(file0.readline().replace('\n',''))))
    data['landing_length'] = np.array(list(eval(file0.readline().replace('\n',''))))
    data['max_altitude']   = np.array(list(eval(file0.readline().replace('\n',''))))
    
    return header, pd.DataFrame.from_dict(data)

def remove_outer_spaces(string):
    # remove the spaces on either side of the name without destroying spaces in the center of the name
    # i.e. South Pole-Aitkens
    lower,upper,length = 0,-1,len(string)
    for i in range(length):
        if string[i]==' ':
            continue
        else:
            lower = i
            break
    for i in range(length):
        if string[-i-1]==' ':
            continue
        else:
            upper = -i
            break
    return string[lower:upper]

def trim(line_values): # reads a line of Basin.txt
    
    # read the ':' separated list
    name = line_values[1]
    North= line_values[2]
    East = line_values[3]
    if line_values[0] == 'B':
        size = line_values[4] 
        order= line_values[5] 
        color= line_values[6]
            
    if line_values[0] == 'B':   
        # returns the name, diameter, and coordinates of the basin
        return remove_outer_spaces(name), float(North), float(East), float(size), -1*float(order), color.replace(' ','')
    if line_values[0] == 'S':   
        # returns the name and coordinates of the sample
        return remove_outer_spaces(name), float(North),float(East)
    
#===================================================================================

def haversine(lat1, lon1, lat2, lon2, Radius): # haversine formula
    dLat = np.radians(lat2 - lat1) # change in latitude
    dLon = np.radians(lon2 - lon1) # change in longitude
    lat1 = np.radians(lat1) # latitude 1
    lat2 = np.radians(lat2) # latitude 2
    a = np.sin(dLat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dLon/2)**2 # intermediate step
    c = 2*np.arcsin(np.sqrt(a)) # interior angle between two coordinates [rad]
    return Radius * c # returns distance between two coordinates

#===================================================================================

def solve(R_planet,R_basin,lam1,lam2):
    numer = 2*(np.sin(R_basin/(2*R_planet)))**2 -1 + np.cos(lam1-lam2) - np.cos(lam1)*np.cos(lam2)
    denom = -np.cos(lam1)*np.cos(lam2)
    return np.arccos(numer/denom)

def process(x_prime,y_prime,phi1):
    xs = [] ; ys = []
    for each in range(2):
        for i in range(len(x_prime)):
            if not each : x = (phi1+x_prime[i]+np.pi)%(2*np.pi) - np.pi
            if     each : x = (phi1-x_prime[-i-1]+np.pi)%(2*np.pi) - np.pi
            if x < -0.99*np.pi or x > 0.99*np.pi:
                xs.append(None) ; ys.append(None)
            xs.append(np.degrees(x))  
            if not each : ys.append(np.degrees(y_prime[i]))
            if     each : ys.append(np.degrees(y_prime[-i-1]))
    xs.append(xs[0]) ; ys.append(ys[0])
    
    return xs,ys

def transform(x,lam_min,lam_max):
    k   = 5
    L   = lam_max - lam_min
    return L * ( 1 + np.exp(-k*x) )**-1 + lam_min
