import numpy as np
import matplotlib.pyplot as plt

def heightcap(r1,r2,d): #intersects if abs(r1-r2) < d < r1+r2 #does not intersect if abs (r1-r2) >= d >= r1+r2 
   #calculates the height of the cap of the sphere
    if r1 <= r2:
        if abs(r1-r2) < d < (r1+r2):
           alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
           h1 = np.cos(alpha)*r1
           a = r1-h1
           return a
        elif (r1-r2) > d > -(r1+r2):
            alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
            h1 = -np.cos(alpha)*r1
            a = r1-h1
            return a
        
        elif (r1-r2) <= d <= abs(r1-r2):
            h1 = -r1
            a = r1-h1
            return a
        
        else:
            h1 = r1
            a = r1-h1
            return a
    elif r1 > r2:
        if (r1-r2) < d < (r1+r2):
           alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
           h1 = np.cos(alpha)*r1
           a = r1-h1
           return a
        elif -(r1+r2) < d < -abs(r1-r2):
            alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
            h1 = -np.cos(alpha)*r1
            a = r1-h1
            return a
        elif (r1-r2) >= d >= -abs(r1-r2):
            h1 = r1
            a = r1-h1
            return a
        else:     
            h1 = r1
            a = r1-h1
            return a

def areacap(r1,h):
    f = 2*np.pi*r1*h
    return abs(f)

def diffarea(r1,h):
    g = 4*np.pi*r1**2-areacap(r1,h)
    return g


def area(r1,r2,d):
    a = 4*np.pi*(r1**2+r2**2)
    b = 2*np.pi*r1*heightcap(r1,r2,d)
    c = 2*np.pi*r2*heightcap(r2,r1,d)
    f = a-b-c
    return f
    
def onefunction(r1,r2,d): #intersects if abs(r1-r2) < d < r1+r2 #does not intersect if abs (r1-r2) >= d >= r1+r2 
   #calculates the area of two intersecting spheres
    if r1 <= r2:
        if abs(r1-r2) < d < (r1+r2):
           alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
           beta = np.arccos((r2**2+d**2-r1**2)/(2*r2*d))
           h1 = np.cos(alpha)*r1
           h2 = np.cos(beta)*r2
           a = r1-h1
           b = r2-h2
           
        elif (r1-r2) > d > -(r1+r2):
            alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
            beta = np.arccos((r2**2+d**2-r1**2)/(2*r2*d))
            h1 = -np.cos(alpha)*r1
            h2 = -np.cos(beta)*r2
            a = r1-h1
            b = r2-h2
            
        
        elif (r1-r2) <= d <= abs(r1-r2):
            h1 = -r1
            h2 = r2
            a = r1-h1
            b = r2-h2
            
        
        else:
            h1 = r1
            h2 = r2
            a = r1-h1
            b = r2-h2
            
    elif r1 > r2:
        if (r1-r2) < d < (r1+r2):
           alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
           beta = np.arccos((r2**2+d**2-r1**2)/(2*r2*d))
           h1 = np.cos(alpha)*r1
           h2 = np.cos(beta)*r2
           a = r1-h1
           b = r2-h2
           
        elif -(r1+r2) < d < -abs(r1-r2):
            alpha = np.arccos((r1**2+d**2-r2**2)/(2*r1*d))
            beta = np.arccos((r2**2+d**2-r1**2)/(2*r2*d))
            h1 = -np.cos(alpha)*r1
            h2 = -np.cos(beta)*r2
            a = r1-h1
            b = r2-h2
            
        elif (r1-r2) >= d >= -abs(r1-r2):
            h1 = r1
            h2 = -r2
            a = r1-h1
            b = r2-h2
            
        else:     
            h1 = r1
            h2 = r2
            a = r1-h1
            b = r2-h2
        
    f = 4*np.pi*(r1**2+r2**2)
    g = 2*np.pi*r1*a
    h = 2*np.pi*r2*b
    return f-g-h
