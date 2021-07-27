
# coding: utf-8

# In[154]:


import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# In[18]:


def vel(h, r, e = None, theta = None, mu = None):
    
    if e == None:
        e = 0   
    if theta == None:
        theta = 0
    if mu == None:
        mu = 398600
        
    theta = theta*np.pi/180
    
    
    v_r = (mu/h)*e*np.sin(theta)
    
    v_t = (h/r)
    
    return v_r, v_t


# In[10]:


def ecc(r_a, r_p):
    
    e = (r_a - r_p)/(r_a + r_p)
    
    return e


# In[5]:


def ang_mom(r_a, r_p, mu = None):
    
    if mu == None:
        mu = 398600
        
    h = ((2*mu)**(1/2))*(((r_a*r_p)/(r_a + r_p))**(1/2))
    
    return h


# In[45]:


def delta_v(r_A, r_B = None, i = None, delta1 = None, mu = None):
    
    if r_B == None:
        r_B = r_A
    if i == None:
        i = 0
    if delta1 == None:
        delta1 = 0 
    if mu == None:
        mu = 398600
    
    i = i*np.pi/180
    delta1 = delta1*np.pi/180
    delta2 = i - delta1
    
    #orbit 1 (circle)
    h_1 = ang_mom(r_A, r_A) #angular momentum of orbit one
    
    v_rA1, v_tA1 = vel(h_1, r_A) #radial and tangential velocity at point A on orbit 1
    
    #orbit 2 (ellipse)
    e_t = ecc(r_A, r_B) #eccentricity of transfer orbit
    
    h_2 = ang_mom(r_A, r_B) #angular momentum of transfer orbit
    
    v_rA2, v_tA2 = vel(h_2, r_A, e_t) #necessary radial and tangential velocity at point A for transfer orbit 
    
    #orbit 3 (circle)
    h_3 = ang_mom(r_B, r_B) #angular momentum for final orbit
    
    v_rB1, v_tB1 = vel(h_2, r_B, e_t, 180) #radial and tangential velocities when craft reaches point B from the 
                                           #transfer orbit
        
    v_rB2, v_tB2 = vel(h_3, r_B) #necessary radial and tangential velocities at point B for GEO
    
    
    #Computing Delta V's            
    delta_vA = (((v_rA2 - v_rA1)**2) + v_tA1**2 + v_tA2**2 - 2*v_tA1*v_tA2*np.cos(delta1))**(1/2) 
    
    delta_vB = (((v_rB2 - v_rB1)**2) + v_tB1**2 + v_tB2**2 - 2*v_tB1*v_tB2*np.cos(delta2))**(1/2)
    
    delta_v_total = delta_vA + delta_vB
    
    #print("The total Delta V is: " + str(round(delta_v_total, 3)))
    
    return delta_v_total


# In[153]:


delta_v_total = delta_v(6678, 42164, 28.5, 0)
print(delta_v_total)


# In[53]:


delta_range = np.arange(28.5,0,-0.5)
print(delta_range)


# In[54]:


list_of_deltas = []

for i in delta_range:
    delta = delta_v(6678, 42164, 28.5, i)
    list_of_deltas.append(delta)
    
print(list_of_deltas)


# In[156]:


angles_and_deltas = pd.DataFrame({
    'Angle at Orbit 1': delta_range,
    'Delta V': list_of_deltas
})

print(angles_and_deltas)


# In[151]:


x = delta_range
y = list_of_deltas

plt.plot(x,y,color='blue')
plt.title('Delta V over Angle')
plt.xlabel('Angle Change from Orbit 1')
plt.ylabel("Total Delta V")
plt.plot(2, 4.23155406, "ro")
plt.show()


# In[150]:


del_v_angle = np.array(list(zip(delta_range,list_of_deltas)))

min_delta = []

for i in del_v_angle:
    if i[1] == del_v_angle[:,1].min():
        min_delta.append(i)
        print(min_delta_v)
    else:
        continue

print("The optimum angle to change in Orbit 1 is " + str(min_delta[0][0]) + " degrees, with the remainder being changed at the insertion point of Orbit2, which provides a minimum 'Delta V' value of " + str(round(min_delta[0][1], 4)) + ".")
print("Adjusting the entire angle from Orbit 1 yields a 'Delta V' value of " + str(round(del_v_angle[0,1], 4)) + ", which is " + str(round(100*((del_v_angle[0,1] - min_delta[0][1])/min_delta[0][1]), 2)) + "% higher than the optimum value.")
print("Adjusting the entire angle from Orbit 2 yields a 'Delta V' value of " + str(round(del_v_angle[len(del_v_angle)-1,1], 4)) + ", which is " + str(round(100*((del_v_angle[len(del_v_angle)-1,1] - min_delta[0][1])/min_delta[0][1]), 2)) + "% higher than the optimum value.")

