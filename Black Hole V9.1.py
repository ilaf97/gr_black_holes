# ======= DEFINE LIBRARIES ======= #

import matplotlib.pyplot as plt      #Graph plotting tools      
import math                          #Mathematical functions
import time as time_                 #Time measurement functions
import numpy as np                   #Basic fucntions toolkit
import sys                           #System commands

# ================================ #


# ======= DEFINE INPUT VARIABLES/CONSTANTS ======= #1.5

c = 299792458               #Speed of light
pi = math.pi                #pi
theta = 1/(2-2**(1./3.))    #Value always used to compute Forrest-Ruth algorithm (constant)
r_s = 1                     #Schwarschild radius
phi_i = 0.0                 #Starting angle
t_i = 0.00                  #Starting time

# ================================================== #



# ======= USER INPUTS TO CHOOSE l VALUE ======= #
 
l_factor = input("l = x√(12)\nEnter value for x (must be 1 < x ≤ 4): ")     #User prompt and input line printed in console
type(l_factor); l_factor = float(l_factor)                                  #Convert from string type to float type    
l = l_factor*math.sqrt(12)                                                  #Use user input to calculate l

if l_factor < 1 or l_factor > 4:                                #Condition to check that user input is valid
    input("Invalid entry - enter value between 1 ≤ x ≤ 4: ")    #User prompt to input different value when user input invalid
    type(l_factor); l_factor = float(l_factor)                  #Convert from string type to float type 
    l = l_factor*math.sqrt(12)                                  #Use user input to calculate l
    
    if l_factor < 1 or l_factor > 4:                            #Condition to check that user input is valid
        print("Program teriminated - please retry")             #Message printed to console if user input invalid
        sys.exit()                                              #Terminate program 


# ======= PLOT EFFECTIVE POTENTIAL AS FUNCTION OF R ======= #
    
    
# ------- Calculate V and create lists for V and distance x ------ #

x = r_s-0.005                                                           #Define initial distance
x_values =[x]                                                           #Create distance values list 
V_values =[-(r_s/(2*x))+((l**2)/(2*(x**2)))-((r_s*(l**2))/(2*(x**3)))]  #Create V values list 

while x < 100*l_factor*r_s:
    V = -(r_s/(2*x))+((l**2)/(2*(x**2)))-((r_s*(l**2))/(2*(x**3)))      #Calculate V value for corresponding x value
    x = x + (0.1*r_s)                                                   #move to next x value 
    x_values.append(x); V_values.append(V)                              #Append new values to lists

# -- End -- #    
#
#
#
# ------- Find turning point - based on code by Malvolio (2013) ------- #

x_at_tp = []; V_at_tp = []                                                  #Create arrays to hold values at turning points
n = 2                                                                       #Loop counter
while n <= (len(x_values))-2:                                               #Loop whilst loop counter less than array length - 2
    if ((V_values[n-2] < V_values[n-1] and V_values[n] < V_values[n-1])     #Condition satisfied at maximum
    or (V_values[n-2] > V_values[n-1] and V_values[n] > V_values[n-1] )):   #Conditiion satisfied at minimum
        x_at_tp.append(x_values[n-1])                                       #Append time at turning points
        V_at_tp.append(V_values[n-1])                                       #Append radius at turning points
    n = n + 1                                                               #Increase loop counter by one
    
# -- End -- #
#
#
#
# ------- Find minimum and maximum values of V ------- #

V_min = min(float(d) for d in V_at_tp)                                  #Find minimum values in V_at_tp
V_max = max(float(d) for d in V_at_tp)                                  #Find maximum values in V_at_tp 
x_min = max(float(d) for d in x_at_tp)                                  #Find minimum values in x_at_tp 
print("Min V =", V_min, ",  Max V =", V_max, ", x at min V =", x_min)   #Print values found

plt.xlabel('Distance [r_s]'); plt.ylabel('V(r)'); plt.title('Effectiv3e Potential')  #Define plot lables
plt.xlim(0, (x_min+0.1*x_min)); plt.ylim((V_min+(0.1*V_min)), (V_max+(0.1*V_max)))  #Scale plot axes to show interesting region
plt.plot(x_values, V_values)                                                        #Define things to plot
plt.show()                                                                          #Show plots


# ============================================================ #



# ======= FIND POSITIONS OF ORBITS ======= #

V_half_max = V_max/2        #Find the half point between zero and V_max
V_quart_min = V_min/4       #Find the half point between V_min and zero

orbit = [0, V_min, V_max, V_quart_min, V_half_max, V_half_max]  #Array of V(r) values corrsponding to different orbits
#In order: NULL, Stable circular, unstable circular, precession, escape, trapped

i = 1                       #Loop counter
V_orbit =[]; r_orbit = []   #Define empty lists to hold values
while i <= 5:               #Condition to cycle through each entry in orbit list
    
    if i == 5:                                                      #Special case for trapped orbits
        j = max(V_values, key=lambda  X: abs(X - orbit[i]))         #Based on code by kennytm (2012)       
        #Find the value to the left of the V(r) peak closest to V_half_max  
        V_orbit.append(j)                                           #Append the V(r) value found in line above      
        r_app = round((x_values[V_values.index(V_orbit[i-1])]), 2)  
        #Find the r value corresponding to the V(r) value found two lines above by finding the index value of the V(r) value in V_values list and looking for corresponding x_values list index value
        
        if r_app < r_s:                         #Condition to check if r value found above is less than Schwarzschild radius
            r_app = r_app + 10*(r_s - r_app)    #Increase r value to ensure it is larger than Schwarzschild radius
            r_orbit.append(r_app)               #Append r value to r_orbit list

        else: r_orbit.append(r_app)             #Append r value to r_orbit list


    else:                                                           #Case for other orbits
        j = min(V_values, key=lambda  X: abs(X - orbit[i]))         #Based on code by kennytm (2012)  
        #Find value to the right of the v(r) peak closest to V value listed in orbit list
        V_orbit.append(j)                                           #Append the V(r) value found in the line above
        r_app = round((x_values[V_values.index(V_orbit[i-1])]), 2)
        #Find the r value corresponding to the V(r) value found two lines above by finding the index value of the V(r) value in V_values list and looking for corresponding x_values list index value

        
        if r_app < r_s:                         #Condition to check if r value found above is less than Schwarzschild radius     
            r_app = r_app + 10*(r_s - r_app)    #Increase r value to ensure it is larger than Schwarzschild radius
            r_orbit.append(r_app)               #Append r value to r_orbit list

        else: r_orbit.append(r_app)             #Append r value to r_orbit list
    
    i = i + 1                                   #Increase loop counter value


# ======================================== #
    
    
    
# ======= USER INPUTS TO CHOOSE E and r ======= #
    
print("--- ORBIT PROFILE EXAMPLE VALUES ---")   #Print new line message to user
print("Stable Circular = ", r_orbit[0], "\nPrecession = ", r_orbit[2], "\nEscape = ", r_orbit[3])   
print("Unstable Circular = ", r_orbit[1], "\nTrapped = ", r_orbit[4])
#Print suggested values for orbits to user 

r_i = input("Enter value for distance: " )  #User prompt to input value for r
r_i = float(r_i)                            #Convert user intput from string type to float type
    
if r_i < r_s:                               #Condition to check if user input value for r is less than Schwarszchild radius
    r_i = input("Value must be larger than 1. PLease enter new value: ")    #Error message and input prompt to user
    r_i = float(r_i)                                                        #Convert user input from string type to float type

    if r_i < r_s:                           #Condition to check if user input value for r is less than Schwarszchild radius
        print("Program terminated - please try again")  #Error message to ell user invalid entry
        sys.exit()                                      #Terminate the program

V_user = -(r_s/(2*r_i))+((l**2)/(2*(r_i**2)))-((r_s*(l**2))/(2*(r_i**3)))   #Calculate V(r) at the position of r chosen by user                                 
V_max_input = V_max + 0.5                                                   #Set max. max input E to be V_max + 0.5

print("\nMax. V(r)", round(V_max, 2))                                       #Print max possible value user can input
print("Enter value for particle energy (", round(V_user, 2), "≤ E <", round(V_max_input, 1), " or V(r)): ") #Show user range of input values for E
E_user = input()                                                            #Read user input from new line

if E_user == 'V(r)' or E_user == 'v' or E_user == 'V' or E_user == round(V_user, 2):    #Read user input to see if it matches any values corresponding V(r)
    E_user = V_user             #Set E_user to be V(r) value from using user defined r
    v_i = 0                     #Particle velocity must be zero as is on V(r) curve (no kinetic energy)
    
else:
    E_user = float(E_user)                                  #Convert user input from string type to float type

    if E_user == round(V_max_input, 1):                     #Condition to check if user has entered max. possible user E
        E_user = V_max_input - 0.01*V_max_input             #Ensure E is less than V_max_input
    
    if E_user > V_max_input or E_user < V_user:             #Condition to check if value is niether greater than max. possible user E or less than V(r) 
        E_user = input("Invalid value! Enter new value: ")  #Message prompting user to input new value
        
        if E_user == 'V(r)' or E_user == 'v' or E_user == 'V':  #Read user input to see if it matches any values corresponding V(r)
            E_user = V_user                                     #Set E_user to be V(r) value from using user defined r
            
        else: E_user = float(E_user)                            #Convert user input from string type to float type
        
        if E_user == round(V_max_input, 1):                     #Condition to check if user has entered max. possible user E
            E_user = V_max_input - 0.01*V_max_input             #Ensure E is less than V_max_input
    
        if E_user > V_max_input or E_user < V_user:     #Condition to check if value is niether greater than max. possible user E or less than V(r)
            print("Program terminated - please retry")  #Error message to ell user invalid entry
            sys.exit()                                  #Terminate the program
            
            
    if E_user > V_user:                                 #Condition to check if user input E is greater than V(r)
        direction = input("Inwards or outwards motion? (Enter + or -)\n")   #User prompt to input direction of particle motion
    
        if direction == "+":                            #Condition when user selects positive particle direction (away from black hole)
            v_i = math.sqrt((2*E_user) - (2*V_user))    #Calculate particle velocity from user inputs and assign positive value
        
        if direction == "-":                            #Condition when user selects negative particle direction (towards from black hole)
            v_i = -math.sqrt((2*E_user) - (2*V_user))   #Calculate particle velocity from user inputs and assign negative value
        

w_i = l/(r_i**2)    #Calculate angular velocity from user inputs                         



# ================================================================================= #



# ======= DEFINE FUNCTIONS ======= #

#These functions are later called in main iterative algorithm to find r, v, phi and w

def dvdt(r):
    return (l**2)/(r**3) - (r_s)/(2*(r**2)) - ((3*r_s*(l**2))/(2*(r**4)))   #Equation (1.8) in report

def dphidt(r):
    return (l/(r**2))               #From equation (1.4) in report

def drdt(v):
    return v                        #From substitution for equation (1.8) in report

def dwdt(r, drdt): 
    return ((-2*l)/(r**3))*drdt     #Equation (1.9) in report


# ================================== #
  


# ======= MAIN FOREST-RUTH ALGORITHM TO CALCULATE NEW VALUES ======= # 

i_max = ((3.34e-9)*c)                                       #Default setting of 1 second (not in metres of time)


def main_algorithm(t, r, v, phi, w, l, i_max):              #Input arguments to function
    
    start_time = time_.time()                               #Calculate start_time based on system time
    run_time = 0                                            #Reset run_time = 0 to prevent errors
    
    
    # ------- Define lists to hold values ------- #
    
    radius = [r]
    rad_vel = [v]
    time = [t]
    ang_vel = [w]
    angle = [phi]
    time_at_tp = []
    radius_at_tp = []
    round_angle = []
    
    
    
    # ------- Function to find optimal step size for r ------- #
    
    def step(r0, v0, h01):                                  
        h1 = h01; j = 1; Error = 1                          #Set inital values for fucntion
        step_r = [r0]                                       #Define an list to hold values later calculated
        
        while abs(Error) > 0.000001:                         #Error to be obtained through convergence
            r1 = r0 + (theta*(h1/2)*drdt(v0))               #Forest-Ruth algorithm to find r for given step h 
            v1 = v0 + (theta*h1*dvdt(r1))                   #...
            r2 = r1 + ((1 - theta)*(h1/2)*drdt(v1))         #...
            v2 = v1 + ((1 - (2*theta))*h1*dvdt(r2))         #...
            r3 = r2 + ((1 - theta)*(h1/2)*drdt(v2))         #...
            v3 = v2 + (theta*h1*dvdt(r3))                   #...
            step_r.append(r3 + (theta*(h1/2)*drdt(v3)))     #Append r values to list
        
            Error = (step_r[j]-step_r[j-1])/step_r[j]       #Calculate difference between subsequent r values
            Ep1 = abs(Error*100)                            #Ouput error as percentage
            h1 = h1/2                                       #Reduce step size and re-calculate r
            j = j+1
        return h1, Ep1                                      #Return program arguments
    
    i = 0
    h0 = (i_max/2)                                          #Define initial 'test' step size (in seconds as based off un-normalised i_max value)
    h, Ep = step(r, v, h0)                                  #Call optimised r step value
    while i <= i_max:                                       #Set maximum number of iterations
         
        
        
        
        # ------- Forest-Ruth to find r and v ------- #
        
        r0 = r 
        v0 = v
    
        r1 = r0 + (theta*(h/2)*drdt(v0))                    
        v1 = v0 + (theta*h*dvdt(r1))
    
        r2 = r1 + ((1 - theta)*(h/2)*drdt(v1))
        v2 = v1 + ((1 - (2*theta))*h*dvdt(r2))
    
        r3 = r2 + ((1 - theta)*(h/2)*drdt(v2))
        v3 = v2 + (theta*h*dvdt(r3))
    
    
        r = r3 + (theta*(h/2)*drdt(v3))                     #Re-set r value 
        v = v3                                              #Re-set v value 
        t = (t + h)                                         #Recalculate time by adding step onto initial value
        
        
        
    
        # ------- Forest-Ruth to find r and v ------- #
        
        phi0 = phi 
        w0 = w
    
        phi1 = phi0 + (theta*(h/2)*dphidt(r0))            
        w1 = w0 + (theta*h*dwdt(r1, drdt(v1)))
    
        phi2 = phi1 + ((1 - theta)*(h/2)*dphidt(r1))
        w2 = w1 + ((1 - (2*theta))*h*dwdt(r2, drdt(v2)))
    
        phi3 = phi2 + ((1 - theta)*(h/2)*dphidt(r2))
        w3 = w2 + (theta*h*dwdt(r3, drdt(v3)))
    
    
        phi = phi3 + (theta*(h/2)*dphidt(r3))               #Re-set phi value 
        w = w3                                              #Re-set v value 
                                                 
        
        
        # ------- Process values from Forest-Ruth algorithms ------- #
        
        # 1. Append new values to lists

        radius.append(r); time.append(t); rad_vel.append(v); ang_vel.append(w); angle.append(phi); round_angle.append(round(phi, 4))
        
    
        # 2. Add to loop counters and calculate run time
    
        i = i + h                               #Add step size to initial i to create next i value
        run_time = time_.time() - start_time    #Calculate the time elapsed since the start of the function


        # 4. Function exit conditions 
                
        if i > 10 and (round(angle[-1], 7) == round(angle[-10], 7)):    #Condition to check if rounded final values in angle list are equal
            i = i_max + 1                                               #Increase loop counter value to exceed i_max to exit loop
            print("Angle constant")                                     #Print message to user
        
        if (angle[-1] < 4*pi and (i+h) >= i_max):       #Condition to check last value in angle array and position in loop 
            i_max = i_max+h                             #Increase i_max to allow for more iteraitons to allow for angle condition to be met
            
        if radius[-1] <= 0.1:                           #Condition to stop code breakdown at r = 0 singularity
            i = i_max+1                                 #Increase loop counter value to exceed i_max to exit loop
            print("Radius too small")                   #Print message to user
                
        if run_time > 10:                               #Condition to check if run time hasn't exceeded 10s
            i = i_max+1                                 #Increase loop counter value to exceed i_max to exit loop
            print("* 10s Processing Time Exceeded *")   #Print message to user explaining run-time is too long
            return radius, time, rad_vel, ang_vel, angle, round_angle, time_at_tp, radius_at_tp    #Return the lists from the function

        
        
        
    # ----- Values to return from function ----- #
    
    return radius, time, rad_vel, ang_vel, angle, round_angle, time_at_tp, radius_at_tp            #Return the lists from the function



radius, time, rad_vel, ang_vel, angle, round_angle, time_at_tp, radius_at_tp = main_algorithm(t_i, r_i, v_i, phi_i, w_i, l, i_max)  #Call main algorithm with argument values 


# =============================================== #


        
# ======= OUTPUT POLAR PLOT ======= #
    
time_s = [x/c for x in time]                #Divide by c as there are n/c seconds in n metres of time

y_max = (max(radius)+((max(radius))*0.1))   #Set axes limits to ensure plot shows most relevant parts of orbits

fig = plt.Figure()                          #Define new plotting object       
ax = plt.subplot(1,1,1, polar=True)         #Create polar subplot 

ax.set(xlim=(0, 2*pi), ylim=(0, y_max))     #Set axes constraints
ax.set_theta_zero_location("N")             #Rotate plot so 0 angle is at top of plot
ax.plot(angle, radius, 'k')                 #Create polar plot and select black line color

circle = plt.Circle((0, 0), 1, transform=ax.transData._b, color='gray') #Create circle vector object to represent black hole
ax.add_artist(circle)                       #Impose circle object onto polar plot
#Code to create circle vector based on code by HYRY (2013)

label_position=ax.get_rlabel_position()     #Call positional argument of label and define new variable
ax.text(np.radians(label_position+20),ax.get_rmax()/2.,'r', rotation=label_position,ha='center',va='center', )  #Reference: tmdavison (2015)
#Adjust position of radial axis label

plt.show()                                  #Show polar plot
                                                     

# ================================= #



# ======= *** END OF PROGRAM *** ======= #

