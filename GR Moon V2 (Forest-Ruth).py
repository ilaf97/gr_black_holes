#-- DEFINE LIBRARIES --
import matplotlib.pyplot as plt
import statistics as sts
import math


#-- DEFINE INPUT VARIABLES/CONSTANTS --
c = 299792458               #Speed of light
r_s = 1                     #Schwarschild radius
pi = math.pi                #pi
G = 6.67408e-11/(c**2)      #Graviational constant (normalised)
M = 1/(2*G)                 #Mass of Earth in r_s = 1 units 
M_earth = 5.972e+24               #Mass of main gravitating body (Earth)
sf = M/M_earth              #Scaling factor
m_moon = 7.3476709e+22
m = sf*7.3476709e+22           #Mass of moon
l = ((2.89805e+34)/(m))/c       #Speicfic angular momentum of Moon (See Newtonian Earth program for reference) (normalised)
theta = 1/(2-2**(1./3.))    #Value always used to compute Forrest-Ruth algorithm (constant)
r = (sf**2)*405.504e+06 
phi = 0            #Radial distance (perigee of moon) - See Newtonian Earth program for references
v = 0.00                    #Radial velocity 
t = 0.00                    #Starting time
w = (l/r**2)                #Agnular velocity (normalised as l already divided by c


#-- DEFINE FUNCTIONS --

def dvdt(r):
    return (l**2)/(r**3) - (r_s)/(2*(r**2)) - ((3*r_s*(l**2))/(2*(r**4))) #In c=1 units

def dphidt(r):
    return (l/(r**2))

def drdt(v):
    return v

def dwdt(r, drdt): #Use v3 to ensure correct value of v used
    return ((-2*l)/(r**3))*drdt


#-- DEFINE LISTS TO HOLD VALUES --
radius = [r]
rad_vel = [v]
time = [t]
ang_vel = [w]
period = []
angle = [phi]
time_at_tp = []
radius_at_tp = []


#-- MAIN ALGORITHM TO CALCULATE NEW VALUES -- see http://physics.ucsc.edu/~peter/242/leapfrog.pdf

error_analysis = 'n' #<<<< INDICATE IF ERROR ANALYSIS REQUIRED Y/N

#-- Function to find optimal step size for r--
def step(r0, v0, h01):
    h1 = h01; j = 1; E = 1
    step_r = [r0]
    
    while abs(E) > 0.001:
        r1 = r0 + (theta*(h1/2)*drdt(v0))
        v1 = v0 + (theta*h1*dvdt(r1))
        r2 = r1 + ((1 - theta)*(h1/2)*drdt(v1))
        v2 = v1 + ((1 - (2*theta))*h1*dvdt(r2))
        r3 = r2 + ((1 - theta)*(h1/2)*drdt(v2))
        v3 = v2 + (theta*h1*dvdt(r3))
        step_r.append(r3 + (theta*(h1/2)*drdt(v3)))
    
        E = (step_r[j]-step_r[j-1])/step_r[j] #Calculate difference between subsequent r valuse
        Ep1 = E*100 
        h1 = h1/2
        j = j+1
    return h1, Ep1


i = 0; i_max = (sf**2)*365*(24*3600)*c #Maximum number of iterations. Note: time is in metres, so need to un-normalise by *c 
h0 = (i_max/2) #Define initial 'test' step size (in seconds as based off un-normalised i_max value)
h, Ep = step(r, v, h0) #Call optimised r step value
print("Step =", h, " --> Converged to", round(Ep, 3), "%") #Print function output

k = 0    #List index loop counter initial value
while i <= i_max:  #Set maximum number of iterations
     
    #-- First stage to find r and values --
    #Rename r and v for simplicity in algorithm
    r0 = r 
    v0 = v

    r1 = r0 + (theta*(h/2)*drdt(v0)) #Theta has no physical meaning - it's just a constant to make the algorithm work 
    v1 = v0 + (theta*h*dvdt(r1))

    r2 = r1 + ((1 - theta)*(h/2)*drdt(v1))
    v2 = v1 + ((1 - (2*theta))*h*dvdt(r2))

    r3 = r2 + ((1 - theta)*(h/2)*drdt(v2))
    v3 = v2 + (theta*h*dvdt(r3))


    r = r3 + (theta*(h/2)*drdt(v3))
    v = v3
    t = (t + h)
    
    
    #-- Second stage to find phi values --
    #Rename r and v for simplicity in algorithm
    phi0 = phi 
    w0 = w

    phi1 = phi0 + (theta*(h/2)*dphidt(r0)) #Theta has no physical meaning - it's just a constant to make the algorithm work 
    w1 = w0 + (theta*h*dwdt(r1, drdt(v1)))

    phi2 = phi1 + ((1 - theta)*(h/2)*dphidt(r1))
    w2 = w1 + ((1 - (2*theta))*h*dwdt(r2, drdt(v2)))

    phi3 = phi2 + ((1 - theta)*(h/2)*dphidt(r2))
    w3 = w2 + (theta*h*dwdt(r3, drdt(v3)))


    phi = phi3 + (theta*(h/2)*dphidt(r3))
    w = w3
    
    T = (2*pi)/(w*((3600*24)*c)) #Calculate period (normalised)

    radius.append(r); time.append(t); rad_vel.append(v); ang_vel.append(w); period.append(T); angle.append(phi) #Append new values to lists
    
    if ((radius[k-2] < radius[k-1] and radius[k] < radius[k-1]) #Maximum
    or (radius[k-2] > radius[k-1] and radius[k] > radius[k-1] )): #Minimum
        time_at_tp.append(time[k-1])    #Append time at turning points to earlier-defined array
        radius_at_tp.append(radius[k-1])  #Append radius at turning points to earlier-defined array
   # Based on code taken from https://stackoverflow.com/questions/19936033/finding-turning-points-of-an-array-in-python


    i = i + h #Add step size to initial i to create next i value
    k = k + 1

#-- End of loop --
     
    
#-- OUTPUT PLOTS --
time_days = [x/((3600*24)*c) for x in time] #Create new array holding all values from time but in units of days (normalised)

plt.xlabel('Time (days)'); plt.ylabel('Orbital Distance (m)'); plt.title('Orbital Distance vs Time') #Define plot lables
plt.plot(time_days, radius) #Define things to plot
plt.show() #Show plots

plt.xlabel('Time (days)'); plt.ylabel('Angle (radians)'); plt.title('Angle vs Time') #Define plot lables
plt.plot(time_days, angle) #Define things to plot
plt.show() #Show plots

plt.polar(angle, radius); plt.title ("Polar Orbital Plot") #Define things to plot on polar plot
plt.show() #Show polar plot

plt.polar(angle, radius); plt.title("Reduced r Axis") #Define things to plot on polar plot
plt.ylim(3.5e+08,4.1e+08)
plt.show() #Show polar plot



#-- ACCURACY ANALYSIS --
if error_analysis == 'y':
    
    times_after_1_period = time_at_tp[1::2] #Takes every even-valued indexed entry in time array, starting with 2nd entry
    
    j = 1   #Initial loop counter value
    computed_periods = []   #Define empty array to append values to next
    while j <= (len(times_after_1_period)-1):   #Define condition that cycles through all available index values in times_after_1_period array
        computed_periods.append(times_after_1_period[j]-times_after_1_period[j-1]) #Calculate periods 
        j = j + 1 #Ammend loop condition
    
    computed_periods = [x/(3600*24*c) for x in computed_periods] #Create new array holding all values from time but in units of days
    
    
    #-- Apply statistical analysis to variables --
    sd_p = sts.stdev(computed_periods)
    sd_maxr = sts.stdev(radius_at_tp[2::2]) #Assuming starting at max. distance from graviating object
    sd_minr = sts.stdev(radius_at_tp[1::2])
    ave_p = sts.mean(computed_periods)
    ave_maxr = sts.mean(radius_at_tp[2::2])
    ave_minr = sts.mean(radius_at_tp[1::2])
    "{:.2e}".format(sd_maxr); "{:.2e}".format(sd_maxr) #Format min and max r into scientific notation

    
    #-- Ouput analysis into console --
    table = [['Period', round(sd_p,2), round(ave_p,2)], ['Max. Dist.', round(sd_maxr,2), "{:.2e}".format(ave_maxr)], ['Min. Dist.', round(sd_minr,2), "{:.2e}".format(ave_minr)]]
    print(table)
    print(' ')
    print('st.dev as % of mean period =', round(((sd_p/ave_p)*100),5))
    print('sd.dev as % of mean max r =', round(((sd_maxr/ave_maxr)*100), 5))
    print('sd.dev as % of mean min r =', round(((sd_minr/ave_minr)*100), 5))
 

