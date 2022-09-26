import math
import csv
from tkinter import N
#####################################################################################
#
#   using bisection method to iteratively approximate first the wet bulb temp       #
#   and subsequently the dry bulb temp.                                             #
#   inputs from measurements are: T-AIR[°C], REL-HUM[%], ETR[g/m³ in KV per second]
#
#####################################################################################


# CONSTANTS
gas_const = 461.52 #J/kgK spedific gas constant for water


# INPUTS - choose either hardcoded or interface by un/commenting

# INTERFACE
#rel_hum = float(input("relative humidity:")) #63.282795 #%
#t_air = float(input("air temp:")) #24.060555 # °C
#ET = float(input("evapotranspiration per second:")) #0.731875837596899 # g/m³KV pro sekunde

# HARDCODED
rel_hum = 55.1136 #%
t_air = 25.4775 # °C
ET = 0.031396 # g/m³KV pro sekunde

#buck saturation vapor pressure at temp [bar]
def buck(temp):
    saturation_vapor_pressure = (0.61121*math.exp((18.678-temp/234.5)*(temp/(257.14+temp))))/100
    return saturation_vapor_pressure


#########################################################
#                                                       #
#       APPROXIMATION OF WET BULB TEMPERATURE           #
#                                                       #
#########################################################

p_vs=buck(t_air)

def error_fct_wb(temp):
    p_vs_wb=buck(temp)#p_vs at twb
    p_v = p_vs_wb-(1.8*(t_air-temp))/2700#apjohn equation for vapor pressure f(twb,tdb) [bar]
    error = (rel_hum * p_vs)/100 - p_v #error fct comparing two approches to solve p_v; apjohn versus p_vs*rel.hum
    return error

# BISECTION METHOD
#low (temp1) and high (temp2) limit for bisection method to calc twb
temp1=-25.0
temp2=30.0
convergence=1.0E-8
max_iteration = 200
iteration = 0

while iteration <= max_iteration:
    temp3 = (temp1+temp2)/2
    error_temp1 = error_fct_wb(temp1)
    error_temp3 = error_fct_wb(temp3)
    #print(error_temp1, error_temp3, temp1, temp2,temp3)
    if (error_temp3 * error_temp1) < 0:
        temp2=temp3
    else:
         temp1=temp3   
    REC = abs(error_temp3/temp3)
    if REC < convergence:
        print("TWB", temp3, "at iteration", iteration, "was found, with convergence of", REC)
        break
    else:
         iteration+=1
print("method ended")

twb = temp3
p_vs_wb=buck(twb) 


#########################################################
#                                                       #
#       APPROXIMATION OF DRY BULB TEMPERATURE           #
#                                                       #
#########################################################


def error_fct_tdb(temp):
    p_vs_db=buck(temp)
    max_hum_ET =(p_vs_db*100000/(gas_const*(temp+273.15)))*1000 # ideal gas law to calc maximum humidity [g/m³]
    max_hum = (p_vs*100000/(gas_const*(t_air+273.15)))*1000
    abs_hum = max_hum * rel_hum /100 #calc for absolute humidity [g/m³]
    abs_hum_ET = abs_hum + ETR # add water from ET
    rel_hum_ET = abs_hum_ET/max_hum_ET #calc for relative humidity inkl ET, NOT in %
    p_v = p_vs_wb-(1.8*(temp-twb))/2700 #apjohn equation for vapor pressure f(twb,tdb) [bar]
    error = (rel_hum_ET * p_vs_db) - p_v #error fct comparing two approches to solve p_v; apjohn versus p_vs*rel.hum
    return error


 # BISECTION METHOD

with open('results_sec.csv', 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['ETR,g/m³s', 'sec', 'tdb,°C', 'temp_depression,°C', 'twb,°C'])
    duration = list(range(200)) #Duration of Evaporation in seconds
    
    for i in duration:
        temp1=-20.0 # low limit twb
        temp2=30.0  # high limit twb
        SEC=duration[i]
        ETR=ET*SEC
        iteration=0
        convergence=1.0E-9
        max_iteration = 100
        print('\n','START: duration', duration[i],'ETR', ETR, 'sec', SEC)
        
        while iteration <= max_iteration:
            temp3 = (temp1+temp2)/2
            error_temp1 = error_fct_tdb(temp1)
            error_temp3 = error_fct_tdb(temp3)
            #print(error_temp1, error_temp3, temp1, temp2,temp3)
	    #print('iteration',iteration,'ETR', ETR, 'SEC', SEC)
            if (error_temp3 * error_temp1) < 0:
                temp2=temp3
            else:
                temp1=temp3
                REC = abs(error_temp3/temp3)
		#print(REC)
            if REC < convergence:
                print("TDB", temp3, "at iteration", iteration, "was found, with convergence of", REC)
                iteration=max_iteration
            else:
                iteration+=1
            iteration+=1
	#print("method ended")
			
        tdb=temp3
        temp_depression = t_air - tdb
        print("temp depression is", round(temp_depression,2))
        
    

#########################################################
#                                                       #
#                   PLAUSIBILITY CHECKS                 #
#                                                       #
#########################################################

        if tdb<twb:
            print("TDB < TWB, not plausible, decrease duration [sec]")
		
        p_vs_db=buck(tdb)
        max_hum_ET =(p_vs_db*100000/(gas_const*(tdb+273.15)))*1000 # ideal gas law to calc maximum humidity [g/m³]
        max_hum = (p_vs*100000/(gas_const*(t_air+273.15)))*1000
        abs_hum = max_hum * rel_hum /100
        abs_hum_ET = abs_hum + ETR

        if abs_hum_ET > max_hum_ET:
            print("absolute humidity exceeds maximum humidity, not plausible, decrease duration [sec]")
        else:
            print("temp depression is", round(temp_depression,2))
            
        writer.writerow([ETR, SEC, tdb, temp_depression, twb])


