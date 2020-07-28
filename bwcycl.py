import math
#######################
## SOURCE SETTINGS
#######################
Bfield=4e12 # in units of G
massns=1.4 # in units of solar mass
radiusns=1 # in units of 10 km
#######################
mdot=1e16 # in units of g/s
#######################
## PHYSICAL CONSTANTS (CGS units) 
#######################
mpro=1.6726e-24  # mass of the proton 
clight=3e10
electronmass=9.109e-28
sigma_th = 6.65e-25
pigreco=3.14
GravConst=6.674e-8 
#######################
## FITTING PARAMETERS
#######################
xi=2.07
delta=6.2
ro=5.5e2 # in cm
kTekeV=4.8 # in keV
#######################
## derived QUANTITIES
#######################
kTe=kTekeV*1.16e7
alpha= 0.2 * massns * xi / radiusns
sigma_par = (math.pi*ro*mpro*clight/(mdot*xi))**2 / sigma_th**2
sigma_ave=(alpha*511)/(3*kTekeV*delta)*(sigma_par)
lum_crit = 2.7e37 * (ro/1e6) * massns / math.sqrt(sigma_par)
temp_mound = 2.32e3 * mdot**(2/5)*ro**(-2/3)
temp_mound_kev=temp_mound/11600e3
rho_mound=4.05e-12 * temp_mound**(7/4) * ro**(-1/2)
ztrap=ro*xi/(2*math.sqrt(sigma_par))
vel_inflow = 7.86e10*mdot*ro**(-3/2)*temp_mound**(-7/4) / clight
tau_therm = 2.64e28 * mdot* (radiusns*1e6) / (massns*2e33*ro**(3/2)*temp_mound**(7/4)*xi)
z_therm = 5.44e15 * mdot* (radiusns*1e6) / (massns*2e33*ro*temp_mound**(7/2)*xi*(sigma_par*sigma_th))
c1=4*GravConst*massns*2e33*ro*xi/(alpha*clight**2*(radiusns*1e6)**2* math.sqrt(sigma_par))
zmax=radiusns*1e6/2 * ( math.sqrt(1+c1) - 1) 
tau_max = sigma_par**(1/4) * math.sqrt(2*zmax/(alpha*xi*ro))
tau_trap = 1/math.sqrt(alpha)
energy_abs= 6.08e12*kTe**(-7/4)*math.sqrt(rho_mound)
jai=mdot/(pigreco*ro**2)
#######################
#######################
print("The local mass accretion rate is {:5.2e} g cm-2 s-1".format(jai))  
print("The critical luminosity for this source is {:5.2e} erg s-1".format(lum_crit))  
print("The altitude of the accretion column is {:5.2e} cm".format(zmax))

print("The longitudinal cross-section area is {:5.2e} sigma_thomson".format(sigma_par))  
print("The angle-averaged cross-section area is {:5.2e} sigma_thomson".format(sigma_ave))  

print("Mound characteristics")  
print("The mound temperature is {:5.2e} keV".format(temp_mound_kev ))  
print("The mound density is {:5.2e} g cm-3".format(rho_mound))  
print("The inflow velocity at the mound top {:5.2e} c".format(vel_inflow))  
print("The optical depth is {:5.2e}".format(tau_therm ))  
print("The altitude zth is {:5.2e} cm".format(z_therm )) 

print("Photons trapping")  
print("The altitude z_trap is {:5.2e} cm".format(ztrap )) 
print("The trapping optical depth is {:5.2e}".format(tau_trap )) 
print("The maximum optical depth of the column is {:5.2e}".format(tau_max)) 
#print("The thermal energy free-free cut-off is {:5.2e} (units of kTe)".format(energy_abs)) 

