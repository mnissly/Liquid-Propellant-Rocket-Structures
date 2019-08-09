# File: EngineSim.py
# Name: Micah Nissly
# Description:
#	This program performs calculations for propulsion system sizing based on 
#	preliminary prupulsion data. Uses a mass model of the system to determine
#	the CG of the system as a function of time. Outputs a series of engine files
#	in a text document for use in OpenRocket. Code was developed from a script 
#	written by the Portland State Aerospace Society.

# FUNCTIONS ------------------------------------------------------------------#
def tankLength(m,r,rho): # determines necessary tank length based on propellant masses
	l = m / (rho*pi*r**2)
	return l

def HeVolume(pProp,vProp,pHe): # determines the initial pressure of Helium tank
	pHef = pProp + FINAL_HE_PRESSURE_ADDON
	vHe = (pProp/(pHe-pHef))**(1/gammaHe)*vProp
	return vHe

def tankThickness(p,r,ys,fs): # determines tank thickness based on tank pressure
	designStress = ys/fs
	thickness = p * r / designStress
	if designStress / p < 10:
		print("Invalid thin-walled pressure vessel assumption!")
	return thickness


def tankMass(l,ro,ri,thickness): # determines tank mass
	volume = pi*(ro**2 - ri**2)*l + 2*pi*ri**2 * thickness
	mass = volume * tank['rho']
	return mass

def mass(t): # determines remaining propellant mass as a function of time
	m = massProp - mdot*t
	m = round(m,3)
	return m

def centerOfMass(t): # calculates the overall CG of the propulsion system as a function of time
	massOx_t = massOx - mdotOx*t
	massFuel_t = massFuel - mdotFuel*t
	cmOx = cmOxTank + lenOx/2.0 - thickness - tankLength(massOx_t,ri,density_o[oxType])/2.0
	cmFuel = cmFuelTank + lenFuel/2.0 - thickness - tankLength(massFuel_t,ri,density_f[fuelType])/2.0
	cmProp = ((massOx_t*cmOx) + (massFuel_t*cmFuel))\
				/(massOx_t + massFuel_t + 0.000001)
	cm = ((dryCM*dryMass) + (cmProp*(massOx_t + massFuel_t))) \
			 / (dryMass + massOx_t + massFuel_t)
	return cm

# SETUP ----------------------------------------------------------------------#
from math import pi, log, sqrt
import os
file_head = []
data_set = []
eng_head = []
# Physics/Math
g0		 =	  9.80665		# m/s^2		Standard Gravity

# Chemistry
density_f	 = {'RP1': 820, 'METH': 423, 'ETH': 790, 'LH2': 71, 'UDMH': 793}
density_o	 = {'LOX': 1142, 'FLOX': 447.6, 'LF': 150, 'N2O4': 1510}
gammaHe     = 1.67
FINAL_HE_PRESSURE_ADDON = 0.517107

# Other (THESE ARE ROUGH ESTIMATES)
massEng   =    2.0			# kg
lenEng    =    0.325		# m
massPlumb =    4.0			# kg
lenPlumb  =    0.300        # m
gaps	  =    0.300		# m
massGaps  =    1.0			# kg

mfg      =  input("Manufacturer: ")

# Tank
ro 		 =	float(input("Propellant Tank Radius (in): "))  #        	Radius of the tanks (ID of the rocket)
rou 	 =	float(input("Helium Tank Radius (in): "))      #        	Radius of the tanks (ID of the rocket)
ro       *=  0.0254								           # 			convert to meters
rou      *=  0.0254

Al		 = {'rho': 2700.0,	# kg/m^3	Density
		   'YS': 276	}	# MPa		Yield Strength
pressure =  float(input("Tank Pressure (MPa): "))
pressure_u = float(input("Helium Tank Pressure (MPa): "))
factorOfSafety = 1.5
thickness = tankThickness(pressure,ro,Al['YS'],factorOfSafety)
thickness_u = tankThickness(pressure_u,rou,Al['YS'],factorOfSafety)
ri       = ro - thickness
riu      = rou - thickness_u
print("Tank thickness:        %0.2f mm" % (thickness*1000))
print("Helium Tank thickness: %0.2f mm" % (thickness_u*1000))

newEng = 'y'
while newEng != 'n':
	code = input("Rocket Code: ")
	
	# Propellants
	print("Oxidizer options:",density_o.keys())
	oxType = input('Oxidizer: ')
	print("Fuel options:",density_f.keys())
	fuelType = input('Fuel: ')

	# Assumptions
	Isp		 =	float(input("Specific Impulse (s): "))	# s 		Specific Impulse
	OF		 =	float(input("O/F Ratio: "))			    #			O/F ratio

	# Variables
	thrust    =    float(input("Thrust (N): "))     # N 		Thrust of engine
	burnTime  =	   input("Burn Time (enter as list separated by spaces): ")	# s Duration of the burn
	burnTime  = burnTime.split()
	tank      =	   Al 							# 			Chosen from tank selection above

	for b in range(len(burnTime)):
		# MASS AND FLOW --------------------------------------------------------------#
		# Total mass flow rate
		mdot = thrust / (g0*Isp)

		# Mass flow rate for each liquid
		mdotOx = mdot / (1 + (1/OF))
		mdotFuel = mdot / (1 + OF)

		# Propellant Mass
		massOx = float(burnTime[b]) * mdotOx
		massFuel = float(burnTime[b]) * mdotFuel
		massProp = massOx + massFuel

		print("Total Propellant Mass: %0.2f kg" % massProp)
		print("Ox mass:               %0.2f kg" % massOx)
		print("Fuel mass:             %0.2f kg" % massFuel)
		print("Mass flow:             %0.2f kg/s" % mdot)
		print("Ox flow:               %0.2f kg/s" % mdotOx)
		print("Fuel flow:             %0.2f kg/s" % mdotFuel)

		# TANK GEOMETRY --------------------------------------------------------------#

		lenOx = tankLength(massOx,ri,density_o[oxType]) + 2*thickness
		lenFuel = tankLength(massFuel,ri,density_f[fuelType]) + 2*thickness
		vHe = HeVolume(pressure,(lenOx+lenFuel)*(pi*ri**2),pressure_u)
		lenHe = (vHe/(pi*riu**2)) + 2*thickness_u		# m
		length = lenHe + lenOx + lenFuel + 2*gaps + lenPlumb + lenEng
		standardLenHe = lenHe * 39.3701
		standardLenOx = lenOx * 39.3701
		standardLenFuel = lenFuel * 39.3701
		standardLength = length * 39.3701
		print()
		print("Helium tank length:    %0.3f in" % standardLenHe)
		print("Ox tank length:        %0.3f in" % standardLenOx)
		print("Fuel tank length:      %0.3f in" % standardLenFuel)
		print("System length:         %0.3f in" % standardLength)

		# TANK MASS ------------------------------------------------------------------#

		massOxTank = tankMass(lenOx,ro,ri,thickness)
		massFuelTank = tankMass(lenFuel,ro,ri,thickness)
		massHelium = tankMass(lenHe,rou,riu,thickness_u)

		print()
		print("Ox tank mass:          %0.2f kg" % massOxTank)
		print("Fuel tank mass:        %0.2f kg" % massFuelTank)
		print("Helium tank mass:	  %0.2f kg" % massHelium)

		# SYSTEM ---------------------------------------------------------------------#

		dryMass = massEng + massPlumb + massOxTank + massFuelTank + massHelium + 2*massGaps

		print()
		print("Dry Mass:              %0.1f kg" % dryMass)

		# lox tank cm is in the center of the tank
		cmHe = lenHe / 2.0
		cmGap1 = lenHe + gaps / 2.0
		cmOxTank = lenHe + gaps + lenOx / 2.0
		cmGap2 = cmOxTank + lenOx / 2.0 + gaps / 2.0
		cmFuelTank = cmGap2 + gaps / 2.0 + lenFuel / 2.0
		cmPlumb = cmFuelTank + lenFuel / 2.0 + lenPlumb / 2.0
		cmEngine = cmPlumb + lenPlumb / 2.0 + lenEng / 2.0
		dryCM = sum([cmHe*massHelium, cmGap1*massGaps, cmOxTank*massOxTank, cmGap2*massGaps, \
					cmFuelTank*massFuelTank, cmPlumb*massPlumb, cmEngine*massEng])
		dryCM /= dryMass
		print("Dry CM:                 %0.3f m" % dryCM)	

		# Nar Letter Code
		impulse = thrust * float(burnTime[b])
		print("Total impulse: %0.0f N.s" % impulse)
		nari = int(log(impulse/2.5)/log(2))
		narPercent = impulse / (2.5*2**(nari+1))
		print('NAR: "%s" (%0.0f%%)' % (chr(66+nari),narPercent*100))

		eng_head.append("""
    		<engine  mfg="{mfg}" code="{code}" Type="Liquid" dia="{diameter}" len="{length}"
    		initWt="{total_mass}" propWt="{M_prop}" delays="0" auto-calc-mass="0" auto-calc-cg="0"
    		avgThrust="{thrust}" peakThrust="{thrust}" throatDia="0." exitDia="0." Itot="{impulse}"
    		burn-time="{burn_time}" massFrac="{m_frac}" Isp="{Isp}" tDiv="10" tStep="-1." tFix="1"
    		FDiv="10" FStep="-1." FFix="1" mDiv="10" mStep="-1." mFix="1" cgDiv="10"
    		cgStep="-1." cgFix="1">
    		<data>
		""".format(**{'diameter': ro*2*1000,
		              'length': length*1000,
		              'total_mass': (massProp+dryMass)*1000,
		              'M_prop': massProp*1000,
		              'thrust': thrust,
		              'burn_time': burnTime[b],
		              'm_frac': dryMass/(dryMass+massProp)*1000,
		              'impulse': impulse,
		              'Isp': Isp,
		              'code': code + "_" + str(burnTime[b]),
		              'mfg': mfg,
		    }))

		data = []
		n = 100
		res = float(burnTime[b])/float(n-1)
		for i in range(n):
		    t = i * res
		    data.append('     <eng-data  t="{t}" f="{thrust}" m="{mass}" cg="{cg}"/>\n'.format(**{
		        't': t,
		        'thrust': thrust,
		        'mass': mass(t)*1000,
		        'cg': centerOfMass(t)*1000,
		    }))
		data_set.append(data)

	newEng = input("Type 'y' to enter a new Engine, type 'n' to finish: ")

eng_tail = """
	</data>
  </engine>"""

file_head = """
<engine-database>
  <engine-list>"""

file_tail = """
</engine-list>
</engine-database>"""

fileName = mfg + ".rse"

prefix = os.path.join(os.getenv("APPDATA"), "OpenRocket/ThrustCurves/")

with open(os.path.join(prefix, fileName), 'w') as eng:
    eng.write(file_head)
    for a in range(len(eng_head)):
    	eng.write(eng_head[a])
    	for d in data_set[a]:
        	eng.write("		")
        	eng.write(d)
    	eng.write(eng_tail)
    eng.write(file_tail)

print("Reopen OpenRocket to run simulation")
