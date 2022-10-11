import numpy as np
from scipy import interpolate



def convert_color2temp(cl, band, feh):


	# Gaia color vs temperature transformation from :	
	# https://arxiv.org/pdf/2106.03882.pdf

	# band can be 'BPRP'

	if band == 'BPRP':
		b = [0.4929, 0.5092, -0.0353, 0.0192, -0.0020, -0.0395]

	elif band == 'BPG': 
		b = [0.5316, 1.2452, -0.4677, 0.0068, -0.0031, -0.0752]
	

	theta = b[0] + b[1]*cl + b[2]*cl**2 +b[3]*feh + b[4]*feh**2 + b[5]*feh*cl

	teff = 5040/theta

	return(teff)

def convert_temp2color(temp, band, feh):

	if band == 'BPRP':
		clmin = 0.39
		clmax = 1.50
	elif band == 'BPG':
		clmin = 0.13
		clmax = 0.69

	cls = np.arange(clmin, clmax, 0.01)
	temps = [ convert_color2temp(cl, band, feh)  for cl in cls ]

	print("Temperatures for the color range %f to %f"%(clmin, clmax))
	print(temps)	

	func = interpolate.interp1d(temps, cls)
	color = func(temp)

	if color<=clmin or color>=clmax:
		print("Color", color, " for temperature", temp, " is out of acceptable range.")	

	return(color)
    
 






#temp = 7500.
#band = 'BPRP'
#feh = -1
#print(convert_temp2color(temp, band, feh))

