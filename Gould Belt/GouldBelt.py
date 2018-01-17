import string
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord


#########################
##  PARSING
#########################

# We define the column labels here for the data file since it does not have a header:
cols = ['a', 'RA', 'DEC', 'Spectral Type', 'V', 'W-B', 'B-V', 'V-R', 'k', 'l']

# Define starts and ends of the data in the fixed-width columns:
colStarts = [0, 12, 19, 28, 38, 44, 52, 60, 68, 71]
colEnds = [10, 19, 26, 38, 44, 52, 60, 68, 71, 74]

# Read in .DAT file using our above settings:
data = ascii.read("WBVR.DAT", format='fixed_width_no_header', col_starts = colStarts, col_ends=colEnds, names=cols)

# Create dictionary keyed from Spectral Type
stars = {'O':{'Color':'Orange', 'V':[], 'W-B':[], 'B-V':[], 'V-R':[], 'RA':[], 'DEC':[], 'SPECTYPE':[], 'DISTANCE':[]},\
         'B':{'Color':'Blue', 'V':[], 'W-B':[], 'B-V':[], 'V-R':[], 'RA':[], 'DEC':[], 'SPECTYPE':[], 'DISTANCE':[]},\
         'A':{'Color':'Red', 'V':[], 'W-B':[], 'B-V':[], 'V-R':[], 'RA':[], 'DEC':[], 'SPECTYPE':[], 'DISTANCE':[]},\
         'F':{'Color':'Cyan', 'V':[], 'W-B':[], 'B-V':[], 'V-R':[], 'RA':[], 'DEC':[], 'SPECTYPE':[], 'DISTANCE':[]},\
         'G':{'Color':'Green', 'V':[], 'W-B':[], 'B-V':[], 'V-R':[], 'RA':[], 'DEC':[], 'SPECTYPE':[], 'DISTANCE':[]},\
         'K':{'Color':'Black', 'V':[], 'W-B':[], 'B-V':[], 'V-R':[], 'RA':[], 'DEC':[], 'SPECTYPE':[], 'DISTANCE':[]},\
         'M':{'Color':'Magenta', 'V':[], 'W-B':[], 'B-V':[], 'V-R':[], 'RA':[], 'DEC':[], 'SPECTYPE':[], 'DISTANCE':[]}}

# Absolute magnitudes for B type stars 
MV = {'B':[-4, -3.25, -2.45, -2.0, -1.6, -1.2, -0.9, -0.6, -0.25, 0.2]}

# Fill our dictionary with information from the ASCII table:
for star in data:
    # Remove lowercase letters from spectral type strings:
    specValue = str(star['Spectral Type']).translate(None, string.ascii_lowercase)[0:2]
    # print(specValue)
    # Skip empty/missing data:
    if '-' in star:
        continue
    
    # If the spectral type is in our dictionary, we add the star to its category:
    if specValue[0] in stars.keys():
        if(len(specValue) > 1 and specValue[1].isdigit()):
            stars[specValue[0]]['V'].append(float(star['V']))
            stars[specValue[0]]['W-B'].append(float(star['W-B']))
            stars[specValue[0]]['B-V'].append(float(star['B-V']))
            stars[specValue[0]]['V-R'].append(float(star['V-R']))
            stars[specValue[0]]['RA'].append(star['RA'])
            stars[specValue[0]]['DEC'].append(star['DEC'])
            stars[specValue[0]]['SPECTYPE'].append(specValue)
            if(specValue[0] == 'B'):
                #Distance Formula:  m - M = 5*log(d/10) 
                stars[specValue[0]]['DISTANCE'].append(((10**(float(star['V']-MV[specValue[0]][int(specValue[1])])/5))/1000)*u.kpc)




#########################
##  PLOTTING
#########################

# W-B VS B-V
for color in stars:
    plt.scatter(stars[color]['B-V'], stars[color]['W-B'],\
                c=stars[color]['Color'], s=1, label=color)

plt.title("W-B vs B-V for WBVR Data")
plt.xlabel("B-V")
plt.ylabel("W-B")
plt.legend()
plt.show()


# B-V VS V-R
for color in stars:
    plt.scatter(stars[color]['V-R'], stars[color]['B-V'],\
                c=stars[color]['Color'], s=1, label=color)

plt.title("B-V vs V-R for WBVR Data")
plt.xlabel("V-R")
plt.ylabel("B-V")
plt.legend()
plt.show()


#Aitoff Projection (All Stars)
plt.subplot(111, projection="aitoff")
plt.title("Aitoff Projection of WBVR Data", y=1.08)
plt.grid(True)
for color in stars:
    coords = SkyCoord(stars[color]['RA'], stars[color]['DEC'], unit=(u.hourangle, u.deg)).galactic
    l_rad = coords.l.wrap_at(180*u.deg).radian
    b_rad = coords.b.radian    
    plt.scatter(l_rad, b_rad, c=stars[color]['Color'], s=0.2)

plt.subplots_adjust(top=0.95, bottom=0.0)
plt.show()


#Aitoff Projection (B Stars)
plt.subplot(111, projection="aitoff")
plt.title("Aitoff Projection of WBVR Data for B Class Stars", y=1.08)
plt.grid(True)
color = 'B'
coords = SkyCoord(stars[color]['RA'], stars[color]['DEC'], unit=(u.hourangle, u.deg)).galactic
l_rad = coords.l.wrap_at(180*u.deg).radian
b_rad = coords.b.radian    
plt.scatter(l_rad, b_rad, c=stars[color]['Color'], s=0.2)

plt.subplots_adjust(top=0.95, bottom=0.0)
plt.show()


#Cartesian Projection (B Stars)
color = 'B'
coords = SkyCoord(ra=stars[color]['RA'], dec=stars[color]['DEC'], distance=stars[color]['DISTANCE'], unit=(u.hourangle, u.deg)).cartesian
x = coords.x
y = coords.y
plt.scatter(x, y, c=stars[color]['Color'], s=0.2)
plt.title("Cartesian Projection of WBVR Data for B Class Stars")
plt.xlabel("V-R")
plt.ylabel("B-V")
plt.show()
