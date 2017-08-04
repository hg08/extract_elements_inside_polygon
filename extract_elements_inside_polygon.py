#=======================================================================================================================================
# Codename: check4.py
# Function: to extract the points inside a polygon and count the number of these points. 
# 2014/12/16 
# Author: Gang Huang (gang@uni-mainz.de)

# Run: python3.2 <codename>  < input_check3   (or use python3.4 instead of python3.2)
# The first line of the input: extral input file 1----transformed position file of Nitrogen atom : eg 'N_990.xyz' 
# The second line of the input: extral input file 2----transformed position file of Gold atom : eg  'au_990.xyz'
# The third line of the input: the number of N atoms
# The forth line of the input: the index of the five Au atoms at the vertices
# Note: MODIFY the input file IF one want to calculate for another snapshort!!

# Result: pol_**.dat, inside_**.xyz, and ts_**.dat 
# pol_**.dat: transformed coordinates (y,z) of the Au atoms at verticesss 
# inside_**.xyz: includes all the Nitrogen atoms whose transformed coordinates(y,z) are inside the polygen (.xyz)
# ts_**.dat: also includes all the Nitrogen atoms whose transformed coordinates(y,z) are inside the polygen (without the element symbols) 

# Drawback: The efficiency is low. check5.py is MUCH BETTER than check4.py!!
#========================================================================================================================================

#=======
#Library
#=======
import math
import linecache

#==========
# CONSTANTS
#==========
thi=((math.pi*2)/float(5))

#================
# Basic questions
#================
# Type in the name of the two snapshort xyz file for Nitrogen atoms and Gold atoms!
# ming and ming2 are temporay varibles
ming =input("What is the position file of Nitrogen atoms (one step):\n")   
ming2=input("What is the position file of Gold atoms (one step):\n")
#print(ming)
#print(ming2)
name = ming
name2= ming2
#
name_part0=name.split(".")[0]
# Type in the number of N atoms in sthe system
num_N=input("How many Nitrogen atoms in the system (int):\n")
n_step=int(num_N)
# Type in the index of the atoms at the five vertices
line=input("What are the index of the five Au atoms at verties (separated by space) :\n")
line_split=line.split(' ')
au_0 =int(line_split[0]) 
au_1 =int(line_split[1]) 
au_2 =int(line_split[2])
au_3 =int(line_split[3])
au_4 =int(line_split[4])

#======================================================
# Read position of the edge of the gold nanoseed: 
# The (y,z) coordinates of five Au atoms at the votices.
#======================================================
y_vertex=[]
z_vertex=[]
y_ave=0
z_ave=0
for i in [au_0 +3,au_1 +3,au_2 +3,au_3 +3,au_4 +3]:      # eg: Au atom  with index 305 conrespond to the 308th line in the xyz file
        au_vertex = linecache.getline(name2, i).split()  # The list au_vertex is a temporay list!
        y_ave+= float(au_vertex[2])
        z_ave+= float(au_vertex[3])
        y_vertex.append(float(au_vertex[2]))
        z_vertex.append(float(au_vertex[3]))
y_ave=y_ave/5
z_ave=z_ave/5
print(y_vertex)
print(z_vertex)
print(y_ave)
print(z_ave)

r=(((y_ave-y_vertex[4])**2)+((z_ave-z_vertex[4])**2))**(1/float(2))
t=math.acos((y_vertex[4]-y_ave)/r)

# Read data from name2
N_Au = int(linecache.getline(name2, 1))  # The Number of Gold atoms in the system
print(N_Au)

#=========
# FUNCTION
#=========
def point_in_poly(x,y,poly):

	n = len(poly)
	inside = False
	p1x,p1y = poly[-1]
	for i in range(n):
		p2x,p2y = poly[i% n]
		if min(p1y,p2y) < y <= max(p1y,p2y) and x <= max(p1x,p2x):
			if p1y != p2y:
				xints = (y-p1y)*(p2x-p1x)/float(p2y-p1y)+p1x
			if p1x == p2x or x <= xints:
				inside = not inside
		p1x,p1y = p2x,p2y
	return inside

#======================================
# Extract the points inside the polygon
#======================================

file1="ts_"+name_part0+".dat"
file2="pol_"+name_part0+".dat"
file3="inside_"+name_part0+".xyz"

fi =open(name,"r")
fo =open(file1,"w")
fo2=open(file2,'w')
fo3=open(file3,'w')
# fo.write('%6s %6s %6s %6s\n' % ("No.","x","y","z"))
ypol=[]
zpol=[]
x=[]
y=[]
z=[]
x1=[]
y1=[]
z1=[]
pol=[]
pol2=[] # For testing, since tuple is not callable

# determining the vertices of the polygon
#========================================
for i in range(5):
        ypol.append(y_ave + r*math.cos(thi*(i+1)+t))
        zpol.append(z_ave + r*math.sin(thi*(i+1)+t))
# for i in range(5):
        pol.append([]) # empty list as an element of the list pol
        pol2.append([])
for i in range(5):
        pol[i].append(ypol[i])
        pol[i].append(zpol[i])
        pol[i]=tuple(pol[i]) # transfer each element to tuple
        pol2[i]=pol[i]
for i in range(5):
        fo2.write('{0:5.6f} {1:5.6f}\n'.format(pol2[i][0],pol2[i][1]))
print('Written in \'',file2,'\' !')

# Read the coordinates from file and print the result into file
#============================================================== 
for i in range(3,n_step+3):      # To skip the first two lines in the xyz file. so start from 3. 
        lis = linecache.getline(name, i).split() 
        x.append(float(lis[1]))      
        y.append(float(lis[2]))
        z.append(float(lis[3]))

p=0
for i in range(n_step):
	if point_in_poly(y[i],z[i],pol)==True: 
		x1.append(x[i])
		y1.append(y[i])
		z1.append(z[i])
		fo.write('{0:5d} {1:5.6f} {2:5.6f} {3:5.6f}\n'.format(p+1,x1[p],y1[p],z1[p]))
		p=p+1
fo3.write('{0:6d}\n'.format(p))
fo3.write('%6s %6s %6s %6s\n' % ("Atom","x","y","z"))
for i in range(p):
	fo3.write('{0:2s} {1:5.6f} {2:5.6f} {3:5.6f}\n'.format('N',x1[i],y1[i],z1[i]))
print(p,' point inside!\n')
print('Written in \'',file1,'\' !')

fi.close()
fo.close()
fo2.close()
fo3.close()

#===========
# Copy files
#===========
with open(name2) as fi2:
	with open("inside_Au_"+name_part0+".xyz", "w") as fo4:
		fo4.write('{0:6d}\n'.format(p+N_Au))
		fo4.write('%6s %6s %6s %6s\n' % ("Atom","x","y","z"))
		for line in fi2:
			if "AU" in line:
				fo4.write(line)
with open(file3) as fi3:
	with open("inside_Au_"+name_part0+".xyz", "a") as fo4:  # The mode 'a' denote 'append'.
		for _ in range(2):                              # Skip two lines
			next(fi3)
		for line in fi3:
			fo4.write(line)
