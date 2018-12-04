import numpy as np
import glob
import re
import os
from math import pi

# def is_number(s):
	# try:
		# float(s)
		# return True
	# except ValueError:
		# pass
 
	# try:
		# import unicodedata
		# unicodedata.numeric(s)
		# return True
	# except (TypeError, ValueError):
		# pass
 
	# return False

	
def is_number(s):
	try:
		p = float(s)
		return True, p
	except ValueError:
		pass
 
	try:
		import unicodedata
		p =unicodedata.numeric(s)
		return True, p
	except (TypeError, ValueError):
		pass 
 
	return False, 0
	
	
buffer = []

for file in glob.glob('./leb/*'):
	with open(file, "r") as f:
		ispolar=0
		iscstw=0
		buffer = []
		for line in f:
			group = re.split("[, ]+", line.strip())
			# print group
			if(len(group)<2):
				continue
			# First check if polar or cartesian
			bool1,val1 = is_number(group[0])
			bool2,val2 = is_number(group[1])
			if(ispolar==0):
				if(bool1 and bool2) :
					if(val1 > 1 or val2 > 1) :
						print "Polar Data found"
						ispolar=1
						# print val1+val2
			
			if(len(group)==3):
				bool3,val3 = is_number(group[2])
				if(bool1 and bool2 and bool3):
					if(ispolar==1):
						x=np.cos(val1*pi/180)*np.sin(val2*pi/180)
						if(abs(x)<=10e-12): x=0
						y=np.sin(val1*pi/180)*np.sin(val2*pi/180)
						if(abs(y)<=10e-12): y=0
						z=np.cos(val2*pi/180)
						if(abs(z)<=10e-12): z=0
						buffer.append( str(x) +' ' + str(y) +' ' + str(z) +' ' +group[2])
					else:
						buffer.append(group[0] + ' ' + group[1] + ' ' + group[2])
						if(iscstw==0 and ispolar==0): 
							iscstw=1
							print file, "Constant Weight found"
				
			
		# print buffer

		# Writing Operation
		if not os.path.exists("./resh"):
			os.makedirs("./resh")
		
		with open("./resh/" + file.split('\\')[-1] +"_resh", "w+") as fp :
			fp.write("$ConstantWeight\n")
			if(iscstw==0): fp.write("0\n")
			else : fp.write("1\n")
			fp.write("$Polar\n")
			if(ispolar==0): fp.write("0\n")
			else : fp.write("1\n")
			fp.write("$Nodes\n")
			fp.write( str(len(buffer)) + "\n")
			for item in buffer:
				fp.write("%s\n" % item)
			fp.write("$EndNodes")

		# thefile.write("%s\n" % item)
			# print buffer
			# m = re.findall(r"[-+]?(\d+\.\d+[eE][-+]\d+)", line)
			
			# print m
		
			
			
			
			
			
	# coord = np.loadtxt(file)
	# print coord
    # coord = np.reshape(coord,(len(coord)/3,3))
    # np.savetxt(file,coord)
