import os
import sys
import numpy as np

if(len(sys.argv)<3):
  print('provide 2 files')
  sys.exit()
  

file1=sys.argv[1]
file2=sys.argv[2]

print("#file1="+str(file1))
print("#file2="+str(file2))

if(not os.path.isfile(file1)):
  print("Error: file "+str(file1)+" does not exist")
  sys.exit()
print("#file1 exists")
  
if(not os.path.isfile(file2)):
  print("Error: file "+str(file2)+" does not exist")
  sys.exit()
print("#file2 exists")

array1 = np.loadtxt(file1,ndmin=2)
array2 = np.loadtxt(file2,ndmin=2)

num_col1 = np.size(array1[0,:]) 
num_row1 = np.size(array1[:,0])  

num_col2 = np.size(array2[0,:]) 
num_row2 = np.size(array2[:,0])  

print("#num_row1,num_col1="+str(num_row1)+","+str(num_col1))
print("#num_row2,num_col2="+str(num_row2)+","+str(num_col2))

if(num_col1!=num_col2):
  print("#num_col1 and num_col2 differ")
  sys.exit()

#check mesh and compute diff when same mesh
if(num_row1==num_row2):
  mesh_error = max(abs(array1[:,0]-array2[:,0]))
  print("#mesh_error="+str(mesh_error))
  if(abs(mesh_error)>1.e-13):
    print("#mesh_error too big")
    sys.exit()
  val=np.zeros((num_col1-1,2))
  err=np.zeros((num_row1,num_col1))
  err=array1-array2   
  for i in range(num_col1-1):
    val[i,0] = max(abs(array1[:,i+1]-array2[:,i+1])) 
    val[i,1] = sum(abs(array1[:,i+1]-array2[:,i+1]))/num_row1 
  if(num_col1==2):
    print "#%1.20g %1.20g" %(val[0,0],val[0,1])
  if(num_col1==3):
    print "#%1.20g %1.20g %1.20g %1.20g" %(val[0,0],val[0,1],val[1,0],val[1,1])
  np.savetxt('err.dat',err)

#check mesh and compute diff when same mesh
if(num_row1!=num_row2):
  factor = (num_row2-1)/(num_row1-1)
  inv_factor = (num_row1-1)/(num_row2-1)
  print("#factor="+str(factor))
  print("#inv_factor="+str(inv_factor))
  if(factor==0):
    if(inv_factor*(num_row2-1)!=num_row1-1):
      print("#bad compatibility")
      sys.exit()
  if(inv_factor==0):
    if(factor*(num_row1-1)!=num_row2-1):
      print("#bad compatibility")
      sys.exit()
  if((factor!=0) and (inv_factor!=0) ):
    print("factor or inv_factor should be 0")
    sys.exit()
  num_row = min(num_row1,num_row2)   
  sol1 = np.zeros((num_row,num_col1))
  sol2 = np.zeros((num_row,num_col1))
  if(factor==0):
    for j in range(num_col1):
      sol1[:,j] = array2[:,j]
      sol2[:,j] = array1[0:num_row1:inv_factor,j]
  if(inv_factor==0):
    for j in range(num_col1):
      sol1[:,j] = array1[:,j]
      sol2[:,j] = array2[0:num_row2:inv_factor,j]
  
  
  #print(sol1[:,0])
  #print(sol2[:,0])
  #sys.exit()
  mesh_error = max(abs(sol1[:,0]-sol2[:,0]))
  print("#mesh_error="+str(mesh_error))
  if(abs(mesh_error)>1.e-13):
    print("#mesh_error too big")
    sys.exit()
  val=np.zeros((num_col1-1,2))
  err=np.zeros((num_row,num_col1))
  err=sol1-sol2   
  for i in range(num_col1-1):
    val[i,0] = max(abs(sol1[:,i+1]-sol2[:,i+1])) 
    val[i,1] = sum(abs(sol1[:,i+1]-sol2[:,i+1]))/num_row1 
  if(num_col1==2):
    print "#%1.20g %1.20g" %(val[0,0],val[0,1])
  if(num_col1==3):
    print "#%1.20g %1.20g %1.20g %1.20g" %(val[0,0],val[0,1],val[1,0],val[1,1])
  np.savetxt('err.dat',err)
    

