import numpy as np
import glob
from numpy import *
import random
import re
import os
from math import pi
# from matplotlib import *
import matplotlib.pyplot as plt

from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib.text import Annotation
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mpl_toolkits.mplot3d as a3




vnode = [ 
1.0, 0, 0, 
-1.0, 0, 0, 
0, 1.0, 0, 
0, -1.0, 0, 
0, 0, 1.0, 
0, 0, -1.0, 
0.57735026919, 0.57735026919, 0.57735026919, 
0.57735026919, 0.57735026919, -0.57735026919, 
0.57735026919, -0.57735026919, 0.57735026919,
0.57735026919, -0.57735026919, -0.57735026919, 
-0.57735026919, 0.57735026919, 0.57735026919, 
-0.57735026919, 0.57735026919, -0.57735026919, 
-0.57735026919, -0.57735026919, 0.57735026919, 
-0.57735026919, -0.57735026919, -0.57735026919
] 

vnormal = [ 1.0, 0, 0, 
-1.0, 0, 0, 
0, 1.0, 0, 
0, -1.0, 0, 
0, 0, 1.0, 
0, 0, -1.0 ]

Nd = len(vnode)/3
Nn = len(vnormal)/3


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class Annotation3D(Annotation):
    '''Annotate the point xyz with text s'''

    def __init__(self, s, xyz, *args, **kwargs):
        Annotation.__init__(self,s, xy=(0,0), *args, **kwargs)
        self._verts3d = xyz        

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.xy=(xs,ys)
        Annotation.draw(self, renderer)
		
def annotate3D(ax, s, *args, **kwargs):
    '''add anotation text s to to Axes3d ax'''

    tag = Annotation3D(s, *args, **kwargs)
    ax.add_artist(tag)	
		

def SpherePlot(pt,nn) :
	#Generate surface of the sphere
	r = 1
	phi, theta = np.mgrid[0.0:np.pi:50j, 0.0:2.0*np.pi:50j]
	x = r*np.sin(phi)*np.cos(theta)
	y = r*np.sin(phi)*np.sin(theta)
	z = r*np.cos(phi)

	#Set colours render and labels
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_xlim([-1,1])
	ax.set_ylim([-1,1])
	ax.set_zlim([-1,1])
	ax.set_aspect("equal")
	ax.set_xlabel('x_values')
	ax.set_ylabel('y_values')
	ax.set_zlabel('z_values')
	ax.scatter(0,0,0,marker='s', c = 'black', s = 90)
	#Draw sphere
	# ax.plot_surface(	x, y, z,  rstride=1, cstride=1, color='b', alpha=0.2, linewidth=0.2)
	
	#Draw all directions and labels on unit sphere
	# for i in range(Nd):
		# if i < 6 :
			# cc = 'b'
		# else : 
			# cc = 'g'
		# ax.scatter(vnode[3*i+0], vnode[3*i+1], vnode[3*i+2],marker='o', c = cc, s = 64)
		# annotate3D(ax, s=str(i), xyz=(vnode[3*i+0], vnode[3*i+1], vnode[3*i+2]), fontsize=10, xytext=(-3,3), textcoords='offset points', ha='right',va='bottom')    
		# a = Arrow3D([0, vnode[3*i+0]], [0, vnode[3*i+1]], [0, vnode[3*i+2]], mutation_scale=20, lw=0.2, arrowstyle="-|>", color="r")
		# ax.add_artist(a)
	# plt.tight_layout()
	# plt.show()
	# return
	#Draw the plane 
	point  = np.array([0, 0, 0])
	normal = np.array([0, 0, 1])
	normal = pt
	dd = -point.dot(normal)
	xx, yy = np.mgrid[-1:1:200j, -1:1:200j]
	zz = (-normal[0] * xx - normal[1] * yy - dd) * 1. /normal[2]
	ax.plot_surface(xx, yy, zz,alpha=0.1, lw=0.2)

		
	#Work on the 3 nearest NN create segment between
	seg = []
	vec1 = [ nn[1][0] - nn[0][0] , nn[1][1] - nn[0][1], nn[1][2] - nn[0][2] ]
	vec2 = [ nn[2][0] - nn[1][0] , nn[2][1] - nn[1][1], nn[2][2] - nn[1][2] ]
	norm = np.cross(vec1,vec2)
	d =  np.dot([nn[0][0] , nn[0][1], nn[0][2] ],norm)/( np.dot([pt[0],pt[1],pt[2]],norm))
	print d
	# a = Arrow3D([0, d*pt[0]], [0, d*pt[1]], [0,d* pt[2]], mutation_scale=20, lw=1, color="black")
	edge_col = Line3DCollection( [( (0,0,0), (d*pt[0],d*pt[1], d*pt[2]) )], lw=2, alpha=0.4,color='black', linestyle='--' )
	ax.add_collection3d(edge_col)
	ax.scatter(d*pt[0],d*pt[1], d*pt[2],marker='x', c = 'k', s = 64)
	# ax.add_artist(a)
	# print norm
	trigl = [ [nn[0][0], nn[0][1], nn[0][2]] , [nn[1][0], nn[1][1], nn[1][2]], [nn[2][0], nn[2][1], nn[2][2]] ]
	# print trigl

	seg.append( ((nn[0][0], nn[0][1], nn[0][2]), (nn[1][0], nn[1][1], nn[1][2]) ))
	seg.append( ((nn[1][0], nn[1][1], nn[1][2]), (nn[2][0], nn[2][1], nn[2][2]) ))
	seg.append( ((nn[2][0], nn[2][1], nn[2][2]), (nn[0][0], nn[0][1], nn[0][2]) ))
	
	seg.append( ((nn[2][0], nn[2][1], nn[2][2]), (nn[0][0], nn[0][1], nn[0][2]) ))
	seg.append( ((nn[2][0], nn[2][1], nn[2][2]), (nn[0][0], nn[0][1], nn[0][2]) ))
	seg.append( ((nn[2][0], nn[2][1], nn[2][2]), (nn[0][0], nn[0][1], nn[0][2]) ))
	seg.append( ((nn[2][0], nn[2][1], nn[2][2]), (nn[0][0], nn[0][1], nn[0][2]) ))
	
	#trouver point NN de la normal a la surface
	#utiliser tableau pour connaitre direct les 7 NN
	#connaitre les relation de corespondance lors que direction K est 
	
	
	
	edge_col = Line3DCollection(seg, lw=1,color='black')
	ax.add_collection3d(edge_col)
	tri = a3.art3d.Poly3DCollection([trigl])
	tri.set_color('red')
	tri.set_alpha(0.3)
	tri.set_edgecolor('k')
	ax.add_collection3d(tri)
	#TODO : draw arc
	
	#Draw elligible directions associated to the 3 NN
	for i in range(9):
		# seg.append( ((nn[i][0], nn[i][1], nn[i][2]), (nn[i+1][0], nn[i+1][1], nn[i+1][2]) ))
		if i == 0 :
			cc = 'green'
		if i == 1:
			cc = 'blue'
		if i ==2 :
			cc = 'red'
		ax.scatter(nn[i][0],nn[i][1], nn[i][2],marker='o', c = cc, s = 64)
		annotate3D(ax, s=str(nn[i][4]), xyz=(nn[i][0],nn[i][1], nn[i][2]), fontsize=10, xytext=(-3,3), textcoords='offset points', ha='right',va='bottom')    
		a = Arrow3D([0, nn[i][0]], [0, nn[i][1]], [0, nn[i][2]], mutation_scale=20, lw=0.7, arrowstyle="-|>", color=cc)
		ax.add_artist(a)

	#Draw the new random vector	
	cc = 'black'
	# ax.scatter(pt[0], pt[1], pt[2],marker='x', c = cc, s = 64)
	annotate3D(ax, s='Normal', xyz=(pt[0], pt[1], pt[2]), fontsize=10, xytext=(-3,3), textcoords='offset points', ha='right',va='bottom')    
	a = Arrow3D([d*pt[0], pt[0]], [d*pt[1], pt[1]], [d*pt[2], pt[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="black")
	ax.add_artist(a)
	plt.tight_layout()
	plt.show()
	return

def RandomPickSphere(n):
	v = [random.gauss(0, 1) for i in range(0, n)]
	inv_len = 1.0 / math.sqrt(sum(coord * coord for coord in v))
	return [coord * inv_len for coord in v]

def Distance(a,b):
	return np.sqrt( (b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) +  (b[2]-a[2])*(b[2]-a[2]))

def DistanceS(a,b):
	return np.arccos(a[0]*b[0] + a[1]*b[1] + a[2]*b[2] )

def FindNN(a):
	dist = 1e5
	arr = np.zeros((14,6))
	print "Random vect : ", a
	for i in range(0,Nd):
		new = Distance(a,vnode[3*i:3*(i+1)])
		arr[i][0] = vnode[3*i+0]
		arr[i][1] = vnode[3*i+1]
		arr[i][2] = vnode[3*i+2]
		arr[i][3] = new
		arr[i][4] = i
		arr[i][5] = vnode[3*i+0]*a[0] + a[1]*vnode[3*i+1] + a[2]*vnode[3*i+2]
		# new = DistanceS(a,vnode[3*i:3*(i+1)])
		# print "Vector : ", vnode[3*i:3*(i+1)], "Distance :", new
		if new < dist:
			dist = new
			idx=i
	return idx,dist,vnode[3*i:3*(i+1)], arr[np.argsort(arr[:, 3])]

		
def ComputeSpec(vnormal, prev_idx):
	vpost = np.zeros(3)
	vpost2 = np.zeros(3)
	vn = vnormal[0]*vnode[3*prev_idx + 0] + vnormal[1]*vnode[3*prev_idx + 1] + vnormal[2]*vnode[3*prev_idx + 2]
	#Collision 
	if (vn > 0) :
		for i in range(0,3):
			vpost[i] = vnode[3*prev_idx + i] - 2*vn*vnormal[i]
			
		#Check  reversed bounceback
		vn2 = -vpost[0]*vnode[3*prev_idx + 0] -vpost[1]*vnode[3*prev_idx + 1] -vpost[2]*vnode[3*prev_idx + 2]
		for i in range(0,3) :
			vpost2[i] = -vpost[i] - 2*vn*vnormal[i]
			assert( -vpost2[i] - vnode[3*prev_idx + i] == 0 )
		
		# print "vnorm = ", vnormal, "\nvprev = ", vnode[3*prev_idx:3*(prev_idx+1)], "\nvpost  = ", vpost
	#No Collision
	else :
		vpost[0] = vnode[3*prev_idx + 0] 
		vpost[1] = vnode[3*prev_idx + 1] 
		vpost[2] = vnode[3*prev_idx + 2] 
		# print "NO Colision case : "
	return vpost
	
	# //6 case
			# // int hsphere6 = {0, 2, 4, 6, 7, 8, 10 };
			# //7 case = 
			# // int hsphere7 = { 0, 2, 5, 6, 7, 9, 11};
			# //8 case = 
			# // int hsphere8 = { 0, 3, 4, 6, 8, 9, 12};
def CheckPostv():
	for i in range(0,Nd) :
		for k in range(0,Nn) :
			vpost =	ComputeSpec(vnormal[3*k:3*(k+1)],i)
			# print  "vnorm = ", vnormal[3*k:3*(k+1)], "vprev = ", vnode[3*i:3*(i+1)], "vpost  = ", vpost   


def TestFindNN():
	LL = []
	LL2 = list(range(0,14))
	for i in range(0,Nd):
		LL.append(FindNN(vnode[3*i:3*(i+1)]))
		assert(LL[i] - LL2[i] == 0)

# TestFindNN()
# a = RandomPickSphere(3)
# a =[1,0,0]
k=1
a= vnode[3*k:3*(k+1)]
idx, dist, vec, nn, = FindNN(a)
# print "Random : ", a
print nn
nnid = np.zeros(Nd)
for i in range(Nd):
	nnid[i] = nn[i][4]

print np.sort(nnid[0:9])
SpherePlot(a,nn)
# exit()

# ComputeSpec(vnormal[0:3],0)
# CheckPostv()
