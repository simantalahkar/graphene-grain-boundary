import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from shapely.geometry import Point, Polygon
#from copy import deepcopy

######################## MAKING PRISTINE GRAPHENE ATOMIC COORDINATES WITH ARMCHAIR ALONG x, 

d=1.42
dy=np.sqrt(3)*0.5*d
xmax=174.66
ymax=174.66
dy=np.sqrt(3)*0.5*d

x=0.25*d
y=0.5*dy
x=x+d
Gpris=np.array([[x,y]])	
Gprisyarm=np.array([[y,x]]) #rotated crystal with armchair along y

dx=3*d

#place=0	#at point a in lattice when the next horizontal jump is 2d equivalent to place is false/0


############################ Making rotated graphene grains

##FOR 16.1 degrees tilt of each grain bw armchair direction and GB periodic distance is 15.359648433; taking the edge length 7 times periodic distance

repeat=15.3596
edge=4*repeat
r=edge/np.sqrt(2)

difsum=0
tol=0.0000001

#fp=open("my_array.npy","r+")

edge1copy=np.loadtxt("edge1coord150k.txt").copy()
print("hi")
print(edge1copy)
vor = Voronoi(edge1copy)
voronoi_plot_2d(vor)
#plt.figure(figsize=(11,5))
plt.axes().set_aspect('equal')
plt.show()
test=0
displace=3*repeat
#while dif>tol:
while test<100000:
	i=0
	difsum=0
	numofgenupdates=0
	copyshift=np.array([[0,0]])
	while i<len(vor.regions):
		i=i+1
		if i>=len(vor.regions):
			break
		if np.prod(vor.regions[i])==1:
			continue
		counter=0
		m=0
		while m<len(vor.regions[i]):
			if vor.regions[i][m]<0:
				counter=1
			m=m+1
		if counter:
			continue
		verti=np.array([[0,0]])
		for j in vor.regions[i]:
			verti=np.append(verti,[vor.vertices[j]],axis=0)
		verti=np.delete(verti,[0],axis=0)
		counter=0
		m=0
		while m<len(verti):
			if (verti[m][0]>=edge+10)|(verti[m][1]>=edge-7.6798)|(verti[m][0]<=edge-10)|(verti[m][1]<=7.6798):
				counter=1
			m=m+1
		if counter:
			continue
		#print(verti)
		tupe1=[(verti[0,0],verti[0,1]),(verti[1,0],verti[1,1])];
		l=2
		while l<len(verti):
			tupe2=[(verti[l,0],verti[l,1])];
			tupe3=tupe1+tupe2;
			tupe1=tupe3;
			l=l+1
		poly1 = Polygon(tupe1)
		poly2 = Polygon(poly1.envelope)
		k=-1
		while k<len(edge1copy):
			k=k+1
			if k>=len(edge1copy):
				break
			#counter=0
			point = Point(edge1copy[k,0],edge1copy[k,1])
			if poly2.contains(point): #update the point to centroid of the vertices(polygon) and store difference
				#print("it is inside")
				temp=edge1copy[k].copy()
				centr=poly2.centroid
				c=np.array([centr.x,centr.y])
				edge1copy[k]=c.copy()
				diff=np.linalg.norm(edge1copy[k]-temp)
				difsum=difsum+diff
				numofgenupdates=numofgenupdates+1
				break
			
			
		#print(i)
	m=0
	indexdel=[]
	while m<len(edge1copy):
			if (edge1copy[m][1]>=edge-7.6798)|(edge1copy[m][1]<=7.6798):
				indexdel.append(m)				
			m=m+1
	edge1copy=np.delete(edge1copy,indexdel,axis=0)
	k=-1
	while k<len(edge1copy):
		k=k+1
		if k>=len(edge1copy):
			break
		if (edge1copy[k][1]<=edge-7.6798)&(edge1copy[k][1]>=edge-repeat):
					copyshift=np.append(copyshift,[[edge1copy[k][0],edge1copy[k][1]-displace]],axis=0)
		elif (edge1copy[k][1]<=repeat)&(edge1copy[k][1]>=7.6798):
					copyshift=np.append(copyshift,[[edge1copy[k][0],edge1copy[k][1]+displace]],axis=0)
	copyshift=np.delete(copyshift,[0],axis=0)
	edge1copy=np.append(edge1copy,copyshift,axis=0)
	vor = Voronoi(edge1copy)
	#voronoi_plot_2d(vor)
	#if difsum<tol:
	#	break
	test=test+1
	print(test)
voronoi_plot_2d(vor)
#plt.figure(figsize=(11,5))
plt.axes().set_aspect('equal')	
plt.show()
print("the number of generator updates={}, and total difference={}ength of edge1copy={}".format(numofgenupdates,difsum,len(edge1copy)))
#f1=open("edge1arr200k","w")
fa = open("edge1coord250k.txt","w")
#np.save('my_array200k', f1)
np.savetxt(fa, edge1copy)


			

	




























