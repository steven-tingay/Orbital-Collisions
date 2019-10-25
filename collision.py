#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:19:12 2019

@author: stingay
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import scipy.constants
import random
from astropy import coordinates 
from mpl_toolkits import mplot3d
import os

# minimum characteristic length
mincl=0.115

# total masses
masses=np.array([100,500,1000])

# output flag
out=0

# fraction of total momentum magntidues as threshold for "zero" net momentum
frac=0.01

# number of time steps (of one second) in movie
T=100

# scaling of axis length relative to maximum distance of object
axsca=0.25

# density of objects
#d=1

# minimum percentage of total mass required in debris cloud
permass=0.6

# heights above Earth (km) - now looped over at end of script
alts=np.array([300,800,1300])

# fragment areas
area300=np.array([0.01,0.1])
area800=np.array([0.05,0.5])
area1300=np.array([0.1,1])

# sampling times
times=np.array([10,30,60])

# figure size
fs=5

# initial speed along x axis (m/s)
#v0=np.sqrt(scipy.constants.G*6e24/(h*1000+6.4e6))

# write a shell script to run the breakup model code with the parameters above, find debris mass greater than a given threshold

for totmass in masses:
    g=open("/Users/stingay/.spyder-py3/collision.sh", "w")
    g.write("#!/bin/sh\n")
    g.write("~stingay/.spyder-py3/collision-new > ~stingay/.spyder-py3/collision.dat << EOF\n")
    g.write("%f\n" % mincl)
    g.write("%i\n" % totmass)
    g.write("%i\n" % out)
    g.write("EOF\n")
    g.close()
    
    tmass=0
    while tmass<permass*totmass:
        os.system('chmod +x ~stingay/.spyder-py3/collision.sh')
        os.system('~stingay/.spyder-py3/collision.sh')
    
        # open the file generated and extract mass and deltav from output, create arrays for mass and deltav
    
        f = open("/Users/stingay/.spyder-py3/collision.dat", "r")
    
        line=f.readline()
    
        nmass=[]
        ndeltav=[]
        narea=[]
        while line:
            if "Cross" in line:
                cs=line.split()
                narea.append(float(cs[4]))
                line=f.readline()
                ml=line.split()
                nmass.append(float(ml[2]))
                line=f.readline()
                dl=line.split()
                ndeltav.append(float(dl[3]))
            line=f.readline()
        f.close()
    
        # calculate total mass
    
        tmass=np.sum(nmass)
        
    print(tmass)
    
    # calculate total magnitude of momentum of objects
    
    dmona=[]
    dmone=[]
    
    tmom=np.sum(np.multiply(nmass,ndeltav))
    #print(tmom)
    
    # loop over assignment of directions until net momentum is below threshold
    
    pt=tmom
    while pt>frac*tmom:
        for x in range(len(nmass)):
            dmona.append(math.radians(random.randint(-90,90)))
            dmone.append(math.radians(random.randint(1,361)))
    
        vx,vy,vz=astropy.coordinates.spherical_to_cartesian(ndeltav,dmona,dmone)
    
        px=np.multiply(vx,nmass)
        py=np.multiply(vy,nmass)
        pz=np.multiply(vz,nmass)
    
        pxt=np.sum(px)
        pyt=np.sum(py)
        pzt=np.sum(pz)
        pt=np.sqrt(pxt*pxt+pyt*pyt+pzt*pzt)
        dmona=[]
        dmone=[]
    
    #print(pxt,pyt,pzt,pt)
    
    # write out the velocity components
    
    np.savetxt('/Users/stingay/.spyder-py3/collision.out',(vx,vy,vz))
    #print(vx,vy,vz)
    
    # plot deltav as function of mass
    
    plt.figure()
    plt.scatter(nmass,ndeltav)
    plt.xlabel("mass (kg)")
    plt.ylabel("delta_v (m/s)")
    plt.close()
    
    # plot histograms of components of momentum
    
    plt.figure()
    plt.hist(px,bins=100)
    plt.hist(py,bins=100)
    plt.hist(pz,bins=100)
    plt.xlabel("p (kg m/s)")
    plt.ylabel("#")
    plt.close()
               
    # find maximum distance of object
               
    XYZ=[np.max(np.abs(T*vx)),np.max(np.abs(T*vy)),np.max(np.abs(T*vz))]
    maxl=axsca*np.max(XYZ)/1000
    
    # generate debris cloud and movie frame at each time step 
    
    perdet=[]          
    for h in alts:
        plt.figure()
        plt.hist(nmass,100,histtype='bar',cumulative=True)
        plt.xlabel('Mass (kg)')
        plt.ylabel('Number')
        #plt.show()
        file="massdist{}-{}.pdf".format(totmass,h)
        plt.savefig(file)
        plt.close()
    
        plt.figure()
        plt.hist(narea,100,histtype='bar',cumulative=True)
        plt.xlabel('Area (sq. m)')
        plt.ylabel('Number')
        #plt.show()
        file="areadist{}-{}.pdf".format(totmass,h)
        plt.savefig(file)
        plt.close()
    
        plt.figure()
        plt.hist(ndeltav,100,histtype='bar',cumulative=True)
        plt.xlabel('Delta-v (m/s)')
        plt.ylabel('Number')
        #plt.show()
        file="deltavdist{}-{}.pdf".format(totmass,h)
        plt.savefig(file)
        plt.close()
        
        if h==300:
            area=area300
        if h==800:
            area=area800
        if h==1300:
            area=area1300
            
        per=[]
        time=[]
        for t in range(1,T):
            x=t*vx/1000
            y=t*vy/1000
            z=t*vz/1000
            file="frame{}.png".format(t+1000)
            fig = plt.figure(figsize=(fs,fs))  
            ax = fig.add_subplot(111, projection='3d')  
#            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.set_zlabel('Z (km)')
        
            # the 3D plot
        
            ax.scatter(x, y, z,c=nmass,cmap='summer') 
            ax.set_xlim(-1*maxl,maxl)
            ax.set_ylim(-1*maxl,maxl)
            ax.set_zlim(-1*maxl,maxl)
            fig.savefig(file)
            plt.close()
        
            # the 2D plot of the sky, as seen from the ground, in angular coordinates
        
            ax=[]
            ay=[]
            d=[]
            index=0
            for i in x:
                ax.append(math.degrees(math.atan2(i,h+z[index])))
                index+=1
            index=0
            for j in y:
                ay.append(math.degrees(math.atan2(j,h+z[index])))
                index+=1
            for i in range(0,len(ay)):
                d.append(np.sqrt(ax[i]*ax[i]+ay[i]*ay[i]))
            per.append(100*sum(i>20 for i in d)/len(nmass))
 #           print(totmass,h,t,per[len(per)-1])
            time.append(t)
            plt.figure(figsize=(fs,fs))
            plt.scatter(ax,ay,c=-1*vz/1000,cmap='Spectral')
            plt.gca().set_aspect('equal', adjustable='box')
            circle=plt.Circle((0,0),20,color='b',fill=False)
            plt.gca().add_artist(circle)
            plt.xlim(40,-40)
            plt.ylim(-40,40)
            plt.xlabel('X (deg)')
            plt.ylabel('Y (deg)')
            filexy="xyframe{}-{}.png".format(t+1000,h)
            plt.savefig(filexy)
            if t in times:
                if t==times[2]:
                    plt.colorbar(fraction=0.046,pad=0.04).set_label('Line-of-sight speed m/s')
                filet="frame{}-{}-{}.pdf".format(totmass,h,t)
                # plt.show()
                plt.savefig(filet)
            plt.close()
            if t in times:
                for ar in area:
                    axr=[]
                    ayr=[]
                    vzr=[]
                    plt.figure(figsize=(fs,fs))
                    alist=np.where(narea>ar)[0]
                    perd=100*(len(alist)/len(ax))
                    if t==times[0]:
                        print(totmass,h,ar,perd)
                    for index in alist:
                        axr.append(ax[index])
                        ayr.append(ay[index])
                        vzr.append(vz[index])
                    plt.scatter(axr,ayr,cmap='Spectral')
                    plt.gca().set_aspect('equal', adjustable='box')
                    circle=plt.Circle((0,0),20,color='b',fill=False)
                    plt.gca().add_artist(circle)
                    plt.xlim(40,-40)
                    plt.ylim(-40,40)
                    plt.xlabel('X (deg)')
                    plt.ylabel('Y (deg)')
                    filet="frame{}-{}-{}-{}.pdf".format(totmass,h,t,ar)
                    # plt.show()
                    plt.savefig(filet)
                    plt.close()
        
        plt.figure()
        plt.plot(time,per)
        plt.xlabel("time (s)")
        plt.ylabel("% fragments outside beam FWHM")
        # plt.show()
        fileper="per-out-beam-{}-{}.pdf".format(totmass,h)
        plt.savefig(fileper)
        plt.close()
    
    # generate the movie and clean up directory
        
    #os.system('for f in frame*.png; do /usr/local/bin/convert $f ${f/f/xyf} +append $f; done')
    os.system('/usr/local/bin/convert frame*.png collision.gif')
    os.system('rm -rf ~stingay/.spyder-py3/*frame*.png')


    