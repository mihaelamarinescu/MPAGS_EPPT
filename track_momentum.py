# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 20:55:53 2021

@author: Mihaela
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import uproot

def fit_track(z_track, x_track):
    
    mean_z = np.mean(z_track)
    mean_x = np.mean(x_track)
    # Total number of values
    n = len(z_track)
    # Using the formula to calculate 'm' and 'c'
    numer = 0
    denom = 0
    chisq = 0
    for i in range(n):
        numer += (z_track[i] - mean_z) * (x_track[i] - mean_x)
        denom += (z_track[i] - mean_z) ** 2
    m = numer / denom
    c = mean_x - (m * mean_z)
    
    return m, c


def chisqfunc(x_tracks, z_tracks, c, m):
    err_x = 100*1e-3 #err in mm
    x = c + m*z_tracks
    chisq = np.sum(((x_tracks - x)/err_x)**2)
    return chisq

def run_tracks(filename, treename):
    
    file = uproot.open(filename+'.root')
    tree = file[treename]
    branches = tree.arrays()
    momenta = []
    
    for entry in range(1000):
        print("Event:", entry)
        dc1_xpos = branches['Dc1HitsVector_x'][entry]
        dc1_ypos = branches['Dc1HitsVector_y'][entry]
        dc1_zpos = branches['Dc1HitsVector_z'][entry]
        dc2_xpos = branches['Dc2HitsVector_x'][entry]
        dc2_ypos = branches['Dc2HitsVector_y'][entry]
        dc2_zpos = branches['Dc2HitsVector_z'][entry]

        # print(dc1_zpos)
        # print(dc1_xpos)
        # print(dc2_zpos)
        # print(dc2_xpos)
        #xcoordinates = np.zeros(len(dc1_xpos)+len(dc2_xpos))
        #zcoordinates = np.zeros(len(dc1_zpos)+len(dc2_zpos))
        xcoordinates1 = []
        zcoordinates1 = []
        
        for i in range(0, len(dc1_zpos)):
            zcoordinates1.append((-6.25 + 0.5*(dc1_zpos[i]))*1e3) #position in mm
            xcoordinates1.append(dc1_xpos[i])
            if dc1_zpos[i] == 4: break # stop at first 4

        xcoordinates2 = []
        zcoordinates2 = []

        for i in range(0, len(dc2_zpos)):
            zcoordinates2.append((2.25 + (dc2_zpos[i])*0.5)*1e3) #position in mm
            xcoordinates2.append(dc2_xpos[i])
            if dc2_zpos[i] == 4: break # stop at first 4

        xcoordinates = np.append(xcoordinates1, xcoordinates2)
        zcoordinates = np.append(zcoordinates1, zcoordinates2)
        
       
        
        # print(xcoordinates1)
        # print(zcoordinates1)
        # print(xcoordinates2)
        # print(zcoordinates2)
        # print(len(xcoordinates1))
        # print(len(zcoordinates1))
        # print(len(xcoordinates2))
        # print(len(zcoordinates2))
        if len(xcoordinates1) < 1 or len(zcoordinates1) < 1 or len(xcoordinates2)< 1 or len(zcoordinates2) < 1: continue
    
        # Initial guess.
        # x0 = np.array([0.0, 0.0, 0.0])
        
        # popt, pcov = opt.curve_fit(chisqfunc, zcoordinates2, xcoordinates2, x0)
        # print(popt)
        
        # coefficients for track 1
        m1, c1 = fit_track(zcoordinates1, xcoordinates1)
         # coefficients for track 2
        m2, c2 = fit_track(zcoordinates2, xcoordinates2)
            
         
        # Printing coefficients
        #print("Coefficients")
        #print(m, c)
 
        # if entry == 78:
        #     print(zcoordinates2)
        #     print(xcoordinates2)
        #     z_fake1 = np.linspace(min(zcoordinates1), zcoordinates2[1], 100)
        #     z_fake2 = np.linspace(zcoordinates1[len(zcoordinates1)-2], max(zcoordinates2), 200)
        #     plt.plot(zcoordinates1, xcoordinates1, 'c+', label='track1')
        #     plt.plot(z_fake1, z_fake1*m1 + c1, 'k-', label = 'fit1')
        #     plt.plot(zcoordinates2, xcoordinates2, 'b+', label='track2')
        #     plt.plot(z_fake2, z_fake2*m2 + c2, 'r-', label = 'fit2')
        #     plt.xlabel('z [mm]')
        #     plt.ylabel('x [mm]')
        #     plt.legend()
        #     plt.show()
        
        z1 = np.sqrt(1/(m1**2)-1)
        x1 = z1*m1 +c1
        #find eq of line perpendicular to track1
        csq = x1+ (1/m1)*z1
        z2 = (csq - c2)/(m2+ 1/m1)
        x2 = z2*m2 + c2
        #print(x1, x2)
        #diff_x = np.sqrt((z1-z2)**2 + (x1-x2)**2)*1e-3
        diff_x = abs((m1*z1+c1)-(m2*z2+c2))*1e-3
        L = 1
        delta_theta = np.arcsin(m1) - np.arcsin(m2)
        momentum = 0.3*np.sqrt(1**2 + diff_x**2)/(2*np.sin(delta_theta/2)) #GeV
        #print(momentum)
        if momentum < 500 and momentum >0: momenta.append(momentum)
    
    print(momenta)
    plt.hist(momenta, bins = 100, label = "Entries: "+str(len(momenta)))
    plt.xlabel(r'pT [GeV]')
    plt.ylabel('Entries')
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    
     run_tracks('case1', 'B5')