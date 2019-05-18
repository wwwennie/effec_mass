#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: Wennie Wang, wwwennie

 helper module for calculating the effective band mass 

 implemented here horiztonal shift for band extrema not 
  at Gamma or high symmetry point

 assumes input file of the form generated from qeigenvals.pl
 i.e., 
   # k_x,y,z[reciprocal] k_x,y,z[2pi/scale] dist[2pi/scale][1/Ang][1/Bohr]
   followed by bands in column format

 Fit to hyperbolic effective mass

     hbar^2 k^2/(2m*) = E (1+ gamma*E)
     alpha k^2 = E(1 + gamma*E)

 where gamma is a fitting parameter and alpha = hbar^2/(2m*)
 Solving for E yields a quadratic equation;
----------------------------------------------
"""

from common import *
from scipy.optimize import curve_fit
from tabulate import tabulate

def run_effmass(filename,bands,init,step,bounds=([100,10,10,100],[0.,0.,-10.,-20.]),iscbm=True,plot=True):
   """ Main function for finding effective mass using hyperbolic fit 

       filename (string): name of band eigenval data, formatted ala above
       bands (list, integers): actual band indices of desired fits
       init, step (integer): start and increment of data truncation for finding
                             best fit to hyperbolic curve, valid near Gamma
    
       iscbm = defaults to finding electron effective mass
       bounds: list of bounds for alpha, gamma, xshift, c; 
              e.g., (pbnds,nbnds) -> ([10,10,10,20],[0,0,0,0])
       plot (boolean): whether to make plot of fits
   """

   print("----> %s" % filename)
   for band in bands:
      kpt_raw, band_raw = read_bands(filename,band)
      label, data_pair = trunc_data(kpt_raw,band_raw,init,step)

      # for each truncation specified by init and step
      popts = []
      for pair in data_pair:
         subkpt = pair[0]
         subband = pair[1]
         if iscbm:
            popt,pcov = fit_curve(energymodel_p,subkpt,subband,bounds)
         else:
            popt,pcov = fit_curve(energymodel_m,subkpt,subband,bounds)
         popts.append(popt)       
      
      table, header = pretty_table(label,popts)
      print("=================== Band %d =====================" % band)
      print("alpha = hbar^2/(2m*), m* in units of free electron mass m0")
      print("gamma = hyperbolicity, xshift = horiz fit param, c = vert fit param, CBM")
      print(tabulate(table,headers=header))
  
      if plot:
         if iscbm:
            plot_fits(data_pair, energymodel_p,popts)
         else:
            plot_fits(data_pair, energymodel_m,popts)    

#======================== Helper functions ============================

def energymodel_p(k, alpha, gamma, xshift,c):
   """ Fit bandstructure along k kpt path to alpha and gamma 
       additional c parameter for vertical shifts
            and xshift for horizontal shifts 
      
       positive solution to quadratic in E
       alpha = hbar^2/(2m*)
   """ 

   energy = (-1 + np.sqrt(1 + 4*alpha*gamma*np.power(k-xshift,2)))/(2*gamma) + c
   return energy

def energymodel_m(k, alpha, gamma, xshift, c):
   """ Fit bandstructure along k kpt path to alpha and gamma 
       additional c parameter for vertical shifts 
           and xshift for vertical shifts
      
       negative solution to quadratic in E
       alpha = hbar^2/(2m*)
   """ 

   energy = (1 - np.sqrt(1 + 4*alpha*gamma*np.power(k-xshift,2)))/(2*gamma) + c
   return energy

def fit_curve(func,kpath,band,bounds):
   """ Wrapper for curve fitting via scipy 
       w/ assumed bounds on alpha, gamma, xshift,c
   """

   pbnds = bounds[0]
   nbnds = bounds[1]
   popt, pcov = curve_fit(func,kpath,band,bounds=(nbnds,pbnds))

   return popt, pcov

def trunc_data(kpath,band,init,step):
    """ Wrapper for truncating data in steps
        returns list of lists containing [kpath, band] in varying lengths of data

        init: integer, at what index to start the truncation
        step: integer, increment amount of data to include
    """
    data_pair = []
    index = np.arange(init,len(kpath),step)

    for i in index:
       data_pair.append([kpath[:i],band[:i]])

    return index, data_pair

def plot_fits(data_pair,func,popts):
   """ Wrapper for plotting all the fits of various truncations in a bunch of subplots 

       data_pair: original data fitted against [[kpts1,band1],[kpts2,band2], ...]
       func; function used to fit
       popts: list of fit parameters 
   """

   assert len(data_pair) == len(popts)
 
   # set up matrix of subplots for each truncation
   ncol = 2
   nrow = len(popts)//ncol + len(popts)%ncol
   print("Creating %d x %d subplot" % (nrow,ncol))

   fig = plt.figure(1,figsize=(nrow*5,ncol*50))
   plot_counter = 1
  
   # largest data set
   data_full = data_pair[-1]
   kpt_full, band_full = data_full[0], data_full[1]   
   xmin,xmax = np.min(kpt_full), np.max(kpt_full)*1.1

   # plot raw data with fit
   for pair in data_pair:
      popt = popts[plot_counter-1] # double-duty index 
      ax = fig.add_subplot(nrow,ncol,plot_counter)

      ax.plot(kpt_full,func(kpt_full,*popt),'b') # fit
      ax.plot(pair[0], pair[1],'ko')   # ref

      ax.set_xlim(xmin,xmax)
      ax.set_xlabel("k [1/Ang]",fontsize=40)
      ax.set_ylabel("Energy [eV]",fontsize=40)
      ax.set_title("%d data pts" % len(pair[0]),fontsize=50)
   
      plot_counter += 1
    
   plt.tight_layout()
   plt.show()

   return popts
   
   
def pretty_table(label, popts):
   """ Make a pretty table of the fits and compute equivalent effective mass """

   # extract parameters
   alpha, gamma, xshift, c = [list(a) for a in zip(*popts)]
   #header = ["# data points", r"$\alpha= \frac{\hbar^2}{2m^*}$ [eV*Ang$^2$]",
   #         r"m$^*$ [$\frac{m^*}{m_0}$",
   #         r"$\gamma$ (hyperbolic fit) [1/eV]",
   #         r"c (fitting) [eV]"]
   header = ["# pts", "alpha [eV*Ang2]", "m* [m0]", "gamma [1/eV]", "xshift [1/Ang]","c [eV]"]
  
   # calculate effective mass in terms of free electron mass
   pre = np.power(hbar,2)/(2*np.power(ang2m,2)*ev2joule*m_e)
   invalpha = np.power(alpha, -1.0)
   effecmass = pre*invalpha
  
   # reshuffle parameters for pretty table printing 
   table = zip(label,alpha,effecmass,gamma,xshift,c)       

   return table, header

def read_bands(filename,bandindex):
   """ Wrapper for reading in input file
       Assumes file format outputte accoridng to qe-eigenvals.sh 
     
      Correct data offset is included in; index bands by actual number in calculation
   """
   
   raw_dat = np.loadtxt(filename, delimiter=" ",skiprows=1)
   kpt_raw = raw_dat[:,8] # 9th column, 1/Ang
   band_raw = raw_dat[:,bandindex+9] # offset by 10 columns

   return kpt_raw,band_raw

        
