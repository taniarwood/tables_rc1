import numpy as np
import scipy as sp
import pickle as pckl
from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class MCEqFluxSpline(object):
    
####################################################
    #LOAD THE DATA
#################################################
    

    def LoadData(self, energylist_file, coszenlist_file, fluxfile2D):
        '''
        takes input of three files, for example energylist_file = "egrid_list.txt", coszenlist_file="cos_zenith_list.txt", fluxfile2D="numu_from_kaions_flux_table.txt"
        '''
        egrid = np.loadtxt(energylist_file)
        #print egrid
        coszenlist = np.loadtxt(coszenlist_file)
        fluxfile2D = np.genfromtxt(fluxfile2D, delimiter=',')
        
        SplinedFlux= RectBivariateSpline(egrid,coszenlist, fluxfile2D , kx=2, ky=2)

        return SplinedFlux
    


    def EvaluateSpline(self, spline_name = None, energy = 0, zenith = -1):
        return self.spline_dict[spline_name].ev(energy, zenith)

    def Plot2Dflux(self, spline):
        x = np.linspace(0., 4, 100)
        y = np.linspace(-1, 1, 100)
        xx, yy = np.meshgrid(x, y)
        xx = 10**xx
        
        values = self.EvaluateSpline(spline, xx, yy)
        plt.pcolor(x, y, values)
        plt.colorbar()
        plt.show()


################################################
    #Make the ID splines (one for numu, one for nue) and flat fluxes
####################################################
    ########## NU MU #################################
    def make_fmu(self,energy_use):
    
       # data = pckl.load(open('/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/tw_neutrino_flux_wsumsDPM_GH.pckl'))
        
	##IDEALLY ONLY DO THIS ONCE

        eedges = np.linspace(1., 1000, 999)
	

        ecenters = (eedges[1:] + eedges[:-1])/2.
        ecenters_new = (eedges[1:] + eedges[:-1])/2.        

	esteps = eedges[1:]-eedges[:-1]

        zedges = np.linspace(-1,1, 101)
        zcenters = (zedges[1:] + zedges[:-1])/2.

    	sum_numu_from_p = np.zeros_like(ecenters)
    	sum_numu_from_k = np.zeros_like(ecenters)
   	sum_numubar_from_p = np.zeros_like(ecenters)
    	sum_numubar_from_k = np.zeros_like(ecenters)
   # 	sum_nue  = np.zeros_like(ecenters)
   # 	sum_nuebar  = np.zeros_like(ecenters)

    	for eii, energyi  in enumerate(ecenters):
        	for zii , czenith  in enumerate(zcenters):
			sum_numu_from_p[eii] +=(self.spline_dict['numu_from_pion'].ev(energyi, czenith))
            		sum_numu_from_p[eii] +=(self.spline_dict['numu_from_pion'].ev(energyi, czenith))
            		sum_numu_from_k[eii] +=(self.spline_dict['numu_from_k'].ev(energyi, czenith))
            		sum_numubar_from_p[eii] +=(self.spline_dict['antinum_from_pion'].ev(energyi, czenith))
            		sum_numubar_from_k[eii] +=(self.spline_dict['antinum_from_k'].ev(energyi, czenith))
            	#	sum_nue[eii]  +=(self.spline_dict['nue'].ev(energyi, czenith)/energyi**3)
            	#	sum_nuebar[eii] +=(self.spline_dict['antinue'].ev(energyi, czenith)/energyi**3)


        numu_tot_use_y =  (sum_numu_from_p + sum_numubar_from_p + sum_numu_from_k + sum_numubar_from_k)
        
        fmu = interp1d(ecenters, numu_tot_use_y, kind='cubic' )
        #pckl.dump(fmu,open('flux_flat_DPMJETIII_GH_numu.pckl', 'w'))
        return fmu(energy_use)
    
    ######### NU E #################################
    def make_fe(self,energy_use):
#        data = pckl.load(open('/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/tw_neutrino_flux_wsumsDPM_GH.pckl'))
     #   eedges = np.logspace(-0.5, 3., 301)

        #eedges = np.linspace(3., 1., 1000)
	eedges = np.linspace(1., 1000, 999)
        ecenters = (eedges[1:] + eedges[:-1])/2.
        ecenters_new = (eedges[1:] + eedges[:-1])/2.

        esteps = eedges[1:]-eedges[:-1]

        zedges = np.linspace(-1,1, 101)
        zcenters = (zedges[1:] + zedges[:-1])/2.

        #sum_numu_from_p = np.zeros_like(ecenters)
       # sum_numu_from_k = np.zeros_like(ecenters)
      #  sum_numubar_from_p = np.zeros_like(ecenters)
     #   sum_numubar_from_k = np.zeros_like(ecenters)
        sum_nue  = np.zeros_like(ecenters)
        sum_nuebar  = np.zeros_like(ecenters)

        for eii, energyi  in enumerate(ecenters):
                for zii , czenith  in enumerate(zcenters):
                       # sum_numu_from_p[eii] +=(self.spline_dict['numu_from_pion'].ev(energyi, czenith)/energyi**3)
                      #  sum_numu_from_p[eii] +=(self.spline_dict['numu_from_pion'].ev(energyi, czenith)/energyi**3)
                     #   sum_numu_from_k[eii] +=(self.spline_dict['numu_from_k'].ev(energyi, czenith)/energyi**3)
                    #    sum_numubar_from_p[eii] +=(self.spline_dict['antinum_from_pion'].ev(energyi, czenith)/energyi**3)
                   #     sum_numubar_from_k[eii] +=(self.spline_dict['antinum_from_k'].ev(energyi, czenith)/energyi**3)
                        sum_nue[eii]  +=(self.spline_dict['nue'].ev(energyi, czenith))
                        sum_nuebar[eii] +=(self.spline_dict['antinue'].ev(energyi, czenith))



       #note nue_tot_use_y is = nu_flux*E^3   (That is what the table is caluclated in)  
	#the confusing issue is I had, as a rule, been returning Flux (not Flux*E^3) 
        nue_tot_use_y =  (sum_nue + sum_nuebar  ) 
        fe = interp1d(ecenters, nue_tot_use_y, kind='cubic' )
        #pckl.dump(fe,open('flux_flat_DPMJETIII_GH_nue.pckl', 'w'))
        return fe(energy_use)
    
    ##################################################
    #  Function to Evaluate/make flat Splines
    ##########################################################
    
    ########## NU E #################################
    def EvaluateSplineEflat(self, spline_name = None, energy = 0, zenith = -1):
        spline_scaling_factorE = 1.0
        return (  (self.spline_dict[spline_name].ev(energy, zenith)) * spline_scaling_factorE/self.make_fe(energy))

    def EvaluateSplineE(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict[spline_name].ev(energy, zenith))) 
    
    ######### NU MU #################################
    def EvaluateSplineMuflat(self, spline_name = None, energy = 0, zenith = -1):
        spline_scaling_factorMu = 4.3
        return (  (self.spline_dict[spline_name].ev(energy, zenith)) * spline_scaling_factorMu/self.make_fmu(energy))
    
    def EvaluateSplineMu(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict[spline_name].ev(energy, zenith))) 

#this looks like an issue here.  the spline (to make flat) is caluclated in Flux*E*3 space .. yikes  maybe switch to doing this in Flux, not flux cubed?  but
# I want this to look flat in flux cubed space in the end. think on it a sec .. make some test plots. 



    #def __call__(self, particle = None, energy = None, zenith = None):
    
    

    #def EvaluateSpline(self, energy, zenith):
    #    result = self.SplinedFlux.ev(energy, zenith)
    #    return result

    def __init__(self):
        '''
        for each flux, 1) load the data , 2) make the spline 3) eval the spline
        '''
        #directory = "/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/"
        #directory = "/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/sybill_2.3_extended_tables/"
#	directory = "/home/trwood/tables/DPMJETIII_GH/"
	directory = "/home/trwood/tables_ICRC/DPMJETIII_GH/"

    #    directory = "/home/trwood/sybill_2.3_extended_tables_package/"
        #total fluxes from all sources for nu/anti nu ratio etc,

        self.spline_dict = {'nue' :self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "nue_totals.txt"),
                            'antinue':self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "antinue_totals.txt"),
                            #'numu':self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "numu_all_Source_table_fm.txt"),
                            #'antinumu':self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "antinumu_allSource.txt"),
                            'numu_from_k': self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "numu_from_kaon.txt"),
                            'antinum_from_k': self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "antinumu_from_kaon.txt"),
                            'numu_from_pion': self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "numu_from_pion.txt"),
                            'antinum_from_pion' : self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid.txt", directory + "antinumu_from_pion.txt")}



             #self.nue = self.EvaluateSpline(self, spline = nue, *kwargs)
