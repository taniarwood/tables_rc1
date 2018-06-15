import numpy as np 
import scipy as sp
import pickle as pckl
from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


################## INFO  #############################
#  This class 
#
#  1) takes input of three files, for example energylist_file = "egrid_list.txt", coszenlist_file="cos_zenith_list.txt", fluxfile2D="numu_from_kaions_flux_table.txt"
#
#  2) Evalues the spline at a given point. The tables are in units of cos(zenith) egrid (energy in GeV for all functions, linear!)
#     the numbers in the tables are the fluxes predicted for a given cos(zeith) and energy. The fluxes are 
#     listed as flux * E^3. We divide the evaluate spline by energy cubed so that we return to the user the flux
#      NOT the Flux*E^3
#
#  3) Contains also a spline to the cos(zenith)integrated flux of the total numu flux and total Nue fluxes.
#	There are then Function to Evaluate/make flat fluxes in E^3 space, while mantaining the zeith shapes of the distributions
#  
#  4) Containss a simpple plotting function to test the basic loading of data and splining for items 1) and 2). For the more in depth varification of item
#     outlined in 3) please see
#    https://wiki.icecube.wisc.edu/index.php/IC86_Low_Energy_Atmospheric_Neutrino_Flux_Analysis_with_IceCube-DeepCore
#
#  Code property of the IceCube Collaboration, pls email Tania Wood: trwood@ualberta.ca for detailed questions 
#  Sept 29 2017
#####################################################


class MCEqFluxSpline(object):

    
####################################################
    #LOAD THE DATA
####################################################
    

    def LoadData(self, energylist_file, coszenlist_file, fluxfile2D):
        '''
        takes input of three files, for example energylist_file = "egrid_list.txt", coszenlist_file="cos_zenith_list.txt", fluxfile2D="numu_from_kaions_flux_table.txt"
        '''
        egrid = np.loadtxt(energylist_file)
        #print egrid
        coszenlist = np.loadtxt(coszenlist_file)
        fluxfile2D = np.genfromtxt(fluxfile2D, delimiter=' ')
        
        SplinedFlux= RectBivariateSpline(egrid,coszenlist, fluxfile2D , kx=2, ky=2)

        return SplinedFlux
    

	#Evalue the spline at a given point. The tables are in units of cos(zenith) egrid (energy in GeV)
	# the numbers in the tables are the fluxes predicted for a given cos(zeith) and energy. The fluxes are 
	#listed as flux * E^3. We divide the evaluate spline by energy cubed so that we return to the user the flux
	# NOT the Flux*E^3
  
	# Deliver the computed neurino Flux at the detector (0- 180 Deg) for a given Energy and zenith
    def EvaluateSpline(self, spline_name = None, energy = 0, zenith = -1):
        return self.spline_dict[spline_name].ev(energy, zenith)/energy**3

    def EvaluateSplinePionsOff(self, spline_name = None, energy = 0, zenith = -1):
        return self.spline_dict_poff[spline_name].ev(energy, zenith)/energy**3


    def Plot2Dflux(self, spline):
        x = np.linspace(0., 4, 100)
        y = np.linspace(-1, 1, 100)
        xx, yy = np.meshgrid(x, y)
        xx = 10**xx
        
        values = self.EvaluateSpline(spline, xx, yy)*(xx**3)
        plt.pcolor(x, y, values)
        plt.colorbar()
        plt.show()


######################################################################
    #Make the ID splines (one for numu, one for nue) and flat fluxes
#####################################################################
    def zenith_integrated(self):
	eedges = np.logspace(-0.5, 3.2, 301)
        ecenters = (eedges[1:] + eedges[:-1])/2.
        ecenters_new = (eedges[1:] + eedges[:-1])/2.

        esteps = eedges[1:]-eedges[:-1]

        zedges = np.linspace(-1,1, 101)
        zcenters = (zedges[1:] + zedges[:-1])/2.

        sum_numu_all = np.zeros_like(ecenters)
        sum_numubar_all = np.zeros_like(ecenters)

        sum_numu_totals_poff = np.zeros_like(ecenters)
        sum_numubar_totals_poff = np.zeros_like(ecenters)

        sum_nue  = np.zeros_like(ecenters)
        sum_nuebar  = np.zeros_like(ecenters)

        sum_nue_poff = np.zeros_like(ecenters)
        sum_nuebar_poff = np.zeros_like(ecenters)


        bin_width = 2.0/len(zcenters)
	#should use this bin_width but did not in orginal flat fluxes, so for today dont' use them either
	#bin_width = 1.0

        for eii, energyi  in enumerate(ecenters):
                for zii , czenith  in enumerate(zcenters):
                	sum_numu_totals_poff[eii] +=(self.spline_dict_poff['numu'].ev(energyi, czenith)/energyi**3)*bin_width 
                 	sum_numubar_totals_poff[eii] +=(self.spline_dict_poff['antinumu'].ev(energyi, czenith)/energyi**3)*bin_width 
                 	
			sum_numu_all[eii] +=(self.spline_dict['numu'].ev(energyi, czenith)/energyi**3)*bin_width 
                        sum_numubar_all[eii] +=(self.spline_dict['antinumu'].ev(energyi, czenith)/energyi**3)*bin_width 
                        
			sum_nue[eii]  +=(self.spline_dict['nue'].ev(energyi, czenith)/energyi**3)*bin_width 
                        sum_nuebar[eii] +=(self.spline_dict['antinue'].ev(energyi, czenith)/energyi**3)*bin_width 

                        sum_nue_poff[eii]  +=(self.spline_dict_poff['nue'].ev(energyi, czenith)/energyi**3)*bin_width
                        sum_nuebar_poff[eii] +=(self.spline_dict_poff['antinue'].ev(energyi, czenith)/energyi**3)*bin_width

	#all the iterations
	#NUMU
	   #from pions = (totals_numu regualar table) - (totals_numu pions off table)
	numu_p_tot_use_y =  (sum_numu_all + sum_numubar_all - sum_numu_totals_poff -  sum_numubar_totals_poff)  * ecenters**3
	   #from kaons = (totals_numu pions off table)
	numu_k_tot_use_y = ( sum_numu_totals_poff +  sum_numubar_totals_poff) * ecenters**3
	   #totals regular table
        numu_tot_use_y =  (sum_numu_all + sum_numubar_all)  * ecenters**3


	#NUE
	  #totals regular table Nue
	nue_tot_use_y = (sum_nue + sum_nuebar)  * ecenters**3
          #from pions = (totals_nue regualar table) - (totals_nue pions off table)
	nue_p_tot_use_y =  ((sum_nue + sum_nuebar) - ( sum_nue_poff + sum_nuebar_poff))  * ecenters**3
	   #from kaons = (totals_numu pions off table)
        nue_k_tot_use_y = ( sum_nue_poff + sum_nuebar_poff) * ecenters**3


	#make the zenith integrated splines to pass to the dictionary
        nue_zenI = interp1d(ecenters, nue_tot_use_y, kind='cubic' )
        numu_zenI = interp1d(ecenters, numu_tot_use_y, kind='cubic' )
        numu_k_zenI = interp1d(ecenters, numu_k_tot_use_y, kind='cubic' )
        numu_p_zenI = interp1d(ecenters, numu_p_tot_use_y, kind='cubic' )

        nue_k_zenI = interp1d(ecenters, nue_k_tot_use_y, kind='cubic' )
        nue_p_zenI = interp1d(ecenters, nue_p_tot_use_y, kind='cubic' )

	####################################
	# Dictionary of splines to return
	######################################
	zen_integrated_dictionary = { 
				'nue_zenI':nue_zenI, 
				'numu_zenI':numu_zenI, 
				'numu_k_zenI':numu_k_zenI,
				'numu_p_zenI':numu_p_zenI,
                                'nue_k_zenI':nue_k_zenI,
                                'nue_p_zenI':nue_p_zenI}
        #pckl.dump((nue_,open('splined_zenith_integ_DPMJETIII_h3a_nue.pckl', 'w'))
        return  zen_integrated_dictionary




####FOR TESTING
    def evalSplineNue_ZenI(self, energy = 10.0):
 #       spline = self.zenith_integrated_dict['nue_zenI']
	return   self.zenith_integrated_dict['nue_zenI'](energy) 

    def EvaluateSplineZen(self, spline_name = None, energy = 10.0):
        return self.zenith_integrated_dict[spline_name](energy)
    
    def adjust_pi_k(self, energy = 10.0, nu_pi_scale=1.0):
        new_numu_flux_total = self.zenith_integrated_dict['numu_k_zenI'](energy) + nu_pi_scale * (self.zenith_integrated_dict['numu_p_zenI'](energy))
	return new_numu_flux_total
	
    def correction_factor(self, energy, nu_pi_scale = 1.0):
	'''
        correction factor fuction added : cor_factor = new total flux /origial total flux
        energy is in GeV (linear)
        '''
	return  ( self.zenith_integrated_dict['numu_k_zenI'](energy) +  (self.zenith_integrated_dict['numu_p_zenI'](energy)) )/ (self.zenith_integrated_dict['numu_k_zenI'](energy) + nu_pi_scale * (self.zenith_integrated_dict['numu_p_zenI'](energy)) )

	
    def correction_factor_jp(self, energy, nu_pi_scale=1.0):
        '''
        Correction factor giving the exact same solution as the one above, but with some algebra to evaluate splines only twice instead of four times as above.       
        energy is in GeV (linear)
        '''
        return 1./((nu_pi_scale-1.)/(self.zenith_integrated_dict['numu_k_zenI'](energy)/self.zenith_integrated_dict['numu_p_zenI'](energy) + 1.) + 1.)

    def correction_factor_jp_e(self, energy, nu_pi_scale=1.0):
        '''
 #       Correction factor giving the exact same solution as the one above, but with some algebra to evaluate splines only twice instead of four times as above.       
 #       energy is in GeV (linear)
 #       '''
        return 1./((nu_pi_scale-1.)/(self.zenith_integrated_dict['nue_k_zenI'](energy)/self.zenith_integrated_dict['nue_p_zenI'](energy) + 1.) + 1.)



#  NOTE!!!  MCEq_rc1 (release candidate 1) has neutrinos from muons on their own.  This means they are missing from 
#           my fluxes! Since They come almost entirely from pions at these energies I will add them to the pion template
#           I will do (numu_total - kaons)  = rest ie "pions" hence forth (is primarily pions or from pion decay)
#

    
    #######################################################################
    #  Function to Evaluate/make 'flat' (flat in E^3 space) fluxes  Splines
    #######################################################################
    
    ########## NU E #################################
    def EvaluateSplineEflat(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict[spline_name].ev(energy, zenith)/energy**3) * self.spline_scaling_factorE/(self.zenith_integrated_dict['nue_zenI'](energy))    )
    def EvaluateSplineE(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict[spline_name].ev(energy, zenith)/energy**3))

    def EvaluateSplineEflat_poff(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict_poff[spline_name].ev(energy, zenith)/energy**3) * self.spline_scaling_factorE/(self.zenith_integrated_dict['nue_zenI'](energy))    )
    def EvaluateSplineE_poff(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict_poff[spline_name].ev(energy, zenith)/energy**3))

 
    
    ######### NU MU #################################
    def EvaluateSplineMuflat(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict[spline_name].ev(energy, zenith)/energy**3) * self.spline_scaling_factorMu/(self.zenith_integrated_dict['numu_zenI'](energy)))
    def EvaluateSplineMu(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict[spline_name].ev(energy, zenith)/energy**3)) 
    
    def EvaluateSplineMuflat_poff(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict_poff[spline_name].ev(energy, zenith)/energy**3) * self.spline_scaling_factorMu/(self.zenith_integrated_dict['numu_zenI'](energy)))
    def EvaluateSplineMu_poff(self, spline_name = None, energy = 0, zenith = -1):
        return (  (self.spline_dict_poff[spline_name].ev(energy, zenith)/energy**3))


#this looks like an issue here, but it is not. It does not matter to the fitter if we are in flux*E^3 space or not appartenly 
#(we are reduning Flux  the spline (to make flat) is caluclated in Flux*E*3 space .. yikes  maybe switch to doing this in Flux, not flux cubed?  but
# I want this to look flat in flux cubed space in the end so although the curernt implementation is not a problem it **might be
#easer to fit in flux*E^3 space instead of flux.. but that will make calculating nu-oscillations difficult.  Ask JP again?
    #def __call__(self, particle = None, energy = None, zenith = None):
    

    def __init__(self, intmodel='h3a', pmodel='DPMJET-III'):
        '''
        for each flux, 1) load the data , 2) make the spline 3) eval the spline
        '''
        intmod =  intmodel
	pmod = pmodel 
                    #/Users/trwood/Downloads/downloaded_notebokos/tables_rc1/DPMJETIII_h3a_rc1
	directory = '/home/trwood/tables_rc1/model_tests/tables_rc1_albrecht_v2/'+ pmod + '_' + intmod + str(".regular/") #'/home/trwood/tables_rc1/DPMJETIII_h3a_rc1/tables/'
	directory_p = '/home/trwood/tables_rc1/model_tests/tables_rc1_albrecht_v2/'+ pmod + '_' + intmod + str("/")
 #       directory_p = '/home/trwood/tables_rc1/model_tests/tables_rc1_albrecht/GH_SIBYLL2.3c/' #'/home/trwood/tables_rc1/h3a_DPMJETII_rc1_PionsOff/'
 #       directory_p = '/Users/trwood/Downloads/downloaded_notebokos/tables_rc1/h3a_DPMJETII_rc1_PionsOff/'
#	directory = '/Users/trwood/Downloads/downloaded_notebokos/tables_rc1/DPMJETIII_h3a_rc1/tables/'
        #total fluxes from all sources for nu/anti nu ratio etc,

	#table grid directory 
	dir_grid ='/home/trwood/tables_rc1/model_tests/tables_rc1_albrecht_v2/communial/' #'cos_zenith_grid_modeltests.txt  egrid.txt
        print directory
        print directory_p	
	#global data members
	self.spline_scaling_factorE = 0.020
        self.spline_scaling_factorMu = 0.100

        self.spline_scaling_factorE_poff = 0.020
        self.spline_scaling_factorMu_poff = 0.100

        self.spline_dict = {'nue' :self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory + "nue_totals.txt"),
                            'antinue':self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory + "antinue_totals.txt"),
                            'numu':self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory + "numu_totals.txt"),
                            'antinumu':self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory + "antinumu_totals.txt")}
                         #   'numu_from_k': self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid_modeltests.txt", directory + "numu_from_kaon.txt"),
                         #   'antinumu_from_k': self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid_modeltests.txt", directory + "antinumu_from_kaon.txt")}
                         #   'nue_from_pion': self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid_modeltests.txt", directory + "numu_from_pion.txt"),
                         #   'antinumu_from_pion' : self.LoadData(directory + "egrid.txt", directory + "cos_zenith_grid_modeltests.txt", directory + "antinumu_from_pion.txt")}


	self.spline_dict_poff = {'nue' :self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory_p + "nue_totals_pions_off.txt"),
                            'antinue':self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory_p + "antinue_totals_pions_off.txt"),
                            'numu':self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory_p + "numu_totals_pions_off.txt"),
                            'antinumu':self.LoadData(dir_grid + "egrid.txt", dir_grid + "cos_zenith_grid_modeltests.txt", directory_p + "antinumu_totals_pions_off.txt")}
                         #   'numu_from_k': self.LoadData(directory_p + "egrid.txt", directory_p + "cos_zenith_grid_modeltests.txt", directory_p + "numu_from_kaon_pions_off.txt"),
                         #   'antinumu_from_k': self.LoadData(directory_p + "egrid.txt", directory_p + "cos_zenith_grid_modeltests.txt", directory_p + "antinumu_from_kaon_pions_off.txt")}
                         #   'numu_from_pion': self.LoadData(directory_p + "egrid.txt", directory_p + "cos_zenith_grid_modeltests.txt", directory_p + "numu_from_pion_pions_off.txt"),
                         #   'antinumu_from_pion' : self.LoadData(directory_p + "egrid.txt", directory_p + "cos_zenith_grid_modeltests.txt", directory_p + "antinumu_from_pion_pions_off.txt")}

        self.zenith_integrated_dict = self.zenith_integrated()
 
