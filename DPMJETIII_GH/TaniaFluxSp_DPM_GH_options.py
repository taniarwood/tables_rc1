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
        return self.spline_dict[spline_name].ev(energy, zenith)/energy**3

    def Plot2Dflux(self, spline):
        x = np.linspace(0., 4, 100)
        y = np.linspace(-1, 1, 100)
        xx, yy = np.meshgrid(x, y)
        xx = 10**xx
        
        values = self.EvaluateSpline(spline, xx, yy)*(xx**3)
        plt.pcolor(x, y, values)
        plt.colorbar()
        plt.show()


################################################
    #Make the ID splines (one for numu, one for nue) and flat fluxes
####################################################
    ########## NU MU #################################
    def make_fmu(self,energy_use):
    
        #data = pckl.load(open('/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/tw_neutrino_flux_wsumsDPM_GH.pckl'))
        data= pckl.load(open('/home/trwood/tables_ICRC/DPMJETIII_GH/tw_neutrino_flux_wsumsDPM_GH.pckl'))
        
        numu_tot_use_y =  (( data['sum_numu_from_p'] + data['sum_numubar_from_p'] + data['sum_numu_from_k'] + data['sum_numubar_from_k']) *data['ecenters']**3)
        
        fmu = interp1d(data['ecenters'], numu_tot_use_y, kind='cubic' )
        #pckl.dump(fmu,open('flux_flat_DPMJETIII_GH_numu.pckl', 'w'))
        return fmu(energy_use)
    
    ######### NU E #################################
    def make_fe(self,energy_use):
        #data = pckl.load(open('/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/tw_neutrino_flux_wsumsDPM_GH.pckl'))
	data= pckl.load(open('/home/trwood/tables_ICRC/DPMJETIII_GH/tw_neutrino_flux_wsumsDPM_GH.pckl'))

        nue_tot_use_y =  (  (  data['sum_nue'] + data['sum_nuebar']  ) *   data['ecenters']**3 )
        fe = interp1d(data['ecenters'], nue_tot_use_y, kind='cubic' )
        #pckl.dump(fe,open('flux_flat_DPMJETIII_GH_nue.pckl', 'w'))
        return fe(energy_use)
    
    ##################################################
    #  Function to Evaluate/make flat Splines
    ##########################################################
    
    ########## NU E #################################
    def EvaluateSplineEflat(self, spline_name = None, energy = 0, zenith = -1):
        
        spline_scaling_factorE = 1.0
        return (  (self.spline_dict[spline_name].ev(energy, zenith)/energy**3) * spline_scaling_factorE/self.make_fe(energy))
    
    
    ######### NU MU #################################
    def EvaluateSplineMuflat(self, spline_name = None, energy = 0, zenith = -1):
        
        spline_scaling_factorMu = 4.3
        return (  (self.spline_dict[spline_name].ev(energy, zenith)/energy**3) * spline_scaling_factorMu/self.make_fmu(energy))
    




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
	directory = "/home/trwood/tables_ICRC/DPMJETIII_GH/"
#	directory = "/Users/trwood/may17/tables_ICRC/DPMJETIII_GH/"

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
