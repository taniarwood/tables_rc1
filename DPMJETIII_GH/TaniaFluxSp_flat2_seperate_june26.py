import numpy as np
import scipy as sp
import pickle 
from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class MCEqFluxSpline(object):


    def LoadData(self, energylist_file, coszenlist_file, fluxfile2D):
        '''
        takes input of three files, for example energylist_file = "egrid_list.txt", coszenlist_file="cos_zenith_list.txt", fluxfile2D="numu_from_kaions_flux_table.txt"
        '''
        egrid = np.loadtxt(energylist_file)
        coszenlist = np.loadtxt(coszenlist_file)
        fluxfile2D = np.genfromtxt(fluxfile2D, delimiter=',')
        
        SplinedFlux= RectBivariateSpline(egrid,coszenlist, fluxfile2D , kx=2, ky=2)

        return SplinedFlux

    def  make_f2(self,energy_use, particle_type):
	
	partype = str(particle_type)
        data = pickle.load(open('/home/trwood/flux_reweight_package/tw_neutrino_flux_wsums.pckl'))

	numu_tot_use_y =  ( (data['sum_numu_from_p'] + data['sum_numubar_from_p'])*data['ecenters']**3 +  (data['sum_numu_from_p'] + data['sum_numubar_from_p'])*data['ecenters']**3 + (data['sum_numu_from_k'] + data['sum_numubar_from_k']) *data['ecenters']**3)
        f2 = interp1d(data['ecenters'], numu_tot_use_y, kind='cubic' )

	if partype == 'nue':
	    numu_tot_use_y = ( (data['sum_numu_from_p'] + data['sum_numubar_from_p'])*data['ecenters']**3)
	if partype == 'numu':
	    print ' numu hello again' 
	    numu_tot_use_y =  ( (data['sum_numu_from_p'] + data['sum_numubar_from_p'])*data['ecenters']**3 + (data['sum_numu_from_k'] + data['sum_numubar_from_k']) *data['ecenters']**3)

	elif partype == None: 
	    numu_tot_use_y =  ( (data['sum_numu_from_p'] + data['sum_numubar_from_p'])*data['ecenters']**3 +  (data['sum_numu_from_p'] + data['sum_numubar_from_p'])*data['ecenters']**3 + (data['sum_numu_from_k'] + data['sum_numubar_from_k']) *data['ecenters']**3)
	f2 = interp1d(data['ecenters'], numu_tot_use_y, kind='cubic' )


	return f2(energy_use)




    def EvaluateSpline(self, spline_name = None, flav = None, energy = 0, zenith = -1):
	#self.make_f2()
	flavor = str(flav)
	print flavor
	spline_scaling_factor = 5.5
	if flavor == 'numu':
	    print 'numu hello'
	    spline_scaling_factor = 4.0
	if flavor == 'nue':
	    spline_scaling_factor = 1.5

	elif flavor == None:
	    spline_scaling_factor = 5.5
        return (  (self.spline_dict[spline_name].ev(energy, zenith)/energy**3) * spline_scaling_factor/self.make_f2(energy,flavor))





    #def EvaluateSpline(self, energy, zenith):
    #    result = self.SplinedFlux.ev(energy, zenith)
    #    return result

    def __init__(self):
        '''
        for each flux, 1) load the data , 2) make the spline 3) eval the spline
        '''
        #directory = "/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/
        directory = "/home/trwood/flux_reweight_package/"

        #total fluxes from all sources for nu/anti nu ratio etc,

        self.spline_dict = {'nue' :self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "not_added_nue_all_table_fine_mesh.txt"),
                            'antinue':self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "antinue_allSource_table_fm.txt"),
                            'numu':self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "numu_all_Source_table_fm.txt"),
                            'antinumu':self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "antinumu_allSource_table_fm.txt"),
                            'numu_from_k': self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "numu_from_kaon_table_fine_mesh.txt"),
                            'antinum_from_k': self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "antinumu_from_kaon_table_fine_mesh.txt"),
                            'numu_from_pion': self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "numu_from_pion_table_fine_mesh.txt"),
                            'antinum_from_pion' : self.LoadData(directory + "egrid.txt", directory + "use_cos_zen_fine_mesh.txt", directory + "antinumu_from_pion_table_fm.txt")}



             #self.nue = self.EvaluateSpline(self, spline = nue, *kwargs)
