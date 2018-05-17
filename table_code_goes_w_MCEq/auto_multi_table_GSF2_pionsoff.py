#!/usr/bin/env python 

import os
#import matplotlib.pyplot as plt
import numpy as np
import math
os.chdir('/home/trwood/MCEq_rc1/')
import sys

import argparse
#from pathlib import Path
#import solver related modules
from MCEq.core import MCEqRun
from mceq_config import config
#import primary model choices
import CRFluxModels as pm

from MCEq.misc import  cornertext
#from MCEq.misc import set_ticks_y, cornertext
#import matplotlib
import cPickle as pickle
from os.path import join
import calendar

#MCEQ PATHS

sys.path.append('/home/trwood/MCEq_rc1/')
sys.path.append('/home/trwood/MCEq_rc1/MCEq')

data_dir = '/home/trwood/MCEq_rc1/data'

#####add function to configure MCQ without some options

def mceq_config_without(key_list):
    r = dict(config)  # make a copy
    for key in key_list:
        del r[key]
    return r

#####add function to automatically create directories if they do not yet exist 

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory) #python 3
        #os.mkdir()
    elif os.path.exists(directory):
        print 'already created'

#FOR system portability, get the path of the home directory + Desktop
desktop = os.path.join(os.path.expanduser("~"),'Desktop')

# This is the modification function, it can explicitely depend on energy
# Particle production matrices will be multiplied with it
def mod_linear(xmat, e_grid, a,*args):
    return (1. + a)*np.ones_like(xmat)

def apply_mod_p(fac_pi_pl, fac_pi_mi, MCEq_obj):

    # Delay re-initialization until all modifications are set
    r = MCEq_obj.set_mod_pprod(2212,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -211, mod_linear, (fac_pi_mi,), delay_init=True)

    r += MCEq_obj.set_mod_pprod(2112,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -211, mod_linear, (fac_pi_mi,), delay_init=True)

    if r > 0:
        MCEq_obj._init_default_matrices(skip_D_matrix=True)


def apply_mod(fac_pi_pl, fac_pi_mi, fac_k_pl, fac_k_mi, fac_k0, fac_k0s, fac_k0l, MCEq_obj):
    #month = str(month)
    # Delay re-initialization until all modifications are set
    r = MCEq_obj.set_mod_pprod(2212,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -211, mod_linear, (fac_pi_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212,  321, mod_linear, (fac_k_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -321, mod_linear, (fac_k_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, 311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, 130, mod_linear, (fac_k0l,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, 310, mod_linear, (fac_k0s,), delay_init=True)


    r += MCEq_obj.set_mod_pprod(2112,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -211, mod_linear, (fac_pi_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112,  321, mod_linear, (fac_k_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -321, mod_linear, (fac_k_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, 311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, 130, mod_linear, (fac_k0l,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, 310, mod_linear, (fac_k0s,), delay_init=True)


    if r > 0:
        MCEq_obj._init_default_matrices(skip_D_matrix=True)

#read in aruguments ##############################################################

parser = argparse.ArgumentParser(description='script for calculating neutrino fluxes')
#parser.add_argument('theta', help='the zenith angle to calculate the density profile with')
parser.add_argument('CRFluxModel', help='CRFlux model to use')
parser.add_argument('hadmodel', help='hadronic model to use')

#parser.add_argument('ofile', help='name of output file for the table')
#parser.add_argument('-s', '--sim', help='turn on extra processing for sim files', action='store_true')

args = parser.parse_args()
crmodel=int(args.CRFluxModel)
intmodel=str(args.hadmodel)
#GlobalSplineFit, None,
if crmodel ==1:
	pmodel = (pm.GlobalSplineFitBeta , None)
elif crmodel == 2:
	pmodel = (pm.HillasGaisser2012,'H3a') 
elif crmodel ==3:
	pmodel = (pm.GaisserHonda,'GH')

print pmodel



cos_theta = np.linspace(-1, 1, 41)
angles = np.arccos(cos_theta )* (180./ np.pi )

#pmodel = (pm.HillasGaisser2012, None)  #  (pm.Thunman, None)   # (pm.GaisserStanevTilav, "3-gen") # (PolyGonato, False)
#pmodel =(pm.HillasGaisser2012, "H3a") 

#pmodel = (pm.GaisserHonda,None)
#intmodel= 'DPMJET-III-2017.1'
mag = 3

#hold_table_modelname=   
#plan to have the table names twice i guess, have to have the pions off version 

modelname_numu = np.zeros(shape=(88,41))
modelname_antinumu = np.zeros(shape=(88,41))
modelname_nue = np.zeros(shape=(88,41))
modelname_antinue =np.zeros(shape=(88,41))

months_short = ['January','March','May','July','September','November']

for zii, zenith in enumerate(angles):
    print 'zii', zii, 
    
    #Initialize empty grid
    flux = {}
    for frac in ['numu_total','antinumu_total','nue_total','antinue_total']:
        flux[frac] = []

    for month in months_short:
            mceq_run = MCEqRun(interaction_model=intmodel, density_model=('MSIS00_IC',('SouthPole',month)),
                               primary_model=pmodel, theta_deg=zenith, **mceq_config_without(['density_model']) )
	    mceq_run.unset_mod_pprod(dont_fill=False) #clear any previous mod
            mod_0 = (-1.0,-1.0,mceq_run)   #set new mod w charged pion first parent fluxes now off (reinteractions too? check)
            apply_mod_p(*mod_0)
            #mceq_run.y.print_mod_pprod()

            mceq_run.solve()

            #print mceq_run.get_solution('total_numu', mag)
            #apparently these want to be initalized... 
            #numus  +=mceq_run.get_solution('total_numu', mag)
            flux['numu_total'].append(mceq_run.get_solution('total_numu', mag))

            flux['antinumu_total'] = (mceq_run.get_solution('total_antinumu', mag))

            flux['nue_total'] = (mceq_run.get_solution('total_nue', mag))

            flux['antinue_total'] = (mceq_run.get_solution('total_antinue', mag))

    #average by the number of months        
    #flux_year_avg_thisdeg['numu_only']= np.array(flux_year_avg_thisdeg['numu_only'])/12.0
    flux['numu_total']   = np.array(flux['numu_total'])
    flux['antinumu_total']= np.array(flux['antinumu_total'])
    flux['nue_total']     = np.array(flux['nue_total'])
    flux['antinue_total'] = np.array(flux['antinue_total'])

    
    flux['numu_total']   = np.sum(flux['numu_total']/6.0, axis = 0)
    flux['antinumu_total']= np.sum(flux['antinumu_total']/6.0, axis = 0)
    flux['nue_total']     = np.sum(flux['nue_total']/6.0, axis = 0)
    flux['antinue_total'] = np.sum(flux['antinue_total']/6.0, axis = 0)
    
       
    modelname_numu[:,zii]      = flux['numu_total'] 
    modelname_antinumu[:,zii]  = flux['antinumu_total']
    modelname_nue[:,zii]       = flux['nue_total']
    modelname_antinue[:,zii]   = flux['antinue_total']
    
    if pmodel[1] == None:
	file_paths = desktop+"/tables_rc1_albrecht/GSF_"+intmodel+"/numu_totals_pions_off.txt"
    else:

        file_paths = desktop+"/tables_rc1_albrecht/"+pmodel[1]+"_"+intmodel+"/numu_totals_pions_off.txt"
#file_paths ='/Users/trwood/my/franks2/file2.txt'
    ensure_dir(file_paths)
 #lib Path does not exist on cobalt.. even though I installed in my env.. (pathlib) 
 #  try:
   #     if not path.parent.exists():
  #          path.parent.mkdir(parents=True)
  #  except OSError:
  #      print "handle error"
    if pmodel[1] == None:
        np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/GSF_"+intmodel, 'numu_totals_pions_off.txt'),'w'),modelname_numu)
    	np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/GSF_"+intmodel, 'antinumu_totals_pions_off.txt'),'w'),modelname_antinumu)
    	np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/GSF_"+intmodel, 'nue_totals_pions_off.txt'),'w'),modelname_nue)
        np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/GSF_"+intmodel, 'antinue_totals_pions_off.txt'),'w'),modelname_antinue)
     
    else:
        np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/"+pmodel[1]+"_"+intmodel, 'numu_totals_pions_off.txt'),'w'),modelname_numu)
        np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/"+pmodel[1]+"_"+intmodel, 'antinumu_totals_pions_off.txt'),'w'),modelname_antinumu)
        np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/"+pmodel[1]+"_"+intmodel, 'nue_totals_pions_off.txt'),'w'),modelname_nue)
        np.savetxt(open(os.path.join(desktop+"/tables_rc1_albrecht/"+pmodel[1]+"_"+intmodel, 'antinue_totals_pions_off.txt'),'w'),modelname_antinue)
