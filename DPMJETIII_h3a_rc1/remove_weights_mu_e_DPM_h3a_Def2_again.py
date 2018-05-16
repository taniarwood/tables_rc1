import pickle
import numpy as np
import os, sys

#indir = '/home/trwood/MSU_contain_removed/flat_tania3/'
#indir = '/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3/'

#indir = '/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_h3a_rc1_plus_epos/'
#indir = '/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_BKBK/flat_tania3_jaspert2345/'
indir =  '/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_BKBK/def2/'
#indir     = '/home/trwood/MSU_contain_removed/flat_tania3'
# Selecting only baseline files
required  = ['600', 'tania']


				#tweight_DMP_h3a_rc1_e_jaspert
e_weights = ['weight_DMP_h3a_rc1_flat_e_k',  'weight_DMP_h3a_rc1_e_k','weight_DMP_h3a_rc1_flat_e_p',  'weight_DMP_h3a_rc1_e_p', 'weight_DMP_h3a_rc1_e','weight_DMP_h3a_rc1_flat_e']
#e_weights = ['tweight_newflat_e', 'tweight_e']

	      #tweight_DMP_h3a_rc1_flat_mu_k_jaspert'
mu_weights = ['weight_DMP_h3a_rc1_flat_mu_k', 'weight_DMP_h3a_rc1_flat_mu_p', 'weight_DMP_h3a_rc1_mu_k', 'weight_DMP_h3a_rc1_mu_p']
#mu_weights = ['weight_mu']
#mu_weights = ['tweight_newflat_mu_k', 'tweight_newflat_mu_p', 'tweight_mu_p', 'tweight_mu_k']


all_files = os.listdir(indir)
for one_file in all_files:
    take_file = True
    for one_req in required:
        if one_req not in one_file:
            take_file = False
    if take_file:
        print one_file

        split_name = one_file.split('tania')

        # This is to remove the muon fluxes
        eflux_only = split_name[0] + 'eflux' + split_name[1]
        print eflux_only
        infile = open(os.path.join(indir, one_file))
        data = pickle.load(infile)
        for one_key in mu_weights:
            data[one_key] = np.zeros_like(data[one_key], dtype=int)
        outfile = open(os.path.join(indir, eflux_only), 'w')
        pickle.dump(data, outfile, protocol=-1)
        outfile.close()
        infile.close()

        # This is to remove the electron fluxes
        muflux_only = split_name[0] + 'muflux' + split_name[1]
        print muflux_only
        infile = open(os.path.join(indir, one_file))
        data = pickle.load(infile)
        for one_key in e_weights:
            data[one_key] = np.zeros_like(data[one_key], dtype=int)
        outfile = open(os.path.join(indir, muflux_only), 'w')
        pickle.dump(data, outfile, protocol=-1)
        outfile.close()
        infile.close()

