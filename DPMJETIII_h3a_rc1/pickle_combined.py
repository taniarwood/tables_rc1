
import os, sys, pickle, shutil
import numpy as np

#import argparse

#sys.path.append('/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks')
#sys.path.append('/home/trwood/flux_reweight_package')
#sys.path.append('/home/trwood/tables_ICRC/DPMJETIII_GH')
sys.path.append('/home/tables_rc1/DPMJETIII_h3a_rc1')

#my spline module
#import TaniaFluxSp_DPM_h3a_rc1_options_update2_jaspert

##### NOW CHANGEING TO this script as this is TOTAL_NUMU instead of pi and k and is mvoing the pi totals to the muons 
import Combined_reg_table_pionsOFFTable_Keep_E3Flat_DPM_h3a_rc1
#import TaniaFluxSp_DPM_h3a_rc1_options_update2_jasper_gil_sumnumu_totals
#import Keep_E3Flat_DPM_h3a_rc1
'''    
parser = argparse.ArgumentParser(description='script for proccessing I3 files')
    parser.add_argument('infile', help='pickle_file to read in')
    parser.add_argument('ofile', help='output pickle file')

filename = args.inflie
outfilename = args.ofile
'''
flux= Combined_reg_table_pionsOFFTable_Keep_E3Flat_DPM_h3a_rc1.MCEqFluxSpline()
t = flux.EvaluateSpline('nue',10, 0.3)
k = flux.EvaluateSplinePionsOff('nue',10, 0.3)
print t, k
print 'try an eval new code'
test_flat =  flux.EvaluateSplineEflat('nue',10, 0.3)
print 'try an eval new code on zeith integrated'
test_zenithinteg =  flux.evalSplineNue_ZenI( 10.)
print test_zenithinteg, 'new dict of zenith inegrated splines'
test_e = flux.EvaluateSplineEflat('nue', 10, 0.3)
#print test
#quit() 
def redoFile(infile_name, outfile_name):

    data = pickle.load(open(infile_name))

    nubool = data['ptype'] > 0
 #   ebool  = data['energy'] < 100000.
    ebool = data['energy'] < 950. 
    factor = 1 
 #set up the weight keys
    data['weight_DMP_h3a_rc1_flat_e_k'] = np.zeros_like(data['weight_e'])
    data['weight_DMP_h3a_rc1_flat_e_p'] = np.zeros_like(data['weight_e'])
    data['weight_DMP_h3a_rc1_flat_mu_k'] = np.zeros_like(data['weight_e'])
    data['weight_DMP_h3a_rc1_flat_mu_p'] = np.zeros_like(data['weight_e'])

 #make the weights
    data['weight_DMP_h3a_rc1_flat_e_k'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSplineEflat_poff('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['weight_DMP_h3a_rc1_flat_e_k'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSplineEflat_poff('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor

    data['weight_DMP_h3a_rc1_flat_e_p'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  (flux.EvaluateSplineEflat('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor - flux.EvaluateSplineEflat_poff('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor)
    data['weight_DMP_h3a_rc1_flat_e_p'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  (flux.EvaluateSplineEflat('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor - flux.EvaluateSplineEflat_poff('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor)

    data['weight_DMP_h3a_rc1_flat_mu_k'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSplineMuflat_poff('numu',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor 
    data['weight_DMP_h3a_rc1_flat_mu_k'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSplineMuflat_poff('antinumu',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor
    
    data['weight_DMP_h3a_rc1_flat_mu_p'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  (flux.EvaluateSplineMuflat('numu',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor - flux.EvaluateSplineMuflat_poff('numu',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor)
    data['weight_DMP_h3a_rc1_flat_mu_p'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  (flux.EvaluateSplineMuflat('antinumu',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor - flux.EvaluateSplineMuflat_poff('antinumu',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor)

   

#NOW PUT IN KEYS WITH THE DMPJETIII h3a Spectrum as is, not flat
        #'weight_DMP_h3a_rc1_e_jaspert'
 #set up the weight keys
    data['weight_DMP_h3a_rc1_e_k'] = np.zeros_like(data['weight_e'])
    data['weight_DMP_h3a_rc1_e_p'] = np.zeros_like(data['weight_e'])

    data['weight_DMP_h3a_rc1_mu_k'] = np.zeros_like(data['weight_e'])
    data['weight_DMP_h3a_rc1_mu_p'] = np.zeros_like(data['weight_e'])

 #make the weights
    data['weight_DMP_h3a_rc1_e_k'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSplinePionsOff('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['weight_DMP_h3a_rc1_e_k'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSplinePionsOff('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor

    data['weight_DMP_h3a_rc1_e_p'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  (flux.EvaluateSpline('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor -  flux.EvaluateSplinePionsOff('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor)
    data['weight_DMP_h3a_rc1_e_p'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  (flux.EvaluateSpline('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor -  flux.EvaluateSplinePionsOff('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor)

    data['weight_DMP_h3a_rc1_mu_k'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSplinePionsOff('numu',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['weight_DMP_h3a_rc1_mu_k'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSplinePionsOff('antinumu',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor
   
    data['weight_DMP_h3a_rc1_mu_p'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  (flux.EvaluateSpline('numu',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor  -  flux.EvaluateSplinePionsOff('numu',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor)
    data['weight_DMP_h3a_rc1_mu_p'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]* (flux.EvaluateSpline('antinumu',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor -  flux.EvaluateSplinePionsOff('antinumu',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor)





    # 'numu_from_pion'
    #    n = 300000
    #    print data['weight_flat_e'][n:n+10]
    #    print data['weight_e'][n:n+10]
    #    print data['weight_e'].max(), data['weight_flat_e'].max()

    print 'Old weights', 
    print 'Stuff', np.sum(data['weight_e']), data['weight_e'].mean(), data['weight_e'].max(), data['weight_e'].min()
    print 'New weights', 
    print 'Stuff', np.sum(data['weight_DMP_h3a_rc1_flat_e_k']), data['weight_DMP_h3a_rc1_flat_e_k'].mean(), data['weight_DMP_h3a_rc1_flat_e_k'].max(), data['weight_DMP_h3a_rc1_flat_e_k'].min()

    bad_apples = data['weight_DMP_h3a_rc1_flat_e_k'] > data['weight_DMP_h3a_rc1_flat_e_k'].mean()

    print 'Energy', data['energy'][bad_apples][:10]
    print 'Zenith', data['zenith'][bad_apples][:10]
    print 'Type', data['ptype'][bad_apples][:10]
    print 'Weight_noflux', data['weight_noflux'][bad_apples][:10]

    print 'Flux', flux.EvaluateSpline('nue', data['energy'][bad_apples][:10], np.cos(data['zenith'][bad_apples][:10]))
    print 'WeightLarson', data['weight_e'][bad_apples][:10]
    #print 'tania' , data,  data['one'][bad_apples][:10] * 
    #print  flux.EvaluateSpline('nue', data['energy'][bad_apples][:10], np.cos(data['zenith'][bad_apples][:10]))
    #print np.sum(data['weight_e'][ebool]), np.sum(data['weight_mu'][ebool])
    #print np.sum(data['weight_flat_e']), np.sum(data['weight_flat_mu_k']) + np.sum(data['weight_flat_mu_k'])

    pickle.dump(data, open(outfile_name,'w'),protocol=-1)    

 #   pickle.dump(data, open(outfile_name,'w'))

if __name__=='__main__':
#    pickle_dir = '/home/trwood/MSU_sample/MSU_sample_sept2016/oscfit/MSU_tania_repickle/oscfitv2_repickle_protocol_minus1'
    #pickle_dir = '/home/trwood/MSU_contain_removed/flat_tania3'
#    pickle_dir = '/gs/project/ngw-282-ac/trwood/MSU_contain_removed/flat_tania3/'
#    pickle_dir = '/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_BKBK/flat_tania3_safe_backup'
#    out_dir    = '/Users/trwood/MSU_sample_sept2016/oscfit/MSU_tania_repickle_flat'

#    out_dir = '/gs/project/ngw-282-ac/trwood/MSU_contain_removed/flat_tania3_oct11/'
    #pickle_dir = '/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_BKBK/flat_tania3_jaspert2'
    pickle_dir  ='/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_BKBK/flat_tania3_jaspert2344'
    out_dir  ='/gs/project/ngw-282-ac/trwood/jasper_home/MSU_contain_removed/flat_tania3_DPM_interm_BAKOct12017_jaspert_BKBK/def2'
 
#   out_dir = '/home/trwood/MSU_sample/MSU_sample_sept2016/oscfit/MSU_tania_repickle_flat_fixed'
    all_files = os.listdir(pickle_dir)

    for one_file in all_files:
        print one_file
        full_filename = os.path.join(pickle_dir, one_file)
        out_filename  = os.path.join(out_dir, one_file)
        # Skip if the file is not a MC file
        if 'muongun' in one_file:
            shutil.copy(full_filename, out_filename)
            continue
        # Write the condition here


        redoFile(full_filename, out_filename)
