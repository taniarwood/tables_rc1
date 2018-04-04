
import os, sys, pickle, shutil
import numpy as np

#import argparse

#sys.path.append('/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks')
#sys.path.append('/home/trwood/flux_reweight_package')
sys.path.append('/home/trwood/tables_ICRC/DPMJETIII_GH')
#my spline module
import TaniaFluxSp_DPM_GH_options_update2_jasper

'''    
parser = argparse.ArgumentParser(description='script for proccessing I3 files')
    parser.add_argument('infile', help='pickle_file to read in')
    parser.add_argument('ofile', help='output pickle file')

filename = args.inflie
outfilename = args.ofile
'''


#flux = TaniaFluxSp.MCEqFluxSpline()
#flux= TaniaFluxSp_flat2.MCEqFluxSpline()
flux = TaniaFluxSp_DPM_GH_options_update2_jasper.MCEqFluxSpline()
test = flux.EvaluateSpline('nue',10, 0.3)
#test_mu = flux.EvaluateSplineMuflat('numu', 10, 0.3)
test_e = flux.EvaluateSplineEflat('nue', 10, 0.3)
#print test

def redoFile(infile_name, outfile_name):

    data = pickle.load(open(infile_name))

    nubool = data['ptype'] > 0
 #   ebool  = data['energy'] < 100000.
    ebool = data['energy'] < 950. 
    factor = 1 

    data['tweight_DMP_GH_flat_e_jasper'] = np.zeros_like(data['weight_e'])
    data['tweight_DMP_GH_flat_mu_k_jasper'] = np.zeros_like(data['weight_e'])
    data['tweight_DMP_GH_flat_mu_p_jasper'] = np.zeros_like(data['weight_e'])


    data['tweight_DMP_GH_flat_e_jasper'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSplineEflat('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['tweight_DMP_GH_flat_e_jasper'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSplineEflat('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor

    data['tweight_DMP_GH_flat_mu_k_jasper'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSplineMuflat('numu_from_k',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['tweight_DMP_GH_flat_mu_k_jasper'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSplineMuflat('antinum_from_k',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor
    # 'tweight_flat_mu_pi' and mubar
    data['tweight_DMP_GH_flat_mu_p_jasper'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSplineMuflat('numu_from_pion',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['tweight_DMP_GH_flat_mu_p_jasper'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSplineMuflat('antinum_from_pion',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor

   

#NOW PUT IN KEYS WITH THE DMPJETIII GH Spectrum as is, not flat

    data['tweight_DMP_GH_e_jasper'] = np.zeros_like(data['weight_e'])
    data['tweight_DMP_GH_mu_k_jasper'] = np.zeros_like(data['weight_e'])
    data['tweight_DMP_GH_mu_p_jasper'] = np.zeros_like(data['weight_e'])

    data['tweight_DMP_GH_e_jasper'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSpline('nue',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['tweight_DMP_GH_e_jasper'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSpline('antinue',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor

    data['tweight_DMP_GH_mu_k_jasper'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSpline('numu_from_k',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['tweight_DMP_GH_mu_k_jasper'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSpline('antinum_from_k',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor
    # 'tweight_flat_mu_pi' and mubar
    data['tweight_DMP_GH_mu_p_jasper'][nubool*ebool] = data['weight_noflux'][nubool*ebool]*  flux.EvaluateSpline('numu_from_pion',data['energy'][nubool*ebool], np.cos(data['zenith'][nubool*ebool]))*factor
    data['tweight_DMP_GH_mu_p_jasper'][~nubool*ebool] = data['weight_noflux'][~nubool*ebool]*  flux.EvaluateSpline('antinum_from_pion',data['energy'][~nubool*ebool], np.cos(data['zenith'][~nubool*ebool]))*factor



    # 'numu_from_pion'
    #    n = 300000
    #    print data['tweight_flat_e'][n:n+10]
    #    print data['weight_e'][n:n+10]
    #    print data['weight_e'].max(), data['tweight_flat_e'].max()

    print 'Old weights', 
    print 'Stuff', np.sum(data['weight_e']), data['weight_e'].mean(), data['weight_e'].max(), data['weight_e'].min()
    print 'New weights', 
    print 'Stuff', np.sum(data['tweight_DMP_GH_flat_e_jasper']), data['tweight_DMP_GH_flat_e_jasper'].mean(), data['tweight_DMP_GH_flat_e_jasper'].max(), data['tweight_DMP_GH_flat_e_jasper'].min()

    bad_apples = data['tweight_DMP_GH_flat_e_jasper'] > data['tweight_DMP_GH_flat_e_jasper'].mean()

    print 'Energy', data['energy'][bad_apples][:10]
    print 'Zenith', data['zenith'][bad_apples][:10]
    print 'Type', data['ptype'][bad_apples][:10]
    print 'Weight_noflux', data['weight_noflux'][bad_apples][:10]

    print 'Flux', flux.EvaluateSpline('nue', data['energy'][bad_apples][:10], np.cos(data['zenith'][bad_apples][:10]))
    print 'WeightLarson', data['weight_e'][bad_apples][:10]
    #print 'tania' , data,  data['one'][bad_apples][:10] * 
    #print  flux.EvaluateSpline('nue', data['energy'][bad_apples][:10], np.cos(data['zenith'][bad_apples][:10]))
    #print np.sum(data['weight_e'][ebool]), np.sum(data['weight_mu'][ebool])
    #print np.sum(data['tweight_flat_e']), np.sum(data['tweight_flat_mu_k']) + np.sum(data['tweight_flat_mu_k'])

    pickle.dump(data, open(outfile_name,'w'),protocol=-1)    

 #   pickle.dump(data, open(outfile_name,'w'))

if __name__=='__main__':
#    pickle_dir = '/home/trwood/MSU_sample/MSU_sample_sept2016/oscfit/MSU_tania_repickle/oscfitv2_repickle_protocol_minus1'
    pickle_dir = '/home/trwood/MSU_contain_removed/flat_tania3/no_Nu'
#    out_dir    = '/Users/trwood/MSU_sample_sept2016/oscfit/MSU_tania_repickle_flat'
    pickle_dir = '/home/trwood/MSU_contain_removed/flat_tania3/no_Nu_jasperin'
#    out_dir = '/home/trwood/MSU_contain_removed/flat_tania3_DPM_interm_BAKaug29_jasper'
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
