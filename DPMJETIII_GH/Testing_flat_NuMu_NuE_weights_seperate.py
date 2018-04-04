
# coding: utf-8

# In[1]:

#testing the flat weights


# In[2]:

import os, sys
# User: set your oscFit installation path here
sys.path.append('/Users/trwood/oscFit3D_v2.0_Tania/modules/')
#sys.path.append('/Users/trwood/oscFit3D_v2.0_Tania/modules')
#import matplotlib.pyplot as plt
#import jp_mpl as jplot
import numpy as np
from scipy import optimize
import pickle
#get_ipython().magic(u'matplotlib inline')


# In[3]:

from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline
import numpy as np
import scipy as sp
import pickle 


# In[ ]:
#DPMJETIII_GH
data = pickle.load(open('/home/trwood/tables_ICRC/DPMJETIII_GH/tw_neutrino_flux_wsumsDPM_GH.pckl'))
#data_flat = pickle.load(open('/Users/trwood/oscFit3D_v2.0_Tania/resources/ipython_notebooks/tw_neutrino_flux_flat.pckl'))
esteps = data['eedges'][1:]-data['eedges'][:-1]


# In[ ]:

sys.path.append('/Users/trwood/may17/tables_ICRC/DPMJETIII_GH')
import TaniaFluxSp_DPM_GH_options


# In[ ]:

#DPM_GH
flux= TaniaFluxSp_DPM_GH_options.MCEqFluxSpline()


# In[ ]:

eedges = np.logspace(-0.5, 3., 301)
ecenters = (eedges[1:] + eedges[:-1])/2.
ecenters_new = (eedges[1:] + eedges[:-1])/2.

print len(eedges),  len(ecenters)

esteps = eedges[1:]-eedges[:-1]

zedges = np.linspace(-1,1, 101)
zcenters = (zedges[1:] + zedges[:-1])/2.

print len(zedges),  len(zcenters)


znumu_from_p = np.zeros([ecenters.size, zcenters.size])
znumu_from_k = np.zeros([ecenters.size, zcenters.size])
znumubar_from_p = np.zeros([ecenters.size, zcenters.size])
znumubar_from_k = np.zeros([ecenters.size, zcenters.size])
ztotal_numu = np.zeros([ecenters.size, zcenters.size])

znue = np.zeros_like(znumu_from_p)
znuebar = np.zeros_like(znumu_from_p)
ztotal_nue = np.zeros_like(znumu_from_p)

znumu =  np.zeros_like(znumu_from_p)
znumubar =  np.zeros_like(znumu_from_p)


# In[ ]:

for ei, energy in enumerate(ecenters):
    for zi, czenith in enumerate(zcenters):
        znumu_from_p[ei,zi]    = flux.EvaluateSplineMuflat('numu_from_pion',energy, czenith)
        znumubar_from_p[ei,zi] = flux.EvaluateSplineMuflat('antinum_from_pion',energy, czenith)
        
        znumu_from_k[ei,zi]    = flux.EvaluateSplineMuflat('numu_from_k',energy, czenith)
        znumubar_from_k[ei,zi] = flux.EvaluateSplineMuflat('antinum_from_k',energy, czenith)
        
        ztotal_numu[ei,zi] =  (flux.EvaluateSplineMuflat('numu_from_pion',energy, czenith) + flux.EvaluateSplineMuflat('antinum_from_pion',energy, czenith)
            + flux.EvaluateSplineMuflat('numu_from_k',energy, czenith) + flux.EvaluateSplineMuflat('antinum_from_k',energy, czenith) )
        
        znue[ei,zi]     = flux.EvaluateSplineEflat('nue',energy, czenith)
        znuebar[ei,zi]  = flux.EvaluateSplineEflat('antinue',energy, czenith)
        
        ztotal_nue[ei,zi] =  flux.EvaluateSplineEflat('nue',energy, czenith) + flux.EvaluateSplineEflat('antinue',energy, czenith)
        
    

        #       znumu[ei,zi] = znumu_from_p[ei,zi] + znumu_from_k[ei,zi] 
        #       znumubar[ei,zi] =  znumubar_from_p[ei,zi] + znumubar_from_k[ei,zi]
        znumu[ei,zi] = flux.EvaluateSplineMuflat('numu_from_pion',energy, czenith) + flux.EvaluateSplineMuflat('numu_from_k',energy, czenith)
        znumubar[ei,zi] = flux.EvaluateSplineMuflat('antinum_from_pion',energy, czenith) +  flux.EvaluateSplineMuflat('antinum_from_k',energy, czenith)

    
sum_numu_from_p = np.zeros_like(data['ecenters'])
sum_numu_from_k = np.zeros_like(data['ecenters'])
sum_numubar_from_p = np.zeros_like(data['ecenters'])
sum_numubar_from_k = np.zeros_like(data['ecenters'])
sum_nue  = np.zeros_like(data['ecenters'])
sum_nuebar  = np.zeros_like(data['ecenters'])

for eii, energy  in enumerate(ecenters):
    for zii , czenith  in enumerate(zcenters):
        sum_numu_from_p[eii] +=flux.EvaluateSplineMuflat('numu_from_pion',energy, czenith)
        sum_numu_from_k[eii] +=flux.EvaluateSplineMuflat('numu_from_k',energy, czenith)
        sum_numubar_from_p[eii] +=flux.EvaluateSplineMuflat('antinum_from_pion',energy, czenith)
        sum_numubar_from_k[eii] +=flux.EvaluateSplineMuflat('antinum_from_k',energy, czenith)
        sum_nue[eii]  +=flux.EvaluateSplineEflat('nue', energy, czenith)
        sum_nuebar[eii] +=flux.EvaluateSplineEflat('antinue', energy, czenith)
        
        
        
pickle.dump({'numu_from_p':znumu_from_p, 'numubar_from_p':znumubar_from_p,
             'numu_from_k':znumu_from_k, 'numubar_from_k':znumubar_from_k,
             'nue_flux':znue, 'nuebar_flux':znuebar,
             'numu_flux':znumu, 'numubar_flux':znumubar,
             'numu_total':ztotal_numu, 'nue_total':ztotal_nue,
             'sum_numu_from_p':sum_numu_from_p,'sum_numu_from_k':sum_numu_from_k,
             'sum_numubar_from_p':sum_numubar_from_p, 'sum_numubar_from_k':sum_numubar_from_k,
             'sum_nue':sum_nue,'sum_nuebar':sum_nuebar,
             'eedges':eedges, 'ecenters':ecenters,
             'zedges':zedges, 'zcenters':zcenters},
            open('/home/trwood/tables_ICRC/DPMJETIII_GH/tw_neutrino_flux_wsums_DH_GH_flat_jun21.pckl','w'))


# In[ ]:



