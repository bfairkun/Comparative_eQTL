#!/usr/bin/env python3
import pandas as pd
import numpy as np

a = pd.read_csv('Dispersion/BoostrapInference.PermutationsCombined.txt', header=0, sep='\t')
Observed=pd.read_csv('../../output/OverdispersionEstimatesFromChimp.txt', header=0, sep='\t')
print('read in files')

# Get distrubtuon of abs difference under null
AbsDifference = np.absolute(np.subtract(a[0:10000], a[10000:20000]))

# Get absolute difference
Actuals = np.absolute(np.subtract(Observed['Chimp.Residual'], Observed['Human.Residual']))
# Pvalue is fraction of null resamples greater than actual value. Add 1 to number of positive resamples, and cap at 1 to ensure range of P-values is (0, 1]
P=np.minimum(1.0, (np.sum(np.greater(AbsDifference, np.expand_dims(Actuals,0))) +1)/10000)
mask=np.isnan(Actuals)
P_masked = pd.DataFrame(np.ma.masked_array(P, mask))
P_masked['gene'] = Observed['gene']
P_masked.to_csv(path_or_buf="../../output/OverdispersionEstimatesFromChimp.txt.Pvals.tab", header=['P', 'gene'], na_rep="NA", sep='\t', index=False)
