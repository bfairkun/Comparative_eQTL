#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm

MatrixFileIn, MatrixFileOut = sys.argv[1:]

def quantile_normalize_using_target(x, target):
    """
    Both `x` and `target` are numpy arrays of equal lengths.
    """

    target_sorted = np.sort(target)

    return target_sorted[x.argsort().argsort()]

#Read file of N rows (genes) and M columns (individuals)
df=pd.read_csv(MatrixFileIn, index_col=0, sep='\t')

#Standardize across individuals to get Z scores.
StandardizedMatrix = df.sub(df.mean(1), axis=0).div(df.std(1), axis=0)

#Get N quantiles of normal distribution
NormalQuantiles = norm.ppf(np.linspace(0,1,df.shape[0]+2))[1:-1]

# Quantile normalize across genes.
QuantileNormalizedMatrix =  StandardizedMatrix.apply(quantile_normalize_using_target, axis=0, result_type='broadcast', target=NormalQuantiles)
QuantileNormalizedMatrix.to_csv(MatrixFileOut)

# Write out file
QuantileNormalizedMatrix.to_csv( MatrixFileOut, sep='\t', index_label="IID")
