#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd

MatrixFileIn, MatrixFileOut = sys.argv[1:]

def quantile_normalize_using_target(x, target):
    """
    Both `x` and `target` are numpy arrays of equal lengths.
    """

    target_sorted = np.sort(target)

    return target_sorted[x.argsort().argsort()]

df=pd.read_csv(MatrixFileIn, index_col=0, sep='\t')
print(df)
