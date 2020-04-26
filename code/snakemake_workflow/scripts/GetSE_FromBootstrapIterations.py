#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# ChimpBootsrapReps = "Dispersion/BoostrapSE.Chimp.PermutationsCombined.txt"
# HumanBootsrapReps = "Dispersion/BoostrapSE.Human.PermutationsCombined.txt"
# ChimpSEOut = "Dispersion/ChimpSE.tab"
# HumanSEOut = "Dispersion/HumanSE.tab"

ChimpBootsrapReps, HumanBootsrapReps, ChimpSEOut, HumanSEOut = sys.argv[1:]

Chimp = pd.read_csv(ChimpBootsrapReps, header=0, sep='\t')
Human = pd.read_csv(HumanBootsrapReps, header=0, sep='\t')

pd.DataFrame(np.std(Chimp)).to_csv(path_or_buf=ChimpSEOut, header=["SE"], index_label="gene", sep='\t', na_rep="NA")
pd.DataFrame(np.std(Human)).to_csv(path_or_buf=HumanSEOut, header=["SE"], index_label="gene", sep='\t', na_rep="NA")
