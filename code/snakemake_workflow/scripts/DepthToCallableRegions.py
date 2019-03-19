#!/usr/bin/env python3

"""
Input from std in the results from samtools depth, and a MinDepth, HengLiConstant,
and GenomeCoverageFile require positional arguments. 

In this case MinDepth is the MinDepth required for all samples (3 is a good number) that will be left
in the output bedfile of callable regions. HengLiConstant is a number (3 to 4 is
suggested) to which all samples must have no more than
SampleGenomeCoverage+HengLiConstant*sqrt(SampleGenomeCoverage) depth at a position to be
output as callable. This is a recommendation from Heng Li's lab, as sites with
too much coverage tend to yield false positive heterozygous calls due to
paralogous mismappings. The GenomeCoverageFile is a tab delimited file with the
first column being the sample name and the second column being the genomewide
average coverage. The sample name is not used, just the order (the order of
the rows must exactly match the sample order given as stdin by the samtools
depth command.

Output to std out each position (bed format) if and only if all samples have depth between MinDepth and HengLiConstant

Usage:
samtools depth [MyBamFile(s)] | python DepthToCallableRegions [MinDepth]
[HengLiConstant] [GenomeCoverageFile] > CallableRegions.bed
"""

import sys
import numpy as np
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

MinDepth=int(sys.argv[1])
HengLiConstant=int(sys.argv[2])
GenomeCoverageFile=sys.argv[3]

GenomeCovOrderedList=[]
fhGenomeCov=open(GenomeCoverageFile, 'rU')
for line in fhGenomeCov:
    GenomeCovOrderedList.append(float(line.strip('\n').split('\t')[1]))
fhGenomeCov.close()

GenomeCovArray = np.asarray(GenomeCovOrderedList)
MaxAllowedCovArray = GenomeCovArray + float(HengLiConstant) * np.sqrt(GenomeCovArray)

for line in sys.stdin:
    chrom, pos, *depths = line.strip('\n').split('\t')
    SiteCoverageArray=np.asarray([int(i) for i in depths])
    if (SiteCoverageArray>=3).all() and (MaxAllowedCovArray - SiteCoverageArray > 0).all():
        print('{}\t{}\t{}'.format(chrom, int(pos)-1, pos))
