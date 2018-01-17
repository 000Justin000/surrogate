#!/bin/bash
#SBATCH -J SOT1001
#SBATCH -p normal
#SBATCH --time=1000:00:00
#SBATCH -n 48
#SBATCH --nodelist=vega[02]

export QCSCRATCH="/scratch/jj585/SOT1001/"
export OMP_NUM_THREADS=48

rm $QCSCRATCH/* -rf 

source activate /home/fs01/jj585/usr/lib/python_envs/chempy36

python  rotSOT.py  di_1001  9  10000    1  2  3  4 \
                                        3  4  5  6 \
                                        4  5  6  7 \
                                        5  6  7  8 \
                                        6  7  8  9 \
                                        7  8  9 10 \
                                        8  9 10 11 \
                                        9 10 11 12 \
                                       11 12 13 14 > output
