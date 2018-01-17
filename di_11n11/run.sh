#!/bin/bash
#SBATCH -J SOT11n11
#SBATCH -p normal
#SBATCH --time=1000:00:00
#SBATCH -n 24
#SBATCH --nodelist=vega[01]

export QCSCRATCH="/scratch/jj585/11n11/"
export OMP_NUM_THREADS=24

rm $QCSCRATCH/* -rf 

# python  rotSOT.py  scan_1001  2   10  1  2  4  8 \
#                                       2  1  3  6

source activate /home/fs01/jj585/usr/lib/python_envs/chempy36

python  rotSOT.py  di_11n11  9  10000    1  2  3  4 \
				 	 3  4  5  6 \
                                         4  5  6  7 \
                                         5  6  7  8 \
                                         6  7  8  9 \
                                         7  8  9 10 \
                                         8  9 10 11 \
                                         9 10 11 12 \
                                        11 12 13 14
