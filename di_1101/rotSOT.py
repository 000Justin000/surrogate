#--------------------------------------------------
import sys, os, math, numpy
import openbabel, pybel, ase
#--------------------------------------------------
from mpi4py import MPI
#--------------------------------------------------
from ase import Atoms
from ase.visualize import view
from ase.calculators.qchem import QChem
from ase.constraints import FixInternals, Hookean
#--------------------------------------------------
import models
import utils
import pySOT
from poap.controller import SerialController, BasicWorkerThread
# from utils import *
#--------------------------------------------------

#################################################
#           main function start here            #
#################################################

#------------------------------------------------
#                 customize area                #
#------------------------------------------------
print(str(sys.argv))
jobname  = str(sys.argv[1])
#------------------------------------------------
numDim   = int(sys.argv[2])
numPts   = int(sys.argv[3])
tsangles = []
#------------------------------------------------
for i in range(0, numDim):
    tsangles.append([int(sys.argv[(i+1)*4+0]), int(sys.argv[(i+1)*4+1]), int(sys.argv[(i+1)*4+2]), int(sys.argv[(i+1)*4+3])])     # dihedral idx
#------------------------------------------------
MMFF       = "mmff94s"
QMFUNC     = 'B3LYP'
DISPERSION = 'd3_op'
QMBASIS    = 'STO-3G'
TASK       = 'optimization'
#------------------------------------------------
MMtol = 1.0e-5
QMtol = 4.5e-3
ertol = 1.0e-10
#------------------------------------------------

#------------------------------------------------
#           initialize mpi parameters           #
#------------------------------------------------
nproc = MPI.COMM_WORLD.Get_size()
iproc = MPI.COMM_WORLD.Get_rank()
#------------------------------------------------

#------------------------------------------------
#         read the molecule with pybel          #
#------------------------------------------------
pybmol = next(pybel.readfile("xyz", jobname+".xyz"))
#------------------------------------------------

#------------------------------------------------
# (1) Optimization Problem
#------------------------------------------------
folding = models.E_QM(pybmol, tsangles)
print(folding.info)
#------------------------------------------------

#------------------------------------------------
# (2) Experimental Design
#------------------------------------------------
design = pySOT.SymmetricLatinHypercube(dim=folding.dim, npts=2*folding.dim+1)
#------------------------------------------------

#------------------------------------------------
# (3) Surrogate Model
#------------------------------------------------
surrogate = pySOT.RBFInterpolant(kernel=pySOT.CubicKernel, tail=pySOT.LinearTail, maxp=numPts)
#------------------------------------------------

#------------------------------------------------
# (4) Adaptive Sampling
#------------------------------------------------
sampling = pySOT.CandidateDYCORS(data=folding, numcand=100*folding.dim)
#------------------------------------------------

#------------------------------------------------
# (5) Controller
#------------------------------------------------
controller = SerialController(folding.objfunction)
#------------------------------------------------

#------------------------------------------------
# (6) Sychronous Strategy
#------------------------------------------------
controller.strategy = pySOT.SyncStrategyNoConstraints(worker_id=0, data=folding, maxeval=numPts, nsamples=1, exp_design=design, response_surface=surrogate, sampling_method=sampling)
#------------------------------------------------

result = controller.run()

fvals = numpy.array([o.value for o in controller.fevals])

# f, ax = plt.subplots()
# ax.plot(np.arange(0,numPts), fvals, 'bo')  # Points
# ax.plot(np.arange(0,numPts), np.minimum.accumulate(fvals), 'r-', linewidth=4.0)  # Best value found
# plt.xlabel('Evaluations')
# plt.ylabel('Function Value')
# plt.title(folding.info)
# plt.show()
# Print the final result

print('Best value found: {0}'.format(result.value))
print('Best solution found: {0}'.format(numpy.array_str(result.params[0], max_line_width=numpy.inf, precision=5, suppress_small=True)))

numpy.savetxt('configure-energy', numpy.concatenate((folding.paramTrace, folding.EETrace), axis=1))

folding.objview()
