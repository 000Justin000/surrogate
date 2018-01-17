import sys
import numpy as np
import math
import openbabel, pybel
import utils
from mpi4py import MPI
from ase.calculators.qchem import QChem

#---------------------------------------------------------------
class E_MM:
    #-----------------------------------------------------------
    def __init__(self, pybmol, tsangles, MMFF="uff", MMtol=1.0e-8):
    #-----------------------------------------------------------
        self.mol        = pybmol.clone
        self.tsangles   = tsangles
        self.MMFF       = MMFF
        self.MMtol      = MMtol
        self.dim        = len(tsangles)
        self.xlow       = -180 * np.ones(self.dim)
        self.xup        =  180 * np.ones(self.dim)
        self.info       = str(self.dim) + " dimensional folding problem"
        self.integer    = []
        self.continuous = np.arange(0, self.dim)
    

    #-----------------------------------------------------------
    def objfunction(self, tsvals):
    #-----------------------------------------------------------

        if (len(tsvals) != self.dim):
            #---------------------------------------------------
            raise ValueError("dimensional mismatch")
            #---------------------------------------------------
        else:
            #---------------------------------------------------
            molr = self.mol.clone
            #---------------------------------------------------
            tsconstraints = []
            for rb, tsval in zip(self.tsangles, tsvals):
                tsconstraints.append([rb, tsval])
            #---------------------------------------------------
            molr, EE = utils.geomOptMM(molr, tsconstraints, self.MMFF, self.MMtol)

            print(str(tsvals) + ": " + str(EE))
        return EE
    

    #-----------------------------------------------------------
    def objview(self):
    #-----------------------------------------------------------
        utils.pybview(self.mol, 0)
        #-------------------------------------------------------



#---------------------------------------------------------------
class E_QM:
    #-----------------------------------------------------------
    def __init__(self, pybmol, tsangles, MMFF="uff", MMtol=1.0e-8, QMFUNC="B3LYP", DISPERSION="d3_op", BASIS="STO-3G", tol_grad=300, tol_disp=1, tol_e=100):
    #-----------------------------------------------------------
        self.mol        = pybmol.clone
        self.tsangles   = tsangles
        self.MMFF       = MMFF
        self.MMtol      = MMtol
        self.QMFUNC     = QMFUNC
        self.DISPERSION = DISPERSION
        self.BASIS      = BASIS
        self.tol_grad   = tol_grad
        self.tol_disp   = tol_disp
        self.tol_e      = tol_e
        self.dim        = len(tsangles)
        self.xlow       = -180 * np.ones(self.dim)
        self.xup        =  180 * np.ones(self.dim)
        self.info       = str(self.dim) + " dimensional folding problem"
        self.integer    = []
        self.continuous = np.arange(0, self.dim)
        self.neval      = 0
        self.EE         = 0
        self.paramTrace = []
        self.EETrace    = []
    

    #-----------------------------------------------------------
    def objfunction(self, tsvals):
    #-----------------------------------------------------------
        nproc = MPI.COMM_WORLD.Get_size()
        iproc = MPI.COMM_WORLD.Get_rank()
	
        if (len(tsvals) != self.dim):
            #---------------------------------------------------
            raise ValueError("dimensional mismatch")
            #---------------------------------------------------
        else:
            #---------------------------------------------------
            molr = self.mol.clone
            #---------------------------------------------------
            tsconstraints = []
            for rb, tsval in zip(self.tsangles, tsvals):
                tsconstraints.append([rb, tsval])
            #---------------------------------------------------
            molr, EE = utils.geomOptMM(molr, tsconstraints, self.MMFF, self.MMtol)
            #---------------------------------------------------

            if (EE <= 1.0e4):
                #----------------------------------------
                asemol = utils.pyb2ase(molr, iproc)
                #----------------------------------------
                prefix = "{:04d}".format(self.neval)
                #----------------------------------------
                calc = QChem(xc=self.QMFUNC,
                             disp=self.DISPERSION,
                             basis=self.BASIS,
                             task="optimization",
                             symmetry=False,
                             tcs=tsconstraints,
                             opt_maxcycle=300,
                             opt_tol_grad=self.tol_grad,
                             opt_tol_disp=self.tol_disp,
                             opt_tol_e   =self.tol_e,
                             thresh=12,
                             scf_convergence=8,
                             maxfile=128,
                             mem_static=400,
                             mem_total=4000,
                             label="tmp_qchem"+"{:04d}".format(iproc)+"/" + prefix)
                asemol, EE = calc.run(asemol)
                #----------------------------------------
                if ((asemol is not None) and (EE is not None)):
                    if (np.isnan(EE)):
                        EE = 0

                    self.EE = EE

                    self.paramTrace.append(tsvals)
                    self.EETrace.append([EE])
                else:
                    EE = self.EE
            else:
                #----------------------------------------
                asemol = utils.pyb2ase(molr, iproc)
                #----------------------------------------
                prefix = "{:04d}".format(self.neval)
                #----------------------------------------
                calc = QChem(xc=self.QMFUNC,
                             disp=self.DISPERSION,
                             basis=self.BASIS,
                             task="single_point",
                             symmetry=False,
                             thresh=12,
                             scf_convergence=8,
                             maxfile=128,
                             mem_static=400,
                             mem_total=4000,
                             label="tmp_qchem"+"{:04d}".format(iproc)+"/" + prefix)
                EE = calc.run(asemol)
                #----------------------------------------
                if (EE is not None):
                    if (np.isnan(EE)):
                        EE = 0

                    self.EE = EE

                    self.paramTrace.append(tsvals)
                    self.EETrace.append([EE])
                else:
                    EE = 0

            sys.stdout.write(str(["{:+10.3f}".format(x) for x in tsvals]) + ": " + str(EE) + "\n")
            sys.stdout.flush()

            self.neval = self.neval + 1
        return EE
    
    #-----------------------------------------------------------
    def objview(self):
    #-----------------------------------------------------------
        utils.pybview(self.mol, 0)
        #-------------------------------------------------------
