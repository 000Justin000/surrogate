#--------------------------------------------------
import sys, os, math, numpy
import openbabel, pybel, ase
#--------------------------------------------------
from mpi4py import MPI
#--------------------------------------------------
from ase import Atoms
from ase.io import read, write
from ase.visualize import view
from ase.calculators.qchem import QChem
from ase.constraints import FixInternals, Hookean
#--------------------------------------------------

#--------------------------------------------------
def pyb2ase(pybmol, pid):
    filename = "tmp"+"{:04d}".format(pid)+".xyz"
    pybmol.write("xyz", filename, overwrite=True)
    asemol = ase.io.read(filename)
    os.remove(filename)
    #--------------
    return asemol

#--------------------------------------------------
def ase2pyb(asemol, pid):
    filename = "tmp"+"{:04d}".format(pid)+".xyz"
    ase.io.write(filename, asemol, comment=filename)
    pybmol = next(pybel.readfile("xyz", filename))
    os.remove(filename)
    #--------------
    return pybmol

#--------------------------------------------------
def geomOptMM(pybmol, tcs, MMFF, tol):
    #----------------------------------------
    constraints = openbabel.OBFFConstraints()
    #----------------------------------------
    dihedrals = []
    #----------------------------------------
    if tcs is not None:
        for tc in tcs:
            pybmol.OBMol.SetTorsion(tc[0][0],tc[0][1],tc[0][2],tc[0][3], tc[1]/360.0*(2*math.pi))
#           constraints.AddTorsionConstraint(tc[0][0],tc[0][1],tc[0][2],tc[0][3], tc[1]/360.0*(2*math.pi))
            dihedrals.append(tc[0])
    #----------------------------------------
    atomsFixed = numpy.unique(sum(dihedrals, []))
    for atom in atomsFixed:
        constraints.AddAtomConstraint(int(atom))
    #----------------------------------------

    # pybview(pybmol, 0)
    #----------------------------------------
    FF = pybel._forcefields[MMFF]
    FF.Setup(pybmol.OBMol, constraints)
    FF.SetConstraints(constraints)
    #----------------------------------------
    EE = FF.Energy()
    dE = EE
    #----------------------------------------
    while(abs(dE / EE) > tol):
        #------------------------------------
    	FF.SteepestDescent(300)
    	dE = FF.Energy() - EE
    	EE = FF.Energy()
    #--------------
    FF.GetCoordinates(pybmol.OBMol)
    #--------------
    return pybmol, EE

#--------------------------------------------------
def pybview(pybmol, pid):
    pybmol.write("xyz", "tmp"+"{:04d}".format(pid)+".xyz", overwrite=True)
    os.system("avogadro tmp"+"{:04d}".format(pid)+".xyz")
    os.remove("tmp"+"{:04d}".format(pid)+".xyz")

#--------------------------------------------------
def getRMSD(pybmol1, pybmol2):
    alg = openbabel.OBAlign(pybmol1.OBMol, pybmol2.OBMol)
    alg.Align()
    return alg.GetRMSD()

#--------------------------------------------------
def getCoords(pybmol):
    coords = []
    for atom in pybmol:
        coords.append( list(atom.coords) )
    return coords

#--------------------------------------------------
def getPybmol(pybmol, coords):
    molr = pybmol.clone
    for atom, coord in zip(molr, coords):
        atom.OBAtom.SetVector(coord[0], coord[1], coord[2])
    return molr

def compareTorsion(pybmol1, pybmol2, rb):
    angle1 = pybmol1.OBMol.GetTorsion(rb[0],rb[1],rb[2],rb[3])/360*(2*math.pi)
    angle2 = pybmol2.OBMol.GetTorsion(rb[0],rb[1],rb[2],rb[3])/360*(2*math.pi)
    return (angle1-angle2) - math.floor((angle1-angle2)/(2*math.pi))*(2*math.pi)

#--------------------------------------------------
def compareFileTorsion(file1, file2, rb):
    file1_ext = os.path.splitext(file1)[1][1:]
    file2_ext = os.path.splitext(file2)[1][1:]
    pybmol1 = next(pybel.readfile(file1_ext, file1))
    pybmol2 = next(pybel.readfile(file2_ext, file2))
    return compareTorsion(pybmol1, pybmol2, rb)
