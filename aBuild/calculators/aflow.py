from os import path
from aBuild import msg # messaging module
from aBuild.utility import chdir
from aBuild.calculators.vasp import POSCAR, INCAR,KPOINTS,POTCAR
import os
import sys
import numpy as np

config = sys.modules["config"]  

class AFLOW:

    def __init__(self, specs,systemSpecies= None, directory=None):
        self.species = systemSpecies
        if isinstance(specs,dict):
            self._init_dict(specs)
        elif isinstance(specs,str):
            if self.species == None:
                msg.error("You can't initialize an AFLOW object using a path without also specifying the system species")
            self.directory = specs
            self._init_path()
               

    # If you initialize the AFLOW object with a dictionary, the dictionary *must* have the following entries:
    #               incar
    #               potcars
    #               kpoints
    #               crystal
    #               system species
    #  This most often happens when you are initialzing off of a yaml file and you augment the dictionary from the yaml
    # with the crystal and the species information.
            
    def _init_dict(self,specs):
        if specs['incar']['build'] != 'auto':
            self.INCAR = INCAR(specs["incar"])
            self.incar_auto_build = False
        else:
            self.INCAR = specs["incar"]
            self.incar_auto_build = True
        if specs['potcars']['build'] != 'auto':
            self.POTCAR = POTCAR(specs["potcars"])
            self.potcar_auto_build = False
        else:
            self.POTCAR = specs["potcars"]
            self.potcar_auto_build = True
            
        if specs['kpoints']['build'] != 'auto':
            self.KPOINTS = KPOINTS(specs["kpoints"])
            self.kpoints_auto_build = False
        else:
            self.KPOINTS = specs["kpoints"]
            self.kpoints_auto_build = True
        self.crystal = specs["crystal"]
        self.species = specs["species"]

    # If you initialize your AFLOW object with the path to the calculation, you must also specify the system's species
    #  This is necessary because often times the POSCAR/POTCAR files will not tell the whole story about what system you are
    # studying.  For example, the sytem may be a ternary, but this particular crystal may only have 2 of the 3 atom type.  If
    # the system species is not specified we have no way of knowing that we are missing an atomt type.
    # This most often happens when you are wanting to suck out the results from a calculation that has finished.
    def _init_path(self):
        from aBuild.database.crystal import Crystal
        from aBuild.calculators.vasp import POTCAR,KPOINTS

        #This routine is called when a path is passed in for initialization.
        # In this case, we want to read all the relevant information from the aflow.in file.
        from os import path

        self.POTCAR = POTCAR(path.join(self.directory,'POTCAR.static.xz'))
        self.KPOINTS = KPOINTS(path.join(self.directory,'KPOINTS.static.xz'))
        self.crystal = Crystal(path.join(self.directory,'POSCAR.static.xz'),self.species,crystalSpecies = self.POTCAR.species)
        
        



    def writeaflowfromposcar(self,filename):
        from os import waitpid
        from subprocess import Popen
        command = " {}  {} < {} > {}".format('aflow', '--poscar2aflowin', filename,'aflow.in')

        child=Popen(command, shell=True, executable="/bin/bash")
        waitpid(child.pid, 0)
        

    def modifyaflowin(self):
        lines = []
        with open('aflow.in','r') as f:
            currentlines = f.readlines()
        for line in currentlines:
            found = False
            for tag in self.INCAR:
                if isinstance(self.INCAR[tag],list) or isinstance(self.INCAR[tag],dict):
                    for subtag in self.INCAR[tag]:
                        if "[" + tag + "]" + subtag in line:
                            lines.append("[" + tag + "]" + subtag + "=" + self.INCAR[tag][subtag] + '\n')
                            found = True
                elif tag in line:
                    lines.append("[" + tag + "]" + self.INCAR[tag] + '\n')
                    found = True
            if not self.kpoints_auto_build:
                if '[VASP_KPOINTS_FILE]' in line:
                    found = True
                if '[VASP_KPOINTS_MODE_IMPLICIT]' in line:
                    lines.append('[VASP_KPOINTS_MODE_EXTERNAL]\n')
                    found = True
            if not found:
                lines.append(line)

        with open('aflow.in','w') as f:
            f.writelines(lines)

    def status(self):
        return 'done'

    def buildFolder(self):
        self.crystal.write('POSCAR.orig')
        self.check_atom_counts_zero()
        self.crystal.write('POSCAR')
        self.writeaflowfromposcar('POSCAR')
        self.modifyaflowin()
        self.KPOINTS.rGP = True
        self.KPOINTS.writeKPOINTS()
        

    def check_atom_counts_zero(self):
        from numpy import array,any
        if any(self.crystal.atom_counts == 0):
            from numpy import  where
            idxKeep = list(where( self.crystal.atom_counts > 0)[0])
            print(idxKeep, 'idx keep')
            print(self.species, 'species')
#            self.POTCAR.species = list(array(self.POTCAR.species)[idxKeep])
            self.crystal.atom_counts = self.crystal.atom_counts[idxKeep]
            self.crystal.species = list(array(self.species)[idxKeep])

    def read_results(self, allElectronic = False, allIonic=False):
        if self.directory is not None and self.status() in ['done','unconverged']:
            with chdir(self.directory):
                self.crystal.results = {}
                self.crystal.results["warning"] = False
                self.crystal.results["energyF"],self.crystal.results["forces"] = self.read_energy_and_forces()
                self.crystal.results["stress"] = self.read_stress()
                #self.POTCAR = POTCAR.from_POTCAR()
                self.crystal.results["species"] = self.POTCAR.species
                self.crystal.results["energypatom"] = self.crystal.results["energyF"]/self.crystal.nAtoms
                if abs(self.crystal.results["energyF"]) > 1000:
                    self.crystal.results["warning"] = True
                if 'pure' not in self.directory:
                    self.crystal.results["fEnth"] = self.formationEnergy
                else:
                    self.crystal.results["fEnth"] = 0
        else:
            print(self.status(), 'not reading results')
            self.crystal.results = None
            msg.info("Unable to extract necessary information from directory! ({})".format(self.directory))

    def read_energy_and_forces(self):
        import lzma
        from numpy import array
        energyZ = None
        energyF = None
        forces = []
#        if allElectronic:
#            energyF = []
 #           energyZ = []
        
        with lzma.open('aflow.qmvasp.out.xz','rt') as f:
            resultslines = f.readlines()

        for idx,line in enumerate(resultslines):
            # Free energy
            if 'TOTAL-FORCE' in line:
                for i in range(idx + 2, idx + 2 + self.crystal.nAtoms + 1):
                    forces.append(array(resultslines[i].split()[-3:]))
            if 'E_cell' in line:
                energyF = float(line.split()[0].split('=')[1])
            # Extrapolated zero point energy
#            if line.startswith('  energy  without entropy'):
#                if allElectronic:
#                    energyZ.append(float(line.split()[-1]))
#                else:
#                    energyZ = float(line.split()[-1])
                        
        return energyF,forces

    def read_stress(self):
        import lzma
        stress = None
        for line in lzma.open('OUTCAR.static.xz','rt'):
            if line.find(' in kB  ') != -1:
                stress = -np.array([float(a) for a in line.split()[2:]])
                stress = stress[[0, 1, 2, 4, 5, 3]] * 1e-1 #* .00624151# * ase.units.GPa.. Gets me to Giga Pascals              
        return stress
