from os import path
from aBuild import msg # messaging module
from aBuild.utility import chdir
from aBuild.calculators.vasp import POSCAR, INCAR,KPOINTS,POTCAR
import os
import sys
import numpy as np

config = sys.modules["config"]  

class AFLOW:

    def __init__(self, INCAR,KPOINTS,POTCAR,crystal,aflowin,directory,autobuild = {"incar":True, "kpoints":True,"potcar":True}):
        self.INCAR = INCAR
        self.KPOINTS = KPOINTS
        self.POTCAR = POTCAR
        self.crystal = crystal
        self.aflowin = aflowin
        self.directory = directory
        self.incar_auto_build = autobuild["incar"]
        self.kpoints_auto_build = autobuild["kpoints"]
        self.potcar_auto_build = autobuild["potcar"]
#        self.species = systemSpecies
#        if isinstance(specs,dict):
#            self._init_dict(specs)
#        elif isinstance(specs,str):
#            if self.species == None:
#                msg.error("You can't initialize an AFLOW object using a path without also specifying the system species")
#            self.directory = specs
#            self._init_path()
               

    # If you initialize the AFLOW object with a dictionary, the dictionary *must* have the following entries:
    #               incar
    #               potcars
    #               kpoints
    #               crystal
    #               system species
    #  This most often happens when you are initialzing off of a yaml file and you augment the dictionary from the yaml
    # with the crystal and the species information.
    @staticmethod
    def from_dictionary(specsDict):
       
        aflowin = specsDict['aflowin']
        auto_build = {}
        if specsDict['incar']['build'] != 'auto':
            incarobj = INCAR(specsDict["incar"])
            auto_build["incar"] = False
        else:
            incarobj = specsDict["incar"]
            auto_build["incar"] = True
        if specsDict['potcars']['build'] != 'auto':
            potcarobj = POTCAR(specsDict["potcars"])
            auto_build["potcar"] = False
        else:
            potcarobj = specsDict["potcars"]
            auto_build["potcar"] = True
            
        if specsDict['kpoints']['build'] != 'auto':
            kpointsobj = KPOINTS(specsDict["kpoints"])
            auto_build["kpoints"] = False
        else:
            kpointsobj = specsDict["kpoints"]
            auto_build["kpoints"] = True

            
        crystal = specsDict["crystal"]
        return AFLOW(incarobj,kpointsobj,potcarobj,crystal,aflowin,autobuild = auto_build)
#        species = specsDict["species"]

    # If you initialize your AFLOW object with the path to the calculation, you must also specify the system's species
    #  This is necessary because often times the POSCAR/POTCAR files will not tell the whole story about what system you are
    # studying.  For example, the sytem may be a ternary, but this particular crystal may only have 2 of the 3 atom type.  If
    # the system species is not specified we have no way of knowing that we are missing an atomt type.
    # This most often happens when you are wanting to suck out the results from a calculation that has finished.

    @staticmethod
    def from_path(directory,species):
        from aBuild.database.crystal import Crystal
        from aBuild.calculators.vasp import POTCAR,KPOINTS

        #This routine is called when a path is passed in for initialization.
        # In this case, we want to read all the relevant information from the aflow.in file.
        from os import path

        potcar = POTCAR.from_path(path.join(directory,'POTCAR.static.xz'))
        kpoints = KPOINTS.from_path(directory)
        incar = INCAR.from_path(path.join(directory,'INCAR.static.xz'))
        crystal = Crystal.from_path(path.join(directory,'POSCAR.static.xz'),species)
        
        return AFLOW(incar,kpoints,potcar,crystal,None,directory)
        
        



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
#            found = False
            needsRemoved = (True in [x in line for x in self.aflowin['remove']]) or (not self.kpoints_auto_build and '[VASP_KPOINTS_FILE]' in line)
            
            if 'add' in self.aflowin:
                needsAdded = True in [x in line for x in self.aflowin['add']]
            else:
                needsAdded = False

            if True in [x in line for x in self.aflowin['replace']]:
                needsReplaced  = True
                replaceTag = line.split()[0]
                replacement = self.aflowin['replace'][replaceTag] + '\n'
            elif not self.kpoints_auto_build and '[VASP_KPOINTS_MODE_IMPLICIT]' in line:
                needsReplaced = True
                replacement = '[VASP_KPOINTS_MODE_EXTERNAL]\n'
            else:
                needsReplaced = False


            if needsRemoved:
                lines.append('#aBuild' + line)
#                continue
            elif needsReplaced:
                lines.append(replacement)
            elif needsAdded:
                print('Added branch triggered')
            else:
                lines.append(line)
#            elif not self.kpoints_auto_build:
#                if '[VASP_KPOINTS_FILE]' in line:
#                    continue
  #                  found = True
#                if '[VASP_KPOINTS_MODE_IMPLICIT]' in line:
#                    lines.append('[VASP_KPOINTS_MODE_EXTERNAL]\n')
  #                  found = True
 #           else:
 #               lines.append(line)

 #           if not found:
 #               lines.append(line)

        with open('aflow.in','w') as f:
            f.writelines(lines)



    def can_extract(self):
        import lzma

        outcar = path.join(self.directory,'OUTCAR.static.xz')
        if not path.isfile(outcar):
            return False
        with lzma.open(outcar,'rt') as f:
            if 'free  energy' in f.read():
                return True
            else:
                return False

    def is_executing(self):

        outcar = path.join(self.directory, "OUTCAR")
        outcars = path.isfile(outcar)
        busy = not self.can_extract()
        return outcars and busy

    def is_setup(self):
        aflow = path.join(self.directory, 'aflow.in')
        return path.isfile(aflow)

    def error(self,searchfile,tag):
        vaspout = path.join(self.directory, searchFile)
        if rgrep(vaspout,tag) is not []:
            return True
        else:
            return False

        
    def has_errors(self):
        from glob import glob
        files = sorted(glob(path.join(self.directory, 'aflow.error') + '*') )
        return files is []

    def is_converged(self):
        from aBuild.utility import rgrep
        import lzma
        oszicar = path.join(self.directory, 'OSZICAR.static.xz')
        electronicIt = 0
        with lzma.open(oszicar,'rt') as f:
            for line in f.readlines():
                if 'DAV:' in line:
                    electronicIt = int(line.split()[1])

        if electronicIt < self.INCAR.tags["nelm"]:
            return True
        else:
            return False
        
    def status(self):
        from os import path
        from time import time
        from aBuild.utility import grep
        import os
        fTagStatic = 'aborting loop because EDIFF is reached'
        fTagRelax = ' writing wavefunctions'
        ctime = time()
        #print('checking directory {}'.format(self.directory))


        if self.can_extract():
            if self.is_converged():
                return 'done'
            else:
                return 'unconverged'
        elif self.is_executing():
            return 'running'
        elif self.has_errors():
            return 'errors'
        elif self.is_setup():
            return 'not_started'

    def buildFolder(self):
        self.crystal.write('POSCAR.orig',keepZeros = True)
        self.crystal.write('POSCAR',keepZeros = False)
        self.writeaflowfromposcar('POSCAR')
        if os.path.getsize('aflow.in') < 10:
            return False
        if not self.kpoints_auto_build:
            # Now build the KPOINTs file
            self.KPOINTS.rGP = True
            success = self.KPOINTS.writeKPOINTS()
            print(success, 'success?')
            if not success:
                currentMethod = self.KPOINTS.specs["method"]
                if currentMethod == 'autogr':
                    self.KPOINTS.specs["method"] = 'mueller'
                    self.KPOINTS.specs["mindistance"] = self.KPOINTS.specs["rmin"]
                    self.KPOINTS.specs["includegamma"] = "True"
                    retrySuccess = self.KPOINTS.writeKPOINTS()
                    self.KPOINTS.specs["method"] = currentMethod
                    if not retrySuccess:
                        msg.info("KPOINTS build failed, reverting to MP grid")
                        self.kpoints_auto_build = True
                else:
                    self.KPOINTS.specs["method"] = 'autogr'
                    self.KPOINTS.specs["rmin"] = self.KPOINTS.specs["mindistance"]
                    self.KPOINTS.specs["eps"] = '1e-12'
                    retrySuccess = self.KPOINTS.writeKPOINTS()
                    self.KPOINTS.specs["method"] = currentMethod
                    if not retrySuccess:
                        msg.info("KPOINTS build failed, reverting to MP grid")
                        self.kpoints_auto_build = True
            self.modifyaflowin()
            return True

    def check_atom_counts_zero(self):
        from numpy import array,any
        if any(self.crystal.atom_counts == 0):
            from numpy import  where
            idxKeep = list(where( self.crystal.atom_counts > 0)[0])
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
