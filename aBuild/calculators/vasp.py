from os import path
from aBuild import msg # messaging module
from aBuild.utility import chdir
import os
import sys
import numpy as np

config = sys.modules["config"]  

class VASP:
    """Class to handle all of the VASP input and output files.
    Args:
        root (str): Path to the calculation folder
        incar (dict): Dictionary containing the INCAR tags to be used
        potcars (dict): Dictionary containing the necessary settings to 
                        find the correct POTCARS.
                        <directory> : where the potcars are located
                        <>
        kpoints (dict): KPOINTS settings
        crystal (CRYSTAL OBJ): Crystal description
    """


    def __init__(self, specs,systemSpecies,directory = None):

        from aBuild.database.crystal import Crystal
        
        if isinstance(specs,dict):
            self.INCAR = INCAR(specs["incar"])
            self.POTCAR = POTCAR(specs["potcar"])
            self.KPOINTS = KPOINTS(specs["kpoints"])
            if isinstance(specs["crystal"],Crystal):
                self.crystal = specs["crystal"]
            else:
                self.crystal = Crystal(specs["crystal"],systemSpecies)
        elif isinstance(specs, str):
            self.POTCAR = POTCAR(path.join(specs,'POTCAR'))
            self.KPOINTS = KPOINTS(path.join(specs,'KPOINTS'))
            self.crystal = Crystal(path.join(specs,'POSCAR'),systemSpecies,crystalSpecies = self.POTCAR.species)
            self.directory = specs
        else:
            msg.fatal("Unable to initialize a VASP object from the data that you passed in:", specs)
        if directory is not None:
            self.directory = directory



        
    def check_atom_counts_zero(self):
        from numpy import array

        if any(self.crystal.atom_counts == 0):
            from numpy import  where
            idxKeep = list(where( self.crystal.atom_counts > 0)[0])
            self.POTCAR.species = list(array(self.POTCAR.species)[idxKeep])
            self.crystal.atom_counts = self.crystal.atom_counts[idxKeep]
            
    def _check_tag_exists(self,file,tag):
        from aBuild.utility import grep

        lines = grep(file,tag)
        if lines == []:
            return False
        else:
            return True



    def _check_file_exists(self,file):
        files = os.listdir('./')
        if file in files:
            return True
        else:
            return False

    def status(self):
        from os import path
        from time import time
        from aBuild.utility import grep
        import os

        ctime = time()
        with chdir(self.directory):
            outcar = self._check_file_exists('OUTCAR')
            incar = self._check_file_exists('INCAR')
            kpoints = self._check_file_exists('KPOINTS')
            potcar = self._check_file_exists('POTCAR')
            poscar = self._check_file_exists('POSCAR')
            output = self._check_file_exists('output')


            inputs = incar and kpoints and potcar and poscar

            ''' Check to see if the input files are present
                if they aren't, no need to proceed, just return
                'not setup'
            '''
            if not inputs:
                return 'not setup'

                ''' If the OUTCAR file is present, we know that we're
                    either running, finished successfully, or finished
                    with errors!
                '''
            elif outcar:
                time = path.getmtime('OUTCAR')
                sgrcon = grep('OUTCAR','SGRCON')
                finalenergyline = grep('OUTCAR','free  energy')

                ''' Check how long since the last file write.  If it was recent
                     then we're probably running.'''
                if (ctime - time) < 60:
                    return 'running'

                    ''' If it's been a while since the last write, we've probably finished
                        the calculation.  Let's proceed to error/warning check. 
                    '''
                else:

                    ''' Let's first check to see if this is a static
                        calculation or a relaxation because the tag 
                        to check for is different.'''
                    if incar:
                        relax = grep('INCAR','IBRION')
                        if '-1' not in relax or relax is []:
                            static = True
                        else:
                            static = False

                    ''' Check finish tag for static calc'''
                    if static and self._check_tag_exists('OUTCAR',
                                        '------------------------ aborting loop because EDIFF is reached ----------------------------------------\n'):
                        folderstat = 'done'
                        ''' Check finish tag for relax calc'''
                    elif self._check_tag_exists('OUTCAR',' writing wavefunctions'):
                        folderstat = 'done'
                    elif finalenergyline == []:
                        return 'error'
                    elif abs( float(finalenergyline[0].split()[-2]) ) > 1000:
                        return 'error'
                    else:
                        return 'idk'

                    if finalenergyline == []:
                        return 'error'
                    if finalenergyline != []  and abs( float(finalenergyline[0].split()[-2]) ) > 1000:
                        return 'error'
            else:
                    folderstat = 'not started'
                    
                        
            if output:
                    warning = grep('output','RRRRR') != [] or grep('output','AAAAAA') != []
                    if warning:
                        return 'warning'



                
        return folderstat
    
    @staticmethod
    def from_file(runpath):

        incar = INCAR.from_file(runpath)
        kpoint = KPOINTS.from_file(runpath)
        potcars = POTCARS.from_file(runpath)
        result = VASP(runpath = runpath,incar=incar,kpoints=kpoints,potcars=potcars)
        return result

    def buildFolder(self,runGetKPoints = True):
        from aBuild.calculators.vasp import POSCAR

        self.KPOINTS.rGP = runGetKPoints
        self.INCAR.writeINCAR()
        print("INCAR built")
        self.crystal.write('POSCAR_orig')
        print("POSCAR_orig built")
        self.check_atom_counts_zero()
        self.crystal.write('POSCAR')
        print("POSCAR built")
        self.KPOINTS.writeKPOINTS()
        
        print("KPOINTS built")
        self.POTCAR.writePOTCAR()
        #print("POTCAR built")

    def read_forces(self,allIonic = True):

        with open('POSCAR','r') as file:
            poslines = file.readlines()
        nAtoms = sum([int(i) for i in poslines[5].split()])
        
        with open('OUTCAR', 'r') as file:
            lines = file.readlines()

        n = 0

        if allIonic:
            forces = []

        n = 0
        for line in lines:
            if line.rfind('TOTAL-FORCE') > -1:
                singleItForces = []
                for i in range(nAtoms):
                    singleItForces.append(np.array([float(f) for f in
                                                lines[n + i + 2].split()[3:6]]))
                msg.info('Found forces for {} atoms.'.format(nAtoms))
                if not '--' in lines[n+nAtoms + 2]:
                    msg.fatal('It appears that there are forces for more atoms than I was expecting!')
                if allIonic:
                    forces.append(singleItForces)
            n+=1
        if not allIonic:
            forces = singleItForces
        if allIonic and len(forces) == 1:
            return forces[0]
        
        return forces

    def read_fermi(self):
        """Method that reads Fermi energy from OUTCAR file"""
        E_f = None
        for line in open('OUTCAR', 'r'):
            if line.rfind('E-fermi') > -1:
                E_f = float(line.split()[2])
        return E_f
                
    def read_nbands(self):
        for line in open('OUTCAR', 'r'):
            line = self.strip_warnings(line)
            if line.rfind('NBANDS') > -1:
                nBands  = int(line.split()[-1])
        return nBands

    def read_energy(self, allElectronic=False):
        energyZ = None
        energyF = None
        if allElectronic:
            energyF = []
            energyZ = []
        for line in open('OUTCAR', 'r'):
            # Free energy
            if line.lower().startswith('  free  energy   toten') or line.lower().startswith('  free energy    toten'):
                if allElectronic:
                    energyF.append(float(line.split()[-2]))
                else:
                    energyF = float(line.split()[-2])
                    # Extrapolated zero point energy
            if line.startswith('  energy  without entropy'):
                if allElectronic:
                    energyZ.append(float(line.split()[-1]))
                else:
                    energyZ = float(line.split()[-1])
                        
        return energyF,energyZ

    def read_stress(self):
        stress = None
        for line in open('OUTCAR'):
            if line.find(' in kB  ') != -1:
                stress = -np.array([float(a) for a in line.split()[2:]])
                stress = stress[[0, 1, 2, 4, 5, 3]] * 1e-1 * .00624151# * ase.units.GPa.. Gets me to Giga Pascals
        return stress

    def read_results(self, allElectronic = False, allIonic=False):
        if self.directory is not None and self.status() is 'done':
            with chdir(self.directory):
                self.crystal.results = {}
                self.crystal.results["warning"] = False
                self.crystal.results["energyF"],self.crystal.results["energyZ"] = self.read_energy(allElectronic=allElectronic)
                self.crystal.results["forces"] = self.read_forces(allIonic=allIonic)
                self.crystal.results["stress"] = self.read_stress()
                #self.POTCAR = POTCAR.from_POTCAR()
                self.crystal.results["species"] = self.POTCAR.species
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
    
    def add_to_results(self,key,item):
        if self.crystal.results is None:
            self.crystal.results = {}
        self.crystal.results[key] = item

    @property
    def formationEnergy(self):
        pures = []
        for i in range(self.crystal.nTypes):
            pureDir = path.join(path.split(self.directory)[0], 'pure' + self.crystal.species[i])
            pureVASP = VASP(pureDir,self.crystal.species)
            pureVASP.read_results()
            pures.append(pureVASP)

        
        formationEnergy = self.crystal.results["energyF"]/self.crystal.nAtoms - sum(   [ pures[i].crystal.results["energyF"]/pures[i].crystal.nAtoms * self.crystal.concentrations[i] for i in range(self.crystal.nTypes)])
        return formationEnergy
        
class POTCAR:

    def __init__(self,specs):

        if isinstance(specs,dict):
            self.srcdirectory = specs["directory"]
            self.xc = specs["xc"]
            self.versions = specs["versions"]
            self.species = list(specs["versions"].keys())
            self.species.sort(reverse=True)
            if sorted(self.species,reverse = True) != self.species:
                msg.fatal('Species are not in reverse alphabetical order... Problem?')
            self.setups = specs["setups"]
        elif isinstance(specs,str):
            self._init_path(specs)



    def _init_path(self,filepath,fileformat = 'POTCAR'):
        import os
        from os import path

              
        if not path.isfile(filepath):
            self.species = None
            self.setups = None
            self.version = None
            self.xc = None
            self.directory = None
            return
        with open(filepath,'r') as f:
            lines = f.readlines()
            
        species = []
        setups = {}
        versions = {}
        xc = []
        for line in lines:
            if 'TITEL' in line.split():
                species.append(line.split()[3].split('_')[0])
                try:
                    setups[line.split()[3].split('_')[0]] = '_' + line.split()[3].split('_')[1]
                except:
                    setups[line.split()[3].split('_')[0]] = ''
                versions[line.split()[3].split('_')[0]] = line.split()[-1]
                xc.append(line.split()[2])

        self.species = species
        self.setups = setups
        self.version = versions
        self.xc = xc
        self.directory = None
        
    # Checks to see that the requested POTCARS are found.  If they are not found, stop and warn the user.
    def _potcarsOK(self):
        from os.path import isfile
        self.species.sort()
        self.species.reverse()
        pots = [path.join(self.srcdirectory,x + self.setups[x],'POTCAR') for x in self.species]
        if all([isfile(x) for x in pots]):
            print('found files')
            potsGood = True
            for pot in pots:
                with open(pot, 'r') as f:
                    lines = f.readlines()
                if lines[0].split()[-1] != self.versions[lines[0].split()[-2].split('_')[0]]:
                    potsGood = False
            
            return potsGood
        else:
            print("can't find files")
            return False
        
    def writePOTCAR(self,filename = 'POTCAR'):

        if not self._potcarsOK():
            ermsg = "Can't find the POTCARS you specified: {}".format(self.versions)
            msg.fatal(ermsg)
        srcpaths = [path.join(self.srcdirectory,x + self.setups[x],'POTCAR') for x in self.species]

        from os import waitpid
        from subprocess import Popen

        command = "cat {} >  {}   ".format(' '.join(srcpaths), filename)
        child=Popen(command, shell=True, executable="/bin/bash")
        waitpid(child.pid, 0)

        
class INCAR:


    def __init__(self,specs):

        if isinstance(specs,dict):
            self.tags = specs
        elif isinstance(specs,str):
            self._init_file(specs)
        else:
            self.setDefaultTags()


        
    def _init_file(path,filename = 'INCAR'):
        with open(path,'r') as f:
            lines = f.readlines()

        self.tags = {}

        for line in lines:
            self.tags[line.split('=')[0]] = line.split('=')[1]


    def setDefaultTags(self):
        self.tags = {}
        self.tags["prec"] = "a"
        self.tags["sigma"] = "0.1"
        self.tags["GGA"] = "PE"
        self.tags["ISMEAR"] = "1"
        self.tags["LWAVE"] = ".FALSE."
        self.tags["LREAL"] = "auto"
        
        #    def processTags(self,tags):
        #for key,val in tags.items():
        ##    self.tags
        #pass
        
    def writeINCAR(self,filename='INCAR'):

        lines = []
        for key,val in self.tags.items():
            lines.append(''.join([key," = ",str(val),'\n']))

        with open(filename,'w') as f:
            f.writelines(lines)


class KPOINTS:

    def __init__(self,specs):


        if isinstance(specs,dict):
            self.method = specs["method"]
            self.density = specs["mindistance"]
            self.includeGamma = True
        elif isinstance(specs,str):
            self._init_file(specs)
            

    def _init_file(self,filepath,filename = 'KPOINTS'):
        from os import path

        if not path.isfile(filepath):
            self.method = None
            self.includeGamma = None
            return
        with open(filepath,'r') as f:
            lines = f.readlines()
        if 'Server' in lines[0].split():
            self.method = 'mueller'
            densityLocation = lines[0].split().index('Angstroms.') - 1
            self.density = float(lines[0].split()[densityLocation])
            if 'Grid includes gamma point' in lines[0]:
                self.includeGamma = True
            else:
                self.includeGamma = False
                #elif 'Automatically generated mesh' in lines[0]:
                #self.method = 'MP'
            
            #            self.
    def writeKPOINTS(self,filename='KPOINTS'):
        
        methodlookup = {'mueller': self.mueller,'equivalent': self.equivalent, "mp": self.monkPack}

        if self.method not in ['mueller','equivalent','mp']:
            mssg.error("I don't recognize the method you have specified: {}".format(self.method))
            
        methodlookup[self.method](filename = filename)



    def mueller(self,filename='KPOINTS',runGetKpoints = False):
        from os import waitpid
        from subprocess import Popen

        self.PRECALC()
        

        if config.GETKPTS is not None:
            print('found GETKPTS')
            if self.rGP:
                command = "{}".format(config.GETKPTS)
                child=Popen(command, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)
            else:
                msg.info("Not running the getKpoints script")
        else:
            msg.fatal("You haven't defined the environment variable: GETKPTS, so I don't know how to generate KPOINT grids ")


    def PRECALC(self):

        lines = []
        lines.append('INCLUDEGAMMA={}\n'.format(self.includeGamma))
        lines.append('MINDISTANCE={}\n'.format(self.density))
        with open('PRECALC','w') as f:
            f.writelines(lines)


    def equivalent(self):
        pass

    def monkPack(self,kvecs,nAtoms): 
        from math import sqrt
        from numpy import dot,array,round,prod,argmin

        magnitudes = [sqrt(dot(k, k)) for k in kvecs]
        sortmags = sorted(magnitudes)
        largest = sortmags[-1]
        relativediff = [magnitudes[x]/largest for x in range(3)]
        possibleDivisions = [round(array(relativediff) * i)  for i in range(1,20) if not any(round(array(relativediff) * i) == 0) ]
        densities = array([prod(i) * nAtoms for i in possibleDivisions])
        uniformMeasure = [sum(abs(round(array(relativediff) * i) - array(relativediff) * i)) for i in range(1,20)]

        locMin = argmin(uniformMeasure)

        bestDensity = argmin(abs(densities- float(self.density)))
        if bestDensity != 0:
            densityChoices = [bestDensity - 1, bestDensity, bestDensity + 1]
            uniformityChoices = uniformMeasure[bestDensity -1: bestDensity + 2]
        else:
            densityChoices = [ bestDensity, bestDensity + 1]
            uniformityChoices = uniformMeasure[bestDensity: bestDensity + 2]

        bestChoice = densityChoices[argmin(uniformityChoices)]
        #        thisOne = densityChoices[argmin(uniformMeasure[argmin(abs(densities- float(targetDensity))) -1: argmin(abs(densities- float(targetDensity))) + 2])]
        return list(map(int,possibleDivisions[bestChoice])) + [0,0,0]

class POSCAR(object):
    """Represents the POSCAR text representation of a crystal. Useful so that
    code can refer to concrete names like 'Lv' or 'Bv' and get their text
    representations instead of referring to some obscure index in the list of
    strings that form the file."""
    def __init__(self, crystal):
        from aBuild.database.crystal import Crystal,Lattice
        """Create the text representation of the crystal in POSCAR format."""
        #The label is the first line in the POSCAR file. ancle uses some
        #naming standards to help organize things.
        if isinstance(crystal, Crystal):
            self._init_crystal(crystal)
        elif isinstance(crystal, Lattice):
            self._init_lattice(crystal)
        elif isinstance(crystal, list):
            self.from_string(crystal)
        else:
            self._init_file(crystal)

    def __str__(self):
        """Returns the string representation of the POSCAR lines to write
        to a file."""
        return '\n'.join(self.lines())
        
    def todict(self):
        """Returns a dictionary containing the POSCAR information with keys
        'Lv', 'basis'=self.atom_counts, 'Bv', 'latpar' and 'coordsys'.
        """
        poscar = {}
        from numpy import array
        poscar["Lv"] = array([map(float, l.split()[0:3]) for l in self.Lv])
        poscar["basis"] = map(int, self.atom_counts.split())
        poscar["Bv"] = array([map(float, l.split()[0:3]) for l in self.Bv])

        #The latpar is sometimes missing and we have a "scale factor" or
        #some other text their instead. If it is missing, just set it to 0.
        try:
            poscar["latpar"] = float(self.latpar)
        except ValueError:
            poscar["latpar"] = 0
        poscar["coordsys"] = self.coordsys
        return poscar

        
    def lines(self, vasp=False):
        """Returns a list of strings for each line in the POSCAR file.

        :arg vasp: when true, the atom_counts line is checked for zeros before
          it is created. Vasp can't handle zero for the number of atoms of a
          certain type; just remove it."""
        result = []
        result.append(self.label)
        result.append(self.latpar)
        result.append(self.Lv_lines)
        result.append(self.atom_counts)
        #        if vasp:
        #    result.append(' '.join([a for a in self.atom_counts if a != '0' and a != ' ']))
        #else:
        #    result.append(' '.join([a for a in self.atom_counts if a != ' ']))
        result.append(self.coordsys)
        result.append(self.Bv_lines)
        return result

    def write(self, filename='POSCAR', vasp=False):
        """Writes the contents of this POSCAR to the specified file."""
        #fullpath = os.path.abspath(filepath)
        with open(filename, 'w') as f:
            if not vasp:
                f.write(self.__str__())
            else:
                f.write('\n'.join(self.lines(vasp)))

    @property
    def Lv_lines(self):
        """Return \n joined lattice vector text lines."""
        return '\n'.join(self.Lv)

    @property
    def Bv_lines(self):
        """Return \n joined basis vector text lines."""
        return '\n'.join(self.Bv)

    def _init_file(self, filename):
        """Initializes the POSCAR lines from a file."""
        with open(os.path.abspath(filename)) as f:
            poscarlines = f.readlines()

        self.from_string(poscarlines)

    def from_string(self, poscarlines):
        """Initializes the POSCAR lines object from a list of strings.

        :arg poscarlines: a list of strings from the POSCAR file.
        """
        self.label = poscarlines[0].strip().split('\n')[0]
        self.latpar = poscarlines[1]
        self.Lv = poscarlines[2:5]
        self.atom_counts = poscarlines[5].strip()
        self.coordsys = poscarlines[6].split('\n')[0]

        nBas = sum(map(int, self.atom_counts.split()))
        self.Bv = poscarlines[7:7+nBas]
        
        if 7 + 2*nBas < len(poscarlines):
            #We could still have concentration information in the POSCAR
            self.concentrations = poscarlines[7+nBas:7+2*nBas]
        else:
            self.concentrations = ""

    def _init_lattice(self, lattice):
        """Initializes the POSCAR lines from a Lattice instance."""
        self.label = "Lattice PosCar"
        self.latpar = str(lattice.latpar)
        self.Lv = ['  '.join([str(i) for i in L]) for L in lattice.Lv ]
        self.atom_counts = ' '.join([str(c) for c in lattice.atom_counts])
        self.coordsys = lattice.coordsys
        self.Bv = ['  '.join([str(i) for i in L]) for L in lattice.Bv ]

    def _init_crystal(self, crystal):
        """Initializes the POSCAR lines from a Crystal object."""
        if 'pure' in str(crystal.title):
            self.label = 'Pure PosCar {}'.format(crystal.strN.split("pure")[1])
        else:
            self.label = crystal.title

        self.latpar = str(crystal.latpar)
        self.Lv = ['  '.join([str(i) for i in L]) for L in crystal.lattice ]
        self.atom_counts = ' '.join([str(c) for c in crystal.atom_counts])
        self.coordsys = crystal.coordsys
        self.Bv = ['  '.join([str(i) for i in L]) for L in crystal.basis ]




    
