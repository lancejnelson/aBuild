from os import path
from aBuild import msg # messaging module
from aBuild.utility import chdir,fileinDir

#from aBuild.calculators import aflow
import os
import sys
import numpy as np
import lzma
config = sys.modules["config"]  

class VASP:
    """Class to handle all of the VASP input and output files.
    Args:
        specs (dict or str):  Either a dictionary containing all of the 
                              necessary settings or a path to a folder that
                              contains all of the files needed.
        systemSpecies (list): List of the species for this system.  Note that this list may not
                              be identical to the species found in the POSCAR or POTCAR.  For example,
                              you may be studying a ternary system, but this particular calculation
                              has 0 of one atom type.  (Optional because if you initialize with a dictionary
                                                        the species list comes in with it)
    """


    def __init__(self, specs,systemSpecies = None, directory = None):

        from aBuild.database.crystal import Crystal

        #Initialize from a dictionary
        if isinstance(specs,dict):
            
            if self._all_present(specs):
                self.POTCAR = POTCAR(specs["potcars"])
                self.KPOINTS = KPOINTS(specs["kpoints"])
                if isinstance(specs["crystal"],Crystal):
                    self.crystal = specs["crystal"]
                else:
                    self.crystal = Crystal(specs["crystal"],specs["species"])
                self.handleSpecialTags(specs)
                    
                self.INCAR = INCAR(specs["incar"])
            else:
                msg.fatal("I don't have all the necessary information to initialize: {}".format(specs.keys()))
        #Initialize from a path
        elif isinstance(specs, str):
            if systemSpecies == None:
                msg.fatal("When initializing a VASP object with a path, you need to supply the system's species list")
            self.directory = specs
            self.getFileNames()
            self.POTCAR = POTCAR(path.join(self.directory,self.potcarName)) if self.potcarName is not None else None
            self.KPOINTS = KPOINTS(path.join(self.directory,self.kpointsName)) if self.kpointsName is not None else None
            self.INCAR = INCAR(path.join(self.directory,self.incarName)) if self.incarName is not None else None
            if self.poscarName is not None:
                if self.potcarName is not None:
                    # Best case scenario: POTCAR is found, systemSpecies has been passed in and there appears to be a POSCAR file
                    self.crystal = Crystal(path.join(self.directory,self.poscarName),systemSpecies,crystalSpecies = self.POTCAR.species)
                else:
                    # If I can't find a POTCAR, I'll let the Crystal object try and determine the species from the POSCAR file
                    self.crystal = Crystal(path.join(self.directory,self.poscarName),systemSpecies)
            else:
                self.crystal = None
        else:
            msg.fatal("Unable to initialize a VASP object from the data that you passed in:", specs)
        if directory is not None:
            self.directory = directory

    def _all_present(self,specs):
        required = ["incar","potcars","kpoints","crystal","species"]
        for tag in required:
            if tag not in specs.keys():
                return False
        return True
    
    def handleSpecialTags(self,specs):
        special = ["AFM","FM"]
        if "FM" in specs.keys():
            specs["incar"]["ispin"] = 2
            specs["incar"]["magmom"] = ''
            
            for idx,species in enumerate(sorted(specs["FM"],reverse = True)):
                specs["incar"]["magmom"] +=  ' '.join(map(str, [ specs["FM"][species] ] * self.crystal.atom_counts[idx]))
                specs["incar"]["magmom"] += ' '
        elif "AFM" in specs.keys():
            if self.crystal.AFMPlanes == None:
                self.crystal.getAFMPlanes([1,0,0])
            if self.crystal.AFMPlanes == None:
                msg.info("You supposedly had an AFM crystal, but I'm not finding the planes")
                return
            specs["incar"]["ispin"] = 2
            #Put in nonzero spin values
            specs["incar"]["magmom"] = ' '.join(map(str,self.crystal.AFMPlanes)) + ' '
            atomsLeft = self.crystal.nAtoms - len(self.crystal.AFMPlanes)
            specs["incar"]["magmom"] += ' '.join(map(str,[0 for x in range(atomsLeft)]))
             
            
    # VASP does not like to have zeros in the atom_counts list
    # but I want to keep track of which atoms are in the crystal.
    # This routine is just here to remove any zeros before I write to
    # the POSCAR file.
    def check_atom_counts_zero(self):
        from numpy import array,any
        if any(self.crystal.atom_counts == 0):
            from numpy import  where
            idxKeep = list(where( self.crystal.atom_counts > 0)[0])
            self.POTCAR.species = list(array(self.POTCAR.species)[idxKeep])
            self.crystal.atom_counts = self.crystal.atom_counts[idxKeep]

            self.crystal.crystalSpecies = self.POTCAR.species
            print(self.crystal.crystalSpecies)
    def _check_tag_exists(self,filename,tag):
        from aBuild.utility import grep
        print('checking that {} tag exists in file {}'.format(tag,filename))
        lines = grep(filename,tag)
        if lines == []:
            print('Not found')
            return False
        else:
            print("Found!")
            return True



    def _check_file_exists(self,file):
        files = os.listdir('./')
        if file in files:
            return True
        else:
            return False


    def getFileNames(self):

        self.incarName = self.fileName('INCAR')
        self.poscarName = self.fileName('POSCAR')
        self.potcarName = self.fileName('POTCAR')
        self.kpointsName = self.fileName('KPOINTS')
        self.vaspOutName = self.fileName('vasp.out')
        self.oszicarName = self.fileName('OSZICAR')
        self.outcarName = self.fileName('OUTCAR')
        self.aflowinName = self.fileName('aflow.in')
        self.aflowendName = self.fileName('aflow.end.out')


    def fileName(self,partName):
        files = fileinDir(partName,self.directory,or_close = True)

        # Let's find the filename of the INCAR we want to initialize
        if True in ['static' in x for x in files]:
            # Must be a static run via AFLOW.  File is INCAR.static.xz probably.
            index = ['static' in x for x in files].index(True)
            if 'xz' in files[index]:
                fileName = partName + '.static.xz'
                filepath = path.join(self.directory,partName + '.static.xz')
            else:
                fileName = partName + '.static'
                filepath = path.join(self.directory,partName+'.static')
        elif True in ['relax' in x for x in files]:
            # Must be a relaxation run via AFLOW
            index = ['relax' in x for x in files].index(True)
            if 'xz' in files[index]:
                fileName = partName + '.relax2.xz'
                filepath = path.join(self.directory,partName+'.relax2.xz')
            else:
                fileName = partName + '.relax2'
                filepath = path.join(self.directory,partName+'.relax2')
        elif files != []:
            fileName = partName
            filepath = path.join(self.directory,partName)
        else:
            fileName = None
            filepath = None
        return fileName
        


    def status(self):
        from os import path
        from time import time
        from aBuild.utility import grep,fileinDir
        import os
        fTagStatic = '------------------------ aborting loop because EDIFF is reached ----------------------------------------\n'
        fTagRelax = ' writing wavefunctions'
        ctime = time()
        print('checking directory {}'.format(self.directory))

        if any([x is None for x in [self.potcarName, self.incarName,self.kpointsName,self.poscarName]]):
            return 'not setup'
        calcType = 'aflow' if True in ['aflow' in x for x in os.listdir(self.directory)] else 'vasp'
        with chdir(self.directory):
#            if calcType == 'aflow':
#                print("Looks like it's and aflow calc")
#                if fileinDir('aflow.end.out','.'):
#                    print("Looks like it's done")
#                    return 'done'
#            outcar = fileinDir('OUTCAR','.', or_close = True)
            outcar = self.outcarName is not None
            incar = self.incarName is not None
            kpoints = self.kpointsName is not None
            potcar = self.potcarName is not None
            poscar = self.poscarName is not None
            output = self.vaspOutName is not None
            oszicar = self.oszicarName is not None
            aflow = self.aflowinName is not None
            aflowend = self.aflowendName is not None
            #if aflowend:
            #    return 'done'
            

            inputs = incar and kpoints and potcar and poscar
            print(inputs,' inputs Found?')
            ''' Check to see if the input files are present
                if they aren't, no need to proceed, just return
                'not setup'
            '''
            
            if not inputs:
                print('inputs not found')
                return 'not setup'

                ''' If the OUTCAR file is present, we know that we're
                    either running, finished successfully, or finished
                    with errors!
                '''
            elif outcar: # OUTCAR present

                sgrcon = grep(self.vaspOutName,'SGRCON')
                tooclose = grep(self.vaspOutName,'HOPE')
                finalenergyline = grep(self.outcarName,'free  energy')
                generalerror = grep(self.vaspOutName,'ERROR')
                # Check to make sure I've converged electonically.
                if grep(self.oszicarName,'DAV:') != []:
                    electronicIteration = int(grep(self.oszicarName,'DAV:')[-1].split()[1])
                else:
                    electronicIteration = 0
                if grep(self.incarName,'nsw') != []:
                    nsw = int(grep(self.incarName,'nsw')[0].split('=')[1])
                    if nsw == 0:
                        nsw = 1
                else:
                    nsw = 1
                if grep(self.oszicarName,'F=') != []:
                    ionicIteration = int(grep(self.oszicarName,'F=')[-1].split()[0])
                else:
                    ionicIteration = 1
                if grep(self.incarName,'nelm') != []:
                    maxelectronic = grep(self.incarName,'nelm')[0].split('=')[1]
                else:
                    maxelectronic = 60
                if ionicIteration == nsw and int(electronicIteration) == int(maxelectronic):
                    return 'unconverged'
                    
                ''' Let's first check to see if this is a static
                calculation or a relaxation because the tag 
                to check for is different.'''
                if incar:
                    relax = grep(self.incarName,'IBRION')
                    if '-1' not in relax or relax is []:
                            static = True
                    else:
                            static = False
                else:
                    return 'not setup'
                    
                ''' Check finish tag for static calc'''
                if static and self._check_tag_exists(self.outcarName, fTagStatic):  #finish tag found
                    if finalenergyline != []:  #Let's double check
                        return 'done'
                    else:  # Apparently not,  why?
                        return 'idk'
                    
                    ''' Check finish tag for relax calc'''
                elif self._check_tag_exists(self.outcarName,fTagRelax): #Looks like it's done
                    if finalenergyline != []:  # Let's double check
                        return 'done'
                    else:  # Apparently not, why?
                        return 'idk'
                else:
                        
                    ''' Check how long since the last file write.  If it was recent
                     then we're probably running.'''
                    time = path.getmtime(self.outcarName)
                    if (ctime - time) < 3600:  # If the OUTCAR was modified in the last hour
                                              # the calculation is probably still running.
                        return 'running'
                    elif sgrcon:
                        return 'sgrcon'
                    elif generalerror:
                        return 'error'
                    elif tooclose:
                        return 'warning'
                    else:
                        return 'killed before done'
                    
            else:
                    return 'not started'
                    
                        
            if output:
                    warning = grep(self.vaspOutName,'RRRRR') != [] or grep(self.vaspOutName,'AAAAAA') != []
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
        print(self.POTCAR.species, 'check here')
        from numpy import  where
        from aBuild.calculators.vasp import POSCAR
        self.KPOINTS.rGP = runGetKPoints
        self.INCAR.writeINCAR()
        print("INCAR built")
        print(self.crystal.crystalSpecies,' species')
        self.crystal.write('POSCAR.orig',keepZeros = True)
        print("POSCAR_orig built")
        #self.check_atom_counts_zero()  # This routine modifies POTCAR.species and crystal.crystalSpecies. Bad idea...?
        self.crystal.write('POSCAR',keepZeros = False)
        print("POSCAR built")
        self.KPOINTS.writeKPOINTS()
        print("KPOINTS built")
        idxKeep = list(where( self.crystal.atom_counts > 0)[0])

        self.POTCAR.writePOTCAR(indices = idxKeep)
        #print("POTCAR built")

        

    def read_forces(self,allIonic = True):
        import lzma
        openDict = {True: lzma.open(self.poscarName,'rt'), False: open(self.poscarName,'r')}

        with openDict['xz' in self.poscarName] as f:
            poslines = f.readlines()
                
        if any(c.isalpha() for c in poslines[5].strip()):  #It's a CONTCAR
            nAtoms = sum([int(i) for i in poslines[6].split()])
        else:

            nAtoms = sum([int(i) for i in poslines[5].split()])

        if 'xz' in self.outcarName:
            import lzma
            with lzma.open(self.outcarName,'rt') as f:
                lines = f.readlines()
        else:
            with open(self.outcarName, 'r') as f:
                lines = f.readlines()

        n = 0

        if allIonic:
            forces = []

        n = 0
        found = False
        for line in lines:
            if line.rfind('TOTAL-FORCE') > -1:
                found = True
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
        if not found:
            msg.info("Couldn't find forces for this calc")
            return None
        if not allIonic:
            forces = singleItForces
        if allIonic and len(forces) == 1:
            return forces[0]
        
        return forces

    def read_fermi(self):
        
        """Method that reads Fermi energy from OUTCAR file"""
        E_f = None

        openDict = {True: lzma.open(self.outcarName,'rt'), False: open(self.outcarName,'r')}
        for line in openDict['xz' in self.outcarName]:
            if line.rfind('E-fermi') > -1:
                E_f = float(line.split()[2])
        return E_f
                
    def read_nbands(self):
        
        openDict = {True: lzma.open(self.outcarName,'rt'), False: open(self.outcarName,'r')}
        for line in openDict['xz' in self.outcarName]:
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
        openDict = {True: lzma.open(self.outcarName,'rt'), False: open(self.outcarName,'r')}
        for line in openDict['xz' in self.outcarName]:
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
        openDict = {True: lzma.open(self.outcarName,'rt'), False: open(self.outcarName,'r')}
        for line in openDict['xz' in self.outcarName]:
            if line.find(' Total    ') != -1:
                try:
                    stress = -np.array([float(a) for a in line.split()[1:]])
                    stress = stress[[0, 1, 2, 4, 5, 3]] #* 1e-1 #* .00624151# * ase.units.GPa.. Gets me to Giga Pascals
                except: 
                    msg.error("Unable to read stress from OUTCAR file")
                    stress = None
        return stress

    def read_results(self, allElectronic = False, allIonic=False,mustHave = ['free energy','stress','species','formation energy']):
        if self.directory is not None and self.status() == 'done':
            with chdir(self.directory):
                self.crystal.results = {}
                self.crystal.results["warning"] = False
                self.crystal.results["energyF"],self.crystal.results["energyZ"] = self.read_energy(allElectronic=allElectronic)
                self.crystal.results["forces"] = self.read_forces(allIonic=allIonic)
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
                msg.info("Succesfully read in results")
        else:
            self.crystal.results = None
            msg.info("Unable to extract necessary information from directory! ({}).  Status: {} ".format(self.directory,self.status))
            
    def add_to_results(self,key,item):
        if self.crystal.results is None:
            self.crystal.results = {}
        self.crystal.results[key] = item

    @property
    def formationEnergy(self):
        from os import path
        pures = []
        for i in range(self.crystal.nTypes):
            pureDir = path.join(path.split(self.directory)[0], 'pure' + self.crystal.systemSpecies[i])
            if path.isdir(pureDir):
                pureVASP = VASP(pureDir,systemSpecies = self.crystal.systemSpecies)
                pureVASP.read_results()
                pures.append(pureVASP)
            else:
                msg.warn("Unable to read pures information, not calculating formation energy")
                return None

        try:
            formationEnergy = self.crystal.results["energyF"]/self.crystal.nAtoms - sum(   [ pures[i].crystal.results["energyF"]/pures[i].crystal.nAtoms * self.crystal.concentrations[i] for i in range(self.crystal.nTypes)])
        except:
            formationEnergy = 10000
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
            if 'build' in specs:
                self.build = specs["build"]
            else:
                self.build = 'manual'
        elif isinstance(specs,str):
            self._init_path(specs)



    def _init_path(self,filepath,fileformat = 'POTCAR'):
        import os
        from os import path
        from aBuild.utility import fileinDir
              
        if not path.isfile(filepath):
            self.species = None
            self.setups = None
            self.version = None
            self.xc = None
            self.directory = None
            return
        if 'xz' in filepath:
            import lzma
            with lzma.open(filepath,'rt') as f:
                lines = f.readlines()
        else:
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
            print('Found POTCAR files')
            potsGood = True
            for pot in pots:
                with open(pot, 'r') as f:
                    lines = f.readlines()
                if lines[0].split()[-1] != self.versions[lines[0].split()[-2].split('_')[0]]:
                    potsGood = False
            
            return potsGood
        else:
            print("can't find POTCAR files")
            return False
        
    def writePOTCAR(self,indices=None,filename = 'POTCAR'):
        from numpy import array
        if self.build == 'aflow':
            return
        if not self._potcarsOK():
            ermsg = "Can't find the POTCARS you specified: {}".format(self.versions)
            msg.fatal(ermsg)
        if indices is not None:
            srcpaths = [path.join(self.srcdirectory,x + self.setups[x],'POTCAR') for x in array(self.species)[indices]]
        else:
            srcpaths = [path.join(self.srcdirectory,x + self.setups[x],'POTCAR') for x in self.species]
        from os import waitpid
        from subprocess import Popen

        command = "cat {} >  {}   ".format(' '.join(srcpaths), filename)
        child=Popen(command, shell=True, executable="/bin/bash")
        waitpid(child.pid, 0)

        
class INCAR:


    def __init__(self,specs):

        if isinstance(specs,dict):
            if 'defaults' in specs and specs["defaults"]:
                self.setDefaultTags()
            else:
                self.tags = specs
        elif isinstance(specs,str):
            self._init_path(specs)
        else:
            self.setDefaultTags()




        
    def _init_path(self,filepath,filename = 'INCAR'):
        if not path.isfile(filepath):
            self.tags = None
            return

        if 'xz' in filepath:
            import lzma
            with lzma.open(filepath,'rt') as f:
                lines = f.readlines()
        else:
            with open(filepath,'r') as f:
                lines = f.readlines()
        self.tags = {}

        for line in lines:
            if line[0] != '#':
                self.tags[line.split('=')[0]] = line.split('=')[1]


    def setDefaultTags(self):
        self.tags = {}
        self.tags["prec"] = "High"
        self.tags["sigma"] = "0.1"
        self.tags["ISMEAR"] = "1" 
        self.tags["ISYM"] = "2"  # Symmetry on
        self.tags["ALGO"] = "Normal"
        self.tags["LWAVE"] = ".FALSE."
        self.tags["LREAL"] = "auto"
        self.tags["LORBIT"] = "10"
        
        #    def processTags(self,tags):
        #for key,val in tags.items():
        ##    self.tags
        #pass
        
    def writeINCAR(self,filename='INCAR'):
        if 'build' in self.tags.keys() and self.tags["build"] == 'aflow':
            return
        lines = []
        for key,val in self.tags.items():
            lines.append(''.join([key," = ",str(val),'\n']))

        with open(filename,'w') as f:
            f.writelines(lines)


class KPOINTS:

    def __init__(self,specs):


        if isinstance(specs,dict):
            self.specs = specs
            #self.method = specs["method"]
            #self.density = specs["mindistance"]
            #self.includeGamma = True
        elif isinstance(specs,str):
            self._init_file(specs)
            

    def _init_file(self,filepath,filename = 'KPOINTS'):
        from os import path

        if not path.isfile(filepath):
            self.method = None
            self.includeGamma = None
            return
        if 'xz' in filepath:
            import lzma
            with lzma.open(filepath,'rt') as f:
                lines = f.readlines()
        else:
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
        
        methodlookup = {'mueller': self.mueller,'equivalent': self.equivalent, "mp": self.monkPack,"autogr": self.autogr}
        print(self.specs["method"],' check here')

        if self.specs["method"] not in ['mueller','equivalent','mp','autogr']:
            msg.error("I don't recognize the method you have specified: {}".format(self.specs["method"]))
            
        return methodlookup[self.specs["method"]](filename = filename)


    def autogr(self,filename='KPOINTS'):
        from os import waitpid,path
        from subprocess import Popen
        self.KPGEN()
        

        if config.AUTOGR is not None:
            print('Found AUTOGR executable')
            if self.rGP:
                command = "{}".format(config.AUTOGR)
                child=Popen(command, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)
            else:
                msg.info("Not running the getKpoints script")
        else:
            msg.fatal("You haven't defined the environment variable: AUTOGR, so I don't know how to generate KPOINT grids ")

        if not path.isfile('KPOINTS'):
            print("Can't find KPOINTS")
            return False
        else:
            print("KPOINTS file FOUND")
            return True


    def mueller(self,filename='KPOINTS',runGetKpoints = False):
        from os import waitpid,path
        from subprocess import Popen

        self.PRECALC()
        

        if config.GETKPTS is not None:
            print('Found GETKPTS executable')
            if self.rGP:
                command = "{}".format(config.GETKPTS)
                child=Popen(command, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)
            else:
                msg.info("Not running the getKpoints script")
        else:
            msg.fatal("You haven't defined the environment variable: GETKPTS, so I don't know how to generate KPOINT grids ")
        if not path.isfile('KPOINTS'):
            return False
        else:
            return True


    def KPGEN(self):

        allowed = ['RMIN','KPDENSITY','KSPACING','KPPRA','SHIFT']
        lines = []
        for opt in self.specs.keys():
            if opt.upper() in allowed:
                lines.append('{}={}\n'.format(opt.upper(),self.specs[opt]))
        #lines.append('SHIFT={}\n'.format(self.includeGamma))
        #lines.append('RMIN={}\n'.format(self.density))
        with open('KPGEN','w') as f:
            f.writelines(lines)

    def PRECALC(self):

        lines = []
        for opt in self.specs.keys():
            lines.append('{}={}\n'.format(opt,self.specs[opt]))
#        lines.append('INCLUDEGAMMA={}\n'.format(self.includeGamma))
#        lines.append('MINDISTANCE={}\n'.format(self.density))
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
        if 'xz' in filename:
            import lzma
            with lzma.open(os.path.abspath(filename),'rt') as f:
                poscarlines = f.readlines()

        else:
            with open(os.path.abspath(filename),'r') as f:
                poscarlines = f.readlines()

        self.from_string(poscarlines)

    def from_string(self, poscarlines):
        """Initializes the POSCAR lines object from a list of strings.

        :arg poscarlines: a list of strings from the POSCAR file.
        """
        self.label = poscarlines[0].strip().split('\n')[0]
        self.latpar = poscarlines[1]
        self.Lv = poscarlines[2:5]

        #CONTCARs have species names in the next line, but typical POSCARs dont.  Let's
        # Figure out which one we have

        if any(c.isalpha() for c in poscarlines[5].strip()):  #Updated styling for
            countsLine = 6
            coordSysLine = 7
            basisStartLine = 8
            self.atom_counts = poscarlines[countsLine].strip()
            nBas = sum(map(int, self.atom_counts.split()))
            self.coordsys = poscarlines[coordSysLine].split('\n')[0]
            self.Bv = [' '.join(x.split()[:3]) for x in poscarlines[basisStartLine:basisStartLine+nBas]]
            #self.Bv = poscarlines[basisStartLine:basisStartLine+nBas]
            self.species = poscarlines[5].split()
            
        else: #It's the older styling
            countsLine = 5
            coordSysLine = 6
            basisStartLine = 7
            self.atom_counts = poscarlines[countsLine].strip()
            nBas = sum(map(int, self.atom_counts.split()))
            self.coordsys = poscarlines[coordSysLine].split('\n')[0]
            self.Bv = [' '.join(x.split()[:3]) for x in poscarlines[basisStartLine:basisStartLine+nBas]]
            if any(c.isalpha() for c in poscarlines[basisStartLine].strip()):
                self.species = []
                for i in poscarlines[basisStartLine:basisStartLine + nBas]:
                    if i.split()[-1] not in self.species:
                        self.species.append(i.split()[-1])
                
            else:
                self.species = None
        
        if basisStartLine + 2*nBas < len(poscarlines):
            #We could still have concentration information in the POSCAR
            self.concentrations = poscarlines[basisStartLine+nBas:basisStartLine+2*nBas]
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




    
