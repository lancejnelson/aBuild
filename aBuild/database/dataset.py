
from aBuild import msg
from aBuild.utility import chdir, _get_reporoot

class dataset:

    def __init__(self,dset,systemSpecies,root=None,calculator = None,lFormat = 'mtpselect'):
        from os import path,makedirs
        from aBuild.database.crystal import Crystal
        from aBuild.calculators.vasp import VASP

        self.calculator = calculator
        self.species = systemSpecies
        
        if isinstance(dset,list):
            
            if isinstance(dset[0], dict):
                self.init_enum(dset,systemSpecies)
            elif isinstance(dset[0], Crystal):
                self.crystals = dset
                self.nCrystals = len(dset)
            elif isinstance(dset[0], VASP):
                self.calcs = dset
                self.nCalcs = len(dset)
            elif isinstance(dset[0],str):
                self.init_paths(dset,systemSpecies)
        elif isinstance(dset, str):
           self.init_file(dset,lFormat)

        self.root = root

    # Used to be called 'buildFoldersFromEnum
    def init_enum(self,enumdicts,systemSpecies,runGetKpoints = True):
        from aBuild.enumeration import Enumerate
        from aBuild.calculators.vasp import VASP
        from aBuild.database.crystal import Crystal
        from aBuild.jobs import Job
        from random import randrange
        from aBuild.utility import chdir
        from os import remove, path

        #        from crystal import Crystal
        from os import path
        import os

        #    if not path.isdir(self.root):
        #    os.mkdir(self.root)
        print("Building database from enumerations")
        self.crystals = []
        #        configIndex = startPoint = self._starting_point
        for eDict in enumdicts:
            enumController = Enumerate(eDict)
            if enumController.nEnumStructs == 0:
                msg.warn('There are no enumerated structures for lattice type {}.  Not building any VASP folders for them.'.format(self.enumDicts[index]["lattice"]))
                enumController.buildInputFile()

                enumController.enumerate()

            # Loop to generate random structures for a given lattice type
            for i in range(eDict["nconfigs"]):
                rStruct = randrange(1,enumController.nEnumStructs)
                print('Adding {} structure # {} to database'.format(eDict["lattice"],rStruct) )
                with open('structNums','a+') as f:
                    f.write(eDict["lattice"] + ' ' + str(rStruct) + '\n')
                    #print("Building VASP folder for {} structure #: {}".format(eDict["lattice"],rStruct))
                enumController.generatePOSCAR(rStruct)

                
                poscarpath = path.join(enumController.root,"poscar.{}.{}".format(eDict["lattice"],rStruct))
                thisCrystal = Crystal(poscarpath, systemSpecies = systemSpecies) #title = ' '.join([self.enumDicts[index]["lattice"]," str #: {}"]).format(rStruct)
                self.crystals.append(thisCrystal)
                delpath = path.join(enumController.root,"poscar.{}.{}".format(eDict["lattice"],rStruct))
                remove(delpath)

    # Sometimes an entire dataset is stored in one file.  I'd like to extract each crystal from the file to 
    # create a list of crystal objects
    def init_file(self,datafile,linesformat):
        from os import path
        handler = {'new_training.cfg':lambda file: self._init_mlpadd(file),'train.cfg': 'mlptrain','structures.in':'ce','dataReport_VASP.txt': lambda file: self._init_dataReport(file)}
        #selectedFile = path.join(self.root,'new_training.cfg')

        handler[path.split(datafile)[-1]](datafile)

    def _init_dataReport(self,datafile):
        with open(datafile,'r') as f:
            lines = f.readlines()

        del lines[:4]
        self.formationenergies = [ float(x.split()[-5]) for x in lines]
        self.concs = [ float(x.split()[-4]) for x in lines]
        
    def _init_mlpadd(self,datafile):
        from aBuild.database.crystal import Crystal
        with open(datafile,'r') as f:
            lines = f.readlines()

        self.crystals = []
        nCrystals = 0
        for index,line in enumerate(lines):
            print(index)
            if line == 'BEGIN_CFG\n':
                #                nCrystals += 1
                #if numOfStructs is not 'all' and (nCrystals < start or nCrystals > start + numOfStructs):
                #    continue
                nAtoms = int(lines[index+2].split()[0])
                structlines = lines[index:index + 18 + nAtoms]
                thisCrystal = Crystal(structlines,self.species,lFormat = 'mlpselect')
                print(thisCrystal.title, "HERE")
                import sys
                sys.exit()
                self.crystals.append(thisCrystal)



    def init_paths(self,paths,systemSpecies):
        from aBuild.database.crystal import Crystal
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.lammps import LAMMPS
        from os import path
        
        self.crystals = []
        for dirpath in paths:
            if self.calculator == 'VASP':
                calc = VASP(dirpath,systemSpecies)
                calc.read_results()

            #Added for LAMMPS compatibility
            if self.calculator == 'LAMMPS':
                calc = LAMMPS(dirpath,systemSpecies)
                calc.read_results()
            if calc.crystal.results is not None:
                self.crystals.append(calc.crystal)

            

    def starting_point(self,folderpath):
        from os import path
        from glob import glob
        from aBuild.utility import chdir

        with chdir(folderpath):
            dirs = glob('E.*')
        prevCalcs = [int(x.split('.')[1])  for x in dirs]
        prevCalcs.sort()
        if prevCalcs != []:
            return prevCalcs[-1] + 1
        else:
            return 1

        #    def write(self):

        
    def buildFolders(self,buildpath,calculator,runGetKpoints = True,foldername = 'E'):
        from os import path
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.lammps import LAMMPS
        from aBuild.calculators.espresso import ESPRESSO
        from aBuild.jobs import Job

        import os
        print("Building folders in {}".format(buildpath))
        if not path.isdir(buildpath):
                os.mkdir(buildpath)
                print('Made path:',buildpath)
        configIndex = startPoint = self.starting_point(buildpath)

        lookupCalc = {'vasp': lambda specs: VASP(specs,self.species),
                  'qe': lambda specs: ESPRESSO(specs,self.species),
                      'lammps': lambda specs: LAMMPS(specs,self.species)}

        lookupSpecs = {'vasp': lambda crystal: {"incar":calculator["vasp"]["incar"],"kpoints":calculator["vasp"]["kpoints"], 'potcar':calculator["vasp"]["potcars"],"crystal":crystal},
                  'qe': lambda crystal : {"crystal":crystal, "pseudopotentials":calculator["qe"]["pseudopotentials"]},
                      'lammps': lambda crystal: {"crystal":crystal, "potential":calculator["lammps"]["potential"]} }

        lookupBuild = {'vasp': lambda obj: obj.buildFolder(runGetKPoints = runGetKpoints),
                  'qe': lambda obj:obj.buildFolder(),
                      'lammps': lambda obj: obj.buildFolder()} 

        for crystal in self.crystals:
            # Initialize the calculation object
            thisCalc = lookupCalc[calculator["active"]](lookupSpecs[calculator["active"]](crystal))
                                  
            # Build the path
            runpath = path.join(buildpath,foldername + ".{}".format(configIndex) )
            if not path.isdir(runpath):
                os.mkdir(runpath)
            else:
                msg.fatal("I'm gonna write over top of a current directory. ({})  I think I'll stop instead.".format(runpath))

            # Change the directory and build the folder
            print("Building folder for structure: {}".format(crystal.title) )
            with chdir(runpath):
                lookupBuild[calculator["active"]](thisCalc)
            configIndex += 1


        # Build the submission script
        exdir = path.join(buildpath,'E.')
        mljob = Job(calculator["execution"],exdir,calculator["execution"]["exec_path"], arrayStart = startPoint,arrayEnd = configIndex - 1)
        with chdir(buildpath):
            print('Building job file')
            mljob.write_jobfile('jobscript_vasp.sh')






        
    


    def build_relax_select_input(self):
        from os import remove,path
        from aBuild.enumeration import Enumerate
        from aBuild.database.crystal import Crystal
        from aBuild.fitting.mtp import MTP
        from aBuild.utility import unpackProtos,getAllPerms
        from glob import glob
        fittingRoot = path.join(self.root,'fitting','mtp')
        
        for ilat  in range(self.nEnums):
            lat = self.enumDicts[ilat]["lattice"]
            enumLattice = Enumerate(self.enumDicts[ilat])

            if lat == 'protos':
                structures = getProtoPaths()
                for struct in structures:
                    scrambleOrder = getAllPerms(self.knary,justCyclic = 'uniqueUnaries' in struct)
                    for scramble in scrambleOrder:
                        thisCrystal = Crystal(struct,species = self.species)
                        thisCrystal.scrambleAtoms(scramble)
                        thisMTP = MTP(fittingRoot,dataSet = [thisCrystal],forRelax=True)
                        with open(path.join(fittingRoot,'to-relax.cfg'),'a+') as f:
                            f.writelines(thisMTP.lines)
                
            else:
                for struct in range(1,enumLattice.nEnumStructs+1):
                    enumLattice.generatePOSCAR(struct)
                    thisCrystal = Crystal.fromPOSCAR(enumLattice.root, self.species,
                                                     filename = "poscar.{}.{}".format(lat,struct),
                                                     title = ' '.join([lat," str #: {}"]).format(struct))
                    thisMTP = MTP(fittingRoot,dataSet = [thisCrystal],forRelax=True)
                    with open(path.join(fittingRoot,'to-relax.cfg'),'a+') as f:
                        f.writelines(thisMTP.lines)

                    delpath = path.join(enumLattice.root,"poscar.{}.{}".format(lat,struct))
                    remove(delpath)
                
        thisMTP.write_relaxin()



    def writeReport(self):
        import datetime
        nAtoms = len(self.crystals[0].species)
        with open('dataReport_' + self.calculator + '.txt', 'w') as f:
            f.write(self.calculator + ' REPORT\n')
            f.write(str(datetime.datetime.now()) + '\n')
            f.write("{:35s} {:14s}{:14s}{:12s}{:10s}{:9s}".format("Title"," T. Energy","F. Energy","Conc.",self.crystals[0].species[0] + "-atoms",self.crystals[0].species[1] + "-atoms\n"))
            f.write('------------------------------------------------------------------------------------------------------------------\n')
            for crystal in self.crystals:
                f.write(crystal.reportline)

    def generateConvexHullPlot(self):
        from scipy.spatial import ConvexHull
        from numpy import array
        from matplotlib import pyplot
        #with open('dataReport_VASP.txt','r') as f:
        #    lines = f.readlines()

        #del lines[0:4]
        #data = [[float(x.split()[-4]),float(x.split()[-5] )] for x in lines]
        #data = [[i.results["fEnth"],i.concentrations[0]] for x in self.crystals]
        data = array([[self.concs[i],self.formationenergies[i]] for i in range(len(self.formationenergies))])
        print(data,'data')
        hull = ConvexHull(data)
        pyplot.plot(self.concs,self.formationenergies,'r+')
        plotConcs = []
        plotEnergies = []
#        pyplot.plot(data[hull.vertices,0], data[hull.vertices,1],'b-',lw = 2)

        for vertex in hull.vertices:
            if self.formationenergies[vertex] <= 0:
                plotConcs.append(data[vertex,0])
                plotEnergies.append(data[vertex,1])
        pyplot.plot(plotConcs,plotEnergies,'k-')
        pyplot.savefig('chull.png')
        









                
