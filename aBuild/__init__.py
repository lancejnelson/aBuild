
from aBuild import config  # Import environment variables that may be needed.
from aBuild import msg # messaging module
from aBuild.utility import chdir, _get_reporoot
import sys

config = sys.modules["config"]  

class Controller(object):
    """ Args:
        config (str): name of the YML file (without the .yml) that
          specifies all information for constructing the set of databases.
        tmpdir (str): path to a temporary directory to use for the
          database. This is for unit testing purposes.
    """
    def __init__(self,inputFile):
        from aBuild.io import read
        from os import path,getcwd

        self.root = getcwd()
        self.inputFile = path.expanduser(path.abspath(inputFile))  # Get full path to the input file.

        # Get the folder and file to the input file
        if path.isabs(inputFile):
            root, inputFile = path.split(inputFile)
        else:
            root,inputFile = path.dirname(self.inputFile),inputFile

        # Read the input file
        self.specs = read(root,inputFile)

        if "directory" not in self.specs["calculator"]["vasp"]["potcars"] or self.specs["calculator"]["vasp"]["potcars"]["directory"] is None:
            print("You did not provide a directory for the POTCARS. Using the environment variable that I found: {}".format(config.POTCAR_DIR))
            self.specs["calculator"]["potcars"]["directory"] = config.POTCAR_DIR

        self.root = path.expanduser(self.specs["root"])  # Set the working directory

        if getcwd().replace("zhome","fslhome") != self.root:
            infoMsg = 'You have specified a working directory  ({}) that is different from your current working directory ({}) .'.format(self.root,getcwd())
            msg.info(infoMsg)

        self.title = self.specs['title']
        self.species = sorted(self.specs["species"],reverse = True)
        self.fpExecution = self.specs.get("execution", {})  # Second argument is the default in case the key doesn't exist.
        self.calculator = self.specs.get("calculator", {}) # Second argument is the default in case the key doesn't exist.
        #        self.nEnumStructuresToSelect = self.specs["trainingset"]["nconfigs"]
        self.knary = len(self.species)
        self.fitting = self.specs.get("fitting",{})


    def setDataSet(self,dataSetType):
        self.dataset = dataSetType
        

    @property
    def nconfigs(self):
        return self.specs[self.dataset]["nconfigs"]
    
    @property
    def enumLattices(self):
        if "lattice" in self.specs[self.dataset]:
            return self.specs[self.dataset]["lattice"]
        else:
            ermsg = 'You have not specified a lattice type'
            msg.fatal(ermsg)  # Print the error message
            

    @property
    def nEnums(self):
        return len(self.enumLattices)

    @property
    def enumConcs(self):
        if "concs" in self.specs[self.dataset]:
            return self.specs[self.dataset]["concs"]
        else:
            return [None for k in range(self.nEnums)]


    @property
    def atomicBasis(self):
        if "basis" in self.specs[self.dataset]:
            return self.specs[self.dataset]["basis"]
        else:
            # Default to [0,0,0] when the enum object gets called up
            return [None for k in range(self.nEnums)]


    @property
    def concRestrictions(self):
        if "concs" in self.specs[self.dataset]:
            return self.specs[self.dataset]["concs"]
        else:
            return [None for k in range(self.nEnums)]


    @property
    def enumSizes(self):
        if "sizes" in self.specs[self.dataset]:
            return self.specs[self.dataset]["sizes"]
        else:
            ermsg = 'You have not specified the enumeration sizes in your input file'
            msg.err(ermsg)  # Print the error message
            import sys
            sys.exit()

    
    @property
    def siteRestrictions(self):
        if "siteRestrictions" in self.specs[self.dataset]:
            return self.specs[self.dataset]["siteRestrictions"]
        else:
            return [None for k in range(self.nEnums)]

    @property
    def enumDicts(self):
        edicts = []
        for i in range(self.nEnums):
            edict = {}
            edict["lattice"] = self.enumLattices[i]
            edict["basis"] = self.atomicBasis[i]
            edict["knary"] = self.knary
            edict["nconfigs"] = self.nconfigs[i]
            edict["sizes"] = self.enumSizes[i]
            edict["site_res"] = self.siteRestrictions[i]
            edict["concs"] = self.concRestrictions[i]
            edict["root"] = self.root
            edict["eps"] = 1e-3
            edicts.append(edict)
        
        return edicts
    
    def enumerate(self,dataset):
        from aBuild.enumeration import Enumerate

        self.dataset = dataset
        for index in range(self.nEnums):
            enumController = Enumerate(self.enumDicts[index])
            enumController.buildInputFile(False)
            enumController.enumerate(False)


    # BUild VASP folders so I can generate training data
    def setup_training_set(self,runGetKpoints = True):
        from os import path
        
        from aBuild.database.dataset import dataset

        self.dataset = "trainingset"
        trainingRoot = path.join(self.root,'training_set')
        trainingSet = dataset(self.enumDicts,self.species)
        trainingSet.buildFolders(trainingRoot,self.calculator,runGetKpoints = runGetKpoints)
        

    # Pull data from VASP calculations and get ready to fit.
    def setup_training_input(self):
        from os import path
        from glob import glob
        from aBuild.calculators.vasp import VASP
        #from aBuild.database.crystal import CrystalsList
        from aBuild.fitting.mtp import MTP
        from aBuild.jobs import Job

        from aBuild.database.dataset import dataset


        # 1. Pull all VASP data and compile into file train.cfg
        # 2. Copy/Build a blank potential file called pot.mtp
        # 3. Build a submission script

        trainingRoot = path.join(self.root, 'training_set')
        with chdir(trainingRoot):
            enumdirs = glob("E.*")
            activedirs = glob("A.*")

        dirs = [path.join(trainingRoot,x) for x in enumdirs + activedirs]
        #        dirs = enumdirs + activedirs

        print('Building dataset')
        trainingSet = dataset(dirs,self.species,calculator='VASP')  
        
        fittingRoot = path.join(self.root,'fitting','mtp')
        thisMTP = MTP(fittingRoot,settings = self.fitting)
        thisMTP.write_blank_pot(self.knary)

        with open(path.join(fittingRoot,'train.cfg'),'a+') as f:
            for crystal in trainingSet.crystals:
                if crystal.results["warning"]:
                    print("Not adding crystal: {} because the energy doesn't look right. Check it out.".format(calc.crystal.title) )
                else:
                    f.writelines('\n'.join(crystal.lines('mtptrain')))


        thisMTP.train(self.fitting["execution"])
#        mlpCommand = 'mlp train pot.mtp train.cfg\n'
#        mljob = Job(self.fitting["execution"],path.join(self.root,"fitting","mtp"),mlpCommand)
#        with chdir(path.join(self.root,"fitting/mtp")):
#            print('Building job file')
#            mljob.write_jobfile()


    # Build files needed to run the relaxation. Lots of structures. Need to write as crystals
    # are generated to save memory
    def setup_relax_input(self,freshStart = False):
        from os import remove,path
        from aBuild.fitting.mtp import MTP
        from glob import glob

        self.dataset = "gss"


        # 1.  If to-relax.cfg does not exit, build it.  This
        #     takes a while because we are putting a bunch of
        #     structures in there.
        # 2.  If to-relax.cfg does exist (iteration > 1), then
        #     copy unrelaxed.cfg from previous iteration to to-relax.cfg
        #     Note that there may be many unrelaxed.cfg_# from the parallel run
        # 3.  Copy Trained.mtp to pot.mtp in preparation for the relaxation
        # 4.  Build a relax.ini from template if it's not already there.
        # 5.  Run calc-grade
        # 6.  Build a submission script.
        fittingRoot = path.join(self.root,'fitting','mtp')
        filePath = path.join(fittingRoot,'to-relax.cfg')
        #unrelaxedfilePath = path.join(fittingRoot,'unrelaxed.cfg')

        torelax = path.isfile(filePath)
        #unrelaxed = path.isfile(unrelaxedfilePath)

        # 1.
        if not torelax:  #This must be the first iteration
            self.build_ToRelax()

        # 2.
        else: # What iteration is it?
            if freshStart:
                remove(filePath)
                self.build_ToRelax()
         #   elif unrelaxed:  #Iteration must be > 1
         ##       unrelaxedFiles = glob(path.join(fittingRoot,"unrelaxed.cfg_*"))
          ##      cat(unrelaxedFiles,path.join(fittingRoot,"to-relax.cfg"))
           #     remove(unrelaxedFiles)
            else: #"to-relax.cfg is present and I don't wanta  fresh start.  Don't do anything."
                msg.info("to-relax.cfg found and you told me you don't want to start fresh.  Proceeding with the to-relax.cfg that's there.")
        #elif path.isfile(filePath):  # Iteration # must be > 1
        #        rename(path.join(self.root,'fitting/mtp')
        

        # 3.
        try:
            print('renaming Trained.mtp to pot.mtp')
            rename(path.join(fittingRoot,'Trained.mtp'),path.join(fittingRoot,'pot.mtp'))
        except:
            if path.isfile(path.join(fittingRoot,'pot.mtp')):
                msg.info('It looks like you have already copied Trained.mtp -> pot.mtp previously')
            else:
                msg.fatal("Can't find a pot.mtp or a Trained.mtp.  Problems")
                          
        fittingRoot = path.join(self.root,'fitting','mtp')
        thisMTP = MTP(fittingRoot,settings = self.fitting)
        # 4.
        thisMTP.write_relaxin()
        # 5.
        thisMTP.calc_grade()
        # 6.
        thisMTP.relax(self.fitting["execution"])



        #        trainingSet = dataset(enumspecs = self.enumDicts,root = self.root,species = self.species,calculator = self.calculator)
        #trainingSet.build_relax_select_input()

    def setup_select_input(self):
        from os import path,remove
        from shutil import copy
        from glob import glob
        from aBuild.utility import cat
        from aBuild.fitting.mtp import MTP
        # 1. Concatenate all of the candidate.cfg_# into one candidate.cfg
        #    and remove all of the other ones.
        # 2. Concatenate all of the selection.log_* files into one file.
        # 3. Concatenate all of the relaxed.cfg_# into one file.  This file should
        #    get bigger and bigger with each iteration (Hopefully), but we don't ever
        #    want to delete one of these files.
        # 4. Build a submission script.

        fittingRoot = path.join(self.root,'fitting','mtp')

        with open(path.join(fittingRoot,'iteration'),'r') as f:
            lines = f.readlines()
        iteration = int(lines[0].split()[0])

        with open(path.join(fittingRoot,'iteration'),'w') as f:
            f.write("{:d}\n".format(iteration + 1))
        
        # 1.
        candFiles = glob(path.join(fittingRoot,"candidate.cfg_*"))
        cat(candFiles,path.join(fittingRoot,"candidate.temp"))
        for file in candFiles:
            remove(file)
        copy( path.join(fittingRoot,'candidate.temp'), path.join(fittingRoot,"candidate_iteration_" + str(iteration) + ".cfg") )
        copy( path.join(fittingRoot,'candidate.temp'), path.join(fittingRoot,"candidate.cfg" ) )
        remove(path.join(fittingRoot,'candidate.temp'))
        
        # 2.
        logFiles = glob(path.join(fittingRoot,"selection.log_*"))
        cat(logFiles,path.join(fittingRoot,"selection.temp"))
        for file in logFiles:
            remove(file)
        copy( path.join(fittingRoot,'selection.temp'), path.join(fittingRoot,"selection_iteration_" + str(iteration) + ".log") )
#        copy( path.join(fittingRoot,'selection.temp'), path.join(fittingRoot,"selection.log" ) )
        remove(path.join(fittingRoot,'selection.temp'))

        
        
        # 3.
        relaxedFiles = glob(path.join(fittingRoot,"relaxed.cfg_*"))
        cat(relaxedFiles,path.join(fittingRoot,"relaxed.temp"))
        for file in relaxedFiles:
            remove(file)
        copy( path.join(fittingRoot,'relaxed.temp'), path.join(fittingRoot,"relaxed_iteration_" + str(iteration) + ".cfg") )
#        copy( path.join(fittingRoot,'relaxed.temp'), path.join(fittingRoot,"relaxed.cfg" ) )
        remove(path.join(fittingRoot,'relaxed.temp'))

            # 3.
        unrelaxedFiles = glob(path.join(fittingRoot,"unrelaxed.cfg_*"))
        cat(unrelaxedFiles,path.join(fittingRoot,"unrelaxed.temp"))
        for file in unrelaxedFiles:
            remove(file)
        copy( path.join(fittingRoot,'unrelaxed.temp'), path.join(fittingRoot,"unrelaxed_iteration_" + str(iteration) + ".cfg") )
#        copy( path.join(fittingRoot,'relaxed.temp'), path.join(fittingRoot,"relaxed.cfg" ) )
        remove(path.join(fittingRoot,'unrelaxed.temp'))

        relaxLogFiles = glob(path.join(fittingRoot,"relax_log.txt_*"))
        cat(relaxLogFiles,path.join(fittingRoot,"relax_log.temp"))
        for file in relaxLogFiles:
            remove(file)
        copy( path.join(fittingRoot,'relax_log.temp'), path.join(fittingRoot,"relax_log_iteration_" + str(iteration) + ".txt") )
#        copy( path.join(fittingRoot,'relaxed.temp'), path.join(fittingRoot,"relaxed.cfg" ) )
        remove(path.join(fittingRoot,'relax_log.temp'))

        # 4.
        thisMTP = MTP(fittingRoot)
        thisMTP.select_add(self.fitting["execution"])
        
    def build_ToRelax(self):
        from aBuild.enumeration import Enumerate
        from aBuild.utility import unpackProtos,getAllPerms,getProtoPaths
        from aBuild.database.crystal import Crystal
        from os import remove,path
        print('Building to-relax.cfg')
        fittingRoot = path.join(self.root,'fitting','mtp')
        for ilat  in range(self.nEnums):
            lat = self.enumDicts[ilat]["lattice"]
            
            if lat == 'protos':
                structures = getProtoPaths(self.knary)
                print(structures)
                #                subdivide = [structures[x:x+100] for x in range() ]
                for struct in structures:
                    print("Proto structure:", struct)
                    scrambleOrder = getAllPerms(self.knary,justCyclic = 'uniqueUnaries' in struct)
                    for scramble in scrambleOrder:
                        thisCrystal = Crystal(struct,self.species)
                        #print("Atom counts before scramble {}".format(thisCrystal.atom_counts))
                        thisCrystal.scrambleAtoms(scramble)
                        #print("Atom counts after scramble {}".format(thisCrystal.atom_counts))
                        with open(path.join(fittingRoot,'to-relax.cfg'),'a+') as f:
                            f.writelines('\n'.join(thisCrystal.lines('mtprelax') ) )
                
            else:
                enumLattice = Enumerate(self.enumDicts[ilat])
                for struct in range(1,enumLattice.nConfigs+1):
                    print("Lattice",lat, "structure:",struct)
                    enumLattice.generatePOSCAR(struct) 
                    thisCrystal = Crystal(path.join(enumLattice.root,"poscar.{}.{}".format(lat,struct)),self.species)
                    with open(path.join(fittingRoot,'to-relax.cfg'),'a+') as f:
                        f.writelines('\n'.join(thisCrystal.lines('mtprelax') ))

                    delpath = path.join(enumLattice.root,"poscar.{}.{}".format(lat,struct))
                    remove(delpath)

    def augmentTraining(self):
        from os import path
        from aBuild.database.dataset import dataset

        newTraining = path.join(self.root,'fitting','mtp','new_training.cfg')
        trainingRoot = path.join(self.root,'training_set')
        dSet = dataset(newTraining,self.species,lFormat = 'mlpselect')
        dSet.buildFolders(trainingRoot,self.calculator,foldername = 'A')
    def statusReport(self):
        from os import path
        from glob import glob
        from aBuild.calculators.vasp import VASP

        trainingRoot = path.join(self.root, 'training_set')
        with chdir(trainingRoot):
            enumdirs = glob("E.*")
            activedirs = glob("A.*")

        dirs = [path.join(trainingRoot,x) for x in enumdirs + activedirs]

        stat = {'done':[],'running':[], 'not started': [], 'error':[], 'not setup':[],'warning':[],'idk':[]}
        for dir in dirs:
            thisVASP = VASP(dir,self.species)
            stat[thisVASP.status()].append(dir.split('/')[-1])
            #            msg.info("Status of directory {} is {} ".format(dir,thisVASP.status()))
        msg.info('Done')
        msg.info(' '.join(stat['done']))
        msg.info('Running')
        msg.info(' '.join(stat['running']))
        msg.info('Not Started')
        msg.info(' '.join(stat['not started']))
        msg.info('Not Setup')
        msg.info(' '.join(stat['not setup']))
        msg.info('Errors')
        msg.info(' '.join(stat['error']))
        msg.info('Warnings')
        msg.info(' '.join(stat['warning']))
        msg.info('Not sure')
        msg.info(' '.join(stat['idk']))

    def gatherResults(self):
        from os import path
        from glob import glob
        from aBuild.database.dataset import dataset

        

        trainingRoot = path.join(self.root, 'training_set')
        with chdir(trainingRoot):
            enumdirs = glob("E.*")
            activedirs = glob("A.*")
            pures = glob("pure*")

        dirs = [path.join(trainingRoot,x) for x in enumdirs + activedirs + pures]
        #        dirs = enumdirs + activedirs

        print('Building dataset')
        if self.calculator["active"].lower() == "lammps":
            trainingSet = dataset(dirs,self.species,calculator = 'LAMMPS')
            print('here')
            trainingSet.writeReport()
        if self.calculator["active"].lower() == 'vasp': 
            trainingSet = dataset(dirs,self.species,calculator = 'VASP')
            print('here')
            trainingSet.writeReport()

    def generateConvexHull(self):
        from os import path
        from aBuild.database.dataset import dataset
        dataFile = path.join(self.root,'dataReport_VASP.txt')
        if not path.isfile(dataFile):
            self.gatherResults()

        data = dataset(dataFile,self.species)
        data.generateConvexHullPlot()
            
