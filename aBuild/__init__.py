
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
    def name(self):
        if "name" in self.specs[self.dataset]:
            savedNames = self.specs[self.dataset]["name"]
            for idx,name in enumerate(savedNames):
                if name is None:
                    savedNames[idx] = self.specs[self.database]["lattice"][idx]
            return savedNames#self.specs[self.dataset]["name"]
        else:
            return [None for k in range(self.nEnums)]
        
    @property
    def coordsys(self):
        if "coordys" in self.specs[self.dataset]:
            return self.specs[self.dataset]["coordys"]
        else:
            return [None for k in range(self.nEnums)]

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
            edict["coordsys"] = self.coordsys[i]
            edict["name"] = self.name[i]
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
        trainingSet = dataset(self.enumDicts,self.species,restrictions = 'AFM' in self.calculator[self.calculator["active"]])
        trainingSet.buildFolders(trainingRoot,self.calculator,runGetKpoints = runGetKpoints)
        
        

    def setupHandler(self,model,tag,start = 1,end = None):
        from aBuild.fitting.mtp import MTP
        from os import path
        
        fittingRoot = path.join(self.root, 'fitting/mtp')
        trainingRoot = path.join(self.root, 'training_set')

        self.dataset = 'gss'
        thisMTP = MTP(fittingRoot,settings=self.fitting)
        handler = {'setup_train':lambda: thisMTP.setup_train(trainingRoot,self.species),
                   'setup_relax':lambda:thisMTP.setup_relax(self.enumDicts,self.species,AFM = 'AFM' in self.calculator[self.calculator["active"]].keys(), start = start,end = end),
                    'setup_select_add':lambda :thisMTP.setup_select()}
        handler[tag]()


    def augmentTraining(self):
        from os import path
        from aBuild.database.dataset import dataset

        newTraining = path.join(self.root,'fitting','mtp','new_training.cfg')
        trainingRoot = path.join(self.root,'training_set')
        dSet = dataset(newTraining,self.species,lFormat = 'mlp')
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
        stat = {'done':[],'running':[], 'not started': [], 'too long':[], 'not setup':[],'warning':[],'idk':[],'unconverged':[],'sgrcon':[],'error':[]}
        for dir in dirs:
            thisVASP = VASP(dir,systemSpecies = self.species)
            stat[thisVASP.status()].append(dir.split('/')[-1])
            #            msg.info("Status of directory {} is {} ".format(dir,thisVASP.status()))
        msg.info('Done (' + str(len(stat['done'])) + ')')
        msg.info(' '.join(stat['done']))
        msg.info('Running (' + str(len(stat['running'])) + ')')
        msg.info(' '.join(stat['running']))
        msg.info('Not Started (' + str(len(stat['not started'])) + ')')
        msg.info(' '.join(stat['not started']))
        msg.info('Not Setup (' + str(len(stat['not setup'])) + ')')
        msg.info(' '.join(stat['not setup']))
        msg.info('Too Long (' + str(len(stat['too long'])) + ')    (last write was > 1 hr ago) ')
        msg.info(' '.join(stat['too long']))
        msg.info('Warnings (' + str(len(stat['warning'])) + ')')
        msg.info(' '.join(stat['warning']))
        msg.info('Not sure (' + str(len(stat['idk'])) + ')  (Found finish tags but couldn''t find final energy.  It''s probably in the process of finishing.)')
        msg.info(' '.join(stat['idk']))
        msg.info('Unconverged (' + str(len(stat['unconverged'])) + ')')
        msg.info(' '.join(stat['unconverged']))
        msg.info('SGRCON error (' + str(len(stat['sgrcon'])) + ')')
        msg.info(' '.join(stat['sgrcon']))
        msg.info('Unknown error (' + str(len(stat['error'])) + ')')
        msg.info(' '.join(stat['error']))

    def gatherResults(self,file=None, folder = None):
        from os import path
        from glob import glob
        from aBuild.database.dataset import dataset
        if file is not None:
            dataSet = dataset(file,self.species,lFormat = 'mlp',)
            dataSet.writeReport('mlp')
        else:
            trainingRoot = path.join(self.root, 'training_set')
            with chdir(trainingRoot):
                enumdirs = glob("E.*")
                activedirs = glob("A.*")
                pures = glob("pure*")

            dirs = [path.join(trainingRoot,x) for x in enumdirs + activedirs + pures]
            trainingSet = dataset(dirs,self.species,calculator = self.calculator["active"].upper())
            trainingSet.writeReport('vasp')


    def generateConvexHull(self,file='dataReport_VASP.txt',plotAll = True):
        from os import path
        from aBuild.database.dataset import dataset
        dataFile = path.join(self.root,file)
        if not path.isfile(dataFile):
            msg.fatal('data file does not exist')

        data = dataset(dataFile,self.species)
        data.generateConvexHullPlot(plotAll = plotAll)
            
    def errorsReport(self,datafile = None, predictFile = None):
        from aBuild.database.dataset import dataset
        from numpy.linalg import norm
        from numpy import average,array

        dataSet = dataset(datafile,self.species,lFormat = 'mlp')
        predictSet = dataset(predictFile,self.species,lFormat = 'mlp')

        diffsEnergy = []
        diffsForces = []
        for i in dataSet.crystals:
            for j in predictSet.crystals:
                if j.title.strip() == i.title.strip():
                    print("Found match")
                    diffsEnergy.append(i.results["energyF"] - j.results["energyF"])
                    diffsForces.append(average(norm(array(i.results["forces"]) - array(j.results["forces"]), axis = 1)))
        from matplotlib import pyplot
        pyplot.subplot(121)
        pyplot.hist(diffsEnergy,bins = 30)
        pyplot.subplot(122)
        pyplot.hist(diffsForces,bins = 30)
        pyplot.savefig('errorPlot.png')
        #pyplot.savefig('errorForces.png')
        
        

    def randomDisplacements(self,POSCAR):
        from os import path
        
        from aBuild.database.crystal import Crystal

        toRelax = path.join(self.root,'fitting','mtp','to_relax.cfg')
        for i in range(1000):
            thisCrystal = Crystal(POSCAR,systemSpecies = self.species)
            thisCrystal.randomDisplace()
            print(thisCrystal.lines('mtprelax'))
            with open(toRelax, 'a+') as f:
                f.writelines('\n'.join(thisCrystal.lines('mtprelax')))
                
        
