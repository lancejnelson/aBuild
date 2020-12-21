
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

        if "srcdirectory" not in self.specs["calculator"]["vasp"]["potcar"] or self.specs["calculator"]["vasp"]["potcar"]["srcdirectory"] is None:
            print("You did not provide a directory for the POTCARS. Using the environment variable that I found: {}".format(config.POTCAR_DIR))
            self.specs["calculator"]["potcar"]["srcdirectory"] = config.POTCAR_DIR

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
    def enums(self):
        from aBuild.enumeration import Enumerate

        nEnums = len(self.specs[self.dataset]["lattice"])

        enums = []
        for index in range(nEnums):
            enums.append(Enumerate.fromYAML(self.specs[self.dataset],index,self.knary,self.root))
        return enums

    def enumerate(self,dataset):
        from aBuild.enumeration import Enumerate

        nEnums = len(self.specs[dataset]["lattice"])

        for index in range(nEnums):
            enumController = Enumerate.fromYAML(self.specs[dataset],index,self.knary,self.root)
            enumController.buildInputFile(False)
            enumController.enumerate(False)


    # BUild VASP folders so I can generate training data
    # This mode is no longer used because with MTP we can
    # start with an empty training set and have MPT pick the
    # first set of training structures
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
                   'setup_relax':lambda:thisMTP.setup_relax(self.enums,self.species,AFM = 'AFM' in self.calculator[self.calculator["active"]].keys(), start = start,end = end),
                    'setup_select_add':lambda :thisMTP.setup_select()}
        handler[tag]()


    def augmentTraining(self):
        from os import path
        from aBuild.database.dataset import dataset

        newTraining = path.join(self.root,'fitting','mtp','new_training.cfg')
        trainingRoot = path.join(self.root,'training_set')
        dSet = dataset.from_file(newTraining,self.species,'mlp')
        dSet.buildFolders(trainingRoot,self.calculator,foldername = 'A')

    def statusReport(self):

        from os import path,listdir,remove
        from glob import glob
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.aflow import AFLOW
        from numpy import argsort
        trainingRoot = path.join(self.root, 'validation_set')
        with chdir(trainingRoot):
            enumdirs = glob("E.*")
            activedirs = glob("A.*")

        delFiles = glob("status.*")
        for i in delFiles:
            remove(i)
        dirs = [path.join(trainingRoot,x) for x in enumdirs + activedirs]
        nDirs = len(dirs)
        numdirs = [int(y.split('.')[-1]) for y in dirs]
        count = 0
        for idx,directory in enumerate([dirs[x] for x in argsort(numdirs)]):
            if count % 100 == 0:
                print("checking directories {} - {}".format(dirs[argsort(numdirs)[count+1]],dirs[argsort(numdirs)[count + 100 if count + 100 < nDirs else nDirs - 1]]))
            count += 1
            if 'aflow.in' in listdir(directory):
                thisAFLOW = AFLOW.from_path(directory,self.species)
                thisStat = thisAFLOW.status()

            else:
                thisVASP = VASP.from_path(directory,self.species)
                thisStat = thisVASP.status()
            if thisStat == 'errors' or thisStat == 'running':
                msg.info("Directory {} is {} ".format(directory, thisStat))

            with open('status.' + thisStat,'a') as f:
                f.write(directory.split('/')[-1] + '  ')

    def gatherResults(self,datafile=None, folder = None):
        from os import path
        from glob import glob
        from aBuild.database.dataset import dataset
        from numpy import argsort

        datafile = path.join(self.root,'fitting','mtp','relaxed_iteration_6.cfg')
        if datafile is not None:
            dataSet = dataset.from_file(datafile,self.species)
            data.generateConvexHull(plot2D= False,plot3D=False,fileOutput = True)  # Find the convex hull
            dataSet.findAllHullDistances()
            dataSet.writeReport('mlptest')
        else:
            trainingRoot = path.join(self.root, 'validation_set')
            with chdir(trainingRoot):
                enumdirs = glob("E.*")
                activedirs = glob("A.*")
                pures = glob("pure*")

            dirs = [path.join(trainingRoot,x) for x in enumdirs + activedirs ]
            numdirs = [int(y.split('.')[-1]) for y in dirs]
            sorteddirs = [dirs[x] for x in argsort(numdirs)]
            trainingSet = dataset.from_paths(sorteddirs,self.species,calculator = self.calculator["active"].upper())
            trainingSet.writeReport('vasp')


    def setupHullPredictions(self):
        from os import path
        from aBuild.database.dataset import dataset
        self.dataset = 'gss'

        dataFile = 'dataReport_mlp.txt'

        data = dataset.from_file(dataFile,self.species,getCrystallographicInfo = True,eDicts = self.enumDicts,onlyCloseToHull = True)
        validationRoot = path.join(self.root,'validation_set')
        data.buildFolders(validationRoot,self.calculator,foldername='A')  # Submit VASP folders for close ones


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



#    def validateGroundStates(self):
