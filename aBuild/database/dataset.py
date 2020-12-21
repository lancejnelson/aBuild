
from aBuild import msg
from aBuild.utility import chdir, _get_reporoot
from aBuild.database.crystal import Crystal



def pointInPolygon(point,polygon):
    # This is the winding number algorithm.
    from numpy.linalg import norm
    from numpy import arccos,array,dot,pi
    refVec = array([1,0])
    angleSum = 0
    for r in range(3):
        vecOne = polygon[r] - point
        vecTwo = polygon[(r+1)%3] - point
        angleOne = arccos(dot(vecOne,refVec)/norm(vecOne)) if vecOne[1] > 0 else 2 * pi - arccos(dot(vecOne,refVec)/norm(vecOne))
        angleTwo = arccos(dot(vecTwo,refVec)/norm(vecTwo)) if vecTwo[1] > 0 else 2 * pi - arccos(dot(vecTwo,refVec)/norm(vecTwo))

        if abs(abs(angleTwo - angleOne) - pi) < 1e-4:
            # If the angle between two adjacent vectors is pi, then we found the plane we need.
            return True
#        theta = arccos(dot(vecOne,vecTwo)/(norm(vecOne) * norm(vecTwo) ))
        if abs(angleTwo - angleOne) < pi:
            if angleTwo - angleOne < 0:
                angleSum += 1
            else:
                angleSum += -1
        else:
            if angleTwo - angleOne < 0:
                angleSum += -1
            else:
                angleSum += 1
    if abs(angleSum) != 3:
        return False
    else:
        return True

class dataset:


    def __init__(self,crystals):
        from os import path,makedirs
        from aBuild.calculators.vasp import VASP

        self.crystals = crystals
        self.species = crystals[0].crystalSpecies


    """
    Method for initializing a dataset object from an enumeration dictionary.
    This method will:
        1- Enumerate derivative crystal structures (if needed)
        2- Generate crystal structures from the enumeration output
        3- Build a list of Crystal objects to initialize with.

    Arg:
        enumDict(dictionary): dictionary containing enumeration information.
        systemSpecies(list of strs): Needed to initialize the Crystal objects
    """
    @staticmethod
    def init_enum(enumdicts,systemSpecies):
        from aBuild.enumeration import Enumerate
        from aBuild.calculators.vasp import VASP
        from aBuild.database.crystal import Crystal
        from aBuild.jobs import Job
        from random import randrange
        from aBuild.utility import chdir
        from numpy import array
        from os import remove, path

        #        from crystal import Crystal
        from os import path
        import os

        print("Building database from enumerations")
        crystals = []
        #        configIndex = startPoint = self._starting_point
        for eDict in enumdicts:
            enumController = Enumerate(eDict)
            if enumController.nEnumStructs == 0:
                msg.warn('There are no enumerated structures for lattice type {}.  Not building any VASP folders for them.'.format(eDict["lattice"]))
                enumController.buildInputFile()

                enumController.enumerate()

            # Loop to generate random structures for a given lattice type
            for i in range(eDict["nconfigs"]):
                print('Adding {} structure # {} to database'.format(eDict["lattice"],rStruct) )
                with open('structNums','a+') as f:
                    f.write(eDict["name"] + ' ' + str(i) + '\n')
                enumController.generatePOSCAR(i)

                poscarpath = path.join(enumController.root,"poscar.{}.{}".format(eDict["name"],rStruct))
                thisCrystal = Crystal.from_poscar(poscarpath, systemSpecies) #title = ' '.join([self.enumDicts[index]["lattice"]," str #: {}"]).format(rStruct)
                crystals.append(thisCrystal)
                delpath = path.join(enumController.root,"poscar.{}.{}".format(eDict["name"],rStruct))
                remove(delpath)
        return dataset(crystals)

    """
    Method for initializing a dataset object from a single data file.  For example, train.cfg
    and new_training.cfg both contain hundreds of crystal structures.
    This method will:
        1- Recognize the file type
        2- Read the file and parse each crystal out.
        3- Construct a list of Crystal objects for initializing the dataset object.

    Arg:
        datafile(str): path to datafile
        species(list of strs): Needed to initialize the Crystal objects
        linesformat(str):  Specifies what type of file it is.
    """
    @staticmethod
    def from_file(datafile,species,eDicts = None,onlyCloseToHull = False, getCrystallographicInfo = False):
        from os import path
        handler = {'new_training.cfg':lambda file: dataset._init_mlp(file),'train.cfg': 'mlptrain','structures.in':'ce',}
        if 'relaxed' in datafile:
            handler[path.split(datafile)[-1]] = lambda file: dataset._init_mlp(file,species)
        if 'dataReport' in datafile:
            handler[path.split(datafile)[-1]] = lambda file: dataset._init_dataReport(file,species,onlyCloseToHull = onlyCloseToHull, getCrystallographicInfo = getCrystallographicInfo, enumDicts = eDicts)
        if 'train' in datafile:
            handler[path.split(datafile)[-1]] = lambda file: dataset._init_mlp(file,species)
        return handler[path.split(datafile)[-1]](datafile)

    @staticmethod
    def _init_dataReport(datafile,species,onlyCloseToHull = False, cutoff = 5e-3,getCrystallographicInfo=False,enums = None):
        if getCrystallographicInfo and enumDicts == None:
            msg.fatal("Can't extract crystallographic information without the enumeration dictionaries")

        with open(datafile,'r') as f:
            lines = f.readlines()

        del lines[:4]
        required = ["lattice","basis", "atom_types","crystalSpecies","latpar","coordsys","systemSpecies"]

        nAtoms = int((len(lines[0].split()) - 9)/2)

        lookup = {}
        for enum in enums:
            lookup[enum.lattice.lattice_name] = enum

        # In case we don't want the full crystallographic information,
        # let's build an empty crystal object that we can fill up with
        # the information that we want.
        templateDict = {}
        for key in required:
            templateDict[key] = None
        templateDict["results"] = {}
#        crystal = Crystal(templateDict)

        count = 0
        crystals = []
        for line in lines:
            formationEnthalpy = float(line.split()[8])
            energyF = float(line.split()[6])
            title = ' '.join(line.split()[:6])
            atomCounts = list(map(int,line.split()[-nAtoms:]))
            lattice = str(line.split()[1].split('_')[1])
            structNum = int(line.split()[5])

            if len(line.split()) == 16:
                distanceToHull = float(line.split()[9])
            else:
                distanceToHull = None
            if onlyCloseToHull:
                if distanceToHull > cutoff:
                    continue

            if getCrystallographicInfo:
                crystal = Crystal.fromEnum(lookup[lattice],structNum,species)
                crystal.results["fEnth"] = formationEnthalpy
                crystal.results["energyF"] = energyF
                crystal.results["distToHull"] = distanceToHull
            else:
                templateDict["fEnth"] = formationEnergy
                templateDict["energyF"] = energyF
                templateDict["distanceToHull"] = distanceToHull
                templateDict["title"] = title
                templateDict["atom_counts"] = atomCounts
                templateDict["crystalSpecies"] = species
                crystal = Crystal(templateDict)


            count += 1
            print("Read in crystal {}".format(count))
            crystals.append(crystal)
        return dataset(crystals)


    @staticmethod
    def _init_mlp(datafile,species):
        from aBuild.database.crystal import Crystal
        import os
        from os import path
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.aflow import AFLOW
        with open(datafile,'r') as f:
            lines = f.readlines()

        crystals = []
        nCrystals = 0
        # Get information for pures so I can calculate formation energies
        root = os.getcwd()
        trainingRoot = path.join(root,'training_set')
        puredirs = [path.join(trainingRoot,'pure' + x) for x in species]
        pures = [AFLOW.from_path(x,species,filesuffix = '.relax2.xz') for x in puredirs]
        for pure in pures:
            pure.read_results()
        # End reading pures

        count = 0
        for index,line in enumerate(lines):
            if 'BEGIN' in line:
                indexStart = index
            elif 'END' in line:
                count += 1
                indexEnd = index
                structlines = lines[indexStart:indexEnd + 1]
                if count %1000 == 0:
                    print("Processed {} crystals".format(nCrystals))
                nCrystals += 1

                thisCrystal = Crystal.from_lines(structlines,species,'mlp')
                # Only add the crystal if mindist is reasonable
                if thisCrystal.minDist > 1.5:
                    # Only calculate formation energy if I have information about the pures


                    if True not in [x.crystal is None or x.crystal.results is None or thisCrystal.results["energyF"] is None for x in pures]:
                        thisCrystal.results["fEnth"] = thisCrystal.results["energyF"]/thisCrystal.nAtoms - sum(   [ pures[i].crystal.results["energyF"]/pures[i].crystal.nAtoms * thisCrystal.concentrations[i] for i in range(thisCrystal.nTypes)])

                    # Otherwise, set it to a ridiculus number
                    else:
                        thisCrystal.results["fEnth"] = 1000
                    thisCrystal.results["distToHull"] = None
                    # Save the crystal for later.
                    crystals.append(thisCrystal)
                else:
                    msg.warn("Mindist is pretty small for this one, so I'm not gonna add it")

        return dataset(crystals)

    @staticmethod
    def from_paths(paths,systemSpecies,calculator='VASP'):
        from aBuild.database.crystal import Crystal
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.aflow import AFLOW

        from aBuild.calculators.lammps import LAMMPS
        from os import path

        # Get pures information
        puredirs = [path.join(path.split(paths[0])[0],'pure' + x) for x in systemSpecies]
        pures = [AFLOW.from_path(x,systemSpecies,filesuffix = '.relax2.xz') for x in puredirs]
        for pure in pures:
            pure.read_results()
        crystals = []

        for dirpath in paths:
            if calculator == 'VASP' or calculator == 'AFLOW':
                if path.isfile(path.join(dirpath,'aflow.in')):
                    print("Initializing AFLOW object from path: {}".format(dirpath))
                    calc = AFLOW.from_path(dirpath,systemSpecies)
                    calc.read_results(pures = pures)
                    if calc.crystal.results is not None:
                        print("Formation energy is {}.".format(calc.crystal.results["fEnth"]))
                else:
                    calc = VASP.from_path(dirpath,systemSpecies)
                    calc.read_results(pures=pures)


            if calc.crystal is not None and calc.crystal.results is not None:
                print("Adding Crystal")
                crystals.append(calc.crystal)
        species = crystals[0].systemSpecies
        return dataset(crystals)


    def starting_point(self,folderpath):
        from os import path
        from glob import glob
        from aBuild.utility import chdir

        with chdir(folderpath):
            dirsE = glob('E.*')
            dirsA = glob('A.*')
        prevCalcs = [int(x.split('.')[1])  for x in dirsE] + [int(x.split('.')[1])  for x in dirsA]
        prevCalcs.sort()
        if prevCalcs != []:
            return prevCalcs[-1] + 1
        else:
            return 1

        #    def write(self):


    def buildFolders(self,buildpath,calculator,runGetKpoints = True,foldername = 'A',onlyCloseToHull = False,distToHull = 5e-3):
        from os import path
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.aflow import AFLOW
        from aBuild.calculators.lammps import LAMMPS
        from aBuild.calculators.espresso import ESPRESSO
        from aBuild.jobs import Job
        from math import floor
        import os

        print("Building folders in {}".format(buildpath))
        if not path.isdir(buildpath):
                os.mkdir(buildpath)
                print('Made path:',buildpath)
        configIndex = startPoint = self.starting_point(buildpath)

        lookupCalc = {'aflow': lambda specs: AFLOW.from_dictionary(specs),
                      'vasp': lambda specs: VASP.from_dictionary(specs),
                      'qe': lambda specs: ESPRESSO(specs,self.species),
                      'lammps': lambda specs: LAMMPS(specs,self.species)}

        lookupBuild = {'aflow': lambda obj: obj.buildFolder(),
                       'vasp': lambda obj: obj.buildFolder(runGetKPoints = runGetKpoints),
                       'qe': lambda obj:obj.buildFolder(),
                      'lammps': lambda obj: obj.buildFolder()}

        for crystal in self.crystals:
            if onlyCloseToHull:
                if crystal.results["distToHull"] is None:
                    msg.fatal("You asked only for cystals that are close to the hull, but I don't have a value for that distance.")
                elif crystal.results["distToHull"] > distToHull:
                    continue
            print("Building crystal {}".format(crystal.title))
            runpath = path.join(buildpath,foldername + ".{}".format(configIndex) )
            #Augment the existing dictionary in preparation for sending it in
            calculator[calculator["active"]]["crystal"] = crystal
            calculator[calculator["active"]]["directory"] = runpath

            # Initialize the calculation object
            print('initializing VASP object')
            thisCalc = lookupCalc[calculator["active"]](calculator[calculator["active"]])

            # Build the path
            if not path.isdir(runpath):
                os.mkdir(runpath)
            else:
                msg.fatal("I'm gonna write over top of a current directory. ({})  I think I'll stop instead.".format(runpath))

                # Change the directory and build the folder
            print("Building folder for structure: {}".format(crystal.title) )
            with chdir(runpath):
                success = lookupBuild[calculator["active"]](thisCalc)
            if not success:
                if calculator["active"] == 'aflow':
                    retryCalc = lookupCalc["vasp"](calculator['vasp'])
                    with chdir(runpath):
                        success = lookupBuild["vasp"](retryCalc)
                    if not success:
                        msg.fatal("I tried building an aflow dir and it failed, then I tried building a VASP dir and it failed too. I give up")
                else:
                    msg.warn("VASP(??) directory build failed, and I'm not sure why")

            configIndex += 1


        # Build the submission script
        exdir = path.join(buildpath,'A.')
        if calculator['active'] == 'aflow':
            calculator["execution"]["exec_path"] = "aflow --run"
        elif calculator["active"] == 'vasp':
            calculator["execution"]["exec_path"] = "vasp6_serial"


        startAdder = int(floor(startPoint/1000)) * 1000
        endAdder = int(floor((configIndex - 1)/1000)) * 1000

        if startAdder == endAdder:  # Don't need to submit two jobs in this case.  Just one, but we might have to add an offset if the numbers are too high.
            msg.info("Building one job submission file")
            calculator["execution"]["offset"] = startAdder
            mljob = Job(calculator["execution"],exdir,calculator["execution"]["exec_path"], arrayStart = startPoint-startAdder,arrayEnd = configIndex - 1 - endAdder)
            with chdir(buildpath):
                print('Building job file')
                mljob.write_jobfile('jobscript_vasp.sh')
        else:  # We're going to have to submit two jobs to span the whole job array.
            msg.info("Building two job submission files")
            #First job..
            calculator["execution"]["offset"] = startAdder
            mljob = Job(calculator["execution"],exdir,calculator["execution"]["exec_path"], arrayStart = startPoint - startAdder,arrayEnd = 999)
            with chdir(buildpath):
                print('Building job file')
                mljob.write_jobfile('jobscript_vasp_1.sh')

            calculator["execution"]["offset"] = endAdder - 1
            mljob = Job(calculator["execution"],exdir,calculator["execution"]["exec_path"], arrayStart = 1,arrayEnd = configIndex - endAdder)
            with chdir(buildpath):
                print('Building job file')
                mljob.write_jobfile('jobscript_vasp_2.sh')




    def writeReport(self,dset):
        import datetime
        nAtoms = len(self.crystals[0].crystalSpecies)
        print(' '.join(self.crystals[0].crystalSpecies))
        print(type(' '.join(self.crystals[0].crystalSpecies)))
        with open('dataReport_' + dset + '.txt', 'w') as f:
            f.write(dset + ' REPORT\n')
            f.write(str(datetime.datetime.now()) + '\n')
            f.write("{:54s} {:14s}{:13s}{:12s}{:27s}{:21s}{:10s}".format("Title"," T. Energy","Enery/Atom","F. Energy","Distance To Hull","Conc.",'      '.join(self.crystals[0].crystalSpecies) + '\n'))
            f.write('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            for crystal in self.crystals:
                f.write(crystal.reportline)


    def ternaryXYcoords(self,concentrations):
        return concentrations[1] * 0.5 + concentrations[2],concentrations[1] * 0.5 * tan(60*pi/180)



    def findAllHullDistances(self):

        index = 1
        for crystal in self.crystals:
            print('crystal number',index)
            crystal.results["distToHull"] = self.distanceToHull(crystal)
            index += 1


    # This is specific to ternary systems!!!!!
    def distanceToHull(self,crystal):
        from numpy import tan, pi,array,sqrt
        from numpy.linalg import norm
        fEnth = crystal.results["fEnth"]

        xCoord = crystal.concentrations[1] * 0.5 + crystal.concentrations[2]
        yCoord = crystal.concentrations[1] * 0.5 * tan(60*pi/180)

        currentPoint = array([xCoord,yCoord])
        if True in [norm(currentPoint - r) < 1e-4 for r in [array([0,0]),array([1,0]),array([0.5,sqrt(3)/2])]]:
            print('Found pure substance-----------------------------------------------------------')
            return fEnth
#        print('Finding hull simplex for point: ({},{})'.format(xCoord,yCoord))
        foundSimplex = False
        index = 1
        for simplex in self.cHull:
            pointOne = array([self.data[simplex[0]][0],self.data[simplex[0]][1],self.data[simplex[0]][2]])
            pointTwo = array([self.data[simplex[1]][0],self.data[simplex[1]][1],self.data[simplex[1]][2]])
            pointThree = array([self.data[simplex[2]][0],self.data[simplex[2]][1],self.data[simplex[2]][2]])

            if pointOne[2] == 0 and pointTwo[2] == 0 and pointThree[2] == 0:
                # We know every point is in this simplex, which is the
                # top of the convex surface.  We only care about the bottom
                # side of the hull.
                print("passing..  This shouldn't happen")
                index += 1
                continue


            if pointInPolygon(currentPoint,[pointOne[:2],pointTwo[:2],pointThree[:2]]):
                print('found simplex')
                print('simplex number',index)
                index += 1
                #print(simplex)
                foundSimplex = True
                break

            index += 1
        if not foundSimplex:
#            return None
            msg.fatal("Can't find simplex for this crystal.")



        xOne = pointOne[0]
        yOne = pointOne[1]
        zOne = pointOne[2]

        xTwo = pointTwo[0]
        yTwo = pointTwo[1]
        zTwo = pointTwo[2]

        xThree = pointThree[0]
        yThree = pointThree[1]
        zThree = pointThree[2]

        planeFunc = lambda x,y: (x * yTwo * zOne - x * yThree * zOne + xOne * y * zTwo - x * yOne *  zTwo + x * yThree * zTwo - xOne * yThree * zTwo + xThree * (y * zOne - yTwo * zOne - y *zTwo + yOne * zTwo) - xOne * y * zThree + x * yOne * zThree - x * yTwo * zThree +xOne * yTwo * zThree + xTwo * (yThree * zOne - yOne * zThree + y * (zThree - zOne)) )/(xThree * (yOne - yTwo) + xOne * (yTwo - yThree) + xTwo * (yThree - yOne))

        distance = fEnth - planeFunc(xCoord,yCoord)
        print("Distance to hull is {} eVs. Crystal f Enth: {}, Hull Value: {}".format(distance,fEnth,planeFunc(xCoord,yCoord)))
        return distance

    def deleteTopOfHull(self):
        from numpy import sqrt,array,tan,pi

        self.cHull = []
        for simplex in self.hull.simplices:
            r1 = array([self.data[simplex[0]][0],self.data[simplex[0]][1],self.data[simplex[0]][2]])
            r2 = array([self.data[simplex[1]][0],self.data[simplex[1]][1],self.data[simplex[1]][2]])
            r3 = array([self.data[simplex[2]][0],self.data[simplex[2]][1],self.data[simplex[2]][2]])
            if all(array([r1[2],r2[2],r3[2]]) == 0):
                # Don't include the top surface
                continue
            if r1[0] * tan(60 * pi/180) - r1[1] < 1e-4 and r2[0] * tan(60 * pi/180) - r2[1] < 1e-4 and (r3[0] * tan(60 * pi/180) - r3[1]) < 1e-4:
                # Don't include vertical plane
                continue
            if all(array([r1[1],r2[1],r3[1]]) == 0):
                # Don't include vertical plane
                continue
            if (-r1[0] * tan(60 * pi/180) + sqrt(3) - r1[1]) < 1e-4 and (-r2[0] * tan(60 * pi/180) + sqrt(3) - r2[1]) < 1e-4 and (-r3[0] * tan(60 * pi/180) + sqrt(3) - r3[1]) < 1e-4:
                # Don't include vertical plane
                continue
            self.cHull.append(simplex)

    def generateConvexHull(self,plot2D = False,plot3D=False,fileOutput = False):
        from scipy.spatial import ConvexHull
        from numpy import array,append,tan,pi,allclose
        from matplotlib import pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib
        import itertools

        fig2d = plt.figure(1)
        fig3d = plt.figure(2)
        ax2d = fig2d.add_subplot(111)
        ax3d = fig3d.add_subplot(111,projection='3d')

        reducedDataSet = [x for x in self.crystals if x.results["fEnth"] <= 0 and x.concentrations not in [[0.0,0.0,1.0],[0.0,1.0,0.0],[1.0,0.0,0.0]]]

        # Only include data for negative formation enthalpies and delete all of the pure crystals.
        self.data = [[x.concentrations[1] * 0.5 + x.concentrations[2],x.concentrations[1] * 0.5 * tan(60*pi/180),x.results["fEnth"]] for x in reducedDataSet]

        # Add the pures back in and set their formation energies to 0
        self.data.append([0.0,0.0,0.0])
        self.data.append([1.0,0.0,0.0])
        self.data.append([0.5,0.5 * tan(60 * pi/180),0.0])

        self.hull = ConvexHull(self.data)

        # Delete top planes and vertical planes
        self.deleteTopOfHull()
        order = 3

        filewrite = open('cHullpts.out','w')
        included = []
        for simplex in self.cHull:
            xlist = []
            ylist = []
            zlist = []
            for point in simplex:
                x = self.data[point][0]
                y = self.data[point][1]
                z = self.data[point][2]
                xlist.append(x)
                ylist.append(y)
                zlist.append(z)
                if point not in included and point < len(reducedDataSet):
                    filewrite.write("{:50s}  {:6.3f}   {:6.3f}   {:6.3f}   {:6.3f}\n".format(reducedDataSet[point].title, reducedDataSet[point].concentrations[0],reducedDataSet[point].concentrations[1],reducedDataSet[point].concentrations[2], reducedDataSet[point].results["fEnth"]) )
                included.append(point)
            xlist.append(xlist[0])
            ylist.append(ylist[0])
            zlist.append(zlist[0])
            if plot3D:
                ax3d.plot(xlist,ylist,zlist)
            if plot2D:
                #ax.scatter(xlist,ylist,zlist,marker = '+')
                ax2d.plot(xlist,ylist,marker = '+')#,[data[i[0]][2],data[i[1]][2],data[i[2]][2]]


        if plot2D:
            print('saving figure')
            cm = plt.cm.get_cmap('autumn')
            plottwoD = ax2d.scatter([float(x[0]) for x in self.data],[float(x[1]) for x in self.data],c=[float(x[2]) for x in self.data],marker= '+',cmap = cm)
            #plt.autumn()
            fig2d.colorbar(plottwoD)
            fig2d.savefig('convex_hull_2d.png')
        if plot3D:
            ax3d.view_init(elev = 20,azim = -85)
           # ax3d.scatter([float(x[0]) for x in self.data],[float(x[1]) for x in self.data],[float(x[2]) for x in self.data],marker= '+')
            fig3d.savefig('convex_hull.png')
        filewrite.close()
      #  import sys
      #  sys.exit()
#        import sys
#        sys.exit()
#        #        pyplot.plot(self.concs,self.formationenergies,'r+')
#        plotConcs = []
#        plotEnergies = []
##        pyplot.plot(data[hull.vertices,0], data[hull.vertices,1],'b-',lw = 2)
#        print(self.hull.vertices, 'verts')
#        print(len(self.data))
#        print(len(self.titles))
#        print([self.titles[x] for x in vertices])
#
#        print(vertices,' verts')
#        pyplot.figure(figsize = (15,10))
#        if plotAll:
#            pyplot.plot([x[0] for x in self.data],[x[1] for x in self.data],'r+')
#        for ivert,vertex in enumerate(vertices):
#            if self.data[vertex,1] <= 0:
#                plotConcs.append(self.data[vertex,0])
#                plotEnergies.append(self.data[vertex,1])
#        pyplot.plot(plotConcs,plotEnergies,'xk-',linewidth=2.3,markersize = 8)
#        font = {'family':'normal',
#                    'weight': 'bold',
#                    'size': 22}
#        matplotlib.rc('font',**font)
#        pyplot.xlabel(' Ag', fontsize = 24)
#        pyplot.ylabel("Formation Energy (eV/atom)", fontsize = 24)
#        pyplot.xticks(fontsize=22)
#        pyplot.yticks(fontsize=22)
#        pyplot.title("Convex Hull Plot")
#        pyplot.savefig('chull.png')
