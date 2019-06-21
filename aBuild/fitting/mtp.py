from aBuild.utility import chdir
from aBuild import msg

class MTP(object):
	


    def __init__(self,root,settings = None):
        from os import makedirs,path
        if not path.isdir(root):
            makedirs(root)

        self.root = root
        self.settings = settings
        #        self.trainingSet = trainingSet  # Should be a list of Crystal objects

    @property
    def relaxDefaults(self):
        relaxDict = {}
        relaxDict["calc_efs"] = 'TRUE'
        relaxDict["fit_setting"] = 'FALSE'
        relaxDict["active_learn"] = 'TRUE'
        relaxDict["site_weight"] = 0.0
        relaxDict["energy_weight"] = 1.0
        relaxDict["force_weight"] = 0.001
        relaxDict["stress_weight"] = 0.0001
        relaxDict["extrap_threshold"] = 5.0
        relaxDict["threshold_break"] = 10.0
        relaxDict["efs_ignore"] = 'FALSE'
        return relaxDict



    def singleStructureLines(self,calc):
        import numpy as np
        if not self.forRelax and calc.results is None:
            return []
        lines = []
        lines.append('BEGIN_CFG\n')
        lines.append('Size\n')
        lines.append(str(calc.crystal.nAtoms) + '\n')
        lines.append('SuperCell\n')
        for lv in calc.crystal.latpar * calc.crystal.lattice:
            lines.append('{:12.6f} {:12.6f} {:12.6f}\n'.format(lv[0],lv[1],lv[2]  ))
        if not self.forRelax:
            lines.append('   AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n')
        else:
            lines.append('   AtomData:  id type       cartes_x      cartes_y      cartes_z\n')
        #counter = crystal.lattice.nTypes - 1
        
        # Took me a few minutes to figure this one out.  Very Pythonic
        #        atomLabels = [ x for sublist in [ [ counter - i for k in range(crystal.lattice.atom_counts[i]) ]   for i in range(counter + 1)] for x in sublist]
        atomLabels = [ x for sublist in [ [  i for k in range(calc.crystal.atom_counts[i]) ]   for i in range(calc.crystal.nTypes)] for x in sublist]
        for i in range(calc.crystal.nAtoms):
            if not self.forRelax:
                forces = calc.results["forces"][i]
            coords = calc.crystal.Bv_cartesian[i]
            if not self.forRelax:
                lines.append('{:16d} {:3d} {:16.6f} {:12.6f} {:12.6f} {:18.6f} {:10.6f} {:10.6f}\n'.format(i+ 1,atomLabels[i], coords[0],coords[1],coords[2], forces[0],forces[1],forces[2]  ))
            else:
                lines.append('{:16d} {:3d} {:16.6f} {:12.6f} {:12.6f}\n'.format(i+ 1,atomLabels[i], coords[0],coords[1],coords[2] ))
            #line +=  ' '.join([map(str,crystal.lattice.Bv_cartesian[i]), crystal.forces[i]])

        if not self.forRelax:
            lines.append('Energy\n')
            lines.append(str(calc.results["energyF"]) + '\n')
        
        
            lines.append(' Stress:   xx          yy           zz            yz           xz           xy\n')
            s = calc.results["stress"]
            stressesline = '{:16.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n'.format(s[0],s[1],s[2],s[3],s[4],s[5])
            lines.append(stressesline)
        lines.append(''.join([' Feature   conf_id ', '  '.join([calc.crystal.symbol,calc.crystal.title,calc.directory]),'\n']))
        lines.append('END_CFG\n')
        
        
        return lines
    @property
    def lines(self):
        lines = []
        for calc in self.dataSet:
            lines += self.singleStructureLines(calc)
        return lines
        

    def write_relaxin(self):

        from os import path
        settings = self.relaxDefaults
        from jinja2 import Environment, PackageLoader  # Package for building files from a template
        env = Environment(loader=PackageLoader('aBuild', 'templates'))
        template = env.get_template("relax.ini")

        target = path.join(self.root, "relax.ini")
        with open(target,'w') as f:
            f.write(template.render(**settings))

    def write_blank_pot(self,nSpecies):
        from os import path
        settings = {}
        settings["n_species"] = nSpecies
        from jinja2 import Environment, PackageLoader  # Package for building files from a template
        env = Environment(loader=PackageLoader('aBuild','templates'))
        print(self.settings)
        template = env.get_template(self.settings['pot']) #"pot.mtp"


        target = path.join(self.root, "pot.mtp")
        if not path.isfile(target):
            with open(target,'w') as f:
                f.write(template.render(**settings))
        else:
            msg.info("File pot.mtp exists already.  You must have copied over Trained.pot_ previously")

    def train(self,executeParams, potential="pot.mtp", tSet="train.cfg",buildJob = True):
        from aBuild.jobs import Job
        from os import path
        if buildJob:
            if executeParams["ntasks"] > 1:
                mlpCommand = 'mpirun -n ' + str(executeParams["ntasks"]) + ' mlp train {} {}'.format(potential,tSet)
            else:
                mlpCommand = 'mlp train {} {}'.format(potential,tSet)
            mljob = Job(executeParams,self.root,mlpCommand)
            with chdir(self.root):
                print('Building job file')
                mljob.write_jobfile('jobscript_train.sh')
        else:
            with chdir(self.root):
                child=Popen(mlpCommand, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)
            


    def calc_grade(self):
        from subprocess import Popen
        from os import waitpid, rename,path
        
        mlpCommand = 'mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg'

        print('Running calc-grade')
        with chdir(self.root):
            child=Popen(mlpCommand, shell=True, executable="/bin/bash")
            waitpid(child.pid, 0)
        print('Done')

    def relax(self,executeParams,buildJob = True):
        from aBuild.jobs import Job
        from subprocess import Popen
        from os import waitpid, rename,path
        
        baseCommand = 'mlp relax relax.ini --cfg-filename=to_relax.cfg --save-relaxed=relaxed.cfg --save-unrelaxed=unrelaxed.cfg --log=relax_log.txt'

        if buildJob:
            if executeParams["ntasks"] > 1:
                mlpCommand = 'mpirun -n ' + str(executeParams["ntasks"]) + ' ' + baseCommand
            else:
                mlpCommand = baseCommand
            mljob = Job(executeParams,self.root,mlpCommand)
            with chdir(self.root):
                print('Building job file')
                mljob.write_jobfile('jobscript_relax.sh')
        else:
            with chdir(self.root):
                child=Popen(mlpCommand, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)

    def select_add(self,executeParams,buildJob = True):
        from subprocess import Popen
        from os import waitpid, rename,path
        from aBuild.jobs import Job
        baseCommand = 'mlp select-add pot.mtp train.cfg candidate.cfg new_training.cfg'

        if buildJob:
            if executeParams["ntasks"] > 1:
                mlpCommand = 'mpirun -n ' + str(executeParams["ntasks"]) + ' ' + baseCommand
            else:
                mlpCommand = baseCommand
            mljob = Job(executeParams,self.root,mlpCommand)
            with chdir(self.root):
                print('Building job file')
                mljob.write_jobfile('jobscript_select.sh')
        else:
            with chdir(self.root):
                child=Popen(mlpCommand, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)
                #    def select_add(self,executeParams



    # This routine pulls all of the VASP data in the training_set folder
    # and puts it into a folder called train.cfg in fitting/mtp
    # Additionaly, the following files are written to the same folder:
    #                   1- a blank MTP potential file, ready to be fit to
    #                   2- a jobscript that can be submitted to train the MTP.
    # Pull data from VASP calculations and get ready to fit.
    def setup_train(self,trainingRoot,species):
        from os import path, remove
        from glob import glob
        from aBuild.calculators.vasp import VASP
        from aBuild.fitting.mtp import MTP
        from aBuild.jobs import Job

        from aBuild.database.dataset import dataset


        # 1. Pull all VASP data and compile into file train.cfg
        # 2. Copy/Build a blank potential file called pot.mtp
        # 3. Build a submission script

        with chdir(trainingRoot):
            enumdirs = glob("E.*")
           # activedirs = []
            activedirs = glob("A.*")
            #puredirs = glob("pure*")

        dirs = [path.join(trainingRoot,x) for x in enumdirs + activedirs]
        if dirs != []:
            
            print('Building dataset')
            trainingSet = dataset(dirs,species,calculator='VASP')  #####

        if path.isfile(path.join(self.root,'train.cfg') ):
            msg.info('train.cfg file found!  Deleting and building from scratch.')
            remove(path.join(self.root,'train.cfg'))
            
        with open(path.join(self.root,'train.cfg'),'a+') as f:
            for crystal in trainingSet.crystals:
                if crystal.results["fEnth"] < 100 and crystal.minDist > 1.5:
                    print("Building: {}".format(crystal.title))
                    f.writelines('\n'.join(crystal.lines('mtptrain')))
                else:
                    print("Not adding structure {}.  Seems like an extreme one.".format(crystal.title))
                    print("Energy: {}".format(crystal.results["energyF"]))
                    print("MinDist: {}".format(crystal.minDist))
                    


        self.write_blank_pot(len(species))
        self.train(self.settings["execution"])

        
    # This routine sets up files needed to perform a relaxation using
    # the MTP that was just built.
    # The file that contains the structures to be relaxed is called to-relax.cfg.
    # If that file isn't present, this routing will build it up (takes a while)
    # from the enumeration and prototypes databases.
    # If it is present, we'll just leave it alone.  (This file should never be disturbed
    # once it's built)
    # Additionally, the following actions are performed:
    #               1- The trained potential (Trained.mtp) is copied to pot.mtp.
    #               2- A relax.ini input file is generated.
    #               3- mlp calc-grade is performed in preparation for the relaxation.
    #                   We just run this interactively because it doesn't take too long.
    #               4- A job submission script is generated.
    def setup_relax(self,enumDicts,species,freshStart = False,AFM = False,start = 1,end = None):
        from os import remove,path,rename
        from aBuild.fitting.mtp import MTP
        from glob import glob
        from shutil import copy
        self.dataset = "gss"


        # 1.  If to-relax.cfg does not exist, build it.  This
        #     takes a while because we are putting a bunch of
        #     structures in there.
        # 2.  If to-relax.cfg does exist (iteration > 1), then
        #     copy unrelaxed.cfg from previous iteration to to-relax.cfg
        #     Note that there may be many unrelaxed.cfg_# from the parallel run
        # 3.  Copy Trained.mtp to pot.mtp in preparation for the relaxation
        # 4.  Build a relax.ini from template if it's not already there.
        # 5.  Run calc-grade
        # 6.  Build a submission script.
        filePath = path.join(self.root,'to_relax.cfg')
        #unrelaxedfilePath = path.join(fittingRoot,'unrelaxed.cfg')

        torelax = path.isfile(filePath)
        #unrelaxed = path.isfile(unrelaxedfilePath)

        # 1.
        if not torelax:  #This must be the first iteration
            print('Building a new to-relax.cfg file')
            self.build_ToRelax(enumDicts,species,AFM = AFM,start = start,end = end)

        # 2.
        else: # What iteration is it?
            if freshStart:
                remove(filePath)
                self.build_ToRelax(enumDicts,species,AFM = AFM,start = start,end = end)
         #   elif unrelaxed:  #Iteration must be > 1
         ##       unrelaxedFiles = glob(path.join(fittingRoot,"unrelaxed.cfg_*"))
          ##      cat(unrelaxedFiles,path.join(fittingRoot,"to-relax.cfg"))
           #     remove(unrelaxedFiles)
            else: #"to-relax.cfg is present and I don't wanta  fresh start.  Don't do anything."
            #    self.build_ToRelax(enumDicts,species)
                msg.info("to-relax.cfg found and you told me to not start fresh.  Proceeding with the to-relax.cfg that's there.")
        #elif path.isfile(filePath):  # Iteration # must be > 1
        #        rename(path.join(self.root,'fitting/mtp')
        

        # 3.
        try:
            print('copying Trained.mtp to pot.mtp')
            copy(path.join(self.root,'Trained.mtp_'),path.join(self.root,'pot.mtp'))
        except:
            if not path.isfile(path.join(self.root,'Trained.mtp_')):
                msg.info("Can't find Trained.mtp_ .You don't appear to have a trained potential ready" )
            else:
                msg.info("Not sure why the copy of Trained.mtp_ -> pot.mtp failed.  You should investigate")
                # if path.isfile(path.join(self.root,'pot.mtp')):
           #     msg.info('It looks like you have already copied Trained.mtp -> pot.mtp previously')
       #     else:
       #         msg.fatal("Can't find a pot.mtp or a Trained.mtp.  Problems")
                          
        # 4.
        self.write_relaxin()
        # 5.
        self.calc_grade()
        # 6.
        self.relax(self.settings["execution"])



        #        trainingSet = dataset(enumspecs = self.enumDicts,root = self.root,species = self.species,calculator = self.calculator)
        #trainingSet.build_relax_select_input()
        


    def build_ToRelax(self,enumDicts,species,AFM = False,start = 1,end = None):
        from aBuild.enumeration import Enumerate
        from aBuild.utility import unpackProtos,getAllPerms,getProtoPaths
        from aBuild.database.crystal import Crystal
        from os import remove,path
        print('Building to-relax.cfg')

        nEnums = len(enumDicts)
        knary = len(species)
        for ilat  in range(nEnums):
            lat = enumDicts[ilat]["name"]
            
            if lat == 'protos':
                structures = getProtoPaths(knary)

                for struct in structures:
                    print("Proto structure:", struct)
                    scrambleOrder = getAllPerms(knary,justCyclic = 'uniqueUnaries' in struct)
                    for scramble in scrambleOrder:
                        thisCrystal = Crystal(struct,species)

                        #print("Atom counts before scramble {}".format(thisCrystal.atom_counts))
                        thisCrystal.scrambleAtoms(scramble)
                        if not thisCrystal.concsOK(concRestrictions = enumDicts[ilat]["concs"]):
                            continue
                        print(thisCrystal.title)
                        mindist = thisCrystal.minDist
                        print(mindist, "actual min dist<-----------------------------------------------------------------------------------------------------")

                        if mindist > 2 and thisCrystal.nAtoms < 60:
                            if not AFM:
                                print('Adding to file')
                                #print("Atom counts after scramble {}".format(thisCrystal.atom_counts))
                                with open(path.join(self.root,'to-relax.cfg'),'a+') as f:
                                    f.writelines('\n'.join(thisCrystal.lines('mtprelax') ) )
                            elif thisCrystal.getAFMPlanes([1,0,0]):
                                print("Original Crystal is AFM compatible")
                                with open(path.join(self.root,'to-relax.cfg_' + str(start)),'a+') as f:
                                    f.writelines('\n'.join(thisCrystal.lines('mtprelax') ))
                              #  break
                            else:
                                print("Checking super-periodics")
                        
                                superCrystal = thisCrystal.superPeriodics(2)
                                if superCrystal != []:
                                    print('Found a super-Periodic that works')
                                    print(superCrystal.basis,'basis')
                                    print(superCrystal.atom_counts, 'atom counts')
                                    with open(path.join(self.root,'to-relax.cfg_' + str(start)),'a+') as f:
                                        f.writelines('\n'.join(superCrystal.lines('mtprelax') ))
                               # break
            else:
                enumLattice = Enumerate(enumDicts[ilat])
                if end == None:
                    end = enumLattice.nConfigs + 1
                    filetag = ''
                else:
                    filetag = '_' + str(start)
                #for struct in range(enumLattice.nConfigs+1):
                for struct in range(start, end):
                    print("Lattice",lat, "structure:",struct)
                    enumLattice.generatePOSCAR(struct) 
                    thisCrystal = Crystal(path.join(enumLattice.root,"poscar.{}.{}".format(lat,struct)),species)
                    if not AFM:
                        with open(path.join(self.root,'to_relax.cfg' + filetag),'a+') as f:
                            f.writelines('\n'.join(thisCrystal.lines('mtprelax') ))
                        
                    elif thisCrystal.getAFMPlanes([1,0,0]):
                        print("Original Crystal is AFM compatible")
                        with open(path.join(self.root,'to_relax.cfg' + filetag),'a+') as f:
                            f.writelines('\n'.join(thisCrystal.lines('mtprelax') ))
                    else:
                        print("Checking super-periodics")
                        
                        superCrystal = thisCrystal.superPeriodics(2)
                        if superCrystal != []:
                            print('Found a super-Periodic that works')
                            with open(path.join(self.root,'to_relax.cfg' + filetag),'a+') as f:
                                f.writelines('\n'.join(superCrystal.lines('mtprelax') ))
                            
                        
                    #                    print(thisCrystal.appMinDist,' approp Min Dist')
 #                   print(thisCrystal.minDist, 'actual min dist')
#                    with open(path.join(self.root,'to-relax.cfg'),'a+') as f:
#                        f.writelines('\n'.join(thisCrystal.lines('mtprelax') ))

                    delpath = path.join(enumLattice.root,"poscar.{}.{}".format(lat,struct))
                    remove(delpath)

    # This routine sets up for selection new structures to add to the
    # training set.  The main thing we need from the relaxation step are
    # the candidate.cfg_* files.  So we will concatenate them all into one
    # file and then remove all of the _* files.  Although we're not sure we'll
    # need them, we'll do the same with the following files:
    #                  1- selection.log_*
    #                  2- relaxed.cfg_*
    #                  3- unrelaxed.cfg_*
    #                  4- relax_log.txt_*
    def setup_select(self):
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

        # Find what iteration we are currently on.
        candIterationFiles = glob(path.join(self.root,"candidate_iteration_*"))
        relaxedIterationFiles = glob(path.join(self.root,"relaxed_iteration_*"))
        unrelaxedIterationFiles = glob(path.join(self.root,"unrelaxed_iteration_*"))
        if len(candIterationFiles) != len(relaxedIterationFiles) or len(candIterationFiles) != len(unrelaxedIterationFiles):
            msg.fatal(" I can't figure out what iteration we're on based on the files present.  Have a look at all of the *_iteration_* files")
        if candIterationFiles != []:
            
            iteration = len(candIterationFiles) + 1
        else:
            iteration = 1

            
        # 1. Bring all of the candidate.cfg_* files together into one
        candFiles = glob(path.join(self.root,"candidate.cfg_*"))
        if candFiles != []:
            cat(candFiles,path.join(self.root,"candidate.temp"),remove=True)
        # Save a record
            copy( path.join(self.root,'candidate.temp'), path.join(self.root,"candidate_iteration_" + str(iteration) + ".cfg") )
        # Prepare to select, which needs this file.
            copy( path.join(self.root,'candidate.temp'), path.join(self.root,"candidate.cfg" ) )
            remove(path.join(self.root,'candidate.temp'))
        else:
            msg.info("Can't find candidate.cfg_* files to concatenate")
        
        # 2. selection.log_*
        logFiles = glob(path.join(self.root,"selection.log_*"))
        if logFiles != []:
            cat(logFiles,path.join(self.root,"selection.temp"),remove=True)
            copy( path.join(self.root,'selection.temp'), path.join(self.root,"selection_iteration_" + str(iteration) + ".log") )
            remove(path.join(self.root,'selection.temp'))
        else:
            msg.info("Can't find selection.log_* files to concatenate")
        
        
        # 3. relaxed.cfg_*
        relaxedFiles = glob(path.join(self.root,"relaxed.cfg_*"))
        if relaxedFiles != []:
            cat(relaxedFiles,path.join(self.root,"relaxed.temp"),remove=True)
            copy( path.join(self.root,'relaxed.temp'), path.join(self.root,"relaxed_iteration_" + str(iteration) + ".cfg") )
            remove(path.join(self.root,'relaxed.temp'))
        else:
            msg.info("Can't find relaxed.cfg_* files to concatenate")

        # 4.  unrelaxed.cfg_*
        unrelaxedFiles = glob(path.join(self.root,"unrelaxed.cfg_*"))
        if unrelaxedFiles != []:
            cat(unrelaxedFiles,path.join(self.root,"unrelaxed.temp"),remove=True)
        # Save a record of this iteration.  May consider copying this file to to-relax.cfg next time we relax
            copy( path.join(self.root,'unrelaxed.temp'), path.join(self.root,"unrelaxed_iteration_" + str(iteration) + ".cfg") )
            remove(path.join(self.root,'unrelaxed.temp'))
        else:
            msg.info("Can't find unrelaxed.cfg_* files to concatenate")

        # 5. relax_log.txt_*
        relaxLogFiles = glob(path.join(self.root,"relax_log.txt_*"))
        if relaxLogFiles != []:
            cat(relaxLogFiles,path.join(self.root,"relax_log.temp"),remove=True)
            copy( path.join(self.root,'relax_log.temp'), path.join(self.root,"relax_log_iteration_" + str(iteration) + ".txt") )
            remove(path.join(self.root,'relax_log.temp'))
        else:
            msg.info("Can't find relax_log.txt_* files to concatenate")


        # 6.  Build job submission script.
        self.select_add(self.settings["execution"])
