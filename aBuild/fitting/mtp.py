from aBuild.utility import chdir


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
        relaxDict["extrap_threshold"] = 2.0
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
        with open(target,'w') as f:
            f.write(template.render(**settings))


    def train(self,executeParams, potential="pot.mtp", tSet="train.cfg",buildJob = True):

        if buildJob:
            if executeParams["ntasks"] > 1:
                mlpCommand = 'mpirun -n ' + executeParams["ntasks"] + ' mlp train {} {}'.format(potential,tSet)
            else:
                mlpCommand = 'mlp train {} {}'.format(potential,tSet)
            mljob = Job(executeParams,path.join(self.root,"fitting","mtp"),mlpCommand)
            with chdir(path.join(self.root,"fitting/mtp")):
                print('Building job file')
                mljob.write_jobfile()
        else:
            with chdir(self.root):
                child=Popen(mlpCommand, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)
            


    def calc_grade(self):
        from subprocess import Popen
        from os import waitpid, rename,path
        
        print('renaming Trained.mtp_ to pot.mtp')
        rename(path.join(self.root,'Trained.mtp_'),path.join(self.root,'pot.mtp'))
        mlpCommand = 'mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg'

        print('Running calc-grade')
        with chdir(self.root):
            child=Popen(mlpCommand, shell=True, executable="/bin/bash")
            waitpid(child.pid, 0)
        print('Done')

    def relax(self,executeParams,buildJob = True):
        from subprocess import Popen
        from os import waitpid, rename,path
        
        mlpCommand = 'mlp relax relax.ini --cfg-filename = to-relax.cfg --save-relaxed=relaxed.cfg --save-unrelaxed=unrelaxed.cfg --log=relax_log.txt'

        if buildJob:
            mljob = Job(executeParams,path.join(self.root,"fitting","mtp"),mlpCommand)
            with chdir(path.join(self.root,"fitting/mtp")):
                print('Building job file')
                mljob.write_jobfile()
        else:
            with chdir(self.root):
                child=Popen(mlpCommand, shell=True, executable="/bin/bash")
                waitpid(child.pid, 0)

                #    def select_add(self,executeParams
