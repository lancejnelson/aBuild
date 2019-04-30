from os import path
from aBuild import msg # messaging module
from aBuild.utility import chdir
import os
import sys
import numpy as np

config = sys.modules["config"]  

class ESPRESSO:
    """Class to handle all of the Quantum Espresso input and output files.
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


    def __init__(self, specs,systemSpecies,directory = None): #crystal,potential,

        from aBuild.database.crystal import Crystal


#        self.crystal = crystal
        self.species = systemSpecies
#        self.potential = potential
        if isinstance(specs,dict):
            self.crystal = specs["crystal"]
            self.pseudopotentials = specs["pseudopotentials"]
        elif isinstance(specs, str):
            self.crystal = Crystal(path.join(specs,'input.in'),systemSpecies)
            
            self.settings = self.parse_lammps_settings(specs)
            self.directory = specs
    def parse_lammps_settings(self,specs):
        filePath = path.join(specs,'input.in')
        with open(filePath, 'r') as f:
            lines = f.readlines()

        self.title = lines[1]
        for line in lines:
            if 'pair_coeff' in line:
                self.potential = line.split()[3]
            if 'units' in line:
                self.units = line.split()[1]
            

        

    def buildFolder(self,runGetKPoints = True):
        import shutil
        from jinja2 import Environment, PackageLoader  # Package for building files from a template
        from os import path
        import aBuild
        from aBuild.calculators.vasp import KPOINTS
                                                                 
        medpath = path.abspath(aBuild.__file__)
        reporoot = path.dirname(path.dirname(medpath))

        settings = {}
        settings["latpar"] = self.crystal.latpar
        settings["nAtoms"] = self.crystal.nAtoms
        settings["knary"] = self.crystal.nTypes
        settings["lVs"] = self.crystal.lattice_lines
        settings["bVs"] = self.crystal.basis_lines_ESPRESSO
        settings["directory"] = self.pseudopotentials["directory"]
        settings["pseudopotentials"] = self.pseudoLines
        kptsdict = {"method":"MP", "mindistance":2000}
        kPts = KPOINTS(kptsdict)
        settings["kPoints"] = ' '.join(map(str,kPts.monkPack(self.crystal.recip_Lv,self.crystal.nAtoms) ) )
        settings["title"] = self.crystal.title
        env = Environment(loader=PackageLoader('aBuild', 'templates'))
        template = env.get_template("qe.in")

        print('here')
        with open("qeinput.in",'w') as f:
            f.write(template.render(**settings))
#        shutil.copy(path.join(reporoot,"aBuild","templates",self.potential),'.')

    @property
    def pseudoLines(self):
        lines = ''
        for i in range(self.crystal.nTypes):
            lines += self.crystal.species[i] + ' 1.0 ' + self.pseudopotentials["versions"][self.crystal.species[i]] + '\n'
        return lines[:-1]
    def read_results(self):
        outFilePath = path.join(self.directory,'lammps.out')
        infile = open(outFilePath,'r')
        self.crystal.results = {}
        for line in infile.readlines():
            if 'Cohesive energy' in line:
                self.crystal.results["energy"] = float(line.split()[4].split(';')[0])

        self.pures = []
        if 'pure' not in self.directory:
            for i in self.crystal.species:
                purepath = path.join(path.split(self.directory)[0],'pure' + i)
                
                pure = LAMMPS(purepath,self.species)
                pure.read_results()
                self.pures.append(pure)

        self.crystal.results["formation enthalpy"] = self.crystal.results["energy"] -   sum([self.pures[x].crystal.results["energy"] *  self.crystal.atom_counts[x] for x in range(len(self.pures))])/sum(self.crystal.atom_counts )
#        from aBuild.calculators.vasp import POSCAR
#
#        self.KPOINTS.rGP = runGetKPoints
#        self.INCAR.writeINCAR()
#        print("INCAR built")
#        self.crystal.write('POSCAR_orig')
#        print("POSCAR_orig built")
#        self.check_atom_counts_zero()
#        self.crystal.write('POSCAR')
#        print("POSCAR built")
#        self.KPOINTS.writeKPOINTS()
#        
#        print("KPOINTS built")
#        self.POTCAR.writePOTCAR()
#        print("POTCAR built")

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
                    print(lines[n+nAtoms + 2])
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

    
        
