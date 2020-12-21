
from six import string_types
from aBuild import msg
from aBuild import config

import sys
config = sys.modules["config"]

class Enumerate:

    def __init__(self,specs,arrows = None,name = "enum",overwrite=False):
        from aBuild.database.crystal import Lattice
        self.root = specs["root"]

        self.lattice = Lattice(specs)
        #        self._get_lattice(specs["lattice"])
        from os import path,makedirs
        if "Enum" not in self.root:
            self.root = path.join(self.root,"Enum")
            if not path.isdir(self.root):
                makedirs(self.root)
                msg.info('Making Enum directory.')
        self.inputFileExists = path.isfile(path.join(self.root, "struct_enum.in." + self.lattice.lattice_name))
        self.enumerationComplete = path.isfile(path.join(self.root, "struct_enum.out." + self.lattice.lattice_name))
        #self._get_basis(specs["basis"])
        self.knary = specs["knary"]
        self.sizeRange = specs["sizes"]
        self.eps = specs["eps"]
        self.nConfigs = specs["nconfigs"]
        self.latticeExpand = specs["latticeExpand"] #In case we want to increase/decrease atomic spacings
        self._processConcRestrictions(specs["concs"])

        self._processSiteRestrictions(specs["site_res"])

        if 'arrows' in specs:
            self.arrow_res = True
            self.arrows = specs["arrows"]
        else:
            self.arrow_res = False




    def setNEnumerated(self,nEnumerated):
        self.nEnumerated = nEnumerated


    @staticmethod
    def fromYAML(specs,index,knary,root):
        edict = {}
        edict["lattice"] = specs["lattice"][index]
        if 'basis' in specs.keys():
            edict["basis"] = specs["basis"][index]
        edict["knary"] = knary
        edict["nconfigs"] = specs["nconfigs"][index]
        edict["sizes"] = specs["sizes"][index]
        if "siteRestrictions" in specs.keys():
            edict["site_res"] = specs["siteRestrictions"][index]
        else:
            edict["site_res"] = None
        if 'coordsys' in specs.keys():
            edict["coordsys"] = specs["coordsys"][index]
        edict["name"] = specs["name"][index]
        if "concs" in specs.keys():
            edict["concs"] = specs["concs"][index]
        else:
            edict["concs"] = None
        edict["root"] = root
        edict["eps"] = 1e-3
        if 'latticeExpand' in specs.keys():
            edict["latticeExpand"] = specs["latticeExpand"][index]
        else:
            edict["latticeExpand"] = 1.0
        return Enumerate(edict)


    @staticmethod
    def from_enum_file(path):
        with open(path,'r') as f:
            lines = f.readlines()

    def _processConcRestrictions(self,concs):
        knarylookup = {2: "binary", 3: "ternary", 4:"quaternary"}
        if concs is not None:
            self.conc_res = True
            self.concRestrictions = concs
            if len(self.concRestrictions) > self.knary - 1:
                msg.warn("You have specified concentration restrictions on {} atom types.\
                For a {} system, you only need to supply restrictions on {} atom types. (<order of system> - 1) "\
                .format(len(self.concRestrictions),knarylookup[self.knary],self.knary - 1 ) )
        else:
            self.conc_res = False

    def _processSiteRestrictions(self,site_res):


        if site_res is None and self.knary is None and self.basis is None:
            self.site_res = None
            self.siteRestrictions = None
            return

        if site_res is not None:
            self.site_res = True
            self.siteRestrictions = site_res
        else:
            self.site_res = False
            self.siteRestrictions = ["/".join([str(i) for i in range(self.knary)]) for j in self.lattice.basis]

        if len(self.siteRestrictions) != self.lattice.nBasis:
            print(self.lattice.nBasis)
            print(self.siteRestrictions)
            msg.fatal("The number of site restrictions is not equal to the number of atomic basis vectors")

        if any( [len(i.split('/') ) != self.knary for i in self.siteRestrictions]):
            print(self.knary, [len(i.split('/') ) for i in self.siteRestrictions])
            msg.fatal("Your site restrictions are not consistent with the system you chose")

#    def _get_lattice(self,lattice):
#        """Gets the lattice vectors for the system.
#
#        Args:
#            lattice (str or list): either a string containing the lattice
#                name or a 3x3 list of the vectors as [a1,a2,a3].
#        """
#        # determine the lattice.
#        import numpy as np
#        if isinstance(lattice,string_types):
#            if lattice.lower() == "fcc":
#                self.lattice = [[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]]
#                self.lattice_name = "fcc"
#            elif lattice.lower() == "bcc":
#                self.lattice = [[0.5,0.5,-0.5],[0.5,-0.5,0.5],[-0.5,0.5,0.5]]
#                self.lattice_name = "bcc"
#            elif lattice.lower() == "sc":
#                self.lattice = [[1,0,0],[0,1,0],[0,0,1]]
#                self.lattice_name = "sc"
#            elif lattice.lower() == "hcp":
#                self.lattice = [[1,0,0],[.5, 0.866025403784439, 0],[0, 0, 1.6329931618554521]]
#                self.lattice_name = "hcp"
#            else: #pragma: no cover
#                msg.err("The lattice type {} is unsupported. Please enter your lattice vectors "
#                        "as a 3x3 matrix with the vectors as rows in the config file "
#                        "(i.e. [a1,a2,a3]).".format(lattice))
#        elif isinstance(lattice,list):
#            if len(lattice) == 3 and len(lattice[0]) == 3 and np.linalg.det(lattice) != 0:
#                self.lattice = lattice
#                self.lattice_name = "custom"
#            else: #pragma: no cover
#                msg.err("The lattice vectors must be a 3x3 matrix with the vectors as rows "
#                        "and the vectors must be linearly independent.")
#        elif lattice is not None: #pragma: no cover
#            msg.err("The lattice vectors must either be a string of 'sc', 'fcc', 'hcp', 'bcc', "
#                    "or a 3x3 matrix with the vectors a rows.")
#        else:
#            self.lattice = None
#            self.lattice_name = None
#
#
#    def _get_basis(self,basis):
#        if basis is not None and isinstance(basis,list):
#            if len(basis[0]) == 3:
#                self.basis = basis
#                self.nBasis = len(basis)
#            else: #pragma: no cover
#                msg.err("The atomic basis must be a list of lists that is nx3 where n is "
#                        "the number of atoms in the basis.")
#        elif basis is None:
#            if self.lattice is not None:
#                if self.lattice_name != "hcp":
#                    self.basis = [[0,0,0]]
#                    self.nBasis = 1
#                else:
#                    self.basis = [[0,0,0],[0.5,0.28867513459,0.81649658093]]
#                    self.nBasis = 2
#            else:
#                self.basis = [[0,0,0]]
#                self.nBasis = 1
#        else: #pragma: no cover
#            msg.err("The atomic basis must be a list of lists that is nx3 where n is "
#                    "the number of atoms in the basis or left blank.")


    def buildInputFile(self,overwrite=False):
        from os import path
        if self.inputFileExists:
            msg.info("The file you are trying to generate already exists")
            if overwrite:
                msg.info("but you said to overwrite it, so I'm moving forward")

            else:
                msg.info("and you told me not to overwrite.  Stopping, no new files generated.!")
                return
        else:
            msg.info("File: struct_enum.in({}) not found! Building it".format(self.lattice.lattice_name))


        from jinja2 import Environment, PackageLoader  # Package for building files from a template

        knaryDict = {2:"binary", 3: "ternary", 4: "quaternary", 5: "quinary"}
        settings = {}
        settings["title"] = self.lattice.lattice_name  + ' ' + knaryDict[self.knary]
        settings["template"] = "struct_enum.in"
        settings["eps"] = self.eps
        settings["size_low"] = self.sizeRange[0]
        settings["size_high"] = self.sizeRange[1]
        if self.conc_res:
            settings["conc_res"] = "T"
            if self.arrow_res:
                temp = []
                for i, a in enumerate(self.arrows):
                    temp.append("{0} {1}".format(" ".join([str(j) for j in self.concs[i]]),a))
            else:
                temp = [" ".join([str(i) for i in j]) for j in self.concRestrictions]

            settings["concentrations"] = temp
        else:
            settings["conc_res"] = "F"
        if self.arrow_res:
            settings["incl_arrows"] = "T"
        else:
            settings["incl_arrows"] = "F"
        settings["lattice"] = [" ".join([str(i) for i in j ]) for j in self.lattice.lattice]
        settings["k_nary"] = self.knary

        settings["atomic_basis"] = [" ".join([str(k) for k in self.lattice.basis[i]]) + " " + self.siteRestrictions[i] for i in range(self.lattice.nBasis)]
        settings["n_basis"] = self.lattice.nBasis


        env = Environment(loader=PackageLoader('aBuild', 'templates'))
        template = env.get_template("struct_enum.in")

        target = path.join(self.root, "struct_enum.in")
        with open(target,'w') as f:
            f.write(template.render(**settings))
        self.inputFileExists = True


    def enumerate(self,overwrite=False):
        import shutil
        from os import path
        if self.enumerationComplete:
            msg.info("It looks like you've already completed an enumeration")
            if overwrite:
                msg.info("but you said to overwrite it, so I'll continue")
            else:
                msg.info("and you said to not overwrite it. Stopping.  No new file generated")
                return
        else:
            msg.info("File: struct_enum.out not found!  Running enumeration code.")

        from os import waitpid
        from subprocess import Popen
        if config.ENUMX is not None:
            command = "cd {}; {}  {} > {}".format(self.root, config.ENUMX, 'struct_enum.in',"output." + self.lattice.lattice_name )
        else:
            msg.fatal("You haven't defined the environment variable: ENUMX, so I don't know how to enumerate")

        child=Popen(command, shell=True, executable="/bin/bash")
        waitpid(child.pid, 0)
        shutil.move(path.join(self.root,"struct_enum.out"),path.join(self.root,"struct_enum.out." + self.lattice.lattice_name))
        shutil.move(path.join(self.root,"struct_enum.in"),path.join(self.root,"struct_enum.in." + self.lattice.lattice_name))
        self.enumerationComplete = True

    @property
    def nEnumStructs(self):
        from os import path
        target = path.join(self.root,'struct_enum.out.' + self.lattice.lattice_name)
        if not path.isfile(target):
            msg.info("No struct_enum.out found.")
            return 0

        with open(target,'r') as f:
            lines = f.readlines()

        try:
            return int(lines[-1].split()[0])
        except:
            return 0

    def generatePOSCAR(self,sNumber):
        from os import waitpid
        from subprocess import Popen
        import subprocess
        if config.MAKESTRX is not None:
            command = "cd {}; {}  {}.{} {} > poscar.{}.{}; cd - ".format(self.root, config.MAKESTRX, 'struct_enum.out',self.lattice.lattice_name, sNumber, self.lattice.lattice_name,sNumber)
        else:
            msg.fatal("You haven't defined the environment variable: MAKESTRX, so I don't know how to generate POSCARs")
        child=Popen(command, shell=True, executable="/bin/bash",stdout = subprocess.PIPE)
        waitpid(child.pid, 0)


    @staticmethod
    def fromOUTPUT(root,filename):
        from os import path

        filepath = path.join(root,filename)
        if not path.isfile(filepath):
            return None
        specs = {}

        beginningLines = 30

        lines = []
        idx = 0
        with open(filepath,'r') as f:
            lines.append(f.readline())
            idx+=1

        specs["root"] = root
        specs["lattice"] = map(float,  [x.split('#')[0].split() for x in lines[3:6] ] )
        nBasis = int(lines[6])
        specs["basis"] = map(float,  [x.split('#')[0].split() for x in lines[7:7 + nBasis]] )
        specs["knary"] = lines[8 + nBasis].split('-')[0]
        specs["sizes"] = map(int,lines[8 + nBasis + 1].split()[:2])
        specs["eps"] = float(lines[8 + nBasis + 2].split()[0])
        specs["site_res"] =   [x.split('#')[1].split()[-1] for x in lines[7:7 + nBasis]]
        specs["concs"] = lines[8 + nBasis + 6].split()

        with open(filepath,'r') as f:
            f.seek(-1,2)
            nEnumerated = int(f.readline().split()[0])
        result = Enumerate(specs)
        result.setNEnumerated(nEnumerated)
        return result

    def fromINPUT(self,root,filename):
        from os import path

        filepath = path.join(root,filename)
        if path.isfile(filepath):

            specs = {}
            with open(filepath,'r') as f:
                lines = f.readlines()

            specs["root"] = root
            specs["lattice"] = map(float,  [x.split() for x in lines[2:5] ] )
            nBasis = int(lines[6])
            specs["basis"] = map(float,  [x.split()[0].split() for x in lines[8:8 + nBasis]] )
            specs["knary"] = lines[5]
            specs["sizes"] = map(int,lines[8 + nBasis].split())
            specs["eps"] = float(lines[8 + nBasis + 1])
            specs["site_res"] =   [x.split('#')[1].split()[-1] for x in lines[7:7 + self.nBasis]]
            specs["concs"] = lines[8 + nBasis + 5].split()

            result = Enumerate(specs)
            return result
        else:
            return None


    def isSame(self,enumObject):
        same = True
        lookup = {"lattice":self.lattice,"sizes":self.sizeRanges}
        if self.lattice != enumObject.lattice:
            return False
        #        if all( array(self.sizeRanges) < array(enumObject.sizeRanges) )
        for i in lookup:
            if enumObject[i] != lookup[i]:
                return False
