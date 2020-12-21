from os import path
from aBuild import msg

import os
import ase

def convert_direct(lattice,vec):
    from numpy import sum as nsum, array,dot
    from numpy.linalg import inv
    from numpy import transpose, array, equal
    inv_lattice = inv(lattice.transpose())
    d_space_vector = dot(inv_lattice, array(vec))

    return d_space_vector

def convert_cartesian(lattice,directVec):
    from numpy import sum as nsum, array
    return nsum([directVec[i]*lattice[i] for i in [0,1,2]], axis=0)


def vec_in_list(vec,veclist):
#    print("Considering if vector {} is in list {}".format(vec,veclist))
#    print([abs(vec - x) for x in veclist], 'diffs')
#    print([ all(abs(vec - x) < 1e-4) for x in veclist], 'check')
    if any([ all(abs(vec - x) < 1e-4) for x in veclist]):
        return True
    else:
        return False


def map_into_cell(vec):
    from math import floor
    from numpy import array
    new_point = []
    for i in vec:
        if i < 0.0 or i > 1.0:
         #   print(i,' i')
         #   print(floor(i),' floor')
         #   print(i - floor(i),' result')
            new_point.append(i - floor(i))
        elif i == 1.0:
            new_point.append(0.0)
        else:
            new_point.append(i)
    return array(new_point)

def orthogonality_defect(lattice):
    from numpy.linalg import norm,det
    from numpy import prod
    return prod([norm(x) for x in lattice])/abs(det(lattice))


def _chop(epsilon, const, i, j):
    """Sets the value of i[j] to exactly 'const' if its value already lies
    within 'epsilon' of 'const'."""
#    print('chopping')
#    print(i[j], 'to be chopped')
#    print(const, 'to go to')
#    print(epsilon, 'epsilon')
    if abs(const-i[j]) <= epsilon:
        i[j] = float(const)

def _chop_all(epsilon, i):
    """Performs chop() on the expected values of +- 0, 0.5, 1 for each
    element of i."""
    for j in range(len(i)):
        _chop(epsilon, 3, i, j)
        _chop(epsilon, -3, i, j)
        _chop(epsilon, 2, i, j)
        _chop(epsilon, -2, i, j)
        _chop(epsilon, 1, i, j)
        _chop(epsilon, -1, i, j)
        _chop(epsilon, 0, i, j)
        _chop(epsilon, 0.5, i, j)
        _chop(epsilon, -0.5, i, j)

    return i





class Lattice:


    def __init__(self,specs,symops=None, clusters=None):

        # If the input is a dictionary, I'm probably passing all of the needed
        # information (lv, bv, etc) in explicitly
        if isinstance(specs["lattice"], list):
            self._init_dict(specs)
        # If the input is a string, I'm just specifying a canonical lattice and I
        # expect the class to infer all the details
        elif isinstance(specs["lattice"], str):
            self._init_string(specs["lattice"])



    def _init_dict(self,enumdict):
        import numpy as np
        necessaryItems = ['lattice','basis','coordsys','name']

        if not all([x in enumdict for x in necessaryItems]):
            msg.fatal('Missing information when initializing Lattice object')

        if len(enumdict["lattice"]) == 3 and len(enumdict["lattice"][0]) == 3 and np.linalg.det(enumdict["lattice"]) != 0:
            self.lattice = enumdict["lattice"]
            self.lattice_name = enumdict["name"]
        else:
            msg.fatal("The lattice vectors must be a 3x3 matrix with the vectors as rows "
                        "and the vectors must be linearly independent.")
        self.basis = enumdict["basis"]
        self.coordsys = enumdict["coordsys"]
        self.lattice_name = enumdict["name"]
        self.nBasis = len(self.basis)

    def _init_string(self,string):
        lVLookupDict = {'sc':[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] ,'fcc': [[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]],'bcc': [[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]],'hcp': [[1,0,0],[.5, 0.866025403784439, 0],[0, 0, 1.6329931618554521]]}
        bVLookupDict = {'sc': [[0.0,0.0,0.0]],'fcc': [[0.0,0.0,0.0]],'bcc': [[0.0,0.0,0.0]],'hcp': [[0,0,0],[0.5,0.28867513459,0.81649658093]]}

        latDict = {}
        self.lattice = lVLookupDict[string]
        self.basis = bVLookupDict[string]
        self.coordsys = 'D'
        self.lattice_name = string
        self.nBasis = len(self.basis)



    @property
    def Bv_cartesian(self):
        from numpy import sum as nsum, array
        if self.coordsys[0].upper() != 'C':
            return [nsum([B[i]*self.Lv[i]*self.latpar for i in [0,1,2]], axis=0) for B in self.Bv]
        else:
            return self.Bv



    @property
    def Bv_direct(self):
        from numpy import sum as nsum, array
        if self.coordsys[0].upper() == 'D':
            return self.Bv

        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.Lv.transpose())
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.Bv ]

        output = []
        for i in d_space_vector:
            if i not in output:
                output.append(_chop_all(epsilon, i))

        return output


    def Lv_strings(self,nolatpar = False):
        if nolatpar:
            return [' '.join(list(map(str,x))) for x in self.Lv * self.latpar]

        else:
            return [' '.join(list(map(str,x))) for x in self.Lv]


class Crystal(object):
    """A crystal object is a lattice + a basis, including the types of atoms in the unit cell.

    Args:
        crystalSpecs (dict or str or list):  Either a dictionary containing all of the
                              necessary settings or a path to a folder that
                              contains all of the files needed.  May also be a list of strings from
                              a file read.  In this case, we will parse the lines one-by-one and read in the
                              crystallographic information for each one.  Note:  We must know what the format
                              of the file is in this case
        systemSpecies (list): List of the species for this system.  Note that this list may not
                              be identical to the species found in the POSCAR or POTCAR.  For example,
                              you may be studying a ternary system, but this particular calculation
                              has 0 of one atom type.  This argument must always be passed in because without it
                              we don't whether we have all of the atom types for the system.
    """

    def __init__(self,crystalSpecs):
        # The system species needs to always be passed in.  The crystals species is optional because we may be able to extract
        # it from input files.

        required = ["lattice","basis", "atom_types","crystalSpecies","atom_counts","latpar","coordsys","title","results","systemSpecies"]
        if True in [x not in crystalSpecs.keys() for x in required]:
            msg.fatal("Not enough information to initialize the Crystal object")
        for spec in required:
            setattr(self,spec,crystalSpecs[spec])
        self.nAtoms = sum(self.atom_counts)
        self.nTypes = len(self.atom_counts)

        if self.basis is not None and self.nAtoms != len(self.basis):
            msg.fatal("We have a problem")

        #  Let's check to see if we have a disagreement between the system species and
        # the crystal species
        self._add_zeros()

        if self.latpar is None and None not in [self.crystalSpecies,self.lattice]:
            self.set_latpar()


    #  Sometimes a crystal object will be instantiated and the number of atom types is not consistent with
    # the order of the system being studied. For example, maybe you are studying a ternary system, but some
    # of the configurations in your training set happened to be binary configurations. (Happens more often at higher order)
    #  In these cases, we need to add zeros to the atom_counts variable to clearly indicate which species are
    # present and which are not.  In other words, our standard for the Crystal object is that the number of species
    # will *always* be equal to the order of the system. No exceptions!
    def _add_zeros(self):
        if len(self.systemSpecies) == self.nTypes:
            #msg.info("No need to add zeros to the atom_types list, the number of types matches the species list")
            return
        elif any(x is None for x in [self.systemSpecies, self.crystalSpecies]):
            self.speciesMismatch = True
            msg.info("I need to know the species list for the system and this crystal to proceed adding zeros in ")
            return

        if len(self.crystalSpecies) != self.nTypes:
            print(self.crystalSpecies, self.nTypes)
            msg.fatal("The number of species that was read in (POTCAR) does not agree with atom_counts (POSCAR)")



        #  This is the case when a crystal *with the atomic species* are read in
        # but it just so happens that the number of species in this crystal does not match
        # the number of species in the system being studied.  In this case, we need to specify which
        # atomic system species are missing from this particular crystal.  We do this by augmenting zeros
        # to atom_counts at the appropriate location.
        if len(self.systemSpecies) != len(self.crystalSpecies):
            from numpy import insert
            lacking = list(set(self.systemSpecies) - set(self.crystalSpecies))
            indices = [self.systemSpecies.index(x) for x in sorted(lacking,reverse = True)]
            for idx,ele in enumerate(indices):
                self.atom_counts = insert(self.atom_counts,ele,0) # + idx ???
            self.nTypes = len(self.atom_counts)
            self.crystalSpecies = self.systemSpecies
            if len(self.crystalSpecies) != self.nTypes:
                msg.fatal('Species list ({})is not consistent with the number of atom types  ({})'.format(self.crystalSpecies,self.nTypes))
        if self.crystalSpecies != self.systemSpecies:
            msg.fatal("Function add_zeros unsuccessful")

    def randomDisplace(self):
        from numpy.random import randn
        msg.warn("Warning:  You are about to move basis atoms around")


        for ibv,bv in enumerate(self.Bv_cartesian):
            disp = 0.05 * randn(3)
            self.basis[ibv] =  bv + disp
        self.coordsys = 'C'


    @staticmethod
    def from_path(filepath,systemSpecies):
        #if self.filename is None:
        #    msg.fatal("Not sure about the format of the file you want me to read")
        #filepath = path.join(self.directory,self.filename)

        if not path.isfile(filepath):
            print('POSCAR not found', filepath)
            return None
        if 'poscar' in filepath.lower():
            return Crystal.from_poscar(filepath,systemSpecies)
        elif 'input.in' in filepath.lower():
            self.from_lammpsin(filepath)
        else:
            msg.fatal("Not sure about the format of the file you want me to read")

    @staticmethod
    def from_lines(lines,species,linesFormat):
        if linesFormat == 'mlp':
            return Crystal.fromMLP(lines,species)

#    @staticmethod
#    def from_dictionary(crystalDict):
#        required = ["lattice","basis", "atom_types","crystalSpecies","atom_counts","latpar","coordsys","title","results","systemSpecies"]
#
#        necessary = ['lattice','basis','atom_counts','coordsys','species']
#
#        if not all([x in crystalDict for x in necessary]):
#            msg.fatal("Some necessary information not set upon initialization of Crystal object")
#
#        from numpy import array
#
#        self.lattice = crystalDict["lattice"]
#        self.basis = crystalDict["basis"]
#        self.atom_counts = crystalDict["atom_counts"]
#        typesList = [[idx] * aCount for idx,aCount in enumerate(self.atom_counts)]
#        self.atom_types = []
#        for i in typesList:
#            self.atom_types += i
#        self.nAtoms = sum(self.atom_counts)
#        self.nTypes = len(self.atom_counts)
#
#        self.coordsys = crystalDict["coordsys"]
#        self.crystalSpecies = crystalDict["species"]
#        if sorted(self.species,reverse = True) != self.species:
#            msg.fatal("The order of your atomic species is not in reverse alphabetical order... OK?")
#
#        if 'title' in crystalDict:
#            self.title = crystalDict["title"]
#        else:
#            self.title = None
#
#        if 'latpar' in crystalDict:
#            self.latpar = crystalDict["latpar"]
#        else:
#            self.latpar = None



            #        self.directory = directory
#        if len(self.species) != self.nTypes:
#            msg.fatal('The number of atom types provided ({})is not consistent with the number of atom types found in the poscar ({})'.format(species,self.lattice.nTypes))
#        try:
#            self.strN = int(self.title.split()[-1])
#        except:
#            self.strN = 0
#        self.calcResults = None



    def __str__(self):
        """Returns the string representation of the POSCAR lines to write
        to a file."""
        return '\n'.join(self.lines())

    # Checks to see if any basis vectors in the crystal
    # are outside of the first unit cell.  If they are, we
    # map them back inside the first cell.
    def validateCrystal(self):
        from numpy import array,any
        from math import floor
       # print(self.Bv_direct, 'D')
        if any(array(self.Bv_direct) > 1) or any(array(self.Bv_direct) < 0):
            basis_inside = []
            for j in self.Bv_direct:
                new_point = []
                for i in j:
                    if i < 0.0 or i > 1.0:
                        new_point.append(i - floor(i))
                    elif i == 1.0:
                        new_point.append(0.0)
                    else:
                        new_point.append(i)
                basis_inside.append(new_point)
            self.basis = basis_inside
            self.coordsys = 'D'
           # msg.info('Crystal is fixed')
        else:
            pass
           # msg.info("Crystal didn't need fixing")

    @property
    def orthogonality_defect(self):
        from numpy.linalg import norm,det
        from numpy import prod

        return prod([norm(x) for x in self.lattice])/abs(det(self.lattice))
    @property
    def cell_volume(self):
        from numpy import dot,cross
        return abs(dot(cross(self.lattice[0],self.lattice[1]),self.lattice[2]))

    # Calculates the distances between all of the atoms and finds
    # the minimum value from all of them.
    @property
    def minDist(self):
        from numpy import array,dot,min,einsum,add,roll,column_stack
        from numpy.linalg import norm
        from itertools import product
        import sys
        import numpy
        numpy.set_printoptions(threshold=sys.maxsize)


        # Need to make sure that all of the atoms are inside the first
        # unit cell before we compile list of distances
        self.validateCrystal()

        #  Calculate all possible shifts.  These are vectors of integers
        # representing the amount of each lattice vector that we are going
        # to add to each basis atom.  We only do combinations of (-1,0,1) because
        # that should be enough to find all possible distances between atoms.
        #  We're not trying to get all atoms out to some cutoff radius, we just want to
        # make sure we get enough atoms in there to find the min separation.
        offsets = array([x for x in product(range(-1,2),repeat = 3)])

        # Now shift every basis atom by every shift previously calculated
        neighborsDirect  = array([self.Bv_direct + array(x) for x in offsets])
        #Convert list of atomic positions to cartesian coordinates
        neighborsCartesian = einsum('abc,cd',neighborsDirect,self.latpar * self.lattice)
        # Flatten the list down to a single list of position vectors
        neighborsCartesian.resize(len(offsets)* self.nAtoms,3)

        # Build a matrix where each row is a shifted version of atomic positions.
        rolledNeighbors = array([roll(neighborsCartesian,x,axis = 0) for x in range(len(neighborsCartesian))])
#        print(neighborsCartesian, 'nc')
#        print(rolledNeighbors,' rn')
#        print(rolledNeighbors.shape,' rn shape')
#        print(neighborsCartesian.shape,' nc shape')
#        import sys
#        sys.exit()
        #Now we can just subtract the first row (unshifted) from all of the other rows
        # and calculate the norm of each vector
        distances = norm(rolledNeighbors[0,:,:] - rolledNeighbors,axis = 2)
        # Return the min, excluding 0 distances.
        from numpy import count_nonzero,nonzero,transpose
        nZeroOccurrences = count_nonzero(distances==0)
        if nZeroOccurrences > self.nAtoms * len(offsets):
            ##import sys
            #import numpy
            #numpy.set_printoptions(threshold=sys.maxsize)
            #print(distances, 'distances')
            #print(transpose(nonzero(distances == 0)))
            msg.fatal("Atoms are on top of each other")
        return min(distances[distances > 1e-5])

#    @property
#    def appMinDist(self):
#        from aBuild.calculators import data
#        return data.nnDistance(self.species,self.atom_counts)

    @property
    def recip_Lv(self):
        if len(self.lattice) == 0:
            raise ValueError("Lattice vectors are required for finding reciprocal"
                             " lattice vectors.")

        from numpy import array,cross,dot
        groupings = [[self.lattice[1], self.lattice[2]], [self.lattice[2], self.lattice[0]],
                     [self.lattice[0], self.lattice[1]]]
        crosses = [cross(v[0], v[1]) for v in groupings]
        dots = [dot(self.lattice[i], crosses[i]) for i in [0,1,2]]
        return [crosses[i]/dots[i] for i in [0,1,2]]


    @property
    def volume(self):

        from numpy import cross, dot
        return dot(cross(self.lattice[0],self.lattice[1]), self.lattice[2])

    @property
    def concentrations(self):
        return [self.atom_counts[x]/sum(self.atom_counts) for x in range(self.nTypes)]

    @property
    def Bv_cartesian(self):
        from numpy import sum as nsum, array
        if self.coordsys[0].upper() != 'C':
            return [nsum([B[i]*self.lattice[i]*self.latpar for i in [0,1,2]], axis=0) for B in self.basis]
        else:
            return self.basis


    def findNeighbors(self,rcut):
        from numpy import array,dot,min,einsum,add,roll,column_stack,append,where,extract,argwhere,set_printoptions,round
        from numpy.linalg import norm
        from itertools import product

        # Get all of the offsets to translate the basis vectors to equivalent
        # positions
        offsets = array([x for x in product(range(-1,2),repeat = 3)])

        # Translate them (in direct coordinates)  Shape: nBasis x nOffset x nCoord (3)
        neighborsDirect  = array([x + offsets for x in self.Bv_direct])
        # Convert coordinates to cartesian
        neighborsCartesian = einsum('abc,cd',neighborsDirect,self.latpar * self.lattice)
        diffs = array([x  - neighborsCartesian for x in self.Bv_cartesian])
        distances = norm(diffs, axis = 3)
        keepIndices = argwhere(distances < rcut)
        self.neighbors = [[] for x in self.atom_types]
        for [centerAtom,neighborAtom,shift] in keepIndices:
            self.neighbors[centerAtom].append([round(neighborsCartesian[neighborAtom,shift],3), self.atom_types[neighborAtom],neighborAtom])



    @property
    def Bv_direct(self):
        from numpy import sum as nsum, array,dot
        if self.coordsys[0].upper() == 'D':
            return self.basis
        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.lattice.transpose()*self.latpar)
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.basis ]

        epsilon = 1e-4
        output = []
        for i in d_space_vector:
            if _chop_all(epsilon, i) not in output:
                output.append(_chop_all(epsilon, i))

        return output


    @property
    def lattice_lines_LAMMPS(self):
        """Return \n joined lattice vector text lines."""
        lines = ' a1 '
        lines += ' '.join(map(str,self.lattice[0]))
        lines += ' a2 '
        lines += ' '.join(map(str,self.lattice[1]))
        lines += ' a3 '
        lines += ' '.join(map(str,self.lattice[2]))

        return lines
#        return zip(['a1','a2','a3'],[' '.join(map(str,x)) for x in self.lattice]) #'\n'.join(list(self.lattice))

    @property
    def lattice_lines(self):
        """Return \n joined lattice vector text lines."""

        return '\n'.join( [' '.join(map(str,x)) for x in self.lattice]) #'\n'.join(list(self.lattice))

    @property
    def lattice_lines_nolatpar(self):
        """Return \n joined lattice vector text lines."""

        return '\n'.join( [' '.join(map(str,x)) for x in self.lattice * self.latpar]) #'\n'.join(list(self.lattice))

    @property
    def basis_lines_LAMMPS(self):
        """Return \n joined lattice vector text lines."""
        lines = ''
        for i in self.basis:
            lines += ' basis '
            lines += ' '.join(map(str,i))

        return lines

    @property
    def basis_lines_ESPRESSO(self):
        """Return \n joined lattice vector text lines."""
        speciesList = []
        for i in range(self.nTypes):
            speciesList += [self.species[i] for x in range(self.atom_counts[i])]

        lines = ''
        for idx,i in enumerate(self.basis):
            lines += speciesList[idx] + ' '
            lines += ' '.join(map(str,i))
            lines += '\n'

        return lines[:-1]

    @property
    def unknown(self):
        """Return \n joined lattice vector text lines."""
        atomLabels = [ x for sublist in [ [  i for k in range(self.atom_counts[i]) ]   for i in range(self.nTypes)] for x in sublist]
        lines = ""
        for i,v in enumerate(atomLabels):
            lines += ' basis '
            lines += str(i+1) + ' ' + str(v+1)

        return lines

    @property
    def basis_lines(self):
        """Return \n joined basis vector text lines."""
        species = []
        for i,n in enumerate(self.atom_counts):
            species += [self.crystalSpecies[i] ] * n
        return '\n'.join( [' '.join(list(map(str,y)) + [species[x]] ) for x,y in enumerate(self.basis)])

    def mtpLines(self,relax = False):
        import numpy as np
        if not relax and self.results is None:
            msg.fatal("You want me to write result information but I don't have any.")
        result = []
        result.append('BEGIN_CFG')
        result.append('Size')
        result.append(str(self.nAtoms))
        result.append('SuperCell')
        for lv in self.latpar * self.lattice:
            result.append('{:12.6f} {:12.6f} {:12.6f}'.format(lv[0],lv[1],lv[2]  ))
        if not relax:
            result.append('   AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz')
        else:
            result.append('   AtomData:  id type       cartes_x      cartes_y      cartes_z')
        #counter = crystal.lattice.nTypes - 1

        # Took me a few minutes to figure this one out.  Very Pythonic
        #        atomLabels = [ x for sublist in [ [ counter - i for k in range(crystal.lattice.atom_counts[i]) ]   for i in range(counter + 1)] for x in sublist]
        atomLabels = [ x for sublist in [ [  i for k in range(self.atom_counts[i]) ]   for i in range(self.nTypes)] for x in sublist]
        for i in range(self.nAtoms):
            if not relax:
                forces = self.results["forces"][i]
            coords = self.Bv_cartesian[i]
#            print(forces,'forces')
            if not relax:
                result.append('{:16d} {:3d} {:16.6f} {:12.6f} {:12.6f} {:18.6f} {:10.6f} {:10.6f}'.format(i+ 1,atomLabels[i], coords[0],coords[1],coords[2], forces[0],forces[1],forces[2]  ))
            else:
                result.append('{:16d} {:3d} {:16.6f} {:12.6f} {:12.6f}'.format(i+ 1,atomLabels[i], coords[0],coords[1],coords[2] ))
            #line +=  ' '.join([map(str,crystal.lattice.Bv_cartesian[i]), crystal.forces[i]])

        if not relax:
            result.append('Energy')
            result.append(str(self.results["energyZ"]) + '')


            result.append(' Stress:   xx          yy           zz            yz           xz           xy')
            s = self.results["stress"]
            stressesline = '{:16.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}'.format(s[0],s[1],s[2],s[3],s[4],s[5])
            result.append(stressesline)
        result.append(''.join([' Feature   conf_id ', '  '.join([self.symbol,self.title]),'']))
        result.append('END_CFG\n')
        return result

    def vasplines(self,keepZeros = False):
        """Returns a list of strings for each line in the POSCAR file.

        :arg vasp: when true, the atom_counts line is checked for zeros before
          it is created. Vasp can't handle zero for the number of atoms of a
          certain type; just remove it."""
        from numpy import where
        result = []
        result.append(self.title)
        result.append(str(self.latpar))
        result.append(self.lattice_lines)

        idxKeep = list(where( self.atom_counts > 0)[0])
        if keepZeros:
            result.append(' '.join(map(str,self.atom_counts)))
        else:
            result.append(' '.join(map(str,self.atom_counts[idxKeep])))

        result.append(self.coordsys)
        result.append(self.basis_lines)
        return result

    def lines(self,fileformat,keepZeros = False):
        if fileformat.lower() == 'vasp':
            return self.vasplines(keepZeros = keepZeros)
        elif fileformat.lower() == 'mtptrain':
            return self.mtpLines()
        elif fileformat.lower() == 'mtprelax':
            return self.mtpLines(relax = True)

    def write(self, filename, fileformat = 'vasp',keepZeros = False):
        """Writes the contents of this POSCAR to the specified file."""
        #fullpath = os.path.abspath(filepath)
#        if not self.speciesMismatch:
        with open(filename, 'w') as f:
            f.write('\n'.join(self.lines(fileformat,keepZeros=keepZeros)))

    def concsOK(self,concRestrictions=None):
        from numpy import all,array
        if concRestrictions == None:
            return True
        lowerlimit = array([x[0]/x[2] for x in concRestrictions])
        upperlimit = array([x[1]/x[2] for x in concRestrictions])
        if all(self.concentrations >= lowerlimit) and all(self.concentrations <= upperlimit):
            return True
        return False
        #  This routine switches A atoms for B atoms and vice versa
        #  For prototype structures, we only provide one crystal structure
        # excluding the concentration brothers.
    def scrambleAtoms(self,scrambleKey):
        from numpy import array
        Bvs = []
        start = 0
        end = self.atom_counts[0]
        Bvs.append(self.basis[start:end])
        for index,aC in enumerate(self.atom_counts[:-1]):
            start = start + aC
            end = end + self.atom_counts[index + 1]
            Bvs.append(self.basis[start:end])

        self.basis = [y for sublist in [Bvs[x] for x in scrambleKey] for y in sublist]
        self.atom_counts = [self.atom_counts[x] for x in scrambleKey]
        typesList = [[idx] * aCount for idx,aCount in enumerate(self.atom_counts)]
        self.atom_types = []
        for i in typesList:
            self.atom_types += i
        self.set_latpar()

    @property
    def symbol(self):
        symbol = ''
        for elem, count in zip(self.systemSpecies,self.atom_counts):
                symbol += elem + '_' + str(count)
        return symbol

    def set_latpar(self,modify = 1.0):
        from aBuild.calculators import data
        #if self.latpar == 1.0 or self.latpar < 0:
            # We must first reverse sorte the species list so we get the right atom in the right place.
        if sorted(self.crystalSpecies,reverse = True) != self.crystalSpecies:
            msg.fatal("Your species are not in reverse alphabetical order... OK?")
        self.latpar = data.vegardsVolume(self.crystalSpecies,self.atom_counts,self.volume)
        self.latpar = self.latpar * modify
    @staticmethod
    def from_poscar(filepath,systemSpecies):
        """Returns an initialized Lattice object using the contents of the
        POSCAR file at the specified filepath.

        :arg strN: an optional structure number. If the label in the POSCAR doesn't
          already include the strN, it will be added to the title.
        """
        crystalDict = {}
        from aBuild.calculators.vasp import POSCAR
        lines = POSCAR.from_path(filepath)
        from numpy import array
        crystalDict["lattice"] =  array([list(map(float, l.strip().split()[0:3])) for l in lines.Lv])
        crystalDict["basis"] = array([list(map(float, b.strip().split()[:3])) for b in lines.Bv])
        if lines.species is not None:
            # Found species inside POSCAR file
            crystalDict["crystalSpecies"] = lines.species
        else:
            # Assume that the passed in species apply to this crystal
            crystalDict["crystalSpecies"] = systemSpecies
        crystalDict["atom_counts"] = array(list(map(int, lines.atom_counts.split() )  ))
        if len(crystalDict["crystalSpecies"]) != len(crystalDict["atom_counts"]):
            msg.fatal("Number of species does not agree with the number of atom types read in")
        typesList = [[idx] * aCount for idx,aCount in enumerate(crystalDict["atom_counts"])]
        atom_types = []
        for i in typesList:
            atom_types += i
        crystalDict["atom_types"] = atom_types
        crystalDict["latpar"] = float(lines.latpar.split()[0])
        crystalDict["coordsys"] = lines.coordsys
        crystalDict["title"] =  lines.label # path.split(filepath)[0].split('/')[-1] + '_' +
        crystalDict["systemSpecies"] = systemSpecies
        crystalDict["results"] = {}

        return Crystal(crystalDict)
    @property
    def reportline(self):
        line = [self.title , self.results["energyF"],self.results["energyF"]/self.nAtoms, self.results["fEnth"] if self.results["fEnth"] is not None else 1000]
        lineWrite = '{:55s} {:9.4f} {:12.4f}{:12.4f}'
        if self.results["distToHull"] is not None:
            line.append(self.results["distToHull"])
            lineWrite += '   {:13.6f}       '
        else:
            line.append('-----')
            lineWrite += '   {:8s}   '

        for i in self.concentrations:
            line.append(i)
            lineWrite += '  {:4.2f}  '
        for i in self.atom_counts:
            line.append(i)
            lineWrite += '  {:4d}  '
        lineWrite += '\n'
            #        for i in range(self.knary):
#            line.append(self.pures[i].results["energy"])
#            lineWrite += '{:8.5f}'
#        for i in range(self.knary):
#            line.append(self.atom_counts[i])
#            lineWrite += '{:2d}'

        print(line)
        return lineWrite.format(*line)


    def from_lammpsin(self,filepath):
        from numpy import array
        with open(filepath,'r') as f:
            lines = f.readlines()

        for idx,line in enumerate(lines):
            if 'lattice' in line:
                self.lattice = array([map(float,line.split()[4:7]),map(float,line.split()[8:11]), map(float,line.strip('&').split()[12:15])])
                nBasis = lines[idx + 1].count('basis')
                self.basis = [map(float,lines[idx+1].strip().strip('basis').split()[i * 4:4 * (i+ 1) - 1]) for i in range(nBasis)]
                self.latpar = float(line.split()[2])
            if 'create_box' in line:
                self.nTypes = int(line.split()[1])

            if 'create_atoms' in line:
                self.atom_counts = [0 for i in range(self.nTypes) ]
                for idx,elem in enumerate(line.split()[5::3]):
                    self.atom_counts[int(elem)-1] += 1
        self.coordsys = 'D'
        self.title = path.split(filepath)[0].split('/')[-1] + lines[1].split('\n')[0]

    @staticmethod
    def fromMLP(lines,species):
        from numpy import array
        import os
        from aBuild.calculators.vasp import VASP
        nAtoms = int(lines[2].split()[0])

        results = {}
        crystalDict = {}
        crystalDict["systemSpecies"] = species
        crystalDict["lattice"] = array([list(map(float,x.split())) for x in lines[4:7]])
        crystalDict["basis"] = array([list(map(float,x.split()[2:5])) for x in lines[8:8 + nAtoms]])
        #LJNself.nAtoms = len(self.basis)
        crystalDict["coordsys"] = 'C'
        crystalDict["atom_types"] = [int(x.split()[1]) for x in lines[8:8 + nAtoms]]
        crystalDict["crystalSpecies"] = [species[x] for x in list(set(crystalDict["atom_types"]))]
        if len(lines[8]) > 5:
            results["forces"] = array([list(map(float,x.split()[5:8])) for x in lines[8:8 + nAtoms]])
        else:
            results["forces"] = None
        crystalDict["atom_counts"] = array([ crystalDict["atom_types"].count(x) for x in set(crystalDict["atom_types"])]) #range(max(atoms)+1)
        #LJNself.nTypes = len(self.atom_counts)
        # THe add_zeros function is here in case you run into a config with few atom types than the
        # system being studied. Like running into a pure config for a binary system.  Or a binary config
        # when studying a ternary system. Commented out because the add_zeros function gets called in
        # the initializer anyway.  No need to call it twice.
        #self.species = [self.systemSpecies[x] for x in list(set(self.atom_types))]
        #self._add_zeros(self.systemSpecies,self.species)
        titleindex = ['conf_id' in x for x in lines].index(True)
        crystalDict["title"] = ' '.join(lines[titleindex].split()[2:])
        if any(['Energy' in x for x in lines] ):
            energyindex = ['Energy' in x for x in lines].index(True)

            results["energyF"] = float(lines[energyindex + 1].strip().split()[0])
#            root = os.getcwd()
#            pures = [VASP(path.join(root,'training_set','pure' + x),systemSpecies = self.species)   for x in self.species]
#            puresDict = {}
#            for ispec,spec in enumerate(self.species):
#                #thispure = VASP(path.join(root,'training_set','pure' + spec,systemSpecies = self.species))
#                pures[ispec].read_results()
##                thispure.read_results()
#                puresDict[spec] = pures[ispec].crystal.results["energypatom"]
#            self.results["fEnth"] = self.results["energyF"]/self.nAtoms - sum(   [ pures[i].crystal.results["energyF"]/pures[i].crystal.nAtoms * self.concentrations[i] for i in range(self.nTypes)])
        else:
            results["energyF"] = None
#        if sum(self.atom_counts) != nAtoms:
#            msg.fatal('atomCounts didn\'t match up with total number of atoms')
#        self.set_latpar()
        crystalDict["latpar"] = 1.0  # MLP files are formatted with no lattice parameter.  It's
                           # already built into the lattice vectors.
#        self.lattice = self.lattice / self.latpar  # I think I did this to ensure that the lattice
                                                    # vectors didn't change but I know that the lattice
                                                    # parameter is just 1.0 for MLP formatting.

        crystalDict["results"] = results
        return Crystal(crystalDict)

    @staticmethod  # Needs fixed!!!
    def fromEnum(enum,structNum,species):
        from aBuild.enumeration import Enumerate

#        enumLattice = Enumerate(enumDict)
        enum.generatePOSCAR(structNum)
        result = Crystal.from_poscar(path.join(enumLattice.root,"poscar.{}.{}".format(enumDict["name"],structNum)),species)
        result.set_latpar()
        return result


#    def setCalcResults(self):
#
#        from aBuild.calculators.vasp import VASP
#        from aBuild.utility import chdir
#
#        thiscalc = VASP(directory = self.directory)
#        #        with chdir(self.directory):
#        thiscalc.read_results(allIonic = False,allElectronic= False)
#        self.calcResults = thiscalc.results
#        if self.calcResults is not None and len(self.calcResults["forces"]) != self.lattice.nAtoms:
#            msg.fatal('The number of forces ({}) is not equal to the number of atoms ({})'.format(self.calcResults["forces"],self.lattice.nAtoms) )


    def mlpLines(self):
        pass

    def setAtypes(self):
        pass

    #    def writePOSCAR(self,filepath):
    #    from aBuild.calculators.vasp import POSCAR

        # Instantiate a POSCAR object from my crystal object so that it intializes to have
        # everything that's inside of crystal
        #    lines = POSCAR(self)
        #lines.write(filepath,vasp=True)


    def superPeriodics(self,size):
        from numpy import dot,array,prod,einsum,any,all,copy
        from numpy.linalg import inv,det,norm
        from numpy import sum as nsum,transpose,matmul,cross
        from itertools import product,combinations

        crystals = []
        # If there are 0 occurrences of the first atom type. Just skip it.
        if self.atom_counts[0] == 0:
            return []
        # Generate all super-periodic lattice vectors
        dimension = len(self.lattice)
        multiplesOf = array([x for x in combinations([x for x in product(range(-2,3),repeat = dimension) if sum(abs(array(x))) !=0],dimension) ])
        lattices = einsum('abc,dc->abd',multiplesOf,transpose(self.lattice))
        #Keep only the lattices with i) positive triple product  ii) specified size increase or less  iii) orthogonality defect
        # less than 2, and sort them by their orthogonality defect.
        keeps = sorted([x  for x in lattices if det(x) > 1e-3 and abs(det(x))/abs(det(self.lattice)) <= size and orthogonality_defect(x)/orthogonality_defect(self.lattice) < 2], key = lambda k: orthogonality_defect(k))

#        sorted(keeps, key = lambda k: orthogonality_defect(k))
        # Done getting super-periodic lattice vectors.

        # Get all atoms by adding multiples of the parent lattice vectors
        combs =  [x for x in product(range(-4,4),repeat = 3)]
        latticeVecCombinations = einsum('ab,bd->ad', combs,self.lattice*self.latpar)
        basisAtoms = [x + latticeVecCombinations for x in self.Bv_cartesian]
        count = 1
        for lattice in keeps:

            count += 1
            basDirect = [map_into_cell(_chop_all(1e-4,map_into_cell(convert_direct(lattice*self.latpar, x)))) for x in self.Bv_cartesian]
            crystalDict = {"lattice":lattice, "basis":basDirect, "coordsys":'D', "atom_counts":[x for x in self.atom_counts],"species":self.systemSpecies,"latpar":self.latpar,"title":self.title + '_super'}
            newCrystal = Crystal(crystalDict,self.sytemSpecies)
            if newCrystal.orthogonality_defect/self.orthogonality_defect > 2:
                print('Cell is too skew, not considering it')
                continue
            newCrystal.validateCrystal() #Maps all the atoms into the first unit cell
            sizeIncrease = round(newCrystal.cell_volume/self.cell_volume)
            if newCrystal.cell_volume < 1e-3:
                print("Found a zero volume cell!! Continuing without considering it")
                continue
            if abs(int(sizeIncrease) - sizeIncrease) > 1e-3:
                print(sizeIncrease,int(sizeIncrease), "Cell volume increase doesn't appear to be an integer multiple of the parent volume")
                import sys
                sys.exit()

            # We don't need to populate the cell if it's size 1 because we know there
            # are the same number of basis atoms as the original
#            if sizeIncrease == 1:
#                continue
            # Now populate the cell with all of the basis atoms.
            print("condisdering lattice {}".format(lattice))
            for idx, origbasis in enumerate(basisAtoms):
                #print(idx,' adding atoms of this type')
                track = 1
                for equiv in origbasis:
                    # First check to see if we have enough atoms of this type already.
                    if track == sizeIncrease:
                        break

                    candidate = map_into_cell(_chop_all(1e-4,map_into_cell(convert_direct(lattice*self.latpar, equiv))))
#                    candidate = _chop_all(1e-4,convert_direct(lattice*self.latpar, equiv))
#                    candidate = _chop_all(1e-4,map_into_cell(_chop_all(1e-4,convert_direct(lattice*self.latpar,equiv))))
                    if not vec_in_list(candidate,newCrystal.basis) and (not any(candidate >= 1) and not any(candidate < 0) ) :
                        newCrystal.basis.append(candidate)
                        atomType = self.atom_types[idx]
                        newCrystal.atom_types.append(atomType)
                        newCrystal.atom_counts[atomType] += 1
                        newCrystal.nAtoms += 1
                        track += 1

            if abs(newCrystal.nAtoms - sizeIncrease * self.nAtoms) > 1e-5:
                print("You don't have enough atoms")
                print(self.nAtoms, 'n primitive atoms')
                print(newCrystal.nAtoms, 'n current atoms')
                print(newCrystal.lattice, 'new lattice')
                print(self.cell_volume, 'size primitive')
                print(newCrystal.cell_volume, 'size current')
                print(det(transpose(lattice)), 'det')
                print(sizeIncrease, 'size increase')
                print(array(newCrystal.basis), 'new basis')
                print(newCrystal.atom_counts,'atom counts')
                print(newCrystal.atom_types,'atom types')
                print(len(crystals))
                import sys
                sys.exit()
            sortKey = sorted(range(newCrystal.nAtoms), key = lambda x: newCrystal.atom_types[x])
            newCrystal.basis = [newCrystal.basis[x] for x in sortKey]
            newCrystal.atom_types = sorted(newCrystal.atom_types)

            print(newCrystal.minDist, 'mindist)')
            with open('supers.out','a+') as f:
                f.writelines("Title\n")
                f.writelines("{}\n".format(sizeIncrease))
                f.writelines(newCrystal.lattice_lines)
                f.writelines("\n{} {}\n".format(newCrystal.atom_counts[0],newCrystal.atom_counts[1]))
                f.writelines(newCrystal.basis_lines)
                f.writelines('\n\n\n')

            if abs(newCrystal.minDist - self.minDist)    > 1e-3:
                msg.fatal("mindist for super-periodic is different than parent")

            print("Got Super-periodic.  Checking for AFM compatible")
            if newCrystal.getAFMPlanes([1,0,0]):
                return newCrystal
            print('not AFM compatible')
            #crystals.append( newCrystal )
        #import sys
        #sys.exit()
        return []
