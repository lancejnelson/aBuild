
from contextlib import contextmanager
from subprocess import Popen,PIPE,run,check_output,STDOUT
from os import waitpid,path
import mmap

@contextmanager
def chdir(target):
    """Context manager for executing some code within a different
    directory after which the current working directory will be set
    back to what it was before.

    Args:
        target (str): path to the directory to change into.
    """
    from os import getcwd, chdir
    current = getcwd()
    try:
        chdir(target)
        yield target
    finally:
        chdir(current)

def head(filename,lineNum,postProcess):
    command = "head -{} {} {}".format(lineNum,filename,postProcess)
    child=Popen(command, shell=True, executable="/bin/bash",stdout=PIPE)
    result = child.communicate()[0]
    #    waitpid(child.pid, 0)
    print(child.stdout.decode('utf-8'))
    return child.stdout

def tail(filename,lineNum,postProcess):
    command = "tail -{} {} {}".format(lineNum,filename,postProcess)
    print('running tail command: {}'.format(command))
    child=Popen(command, shell=True, executable="/bin/bash",stdout=PIPE)
    result = child.communicate()[0]
    #    waitpid(child.pid, 0)
#    print(child.stdout.decode('utf-8'))
    return child.stdout

#LJNdef cat(filename,postProcess):
#LJN    command = "cat {} {}".format(filename,postProcess)
#LJN    child=Popen(command, shell=True, executable="/bin/bash",stdout=PIPE)
#LJN    result = child.communicate()[0]
#LJN    #    waitpid(child.pid, 0)
#LJN    print(result.stdout.decode('utf-8').split('\n')[:-1])
#LJN    return result.stdout.decode('utf-8').split('\n')[:-1]


def sed(options,string,fileName,postProcess):
    command = "sed {} '{}' {} {}".format(options,string,fileName,postProcess)
    child=Popen(command, shell=True, executable="/bin/bash",stdout=PIPE)
    result = child.communicate()[0]
    #    waitpid(child.pid, 0)
    return result.decode('utf-8').split('\n')[:-1]


def rgrep(filename,tag):
    if not path.isfile(filename):
        return None
    lines = []
    with open(filename,'r') as f:
        m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
        i = m.rfind(str.encode(tag))
        print(i, 'i')
        if i > 0:
            m.seek(i)
            return m.readline()
        
def grep(filename,tag,options='',postprocess = ''):

#    from os import waitpid
#    from subprocess import Popen,PIPE
#
#    if 'xz' in filename:
#        import lzma
#        with lzma.open(filename,'rt') as f:
#            lines = f.readlines()
#    else:
#        with open(filename,'r') as f:
#            lines = f.readlines()
#    matchedlines = []
#    for line in lines:
#        if tag in line:
#            matchedlines.append(line)
#    return matchedlines

#LjN    print("Grepping for {} in file {}".format(tag,filename))
#LjN    if not path.isfile(filename):
#LjN        return None
#LjN    lines = []
#LjN    with open(filename,'r') as f:
#LjN        m = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
#LjN        i = m.find(str.encode(tag))
#LjN        print(i,'i')
#LjN        while i > 0:
#LjN            m.seek(i)
#LjN            lines.append(m.readline())
#LjN            i = m.find(tag)
#LjN    return lines
    if 'xz' in filename:
        command = 'xzgrep {} "{}" {} {}'.format(options, tag,filename,postprocess)
    else:
        command = 'grep {} "{}" {} {}'.format(options, tag,filename,postprocess)
    print('greping file {} for tag {} with command: {}'.format(filename,tag,command))
    child=Popen(command, shell=True, executable="/bin/bash",stdout=PIPE)
    result = child.communicate()[0]
#    print('Result',result.decode('utf-8').split('\n')[:-1])
    return result.decode('utf-8').split('\n')[:-1]#child.stdout
    #if len(result.decode('utf-8').split('\n')) == 1 and result.decode('utf-8').split('\n')[0] == '':
    #    return []
    #elif 
    #else:
    #    waitpid(child.pid, 0)
    #    return result.decode('utf-8').split('\n')[:-1]#child.stdout

def _get_reporoot():
    """Returns the absolute path to the repo root directory on the current
    system.
    """
    from os import path
    import aBuild
    #    global reporoot
    #if reporoot is None:
    medpath = path.abspath(aBuild.__file__)
    reporoot = path.dirname(path.dirname(medpath))

    return reporoot


def unpackProtos():
    from os import path
    template_root = path.join(_get_reporoot(), "aBuild", "templates")
    if not path.isdir(path.join(template_root, "uniqueUnaries")):
        with chdir(template_root):
            tarf = "prototypes.tar.gz"
            tar = tarfile.open(tarf, "r:gz")
            tar.extractall()
            tar.close()


def cat(files, target,remove = False):
    """Combines the specified list of files into a single file.

    Args:
        files (list): of `str` file paths to combine.
        target (str): name/path of the output file that will include all of the
          combined files.
    """
    from os import remove
    if files == []:
        return
    
    with open(target, 'w') as outfile:
        for fname in files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    if remove:
        for file in files:
            remove(file)
    

def getAllPerms(knary,justCyclic=False):
    if justCyclic:
        start = list(range(knary))
        perms = [start]
        for i in range(1,knary):
            perms.append(start[i:] + start[:i])
    else:
        from itertools import permutations
        perms = list(permutations(range(knary)))

    return perms

    
def getProtoPaths(knary):
    from os import path
    from glob import glob

    templatePath = path.join(_get_reporoot(),'aBuild','templates')
    unpackProtos()
    unaries = glob('{}/*'.format(path.join(templatePath,'uniqueUnaries')))
    binaries = glob('{}/*'.format(path.join(templatePath,'uniqueBinaries')))
    ternaries = glob('{}/*'.format(path.join(templatePath,'uniqueTernaries')))
    if knary == 1:
        return unaries
    elif knary == 2:
        return unaries + binaries
    else:
        return unaries + binaries + ternaries
    
def convert_direct(lattice,vec):
    from numpy import sum as nsum, array,dot
    from numpy.linalg import inv
    from numpy import transpose, array, equal
    inv_lattice = inv(lattice.transpose())
    d_space_vector = dot(inv_lattice, array(vec))

    return d_space_vector


def vec_in_list(vec,veclist):
#    print("Considering if vector {} is in list {}".format(vec,list))
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


def _chop(epsilon, const, i, j):
    """Sets the value of i[j] to exactly 'const' if its value already lies
    within 'epsilon' of 'const'."""
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


    
#def fileinDir(searchfile,directory, or_close=False,returnFiles = True):
#    from os import listdir
#    if not or_close:
#        return True in [searchfile == x for x  in listdir(directory)]
#    else:
#        if returnFiles:
#            return [x for x in listdir(directory) if searchfile in x]
#        else:
#            return True in [searchfile in x for x  in listdir(directory)]
