
from contextlib import contextmanager


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

def grep(file,tag):

    from os import waitpid
    from subprocess import Popen,PIPE

    with open(file,'r') as f:
        lines = f.readlines()
    matchedlines = []
    for line in lines:
        if tag in line:
            matchedlines.append(line)
    return matchedlines
#    command = "grep {} {}".format(tag,file)
#    child=Popen(command, shell=True, executable="/bin/bash",stdout=PIPE)
#    result = child.communicate()[0]
#    #    waitpid(child.pid, 0)
#    print(child.stdout.decode('utf-8'))
#    return child.stdout

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
    
