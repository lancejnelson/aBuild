


from aBuild.utility import chdir
from os import path
import yaml
def read(context, yfile):
    """Reads in the specified YAML file, following any additional file
    directives to compile a full representation of the template hierarchy for
    the root file.

    Args:
        context (str): path to the root folder where the yaml file is
          located. Needed for relative paths of file links.
        yfile (str): name of the template YAML file *relative* to
        `context`. Should *not* include the `.yaml` or `.yml` extension.
    """
    with chdir(context):
        if yfile[0] == ":":
            root = path.abspath(yfile[1:])
        else:
            root = path.abspath(yfile)

    if path.isfile(root + ".yml"):
        target = root + ".yml"
    else:
        emsg = ("The specified template file '{}' was not found relative "
                "to the given context directory ('{}'). Note that all files"
                " should use the `.yml` extension, *not* `.yaml`.")
        raise ValueError(emsg.format(yfile, context))
        
    with open(target, 'r') as stream:
        result = yaml.load(stream)

    #Determine the new context for recursive file links within the values of
    #this file.
    ncontext = path.dirname(target)

    #The specification allows for a "local" context that describes folder
    #locations for specific items within the template.
    lcontext = None
    if isinstance(result, dict) and "context" in result:
        lcontext = result["context"]
        del result["context"]
    
    #The unpacking command will mutate the values in result so that file links
    #are expanded to be full-fledged python objects from their YAML files.
    _unpack_obj(ncontext, result, lcontext)
    return result

def is_link(obj):
    """Determines whether the specified object is a link according to the `aBuild`
    templating specification.
    """
    result = False
    if isinstance(obj, str):
        if len(obj) > 0:
            result = obj[0] == ":"
    return result

def _unpack_obj(context, obj, lcontext=None):
    """Unpacks each item of the specified object recursively so that all
    dictionary values are visited and all list items are also visited.

    .. warning:: `obj` will be mutated if any value it considers turns out to be
      a link (according to :func:`is_link`). In that case, the file descriptor
      will be placed by the actual contents of the YAML file that the link
      points to.

    Args:
        context (str): path to the root folder where the yaml file is
          located. Needed for relative paths of file links.
        lcontext (dict): local context for the items in `obj`. Keys are the
          names of keys in `obj`; values are relative folder paths that should
          be used as the context for reads within that item.
    """
    if isinstance(obj, dict):
        result = obj
        for k, o in obj.items():
            ncontext = context
            #If the template specifies a relative context for this item,
            #then switch out the context for all of its children.
            if lcontext is not None and k in lcontext:
                with chdir(context):
                    ncontext = path.abspath(lcontext[k])
            
            if is_link(o):
                result[k] = read(ncontext, o)
            else:
                result[k] = _unpack_obj(ncontext, o)
    elif isinstance(obj, (list, set, tuple)):
        result = []
        for o in obj:
            if is_link(o):
                result.append(read(context, o))
            else:
                result.append(_unpack_obj(context, o))
    else:
        result = obj
                
    return result

