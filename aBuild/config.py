#This config module handles retrieval of global variable values for pyc from environment
#variables or a config XML file.
import types
from sys import modules
from os import path as ppath

class _config(types.ModuleType):
    def __init__(self):
        """A configuration class to store global variables under a configuration
        module name. Exposes values as properties to allow arbitrary variable
        storage via XML as well as static oft-used variables."""
        self._vardict = {}
        self._initialized = False
        self.getenvar("MAKESTRX")
        self.getenvar("GETKPTS")
        self.getenvar("POTCAR_DIR")
        self.getenvar("INPUT_DIR")
        self.getenvar("RUN_DIR")
        self.getenvar("STORE_DIR")
        self.getenvar("TEMPL_DIR")
        self.getenvar("UNCLEX")
        self.getenvar("ENUMX")
        self.getenvar("PWX")
        self._initialized = True             

#    def UNCLE_static(self, filename):
#        """Returns the full path to the specified static library file.
#
#        :arg filename: the name of the *.so static library file.
#        """
#        uncle = self.property_get("UNCLE_LIB")
#        if uncle is not None:
#            jpath = ppath.join(uncle, filename)
#            return ppath.expanduser(jpath)
#        else:
#            return filename
#
#    def static_symbol(self, module, method, lib=None):
#        """Returns the symbol for the specified *fortran* module and subroutine/function
#        that has been compiled into a shared library with the compiler listed in the
#        'COMPILER' config variable.
#
#        :arg lib: the loaded ctypes library python object obtained by calling
#          ctypes.cdll.LoadLibrary('shared.so')
#        """
#        compiler = self.property_get("COMPILER", "gfortran")
#        identifier = "{}.{}".format(module, method)
#        if compiler == "gfortran":
#            symbol = "__" + identifier.lower().replace(".", "_MOD_")
#        else:
#            #We assume the only other option is ifort.
#            symbol = identifier.lower().replace(".", "_mp_") + "_"
#
#        if lib is None:
#            return symbol
#        elif hasattr(lib, symbol):
#            return getattr(lib, symbol)
#        else:
#            return None
#
#    @property
#    def implicit_XML(self):
#        """Returns the full path to an XML file that contains path information for ANCLE."""
#        return self.property_get("ANCLE_XML")
#
#    @property
#    def gateway(self):
#        """Returns the SMTP gateway for sending email reports."""
#        return self.property_get("GATEWAY")

    @property
    def MAKESTRX(self):
        """Returns the path to the makestr.x executable."""
        return self.property_get("MAKESTRX")

    @property
    def GETKPTS(self):
        """Returns the path to the makestr.x executable."""
        return self.property_get("GETKPTS")

    @property
    def UNCLEX(self):
        """Returns the path to the compiled UNCLE executable."""
        return self.property_get("UNCLEX")
        
    @property
    def PWX(self):
        """Returns the path to the compiled UNCLE executable."""
        return self.property_get("PWX")

    @property 
    def POTCAR_DIR(self):
        """Directory for all the POTCARs; can be relative user path."""
        return self.property_get("POTCAR_DIR")

    @property
    def INPUT_DIR(self):
        """Directory that contains all the input file templates referenced in the
        ANCE modules."""
        from os import path
        return path.expanduser(self.property_get("INPUT_DIR"))

    @property
    def TEMPL_DIR(self):
        """Directory containing XML conversion templates for input files."""
        return self.property_get("TEMPL_DIR")

    @property
    def RUN_DIR(self):
        """The folder on compute with faster HDD to speed up computation times."""
        return self.property_get("RUN_DIR")

    @property
    def STORE_DIR(self):
        """The directory on the archival storage to move run results to."""
        return self.property_get("STORE_DIR")

    @property
    def ENUMX(self):
        """The directory on the archival storage to move run results to."""
        return self.property_get("ENUMX")

    @property
    def latpar_path(self):
        """Returns the full path to the latpars input file for an ANCLE vasp build."""
        from os import path
        return path.expanduser(path.join(self.INPUT_DIR, "latpars"))

    def property_get(self, key, default=None):
        if key in self._vardict:
            return self._vardict[key]
        else:
            return default

    def getenvar(self, envar):
        from os import getenv
        """Retrieves the value of an environment variable if it exists."""
        if getenv(envar) is not None:
            self._vardict[envar] = getenv(envar)
        else:
            self._vardict[envar] = None
            

#    def load_xml(self, filepath):
#        """Loads the values of the configuration variables from an XML path."""
#        from os import path
#        import xml.etree.ElementTree as ET
#        #Make sure the file exists and then import it as XML and read the values out.
#        uxpath = path.expanduser(filepath)
#        if path.isfile(uxpath):
#            tree = ET.parse(uxpath)
#            root = tree.getroot()
#
#            for child in root:
#                if child.tag == "var":
#                    self._vardict[child.attrib["name"]] = child.attrib["value"]
      
modules["config"] = _config()
