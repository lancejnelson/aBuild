

class Job:

    def __init__(self,settings,exdir,command,arrayStart=None,arrayEnd=None):
        self.jobSettings = settings
        self.jobSettings["execution_path"] = exdir
        self.jobSettings["exec_path"] = command
        if arrayStart is not None:
            self.jobSettings["array_start"] = arrayStart
            self.jobSettings["array_end"] = arrayEnd
            self.arrayJob = True
        else:
            self.arrayJob = False


    def write_jobfile(self,target):


        from jinja2 import Environment, PackageLoader  # Package for building files from a template
        from os import path
        env = Environment(loader=PackageLoader('aBuild', 'templates'))
        template = env.get_template(self.jobSettings["template"])
#        if self.arrayJob:
#            template = env.get_template("run_array_ml.sh")
#        else:
#            template = env.get_template("run_mtp_ml.sh")

#        target = 'jobscript.sh'
        with open(target,'w') as f:
            f.write(template.render(**self.jobSettings))

        
