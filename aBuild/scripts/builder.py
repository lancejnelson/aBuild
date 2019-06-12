

from aBuild.io import read

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-POSCAR": {"help": "When converting a POSCAR to a cif, you can supply the path to the POSCAR here.",'type':str},
    "-enum": {"action": "store_true",
           "help": "Enumerate derivative superstructures."},
    "-write": {"action": "store_true",
           "help": "Build VASP folders."},
    "-setup_train": {"action": "store_true",
           "help": "Build input files for fitting."},
    "-report": {"action": "store_true",
           "help": "Build a report of the results of the calculations."},
    "-cif": {"action": "store_true",
           "help": "Build a cif file for the POSCAR specified by the -POSCAR tag"},
    "-rgk": {"action": "store_true",
           "help": "Do you want to run the getKpoint script?."},
    "-start": {"type": int, "default": 1,
           "help": "Which structure do you want to start with"},
    "-end": {"type": int, "default": None,
           "help": "Which structure do you want to start with"},
    "-s": {"action": "store_true",
           "help": "Run the setup method for each database."},
    "-x": {"action": "store_true",
           "help": ("Submit the job array file for each database that has "
                    "folders ready to run. If --recover is specified, then "
                    "the recovery jobfile is submitted instead.")},
    "-e": {"action": "store_true",
           "help": ("Process the databases that have completed execution "
                    "so that results can be extracted.")},
    "-status": {"action": "store_true",
                "help": ("Determines status of the databases "
                         "based on presence of completed VASP "
                         "directories. Sanity check before `-x`.")},
    "-setup_relax": {"action": "store_true",
                "help": ("Setup input files for an mlp relax run "
                         "based on presence of completed VASP "
                         "directories. Sanity check before `-x`.")},
    "-setup_select_add": {"action": "store_true",
                "help": ("Setup input files for an mlp select-add run "
                             "directories. Sanity check before `-x`.")},
    "-add": {"action": "store_true",
                "help": ("Setup input files for an mlp relax run "
                         "based on presence of completed VASP "
                         "directories. Sanity check before `-x`.")},
    
    "--rerun": {"type": int, "default": 0,
                "help": ("Re-run the specified option, even if it has already "
                         "been done before. Higher values re-run at a deeper level.")},
    "--dfilter": {"nargs": "+",
                  "help": ("Specify a list of patterns to match against _database_ "
                           "names that should be *included*.")},
    "--clean": {"default": "default",
                "type": str,
                "help": ("Specify the cleanup level for the database, 'aggressive', "
                         "'default', or 'light'.")},
    "--busy": {"action": "store_true",
               "help": ("Display a list of configurations that haven't "
                        "finished running in DFT yet.")},
    "--recover": {"action": "store_true",
                  "help": ("Creates a jobfile for those configs that didn't "
                           "finish computing so that they can be re-run.")},
    "--dryrun": {"action": "store_true",
                  "help": ("For execution, doesn't actually submit any jobs, "
                           "just prints what would be done.")}
    }

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse
    import sys
    #    from matdb import base
    pdescr = "aBuild Constructor"
    parser = argparse.ArgumentParser(description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)
        
        #    args = base.exhandler(examples, parser)
        #if args is None:
        #return

    return parser.parse_args()


def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    if args is None:
        return

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    from aBuild import Controller
    cdb = Controller(args.dbspec)
    print(args.POSCAR,'POSCAR')
    if args.enum:
        cdb.enumerate('trainingset')
    if args.write:
        cdb.setup_training_set(runGetKpoints = args.rgk)

    if args.setup_train:
        cdb.setupHandler('mtp',"setup_train")
    if args.setup_relax:
        cdb.setupHandler('mtp',"setup_relax",start = args.start,end = args.end)
    if args.setup_select_add:
        cdb.setupHandler('mtp',"setup_select_add")

#    if args.setup_train:
#        cdb.setup_training_input()
    if args.status:
        cdb.statusReport()
#########################
#THIS IS WHAT I'VE ADDED
    if args.cif:
        cdb.generate_cif(args.POSCAR)
#########################
#    if args.setup_relax:
#        cdb.setup_relax_input()
#    if args.setup_select_add:
#        cdb.setup_select_input()
    if args.add:
        cdb.augmentTraining()
    if args.report:
        cdb.gatherResults()
#    if args["s"]:
#        cdb.setup(args["rerun"], args["dfilter"])
#    if args["x"]:
#        cdb.execute(args["recover"], args["dfilter"], dryrun=args["dryrun"])
#    if args["e"]:
#        cdb.extract(args["dfilter"], cleanup=args["clean"])    
#
#    if args["recover"] and not args["x"]:
#        cdb.recover(args["rerun"], args["dfilter"])
#        
#    if args["status"]:
#        cdb.status(args["busy"], args["dfilter"])
        
if __name__ == '__main__': # pragma: no cover
    parser = _parser_options()
    run(parser)


#path = './'
#file = 'AgPt'
#
##results = read(path,file)
#
#inputFile = 'AgPt'
#myController = Controller(inputFile)
#
#myController.enumerate()




