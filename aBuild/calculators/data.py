"""Data about elemental and alloy systems that doesn't change with time.
"""
latpars ={"H": 3.75,"He": 3.57,"Li": 3.49,"Be": 2.29,"B": 8.73,"C": 3.57,"N": 4.039,
          "O": 6.83,"Ne": 4.43,"Na": 4.23,"Mg": 3.21,"Al": 4.05,"Si": 5.43,"P": 7.17,
          "S": 10.47,"Cl": 6.24,"Ar": 5.26,"K": 5.23,"Ca": 5.58,"Sc": 3.31,"Ti": 2.95,
          "V": 3.02,"Cr": 2.88,"Mn": 8.89,"Fe": 2.87,"Co": 2.51,"Ni": 3.52,"Cu": 3.61,
          "Zn": 2.66,"Ga": 4.51,"Ge": 5.66,"As": 4.13,"Se": 4.36,"Br": 6.67,"Kr": 5.72,
          "Rb": 5.59,"Sr": 6.08,"Y": 3.65,"Zr": 3.23,"Nb": 3.3,"Mo": 3.15,"Tc": 2.74,
          "Ru": 2.7,"Rh": 3.8,"Pd": 3.89,"Ag": 4.09,"Cd": 2.98,"In": 4.59,"Sn": 5.82,
          "Sb": 4.51,"Te": 4.45,"I": 7.27,"Xe": 6.2,"Cs": 6.05,"Ba": 5.02,"Hf": 3.2,
          "Ta": 3.31,"W": 3.16,"Re": 2.76,"Os": 2.64,"Ir": 3.84,"Pt": 3.92,"Au": 4.08,
          "Hg": 2.99,"Tl": 3.46,"Pb": 4.95,"Bi": 4.75}

element_volume ={"H":37.2958,"He":32.1789,"Li":21.2543,"Be":8.49323,"B":7.24205,"C":5.68741,
                 "N":46.6002,"O":22.2802,"F":17.0258,"Ne":21.7346,"Na":23.2596,"Mg":23.3928,
                 "Al":16.6075,"Si":7.8511,"P":9.1459,"S":17.1672,"Cl":35.2074,"Ar":36.3829,
                 "K":71.5278,"Ca":43.4353,"Sc":25.6478,"Ti":18.1565,"V":13.7718,"Cr":11.9439,
                 "Mn":19.3207,"Fe":11.82,"Co":11.1838,"Ni":10.9036,"Cu":11.7615,"Zn":13.311,
                 "Ga":18.4496,"Ge":19.3638,"As":20.4270,"Se":58.6173,"Br":33.3170,"Kr":46.7873,
                 "Rb":87.3384,"Sr":56.1889,"Y":33.0792,"Zr":23.8327,"Nb":17.9685,"Mo":15.6279,
                 "Tc":14.5458,"Ru":13.9206,"Rh":13.718,"Pd":14.716,"Ag":17.1045,"Cd":18.7161,
                 "In":26.6861,"Sn":29.3238,"Sb":27.1733,"Te":62.3227,"I":24.3807,"Xe":59.582,
                 "Cs":110.723,"Ba":63.253,"Hf":23.1748,"Ta":18.1323,"W":15.7772,"Re":14.8694,
                 "Os":14.5485,"Ir":14.1558,"Pt":15.0591,"Au":16.9793,"Hg":27.6914,"Tl":29.2949,
                 "Pd":30.3218,"Bi":31.2849}

volumeperatom = {"Au": latpars["Au"]**3/4, "Ag": latpars["Ag"]**3/4, "Cu": latpars["Cu"]**3/4}


# This routine is to calculate the appropriate lattice
# parameter for the crystal.  We do this by ensuring that
# the volume per atom is right.  There are two ways: One is
# the way I thought up, and the other is Wiley's way.  Wiley's
# way is probably right because he has thought of it for longer.
# My way:
#     1- Find the concentration-weighted average of the pure atomic
#        volumes.  This is the target volume for the crystal.
#     2- The lattice parameter is then found by dividing the target
#        atomic volume by the current atomic volume and then taking
#        cubed root.
# Wiley's way:
#     1- Find the lattice paramter that the crystal would have if it was
#        pure A(B,C,etc).
#     2- Perform a concentration-weighted average over the lattice parameters.
#
#     The two methods are mathematically very similar but if the
#     mismatch is large, the difference can get larger.  
def vegardsVolume(elements,atom_counts,volume):
    nAtoms = sum(atom_counts)
    nTypes = len(atom_counts)
    concentrations = [x/nAtoms for x in atom_counts]

    #  The next two lines are Wiley's way.
    # Find lat pars if crystal was pure.
    purelatPars = [( element_volume[elements[x]]/(volume/nAtoms) )**(1./3.) for x in range(nTypes)]
    # Concentration-weighted average over lattice parameters
    wiley = sum([purelatPars[x] * concentrations[x] for x in range(nTypes)])
    

    # The next two lines are my way.
    # concentration-weighted average of volumes
    avgVolume = sum([element_volume[elements[x]] *concentrations[x] for x in range(nTypes)] )
    # Calculate lattice parameter
    mine = ( avgVolume/(volume/nAtoms) )**(1/3.)
#    print("My approach: {}.  Wiley's approach: {}".format(mine,wiley))
    return wiley


# Should not be used! Was foolish to ever use in the first place.  Really only works
# when the constituents belong to the same lattice as the crystal.
#def vegard(elements, concs):
#    """Interpolates the lattice parameters for the given alloy system using
#    vegard's law.
#
#    Args:
#        elements (list): of `str` chemical symbols; keys in :data:`latpars`.
#        concs (list): of concentrations for each element in `elements`.
#    """
#    conctot = sum(concs)
#    lat = [latpars[e] for e in elements]
#    return sum(a*b/conctot for a, b in zip(lat, concs))
