
def resName2Type(resName=None):
    """
    definition of which parameter-group each residues belongs to
    """

    resType = {'C- ': "COO",
               'ASP': "COO",
               'GLU': "COO",
               'HIS': "HIS",
               'CYS': "CYS",
               'TYR': "TYR",
               'LYS': "LYS",
               'ARG': "ARG",
               'N+ ': "N+ ",
               'SER': "ROH",
               'THR': "ROH",
               'ASN': "AMD",
               'GLN': "AMD",
               'TRP': "TRP"}

    if resName in resType:
        return resType[resName]
    else:
        return resType



def getQs(resName=None):
    """
    Returns a dictionary with residue charges
    """
    Q  =      {'COO': -1.0,
               'ASP': -1.0,
               'GLU': -1.0,
               'C- ': -1.0,
               'TYR': -1.0,
               'CYS': -1.0,
               'HIS':  1.0,
               'LYS':  1.0,
               'ARG':  1.0,
               'N+ ':  1.0}

    if resName in Q:
        return Q[resName]
    else:
        return Q


def pKa_mod(resName=None):
    """
    returns a dictionary with model/water pKa values
    """
    pKa_mod = {'C- ':      3.20,
               'ASP':      3.80,
               'GLU':      4.50,
               'HIS':      6.50,
               'CYS':      9.00,
               'TYR':     10.00,
               'LYS':     10.50,
               'ARG':     12.50,
               'N+ ':      8.00,
               'default': 20.00}

    if   resName == None:
        return  pKa_mod
    elif resName in pKa_mod:
        return  pKa_mod[resName]
    else:
        # generic value for 'uninteresting' residues, e.g. ASN, GLN
        return  20.00


def getInteraction():
    """
      matrix for propka interactions; Note that only the LOWER part of the matrix is used!
    
      'N'   non-iterative interaction
      'I'   iterative interaction
      '-'   no interaction
    """

    #                     COO  CYS  TYR  HIS   N+  LYS  ARG
    side_chain = {'COO': ["I", "I", "N", "I", "N", "N", "N"],
                  'CYS': ["I", "I", "N", "I", "N", "N", "N"],
                  'TYR': ["N", "N", "I", "I", "I", "I", "N"],
                  'HIS': ["I", "I", "I", "I", "N", "N", "N"],
                  'N+ ': ["N", "N", "I", "N", "I", "N", "N"],
                  'LYS': ["N", "N", "I", "N", "N", "I", "N"],
                  'ARG': ["N", "N", "N", "N", "N", "N", "I"],
                  'ROH': ["N", "N", "N", "-", "-", "-", "-"],
                  'AMD': ["N", "N", "N", "N", "N", "N", "N"],
                  'TRP': ["N", "N", "N", "-", "-", "-", "-"]}

    return  side_chain



# ------- Coulomb parameters --------- #


def getCoulombParameters(label=None):
    """
    storage of Coulomb default parameters
    """
    CoulombParameters = {}
    CoulombParameters['Linear']                 = {'cutoff':             [4.0,  7.0],
                                                   'max_dpka':                  2.40,
                                                   'scaled':                    True,
                                                  }
    CoulombParameters['Coulomb']                = {'cutoff':             [4.0, 10.0],
                                                   'diel':                     80.00,
                                                   'scaled':                    True,
                                                  }
    CoulombParameters['DistanceScaledCoulomb']  = {'cutoff':             [4.0, 10.0],
                                                   'diel':             [30.0, 160.0],
                                                   'scaled':                   False,
                                                  }

    if label in CoulombParameters:
      return CoulombParameters[label]
    else:
      return CoulombParameters



# ------- Desolvation parameters --------- #


def getDesolvationParameters(label=None):
    """
    storage of desolvation default parameters
    """
    DesolvationParameters = {}
    DesolvationParameters['propka2']            = {'allowance':               400.00,
                                                   'prefactor':                -0.01,
                                                   'local':                    -0.07,
                                                   'radii':          getLocalRadii(),
                                                  }
    DesolvationParameters['ContactModel']       = {'allowance':               400.00,
                                                   'prefactor':                -0.01,
                                                   'local':                    -0.07,
                                                   'radii':          getLocalRadii(),
                                                  }
    DesolvationParameters['VolumeModel']        = {'allowance':                 0.00,
                                                   'prefactor':               -13.50,
                                                   'surface':                   0.25,
                                                   'volume': getVanDerWaalsVolumes(),
                                                  }
    DesolvationParameters['ScaledVolumeModel']  = {'allowance':                 0.00,
                                                   'prefactor':               -13.50,
                                                   'surface':                   0.00,
                                                   'volume': getVanDerWaalsVolumes(),
                                                  }

    if label in DesolvationParameters:
      return DesolvationParameters[label]
    else:
      return DesolvationParameters


def getVanDerWaalsVolumes():
    """
    storing relative Van der Waals volumes for volume desolvation models
    """
    #                         relative     volume  radius
    VanDerWaalsVolume = {'C':    1.40,    # 20.58   1.70   all 'C' and 'CA' atoms
                         'C4':   2.64,    # 38.79   2.10   hydrodphobic carbon atoms + unidentified atoms
                         'N':    1.06,    # 15.60   1.55   all nitrogen atoms
                         'O':    1.00,    # 14.71   1.52   all oxygen atoms
                         'S':    1.66,    # 24.43   1.80   all sulphur atoms
                        }

    return VanDerWaalsVolume


def getLocalRadii():
    """
    local radii used in the 'propka2' and 'contact' desolvation models
    """
    local_radius = {'ASP': 4.5,
                    'GLU': 4.5,
                    'HIS': 4.5,
                    'CYS': 3.5,
                    'TYR': 3.5,
                    'LYS': 4.5,
                    'ARG': 5.0,
                    'C- ': 4.5,
                    'N+ ': 4.5}

    return local_radius





# ------- hydrogen-bond parameters --------- #


def getHydrogenBondParameters(type=None):
    """ 
    definitions of default back-bone or side-chain interaction parameters
    IMPORTANT: parameters with assigned to 'None' are given by the reverse
               (e.g. CYS-COO is given by COO-CYS) generated at the end.
    """
    if type == "back-bone":

            # --- new back-bone parameter set ---
            # parameters determining the interaction with back-bone NH or CO groups

            parameters        = {"COO": [-0.80, [2.00, 3.00]],
                                 "CYS": [-0.80, [3.00, 4.00]],
                                 "TYR": [-1.20, [2.20, 3.20]],
                                 "HIS": [ 0.80, [2.00, 3.00]],
                                 "N+ ": [ 0.80, [2.80, 3.80]],
                                 "LYS": [ 0.80, [2.80, 3.80]],
                                 "ARG": [ 0.80, [2.00, 3.00]]}
      
    elif type == "side-chain":


            # --- new side-chain parameter set ---
            # parameters determining the interaction with side-chain NH or CO groups
            # IMPORTANT: parameters with assigned to 'None' are given by the reverse
            # (e.g. CYS-COO is given by COO-CYS) generated at the end.

            parameters = {}
            parameters["COO"] = {"COO": [-0.80, [ 2.50, 3.50]],
                                 "CYS": [-0.80, [ 3.00, 4.00]],
                                 "TYR": [-0.80, [ 2.65, 3.65]],
                                 "HIS": [-0.80, [ 2.00, 3.00]],
                                 "N+ ": [-0.80, [ 2.85, 3.85]],
                                 "LYS": [-0.80, [ 2.85, 3.85]],
                                 "ARG": [-0.80, [ 1.85, 2.85]],
                                 "ROH": [-0.80, [ 2.65, 3.65]],
                                 "AMD": [-0.80, [ 2.00, 3.00]],
                                 "TRP": [-0.80, [ 2.00, 3.00]]}
            parameters["CYS"] = {"COO": None,
                                 "CYS": [-1.60, [ 3.00, 5.00]],
                                 "TYR": [-0.80, [ 3.50, 4.50]],
                                 "HIS": [-1.60, [ 3.00, 4.00]],
                                 "N+ ": [-2.40, [ 3.00, 4.50]],
                                 "LYS": [-1.60, [ 3.00, 4.00]],
                                 "ARG": [-1.60, [ 2.50, 4.00]],
                                 "ROH": [-1.60, [ 3.50, 4.50]],
                                 "AMD": [-1.60, [ 2.50, 3.50]],
                                 "TRP": [-1.60, [ 2.50, 3.50]]}
            parameters["TYR"] = {"COO": None,
                                 "CYS": None,
                                 "TYR": [ 0.80, [ 3.50, 4.50]],
                                 "HIS": [-0.80, [ 2.00, 3.00]],
                                 "N+ ": [-1.20, [ 3.00, 4.50]],
                                 "LYS": [-0.80, [ 3.00, 4.00]],
                                 "ARG": [-0.80, [ 2.50, 4.00]],
                                 "ROH": [-0.80, [ 3.50, 4.50]],
                                 "AMD": [-0.80, [ 2.50, 3.50]],
                                 "TRP": [-0.80, [ 2.50, 3.50]]}
            parameters["HIS"] = {"COO": None,
                                 "CYS": None,
                                 "TYR": None,
                                 "HIS": [ 0.00, [ 0.00, 0.00]],
                                 "N+ ": [ 0.00, [ 0.00, 0.00]],
                                 "LYS": [ 0.00, [ 0.00, 0.00]],
                                 "ARG": [ 0.00, [ 0.00, 0.00]],
                                 "ROH": [ 0.00, [ 0.00, 0.00]],
                                 "AMD": [ 0.80, [ 2.00, 3.00]],
                                 "TRP": [ 0.00, [ 0.00, 0.00]]}
            parameters["N+ "] = {"COO": None,
                                 "CYS": None,
                                 "TYR": None,
                                 "HIS": None,
                                 "N+ ": [ 0.00, [ 0.00, 0.00]],
                                 "LYS": [ 0.00, [ 0.00, 0.00]],
                                 "ARG": [ 0.00, [ 0.00, 0.00]],
                                 "ROH": [ 0.00, [ 0.00, 0.00]],
                                 "AMD": [ 0.00, [ 0.00, 0.00]],
                                 "TRP": [ 0.00, [ 0.00, 0.00]]}
            parameters["LYS"] = {"COO": None,
                                 "CYS": None,
                                 "TYR": None,
                                 "HIS": None,
                                 "N+ ": None,
                                 "LYS": [ 0.00, [ 0.00, 0.00]],
                                 "ARG": [ 0.00, [ 0.00, 0.00]],
                                 "ROH": [ 0.00, [ 0.00, 0.00]],
                                 "AMD": [ 0.00, [ 0.00, 0.00]],
                                 "TRP": [ 0.00, [ 0.00, 0.00]]}
            parameters["ARG"] = {"COO": None,
                                 "CYS": None,
                                 "TYR": None,
                                 "HIS": None,
                                 "N+ ": None,
                                 "LYS": None,
                                 "ARG": [ 0.00, [ 0.00, 0.00]],
                                 "ROH": [ 0.00, [ 0.00, 0.00]],
                                 "AMD": [ 0.00, [ 0.00, 0.00]],
                                 "TRP": [ 0.00, [ 0.00, 0.00]]}
            parameters["ROH"] = {"COO": None,
                                 "CYS": None,
                                 "TYR": None,
                                 "HIS": None,
                                 "N+ ": None,
                                 "LYS": None,
                                 "ARG": None,
                                 "ROH": [ 0.00, [ 0.00, 0.00]],
                                 "AMD": [ 0.00, [ 0.00, 0.00]],
                                 "TRP": [ 0.00, [ 0.00, 0.00]]}
            parameters["AMD"] = {"COO": None,
                                 "CYS": None,
                                 "TYR": None,
                                 "HIS": None,
                                 "N+ ": None,
                                 "LYS": None,
                                 "ARG": None,
                                 "ROH": None,
                                 "AMD": [ 0.00, [ 0.00, 0.00]],
                                 "TRP": [ 0.00, [ 0.00, 0.00]]}
            parameters["TRP"] = {"COO": None,
                                 "CYS": None,
                                 "TYR": None,
                                 "HIS": None,
                                 "N+ ": None,
                                 "LYS": None,
                                 "ARG": None,
                                 "ROH": None,
                                 "AMD": None,
                                 "TRP": [ 0.00, [ 0.00, 0.00]]}


            # updating parameter matrix to full matrix
            keys = parameters.keys()
            for key1 in keys:
              for key2 in keys:
                if key2 not in parameters[key1]:
                  parameters[key1][key2] == [ 0.00, [ 0.00, 0.00]]
                elif parameters[key1][key2] == None:
                  parameters[key1][key2] = parameters[key2][key1]


    else:
      print("cannot determine what type of hydrogen-bonding interactions you want type=\"%s\" ['back-bone', 'side-chain']" % (label))
      sys.exit(9)

    return parameters


