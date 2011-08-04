import string
import Source.lib as lib
import Source.calculator as calculate


def interactionMatrix(interaction):
    """
    printing out all information in resInfo
    """
    keys = ["COO", "CYS", "TYR", "HIS", "N+ ", "LYS", "ARG", "ROH", "AMD", "TRP"]

    print("interaction matrix:")
    for key1 in keys:
      str = "%6s:" % (key1)
      for key2 in keys:
        do_pair, iterative = interaction[key1][key2]
        if   do_pair == True and iterative == True:
          str += "%3s" % ("I")
        elif do_pair == True and iterative == False:
          str += "%3s" % ("N")
        else:
          str += "%3s" % ("-")

      print(str)


def printResInfo(resInfo):
    """
    printing out all information in resInfo
    """
    print("in resInfo:")
    for key1 in resInfo.keys():
      print(" --- %s ---" % (key1))
      for key2 in resInfo[key1].keys():
        print(key2, resInfo[key1][key2])


def printCooArgAtomDistances(residue_coo, residue_arg):
    """
    printing out all COO and ARG distances for debugging, picking closest + runner-up
    """
    for atom_coo in residue_coo.makeDeterminantAtomList(residue_arg.resName):
      for atom_arg in residue_arg.makeDeterminantAtomList(residue_coo.resName):
        distance = calculate.InterAtomDistance(atom_coo, atom_arg)
        print("%3s %3s %6.2lf" % (atom_coo.name, atom_arg.name, distance))


def printBackBoneAtoms(list):
    """
    Prints out determinant information for debugging
    """
    print("  --- debug back-bone atom list --- ")
    for atoms in list:
        label1 = "%s%4d%2s" % (atoms[0].resName, atoms[0].resNumb, atoms[0].chainID)
        label2 = "%s%4d%2s" % (atoms[1].resName, atoms[1].resNumb, atoms[1].chainID)
        print("%s - %s" % (label1, label2))


def printResidues(residue_list):
    """
    Prints out determinant information for debugging
    """
    print("  --- debug residue list --- ")
    for residue in residue_list:
        residue.printLabel()


def printIterativeDeterminants(all_determinants):
    """
    Prints out determinant information for debugging
    """

    if True:
      print(" --- Iterative determinants ---")
      print("%28s%8s" % ("H-bond", "Coulomb") )
      for determinant in all_determinants:
        pair   =  determinant[0]
        values =  determinant[1]
        annihilation = determinant[2]
        residue1 = pair[0]
        residue2 = pair[1]
        str  = ""
        str += " %s %s" % (residue1.label, residue2.label)
        str += " %6.2lf %6.2lf" % (values[0], values[1])
        print(str)
    else:
      # priting out ALL types of determinants - mainly debugging
      print("  --- debug iterative determinants --- ")
      print("len(sidechain_determinants) = %d" % (len(all_determinants[0])))
      print("len(backbone_determinants)  = %d" % (len(all_determinants[1])))
      print("len(coulomb_determinants)   = %d" % (len(all_determinants[2])))
    
      sidechain_determinants = all_determinants[0]
      backbone_determinants  = all_determinants[1]
      coulomb_determinants   = all_determinants[2]
    
      if True:
        for type in range(0,3):
          determinants = all_determinants[type]
          if   type == 0:
            print("Iterative side-chain interactions:")
          elif type == 1:
            print("Iterative back-bone  interactions:")
          elif type == 2:
            print("Iterative Coulomb    interactions:")
          for determinant in determinants:
            pair   =  determinant[0]
            values =  determinant[1]
    
            pair   = determinant[0]
            value  = determinant[1]
            residue1 = pair[0]
            residue2 = pair[1]
            str  = ""
            str += "'%s'  '%s'" % (residue1.label, residue2.label)
            str += "   "
            str += "%6.2lf  %6.2lf" % (value, value)
            print(str)


def printAlignment(alignment):
    """
    Prints out alignment information for debugging
    """
    for key in alignment.keys():
      print( " --- %s ---" % (key) )
      for key2 in alignment[key].keys():
        print("%s %5d%2s" % (alignment[key][key2]["name"], alignment[key][key2]["resNumb"], alignment[key][key2]["chainID"]))
        print("%s\n" % (alignment[key][key2]["sequence"]))


def printAllAtoms(protein):
    """
    Prints out determinant information for debugging
    """
    print("Not implemented yet")
