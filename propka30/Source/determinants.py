import math, time

import Source.iterative as iterative
import Source.lib as lib
#import Source.debug as debug
import Source.calculator as calculate
from   Source.determinant import Determinant


def setDeterminants(propka_residues, version=None, options=None):
    """
    adding side-chain and coulomb determinants/perturbations to all residues - note, backbone determinants are set separately
    """
    #debug.printResidues(propka_residues)
    iterative_interactions = []
    # --- NonIterative section ---#
    for residue1 in propka_residues:
      for residue2 in propka_residues:

        if residue1 == residue2:
            break

        distance = calculate.InterResidueDistance(residue1, residue2)

        if distance < version.sidechain_cutoff or distance < version.coulomb_cutoff[1]:
          do_pair, iterative_interaction = version.interaction[residue1.resType][residue2.resType]
         
          if do_pair == True:
            if iterative_interaction == True:
              iterative.addtoDeterminantList(residue1, residue2, distance, iterative_interactions, version=version)
              #print "%s - %s I" % (residue1.label, residue2.label)
            else:
              addDeterminants(residue1, residue2, distance, version=version)
              #print "%s - %s" % (residue1.label, residue2.label)
          else:
            """ False - don't do this at home folks """

    # --- Iterative section ---#
    #debug.printIterativeDeterminants(iterative_interactions)
    iterative.addDeterminants(iterative_interactions, version, options=options)


def addDeterminants(residue1, residue2, distance, version=None):
    """
    adding determinants/perturbations, distance(R1, R2) < coulomb_cutoff always
    """

    # side-chain determinant
    if distance < version.sidechain_cutoff:
        # Currently we don't want any hydrogen bonds to ligands
        if 'ligand' not in [residue1.type, residue2.type]:
            addSidechainDeterminants(residue1, residue2, version)

    do_coulomb = version.checkCoulombPair(residue1, residue2, distance)

    # Coulomb determinant
    if do_coulomb == True:
        addCoulombDeterminants(residue1, residue2, distance, version)


def addSidechainDeterminants(residue1, residue2, version=None):
    """
    adding side-chain determinants/perturbations
    Note, resNumb1 > resNumb2
    """
    distance = 999.0
    closest_atom1 = None
    closest_atom2 = None
    atoms1 = residue1.makeDeterminantAtomList(residue2.resName, version=version)
    atoms2 = residue2.makeDeterminantAtomList(residue1.resName, version=version)
    for atom1 in atoms1:
        for atom2 in atoms2:
            # select the smallest inter-atom distance
            current_distance = calculate.InterAtomDistance(atom1, atom2)
            if current_distance < distance:
                closest_atom1 = atom1
                closest_atom2 = atom2
                distance = current_distance

    dpka_max, cutoff = version.SideChainParameters[residue1.resType][residue2.resType]
    if distance < cutoff[1]:
        if   residue2.resType in version.angularDependentSideChainInteractions:
          atom3 = residue2.getThirdAtomInAngle(closest_atom2)
          distance, f_angle, nada = calculate.AngleFactorX(closest_atom1, closest_atom2, atom3)
        elif residue1.resType in version.angularDependentSideChainInteractions:
          atom3 = residue1.getThirdAtomInAngle(closest_atom1)
          distance, f_angle, nada = calculate.AngleFactorX(closest_atom2, closest_atom1, atom3)
        else:
          # i.e. no angular dependence
          f_angle = 1.0

        weight = version.calculatePairWeight(residue1.Nmass, residue2.Nmass)
        exception, value = version.checkExceptions(residue1, residue2)
        #exception = False # circumventing exception
        if exception == True:
            """ do nothing, value should have been assigned """
            #print(" exception for %s %s %6.2lf" % (residue1.label, residue2.label, value))
        else:
            value = version.calculateSideChainEnergy(distance, dpka_max, cutoff, weight, f_angle)
        if residue1.Q == residue2.Q:
          # acid pair or base pair
          if residue1.pKa_mod < residue2.pKa_mod:
            newDeterminant1 = Determinant(residue2.label, -value)
            newDeterminant2 = Determinant(residue1.label,  value)
          else:
            newDeterminant1 = Determinant(residue2.label,  value)
            newDeterminant2 = Determinant(residue1.label, -value)
        else:
          newDeterminant1 = Determinant(residue2.label, value*residue1.Q)
          newDeterminant2 = Determinant(residue1.label, value*residue2.Q)
        if residue1.resName not in version.exclude_sidechain_interactions:
          residue1.determinants[0].append(newDeterminant1)
        if residue2.resName not in version.exclude_sidechain_interactions:
          residue2.determinants[0].append(newDeterminant2)


def addCoulombDeterminants(residue1, residue2, distance, version):
    """
    adding NonIterative Coulomb determinants/perturbations
    """
    weight = version.calculatePairWeight(residue1.Nmass, residue2.Nmass)
    value  = version.calculateCoulombEnergy(distance, weight)
    Q1 = residue1.Q
    Q2 = residue2.Q

    # assigning the Coulombic interaction
    if   Q1 < 0.0 and Q2 < 0.0:
        """ both are acids """
        addCoulombAcidPair(residue1, residue2, value)
    elif Q1 > 0.0 and Q2 > 0.0:
        """ both are bases """
        addCoulombBasePair(residue1, residue2, value)
    else:
        """ one of each """
        addCoulombIonPair(residue1, residue2, value)


def addCoulombAcidPair(object1, object2, value):
    """
    Adding the Coulomb interaction (an acid pair):
    the higher pKa is raised
    """
    label1 = object1.label
    label2 = object2.label
    if object1.pKa_mod > object2.pKa_mod:
        newDeterminant = Determinant(label2, value)
        object1.determinants[2].append(newDeterminant)
    else:
        newDeterminant = Determinant(label1, value)
        object2.determinants[2].append(newDeterminant)


def addCoulombBasePair(object1, object2, value):
    """
    Adding the Coulomb interaction (a base pair):
    the lower pKa is lowered
    """
    label1 = object1.label
    label2 = object2.label
    if object1.pKa_mod < object2.pKa_mod:
        newDeterminant = Determinant(label2, -value)
        object1.determinants[2].append(newDeterminant)
    else:
        newDeterminant = Determinant(label1, -value)
        object2.determinants[2].append(newDeterminant)


def addCoulombIonPair(object1, object2, value):
    """
    Adding the Coulomb interaction (an acid-base pair):
    the pKa of the acid is lowered & the pKa of the base is raised
    """
    label1 = object1.label
    label2 = object2.label

    # residue1
    Q1 = object1.Q
    newDeterminant = Determinant(label2, Q1*value)
    object1.determinants[2].append(newDeterminant)

    # residue2
    Q2 = object2.Q
    newDeterminant = Determinant(label1, Q2*value)
    object2.determinants[2].append(newDeterminant)




def setIonDeterminants(protein, version=None):
    """
    adding ion determinants/perturbations
    """
    ionizable_residues = lib.residueList("propka1")
    for residue in protein.propka_residues:
      if residue.resName in ionizable_residues:
        for ion in protein.residue_dictionary["ION"]:
          distance = calculate.InterResidueDistance(residue, ion)
          if distance < version.coulomb_cutoff[1]:
            label  = "%s%4d%2s" % (ion.resName, ion.resNumb, ion.chainID)
            weight = version.calculatePairWeight(residue.Nmass, ion.Nmass)
            # the pKa of both acids and bases are shifted up by negative ions (and vice versa)
            value  =  (-ion.Q) * version.calculateCoulombEnergy(distance, weight)
            newDeterminant = Determinant(label, value)
            residue.determinants[2].append(newDeterminant)


def setBackBoneDeterminants(backbone_interactions, version=None):
    """
    adding back-bone determinants/perturbations
    Angle: atom1 -- atom2-atom3
    backbone_interactions = [[acids, NH], [bases, CO]]
    changing the code with minimum effect of method calls
    """
    setBackBoneAcidDeterminants(backbone_interactions[0], version=version)
    setBackBoneBaseDeterminants(backbone_interactions[1], version=version)


def setBackBoneAcidDeterminants(data_clump, version=None):
    """
    adding back-bone determinants/perturbations for acids:
    Angle: atom1 -- atom2-atom3, i.e. COO -- H-N
    data_clump = [acids, NH]
    """
    residues, interactions = data_clump
    for residue in residues:
      if residue.location != "BONDED":
        dpKa_max, cutoff = version.BackBoneParameters[residue.resType]
        for interaction in interactions:
          atom2 = interaction[1]
          atom3 = interaction[0]
          atoms = residue.makeDeterminantAtomList("back-bone", version=version)
          shortest_distance = 999.
          for atom in atoms:
              distance = calculate.InterAtomDistance(atom, atom2)
              if distance < shortest_distance:
                  shortest_distance = distance
                  atom1 = atom
          distance, f_angle, nada = calculate.AngleFactorX(atom1, atom2, atom3)
          if distance < cutoff[1] and f_angle > 0.001:
              label = "%s%4d%2s" % (atom2.resName, atom2.resNumb, atom2.chainID)
              value = residue.Q * calculate.HydrogenBondEnergy(distance, dpKa_max, cutoff, f_angle)
              newDeterminant = Determinant(label, value)
              residue.determinants[1].append(newDeterminant)

        
def setBackBoneBaseDeterminants(data_clump, version=None):
    """
    adding back-bone determinants/perturbations for bases:
    Angle: atom1 -- atom2-atom3, i.e. C=O -- H-N(HIS)
    data_clump = [bases, CO]
    """
    residues, interactions = data_clump
    for residue in residues:
      if residue.location != "BONDED":
          dpKa_max, cutoff = version.BackBoneParameters[residue.resType]
          for interaction in interactions:
            distance = 999.
            atom1 = interaction[1]
            atoms = residue.makeDeterminantAtomList("back-bone", version=version)
            for atom in atoms:
              current_distance = calculate.InterAtomDistance(atom1, atom)
              if current_distance < distance:
                  atom2 = atom
                  distance = current_distance
            if distance < cutoff[1]:
              if residue.resType in version.angularDependentSideChainInteractions:
                atom3 = residue.getThirdAtomInAngle(atom2)
                distance, f_angle, nada = calculate.AngleFactorX(atom1, atom2, atom3)
              else:
                f_angle = 1.0
              if f_angle > 0.001:
                # add determinant
                label = "%s%4d%2s" % (atom2.resName, atom2.resNumb, atom2.chainID)
                value = residue.Q * calculate.HydrogenBondEnergy(distance, dpKa_max, cutoff, f_angle)
                newDeterminant = Determinant(label, value)
                residue.determinants[1].append(newDeterminant)

