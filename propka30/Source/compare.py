import  sys, math, string
from lib import int2roman, convertResidueCode


def compareFoldingContributions(target=None, template=None, options=None):
    """
    1. check if pKa values are calculated
    2. calculate contributions to dG_fold
    3. read alignment
    4. calculate the difference
    5. printout sorted result
    """
    from mutate import readAlignmentFiles

    # checking that pKa values are available
    checkDonePKA(target, template)

    # reading alignment files
    alignment = readAlignmentFiles(filenames=options.alignment, mesophile=target.name, options=options)
    names = makeNameList(name=target.name, alignment=alignment, options=options)
    print(names)
    printAlignment(names=names, alignment=alignment)

    # setup contribution differences
    positions = makeFoldingEnergyDifferences(target, template, alignment, names=names, options=options)

    # print out contribution differences
    printFoldingEnergyDifferences(positions, names=names, template=template, options=options)

    # print out suggested mutations
    suggestMutations(positions, names=names, template=template, options=options)

    return  None



def makeNameList(name=None, alignment=None, options=None):
    """
    make a list of pdbcodes in 'alignment'
    """
    names = [name]
    for key in alignment.keys():
      if key not in names:
        names.append(key)

    return  names


def checkDonePKA(target=None, template=None, alignment=None, options=None):
    """
    check if pKa values have been calculated for target and template
    """
    if target.status['done pka'] == False:
      print("please calculate pKa values for target %s before comparing" % (target.name))
    if template.status['done pka'] == False:
      print("please calculate pKa values for template %s before comparing" % (template.name))

    if target.status['done pka'] == True and template.status['done pka'] == True:
      if False:
        print("target = %s, template = %s : pKas done" % (target.name, template.name))
    else:
      sys.exit(8)

    return


def printAlignment(names=None, alignment=None):
    """
    printing out alignment
    """
    str = \
"""
sequence alignment:
                        1         2         3         4         5         6         7         8         9        10
               1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
               .........|.........|.........|.........|.........|.........|.........|.........|.........|.........|
"""
    str = str[:-1]
    print(str)
    for key1 in alignment.keys():
      for key2 in alignment[key1].keys():
        str = "    %s 100 %s 0" % (key2, "%")
        index = 0
        for code in alignment[key1][key2]['sequence']:
          index += 1
          if index%100 == 0:
            str += "\n               "
          str += "%s" % (code)
        print(str)

      #print(alignment[key])

    return


def makeFoldingEnergyDifferences(target, template, alignment, names=None, options=None):
    """
    making an array with folding energy differences and related information
    """
    number_of_positions = len(alignment[template.name][template.name]['sequence'])

    # setting up the list of sequence positions, initiated with dictionary
    positions = [{target.name: {}, template.name: {}} for position in range(number_of_positions)]

    # setting residue labels and contributions to sequence positions
    for protein in [target, template]:
      i_position = 0
      for chain in protein.chains:
        for residue in chain.residues:
          while alignment[template.name][protein.name]['sequence'][i_position] in ["-", "?"]:
            position = positions[i_position]
            position[protein.name]['label']        = "   gap   "
            position[protein.name]['contribution'] = 0.00
            i_position += 1
          if residue.resName == "N+ ":
            position = positions[0]
          else:
            position = positions[i_position]
            i_position += 1
          position[protein.name]['label']        = residue.label
          position[protein.name]['contribution'] = residue.calculateFoldingEnergy(options=options)

    # filling position differences
    print("\n unsorted contributions to the folding energy: (kcal/mol)")
    print("-"*64)
    i_position = 0
    for position in positions:
      str  = "%5d  " % (i_position)
      for name in names:
        position[name]['difference'] = position[name]['contribution'] - position[names[0]]['contribution']
        str += "  %s" % (position[name]['label'])
        str += " %6.2lf" % (position[name]['contribution'])
        str += "    "
      print(str)
      i_position += 1


    return  positions


def sortAccordingToMin(positions, key=None, options=None):
    """
    making an array with references in order
    """
    if key == None:
      key = '2vuj'
    min = []
    for position in positions:
      inserted = False
      for i in range(len(min)):
        if position[key]['difference']  <  min[i][key]['difference']:
          min.insert(i, position)
          inserted = True
          break
      if inserted == False:
        min.append(position)

    return  min


def printFoldingEnergyDifferences(positions, names=None, template=None, options=None):
    """
    making an array with folding energy differences and related information
    """
    print("\n the most stabilizing differences: (kcal/mol)")
    print("-"*64)
    sorted_positions = sortAccordingToMin(positions, key=names[1])

    i_position = 0
    for position in sorted_positions:
      target_label   = position[names[0]]['label']
      template_label = position[names[1]]['label']
      difference     = position[names[1]]['difference']
      str  = "%5d    %s   ->   %s  %6.2lf   " % (i_position, target_label, template_label, difference)
      if difference < -0.50:
        str += suggestMutation(positions=positions, label=template_label, names=names, template=template, options=options)
      print(str)
      i_position += 1

    return


def suggestMutation(positions=None, label=None, names=None, template=None, options=None):
    """
    making a suggestion to this mutation
    """
    target_label   = getCorrespondingResidueLabel(positions, key=names[1], label=label)
    template_label = label
    if   target_label == "   gap   " or template_label == "   gap   ":
      return  ""
    elif target_label[:3] in ["C- ", "N+ "] or template_label[:3] in ["C- ", "N+ "]:
      return  ""
    else:
      determinant_labels = []
      # set self-mutation
      mutation = makeMutationAddendum(target=target_label, template=template_label)
      residue  = template.getResidue(label=template_label)
      # get determinants
      for determinants in [residue.determinants[0], residue.determinants[2]]:
        for determinant in determinants:
          if determinant.label not in determinant_labels:
            determinant_labels.append(determinant.label)
      # set together mutation
      for determinant_label in determinant_labels:
        target_label = getCorrespondingResidueLabel(positions, key=names[1], label=determinant_label)
        if target_label != "   gap   " and target_label[:3] not in ["C- ", "N+ "]:
          mutation += "/%s" % (makeMutationAddendum(target=target_label, template=determinant_label))

      return  mutation


def suggestMutations(positions, names=None, template=None, options=None):
    """
    making an array with folding energy differences and related information
    """

    print("\n suggesting mutations")
    print("-"*64)
    sorted_positions = sortAccordingToMin(positions, key=names[1])

    number = 0; weight = None
    for position in sorted_positions:
      target_label   = position[names[0]]['label']
      template_label = position[names[1]]['label']
      if   position[names[1]]['difference'] > -0.5:
        break
      elif target_label == "   gap   " or template_label == "   gap   ":
        """ do nothing """
      elif target_label[:3] in ["C- ", "N+ "] or template_label[:3] in ["C- ", "N+ "]:
        """ do nothing """
      else:
        number +=1; roman = int2roman(number)
        determinant_labels = []
        mutation = makeMutationAddendum(target=target_label, template=template_label)
        #1. get residue
        residue = template.getResidue(label=template_label)
        #2. get determinants
        for determinants in [residue.determinants[0], residue.determinants[2]]:
          for determinant in determinants:
            if determinant.label not in determinant_labels:
              determinant_labels.append(determinant.label)
        #3. set together mutation
        weight = residue.buried
        for determinant_label in determinant_labels:
          target_label = getCorrespondingResidueLabel(positions, key=names[1], label=determinant_label)
          if target_label != "   gap   " and target_label[:3] not in ["C- ", "N+ "]:
            mutation += "/%s" % (makeMutationAddendum(target=target_label, template=determinant_label))
            residue = template.getResidue(label=determinant_label)
            weight += residue.buried
        weight = int( 100.*weight/(len(determinant_labels) + 1) )
        print(" %-5s %3d%2s   %s" % (roman, weight, "%", mutation))

    return


def getCorrespondingResidueLabel(positions=None, key=None, label=None):
    """
    going through positions and return the target label corresponding to the template label 'label'
    """
    for position in positions:
      if position[key]['label'] == label:
        return  position["1xnb"]['label']


def makeMutationAddendum(target=None, template=None, options=None):
    """
    converting two residue labels to a mutation tag
    """
    if target == "   gap   " or template == "   gap   ":
      return  None
    else:
      code1, resName = convertResidueCode(resName=target[:3])
      code2, resName = convertResidueCode(resName=template[:3])
      resNumb = int(target[3:7])
    
      return  "%s%d%s" % (code1, resNumb, code2)


