import math
#import data


if False:
    def compareWithExperiment(protein, file=None, list=None, set=None, label="ALL", **argv):
        """
        Compares the calculated pKa values with experiment if it is found in 'experiment'
        """
        str = "# %s\n" % (protein.name)
        file.write(str)
        experimental = data.getExperiment(set=set, name=protein.name)
        for chain in protein.chains:
          for residue in chain.residues:
            key = residue.label
            if label == key or label == "ALL":
              if key in experimental:
                pka_exp  = experimental[key]
                dpka_exp = pka_exp - residue.pKa_mod
                dpka_clc = residue.pKa_pro - residue.pKa_mod
                diff     = residue.pKa_pro - pka_exp
                print( "compare:%s%6.2lf%6.2lf %6.2lf" % (key, residue.pKa_pro, pka_exp, diff) )
                # adding to 'compare' list
                list.append(diff)
                # writing to 'comp.dat' file
                str = "%8.2lf%6.2lf\n" % (dpka_exp, dpka_clc)
                file.write(str)


def makeErrorPlot(points, filename=None):
    """
    do patented error graph data
    """
    # open file if needed
    if filename != None:
      file = open(filename, 'w')

    # setup error-list
    error_list = []
    for i in range(0, 301):
      error_list.append(0.00)

    # loop over points and error bins
    for point in points:
      abs_diff = abs(point)
      i = 0
      while abs_diff > float(i)/100.0 and i < 301:
        error_list[i] += 1.0
        i += 1

    number_of_points = float(len(points))

    # print/write out error fraction
    for i in range(len(error_list)):
      error = float(i)/100.0
      fraction = error_list[i]/number_of_points
      str = "%6.2lf%6.3lf" % (error, fraction)
      if filename == None:
        print(str)
      else:
        str += "\n"; file.write(str)

    if filename != None:
      file = open(filename, 'w')


def calculateMUE(points):
    """
    do mean unsigned error of all points
    """
    number_of_points = len(points)
    abs_sum = 0.00
    for point in points:
        abs_sum += abs(point)
    mue  = abs_sum/number_of_points

    print("mue  =     %6.4lf" % (mue))


def calculateRMSD(points):
    """
    do rmsd of all points
    """
    number_of_points = len(points)
    sqr_sum = 0.00
    for point in points:
        sqr_sum += pow(point, 2)
    rmsd = pow(sqr_sum/number_of_points, 0.5)

    print("rmsd =     %6.4lf (%d)" % (rmsd, number_of_points))


def calculateShift(points):
    """
    calculate an average shift dG_clc - dG_exp
    """
    number_of_points = len(points)
    sum_diff = 0.00
    for pH, experiment, calculation in points:
      sum_diff += experiment - calculation

    return  sum_diff/number_of_points


def printDesolvationDistances(protein):
    """
    printing out wanted distances
    """
    label = "GLU  66 A"
    sum = 0.00; enumerator = 130.0
    center_residue = protein.getResidue(label=label)
    for chain in protein.chains:
      for residue in chain.residues:
        distance = 999.; closest = None
        for atom in residue.atoms:
          dX = atom.x - center_residue.x
          dY = atom.y - center_residue.y
          dZ = atom.z - center_residue.z
          test_distance = math.sqrt(dX*dX + dY*dY + dZ*dZ)
          if residue.label in ["THR  62 A", "VAL  23 A", "LEU  14 A", "VAL  99 A", "ILE  92 A", "LEU  36 A"]:
            if atom.element == "C" and atom.name not in ["C", "CA"]:
              contribution = enumerator/pow(test_distance, 4)
              sum += contribution
              print("%s  %-3s %6.2lf %6.2lf" % (residue.label, atom.name, test_distance, contribution))
          if test_distance < distance and atom.element == "C":
            distance = test_distance
            closest = atom
        if distance <  8.0 and False:
          print("%s %6.2lf  %s" % (residue.label, distance, closest.name))
    print("sum:%8.4lf" % (sum))


def regression(points):
    """
    do linear regression of all points
    """
    sum_x  = 0.00
    sum_y  = 0.00
    sum_xx = 0.00
    sum_xy = 0.00
    sum_yy = 0.00
    number_of_points = len(points)
    for x, y in points:
      print("%6.2lf %6.2lf" % (x, y))
      sum_x  += x
      sum_y  += y
      sum_xx += x*x
      sum_xy += x*y
      sum_yy += y*y

    slope       = (number_of_points*sum_xy - sum_x*sum_y) / (number_of_points*sum_xx - sum_x*sum_x)
    intercept   = (sum_y*sum_xx - sum_x*sum_xy) / (number_of_points*sum_xx - sum_x*sum_x)
    correlation = (number_of_points*sum_xy - sum_x*sum_y) / math.sqrt((number_of_points*sum_xx - sum_x*sum_x)*(number_of_points*sum_yy - sum_y*sum_y))

    print("%6.2lf + %8.4lf*x, r=%6.3lf" % (intercept, slope, correlation))


