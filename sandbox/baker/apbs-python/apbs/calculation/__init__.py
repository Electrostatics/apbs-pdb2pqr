""" The main APBS calculation module """
import logging
import unittest
import queue
import parser
_LOGGER = logging.getLogger(__name__)

class APBSCalculation(queue.Queue):
    """ Stores APBS calculation objects and handles execution """
    # TODO - The use of a queue is intended to make future PDB2PQR execution more practical; however, it creates the problem
    def __init__(self, apbs_input, *args, **kwargs):
        """ The constructor fills the queue full of tasks """
        super(APBSCalculation, self).__init__(*args, **kwargs)
        self.apbs_input = apbs_input
        for read in self.apbs_input.reads:
            self.put(read)
        for calc in self.apbs_input.calcs:
            self.put(calc)
        for print_ in self.apbs_input.prints:
            self.put(print_)
    def process(self):
        """ Process the items in the queue """
        while not self.empty():
            item = self.get()
            print(item)

class _TestSetup(unittest.TestCase):
    """ Test ability to set up APBS calculations """
    def setUp(self):
        from parser import read_parser
        read = read_parser.Read()
        read.mol = [read_parser.Mol(format_="pqr", path="./examples/fas2.pqr")]
        read.validate()
        _LOGGER.debug(read)
        from parser import elec_parser
        elec_solv = elec_parser.Elec()
        elec_solv.name = elec_parser.Name("elec_solv")
        elec_solv.solvtype = elec_parser.Solvtype("mg-auto")
        elec_solv.bcfl = elec_parser.Bcfl("mdh")
        elec_solv.chgm = elec_parser.Chgm("spl0")
        # TODO - with this new input file parser we should be able to give molecules names that we can refer to them by
        elec_solv.cgcent = elec_parser.Cgcent("mol 1")
        elec_solv.cglen = elec_parser.Cglen((120., 120., 120.))
        elec_solv.dime = elec_parser.Dime((97, 97, 97))
        elec_solv.fgcent = elec_parser.Fgcent("mol 1")
        elec_solv.fglen = elec_parser.Fglen((60., 60., 60.))
        elec_solv.ion = [elec_parser.Ion(charge=+1.0, conc=0.1, radius=2.0),
                         elec_parser.Ion(charge=-1.0, conc=0.1, radius=2.0)]
        elec_solv.eqntype = elec_parser.Eqntype("lpbe")
        elec_solv.mol = elec_parser.Mol(1)
        elec_solv.pdie = elec_parser.Pdie(2.0)
        elec_solv.sdens = elec_parser.Sdens(10.0)
        elec_solv.sdie = elec_parser.Sdie(78.54)
        elec_solv.srad = elec_parser.Srad(1.4)
        elec_solv.srfm = elec_parser.Srfm("smol")
        elec_solv.swin = elec_parser.Swin(0.2)
        elec_solv.temp = elec_parser.Temp(298.15)
        elec_solv.calcenergy = elec_parser.Calcenergy()
        import copy
        elec_vac = copy.deepcopy(elec_solv)
        elec_vac.ion = []
        elec_vac.sdie.set_value(1.0)
        from parser import print_parser
        print_ = print_parser.Print(opstring="elecEnergy 1 - 2")
        self.apbs_input = parser.CalcInput()
        self.apbs_input.reads = [read]
        self.apbs_input.calcs = [elec_solv, elec_vac]
        self.apbs_input.prints = [print_]
    def test_load_queue(self):
        """ Test loading the calculation queue """
        _LOGGER.info("Testing APBS calculation queue loading...")
        calc_queue = APBSCalculation(self.apbs_input)
        while True:
            try:
                item = calc_queue.get(block=False)
            except queue.Empty:
                break
            _LOGGER.info("The queue contains a %s" % type(item))
        _LOGGER.info("I hope that's what you expected because that is all this unit test does.")
