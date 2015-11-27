""" The main APBS calculation module """
import logging
import unittest
import queue
import parser
_LOGGER = logging.getLogger("calculation")

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
    INPUT_PATH = "./examples/uber-input.in"
    def setUp(self):
        with open(INPUT_PATH, "rt") as input_file:
            self.apbs_from_text = TextDecoder().decode(input_file.read())
