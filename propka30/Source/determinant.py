
class Determinant:
    """
        Determinant class - set up for later structurization
    """

    def __init__(self, label, value):
        """
        Contructer of determinant object - simple, but helps in creating structure!
        """
        self.label  = label
        self.value  = value


    def add(self, value):
        """
        adding a value to determinant
        """
        self.value += value

