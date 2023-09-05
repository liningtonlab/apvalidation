# Import specific classes from the module files
from .varian import Varian
from .bruker import Bruker
from .jeol import JEOL
from .jcampdx import Jcampdx


# You can also define what gets imported when using 'from apvalidation.extract import *'
__all__ = ["Varian", "Bruker", "Jcampdx", "JEOL"]