

from msmbuilder import Trajectory
from msmbuilder import arglib
import numpy as np
from montecarlo import MonteCarlo

parser = arglib.ArgumentParser(get_metric=True)

parser.add_argument('generators', help='Generators filename')
