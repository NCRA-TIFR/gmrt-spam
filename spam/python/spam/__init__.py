###############################################################################

import atexit

###############################################################################

# import package modules
from files import *
from aips import *
from plot import *
from acalc import *
from sphere import *
from parameter import *
from skymodel import *
from image import *
from peel import *
from mpfit import *
from ionosphere import *
from solutions import *
from extraction import *
from instrumental import *
from pointing import *
from error import *
from unwrap import *
from instrumental import *
from flag import *
from calibrate import *
from gmrt import *
from pipeline import *
from tgss import *

###############################################################################

# initialize matplotlib
#matplotlib.rc( 'text', usetex = True )
matplotlib.rc( 'text', usetex = False )

###############################################################################

#aips_id = allocate_aips_id()
#print 'Allocated AIPS ID = %d (%s)' % ( aips_id, decimal_to_ehex( aips_id ) )
#atexit.register( free_aips_id )

###############################################################################

# initialize numpy
import numpy
numpy.seterr( all = 'raise' )

###############################################################################

