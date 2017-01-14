#! /usr/bin/python

from sys import *
from datetime import *

dt = datetime.now()
print '%02d%02d%02d' % ( dt.hour, dt.minute, dt.second )

