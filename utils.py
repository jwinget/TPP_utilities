#-- Utilities
#
# Some general-purpose tools
# Author: Jason Winget
# Version: 0.1
#
#--

import time

def Timer(fn):
	''' Decorator to time a function's execution '''
	def timed(*args, **kw):
		start = time.time()
		result = fn(*args, **kw)
		elapsed = str(round(time.time() - start, 2))
		print fn.__name__+' completed in '+elapsed+' s'
		return result
	return timed
