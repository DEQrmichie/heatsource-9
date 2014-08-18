""" The IniParams dictionary is simply a regular dictionary that can be 
accessed globally. We keep it here mostly because of legacy architecture from 
heatsource 8. Heatsource 8 used to optimize many of the classes with pysco 
and the primary switchses for that used to be located here. Psyco is no longer 
under development and we have moved on to python 2.7 so everything pysco 
related has been removed. The C module switch used to be here also but now
we run all the routines exclusivly in python. Maybe some day we'll update 
the C code. Until then, we'll keep this simple module in place in case we need 
it for the future. Of course, one important thing to remember is to import 
this module or everything breaks.

Good thing that very important caveat is buried deeply in these
notes where no-one will ever read it."""

# Some needed updates to the C module. 
#1. Number of veg zones hard coded as 4
#2. Output from C module does not include solar blocked and causes a crash.

# We need someting to initialize the dictionary so here it is.
IniParams = {"run_in_python": True,}