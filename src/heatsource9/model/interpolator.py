"""Interpolation utilities used by model setup and during model runs."""


from bisect import bisect_left, bisect_right
from collections import defaultdict
from time import gmtime


class Interpolator(defaultdict):
    """
    Linearly interpolated dictionary class for numeric key/value pairs.
    This class will store a timeseries and use linear interpolation 
    between all the values for each model time step filling the dictionary 
    with the results.
    """

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.sortedkeys = None

    def __missing__(self, key):
        """
        Fill missing input timestamps using linearly interpolation. 
        Interploates between the surrounding timestamps
        """
        # First time the dictionary is in actual use, so we do some setup
        if not self.sortedkeys:
            if not len(self.keys()):
                self.__missing__ = lambda x: (0.0,)
                return (0.0,)
            self.sortedkeys = sorted(self.keys())
        
         # Find the previous and next available keys
        ind = bisect_right(self.sortedkeys, key) - 1
        x0 = int(self.sortedkeys[ind])
        x1 = int(self.sortedkeys[ind + 1])

        y0 = self[x0]
        y1 = self[x1]

        if isinstance(y0, tuple):
            if not len(y0):
                # We have nothing in the tuple, 
                # so return a blank tuple (not 'val', which is None)
                return ()
            vals = []
            for i in range(len(y0)):
                vals.append(y0[i] + ((y1[i] - y0[i]) * (key - x0)) / (x1 - x0))
            return tuple(vals)

        return y0 + ((y1 - y0) * (key - x0)) / (x1 - x0)

    def view(self, minkey, maxkey, fore=None, aft=None):
        """"
        Return dictionary subset

        Return a subset of the current dictionary containing items
        with keys between minkey and maxkey. i.e. compares the UTC 
        date (year, month, day) of minkey/maxkey to the dataset bounds, 
        and returns self unchanged if they match exactly. Otherwise builds 
        a new Interpolator containing keys in the requested range.
        
        If either or both are fore and/or aft are anything but None, then the returned
        dictionary will also contain the next element before or
        after minkey and maxkey, respectively.

        
        """
        keys = sorted(self.keys())
        if not keys:
            return self
        
        # If our subset includes all values (i.e. we're not really
        # subsetting), just return ourself.
        start_new = gmtime(minkey)[0:3]
        start_old = gmtime(min(keys))[0:3]
        stop_new = gmtime(maxkey)[0:3]
        stop_old = gmtime(max(keys))[0:3]
        if (start_new == start_old) and (stop_new == stop_old):
            return self
        
        # Get the minimum and maximum indices, including the one
        # before and one after if fore or aft are anything but None
        newmin = bisect_left(keys, minkey) - (fore is not None)
        newmax = bisect_right(keys, maxkey) + (aft is not None)

        d = Interpolator()
        for k in keys[newmin:newmax]:
            d[k] = self[k]
        return d
