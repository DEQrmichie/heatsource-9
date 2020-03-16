# Heat Source, Copyright (C) 2000-2019, 
# Oregon Department of Environmental Quality

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>..

from __future__ import division
import logging
import sys

logger = logging.getLogger(__name__)

class Printer(object):
    """Class to print messages to the console"""
    def __init__(self, msg, progress_bar=False, current=100, total=100):
        if progress_bar is True:
            if int(current / total * 100) < 100:
                sys.stdout.write(self.progress_bar(msg, current, total) + '\r')
            else:
                sys.stdout.write(self.progress_bar(msg, current, total) + '\n')
        else:
            sys.stdout.write(msg + '\n')            
        sys.stdout.flush()
        #time.sleep(0.05)        
                
    def progress_bar(self, msg, current, total):
        
        percent = int(current / total * 100)
        
        prefix = '  {0} {1} {2}'.format(msg, current, total)
        bar_start = ' ['
        bar_end = '] {0}% '.format(percent)
    
        bar_size = 79 - len(prefix + bar_start + bar_end)
        amount = int(current / (total / float(bar_size)))
        remain = bar_size - amount
    
        bar = '.' * amount + ' ' * remain
        return prefix + bar_start + bar + bar_end