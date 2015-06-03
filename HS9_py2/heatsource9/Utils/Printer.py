
from __future__ import division
import sys

class Printer(object):
    def __init__(self, msg, current=0, total=100):
 
        sys.stdout.write(self.progress_bar(msg, current, total) + '\r')
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
    
        bar = '>' * amount + ' ' * remain
        return prefix + bar_start + bar + bar_end