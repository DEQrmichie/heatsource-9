# Heat Source, Copyright (C) 2000-2015, Oregon Department of Environmental Quality

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import time

class LoggerDiety(object):
    def __init__(self):
        self._last = None
        self._file = None
    def __del__(self):
        try:
            self._file.close()
        except AttributeError:
            pass
    def SetFile(self, filename):
        try:
            self._file = open(filename,"w")
        except IOError:
            raise IOError("Opening output directory failed. Make sure directory exists before running.")
    def __call__(self, message, n=None,t=None): self.write(message)
    def write(self, message):
        if message != self._last:
            t = time.strftime("%H:%M:%S",time.gmtime(time.time()))
            self._file.write("%s-> %s\n"% (t,message))
            self._last = message
    def progress(self): self._file.write(".")
Logger = LoggerDiety()
