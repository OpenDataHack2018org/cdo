import os,re,subprocess

# Copyright (C) 2011-2012 Ralf Mueller, ralf.mueller@zmaw.de
# See COPYING file for copying and redistribution conditions.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

class Cdo(object):

    def __init__(self):
        if os.environ.has_key('CDO'):
          self.CDO = os.environ['CDO']
        else:
          self.CDO = 'cdo'

        self.operators   = self.getOperators()
        self.returnArray = False

        self.debug = False

    def __getattr__(self, method_name):
        def get(self, *args,**kwargs):
            operator = [method_name]
            if args.__len__() != 0:
              for arg in args:
                operator.append(arg)
              print "args:"
              print args

            io = []
            if kwargs.__contains__("input"):
              io.append(kwargs["input"])

            if kwargs.__contains__("output"):
              io.append(kwargs["output"])

            if not kwargs.__contains__("options"):
              kwargs["options"] = ""

            if not kwargs.__contains__("returnArray"):
              kwargs["returnArray"] = False

            call = [self.CDO,kwargs["options"],','.join(operator),' '.join(io)]

            if self.debug:
              print ' '.join(call)

            proc = subprocess.Popen(' '.join(call),
                                    shell  = True,
                                    stderr = subprocess.PIPE,
                                    stdout = subprocess.PIPE)
            retvals = proc.communicate()

            if self.debug:
              print retvals[0]
              print retvals[1]

            if re.search('(info|show|griddes)',method_name):
              return retvals[1].split('\n')
            else:
              if self.returnArray or kwargs["returnArray"]:
                if not self.returnArray:
                  self.loadCdf()

                return self.cdf(kwargs["output"])
              else:
                return kwargs["output"]

          

        
        if ((method_name in self.__dict__) or (method_name in self.operators)):
          if self.debug:
            print("Found method:" + method_name)

          return get.__get__(self)
        else:
          # If the method isn't in our dictionary, act normal.
          print("#=====================================================")
          print("Cannot find method:" + method_name)
          raise AttributeError, method_name

    def getOperators(self):
        proc = subprocess.Popen([self.CDO,'-h'],stderr = subprocess.PIPE,stdout = subprocess.PIPE)
        ret  = proc.communicate()
        l    = ret[1].find("Operators:")
        ops  = ret[1][l:-1].split("\n")[1:-1]
        endI = ops.index('')
        s    = ' '.join(ops[:endI]).strip()
        s    = re.sub("\s+" , " ", s)
        return s.split(" ")

    def loadCdf(self):
      try:
        import pycdf as cdf
        self.returnArray = True
        self.cdf         = cdf.CDF
      except ImportError:
        raise ImportError,"Module pycdf is required to return numpy arrays."

    def setReturnArray(self,value=True):
      self.returnArray = value
      if value:
        self.loadCdf()


    def unsetReturnArray(self):
      self.setReturnArray(False)
