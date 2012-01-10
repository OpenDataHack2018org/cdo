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

        self.operators = self.getOperators()
    
    def __getattr__(self, method_name):
        """
            This is called every time a class method or property 
            is checked and/or called.
            
            In here we'll return a new function to handle what we
            want to do.
        """
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

            call = [self.CDO,kwargs["options"],','.join(operator),' '.join(io)]
            print ' '.join(call)
            print 'before'

            proc = subprocess.Popen(' '.join(call),
                shell  = True,
                stderr = subprocess.PIPE,
                stdout = subprocess.PIPE)
            retvals = proc.communicate()
            print retvals[0]
            print retvals[1]

          

        
        if method_name in self.operators:
            return get.__get__(self)
        else:
            # If the method isn't in our dictionary, act normal.
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


if __name__ == '__main__':
   cdo = Cdo()
   print(cdo.CDO)
   s = cdo.operators
   # Call an arbitrary operator of CDO
   cdo.sinfov(input="-topo",options="-f nc")
   print "#=============================================================="
   cdo.sinfov(input="-remapnn,r36x18 -topo",options="-f nc")
   print "#==========================================================================================="
   f = 'ofile.nc'
   cdo.expr("'z=log(abs(topo+1))*9.81'",input="-topo",output = f,options="-f nc")
   cdo.infov(input=f)
