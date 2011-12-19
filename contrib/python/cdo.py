from os import system,popen
import re
import subprocess

method_dictionary = ["sinfov","remap"]

class Cdo(object):

    def getOperators(self):
        proc = subprocess.Popen(['cdo','-h'],stderr = subprocess.PIPE)
        ret  = proc.communicate()
        l    = ret[1].find("Operators:")
        ops  = ret[1][l:-1].split("\n")[1:-1]
        endI = ops.index('')
        s    = ' '.join(ops[:endI]).strip()
        s    = re.sub("\s+" , " ", s)
        return s.split(" ")

    def __init__(self):
        """
            Store params and junk. Ordinarily more verbose.
        """
        self.operators = self.getOperators()
    
    def __getattr__(self, method_name):
        """
            This is called every time a class method or property 
            is checked and/or called.
            
            In here we'll return a new function to handle what we
            want to do.
        """
        def get(self, *args):
            # Make our API calls, return data, etc
            print method_name
            print ' '.join(['cdo', ','.join([method_name, ','.join(args)])])

            return
        
        if method_name in method_dictionary:
            return get.__get__(self)
        else:
            # If the method isn't in our dictionary, act normal.
            raise AttributeError, method_name


if __name__ == '__main__':
    cdo = Cdo()
    s = cdo.operators
    print "#=============================================================="
    for ss in s:
      print ss

    # Call an arbitrary method.
    cdo.sinfov('1','str')
    cdo.remap('gridFile','weightfile')
