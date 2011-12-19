method_dictionary = ["sinfov"]

class Cdo(object):
    def __init__(self):
        """
            Store params and junk. Ordinarily more verbose.
        """
    
    def __getattr__(self, method_name):
        """
            This is called every time a class method or property 
            is checked and/or called.
            
            In here we'll return a new function to handle what we
            want to do.
        """
        def get(self, *args):
            # Make our API calls, return data, etc
            for arg in args:
              print arg

            return
        
        if method_name in method_dictionary:
            return get.__get__(self)
        else:
            # If the method isn't in our dictionary, act normal.
            raise AttributeError, method_name

# Instantiate...
cdo = Cdo()

# Call an arbitrary method.
cdo.sinfov(1,'str')
