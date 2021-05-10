
class SlawStepException(Exception):
    def __init__(self,message):
        super().__init__(message)

class ParametersException(Exception):
    def __init__ (self,arg,type,value):
        message = "Parameters: "+arg+" value:"+value+"should be of type"+type
        self.message = message


class OptimizationException(Exception):
    def __init__(self,message):
        self.message = message