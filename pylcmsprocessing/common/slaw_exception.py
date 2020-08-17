
class SlawStepException(Exception):
    def __init__(self,message):
        super(SlawStepException).__init__(message)


class ParametersException(Exception):
    def __init__ (self,arg,type,value):
        message = "Parameters: "+arg+" value:"+value+"should be of type"+type
        super(ParametersException).__init__(message)
