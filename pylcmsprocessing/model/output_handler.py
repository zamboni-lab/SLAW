import os


class OutputDirectory:
    def __init__(self,directory):
        if not os.path.isabs(directory):
            raise FileNotFoundError("Directory path"+directory+" is not an absoult path p[lease input an absolute path")
        self.dir = directory

    ###Given a relative or absoulte path return an absolute path
    def convertPath(self,path):
        if os.path.isabs(path):
            return os.path.normpath(path)
        else :
            return os.path.normpath(os.path.abspath(os.path.join(self.dir, path)))
    def getRoot(self):
        return(self.dir)

    def getDir(self,path):
        val = self.convertPath(path)
        if not os.path.isdir(val):
            os.makedirs(val)
        return val

    def getFile(self,path):
        val = self.convertPath(path)
        if not os.path.isdir(os.path.dirname(val)):
            os.makedirs(os.path.dirname(val))
        return val

    def exists(self,path):
        path = self.convertPath(path)
        return os.path.exists(path)


