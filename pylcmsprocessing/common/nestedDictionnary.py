class nestedDictionnary:
    def __init__(self,dic):
        self.dic = dic

    def recurr_access_dic(self,dic,path):
        if len(path)==1:
            return(dic[path[0]])
        return self.recurr_access_dic(self,dic[path[0]],path[1:])

    def access_nested_dic(self,path,sep="/"):
        if isinstance(path,str):
            path = path.split(sep)
        return self.recurr_access_dic(self,self.dic,path)

    def recurr_set_dic(self,dic,path,val):
        if len(path)==1:
            dic[path[0]]=val
            return
        self.recurr_set_dic(self,dic[path[0]],path[1:],val)

    def set_nested_dic(self,dic,path,val,sep="/"):
        if isinstance(path,str):
            path = path.split(sep)
        self.recurr_set_dic(self,dic,path,val)


