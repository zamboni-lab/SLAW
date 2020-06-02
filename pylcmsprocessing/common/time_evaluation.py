import time
import timeit
import os
class Timer:
    def __init__(self):
        self.names= []
        self.measures = []

    def store_point(self,name):
        self.measures.append(timeit.default_timer())
        self.names.append(name)

    def __str__(self):
        to_join = [self.names[idx]+" "+str(self.measures[idx]-self.measures[0]) for idx in range(len(self.names))]
        return "\n".join(to_join)

    def print_point(self,name):
        pindex = self.names.index(name)
        tval = self.measures[pindex]-self.measures[0]
        last_interval = self.measures[pindex]-self.measures[pindex-1]
        print("STEP:",name,"TOTAL_TIME:","%0.2f" % tval+"s"," LAST_STEP:","%0.2f" % last_interval+"s")


if __name__=="__main__":
    trt = Timer()
    trt.store_point("Init")
    ##Do some stuff
    time.sleep(2)
    trt.store_point("after_wait")
    print(trt)
    trt.print_point("after_wait")