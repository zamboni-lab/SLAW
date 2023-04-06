import time
import timeit
import os
import logging
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
        tval = (self.measures[pindex]-self.measures[0])/60
        last_interval = (self.measures[pindex]-self.measures[pindex-1])/60
        logging.info("STEP: " + name +"   DURATION: " + "%0.1f" % last_interval +" min   TOTAL:"+"%0.1f" % tval +" min")


if __name__=="__main__":
    trt = Timer()
    trt.store_point("Init")
    ##Do some stuff
    time.sleep(2)
    trt.store_point("after_wait")
    print(trt)
    trt.print_point("after_wait")