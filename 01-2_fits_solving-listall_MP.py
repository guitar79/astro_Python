# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

No module named 'ccdproc'
conda install -c conda-forge ccdproc
"""
import os
import shutil 
import Python_utilities
import astro_utilities

#########################################
from multiprocessing import Process, Queue

class Multiprocessor():
    def __init__(self):
        self.processes = []
        self.queue = Queue()

    @staticmethod
    def _wrapper(func, queue, args, kwargs):
        ret = func(*args, **kwargs)
        queue.put(ret)

    def restart(self):
        self.processes = []
        self.queue = Queue()

    def run(self, func, *args, **kwargs):
        args2 = [func, self.queue, args, kwargs]
        p = Process(target=self._wrapper, args=args2)
        self.processes.append(p)
        p.start()

    def wait(self):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        for p in self.processes:
            p.join()
        return rets

#########################################

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

base_dir = "../CCD_new_files/"

fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))


#########################################
myMP = Multiprocessor()
num_cpu = 6
values = []
num_batches = len(fullnames) // num_cpu + 1

for batch in range(num_batches):
    myMP.restart()
    for fullname in fullnames[batch*num_batches:(batch+1)*num_batches]:
        myMP.run(astro_utilities.KevinSolver, fullname)

    print("Batch " + str(batch))
    myMP.wait()
    values.append(myMP.wait())
    print("OK batch" + str(batch))


#%%
#############################################################################
#Check existence tmp file and rename ...
#############################################################################
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

fullnames_tmp = [w for w in fullnames if ".tmp" in w]

#%%
n = 0
for fullname in fullnames_tmp[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames_tmp))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    try:
        shutil.move(r"{}".format(fullname), \
                        r"{}.fit".format(fullname[:-4]))

    except Exception as err:
        Python_utilities.write_log(err_log_file,
                    '{2} ::: {0} There is no {1} '.format(err, fullname, datetime.now()))      
