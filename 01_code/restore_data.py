import numpy as np
import sys
import os
from pathlib import Path

start_t = 0.0
end_t = float(sys.argv[1])
print('End time: '+str(end_t))
dt = 0.05
B = float(sys.argv[2])
J = float(sys.argv[3])
De = float(sys.argv[4])

cwd = os.getcwd()
fname = cwd+'/restore_data'
my_file = Path(fname)
if not my_file.is_file():
  print(fname)
  cmd = 'qcc -O2 -Wall restore_data.c -o restore_data -lm'
  os.system(cmd)
  input("Press Enter to continue")

for i in np.arange(start_t,end_t,dt):
    nameIn = f'intermediate_f/snapshot-{i:.4f}'
    nameOut = f'intermediate_f/data-{i:.4f}'
    cmd = "./restore_data %s %s %s %s %s" % (nameIn,nameOut,str(B),str(J),str(De))
    print(cmd)
    os.system(cmd)
    print('File: '+nameIn+' t = '+str(i))
