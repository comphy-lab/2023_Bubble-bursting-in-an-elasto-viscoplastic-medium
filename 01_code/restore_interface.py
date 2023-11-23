import numpy as np
import sys
import os
from pathlib import Path

start_t = float(sys.argv[1])
end_t = float(sys.argv[2])
B = float(sys.argv[3])
J = float(sys.argv[4])
De = float(sys.argv[5])

dt = 0.01

cwd = os.getcwd()
fname = cwd+'/restore_interface'
my_file = Path(fname)
if not my_file.is_file():
#  print(fname)
  cmd = 'qcc -O2 -Wall restore_interface.c -o restore_interface -lm'
  os.system(cmd)
#  input("Press Enter to continue")

for i in np.arange(start_t,end_t,dt):
    nameIn = f'intermediate_f/snapshot-{i:.4f}'
    nameOut = f'intermediate_f/interface-{i:.4f}.dat'
    cmd = "./restore_interface %s %s %s %s %s" % (nameIn,nameOut,str(B),str(J),str(De))
    os.system(cmd)
    print('File: '+nameIn+' t = '+str(i))
