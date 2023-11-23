import numpy as np
import sys
import os
from pathlib import Path

start_t = 0.0
end_t = float(sys.argv[1])
dt = 0.01
B = float(sys.argv[2])
J = float(sys.argv[3])
De = float(sys.argv[4])

cwd = os.getcwd()
fname = cwd+'/restore_ke'
my_file = Path(fname)
if not my_file.is_file():
  print(fname)
  cmd = 'qcc -O2 -Wall restore_ke.c -o restore_ke -lm'
  os.system(cmd)
#  input("Press Enter to continue")

for i in np.arange(start_t,end_t,dt):
    nameIn = f'intermediate_f/snapshot-{i:.4f}'
    cmd = f"./restore_ke %s {i:.4f} %s %s %s" % (nameIn, str(B), str(J), str(De))
    os.system(cmd)
    print('File: '+nameIn+' t = '+str(i))
