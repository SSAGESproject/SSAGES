# Generates patch for gromacs 
import subprocess

gmx_orig = '/home/hsidky/Code/gromacs-original/gromacs-5.1.3'
gmx_mod = '/home/hsidky/Code/gromacs-5.1.3'

output = subprocess.Popen(['diff', '-ruN', '{0}/src/'.format(gmx_orig),'{0}/src/'.format(gmx_mod)], stdout=subprocess.PIPE).communicate()[0]
output = str.replace(output, gmx_orig, '/gromacs-original')
output = str.replace(output, gmx_mod, '/gromacs-ssages')
print(output)