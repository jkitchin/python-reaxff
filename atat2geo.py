from geo import *
'''
when run as a script, this module creates a geo file from
directories created by ATAT. It is assumed for now that vasp was
used to run the calculations.
'''
from optparse import OptionParser
import sys

parser = OptionParser(usage='geo.py [-f geofile]')

parser.add_option('-f',
                  nargs = 1,
                  help='Specify the geofile to add entry to. default = geo')

options, args = parser.parse_args()

if options.f is None:
    geofile = 'geo'

from jasp import *
import glob, os

# first check if we are up to date
if os.path.exists('geo'):
    geotime = os.path.getmtime('geo')
else:
    geotime = None

dirs = glob.glob('[0-9]*')

if newer_files_exist():
    if os.path.exists('geo'):
        os.unlink('geo')
else:
    print 'No newer files found. geo is up to date.'
    sys.exit()

count = 0
for d in dirs:
    if (os.path.exists(os.path.join(d,'energy'))
        and not os.path.exists(os.path.join(d,'error'))):
        with jasp(d) as calc:
            atoms = calc.get_atoms()
        write_bgf(atoms,'geo')

        # now we add the equation of state data points
        if os.path.exists(os.path.join(d,'eos-exp')):
            # get directories in eos-exp
            eos_d = os.path.join(d,'eos-exp')
            for d in os.listdir(eos_d):
                if (os.path.isdir(os.path.join(eos_d,d))
                    and os.path.exists(os.path.join(eos_d,d,'energy'))
                    and not os.path.exists(os.path.join(eos_d,d,'error'))):
                    with jasp(os.path.join(eos_d,d)) as calc:
                        atoms = calc.get_atoms()
                    write_bgf(atoms, geofile)
                    count += 1
print 'Wrote {0:d} entries'.format(count)
