'''
when run as a script, this module creates a trainset.in file from
directories created by ATAT. It is assumed for now that vasp was
used to run the calculations.

The script only works up to a ternary system.
'''
import sys
from optparse import OptionParser

parser = OptionParser(usage='trainset.py [-f trainsetfile]')
parser.add_option('-f',
                  nargs = 1,
                  help='Specify the trainset.in to create. default = trainset.in ')

options, ncfiles = parser.parse_args()

if options.f is None:
    trainset_file = 'trainset.in'

from jasp import *
import glob, os
import numpy as np
from geo import newer_files_exist

if newer_files_exist(trainset_file):
    if os.path.exists(trainset_file):
        os.unlink(trainset_file)
else:
    print 'No newer files found. trainset.in is up to date.'
    sys.exit()

cells = parse_trainset()

dirs = glob.glob('[0-9]*')

if os.path.exists('atoms.out'):
    mmaps = True
else:
    mmaps = False

# reference endpoints
with jasp('0') as c0:
    atoms0 = c0.get_atoms()
    e0 = atoms0.get_potential_energy()/len(atoms0)
    d0 = reax_hash(atoms0)

with jasp('1') as c1:
    atoms1 = c1.get_atoms()
    e1 = atoms1.get_potential_energy()/len(atoms1)
    d1 = reax_hash(atoms1)

if mmaps:
    with jasp('2') as c2:
         atoms2 = c2.get_atoms()
         e2 = atoms2.get_potential_energy()/len(atoms2)
         d2 = reax_hash(atoms2)

#reference energies in eV/atom
reference_energies = {}
reference_energies[atoms0[0].symbol] = (d0, e0)  # (description, energy)
reference_energies[atoms1[0].symbol] = (d1, e1)
if mmaps:
    reference_energies[atoms2[0].symbol] = (d2, e2)

'''
this is how Adri likes the trainset.in
ENERGY
# Volume [phase_name1]
[Equation of state data]
# Volume [phase_name2]
[Equation of state data]
# Heats_of_formation
[All heat of formation data[
ENDENERGY
'''
# now we add formation energies and eos
counter = 0
eos_counter = 0
for d in dirs:
    if (os.path.exists(os.path.join(d,'energy'))
        and not os.path.exists(os.path.join(d,'error'))):
        # geometry of cell
        with jasp(d) as calc:
            atoms = calc.get_atoms()
        cells = cell_parameters(atoms, weight=1.0, cells=cells)

        # equation of state Adri suggested we put a higher weight on
        # points close to the minimum
        eosdir = os.path.join(d,'eos-exp')
        if os.path.isdir(eosdir):
            edirs = glob.glob('%s/*' % eosdir)
            cells['ENERGY'].append('# Volume %s' % reax_hash(atoms))
            for ed in edirs:
                if (os.path.exists(os.path.join(ed,'energy'))
                    and not os.path.exists(os.path.join(ed,'error'))):
                    with jasp(ed) as calc2:
                        atoms2 = calc2.get_atoms()
                        e1 = atoms.get_potential_energy()
                        e2 = atoms2.get_potential_energy()
                        weight = 0.1 + 0.25*np.abs((e2 - e1))*23.061 #kcal

                        cells = eos(atoms,
                                    atoms2,
                                    weight=weight,
                                    cells=cells)
                        eos_counter += 1

# now heats of formation go last
cells['ENERGY'].append('# Heats_of_formation')
for d in dirs:
    if (os.path.exists(os.path.join(d,'energy'))
        and not os.path.exists(os.path.join(d,'error'))):
        # geometry of cell
        with jasp(d) as calc:
            atoms = calc.get_atoms()

        # This is a heat of formation
        cells = energy(atoms,
                       reference_energies,
                       weight=1.0,
                       cells=cells)
        counter += 1

write_trainset(trainset_file, cells=cells)
print 'wrote {0:d} heats of formation'.format(counter)
print 'wrote {0:d} EOS points'.format(eos_counter)
