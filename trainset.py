#!/usr/bin/env python
'''
trainset.in contains input training data for reax.

Here is the basic format of these sections.
I am not sure what a Heat of formation in reax is. those appear to be captured in the ENERGY section.

HEATFO
descriptor weight reference_value
ENDHEATFO

CELL PARAMETERS
#descriptor weight type reference value
fcc0 1 a 4.0
fcc0 1 b 4.0
fcc0 1 c 4.0
fcc0 1 alpha 90
fcc0 1 beta 90
fcc0 1 gamma 90
# comment
fccm1 1 a 3.9
fccm1 1 b 3.9
fccm1 1 c 3.9
fccm1 1 alpha 90
fccm1 1 beta 90
fccm1 1 gamma 90
#
fccp1 1 a 4.1
fccp1 1 b 4.1
fccp1 1 c 4.1
fccp1 1 alpha 90
fccp1 1 beta 90
fccp1 1 gamma 90
END CELL PARAMETERS

ENERGY
weight + descriptor1 /a +- descriptor2 /b +- descriptor3 /c reference_value
ENDENERGY
'''
from Cheetah.Template import Template
from Scientific.IO.FortranFormat import FortranFormat, FortranLine
from Scientific.Geometry import Vector
import  math, pickle, os
from reax import *

def parse_trainset(trainset_file='trainset.in'):
    '''
    reads all the lines *between* the section indicators in
    trainset.in and stores them in a cells dictionary.

    this is so I can append new lines to existing sections, and later
    write them back out. Each line is stripped.
    '''

    cells = {}
    cells['ENERGY'] = []
    cells['CELL PARAMETERS'] = []
    cells['HEATFO'] = []
    if not os.path.exists(trainset_file):
        return cells

    EFLAG, HFLAG, CFLAG = False, False, False
    f = open(trainset_file,'r')
    for line in f:
        if line.startswith('ENERGY'):
                EFLAG = True
                continue
        elif line.startswith('ENDENERGY'):
                EFLAG = False

        elif line.startswith('HEATFO'):
                HFLAG = True
                continue
        elif line.startswith('ENDHEATFO'):
                HFLAG = False

        elif line.startswith('CELL PARAMETERS'):
                CFLAG = True
                continue
        elif line.startswith('END CELL PARAMETERS'):
                CFLAG = False

        if EFLAG:
            cells['ENERGY'].append(line.strip())
        elif CFLAG:
            cells['CELL PARAMETERS'].append(line.strip())
        elif HFLAG:
            cells['HEATFO'].append(line.strip())
    f.close()
    return cells

def write_trainset(trainset_file='trainset.in', cells=None):
    '''
    writes all entries contained in the cells dictionary to the
    trainset file.
    '''

    f = open(trainset_file, 'w')
    for key in cells:
        if cells[key] == []:
            continue
        f.write('%s\n' % key)
        for line in cells[key]:
            f.write(line+'\n')
        f.write('END%s\n' % key)
    f.close()

#######################################################
# functions to generate trainset entries

def cell_parameters(atoms, weight=1.0, cells=None):
    '''
    adds entries to CELL PARAMETERS for the unit cell geometry of the
    atoms object.
    '''

    description = reax_hash(atoms)

    # get all the parameters needed for a cell parameters block
    uc = atoms.get_cell()
    a,b,c = [Vector(x) for x in uc]
    A = a.length()
    B = b.length()
    C = c.length()

    rad2deg = 360./(2*math.pi)
    alpha = b.angle(c)*rad2deg
    beta  = a.angle(c)*rad2deg
    gamma = a.angle(b)*rad2deg

    cells['CELL PARAMETERS'].append('# created by reax.trainset')
    # this is an org-link to the geo file
    #cells['CELL PARAMETERS'].append('# file:geo::%s' % description)
    for key,val in [('a',A),
                    ('b',B),
                    ('c',C),
                    ('alpha',alpha),
                    ('beta', beta),
                    ('gamma', gamma)]:
        s = '%(description)s %(weight)f %(key)s %(val)s' % locals()
        if s not in cells['CELL PARAMETERS']:
            cells['CELL PARAMETERS'].append(s)

    return cells

def energy(atoms, reference_energies,
           comment=None,
           weight=1.0, cells=None):
    '''generate lines for formation energies

    atoms is an ase.Atoms object
    reference_energies is a dictionary:
    chemical symbol: (description, energy)

    cells is a dictionary that you are appending these entries too.
    '''
    description = reax_hash(atoms)
    # add org-link
    #cells['ENERGY'].append('# file:geo::%s ' % description)

    e = atoms.get_potential_energy()
    EE = e

    c = {}
    for atom in atoms:
        if atom.symbol in c:
            c[atom.symbol] += 1
        else:
            c[atom.symbol] = 1

    x = {}
    for key in c:
        x[key] = c[key]/float(len(atoms)) # mole fractions

    # compute formation energy
    for atom in atoms:
        e -= reference_energies[atom.symbol][1]

    # construct the string for reax to compute formation energy
    s = '%1.2f + %s /1' % (weight, description)
    for entry in reference_energies:
        desc, ref_energy = reference_energies[entry]
        xi = x.get(entry,None)
        if xi:
            s += ' - %s / %f ' % (desc, 1./(float(xi)*len(atoms)))
            EE -= ref_energy / (1./(float(xi)*len(atoms)))
    s += '  %f' % (e*23.061)

    cells['ENERGY'].append(s)
    return cells

def eos(atoms1, atoms2, weight=1.0, cells=None):
    '''
    add equation of state data to the trainset.in
    '''
    d1 = reax_hash(atoms1)
    d2 = reax_hash(atoms2)
    e = atoms1.get_potential_energy() - atoms2.get_potential_energy()
    s = '%1.2f + %s /1 - %s /1 %f' % (weight,
                                      d1, d2,
                                      e*23.061)
    cells['ENERGY'].append(s)
    return cells

if __name__ in ['__main__']:
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
