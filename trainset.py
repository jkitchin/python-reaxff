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

