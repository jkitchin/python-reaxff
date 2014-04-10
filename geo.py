#!/usr/bin/env python
'''
Converts an atoms object into an entry in the reaxff geo file.

This file is a module and script.

commandline usage: geo.py -f geo ncfile

python usage:
from reaxx.geo import write_bgf
write_bgf(atoms)
'''
from reax import *
from Cheetah.Template import Template
from Scientific.IO.FortranFormat import FortranFormat, FortranLine
from Scientific.Geometry import Vector
import math, os, pickle

bgf_template = '''
# file:$directory
XTLGRF 200
DESCRP $description
REMARK created by reax.geo
$crstx
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
#for $s in $atom_lines
$s
#end for
FORMAT CONECT (a6,12i6)
UNIT ENERGY   kcal
ENERGY        $energy
END
'''
def reax_coordinates(atoms):
	'''
	returns an array of reax coordinates
		'''
	import numpy as np
	pi=3.14159265
	rdndgr=180.0/pi
	dgrrdn=1.0/rdndgr

	uc = atoms.get_cell()
	vscale = 1.0 # I think since i get the uc from the atoms,
				 # i dont' need vscale like adri did

	tm11, tm21, tm31 = uc[0,:]
	tm12, tm22, tm32 = uc[1,:]
	tm13, tm23, tm33 = uc[2,:]

	fc = atoms.get_scaled_positions()

	## get a,b,c, alpha, beta, gamma
	aaxis = np.sqrt(np.dot(uc[0,:],uc[0,:]))
	baxis = np.sqrt(np.dot(uc[1,:],uc[1,:]))
	caxis = np.sqrt(np.dot(uc[2,:],uc[2,:]))

	adotb = np.dot(uc[0,:],uc[1,:])
	adotc = np.dot(uc[0,:],uc[2,:])
	bdotc = np.dot(uc[1,:],uc[2,:])

	hgamma=rdndgr*np.arccos(adotb/(aaxis*baxis))
	hbeta=rdndgr*np.arccos(adotc/(aaxis*caxis))
	halfa=rdndgr*np.arccos(bdotc/(baxis*caxis))

	halfa2=halfa*dgrrdn                          #Convert to ReaxFF cell matrix
	hbeta2=hbeta*dgrrdn
	hgamma2=hgamma*dgrrdn
	sinalf=np.sin(halfa2)
	cosalf=np.cos(halfa2)
	sinbet=np.sin(hbeta2)
	cosbet=np.cos(hbeta2)
	cosphi=(np.cos(hgamma2)-cosalf*cosbet)/(sinalf*sinbet)
	if (cosphi > 1.0):
		cosphi=1.0

	sinphi=np.sqrt(1.0-cosphi*cosphi)
	tm11r=vscale*aaxis*sinbet*sinphi
	tm21r=vscale*aaxis*sinbet*cosphi
	tm31r=vscale*aaxis*cosbet
	tm22r=vscale*baxis*sinalf
	tm32r=vscale*baxis*cosalf
	tm33r=vscale*caxis

	nat = len(atoms)

	c = []
	for i in range(nat):
		ct1 = fc[i,0]*tm11r
		ct2 = fc[i,0]*tm21r + fc[i,1]*tm22r
		ct3 = fc[i,0]*tm31r + fc[i,1]*tm32r + fc[i,2]*tm33r
		c.append([ct1, ct2, ct3])

	return c

def write_bgf(atoms, geofile=None):
	'''
	prepares a bgf entry from an atoms object and writes it to the geofile
	'''
	if geofile is None:
		geofile = 'geo'

	# create the geofile if it does not already exist
	if not os.path.exists(geofile):
		f = open(geofile,'w')
		f.close()

	# get descriptions to make sure they are all unique.
	f = open(geofile,'r')
	DESCRIPTIONS = []
	for line in f.readlines():
		if line.startswith('DESCRP'):
			DESCRIPTIONS.append(line[6:].strip())
	f.close()

	# we need a hash of something to check for uniqueness of the atoms.
	description = reax_hash(atoms)

	if description in DESCRIPTIONS:
		raise Exception, 'Non-unique description "%s" in geofile' % description

	# now prepare the unit cell information for the CRYSTX line
	uc = atoms.get_cell()
	a,b,c = [Vector(x) for x in uc]
	A = a.length()
	B = b.length()
	C = c.length()

	rad2deg = 360./(2*math.pi)
	alpha = b.angle(c)*rad2deg
	beta  = a.angle(c)*rad2deg
	gamma = a.angle(b)*rad2deg

	crstx = str(FortranLine(('CRYSTX',
							 A,
							 B,
							 C,
							 alpha,
							 beta,
							 gamma),
							 FortranFormat('a6,1x,6f11.5')))

	# now we prepare each atom line. These are printed in the template.
	atom_lines = []
	reax_coords = reax_coordinates(atoms)
	fmt = FortranFormat('a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5')
	for i,atom in enumerate(atoms):
		x,y,z = reax_coords[i]
		line = FortranLine(('HETATM', i, atom.symbol,
							'', '', '',
							x,y,z,
							atom.symbol,
							1, 1, 0.0),
							fmt)
		atom_lines.append(line)

	# energy in kcal
	calc = atoms.get_calculator()
	if calc is not None:
		directory = os.path.abspath(calc.vaspdir)
		energy = atoms.get_potential_energy()*23.061
	else:
		directory = 'None'
		energy = 'None'


	bgf = Template(bgf_template, searchList=[locals()])

	entry = bgf.respond()

	if geofile is None:
		geofile = 'geo'

	f = open(geofile, 'a')
	f.write(entry)
	f.close()

def newer_files_exist(reference_file='geo'):
    import glob, os
    if not os.path.exists(reference_file):
        # no reference file, so there are newer files
        return True

    geotime = os.path.getmtime(reference_file)
    dirs = glob.glob('[0-9]*')
    for d in dirs:
        energy = os.path.join(d,'energy')
        # first check energy files
        if os.path.exists(energy):
            if os.path.getmtime(energy) > geotime:
                return True

        # now check EOS files
        if os.path.exists(os.path.join(d,'eos-exp')):

            eos_d = os.path.join(d,'eos-exp')
            for d in os.listdir(eos_d):
                nextdir = os.path.join(eos_d,d)
                if (os.path.isdir(nextdir)):
                    energy = os.path.join(nextdir,'energy')
                    if os.path.exists(energy):
                        if os.path.getmtime(energy) > geotime:
                            return True
    # getting here means no newer files found.
    return False

