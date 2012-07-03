from Cheetah.Template import Template
from Scientific.IO.FortranFormat import FortranFormat, FortranLine
from Scientific.Geometry import Vector
import math

bgf_template = '''
XTLGRF 200                                                                      
DESCRP $description                                                                
REMARK                                                                          
RUTYPE CELL OPT     0                                                           
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
        
def write_bgf(atoms,description=None):

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

    atom_lines = []
    fmt = FortranFormat('a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5')
    for i,atom in enumerate(atoms):
        
        line = str(FortranLine(('HETATM',
                               i,
                               atom.symbol,
                               '',
                               '',
                               '',
                               atom.x,
                               atom.y,
                               atom.z,
                               atom.symbol,
                               1,
                               1,
                               0.0),
                               fmt))
        atom_lines.append(line)

    # energy in kcal
    energy = atoms.get_potential_energy()*23.061


    bgf = Template(bgf_template, searchList=[locals()])

    return bgf.respond()

trainset_template = '''
ENERGY
# heats of formation
#for ($weight, $label1, $label2, $m, $label3, $n, $dft_energy) in $lines
$weight + $label1 /1 - $label2 /$m - $label3 /$n $dft_energy
#end for
ENDENERGY
'''

def write_trainset(list_of_atoms, weight=1.0):
    '''
    create the trainset.in function based on heats of formations. This
    file is limited in format, e.g. one cannot have arbitrary algebra
    in it, only + and - and division.i will have to be clever about
    computing formation energies this way. e.g assume a1 is a Ru3Pt
    alloy, r1 is a 4 atom Ru structure, and p1 is a 4 atom Pt
    structure. then this is the line that goes in trainset

    weight + a1 /1 - r1 /2 - r1/4 -p1 /4   dft_energy

    furthermore than cannot be more than 5

    TODO: this function prepares the heat of formation part of the
    trainset. so, it should take a list of atoms objects.
    
    '''
    from collections import defaultdict

    lines = []
    for atoms in list_of_atoms:
        # get composition
        syms = atoms.get_chemical_symbols()
        d = defaultdict(int)
        for s in syms:
            d[s] += 1

        items = d.items()

        natoms = len(atoms)

        calc = atoms.get_calculator()

        # TODO these need to be connected to the labels in geo
        label1 = calc.get_nc()
        label2 = items[0][0] + str(items[0][1])
        label3 = items[1][0] + str(items[1][1])

        m = float(items[0][1])/natoms
        n = float(items[1][1])/natoms

        # TODO this should be a formation energy
        dft_energy = atoms.get_potential_energy()

        lines.append((weight, label1, label2, m, label3, n, dft_energy)) 

    tmpl = Template(trainset_template, searchList=[locals()])
    return tmpl.respond()
    
if __name__ in['__main__']:
    from ase.calculators.jacapo import *
    calc = Jacapo('fe-al/36/out.nc')
    atoms = calc.get_atoms()
    #print write_bgf(atoms)
    print write_trainset([atoms, atoms])
