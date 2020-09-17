# -*- coding: utf-8 -*-


"""Molecule geometry optimization with MOPAC."""


import os

from openbabel import pybel

from pychem import vector3d


class Atom:
    """Wrapper for atomic properties."""

    def __init__(self, Coordinates: List[float]):
        """Initialize an Atom object."""
        self.pos = vector3d.Vector3D()
        self.radius = 0.0
        self.Coordinates = Coordinates
        self.Element = ''

    def SetCoordinates(self):
        """Parse raw ARC coordinates."""
        temp = self.Coordinates
        self.pos.x = float(temp[1])
        self.pos.y = float(temp[2])
        self.pos.z = float(temp[3])

    def GetCoordinates(self):
        """Get coordinates of the atom."""
        self.SetCoordinates()
        return self.pos

    def SetElement(self):
        """Set element from raw ARC coordinates."""
        temp = self.Coordinates
        self.Element = temp[0]

    def GetElement(self):
        """Get element."""
        self.SetElement()
        return self.Element

    def SetRadius(self):
        """Set radius."""
        radii = {'H': 1.20, 'N': 1.55, 'Na': 2.27, 'Cu': 1.40, 'Cl': 1.75, 'C': 1.70,
                 'O': 1.52, 'I': 1.98, 'P': 1.80, 'B': 1.85, 'Br': 1.85, 'S': 1.80, 'Se': 1.90,
                 'F': 1.47, 'Fe': 1.80, 'K': 2.75, 'Mn': 1.73, 'Mg': 1.73, 'Zn': 1.39, 'Hg': 1.8,
                 'Li': 1.8, '.': 1.8}
        temp = self.GetElement()
        if temp in radii.keys():
            self.radius = radii[temp]
        else:
            self.radius = radii['.']

    def GetRadius(self):
        """Get the radius of the atom."""
        self.SetRadius()
        return self.radius


def GetAtomClassList(Coordinates: List[float]) -> List[Atom]:
    """Get a list of atoms from a list of raw ARC coordinates.

    :param Coordinates: raw ARC coordinates as returned by _ReadCoordinates
    """
    Atoms = []
    for i in Coordinates:
        atom = Atom(i)
        atom.SetCoordinates()
        atom.SetElement()
        atom.SetRadius()
        Atoms.append(atom)
    return Atoms


def _ReadCoordinates(arc_file):
    """Read coordinates and charges of atoms from a MOPAC ARC file.

    :param arc_file: Path to MOPAC .arc file
    """
    res = []
    with open(arc_file, 'r') as f:
        templine = f.readlines()
    for line in range(len(templine)):
        if templine[line][-7: -1] == "CHARGE":
            k = line
            break
    for i in templine[k + 4: len(templine) - 1]:
        temp = i.split()
        ElementCoordinate = [temp[0].strip(), temp[1].strip(),
                             temp[3].strip(), temp[5].strip(),
                             temp[-1].strip()]
        res.append(ElementCoordinate)
    return res


def FormatConversion(inputmol: Chem.Mol) -> Tuple[str, str]:
    """Prepare a molecule to be optimized by MOPAC."""
    #inputmol.removeh()
    inputmol.addh()
    inputmol.make3D(forcefield='mmff94', steps=50)    ##Gemetrical optimization
    ##forcefields = ['uff', 'mmff94', 'ghemical']
    #make3D(self, forcefield = "mmff94", steps = 50)
    ##inputmol.localopt(forcefield='mmff94',steps=50)
    outputmol = pybel.Outputfile('mop', "temp.dat", overwrite=True)
    outputmol.write(inputmol)
    outputmol.close()
    f = file('temp.dat','r+')
    f.write('AM1              ')
    f.close()


def RunMOPAC(filename: str) -> int:
    """Run the MOPAC on a well prepared input file."""
    
    itest=os.system("run_mopac7"+" "+filename)
    #time.sleep(1)
    return itest

############################################################################ 
def GetARCFile(inputmol: Chem.Mol) -> None:
    """
    #################################################################
    Get ARC file for each molecule
    #################################################################
    """
    
    FormatConversion(inputmol)
    itest = RunMOPAC('temp')
    if not itest:
        print(itest,'\t', 'Finshed successfully!')
    else:
        print(itest,'\t', 'Failed Finished........')
    os.remove('temp.dat')
    os.remove('temp.log')
    os.remove('temp.OUT')    
    #os.remove('temp.arc')
    oldpath = os.getcwd() + '/temp.arc'
    newpath = os.getcwd() + '/temp'
    os.rename(oldpath, newpath)


##############################################################################
# if __name__=="__main__":
#     mol='C1C=CCS1'
#     mol='SCCC(=O)N1[C@@H](CCC1)C(=O)OCC'
#     inputmol=pybel.readstring('smi',mol)
#     dir = GetARCFile(inputmol)
#     res=_ReadCoordinates(dir)
#     print(res)
#     Dispose(dir)
