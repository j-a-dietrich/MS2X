import math

import numpy

from rdkit import Chem
from rdkit.Chem import Graphs, rdchem, rdMolDescriptors
from rdkit.ML.InfoTheory import entropy

    
def _AssignSymmetryClasses(mol, vdList, bdMat, forceBDMat, numAtoms, cutoff):
  """
     Used by BertzCT

     vdList: the number of neighbors each atom has
     bdMat: "balaban" distance matrix

  """
  if forceBDMat:
    bdMat = Chem.GetDistanceMatrix(mol, useBO=1, useAtomWts=0, force=1, prefix="Balaban")
    mol._balabanMat = bdMat

  keysSeen = []
  symList = [0] * numAtoms
  for i in range(numAtoms):
    tmpList = bdMat[i].tolist()
    tmpList.sort()
    theKey = tuple(['%.4f' % x for x in tmpList[:cutoff]])
    try:
      idx = keysSeen.index(theKey)
    except ValueError:
      idx = len(keysSeen)
      keysSeen.append(theKey)
    symList[i] = idx + 1
  return tuple(symList)


def _LookUpBondOrder(atom1Id, atom2Id, bondDic, numAtoms):
  """
     Used by BertzCT
  """
  if atom1Id < atom2Id:
    theKey = (atom1Id, atom2Id)
  else:
    theKey = (atom2Id, atom1Id)
  tmp = bondDic[theKey]
  if tmp == Chem.BondType.AROMATIC:  
    tmp = 2.2
  else:
    tmp = float(tmp)
  # tmp = int(tmp)
  return tmp


def _CalculateEntropies(connectionDict, vdList, bondOrders, symmetryClasses):
  """
     Used by BertzCT
  """

  degreeIE = entropy.InfoEntropy(numpy.array(vdList))
  bondIE = entropy.InfoEntropy(numpy.array(bondOrders))
  symmetryIE = entropy.InfoEntropy(numpy.array(symmetryClasses))
    
  connectionList = list(connectionDict.values())
  totConnections = sum(connectionList)
  connectionIE = (totConnections) * (entropy.InfoEntropy(numpy.array(connectionList, float)) +
                                   math.log(totConnections) / _log2val)

  return connectionIE

_log2val = math.log(2)

def _CreateBondDictEtc(mol, numAtoms):
  """ _Internal Use Only_
     Used by BertzCT

  """
  bondDict = {}
  nList = [None] * numAtoms
  vdList = [0] * numAtoms
  for aBond in mol.GetBonds():
    atom1 = aBond.GetBeginAtomIdx()
    atom2 = aBond.GetEndAtomIdx()
    if atom1 > atom2:
      atom2, atom1 = atom1, atom2
    if not aBond.GetIsAromatic():
      bondDict[(atom1, atom2)] = aBond.GetBondType()
    else:
      # mark Kekulized systems as aromatic
      bondDict[(atom1, atom2)] = Chem.BondType.AROMATIC
    if nList[atom1] is None:
      nList[atom1] = [atom2]
    elif atom2 not in nList[atom1]:
      nList[atom1].append(atom2)
    if nList[atom2] is None:
      nList[atom2] = [atom1]
    elif atom1 not in nList[atom2]:
      nList[atom2].append(atom1)

  for i, element in enumerate(nList):
    try:
      element.sort()
      vdList[i] = len(element)
    except Exception:
      vdList[i] = 0
  return bondDict, nList, vdList


def BertzCT(mol, cutoff=100, dMat=None, forceDMat=1):
  """ A topological index meant to quantify "complexity" of molecules.

     Consists of a sum of two terms, one representing the complexity
     of the bonding, the other representing the complexity of the
     distribution of heteroatoms.

     From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)

     "cutoff" is an integer value used to limit the computational
     expense.  A cutoff value tells the program to consider vertices
     topologically identical if their distance vectors (sets of
     distances to all other vertices) are equal out to the "cutoff"th
     nearest-neighbor.

     **NOTE**  The original implementation had the following comment:
         > this implementation treats aromatic rings as the
         > corresponding Kekule structure with alternating bonds,
         > for purposes of counting "connections".
       Upon further thought, this is the WRONG thing to do.  It
        results in the possibility of a molecule giving two different
        CT values depending on the kekulization.  For example, in the
        old implementation, these two SMILES:
           CC2=CN=C1C3=C(C(C)=C(C=N3)C)C=CC1=C2C
           CC3=CN=C2C1=NC=C(C)C(C)=C1C=CC2=C3C
        which correspond to differentk kekule forms, yield different
        values.
       The new implementation uses consistent (aromatic) bond orders
        for aromatic bonds.

       THIS MEANS THAT THIS IMPLEMENTATION IS NOT BACKWARDS COMPATIBLE.

       Any molecule containing aromatic rings will yield different
       values with this implementation.  The new behavior is the correct
       one, so we're going to live with the breakage.

     **NOTE** this barfs if the molecule contains a second (or
       nth) fragment that is one atom.

  """
  atomTypeDict = {}
  connectionDict = {}
  numAtoms = mol.GetNumAtoms()
  if forceDMat or dMat is None:
    if forceDMat:
      # nope, gotta calculate one
      dMat = Chem.GetDistanceMatrix(mol, useBO=0, useAtomWts=0, force=1)
      mol._adjMat = dMat
    else:
      try:
        dMat = mol._adjMat
      except AttributeError:
        dMat = Chem.GetDistanceMatrix(mol, useBO=0, useAtomWts=0, force=1)
        mol._adjMat = dMat

  if numAtoms < 2:
    return 0

  bondOrders = []
  bondDict, neighborList, vdList = _CreateBondDictEtc(mol, numAtoms)
  symmetryClasses = _AssignSymmetryClasses(mol, vdList, dMat, forceDMat, numAtoms, cutoff)
  #print('Symmm Classes:',symmetryClasses)
  for atomIdx in range(numAtoms):
    hingeAtomNumber = mol.GetAtomWithIdx(atomIdx).GetAtomicNum()
    atomTypeDict[hingeAtomNumber] = atomTypeDict.get(hingeAtomNumber, 0) + 1

    hingeAtomClass = symmetryClasses[atomIdx]
    numNeighbors = vdList[atomIdx]
    
    for i in range(numNeighbors):
      neighbor_iIdx = neighborList[atomIdx][i]
      NiClass = symmetryClasses[neighbor_iIdx]
      bond_i_order = _LookUpBondOrder(atomIdx, neighbor_iIdx, bondDict, numAtoms) 
      #print('\t',atomIdx,i,hingeAtomClass,NiClass,bond_i_order)
      bondOrders.append(bond_i_order)
      if (bond_i_order > 1) and (neighbor_iIdx > atomIdx):
        numConnections = bond_i_order * (bond_i_order - 1) / 2
        connectionKey = (min(hingeAtomClass, NiClass), max(hingeAtomClass, NiClass))
        connectionDict[connectionKey] = connectionDict.get(connectionKey, 0) + numConnections

      for j in range(i + 1, numNeighbors):
        neighbor_jIdx = neighborList[atomIdx][j]
        NjClass = symmetryClasses[neighbor_jIdx]
        bond_j_order = _LookUpBondOrder(atomIdx, neighbor_jIdx, bondDict, numAtoms)
        numConnections = bond_i_order * bond_j_order
        connectionKey = (min(NiClass, NjClass), hingeAtomClass, max(NiClass, NjClass))
        connectionDict[connectionKey] = connectionDict.get(connectionKey, 0) + numConnections
        

  if not connectionDict:
    connectionDict = {'a': 1}

  return _CalculateEntropies(connectionDict, vdList, bondOrders, symmetryClasses)