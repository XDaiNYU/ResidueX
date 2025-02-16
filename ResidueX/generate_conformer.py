import sys, os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from sklearn.cluster import DBSCAN
import rmsd #rmsd-1.4

def CoreAddConformer(core, refMol=None, mol = None, confId=-1):
    '''
    Add conformer to newCore
    '''
    ed_mol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    ref_core_Idxs = list(refMol.GetSubstructMatch(core))
    mol_core_Idxs = list(mol.GetSubstructMatch(core))
    
    mapOldIdxToNewIdx = {}
    
    conf_old = refMol.GetConformer(confId)
    conf_new = Chem.rdchem.Conformer()   
    
    for newAtomIdx,(oldAtomIdx1,oldAtomIdx2, oldAtom) in enumerate(zip(ref_core_Idxs, mol_core_Idxs, core.GetAtoms())):
        oldAtom2 = mol.GetAtomWithIdx(oldAtomIdx2)
        newAtom = Chem.rdchem.Atom(oldAtom2.GetAtomicNum())
        newAtom.SetFormalCharge(oldAtom2.GetFormalCharge())
        newAtom.SetHybridization(oldAtom2.GetHybridization())
        newAtomIdx = ed_mol.AddAtom(newAtom)
        point3D = conf_old.GetAtomPosition(oldAtomIdx1)
        conf_new.SetAtomPosition(newAtomIdx, point3D)
        mapOldIdxToNewIdx.update({oldAtomIdx2: newAtomIdx})
    
    for bond in mol.GetBonds():
        bIdx = bond.GetBeginAtomIdx()
        eIdx = bond.GetEndAtomIdx()
        if bIdx in mol_core_Idxs and eIdx in mol_core_Idxs:
            ed_mol.AddBond(mapOldIdxToNewIdx[bIdx], mapOldIdxToNewIdx[eIdx],order = bond.GetBondType())
            
    core1 = ed_mol.GetMol()
    core1.AddConformer(conf_new)
    return core1

def FindCore(mol1, mol2):
    MCS = rdFMCS.FindMCS([mol1,mol2],timeout=20,
                         atomCompare=Chem.rdFMCS.AtomCompare.CompareElements,
                         bondCompare=Chem.rdFMCS.BondCompare.CompareOrderExact)#,ringMatchesRingOnly=True) #, completeRingsOnly=True)
    
    #CompareElements
    core = Chem.MolFromSmarts(MCS.smartsString)


    return core

def GetCoordMap(refMol, mol, core):
    refMol_Idxs = list(refMol.GetSubstructMatch(core))
    mol_Idxs = list(mol.GetSubstructMatch(core))
    coordMap = {}
    
    for atomIdx1, atomIdx2 in zip(refMol_Idxs, mol_Idxs):
        pos = refMol.GetConformer().GetAtomPosition(atomIdx1)
        coordMap[atomIdx2] = pos
    # print(coordMap)
    # print(list(zip(mol_Idxs, refMol_Idxs)))
    # print(mol_Idxs)
    return coordMap, list(zip(mol_Idxs, refMol_Idxs)), mol_Idxs

def calc_energy(mol, corePart_nonH_Idxs, conformerId=-1, minimizeIts=5000):
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    mp.SetMMFFDielectricConstant(dielConst=80)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conformerId)
    for i in corePart_nonH_Idxs:
        ff.MMFFAddPositionConstraint(i, 0.5,10.0)

    ff.Minimize(maxIts=minimizeIts)
    results = ff.CalcEnergy()
    return results

def GetPosition(mol, atom_ids, confId=0):
    '''
    Get the position (x,y,z) matrix of the selected atoms
    '''
    mat = []
    for i in atom_ids:
        pos = mol.GetConformer(confId).GetAtomPosition(i)
        mat.append([pos.x, pos.y, pos.z])
    mat = np.array(mat)
    return mat

def cluster_conformers(mol, atom_ids, confIds, eps=0.5):
    """
    DBSCAN clustering the conformers of selected frag in molecule
    """    
    RMS_mat = []
    for i in confIds:
        vec_i = GetPosition(mol, atom_ids, i)
        mat = []
        for j in confIds:
            vec_j = GetPosition(mol, atom_ids, j)
            mat.append(rmsd.rmsd(vec_i, vec_j))
        RMS_mat.append(mat)
    RMS_mat = np.array(RMS_mat)
    ## clustering
    clustering = DBSCAN(eps, min_samples=1, metric = 'precomputed').fit(RMS_mat)
    _,sele_Ids = np.unique(clustering.labels_, return_index=True)
    sele_confIds = [confIds[x] for x in sele_Ids]
    
    return sele_confIds

def GenerateConformer(mol, refMol, path, sidechainIds, optimize=True):
    """
    mol --- molecule conformers need to be generated
    refMol --- reference molecule with available conformer
    path --- location where the conformers stored in sdf file
    sidechainId --- selected fragment for conformation clustering
    """
    #core = FindCore(refMol, mol)
    #new_core = CoreAddConformer(core, refMol, mol)
    core=Chem.MolFromSmarts('NCC(=O)')
    coordMap, atommap, corePart_Idxs = GetCoordMap(refMol, mol, core) #changed in 2023 12 11
    print(coordMap, atommap, corePart_Idxs)

    mol = Chem.AddHs(mol)

    conformerIds = AllChem.EmbedMultipleConfs(mol, maxAttempts=10,
                                              numConfs=200,
                                              pruneRmsThresh=0.5,
                                              coordMap=coordMap,
                                              ignoreSmoothingFailures=True,
                                              useRandomCoords=True,
                                              useExpTorsionAnglePrefs=False)
 
    print ('%d conformers generated in initial step'%len(conformerIds))

    if len(conformerIds) > 0:
        corePart_nonH_Idxs = [i for i in corePart_Idxs if mol.GetAtomWithIdx(i).GetAtomicNum() != 1]
        if optimize:
            ## conformer optimization
            for count, i in enumerate(conformerIds, 1):
                calc_energy(mol, corePart_nonH_Idxs, i)

        ## let's calculate the RMSD between cores
        rms_list = []
        for ID in conformerIds:
            rms, transForm = AllChem.GetAlignmentTransform(mol, refMol, prbCid=ID, refCid=-1, atomMap=atommap)
            rms_list.append(rms)
            AllChem.TransformMol(mol, transForm, confId=ID, keepConfs=True)
        print ("Min RMSD between mol's core and reference core: %.2f"%min(rms_list))

        if min(rms_list) < 1.0:
            os.makedirs(path, exist_ok=True)            
            id_list = [x for x in conformerIds if rms_list[x] < 1.0]
            print ("Get %d conformers with core RMSD < 1.0."%len(id_list))
            ## Clustering the conformations of selected fragment in molecule
            sele_confIds = cluster_conformers(mol, sidechainIds, id_list, eps=0.5)
            print ("Finally, %d conformers generated after sidechain clustering."%len(sele_confIds))
            ## write down the molecule in sdf files.
            for c,i in enumerate(sele_confIds,1):
                Chem.MolToMolFile(mol, '%s/%02d.sdf'%(path,c),confId=i)
        
    return core,coordMap, atommap, corePart_Idxs #No new core any more