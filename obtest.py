import openbabel

def btd_info(fileformat,input_name):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(fileformat, "report")
#    input_name="NaT2_2.txyz"
#output_name="NaT2_2.rep"
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_name)   # Open Babel will uncompress automatically
#obConversion.WriteFile(mol,output_name)   # Open Babel will uncompress automatically
    outMDL = obConversion.WriteString(mol)
#    print outMDL
    output = outMDL.split('\n')
#print len(output)
    ba_flag = 0
#    ba_ct = 0
    bangle1=[]
    bangle2=[]
    bangle3=[]
    ta_flag = 0
#    ta_ct = 0
    tangle1=[]
    tangle2=[]
    tangle3=[]
    tangle4=[]
    for line in output:
        mystring = line
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        if (ba_flag == -1 and len(mylist) > 6):
#            ba_ct += 1
            bangle1.append(int(mylist[0]))
            bangle2.append(int(mylist[1]))
            bangle3.append(int(mylist[2]))
        elif (ba_flag == -1):
            ba_flag == 0
        if (mylist[0] == 'BOND' and mylist[1] == 'ANGLES'):
            ba_flag = -1
        if (ta_flag == -1 and len(mylist) > 3):
#            ta_ct += 1
            tangle1.append(int(mylist[0]))
            tangle2.append(int(mylist[1]))
            tangle3.append(int(mylist[2]))
            tangle4.append(int(mylist[3]))
        elif (ta_flag == -1):
            ta_flag == 0
        if (mylist[0] == 'TORSION' and mylist[1] == 'ANGLES'):
            ta_flag = -1
#print mol.GetBond(1,2)
#print mol.NumAtoms()
#print mol.NumBonds()
#print mol.NumResidues()
#print len(outMDL)
    obConversion.SetInAndOutFormats(fileformat, "molreport")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_name)   # Open Babel will uncompress automatically
#obConversion.WriteFile(mol,output_name)   # Open Babel will uncompress automatically
    outMDL = obConversion.WriteString(mol)
#    print outMDL
    output = outMDL.split('\n')
#    pair_ct=0
    pair1=[]
    pair2=[]
    for line in output:
        mystring = line
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        if (mylist[0] == "BOND:"):
#            pair_ct += 1
            pair1.append(int(mylist[3]))
            pair2.append(int(mylist[5]))
#    for ring in openbabel.OB_RINGTYPES_MOL(mol)        
###   mol.Clear();
# Can't clear yet
    return pair1,pair2,bangle1,bangle2,bangle3,tangle1,tangle2,tangle3,tangle4,mol

def openbabel_count_aromatic_rings(mol):
    count = 0
    for ring in mol.GetSSSR():
# Note: the OB implementation is wrong. It assumes that if all
# atoms in the ring are aromatic then the ring itself must be
# aromatic. This is not necessarily true.
        a = int(ring.Size())
        b = [ring.IsAromatic()]
        c = [ring.GetType()]
        print a,b[0],c[0]
        print ring._path
        print ring.Size(), ring.IsAromatic(), ring.GetType(), 
        if ring.IsAromatic():
            count += 1
    return count
# 
def openbabel_count_hetero_rings(mol):
    count = 0
    for ring in mol.GetSSSR():
        if (ring.IsAromatic() and
               any(mol.GetAtom(atom_id).GetAtomicNum() != 6 for atom_id in ring._path)):
            test = [any(mol.GetAtom(atom_id) for atom_id in ring._path )]
            count += 1
    print test,[mol.GetAtom(atom_id)for atom_id in ring._path]   
    return count
"""
#
# test code
#

pair1,pair2,bangle1,bangle2,bangle3,tangle1,tangle2,tangle3,tangle4,mol = btd_info("txyz","/home/winokur/hoomd_test/NaT2_2.txyz")
print openbabel_count_aromatic_rings(mol)
for i in range(len(pair1)):
    print i,pair1[i],pair2[i]
for i in range(len(bangle1)):
    print i,bangle1[i],bangle2[i],bangle3[i]
for i in range(len(tangle1)):
    print i,tangle1[i],tangle2[i],tangle3[i],tangle4[4]
"""