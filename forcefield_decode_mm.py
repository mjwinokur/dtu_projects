#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:01:05 2017

@author: winokur
"""
import re
import numpy as np
def mm3_decode(fffile,ptypes,aatypes,bdtypes,angtypes,tortypes):
    ic=len(ptypes)
    amass = np.empty(ic,dtype=float)
    vdw1 = np.empty(ic, dtype=float)
    vdw2 = np.empty(ic, dtype=float)
    vdwpr1 = np.empty(ic,dtype=float)
    vdwpr2 = np.empty(ic,dtype=float)
    bond1 = []
    bond2 = []
    bdtypes_tot = []
    ang1 = []
    ang2 =  []
    angtypes_tot = []
    tor1 = []
    tor2 = []
    tor3 = []
    tortypes_tot = []

# Because there are multiply type of pair and 
# angles type bonds the arrays may need to be appended
    with open(fffile, "r") as file:
        text = file.readlines()
        for line in text:
#        print line
            mystring = line
            newline = ' '.join(mystring.split())
            mylist = newline.split(" ")
#            alen = len(mylist)
#            print alen,mylist
            if (mylist[0]=='atom'):
# Need to resplit with the "xxx yyy" stuff accounted for
#            lcut = re.findall('".+"',line)
                mystring = re.sub('".+"',' ',line)
                newline = ' '.join(mystring.split())
                mylist = newline.split(" ")
                iatm=int(mylist[1])
                if iatm in ptypes:
                    ic = ptypes.index(iatm)
                    amass[ic]=(float(mylist[4]))
#                    print 'atom ',iatm,mylist,ic
            elif (mylist[0]=='vdw'):
                iatm=int(mylist[1])
                if iatm in ptypes:
                    ic = ptypes.index(iatm)
                    vdw1[ic]=(float(mylist[2]))
                    vdw2[ic]=(float(mylist[3]))
                    print 'vdw ',iatm,mylist,ic
            elif (mylist[0]=='vdwpr'):
                iatm=int(mylist[1])
                iatm2=int(mylist[2])
                if ((iatm in ptypes) and (iatm2 in ptypes)):
                    ic = ptypes.index(iatm)
                    vdwpr1[ic]=(float(mylist[3]))
                    vdwpr2[ic]=(float(mylist[4]))
                    ic = ptypes.index(iatm2)
                    vdwpr1[ic]=(float(mylist[3]))
                    vdwpr2[ic]=(float(mylist[4]))
                    print 'vdwpr ',iatm,iatm2,mylist,ic
            elif (mylist[0]=='bond' or mylist[0]=='bond5'):
                v1=int(mylist[1])
                v2=int(mylist[2])
                if ((v1 in ptypes) and (v2 in ptypes)):
                    t1 = aatypes[ptypes.index(v1)]+mylist[1]
                    t2 = aatypes[ptypes.index(v2)]+mylist[2]
                    tf = t1+t2
                    tr = t2+t1  # Reverse sequence is needed MM3 just goes from low to high
                    if ((tf == tr) and (tf in angtypes)):
                        if ( (mylist[0]+'_'+tf) not in bdtypes_tot):
                            bond1.append(float(mylist[3]))
                            bond2.append(float(mylist[4]))
                            bdtypes_tot.append(mylist[0]+'_'+tf)
                    elif ((tf in bdtypes) or (tr in bdtypes) ):
                        if ( ((mylist[0]+'_'+tf) not in bdtypes_tot) and ((mylist[0]+'_'+tr) not in bdtypes_tot)):
                            bond1.append(float(mylist[3]))
                            bond2.append(float(mylist[4]))
                            bdtypes_tot.append(mylist[0]+'_'+tf+'_'+tr)
#                    print mylist
#                    print mylist[0],v1,v2,'index:',idx,tf,bond1[idx],bond2[idx],'\n'
            elif (mylist[0]=='hbond'):
                v1=int(mylist[1])
                v2=int(mylist[2])
                if ((v1 in ptypes) and (v2 in ptypes)):
                    print mylist[0],"Not done"
            elif (mylist[0]=='angle' or mylist[0]=='angle5'):
                v1=int(mylist[1])
                v2=int(mylist[2])
                v3=int(mylist[3])
                if ((v1 in ptypes) and (v2 in ptypes) and (v3 in ptypes)):
                    t1 = aatypes[ptypes.index(v1)]+mylist[1]
                    t2 = aatypes[ptypes.index(v2)]+mylist[2]
                    t3 = aatypes[ptypes.index(v3)]+mylist[3]
                    tf = t1+t2+t3
                    tr = t3+t2+t1
                    if ((tr == tf) and (tf in angtypes)): # no duplicates
                        if ((mylist[0]+'_'+tf) not in angtypes_tot): # don't use if not needed
                            ang1.append(float(mylist[4]))
                            ang2.append(float(mylist[5]))
                            angtypes_tot.append(mylist[0]+'_'+tf)
                    elif ((tf in angtypes) or (tr in angtypes)) : # no duplicates
                        if (((mylist[0]+'_'+tr+'_'+tf) not in angtypes_tot) and ((mylist[0]+'_'+tf+'_'+tr) not in angtypes_tot)): # don't use if not needed
                            ang1.append(float(mylist[4]))
                            ang2.append(float(mylist[5]))
                            angtypes_tot.append(mylist[0]+'_'+tf+'_'+tr)
#                    print mylist
#                    print mylist[0],v1,v2,v3,'index:',idx,tf,angle1[idx],angle2[idx],'\n'
# Need logic to increase force type arrays if needed
# A list of used force ids and then add multiples                        
            elif (mylist[0]=='strbnd'):
                v1=int(mylist[1])
                v2=int(mylist[2])
                v3=int(mylist[3])
                if ((v1 in ptypes) and (v2 in ptypes) and (v3 in ptypes)):
                    print mylist[0],"Not done" 
            elif (mylist[0]=='angang'):
                iatm=int(mylist[1])
                if ((iatm in ptypes)):
                    print mylist[0],"Not done"
            elif (mylist[0]=='opbend'):
                iatm=int(mylist[1])
                iatm2=int(mylist[2])
                if ((iatm in ptypes) and (iatm2 in ptypes)):
                    print mylist[0],"Not done"
# MM3looks to use OPLS style
            elif ( (mylist[0]=='torsion') or (mylist[0]=='torsion5') or (mylist[0]=='torsion4') ):
                v1=int(mylist[1])
                v2=int(mylist[2])
                v3=int(mylist[3])
                v4=int(mylist[4])
                if ((v1 in ptypes) and (v2 in ptypes) and (v3 in ptypes) and (v4 in ptypes)): # Start to sift and winnow
                    t1 = aatypes[ptypes.index(v1)]+mylist[1]
                    t2 = aatypes[ptypes.index(v2)]+mylist[2]
                    t3 = aatypes[ptypes.index(v3)]+mylist[3]
                    t4 = aatypes[ptypes.index(v4)]+mylist[4]
                    tf = t1+t2+t3+t4
                    tr = t4+t3+t2+t1
                    prefix = mylist[0]+'_'
                    if ((tf == tr) and (tf in tortypes)):
                        if (((prefix+tr+'_'+tf) not in tortypes_tot) and ((prefix+tf+'_'+tr) not in tortypes_tot)):
                            tor1.append(float(mylist[5]))
                            tor2.append(float(mylist[8]))
                            tor3.append(float(mylist[11]))
                            tortypes_tot.append(prefix+tf)
#
                    elif (((tf in tortypes)) or (tr in tortypes)):
                        if (((prefix+tr+'_'+tf) not in tortypes_tot) and ((prefix+tf+'_'+tr) not in tortypes_tot)):
                            tor1.append(float(mylist[5]))
                            tor2.append(float(mylist[8]))
                            tor3.append(float(mylist[11]))
                            tortypes_tot.append(prefix+tf+'_'+tr)
            elif (mylist[0]=='strtors'):
                iatm2=int(mylist[2])
                iatm3=int(mylist[3])
                if ((iatm2 in ptypes) and (iatm3 in ptypes)):
                    print mylist[0],"Not done"
#            else:
#                print "Not on list"
              
#print amass
# but the order is wrong
    file.close()
    return amass,bdtypes_tot,bond1,bond2,angtypes_tot,ang1,ang2,tortypes_tot,tor1,tor2,tor3,vdw1,vdw2,vdwpr1,vdwpr2
