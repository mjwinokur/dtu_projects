#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 13:17:30 2017

@author: winokur
"""
#import math
from mw_library import read_txyz
from mw_library import iofiles
from mw_library import write_xyz
import sys
txyz_file_in,xyz_file_out,text = iofiles(sys.argv[1:])
if (len(txyz_file_in)== 0):
    print ('No input txyz file:')
    raw_input()
if (len(xyz_file_out)== 0):
    print ('No output xyz file:')
    raw_input()
print 'Input txyz file is: ',txyz_file_in
print 'Output xyz file is: ',xyz_file_out
l2_num, aax,aay,aaz,atype, astring, abt, uc, header = read_txyz(txyz_file_in)
write_xyz(xyz_file_out,atype,aax,aay,aaz)
print 'Done'