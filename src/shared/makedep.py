#!/usr/bin/env python

# usage: makedep.py $(SRC_DIR) $(OBJ_DIR) $(DEP_FILE)

# Read in every file in $SRC_DIR (argument 1)
# Generate $DEP_FILE in $OBJ_DIR (argument 3 and 2, respectively)
# If no $DEP_FILE, generates depends.d

import os, sys, re

try:
  src_dir = sys.argv[1]
except:
  src_dir = '.'

try:
  obj_dir = sys.argv[2]
except:
  obj_dir = '.'

try:
  dep_file = sys.argv[3]
except:
  dep_file = "depends.d"

fout = open(dep_file, 'w')
files_in_src_dir = os.listdir(src_dir)
for src_file in files_in_src_dir:
  file_name, file_ext = os.path.splitext(src_file)
  if file_ext == '.F90':
    fin = open(src_dir+'/'+src_file,"r")
    for line in fin:
      if re.match('^ *[Uu][Ss][Ee]',line):
        line_array = line.split()
        # statements are usually "use module, only : subroutine"
        # so we need to strip away the , to get the module name
        file_used = line_array[1].split(',')[0]
        print file_name+'.o depends on '+file_used+'.o'
        fout.write(obj_dir+'/'+file_name+'.o: '+obj_dir+'/'+file_used+'.o\n')
    fin.close


