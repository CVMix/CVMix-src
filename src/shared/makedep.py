#!/usr/bin/env python

# usage: makedep.py $(DEP_FILE) $(OBJ_DIR) $(SRC_DIR) [$(SRC_DIR2)]

# Generate $DEP_FILE in $OBJ_DIR (arguments 1 and 2, respectively)
# Read in every file in $SRC_DIR and $SRC_DIR2 (arguments 3 and 4)
# Only depend on modules located in $SRC_DIR or $SRC_DIR2

import os, sys, re

try:
  dep_file = sys.argv[1]
except:
  dep_file = "depends.d"

try:
  obj_dir = sys.argv[2]
except:
  obj_dir = '.'

try:
  src_dir = sys.argv[3]
except:
  src_dir = '.'

try:
  src_dir2 = sys.argv[4]
except:
  src_dir2 = src_dir

fout = open(dep_file, 'w')
files_in_src_dir =  os.listdir(src_dir)
if src_dir != src_dir2:
  files_in_src_dir.extend(os.listdir(src_dir2))

for src_file in files_in_src_dir:
  file_name, file_ext = os.path.splitext(src_file)
  if file_ext == '.F90':
    try:
      fin = open(src_dir+'/'+src_file,"r")
    except:
      fin = open(src_dir2+'/'+src_file,"r")
    for line in fin:
      if re.match('^ *[Uu][Ss][Ee]',line):
        line_array = line.split()
        # statements are usually "use module, only : subroutine"
        # so we need to strip away the , to get the module name
        file_used = line_array[1].split(',')[0]
        if file_used+'.F90' in files_in_src_dir:
          print file_name+'.o depends on '+file_used+'.o'
          fout.write(obj_dir+'/'+file_name+'.o: '+obj_dir+'/'+file_used+'.o\n')
    fin.close
