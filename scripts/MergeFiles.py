#! /usr/bin/env python

import sys
import os
import re
import glob
from subprocess import Popen

# Sort the given list in the way that humans expect.
def sort_nicely( l ):
  convert = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
  l.sort( key=alphanum_key )

# Determine if file is an executable
def is_exe(fpath):
  return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

# Determine if executable in path
def which(program): 
  # Determine if executable in path
  fpath, fname = os.path.split(program)
  if fpath:
    if is_exe(program):
      return program
  else:
    for path in os.environ["PATH"].split(os.pathsep):
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file
  return None 


if __name__ == "__main__":

    if len(sys.argv) < 2:
      print("\n----------------Peridigm Exodus Outfile Merger----------------\n")
      print("Usage:  MergeFiles.py   <File Base Name>    <# of Processors>\n")
      sys.exit(1)

    path = sys.argv[0];
    base_name = sys.argv[1]
    num_proc = sys.argv[2]

    result = 0
    
    # Define path to location of MergeFiles to call epu and conjoin later
    start = 0
    end   = path.find('MergeFiles')
    newfile = path[start:end]
    path = newfile

    # remove the output .e file if it exists
    if os.path.exists(base_name+".e"):
      os.remove(base_name+".e")

    # open log file
    log_file_name = base_name + ".log"
    if os.path.exists(log_file_name):
      os.remove(log_file_name)
    logfile = open(log_file_name, 'w')

    # Generate listing of files <basename>-s* for merging
    files = glob.glob('*.e*')
    files_to_join = []
    for file in files:
      start = 0
      end   = file.find('.e')
      newfile = file[start:end];
      files_to_join.append(newfile)

    # Remove duplicates
    files_to_join = list(set(files_to_join))
    files_to_join.sort()

    # Delete empty exodus databases and rename with new number of processors
    ratio=100
    for file in files_to_join:
      maxsize=0
      for i in range(0,int(num_proc)):
        filename=file+".e."+num_proc+"."+str(i).zfill(len(str(num_proc)))
        if not os.path.exists(filename):
          break
        statinfo = os.stat(filename)
        filesize=statinfo.st_size
        if filesize>maxsize:
          maxsize=filesize

      new_num_proc=0
      for i in range(0,int(num_proc)):
        filename=file+".e."+num_proc+"."+str(i).zfill(len(str(num_proc)))
        if not os.path.exists(filename):
          break
        statinfo = os.stat(filename)
        filesize=statinfo.st_size
        if filesize>=maxsize/ratio:
          new_num_proc=new_num_proc+1
      new_num_proc=str(new_num_proc)

      j=0
      if int(new_num_proc)<int(num_proc):
        for i in range(0,int(num_proc)):
          filename=file+".e."+num_proc+"."+str(i).zfill(len(str(num_proc)))
          if not os.path.exists(filename):
            print filename + " does not exist."
            break
          statinfo = os.stat(filename)
          filesize=statinfo.st_size
          if filesize>=maxsize/ratio:
            os.rename(filename,file+".e."+new_num_proc+"."+str(j)).zfill(len(str(new_num_proc)))
            j=j+1
          else:
            os.remove(filename)

    # First merge all distributed exodus databases for each time stamp

      if is_exe(path+"epu"):
        command = [path+"epu", "-p", new_num_proc, file]
      elif which("epu") != None:
        command = ["epu", "-p", new_num_proc, file]
      else:
        print("Error: epu not found! Please execute this script from 'scripts' in your Peridigm build directory or put 'epu' in your path. 'epu' can be found in your Trilinos install 'bin' directory.")
        sys.exit(-1)
      p = Popen(command, stdout=logfile, stderr=logfile)
      return_code = p.wait()
      if return_code != 0:
          result = return_code

    # Check for any "-s" files
    check = glob.glob('*-s*')
    if len(check) != 0:
      # Add .e to list of files for conjoin
      files_to_conjoin = []
      for file in files_to_join:
        newfile = file + ".e"
        files_to_conjoin.append(newfile)
      sort_nicely(files_to_conjoin)
      # Now combine time series from all databases
      if is_exe(path+"conjoin"):
        command = [path+"conjoin", "-output", base_name+".e"]
      elif which("conjoin") != None:
        command = ["conjoin", "-output", base_name+".e"]
      else:
        print("Error: conjoin not found! Please execute this script from 'scripts' in your Peridigm build directory or put 'conjoin' in your path. 'conjoin' can be found in your Trilinos install 'bin' directory.")
        sys.exit(-1)
      for file in files_to_conjoin:
        command.append(file)
      p = Popen(command, stdout=logfile, stderr=logfile)
      return_code = p.wait()
      if return_code != 0:
        result = return_code
    
      # If error print message to screen
      if result != 0:
        print("\n Error when merging files! Check .log file for error details.\n")

    sys.exit(result)
