#!/usr/bin/env python
"""Search for inlined CPP macros in ABINIT src files"""

__author__  = "M. Giantomassi"

import os
import re
import sys

# Files that will be checked.
re_srcfile = re.compile("\.([Ff]|[Ff]90|finc)$")

def is_srcfile(dirpath, fname):
  return re_srcfile.search(fname)

# List of macros that should not be inlined with if statements or other instructions.
# These macros are defined in src/incs/abi_common.h
# 
MACRO_NAMES = [
"ABI_ALLOCATE",
"ABI_DEALLOCATE",
"ABI_DATATYPE_ALLOCATE",
"ABI_DATATYPE_DEALLOCATE",
"ABI_MALLOC",
"ABI_FREE",
"ABI_CHECK_ALLOC",
]

# Regular expressions for each macro.
regexps = dict()

for name in MACRO_NAMES:
  regexps[name] = re.compile(name + " ?\(.*\)")
  #regexps[name] = re.compile(name + "\(.*\)") blanks between macro name and () are not permitted


def wrong_string(string):
  "Return an empty string if input string does not contain inlined macros."
  string = string.strip()
  for macro_name in MACRO_NAMES:
    pattern = regexps[macro_name]
    if macro_name in string and not string.startswith("!"): 
      s = re.sub(pattern, "", string, count=0).strip()
      if not s.startswith("!"): return s
  else:
    return ""

def abinit_test_generator():
  def test_func(abenv):
     """Search for inlined CPP macros in ABINIT src files"""
     top = abenv.apath_of("src")
     try:
       return main(top)
     except Exception:
       import sys
       raise sys.exc_info()[1] # Reraise current exception (py2.4 compliant)
  return {"test_func" : test_func}

def main(top):
  exit_status = 0
  for dirpath, dirnames, files in os.walk(top):
    for src in files:
      if is_srcfile(dirpath, src):
        fpath = os.path.join(dirpath,src)
        for lno, line in enumerate(file(fpath)):
          #
          s = wrong_string(line)
          if s:
            print "(INLINED MACRO at %s:%d):  %s " % (src, lno+1, line)
            exit_status += 1

  if exit_status > 0:
    err_msg = """
  Please, avoid instructions like:

    if (allocated(arr)) ABI_DEALLOCATE(arr)

  When the code is compiled in profile mode, indeed, ABI_DEALLOCATE expands to 
  the set of Fortran instructions:

    deallocate(arr)
    call memocc_abi()

  These instructions MUST be placed inside an "if then" "end if" block.
  This limitation can be lifted, but we need support at the level of the build system.
  For the time being, one has to use the more verbose form:

    if (allocated(arr)) then 
      ABI_DEALLOCATE(arr)
    end if

  This is the list of macros that cannot be inlined:

  %(MACRO_NAMES)s
""" % globals()
    print err_msg

  return exit_status

if __name__ == "__main__":

  if len(sys.argv) == 1: 
    top = "../../../src"
    print "-------------------------------------------------------"
    print " Searching for inlined CPP macros in ABINIT src files  "
    print "-------------------------------------------------------"
  else:
    top = sys.argv[1] 

  exit_status = main(top)
  sys.exit(exit_status)
