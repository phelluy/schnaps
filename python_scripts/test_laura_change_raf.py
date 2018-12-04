"""
This script was made with the objective of changing the refinement
of a SCHNAPS simulation and puts the given one.
We suppose that:
    1. The refinement is defined as such:   "int theraf = 16;"
    2. You want to replace the '16' by any other number
    3. The script is being called from the build file, and the source file is
       one directory above and in the example directory.
    4. The line where the refinement is defined is at line 76 of the source code

To use the script use:

Author: Laura S. Mendoza
Modified date: 17/07/2017
"""
__all__ = ['find_and_replace','main']
__docformat__ = 'reStructuredText'

from tempfile import mkstemp
from shutil import move
from os import fdopen, remove


#==============================================================================
# MAIN FUNCTION (REPLACE FUNCTION)
#==============================================================================

# --------------------------------------------------------------------
# Function that goes through file and replaces the refinement file
# Arguments:
#   file_path : string path to the original file
#   subst     : string new line to replace the old line
#
def find_and_replace(file_path, subst):
    # We open the file
    with open(file_path, 'r') as file:
        # read all lines of file
        data = file.readlines()

    # now change the 76th line of the file
    newline = "  int theraf = " + subst + ";\n"
    data[75] = newline

    # and write everything back
    with open(file_path, 'w') as file:
        file.writelines( data )

def parse_input():

  import argparse, sys

  parser = argparse.ArgumentParser (
      prog        = 'python3 '+ sys.argv[0],
      description = 'Redefines the refinement of a SCHNAPS simulation.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument( metavar = 'ROOT',
                       dest    = 'root',
                       help    = 'relative path of the test file' )

  parser.add_argument( metavar = 'NEW_NAME',
                       dest    = 'new_name',
                       help    = '(integer) new refinement' )

  return parser.parse_args()

#==============================================================================
# SCRIPT FUNCTIONALITY
#==============================================================================

def main():

    # Parse input arguments
    print('')
    args = parse_input()
    print(args)
    print('')

    # Walk directory tree and change names of all library modules
    find_and_replace( args.root, args.new_name )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()
