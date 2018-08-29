#!/usr/bin/env/ python

#-----------------------------------------------------------------------------
#							RF matrix calculator 
#-----------------------------------------------------------------------------

# import required libraries
import numpy as np
import codecs, json, sys
from tqdm import tqdm
import ete3 as ete
import os
import pandas as pd
try: from optparse import OptionParser
except ImportError:
    print "\n\tError: OptionParser (optparse) is not installed/loaded"
    sys.exit()

def listabs(directory):
    '''Returns the absolute path of all items in a directory
    
    Parameters
    ----------
    directory : path to directory'''
    return [os.path.join(directory, filename) for filename in os.listdir(directory)]

def constructmatrix(directory):
	print("Hello")
	trees = [filename for filename in listabs(directory) if filename.endswith('.tre')]

	eteRF = np.zeros((len(trees), len(trees)))
	i=0
	j=0
	for t1 in tqdm(trees):
	    with open(t1, 'r') as t: tree1=ete.Tree(newick=t.readline())
	    for t2 in trees:
	        with open(t2, 'r') as t: tree2=ete.Tree(newick=t.readline())
	        
	        tdict = tree1.compare(ref_tree=tree2, unrooted=True)
	        try: eteRF[i,j] = tdict['rf']
	        except ValueError:
	            eteRF[i,j] = len(tree1.get_descendants())+len(tree2.get_descendants())
	        j+=1
	    j=0
	    i+=1

	return eteRF

def main():
    
    print('\n++++++++ RF matrix calculator ++++++++\n')
    
    parser = OptionParser(prog="RFmatcalc", usage="%prog [options]", version="%prog 1.0")
    parser.add_option("-d", "--directory",
                      action="store",
                      dest="directory",
                      help="path to directory containing newick trees")
    parser.add_option("-o", "--outfile",
                      action="store",
                      dest="outfile",
                      help="outfile name")
    
    (options, args) = parser.parse_args()
    if not all((options.directory, options.outfile)):
        print "\n\tmust include all options -d and -o"
        sys.exit()
        
    eteRF = constructmatrix(directory=options.directory)

    np.savetxt(outfile, eteRF, delimiter=',')

#+------------------------------------------------------------------------
if __name__ == '__main__':
    main()







