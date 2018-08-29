#!/usr/bin/env python

#+------------------------------------------------------------------------
#                 ++++++++++++++++ G2Genes ++++++++++++++++                
#
# This script reads in a .gphocs file from pyrad/ipyrad output and splits 
# it into individual .nex and .phy files. 
#                                                         30 November 2017
#                                                            Ian S. Gilman 
#+------------------------------------------------------------------------

try: import sys
except ImportError:
    print "\n\tError: sys is not installed/loaded"
    sys.exit()
try: import os
except ImportError:
    print "\n\tError: os is not installed/loaded"
    sys.exit()
try: from optparse import OptionParser
except ImportError:
    print "\n\tError: OptionParser (optparse) is not installed/loaded"
    sys.exit()
try: from Bio import AlignIO
except ImportError:
    print "\n\tError: Bio is not installed/loaded"
    sys.exit()
try: from tqdm import tqdm
except ImportError:
    print "\n\tError: tqdm is not installed/loaded."
    print "\n\tTo install, try 'pip install tqdm' or 'conda install -c conda-forge tqdm'"
    sys.exit()

#+------------------------------------------------------------------------
def trim_path(path, n=1, sep='/'):
    '''The function trim_path takes a string and lops off the final substring that has been delimited
    
    Parameters
    ----------
    path: string (usually a file path) to be trimmed
    n: how many substrings to lop off (default=1)'''

    try: int(n)
    except TypeError:
        print('\nn must be an integer')
        sys.exit()      
    for i in range(n):
        if len(path.split(str(sep)))<=1:
            path=os.getcwd()
        else:
            path = str(sep).join(path.split(str(sep))[:-1])
    return path

#+------------------------------------------------------------------------
class Vividict(dict):
    '''An infinitely nestable dictionary created by Aaron Hall. See
    https://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries.'''
    def __missing__(self, key):
        value = self[key] = type(self)() # retain local pointer to value
        return value                     # faster to return than dict lookup

#+------------------------------------------------------------------------
def gphocs2dict(gphocs_file):
    ''' gphocs2dict reads in a gphocs file and parses it into a dict for
    downstream separation into nexus and phylip files
    Parameters
    ----------
    gphocs_file: a gphocs formatted file, usually from pyRAD output
    '''
    with open(gphocs_file, 'r') as f:
        data = f.readlines()

    gphocs_dict = Vividict()
    for i, line in enumerate(data):
        elements = len(line.split())
        if elements not in [2,3]:
            continue
        elif len(line.split())==3:
            locus = line.split()[0]
            n_sites = line.split()[2]
            gphocs_dict[locus]['n_taxa'] = line.split()[1]
            gphocs_dict[locus]['n_sites'] = line.split()[2]
            gphocs_dict[locus]['seqs'] = []
            gphocs_dict[locus]['missing'] = 'N'
            gphocs_dict[locus]['gap'] = '-'
            gphocs_dict[locus]['datatype'] = 'DNA'
        else:
            gphocs_dict[locus]['seqs'].append(line)
            # Replace ambiguity character
            if '?' in line:
                line.replace('?', 'N')
            
    return gphocs_dict

#+------------------------------------------------------------------------
def gphocs2nexphy(gphocs_file):
    '''gpohcs2nexphy converts a gphocs file to nexus and phylip formats via
    gphocs2dict, which parses a gpohcs file into a dictionary so that it 
    can be reassembled in another format.

        Begin DATA; [locus1]
        Dimensions NTAX= 2 NCHAR=28;
        Format MISSING=N GAP=- DATATYPE=DNA;
            Matrix

        Sample_1
        ACGTACGTACGTACGTACGTACGTACGT
        Sample_2
        ACGTACGTACGTACGTACGTACGTACGT
            ;
        END;

    Parameters
    ----------
    gphocs_file: a gphocs formatted file, usually from pyRAD output'''
    
    nexus_dir = os.path.join(trim_path(gphocs_file), 'nexus_gfiles')
    try: os.mkdir(nexus_dir)
    except OSError: print('Directory %s already exists') % (nexus_dir)
        
    phylip_dir = os.path.join(trim_path(gphocs_file), 'phylip_gfiles')
    try: os.mkdir(phylip_dir)
    except OSError: print('Directory %s already exists') % (phylip_dir)
    
    gdict = gphocs2dict(gphocs_file=gphocs_file)
    i=0
    for key in tqdm(gdict.keys()):
        nex_out = os.path.join(nexus_dir, key+'.nexus')
        phy_out = os.path.join(phylip_dir, key+'.phylip')
        
        # Write PAUP block
        with open(nex_out, 'w') as f:
            nexus_header = ('#NEXUS\nBegin DATA; [%s]\n\tDimensions NTAX=%s NCHAR=%s;\n\tFormat MISSING=%s GAP=%s DATATYPE=%s;\n\tMatrix\n') %(key, gdict[key]['n_taxa'], gdict[key]['n_sites'], gdict[key]['missing'], gdict[key]['gap'], gdict[key]['datatype'])
            f.write(nexus_header)
            f.write(''.join(gdict[key]['seqs']))
            nexus_footer = '\t;\nEND;'
            f.write(nexus_footer)
        
        input_handle = open(nex_out, "rU")
        output_handle = open(phy_out, "w")
        
        alignments = AlignIO.parse(input_handle, "nexus")
        AlignIO.write(alignments, output_handle, "phylip-relaxed")

        output_handle.close()
        input_handle.close()
        
        i+=1
    print('\nWrote %d nexus to %s') % (i, nexus_dir)
    print('Wrote %d phylip to %s') % (i, phylip_dir)

#+------------------------------------------------------------------------
def main():
    print('\n++++++++ G2Genes ++++++++\n')
    
    parser = OptionParser(prog="G2Genes", usage="%prog [options]", version="%prog 1.0")
    parser.add_option("-i", "--infile",
                      action="store",
                      dest="infile",
                      help="path to gphocs file")
    
    (options, args) = parser.parse_args()

    if not all(options.infile):
        print "\n\tmust include option -i (input gphocs file)"
        sys.exit()
    
    gphocs2nexphy(gphocs_file=options.infile)

#+------------------------------------------------------------------------
if __name__ == "__main__":
    main()

