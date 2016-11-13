#!/usr/bin/env python

# ============================================================================ #
# poremapstats.py                                                              #
# ============================================================================ #
'''
map and generate mapping statistics for Oxford Nanopore reads

Given a BAM file of primary and/or secondary alignments for a set of reads,
prints a file of statistics for the alignments, for each read, and for the run.

It was originally designed for generating statistics for 2D reads from one
run of an Oxford Nanopore MinION flow cell but it could be used for generating
statistics from any BAM file.

Input:
- entire set of mapped Oxford Nanopore 2D reads in BAM format
- target and control reference sequence(s) in FASTA format
- the mapping program and parameters to use for the mapping and random
  alignment classification steps

Output:
- outdir/outprefix_aligstats.txt
- outdir/outprefix_readstats.txt
- outdir/outprefix_runstats.txt
- outdir/outprefix_readid_flag_refcontigid_pos_mapq[_rand]_aln.fasta [optional]
'''
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# April 2015                                                                   #
# ============================================================================ #

# ============================================================================ #
# Import system-level Python modules                                           #
# ============================================================================ #

import argparse, itertools, numpy as np, os, shlex, subprocess as sp, sys
from subprocess import PIPE
from Bio import SeqIO

# ============================================================================ #
# Set _bindir and _profilepath, then import onttools                           #
#                                                                              #
# This complicated pre-onttools command-line parsing was introduced so that    #
# each Ont program runs independently using auxiliary config and programs      #
# located in the same 'bin' directory, and can be run using qsub using         #
# absolute pathnames (otherwise it doesn't work).                              #
# ============================================================================ #

_errorargsinvalid = 3

_bindir = None
if len(sys.argv) > 1:
    _bindir = None
    try:
        for i in range(1, len(sys.argv)):
            if sys.argv[i].startswith('--bindir'):
                if '=' in sys.argv[i]:
                    _bindir = os.path.realpath(os.path.expandvars('='.join(sys.argv[i].split('=')[1:])))
                else:
                    _bindir = os.path.realpath(os.path.expandvars(sys.argv[sys.arvg[i+1]]))
    except:
        _bindir = os.path.dirname(os.path.realpath(os.path.expandvars(sys.argv[0])))
    if not _bindir:
        _bindir = os.path.dirname(os.path.realpath(os.path.expandvars(sys.argv[0])))
    if _bindir and os.path.exists(_bindir):
        try:
            sys.path.insert(0, _bindir)
            from poremaplib.poretypes import Poretypes
            from poremaplib.poreerr import Poreerr
            from poremaplib.porelog import Porelog
            from poremaplib.poreargs import Poreargs
        except:
            sys.stderr.write('Erro: Failed to import onttools module\n')
            sys.exit(_errorargsinvalid)
    else:
        sys.path.insert(0, '/'.join(os.path.dirname(os.path.realpath(sys.argv[0])).split('/')[:-2] + ['bin']))
        sys.path.insert(0, os.path.dirname(os.path.realpath(sys.argv[0])))
        from poremaplib.poretypes import Poretypes
        from poremaplib.poreerr import Poreerr
        from poremaplib.porelog import Porelog
        from poremaplib.poreargs import Poreargs

    _profilepath = None
    try:
        for i in range(1, len(sys.argv)):
            if sys.argv[i].startswith('--profilepath'):
                if '=' in sys.argv[i]:
                    _profilepath = os.path.realpath(os.path.expandvars('='.join(sys.argv[i].split('=')[1:])))
                else:
                    _profilepath = os.path.realpath(os.path.expandvars(sys.argv[sys.arvg[i+1]]))
    except:
        if os.path.exists(_bindir):
            _profilepath = os.path.join(_bindir, 'ont.profile')
        else:
            _profilepath = 'ont.profile'
    if not _profilepath:
        _profilepath = os.path.join(_bindir, 'ont.profile')
    if os.path.exists(_profilepath):
        with open (_profilepath, 'r') as in_fp:
            for line in in_fp:
                if not line.startswith('export'):
                    continue
                info = line.strip().split(' ')[1].split('=')
                if len(info) == 2:
                    var, val = info
                    print 'Info: Setting {0}={1}'.format(var, val)
                    os.environ[var] = val
else:
    sys.path.insert(0, '/'.join(os.path.dirname(os.path.realpath(sys.argv[0])).split('/')[:-2] + ['bin']))
    sys.path.insert(0, os.path.dirname(os.path.realpath(sys.argv[0])))
    from poremaplib.poretypes import Poretypes
    from poremaplib.poreerr import Poreerr
    from poremaplib.porelog import Porelog
    from poremaplib.poreargs import Poreargs

# ============================================================================ #
# Global variables                                                             #
# ============================================================================ #

_version = '0.1.1'
#_progdesc = 'Alignment Statistics Program version {0}\n\nPrint read and run statistics for a mapped BAM file.'.format(_version)
_progdesc = 'poremapstats program v{0}'.format(_version)

#_ErrorSysCall = 0
#_ErrorDirCreate = 13
#_ErrorInvalidData = 18
#_ErrorFileMissing = 23
#_ErrorOutfileExists = 31

#_linewidth = 100

_C = {}		# Dictionary of isnonrandomaln classifications
_A = {}		# Dictionary of alignment stats
_E = {}		# Dictionary of read stats
_R = {}		# Dictionary of run stats

_ref = {}	# Dictionary of reference contigs
_targetrefidL = []
_controlrefidL = []

# ============================================================================ #
# Column information                                                           #
# ============================================================================ #

def Print_Column_Information(progname):
    'Print tab-separated text file containing column information.'

    print '# {0} version {1}'.format(progname, _version)
    print '# initstats.txt and alignstats.txt'
    print 'column\tdescription'
    for elt in _PT.alignstatsdesc:
        field, desc = elt
        print '{0}\t{1}'.format(field, desc)
    print '# readstats.txt'
    print 'column\tdescription'
    for elt in _PT.readstatsdesc:
        field, desc = elt
        print '{0}\t{1}'.format(field, desc)
    print '# runstats.txt'
    print 'column\tdescription'
    for elt in _PT.runstatsdesc:
        field, desc = elt
        print '{0}\t{1}'.format(field, desc)

# ============================================================================ #
# Program usage                                                                #
# ============================================================================ #

def Initialise():
    'Parse and validate command-line arguments, set up global variables.'

    # Set up globals
    global _progdir, _progname
    _progdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    _progname = os.path.basename(os.path.realpath(sys.argv[0]))
    progexamples = []

    global _PA, _PT, _PL, _PE
    _progname = os.path.basename(os.path.realpath(sys.argv[0]))
    progexamples = []
    _PT = Poretypes()
    _PE = Poreerr()
    _PL = Porelog(progdoc=__doc__, bindir=_bindir, ini_version=_version)
    _PA = Poreargs(progdoc=__doc__, progexamples=progexamples, bindir=_bindir)

    # Print column info and exit
    if 'printinfo' in ' '.join(sys.argv):
        Print_Column_Information(_progname)
        sys.exit(0)

    # Argument parsing
    _PA.arg_addstd( '--bindir', required=True )
    _PA.arg_addstd( '--profilepath', required=False )
    _PA.arg_addstd( '--runid', required=True )
    _PA.arg_addstd( '--readtype', required=True )
    _PA.arg_addstd( '--readclass', required=True )
    _PA.arg_addstd( '--datatype', required=True )
    _PA.arg_addstd( '--mapprog', required=True )
    _PA.arg_addstd( '--mapparams', required=True )
    _PA.arg_addstd( '--alignclasspath', required=False )
    _PA.arg_addstd( '--readsbam', required=False )
    _PA.arg_addstd( '--targetrefpath', required=True )
    _PA.arg_addstd( '--controlrefpath', required=True )
    _PA.arg_addstd( '--outdir', required=True )
    _PA.arg_addstd( '--outprefix', required=True )
    _PA.arg_addstd( '--savealignments', required=False, default=False )
    _PA.arg_addstd( '--fastalinewidth', required=False )
    _PA.arg_addstd( '--printinfo', required=False, default=False )
    _PA.arg_addstd( '--overwrite', required=False, default=False )
    _PA.arg_addstd( '--verbose', required=False, default=False )
    _PA.arg_addstd( '--logindentwidth', required=False )
    _PA.arg_parseall()

#    parser = argparse.ArgumentParser(description=_progdesc,
#        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument('--bindir', metavar='str', dest='bindir',
#        type=str, default=None, help='Absolute path to this program and all auxiliary Python files', required=True)
#    parser.add_argument('--runid', metavar='str', dest='runid',
#        type=str, default=None, help='Run identifier, must be unique across all ONT MinION runs', required=True)
#    parser.add_argument('--readtype', metavar='str', dest='readtype', choices=['2d', 'temp', 'comp', 'mixed', 'unknown'],
#        type=str, default='2D', help='Type of ONT reads: 2d, temp, comp, mixed, unknown', required=True)
#    parser.add_argument('--readclass', metavar='str', dest='readclass', choices=['all', 'pass', 'fail'],
#        type=str, default='all', help='Class of ONT reads: all, pass, fail', required=True)
#    parser.add_argument('--datatype', metavar='str', dest='datatype', choices=['minion', 'rand', 'mixed', 'unknown'],
#        type=str, default='minion', help='Type of run data: minion, random, mixed, unknown', required=True)
#    parser.add_argument('--mapprog', metavar='str', dest='mapprog',
#        type=str, default=None, help='Name of the mapping program used, in lower case, surrounded by quotes if it contains spaces (e.g., bwa)', required=True)
#    parser.add_argument('--mapparams', metavar='str', dest='mapparams',
#        type=str, default=None, help='The arguments passed to --mapprog, surrounded by quotes (e.g., "mem -x ont2d -M")', required=True)
#    parser.add_argument('--alignclass', metavar='str', dest='alignclass',
#        type=str, default=None, help='If None, output the initstats.txt file. If an alignclass.txt file containing the isnonrandomaln information, output the alignstats.txt, readstats.txt and runstats.txt files as output', required=False )
#    parser.add_argument('--readsbam', metavar='str', dest='readsbam',
#        type=str, default=None, help='Path to uniquely named, single-end reads file in BAM format', required=True)
#    parser.add_argument('--targetrefpath', metavar='str', dest='targetrefpath',
#        type=str, default=None, help='Path to FASTA file for extracting reads of target genome', required=True)
#    parser.add_argument('--controlrefpath', metavar='str', dest='controlrefpath',
#        type=str, default=None, help='Path to FASTA file for extracting reads of control genome', required=True)

#    parser.add_argument('--outdir', metavar='str', dest='outdir',
#        type=str, default=None, help='Output directory', required=True)
#    parser.add_argument('--outprefix', metavar='str', dest='outprefix',
#        type=str, default=None, help='Output file prefix', required=True)
#    parser.add_argument('--savealignments', action='store_true', dest='savealignments',
#        default=False, help='Output the alignment for each read as a FASTA file', required=False)
#    parser.add_argument('--printinfo', action='store_true', dest='printinfo',
#        default=False, help='Print a description of each of the output file fields', required=False)
#    parser.add_argument('--overwrite', action='store_true', dest='overwrite',
#        default=False, help='Print warning and overwrite any output files that may already exist', required=False)
#    parser.add_argument('--verbose', action='store_true', dest='verbose',
#        default=False, help='Print verbose output', required=False)
#    _args = parser.parse_args()

#    if _args.verbose:
#        sys.stdout.write('Info: Initialise()\n')

  # Initialise the log output
    _PL.log_2out('Info: Started')
    _PL.log_2out('{0}'.format(_PL.log_runinfo()))
    _PL.log_2out('Info:\n\n{valuelines}'.format(
        valuelines=_PA.arg_valuelines(indentwidth=int(_PA.args.logindentwidth))))

    # Set up global variables for out paths
    global _outalignstatspath, _outreadstatspath, _outrunstatspath
    if _PA.args.alignclasspath is None or _PA.args.alignclasspath.lower() == 'none':
        _outalignstatspath = os.path.join(_PA.args.outdir, "{0}_initstats.txt".format(_PA.args.outprefix))
        _PA.args.alignclasspath = None
    else:
        _outalignstatspath = os.path.join(_PA.args.outdir, "{0}_alignstats.txt".format(_PA.args.outprefix))
    _outreadstatspath = os.path.join(_PA.args.outdir, "{0}_readstats.txt".format(_PA.args.outprefix))
    _outrunstatspath = os.path.join(_PA.args.outdir, "{0}_runstats.txt".format(_PA.args.outprefix))

    # Fix type of non-essential input paths
    if _PA.args.targetrefpath.lower() == 'none':
        _PA.args.targetrefpath = None
    if _PA.args.controlrefpath.lower() == 'none':
        _PA.args.controlrefpath = None

    # Check input and exit on error
    if not os.path.exists(os.path.expandvars(_PA.args.readsbam)):
        _PL.log_2err('Erro: --readsbam does not exist ({0})'.format(_PA.args.readsbam))
        Exit('ErrorFileMissing')
    if _PA.args.targetrefpath is None and _PA.args.controlrefpath is None:
        _PL.log_2err('Erro: At least one of the following must be specified (--targetrefpath and --controlrefpath)')
        Exit('ErrorFileMissing')
    if _PA.args.targetrefpath is not None and not os.path.exists(os.path.expandvars(_PA.args.targetrefpath)):
        _PL.log_2err('Erro: --targetrefpath path not exist ({0})'.format(_PA.args.targetrefpath))
        Exit('ErrorFileMissing')
    if _PA.args.controlrefpath is not None and not os.path.exists(os.path.expandvars(_PA.args.controlrefpath)):
        _PL.log_2err('Erro: --controlrefpath path not exist ({0})'.format(_PA.args.controlrefpath))
        Exit('ErrorFileMissing')
    if _PA.args.mapprog is None:
        _PL.log_2out('Warn: --mapprog has not been specified ({0})'.format(_PA.args.mapprog))
    if _PA.args.mapparams is None:
        _PL.log_2out('Warn: --mapparams has not been specified ({0})'.format(_PA.args.mapparams))
    if not os.path.exists(os.path.expandvars(_PA.args.outdir)):
        try:
            os.makedirs(os.path.expandvars(_PA.args.outdir))
        except:
            _PL.log_2out('Erro: Failed to create non-existent --outdir ({0})\n'.format(_PA.args.outdir))
            Exit('ErrorDirCreate')
        sys.stdout.write('Info: Non-existent --outdir has been created ({0})\n'.format(_PA.args.outdir))
    if os.path.exists(_outalignstatspath):
        if _PA.args.overwrite:
            sys.stdout.write('Warn: Output alignstats file will be overwritten ({0})\n'.format(_outalignstatspath))
        else:
            sys.stdout.write('Warn: Output alignstats file already exists and overwrite is False ({0})\n'.format(_outalignstatspath))
            Exit('ErrorOutfileExists')
    if os.path.exists(_outreadstatspath):
        if _PA.args.overwrite:
            sys.stdout.write('Warn: Output readstats file will be overwritten ({0})\n'.format(_outreadstatspath))
        else:
            sys.stdout.write('Warn: Output readstats file already exists and overwrite is False ({0})\n'.format(_outreadstatspath))
            Exit('ErrorOutfileExists')
    if os.path.exists(_outrunstatspath):
        if _PA.args.overwrite:
            sys.stdout.write('Warn: Output runstats file will be overwritten ({0})\n'.format(_outrunstatspath))
        else:
            sys.stdout.write('Warn: Output runstats file alruny exists and overwrite is False ({0})\n'.format(_outrunstatspath))
            Exit('ErrorOutfileExists')

def Exit(retkey):
    'Pring diagnostic messages and exit the program with a return value and error message from the poreerr class.'
    _PL.log_exit(_PE.err_retkey(retkey), _PE.err_text(retkey), print2out=True, print2dump=True, exitprogram=True)

# ============================================================================ #
# Cigar string manipulation                                                    #
# ============================================================================ #

def CigarOperation(alignmentColumn):
    'Return a cigar operation. Copied from cigarCategory in maf-convert.py'
    x, y = alignmentColumn
    if x == "-":
        if y == "-": return "P"
        else: return "I"
    else:
        if y == "-": return "D"
        else: return "M"

def CigarOperationWithCorrectBases(alignmentColumn):
    'Return a pseudo-cigar operation with C denoting same-as-ref matches and M for other matches.'
    x, y = alignmentColumn
    if x == "-":
        if y == "-": return "P"
        else: return "I"
    else:
        if y == "-": return "D"
        else:
            if x == y:
                return "C"
            else:
                return "M"

def CigarParts(beg, alignmentcolumns, end):
    '''
    Split a list of column tuples into a list of cigar parts of format "nC"
    where n is a number and C is a cigar operation. Copied from maf-convert.py.
    '''
    if beg:
        yield str(beg) + "S"
    # (doesn't handle translated alignments)
    for k, v in itertools.groupby(alignmentcolumns, CigarOperation):
        yield str(sum(1 for _ in v)) + k
    if end:
        yield str(end) + "S"

def CigarPartsWithCorrectBases(beg, alignmentcolumns, end):
    '''
    Split a list of column tuples into a list of cigar parts of format "nC"
    where n is a number and C is a cigar operation. Copied from maf-convert.py.
    '''
    if beg:
        yield str(beg) + "S"
    # (doesn't handle translated alignments)
    for k, v in itertools.groupby(alignmentcolumns, CigarOperationWithCorrectBases):
        yield str(sum(1 for _ in v)) + k
    if end:
        yield str(end) + "S"

def Cigar_StoL(cigar):
    '''
    Transform a cigar string into a list of the nC components.
    Note that the 'C' cigar operation is my own that denotes "same as reference"
    so that I can more easily count stretches of matches that do not contain SNPs.
    '''
    replacements = { 'M':'M ', 'I':'I ', 'D':'D ', 'H':'H ', 'S':'S ', 'C':'C ', 'P':'P ' }
    cigarL = "".join([replacements.get(c,c) for c in cigar]).strip().split(' ')
    return cigarL

def Cigar_LtoAlignment(readid, refseq, refstartpos, readseq, cigarL):
    'Use a cigar string to reconstruct the alignment between the reference and read sequences.'
    refaln = []
    readaln = []
    refidx = refstartpos-1
    readidx = 0
    for e in cigarL:
        n = int(e[0:-1])
        C = e[-1]
        if e.endswith('H'):
            pass
        elif e.endswith('M'):
            refaln += list(refseq[refidx:refidx+n])
            readaln += list(readseq[readidx:readidx+n])
            refidx += n
            readidx += n
        elif e.endswith('I'):
            refaln += ['-']*n
            readaln += list(readseq[readidx:readidx+n])
            readidx += n
        elif e.endswith('D'):
            refaln += list(refseq[refidx:refidx+n])
            readaln += ['-']*n
            refidx += n
        elif e.endswith('S'):
            refaln += ['-']*n
            readaln += ['-']*n
            readidx += n
        else:
            _PL.log_2out('Erro: Unrecognised cigar instruction character ({0}:{1})\n'.format(readid, e))
            Exit('ErrorInvalidData')
    return refaln, readaln

def Cigar_LtoAlignmentPrintable(readid, refseq, refstartpos, readseq, cigarL):
    'Use a cigar string to reconstruct the alignment between the reference and read sequences.'
    refaln = []
    readaln = []
    refidx = refstartpos-1
    readidx = 0
    for e in cigarL:
        n = int(e[0:-1])
        C = e[-1]
        if e.endswith('H'):
            pass
        elif e.endswith('M'):
            refaln += list(refseq[refidx:refidx+n])
            readaln += list(readseq[readidx:readidx+n])
            refidx += n
            readidx += n
        elif e.endswith('I'):
            refaln += ['-']*n
            readaln += list(readseq[readidx:readidx+n])
            readidx += n
        elif e.endswith('D'):
            refaln += list(refseq[refidx:refidx+n])
            readaln += ['-']*n
            refidx += n
        elif e.endswith('S'):
            refaln += ['-']*n
            readaln += list(readseq[readidx:readidx+n])
            readidx += n
        else:
            _PL.log_2out('Erro: Unrecognised cigar instruction character ({0}:{1})\n'.format(readid, e))
            Exit('ErrorInvalidData')
    return refaln, readaln

# ============================================================================ #
# Reference                                                                    #
# ============================================================================ #

def Reference_Parse():
    'Read the reference contigs into a dictionary.'
    global _ref, _targetrefidL, _controlrefidL
    _ref = {}
    _targetrefidL = []
    _controlrefidL = []
    if _PA.args.targetrefpath is not None:
        for record in SeqIO.parse(_PA.args.targetrefpath, 'fasta'):
            _ref[record.id] = record.seq
            _targetrefidL.append(record.id)
    if _PA.args.controlrefpath is not None:
        for record in SeqIO.parse(_PA.args.controlrefpath, 'fasta'):
            _ref[record.id] = record.seq
            _controlrefidL.append(record.id)

# ============================================================================ #
# Alignment statistics                                                         #
# ============================================================================ #

def Parse_FastqReads():
    'Read the contents of the reads files in FASTQ format into flobal variable _F.'
    if _PA.args.verbose:
        sys.stdout.write('Info: ReadFastqReads()\n')

    global _F
    _F = {}
    in_fp = open(os.path.expandvars(_PA.args.readsfastq), 'r')
    if not in_fp:
        _PL.log_2out('Error: failed to open readsfastq file ({0})'.format(_PA.args.readsfastq))
        Exit('ErrorFileMissing')
    for record in SeqIO.parse(in_fp, 'fastq'):
        if not _F.has_key(record.id):
            _F[record.id] = [str(record.seq), ''.join([chr(x+33) for x in record.letter_annotations.values()[0]])]
        else:
            _PL.log_2out('Warn: more than one record with same id - overwriting previous bq value ({0})\n'.format(record.id))
    in_fp.close()

def Convert_BamToSam(bampath, sampath):
    #cmd = 'samtools view {bampath} > {sampath}'.format(bampath=bampath, sampath=sampath)
    cmd = 'samtools view {bampath}'.format(bampath=bampath)
    #print cmd
    rh = sp.Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
    ro, re = rh.communicate()
    rc = rh.returncode
    #os.system(cmd)
    #if not os.path.exists(sampath) or os.stat(sampath).st_size == 0:
    if rc != 0:
        _PL.log_2err('Erro: Failed to convert BAM to SAM before ({0})'.format(cmd))
        Exit('ErrorSysCall')
    with open(sampath, 'w') as out_fp:
        out_fp.write('{0}\n'.format(ro))
    return 0

def Parse_SamLine(line):
    L = line.strip().split('\t')
    L[1] = int(L[1])
    L[3] = int(L[3])
    L[4] = int(L[4])
    L[7] = int(L[7])
    L[8] = int(L[8])
    return L

def AlignStats_Compute_One(samlineL, savealignments):
    'Return a list containing the mapping statistics for this read.'

    result = {}		# Dictionary keyed by all the variables in alignstatsheader column names

    # Parse the SAM record    
    readid, samflag, refcontigid, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = samlineL[0:11]

    # Set output values that are relevant for mapped and unmapped reads
    result['runid'] = _PA.args.runid
    result['readid'] = readid
    result['readtype'] = _PA.args.readtype
    result['readclass'] = _PA.args.readclass
    result['datatype'] = _PA.args.datatype
    result['mapprog'] = _PA.args.mapprog
    result['mapparams'] = _PA.args.mapparams
    result['samflag'] = samflag
    result['refcontigid'] = refcontigid
    result['refcontigpos1'] = pos
    result['mapq'] = mapq
    key = (_PA.args.runid, readid, _PA.args.readtype, _PA.args.readclass, _PA.args.datatype, _PA.args.mapprog, _PA.args.mapparams, samflag, refcontigid, pos, mapq)

    # Set values for unmapped reads
    if refcontigid == '*' or mapq == 0:
        if len(cigar) and cigar != '*':
            cigarL = Cigar_StoL(cigar)
            if len(cigarL):
                result['readbp'] = sum([int(x[0:-1]) for x in cigarL if x[-1]!='D'])
            else:
                result['readbp'] = len(seq)
        else:
            result['readbp'] = len(seq)
        result['isprimaryaln'] = 0
        result['isnonrandomaln'] = 1
        result['alnlen'] = 0
        result['alnrefbp'], result['alnrefstartpos1'], result['alnrefendpos1'] = [0]*3
        result['alnreadbp'], result['alnreadstartpos1'], result['alnreadendpos1'] = [0]*3
        result['numC'], result['numM'], result['numI'], result['numD'] = [0]*4
        result['freqCperbp'], result['freqMperbp'], result['freqIperbp'], result['freqDperbp'] = [0.0]*4
        result['lenS'], result['lenC'], result['lenM'], result['lenI'], result['lenD'] = [0]*5
        result['pctS'], result['pctC'], result['pctM'], result['pctI'], result['pctD'] = [0.0]*5
        result['meanrunlenC'], result['meanrunlenM'], result['meanrunlenI'], result['meanrunlenD'] = [0.0]*4
        result['meanbq'], result['meanbqC'], result['meanbqM'], result['meanbqI'] = [0.0]*4
 
    # Set values for mapped reads
    else:
        cigarL = Cigar_StoL(cigar)
        # Note: len(refaln) == len(readaln) == len(elts in entire cigarL including leading/trailing soft/hard clipped bases)
        refaln, readaln = Cigar_LtoAlignment(readid, _ref[refcontigid], pos, seq, cigarL)
        cigarC = ''.join(CigarPartsWithCorrectBases(0, zip(refaln, readaln), 0))
        cigarCL = Cigar_StoL(cigarC)
        cigarLaln = [cigarL[i] for i in range(0 if cigarL[0][-1] not in 'HS' else 1, len(cigarL) if cigarL[-1][-1] not in 'HS' else len(cigarL)-1)]

        result['readbp'] = sum([int(x[0:-1]) for x in cigarL if x[-1]!='D'])
        if readid == '1253_1_channel_62_read_51_2D':
            pass
        result['isprimaryaln'] = int((not (samflag & 128)) and (not (samflag & 256)) and (not (samflag & 2048)))
        result['isnonrandomaln'] = -1 if not _C.has_key(key) else _C[key]
        result['alnlen'] = sum([int(x[0:-1]) for x in cigarLaln])

        result['alnrefbp'] = sum([int(x[0:-1]) for x in cigarLaln if x[-1] in 'MD'])
        alnrefbp_check = len([x for x in refaln if x != '-'])
        if result['alnrefbp'] != alnrefbp_check:
            sys.stdout.write('Warn: readid={0}, alnrefbp={1} != alnrefbp_check={2}\n'.format(
                readid, result['alnrefbp'], alnrefbp_check))
        result['alnrefstartpos1'] = pos if cigarL[0][-1]=='H' else (pos+int(cigarL[0][0:-1]) if cigarL[0][-1]=='S' else pos)
        result['alnrefendpos1'] = result['alnrefstartpos1'] + result['alnrefbp'] - 1

        result['alnreadbp'] = sum([int(x[0:-1]) for x in cigarLaln if x[-1] in 'MI'])
        alnreadbp_check = len([x for x in readaln if x != '-'])
        if result['alnreadbp'] != alnreadbp_check:
            sys.stdout.write('Warn: readid={0}, alnreadbp={1} != alnreadbp_check={2}\n'.format(
                readid, result['alnreadbp'], alnreadbp_check))
        result['alnreadstartpos1'] = 1 if cigarL[0][-1]=='H' else (1+int(cigarL[0][0:-1]) if cigarL[0][-1]=='S' else 1)
        result['alnreadendpos1'] = result['alnreadstartpos1'] + result['alnreadbp'] -1
        if result['alnreadstartpos1'] > result['readbp']:
            sys.stdout.write('Warn: readid {0}: alnreadstartpos {1} > read length {2}\n'.format(
                result['readid'], result['alnreadstartpos1'], result['readbp']))
        if result['alnreadendpos1'] > result['readbp']:
            sys.stdout.write('Warn: readid {0}: alnreadendpos {1} > read length {2}\n'.format(
                result['readid'], result['alnreadendpos1'], result['readbp']))
        result['numC'] = sum([1 for x in cigarCL if x.endswith('C')])
        result['numM'] = sum([1 for x in cigarL if x.endswith('M')])
        result['numI'] = sum([1 for x in cigarL if x.endswith('I')])
        result['numD'] = sum([1 for x in cigarL if x.endswith('D')])
        result['freqCperbp'] = round(result['numC'] / float(result['alnlen']), 6) if result['alnlen'] else 0.0
        result['freqMperbp'] = round(result['numM'] / float(result['alnlen']), 6) if result['alnlen'] else 0.0
        result['freqIperbp'] = round(result['numI'] / float(result['alnlen']), 6) if result['alnlen'] else 0.0
        result['freqDperbp'] = round(result['numD'] / float(result['alnlen']), 6) if result['alnlen'] else 0.0

        result['lenS'] = sum([1 for x in itertools.izip(refaln, readaln) if x[0] != x[1] and x[0] != '-' and x[1] != '-'])
        result['lenC'] = sum([int(x[0:-1]) for x in cigarCL if x.endswith('C')])
        result['lenM'] = sum([int(x[0:-1]) for x in cigarL if x.endswith('M')])
        result['lenI'] = sum([int(x[0:-1]) for x in cigarL if x.endswith('I')])
        result['lenD'] = sum([int(x[0:-1]) for x in cigarL if x.endswith('D')])
        result['pctS'] = round(result['lenS'] / float(result['alnlen']) * 100.0, 6) if result['alnlen'] else 0.0
        result['pctC'] = round(result['lenC'] / float(result['alnlen']) * 100.0, 6) if result['alnlen'] else 0.0
        result['pctM'] = round(result['lenM'] / float(result['alnlen']) * 100.0, 6) if result['alnlen'] else 0.0
        result['pctI'] = round(result['lenI'] / float(result['alnlen']) * 100.0, 6) if result['alnlen'] else 0.0
        result['pctD'] = round(result['lenD'] / float(result['alnlen']) * 100.0, 6) if result['alnlen'] else 0.0

        result['meanrunlenC'] = round(result['lenC'] / float(result['numC']), 3) if result['numC'] else 0.0
        result['meanrunlenM'] = round(result['lenM'] / float(result['numM']), 3) if result['numM'] else 0.0
        result['meanrunlenI'] = round(result['lenI'] / float(result['numI']), 3) if result['numI'] else 0.0
        result['meanrunlenD'] = round(result['lenD'] / float(result['numD']), 3) if result['numD'] else 0.0

        result['meanbq'] = round(sum([ord(x)-33 for x in qual]) / float(len(qual)), 3)

        mask = [b=='C' for b in ''.join([x[-1]*int(x[0:-1]) for x in cigarCL if x[-1] in 'CMISP'])]
        bq = [ord(a)-33 for (a, isgood) in zip([x for x in qual], mask) if isgood]
        if len(mask) != len(qual):
            _PL.log_2out('Warn: readid={0}, len(Cmask)={1} != len(qual)={2}\n'.format(
                readid, len(mask), len(qual)))
        result['meanbqC'] = round(sum(bq) / float(len(bq)), 3) if len(bq) else 0.0

        mask = [b=='M' for b in ''.join([x[-1]*int(x[0:-1]) for x in cigarL if x[-1] in 'MIS'])]
        bq = [ord(a)-33 for (a, isgood) in zip([x for x in qual], mask) if isgood]
        if len(mask) != len(qual):
            _PL.log_2out('Warn: readid={0}, len(Mmask)={1} != len(qual)={2}\n'.format(
                readid, len(mask), len(qual)))
        result['meanbqM'] = round(sum(bq) / float(len(bq)), 3) if len(bq) else 0.0

        mask = [b=='I' for b in ''.join([x[-1]*int(x[0:-1]) for x in cigarL if x[-1] in 'MIS'])]
        bq = [ord(a)-33 for (a, isgood) in zip([x for x in qual], mask) if isgood]
        if len(mask) != len(qual):
            _PL.log_2out('Warn: readid={0}, len(Imask)={1} != len(qual)={2}\n'.format(
                readid, len(mask), len(qual)))
        result['meanbqI'] = round(sum(bq) / float(len(bq)), 3) if len(bq) else 0.0

        # Save the alignment to a fasta file for debugging purposes
        if savealignments:
            refalnp, readalnp = Cigar_LtoAlignmentPrintable(readid, _ref[refcontigid], pos, seq, cigarL)
            fasta_path = os.path.join(_PA.args.outdir,
                '{outprefix}_{readid}_{samflag}_{refcontigid}_{pos}_{mapq}_aln.fasta'.format(
                outprefix=_PA.args.outprefix, readid=readid, samflag=samflag, refcontigid=refcontigid,
                pos=pos, mapq=mapq, datatype=_PA.args.datatype))
            with open(fasta_path, 'w') as out_fp:
                seq = ''.join(refalnp)
                out_fp.write('>{header}\n{seq}\n'.format(
                    header='ref',
                    seq="\n".join([seq[i:i+int(_PA.args.fastalinewidth)] for i in range(0, len(seq), int(_PA.args.fastalinewidth))])))
                seq = ''.join(readalnp)
                out_fp.write('>{header}\n{seq}\n'.format(
                    header='read',
                    seq="\n".join([seq[i:i+int(_PA.args.fastalinewidth)] for i in range(0, len(seq), int(_PA.args.fastalinewidth))])))

    return result

def AlignStats_Compute(savealignments):
    '''
    For every read in the BAM file, use the read sequence, CIGAR string and the corresponding
    segment of the refcontig sequence to infer the statistics for this read. These are stored
    in the global _A dictionary.
    '''
    tmpsampath = os.path.join(os.path.dirname(_PA.args.readsbam), os.path.basename(_PA.args.readsbam).replace('.bam', '.tmp.sam'))
    Convert_BamToSam(_PA.args.readsbam, tmpsampath)

    # Read in the --alignclasspath file and use the isnonrandomaln values for each alignment, otherwise set them to the default values.
    global _C
    _C = {}
    if _PA.args.alignclasspath is not None and os.path.exists(_PA.args.alignclasspath):
        data = np.loadtxt(_PA.args.alignclasspath, dtype=_PT.alignclassdtype, delimiter='\t', skiprows=1)
        for rowT in data:
            rowL = list(rowT)
	    key = tuple(rowL[0:11])
            val = rowL[11]
            _C[key] = val

    global _A
    _A = {}
    with open(tmpsampath, 'r') as in_fp:
        cnt = 0
        for line in in_fp:
            cnt += 1
            if len(line.strip()) == 0 or line.startswith('@'):
                continue
            samlineL = Parse_SamLine(line)
            readid, samflag, refcontigid, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = samlineL[0:11]
            key = (readid, samflag, refcontigid, pos, mapq)
            _A[key] = AlignStats_Compute_One(samlineL, savealignments)

    if os.path.exists(tmpsampath):
        os.remove(tmpsampath)

def AlignStats_Save():
    'Save the pre-formatted alignment statistics in _A to a file.'
    keyL = _A.keys()
    keyL.sort()
    with open(_outalignstatspath, 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join([var for var in _PT.alignstatsheader])))
        for k in keyL:
            out_fp.write('{0}\n'.format('\t'.join([str(_A[k][var]) for var in _PT.alignstatsheader])))

def AlignStats_Print(savealignments):
    AlignStats_Compute(savealignments)
    AlignStats_Save()

# ============================================================================ #
# Read statistics                                                              #
# ============================================================================ #

def NumAlignedReadBases(alignstats):
    ranges = [[x['alnreadstartpos1'], x['alnreadendpos1']] for x in alignstats]
    basecoordsinranges = []
    for elt in ranges:
        basecoordsinranges += range(elt[0], elt[1]+1)
    numunique = len(set(basecoordsinranges))
    return numunique

def ReadStats_Compute_AlnSubStats(readid, totalreadbp, alnstats):
    '''
    Compute the 16 statistics for the alnstats entries for this read.
    The calling function must choose a subset of these if necessary.
    Decided to treat the case where there is only one alignment in alnstats
    as a special case to save time on doing all those list comprehensions.
    '''
    if len(alnstats) == 1:
        A = alnstats[0]
        alnlen = A['alnlen'] if A['isnonrandomaln']!=0 else 0
        readbp = A['alnreadbp'] if A['isnonrandomaln']!=0 else 0
        refbp = A['alnrefbp'] if A['isnonrandomaln']!=0 else 0
        meanfreqC = A['freqCperbp']
        meanfreqM = A['freqMperbp']
        meanfreqI = A['freqIperbp']
        meanfreqD = A['freqDperbp']
        pctS = A['pctS']
        pctC = A['pctC']
        pctM = A['pctM']
        pctI = A['pctI']
        pctD = A['pctD']
        meanrunlenC = A['meanrunlenC']
        meanrunlenM = A['meanrunlenM']
        meanrunlenI = A['meanrunlenI']
        meanrunlenD = A['meanrunlenD']
        meanbq = A['meanbq']
        meanbqC = A['meanbqC']
        meanbqM = A['meanbqM']
        meanbqI = A['meanbqI']
        targetaligncnt = 1 if A['refcontigid'] in _targetrefidL and A['isnonrandomaln']!=0 else 0
        targetalignbp = A['alnreadbp'] if A['refcontigid'] in _targetrefidL and A['isnonrandomaln']!=0 else 0
        controlaligncnt = 1 if A['refcontigid'] in _controlrefidL and A['isnonrandomaln']!=0 else 0
        controlalignbp = A['alnreadbp'] if A['refcontigid'] in _controlrefidL and A['isnonrandomaln']!=0 else 0

        randalncnt = 1 if A['isnonrandomaln']==0 else 0
        randalncntpct = 100.0 if A['isnonrandomaln']==0 else 0.0
        randalnbp = A['alnreadbp'] if A['isnonrandomaln']==0 else 0
        randalnbppct = 100.0 if A['isnonrandomaln']==0 else 0.0
        randalnmeanbp = A['meanbq'] if A['isnonrandomaln']==0 else 0.0

    else:
        alnnum = len([x for x in alnstats if x['isnonrandomaln']!=0])
        alnlen = sum([x['alnlen'] for x in alnstats if x['isnonrandomaln']!=0])
        readbp = NumAlignedReadBases([x for x in alnstats if x['isnonrandomaln']!=0])
        if readbp > totalreadbp:
            sys.stdout.write('Warn: readid {0}: sum of aligned read bases {1} > read length {2}\n'.format(readid, readbp, totalreadbp))
        refbp = sum([x['alnrefbp'] for x in alnstats if x['isnonrandomaln']!=0])
        meanfreqC = round(sum([x['numC'] for x in alnstats if x['isnonrandomaln']!=0 and x['numC']>0]) / float(alnlen), 6) if alnlen else 0.0
        meanfreqM = round(sum([x['numM'] for x in alnstats if x['isnonrandomaln']!=0 and x['numM']>0]) / float(alnlen), 6) if alnlen else 0.0
        meanfreqI = round(sum([x['numI'] for x in alnstats if x['isnonrandomaln']!=0 and x['numI']>0]) / float(alnlen), 6) if alnlen else 0.0
        meanfreqD = round(sum([x['numD'] for x in alnstats if x['isnonrandomaln']!=0 and x['numD']>0]) / float(alnlen), 6) if alnlen else 0.0
        pctS = round(sum([x['lenS'] for x in alnstats if x['isnonrandomaln']!=0 and x['lenS']>0]) / float(alnlen) * 100.0, 6) if alnlen else 0.0
        pctC = round(sum([x['lenC'] for x in alnstats if x['isnonrandomaln']!=0 and x['lenC']>0]) / float(alnlen) * 100.0, 6) if alnlen else 0.0
        pctM = round(sum([x['lenM'] for x in alnstats if x['isnonrandomaln']!=0 and x['lenM']>0]) / float(alnlen) * 100.0, 6) if alnlen else 0.0
        pctI = round(sum([x['lenI'] for x in alnstats if x['isnonrandomaln']!=0 and x['lenI']>0]) / float(alnlen) * 100.0, 6) if alnlen else 0.0
        pctD = round(sum([x['lenD'] for x in alnstats if x['isnonrandomaln']!=0 and x['lenD']>0]) / float(alnlen) * 100.0, 6) if alnlen else 0.0
        meanrunlenC = round(sum([x['meanrunlenC'] for x in alnstats if x['isnonrandomaln']!=0 and x['meanrunlenC']>0]) / \
            float(sum([x['numC'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['numC'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        meanrunlenM = round(sum([x['meanrunlenM'] for x in alnstats if x['isnonrandomaln']!=0 and x['meanrunlenM']>0]) / \
            float(sum([x['numM'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['numM'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        meanrunlenI = round(sum([x['meanrunlenI'] for x in alnstats if x['isnonrandomaln']!=0 and x['meanrunlenI']>0]) / \
            float(sum([x['numI'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['numI'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        meanrunlenD = round(sum([x['meanrunlenD'] for x in alnstats if x['isnonrandomaln']!=0 and x['meanrunlenD']>0]) / \
            float(sum([x['numD'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['numD'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        meanbq = round(sum([x['meanbq']*x['readbp'] for x in alnstats if x['isnonrandomaln']!=0]) / \
            float(sum([x['readbp'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['readbp'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        meanbqC = round(sum([x['meanbqC']*x['numC'] for x in alnstats if x['isnonrandomaln']!=0]) / \
            float(sum([x['numC'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['numC'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        meanbqM = round(sum([x['meanbqM']*x['numM'] for x in alnstats if x['isnonrandomaln']!=0]) / \
            float(sum([x['numM'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['numM'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        meanbqI = round(sum([x['meanbqI']*x['numI'] for x in alnstats if x['isnonrandomaln']!=0]) / \
            float(sum([x['numI'] for x in alnstats if x['isnonrandomaln']!=0])), 3) \
            if sum([x['numI'] for x in alnstats if x['isnonrandomaln']!=0]) else 0.0
        targetaligncnt = sum([1 for x in alnstats if x['refcontigid'] in _targetrefidL and x['isnonrandomaln']!=0])
        targetalignbp = sum([x['alnreadbp'] for x in alnstats if  x['refcontigid'] in _targetrefidL and x['isnonrandomaln']!=0])
        controlaligncnt = sum([1 for x in alnstats if x['refcontigid'] in _controlrefidL and x['isnonrandomaln']!=0])
        controlalignbp = sum([x['alnreadbp'] for x in alnstats if  x['refcontigid'] in _controlrefidL and x['isnonrandomaln']!=0])

        randalncnt = sum([1 for x in alnstats if x['isnonrandomaln']==0])
        randalncntpct = round(randalncnt / float(alnnum) * 100.0, 3) if alnnum else 0.0
        randalnbp = sum([x['alnreadbp'] for x in alnstats if x['isnonrandomaln']==0])
        randalnbppct = round(randalnbp / float(alnlen) * 100.0, 3) if alnlen else 0.0
        randalnmeanbp = round(sum([x['meanbq']*x['readbp'] for x in alnstats if x['isnonrandomaln']==0]) / \
            float(sum([x['readbp'] for x in alnstats if x['isnonrandomaln']==0])), 3) \
            if sum([x['readbp'] for x in alnstats if x['isnonrandomaln']==0]) else 0.0

    result =  [
        alnlen, readbp, refbp,
        meanfreqC, meanfreqM, meanfreqI, meanfreqD,
        pctS, pctC, pctM, pctI, pctD,
        meanrunlenC, meanrunlenM, meanrunlenI, meanrunlenD,
        meanbq, meanbqC, meanbqM, meanbqI,
        targetaligncnt, targetalignbp, controlaligncnt, controlalignbp,
        randalncnt, randalncntpct, randalnbp, randalnbppct, randalnmeanbp
    ]
    return result

def ReadStats_Compute_One(readid, alnstats):
    'Return dictionary of statistics keyed on variable names in readstatsheader.'

    if not alnstats or not len(alnstats):
        _PL.log_2out('Erro: No alignment statistics for readid ({0})\n'.format(readid))
        Exit('ErrorInvalidData')

    result = {}
    result['runid'] = _PA.args.runid
    result['readid'] = readid
    result['readtype'] = _PA.args.readtype
    result['readclass'] = _PA.args.readclass
    result['datatype'] = _PA.args.datatype
    result['mapprog'] = _PA.args.mapprog
    result['mapparams'] = _PA.args.mapparams

    result['readbp'] = alnstats[0]['readbp']
    result['ismapped'] = int(len([x for x in alnstats if x['isnonrandomaln']!=0 and x['refcontigid']!='*'])>0)    # At least 1 non-random (or unknown), named refcontig
    result['numalignments'] = len([x for x in alnstats if x['isnonrandomaln']!=0 and x['refcontigid']!='*'])
    result['numrefcontigs'] = len(set([x['refcontigid'] for x in alnstats if x['isnonrandomaln']!=0 and x['refcontigid']!='*'])) # Num named refcontigs involved in non-random (or unknown) alignments

    result['allalnlen'], result['allalnreadbp'], result['allalnrefbp'], \
    result['allalnfreqCperbp'], result['allalnfreqMperbp'], result['allalnfreqIperbp'], result['allalnfreqDperbp'], \
    result['allalnpctS'], result['allalnpctC'], result['allalnpctM'], result['allalnpctI'], result['allalnpctD'], \
    result['allalnmeanrunlenC'], result['allalnmeanrunlenM'], result['allalnmeanrunlenI'], result['allalnmeanrunlenD'], \
    result['allalnmeanbq'], result['allalnmeanbqC'], result['allalnmeanbqM'], result['allalnmeanbqI'], \
    result['alltargetaligncnt'], result['alltargetalignbp'], result['allcontrolaligncnt'], result['allcontrolalignbp'], \
    result['allrandalncnt'], result['allrandalncntpct'], result['allrandalnbp'], result['allrandalnbppct'], result['allrandalnmeanbp'] \
        = ReadStats_Compute_AlnSubStats(readid, result['readbp'], alnstats)

    prialnstats = [x for x in alnstats if x['isprimaryaln']==1]
    result['prialnlen'], result['prialnreadbp'], result['prialnrefbp'], \
    result['prialnfreqCperbp'], result['prialnfreqMperbp'], result['prialnfreqIperbp'], result['prialnfreqDperbp'], \
    result['prialnpctS'], result['prialnpctC'], result['prialnpctM'], result['prialnpctI'], result['prialnpctD'], \
    result['prialnmeanrunlenC'], result['prialnmeanrunlenM'], result['prialnmeanrunlenI'], result['prialnmeanrunlenD'], \
    result['prialnmeanbq'], result['prialnmeanbqC'], result['prialnmeanbqM'], result['prialnmeanbqI'], \
    result['pritargetaligncnt'], result['pritargetalignbp'], result['pricontrolaligncnt'], result['pricontrolalignbp'], \
    result['prirandalncnt'], result['prirandalncntpct'], result['prirandalnbp'], result['prirandalnbppct'], result['prirandalnmeanbp'] \
        = ReadStats_Compute_AlnSubStats(readid, result['readbp'], prialnstats)
    result['prialncnt'] = len([x for x in alnstats if x['refcontigid']!='*' and x['mapq']!=0 and x['isprimaryaln']==1 and x['isnonrandomaln']!=0])

    return result

def ReadStats_Compute():
    'Compute the read statistics from the alignments in _A, and save in _E.'

    readidL = list(set([k[0] for k in _A.keys()]))
    readidL.sort()

    global _E
    _E = {}
    for readid in readidL:
        alignments = [_A[k] for k in _A.keys() if k[0] == readid]
        key = (_PA.args.runid, readid, _PA.args.readtype, _PA.args.readclass, _PA.args.datatype, _PA.args.mapprog, _PA.args.mapparams)
        _E[key] = ReadStats_Compute_One(readid, alignments)
        pass
    pass

def ReadStats_Save():
    'Save the pre-formatted read statistics in _E to a file.'
    keyL = _E.keys()
    keyL.sort()
    with open(_outreadstatspath, 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join([var for var in _PT.readstatsheader])))
        for k in keyL:
            out_fp.write('{0}\n'.format('\t'.join([str(_E[k][var]) for var in _PT.readstatsheader])))

def ReadStats_Print():
    ReadStats_Compute()
    ReadStats_Save()

# ============================================================================ #
# Run statistics                                                               #
# ============================================================================ #

def RunStats_Compute_FromSubset(readstats, totalreadcnt, totalreadbp, readfields):
    'Return the run statistics for the list of readstats provided.'

    if readfields not in ['all', 'pri']:
        _PL.log_2out('Erro: Invalid readfields parameter passed to RunStats_Compute_FromSubset ({0})\n'.format(readfields))
        Exit('ErrorInvalidData')

    n = len(readstats)
    nf = float(n)

    mapreadcnt = len([x for x in readstats if x['ismapped']==1])		# For mapped reads only
    mapreadcntpct = round(mapreadcnt / float(totalreadcnt) * 100.0, 3) if totalreadcnt else 0.0

    mapreadbp = sum([x['readbp'] for x in readstats if x['ismapped']==1])	# For mapped reads only
    mapreadbppct = round(mapreadbp / float(totalreadbp) * 100.0, 3) if totalreadbp else 0.0
    mapmeanreadbp = round(mapreadbp / float(mapreadcnt), 3) if mapreadcnt else 0.0

    if readfields == 'all':
        alncnt = sum([x['numalignments'] for x in readstats])
        alnbp = sum([x['allalnreadbp'] for x in readstats])
        alnbppct = round(alnbp / float(totalreadbp) * 100.0, 6) if totalreadbp else 0.0
        alnmeanbp = round(alnbp / nf, 3) if n else 0.0
        Spct = round(sum([x['allalnpctS'] for x in readstats]) / nf, 6) if n else 0.0
        Cpct = round(sum([x['allalnpctC'] for x in readstats]) / nf, 6) if n else 0.0
        Mpct = round(sum([x['allalnpctM'] for x in readstats]) / nf, 6) if n else 0.0
        Ipct = round(sum([x['allalnpctI'] for x in readstats]) / nf, 6) if n else 0.0
        Dpct = round(sum([x['allalnpctD'] for x in readstats]) / nf, 6) if n else 0.0

        Cmeandistbp = round(1 / (sum([x['allalnfreqCperbp'] for x in readstats]) / nf), 6) if (n and sum([x['allalnfreqCperbp'] for x in readstats])) else 0.0
        Mmeandistbp = round(1 / (sum([x['allalnfreqMperbp'] for x in readstats]) / nf), 6) if (n and sum([x['allalnfreqMperbp'] for x in readstats])) else 0.0
        Imeandistbp = round(1 / (sum([x['allalnfreqIperbp'] for x in readstats]) / nf), 6) if (n and sum([x['allalnfreqIperbp'] for x in readstats])) else 0.0
        Dmeandistbp = round(1 / (sum([x['allalnfreqDperbp'] for x in readstats]) / nf), 6) if (n and sum([x['allalnfreqDperbp'] for x in readstats])) else 0.0

        meanrunlenC = round(sum([x['allalnmeanrunlenC'] for x in readstats]) / nf, 3) if n else 0.0
        meanrunlenM = round(sum([x['allalnmeanrunlenM'] for x in readstats]) / nf, 3) if n else 0.0
        meanrunlenI = round(sum([x['allalnmeanrunlenI'] for x in readstats]) / nf, 3) if n else 0.0
        meanrunlenD = round(sum([x['allalnmeanrunlenD'] for x in readstats]) / nf, 3) if n else 0.0

        targetaligncnt = sum([x['alltargetaligncnt'] for x in readstats])
        targetalignbp = sum([x['alltargetalignbp'] for x in readstats])
        controlaligncnt = sum([x['allcontrolaligncnt'] for x in readstats])
        controlalignbp = sum([x['allcontrolalignbp'] for x in readstats])

        randalncnt = sum([x['allrandalncnt'] for x in readstats])
        randalncntpct = randalncnt / float(mapreadcnt) * 100.0 if mapreadcnt else 0.0
        randalnbp = sum([x['allrandalnbp'] for x in readstats])
        randalnbppct = randalnbp / float(mapreadbp) * 100.0 if mapreadbp else 0.0
        randalnmeanbp = round(sum([x['allrandalnmeanbp'] for x in readstats]) / nf, 3) if n else 0.0

    elif readfields == 'pri':
        alncnt = sum([x['prialncnt'] for x in readstats])
        alnbp = sum([x['prialnreadbp'] for x in readstats if x['ismapped']==1])
        alnbppct = round(alnbp / float(totalreadbp) * 100.0, 6) if totalreadbp else 0.0
        alnmeanbp = round(alnbp / nf, 3) if n else 0.0
        Spct = round(sum([x['prialnpctS'] for x in readstats]) / nf, 6) if n else 0.0
        Cpct = round(sum([x['prialnpctC'] for x in readstats]) / nf, 6) if n else 0.0
        Mpct = round(sum([x['prialnpctM'] for x in readstats]) / nf, 6) if n else 0.0
        Ipct = round(sum([x['prialnpctI'] for x in readstats]) / nf, 6) if n else 0.0
        Dpct = round(sum([x['prialnpctD'] for x in readstats]) / nf, 6) if n else 0.0

        Cmeandistbp = round(1 / (sum([x['prialnfreqCperbp'] for x in readstats]) / nf), 6) if (n and sum([x['prialnfreqCperbp'] for x in readstats])) else 0.0
        Mmeandistbp = round(1 / (sum([x['prialnfreqMperbp'] for x in readstats]) / nf), 6) if (n and sum([x['prialnfreqMperbp'] for x in readstats])) else 0.0
        Imeandistbp = round(1 / (sum([x['prialnfreqIperbp'] for x in readstats]) / nf), 6) if (n and sum([x['prialnfreqIperbp'] for x in readstats])) else 0.0
        Dmeandistbp = round(1 / (sum([x['prialnfreqDperbp'] for x in readstats]) / nf), 6) if (n and sum([x['prialnfreqDperbp'] for x in readstats])) else 0.0

        meanrunlenC = round(sum([x['prialnmeanrunlenC'] for x in readstats]) / nf, 3) if n else 0.0
        meanrunlenM = round(sum([x['prialnmeanrunlenM'] for x in readstats]) / nf, 3) if n else 0.0
        meanrunlenI = round(sum([x['prialnmeanrunlenI'] for x in readstats]) / nf, 3) if n else 0.0
        meanrunlenD = round(sum([x['prialnmeanrunlenD'] for x in readstats]) / nf, 3) if n else 0.0

        targetaligncnt = sum([x['pritargetaligncnt'] for x in readstats])
        targetalignbp = sum([x['pritargetalignbp'] for x in readstats])
        controlaligncnt = sum([x['pricontrolaligncnt'] for x in readstats])
        controlalignbp = sum([x['pricontrolalignbp'] for x in readstats])

        randalncnt = sum([x['prirandalncnt'] for x in readstats])
        randalncntpct = randalncnt / float(mapreadcnt) * 100.0 if mapreadcnt else 0.0
        randalnbp = sum([x['prirandalnbp'] for x in readstats])
        randalnbppct = randalnbp / float(mapreadbp) * 100.0 if mapreadbp else 0.0
        randalnmeanbp = round(sum([x['prirandalnmeanbp'] for x in readstats]) / nf, 3) if n else 0.0

    result = [
        mapreadcnt, mapreadcntpct, mapreadbp, mapreadbppct, mapmeanreadbp,
        alncnt, alnbp, alnbppct, alnmeanbp,
        Spct, Cpct, Mpct, Ipct, Dpct,
        Cmeandistbp, Mmeandistbp, Imeandistbp, Dmeandistbp,
        meanrunlenC, meanrunlenM, meanrunlenI, meanrunlenD,
        targetaligncnt, targetalignbp, controlaligncnt, controlalignbp,
        randalncnt, randalncntpct, randalnbp, randalnbppct, randalnmeanbp
    ]
    return result

def RunStats_Compute_One(readstats, readfields, runstattype):
    'Compute the run statistics from the contents of the read statstics in _A.'

    result = {}
    result['runid'] = _PA.args.runid
    result['readtype'] = _PA.args.readtype
    result['readclass'] = _PA.args.readclass
    result['datatype'] = _PA.args.datatype
    result['mapprog'] = _PA.args.mapprog
    result['mapparams'] = _PA.args.mapparams
    result['runstattype'] = runstattype

    result['runreadcnt'] = len(_E.keys())
    result['runreadbp'] = sum([_E[k]['readbp'] for k in _E.keys()])
    result['runmeanreadbp'] = round(result['runreadbp'] / float(result['runreadcnt']), 3) if result['runreadcnt'] else 0.0

    result['umapreadcnt'] = len([x for x in readstats if x['ismapped']==0])
    result['umapreadcntpct'] = round(result['umapreadcnt'] / float(result['runreadcnt']) * 100.0, 3) if result['runreadcnt'] else 0.0
    result['umapreadbp'] = sum([x['readbp'] for x in readstats if x['ismapped']==0])
    result['umapreadbppct'] = round(result['umapreadbp'] / float(result['runreadbp']) * 100.0, 3) if result['runreadbp'] else 0.0
    result['umapmeanreadbp'] = round(result['umapreadbp'] / float(result['umapreadcnt']), 3) if result['umapreadcnt'] else 0.0

    result['mapreadcnt'], result['mapreadcntpct'], result['mapreadbp'], result['mapreadbppct'], result['mapmeanreadbp'], \
    result['alncnt'], result['alnbp'], result['alnbppct'], result['alnmeanbp'], \
    result['Spct'], result['Cpct'], result['Mpct'], result['Ipct'], result['Dpct'], \
    result['Cmeandistbp'], result['Mmeandistbp'], result['Imeandistbp'], result['Dmeandistbp'], \
    result['Cmeanlen'], result['Mmeanlen'], result['Imeanlen'], result['Dmeanlen'], \
    result['targetaligncnt'], result['targetalignbp'], result['controlaligncnt'], result['controlalignbp'], \
    result['randalncnt'], result['randalncntpct'], result['randalnbp'], result['randalnbppct'], result['randalnmeanbp'] \
        = RunStats_Compute_FromSubset(readstats, result['runreadcnt'], result['runreadbp'], readfields)

    return result

def RunStats_Compute():
    '''
    Compute the read statistics from the readstats in _E.
    There should be 4 sets of statistics, for:
    allaln : all primary and secondary alignments
    prialn : all primary alignments
    oneref : all primary and secondary alignments for reads that match only one refcontig, and
    mulref : all primary and secondary alignments for reads that match more than one refcontig.
    '''

    allaln_readstats = [_E[k] for k in _E.keys()]
    oneref_readstats = [_E[k] for k in _E.keys() if _E[k]['numrefcontigs']<=1]      # aligned to 1 refcontig or unaligned
    mulref_readstats = [_E[k] for k in _E.keys() if _E[k]['numrefcontigs']!=1]      # aligned to >1 refcontig or unaligned
    global _R
    _R = {}
    _R['allaln'] = RunStats_Compute_One(allaln_readstats, 'all', 'allaln')
    _R['prialn'] = RunStats_Compute_One(allaln_readstats, 'pri', 'prialn')
    _R['oneref'] = RunStats_Compute_One(oneref_readstats, 'all', 'oneref')
    _R['mulref'] = RunStats_Compute_One(mulref_readstats, 'all', 'mulref')

def RunStats_Save():
    'Print the run statistics to a file.'
    keyL = ['allaln', 'prialn', 'oneref', 'mulref']
    with open(_outrunstatspath, 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join([var for var in _PT.runstatsheader])))
        for k in keyL:
            out_fp.write('{0}\n'.format('\t'.join([str(_R[k][var]) for var in _PT.runstatsheader])))

def RunStats_Print():
    RunStats_Compute()
    RunStats_Save()

# ============================================================================ #
# Main                                                                         #
# ============================================================================ #

if __name__ == '__main__':

    Initialise()
    Reference_Parse()
    AlignStats_Print(_PA.args.savealignments)
    ReadStats_Print()
    RunStats_Print()
    sys.stdout.write('Info: Finished poremapstats successfully\n')
    sys.exit(0)

# ============================================================================ #
