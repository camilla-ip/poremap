#!/usr/bin/env python

# ============================================================================ #
# poremap.py                                                                   #
# ============================================================================ #
'''
pipeline for nanopore read alignment and statistics

A mapping-based pipeline to report:
- genuine mappings of reads to a reference
- quality statistics for reads
- yield statistics for the run

Should work for data from:
- single-chromosome, pure bacterial or viral samples
- bacterial samples with multiple chromosome(s) and/or plasmid(s)
  where there may be regions of significant homology between parts
  of the chromosomes and plasmids
- eukaryotic chromosomes
- chimeric reads

Pipeline steps:
1.  Create a single FASTA with one contig per chromosome/plasmid/spike-in
2.  Align all reads to the set of target and control references using "bwa mem"
    and generate alignment statistics.
3.  Let the subset of reads that map to the control reference(s) be
    part of the training set for class 'real' and generate alignment statistics.
4.  Let the subset of reads that map to the control reference(s) be
    randomised in base order and quality to be the training set for class 'random'
    and generate alignment statistics.
5.  Let the subset of reads that map to the target reference(s)
    be the set of alignments to be classified.
6.  Classify the target alignments based on whether a linear regression
    of the principal components of some discriminating alignment statistics
    deem the alignment close to random alignments or non-random alignments.
7.  Generate the final set of alignment, read and run statistics based on
    non-random alignments.
8.  Save unique, chimeric and unmapped reads in BAM format (NOT IMPLEMENTED)
9.  Save statistics report in HTML format (NOT IMPLEMENTED)

Input:
- reads in FASTQ format
- reference(s) for target sample(s) in a single FASTA file
- reference(s) for control sample(s) in a single FASTA file

Output:
- outdir/outprefix_alignstats.txt
- outdir/outprefix_readstats.txt
- outdir/outprefix_runstats.txt

Pre-conditions:
- bindir must appear on the command-line to define where
  to find auxiliary programs and python packages.
- profilepath is useful for storing BASH commands for changing the environment
  before running it SGE scheduling.
'''
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# April 2015                                                                   #
# ============================================================================ #

# ============================================================================ #
# Import system-level Python modules                                           #
# ============================================================================ #

import subprocess as sp
import itertools, numpy as np, os, random, re, shlex, shutil, stat, sys, time
from Bio import SeqIO
from subprocess import PIPE
from ruffus import *

# ============================================================================ #
# Set _bindir and _profilepath, then import the classes in poremaplib/         #
#                                                                              #
# This complicated command-line parsing was introduced so the program could be #
# run using qsub using absolute pathnames (otherwise it doesn't work).         # 
# ============================================================================ #

_errorargsinvalid = 3

_bindir=None
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
            sys.stderr.write('Erro: Failed to import poremaplib classes\n')
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

_progdir = None
_progname = None
_intdir = None
_version = '0.1.1'
_progdesc = 'poremap mapping pipeline v{0}'.format(_version)

_mapper = {
    'bwa' : { 'params' : 'mem -x ont2d -M' },
    'graphmap' : { 'params' : '-x nanopore -C' }
}

_R = {}	# _R[refcontigid] = [reftype, refseq] where reftype is either 'target' or 'control'

_C = {} # Contents of alignclass path. _C[key] = isnonrandomaln, where key is a tuple of 11 fields

# ============================================================================ #
# Program usage                                                                #
# ============================================================================ #

def Initialise():
    'Parse and validate command-line arguments, set up global variables.'

    global _progname, _PA, _PT, _PL, _PE

    _progname = os.path.basename(os.path.realpath(sys.argv[0]))
    progexamples = []
    _PT = Poretypes()
    _PE = Poreerr()
    _PL = Porelog(progdoc=__doc__, bindir=_bindir, ini_version=_version)
    _PA = Poreargs(progdoc=__doc__, progexamples=progexamples, bindir=_bindir, version=_version)

  # General purpose arguments
    _PA.arg_addstd( '--bindir', required=True )
    _PA.arg_addstd( '--profilepath', required=False )
    _PA.arg_addstd( '--runid', required=True )
    _PA.arg_addstd( '--readtype', required=True )
    _PA.arg_addstd( '--readclass', required=True )
    _PA.arg_addstd( '--infastqpath', required=True )
    _PA.arg_addstd( '--targetrefpath', required=False )
    _PA.arg_addstd( '--controlrefpath', required=False )
    _PA.arg_addstd( '--outdir', required=True )
    _PA.arg_addstd( '--outprefix', required=True )
    _PA.arg_addstd( '--useintdir', required=False )
    _PA.arg_addstd( '--threads', required=False )
    _PA.arg_addstd( '--realmapprog', required=True, default='bwa' )
    _PA.arg_addstd( '--realmapparams', required=False, default=_mapper['bwa']['params'] )
    _PA.arg_addstd( '--randmapprog', required=True, default='graphmap' )
    _PA.arg_addstd( '--randmapparams', required=False, default=_mapper['graphmap']['params'] )
    _PA.arg_addstd( '--overwrite', required=False, default=False )
    _PA.arg_addstd( '--savealignments', required=False, default=False )
    _PA.arg_addstd( '--dryrun', required=False )
    _PA.arg_addstd( '--touch', required=False )
    _PA.arg_addstd( '--ruffusgraph', required=False )
    _PA.arg_addstd( '--verbose', required=False, default=False )

  # Formatting and misc
    _PA.arg_addstd( '--deleteintfiles', required=False )
    _PA.arg_addstd( '--logindentwidth', required=False )
    _PA.arg_addstd( '--covplotwinsz', required=False )
    _PA.arg_addstd( '--usecbcolours', required=False )

  # Third-party software
    _PA.arg_addstd( '--prog_bwa', required=False )
    _PA.arg_addstd( '--prog_graphmap', required=False )
    _PA.arg_addstd( '--prog_samtools', required=False )
    _PA.arg_addstd( '--prog_rscript', required=False )

  # This is the line that causes the program to exit if the arguments supplied are not sufficient or valid.
    _PA.arg_parseall()

  # Initialise the log output
    _PL.log_2out('Info: Started')
    _PL.log_2out('{0}'.format(_PL.log_runinfo()))
    _PL.log_2out('Info:\n\n{valuelines}'.format(
        valuelines=_PA.arg_valuelines(indentwidth=int(_PA.args.logindentwidth))))

def Exit(retkey):
    'Pring diagnostic messages and exit the program with a return value and error message from the poreerr class.'
    _PL.log_exit(_PE.err_retkey(retkey), _PE.err_text(retkey), print2out=True, print2dump=True, exitprogram=True)

def CheckInput():
    'If any errors found in input, print an error message and exit with the appropriate return value.'

    global _PA
    if _PA.args.targetrefpath is not None and _PA.args.targetrefpath.lower() == 'none':
        _PA.args.targetrefpath = None
    if _PA.args.controlrefpath is not None and _PA.args.controlrefpath.lower() == 'none':
        _PA.args.controlrefpath = None
    if _PA.args.profilepath is not None and _PA.args.profilepath.lower() == 'none':
        _PA.args.profilepath = None

    if not os.path.exists(os.path.expandvars(_PA.args.bindir)):
        _PL.log_2out('Erro: --bindir does not exist *{0}*'.format(_PA.args.bindir))
        Exit('ErrorDirMissing')
    if _PA.args.profilepath is not None and not os.path.exists(os.path.expandvars(_PA.args.profilepath)):
        _PL.log_2out('Erro: --profilepath does not exist *{0}*'.format(_PA.args.profilepath))
        Exit('ErrorFileMissing')
    if not os.path.exists(os.path.expandvars(_PA.args.infastqpath)):
        _PL.log_2out('Erro: --infastqpath does not exist *{0}*'.format(_PA.args.infastqpath))
        Exit('ErrorFileMissing')
    if _PA.args.targetrefpath is None and _PA.args.controlrefpath is None:
        _PL.log_2out('Erro: Both --targetrefpath and --controlrefpath are None - nothing to do')
        Exit('ErrorFileMissing')
    if _PA.args.targetrefpath is not None and not os.path.exists(os.path.expandvars(_PA.args.targetrefpath)):
        _PL.log_2out('Erro: --targetrefpath does not exist *{0}*'.format(_PA.args.targetrefpath))
        Exit('ErrorFileMissing')
    if _PA.args.controlrefpath is not None and not os.path.exists(os.path.expandvars(_PA.args.controlrefpath)):
        _PL.log_2out('Erro: --controlrefpath does not exist *{0}*'.format(_PA.args.controlrefpath))
        Exit('ErrorFileMissing')
    if not os.path.exists(os.path.expandvars(_PA.args.outdir)):
        try:
            os.makedirs(_PA.args.outdir)
        except:
            _PL.log_2out('Erro: --outdir does not exist *{0}*'.format(_PA.args.outdir))
            Exit('ErrorDirMissing')

    global prog_poremapstats
    prog_poremapstats = os.path.join(os.path.expandvars(_PA.args.bindir), "poremapstats.py")
    if not os.path.exists(prog_poremapstats):
        _PL.log_2out('Erro: poremapstats program does not exist *{0}*'.format(prog_poremapstats))
        Exit('ErrorFileMissing')

    global prog_poremapclassifier
    prog_poremapclassifier = os.path.join(os.path.expandvars(_PA.args.bindir), "poremapclassifier.R")
    if not os.path.exists(prog_poremapclassifier):
        _PL.log_2out('Erro: poremapclassifier program does not exist *{0}*'.format(prog_poremapclassifier))
        Exit('ErrorFileMissing')

    if not os.path.exists(os.path.expandvars(_PA.args.prog_bwa)):
        _PL.log_2out('Erro: bwa program does not exist ({0})'.format(_PA.args.prog_bwa))
        Exit('ErrorFileMissing')
    if not os.path.exists(os.path.expandvars(_PA.args.prog_graphmap)):
        _PL.log_2out('Erro: graphmap program does not exist ({0})'.format(_PA.args.prog_graphmap))
        Exit('ErrorFileMissing')
    if not os.path.exists(os.path.expandvars(_PA.args.prog_rscript)):
        _PL.log_2out('Erro: Rscript program does not exist ({0})'.format(_PA.args.prog_rscript))
        Exit('ErrorFileMissing')
    if not os.path.exists(os.path.expandvars(_PA.args.prog_samtools)):
        _PL.log_2out('Erro: samtools program does not exist ({0})'.format(_PA.args.prog_samtools))
        Exit('ErrorFileMissing')

def Step_Setup():
    'Do pre-workflow setup, like creating a sub-directory for intermediate files.'
    global _intdir, _logpath, _tmpdir

  # Create directory for intermediate files
    if _PA.args.useintdir:
        _intdir = os.path.join(_PA.args.outdir, "int")
        _tmpdir = os.path.join(_PA.args.outdir, "tmp")
    else:
        _intdir = _PA.args.outdir
        _tmpdir = os.path.join(_PA.args.outdir, "tmp")
    if not os.path.exists(_intdir):
        _PA.sys_makedir(_intdir)
    if not os.path.exists(_intdir):
        Exit('ErrorDirCreate')
    if not os.path.exists(_tmpdir):
        _PA.sys_makedir(_tmpdir)
    if not os.path.exists(_tmpdir):
        Exit('ErrorDirCreate')

  # Set the log file path to be outdir/outprefix_poremap_num.log where we try num from 1 to 1000
    for i in range(1, 1001):
        _logpath = os.path.join(_PA.args.outdir, "{outprefix}_poremap_{num}.log".format(outprefix=_PA.args.outprefix, num=i))
        logpathgz = "{0}.gz".format(_logpath)
        if not os.path.exists(_logpath) and not os.path.exists(logpathgz):
            break
    _PL.log_setoutpath(logoutpath=_logpath)
    _PL.log_dumpout('Info: Started')
    _PL.log_dumpout('{0}'.format(_PL.log_runinfo()))
    _PL.log_dump('Info:\n\n{valuelines}'.format(valuelines=_PA.arg_valuelines(indentwidth=int(_PA.args.logindentwidth))))

# ============================================================================ #
# Run the initialisation steps                                                 #
# ============================================================================ #

if __name__ == "__main__":
    # From here, all diagnostic messages written to stdout
    Initialise()
    CheckInput()
    Step_Setup()

# ============================================================================ #
# Generic utility functions                                                    #
# ============================================================================ #

#def sys_exec(cmd):
#    '''
#    Execute a command using the subprocess module to trap the
#    return value, the stdout string and stderr string.
#    '''
#    proc_handle = sp.Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
#    proc_stdout, proc_stderr = proc_handle.communicate()
#    proc_returncode = proc_handle.returncode
#    return [proc_returncode, proc_stdout, proc_stderr]

# ============================================================================ #
# Workflow                                                                     #
# ============================================================================ #

def Run_Shell_Script(stepname, cmd, dump_stdout=True):
    '''
    Run the command using the sub-process library, printing stepname in any
    diagnostic messages, print any stdout or stderr text to the log file and 
    return the retcode.
    '''
    _PL.log_dumpout('Info: Step {stepname} start: {cmd}'.format(stepname=stepname, cmd=cmd))
    retcode, outstr, errstr = _PA.sys_exec('{cmd}'.format(cmd=cmd))
    if dump_stdout and outstr and len(outstr.strip()):
        _PL.log_dumpout('Info: Step {stepname} after: {progname} script stdout:'.format(
            stepname=stepname, progname=_progname))
        _PL.log_dumpout(_PL.str_indented('\n'+outstr, indentwidth=int(_PA.args.logindentwidth)))
    if errstr and len(errstr):
        _PL.log_dumpout('Info: Step {stepname} after: {progname} script stderr:'.format(
            stepname=stepname, progname=_progname))
        _PL.log_dumpout(_PL.str_indented('\n'+errstr, indentwidth=int(_PA.args.logindentwidth)))
    _PL.log_dumpout('Info: Step {stepname} end: {progname} script retcode={retcode}'.format(
        stepname=stepname, progname=_progname, retcode=retcode))
    return retcode, outstr, errstr

def Step_Setup():
    'Do pre-workflow setup, like creating a sub-directory for intermediate files.'

# ============================================================================ #
# References                                                                   #
# ============================================================================ #

_pipeline_input = [ os.path.expandvars(_PA.args.infastqpath) ]

def Read_References():
    'Read both files into a single global variable _R.'
    global _R
    if not _R:
        if _PA.args.targetrefpath is not None:
            targetrefpath = os.path.expandvars(_PA.args.targetrefpath)
            for record in SeqIO.parse(targetrefpath, 'fasta'):
                _R[record.id] = ['target', str(record.seq)]
        if _PA.args.controlrefpath is not None:
            controlrefpath = os.path.expandvars(_PA.args.controlrefpath)
            for record in SeqIO.parse(controlrefpath, 'fasta'):
                _R[record.id] = ['control', str(record.seq)]

def Read_AlignClassFile():
    'Read alignclasspath into global dict _C.'
    global _C
    if len(_C):
        return
    alignclass_path = os.path.join(_intdir, _PA.args.outprefix+'_realmapped_alignclass.txt')
    data = np.loadtxt(alignclass_path, dtype=_PT.alignclassdtype, delimiter='\t', skiprows=1)
    for rowT in data:
        rowL = list(rowT)
        key = tuple(rowL[0:11])
        val = rowL[11]
        _C[key] = val

# ============================================================================ #
# Mapping                                                                      #
# ============================================================================ #

def Run_BwaMem(jobdesc, infastqpath, reffastapath, mapparams, outbampath):
    'Map the reads in the input BAM file with bwa mem, assuming they are ONT reads, saving secondary alignments.'
    # Set up file paths
    outsampath = os.path.join(_intdir, os.path.basename(outbampath).replace('.bam', '.sam'))
    bwalogpath = os.path.join(_intdir, os.path.basename(outbampath).replace('.bam', '.log'))
    outbamprefix = os.path.join(_intdir, os.path.basename(outbampath).replace('.bam', ''))
    tmpbampath = outbampath.replace('.bam', '.tmp.bam')
    # Run bwa mem
    if 1:
        cmd = '{prog} {params} {ref} {reads}'.format( \
            prog=_PA.args.prog_bwa, params=mapparams, ref=reffastapath, reads=infastqpath)
        rc, ro, re = Run_Shell_Script(jobdesc, cmd, dump_stdout=False)
        try:
            ro = ro.rstrip('\n')
            with open(outsampath, 'w') as sam_fp:
                sam_fp.write('{0}\n'.format(ro))
            with open(bwalogpath, 'w') as log_fp:
                log_fp.write('{0}\n'.format(re))
        except:
            pass
        if rc != 0:
             Exit('ErrorSysCall')
    else:
        cmd = '{prog} {params} {ref} {reads} > {sampath} 2> {logpath}'.format( \
            prog=_PA.args.prog_bwa, params=mapparams, ref=reffastapath, reads=infastqpath, sampath=outsampath, logpath=bwalogpath)
        os.system(cmd)
    # Convert SAM to BAM
    cmd = '{prog} view -Sb -o {bampath} {sampath}'.format(prog=_PA.args.prog_samtools, bampath=tmpbampath, sampath=outsampath)
    rc, ro, re = Run_Shell_Script(jobdesc, cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Sort BAM
    cmd = '{prog} sort {inbam} {outbamprefix}'.format(prog=_PA.args.prog_samtools, inbam=tmpbampath, outbamprefix=outbamprefix)
    rc, ro, re = Run_Shell_Script(jobdesc, cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Create BAM index
    cmd = '{prog} index {bampath}'.format(prog=_PA.args.prog_samtools, bampath=outbampath)
    rc, ro, re = Run_Shell_Script(jobdesc, cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Remove intermediate files
    if os.path.exists(outsampath):
        os.remove(outsampath)
    if os.path.exists(tmpbampath):
        os.remove(tmpbampath)
    # All ok
    return 0

def CleanSamIDs(tmpsam1, tmpsam2):
    '''
    Clean GraphMap output so that the readid is only the first space-separated word
    and the chromosome is only the first space-separated word.
    '''
    # Clean the readids and chromnames, save as SAM
    with open(tmpsam1, 'r') as in_fp:
        with open(tmpsam2, 'w') as out_fp:
            for line in in_fp:
                if line.startswith('@'):
                    if line.startswith('@SQ'):
                        L = line.strip().split('\t')
                        L[1] = L[1].split(' ')[0]
                        out_fp.write('{0}\n'.format('\t'.join(L)))
                    else:
                        out_fp.write(line)
                else:
                    L = line.strip().split('\t')
                    L[0] = L[0].split(' ')[0]
                    L[2] = L[2].split(' ')[0]
                    if L[2] == '*' and L[3] != '0':
                        L[3] = '0'
                    out_fp.write('{0}\n'.format('\t'.join(L)))
 
def Run_GraphMap(jobdesc, infastqpath, reffastapath, outbampath):
    'Map the fastq reads with GraphMap.'
    # Tmp files
    intsampath1 = outbampath.replace('.bam', '.withlongids.sam')
    intsampath2 = outbampath.replace('.bam', '.shortids.sam')
    outbamprefix = os.path.join(_intdir, os.path.basename(outbampath).replace('.bam', ''))
    tmpbampath = outbampath.replace('.bam', '.tmp1.bam')
    outpath = outbampath.replace('.bam', '.out')
    errpath = outbampath.replace('.bam', '.err')
    # Run GraphMap
    cmd = '{prog} {params} -r {reffasta} -d {infastq} -o {outsam} -t {threads}'.format(
        prog=_PA.args.prog_graphmap,
        params=_PA.args.randmapparams,
        reffasta=reffastapath,
        infastq=infastqpath,
        outsam=intsampath1,
        threads=_PA.args.threads
    )
    rc, ro, re = Run_Shell_Script(jobdesc, cmd, dump_stdout=False)
    try:
        if len(ro):
            with open(outpath, 'w') as out_fp:
                out_fp.write('{0}\n'.format(ro))
        if len(re):
            with open(errpath, 'w') as err_fp:
                err_fp.write('{0}\n'.format(re))
    except:
        pass
    if rc != 0:
         Exit('ErrorSysCall')
    # Clean the readids and chromname in the GraphMap output
    CleanSamIDs(intsampath1, intsampath2)
    # Convert SAM to BAM
    cmd = '{prog} view -bS -o {outbam} {outsam}'.format(prog=_PA.args.prog_samtools, outbam=tmpbampath, outsam=intsampath2)
    print cmd
    rc, ro, re = Run_Shell_Script(jobdesc, cmd)
    if rc != 0:
         Exit('ErrorSysCall')
    # Sort BAM
    cmd = '{prog} sort {inbam} {outbamprefix}'.format(prog=_PA.args.prog_samtools, inbam=tmpbampath, outbamprefix=outbamprefix)
    rc, ro, re = Run_Shell_Script(jobdesc, cmd)
    if rc != 0:
         Exit('ErrorSysCall')
    # Create BAM index
    cmd = '{prog} index {bampath}'.format(prog=_PA.args.prog_samtools, bampath=outbampath)
    print cmd
    rc, ro, re = Run_Shell_Script(jobdesc, cmd)
    if rc != 0:
         Exit('ErrorSysCall')
    # Remove intermediate files
    if os.path.exists(intsampath1):
        os.remove(intsampath1)
    if os.path.exists(intsampath2):
        os.remove(intsampath2)
    if os.path.exists(tmpbampath):
        os.remove(tmpbampath)
    # All ok
    return 0

def Count_BamReads(inreadsbam):
    'Return number of reads in the bam file.'
    cmd  = '{samtools} view -c {bampath}'.format(samtools=_PA.args.prog_samtools, bampath=inreadsbam)
    print cmd
    rc, ro, re = Run_Shell_Script('Count_BamReads', cmd, dump_stdout=False)
    if rc != 0:
        Exit('ErrorSysCall')
    ro = ro.strip()
    readcnt = int(ro) if len(ro) else 0
    return readcnt

def Run_Poremapstats(jobdesc, runid, readtype, readclass, datatype, mapprog, mapparams, alignclasspath, inreadsbam, inreffasta, outdir, outprefix):
    'Run the poremapstats.py program to generate statistics on each read, or just return if there are no reads in the BAM file.'
    # Run poremapstats on non-empty BAM file
    cmd ='{0}'.format(prog_poremapstats)
    cmd += ' --bindir {0}'.format(_PA.args.bindir)
    cmd += ' --runid {0}'.format(runid)
    cmd += ' --readtype {0}'.format(readtype)
    cmd += ' --readclass {0}'.format(readclass)
    cmd += ' --datatype {0}'.format(datatype)
    cmd += ' --mapprog {0}'.format(mapprog)
    cmd += ' --mapparams "{0}"'.format(mapparams)
    cmd += ' --alignclasspath {0}'.format(alignclasspath)
    cmd += ' --readsbam {0}'.format(inreadsbam)
    cmd += ' --targetrefpath {0}'.format(_PA.args.targetrefpath)
    cmd += ' --controlrefpath {0}'.format(_PA.args.controlrefpath)
    cmd += ' --outdir {0}'.format(outdir)
    cmd += ' --outprefix {0}'.format(outprefix)
    cmd += ' --overwrite {0}'.format(_PA.args.overwrite)
    if _PA.args.savealignments:
        cmd += ' --savealignments True'
    print cmd
    rc, ro, re = Run_Shell_Script(jobdesc, cmd)
    if rc != 0:
         Exit('ErrorSysCall')
    return 0

# ============================================================================ #
# New pipeline : alignment
# ============================================================================ #

def Align_Reads(jobdesc, readspath, refpath, mapprog, mapparams, bampath):
    'Map the reads to the reference contigs using the program specified, saving the output in BAM format.'

    if mapprog == 'bwa':
        Run_BwaMem(jobdesc, readspath, refpath, mapparams, bampath)
    elif mapprog == 'graphmap':
        Run_GraphMap()
    else:
        _PL.log_dumpout('Erro: Invalid alignment program specified ({0})'.format(mapprog))
        Exit('ErrorArgsInvalid')

def Convert_InitstatsToAlignclass(alignstats, alignclass):
    'Convert the initstats files into an alignclass file.'
    with open(alignclass, 'w') as out_fp:
        with open(alignstats, 'r') as in_fp:
            for line in in_fp:
                L = line[:-1].split('\t')
                N = L[0:11] + [L[13]]
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in N])))

def Run_AlignClassifier(runid, trainingalignstats, testingalignstats, outdir, outprefix):
    'Run the poremapclassifier program to get the isnonrandomaln field for each target alignment.'
    _PL.log_dumpout('Info: Run_AlignClassifier: runid={0}, trainingalignstats={1}, testingalignstats={2}, outdir={3}, outprefix={4}'.format(runid, trainingalignstats, testingalignstats, outdir, outprefix))

    if _PA.args.targetrefpath is None:
        try:
            testingalignclass = os.path.join(outdir, outprefix+'_alignclass.txt')
            Convert_InitstatsToAlignclass(testingalignstats, testingalignclass)
            return
        except:
            _PL.log_dumpout('Erro: Failed to convert target_initstats to target_alignclass file when target is None ({0}, {1})'.format(
                testingalignstats, testingalignclass))
            Exit('ErrorSysCall')

    data = open(trainingalignstats, 'r').read().strip().split('\n')
    if len(data) == 1 and data[0].startswith('runid'):
        try:
            testingalignclass = os.path.join(outdir, outprefix+'_alignclass.txt')
            Convert_InitstatsToAlignclass(testingalignstats, testingalignclass)
            return
        except:
            _PL.log_dumpout('Erro: Failed to convert target_initstats to target_alignclass file when no control reads ({0}, {1})'.format(
                testingalignstats, testingalignclass))
            Exit('ErrorSysCall')

    data = open(testingalignstats, 'r').read().strip().split('\n')
    if len(data) == 1 and data[0].startswith('runid'):
        try:
            testingalignclass = os.path.join(outdir, outprefix+'_alignclass.txt')
            Convert_InitstatsToAlignclass(testingalignstats, testingalignclass)
            return
        except:
            _PL.log_dumpout('Erro: Failed to convert target_initstats to target_alignclass file when no target reads ({0}, {1})'.format(
                testingalignstats, testingalignclass))
            Exit('ErrorSysCall')

    cmd = '{rscript} {prog} {runid} {trainingdata} {testingdata} {outdir} {outprefix}'.format(
        runid=runid,
        rscript=_PA.args.prog_rscript,
        prog=prog_poremapclassifier,
        trainingdata=trainingalignstats,
        testingdata=testingalignstats,
        outdir=outdir,
        outprefix=outprefix)
    rc, ro, re = Run_Shell_Script('Run_AlignClassifier', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.strip()
    if len(ro):
        log_path = os.path.join(outdir, outprefix+'_alignclass.log')
        with open(log_path, 'w') as out_fp:
            out_fp.write('{0}\n'.format(ro))

# ============================================================================ #
# New pipeline :steps
# ============================================================================ #

@transform( \
    _pipeline_input, \
    regex(r"(\S+)/(\S+).fastq"), \
    [r'{intdir}/{outprefix}_targetandcontrol.fasta'.format(intdir=_intdir, outprefix=_PA.args.outprefix), \
     r'{intdir}/{outprefix}_targetandcontrol.fasta.bwt'.format(intdir=_intdir, outprefix=_PA.args.outprefix)])
def References_Setup(inpaths, outpaths):
    '''
    Concatenate the reference files into a single .fasta file.
    Currently can't get the suffix specifier to match all the the .fasta files properly
    so have done it the kludgy way for now.
    '''
    _PL.log_dumpout('Info: Step References_Setup: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))
    # Merge both references files into a single FASTA file
    mergedrefpath = os.path.join(_intdir, _PA.args.outprefix+'_targetandcontrol.fasta')
    out_fp = open(mergedrefpath,'wb')
    if _PA.args.targetrefpath is not None:
        targetrefpath = os.path.expandvars(_PA.args.targetrefpath)
        shutil.copyfileobj(open(targetrefpath,'rb'), out_fp)
    if _PA.args.controlrefpath is not None:
        controlrefpath = os.path.expandvars(_PA.args.controlrefpath)
        shutil.copyfileobj(open(controlrefpath,'rb'), out_fp)
    out_fp.close()
    cmd = '{prog} index {fasta}'.format(prog=_PA.args.prog_bwa, fasta=mergedrefpath)
    rc, ro, re = Run_Shell_Script('References_Setup', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    #time.sleep(2)
    return 0

@transform( \
    References_Setup, \
    regex(r"(\S+)/(\S+)_targetandcontrol.fasta"), \
    [r'{intdir}/{outprefix}_realmapped.bam'.format(intdir=_intdir, outprefix=_PA.args.outprefix), \
     r'{intdir}/{outprefix}_realmapped.bam.bai'.format(intdir=_intdir, outprefix=_PA.args.outprefix)])
def Get_Alignments_Initial(inpaths, outpaths):
    'Map all the input fastq reads with our trusted read mapper.'
    _PL.log_dumpout('Info: Get_Alignments_Initial: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))
    jobdesc = 'Mapping all real reads with {0}'.format(_PA.args.realmapprog)
    readspath = _pipeline_input[0]
    mapprog = _PA.args.realmapprog
    mapparams = _PA.args.realmapparams
    refpath = os.path.join(_intdir, _PA.args.outprefix+'_targetandcontrol.fasta')
    bampath = outpaths[0]
    Align_Reads(jobdesc, readspath, refpath, mapprog, mapparams, bampath)
    return 0

@transform( \
    Get_Alignments_Initial, \
    regex(r"(\S+)/(\S+)_realmapped.bam"), \
    [r'{intdir}/\2_realmapped_control.bam'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_control.bam.bai'.format(intdir=_intdir)])
def Get_Alignments_Training_ControlReal(inpaths, outpaths):
    'Extract the primary alignments to the control references from _realmapped.bam into _realmappedcontrol.bam.'
    _PL.log_dumpout('Info: Get_Alignments_Training_ControlReal: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))
    # Set up file paths
    inbam = [x for x in inpaths if x.endswith('_realmapped.bam')][0]
    tmpsam = outpaths[0].replace('.bam', '.sam')
    tmpbam = outpaths[0].replace('.bam', '.tmp.bam')
    outbam = outpaths[0]
    outbamprefix = outpaths[0].replace('.bam', '')
    # Get SAM records : inbam -> tmpsam
    cmd = '{samtools} view -H {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Training_ControlReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.rstrip('\n')
    if len(ro):
        with open(tmpsam, 'w') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    Read_References()
    controlrefcontigidL = [refcontigid for refcontigid in _R.keys() if _R[refcontigid][0]=='control']
    if len(controlrefcontigidL):
        cmd1 = '{samtools} view -F2304 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        cmd2 = 'awk \'{0}\''.format(' || '.join(["$3 == \"{0}\"".format(x) for x in controlrefcontigidL]))
        _PL.log_dump('Info: Get_Alignments_Training_ControlReal: {0} | {1}'.format(cmd1, cmd2))
        cmd1_h = sp.Popen(shlex.split(cmd1), stdout=PIPE)
        cmd2_h = sp.Popen(shlex.split(cmd2), stdin=cmd1_h.stdout, stdout=PIPE)
        ro, re = cmd2_h.communicate()
        rc = cmd2_h.returncode
    else:
        cmd = '{samtools} view -F2304 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        _PL.log_dump('Info: Get_Alignments_Training_ControlReal: {0}'.format(cmd))
        rc, ro, re = Run_Shell_Script('Get_Alignments_Training_ControlReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.rstrip('\n')
    if len(ro):
        with open(tmpsam, 'a') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    # Convert SAM to BAM : tmpsam -> tmpbam
    cmd = '{prog} view -Sb -o {bampath} {sampath}'.format(prog=_PA.args.prog_samtools, bampath=tmpbam, sampath=tmpsam)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Training_ControlReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Sort BAM : tmpbam -> outbam
    cmd = '{prog} sort {inbam} {outbamprefix}'.format(prog=_PA.args.prog_samtools, inbam=tmpbam, outbamprefix=outbamprefix)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Training_ControlReal', cmd)
    if rc != 0:
         Exit('ErrorSysCall')
    # Create BAM index : outbam -> outbam.bai
    cmd = '{prog} index {bampath}'.format(prog=_PA.args.prog_samtools, bampath=outbam)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Training_ControlReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Remove intermediate files
    if os.path.exists(tmpsam):
        os.remove(tmpsam)
    if os.path.exists(tmpbam):
        os.remove(tmpbam)
    # All ok
    return 0

def ReadidsFromBam(bampath):
    'Return the list of unique readids from a BAM file.'
    cmd1 = '{samtools} view {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=bampath)
    cmd2 = 'cut -f1'
    cmd3 = 'sort -u'
    cmd1_h = sp.Popen(shlex.split(cmd1), stdout=PIPE)
    cmd2_h = sp.Popen(shlex.split(cmd2), stdin=cmd1_h.stdout, stdout=PIPE)
    cmd3_h = sp.Popen(shlex.split(cmd3), stdin=cmd2_h.stdout, stdout=PIPE)
    ro, re = cmd3_h.communicate()
    rc = cmd3_h.returncode
    if rc != 0:
         Exit('ErrorSysCall')
    readidL = ro.strip().split('\n')
    return readidL

@transform(
    Get_Alignments_Training_ControlReal, \
    regex(r"(\S+)/(\S+)_realmapped_control.bam"), \
    [r'{intdir}/\2_randmapped_control.bam'.format(intdir=_intdir), \
     r'{intdir}/\2_randmapped_control.bam.bai'.format(intdir=_intdir)])
def Get_Alignments_Training_ControlRand(inpaths, outpaths):
    '''
    Create pseudo-random control reads by shuffling the order of the bases and base
    qualities of the reads mentioned in _realmapped_control.bam. Then map these reads
    to all the target and control references, saving the output alignments as a
    sorted bam file with index.
    '''
    _PL.log_dumpout('Info: Get_Alignments_Training_ControlRand: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))
    # File paths
    infastq = _pipeline_input[0]
    randcontrolfastq = outpaths[0].replace('_randmapped_control.bam', '_rand_control.fastq')
    inbam = [x for x in inpaths if x.endswith('.bam')][0]
    outbam = outpaths[0]
    reffastapath = os.path.join(_intdir, _PA.args.outprefix+'_targetandcontrol.fasta')
    # Set up random control reads in FASTQ format
    if Count_BamReads(inbam):
        readidL = ReadidsFromBam(inbam)
    else:
        readidL = []
    with open(randcontrolfastq, 'w') as out_fp:
        for record in SeqIO.parse(infastq, 'fastq'):
            id = record.id
            if id in readidL:
                tmpB = list(str(record.seq))
                random.shuffle(tmpB)
                seq = ''.join(tmpB)
                tmpQ = record.letter_annotations['phred_quality']
                random.shuffle(tmpQ)
                bq = ''.join([chr(n+33) for n in tmpQ])
                out_fp.write('@{id}\n{seq}\n+\n{bq}\n'.format(id=id, seq=seq, bq=bq))
    # Map the random control reads to all the target and control reference contigs.
    Run_GraphMap('Get_Alignments_Training_ControlRand', randcontrolfastq, reffastapath, outbam)
    # Remove intermediate files.
    if os.path.exists(randcontrolfastq):
        os.remove(randcontrolfastq)
    return 0

@transform( \
    Get_Alignments_Training_ControlRand, \
    regex(r"(\S+)/(\S+)_randmapped_control.bam"), \
    [r'{intdir}/\2_realmapped_target.bam'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_target.bam.bai'.format(intdir=_intdir)])
def Get_Alignments_Testing_TargetReal(inpaths, outpaths):
    '''
    From a statistical point of view, it is not valid to have any overlap between
    the training and test data. Hence, the _realmapped_target.bam file of alignments
    to be classified consists of all alignments that are not in the _realmapped_control.bam file.
    '''
    _PL.log_dumpout('Info: Get_Alignments_Testing_TargetReal: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))
    # Set up file paths
    inbam = [x for x in inpaths if x.endswith('_randmapped_control.bam')][0].replace('_randmapped_control.bam', '_realmapped.bam')
    tmpsam = outpaths[0].replace('.bam', '.sam')
    tmpbam = outpaths[0].replace('.bam', '.tmp.bam')
    outbam = outpaths[0]
    outbamprefix = outpaths[0].replace('.bam', '')
    # Get SAM records : inbam -> tmpsam
    # - Get header
    cmd = '{samtools} view -H {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Testing_TargetReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.rstrip('\n')
    if len(ro):
        with open(tmpsam, 'w') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    # - Get non-primary alignments to control refcontig(s)
    Read_References()
    controlrefcontigidL = [refcontigid for refcontigid in _R.keys() if _R[refcontigid][0]=='control']
    if len(controlrefcontigidL):
        cmd1 = '{samtools} view -f256 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        cmd2 = 'awk \'{0}\''.format(' || '.join(["$3 == \"{0}\"".format(x) for x in controlrefcontigidL]))
        _PL.log_dump('Info: Get_Alignments_Testing_TargetReal: {0} | {1}'.format(cmd1, cmd2))
        cmd1_h = sp.Popen(shlex.split(cmd1), stdout=PIPE)
        cmd2_h = sp.Popen(shlex.split(cmd2), stdin=cmd1_h.stdout, stdout=PIPE)
        ro, re = cmd2_h.communicate()
        rc = cmd2_h.returncode
    else:
        cmd = '{samtools} view -f256 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        rc, ro, re = Run_Shell_Script('Get_Alignments_Testing_TargetReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.rstrip('\n')
    if len(ro):
        with open(tmpsam, 'a') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    # - Get supplementary alignments to control refcontig(s)
    if len(controlrefcontigidL):
        cmd1 = '{samtools} view -f2048 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        cmd2 = 'awk \'{0}\''.format(' || '.join(["$3 == \"{0}\"".format(x) for x in controlrefcontigidL]))
        _PL.log_dump('Info: Get_Alignments_Testing_TargetReal: {0} | {1}'.format(cmd1, cmd2))
        cmd1_h = sp.Popen(shlex.split(cmd1), stdout=PIPE)
        cmd2_h = sp.Popen(shlex.split(cmd2), stdin=cmd1_h.stdout, stdout=PIPE)
        ro, re = cmd2_h.communicate()
        rc = cmd2_h.returncode
    else:
        cmd = '{samtools} view -f2048 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        rc, ro, re = Run_Shell_Script('Get_Alignments_Testing_TargetReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.rstrip('\n')
    if len(ro):
        with open(tmpsam, 'a') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    # - Get all alignments to target (i.e., non-control) refcontig(s)
    if len(controlrefcontigidL):
        cmd1 = '{samtools} view -F4 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        cmd2 = 'awk \'{0}\''.format(' && '.join(["$3 != \"{0}\"".format(x) for x in controlrefcontigidL]))
        _PL.log_dump('Info: Get_Alignments_Testing_TargetReal: {0} | {1}'.format(cmd1, cmd2))
        cmd1_h = sp.Popen(shlex.split(cmd1), stdout=PIPE)
        cmd2_h = sp.Popen(shlex.split(cmd2), stdin=cmd1_h.stdout, stdout=PIPE)
        ro, re = cmd2_h.communicate()
        rc = cmd2_h.returncode
    else:
        cmd = '{samtools} view -F4 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
        rc, ro, re = Run_Shell_Script('Get_Alignments_Testing_TargetReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.rstrip('\n')
    if len(ro):
        with open(tmpsam, 'a') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    # Convert SAM to BAM : tmpsam -> tmpbam
    cmd = '{prog} view -Sb -o {bampath} {sampath}'.format(prog=_PA.args.prog_samtools, bampath=tmpbam, sampath=tmpsam)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Testing_TargetReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Sort BAM : tmpbam -> outbam
    cmd = '{prog} sort {inbam} {outbamprefix}'.format(prog=_PA.args.prog_samtools, inbam=tmpbam, outbamprefix=outbamprefix)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Testing_TargetReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Create BAM index : outbam -> outbam.bai
    cmd = '{prog} index {bampath}'.format(prog=_PA.args.prog_samtools, bampath=outbam)
    rc, ro, re = Run_Shell_Script('Get_Alignments_Testing_TargetReal', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    # Remove intermediate files
    if os.path.exists(tmpsam):
        os.remove(tmpsam)
    if os.path.exists(tmpbam):
        os.remove(tmpbam)
    # All ok
    return 0

@transform( \
    Get_Alignments_Testing_TargetReal, \
    regex(r"(\S+)/(\S+)_realmapped_target.bam"), \
    [r'{intdir}/\2_target_initstats.txt'.format(intdir=_intdir), \
     r'{intdir}/\2_control_initstats.txt'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_target_alignclass.txt'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_alignclass.txt'.format(intdir=_intdir)])
def Get_Alignment_TargetReal_Classifications(inpaths, outpaths):
    '''
    Let the training data be the primary alignments to the control refcontig(s)
    plus the randomised versions of any readids in the primary control ref alignments.
    Let the testing data be all alignments of "real" reads not in the training set.
    Generate unclassified alignstats.txt files for the training and testing sets.
    Run the classifier to output a classified alignstats.txt file for the testing set. 
    '''
    _PL.log_dumpout('Info: Get_Alignment_TargetReal_Classifications: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))
    # D <- data[c("readtype", "Iperalnbase", "meanrunlenC", "Mperalnbase", "proplenC", "proplenI")]
    # File paths
    outdir = os.path.dirname(inpaths[0])
    outprefix = [os.path.basename(x) for x in inpaths if x.endswith('_realmapped_target.bam')][0].replace('_realmapped_target.bam', '')
    reffasta = os.path.join(_intdir, _PA.args.outprefix+'_targetandcontrol.fasta')
    target_real_bam = os.path.join(outdir, outprefix+'_realmapped_target.bam')
    control_real_bam = os.path.join(outdir, outprefix + '_realmapped_control.bam')
    control_rand_bam = os.path.join(outdir,outprefix + '_randmapped_control.bam')
    control_real_initstatspath = os.path.join(outdir, outprefix + '_realmapped_control_initstats.txt')
    control_rand_initstatspath = os.path.join(outdir, outprefix + '_randmapped_control_initstats.txt')
    control_initstatspath = os.path.join(outdir, outprefix + '_control_initstats.txt')
    target_initstatspath = os.path.join(outdir, outprefix + '_target_initstats.txt')
    target_alignclasspath = os.path.join(outdir, outprefix + '_realmapped_target_alignclass.txt')
    final_alignclasspath = os.path.join(outdir, outprefix + '_realmapped_alignclass.txt')

    #target_alignstatspath = outprefix + '_target_alignstats.txt'
    #final_alignstatspath = outprefix + '_alignstats.txt'

    # Count the number of reads in each of the BAM files
    target_readcnt = Count_BamReads(target_real_bam)
    control_readcnt = Count_BamReads(control_real_bam)

    if target_readcnt == 0:
        #open(target_initstatspath, 'a').close()
        with open(target_initstatspath, 'w') as out_fp:
            out_fp.write('{0}\n'.format('\t'.join(_PT.alignstatsheader)))
    else:
        Run_Poremapstats('Run_Poremapstats real target',
            _PA.args.runid, _PA.args.readtype, _PA.args.readclass, 'minion',
            _PA.args.realmapprog, _PA.args.realmapparams, 'None',
            target_real_bam, reffasta,
            outdir, outprefix+'_realmapped_target')
        os.rename(os.path.join(outdir, outprefix+'_realmapped_target_initstats.txt'), target_initstatspath)

    # Collate the initstats.txt files for the training and testing sets of alignments
    if control_readcnt == 0:
        #open(control_initstatspath, 'a').close()
        #open(control_rand_initstatspath, 'a').close()
        with open(control_initstatspath, 'w') as out_fp:
            out_fp.write('{0}\n'.format('\t'.join(_PT.alignstatsheader)))
        with open(control_real_initstatspath, 'w') as out_fp:
            out_fp.write('{0}\n'.format('\t'.join(_PT.alignclassheader)))
        with open(control_rand_initstatspath, 'w') as out_fp:
            out_fp.write('{0}\n'.format('\t'.join(_PT.alignclassheader)))
    else:
        Run_Poremapstats('Run_Poremapstats real control',
            _PA.args.runid, _PA.args.readtype, _PA.args.readclass, 'minion',
            _PA.args.realmapprog, _PA.args.realmapparams, 'None',
            control_real_bam, reffasta,
            outdir, outprefix+'_realmapped_control')
        Run_Poremapstats('Run_Poremapstats rand control',
            _PA.args.runid, _PA.args.readtype, _PA.args.readclass, 'rand',
            _PA.args.randmapprog, _PA.args.randmapparams, 'None',
            control_rand_bam, reffasta,
            outdir, outprefix+'_randmapped_control')

    if os.path.exists(control_real_initstatspath) and os.stat(control_real_initstatspath).st_size > 0:
        cmd = 'head -n 1 {0}'.format(control_real_initstatspath)
    elif os.path.exists(target_initstatspath) and os.stat(target_initstatspath).st_size > 0:
        cmd = 'head -n 1 {0}'.format(target_initstatspath)
    rc, ro, re = Run_Shell_Script('Get_Alignment_TargetReal_Classifications', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.strip()
    if len(ro):
        with open(control_initstatspath, 'w') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    cmd = 'tail -n +2 {0}'.format(control_real_initstatspath)
    rc, ro, re = Run_Shell_Script('Get_Alignment_TargetReal_Classifications', cmd, dump_stdout=False)
    #if rc != 0:
    #     Exit('ErrorSysCall')
    ro = ro.strip()
    if len(ro):
        with open(control_initstatspath, 'a') as out_fp:
            out_fp.write('{0}\n'.format(ro))
    cmd = 'tail -n +2 {0}'.format(control_rand_initstatspath)
    rc, ro, re = Run_Shell_Script('Get_Alignment_TargetReal_Classifications', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')
    ro = ro.strip()
    if len(ro):
        with open(control_initstatspath, 'a') as out_fp:
            out_fp.write('{0}\n'.format(ro))

    # Run the classifier, which outputs the alignstats.txt file for the testing set of alignments.
    outclassifierprefix = os.path.basename(outprefix+'_realmapped_target')
    Run_AlignClassifier(_PA.args.runid, control_initstatspath, target_initstatspath, outdir, outclassifierprefix)

    # Create the final initstats.txt file with:
    # - all the records from _realmapped_target_initstats.txt,
    # - the realmapped_control_initstats.txt records with isnonrandomaln set to 1
    # - the unmapped reads with isnonrandomaln set to -1
    shutil.copyfile(target_alignclasspath, final_alignclasspath)
    with open(final_alignclasspath, 'a') as out_fp:
        with open(control_real_initstatspath, 'r') as in_fp:
            # real mapped control alignments
            linecnt = 0
            for line in in_fp:
                linecnt += 1
                if linecnt == 1:
                    continue
                L = line.strip().split('\t')
                newL = L[0:11] + ['1']
                out_fp.write('{0}\n'.format('\t'.join(newL)))
            # unmapped reads
            inbam=inpaths[0].replace('_realmapped_target.bam', '_realmapped.bam')
            cmd1 = '{samtools} view -f4 {bam}'.format(samtools=_PA.args.prog_samtools, bam=inbam)
            cmd2 = 'cut -f1-5'
            _PL.log_dump('Info: Get_Alignment_TargetReal_Classifications: {0} | {1}'.format(cmd1, cmd2))
            cmd1_h = sp.Popen(shlex.split(cmd1), stdout=PIPE)
            cmd2_h = sp.Popen(shlex.split(cmd2), stdin=cmd1_h.stdout, stdout=PIPE)
            ro, re = cmd2_h.communicate()
            rc = cmd2_h.returncode
            if rc != 0:
                 Exit('ErrorSysCall')
            if len(ro) and ro[-2] == '\n':
                ro = ro[:-1]
            recordL = ro.split('\n')
            for line in recordL:
                line = line.strip()
                if not len(line):
                    continue
                readid, samflag, refcontigid, refcontigpos1, mapq = line.split('\t')
                datatype = 'minion'
                isnonrandomaln = '-1'
                newL = [_PA.args.runid, readid, _PA.args.readtype, _PA.args.readclass, datatype,
                    _PA.args.realmapprog, _PA.args.realmapparams, samflag, refcontigid, refcontigpos1,
                    mapq, isnonrandomaln]
                out_fp.write('{0}\n'.format('\t'.join(newL)))

    # Remove intermediate files
    return 0

@transform( \
    Get_Alignment_TargetReal_Classifications, \
    regex(r"(\S+)/(\S+)_target_initstats.txt"), \
    [r'{intdir}/\2_realmapped_alignclass.txt'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_alignstats.txt'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_readstats.txt'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_runstats.txt'.format(intdir=_intdir)])
def Get_Alignment_Statistics(inpaths, outpaths):
    ''
    _PL.log_dumpout('Info: Get_Alignment_Statistics: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))
    # File paths
    alignclasspath = [x for x in inpaths if x.endswith('_realmapped_alignclass.txt')][0]
    realmapped_bam = alignclasspath.replace('_realmapped_alignclass.txt', '_realmapped.bam')
    reffasta = os.path.join(_intdir, _PA.args.outprefix+'_targetandcontrol.fasta')
    outdir = os.path.dirname(inpaths[0])
    outprefix = os.path.basename(alignclasspath).replace('_realmapped_alignclass.txt', '')

    # Generate the final statistics files
    Run_Poremapstats('Run_Poremapstats all alignments',
        _PA.args.runid, _PA.args.readtype, _PA.args.readclass, 'minion',
        _PA.args.realmapprog, _PA.args.realmapparams, alignclasspath,
        realmapped_bam, reffasta,
        outdir, outprefix+'_realmapped')

    return 0

def SAMtoBAM(jobdesc, insam, outbam):
    cmd = '{samtools} view -Sb -o {outbam} {insam}'.format(
        samtools=_PA.args.prog_samtools,
        outbam=outbam,
        insam=insam)
    rc, ro, re = Run_Shell_Script(jobdesc, cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')

def BAMindex(jobdesc, inbam):
    cmd = 'samtools index {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=inbam)
    rc, ro, re = Run_Shell_Script(jobdesc, cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')

@transform( \
    Get_Alignment_Statistics, \
    regex(r"(\S+)/(\S+)_realmapped_alignclass.txt"), \
    [r'{intdir}/\2_realmapped_targetrand.bam'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_targetnonrand.bam'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_controlrand.bam'.format(intdir=_intdir), \
     r'{intdir}/\2_realmapped_controlnonrand.bam'.format(intdir=_intdir)])
def Get_SplitBams(inpaths, outpaths):
    '''
    Split the 2 prefix_realmaped_[control|target].bam files into
    4 prefix_realmaped_[control|target]_[nonrand|rand].bam files
    based on the classification in the alignclass.txt file.
    '''
    _PL.log_dumpout('Info: Get_SplitBams: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))

    # Pre-requisites
    Read_References()

    # IO paths
    alignclass_path = [x for x in inpaths if x.endswith('_alignclass.txt')][0]
    realmapped_path = alignclass_path.replace('_alignclass.txt', '.bam')
    realmappedtargetrand_path = realmapped_path.replace('.bam', '_targetrand.bam')
    realmappedtargetnonrand_path = realmapped_path.replace('.bam', '_targetnonrand.bam')
    realmappedcontrolrand_path = realmapped_path.replace('.bam', '_controlrand.bam')
    realmappedcontrolnonrand_path = realmapped_path.replace('.bam', '_controlnonrand.bam')

    # Intermediate paths
    realmappedsam_path = realmapped_path.replace('.bam', '.tmp.sam')
    realmappedtargetrandsam_path = realmapped_path.replace('.bam', '_targetrand.tmp.sam')
    realmappedtargetnonrandsam_path = realmapped_path.replace('.bam', '_targetnonrand.tmp.sam')
    realmappedcontrolrandsam_path = realmapped_path.replace('.bam', '_controlrand.tmp.sam')
    realmappedcontrolnonrandsam_path = realmapped_path.replace('.bam', '_controlnonrand.tmp.sam')

    # Read in alignclass file
    Read_AlignClassFile()
    #global _C
    #data = np.loadtxt(alignclass_path, dtype=_PT.alignclassdtype, delimiter='\t', skiprows=1)
    #_C = {}
    #for rowT in data:
    #    rowL = list(rowT)
    #    key = tuple(rowL[0:11])
    #    val = rowL[11]
    #    _C[key] = val

    # Convert input BAM to SAM
    cmd = '{samtools} view -h -o {outsam} {inbam}'.format(
        samtools=_PA.args.prog_samtools,
        outsam=realmappedsam_path,
        inbam=realmapped_path)
    rc, ro, re = Run_Shell_Script('realmapped BAM to SAM', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')

    # Iterate through SAM file, outputting records to appropriate output SAM file
    out_fp = {}
    out_fp['target'] = {}
    out_fp['control'] = {}
    out_fp['target']['rand'] = open(realmappedtargetrandsam_path, 'w')
    out_fp['target']['nonrand'] = open(realmappedtargetnonrandsam_path, 'w')
    out_fp['control']['rand'] = open(realmappedcontrolrandsam_path, 'w')
    out_fp['control']['nonrand'] = open(realmappedcontrolnonrandsam_path, 'w')
    with open(realmappedsam_path, 'r') as in_fp:
        for line in in_fp:
            if line.startswith('@'):
                for readtype in ['target', 'control']:
                    for aligntype in ['rand', 'nonrand']:
                        out_fp[readtype][aligntype].write(line)
            else:
                L = line.split('\t')
                readid, samflag, refcontigid, refcontigpos1, mapq = L[0:5]
                if refcontigid not in _R.keys():
                    continue
                datatype = 'minion'
                key = (_PA.args.runid, readid, _PA.args.readtype, _PA.args.readclass, datatype,
                    _PA.args.realmapprog, _PA.args.realmapparams, int(samflag), refcontigid, int(refcontigpos1), int(mapq))
                readtype = _R[refcontigid][0]
                aligntype = 'rand' if (_C[key]==0) else 'nonrand'
                out_fp[readtype][aligntype].write(line)
    for readtype in ['target', 'control']:
        for aligntype in ['rand', 'nonrand']:
            if out_fp[readtype][aligntype]:
                out_fp[readtype][aligntype].close()

    # Convert output SAM to BAM files
    SAMtoBAM('realmapped_target_rand BAM to SAM', realmappedtargetrandsam_path, realmappedtargetrand_path)
    SAMtoBAM('realmapped_target_nonrand BAM to SAM', realmappedtargetnonrandsam_path, realmappedtargetnonrand_path)
    SAMtoBAM('realmapped_control_rand BAM to SAM', realmappedcontrolrandsam_path, realmappedcontrolrand_path)
    SAMtoBAM('realmapped_control_nonrand BAM to SAM', realmappedcontrolnonrandsam_path, realmappedcontrolnonrand_path)

    # Create output BAM index files
    BAMindex('realmapped_target_rand BAM index', realmappedtargetrand_path)
    BAMindex('realmapped_target_nonrand BAM index', realmappedtargetnonrand_path)
    BAMindex('realmapped_control_rand BAM index', realmappedcontrolrand_path)
    BAMindex('realmapped_control_nonrand BAM index', realmappedcontrolnonrand_path)

    # Remove intermediate files
    for path in [realmappedsam_path, realmappedtargetrandsam_path, realmappedtargetnonrandsam_path,
        realmappedcontrolrandsam_path, realmappedcontrolnonrandsam_path]:
        if os.path.exists(path):
            os.remove(path)

#def Extract_ReadDepth(inbam, refcontigid, outdepth, jobdesc):
#    cmd1 = '{samtools} depth -l 0 -q 0 -Q 0 {inbam}'.format(samtools=_PA.args.prog_samtools, inbam=controlnonrand_bam)
#    cmd2 = 'grep {refcontigid}'.format(refcontigid=refcontigid)
#    cmd1_h = sp.Popen(shlex.split(cmd1), stdout=PIPE)
#    cmd2_h = sp.Popen(shlex.split(cmd2), stdin=cmd1_h.stdout, stdout=PIPE)
#    ro, re = cmd2_h.communicate()
#    rc = cmd2_h.returncode
#    if rc != 0:
#         Exit('ErrorSysCall')
#    if len(ro):
#        if ro[-1]=='\n':
#            ro = ro[:-1]
#        with open(outdepth, 'w') as out_fp:
#            out_fp.write('{0}\n'.format(ro))
#    else:
#        open(outdepth, 'a').close()

def Extract_ReadDepth(goodbam_path, randbam_path, refcontigid, refcontiglen, depth_path):
    'Run samtools to extract the read depth for alignments to the specified refcontigid.'

    # Intermediate files
    gooddepth_path = depth_path + '.goodraw'
    randdepth_path = depth_path + '.randraw'

    # Read non-zero-depth sites from the good BAM file
    cmd = '{samtools} depth -l 0 -q 0 -Q 0 {inbam} | grep "{refcontigid}" > {outdepth}'.format(
        samtools=_PA.args.prog_samtools,
        inbam=goodbam_path,
        refcontigid=refcontigid,
        outdepth=gooddepth_path)
    os.system(cmd)
    if not os.path.exists(gooddepth_path):
        Exit('ErrorSysCall')
    gooddepth_data = np.loadtxt(gooddepth_path, dtype=[('refcontigid', '|S20'), ('pos', 'int'), ('depth', 'int')])

    # Read non-zero-depth sites from the rand BAM file
    cmd = '{samtools} depth -l 0 -q 0 -Q 0 {inbam} | grep "{refcontigid}" > {outdepth}'.format(
        samtools=_PA.args.prog_samtools,
        inbam=randbam_path,
        refcontigid=refcontigid,
        outdepth=randdepth_path)
    os.system(cmd)
    if not os.path.exists(gooddepth_path):
        Exit('ErrorSysCall')
    randdepth_data = np.loadtxt(randdepth_path, dtype=[('refcontigid', '|S20'), ('pos', 'int'), ('depth', 'int')])

    # Final table has columns: pos, rawgood, rawrand, rawboth, wingood, winrand, winboth
    depth_data = np.zeros(refcontiglen+1, dtype=[('pos', 'int'),
        ('rawgood', 'int'), ('rawrand', 'int'), ('rawboth', 'int'),
        ('wingood', 'float'), ('winrand', 'float'), ('winboth', 'float')])
    for i in range(1, refcontiglen+1):	        # pos
        depth_data[i][0] = i
    for i in range(0, len(gooddepth_data)):	# rawgood
        depth_data[gooddepth_data[i][1]][1] = gooddepth_data[i][2]
    for i in range(0, len(randdepth_data)):	# rawrand
        depth_data[randdepth_data[i][1]][2] = randdepth_data[i][2]
    for i in range(0, refcontiglen+1):		# rawboth
        depth_data[i][3] = depth_data[i][1] + depth_data[i][2]
    for i in range(0, len(gooddepth_data)):     # wingood
        depth_data[gooddepth_data[i][1]][4] = round(gooddepth_data[i][2], 3)
    for i in range(0, len(randdepth_data)):     # winrand
        depth_data[randdepth_data[i][1]][5] = round(randdepth_data[i][2], 3)
    for i in range(0, refcontiglen+1):          # winboth
        depth_data[i][6] = depth_data[i][4] + depth_data[i][5]

    #Save copy with missing zeros to a file.
    with open(depth_path, 'w') as out_fp:
        L = ['pos', 'rawgood', 'rawrand', 'rawboth', 'wingood', 'winrand', 'winboth']
        out_fp.write('{0}\n'.format('\t'.join(L)))
        for i in range(1, refcontiglen+1):
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in depth_data[i]])))

def Plot_Refcontig_ReadDepth(nonrandbam, randbam, refcontigid, refcontiglen, reftype, outdir, outprefix):
    '''
    Create one outprefix_reftype_refcontigid_alntype_readdepth.png file where
    reftype is "target" or "control".
    '''
    # Extract the depth information, one per 1-based position, with zeros,
    # with 7 columns: pos, rawgood, rawrand, rawboth, wingood, winrand, winboth
    depth_path = os.path.join(outdir,
        '{outprefix}_realmapped_{reftype}_{refcontigid}.depth'.format(
        outprefix=outprefix, reftype=reftype, refcontigid=refcontigid))
    Extract_ReadDepth(nonrandbam, randbam, refcontigid, refcontiglen, depth_path)

    # Generate the read depth plots
    outprefix = '{outprefix}_{reftype}_{refcontigid}'.format(
        outprefix=_PA.args.outprefix, reftype=reftype, refcontigid=refcontigid)
    cmd = '{rscript} {poremapcoverage} {runid} {refcontigid} {refcontiglen} {depthpath} ' \
        '{outdir} {outprefix}'.format(
        rscript=_PA.args.prog_rscript,
        poremapcoverage=os.path.join(_PA.args.bindir, 'poremapcoverage.R'),
        runid=_PA.args.runid,
        refcontigid=refcontigid,
        refcontiglen=len(_R[refcontigid][1]),
        depthpath=depth_path,
        outdir=_PA.args.outdir,
        outprefix=outprefix)
    rc, ro, re = Run_Shell_Script('poremapcoverage', cmd, dump_stdout=False)
    if rc != 0:
         Exit('ErrorSysCall')

@transform( \
    Get_SplitBams, \
    regex(r"(\S+)/(\S+)_realmapped_targetrand.bam"), \
    [r'{intdir}/\2_coverage_sentinel.txt'.format(intdir=_intdir)])
def Get_ReadDepth_Plots(inpaths, outpaths):
    '''
    Use samtools to get the read depth from the [target|control]_[nonrand|rand].bam files and
    output _[target|control]_REFCONTIGID_[nonrand|rand].png files.
    '''
    _PL.log_dumpout('Info: Get_ReadDepth_Plots: inpaths={0}, outpaths={1}'.format(inpaths, outpaths))

    # Read in alignclass file
    Read_References()
    Read_AlignClassFile()

    # IO paths
    targetnonrand_bam = os.path.join(_intdir, _PA.args.outprefix+'_realmapped_targetnonrand.bam')
    targetrand_bam = os.path.join(_intdir, _PA.args.outprefix+'_realmapped_targetrand.bam')
    controlnonrand_bam = os.path.join(_intdir, _PA.args.outprefix+'_realmapped_controlnonrand.bam')
    controlrand_bam = os.path.join(_intdir, _PA.args.outprefix+'_realmapped_controlrand.bam')

    # Extract a _nonrand.depth and _rand.depth file for each refcontigid
    # Create one PNG using for each _[nonrand|rand].depth pair using poremapcoverage.R
    target_refcontigidL = [refcontigid for refcontigid in _R.keys() if _R[refcontigid][0] == 'target']
    for refcontigid in target_refcontigidL:
        refcontiglen = len(_R[refcontigid][1])
        Plot_Refcontig_ReadDepth(targetnonrand_bam, targetrand_bam, refcontigid, refcontiglen, 'target', _intdir, _PA.args.outprefix)
    control_refcontigidL = [refcontigid for refcontigid in _R.keys() if _R[refcontigid][0] == 'control']
    for refcontigid in control_refcontigidL:
        refcontiglen = len(_R[refcontigid][1])
        Plot_Refcontig_ReadDepth(controlnonrand_bam, controlrand_bam, refcontigid, refcontiglen, 'control', _intdir, _PA.args.outprefix)

    # Create sentinel file
    sentinel_path = os.path.join(_intdir, _PA.args.outprefix+'_coverage_sentinel.txt')
    open(sentinel_path, 'a').close()

# ============================================================================ #
# Main                                                                         #
# ============================================================================ #

if __name__ == '__main__':

    task_list = [
        References_Setup,
        Get_Alignments_Initial,
        Get_Alignments_Training_ControlReal,
        Get_Alignments_Training_ControlRand,
        Get_Alignments_Testing_TargetReal,
        Get_Alignment_TargetReal_Classifications,
        Get_Alignment_Statistics,
        Get_SplitBams,
        Get_ReadDepth_Plots
    ]
    if _PA.args.ruffusgraph:
        pipeline_printout_graph(os.path.join(_PA.args.outdir, "ruffus_graph.png"), "png", task_list)
    elif _PA.args.touch:
        pipeline_run(task_list, touch_files_only=True)
    elif _PA.args.dryrun:
        pipeline_printout(sys.stdout, task_list, verbose=6)
    elif _PA.args.overwrite:
        pipeline_run(task_list, forcedtorun_tasks=task_list)
    else:
        pipeline_run(task_list)
    Exit('SuccessReturn')

# ============================================================================ #
