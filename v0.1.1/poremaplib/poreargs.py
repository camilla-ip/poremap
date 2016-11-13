#!/usr/bin/env python

# ============================================================================ #
# poreargs.py                                                                  #
# ============================================================================ #
'''
Argument parsing class shared by poremap programs
'''
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# April 2015                                                                   #
# ============================================================================ #

from subprocess import PIPE
import argparse, getpass, os, shlex, socket, stat, subprocess as sp, sys, time

class Poreargs(object):
    'A class for poremap program argument parsing.'

    def __init__(self, progdoc='', progexamples=[], progaddhelp=True, progusagemsg='', bindir=None, version=None):

      # ============================================================================ #
      # Data                                                                         #
      # ============================================================================ #

      # Ont ontlib data - that does not depend on anything else
        self._progdoc = '\n'.join(progdoc.strip().split('\n')[1:]).strip()
        self._progdesc = progdoc.strip().split('\n')[0].strip()
        self._exec_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
        self._hostname = socket.gethostname()
        self._userid = getpass.getuser()
        self._bindir = bindir if bindir is not None else './'
        self._version = version if version is not None else ''

      # Ont src data
        self._filepath = os.path.expandvars(sys.argv[0])
        self._repo_info = self._src_info()


        self.args = None

        # Call all the base class functions to set up data for Ont.
        # Do it here instead of in the base class so that the base classes have no dependencies.
        #  1. flag        2. metavar     3. action      4. dest       5. type
        #  6. default     7. nargs       8. choices     9. help
        self._arg_stdarglist =  [

        [ '--bindir', 'DIR', 'store', 'bindir', dir, None, None, '',
            'Absolute path of directory containing this program and .ini and .profile config files' ],
        [ '--profilepath', 'STR', 'store', 'profilepath', str, None, None, '',
            'Path to BASH environment variables path' ],
        [ '--runid', 'STR', 'store', 'runid', str, None, None, '',
             'RunID, which should not contain any spaces or unusual punctuation' ],
        [ '--readtype', 'STR', 'store', 'readtype', str, None, None, ['2d', 'temp', 'comp', 'mixed', 'unknown'],
             'Type of ONT reads' ],
        [ '--readclass', 'STR', 'store', 'readclass', str, None, None, ['all', 'pass', 'fail'],
             'Class of ONT reads' ],
        [ '--datatype', 'STR', 'store', 'datatype', str, 'minion', None, ['minion', 'rand', 'mixed', 'unknown'],
             'Type of run data' ],
        [ '--infastqpath', 'FILE', 'store', 'infastqpath', file, None, None, '', 'Input reads path (FASTQ)' ],
        [ '--readsbam', 'FILE', 'store', 'readsbam', file, None, None, '', 'Uniquely-named, single-end reads (BAM)' ],
        [ '--targetrefpath', 'FILE', 'store', 'targetrefpath', file, None, None, '',
            'Path to FASTA file for extracting reads of target genome' ],
        [ '--controlrefpath', 'FILE', 'store', 'controlrefpath', file, None, None, '',
            'Reference control sample(s) (FASTA)' ],
        [ '--alignclasspath', 'FILE', 'store', 'alignclasspath', file, None, None, '',
            'Path to the alignclass file containing the isnonrandomaln information' ],
        [ '--outdir', 'DIR', 'store', 'outdir', dir, None, None, '', 'Output directory' ],
        [ '--outprefix', 'STR', 'store', 'outprefix', str, None, None, '', 'Prefix for output files' ],
        [ '--useintdir', 'BOOL', 'store', 'useintdir', str, False, None, '',
            'Run pipeline to create output files in an intermediate subdir first, then move them to the ' \
            'final directory if the pipeline completes without error' ],
        [ '--threads', 'INT', 'store', 'threads', str, 1, None, '', 'Maximum number of threads to use' ],
        [ '--mapprog', 'STR', 'store', 'mapprog', str, None, None, ['bwa', 'graphmap'],
             'Mapping program to use for the target and/or control reads: bwa, graphmap' ],
        [ '--mapparams', 'STR', 'store', 'mapparams', str, None, None, '',
             'Additional parameters to pass to --mapprog surrounded by double-quotes if it contains spaces' ],
        [ '--realmapprog', 'STR', 'store', 'realmapprog', str, None, None, ['bwa', 'graphmap'],
             'Mapping program to use for the real ONT target and/or control reads: bwa, graphmap' ],
        [ '--realmapparams', 'STR', 'store', 'realmapparams', str, None, None, '',
             'Additional parameters to pass to --realmapprog surrounded by double-quotes if it contains spaces' ],
        [ '--randmapprog', 'STR', 'store', 'randmapprog', str, None, None, ['bwa', 'graphmap'],
             'Mapping program to use for the randomised ONT control reads: bwa, graphmap' ],
        [ '--randmapparams', 'STR', 'store', 'randmapparams', str, None, None, '',
             'Additional parameters to pass to --randmapprog surrounded by double-quotes if it contains spaces' ],
        [ '--savealignments', 'BOOL', 'store', 'savealignments', str, False, None, '',
             'Save each alignment to a separate FASTA file for diagnostic purposes' ],
        [ '--dryrun', 'BOOL', 'store', 'dryrun', str, False, None, '',
            'Run without creating any output or changing any states' ],
        [ '--touch', 'BOOL', 'store', 'touch', str, False, None, '',
             'Force Ruffus to touch all files only, making them seem up to date' ],
        [ '--ruffusgraph', 'BOOL', 'store', 'ruffusgraph', str, False, None, '',
             'Produce Ruffus dependency graph only, without running pipeline' ],
        [ '--printinfo', 'BOOL', 'store', 'printinfo', str, False, None, '',
            'Print a description of each of the output file fields' ],
        [ '--overwrite', 'BOOL', 'store', 'overwrite', str, False, None, '',
             'Print warning and overwrite any output files that may already exist' ],
        [ '--verbose', 'BOOL', 'store', 'verbose', str, False, None, '',
            'Verbose diagnostic output' ],
        [ '--deleteintfiles', 'BOOL', 'store', 'deleteintfiles', str, False, None, '',
            'Delete intermediate output files before exiting' ],
        [ '--logindentwidth', 'INT', 'store', 'logindentwidth', int, 4, None, '',
            'Print log messages from each new spawned process indented by this many spaces' ],
        [ '--covplotwinsz', 'INT', 'store', 'covplotwinsz', int, 100, None, '',
            'Print log messages from each new spawned process indented by this many spaces' ],
        [ '--usecbcolours', 'BOOL', 'store', 'usecbcolours', str, False, None, '',
            'Use colour-blind-friendly colours' ],
        [ '--fastalinewidth', 'INT', 'store', 'fastalinewidth', int, 100, None, '',
            'Number of characters in each line of the output FASTA file(s)' ],
        [ '--prog_bwa', 'STR', 'store', 'prog_bwa', str, 'bwa', None, '',
            'BWA program (absolute pathname)' ],
        [ '--prog_graphmap', 'STR', 'store', 'prog_graphmap', str, 'graphmap', None, '',
            'GraphMap program (absolute pathname)' ],
        [ '--prog_samtools', 'STR', 'store', 'prog_samtools', str, 'samtools', None, '',
            'SAMtools program (absolute pathname)' ],
        [ '--prog_rscript', 'STR', 'store', 'prog_rscript', str, 'Rscript', None, '',
            'Rscript program (absolute pathname)' ]
        ] # END

      # Ont arg data, passed to Ont class
        self._prog = os.path.basename(sys.argv[0])
        self._progexamples = progexamples

      # Ont arg data
        # formatter_class choices are: RawDescriptionHelpFormatter, RawText, ArgumentDefaults, MetavarType
        if len(progusagemsg):
            self._parser = argparse.ArgumentParser(
               #description=self._progdoc,
               #epilog='\n'.join(['examples:'] + ['  '+x for x in progexamples]),
                add_help=progaddhelp,
                prog=self._prog,
                usage=progusagemsg,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        else:
            self._parser = argparse.ArgumentParser(
               #description=self._progdoc,
               #epilog='\n'.join(['examples:'] + ['  '+x for x in progexamples]),
                add_help=progaddhelp,
                prog=self._prog,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self._arg_orderedlist = []      # _arg_orderedlist = [[optionname, variablename, metavar], ...]
                                        # e.g., ['--prog_annotvcf', 'ont_prog_annotvcf' ]
                                        # It is useful to know the order in which these were
                                        # specified in the calling program so we can iterate
                                        # through them in order later.

    # ============================================================================ #
    # Methods                                                                      #
    # ============================================================================ #

    def sys_makedir(self, dir):
        """
        Create dir with file permissions from ontdirmode in the ini file, if necessary.
        """
        dirpath = os.path.expandvars(dir)
        if dirpath and not os.path.exists(dirpath):
            try:
                os.makedirs(dirpath, mode=stat.S_IRWXU|stat.S_IRGRP|stat.S_IROTH)
            except:
                return False
        if not os.path.exists(dirpath):
            return False
        return True

    def sys_exec(self, cmd):
        """
        Execute a command using the subprocess module to trap the
        return value, the stdout string and stderr string.
        """
        proc_handle = sp.Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
        proc_stdout, proc_stderr = proc_handle.communicate()
        proc_returncode = proc_handle.returncode
        return [proc_returncode, proc_stdout, proc_stderr]

    def _src_info(self):
        '''
        Return all the information from the "svn info ont.ini:mmmrep" command as a dictionary.
        If call to svn repository does not work, then try to read the repinfo.ini file
        from the same directory as the ont.ini file (which may have been set by the
        make install program). If that fails, return the default "notfound" strings
        for each of the fields.
        '''
        info = {
            'Path' : 'notfound',
            'URL' : 'notfound',
            'RepositoryRoot' : 'notfound',
            'RepositoryUUID' : 'notfound',
            'Revision' : 'NOTFOUND',
            'NodeKind' : 'notfound',
            'LastChangedAuthor' : 'notfound',
            'LastChangedRev' : 'notfound',
            'LastChangedDate' : '0000-00-00 00:00:00'
        }
        found = False

        svninfopath = os.path.join(self._bindir, 'repinfo.ini')
        if os.path.exists(svninfopath):
            try:
                cmd = 'cat {svninfopath}'.format(svninfopath=svninfopath)
                proc_returncode, proc_stdout, proc_stderr = self.sys_exec(cmd)
                if proc_returncode != 0:
                    return info
                L = [x for x in proc_stdout.split('\n') if len(x)]
                D = {}
                for elt in L:
                    var = elt.split(':')[0].replace(' ', '')
                    val = ':'.join(elt.split(':')[1:]).strip()
                    D[var] = val
                info = copy.deepcopy(D)
                found = True
            except:
                pass

        #if not found:
        #    cmd = '{prog} info {url}'.format('svn', url=os.path.join(self.ini_val('mmmrep'), 'src'))
        #    proc_returncode, proc_stdout, proc_stderr = self.sys_exec(cmd)
        #    if proc_returncode == 0:
        #        try:
        #            L = [x for x in proc_stdout.split('\n') if len(x)]
        #            D = {}
        #            for elt in L:
        #                var = elt.split(':')[0].replace(' ', '')
        #                val = ':'.join(elt.split(':')[1:]).strip()
        #                D[var] = val
        #            info = copy.deepcopy(D)
        #            found = True
        #        except:
        #            pass

        return info

    def src_revision(self):
        '''Return just the rN part of the revision information retrieved from svn log.'''
        result = ''
        if len(self._repo_info['Revision']) and self._repo_info['Revision'].upper() != 'NOTFOUND':
            result = 'r{revisionnumber}'.format(revisionnumber=self._repo_info['Revision'])
        return result

    def src_revisiondatetime(self):
        '''Return a string of the format (rN YYYY-MM-DD HH:MM:SS) indicating the revision and last edit date+time.'''
        return '{revision} {lastchangeddate}'.format(
            revision=self.src_revision(),
            lastchangeddate=' '.join(self._repo_info['LastChangedDate'].split(' ')[0:2]))

    def arg_setstdlist(self, stdarglist):
        'Set the class variable storing the set of standard arguments.'
        self._arg_stdarglist = stdarglist

    def arg_addnew(self, name, **kwargs):
        '''
        Wrapper to the argparse add_argument function.

        The value flags are:

        name or flags - Either a name or a list of option strings, e.g. foo or -f, --foo.
        action - The basic type of action to be taken when this arg is encountered at the command line.
        nargs - The number of command-line arguments that should be consumed.
        const - A constant value required by some action and nargs selections.
        default - The value produced if the argument is absent from the command line.
        type - The type to which the command-line argument should be converted.
        choices - A container of the allowable values for the argument.
        required - Whether or not the command-line option may be omitted (optionals only).
        help - A brief description of what the argument does.
        metavar - A name for the argument in usage messages.
        dest - The name of the attribute to be added to the object returned by parse_args().
        '''
        self._parser.add_argument(name, **kwargs)
        metavar = kwargs['metavar'] if kwargs.has_key('metavar') else 'STR'
        self._arg_orderedlist.append([name, kwargs['dest'], metavar])

    #  1. name or flag 2. metavar     3. action      4. dest       5. type
    #  6. default      7. nargs       8. choices     9. help
    def arg_addstd(self, name, **kwargs):
        'Add a standard argument with consistent help message and default.'
        elt = [x for x in self._arg_stdarglist if x[0] == name]
        info = elt[0] if len(elt) else []
        if not len(elt) or not len(info):
            sys.stderr.write('Error: Onttools::arg_addstd(): Unrecognised standard argument in calling program *{0}*\n'.format(name))
            sys.exit(3) # ErrorArgsInvalid
        params = {}
        params['required'] = kwargs['required'] if kwargs.has_key('required') else False
        if (kwargs.has_key('choices') and len(kwargs['choices'])):
            params['metavar'] = '{0}{1}{2}'.format('{', ','.join(kwargs['choices']), '}')
        elif (info[7]):
            params['metavar'] = '{0}{1}{2}'.format('{', ','.join(info[7]), '}')
        elif kwargs.has_key('metavar') and kwargs.has_key('metavar'):
            params['metavar'] = '{0}{1}{2}'.format('{', ','.join(kwargs['metavar']), '}')
        elif info[1]:
            params['metavar'] = info[1]
        metavartype = kwargs['metavar'] if kwargs.has_key('metavar') else info[1]
        params['action'] = kwargs['action'] if kwargs.has_key('action') else info[2]
        params['dest'] = kwargs['dest'] if kwargs.has_key('dest') else info[3]
        params['default'] = kwargs['default'] if kwargs.has_key('default') else info[5]
        params['nargs'] = kwargs['nargs'] if kwargs.has_key('nargs') else info[6]
        if kwargs.has_key('choices'):
            params['choices'] = kwargs['choices']
        elif info[7] and len(info[7]):
            params['choices'] = info[7]
       #params['choices'] = kwargs['choices'] if kwargs.has_key('choices') else info[7]
        params['help'] = kwargs['help'] if kwargs.has_key('help') else info[8]
        self._parser.add_argument(name, **params)
        self._arg_orderedlist.append([name, params['dest'], metavartype])
        return True

    def arg_parse(self):
        'Function that prints usage, help or version messages as required.'
        # Wrapper to the parse_args fn from the Python arg module.'
        # This function should not be called from a program. Should use
        # Onttools.parseall() instead.
        self.args = self._parser.parse_args()
        # Now have to fix all the variables of type bool to pick out the first element of this type

    def arg_printhelp(self, fstream=sys.stdout):
        'Print a longer help message for this program.'
       #try:
       #    helpwidth = int(os.environ['COLUMNS'])
       #except:
       #    helpwidth = 80
       #indentwidth = 2
       #helpwidth -= indentwidth
      # Get the original help text from the Python argparse module
        msg = self._progsummaryline() + '\n' + self._parser.format_help()
      # Replace the description with the docstring with the original line breaks.
      # desc_stripped = self._parser.description.strip()
      # desc_line1 = desc_stripped.split('\n')[0].strip()
      # desc_lineN = desc_stripped.split('\n')[-1].strip()
      # newdesc = '  ' + '\n  '.join(self._parser.description.strip().split('\n'))
      # sttidx = msg.find(desc_line1)
      # endidx = msg.find(desc_lineN)+len(desc_lineN)
      # if (sttidx <= endidx):
      #     msg = msg.replace(msg[sttidx:endidx], newdesc)
      # Fix up minor formatting issues
        msg = msg.replace('program:', 'Program: ')
        msg = msg.replace('version:', 'Version: ')
        msg = msg.replace('usage:', '\nUsage:\n')
        msg = msg.replace('positional arguments:', 'Positional Arguments:\n')
        msg = msg.replace('optional arguments:', 'Optional Arguments:\n')
       #msg = msg.replace('examples:\n', '')
        msg = msg.replace(' {0}'.format(self._prog), '\n  {0}'.format(self._prog))
        msg += '\n'
      # Add the examples lines
       #msg += 'Examples:\n'
       #msg = msg.replace('examples:', 'Examples:\n')
        msg += '\n'.join(['Examples:\n'] + ['  '+x for x in self._progexamples])
        msg += '\n'
      # Add all the description paragraphs
        msg += '\nDescription:\n\n'
        msg += '\n'.join([x if x.endswith(':') else ('  {0}'.format(x)) for x in self._progdoc.split('\n')])
        msg += '\n\n'
      # Print the whole msg
        fstream.write(msg)

    def arg_printusage(self, fstream=sys.stdout):
        'Print a brief usage message for this program.'
        msg = self._progsummaryline() + '\n\n' + self._parser.format_usage() + '\n'
        msg = msg.replace('usage: ', 'Usage:\n\n  ')
        fstream.write(msg)

    def arg_validate(self):
        'Do additional validation checks on the argument values.'

    def arg_progusage(self):
        '''Return a program usage text which is split over multiple lines.'''
       #return '\n'.join([self._progdescline(), self._progversionline(), "", self._progusageline(), ""])
        return self._progsummaryline()

    def _progsummaryline(self):
        '''Return the program name, description, version and revision string.'''
        result = '{progname}'.format(progname=self._prog)
        if len(self._progdesc):
            result += ' : {progdesc}'.format(progdesc=self._progdesc)
        if len(self.src_revision()):
            result += ' ({revision})'.format(revision=self.src_revision())
        return result

    def arg_progusage(self):
        '''Return a program usage text which is split over multiple lines.'''
       #return '\n'.join([self._progdescline(), self._progversionline(), "", self._progusageline(), ""])
        return self._progsummaryline()

    def str2bool(self, s):
        """Return the Python boolean value True or False dependingon what the string is."""
        return (s is not None and s.lower() in ['1', 't', 'true', 'y', 'yes'])

    def arg_fixtypes(self):
        '''
        Iterate through self.args.__dict__ and if the key has metavar 'BOOL' in self._arg_orderedlist
        then use str2bool() to change the type of the value to boolean.
        '''
        for varname in self.args.__dict__.keys():
            try:
                varinfo = [x for x in self._arg_orderedlist if x[1] == varname]
                #if varinfo and varinfo[0][2].upper() == 'BOOL' and type(self.args.__dict__[varname]) is str:
                if varinfo and varinfo[0][2].upper() == 'BOOL':
                    self.args.__dict__[varname] = self.str2bool(self.args.__dict__[varname])
            except:
                #sys.stderr.write('Erro: in ontlib.py arg_fixtypes() when trying to fix varname=*{varname}*\n'.format(varname=varname))
                pass

    def _progdescline(self):
        result = "Program: {progname}".format(progname=os.path.basename(self._prog))
        if (len(self._progdesc)):
            result += " ({progdesc})".format(progdesc=self._progdesc)
        return result

    def _progusageline(self):
        """Return a string containing the usage pattern for the calling program."""
        result = "Usage:   {0}".format(os.path.basename(sys.argv[0]))
        if (len(self._progusagepattern)):
            result += " " + self._progusagepattern
        return result

    def _progversionline(self):
        """Return a string listing the version (SVN revision) of the calling program."""
        return "Version: {0} ({1})".format( self.ini_version(), self.src_revision())

    def progversionrevision(self):
        """This is the program basename and the current version from the ont.ini file."""
        return '{0} version {1}'.format(os.path.basename(sys.argv[0]), self._version)

    def arg_parseall(self):
        '''Function that should be called to parse all the command-line arguments.'''
        usagestrings = set([ 'usage', '-usage', '--usage' ])
        helpstrings = set([ 'help', '-h', '-help', '--help', '?' ])
        versionstrings = set([ '-v', '--version', '-about', '--about' ])
        lowerargs = set([x.lower() for x in sys.argv])
        if (len(sys.argv) == 1 or not lowerargs.isdisjoint(usagestrings)):
            self.arg_printusage()
            sys.exit(11) # SuccessUsage
        if (not lowerargs.isdisjoint(helpstrings)):
            self.arg_printhelp()
            sys.exit(1) # SuccessHelp
        if (not lowerargs.isdisjoint(versionstrings)):
            print self.progversionrevision()
            sys.exit(12) # SuccessVersion
        try:
            self.arg_parse()
            self.arg_fixtypes()
        except:
            sys.exit(4) # ErrorArgsParseFailure

    def arg_valuelines(self, indentwidth=8):
        '''Return a string of key=value pairs, one per line, of the program arguments after parsing.'''
        indent = ' ' * indentwidth
        s = '{indent}Options:\n'.format(indent=indent)
        for pair in self._arg_orderedlist:
            s += '{indent}{key} {val}\n'.format(indent=indent, key=pair[0], val=self.args.__dict__[pair[1]])
        return s

    def arg_options(self):
        '''Return the list of the [key,value] pairs of the command-line arguments.'''
        L = []
        for pair in self._arg_orderedlist:
            L.append([pair[0], self.args.__dict__[pair[1]]])
        return L

# ============================================================================ #
# Class tests                                                                  #
# ============================================================================ #

if __name__ == "__main__":

    pass

# ============================================================================ #
