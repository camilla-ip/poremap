#!/usr/bin/env python

# ============================================================================ #
# porelog.py                                                                   #
# ============================================================================ #
'''
Logging functions for poremap programs
'''
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# April 2015                                                                   #
# ============================================================================ #

import os, sys, time

class Porelog(object):
    'A class for poremap program diagnostic and error message logging.'

  # ============================================================================ #
  # Data                                                                         #
  # ============================================================================ #

    def __init__(self, progdoc='', logoutdir=None, logoutpath=None, logoutfileend=None, logversion=None, bindir=None, ini_version=None):

        self._prog = os.path.basename(sys.argv[0])
        self._progdesc = progdoc.strip().split('\n')[0].strip()
        self._bindir = bindir if bindir is not None else './'
        self._ini_version = ini_version

        self._repo_info = self._src_info()

        self._starttime = time.localtime()
        if (logoutdir is None and logoutpath is None and logoutfileend is None):
            path = os.path.join('$ONT', 'log', time.strftime('%Y', self._starttime),
                time.strftime('%m', self._starttime),
                '{0}.log'.format(time.strftime('%Y%m%d%H%M%S', self._starttime)))
        elif (logoutdir is not None):
            path = os.path.join(logoutdir,
                '{0}.log'.format(time.strftime('%Y%m%d%H%M%S', self._starttime)))
        elif (logoutpath is not None):
            path = logoutpath
        elif (logoutfileend is not None):
            path = os.path.join('$ONT', 'log', time.strftime('%Y', self._starttime),
                time.strftime('%m', self._starttime),
                '{0}_{1}.log'.format(time.strftime('%Y%m%d%H%M%S', self._starttime), logoutfileend))
        self._fpath = os.path.expandvars(path)
        self._fpth = None
        self._fpathgz = os.path.expandvars(path + '.gz')
        self._filehasbeenopened = False

        self._fpath = None
        self._fpth = None
        self._fpathgz = None
        self._filehasbeenopened = False
        self.log_setoutpath(logoutdir, logoutpath, logoutfileend)

  # ============================================================================ #
  # Methods                                                                      #
  # ============================================================================ #

    def _src_info(self):
        """
        Return all the information from the "svn info ont.ini:mmmrep" command as a dictionary.
        If call to svn repository does not work, then try to read the repinfo.ini file
        from the same directory as the ont.ini file (which may have been set by the
        make install program). If that fails, return the default "notfound" strings
        for each of the fields.
        """
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
        #    cmd = '{prog} info {url}'.format(prog=self.ini_val('prog_svn'), url=os.path.join(self.ini_val('mmmrep'), 'src'))
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
        """Return just the rN part of the revision information retrieved from svn log."""
        return 'r{revisionnumber}'.format(revisionnumber=self._repo_info['Revision'])

    def src_revisiondatetime(self):
        """Return a string of the format (rN YYYY-MM-DD HH:MM:SS) indicating the revision and last edit date+time."""
        return '{revision} {lastchangeddate}'.format(
            revision=self.src_revision(),
            lastchangeddate=' '.join(self._repo_info['LastChangedDate'].split(' ')[0:2]))

    def _progdescline(self):
        result = "Program: {progname}".format(progname=os.path.basename(self._prog))
        if (len(self._progdesc)):
            result += " ({progdesc})".format(progdesc=self._progdesc)
        return result

    def _progversionline(self):
        """Return a string listing the version (SVN revision) of the calling program."""
        if self._repo_info['Revision'] != 'NOTFOUND':
            return "Version: {0} ({1})".format(self._ini_version, self.src_revision())
        else:
            return "Version: {0}".format(self._ini_version)

    def _msgprefix(self):
        'Assemble prefix for all messages of format "YYYYMMDD-HHMMSS (progname): ".'
        return '{datetime} [{progname}]: '.format(datetime=time.strftime('%Y%m%d-%H%M%S', time.localtime()), progname=self._prog)

    def _createdirsandopenfile(self):
        'The first time a dump has been called, create intermediate dirs and open the output file.'
        if (not os.path.exists(self._fpath) and not os.path.exists(self._fpathgz)):
            outdir = os.path.dirname(self._fpath)
            if (not os.path.exists(outdir)):
                try:
                    os.makedirs(os.path.dirname(self._fpath), mode=0775)
                except OSError, e:
                   #sys.stderr.write('Error: Cannot create intermediate dirs for log file [{0}: {1}]\n'.format(
                   #    e.code, e.errorcode[e.errno]))
                   sys.stderr.write('Error: Cannot create intermediate dirs for log file [{0}: {1}]\n'.format(
                       e.errno, e.strerror))
            try:
                self._fptr = open(self._fpath, 'a') if not self._fpath.endswith('.gz') else gzip.open(self._fpath, 'ab')
                if (self._fptr):
                    self._filehasbeenopened = True
                self._fptr.write('{0}Starting {1}\n'.format(self._msgprefix(), ' '.join(sys.argv)))
            except OSError, e:
               #sys.stderr.write('Error: Cannot open log file in write mode [{0}: {1}]\n'.format(
               #    e.code, e.errorcode[e.errno]))
               sys.stderr.write('Error: Cannot open log file in write mode [{0}: {1}]\n'.format(
                   e.errno, e.strerror))
        else:
            sys.stderr.write('Error: Should not be writing to existing log file [{0}(.gz)]\n'.format(self._fpath))
            self.err_tidy('ErrorExistingLogFile', exitprogram=True)

    def str_indented(self, s, indentwidth=8):
        'Return the multi-line string with indentwidth leading spaces.'
        return '\n'.join(['{0}{1}'.format(' '*indentwidth, x) for x in s.split('\n')])

    def log_header(self):
        'A nicely formatted string that would go well at the top of a log file.'
        pline = "{datetime} [{progname}]: {text}".format(
            datetime=time.strftime('%Y%m%d-%H%M%S', time.localtime()),
            progname=self._prog,
            text=self._progdescline())
        vline = "{datetime} [{progname}]: {text}".format(
            datetime=time.strftime('%Y%m%d-%H%M%S', time.localtime()),
            progname=self._prog,
            text=self._progversionline())
        cline = "{datetime} [{progname}]: Command: {text}".format(
            datetime=time.strftime('%Y%m%d-%H%M%S', time.localtime()),
            progname=self._prog,
            text=' '.join([os.path.expandvars(sys.argv[0])] + sys.argv[1:]))
        return '\n'.join([pline, vline, cline])

    def log_runinfo(self):
        '''
        Return the multi-line string that can be printed to the top of a log file
        using something like _GT.log_2out(_GT.log_runinfo()).
        '''
        return '\n{0}'.format('\n'.join(['{0}'.format(x) for x in self.log_header().split('\n')]))

    def log_setoutpath(self, logoutdir=None, logoutpath=None, logoutfileend=None):
        'Set the ontlib private class variables for the output dir for the log_ functions.'
        if (self._filehasbeenopened):
            try:
                self._fpth.close()
                self._fpath = None
                self._fpth = None
                self._fpathgz = None
                self._filehasbeenopened = False
            except:
                pass
        if (logoutdir is None and logoutpath is None and logoutfileend is None):
            path = os.path.join('$ONT', 'log', time.strftime('%Y', self._starttime),
                time.strftime('%m', self._starttime),
                '{0}.log'.format(time.strftime('%Y%m%d%H%M%S', self._starttime)))
        elif (logoutdir is not None):
            path = os.path.join(logoutdir,
                '{0}.log'.format(time.strftime('%Y%m%d%H%M%S', self._starttime)))
        elif (logoutpath is not None):
            path = logoutpath
        elif (logoutfileend is not None):
            path = os.path.join('$ONT', 'log', time.strftime('%Y', self._starttime),
                time.strftime('%m', self._starttime),
                '{0}_{1}.log'.format(time.strftime('%Y%m%d%H%M%S', self._starttime), logoutfileend))
        self._fpath = os.path.expandvars(path)
        self._fpth = None
        self._fpathgz = os.path.expandvars(path + '.gz')
        self._filehasbeenopened = False

    def log_dump(self, str):
        'Print the text message to the current logfile.'
        if (not self._filehasbeenopened):
            self._createdirsandopenfile()
        if self._fptr:
            self._fptr.write('{0}{1}\n'.format(self._msgprefix(), str))
            self._fptr.flush()

    def log_dumpout(self, str):
        'Print the text message to the current logfile and to stdout.'
        if (not self._filehasbeenopened):
            self._createdirsandopenfile()
        self.log_dump(str)
        sys.stdout.write('{0}{1}\n'.format(self._msgprefix(), str))
        sys.stdout.flush()

    def log_dumperr(self, str):
        'Print the text message to the current logfile and to stdout.'
        if (not self._filehasbeenopened):
            self._createdirsandopenfile()
        self.log_dump(str)
        sys.stderr.write('{0}{1}\n'.format(self._msgprefix(), str))
        sys.stderr.flush()

    def log_2out(self, str):
        'Print the text message to the current logfile and to stdout.'
        sys.stdout.write('{0}{1}\n'.format(self._msgprefix(), str))
        sys.stdout.flush()

    def log_2err(self, str):
        'Print the text message to the current logfile and to stdout.'
        sys.stderr.write('{0}{1}\n'.format(self._msgprefix(), str))
        sys.stderr.flush()

    def log_tidy(self):
        'Print final log message, close the log file and gzip it.'
        if self._filehasbeenopened:
            try:
                self._fptr.write('{0}Finishing {1}\n'.format(self._msgprefix(), ' '.join(sys.argv)))
                self._fptr.flush()
                self._fptr.close()
            except OSError, err:
                sys.stderr.write('Error: Failed to close log file [{0}: {1}]\n'.format(
                    e.code, e.errorcode[e.errno]))
        self._fptr = None
        if os.path.exists(self._fpath) and not self._fpath.endswith('.gz'):
            gzipcmd = 'gzip {0}'.format(self._fpath)
            os.system(gzipcmd)

    def log_exit(self, ret_code, err_text, print2out=True, print2err=False, print2dump=False, exitprogram=False):
        """
        Function to call just before exiting program.
        """
        if print2out:
            self.log_2out(err_text)
        if print2err:
            self.log_2out(err_text)
        if print2dump:
            self.log_dump(err_text)
        self.log_tidy()
        if exitprogram:
            sys.exit(ret_code)

# ============================================================================ #
