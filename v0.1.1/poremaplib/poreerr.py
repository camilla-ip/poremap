#!/usr/bin/env python

# ============================================================================ #
# poreerr.py                                                                   #
# ============================================================================ #
'''
Return value functions for poremap programs
'''
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# April 2015                                                                   #
# ============================================================================ #

import os, sys, time

class Poreerr(object):
    'A class for poremap program for consistent return value codes.'

  # ============================================================================ #
  # Data                                                                         #
  # ============================================================================ #

    def __init__(self):

        self.err = {
        'SuccessReturn'         : [  0, 'Info', 'Successfully completed' ],
        'SuccessHelp'           : [  1, 'Info', 'Successfully printed help message' ],
        'ErrorArgsMissing'      : [  2, 'Erro', 'Failed to retrieve expected command-line argument(s)' ],
        'ErrorArgsInvalid'      : [  3, 'Erro', 'Failed due to invalid command-line argument(s)' ],
        'ErrorArgsParseFailure' : [  4, 'Erro', 'Failed to parse command-line argument(s)' ],
        'ErrorSagaConnect'      : [  5, 'Erro', 'Failed to connect to SCAMPI db' ],
        'ErrorSagaDisconnect'   : [  6, 'Erro', 'Failed to disconnect from SCAMPI db' ],
        'ErrorFileOpen'         : [  7, 'Erro', 'Failed to open file' ],
        'ErrorSysCall'          : [  8, 'Erro', 'Failed system call' ],
        'ErrorFork1Failure'     : [  9, 'Erro', 'Failed to do first process fork' ],
        'ErrorFork2Failure'     : [ 10, 'Erro', 'Failed to do second process fork' ],
        'SuccessUsage'          : [ 11, 'Info', 'Successfully printed usage message' ],
        'SuccessVersion'        : [ 12, 'Info', 'Successfully printed version message' ],
        'ErrorDirCreate'        : [ 13, 'Erro', 'Failed to create directory' ],
        'ErrorExistingLogFile'  : [ 14, 'Erro', 'Failed because log file already exists' ],
        'ErrorDirDoesNotExist'  : [ 15, 'Erro', 'Failed because output directory does not exist' ],
        'SuccessNothingToDo'    : [ 16, 'Info', 'Successfully terminated because nothing to do' ],
        'ErrorIniFileMissing'   : [ 17, 'Erro', 'Failed to find ont.ini file' ],
        'ErrorInvalidData'      : [ 18, 'Erro', 'Failed due to invalid input data' ],
        'ErrorExternalSoftware' : [ 19, 'Erro', 'Failed in call to external software' ],
        'SuccessDryRun'         : [ 20, 'Info', 'Successfully terminated in dryrun mode' ],
        'SuccessKilled'         : [ 21, 'Info', 'Successfully terminated by KILL signal' ],
        'ErrorKilled'           : [ 22, 'Erro', 'Failed because terminated by KILL signal' ],
        'ErrorFileMissing'      : [ 23, 'Erro', 'Failed because file is missing' ],
        'ErrorChmod'            : [ 24, 'Erro', 'Failed to change mode of file or directory' ],
        'ErrorScriptCompute'    : [ 25, 'Erro', 'Failed in spawned script' ],
        'ErrorDeletingFile'     : [ 26, 'Erro', 'Failed to delete file' ],
        'ErrorEmptyFile'        : [ 27, 'Erro', 'Failed due to file empty' ],
        'ErrorOutputIncomplete' : [ 28, 'Erro', 'Failed due to incomplete output data' ],
        'ErrorInvalidRetKey'    : [ 29, 'Erro', 'Failed due to invalid return code key in program' ],
        'ErrorSagaDataGet'      : [ 30, 'Erro', 'Failed to retrieve data from SCAMPI' ],
        'ErrorOutfileExists'    : [ 31, 'Erro', 'Failed because output file already exists' ],
        'ErrorDirRename'        : [ 32, 'Erro', 'Failed to rename dir' ],
        'ErrorFileRename'       : [ 33, 'Erro', 'Failed to rename file' ],
        'ErrorReadingData'      : [ 34, 'Erro', 'Failed to read data' ],
        'ErrorFileCopy'         : [ 35, 'Erro', 'Failed to copy file' ],
        'ErrorFileNotLink'      : [ 36, 'Erro', 'Have a file not a symlink' ]
        }

  # ============================================================================ #
  # Methods                                                                      #
  # ============================================================================ #

    def _msgprefix(self):
        'Assemble prefix for all messages of format "YYYYMMDD-HHMMSS (progname): ".'
        return '{datetime} [{progname}]: '.format(datetime=time.strftime('%Y%m%d-%H%M%S', time.localtime()), progname=self._prog)

    def err_code(self, retkey):
        'Return the numerical return code value corresponding to the return code descriptor.'
        return self.err[retkey][0]

    def err_retkey(self, retcode):
        'Return the return code descriptor string corresponding to the numerical return code.'
        return [key for key in self.err.keys() if self.err[key][0] == retcode]

    def err_text(self, retkey):
        'Return a string reporting the return code and the standardised message for that return code.'
        if self.err.has_key(retkey):
            return '{type}: Program finished with Ont retcode={code} - {msg}\n'.format(
                type=self.err[retkey][1], code=self.err[retkey][0], msg=self.err[retkey][2])
        return ''

    def err_exit(self, retkey, exitprogram=False):
        'Function to report exit code and exit program.'
        if exitprogram:
            sys.exit(self.err[retkey][0])

    #def err_tidy(self, retkey, print2out=True, print2err=False, print2dump=False, exitprogram=False):
    #    'Function to call just before exiting program.'
    #    if print2out:
    #        sys.stdout.write('{0}{1}\n'.format(self._msgprefix(), self.err_text(retkey)))
    #        sys.stdout.flush()
    #    if print2err:
    #        sys.stderr.write('{0}{1}\n'.format(self._msgprefix(), self.err_text(retkey)))
    #        sys.stderr.flush()
    #    if print2dump:
    #        self.log_dump(self.err_text(retkey))
    #    self.log_tidy()
    #    if exitprogram:
    #        self.err_exit(retkey, exitprogram=exitprogram)

# ============================================================================ #
