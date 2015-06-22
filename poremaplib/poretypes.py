#!/usr/bin/env python

# ============================================================================ #
# poretypes.py                                                                 #
# ============================================================================ #
'''
Data types and methods shared across the poremap programs
'''
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# April 2015                                                                   #
# ============================================================================ #

import numpy as np

class Poretypes(object):
    'A class for poremap program data types and methods.'

    def __init__(self):

        self.alignstatsdef = [ # Primary key : 11 fields
            [ 'runid', 'S50',  'Run identifier' ], # 0
            [ 'readid', 'S100',  'Read identifier' ],
            [ 'readtype', 'S10',  'Read type : 2d, temp, comp, mixed, unknown' ],
            [ 'readclass', 'S10',  'Read class : all, pass, fail' ],
            [ 'datatype', 'S10',  'Data type : minion, random, mixed, unknown' ],
            [ 'mapprog', 'S20',  'Mapping program' ], # 5
            [ 'mapparams', 'S500', 'Mapping program parameters' ],
            [ 'samflag', np.int,  'SAM bitwise flag (decimal)' ],
            [ 'refcontigid', 'S20',  'Refcontigid of match ' ],
            [ 'refcontigpos1', np.int, '1-based start of alignment in refcontigid' ],
            [ 'mapq', np.int, 'Mapping quality score of alignment' ], # 10
        
            [ 'readbp', np.int, 'Total read length' ],
            [ 'isprimaryaln', np.int, '1 if it is the primary alignment, otherwise 0' ],
            [ 'isnonrandomaln', np.int, '1 if it has characteristics of a non-random alignment, ' \
                                        '0 if it looks random, -1 if it has not been implemented yet.' ],
            [ 'alnlen', np.int, 'Length of the alignment between the read and refcontig' ],
            [ 'alnrefbp', np.int, 'Number of ref bases involved in the alignment' ], # 15
            [ 'alnrefstartpos1', np.int, '1-based start position of alignment in refcontigid (inclusive)' ],
            [ 'alnrefendpos1', np.int, '1-based end position of alignment in refcontigid (inclusive)' ],
            [ 'alnreadbp', np.int, 'Number of read bases involved in the alignment' ],
            [ 'alnreadstartpos1', np.int, '1-based start position of alignment in read (inclusive)' ],
            [ 'alnreadendpos1', np.int, '1-based end position of alignment in read (inclusive)' ], # 20
            [ 'numC', np.int, 'Number of same-as-ref segments in alignment' ],
            [ 'numM', np.int, 'Number of matching segments' ],
            [ 'numI', np.int, 'Number of insertions (correct irrespective of leading/trailing Hard/Soft clipped segments)' ],
            [ 'numD', np.int, 'Number of deletions (correct irrespective of leading/trailing Hard/Soft clipped segments)' ],
            [ 'freqCperbp', np.float, 'C freq : Number of same-as-ref segments per aligned base' ], # 25
            [ 'freqMperbp', np.float, 'M freq : Number of match segments per aligned base' ],
            [ 'freqIperbp', np.float, 'I freq : Number of insertion segments per aligned base' ],
            [ 'freqDperbp', np.float, 'D freq : Number of deletion segments per aligned base' ],
            [ 'lenS', np.int, 'Total number of diff-to-ref bases in aligned read' ],
            [ 'lenC', np.int, 'Total number of same-as-ref bases in aligned read' ], # 30
            [ 'lenM', np.int, 'Total number of matching bases in aligned read' ],
            [ 'lenI', np.int, 'Total number of insertion bases in aligned read' ],
            [ 'lenD', np.int, 'Total number of deleted bases in aligned read' ],
            [ 'pctS', np.float, 'Pct bases with diff-to-ref call in read' ],
            [ 'pctC', np.float, 'Pct bases with same-as-ref call in read' ], # 35
            [ 'pctM', np.float, 'Pct bases in matching segment of read' ],
            [ 'pctI', np.float, 'Pct bases in insertions in read' ],
            [ 'pctD', np.float, 'Pct bases in deletions in read' ],
            [ 'meanrunlenC', np.float, 'Mean length of same-as-ref segments in alignment' ],
            [ 'meanrunlenM', np.float, 'Mean length of match segments in alignment' ], # 40
            [ 'meanrunlenI', np.float, 'Mean length of insertion segments in alignment' ],
            [ 'meanrunlenD', np.float, 'Mean length of deletion segments in alignment' ],
            [ 'meanbq', np.float, 'Mean base quality of entire read' ],
            [ 'meanbqC', np.float, 'Mean base quality of same-as-ref calls' ],
            [ 'meanbqM', np.float, 'Mean base quality of matching bases' ], # 45
            [ 'meanbqI', np.float, 'Mean base quality of insertion bases' ]
        ]
        self.alignstatsheader = [x[0] for x in self.alignstatsdef]
        self.alignstatsdtype = [(x[0], x[1]) for x in self.alignstatsdef]
        self.alignstatsdesc = [(x[0], x[2]) for x in self.alignstatsdef]

        self.alignclassdef = [ # Primary key : 11 fields
            [ 'runid', 'S50',  'Run identifier' ],
            [ 'readid', 'S100',  'Read identifier' ],
            [ 'readtype', 'S10',  'Read type : 2d, temp, comp, mixed, unknown' ],
            [ 'readclass', 'S10',  'Read class : all, pass, fail' ],
            [ 'datatype', 'S10',  'Data type : minion, random, mixed, unknown' ],
            [ 'mapprog', 'S20',  'Mapping program' ],
            [ 'mapparams', 'S500', 'Mapping program parameters' ],
            [ 'samflag', np.int,  'SAM bitwise flag (decimal)' ],
            [ 'refcontigid', 'S20',  'Refcontigid of match ' ],
            [ 'refcontigpos1', np.int, '1-based start of alignment in refcontigid' ],
            [ 'mapq', np.int, 'Mapping quality score of alignment' ],
            [ 'isnonrandomaln', np.int, '1 if it has characteristics of a non-random alignment, ' \
                                        '0 if it looks random, -1 if it has not been implemented yet.' ]
        ]
        self.alignclassheader = [x[0] for x in self.alignclassdef]
        self.alignclassdtype = [(x[0], x[1]) for x in self.alignclassdef]
        self.alignclassdesc = [(x[0], x[2]) for x in self.alignclassdef]

        self.readstatsdef = [ # Primary key : 7 fields
            [ 'runid', 'S50', 'Run identifier' ], # 0
            [ 'readid', 'S100', 'Read identifier' ],
            [ 'readtype', 'S10', 'Read type : 2d, temp, comp, mixed, unknown' ],
            [ 'readclass', 'S10', 'Read class : all, pass, fail' ],
            [ 'datatype', 'S10', 'Data type : minion, random, mixed, unknown' ],
            [ 'mapprog', 'S20', 'Mapping program' ], # 5
            [ 'mapparams', 'S500', 'Mapping program parameters' ],
        
            [ 'readbp', np.int, 'Total read length' ],
            [ 'ismapped', np.int, '1 if mapped to any reference supplied, otherwise 0' ],
            [ 'numalignments', np.int, 'The number of non-random alignments for this read' ],
            [ 'numrefcontigs', np.int, 'Number of refcontigs with genuine alignments' ], # 10
            [ 'allalnlen', np.int, 'Sum of length of alignment(s) between read and refcontig(s)' ],
            [ 'allalnreadbp', np.int, 'Total number of read bases involved in alignment(s)' ],
            [ 'allalnrefbp', np.int, 'Total number of ref bases involved in alignment(s)' ],
            [ 'allalnfreqCperbp', np.float, 'C freq : Number of same-as-ref segments per aligned base' ],
            [ 'allalnfreqMperbp', np.float, 'M freq : Number of match segments per aligned base' ], # 15
            [ 'allalnfreqIperbp', np.float,'I freq : Number of insertion segments per aligned base' ],
            [ 'allalnfreqDperbp', np.float, 'D freq : Number of deletion segments per aligned base' ],
            [ 'allalnpctS', np.float, 'Pct bases with diff-to-ref call in alignment(s)' ],
            [ 'allalnpctC', np.float, 'Pct bases with same-as-ref call in alignment(s)' ],
            [ 'allalnpctM', np.float, 'Pct bases in matches in alignment(s)' ], # 20
            [ 'allalnpctI', np.float, 'Pct bases in insertions in alignment(s)' ],
            [ 'allalnpctD', np.float, 'Pct bases in deletions in alignment(s)' ],
            [ 'allalnmeanrunlenC', np.float, 'Mean length of same-as-ref matches in alignment(s)' ],
            [ 'allalnmeanrunlenM', np.float, 'Mean length of matches in alignment(s)' ],
            [ 'allalnmeanrunlenI', np.float, 'Mean length of insertions in alignment(s)' ], # 25
            [ 'allalnmeanrunlenD', np.float, 'Mean length of deletions in alignment(s)' ],
            [ 'allalnmeanbq', np.float, 'Mean base quality for all bases in alignment(s)' ],
            [ 'allalnmeanbqC', np.float, 'Mean base quality for all same-as-ref calls in alignment(s)' ],
            [ 'allalnmeanbqM', np.float, 'Mean base quality for all matching bases in alignment(s)' ],
            [ 'allalnmeanbqI', np.float, 'Mean base quality for insertions in alignment(s)' ], # 30
            [ 'alltargetaligncnt', np.int, 'Number of reads that align to a target refcontig' ],
            [ 'alltargetalignbp', np.int, 'Total bases in alignments to a target refcontig' ],
            [ 'allcontrolaligncnt', np.int, 'Total bases in reads that align to a control refcontig' ],
            [ 'allcontrolalignbp', np.int, 'Total bases in alignments to a control refcontig' ],
            [ 'allrandalncnt', np.int, 'Number of alignments rejected for being more similar to random alignments' ], # 35
            [ 'allrandalncntpct', np.float, 'Percent of alignments rejected for being more similar to random alignments' ],
            [ 'allrandalnbp', np.int, 'Bases in random alignments' ],
            [ 'allrandalnbppct', np.float, 'Percent bases in random alignments' ],
            [ 'allrandalnmeanbp', np.float, 'Mean length of random alignments' ],
            [ 'prialnlen', np.int, 'Length of primary alignment between read and refcontig' ], # 40
            [ 'prialnreadbp', np.int, 'Number of read bases involved in primary alignment' ],
            [ 'prialnrefbp', np.int, 'Number of ref bases involved in primary alignment' ],
            [ 'prialnfreqCperbp', np.float, 'C freq : Number of same-as-ref segments per aligned base' ],
            [ 'prialnfreqMperbp', np.float, 'M freq : Number of match segments per aligned base' ],
            [ 'prialnfreqIperbp', np.float, 'I freq : Number of insertion segments per aligned base' ], # 45
            [ 'prialnfreqDperbp', np.float, 'D freq : Number of deletion segments per aligned base' ],
            [ 'prialnpctS', np.float, 'Pct bases with diff-to-ref call in primary alignment' ],
            [ 'prialnpctC', np.float, 'Pct bases with same-as-ref call in primary alignment' ],
            [ 'prialnpctM', np.float, 'Pct bases in matches in primary alignment' ],
            [ 'prialnpctI', np.float, 'Pct bases in insertions in primary alignment' ], # 50
            [ 'prialnpctD', np.float, 'Pct bases in deletions in primary alignment' ],
            [ 'prialnmeanrunlenC', np.float, 'Mean length of same-as-ref segments in primary alignment' ],
            [ 'prialnmeanrunlenM', np.float, 'Mean length of matches in primary alignment' ],
            [ 'prialnmeanrunlenI', np.float, 'Mean length of insertions in primary alignment' ],
            [ 'prialnmeanrunlenD', np.float, 'Mean length of deletions in primary alignment' ], # 55
            [ 'prialnmeanbq', np.float, 'Mean base quality for all bases in primary alignment' ],
            [ 'prialnmeanbqC', np.float, 'Mean base quality for all same-as-ref calls in primary alignment' ],
            [ 'prialnmeanbqM', np.float, 'Mean base quality for all matching bases in primary alignment' ],
            [ 'prialnmeanbqI', np.float, 'Mean base quality for insertions in primary alignment' ],
            [ 'pritargetaligncnt', np.int, 'Number of reads that align to a target refcontig' ], # 60
            [ 'pritargetalignbp', np.int, 'Total bases in alignments to a target refcontig' ],
            [ 'pricontrolaligncnt', np.int, 'Total bases in reads that align to a control refcontig' ],
            [ 'pricontrolalignbp', np.int, 'Total bases in alignments to a control refcontig' ],
            [ 'prirandalncnt', np.int, 'Number of alignments rejected for being more similar to random alignments' ],
            [ 'prirandalncntpct', np.float, 'Percent of alignments rejected for being more similar to random alignments' ], # 65
            [ 'prirandalnbp', np.int, 'Bases in random alignments' ],
            [ 'prirandalnbppct', np.float, 'Percent bases in random alignments' ],
            [ 'prirandalnmeanbp', np.float, 'Mean length of random alignments' ],
            [ 'prialncnt', np.int, 'Number of non-random primary alignments' ]
        ]
        self.readstatsheader = [x[0] for x in self.readstatsdef]
        self.readstatsdtype = [(x[0], x[1]) for x in self.readstatsdef]
        self.readstatsdesc = [(x[0], x[2]) for x in self.readstatsdef]

        self.runstatsdef = [ # Primary key : 7 fields
            [ 'runid', 'S50', 'Run identifier' ], # 0
            [ 'readtype', 'S10', 'Read type : 2d, temp, comp, mixed, unknown' ],
            [ 'readclass', 'S10', 'Read class : all, pass, fail' ],
            [ 'datatype', 'S10', 'Data type : minion, rand, mixed or unknown' ],
            [ 'mapprog', 'S20', 'Mapping program' ],
            [ 'mapparams', 'S500', 'Mapping program parameters' ], # 5
            [ 'runstattype', 'S10', 'Run statistic type: allaln, prialn, oneref, mulref' ],
        
            [ 'runreadcnt', np.int, 'Total number of reads' ],
            [ 'runreadbp', np.int, 'Total bases in all reads' ],
            [ 'runmeanreadbp', np.int, 'Mean length of all reads' ],
            [ 'umapreadcnt', np.int, 'Number of reads that do not map to any refcontig' ], # 10
            [ 'umapreadcntpct', np.float, '% reads that do not map to any refcontig' ],
            [ 'umapreadbp', np.int, 'Total bases in unmapped reads' ],
            [ 'umapreadbppct', np.float, '% bases in unmapped reads' ],
            [ 'umapmeanreadbp', np.int, 'Mean length of unmapped reads' ],
            [ 'mapreadcnt', np.int, 'Number of reads' ], # 15
            [ 'mapreadcntpct', np.float, 'Number of reads' ],
            [ 'mapreadbp', np.int, 'Total bases' ],
            [ 'mapreadbppct', np.float, 'Total bases' ],
            [ 'mapmeanreadbp', np.int, 'Mean length of mapped reads' ],

            [ 'alncnt', np.int, 'Number of alignments' ], # 20
            [ 'alnbp', np.int, 'Aligned bases' ],
            [ 'alnbppct', np.float, 'Yield as percent of total bases' ],
            [ 'alnmeanbp', np.int, 'Mean length of aligned segments' ],
            [ 'Spct', np.float, 'Mean pct bases with diff-to-ref call' ],
            [ 'Cpct', np.float, 'Mean pct bases with same-as-ref call' ], # 25
            [ 'Mpct', np.float, 'Mean pct bases in aligned segments' ],
            [ 'Ipct', np.float, 'Mean pct bases in insertions' ],
            [ 'Dpct', np.float, 'Mean pct bases in deletions' ],
            [ 'Cmeandistbp', np.float, 'Mean number of bases separating segments of same-as-ref segments' ],
            [ 'Mmeandistbp', np.float, 'Mean number of bases separating matching segments' ], # 30
            [ 'Imeandistbp', np.float, 'Mean number of bases separating insertion segments' ],
            [ 'Dmeandistbp', np.float, 'Mean number of bases separating deletion segments' ],
            [ 'Cmeanlen', np.float, 'Mean length of same-as-ref matches' ],
            [ 'Mmeanlen', np.float, 'Mean length of matches' ],
            [ 'Imeanlen', np.float, 'Mean length of insertions' ], # 35
            [ 'Dmeanlen', np.float, 'Mean length of deletions' ],
            [ 'targetaligncnt', np.int, 'Number of reads that align to a target refcontig' ],
            [ 'targetalignbp', np.int, 'Total bases in alignments to a target refcontig' ],
            [ 'controlaligncnt', np.int, 'Total bases in reads that align to a control refcontig' ],
            [ 'controlalignbp', np.int, 'Total bases in alignments to a control refcontig' ], # 40
            [ 'randalncnt', np.int, 'Number of alignments rejected for being more similar to random alignments' ],
            [ 'randalncntpct', np.float, 'Percent of alignments rejected for being more similar to random alignments' ],
            [ 'randalnbp', np.int, 'Bases in random alignments' ],
            [ 'randalnbppct', np.float, 'Percent bases in random alignments' ],
            [ 'randalnmeanbp', np.float, 'Mean length of random alignments' ] # 45
        ]
        self.runstatsheader = [x[0] for x in self.runstatsdef]
        self.runstatsdtype = [(x[0], x[1]) for x in self.runstatsdef]
        self.runstatsdesc = [(x[0], x[2]) for x in self.runstatsdef]

    # ============================================================================ #
    # Methods                                                                      #
    # ============================================================================ #

# ============================================================================ #
# Class tests                                                                  #
# ============================================================================ #

if __name__ == "__main__":

    pass

# ============================================================================ #
