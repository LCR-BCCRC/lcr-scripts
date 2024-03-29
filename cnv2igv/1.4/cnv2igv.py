#!/usr/bin/env python
import argparse
import re
import math
from abc import ABC, abstractmethod

# Globals ----------------------------------------------------------------

MODES = ['sclust', 'sequenza', 'titan', 'cnvkit','battenberg', 'controlfreec'] # add new seg filetypes here
LOH_TYPES = ['neutral', 'deletion', 'any']

# Classes ----------------------------------------------------------------

class Parser:
    ''' Extend this class and implement is_header(), parse_segment() methods.
        get_loh_flag() is optional.
    '''
    def __init__(self, stream, sample, loh_type):
        self.stream   = open(stream, 'r')
        self.sample   = sample
        self.loh_type = loh_type

    @abstractmethod
    def is_header(self, line):
        ''' Return true if line is part of header '''
        pass

    @abstractmethod
    def parse_segment(self, line):
        ''' Return Segment object '''
        pass

    def get_loh_flag(self):
        ''' Return string indicating if segment is LOH or not. '''
        pass

    # battenberg needs its own function because it outputs CN for each subclone separately
    def calculate_logratio_battenberg(self, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A):
        logr = None
        # replace NAs with zero, othervise convert string to floats for future calculations
        nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A = [0.0 if "NA" in value else float(value) for value in [nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A]]
        # calculate CN
        # Determine whether it's the major or minor allele represented by two states
        is_subclonal_maj = abs(nMaj1_A - nMaj2_A) > 0
        if is_subclonal_maj:
          segment_states_min = nMin1_A * 1
          segment_states_maj = nMaj1_A * frac1_A  + nMaj2_A * frac2_A
        # If it's not major, then it is the minor allele represented by two states
        else:
          segment_states_min = nMin1_A * frac1_A  + nMin2_A * frac2_A
          segment_states_maj = nMaj1_A*1
        cn = segment_states_min + segment_states_maj
        if cn == 0:
            logr = '-10'
        else:
            logr = math.log(cn, 2) - 1
        return(str(logr))

    def calculate_logratio(self, cn):
        cn   = float(cn)
        logr = None
        if cn == 0:
            logr = '-10'
        else:
            logr = math.log(cn, 2) - 1
        return(str(logr))

class CNVKitParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        chrm = line.split('\t', 1)[0]
        return True if chrm.startswith("chr") else False

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end = _line[0:3]
        cn = _line[5]
        logr = self.calculate_logratio(cn)

        return(Segment(chrm, start, end, cn, logr, self.sample))

class SequenzaParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        chrm = line.split('\t', 1)[0]

        if chrm.startswith('"'):
            chrm = chrm[1:]
        if chrm.endswith('"'):
            chrm = chrm[0:-1]

        if chrm == "chromosome":
            return True
        else:
            return False

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end = _line[0:3]
        if chrm.startswith('"'):
            chrm = chrm[1:]
        if chrm.endswith('"'):
            chrm = chrm[0:-1]
        cn, a, b = _line[9:12]
        if cn == 'NA':
            cn = 2 # assume it's diploid
        if a == 'NA':
            a = 1
        if b == 'NA':
            b =1
        loh_flag = self.get_loh_flag(cn, a, b)
        logr = self.calculate_logratio(cn)
        start = str(int(float(start)))
        end = str(int(float(end)))

        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))

    def get_loh_flag(self, cn, a, b):
        cn = int(cn)
        loh_flag = '0'
        if self.loh_type == 'neutral':
            if a == 'NA':
                loh_flag = 'NA'
            elif cn == 2 and int(a) == 2:
                loh_flag = '1'
        elif self.loh_type == 'deletion':
            if a == 'NA' or b == 'NA':
                loh_flag = 'NA'
            elif cn == 1 and (int(a) + int(b)) == 1:
                loh_flag = '1'
        elif self.loh_type == 'any':
            if b == 'NA':
                loh_flag = 'NA'
            elif cn <= 2 and int(b) == 0:
                loh_flag = '1'
        return(loh_flag)

class BattenbergParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        chrm = line.split('\t', 1)[0]
        return True if chrm == "chr" else False

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end, BAF, pval, orig_logr, cn, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A = _line[0:13]
        if not chrm.startswith("chr"):
            chrm = "chr"+ str(chrm)
        loh_flag = self.get_loh_flag(nMaj1_A, nMin1_A, nMin2_A, frac1_A, frac2_A)
        logr = self.calculate_logratio_battenberg(nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A)
        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))

    #actually use the LOH information from Battenberg, The column is nMin1_A (if < 1, LOH)
    def get_loh_flag(self, nMaj1_A, nMin1_A, nMin2_A, frac1_A, frac2_A):
        loh_flag = '0'
        if self.loh_type == 'neutral':
            if int(nMaj1_A) == 2 and int(nMin1_A) == 0:
                loh_flag = '1'
        elif self.loh_type == 'deletion':
           if (int(nMaj1_A) + int(nMin1_A)) == 1:
                loh_flag = '1'
        elif self.loh_type == 'any':
            # segments with no information of their copy number state:
            if "NA" in str(nMin1_A):
                loh_flag = '0'
            # events with no subclones and no minor allele are set with loh flag
            elif "NA" in str(frac2_A) and float(nMin1_A) == 0:
                loh_flag = '1'
            # events with no subclones, but where minor allele is present, are not assigned loh flag
            elif "NA" in str(frac2_A) and not float(nMin1_A) == 0:
                loh_flag = '0'
            # the rest of events have subclones
            else:
                # for events where major clone is prevalent
                if float(frac1_A) > float(frac2_A):
                    if float(nMin1_A) > 0:
                        loh_flag = '0'
                    # set loh flag if there is no minor allele
                    elif float(nMin1_A) == 0:
                        loh_flag = '1'
                # for events where subclone is prevalent
                else:
                    if float(nMin2_A) > 0:
                        loh_flag = '0'
                    # set loh flag if there is no minor allele
                    elif float(nMin2_A) == 0:
                        loh_flag = '1'
        return(loh_flag)

class SClustParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        return(False)

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end, x, cn = _line[1:] # not sure what the x column represents
        loh_flag = self.get_loh_flag()
        logr = self.calculate_logratio(cn)
        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))

    def get_loh_flag(self):
        return('NA')

class TitanParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        p = re.compile('start', re.I)
        m = re.search(p, line)
        return(m)

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end = _line[1:4]
        logr = _line[6]
        call = _line[8]
        cn = _line[9]
        loh_flag = self.get_loh_flag(call)

        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))

    def get_loh_flag(self, call):
        loh_flag = '0'
        if self.loh_type == 'neutral':
            if call == 'NLOH':
                loh_flag = '1'
        elif self.loh_type == 'deletion':
            if call == 'DLOH':
                loh_flag = '1'
        elif self.loh_type == 'any':
            if 'LOH' in call:
                loh_flag = '1'
        return(loh_flag)


class ControlfreecParser(Parser):
    def __init__(self, stream, sample, loh_type):
        super().__init__(stream, sample, loh_type)

    def is_header(self, line):
        chrm = line.split('\t', 1)[0]
        return True if chrm == "chr" else False

    def parse_segment(self, line):
        _line = line.split('\t')
        chrm, start, end, cn, status, genotype, uncert, somgerm, pgerml, wilk, pval = _line[0:11]
        if not chrm.startswith("chr"):
            chrm = "chr"+ str(chrm)
        # if both alleles are not present in genotype, it indicates LOH event
        alleles = ["A", "B"]
        if all(x in genotype for x in alleles):
            loh_flag = str(0)
        else:
            loh_flag = str(1)
        # calsulate logratio for somatic events that pass significance threshold
        if somgerm == "somatic" and not "NA" in str(pval) and float(pval) <= 0.1:
            logr = self.calculate_logratio(cn)
        else:
            logr = str(0.0)
        return(Segment(chrm, start, end, cn, logr, self.sample, loh_flag))


class Segment:
    def __init__(self, chrm, start, end, cn, logr, sample, loh_flag = 'NA'):
        self.chrm     = chrm
        self.start    = start
        self.end      = end
        self.loh      = loh_flag
        self.cn       = cn
        self.cn_state = self.get_cnv_state()
        self.logr     = logr
        self.sample   = sample

    def get_cnv_state(self):
        cn_state = 'NEUT'
        cn = float(self.cn)

        if cn > 3:
            cn_state = 'AMP'
        elif cn == 3:
            cn_state = 'GAIN'
        elif cn == 1:
            cn_state = 'HETD'
        elif cn == 0:
            cn_state = 'HOMD'
        return(cn_state)


    def to_igv(self, prepend):
        if prepend:
            self.chrm = "chr"+self.chrm
        return('\t'.join([self.sample, self.chrm, self.start, \
                          self.end, self.loh, self.logr]))

    def to_oncocircos(self, prepend):
        if prepend:
            self.chrm = "chr"+self.chrm
        return('\t'.join([self.sample, self.chrm, self.start, \
                          self.end, self.loh, self.logr, self.cn_state]))

# Argument Parser --------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--include_loh', '--loh_type', help = 'Type of LOH events to consider for LOH flag column. Default: any.',
                        choices = LOH_TYPES, default = 'any', const = 'any', nargs = '?', dest = 'loh_type')
    parser.add_argument('--prepend_chr', help = 'Prepend "chr" to chromosome column.', action = 'store_true')
    parser.add_argument('--oncocircos', help = 'Output OncoCircos compatible output.', action = 'store_true')
    parser.add_argument('--gistic', help = 'Output GISTIC2.0 compatible output.', action = 'store_true')
    parser.add_argument('--mode', help = 'Input file type.', choices = MODES, required = True)
    parser.add_argument('--sample', '--sequenza_sample', help = 'Sample name for IGV ID column.', dest = 'sample', required = True)
    parser.add_argument('seg_file', help = 'Segmentation file to convert to IGV format.')
    return(parser.parse_args())

# Main -------------------------------------------------------------------

def main():
    args       = parse_arguments()
    mode       = args.mode
    sample     = args.sample
    seg_file   = args.seg_file
    prepend    = args.prepend_chr
    loh_type   = args.loh_type
    oncocircos = args.oncocircos
    gistic     = args.gistic
    parser     = None

    if gistic:
        prepend = True # backwards compatibility

    # - implement your own parser for any new segment file types by extending Parser
    # - add filetype to MODES list above
    # - add condition for new filetype here and instantiate your own Parser as parser
    if mode == 'sequenza':
        parser = SequenzaParser(seg_file, sample, loh_type)
    elif mode == 'sclust':
        parser = SClustParser(seg_file, sample, loh_type)
    elif mode == 'titan':
        parser = TitanParser(seg_file, sample, loh_type)
    elif mode == 'cnvkit':
        parser = CNVKitParser(seg_file, sample, loh_type)
    elif mode == "battenberg":
        parser = BattenbergParser(seg_file, sample, loh_type)
    elif mode == 'controlfreec':
        parser = ControlfreecParser(seg_file, sample, loh_type)

    # header
    if not oncocircos:
        print('ID\tchrom\tstart\tend\tLOH_flag\tlog.ratio')

    for line in parser.stream:
        if parser.is_header(line):
            continue
        seg = parser.parse_segment(line)

        if oncocircos:
            print(seg.to_oncocircos(prepend))
        else:
            print(seg.to_igv(prepend))

if __name__ == '__main__':
    main()
