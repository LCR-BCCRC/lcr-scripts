#!/usr/bin/env python
import argparse
import re
import os
import math
from abc import ABC, abstractmethod

# Globals ----------------------------------------------------------------

MODES = ['sclust', 'sequenza', 'titan', 'cnvkit', 'purecn', 'purecn_cnvkit', 'battenberg', 'controlfreec'] # add new seg filetypes here
LOH_TYPES = ['neutral', 'deletion', 'any']

# Classes ----------------------------------------------------------------

class Parser:
    ''' Extend this class and implement is_header(), parse_segment() methods.
        get_loh_flag() is optional.
    '''
    def __init__(self, stream, sample, mode, loh_type, logr_type = "corrected"):
        self.stream   = open(stream, 'r')
        self.filename = stream
        self.sample   = sample
        self.mode = mode
        self.loh_type = loh_type
        self.logr_type = logr_type

        # If the user provided --sample (not default), prefer it verbatim
        self._prefer_arg_sample = (sample is not None and sample != 'SAMPLE')

    @abstractmethod
    def is_header(self, line):
        ''' Return true if line is part of header '''
        pass

    @abstractmethod
    def parse_segment(self, line, logr_type, mode):
        ''' Return Segment object '''
        pass

    def resolve_sample(self, line_sample=None):
        """Return the ID to use for this row without mutating self.sample."""
        if self._prefer_arg_sample:
            return self.sample
        # fallbacks if user did not provide --sample
        if line_sample and line_sample.strip():
            return line_sample.strip()
        return self.sample  # last resort

    def get_loh_flag(self):
        ''' Return string indicating if segment is LOH or not. '''
        pass

    # battenberg needs its own function because it outputs CN for each subclone separately
    def calculate_cn_battenberg(self, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A):
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
        return(str(cn))

    def calculate_logratio(self, cn):
        cn   = float(cn)
        logr = None
        if cn == 0:
            logr = '-10'
        else:
            logr = math.log(cn, 2) - 1
        return(str(logr))

class PurecnParser(Parser):
    def __init__(self, stream, sample, mode, loh_type, logr_type):
        super().__init__(stream, sample, mode, loh_type, logr_type)

    def is_header(self, line):
        toks = line.rstrip('\n').split('\t')
        # expect header like: ID    chromosome  start  end ...
        return len(toks) > 1 and (toks[0].lower() in ('id', 'sample')
                                  or toks[1].lower().startswith('chrom'))

    def parse_segment(self, line, logr_type, mode):
        t = line.rstrip('\n').split('\t')
        line_sample, chrm, start, end, num_markers = t[0:5]
        cn = t[6]
        logr = self.calculate_logratio(cn) if logr_type == "corrected" else t[5]
        sample_id = self.resolve_sample(line_sample)
        return Segment(chrm, start, end, cn, num_markers, logr, sample_id, mode)

class CNVKitParser(Parser):
    def __init__(self, stream, sample, mode, loh_type, logr_type):
        super().__init__(stream, sample, mode, loh_type, logr_type)

    def is_header(self, line):
        # expect header like: chromosome, start, end, gene, log2, baf, cn, ...
        chrm = line.split('\t', 1)[0]
        if chrm.startswith('"'):
            chrm = chrm[1:]
        if chrm.endswith('"'):
            chrm = chrm[0:-1]
        if chrm == "chromosome":
            return True
        else:
            return False

    def parse_segment(self, line, logr_type, mode):
        _line = line.split('\t')
        chrm, start, end = _line[0:3]
        if chrm.startswith('"'):
            chrm = chrm[1:]
        if chrm.endswith('"'):
            chrm = chrm[0:-1]
        cn = _line[6]
        cn1 = None if _line[7] == '' else _line[7]
        cn2 =  None if _line[8] == '' else _line[8]
        loh_flag = self.get_loh_flag(cn, cn1, cn2)
        logr = self.calculate_logratio(cn) if logr_type == "corrected" else _line[4]
        num_markers = _line[10]
        return(Segment(chrm, start, end, cn, num_markers, logr, self.sample, mode, loh_flag))

    def get_loh_flag(self, cn, cn1, cn2):
        cn = int(cn)
        loh_flag = '0'
        if self.loh_type == 'neutral':
            if cn1 == None:
                loh_flag = 'NA'
            elif cn == 2 and int(cn2) == 0:
                loh_flag = '1'
        elif self.loh_type == 'deletion':
            if cn1 == None:
                loh_flag = 'NA'
            elif cn == 1 and (int(cn1) + int(cn2)) == 1:
                loh_flag = '1'
        elif self.loh_type == 'any':
            if cn1 == None:
                loh_flag = 'NA'
            elif int(cn2) == 0 and int(cn1) + int(cn2) > 0:
                loh_flag = '1'
        return(loh_flag)

class SequenzaParser(Parser):
    def __init__(self, stream, sample, mode, loh_type, logr_type):
        super().__init__(stream, sample, mode, loh_type, logr_type)

    def is_header(self, line):
        # expect header like: chromosome, start.pos, end.pos, ...
        chrm = line.split('\t', 1)[0]
        if chrm.startswith('"'):
            chrm = chrm[1:]
        if chrm.endswith('"'):
            chrm = chrm[0:-1]
        if chrm == "chromosome":
            return True
        else:
            return False

    def parse_segment(self, line, logr_type, mode):
        _line = line.split('\t')
        chrm, start, end = _line[0:3]
        if chrm.startswith('"'):
            chrm = chrm[1:]
        if chrm.endswith('"'):
            chrm = chrm[0:-1]
        cn, a, b = _line[9:12]
        depth_ratio = float(_line[6])
        num_markers = _line[4]
        # Follow what is done for battenberg: replace NAs with zero
        # will set logr = -10 in "corrected mode"
        if cn == 'NA':
            cn = 0
        if a == 'NA':
            a = 0
        if b == 'NA':
            b = 0
        loh_flag = self.get_loh_flag(cn, a, b)
        if logr_type == "corrected":
            logr = self.calculate_logratio(cn)
        else:
            logr = str(math.log(depth_ratio, 2))
        return(Segment(chrm, start, end, cn, num_markers, logr, self.sample, mode, loh_flag))

    def get_loh_flag(self, cn, a, b):
        cn = int(cn)
        loh_flag = '0'
        if self.loh_type == 'neutral':
            if cn == 0 :
                loh_flag = 'NA'
            elif cn == 2 and int(a) == 2:
                loh_flag = '1'
        elif self.loh_type == 'deletion':
            if cn == 0:
                loh_flag = 'NA'
            elif cn == 1 and (int(a) + int(b)) == 1:
                loh_flag = '1'
        elif self.loh_type == 'any':
            if cn == 0:
                loh_flag = 'NA'
            elif int(b) == 0 and int(a) + int(b) > 0:
                loh_flag = '1'
        return(loh_flag)

class BattenbergParser(Parser):
    def __init__(self, stream, sample, mode, loh_type, logr_type):
        super().__init__(stream, sample, mode, loh_type, logr_type)

    def is_header(self, line):
        chrm = line.split('\t', 1)[0]
        return True if chrm == "chr" else False

    def parse_segment(self, line, logr_type, mode):
        _line = line.split('\t')
        chrm, start, end, BAF, pval, num_markers, orig_logr, cn_orig, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A = _line[0:14]
        loh_flag = self.get_loh_flag(nMaj1_A, nMin1_A, nMin2_A, frac1_A, frac2_A)
        cn = self.calculate_cn_battenberg(nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A)
        logr = self.calculate_logratio(cn) if logr_type == "corrected" else orig_logr
        return(Segment(chrm, start, end, cn, num_markers, logr, self.sample, mode, loh_flag))

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
                loh_flag = 'NA'
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
    def __init__(self, stream, sample, mode, loh_type, logr_type):
        super().__init__(stream, sample, mode, loh_type, logr_type)

    def is_header(self, line):
        chrm = line.split('\t', 1)[0]
        return True if chrm == "chr" else False

    def parse_segment(self, line, mode, logr_type):
        _line = line.rstrip('\n').split('\t')
        chrm, start, end, cn, status, genotype, uncert, somgerm, pgerml, wilk, pval, num_markers, log2 = _line[0:13]
        # if both alleles are not present in genotype, it indicates LOH event
        alleles = ["A", "B"]
        if all(x in genotype for x in alleles):
            loh_flag = str(0)
        else:
            loh_flag = str(1)
        # calculate logratio for somatic events that pass significance threshold
        if logr_type == "corrected" and somgerm == "somatic" and not "NA" in str(pval) and float(pval) <= 0.1:
            logr = self.calculate_logratio(cn)
        elif logr_type == "raw" and somgerm == "somatic" and not "NA" in str(pval) and float(pval) <= 0.1:
            logr = str(log2)
        else:
            logr = str(0.0)
        return(Segment(chrm, start, end, cn, num_markers, logr, self.sample, mode, loh_flag))


class Segment:
    def __init__(self, chrm, start, end, cn, num_markers, logr, sample, mode, loh_flag = 'NA'):
        self.chrm     = str(chrm)
        self.start    = str(int(float(start)))
        self.end      = str(int(float(end)))
        self.loh      = str(loh_flag)
        self.cn       = str(cn)
        self.cn_state = self.get_cnv_state()
        self.num_markers = str(num_markers)
        self.logr     = str(logr)
        self.sample   = str(sample)
        self.mode   = str(mode)

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

    @staticmethod
    def _maybe_chr(prepend, chrm):
        if not prepend:
            return chrm
        return chrm if str(chrm).startswith("chr") else "chr" + str(chrm)

    def to_igv(self, prepend):
        ch = self._maybe_chr(prepend, self.chrm)
        return '\t'.join([self.sample, ch, str(self.start), str(self.end),
                          self.mode, self.loh, self.cn, self.num_markers, self.logr])

    def to_oncocircos(self, prepend):
        ch = self._maybe_chr(prepend, self.chrm)
        return '\t'.join([self.sample, ch, str(self.start), str(self.end),
                          self.loh, self.logr, self.cn_state])

# Argument Parser --------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--include_loh', '--loh_type', help = 'Type of LOH events to consider for LOH flag column. Default: any.',
                        choices = LOH_TYPES, default = 'any', const = 'any', nargs = '?', dest = 'loh_type')
    parser.add_argument('--prepend_chr', help = 'Prepend "chr" to chromosome column.', action = 'store_true')
    parser.add_argument('--oncocircos', help = 'Output OncoCircos compatible output.', action = 'store_true')
    parser.add_argument('--gistic', help = 'Output GISTIC2.0 compatible output.', action = 'store_true')
    parser.add_argument('--preserve_log_ratio', help = "Preserve the log ratio and report absolute CN as a separate column", action = 'store_true')
    parser.add_argument('--mode', help = 'Input file type.', choices = MODES, required = True)
    parser.add_argument('--sample', '--sequenza_sample', help = 'Sample name for IGV ID column.', dest = 'sample', default = 'SAMPLE', required = False)
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
    preserve   = args.preserve_log_ratio
    parser     = None
    if preserve:
        logr_type = "raw"
    else:
        logr_type = "corrected"
    if gistic:
        prepend = True # backwards compatibility

    # - implement your own parser for any new segment file types by extending Parser
    # - add filetype to MODES list above
    # - add condition for new filetype here and instantiate your own Parser as parser
    if mode == 'sequenza':
        parser = SequenzaParser(seg_file, sample, mode, loh_type, logr_type)
    elif mode == 'purecn' or mode == 'purecn_cnvkit':
        parser = PurecnParser(seg_file, sample, mode, loh_type, logr_type)
    elif mode == 'sclust':
        parser = SClustParser(seg_file, sample, loh_type)
    elif mode == 'titan':
        parser = TitanParser(seg_file, sample, loh_type)
    elif mode == 'cnvkit':
        parser = CNVKitParser(seg_file, sample, mode, loh_type, logr_type)
    elif mode == "battenberg":
        parser = BattenbergParser(seg_file, sample, mode, loh_type, logr_type)
    elif mode == 'controlfreec':
        parser = ControlfreecParser(seg_file, sample, mode, loh_type, logr_type)

    header = 'ID\tchrom\tstart\tend\tmodule\tLOH_flag\tCN\tnum_markers\tlog.ratio'
    if not oncocircos:
        print(header)

    for line in parser.stream:
        if parser.is_header(line):
            continue
        seg = parser.parse_segment(line, logr_type = logr_type, mode = mode)

        if oncocircos:
            print(seg.to_oncocircos(prepend))
        else:
            print(seg.to_igv(prepend))

if __name__ == '__main__':
    main()
