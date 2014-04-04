#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import re

import string
import pylab 

class Table:

    def __init__(self, fn, min_mapq):
        ncol = -1

        # read in the data
        data = list()
        fh = open(fn, 'r')
        for line in fh:
            line = line.rstrip("\r\n")
            if '#' == line[0]:
                continue
            else:
                table = self.__parse(line)
                if table[0] < min_mapq:
                    continue
                data.append(table)
                if -1 == ncol:
                    ncol = len(table)
                elif ncol != len(table):
                    sys.stderr.write("Error: table was misformatted\n")
                    sys.exit(1)
                    
        fh.close()

        # re-arrange the data
        self.matrix = [[data[i][j] for j in xrange(ncol)] for i in xrange(len(data))]

    # NB: assumes sorted descending
    def __parse(self, line):
        tokens = re.findall(r'[\w+-.]+', line)
        for i in range(0, 13):
            tokens[i] = int(tokens[i])
        for i in range(13, len(tokens)):
            tokens[i] = float(tokens[i])
        return tokens

def main(options):
    fmts = ['-', '--', '-.', ':', '.']

    if None == options.fns or 0 == len(options.fns):
        return

    # read in through the files
    fmts_i = 0
    for fn in options.fns:
        table = Table(fn, options.min_mapq)
        # extract sensitivity 
        x = [table.matrix[i][17] for i in range(len(table.matrix))]
        # extract specificity
        y = [table.matrix[i][16] for i in range(len(table.matrix))]
        pylab.plot(x, y, fmts[fmts_i])
        fmts_i = fmts_i + 1
        if len(fmts) <= fmts_i:
            fmts_i = 0
    # plot values
    pylab.title('DWGSIM ROC')
    pylab.xlabel('specificity')
    pylab.ylabel('sensitivity')
    # xlim
    if None != options.xlim:
        nums = re.findall(r'\d+\.?\d*', options.xlim)
        if 2 == len(nums):
            pylab.xlim(xmin=float(nums[0]), xmax=float(nums[1]))
    # ylim
    if None != options.ylim:
        nums = re.findall(r'\d+\.?\d*', options.ylim)
        if 2 == len(nums):
            pylab.ylim(ymin=float(nums[0]), ymax=float(nums[1]))
    # legend
    if not options.no_legend:
        if None != options.ids and len(options.fns) == len(options.ids):
            pylab.legend(options.ids, title='IDs', loc='lower left')
        elif options.infer_ids:
            ids = list()
            for fn in options.fns:
                i = re.sub(r'^.*\/', '', fn)
                i = re.sub(r'\.sam\.txt$', '', i)
                i = re.sub(r'^out', '', i)
                i = re.sub(r'^dwgsim_eval', '', i)
                i = re.sub(r'^\.', '', i)
                ids.append(i)
            pylab.legend(ids, title='IDs', loc='lower left')
    # show the plot
    if None == options.out:
        pylab.show()
    else:
        # check if the file extension is present, if not, add it
        if len(options.out) <= len(options.outtype) or options.outtype != options.out[len(options.out)-len(options.outtype):]:
            pylab.savefig("%s.%s" % (options.out, options.outtype), format=options.outtype)
        else:
            pylab.savefig(options.out, format=options.outtype)
        pylab.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--fn', help="an output file from dwgsim_eval", action="append", dest="fns")
    parser.add_option('--id', help="the corresponding ID to add to the plot legened", action="append", dest="ids")
    parser.add_option('--infer-ids', help="infer the ids from the file names", action="store_true", default=True, dest="infer_ids")
    parser.add_option('--xlim', help="the x-axis range", dest="xlim")
    parser.add_option('--ylim', help="the y-axis range", dest="ylim")
    parser.add_option('--out', help="the output file", dest="out")
    parser.add_option('--outtype', help="the output file type (ex. pdf, png)", dest="outtype", default="pdf")
    parser.add_option('--min-mapq', help="minimum mapping quality (inclusive)", dest="min_mapq", type=int, default=0)
    parser.add_option('--no-legend', help="do not add a legend to the plot", dest="no_legend", action="store_true", default=False)
    if len(sys.argv[1:]) < 1:
        parser.print_help()
    else:
        options, args = parser.parse_args()
        main(options)
