#!/usr/bin/env python

"""
Runs DWGSIM

usage: dwgsim_wrapper.py [options]
    -e,--eRateFirstread=e: base/color error rate of the first read [0.020]
    -E,--eRateSecondread=E: base/color error rate of the second read [0.020]
    -d,--innerDist=d: inner distance between the two ends [500]
    -s,--stdev=s: standard deviation [50]
    -N,--numReads=N: number of read pairs [1000000]
    -1,--loFirstread=1: length of the first read [70]
    -2,--loSecondread=2: length of the second read [70]
    -r,--mutRate=r: rate of mutations [0.0010]
    -R,--fracIndels=R: fraction of mutations that are indels [0.10]
    -X,--indelExt=X: probability an indel is extended [0.30]
    -y,--randProb=y: probability of a random DNA read [0.10]
    -n,--maxN=n: maximum number of Ns allowed in a given read [0]
    -c,--platformType=c: generate reads for Illumina (nuc space), for SOLiD (in color space) or Ion Torrent 
    -S,--strand=S: strand 0: default, 1: same strand, 2: opposite strand
    -f,--flowOrder=f: the flow order for Ion Torrent data [(null)]
    -O,--splitOutput=O: which output needs to be reported, one: one single; two: two files; all: all output options
    -H,--haploid=H: haploid mode
    -i,--input=i: the reference genome FASTA
    -3,--outputBFAST=3: the BFAST output FASTQ
    -4,--outputBWA1=4: the first BWA output FASTQ
    -5,--outputBWA2=5: the second BWA output FASTQ
    -6,--outputMutations=6: the output mutations TXT
"""

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def run_process ( cmd, name, tmp_dir, buffsize ):
    try:
        tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
        tmp_stderr = open( tmp, 'wb' )
        proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp, 'rb' )
        stderr = ''
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        raise Exception, 'Error in \'' + name + '\'. \n' + str( e )

def check_output ( output, canBeEmpty ):
    if 0 < os.path.getsize( output ):
        return True
    elif False == canBeEmpty:
        raise Exception, 'The output file is empty:' + output

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-e', '--eRateFirstread', dest='eRateFirstread', type='float', help='base/color error rate of the first read' )
    parser.add_option( '-E', '--eRateSecondread', dest='eRateSecondread', type='float', help='base/color error rate of the second read' )
    parser.add_option( '-d', '--innerDist', dest='innerDist', type='int', help='inner distance between the two ends' )
    parser.add_option( '-s', '--stdev', dest='stdev', type='float', help='standard deviation' )
    parser.add_option( '-N', '--numReads', dest='numReads', type='int', help='number of read pairs' )
    parser.add_option( '-r', '--mutRate', dest='mutRate', type='float', help='rate of mutations' )
    parser.add_option( '-R', '--fracIndels', dest='fracIndels', type='float', help='fraction of mutations that are indels' )
    parser.add_option( '-X', '--indelExt', dest='indelExt', type='float', help='probability an indel is extended' )
    parser.add_option( '-y', '--randProb', dest='randProb', type='float', help='probability of a random DNA read' )
    parser.add_option( '-n', '--maxN', dest='maxN', type='int', help='maximum number of Ns allowed in a given read' )
    parser.add_option( '-c', '--platformType', dest='platformType', type='int', help='Platform to simulate: 0: Illumina; 1: solid; 2: Ion Torrent' )
    parser.add_option( '-S', '--strand', dest='strand', type='choice', default='0', choices=('0', '1', '2'), help='strand 0: default, 1: same strand, 2: opposite strand' )
    parser.add_option( '-f', '--flowOrder', dest='flowOrder', type='int', default='2')
    parser.add_option( '-H', '--haploid', action='store_true', dest='haploid', default=False, help='haploid mode' )
    parser.add_option( '-i', '--input', dest='input', help='The reference genome fasta' )
    parser.add_option( '-O', dest='splitoutput', type='choice', default='two',choices=('two', 'one', 'all'), help='one: one single; two: two files; all: all output options' )
    parser.add_option( '-1', '--loFirstread', dest='loFirstread', type='int', help='length of the first read' )
    parser.add_option( '-2', '--loSecondread', dest='loSecondread', type='int', help='length of the second read' )
    parser.add_option( '-6', '--outputMutations', dest='outputMutations', help='the output mutations TXT' )
    parser.add_option( '-4', '--outputBWA1', dest='outputBWA1', help='the first BWA output FASTQ' )
    parser.add_option( '-5', '--outputBWA2', dest='outputBWA2', help='the second BWA output FASTQ' )
    parser.add_option( '-3', '--outputBFAST', dest='outputBFAST', help='the BFAST output FASTQ' )
							
    (options, args) = parser.parse_args()

    # output version # of tool
    try:
        tmp = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp, 'wb' )
        proc = subprocess.Popen( args='dwgsim 2>&1', shell=True, stdout=tmp_stdout )
        tmp_stdout.close()
        returncode = proc.wait()
        stdout = None
        for line in open( tmp_stdout.name, 'rb' ):
            if line.lower().find( 'version' ) >= 0:
                stdout = line.strip()
                break
        if stdout:
            sys.stdout.write( '%s\n' % stdout )
        else:
            raise Exception
    except:
        sys.stdout.write( 'Could not determine DWGSIM version\n' )

    buffsize = 1048576

    # make temp directory for dwgsim, requires trailing slash
    tmp_dir = '%s/' % tempfile.mkdtemp()

    #'generic' options used in all dwgsim commands here

    try:
        reference_filepath = tempfile.NamedTemporaryFile( dir=tmp_dir, suffix='.fa' ).name
        os.symlink( options.input, reference_filepath )
        assert reference_filepath and os.path.exists( reference_filepath ), 'A valid genome reference was not provided.'
        tmp_dir = '%s/' % tempfile.mkdtemp()
        dwgsim_output_prefix = "dwgsim_output_prefix"
        dwgsim_cmd = 'dwgsim -e %s -E %s -d %s -s %s -N %s -1 %s -2 %s -r %s -R %s -X %s -y %s -c %s -n %s' % \
                (options.eRateFirstread, \
                options.eRateSecondread, \
                options.innerDist, \
                options.stdev, \
                options.numReads, \
                options.loFirstread, \
                options.loSecondread, \
                options.mutRate, \
                options.fracIndels, \
                options.indelExt, \
                options.randProb, \
                options.platformType, \
                options.maxN)
        if options.haploid:
            dwgsim_cmd = dwgsim_cmd + ' -H'
        if options.platformType == "2":
            dwgsim_cmd = dwgsim_cmd + options.flowOrder
        dwgsim_cmd = dwgsim_cmd + ' ' + options.input
        dwgsim_cmd = dwgsim_cmd + ' ' + tmp_dir + dwgsim_output_prefix
        
        # temporary file to check the command being launched
        tmpfile = open('/tmp/dwgsim.log', 'w')
        tmpfile.write("%s" % dwgsim_cmd)

        # need to nest try-except in try-finally to handle 2.4
        try:
			# launch dwgsim + accept mutations.txt output
			run_process ( dwgsim_cmd, 'dwgsim', tmp_dir, buffsize )
			cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.mutations.txt' + ' ' + options.outputMutations
			run_process ( cmd, 'mv #1', tmp_dir, buffsize )
			check_output ( options.outputMutations, True )
			# check which output required and move files accordingly
			tmpfile.write("\nsplitoutput: %s \ncondition one: %s \ncondition two: %s \ncondition all: %s\n" % (options.splitoutput, str((options.splitoutput == 'one')), str((options.splitoutput == 'two')), str((options.splitoutput == 'all'))))
			if  any( [options.splitoutput == 'two' , options.splitoutput == 'all'] ) :	# two files
				cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.bwa.read1.fastq' + ' ' + options.outputBWA1
				run_process ( cmd, 'mv #3', tmp_dir, buffsize )
				cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.bwa.read2.fastq' + ' ' + options.outputBWA2
				run_process ( cmd, 'mv #4', tmp_dir, buffsize )
				tmpfile.write("\nmoved!!: %s" % options.splitoutput)
				# check that there are results in the output file
				check_output ( options.outputBWA1, False )
				check_output ( options.outputBWA2, False )
			if any( [options.splitoutput == 'one' , options.splitoutput == 'all'] ):
				cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.bfast.fastq' + ' ' + options.outputBFAST
				run_process ( cmd, 'mv #2', tmp_dir, buffsize )
				# check that there are results in the output file
				check_output ( options.outputBFAST, False )
			sys.stdout.write( 'DWGSIM successful' )
        except Exception, e:
			stop_err( 'DWGSIM failed.\n' + str( e ) )
    finally:
        # clean up temp dir
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )

if __name__=="__main__": __main__()
