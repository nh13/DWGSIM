#!/usr/bin/env python

"""
Runs DWGSIM

usage: dwgsim_wrapper.py [options]
    -e,--errorOne=e: base/color error rate of the first read [0.020]
    -E,--errorTwo=E: base/color error rate of the second read [0.020]
    -d,--innertDist=d: inner distance between the two ends [500]
    -s,--stdev=s: standard deviation [50]
    -N,--numReads=N: number of read pairs [1000000]
    -1,--lengthOne=1: length of the first read [70]
    -2,--lengthTwo=2: length of the second read [70]
    -r,--mutRate=r: rate of mutations [0.0010]
    -R,--fracIndels=R: fraction of mutations that are indels [0.10]
    -X,--indelExt=X: probability an indel is extended [0.30]
    -y,--randProb=y: probability of a random DNA read [0.10]
    -n,--maxN=n: maximum number of Ns allowed in a given read [0]
    -c,--colorSpace=c: generate reads in color space (SOLiD reads)
    -S,--strand=S: strand 0: default, 1: same strand, 2: opposite strand
    -H,--haploid=H: haploid mode
    -f,--fasta=f: the reference genome FASTA
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
    parser.add_option( '-e', '--errorOne', dest='errorOne', type='float', help='base/color error rate of the first read' )
    parser.add_option( '-E', '--errorTwo', dest='errorTwo', type='float', help='base/color error rate of the second read' )
    parser.add_option( '-d', '--innertDist', dest='innertDist', type='int', help='inner distance between the two ends' )
    parser.add_option( '-s', '--stdev', dest='stdev', type='float', help='standard deviation' )
    parser.add_option( '-N', '--numReads', dest='numReads', type='int', help='number of read pairs' )
    parser.add_option( '-1', '--lengthOne', dest='lengthOne', type='int', help='length of the first read' )
    parser.add_option( '-2', '--lengthTwo', dest='lengthTwo', type='int', help='length of the second read' )
    parser.add_option( '-r', '--mutRate', dest='mutRate', type='float', help='rate of mutations' )
    parser.add_option( '-R', '--fracIndels', dest='fracIndels', type='float', help='fraction of mutations that are indels' )
    parser.add_option( '-X', '--indelExt', dest='indelExt', type='float', help='probability an indel is extended' )
    parser.add_option( '-y', '--randProb', dest='randProb', type='float', help='probability of a random DNA read' )
    parser.add_option( '-n', '--maxN', dest='maxN', type='int', help='maximum number of Ns allowed in a given read' )
    parser.add_option( '-c', '--colorSpace', action='store_true', dest='colorSpace', default=False, help='generate reads in color space (SOLiD reads)' )
    parser.add_option( '-S', '--strand', dest='strand', type='choice', default='0', choices=('0', '1', '2'), help='strand 0: default, 1: same strand, 2: opposite strand' )
    parser.add_option( '-H', '--haploid', action='store_true', dest='haploid', default=False, help='haploid mode' )
    parser.add_option( '-f', '--fasta', dest='fasta', help='The reference genome fasta' )
    parser.add_option( '-3', '--outputBFAST', dest='outputBFAST', help='the BFAST output FASTQ' )
    parser.add_option( '-4', '--outputBWA1', dest='outputBWA1', help='the first BWA output FASTQ' )
    parser.add_option( '-5', '--outputBWA2', dest='outputBWA2', help='the second BWA output FASTQ' )
    parser.add_option( '-6', '--outputMutations', dest='outputMutations', help='the output mutations TXT' )
    
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
        os.symlink( options.fasta, reference_filepath )
        assert reference_filepath and os.path.exists( reference_filepath ), 'A valid genome reference was not provided.'
        tmp_dir = '%s/' % tempfile.mkdtemp()
        dwgsim_output_prefix = "dwgsim_output_prefix"
        dwgsim_cmd = 'dwgsim -e %s -E %s -d %s -s %s -N %s -1 %s -2 %s -r %s -R %s -X %s -y %s -n %s' % \
                (options.errorOne, \
                options.errorTwo, \
                options.innertDist, \
                options.stdev, \
                options.numReads, \
                options.lengthOne, \
                options.lengthTwo, \
                options.mutRate, \
                options.fracIndels, \
                options.indelExt, \
                options.randProb, \
                options.maxN)
        if options.colorSpace:
            dwgsim_cmd = dwgsim_cmd + ' -c'
        if options.haploid:
            dwgsim_cmd = dwgsim_cmd + ' -H'
        dwgsim_cmd = dwgsim_cmd + ' ' + options.fasta
        dwgsim_cmd = dwgsim_cmd + ' ' + tmp_dir + '/' + dwgsim_output_prefix

        # need to nest try-except in try-finally to handle 2.4
        try:
            # dwgsim
            run_process ( dwgsim_cmd, 'dwgsim', tmp_dir, buffsize )
            # Move files
            cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.mutations.txt' + ' ' + options.outputMutations
            run_process ( cmd, 'mv #1', tmp_dir, buffsize )
            cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.bfast.fastq' + ' ' + options.outputBFAST
            run_process ( cmd, 'mv #2', tmp_dir, buffsize )
            cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.bwa.read1.fastq' + ' ' + options.outputBWA1
            run_process ( cmd, 'mv #3', tmp_dir, buffsize )
            cmd = 'mv ' + tmp_dir + '/' + dwgsim_output_prefix + '.bwa.read2.fastq' + ' ' + options.outputBWA2
            run_process ( cmd, 'mv #4', tmp_dir, buffsize )
            # check that there are results in the output file
            check_output ( options.outputMutations, True )
            check_output ( options.outputBFAST, False )
            check_output ( options.outputBWA1, False )
            check_output ( options.outputBWA2, False )
            sys.stdout.write( 'DWGSIM successful' )
        except Exception, e:
            stop_err( 'DWGSIM failed.\n' + str( e ) )
    finally:
        # clean up temp dir
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )

if __name__=="__main__": __main__()
