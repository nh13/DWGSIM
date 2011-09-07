#!/usr/bin/env python

"""
Runs DWGSIM_EVAL

usage: dwgsim_eval_wrapper.py [options]
    -a,--alignmentScore=a: split alignments by alignment score instead of mapping quality
    -b,--bwa=b: alignments are from BWA (SOLiD only)
    -c,--colorSpace=c: color space alignments
    -d,--scoreFactor=d: divide quality/alignment score by this factor
    -g,--wiggle=g: gap "wiggle"
    -n,--numReads=n: number of raw input paired-end reads (otherwise, inferred from all SAM records present).
    -q,--minMapq=q: consider only alignments with this mapping quality or greater.
    -z,--singleEnd=z: input contains only single end reads
    -S,--sam=s: input SAM file
    -B,--bam=B: input BAM file
    -p,--printIncorrect=p: print incorrect alignments
    -s,--numSnps=s: consider only alignments with the number of specified SNPs
    -e,--numErrors=e: consider only alignments with the number of specified errors
    -i,--indels=i: consider only alignments with indels
    -o,--output=u: The file to save the output 
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

    parser.add_option( '-a', '--alignmentScore', action='store_true', dest='alignmentScore', default=False, help='split alignments by alignment score instead of mapping quality' )
    parser.add_option( '-b', '--bwa', action='store_true', dest='bwa', default=False, help='alignments are from BWA (SOLiD only)' )
    parser.add_option( '-c', '--colorSpace', action='store_true', dest='colorSpace', default=False, help='generate reads in color space (SOLiD reads)' )
    parser.add_option( '-d', '--scoreFactor', dest='scoreFactor', type='int', help='divide quality/alignment score by this factor' )
    parser.add_option( '-g', '--wiggle', dest='wiggle', type='int', help='gap "wiggle"' )
    parser.add_option( '-n', '--numReads', dest='numReads', type='int', help='number of raw input paired-end reads (otherwise, inferred from all SAM records present).' )
    parser.add_option( '-q', '--minMapq', dest='minMapq', type='int', help='consider only alignments with this mapping quality or greater.' )
    parser.add_option( '-z', '--singleEnd', action='store_true', dest='singleEnd', default=False, help='input contains only single end reads' )
    parser.add_option( '-S', '--sam', dest='sam', default=None, help='input SAM' )
    parser.add_option( '-B', '--bam', dest='bam', default=None, help='input BAM' )
    parser.add_option( '-p', '--printIncorrect', action='store_true', dest='printIncorrect', default=False, help='print incorrect alignments' )
    parser.add_option( '-s', '--numSnps', dest='numSnps', type="int", help='consider only alignments with the number of specified SNPs' )
    parser.add_option( '-e', '--numErrors', dest='numErrors', type="int", default=False, help='consider only alignments with the number of specified errors' )
    parser.add_option( '-i', '--indels', action='store_true', dest='indels', default=False, help='consider only alignments with indels' )
    parser.add_option( '-o', '--output', dest='output', help='The file to save the output' )
    
    (options, args) = parser.parse_args()

    # output version # of tool
    try:
        tmp = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp, 'wb' )
        proc = subprocess.Popen( args='dwgsim_eval 2>&1', shell=True, stdout=tmp_stdout )
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
        sys.stdout.write( 'Could not determine DWGSIM_EVAL version\n' )

    buffsize = 1048576

    # make temp directory for dwgsim, requires trailing slash
    tmp_dir = '%s/' % tempfile.mkdtemp()

    #'generic' options used in all dwgsim commands here

    try:
        tmp_dir = '%s/' % tempfile.mkdtemp()
        dwgsim_eval_cmd = 'dwgsim_eval'
        if True == options.alignmentScore:
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' -a'
        if True == options.bwa:
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' -b'
        if True == options.colorSpace:
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' -c'
        use_p = False
        if 0 <= options.numSnps:
            use_p = True
            dwgsim_eval_cmd = dwgsim_eval_cmd + (' -s %s' % options.numSnps)
        if 0 <= options.numErrors:
            use_p = True
            dwgsim_eval_cmd = dwgsim_eval_cmd + (' -e %s' % options.numErrors)
        if True == options.indels:
            use_p = True
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' -i'
        if True == use_p or True == options.printIncorrect:
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' -p'
        if True == options.singleEnd:
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' -z'
        dwgsim_eval_cmd = '%s -d %s -g %s -n %s -q %s' % (dwgsim_eval_cmd, \
                options.scoreFactor, \
                options.wiggle, \
                options.numReads, \
                options.minMapq)
        if None != options.sam:
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' -S ' + options.sam
        elif None != options.bam:
            dwgsim_eval_cmd = dwgsim_eval_cmd + ' ' + options.bam
        else:
            raise Exception, 'Input file was neither a SAM nor BAM'
        dwgsim_eval_cmd = dwgsim_eval_cmd + ' > ' + options.output

        # need to nest try-except in try-finally to handle 2.4
        try:
            # dwgsim
            run_process ( dwgsim_eval_cmd, 'dwgsim', tmp_dir, buffsize )
            # check that there are results in the output file
            check_output ( options.output, False )
            sys.stdout.write( 'DWGSIM_EVAL successful' )
        except Exception, e:
            stop_err( 'DWGSIM_EVAL failed.\n' + str( e ) )
    finally:
        # clean up temp dir
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )

if __name__=="__main__": __main__()
