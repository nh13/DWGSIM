#!/usr/bin/perl 
# Wrapper script: Joachim.Jacob@vib.be
# dwgsim: Nils Homer
 
 
use strict;
 
# Initiating
##########################
open(TMP,">/tmp/dwgsim.log");
my $passed_string = join (" ",@ARGV);
print TMP $passed_string,"\n\n";
 
# Parsing parameters
#########################
my $tempDir=pop(@ARGV);
my $noOutputFile=pop @ARGV;
print TMP "Nooutputfiles: ",$noOutputFile,"\n\n";
 
my @outputfiles;
for (my $i=0;$i<4;$i++){
	my $tmp = pop(@ARGV);
	push(@outputfiles,$tmp);
	print TMP "Outputfile: $tmp\n\n";
}
 
# Prepare
#########################
my $cmd_string = join (" ",@ARGV);
my @letters=('a'..'z','A'..'Z','0'..'9');
my $noLetters = @letters;
my $rand;
for (my $j=0;$j<6;$j++){ $rand.= $letters[rand $noLetters];}
my $workDir = $tempDir."/tmp".$rand;
print TMP "Curdir=$workDir\n\n";
 
# Execute dwgsim
#########################
$cmd_string .= " dwgsim_galaxy";  # dwgsim_galaxy is the prefix for the output
print TMP $cmd_string,"\n";
system("mkdir $workDir; cd $workDir; /home/galaxy/galaxy-dist/tools/bitstools/discard_stderr_wrapper.sh dwgsim $cmd_string");
# /home/galaxy/galaxy-dist/tools/bitstools/discard_stderr_wrapper.sh 
 
# Reformat results
#########################
my $file = pop(@outputfiles);
my $mvcmd= "mv -f ".$workDir."/dwgsim_galaxy.bfast.fastq ".$file;
print TMP "\nMove files: $mvcmd\n\n";
system("$mvcmd");
$file = pop(@outputfiles);
$mvcmd="mv -f ".$workDir."/dwgsim_galaxy.bwa.read1.fastq ".$file;
system("$mvcmd");
$file = pop(@outputfiles);
$mvcmd="mv -f ".$workDir."/dwgsim_galaxy.bwa.read2.fastq ".$file;
system("$mvcmd");
$file = pop(@outputfiles);
$mvcmd="mv -f ".$workDir."/dwgsim_galaxy.mutations.txt ".$file;
system("$mvcmd");
 
close(TMP);
# exit
########
exit 0;
