#!/bin/perl

use strict;
use warnings;

# TODO:
# getopt
# better messaging
# usage
# pod...

my $fn = shift || die "No input file";
my %to_alt = (
	"AR" => "G", "GR" => "A",
	"CY" => "T", "TY" => "C",
	"CS" => "G", "GS" => "C",
	"AW" => "T", "TW" => "A",
	"GK" => "T", "TK" => "G",
	"AM" => "C", "CM" => "A"
);

open(FH, "$fn") || die "Could not open '$fn' for reading.\n";

print STDOUT "##fileformat=VCFv4.1\n";
print STDOUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n";

my $type = -1; # 1 - SUB, 2 - INS, 3 - DEL
my $prev_name = -1;
my $prev_pos = -1;
my $del_n = 0;
my $del_name = "";
my $del_pos = -1;
my $del_ref = "";
my $del_s = -1;
my $line_n = 0;
while(defined(my $line = <FH>)) {
	if($line =~ m/^(.+)\t(\d+)\t(.+)\t(.+)\t(\d+)$/) {
        $line_n++;
		my $name = $1;
		my $pos = $2;
		my $ref = $3;
		my $alt = $4;
		my $s = $5;
		if($prev_name ne $name) {
			if(0 < $del_n) {
				printf(STDOUT "%s\t%d\t.\t%s\t.\t100\tPASS\tpl=%d\n",
					$del_name, $del_pos, $del_ref, $del_s);
			}
			# reset
			$prev_pos = -1; $del_n = 0; $del_s = -1;
			$del_name = ""; $del_pos = -1; $del_ref = "";
		}
		if("-" eq $alt) { # DEL
			if(0 < $del_n && $del_pos + $del_n == $pos && $s == $del_s) {
				$del_n++;
				$del_ref .= $ref;
			}
			else {
				if(0 < $del_n) {
					printf(STDOUT "%s\t%d\t.\t%s\t.\t100\tPASS\tpl=%d\n",
						$del_name, $del_pos, $del_ref, $del_s);
				}
				$del_n = 1;
				$del_name = $name;
				$del_pos = $pos;
				$del_ref = $ref;
				$del_s = $s;
			}
		}
		else {
			# print deletion, if it exists
			if(0 < $del_n) {
				printf(STDOUT "%s\t%d\t.\t%s\t.\t100\tPASS\tpl=%d\n",
					$del_name, $del_pos, $del_ref, $s);
				# reset
				$prev_pos = -1; $del_n = 0; $del_s = -1;
				$del_name = ""; $del_pos = -1; $del_ref = "";
			}
			if("-" eq $ref) { # INS
				printf(STDOUT "%s\t%d\t.\t.\t%s\t100\tPASS\tpl=%d\n",
					$name, $pos, $alt, $s);
			}
			else { # SUB
				if(3 != $s) {
					if(!defined($to_alt{"$ref"."$alt"})) {
						die "Could not convert $ref$alt"." on $line_n.\n";
					}
					$alt = $to_alt{"$ref"."$alt"};
				}
				printf(STDOUT "%s\t%d\t.\t%s\t%s\t100\tPASS\tpl=%d\n",
					$name, $pos, $ref, $alt, $s);
			}
		}
		$prev_name = $name;
		$prev_pos = $pos;
	}
	else {
		die "Input file is not in the proper format.\n";
	}
}
# print deletion, if it exists
if(0 < $del_n) {
	printf(STDOUT "%s\t%d\t.\t%s\t.\t100\tPASS\tpl=%d\n",
		$del_name, $del_pos, $del_ref, $del_s);
}
close(FH);
