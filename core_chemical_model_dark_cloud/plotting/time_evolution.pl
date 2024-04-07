#!/usr/bin/perl

# Author: Catherine Walsh 05/2007
# Adapted from graph4.pl by David Tideswell et al for 2D disk chemistry
# Program to extact abundances of several species at a particular grid point
# Abundances returned as a function of time
# Produces GNUPLOT friendly data and GNUPLOT script
# Produces LOGARITHMIC x axis vs LOGARITHMIC y axis graph

use strict;

use POSIX;

#Check for correct number of arguments

if(@ARGV<5)
{
    	die ("Usage: time_evolution.pl <results.dat> <output.dat> <radius> <height> <species list>
    	e.g. time_evolution.pl disk_output.dat graph.dat 1.0 0.5 H HCN CO H2O\n");
}

# Set variable $infile to first command line argument - file to extract data from
# Set variable $outfile to second command line argument - file to write data to
# Set variable $radius to third command line argument

chomp(my $infile = shift @ARGV);
chomp(my $outfile = shift @ARGV);
chomp(my $radius = shift @ARGV);
chomp(my $height = shift @ARGV);


# Cycle through remaining arguments and put into a species array
# $smax counts number of species

my @species;
my $smax = 0;
my $s;

while(<@ARGV>)
{
   	push @species,$_;
        $smax = $smax + 1;
}

$smax = $smax - 1;

# Append GNUPLOT script file name with .gp

#my $outfile1 = $outfile. ".gp";

# Append data output file name with .dat

$outfile = $outfile. ".dat";

# Open files for reading and writing

open(FILE1,"<",$infile) or die("Unable to open file - check filename!");
open(FILE2,">",$outfile);
#open(FILE3,">",$outfile1);

# Read input file - remove newlines split lines into strings and store in an array

my @r;
my @z;
my @t;
my $imax = 0;
my $i;
my $jmax = 0;
my $j;
my @ab;
my @abundance;
my @he;
my @h2;


#while (<FILE1>)
#{
#	chomp;
#        my @line = split /\s+/, $_;
#        
#        if ($line[1] eq "TIME")
#        {
#        	shift @line;
#                shift @line;
#                @t = @line;
#                $jmax = @t;
#        }
#}

$jmax = $jmax - 1;

while (<FILE1>)
{
	chomp;
        my @line = split /\s+/, $_;
        
        if ($line[1] eq "RADIUS")
        {
        	push @r, $line[3];
	        $imax = $imax + 1; 
        }

	elsif ($line[1] eq "HEIGHT")
        {
        	push @z, $line[3];
        }
        
        elsif ($line[1] eq "TIME")
        {
                shift @line;
                shift @line;
                @t = @line;
                $jmax = @t;
        }
        
        if ($line[1] eq 'H2')
        {
		shift @line;
        	shift @line;
                
         	for $j (0..$jmax)
         	{
           	    	$h2[$imax-1][$j] = $line[$j];
		}
                	
        }
                
        if ($line[1] eq 'He')
        {
		shift @line;
        	shift @line;
                
         	for $j (0..$jmax)
         	{
           	    	 $he[$imax-1][$j] = $line[$j];
			}


        }
        
        
        for $s (0..$smax)
        {
        	
                if ($line[1] eq $species[$s])
                {
                        @ab = @line;
                        shift @ab;
                        shift @ab;
                        
                        for $j (0..$jmax)
                        {
                        	$abundance[$j][$s][$imax-1] = $ab[$j];
                	}
                }
        }
}

$imax = $imax - 1;
$jmax = $jmax - 1;

close FILE1;

# Set input radius and height to closest value in output file

my $rad = 0;
my $h = 0;

for $i (0..$imax)
{
	                
        if ($r[$i] <= $radius)
        {
        	$rad = $r[$i];    
        }
	        
        if (($r[$i]==$rad) && ($z[$i] <= $height))
        {
        	$h = $z[$i];
                
        } 
        
        for $j (0..$jmax)
        {
#        	my $ratio = 1-$h2[$i][$j]+$he[$i][$j];
                my $ratio = $h2[$i][$j];

        	for $s (0..$smax)
                {
                	$abundance[$j][$s][$i] = $abundance[$j][$s][$i]/$ratio;
                } 
                
                $h2[$i][$j] = $h2[$i][$j]/$ratio;  
                $he[$i][$j] = $he[$i][$j]/$ratio;       
     
        }
        
}

# Print data file

print FILE2 "# Disk Chemistry Output \n \n";
print FILE2 "# $outfile in GNUPLOT format produced by $0 \n \n";
print FILE2 "# Time evolution of abundances for ";

foreach (@species)
{
	print FILE2 "$_ ";
}

print FILE2 "at radius = $rad AU and height = $h AU \n \n";

print FILE2 "# TIME (YRS) \t";
print FILE2 "H2 \t\t";
print FILE2 "He \t\t";

foreach (@species)
{
	print FILE2 "$_ \t \t";
}

print FILE2 "\n \n";


for $i (0..$imax)
{        
        if (($r[$i] == $rad) && ($z[$i] == $h))
        {
                for $j (0..$jmax)
                {
                	print FILE2 "$t[$j]\t";
                        
                        $h2[$i][$j] = sprintf('%9.3e',$h2[$i][$j]);
                        $he[$i][$j] = sprintf('%9.3e',$he[$i][$j]);
                        
                        print FILE2 "$h2[$i][$j]\t";
                        print FILE2 "$he[$i][$j]\t";
                        
                        for $s (0..$smax)
                        {	                        
                        	$abundance[$j][$s][$i] = sprintf('%9.3e',$abundance[$j][$s][$i]);
                                print FILE2 "$abundance[$j][$s][$i]\t";
                        }
                        print FILE2 "\n";
                }

        }        
}


close FILE2;

# 
# 
# # Now print output for GNUPLOT script
# # Print to FILE3 -  the .gp file
# 
# 
# # Values for zmin and zmax and amin and amax - start by initilising
# 
# 
# my $tmin = 1;
# my $tmax = $t[$jmax];
# my $amin = $abundance[0][0][0];
# my $amax = $abundance[0][0][0];
# 
# for $i (0..$imax)
# {
# 	if (($r[$i] == $rad) and ($z[$i] == $h))
#         {
# 	
#         	for $j (0..$jmax)
#                 {
#         		for $s (0..$smax)
#         		{
# 				if ($abundance[$j][$s][$i] > $amax)
#         	        	{
# 					$amax = $abundance[$j][$s][$i];
# 				}
#         		
# 				if (($abundance[$j][$s][$i] < $amin) && ($abundance[$j][$s][$i] > 0))
# 				{
# 					$amin = $abundance[$j][$s][$i];
#             			}        
#                        }         
#         	}	
#         }
# }
# 
# 
# # Round the limits to plot logscale
# # No log_10 in perl so use log_e with identity
# 
# $tmin = 10.0**(floor(log($tmin)/log(10.0)));
# $tmax = 10.0**(ceil(log($tmax)/log(10.0)));
# $amin = 10.0**(floor(log($amin)/log(10.0)));
# $amax = 10.0**(ceil(log($amax)/log(10.0)));
# 
# if ($amin < 1.0e-13)
# {
# 	$amin = 1.0e-13;
# }
# 
# Print GNUPLOT commands

#&GNUPLOT;

#print FILE3 "plot \"$outfile\" u 1:4 t \"$species[0]\" w lines lw 4, \\\n";
		
#for $s (1..$smax)
#{
#	my $p = $s + 4;
	
#        if ($s < $smax) 
#        {
#        	print FILE3 "     \"$outfile\" u 1:$p t \"$species[$s]\" w lines lw 4, \\\n";
#        }#
#	
#        elsif ($s == $smax) #
#	{
#		print FILE3 "     \"$outfile\" u 1:$p t \"$species[$s]\" w lines lw 4 \n";
#        }
#}

#close FILE3;

#print "GNUPLOT data saved to $outfile \n";
#print "GNUPLOT script saved as $outfile1 \n";

# Call GNUPLOT

#`gnuplot $outfile1`;

#sub GNUPLOT 
#{
#print FILE3 <<_EOT_

#\# GNUPLOT $outfile1

#set title \" Radius = $rad AU and Height = $h AU \"
#set term postscript enhanced eps colour 20 \"Helvetica\"
#set out \"$outfile.eps\"

#set logscale x
#set xrange [$tmin:$tmax]
#set format x "10^{%T}" 
#set mxtics 10

#set logscale y
#set yrange [$amin:$amax]
#set format y "10^{%T}" 
#set ytics mirror $amin,10,$amax
#set mytics 10

#set key top left

#set xlabel \"Time (Yrs)\"
#set ylabel \"Fractional Abundance\"

#_EOT_
#}








