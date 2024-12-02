#! /GWSPH/groups/liu_price_lab/tools/anaconda3/bin/perl
#Author : Maliha Aziz
#Date 04/03/2021
use strict;
use List::MoreUtils qw(first_index);
use List::Util qw(all);

#RDP default output + species info
my $input = shift;
my $fh;


open FILE, "$input" or die "cannot open $input for reading \n";
open ($fh, '>', "$input.c80_BIOMversion.tsv");
print $fh "#OTU ID\ttaxonomy\n";

while (<FILE>) {
	chomp;
	my @out = split ("\t", $_);
	my @new_taxa_arr;
	#print "this is the array\n";
	#print join("\n",@out),"\n";
	#print join("\n",@new_taxa_arr),"\n";

	my $dmnidx = first_index { $_ eq "domain" } @out;
	if (defined $dmnidx && $dmnidx ne '' && $dmnidx ne -1)
	{
		if ($out[$dmnidx+1] > 0.8 )
		{ 
			push @new_taxa_arr, $out[$dmnidx-1];
			push @new_taxa_arr, "Domain";
			push @new_taxa_arr, $out[$dmnidx+1];
		}
		else
		{
			print $fh "$out[0]\tUnclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified\n";
		#	print "$out[0]\tUnclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified\n";
          		next;
		}
	}
	else
	{
		push @new_taxa_arr, "Unclassified";
		push @new_taxa_arr, "domain";
		push @new_taxa_arr, "NULL";
	}
	
	#print "domain:\n";
	#print join("\n",@new_taxa_arr),"\n";
	#print("\n");

	my $phylaidx=first_index { $_ eq "phylum" } @out;
	if (defined $phylaidx && $phylaidx ne '' && $phylaidx ne -1)
	{
		if ($out[$phylaidx +1] > 0.8 )
		{ 
			push @new_taxa_arr, $out[$phylaidx-1];
			push @new_taxa_arr, "Phyla";
			push @new_taxa_arr, $out[$phylaidx +1];
		}
		else
		{
			
			print $fh "$out[0]\t$new_taxa_arr[0];Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified\n";
		#	print "$out[0]\t$new_taxa_arr[0];Unclassified;Unclassified;Unclassified;Unclassified;Unclassified\n";
          		next;
		}
	}
	else
	{
		push @new_taxa_arr, "Unclassified";
		push @new_taxa_arr, "phylum";
		push @new_taxa_arr, "NULL";
	}

	#print "phylum:\n";
	#print join("\n",@new_taxa_arr),"\n";
	#print("\n");

	my $classidx=first_index { $_ eq "class" } @out;
	if (defined $classidx && $classidx ne '' && $classidx ne -1)
	{
		if ($out[$classidx +1] > 0.8 )
		{ 
			push @new_taxa_arr, $out[$classidx-1];
			push @new_taxa_arr, "Class";
			push @new_taxa_arr, $out[$classidx +1];
		}
		else
		{
			
			print $fh "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];Unclassified;Unclassified;Unclassified;Unclassified;Unclassified\n";
		#	print "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];Unclassified;Unclassified;Unclassified;Unclassified\n";
          		next;
		}
	}
	else
	{
		push @new_taxa_arr, "Unclassified";
		push @new_taxa_arr, "class";
		push @new_taxa_arr, "NULL";
	}
	
	#print "class:\n";
	#print join("\n",@new_taxa_arr),"\n";
	#print("\n");

	my $orderidx=first_index { $_ eq "order" } @out;
	if (defined $orderidx && $orderidx ne '' && $orderidx ne -1)
	{
		if ($out[$orderidx +1] > 0.8 )
		{ 
			push @new_taxa_arr, $out[$orderidx-1];
			push @new_taxa_arr, "Order";
			push @new_taxa_arr, $out[$orderidx +1];
		}
		else
		{
			
			print $fh "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];Unclassified;Unclassified;Unclassified;Unclassified\n";
		#	print  "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];Unclassified;Unclassified;Unclassified\n";
          		next;
		}
	}
	else
	{
		push @new_taxa_arr, "Unclassified";
		push @new_taxa_arr, "order";
		push @new_taxa_arr, "NULL";
	}

	#print "order:\n";
	#print join("\n",@new_taxa_arr),"\n";
	#print("\n");

	my $familyidx=first_index { $_ eq "family" } @out;
	#printf "family: ".$familyidx;
	if (defined $familyidx && $familyidx ne '' && $familyidx ne -1)
	{
		if ($out[$familyidx +1] > 0.8 )
		{ 
			push @new_taxa_arr, $out[$familyidx-1];
			push @new_taxa_arr, "Family";
			push @new_taxa_arr, $out[$familyidx +1];
		}
		else
		{
			
			print $fh "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];Unclassified;Unclassified;Unclassified\n";
			#print  "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];Unclassified;Unclassified\n";
          		next;
		}
	}
	else
	{
		push @new_taxa_arr, "Unclassified";
		push @new_taxa_arr, "family";
		push @new_taxa_arr, "NULL";
	}

	#print "family:\n";
	#print join("\n",@new_taxa_arr),"\n";
	#print("\n");

	my $genusidx=first_index { $_ eq "genus" } @out;
	if (defined $genusidx && $genusidx ne '' && $genusidx ne -1)
	{
		if ($out[$genusidx +1] > 0.8 )
		{ 
			push @new_taxa_arr, $out[$genusidx-1];
			push @new_taxa_arr, "Genus";
			push @new_taxa_arr, $out[$genusidx +1];
		}
		else
		{
			
			print $fh "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];$new_taxa_arr[12];Unclassified;Unclassified\n";
		#	print "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];$new_taxa_arr[12];Unclassified\n";
          		next;
		}
	}
	else
	{
		push @new_taxa_arr, "Unclassified";
		push @new_taxa_arr, "genus";
		push @new_taxa_arr, "NULL";
	}

	#print "genus:\n";
	#print join("\n",@new_taxa_arr),"\n";
	#print("\n");
	
	my $speciesidx=first_index { $_ eq "species" } @out;
	if (defined $speciesidx && $speciesidx ne '' && $speciesidx ne -1)
	{
		if ($out[$speciesidx +1] > 0.8 )
		{ 
			push @new_taxa_arr, $out[$speciesidx-1];
			push @new_taxa_arr, "Species";
			push @new_taxa_arr, $out[$speciesidx +1];
			#print $fh "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];$new_taxa_arr[12];$new_taxa_arr[15];$out[$speciesidx-1]\n";
		}
		else
		{
			
			print $fh "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];$new_taxa_arr[12];$new_taxa_arr[15];$out[$speciesidx-1]<0.8\n";
			#print "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];$new_taxa_arr[12];$new_taxa_arr[15];Unclassified\n";
          		next;
		}
	}
	else
	{
		push @new_taxa_arr, "Unclassified";
		push @new_taxa_arr, "species";
		push @new_taxa_arr, "NULL";
	}

	#print "species";
	#print join("\n",@new_taxa_arr),"\n";
	#print("\n");
	print $fh "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];$new_taxa_arr[12];$new_taxa_arr[15];$new_taxa_arr[18]\n";
	#print "$out[0]\t$new_taxa_arr[0];$new_taxa_arr[3];$new_taxa_arr[6];$new_taxa_arr[9];$new_taxa_arr[12];$new_taxa_arr[15];$new_taxa_arr[18]\n";

}
