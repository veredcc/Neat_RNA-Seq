#Author: Vered Caspi
#Date:   31.5.2023
#This script maps between Ensembl IDs and KEGG IDs:
#A. Directly, using a conversion file Ensembl-KEGG from KEGG REST API
#B. Indirectly through uniprot, using a conversion file Ensembl2uniprot
#   from Uniprot web site and a conversion file uniprot2kegg from KEGG REST API
#Rules are applied as follows:
#1. Start from all KEGG gene IDs from KEGG info file (per organism) retrieved from KEGG REST API
#2. If there is a direct conversion to Ensembl (from A above) take all Ensembl IDs.
#   If among them are Ensembl IDs with the same gene symbol as the KEGG ID, take only them.
#   If none of them have the same gene symbl as the KEGG ID, take all of them.
#3. If there is no Ensembl assigned to them in the direct file (from A above),
#   look at the indirect conversion via uniprot (B above).
#   If among them are Ensembl IDs with the same gene symbol as the KEGG ID, take only them.
#   If none of them have the same gene symbl as the KEGG ID, take all of them.
#   For the retrieved Ensembls, take only those ensembls which are not assigned to another
#   KEGG ID in the direct assignment from KEGG.

#The script produces two result files:
#1. A detailed file for QA purposes
#2. A 2-column conversion file with KEGG ID and Ensembl ID
#In addition, it prepares the necessary files for KEGG enrichment analysis in R

#The conversion file is used for:
#1. preparing files for upload to KEGG mapping-color tool
#2. counting the no. of genes per cluster per pathway
#3. enrichment testing in R 

#Instructions how to retrieve the script's input files from KEGG REST API:
#Get gene to pathway in mouse:
#https://rest.kegg.jp/link/pathway/mmu
#
#Get pathway to name in mouse:
#https://rest.kegg.jp/list/pathway/mmu
#
#Get all mouse genes with their type, name and chromosomal location
#https://rest.kegg.jp/list/mmu
#
#Get all KEGG mouse genes and their ensembl ID
#http://rest.genome.jp/link/ensembl/mmu
#
#KEGG Rest instructions:
#https://www.genome.jp/linkdb/linkdb_api.html
#or use the web page:
#https://www.genome.jp/linkdb/
#
#Get all KEGG mouse genes and their uniprot
#https://rest.kegg.jp/conv/uniprot/mmu

use strict;
use warnings;

#input files
my $gene2pathway_file              = "From_KEGGREST_API/Mouse_gene2pathway.txt";
my $pathway2name_file              = "From_KEGGREST_API/Mouse_pathway2name.txt";
my $keggGene_info_file             = "From_KEGGREST_API/Mouse_genes_info.txt";
my $kegg2ensembl_file_from_KEGG    = "From_KEGGREST_API/Mouse_gene2ensembl.txt";
my $ensembl_from_biomart_file      = "From_BioMart/mart_export_21.5.2023.txt";
my $enst2ensg_file_from_ensmart    = "From_BioMart/mart_export_ensg2enst_28.5.2023.txt";
my $uniprot2enst_from_uniprot_file = "From_Uniprot/uniprot2ensembl_from_uniprot_28.5.2023.tsv";
my $uniprot2kegg_file_from_kegg    = "From_KEGGREST_API/Mouse_genes2uniprot.txt";

#output dir
my $output_dir = "Functional_annotation_results/";

#result files

my $all_info_for_QA_file      = $output_dir . "KEGG_gene2all_info.txt";     #for QA
my $kegg_gene_to_ensembl_file = $output_dir . "KEGG_gene2ensembl_gene.txt"; #for perl script to prepare files for KEGG map and color
my $kegg_path_to_ensembl_file = $output_dir . "KEGG_pathway2gene.tab";   #for R enrichment analysis
my $kegg_path_to_name_file    = $output_dir . "KEGG_pathway2name.tab";   #for R enrichment analysis

#general info
#############

#open general files
open (GENE2PATH,     $gene2pathway_file) or die $!;
open (PATH2NAME,     $pathway2name_file) or die $!;
open (GENE2NAME,     $keggGene_info_file) or die $!;
open (KEGG2ENS,      $kegg2ensembl_file_from_KEGG) or die $!;
open (MART,          $ensembl_from_biomart_file) or die $!;
open (ENST2ENSG,     $enst2ensg_file_from_ensmart) or die $!;
open (UP2ENST,       $uniprot2enst_from_uniprot_file) or die $!;
open (UP2KEGG,       $uniprot2kegg_file_from_kegg) or die $!;

#data structures
my $gene2path;
my $path2genes;
my $path2name;
my $KeggID2info;
my $EnsID2infoBM;  #BM for BioMart

my %enst2ensg; #I checked and the relation is 1 to many
my $up2ens;
my $kegg2up;

my $kegg2ens_direct;
my $ens2kegg_direct;
my $kegg2ens_via_up;
my $kegg2ens_combined;
my $path2ens;

#read and store general info
############################

while (<GENE2PATH>) {
	chomp;
	my ($kegg_gene_id, $kegg_path) = split (/\t/);
	$kegg_path =~ s/^path://;
	$gene2path->{$kegg_gene_id}->{$kegg_path} = 1;
	$path2genes->{$kegg_path}->{$kegg_gene_id} = 1;
}
while (<PATH2NAME>) {
	chomp;
	my ($kegg_path, $path_description) = split (/\t/);
	$path_description =~ s/ - Mus musculus (house mouse)$//;
	$path2name->{$kegg_path} = $path_description;
}
while (<GENE2NAME>) {
	chomp;
	my ($kegg_gene_id, $gene_type, $gene_loc, $gene_desc) = split (/\t/);
	my ($gene_symbols, $gene_title);
	if ($gene_desc =~ /;/) {
		($gene_symbols, $gene_title) = split (/; /, $gene_desc);
	} else {
		$gene_symbols = "";
		$gene_title = $gene_desc;
	}
	$KeggID2info->{$kegg_gene_id}->{"symbols"} = $gene_symbols;
	$KeggID2info->{$kegg_gene_id}->{"title"} = $gene_title;	
}
while (<KEGG2ENS>) {
	chomp;
	my ($kegg_id, $ens_id, $orig) = split (/\t/);
	$ens_id =~ s/^ensembl://;
	$kegg2ens_direct->{$kegg_id}->{$ens_id}=1;
	$ens2kegg_direct->{$ens_id}->{$kegg_id} = 1;
}
my $header1 = <MART>;
while(<MART>) {
	chomp;
	my ($Gene_stable_ID, $Gene_name, $Gene_description, $Gene_type) = split (/\t/);
	$EnsID2infoBM->{$Gene_stable_ID}->{"symbol"} = $Gene_name;
	$EnsID2infoBM->{$Gene_stable_ID}->{"title"} = $Gene_description;
}

#create ensembl to kegg via uniprot
###################################

#read files

while (<UP2KEGG>) {
	chomp;
	my ($up, $kegg) = split (/\t/);
	$up =~ s/^up://;
	$kegg2up->{$kegg}->{$up}=1
}

my $header2 = <ENST2ENSG>;
while (<ENST2ENSG>) {
	chomp;
	my ($ensg, $ensg_version, $enst, $enst_version) = split (/\t/);
	$enst2ensg{$enst} = $ensg;
}

my $header3 = <UP2ENST>;
while (<UP2ENST>) {
	chomp;
	my ($up, $Reviewed, $Entry_Name, $enst) = split (/\t/);
	
	if ($enst) {
		$enst =~ s/;$//;
		my @ensts = split (/;/, $enst);	
		my $enst;
		foreach $enst (@ensts) {
			$enst =~ s/ \[.*$//;
			$enst =~ s/\.\d+$//;
			if ($enst2ensg{$enst}) {
				my $ensg = $enst2ensg{$enst};		
				$up2ens->{$up}->{$ensg} = 1;				
			}	
		}		
	}	
}

#store in data structure

my $kegg;
foreach $kegg (sort keys %{$KeggID2info}) {
	if ($kegg2up->{$kegg}) {
		my $up;
		foreach $up (sort keys (%{$kegg2up->{$kegg}})) {
			if ($up2ens->{$up}) {
				my $ens;
				foreach $ens (sort keys (%{$up2ens->{$up}})) {
					$kegg2ens_via_up->{$kegg}->{$ens} = 1;
				}
			}
		}
	}
}

#close general files
close (GENE2PATH);
close (PATH2NAME);
close (MART);
close (KEGG2ENS);
close (UP2ENST);
close (ENST2ENSG);
close (UP2KEGG);

#Main combined
##############

open (OUT, ">$all_info_for_QA_file") or die $!;
print OUT join ("\t", "KEGG ID", "Orig No. EnsIDs", "Final No. EnsIDs", "Source", "KEGG Symbol", "KEGG Names",  "KEGG Title", "Ensembl ID", "EnsMart Symbol", "EnsMart Title"), "\n";

open (OUT_SHORT, ">$kegg_gene_to_ensembl_file") or die $!;

my $keggid;
foreach $keggid (sort keys %{$KeggID2info}) {
	my $symbols = $KeggID2info->{$keggid}->{"symbols"};
	my $kegg_symbol;
	if ($symbols =~ /,/) {
		($kegg_symbol) = $symbols =~ /^(.*?),/;
	} else {
		$kegg_symbol = $symbols;
	}
	my $keggid_title = "";
	$keggid_title = $KeggID2info->{$keggid}->{"title"} if ($KeggID2info->{$keggid}->{"title"});
	
	my $source;
	my $nr_ens_ids_direct = scalar (keys %{$kegg2ens_direct->{$keggid}});
	my $nr_ens_ids_via_up = scalar (keys %{$kegg2ens_via_up->{$keggid}});
	if ($nr_ens_ids_direct > 0) {  #if there is a direct ensembl (from KEGG), use it
		$source = "direct";
		my @all_ensids = sort keys %{$kegg2ens_direct->{$keggid}};
		my @filtered_ensids = ();
		my $ensid;
		foreach $ensid (@all_ensids) {
			my $ensmart_symbol = $EnsID2infoBM->{$ensid}->{"symbol"};
			if ($ensmart_symbol eq $kegg_symbol) {
				push (@filtered_ensids, $ensid);
			}
		}
		my @final_ensids;
		if (@filtered_ensids) {
			@final_ensids = @filtered_ensids;
		} else {
			@final_ensids = @all_ensids;
		}
		my $final_no_ensids1 = scalar (@final_ensids);
		foreach $ensid (@final_ensids) {
			my $ensmart_symbol = $EnsID2infoBM->{$ensid}->{"symbol"};
			my $ensmart_title  = $EnsID2infoBM->{$ensid}->{"title"};
			print OUT join ("\t", $keggid, $nr_ens_ids_direct, $final_no_ensids1, $source, $kegg_symbol, $symbols, $keggid_title, $ensid, $ensmart_symbol, $ensmart_title), "\n";
			print OUT_SHORT join ("\t", $keggid, $ensid), "\n";
			$kegg2ens_combined->{$keggid}->{$ensid} = 1;
		}					
	} elsif ($nr_ens_ids_via_up > 0) { #else, take it from the ensembl assignment through uniprot, but only if that ensembl is not assigned to another kegg
		$source = "via up";
		my @all_ensids = sort keys %{$kegg2ens_via_up->{$keggid}};
		my @filtered_ensids = ();
		my $ensid;
		foreach $ensid (@all_ensids) {
			my $ensmart_symbol = "";
			$ensmart_symbol = $EnsID2infoBM->{$ensid}->{"symbol"} if ($EnsID2infoBM->{$ensid}->{"symbol"});
			if ($ensmart_symbol eq $kegg_symbol) {
				push (@filtered_ensids, $ensid);
			}
		}
		my @final_ensids;
		if (@filtered_ensids) {
			@final_ensids = @filtered_ensids;
		} else {
			@final_ensids = @all_ensids;
		}
		my $final_no_ensids1 = scalar (@final_ensids);
		foreach $ensid (@final_ensids) {
			my $ensmart_symbol = "";
			my $ensmart_title = "";
			$ensmart_symbol = $EnsID2infoBM->{$ensid}->{"symbol"} if ($EnsID2infoBM->{$ensid}->{"symbol"});;
			$ensmart_title  = $EnsID2infoBM->{$ensid}->{"title"}  if ($EnsID2infoBM->{$ensid}->{"title"});
			my @other_keggs = sort keys (%{$ens2kegg_direct->{$ensid}}); #unnecessary but without it there is an error in the next command
			my %all_keggs_per_ens = %{$ens2kegg_direct->{$ensid}};
			unless (%all_keggs_per_ens) {  #take only those ensembls which are not assigned to another kegg id in the direct assignment from kegg
				print OUT join ("\t", $keggid, $nr_ens_ids_via_up, $final_no_ensids1, $source, $kegg_symbol, $symbols, $keggid_title, $ensid, $ensmart_symbol, $ensmart_title), "\n";
				print OUT_SHORT join ("\t", $keggid, $ensid), "\n";
				$kegg2ens_combined->{$keggid}->{$ensid} = 1;
			}
		}			
	} 
}

close (OUT);
close (OUT_SHORT);

#create files for R enrichment analysis
#######################################

#create kegg path 2 ensembl genes

open (OUT, ">$kegg_path_to_ensembl_file") or die $!;
print OUT join ("\t", "v1", "index"), "\n";

my $path;
foreach $path (sort keys %{$path2genes}) {
	my $kegg;
	foreach $kegg (sort keys %{$path2genes->{$path}}) {
		my $ens;
		foreach $ens (sort keys %{$kegg2ens_combined->{$kegg}}) {  #here is the bug
			$path2ens->{$path}->{$ens} = 1;
		}
	}
	my $ens;
	foreach $ens (sort keys %{$path2ens->{$path}}) {
		print OUT join ("\t", "path:" . $path, $ens), "\n";
	}		
}

close (OUT);

#create path 2 name
###################

open (OUT, ">$kegg_path_to_name_file") or die $!;
print OUT join ("\t", "pathway", "info"), "\n";

my $pathid;
foreach $pathid (sort keys %{$path2name}) {
	my $name = $path2name->{$pathid};
	my $path = "path:" . $pathid;
	print OUT join ("\t", $path, $name), "\n";
}
close (OUT);

print "Done\n";
