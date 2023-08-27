#This script does two things:
#creates the files for upload to KEGG Mapper – Color” tool
#counts the no. of genes per pathway per cluster

#Regarding creating files for KEGG mapper-color tool
#this script is less comprehensive than the previous ones, and only colors by cluster
#unlike the previous ones, it uses the results of Vered's script and not Liron's DeSeq2 result files

#Regarding the counting of the no. of genes per pathway in each cluster:
#The result files contain columns for the 8 clusters, 4 unified clusters, and “up” (clusters 1-4) and “down” (clusters 5-8) clusters.
#Pathways are in the rows and values are either:
#1. no. of unique ensembl genes
#2. List of unique ensemble gene IDs
#3. List of unique gene symbols from ensemble annotation
#4. no. of unique KEGG genes
#5. List of unique KEGG gene IDs
#6. List of unique gene symbols from KEGG annotation
#(one file for each value type)

use strict;
use warnings;

#input files - annotation

my $gene2pathway_file           = "From_KEGGREST_API/Mouse_gene2pathway.txt";
my $pathway2name_file           = "From_KEGGREST_API/Mouse_pathway2name.txt";
my $keggGene_info_file          = "From_KEGGREST_API/Mouse_genes_info.txt";
my $ensembl_from_biomart_file   = "From_BioMart/mart_export_21.5.2023.txt";
my $KEGG_gene2ensembl_gene_file = "Functional_annotation_results/KEGG_gene2ensembl_gene.txt"; #result file from the perl script "Retrieve_KEGG_annotation.pl"

#input files - gene2cluster (from RNA-Seq analysis results)

my @datasets = ("DatasetA", "DatasetB", "DatasetC");  #modify to reflect your dataset (e.g. analysis round, tissue etc.) names

my $input_dir = "From_RNA-Seq_analysis/";  #files from RNA-Seq analysis that map gene 2 cluster

my %dataset_files;
$dataset_files{"DatasetA"} = $input_dir . "Gene2Cluster_DatasetA.txt";
$dataset_files{"DatasetB"} = $input_dir . "Gene2Cluster_DatasetB.txt";
$dataset_files{"DatasetC"} = $input_dir . "Gene2Cluster_DatasetC.txt";

#output dir

my $output_dir = "KEGG_pathways_in_clustering_results/";
my $for_kegg_maps_dir  = $output_dir . "for_KEGG_maps/";
my $summary_stats_dir  = $output_dir . "summary_genes_in_clusters/";

#define colors

#see color examples in the following links
#(https://www.w3schools.com/tags/ref_colornames.asp)
#(https://www.w3schools.com/colors/colors_groups.asp)
#(https://www.htmlcsscolor.com/hex/5CD6FF)

#8 clusters
my %color  = ("1" => "#FFBDBD", # medium red
              "2" => "#FFDDDD", # light red
              "3" => "#FF7D7D", # dark red
              "4" => "#B88AA1", # violet
              "5" => "#97FF97", # medium green
              "6" => "#CAFFCA", # light green
              "7" => "#33FF33", # dark green
              "8" => "#91A486");# olive green

my %unified_clusters = ("1" => "persistent",
                        "5" => "persistent",
                        "2" => "reversed",
                        "6" => "reversed",
                        "3" => "aggravated",
                        "7" => "aggravated",
                        "4" => "reverse_effect",
                        "8" => "reverse_effect");

my %up_down_clusters = ("1" => "up", 
                        "2" => "up", 
                        "3" => "up", 
                        "4" => "up", 
                        "5" => "down", 
                        "6" => "down", 
                        "7" => "down", 
                        "8" => "down");

my @cluster_order = (1..8, "persistent", "reversed", "aggravated", "reverse_effect", "up", "down");

#data structures

my $ens2kegg;
my $ens2path;
my $kegg2path;
my $path2kegg;
my $path2name;
my $kegg2info;
my $ens2info;
my $struct;

#general info
#############

#open general files

open (KEGG2ENS,  $KEGG_gene2ensembl_gene_file) or die $!;
open (GENE2PATH, $gene2pathway_file) or die $!;
open (PATH2NAME, $pathway2name_file) or die $!;
open (KEGG2INFO, $keggGene_info_file) or die $!;
open (ENS2INFO,  $ensembl_from_biomart_file) or die $!;

while (<KEGG2ENS>) {
	chomp;
	my ($kegg_id, $ens_id) = split (/\t/);
	$ens2kegg->{$ens_id}->{$kegg_id}=1;
}

while (<GENE2PATH>) {
	chomp;
	my ($kegg_gene_id, $kegg_path) = split (/\t/);
	$kegg_path =~ s/^path://;
	$kegg2path->{$kegg_gene_id}->{$kegg_path} = 1;
	$path2kegg->{$kegg_path}->{$kegg_gene_id} = 1;
}
while (<PATH2NAME>) {
	chomp;
	my ($kegg_path, $path_description) = split (/\t/);
	$path_description =~ s/ \- Mus musculus \(house mouse\)$//;
	$path2name->{$kegg_path} = $path_description;
}
while (<KEGG2INFO>) {
	chomp;
	my ($kegg_gene_id, $gene_type, $gene_loc, $gene_desc) = split (/\t/);
	my ($gene_symbols, $gene_title);
	if ($gene_desc =~ /;/) {
		($gene_symbols, $gene_title) = split (/; /, $gene_desc);
	} else {
		$gene_symbols = "";
		$gene_title = $gene_desc;
	}
	$kegg2info->{$kegg_gene_id}->{"symbols"} = $gene_symbols;
	$kegg2info->{$kegg_gene_id}->{"title"} = $gene_title;	
}

my $header1 = <ENS2INFO>;
while(<ENS2INFO>) {
	chomp;
	my ($Gene_stable_ID, $Gene_name, $Gene_description, $Gene_type) = split (/\t/);
	$ens2info->{$Gene_stable_ID}->{"symbol"} = $Gene_name;
	$ens2info->{$Gene_stable_ID}->{"title"} = $Gene_description;
}

close (KEGG2ENS);
close (GENE2PATH);
close (PATH2NAME);
close (KEGG2INFO);
close (ENS2INFO);

#clusters info
##############

my $dataset;
foreach $dataset (@datasets) {
	print $dataset, "\n";
	my $dataset_file = $dataset_files{$dataset};
	my $result_file  = $for_kegg_maps_dir . $dataset . "_for_KEGG_by_cluster.txt";
	
	#data structure
	my $kegg2cluster;
	
	open (IN, $dataset_file);
	my $header2 = <IN>; 
	while (<IN>) {
		chomp;
		my ($ensembl_id, $cluster, @rest) = split (/\t/);
		my $unified_cluster = $unified_clusters{$cluster};
		my $up_down_cluster = $up_down_clusters{$cluster};
		
		my (@keggids) = keys (%{$ens2kegg->{$ensembl_id}});
		my $keggid;
		foreach $keggid (@keggids) {
			
			#for colored kegg maps
			$kegg2cluster->{$keggid}->{$cluster}=1;  
			
			#for summary
			my @paths = sort keys (%{$kegg2path->{$keggid}});
			my $path;
			foreach $path (@paths) {
				$struct->{$dataset}->{"ens"}->{$path}->{$cluster}->{$ensembl_id}=1;
				$struct->{$dataset}->{"ens"}->{$path}->{$unified_cluster}->{$ensembl_id}=1;
				$struct->{$dataset}->{"ens"}->{$path}->{$up_down_cluster}->{$ensembl_id}=1;
				$struct->{$dataset}->{"kegg"}->{$path}->{$cluster}->{$keggid}=1;
				$struct->{$dataset}->{"kegg"}->{$path}->{$unified_cluster}->{$keggid}=1;
				$struct->{$dataset}->{"kegg"}->{$path}->{$up_down_cluster}->{$keggid}=1;
			}
		}		  
	}
	close (IN);
	
	#print files for "kegg map - color" tool
	
	open (OUT, ">$result_file") or die $!;
	my $keggid;
	foreach $keggid (keys %{$kegg2cluster}) {
		my @clusters = keys %{$kegg2cluster->{$keggid}};
		my $cluster;
		foreach $cluster (@clusters) {
			print OUT join ("\t", $keggid, $color{$cluster}), "\n";
		}
	}
	print "Done printing files for kegg colored maps\n";
	close (OUT);	
	
	#print summary files
	
	my $id_type;
	foreach $id_type ("ens", "kegg") {
		print "Printing summary files for $id_type\n";
		my $file_path_start = $summary_stats_dir . $dataset . "_" . $id_type;
		my $summary_file          = $file_path_start . "_summary.txt";
		my $summary_file_gene_ids = $file_path_start . "_gene_IDs.txt";
		my $summary_file_symbols  = $file_path_start . "_gene_symbols.txt";
		open (SUMMARY, ">$summary_file") or die $!;
		open (GENE_IDS, ">$summary_file_gene_ids") or die $!;
		open (SYMBOLS, ">$summary_file_symbols") or die $!;
		my $summary_header = join ("\t", "Pathway", "Cluster", @cluster_order) . "\n";
		print SUMMARY  $summary_header;
		print GENE_IDS $summary_header;  
		print SYMBOLS  $summary_header;
		my $path;
		foreach $path (sort keys (%{$path2name})) {
			my $path_name = $path2name->{$path};
			my @items_summary = ();
			my @items_gene_IDs = ();
			my @items_symbols = ();
			if ($struct->{$dataset}->{$id_type}->{$path}) {
				my $cluster;
				foreach $cluster (@cluster_order) {
					my @geneIDs = sort keys (%{$struct->{$dataset}->{$id_type}->{$path}->{$cluster}});
					my %gene_symbols;
					my $symbol;
					my $geneID;
					foreach $geneID (@geneIDs) {
						if ($id_type eq "ens") {
							if ($ens2info->{$geneID}->{"symbol"}) {
								$symbol = $ens2info->{$geneID}->{"symbol"};
							} else {
								$symbol = $geneID;
							}
							
							$gene_symbols{$symbol}=1;
						} elsif ($id_type eq "kegg") {
							if ($kegg2info->{$geneID}->{"symbols"}) {
								my $symbols = $kegg2info->{$geneID}->{"symbols"};
								if ($symbols =~ /,/) {
									($symbol) = $symbols =~ /^(.*?),/;
								} else {
									$symbol = $symbols;
								}
								
							} else {
								$symbol = $geneID;
							}
							$gene_symbols{$symbol}=1;
						}
					}
					push (@items_summary, scalar (@geneIDs));
					push (@items_gene_IDs, join (",", @geneIDs));
					push (@items_symbols,  join (",", sort keys (%gene_symbols)));
				}
			}
			my $path_to_print = "path:" . $path;
			print SUMMARY  join ("\t", $path_to_print, $path_name, @items_summary), "\n";
			print GENE_IDS join ("\t", $path_to_print, $path_name, @items_gene_IDs), "\n";
			print SYMBOLS  join ("\t", $path_to_print, $path_name, @items_symbols), "\n";
		}
		
		close (SUMMARY);
		close (GENE_IDS);
		close (SYMBOLS);
	}
}

print "End\n";
