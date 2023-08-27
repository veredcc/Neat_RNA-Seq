use strict;
use warnings;

#KEGG annot for pathway enrichment analysis in Vered's R script
#KAAS was run on the representative sequences from the transcriptome (i.e. one representative sequence per gene)
#KEGG REST API:
#Get pathway to name (for reference pathways): https://rest.kegg.jp/list/pathway
#Get kegg gene to pathway id (for ko genes): https://rest.kegg.jp/link/pathway/ko
#optional:
#Get all kegg ko genes with their short and long names: https://rest.kegg.jp/list/ko

#In addition, the script can prepare an annotation file given a list of DE genes from the RNA-Seq analysis (or actually any list of genes)

#in files
my $transcript2ko_file = "From_KEGG_KAAS/query.ko.txt";
my $ko2name_file       = "From_KEGGREST_API/ko_to_name.txt";
my $ko2path_file       = "From_KEGGREST_API/ko_to_path.txt";
my $path2name_file     = "From_KEGGREST_API/pathway_to_name.txt";
my $DE_genes_list_file = "From_RNA-Seq_analysis/DE_genes_list.txt";

#out dir
my $out_dir = "Results";

#out files
my $path2name_file_for_R  = "$out_dir/KEGG_pathway2name.tab";
my $path2genes_file_for_R = "$out_dir/KEGG_pathway2gene.tab";
my $DE_annot_result_file  = "$out_dir/KEGG_annot_DE_genes.txt";

#data structures
my $contig2ko;
my $ko2info;
my $ko2path;
my $path2name;

#read and store contig to KO mapping data

open (CONTIG2KO, $transcript2ko_file) or die $!;
while (<CONTIG2KO>) {
	chomp;
	my ($contig, $ko) = split (/\t/);
	if ($ko) {
		$contig =~ s/_i\d+$//;  #remove isoform indication
		$contig2ko->{$contig}->{$ko} = 1;
	}
}
close (CONTIG2KO);

#read and store KEGG data
#print pathway to name file for R

#open files
open (KO2NAME, $ko2name_file) or die $!;
open (KO2PATH, $ko2path_file) or die$!;
open (PATH2NAME, $path2name_file) or die $!;
open (PATH2NAME4R, ">$path2name_file_for_R") or die $!;

while (<KO2NAME>) {
	chomp;
	my ($ko, $info) = split (/\t/);
	my ($names, $title);
	if ($info =~ /\;/) {
		($names, $title) = split (/\; /, $info);
	} else {
		$names = "";
		$title = $info;
	}
	my $ec_nr = "";
	my $title1 = $title;
	if ($title =~ /\[EC\:.*\]/) {
		($title1, $ec_nr) = $title =~ /(.*) \[(EC\:.*)\]/;
	}
	$ko2info->{$ko}->{"names"} = $names;
	$ko2info->{$ko}->{"title"} = $title1;
	$ko2info->{$ko}->{"ec"}    = $ec_nr;
}

while (<KO2PATH>) {
	chomp;
	my ($ko, $path) = split (/\t/);
	$ko =~ s/^ko\://;
	$path =~ s/^path\://;
	if ($path =~ /^map/) {
		$ko2path->{$ko}->{$path} = 1;
	}
}

print PATH2NAME4R "pathway\tinfo\n";
while (<PATH2NAME>) {
	chomp;
	my ($path, $name) = split (/\t/);
	$path2name->{$path} = $name;
	print PATH2NAME4R join ("\t", $path, $name), "\n";
}

close (KO2NAME);
close (KO2PATH);
close (PATH2NAME);
close (PATH2NAME4R);

#read DE genes list and write DE annot files

open (DE, $DE_genes_list_file) or die $!;
open (DE_ANNOT, ">$DE_annot_result_file") or die $!;

print DE_ANNOT join ("\t", "Gene", "KEGG ID", "KEGG name(s)", "KEGG description", "EC number", "Pathway ID(s)", "Pathway name(s)"), "\n";
while (<DE>) {
	chomp;
	my $gene = $_;
	my ($ko);
	my @kos;
	my @names;
	my @titles;
	my @EC_nrs;
	my @paths;
	my @path_names;
	foreach $ko (keys (%{$contig2ko->{$gene}})) {
		push (@kos, $ko);
		push (@names, $ko2info->{$ko}->{"names"});
		push (@titles, $ko2info->{$ko}->{"title"});
		push (@EC_nrs, $ko2info->{$ko}->{"ec"});
		my $path;
		foreach $path (keys (%{$ko2path->{$ko}})) {
			my $path_name = $path2name->{$path};
			push (@paths, $path);
			push (@path_names, $path_name);
		}
	} 
	print DE_ANNOT join ("\t", 
	                     $gene,
	                     join (" | ", @kos),
	                     join (" | ", @names), 
	                     join (" | ", @titles),
	                     join (" | ", @EC_nrs),
	                     join (" | ", @paths),
	                     join (" | ", @path_names)
	                     ),
	                     "\n";
	                          
}	

close (DE);
close (DE_ANNOT);


#write path to genes files for R
	
#data struct
my $path2contigs;  

my ($contig);
my ($ko);
my ($path);

foreach $contig (keys (%{$contig2ko})) {
	foreach $ko (keys (%{$contig2ko->{$contig}})) {
		foreach $path (keys %{$ko2path->{$ko}}) {
			$path2contigs->{$path}->{$contig} = 1
		}
	}
}


open (PATH2GENE, ">$path2genes_file_for_R") or die $!;

print PATH2GENE "v1\tindex\n";
foreach $path (sort keys (%{$path2contigs})) {
	my $contig;
	foreach $contig (sort keys (%{$path2contigs->{$path}})) {
		print PATH2GENE join ("\t", $path, $contig), "\n";
	}
}

close (PATH2GENE);


print "Done\n";
