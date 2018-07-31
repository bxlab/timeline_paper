# This is a comprehensive commands list for the basic workflow for QIIME. For this to work, the working folder needs to contain the folder fastq_files, which would containt he demultiplexed read files. In the above directory, you will need params.txt, qiime_paremeters.txt, make_mapping_file.py (custom script), Silva_Database.fa, and taxonomy_7_levels.txt. Once everything is in place, feel free to copy and paste the entire file into the command line! -German :)


#Read pre-processing:
multiple_join_paired_ends.py -i fastq_files --read1_indicator _R1_001 --read2_indicator _R2_001  -o joint_output
multiple_split_libraries_fastq.py -i joint_output  -o jointsplit_output --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name


#Picking OTUs:
cp jointsplit_output/seqs.fna ./
pick_open_reference_otus.py -i seqs.fna -r ~/Desktop/16S_QIIME/Silva_Database.fa -o otu_output --suppress_step4 -m usearch61 -p ~/Desktop/16S_QIIME/params.txt -v
filter_otus_from_otu_table.py -i otu_output/otu_table_mc2_w_tax_no_pynast_failures.biom -o otu_output/otu_table_no_singletons.biom -n 2
biom convert -i otu_output/otu_table_no_singletons.biom -o otu_table.txt --to-tsv --header-key taxonomy
cp otu_output/otu_table_no_singletons.biom otu_table.biom

#Summarize the OTUs with basic figures:
summarize_taxa.py -i otu_table.biom -L 2 -o ./phylum
plot_taxa_summary.py -i phylum/otu_table_L2.txt -l phylum -c pie,bar,area -o phylum_charts
#open phylum_charts/bar_charts.html
summarize_taxa_through_plots.py -s -i otu_table.biom -o figures -f
open figures/taxa_summary_plots/bar_charts.html


#Alpha diversity metrics:
alpha_rarefaction.py -i otu_table.biom -o alpha_diversity -p ../qiime_paremeters.txt -m mapping_file.txt -f
#open alpha_diversity/alpha_rarefaction_plots/rarefaction_plots.html

#Beta diversity plots
normalize_table.py -i otu_table.biom -a CSS -o CSS_normalized_otu_table.biom
beta_diversity.py -i CSS_normalized_otu_table.biom -m bray_curtis,unweighted_unifrac,weighted_unifrac -t ./otu_output/rep_set.tre -o beta_div
#principal_coordinates.py -i ./beta_div/bray_curtis_CSS_normalized_otu_table.txt -o ./beta_div_coords.txt
principal_coordinates.py -i ./beta_div/weighted_unifrac_CSS_normalized_otu_table.txt -o ./beta_div_coords.txt
make_2d_plots.py -i beta_div_coords.txt -m metadata.txt
open beta_div_coords_2D_PCoA_plots.html


#Determine if two groups of samples are signifficantly different:
compare_categories.py --method anosim -i ./beta_div/bray_curtis_CSS_normalized_otu_table.txt -m mapping_file.txt -c Description -o anosim_out
cat anosim_out/anosim_results.txt
compare_categories.py --method adonis -i ./beta_div/bray_curtis_CSS_normalized_otu_table.txt -m mapping_file.txt -c Description -o adonis_out
cat adonis_out/adonis_results.txt

#make kronagram from the otu table
../krona_from_otu_table.py mapping_file.txt otu_table.txt
ktImportText -l -o krona_summary.html krona_summary/*
open krona_summary.html


#heatmap
upgma_cluster.py -i beta_div -o sample_cluster
../shorten_otu_table.py otu_table.txt 200 > otu_table_short.txt
make_otu_heatmap.py -i otu_table_short.txt -o heatmap.png -g png --dpi 500 -m mapping_file.txt -s sample_cluster/upgma_bray_curtis_CSS_normalized_otu_table.tre
open heatmap.png
