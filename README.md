New York Genome Center Macmillian Summer 2024 intern Daniel Cha

Contact: kwonryuc@andrew.cmu.edu
Advisor: Nyasha Chambwe

Directory Heirarchy

```
APGmapping/
├── alignmentCSVs ***(Contains outputs of "generateRefCSV()")***
│   ├── GRCh38.csv
│   ├── GRCh38_withoutAlts.csv
│   ├── HG002.csv
│   ├── HG00438.csv
...
│   ├── HG03540.csv
│   ├── HG03579.csv
│   └── T2T.csv
├── alignmentTables ***(Contains outputs of "create_alignment_table(), used for plotting chromosome ideogram")***
│   ├── alignments_with_names.xlsx
│   └── alignments.xlsx
├── analysis_code ***(Main Source Code)***
│   ├── analyse.py **(Source Code)**
│   ├── getPaths.sh **(generates txt file containing file paths)**
│   ├── csvFilePaths.txt
│   └── run.sh 
├── figures (all figures)
│   ├── aligns
│   │   ├── cumul_reasonable_chart.png
│   │   └── single_reasonable_chart.png
│   ├── barcharts
│   │   ├── 0.5_0.8CumulPassed.csv
│   │   ├── 0.5_0.8SinglePassed.csv
│   │   ├── 0.8_0.9CumulPassed.csv
│   │   ├── 0.8_0.9SinglePassed.csv
│   │   ├── CumulPassedchart.png
│   │   └── SinglePassedchart.png
│   ├── chromosomeAnnotote
│   │   ├── chromosome_file.txt
│   │   ├── chromosome_ideogram_with_density_heatmap_unique_regions3.png
│   │   ├── chromosome_ideogram_with_density_heatmap_unique_regions_bins.png
│   │   ├── chromosome_ideogram_with_density_heatmap_unique_regions.png
│   │   ├── chromosome_plot_with_density_heatmap_unique_regions.png
│   │   ├── plots.py **(Plotting code to create chromosome ideogram)**
│   │   └── run_script.sh
│   ├── histograms
│   │   ├── coverage_histo.png
│   │   ├── identity_histo.png
│   │   └── unmasked_CpGTracks
│   │       ├── unmasked_coverage_histo.png
│   │       ├── unmasked_CpG_coverage_histo.png
│   │       ├── unmasked_CpG_identity_histo.png
│   │       └── unmasked_identity_histo.png
│   └── scatter
│       ├── Heat_GRCh38_Contig_Scatter.png
│       ├── Heat_T2T_Contig_Scatter.png
│       └── scatter_Data
│           ├── GRCh38scatterData.csv
│           └── T2TscatterData.csv
└── logfiles
    ├── APG_T2T.sh.e7610171
    ├── APG_T2T.sh.o7610171
    ├── download.sh.e7606030
    └── download.sh.o7606030

APG_GRCh38 ***(sam files and alignment scripts for GRCh38)***
├── APG_GRCh38alignment.sam
├── APG_GRCh38.sh
├── APG_GRCh38.sh.e7619935
└── APG_GRCh38.sh.o7619935
APG_HPRC ***(sam files and alignment scripts for HPRC assemblies)***
├── HPRCSamPaths.txt
├── samFiles
│   ├── APG_HG002alignment.sam
│   ├── APG_HG00438alignment.sam
│   ├── APG_HG01891alignment.sam
...
└── setup
    ├── alignParallel.sh
    ├── downloadLinks.tsv
    ├── logs
    ├── progress.log
    ├── refPaths.txt
    └── setupData.sh
APG_T2T ***(sam files and alignment scripts for T2T-CHM13)***
├── APG_T2Talignment.bam
├── APG_T2Talignment.sam
└── APG_T2T.sh
```
