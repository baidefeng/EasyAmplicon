
1. Initialize the working environment and create necessary directories.

  wd=/d/Download/Easyamplicon2
  db=/d/Download/Easymicrobiome
  PATH=$PATH:${db}/win
  cd ${wd}
  mkdir -p seq result temp 

2. Generate the feature table and representative sequences for amplicon data analysis.
2.1 short-read amplicon sequencing: prepare seqeucning data in seq/*.fq.gz and metadata.txt in result/

  bash easyamplicon.sh feature_table -d ${db} \
      --raw-reads-dir seq \
      --metadata-file result/metadata.txt \
      --left-strip 29 --right-strip 18 \
      --chimera-ref-db ${db}/usearch/rdp_16s_v18.fa \
      --sintax-db ${db}/usearch/rdp_16s_v18.fa \
      --min-unique-size 10 \
      --output-dir result/ \
      --otutab-filename otutab.txt \
      --otus-filename otus.fa 

2.2 Pacbio amplicon sequencing 

  bash easyamplicon.sh feature_pacbio --dbdir /mnt/d/EasyMicrobiome \
    --raw-reads-dir Pacbi/seq \
    --metadata-file Pacbio/result/metadata.txt \
    --temp-dir Pacbio/temp \
    --output-dir Pacbio/result \
    --fastq-qmax 93 \
    --sintax-db /mnt/c/EasyMicrobiome/usearch/sintax_defalut_emu_database.fasta 
    
2.3 Nanopore amplicon sequencing

  bash easyamplicon.sh feature_nanopore --dbdir /mnt/d/EasyMicrobiome \
    --raw-reads-dir Nanopore/seq \
    --metadata-file Nanopore/result/metadata.txt \
    --temp-dir Nanopore/temp \
    --output-dir Nanopore/result \
    --fastq-qmax 93 \
    --sintax-db /mnt/c/EasyMicrobiome/usearch/sintax_defalut_emu_database.fasta

3. Filter and normalize the generated feature table

  # 过滤不需要otus.fa吧，过滤值也需要根据前面结果推荐 
  bash easyamplicon.sh feature_filter -d ${db} \
      --input-otutab-file result/otutab.txt \
      --input-otus-file result/otus.fa \
      --depth 10000 \
      --output-otutab-stat-file result/otutab.stat \
      --output-rarefied-otutab-file result/otutab_rare.txt

4. Perform taxonomic annotation to classify the representative sequences.

  bash easyamplicon.sh taxonomy -d ${db} \
    --sintax-input result/otus.sintax \
    --tax-output-dir result/tax


5. Calculate alpha diversity indices to assess within-sample diversity.

  bash easyamplicon.sh diversity_alpha -d ${db} \
    --rarefied-otu-table result/otutab_rare.txt \
    --feature-seqs result/otus.fa \
    --alpha-diversity-output result/alpha/alpha.txt \
    --alpha-rare-output result/alpha/alpha_rare.txt

6. Analyze beta diversity to compare microbial community structures between samples.

  bash easyamplicon.sh diversity_beta -d ${db} \
    --rarefied-otu-table result/otutab_rare.txt \
    --otus-tree result/otus.tree \
    --beta-output-dir result/beta
    
7. Differential analysis to identify statistically significant differences between groups.
 
  bash easyamplicon.sh compare -d ${db} \
    --feature-table result/otutab.txt \
    --taxonomy-file result/taxonomy.txt \
    --tax-summary-dir result/tax \
    --group "Group" \
    --compare "KO-WT" \
    --compare-output-dir result/compare \
    --stamp-output-dir result/stamp \
    --lefse-output-dir result/lefse


8. Predict functional profiles of the microbial communities.
 
  bash easyamplicon.sh function -d /mnt/mnt/d/Easymicrobiome \
    --input-filtered-fa temp/filtered.fa \
    --group "Group" \
    --output-bugbase-dir result/bugbase/ 
 
 
  


Modular pipeline minimal dependencies (easyamplicon.sh)

If you only want to run the modular scripts (`easyamplicon.sh` + `pipeline_modules/`), follow:

  EASYAMPLICON_MODULAR_INSTALL.md

This distilled guide covers the minimum required tools (recommended: WSL + conda),
EasyMicrobiome path setup (`EASY_MICROBIOME` / `--dbdir`), and quick verification commands.