######################  DrLink pipeline (documentation)                      ###########################################
######################  under development by Hequan Sun                      ###########################################
######################  at Schneeberger Lab, MPIPZ (Koeln, Germany)          ###########################################
####       given longranger read alignments (molecule recovery), read-barcodes in vcf (barcode extraction) #############
####       and parental SNP markers, DrLink predicts CO breakpoints                                    #################
#
# step 0. compile DrLink under g++ (Debian 6.3.0-18+deb9u1) (supposing longranger > version 2.1.6 has been installed): check INSTALL file
#
# step 1. read alignment and variant calling (/barcode-assignment) using longranger (RAM 96GB in days to 2 weeks in Athal)
#
sample=sample_name
readpath=/your_read_path/
longranger wgs --id=${sample}_align_new --reference=/your_reference_index_folder/ --fastqs=${readpath} --sample=${sample}_run569_XXXXXX --library=pBlibloadXXng --jobmode=local --localcores=16 --localmem=96 --noloupe --sex=male --vcmode=freebayes >longranger_${sample}_new.log
#
# step 2: three sub processes regarding molecule recovery
#
export PATH=/your_path_to/DrLink_src/:$PATH
sample=sample_name  #
RD=25000            # maximum inter-read-distance-allowed
minMoleSize=1000    #
#
# 2.1 preprocess
samtools view phased_possorted_bam.bam | DrLink preprocess - out_pool_${sample} > out_pool_${sample}_drlin_preprocess.log
#
# 2.2 find statistics of molecules ~ 15 Gb in athal, depending sequencing depth, it may require more.
DrLink molecule out_pool_${sample}_bamInfo_trashme.gz ${RD} 1000 out_pool_${sample}_${RD} > out_pool_${sample}_drlin_molecule_${RD}.log
#
# 2.3 visualize statistics
Rscript /your_path_to/DrLink_src/Rplot_src/plotMoleculeCovLen.R out_pool_${sample}_${RD}_min${minMoleSize}bp_moleCov_stat.txt  out_pool_${sample}_${RD}_min${minMoleSize}bp_moleLen_stat.txt  out_pool_${sample}_${RD}_min${minMoleSize}bp_moleNumPerBarc_stat.txt out_pool_${sample}_${RD}_min${minMoleSize}bp_readNum_stat.txt out_pool_${sample}_${RD}_min${minMoleSize}bp_readDist_stat.txt out_pool_${sample}_${RD}
#
## step 2.4: recombination breakpoints prediction with parent1 genome as reference: ~5 Gb RAM, it may require.
##           note: one of the marker files parent1.txt and parent2.txt can be null if you only have markers defined using one parent against the reference genome 
##                 each row representing for a snp marker should be in format "label chr pos ref-base alt-base" (tab-separated), for example, "p1  1   71348   C   T"
##           run "DrLink recombis" to see details of options, set them with your own molecule/read statistically (need to consider specific molecule coverage etc; specifically for maximum settings, according to the distributions from 2.3, 1.5x the corresponding x showing the peak). 
##           result CO file like: out_RD${RD}_ks_acd${acd}_icd${icd}_${min_marker}_****_raw_BP_sorted.txt, with details (on reads and marker alleles) given in folder "out_RD${RD}_ks_acd${acd}_icd${icd}_${min_marker}__aMLxxx_aRNxx_iSCxx_iMKx_iSPxxx_iVCxx_aVCxxxx_iAFxxx_aAFx_irdx" for each CO
#
sample=sample_name
mkdir /your_path_to/${sample}_align_new/outs/splitted_vcf_mole/
mkdir /your_path_to/${sample}_align_new/outs/systematic/
cd /your_path_to/${sample}_align_new/outs/systematic/
for RD in 25000; do
    for min_score in 0.80; do
        for min_span in 1000; do
            for ML in 55000; do 
                for RN in 80; do
                    cd /your_path_to/${sample}_align_new/outs/systematic/
                    iaf=0.200
                    aaf=1.000          
                    ivc=20
                    avc=1000
                    ird=2
                    mkdir RD${RD}_SCO${min_score}_span${min_span}_ML${ML}_RN${RN}_${iaf}to${aaf}_${ivc}to${avc}_MSet2
                    cd RD${RD}_SCO${min_score}_span${min_span}_ML${ML}_RN${RN}_${iaf}to${aaf}_${ivc}to${avc}_MSet2
                    for acd in 5000; do
                        for icd in 600; do
                            for min_marker in 3; do
                                minMoleSize=1000 
                                outprefix=out_RD${RD}_ks_acd${acd}_icd${icd}_${min_marker}
                                DrLink recombis --split ../../splitted_vcf_mole/span${RD} --var ../../phased_variants.vcf.gz --mol ../../out_pool_${sample}_${RD}_min${minMoleSize}bp_molecule_table_trashme.txt.gz --str ../../smoothed_out_pool_${sample}_${RD}_min${minMoleSize}bp_readNum_stat.txt --stl ../../out_pool_${sample}_${RD}_min${minMoleSize}bp_moleLen_stat.txt --spa /your_path_to/marker_set1/parent1.txt --sma /your_path_to/marker_set2/parent2.txt --ims ${min_score} --imn ${min_marker} --imr ${min_span} --out ${outprefix} --ird ${ird} --ivc ${ivc} --avc ${avc} --iaf ${iaf} --aaf ${aaf} --nur ${RN} --nul ${ML} --acd ${acd} --icd ${icd} --aco 200000 > out_RD${RD}_ks_acd${acd}_icd${icd}_${min_marker}_check_DrLink.log
                            done
                        done
                    done
                done
            done
        done    
    done
done
#
################################ till now, you have one set of COs predicted based on parental genome as the reference ##############
#
#
# false positives can come from duplication/triplications in one parental genome, as compared the other. To remove such fps, you can repeat the above processes with the other parental genome as reference, leading to another set of CO prediction. By intersecting two sets of COs, false positives due to differences in repetitve regions between two parental genomes can be removed.
# intersect two CO sets based on two reference genomes as read alignment and variant calling
#
# sample
sample=sample_name
#
cpath=/your_path_to_parent1_genome_based_prediction/${sample}_align_new/outs/systematic                                    # a path to CO prediction based on reference genome 1
cco=${cpath}/out_RD25000_xx_acdxxx_icdxxx_minMarkerxxx_****_raw_BP_sorted.txt                                              # a file listing CO predictions 
cread=${cpath}/out_RD25000_mkrxx_acdxxx_icdxxx_mkrnumxxx_aMLxx_aRNxx_iSCxxx_iMKxxx_iSPxxx_iVCxx_aVCxxxx_iAFxx_aAFxx_irdxx/ # a folder lists read alignment info for each CO prediction
#
lpath=/your_path_to_parent2_genome_based_prediction/${sample}_align_new/outs/systematic                                    # a path to CO prediction based on reference genome 2
lco=${lpath}/out_RD25000_xx_acdxxx_icdxxx_minMarkerxxx_****_raw_BP_sorted.txt                                              # a file listing CO predictions 
lread=${lpath}/out_RD25000_mkrxx_acdxxx_icdxxx_mkrnumxxx_aMLxx_aRNxx_iSCxxx_iMKxxx_iSPxxx_iVCxx_aVCxxxx_iAFxx_aAFxx_irdxx/ # a folder lists read alignment info for each CO prediction
#
DrLink intersec ${cco} ${lco} ${cread} ${lread} outfolder_flagstr
#
#
#
#
#
#
#
#
#
#
#
