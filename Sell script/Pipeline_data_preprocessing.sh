#!/bin/bash                                                                       

INPUT_DIR=""
TRIM_DIR=""
MERGE_DIR=""
FILTER_DIR=""
STAT_DIR=""
CHIM_DIR=""
TAX_DIR=""

MINLEN=350
MAXLEN=550
arg_s=1
arg_e=12
arg_ee=12
THREADS=4
CHIMERA_REF="/data/sharedcode/scpark/default/gold.fasta"
PERL=$(which perl)
VSEARCH=$(which vsearch)
CUTADAPT=$(which cutadapt)
CASPER="casper"

MATCH_OTU="/data/sharedcode/scpark/python/match_otu_to_tax.py"
PERL_SCRIPT="/data/sharedcode/scpark/perl/map.pl"
BOKULICH="/data/sharedcode/scpark/python/bokulich.py"
PARALLEL_R="/data/sharedcode/kjkim/R/go.R"

PRIMER_FOR=""
PRIMER_REV=""

FILE_IDENTIFER_F="_1"
FILE_IDENTIFER_R="_2"

METHOD="ez"


DB_EZ="/data/MICROBIOME/MDHC/eztaxon/eztaxon_qiime_full.fasta"
TAX_ID_EZ="/data/MICROBIOME/MDHC/eztaxon/eztaxon_id_taxonomy.txt"
TAX_ID_EZ_CODE="/data/sharedcode/kjkim/taxonomy_code/eztaxon_id_taxonomy_togo.txt"

DB_SILVA="/data/MICROBIOME/MDHC/SILVA_128_QIIME_release/rep_set/rep_set_all/97/97_otus.fasta"
TAX_ID_SILVA="/data/MICROBIOME/MDHC/SILVA_128_QIIME_release/taxonomy/taxonomy_all/97/taxonomy_7_levels.txt"
TAX_ID_SILVA_CODE="/data/sharedcode/kjkim/taxonomy_code/taxonomy_7_levels_code_togo.txt"

DB_GG="/data/software/miniconda2/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta"
TAX_ID_GG="/data/software/miniconda2/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt"


setDirectories(){
    INPUT_DIR=$(readlink -f $INPUT_DIR)
    TRIM_DIR=$(dirname $INPUT_DIR)/"trimed" 
    MERGE_DIR=$(dirname $INPUT_DIR)/"merged"    
    FILTER_DIR=$(dirname $INPUT_DIR)/"filter"   
    STAT_DIR=$(dirname $INPUT_DIR)/"stat"   
    DEREP_DIR=$(dirname $INPUT_DIR)/"derep"  
    CHIM_DIR=$(dirname $INPUT_DIR)/"chimera"    
    OTU_DIR=$(dirname $INPUT_DIR)/"otu"
    if [ "$METHOD" == "ez" ]; then
        TAX_DIR=$OTU_DIR/"tax_ez"
        TAX_ID_GO=$TAX_ID_EZ
    elif [ "$METHOD" == "silva" ];  then
        TAX_DIR=$OTU_DIR/"tax_silva"
	TAX_ID_GO=$TAX_ID_SILVA
    fi  
    echo "MIN MAXLEN : $MINLEN - $MAXLEN"
    echo "TAX_DIR : $TAX_DIR"
}

unzip(){
    echo Unzip Fastq Files for primer remove
    echo Decompressing...
    gzip -d $INPUT_DIR/*.gz
    echo
}

unzip2(){
    echo Unzip Fastq Files for primer remove
    echo Decompressing...
    cd $INPUT_DIR
    rm $INPUT_DIR/".."/"lala.txt"
    for file in $(find ./ -maxdepth 1 -type f)
    do
            echo 'gzip -d '$INPUT_DIR/$file >> $INPUT_DIR/".."/"lala.txt"
    done
    Rscript $PARALLEL_R $THREADS $INPUT_DIR/".."/
    rm $INPUT_DIR/".."/"lala.txt"
    echo
}

obtainChimeraRef(){
    echo Obtaining Gold reference database for chimera detection

    if [ ! -e $CHIMERA_REF ]; then
        cd $(dirname $CHIMERA_REF)
        if [ ! -e Silva.gold.bacteria.zip ]; then
            wget https://www.mothur.org/w/images/f/f1/Silva.gold.bacteria.zip
    fi

    echo Decompressing and reformatting...

    echo "unzip -p Silva.gold.bacteria.zip silva.gold.align | sed -e "s/[.-]//g" > $(basename $CHIMERA_REF)"
    unzip -p Silva.gold.bacteria.zip silva.gold.align | sed -e "s/[.-]//g" > $(basename $CHIMERA_REF)
    fi
}

checkFQformat(){
    cd $INPUT_DIR
    echo Checking FASTQ format version for one file
    $VSEARCH --threads $THREADS --fastq_chars $(ls -1 *.fastq | head -1)
}

checkPair(){
    cd $INPUT_DIR
    for file in $(find ./ -maxdepth 1 -type f)
    do
        if [[ "$file" == *$FILE_IDENTIFER_F* ]]; then
            local str=$(basename $file)
            local TMP_ID=${str//$FILE_IDENTIFER_F*/""}
            local FILE_EXT=${str//*$FILE_IDENTIFER_F/""}
            local REV_FILE=./$TMP_ID$FILE_IDENTIFER_R$FILE_EXT

          if [ ! -e $REV_FILE ]; then
            echo "$REV_FILE not found."         
        fi
       fi
    done
    echo
}

trimAdapter(){

    echo Trim adapter with cutadapt

    cd $INPUT_DIR
    local OUT_DIR=$TRIM_DIR
    mkdir -p $OUT_DIR

    for file in $(find ./ -maxdepth 1 -type f)
    do
	echo $file
	echo $arg_T	
        local OUTFILE=$OUT_DIR/$file
	if [ "$arg_T" == "merged" ]; then
            $CUTADAPT -a $PRIMER_FOR...$PRIMER_REV -o $OUTFILE $file -O 11 -e 0.15 -m 10
        else
		if [[ $file == *$FILE_IDENTIFER_F* ]]; then

            	local str=$(basename $file)
            	local TMP_ID=${str//$FILE_IDENTIFER_F*/""}
            	local FILE_EXT=${str//*$FILE_IDENTIFER_F/""}
            	local REV_FILE=./$TMP_ID$FILE_IDENTIFER_R$FILE_EXT

	        local OUTFILE_FWD=$OUT_DIR/$TMP_ID$FILE_IDENTIFER_F$FILE_EXT
                local OUTFILE_REV=$OUT_DIR/$TMP_ID$FILE_IDENTIFER_R$FILE_EXT
	        fi
        	    $CUTADAPT -g $PRIMER_FOR -G $PRIMER_REV -o $OUTFILE_FWD -p $OUTFILE_REV $file $REV_FILE -O 11 -e 0.15 -m 10
	fi
    done

    echo
}

merge(){
    echo Merge sequence

    cd $TRIM_DIR
    local OUT_DIR=$MERGE_DIR
    mkdir -p $OUT_DIR
    rm $TRIM_DIR/"lala.txt"
    for file in $(find ./ -maxdepth 1 -type f)
    do
        if [[ $file == *$FILE_IDENTIFER_F* ]]; then
            local str=$(basename $file)
            local TMP_ID=${str//$FILE_IDENTIFER_F*/""}
            local FILE_EXT=${str//*$FILE_IDENTIFER_F/""}
            local REV_FILE=./$TMP_ID$FILE_IDENTIFER_R$FILE_EXT

            local OUT_FILE=$OUT_DIR/$TMP_ID
            echo $CASPER $file $REV_FILE' -t '$THREADS' -g 0.27 -o '$OUT_FILE >> $TRIM_DIR/"lala.txt"
        fi
    done
    echo
    Rscript $PARALLEL_R $THREADS $TRIM_DIR
    rm $TRIM_DIR/"lala.txt"

}

filter(){
    echo Filtering sequence by BOKULICH

    cd $MERGE_DIR
    local OUT_DIR=$FILTER_DIR
    mkdir -p $OUT_DIR

    source deactivate qiime1
    rm $MERGE_DIR/"lala.txt"
    for file in $(find ./ -maxdepth 1 -type f)
    do
        local str=$(basename $file)
        local OUT_FILE=$OUT_DIR/$str
        #python $BOKULICH -i $file -o $OUT_FILE -m $MINLEN -M $MAXLEN
	echo 'python '$BOKULICH' -i '$file' -o '$OUT_FILE' -m '$MINLEN' -M '$MAXLEN >> $MERGE_DIR/"lala.txt"
    done
    echo
    Rscript $PARALLEL_R $THREADS $MERGE_DIR
    rm $MERGE_DIR/"lala.txt" 
    source activate qiime1
}

dereplication(){
    echo ====================================
    echo Processing sample $s
    cd $FILTER_DIR
    local OUT_DIR=$DEREP_DIR
    mkdir -p $OUT_DIR
    for file in $(find ./ -maxdepth 1 -type f)
    do
        local TMP_ID=$(basename "$file" | cut -d. -f1)
        echo $TMP_ID 
        $VSEARCH --threads $THREADS \
        --derep_fulllength $file \
        --strand plus \
        --output $OUT_DIR/$TMP_ID.derep.fasta  \
        --sizeout \
        --uc $OUT_DIR/$TMP_ID.derep.uc \
        --relabel $TMP_ID. \
        --fasta_width 0
    done
    echo
    echo Sum of unique sequences in each sample: $(cat  $OUT_DIR/*.derep.fasta | grep -c "^>")

#   rm -f $OUT_DIR/all.derep.fasta $OUT_DIR/all.nonchimeras.derep.fasta
#   cat $OUT_DIR/*.derep.fasta > $OUT_DIR/all.fasta
}

removeChiRef(){

    echo Reference chimera detection

    #obtainChimeraRef

    cd $DEREP_DIR
    local OUT_DIR=$CHIM_DIR
    mkdir -p $OUT_DIR

    for file in $(find ./ -maxdepth 1 -type f)
    do
        if [[ $file == *"fasta"* ]]; then
            
            local TMP_ID=$(basename "$file" | cut -d. -f1)
            echo $file
            $VSEARCH --threads $THREADS \
                --uchime_ref $file \
                --db $CHIMERA_REF \
                --sizein \
                --sizeout \
                --fasta_width 0 \
                --nonchimeras $OUT_DIR/$TMP_ID.derep.chim.fasta
            fi
    done

    rm -f all.nonchimeras.derep.fasta
    cat $OUT_DIR/*.derep.chim.fasta > $OUT_DIR/all.nonchimeras.derep.fasta
    echo
}

cluster(){
    cd $CHIM_DIR
    local OUT_DIR=$OTU_DIR
    mkdir -p $OUT_DIR
    # echo "$VSEARCH --threads $THREADS \
    #     --cluster_size all.nonchimeras.derep.fasta \
    #     --id 0.97 \
    #     --strand plus \
    #     --sizein \
    #     --sizeout \
    #     --fasta_width 0 \
    #     --uc $OUT_DIR/all.clustered.uc \
    #     --relabel OTU_ \
    #     --centroids $OUT_DIR/all.otu.centroids.fasta \
    #     --otutabout $OUT_DIR/all.otutab.txt"
    $VSEARCH --threads $THREADS \
        --cluster_size all.nonchimeras.derep.fasta \
        --id 0.97 \
        --strand plus \
        --sizein \
        --fasta_width 0 \
        --uc $OUT_DIR/all.clustered.uc \
        --relabel OTU_ \
        --centroids $OUT_DIR/all.otu.centroids.fasta \
        --otutabout $OUT_DIR/all.otutab.txt 
}

assign_tax(){
    local DB=""
    local TAX=""

    if [ "$METHOD" == "ez" ]; then
        echo database is ez
        DB=$DB_EZ
        TAX=$TAX_ID_EZ_CODE
        TAX_DIR=$OTU_DIR/"tax_ez"
    elif [ "$METHOD" == "silva" ];  then
        echo database is silva
        DB=$DB_SILVA
        TAX=$TAX_ID_SILVA_CODE
        TAX_DIR=$OTU_DIR/"tax_silva"
    fi
    source activate qiime1
    /usr/bin/time -p -f "%E %K" parallel_assign_taxonomy_uclust.py \
                                -i $OTU_DIR/all.otu.centroids.fasta \
                                -o $TAX_DIR \
                                --reference_seqs_fp $DB \
                                --id_to_taxonomy_fp $TAX \
                                -O $THREADS \
                                --uclust_max_accepts 1 \
                                --min_consensus_fraction 1
}


post_otu(){
   
# cd $OTU_DIR/..
 #   python /data/sharedcode/kjkim/python/taxonomy_assign2.py -d $METHOD    
    
    #local FILE_ID="all.otu.centroids"
    #source activate qiime1    
    #echo "make otu able"
    #/usr/bin/time -p -f "%E" $MATCH_OTU -t $TAX_DIR/$FILE_ID"_tax_assignments.txt" -m all.otutab.txt  -o otu_comp.txt
    # /usr/bin/time -p -f "%E" make_otu_table.py -i all.otutab.txt  -o otu_table_no_tax.biom 
    # /usr/bin/time -p -f "%E" filter_otus_from_otu_table.py -i otu_table_no_tax.biom -o otu_table_mc2_no_tax.biom -n 2

    # echo "biom add taxonomy"
    # /usr/bin/time -p -f "%E" python $MATCH_OTU -t all.otutab.txt -m TAX_DIR/"$FILE_ID"_tax_assignments.txt -o otu_table_comp_"$METHOD"_tax.txt


    #echo "parallel align seqs_pynast"
    #/usr/bin/time -p -f "%E" parallel_align_seqs_pynast.py -i all.otu.centroids.fasta -o pynast_aligned_seqs -T --jobs_to_start $THREADS

    #echo "make phylogeny"
    #qiime2
    #export OMP_NUM_THREADS=$THREADS
    #FastTreeMP -fasetest -nt pynast_aligned_seqs/"$FILE_ID"_aligned.fasta > rep_set.tre    
    cd $OTU_DIR
    local FILE_ID="all.otu.centroids"

    echo "make otu able"
    /usr/bin/time -p -f "%E" biom convert -i all.otutab.txt -o otu_table.biom --table-type="OTU table" --to-hdf5
    /usr/bin/time -p -f "%E" filter_otus_from_otu_table.py -i otu_table.biom -o otu_table_mc2.biom -n 2
    /usr/bin/time -p -f "%E" biom convert -i otu_table_mc2.biom -o all.otutab_mc2.txt --to-tsv
    /usr/bin/time -p -f "%E" python $MATCH_OTU -t $TAX_DIR/"$FILE_ID"_tax_assignments.txt -m all.otutab_mc2.txt  -o otu_comp_$METHOD\.txt

}

post2(){
	cd $OTU_DIR
	Rscript /data/sharedcode/kjkim/R/taxonomy_attach_DB.R ./otu_comp_$METHOD\.txt $METHOD $TAX_ID_GO $THREADS 
}

run1(){
        unzip2
	#mkdir -p /tmp/jobs
    	#chmod 777 -R /tmp/jobs
}
run2(){
        checkFQformat
}
run3(){
        checkPair
}
run4(){
        trimAdapter
}   
run5(){
        merge
}
run6(){
        filter
}
run7(){
        dereplication
}
run8(){
        removeChiRef
}
run9(){
        cluster
}
run10(){
        assign_tax
}
run11(){
        post_otu
}
run12(){
	post2
}

run(){
    #unzip2
    #checkFQformat
    #checkPair
    #trimAdapter
    #merge
    #filter
    #dereplication
    #removeChiRef
    #cluster
    #assign_tax
    #post_otu
    mkdir -p /tmp/jobs
    chmod 777 -R /tmp/jobs
    i=0

    while [ $i -le $arg_ee ]

    do
        i=$(($i+1))
        if [ $i -ge $arg_s ] && [ $i -le $arg_e ]; then
            echo $i
            run$i
        fi
    done

}

help() {
    echo "process [OPTIONS] FILE"
    echo "    -h         help"
    echo "    -d ARG     input dirname"
    echo "    -f ARG     forward primer"
    echo "    -r ARG     reverse primer"
    echo "    -F ARG     identifier forward seq"
    echo "    -R ARG     identifier reverse seq"
    echo "    -t ARG     threads"
    echo "    -z ARG     ez, silva" 
    echo "    -m ARG     min_len : default 350"
    echo "    -M ARG     max_len : default 550"
    echo "    -s ARG     start step"
    echo "    -e ARG     end step"
    echo "    -T ARG	 sequence type"
    echo ""
    echo "    #Step description: "
    echo "    1-unzip 2-checkFQformat 3-checkPair"
    echo "    4-trimAdapter 5-merge 6-filter"
    echo "    7-dereplication 8-removeChiRef 9-cluster"
    echo "    10-assign_tax 11-post_otu(attach unique otu code) 12-post2(attach taxonomy)"
    echo ""
        exit 0
}


while getopts "d:f:r:F:R:t:z:m:M:s:e:T:" opt
do
    case $opt in
        d) arg_d=$OPTARG
            INPUT_DIR=$arg_d
            ;;
        f) arg_f=$OPTARG
            PRIMER_FOR=$arg_f
            ;;
        r) arg_r=$OPTARG
            PRIMER_REV=$arg_r
            ;;
        F) arg_F=$OPTARG
            FILE_IDENTIFER_F=$arg_F
            ;;
        R) arg_R=$OPTARG
            FILE_IDENTIFER_R=$arg_R
            ;;    
        t) arg_t=$OPTARG
            THREADS=$arg_t
            ;;
        z) arg_z=$OPTARG
            METHOD=$arg_z
            ;;
        m) arg_m=$OPTARG
            MINLEN=$arg_m
            ;;
        M) arg_M=$OPTARG
            MAXLEN=$arg_M
            ;;
        s) arg_s=$OPTARG
           ;;
        e) arg_e=$OPTARG
           ;;
	T) arg_T=$OPTARG
           ;;

        h) help ;;
        ?) help ;;
    esac
done

echo $arg_z
if [ -z "$arg_d" ]  || [ -z "$arg_z" ]; then
    help
fi

#source activate qiime1

setDirectories
run
 
# # getopt 부분 끝나고 난 후의 인자(FILE) 읽기
# shift $(( $OPTIND - 1))
# file=$1
# echo "$file"



