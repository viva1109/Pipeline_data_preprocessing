#!/bin/bash               

arg_s=1
arg_e=3
arg_ee=3
arg_C=1
Set_Data="/data/sharedcode/kjkim/RFunctions/get_Set.R"
Get_Pvalue="/data/sharedcode/kjkim/RFunctions/Get_pval.R"
Get_Plots="/data/sharedcode/kjkim/RFunctions/Get_plots.R"
arg_r=5
arg_M="TMAT15_OMiAT_wilcoxon_ANCOM"
arg_n="Control_Case"
Tree_EZ="/data/sharedcode/kjkim/ez/ez_sina_fastmp.tre"
Tree_Silva="/data/sharedcode/kjkim/silva128/97_otus.tre"
arg_q=0.05
setDirectories(){
    INPUT_DIR=$(readlink -f $INPUT_DIR)
    OUTPUT_DIR=$(readlink -f $OUTPUT_DIR)
}

help() {
    echo "process [OPTIONS] FILE"
    echo "    -h         Help"
    echo "    -d ARG     Input dirname (S1)"
    echo "    -D ARG     Output dirname(S1-3)"
    echo "    -t ARG     Threads (S2)"
    echo "    -z ARG     Database choice: ez, silva (S2-3, for S3 you can put it like this: ez_silva)" 
    echo "    -s ARG     Start step"
    echo "    -e ARG     End step"
    echo "    -M ARG     Methods choice. Default: TMAT15_OMiAT_wilcoxon_ANCOM (S1-3, for S3, you need to put three or four methods)"
    echo "    -p ARG     Number of Permutation (S2)"
    echo "    -c ARG     Number of Columns of Plots (S2-3)"
    echo "    -P ARG     Phenotype file name. Just a name. Not a path. (S1)"
    echo "    -C ARG     Is Phenotype Continuous or Ordinal? If yes, 1. If not, 0. (S2)"
    echo "    -r ARG     Toxonomy rank. 1: phylum 2: Class 3: Order 4: Family 5: Genus 6: Species"
    echo "    -n ARG     Names of phenotype from 0 to k k=Number of levels of Phenotype
                         ex) -n Control_Case => 0:Control, 1:Case (S3)"
    echo "    -q ARG     FDR-q value for ANCOM"
    echo ""
    echo "    #Step description: "
    echo "    1-Set data 2-Get p-value 3-Get plots"
    echo ""
        exit 0
}

while getopts "d:D:t:z:s:e:M:p:c:P:C:n:q:r:" opt
do
    case $opt in
        d) arg_d=$OPTARG
            INPUT_DIR=$arg_d
            ;;
	r) arg_r=$OPTARG
	    ;;
        D) arg_D=$OPTARG
	    OUTPUT_DIR=$arg_D
	    ;;
        t) arg_t=$OPTARG
            THREADS=$arg_t
            ;;
        z) arg_z=$OPTARG
            METHOD=$arg_z
            ;;
        s) arg_s=$OPTARG
           ;;
        e) arg_e=$OPTARG
           ;;
	M) arg_M=$OPTARG
           ;;
	p) arg_p=$OPTARG
           ;;
        c) arg_c=$OPTARG
           ;;
	P) arg_P=$OPTARG
	   ;;
	C) arg_C=$OPTARG
	   ;;
	n) arg_n=$OPTARG
	   ;;
	q) arg_q=$OPTARG
	   ;;
        h) help ;;
        ?) help ;;
    esac
done


run1(){
        Rscript $Set_Data $OUTPUT_DIR $INPUT_DIR $arg_M $Tree_EZ $Tree_Silva $arg_P $arg_z $arg_r
}
run2(){
        Rscript $Get_Pvalue $OUTPUT_DIR $arg_z $arg_p $arg_M $arg_c $THREADS $arg_C $arg_q
}
run3(){
        Rscript $Get_Plots $OUTPUT_DIR $arg_n $arg_c $arg_z $arg_M
}

run(){
	bash /data/sharedcode/kjkim/bash/Updata_RLibraries.sh	
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

mkdir -p $OUTPUT_DIR
setDirectories
run


