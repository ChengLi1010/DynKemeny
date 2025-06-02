#!/bin/sh
g_arr=(powergrid) #powergrid astroph hepth emailenron amazon dblp youtube roadnetPA roadnetCA orkut)
alg_arr=(ttf spantree forestmc) 

make
# make preprocess

for g in "${g_arr[@]}"; do
    # for alg in "${alg_arr[@]}"; do
    #     ./main -d ${g} -a ${alg} -dynamic -log datasets/${g}/res_dynamic.txt
    # done
    ./main -d ${g} -a ttf -dynamic -basic -log datasets/${g}/res_dynamic.txt
    ./main -d ${g} -a ttf -dynamic -imp -log datasets/${g}/res_dynamic.txt
done