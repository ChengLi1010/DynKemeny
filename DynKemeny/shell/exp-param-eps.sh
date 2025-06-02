#!/bin/sh
# g_arr=(powergrid astroph hepth emailenron amazon dblp youtube roadnetPA roadnetCA orkut)
g_arr=(powergrid)
eps_arr=(0.5 0.4 0.3 0.2 0.1)
alg_arr=(ttf spantree forestmc) 

# make -B
make test

for g in "${g_arr[@]}"; do
    for alg in "${alg_arr[@]}"; do
        for eps in "${eps_arr[@]}"; do
            ./main -d ${g} -a ${alg} -eps ${eps} -static -log datasets/${g}/eps_${alg}.txt
        done
    done
done


