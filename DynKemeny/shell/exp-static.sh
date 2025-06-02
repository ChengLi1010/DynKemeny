#!/bin/sh
g_arr=(powergrid hepth) #powergrid astroph hepth emailenron amazon dblp youtube roadnetPA roadnetCA orkut)
alg_arr=(ttf spantree forestmc) 
make

for g in "${g_arr[@]}"; do
    for alg in "${alg_arr[@]}"; do
        ./main -d ${g} -a ${alg} -static -log datasets/${g}/res_static.txt
    done
done


