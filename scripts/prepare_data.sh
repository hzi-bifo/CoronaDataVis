#!/usr/bin/env bash


wget https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv -P ../data
wget https://raw.githubusercontent.com/nextstrain/ncov/master/data/metadata.tsv -P ../data


awk -F"\t" 'NR==1{print}$2=="ncov"&&$5~/20[12][90]-[0-1][0-9]-[0-9][0-9]/&&$13=="genome"&&$14>=28000&&$15=="Human"' ../data/metadata.tsv > ../data/metadata_fullgenome_filtered.tsv