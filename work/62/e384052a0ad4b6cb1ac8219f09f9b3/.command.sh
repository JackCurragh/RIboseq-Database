#!/bin/bash -ue
sleep 5
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185458/miniml/GSE185458_family.xml.tgz
tar -zxvf GSE185458_family.xml.tgz
