#!/bin/bash -ue
sleep 9
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136940/miniml/GSE136940_family.xml.tgz
tar -zxvf GSE136940_family.xml.tgz
