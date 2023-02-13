#!/bin/bash -ue
sleep 9
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158141/miniml/GSE158141_family.xml.tgz
tar -zxvf GSE158141_family.xml.tgz
