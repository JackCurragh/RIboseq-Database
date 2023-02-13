#!/bin/bash -ue
sleep 5
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173856/miniml/GSE173856_family.xml.tgz
tar -zxvf GSE173856_family.xml.tgz
