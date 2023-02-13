#!/bin/bash -ue
sleep 4
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157361/miniml/GSE157361_family.xml.tgz
tar -zxvf GSE157361_family.xml.tgz
