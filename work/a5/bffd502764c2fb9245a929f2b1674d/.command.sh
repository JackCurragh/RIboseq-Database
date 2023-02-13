#!/bin/bash -ue
sleep 7
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157423/miniml/GSE157423_family.xml.tgz
tar -zxvf GSE157423_family.xml.tgz
