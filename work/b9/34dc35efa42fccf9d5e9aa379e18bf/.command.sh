#!/bin/bash -ue
sleep 6
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180669/miniml/GSE180669_family.xml.tgz
tar -zxvf GSE180669_family.xml.tgz
