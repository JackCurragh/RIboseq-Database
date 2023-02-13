#!/bin/bash -ue
sleep 4
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185286/miniml/GSE185286_family.xml.tgz
tar -zxvf GSE185286_family.xml.tgz
