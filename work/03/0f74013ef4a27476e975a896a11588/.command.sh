#!/bin/bash -ue
sleep 6
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130465/miniml/GSE130465_family.xml.tgz
tar -zxvf GSE130465_family.xml.tgz
