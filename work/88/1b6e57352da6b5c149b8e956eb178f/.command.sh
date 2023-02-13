#!/bin/bash -ue
sleep 6
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156796/miniml/GSE156796_family.xml.tgz
tar -zxvf GSE156796_family.xml.tgz
