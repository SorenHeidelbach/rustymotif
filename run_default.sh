#!/bin/bash

REFERENCE="/home/bio.aau.dk/lx38ll/dark-science/motif-identification/analysis/monocultures/m_ruber/medaka/assembly.polished.fasta"
PILEUP="/home/bio.aau.dk/lx38ll/dark-science/motif-identification/data/local/m_ruber/modkit.pileup.bed"

cargo run --release -- $REFERENCE  $PILEUP





###

###