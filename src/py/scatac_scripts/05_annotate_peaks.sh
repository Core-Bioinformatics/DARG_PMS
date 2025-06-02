#!/usr/bin/bash

peaks_folder=
input_file="${peaks_folder}/aggregate_peaks.bed"
output_file="${peaks_folder}/aggregate_annotated_peaks.txt"
homer_bin="../homer/bin"

annotatePeaks.pl ${input_file} hg38 -annStats ${output_file}
