# Intro
This repository contains the source codes used in the `Experimental Results` section of the research paper [Restricted Hamming-Huffman trees](https://www.rairo-ro.org/component/article?access=doi&doi=10.1051/ro/2022066)
published in [RAIRO - Operations Research](https://www.rairo-ro.org/).


The relation among each source code and experiment is as follows: 

1.`experiment_1.cpp`: evaluation of the costs described  in section `5.1. uniform [3]-HHT optimality hypothesis`.

2.`experiment_2.cpp`: evaluation of the costs described in section `5.2. [2]-HHTs efficiency`.

3.`experiment_3.cpp`: evaluation of the error detection capabilities described in section `5.2. [2]-HHTs efficiency`.

# Running

To run the programs in a Linux environment use:

`g++ {file_name} -o {desired_exe_name} && ./{desired_exe_name}`

For instance, to run `experiment_1.cpp`, one may use:

`g++ {experiment_1.cpp} -o {exp1} && ./{exp1}`
