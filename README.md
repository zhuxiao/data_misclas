# misClus Experiments

## Prerequisites
#### Installation

```bash
# Debian / Ubuntu:
$ sudo apt update  # Ensure the package list is up to date
$ sudo apt install mlpack-bin g++ python3 autoconf automake libtool samtools openjdk-11-jdk
```

```bash
# RedHat / CentOS:
$ sudo yum install gcc-c++ python3 autoconf automake libtool samtools java-11-openjdk-devel
```

```bash
# htslib
$ wget -c https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
$ tar -xvjf htslib-1.19.tar.bz2
$ cd htslib-1.19
$ autoreconf -i
$ ./configure
$ make 
$ make install

# bwa
$ wget -c https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
$ tar -xvjf bwa-0.7.17.tar.bz2
$ cd bwa
$ make
```

 htslib (1.9 version or later)  and  bwa be installed from source files (it cannot be installed by `apt install` or `yum install`)

```bash
# misasm
$ git clone https://github.com/zhuxiao/misasm.git
$ cd misasm/src/
$ make
```

And the binary file `misasm` will be output into the folder bin in this package directory.

```bash
# misfinder 
$ git clone https://github.com/zhuxiao/misFinder.git
$ cd misFinder
$ ./autogen.sh
```

And the binary file `mf` will be output into the folder bin in this package directory.

```bash
# Pilon 
$ wget -c https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
```

And the jar file `pilon` will be download in this package directory.

```bash
#Quast 
$ wget -c https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
$ tar -zvxf quast-5.2.0.tar.gz
$ cd quast-5.2.0
```

And the python file `quast.py` will be put in this package directory.

## Experiments

####  Schizosaccharomyces pombe

```bash
# Download  
$ prefetch SRR16381171
# Assemble paired-end reads to contigs
$ fasterq-dump --split-files SRR16381171.sra
$ fq2fa --merge SRR16381171_1.fq SRR16381171_2.fq SRR16381171.fa
$ idba_ud --mink 20 --maxk 120 --min_contig 1000 --num_threads 8 -r SRR16381171.fa -o S.pombe_out
# Align paired-end reads onto contigs
$ cd S.pombe_out
$ bwa index scaffold.fa
$ bwa mem -t 8 scaffold.fa SRR16381171_1.fastq SRR16381171_2.fastq > S.pombe.sam
$ samtools view -b -S -@ 8 S.pombe.sam > S.pombe.bam
$ samtools sort -@ 8 S.pombe.bam -o S.pombe_sorted.bam
$ samtools index -@ 8 S.pombe_sorted.bam S.pombe_sorted.bai
```

We used [misasm](https://github.com/zhuxiao/misasm) to detect misassemblies based on paired-end reads and used `misclus` to  analyze assembly errors in the misassembly regions.

```bash
# misasm candidate misassemblies
$ misasm scaffold.fa S.pombe_sorted.bam
$ cd output
# Classify misassemblies
$ cat genome_Indel genome_Misjoin | awk '{split($1,a,":"); split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > MisasmRaw
$ generateRandReg.py MisasmRaw MisasmReg
#Misassembly clustering
$ misclus MisasmReg scaffold.fa S.pombe_sorted.bam
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.

We used [misfinder](https://github.com/zhuxiao/misFinder) to identify mis-assemblies in an unbiased manner using reference and paired-end reads and used `misclus` to  analyze assembly errors in the misassembly regions.

```bash
# misfinder candidate misassemblies  
$ mf all -t 8 -o output config 
$ cd output
$ cat misassemblies_misFinder | awk '{print $1"\t"$3"\t"$4}' > misFindeRaw
$ generateRandReg.py MisasmRaw MisasmReg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa S.pombe_sorted.bam
$ cat result_errors | awk '{print $1"\t"$3"\t"$4}' > misFindeRaw
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.

We used [Pilon](https://github.com/broadinstitute/pilon) to read alignment analysis to identify inconsistencies between the input genome and the evidence in the reads.We then clustered the misassembly regions using misclus  to determine if the intervals included assembly errors.

```bash
# Pilon candidate misassemblies
$ java -jar pilon --genome scaffold.fa --frags S.pombe_sorted.bam --changes > misassemblies_Pilon
$ cat misassemblies_Pilon | grep "fix" | grep ":" | cut -d ":" -f 2,3 | awk '{split($1,a,":");split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > PilonRaw
$ generateRandReg.py MisasmRaw MisasmReg
# Misassembly clustering
$ misclus PilonReg scaffold.fa S.pombe_sorted.bam
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.

[Quast](https://github.com/ablab/quast) is used to evaluates genome/metagenome assemblies by computing various metrics.We then clustered the misassembly regions using misclus  to determine if the intervals included assembly errors.

```bash
# Quast candidate misassemblies
$ quast.py -r GCF_000002945.1_ASM294v2_genomic.fa --extensive-mis-size 200 scaffold.fa -m 1000
$ cd ./quast_result/latest
$ extractQUASTreg.py misassemblies_QUAST > QUASTraw
$ generateRandReg.py MisasmRaw MisasmReg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa S.pombe_sorted.bam
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.

#### Simulate Schizosaccharomyces pombe data

To benchmark the evaluation accuracy of misclus, we used a simulated Schizosaccharomyces pombe assembly containing both structural and small-scale assembly errors.

```bash
$ /home/jlgao/hxliu/tools/mason2-2.0.9-Linux-x86_64/bin/mason_variator 
        -ir GCF_000002945.1_ASM294v2_genomic.fa 
        -ov simulate_Spombe.vcf \
        --snp-rate 0 \
        --small-indel-rate 0.00001 \
        --min-small-indel-size 20 \
        --max-small-indel-size 50 \
        --sv-indel-rate 0.0000001 \
        --sv-inversion-rate 0 \
        --sv-translocation-rate 0.000001 \
        --sv-duplication-rate 0
```

We used Mason2 to generate the answer VCF file based on the Schizosaccharomyces pombe reference data and used this answer VCF file to create simulated Schizosaccharomyces pombe assembly data that includes structural and small-scale assembly errors.

```bash
$ mason_simulator -ir GCF_000002945.1_ASM294v2_genomic.fa \
        -n 2490000 \
        -iv simulate_Spombe.vcf \
        -o left_reads.fq \
        -or right_reads.fq \
        --fragment-size-model normal \
        --fragment-mean-size 450 \  
        --fragment-size-std-dev 50 \ 
        --seq-technology illumina \
        --seq-mate-orientation FR \
        --illumina-read-length 150 \  
```

We align paired-end simulated reads onto contigs ,the assembled contigs were simulated by the simulated human genome .

```bash
$ bwa mem -t 8 GCF_000002945.1_ASM294v2_genomic.fa left_reads.fq right_reads.fq > simulated_Spombe.sam
$ samtools view -b -S -@ 8 simulated_Spombe.sam > simulated_Spombe.bam
$ samtools sort simulated_Spombe.bam -o simulated_Spombe_sorted.bam
$ samtools index -@ 8 simulated_Spombe_sorted.bam simulated_Spombe_sorted.bai
```

We randomly select a subset of loci from the VCF file to create a region file for evaluating misclus.

```bash
#misassemble clustering
$ misclus simulateReg ref_Spombe.fa simulate_Spombe_sorted.bam
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.

####  HG002 

```bash
# Download  HG002-NYGC-NovaSeq-2x250.
$ prefetch SRR11321732
# Assemble paired-end reads to contigs
$ fasterq-dump --split-files SRR11321732.sra
$ masurca -t 32 -i SRR11321732_1.fastq,SRR11321732_2.fastq -o SRR11321732_masurca
# Align paired-end reads onto contigs
$ cd SRR11321732_masurca
$ bwa index scaffold.fa
$ bwa mem -t 8 scaffold.fa SRR11321732_1.fastq SRR11321732_2.fastq > masurca.hg002.sam
$ samtools view -b -S -@ 8 hg002.sam > hg002.bam
$ samtools sort -@ 8 hg002.bam -o hg002_sorted.bam
$ samtools index -@ 8 hg002_sorted.bam hg002_sorted.bai
```

We used [misasm](https://github.com/zhuxiao/misasm) to detect misassemblies based on paired-end reads and used `misclus` to  analyze assembly errors in the misassembly regions.

```bash
# misasm candidate misassemblies
$ misasm scaffold.fa hg002_sorted.bam
# Classify misassemblies
$ cat genome_Indel genome_Misjoin | awk '{split($1,a,":"); split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > MisasmRaw
$ filterOverlapReg.py MisasmRaw scaffold.fa.fai MisasmReg
# Misassembly clustering
$ misclus MisasmReg scaffold.fa $ misclus MisasmReg scaffold.fa hg002_sorted.bam 
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.

We used [Pilon](https://github.com/broadinstitute/pilon) to read alignment analysis to identify inconsistencies between the input genome and the evidence in the reads.We then clustered the misassembly regions using misclus  to determine if the intervals included assembly errors.

```bash
# Pilon candidate misassemblies ？
$ java -jar pilon --genome scaffold.fa --frags S.pombe_sorted.bam --changes > misassemblies_Pilon
$ cat misassemblies_Pilon | grep "fix" | grep ":" | cut -d ":" -f 2,3 | awk '{split($1,a,":");split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > PilonRaw
$ filterOverlapReg.py PilonRaw scaffold.fa.fai PilonReg
# Misassembly clustering
$ misclus PilonReg scaffold.fa hg002_sorted.bam
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.

[Quast](https://github.com/ablab/quast) is used to evaluates genome/metagenome assemblies by computing various metrics.We then clustered the misassembly regions using misclus  to determine if the intervals included assembly errors.

```bash
# Quast candidate misassemblies
$ quast.py -r hs38.fa --extensive-mis-size 200 scaffold.fa -m 1000
$ cd ./quast_result/latest
$ extractQUASTreg.py misassemblies_QUAST > QUASTraw
$ filterOverlapReg.py QUASTraw scaffold.fa.fai QUASTreg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa hg002_sorted.bam
```

#### Simulate GRCh38 data

To benchmark the evaluation accuracy of misclus, we used a simulated human whole-genome assembly containing both structural and small-scale assembly errors.

```bash
$ /home/jlgao/hxliu/tools/mason2-2.0.9-Linux-x86_64/bin/mason_variator \
        -ir hg38.fa \
        -ov simulate_hg38.vcf \
        --snp-rate 0 \
        --small-indel-rate 0.000005 \
        --min-small-indel-size 20 \
        --max-small-indel-size 50 \
        --sv-indel-rate 0.0000001 \
        --sv-inversion-rate 0 \
        --sv-translocation-rate 0.000001 \
        --sv-duplication-rate 0
```

We used Mason2 to generate the answer VCF file based on the Homosapiens chromosome 14 reference data and used this answer VCF file to create simulated Homosapiens chromosome assembly data that includes structural and small-scale assembly errors.

```bash
$ mason_simulator --num-threads 16 -ir hg38.fa \
                -n 320000000 \
                -iv simulate_hg38.vcf \
                -o left_reads_2.fq \
                -or right_reads_2.fq \
                --fragment-size-model normal \
                --fragment-mean-size 450 \
                --fragment-size-std-dev 50 \
                --seq-technology illumina \
                --seq-mate-orientation FR \
                --illumina-read-length 150 \	
```

We align paired-end simulated reads onto contigs ,the assembled contigs were simulated by the simulated human genome.

```bash
$ bwa mem -t 8 hg38.fa left_reads.fq right_reads.fq > simulated_hg38.sam
$ samtools view -b -S -@ 8 simulated_hg38.sam > simulated_hg38.bam
$ samtools sort simulated_hg38.bam -o simulated_hg38_sorted.bam
$ samtools index -@ 8 simulated_hg38_sorted.bam simulated_hg38_sorted.bai
```

We randomly select a subset of loci from the VCF file to create a region file for evaluating misclus.

```bash
# misassemble clustering
$ misclus simulateReg ref_Spombe.fa simulate_Spombe_sorted.bam
```

The analysis results are saved in the `final_result` file in the `cluster_result` folder.



## Results

Evaluation results are saved within this repository with each tool's results saved in the subfolder named after the respective tool. 

## Contact ##

If you have problems or some suggestions, please contact: [xzhu@ytu.edu.cn](xzhu@ytu.edu.cn) without hesitation. 

---- Enjoy !!! -----
