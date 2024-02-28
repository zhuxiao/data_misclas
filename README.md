# Misassembly clustering Results for misClus Experiments

## Prerequisites
#### Install tools:

```bash
# Debian / Ubuntu:
$ sudo apt update  # Ensure the package list is up to date
$ sudo apt install mlpack-bin g++ python3 autoconf automake libtool samtools openjdk-11-jdk
```

```bash
# RedHat / CentOS:
$ sudo yum install gcc-c++ python3 autoconf automake libtool samtools java-11-openjkd-devel
```

 htslib (1.9 version or later)  and  bwa be installed from source files (it cannot be install by apt install or yum install)

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

misassembly detection tools 

```bash
# misasm
$ git clone https://github.com/zhuxiao/misasm.git
$ cd misasm/src/
$ make

# misfinder 
$ git clone https://github.com/zhuxiao/misFinder.git
$ cd misFinder
$ ./autogen.sh

# Pilon 
$ wget -c https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar

#Quast 
$ wget -c https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
$ tar -zvxf quast-5.2.0.tar.gz
$ cd quast-5.2.0
$ ./setup.py 
```

## Experiments

#### Escherichia coli K-12 MG1655

```bash
# Download E. coli K-12 MG1655  
$ prefetch SRR12717904
# Assemble paired-end reads to contigs
$ fasterq-dump --split-files SRR12717904.sra
$ fq2fa --merge SRR12717904_1.fq SRR12717904_2.fq SRR12717904.fa
$ idba_ud --mink 20 --maxk 120 --min_contig 1000 --num_threads 8 -r SRR12717904.fa -o ecoli_idba_out
spades
# Align paired-end reads onto contigs
$ cd ecoli_idba_out
$ bwa index scaffold.fa
$ bwa mem -t 8 scaffold.fa SRR12717904_1.fq SRR12717904_2.fq > E.coli.sam
$ samtools view -b -S -@ 8 E.coli.sam > E.coli.bam
$ samtools sort -@ 8 E.coli.bam -o E.coli_sorted.bam
$ samtools index -@ 8 E.coli_sorted.bam E.coli_sorted.bai
```

```bash
# misasm candidate misassemblies 
$ misasm scaffold.fa E.coli_sorted.bam
$ cd output
$ cat genome_Indel genome_Misjoin | awk '{split($1,a,":"); split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > MisasmRaw
$ filterOverlapReg.py MisasmRaw scaffold.fa.fai MisasmReg
#Misassembly clustering
$ misclus MisasmReg scaffold.fa E.coli_sorted.bam
```

```bash
# misfinder candidate misassemblies  
$ mf all -t 8 -o output config 
$ cd output
$ cat misassemblies_misFinder | awk '{print $1"\t"$3"\t"$4}' > misFindeRaw
$ filterOverlapReg.py misFindeRaw scaffold.fa.fai misFindeReg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa E.coli_sorted.bam
```

```bash
# Pilon candidate misassemblies
$ java -jar pilon --genome scaffold.fa --frags S.pombe_sorted.bam --changes > misassemblies_Pilon
$ cat misassemblies_Pilon | grep "fix" | grep ":" | cut -d ":" -f 2,3 | awk '{split($1,a,":");split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > PilonRaw
$ filterOverlapReg.py PilonRaw scaffold.fa.fai PilonReg
# Misassembly clustering
$ misclus PilonReg scaffold.fa E.coli_sorted.bam
```

```bash
# Quast candidate misassemblies
$ quast.py -r NC_000913.3.fa --extensive-mis-size 200 scaffold.fa -m 1000
$ cd ./quast_result/latest
$ extractQUASTreg.py contigs_report_scaffold.stdout > QUASTraw
$ filterOverlapReg.py QUASTraw scaffold.fa.fai QUASTreg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa E.coli_sorted.bam
```



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

```bash
# misasm candidate misassemblies
$ misasm scaffold.fa S.pombe_sorted.bam
$ cd output
# Classify misassemblies
$ cat genome_Indel genome_Misjoin | awk '{split($1,a,":"); split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > MisasmRaw
$ filterOverlapReg.py MisasmRaw scaffold.fa.fai MisasmReg
#Misassembly clustering
$ misclus MisasmReg scaffold.fa S.pombe_sorted.bam
```

```bash
# misfinder candidate misassemblies  
$ mf all -t 8 -o output config 
$ cd output
$ cat misassemblies_misFinder | awk '{print $1"\t"$3"\t"$4}' > misFindeRaw
$ filterOverlapReg.py misFindeRaw scaffold.fa.fai misFindeReg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa S.pombe_sorted.bam
$ cat result_errors | awk '{print $1"\t"$3"\t"$4}' > misFindeRaw
```

```bash
# Pilon candidate misassemblies
$ java -jar pilon --genome scaffold.fa --frags S.pombe_sorted.bam --changes > misassemblies_Pilon
$ cat misassemblies_Pilon | grep "fix" | grep ":" | cut -d ":" -f 2,3 | awk '{split($1,a,":");split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > PilonRaw
$ filterOverlapReg.py PilonRaw scaffold.fa.fai PilonReg
# Misassembly clustering
$ misclus PilonReg scaffold.fa S.pombe_sorted.bam
```

```bash
# Quast candidate misassemblies
$ quast.py -r GCF_000002945.1_ASM294v2_genomic.fa --extensive-mis-size 200 scaffold.fa -m 1000
$ cd ./quast_result/latest
$ extractQUASTreg.py misassemblies_QUAST > QUASTraw
$ filterOverlapReg.py QUASTraw scaffold.fa.fai QUASTreg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa S.pombe_sorted.bam
```



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

```bash
# misasm candidate misassemblies
$ misasm scaffold.fa hg002_sorted.bam
# Classify misassemblies
$ cat genome_Indel genome_Misjoin | awk '{split($1,a,":"); split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > MisasmRaw
$ filterOverlapReg.py MisasmRaw scaffold.fa.fai MisasmReg
#Misassembly clustering
$ misclus MisasmReg scaffold.fa $ misclus MisasmReg scaffold.fa hg002_sorted.bam 
```

```bash
# Pilon candidate misassemblies ？
$ java -jar pilon --genome scaffold.fa --frags S.pombe_sorted.bam --changes > misassemblies_Pilon
$ cat misassemblies_Pilon | grep "fix" | grep ":" | cut -d ":" -f 2,3 | awk '{split($1,a,":");split(a[2],b,"-");print a[1]"\t"b[1]"\t"b[2]}' > PilonRaw
$ filterOverlapReg.py PilonRaw scaffold.fa.fai PilonReg
# Misassembly clustering
$ misclus PilonReg scaffold.fa hg002_sorted.bam
```

```bash
# Quast candidate misassemblies ？
$ quast.py -r hs38.fa --extensive-mis-size 200 scaffold.fa -m 1000
$ cd ./quast_result/latest
$ extractQUASTreg.py misassemblies_QUAST > QUASTraw
$ filterOverlapReg.py QUASTraw scaffold.fa.fai QUASTreg
# Misassembly clustering
$ misclus misFindeReg scaffold.fa hg002_sorted.bam
```



## Result

Typically, evaluation results are saved within [zhuxiao/misclus-experiments: Data used in misClus experiments (github.com)](https://github.com/zhuxiao/misclus-experiments), with each tool's results saved in the subfolder named after the respective tool. 

## Contact ##

If you have problems or some suggestions, please contact: [xzhu@ytu.edu.cn](xzhu@ytu.edu.cn) without hesitation. 

---- Enjoy !!! -----
