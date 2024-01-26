# Misassemly Classification Results for misClas Experiments

Misassemly Classification Results for misClas Experiments

## Prerequisites
### Install tools:

#### xxx

```bash
Debian / Ubuntu:

$ sudo apt update  # Ensure the package list is up to date
$ sudo apt install mlpack-bin g++ python3 autoconf automake libtool

RedHat / CentOS:

$ sudo yum install g++ python3 autoconf automake libtool

# htslib (1.9 version or later) be installed from source files (it cannot be install by apt install or yum install)
$ wget -c https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
$ tar -xvjf htslib-1.19.tar.bz2
$ cd htslib-1.14
$ autoreconf -i
$ ./configure
$ make 
$ make install
```

#### dokcer

```bash
$  docker pull xzhu2/misclus:0.5.0
$ docker run it --name your_contaniner_name misclus_iamges:0.5.0
$ cd /opt/misclus/
$ ./autogen.sh
```

### Download reference and annotations:

```
# Reference


```

