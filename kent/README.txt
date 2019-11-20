
This is a trimmed down version of the Kent source code, which has fewer 
dependencies. It will only compile selected tools for generating or working
with genome alignment chains and nets. The code of these tools remains unchanged
when compared with the original UCSC version.

Original Kent source code was cloned from (on Aug 26, 2019):

https://github.com/ucscGenomeBrowser/kent.git
[ commit de59f1cfae7ff48a75874ff8b19b8f0721222ac9 ]

'README_ucsc' is the original README file provided by the UCSC. 
Please, do NOT follow the install instructions given in 'README_ucsc'.

<Tue Sep 10 14:57:10 CEST 2019, Nikolai Hecker, hecker@mpi-cbg.de>



Dependencies
////////////

* GCC including development environment; tested with: gcc (GCC) 6.2.0
- math library: libm

* openSSL libraries including development environment (e.g. libssl-dev): libssl and libcrypto 

* compression library zlib including development environment (e.g. zlib1g-dev): libz 

* uuid library including development environment (e.g. uuid-dev): libuuid

* pthreads library: lpthreads





INSTALLATION
////////////

Let us assume you cloned this trimmed down version of the Kent source code,
for example, to: /home/yourname/source/kent/

1) change into the source directory, e.g.:

$ cd /home/yourname/source/kent/
$ cd src


2) compile binaries using 'make'
$ make

binaries will be inside the 'bin' directory in the Kent parent folder, 
e.g. /home/yourname/source/kent/bin/


3) (optional) export binaries to make them globally accessible
If you are using BASH, for example, add following lines to your ~/.bashrc:

# Kent binaries
PATH=$PATH:/home/yourname/source/kent/bin;export PATH

and source your .bashrc:
$ source ~/.bashrc

