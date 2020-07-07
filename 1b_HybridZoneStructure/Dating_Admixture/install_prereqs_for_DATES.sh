# install gsl 
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz; tar -xvf gsl-2.5.tar.gz
./configure --prefix=/data/tunglab/tpv/Programs/gsl-2.5/; make; make install 
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.3.tar.gz; tar -xvf gsl-2.3.tar.gz
./configure --prefix=/data/tunglab/tpv/Programs/gsl-2.3/; make; make install 
# install fftw 
wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.5.tar.gz; tar -xvf fftw-3.3.5.tar.gz
cd fftw-3.3.5; ./configure --prefix=/data/tunglab/tpv/Programs/fftw-3.3.5/; make; make install 
wget http://www.fftw.org/fftw-3.3.3.tar.gz; tar -xzf fftw-3.3.3.tar.gz
cd fftw-3.3.3; ./configure --enable-shared --prefix=/data/tunglab/tpv/Programs/fftw-3.3.3/; make; make install 
# Install DATES 
cd /data/tunglab/tpv/dating_admixture/DATES-master/src/
vi Makefile ## Edits to the makefile
		override CFLAGS += -I/nfs/software/helmod/apps/Core/OpenBLAS/0.2.20-gcb01/include -I/data/tunglab/tpv/Programs/gsl-2.3/include -I/data/tunglab/tpv/Programs/fftw-3.3.3/include
		override LDFLAGS += -L/nfs/software/helmod/apps/Core/OpenBLAS/0.2.20-gcb01/lib -L/data/tunglab/tpv/Programs/gsl-2.3/lib/  -L/data/tunglab/tpv/Programs/fftw-3.3.3/lib
make install 

