Note: Please update the "LIB" line with your own install directory for ARPACK!

How to compile + install ARPACK:

Inside a useful directory, run:

mkdir build    # where build will happen
mkdir install  # where install will go
git clone https://github.com/pv/arpack-ng.git
cd build
./configure --prefix=[install directory, absolute path] {--with-blas=[absolute path to blas lib directory] --with-lapack=[absolute path to lapack lib directory]}
make
make install