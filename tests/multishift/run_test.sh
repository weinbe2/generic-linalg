#! /bin/bash

# Run the bicgstab-l test, produce output files, make plots, etc.

make multishift

# Run the multishift test, save to file.
./multishift > ./test_output/mshift.out

# Isolate the relative residual history for each method.
grep "\[CG-M\]:" ./test_output/mshift.out | grep -v "Success" | awk ' { print $8 } ' > ./test_output/cgm_out.dat
grep "\[CG_mass" ./test_output/mshift.out | grep -v "Success" | awk ' { print $8 } ' > ./test_output/cg_out.dat

cd test_output

# Prepare plots.
gnuplot -e "load \"plots.plt\""

# Build plots.
function buildplot {
  epstopdf ${1}-inc.eps
  pdflatex ${1}.tex
  pdflatex ${1}.tex
  pdflatex ${1}.tex
}

buildplot free_laplace_cg

# Clean up plot byproducts.

rm *.eps
rm *.tex
rm *.log
rm *.aux
rm *-inc.pdf

cd ..
