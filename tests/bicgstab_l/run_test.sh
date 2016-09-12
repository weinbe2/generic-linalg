#! /bin/bash

# Run the bicgstab-l test, produce output files, make plots, etc.

make bicgstab_l

# Run a unit gauge field test with mass 1e-3. The beta = 7.0 tricks the program into using a unit gauge field... sneaky, I know.
./bicgstab_l --mass 1e-3 --beta 7.0 --lattice-size 64 > ./test_output/unit_output.out

# Run a beta = 6.0 gauge field test with mass 1e-3.
./bicgstab_l --mass 1e-3 --beta 6.0 --lattice-size 64 > ./test_output/b60_output.out

# Sanitize the output quickly.
./test_output/postappend_time.sh ./test_output/unit_output.out ./test_output/unit_output_fix.out
./test_output/postappend_time.sh ./test_output/b60_output.out ./test_output/b60_output_fix.out

# Isolate just the iterations and timings.
grep "Time" ./test_output/b60_output_fix.out > ./test_output/b60_summary.dat
grep "Time" ./test_output/unit_output_fix.out > ./test_output/unit_summary.dat

# Isolate the relative residual history for each solve separately.
grep "\[BICGSTAB\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab.dat
grep "\[BICGSTAB-1\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab1.dat
grep "\[BICGSTAB-2\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab2.dat
grep "\[BICGSTAB-4\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab4.dat
grep "\[BICGSTAB-8\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab8.dat
grep "\[BICGSTAB-16\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab16.dat

grep "\[BICGSTAB\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } ' > ./test_output/unit_detail_bicgstab.dat
grep "\[BICGSTAB-1\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/unit_detail_bicgstab1.dat
grep "\[BICGSTAB-2\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/unit_detail_bicgstab2.dat
grep "\[BICGSTAB-4\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/unit_detail_bicgstab4.dat
grep "\[BICGSTAB-8\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/unit_detail_bicgstab8.dat
grep "\[BICGSTAB-16\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/unit_detail_bicgstab16.dat

# Prepare plots.
gnuplot -e "load \"plots.plt\""

# Build plots.
function buildplot {
  epstopdf ${1}-inc.eps
  pdflatex ${1}.tex
  pdflatex ${1}.tex
  pdflatex ${1}.tex
}

buildplot unit_test
buildplot b60_test

# Clean up plot byproducts.

rm *.eps
rm *.tex
rm *.log
rm *.aux
rm *-inc.pdf
