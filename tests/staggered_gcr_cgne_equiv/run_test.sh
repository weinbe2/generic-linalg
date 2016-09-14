#! /bin/bash

# Run the bicgstab-l test, produce output files, make plots, etc.

make clean; make

# Run a unit gauge field test with mass 1e-3. The beta = 7.0 tricks the program into using a unit gauge field... sneaky, I know.
./gcr_cgne_equiv --mass 1e-3 --beta 7.0 --lattice-size 64 > ./test_output/unit_output.out

# Run a beta = 6.0 gauge field test with mass 1e-3.
#./bicgstab_l --mass 1e-3 --beta 6.0 --lattice-size 64 > ./test_output/b60_output.out

# Sanitize the output quickly.
./test_output/postappend_time.sh ./test_output/unit_output.out ./test_output/unit_output_fix.out
#./test_output/postappend_time.sh ./test_output/b60_output.out ./test_output/b60_output_fix.out

# Isolate just the iterations and timings.
#grep "Time" ./test_output/b60_output_fix.out > ./test_output/b60_summary.dat
grep "Time" ./test_output/unit_output_fix.out > ./test_output/unit_summary.dat

# Put some outputs from this into the bare plots file.

utcg=$(grep "\[CG_NORMAL\]:" ./test_output/unit_summary.dat | awk ' { printf "%3.2f", $12 } ')
utcgnd=$(grep "\[CG_NORMAL_NODDAG\]:" ./test_output/unit_summary.dat | awk ' { printf "%3.2f", $12 } ')
utgcr=$(grep "\[GCR\]:" ./test_output/unit_summary.dat | awk ' { printf "%3.2f", $12 } ')
utgmres=$(grep "\[GMRES\]:" ./test_output/unit_summary.dat | awk ' { printf "%3.2f", $12 } ')

#btcg=$(grep "\[CG_NORMAL\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')
#bt0=$(grep "\[BICGSTAB\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')
#bt1=$(grep "\[BICGSTAB-1\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')
#bt2=$(grep "\[BICGSTAB-2\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')
#bt4=$(grep "\[BICGSTAB-4\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')
#bt8=$(grep "\[BICGSTAB-8\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')
#bt16=$(grep "\[BICGSTAB-16\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')
#btgcr=$(grep "\[GCR\]:" ./test_output/b60_summary.dat | awk ' { printf "%3.2f", $12 } ')

sed -e "s/:UNITGCR:/${utgcr}/g" -e "s/:UNITCG:/${utcg}/g" -e "s/:UNITCGND:/${utcgnd}/g" -e "s/:UNITGMRES:/${utgmres}/g" ./test_output/plots_bare.plt > ./test_output/plots.plt

# Isolate the relative residual history for each solve separately.
#grep "\[CG_NORMAL\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_cg.dat
#grep "\[GCR\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_gcr.dat
#grep "\[BICGSTAB\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab.dat
#grep "\[BICGSTAB-1\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab1.dat
#grep "\[BICGSTAB-2\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab2.dat
#grep "\[BICGSTAB-4\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab4.dat
#grep "\[BICGSTAB-8\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab8.dat
#grep "\[BICGSTAB-16\]:" ./test_output/b60_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } '  > ./test_output/b60_detail_bicgstab16.dat

grep "\[CG_NORMAL\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } ' > ./test_output/unit_detail_cg.dat
grep "\[GCR\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } ' > ./test_output/unit_detail_gcr.dat
grep "\[GMRES\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } ' > ./test_output/unit_detail_gmres.dat
grep "\[CG_NORMAL_NODDAG\]:" ./test_output/unit_output_fix.out | grep -v "Success" | awk ' { print $4,$6,$8 } ' > ./test_output/unit_detail_cgnd.dat

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

buildplot unit_test
#buildplot b60_test

# Clean up plot byproducts.

rm *.eps
rm *.tex
rm *.log
rm *.aux
rm *-inc.pdf
rm plots.plt

cd ..
