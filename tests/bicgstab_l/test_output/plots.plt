set style line 1 lt 1 lc rgb '#000000' lw 2 pt 11 ps 0.5
set style line 2 lt 2 lc rgb '#B8860B' lw 2 pt 8 ps 0.5
set style line 3 lt 3 lc rgb '#CC0000' lw 2 pt 6 ps 0.5
set style line 4 lt 4 lc rgb '#0000CC' lw 2 pt 4 ps 0.5
set style line 5 lt 5 lc rgb '#009900' lw 2 pt 12 ps 0.5
set style line 11 lt 6 lc rgb '#555555' lw 2 pt 10 ps 0.5
set style line 12 lt 1 lc rgb '#D0A015' lw 2 pt 8 #ps 0.5
set style line 13 lt 1 lc rgb '#885555' lw 2 pt 7 #ps 0.5
set style line 14 lt 1 lc rgb '#5555FF' lw 2 pt 5 #ps 0.5
set style line 15 lt 1 lc rgb '#558855' lw 2 pt 13 #ps 0.5

set key top right Left reverse
set terminal epslatex color standalone

set log y

# Unit tests.

set xlabel "Iteration"
set ylabel "Relative Residual"
set yrange [1e-11:1e3]
set title "Unit gauge, $64^2$, $m = 10^{-3}$, RelRes $10^{-10}$"
set output "unit_test.tex"
plot "./unit_detail_bicgstab.dat" using 1:3 ls 1 w l title "BiCGStab", "./unit_detail_bicgstab1.dat" using 1:3 ls 2 w l title "BiCGStab-1", "./unit_detail_bicgstab2.dat" using 1:3 ls 3 w l title "BiCGStab-2", "./unit_detail_bicgstab4.dat" using 1:3 ls 4 w l title "BiCGStab-4", "./unit_detail_bicgstab8.dat" using 1:3 ls 5 w l title "BiCGStab-8", "./unit_detail_bicgstab16.dat" using 1:3 ls 11 w l title "BiCGStab-16"
set output

# beta 6.0 tests.


set xlabel "Iteration"
set ylabel "Relative Residual"
set yrange [1e-11:1e3]
set title "$\\beta = 6.0$, $64^2$, $m = 10^{-3}$, RelRes $10^{-10}$"
set output "b60_test.tex"
plot "./b60_detail_bicgstab.dat" using 1:3 ls 1 w l title "BiCGStab", "./b60_detail_bicgstab1.dat" using 1:3 ls 2 w l title "BiCGStab-1", "./b60_detail_bicgstab2.dat" using 1:3 ls 3 w l title "BiCGStab-2", "./b60_detail_bicgstab4.dat" using 1:3 ls 4 w l title "BiCGStab-4", "./b60_detail_bicgstab8.dat" using 1:3 ls 5 w l title "BiCGStab-8", "./b60_detail_bicgstab16.dat" using 1:3 ls 11 w l title "BiCGStab-16"
set output
