set style line 1 lt 1 lc rgb '#000000' lw 2 pt 11 ps 0.5
set style line 2 lt 2 lc rgb '#B8860B' lw 2 pt 8 ps 0.5
set style line 3 lt 3 lc rgb '#CC0000' lw 2 pt 6 ps 0.5
set style line 4 lt 4 lc rgb '#0000CC' lw 2 pt 4 ps 0.5
set style line 5 lt 5 lc rgb '#009900' lw 2 pt 12 ps 0.5
set style line 11 lt 6 lc rgb '#555555' lw 2 pt 10 ps 0.5
set style line 12 lt 7 lc rgb '#D0A015' lw 2 pt 8 ps 0.5
set style line 13 lt 8 lc rgb '#885555' lw 2 pt 7 ps 0.5
set style line 14 lt 9 lc rgb '#5555FF' lw 2 pt 5 ps 0.5
set style line 15 lt 1 lc rgb '#558855' lw 2 pt 13 #ps 0.5

set key top right Left reverse font ",8"
set terminal epslatex color standalone

set log y

# Unit tests.

set xlabel "Iteration"
set ylabel "Relative Residual"
set yrange [1e-11:1e2]
set title "Free 2D Square Laplace, $64^3$,\n$m = 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1$"
set output "free_laplace_cg.tex"
plot "./cgm_out.dat" using 1 ls 1 w l title "CG-M", "./cg_out.dat" using 1 ls 14 w l title "Sequential CG"
set output
