reset
set term gif animate
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb"#ccffcc" behind 
set output "animate4.gif"
set style data dots
set xrange [-1.0:1.0]
set yrange [-1.0:1.0]
set size square
plot "vortices0.txt"
plot "vortices1.txt"
plot "vortices2.txt"
plot "vortices3.txt"
plot "vortices4.txt"
plot "vortices5.txt"
plot "vortices6.txt"
plot "vortices7.txt"
plot "vortices8.txt"
plot "vortices9.txt"
