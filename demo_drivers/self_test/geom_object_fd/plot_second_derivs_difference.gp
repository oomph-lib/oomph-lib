# This script plots the finite-difference derivative errors against the FD step size
# in order to help identify the optimal step size for the derivative approximation.
#
# It is used for the second derivative. The corresponding section of the code 
# (in geom_object_fd.cc) has been commented out. You must uncomment that section to 
# generate the data files.

reset

# Set the output terminal and image size
set terminal pngcairo size 800,600

# Uncomment this for interactive window display
#set terminal qt

# Specify the output file
set output sprintf("build/RESLT/plot_d2.png")

# Read data file statistics to determine number of columns
stats "build/RESLT/second_derivs_difference.dat" nooutput
n = STATS_columns

# Set axis labels
set xlabel "dzeta"
set ylabel "difference"

# Use logarithmic scale for both axes
set log x
set log y

# Format x-axis tics in scientific notation using powers 
# of 10 (for log scale)
set format x "1×10^{%L}"

# Plot absolute differences between paired columns
# Each pair (2*i+2, 2*i+3) represents FD vs exact values
# taking |FD - exact| to show magnitude of numerical error
plot for [i=0:(n-3)/2] \
    "build/RESLT/second_derivs_difference.dat" \
    using 1:(abs(column(2*i+2) - column(2*i+3))) \
    with lines lw 2 title sprintf("col %d", i+1)
