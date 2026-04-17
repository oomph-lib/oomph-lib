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

# Plot all columns (from column 2 to n) against column 1
# Each column is plotted as a separate line
plot for [i=2:n] "build/RESLT/second_derivs_difference.dat" using 1:i with lines title sprintf("col %d", i-1)
