#! /bin/bash

#Run or plot?
plot_list="0 1"

# Run plot_it.bash to re-generate data?
run_plot_it=1

# Regenerate pngs (i.e. pvsm or py file changed)?
regenerate_png=1

# Parameter study:
lambda_sq_list="0.5 1.0" 
alpha_list="0.5" 
density_ratio_list="1.0" 

# Number of pics in row
n_pics_in_row=2

# Number of rows per page
n_rows_per_page=2


start_dir=`pwd`
for plot in `echo $plot_list`; do

# Directory to store data
dir=./NEW_RUNS

if [ $plot -ne 1 ]; then
    
    if [ -e `echo $dir` ];  then
        echo " "
        echo "ERROR: Please delete directory $dir and try again"
        echo " "
        exit
    fi
    mkdir $dir
    
    
    # Rebuild executable
    make adaptive_unstructured_two_d_curved
    
    # Move to working directory
    cp adaptive_unstructured_two_d_curved       $dir
    cp unstructured_two_d_curved.cc    $dir

else

    cp plot_it.bash                    $dir
    cp animate.pvsm                    $dir
    cp animate.py                      $dir
    
fi

cd $dir

if [ $plot -eq 1 ]; then
    

    echo "\documentclass{beamer} " > compare.tex
    echo "\usepackage{beamerthemesplit} " >> compare.tex
    echo "\usepackage{ifthen} " >> compare.tex
    echo "\usepackage{multimedia} " >> compare.tex
    echo "\usepackage{xmpmulti} " >> compare.tex
    echo "\usepackage{ulem} " >> compare.tex
    echo "\usepackage{graphicx} " >> compare.tex
    echo "\usepackage{amsbsy} " >> compare.tex
    echo "\usepackage{amssymb} " >> compare.tex
    echo "\usepackage{rotating} " >> compare.tex
    echo "\usepackage{color} " >> compare.tex
    echo "\usepackage{animate} " >> compare.tex
    echo "\providecommand\thispdfpagelabel[1]{} " >> compare.tex
    echo "\begin{document} " >> compare.tex
    echo "\frame[plain]{" >> compare.tex
    echo "\begin{overlayarea}{\textwidth}{\textheight} " >> compare.tex 
    echo "\vspace{-0.3cm} " >> compare.tex


fi

mbox_open=0
col_count=0
row_count=0

for lambda_sq in `echo $lambda_sq_list`; do
    for alpha in `echo $alpha_list`; do
        for density_ratio in `echo $density_ratio_list`; do
            
            dir=R_lam_sq`echo $lambda_sq`_alp`echo $alpha`_den_rat`echo $density_ratio`


            if [ $plot -eq 0 ]; then
                mkdir $dir
                ./adaptive_unstructured_two_d_curved --t_tanh 1.5 --alpha_tanh 5.0 --n_steps 2000 --dt 0.01 --lambda_sq $lambda_sq --alpha $alpha --density_ratio $density_ratio  --dir $dir > `echo $dir`/OUTPUT &
                
            else
                cd $dir
                
                
                if [ $run_plot_it -eq 1 ]; then
                    ../plot_it.bash
                fi
                if [ $regenerate_png -eq 1 ]; then
                    cp ../animate.pvsm .
                    pvbatch ../animate.py
                fi
                output_label=1
                count=0
                rm -f frame*png
                while [  $count -lt 1000 ]; do
                    if [ $count -lt 10 ]; then
                        file=animation.000`echo $count`.png
                        if [ -e $file ]; then
                            ln -sf $file frame`echo $output_label`.png
                            let output_label=output_label+1
                        fi
                    else 
                        if [  $count -lt 100 ]; then
                            file=animation.00`echo $count`.png
                            if [ -e $file ]; then
                                ln -sf $file frame`echo $output_label`.png
                                let output_label=output_label+1
                            fi        
                        else 
                            file=animation.0`echo $count`.png
                            if [ -e $file ]; then
                                ln -sf $file frame`echo $output_label`.png
                                let output_label=output_label+1
                            fi
                        fi
                    fi
                    let count=count+1
                done
                
                let total_number=output_label-1
                subdir=`pwd`
                
                
                echo "col count: " $col_count
                if [ $col_count -eq 0 ]; then
                    echo "doig it...: " 
                    echo "\mbox{\hspace{-1cm} " >> ../compare.tex
                    mbox_open=1
                fi
                echo "\begin{minipage}[h!]{0.36\textwidth}" >> ../compare.tex
                echo "\scalebox{0.13}{\animategraphics[autoplay,loop,controls]{10}{"$subdir"/frame}{1}{"$total_number"}}" >> ../compare.tex 
                
                
                echo "\scalebox{0.7}{\mbox{\tiny $\Lambda^2="$lambda_sq", \alpha="$alpha", \frac{\rho_{\rm f}}{\rho_{\rm s}}="$density_ratio"$}}" >> ../compare.tex
                echo "\end{minipage}" >> ../compare.tex
                
                let col_count=col_count+1
                if [ $col_count -eq  $n_pics_in_row ]; then
                    col_count=0
                    echo "} " >> ../compare.tex
                    mbox_open=0
                    let row_count=row_count+1
                    
                    echo "ROW COUNT "$row_count " " $n_rows_per_page
                    if [ $row_count -eq $n_rows_per_page ]; then
                        echo "\end{overlayarea} " >> ../compare.tex
                        echo "} " >> ../compare.tex
                        echo "\frame[plain]{" >> ../compare.tex
                        echo "\begin{overlayarea}{\textwidth}{\textheight} " >> ../compare.tex 
                        echo "\vspace{-0.3cm} " >> ../compare.tex
                        row_count=0
                    fi
                fi
                cd ..  
            fi
        done
    done
done

if [ $plot -eq 1 ]; then
    
    if [ $mbox_open -eq 1 ]; then
        echo "} " >> compare.tex
    fi
    echo "\end{overlayarea} " >> compare.tex
    echo "} " >> compare.tex
    echo "\end{document} " >> compare.tex

    pdflatex compare
fi


echo "WAITING FOR BG JOBS TO FINISH"
wait
echo "CONTINUING"

cd $start_dir

done

