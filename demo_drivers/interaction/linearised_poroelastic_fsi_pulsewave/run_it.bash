#! /bin/bash



#Number of timesteps
nstep=1000

#Run or plot?
plot_list="0 1"

# PVMS file
pvsm_file=animate.pvsm

# Run plot_it.bash to re-generate data?
run_plot_it=1

# Regenerate pngs (i.e. pvsm or py file changed)?
regenerate_png=1

# Permeability variations
perm_list="0.005 0.05 0.5 5.0" # 

# Pin darcy?
pin_darcy_list="0"

# Number of pics in row
n_pics_in_row=2

# Number of rows per page
n_rows_per_page=2

start_dir=`pwd`
for plot in `echo $plot_list`; do

# Directory to store data
dir=./PermRangeBeforeCheckIn/NEW_RUNS #./NEW_RUNS

if [ $plot -ne 1 ]; then
    
    if [ -e `echo $dir` ];  then
        echo " "
        echo "ERROR: Please delete directory $dir and try again"
        echo " "
        exit
    fi
    mkdir $dir
    
    
    # Rebuild executable
    make linearised_poroelastic_fsi_pulsewave
    
    # Move to working directory
    cp linearised_poroelastic_fsi_pulsewave       $dir
    cp linearised_poroelastic_fsi_pulsewave.cc    $dir

else

    cp plot_it.bash                    $dir
    cp $pvsm_file                      $dir/fpsi.pvsm
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



for perm in `echo $perm_list`; do
for pin_darcy in `echo $pin_darcy_list`; do

    dir=RESLT_Perm`echo $perm`

    pin_darcy_flag=""
    pin_darcy_caption=" [Darcy switched on] "
    if [ $pin_darcy -eq 1 ]; then
        dir=`echo $dir`"_pin_darcy"
        pin_darcy_flag=" --pin_darcy "
        pin_darcy_caption=" [Darcy switched off] "
    else
        dir=`echo $dir`"_no_pin_darcy"
    fi


    if [ $plot -eq 0 ]; then
        mkdir $dir
        ./linearised_poroelastic_fsi_pulsewave $pin_darcy_flag  --nstep $nstep --permeability $perm --dir $dir > `echo $dir`/OUTPUT &
        
    else
        cd $dir
        
        
        if [ $run_plot_it -eq 1 ]; then
            ../plot_it.bash
        fi
        if [ $regenerate_png -eq 1 ]; then
            cp ../fpsi.pvsm .
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
        
        
        echo "\scalebox{0.7}{\mbox{\tiny Permeability="$perm"; "$pin_darcy_caption"}}" >> ../compare.tex
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

