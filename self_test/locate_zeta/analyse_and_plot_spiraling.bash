#! /bin/bash

#------------------------------------------------------
# Run this to post-process locate_zeta self-tests
#------------------------------------------------------

if [ ! -e Validation ]; then
    echo "Please run self-tests first"
    exit
fi

cd Validation

cd RESLT_2D
oomph-convert ref_bin.dat
oomph-convert -p2 sample_points_2d_ref_bin.dat
oomph-convert non_ref_bin.dat 
oomph-convert -p2 sample_points_2d_nonref_bin.dat 
oomph-convert -p2 sample_points_2d_cgal.dat
cp ../../spiraling_2D.pvsm .
cp ../../distance_2d.lay .
echo " "
echo "================================================================="
echo "Post-process with:"
echo " "
echo "     cd "`pwd`
echo "     paraview --state=spiraling_2D.pvsm "
echo "     tecplot distance_2d.lay"
echo " "
echo "================================================================="
echo " "
cd ..

cd RESLT_triangle
oomph-convert ref_bin.dat
oomph-convert -p2 sample_points_2d_ref_bin.dat
oomph-convert non_ref_bin.dat 
oomph-convert -p2 sample_points_2d_nonref_bin.dat 
oomph-convert -p2 sample_points_2d_cgal.dat

cp ../../spiraling_2D.pvsm .
cp ../../distance_2d.lay .
echo " "
echo "================================================================="
echo "Post-process with:"
echo " "
echo "     cd "`pwd`
echo "     paraview --state=spiraling_2D.pvsm "
echo "     tecplot distance_2d.lay"
echo " "
echo "================================================================="
echo " "
cd ..

cd RESLT_3D
oomph-convert ref_bin.dat
oomph-convert -p3 sample_points_3d_ref_bin.dat
oomph-convert non_ref_bin.dat 
oomph-convert -p3 sample_points_3d_nonref_bin.dat 
oomph-convert -p3 sample_points_3d_cgal.dat
cp ../../spiraling_3D.pvsm .
cp ../../distance_3d.lay .
echo " "
echo "================================================================="
echo "Post-process with:"
echo " "
echo "     cd "`pwd`
echo "     paraview --state=spiraling_3D.pvsm "
echo "     tecplot distance_3d.lay"
echo " "
echo "================================================================="
echo " "
cd ..


cd RESLT_tetgen
oomph-convert ref_bin.dat
oomph-convert -p3 sample_points_3d_ref_bin.dat
oomph-convert non_ref_bin.dat 
oomph-convert -p3 sample_points_3d_nonref_bin.dat 
oomph-convert -p3 sample_points_3d_cgal.dat
cp ../../spiraling_3D.pvsm .
cp ../../distance_3d.lay .
echo " "
echo "================================================================="
echo "Post-process with:"
echo " "
echo "     cd "`pwd`
echo "     paraview --state=spiraling_3D.pvsm "
echo "     tecplot distance_3d.lay"
echo " "
echo "================================================================="
echo " "
cd ..
