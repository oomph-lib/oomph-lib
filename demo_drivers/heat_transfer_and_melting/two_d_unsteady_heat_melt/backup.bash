dirname="/home/mheil/backup_"`date | awk '{print $4}'`

echo "Making backup in $dirname"
mkdir $dirname
cp *.bash *.pvsm *.h *.cc *.map $dirname


