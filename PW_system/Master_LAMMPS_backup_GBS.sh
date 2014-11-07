#!/bin/sh
###################################################################
#                                                                 #
# BACKUP script of lammps-generated files for masterproject       #
#                                                                 #
# Masterproject in molecular dynamics for Goeran Brekke Svaland   #
#                                                                 #
###################################################################

# What to back up
backup_files="/home/goran/lammps-28Jun14/examples/water_portlandite_system"

# Where to back up to
dest="/media/goran/9196bccb-acd5-42ba-ae64-4231e568fd1d/home/goran/LAMMPS_backup"

# Create archive filename
day=$(date +%F)
hostname=$(hostname -s)
archive_file="$hostname-$day.tgz"

# Print start status message
echo "Backing up $backup_files to $dest/$archive_file"
date
echo

# Back up the files using tar. c-create v-verbose z-gzip f-archive (create archive using gzip verbosly)
tar cvzf $dest/$archive_file $backup_files

# Print end status message
echo
echo "Backup finished"
date

# Long listing of files in $dest to check file sizes
ls -lh $dest

# we should think about removing older archive files..
# use newest archive to process information and generate plots in ipython notebook

