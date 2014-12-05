#!/bin/bash
################################################################
#                                                              #
#                     BACKUP script                            #
#                                                              #
# When ran, it takes backup of the input filelocation in the   #
# desired folder location.                                     #
# Be wared! If a too large file-location is choosen, this will #
# take a fucking loooong time!                                 #
# When appropriate execute (chmod) access is set, the script   #
# is to be executed like this:                                 #
#                                                              #
# >> ./Abel_backup.sh path_to_folder_that_you_wish_to_back_up  #
#                                                              #
# default path is:                                             #
# ~/lammps-28jun14/examples/Abel_runs/PW_system/flat_system    #
################################################################

default_path="/home/goran/lammps-28jun14/examples/Abel_runs"
dest="/media/goran/9196bccb-acd5-42ba-ae64-4231e568fd1d/home/goran/LAMMPS_backup"

if [[ $1 == /home/goran* ]]
then
	echo "$1 ...seems to be a valid adress... proceeding"
	backup_files=$1
else
	echo "argument given is not a valid path. path must start with /home/goran"
	echo "agument given on commandline was:"
	echo "$1"
	echo "using default path = $default_path"
	backup_files=$default_path
fi

## Destination of archive file:
day=$(date +%F)
hostname=$(hostname -s)         ## name of the computer
base=${backup_files##*/}        ## extracts the part of the string after the last backslash

archive_file="$base-$day.tgz"   ## name of archive file


## Printing start status message...
echo "Backing up $backup_files to $dest/$archive_file"
date
echo

## Back up the files using tar. c-create, v-verbose, z-gzip, f-archivefile
## creates archive using gzip verbosly

tar cvzf $dest/$archive_file $backup_files

## print end status message

echo 
echo "Backup finished.."
date

## Long listing of files in $dest to check file sizes:
ls -lh $dest


