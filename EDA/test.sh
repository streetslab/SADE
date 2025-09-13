#!/bin/bash

THREDS=100

while getopts "e:f:" opt; do
	case $opt in 
		e) threds="$OPTARG"
			;;
		f) file="$OPTARG"
			;;
	        \?) echo "Invalid option: -$OPTARG" >&2
		    exit 1
			;;
	esac
done


if [[ -z ${threds+x} ]] ; then
	threds=$THREDS
fi

echo "file is ${file}"
echo "threds is ${threds}"

