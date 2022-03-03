#! /bin/bash
for folder in /data/rewrite/*; do
    f_name="$(basename -- $folder)"
    input=/data/rewrite/$f_name
    output=/data/zipped_1/$f_name
    outstore=/data/feh230_updated_1
    if [ ! -f $outstore/$f_name.tar ] 
    then
    	echo $outstore/$f_name.tar
        mkdir $output

        cd $input

        for fname in $input/*; do

            gzip $fname
    
            mv $fname.gz $output
    
        done

        cd $output

        tar -cf $f_name.tar *.gz

        mv $f_name.tar $outstore
    fi
done
