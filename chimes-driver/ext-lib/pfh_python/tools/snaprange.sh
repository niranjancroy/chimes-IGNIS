### this function uses grep to search through a set of logfiles and 
### find the first and last instance of a snapshot being finished in each logfile
### then uses those instances to construct a set of snapshot numbers to skip
### since a finished snapshot doesn't need to be touched again


### this function sometimes mangles ${fname} for some reason in certain
### shells.  I can't for the life of me figure out why, and after tearing 
### my hair out over it I'm ready to give up.

snaprange ()
{
    echo "This function is deprecated and may not work properly"
    echo "Recommend using ${PYDIR}/tools/identify_missing_frames.py to build a frame_skip_list"
    files=$(ls -v "$1"*);
    arangestring="skiplist = np.r_[ ";
    for fname in $files;    
    do
        low=$(fgrep "finished with snapshot " ${fname} | head -n 1);   #grep for snapshot file loaded and find the first
        low=${low/"...... finished with snapshot "/""};   #trim off first part the string
        lownum=${low% *}; #trim off parts after the space

        high=$(fgrep "...... finished with snapshot " ${fname} | tail -n 1); #do it again for the last snapshot loaded
        high=${high/"...... finished with snapshot "/""};
        highnum=${high% *};

        delta=$((highnum-lownum));
        if [ $delta -le 2 ]; then
            echo "${fname}:  ${lownum} -- ${highnum}; delta = ${delta} -- !!!";
        else
            echo "${fname}:  ${lownum} -- ${highnum}; delta = ${delta}";
            arangestring="$arangestring np.arange($lownum, $highnum),";
        fi;
    done;
    chrlen=${#arangestring};
    arangestring=${arangestring:0:$chrlen-1}" ]";
    echo "Copy and paste this arange to cover all the above with a delta > 2:";
    echo $arangestring
}
