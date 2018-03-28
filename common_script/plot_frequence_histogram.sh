if [ $# != 2 ]; then
    echo "Usage: $0 data.file scale"
    exit 0
fi

cat $1 | awk -v scale=$2 '{b = int($0/scale); if(b > 50) b = 50; if (b!= 50) a[int(b)] += 1} END{ for(i in a) print i*scale" "a[i]}' > $1.data

echo "
    set term jpeg 
    set output \"${1}.jpg\"
    set size 1, 1
    set xlabel \"coverage\"
    set ylabel \"times\"
    set key right bottom
    plot \"${1}.data\" title \"dis&franquce\" with impulse 
    " | gnuplot 
#rm a.data
