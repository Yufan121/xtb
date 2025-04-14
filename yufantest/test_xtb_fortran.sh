# Set the XTBPATH environment variable
export XTBPATH=$2

# Obtain name without extension for $1
filename=$(basename -- "$1")
# Remove filename.engrad if it exists
rm -f $filename.engrad
rm -f vibspectrum
exec=$5

# Set default parallel number to 4 if not provided
parallel=${6:-4}

# Check if the flag output is enabled
if [ "$3" == "--output" ]; then
    if [ "$4" == "--fe" ]; then
        $5 $1 --sp --grad --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000
    fi
    if [ "$4" == "--vib" ]; then
        $5 $1 --hess --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000
    fi
    # if [ "$4" == "--nci" ]; then
    # fi
else
    if [ "$4" == "--fe" ]; then
        $5 $1 --sp --grad --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000 > /dev/null 2>&1
    fi
    if [ "$4" == "--vib" ]; then
        $5 $1 --hess --gfn 2 --parallel $parallel --norestart --acc 0.0001 --iterations 1000 > /dev/null 2>&1
    fi
    # if [ "$4" == "--nci" ]; then
    # fi
fi

