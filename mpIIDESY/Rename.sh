# Rename all *.root to *.bin 
for f in $1*.root; do 
    mv -- "$f" "${f%.root}.bin"
done