# Rename all *.root to *.bin 
for f in $1*.root; do 
    echo "$f"
    echo "${f%.root}.bin"
    mv -f "$f" "${f%.root}.bin"
done