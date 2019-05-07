if [ ! -f FileList.txt ]; then
    echo "FileList.txt not found"
    return
fi

if [ -z $GM2_VERSION ]; then
    echo "GM2 software not setup"
    return 
fi

for file in `cat FileList.txt`; do echo $file; done | xargs -i --max-procs 8 bash -c ". runJob2.sh {}"