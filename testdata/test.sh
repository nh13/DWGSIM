#!/bin/bash

set -e

TESTDATADIR=$(dirname ${0});

# Decompress the test data
pushd $TESTDATADIR;
for FILE in $(ls -1 *gz)
do 
	gunzip -c ${FILE} > $(basename ${FILE} .gz);
done
popd

mkdir tmp

# Generate the new test data
./dwgsim -z 13 -N 10000 samtools/examples/ex1.fa tmp/ex1.test

# Test the differences
for GZFILE in $(ls -1 tmp/ex1.test*gz)
do 
	gunzip $GZFILE;
    FILE=$(basename $GZFILE .gz);
	diff -q tmp/${FILE} ${TESTDATADIR}/${FILE}
done

# Clean up the testdata
find ${TESTDATADIR} \! -name "*gz" -type f | grep -v sh$ | xargs rm; 

rm -r tmp
