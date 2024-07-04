#!/bin/bash
testdir=`pwd`

# executable
var=test_siem_mesh

# title
echo >> $testdir/results.log
echo "test: $var" >> $testdir/results.log
echo >> $testdir/results.log

echo "directory: `pwd`" >> $testdir/results.log

# clean
mkdir -p bin
rm -f ./bin/$var

# single compilation
echo "compilation: $var" >> $testdir/results.log

make -j8 xmeshfem3D >> $testdir/results.log 2>&1

echo "" >> $testdir/results.log

# check
if [ ! -e ./bin/xmeshfem3D ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
else
    echo "WARNING: TEST IS SIMPLY BUILDING MESHFEM AT THE MOMENT"
fi

# runs test

if [ ! -d  $testdir/DATABASES_MPI ]; then
  mkdir DATABASES_MPI
fi


echo "run: `date`" >> $testdir/results.log
mpirun -np 24 ./bin/xmeshfem3D >> $testdir/results.log 2>$testdir/error.log

# checks exit code
if [[ $? -ne 0 ]]; then
  echo "test failed"; echo "error log:"; cat $testdir/error.log; echo ""
  exit 1
fi

# checks error output (note: fortran stop returns with a zero-exit code)
if [[ -s $testdir/error.log ]]; then
  echo "returned ERROR output:" >> $testdir/results.log
  cat $testdir/error.log >> $testdir/results.log
  exit 1
fi
rm -f $testdir/error.log

#cleanup
rm -f bin/$var
# done
echo "successfully tested: `date`" >> $testdir/results.log
