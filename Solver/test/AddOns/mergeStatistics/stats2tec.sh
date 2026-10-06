#!/bin/bash
#
#  Converts a statistics file to Tecplot with horses2plt (control file mode),
#  so that it can be compared with compare_tec.py.
#
#  Usage: ./stats2tec.sh MESH.hmesh FILE.stats.hsol [--has-gradients]
#
#  The output is FILE.stats.tec, next to the statistics file.
#
mesh=$1
stats=$2
variables="Q,Vmean,Sij"
hasGradients=.false.
if [ "$3" == "--has-gradients" ]; then
   variables="$variables,gradV"
   hasGradients=.true.
fi

[ -s "$stats" ] || { echo "FAILED: $stats does not exist or is empty"; exit 1; }

convertFile=$(mktemp ./stats2tec_XXXXXX.convert)
cat > $convertFile << EOF
hmesh file       = $mesh
hsol file        = $stats
has gradients    = $hasGradients
output variables = $variables
EOF

rm -f "${stats%.hsol}.tec"
./horses2plt $convertFile > /dev/null
rc=$?
rm -f $convertFile

[ $rc -eq 0 ] || { echo "FAILED: horses2plt exited with code $rc for $stats"; exit 1; }
[ -s "${stats%.hsol}.tec" ] || { echo "FAILED: ${stats%.hsol}.tec was not generated"; exit 1; }
echo "OK: ${stats%.hsol}.tec"
