#!/usr/bin/env python3
#
#  Builds a legacy statistics file from a current one, to test
#  "legacy stats = .true." in horses2plt.
#
#  Current .stats.hsol files store, for each element, the 9 statistics
#  (u, v, w, uu, vv, ww, uv, uw, vw) followed by the averaged conserved
#  variables (and their gradients, if saved). Legacy files only store the
#  9 statistics. The header (file name, type, node type, number of
#  elements, iteration, time, reference values) is the same.
#
#  Usage: ./make_legacy_stats.py CURRENT.stats.hsol LEGACY.stats.hsol
#
import os
import struct
import sys

SOLFILE_STR_LEN = 128
SIZEOF_INT = 4
SIZEOF_RP = 8
NO_OF_SAVED_REFS = 6
NO_OF_STATS = 9

# Byte offset of the first element (POS_INIT_DATA - 1 in SolutionFile.f90)
HEADER_SIZE = SOLFILE_STR_LEN + 4 * SIZEOF_INT + SIZEOF_RP + NO_OF_SAVED_REFS * SIZEOF_RP + SIZEOF_INT
POS_NOOFELEMENTS = SOLFILE_STR_LEN + 2 * SIZEOF_INT

inFile, outFile = sys.argv[1], sys.argv[2]

with open(inFile, "rb") as f:
    data = f.read()

header = data[:HEADER_SIZE]
noOfElements = struct.unpack_from("<i", data, POS_NOOFELEMENTS)[0]

# Each element: rank, 4 dimensions, then (9 + extra) values per point. The number of
# extra values per point is the same for all elements; find it by walking the file.
# The file ends with a 4-byte trailer, which is copied as is.
pos = HEADER_SIZE
dims = struct.unpack_from("<5i", data, pos)
if dims[0] != 4 or dims[1] != NO_OF_STATS:
    sys.exit(f"FAILED: unexpected array header {dims} in the first element")

# Walk the file for every possible number of values per point until it ends at the trailer
for valuesPerPoint in range(NO_OF_STATS, 200):
    pos, blocks = HEADER_SIZE, []
    ok = True
    for eID in range(noOfElements):
        if pos + 5 * SIZEOF_INT > len(data):
            ok = False
            break
        rank, n1, n2, n3, n4 = struct.unpack_from("<5i", data, pos)
        if rank != 4 or n1 != NO_OF_STATS:
            ok = False
            break
        npoints = n2 * n3 * n4
        blocks.append((pos, npoints))
        pos += 5 * SIZEOF_INT + valuesPerPoint * npoints * SIZEOF_RP
    if ok and len(data) - pos == SIZEOF_INT:
        break
else:
    sys.exit(f"FAILED: could not find the element layout of {inFile}")

with open(outFile, "wb") as f:
    f.write(header)
    for start, npoints in blocks:
        f.write(data[start:start + 5 * SIZEOF_INT + NO_OF_STATS * npoints * SIZEOF_RP])
    f.write(data[pos:])

print(f"OK: {outFile} ({noOfElements} elements, dropped {valuesPerPoint - NO_OF_STATS} "
      f"values per point, {os.path.getsize(outFile)} bytes)")
