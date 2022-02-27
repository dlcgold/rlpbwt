# rlpbwt

## Type of RLPBWT
In this repository we have 4 possible versions of RLPBWT:
- `rlpbwt_naive`
- `rlpbwt_bv` (practically the previous one with bitvectors for runs and u/v 
   values)
- `rlpbwt_ms` with the panel as bitvectors
- `rlpbwt_ms` with the panel as an SLP

The first two can be queried to obtain the number of rows that are matching at 
column k while the last two (if extended with the phi data structure) could also
report the exact rows that are matching.

## Build
to build:
```
cmake -S . -B build -D BUILD_TESTS=OFF
cmake --build build
cd build
./rlpbwt
```
usage:
```Usage: RLPBWT [options]

Options:
  -i, --input_file <path>
  -s, --input_slp <path>
  -o, --output <path>
  -q, --query <path>
  -N, --Naive: naive RLPBWT (only one mode allowed)
  -B, --Bv: bitvectors RLPBWT (only one mode allowed)
  -S, --Slp: slp RLPBWT (only one mode allowed | slp file required)
  -P, --Panel: panel RLPBWT (only one mode allowed)
  -t, --thresholds: enable thresholds (slp/panel mode only)
  -e, --extend: extend matches (slp/panel mode only)
  -v, --verbose: extra prints
  -h, --help: show this help message and exit
```
example using RLPBWT with slp with lce queries and matches extended:
```./rlpbwt -i <input matrix> -s <slp file> -o <output file> -q <query file> -S -e```
## Test
to build tests:
```
cmake -S . -B build_test -D BUILD_TESTS=ON
cmake --build build_test
cd build_test
./rlpbwt_test
```
