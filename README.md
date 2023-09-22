# rlpbwt

Other run-length encoded PBWT implementations can be found at [koeppl/pbwt](https://github.com/koeppl/pbwt).
## Type of RLPBWT

In this repository we have various  possible versions of RLPBWT:

- `rlpbwt_naive`
- `rlpbwt_bv` (practically the previous one with bitvectors for runs and u/v
  values)
- `rlpbwt_ms` with the panel as bitvectors
- `rlpbwt_ms` with the panel as an SLP (for random access or LCE query)

The first two can be queried to obtain the number of rows that are matching at
column k while the last two (if extended with the phi data structure) could also
report the exact rows that are matching.

## Build

to build:

```
cmake -S . -B build -D BUILD_TESTS=OFF
cmake --build build
```
to execute:
```
cd build
./rlpbwr [options]
```
usage:

```
Usage: RLPBWT [options]

Options:
  -i, --input_file <path>	 macs file for panel
  -m, --memorize <path>	  path to save serialization 
  -l, --load <path>	 path to load serialization
  -s, --input_slp <path>	 path to SLP file
  -o, --output <path>	 path to query output
  -q, --query <path>	 path to macs query file
  -N, --Naive	 naive RLPBWT (only one mode allowed)
  -B, --Bv	 bitvectors RLPBWT (only one mode allowed)
  -S, --Slp	 slp RLPBWT (only one mode allowed | slp file required)
  -P, --Panel	 panel RLPBWT (only one mode allowed)
  -t, --thresholds	 enable thresholds (slp/panel mode only)
  -e, --extend	 extend matches (slp/panel mode only)
  -r, --raw	use compressed intvectors instead of sparse bitvectors
  -v, --verbose	 extra prints
  -V, --fverbose	 extra prints for functions (cautions)
  -h, --help	 show this help message and exit

```

example using RLPBWT (slp mode) with lce queries and matches extended:

```./rlpbwt -i <input matrix> -s <slp file> -o <output file> -q <query file> -S -e```

## Test

to build tests:

```
cmake -S . -B build_test -D BUILD_TESTS=ON
cmake --build build_test
```

to execute:
```
mkdir output
cd build_test
./rlpbwt_test
```
