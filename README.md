# rlpbwt

## Type of RLPBWT
In this repository we have 4 possible versions of RLPBWT:
- `rlpbwt_naive`
- `rlpbwt_bv` (practically the previous one with bitvectors for runs and u/v 
   values)
- `rlpbwt_ra` with the panel as bitvectors
- `rlpbwt_ra` with the panel as an SLP

The first two can be queried to obtain the number of rows that are matching at 
column k while the last two (if extended with the phi data structure) could also
report the exact rows that are matching.

## Build
Normal build:
```
cmake -S . -B build -D BUILD_TESTS=OFF
cmake --build build
cd build
./rlpbwt
```

test:
```
cmake -S . -B build_test -D BUILD_TESTS=ON
cmake --build build_test
cd build_test
./rlpbwt_test
```
