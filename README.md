# rlpbwt
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
