# Rational function reconstruction benchmarks

This repository contains benchmarks to compare different algorithms
for rational function reconstruction. See
[data/Readme.md](data/Readme.md) for more information on the sample
rational functions.

## How to run

Note that a run will typically take several minutes.

### [rare](1)

If [Rust and Cargo](https://www.rust-lang.org/) are installed on your
system, run

    cargo run --release

### [FireFly](2)

First, install [FireFly](2). Then compile `FireFly/rec.cpp`, e.g. with

```
g++ -o rec FireFly/rec.cpp -lfirefly -lgmp -lgmpxx -lflint
```

and run the resulting executable.

### [FiniteFlow](3)

For this, Wolfram Mathematica is required. After installing
[FiniteFlow](3) run

```
FiniteFlow/rec.wls
```

[1]: https://github.com/a-maier/rare
[2]: https://github.com/jklappert/FireFly/
[3]: https://github.com/peraro/finiteflow
