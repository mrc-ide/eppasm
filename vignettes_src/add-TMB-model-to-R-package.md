# Add TMB model to R package with a `Makefile`

Steps to add TMB model `tmbmodel` to R package `rpackage`

## From `tmbmodel`'s side

Write and test `tmbmodel.cpp` as usual

## From `rpackage`'s side

The `rpackage`'s directory tree

```
├── DESCRIPTION
├── NAMESPACE
├── R
│   └── kkk.R
├── src
│   ├── *.hpp
│   ├── Makefile
│   ├── *.cpp
├── tmb
    ├── tmbmodel.cpp
```

1. Add to `DESCRIPTION` TMB and dependencies

```
...
Imports: 
    TMB,
    RcppEigen
    ...
```

2. To load compiled TMB DLLs into `rpackage`'s name space, add TMB DLL names to `NAMESPACE`, e.g.,

```
...
useDynLib(tmbmodel)
...
```

or using `roxygen2`'s approach (see `?roxygenise`).

3. Finally, add to `Makefile`

```
all: buildr tmbi

buildr:
    R CMD COMPILE *.cpp *.c;\
    R CMD SHLIB -o rpackage.so *.o
tmbi:
    cd ../tmb;\
    Rscript --vanilla -e "TMB::compile('tmbmodel.cpp')";\
    mv *.so ../src/;\
clean:
    rm *.o; rm ../tmb/*.o
```

If we have other TMB models, add `.cpp` files to `tmb`, add their names to `NAMESPACE` and `Makefile` (section `tmbi`).