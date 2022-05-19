# HDX-MD
Gromacs plugin which helps you to calculate distance of three nearest water molecules from an amide nitrogen. This allows identification of MD frames which can go HDX as per the criterion proposed by Halle and Persson https://www.pnas.org/doi/10.1073/pnas.1506079112

This plugin is for Gromacs 2019.4. The compilation steps are the same as for standard gromacs:

```
cd gromacs-2019.4
tar xvf  mindist2_patch.tar
mkdir build
cd build
cmake -DGMX_BUILD_OWN_FFTW=ON ..
make -j4
```

Usage:

1. Don't forget to save waters.
2. Try ``` mindist2 -h ```

Contact: Pär Söderhjelm, https://www.cmps.lu.se/bpc/people/par-soderhjelm/
