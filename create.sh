#!/bin/bash
CWD=$PWD
module --force purge
module load StdEnv/2020  gcc/9.3.0 openblas/0.3.17 gsl sqlite/3.36 llvm/9.0.1 r-bundle-bioconductor/3.12 python/3
module save epiclomalstack
virtualenv --no-download ec-env && source ec-env/bin/activate
sed -i '3 i module restore epiclomalstack' "${VIRTUAL_ENV}/bin/activate"
mkdir -p "${VIRTUAL_ENV}/R"
cd "${VIRTUAL_ENV}"
git clone git@github.com:molonc/Epiclomal.git
sed -i "s/datrie=0.8/datrie=0.8.2/g" Epiclomal/conda_packages.txt
sed -i "s/llvmlite=0.29.0/llvmlite=0.32.1/g" Epiclomal/conda_packages.txt
sed -i "s/pandas=0.25.0/pandas=0.25.3/g" Epiclomal/conda_packages.txt
wget https://gitlab.freedesktop.org/pixman/pixman/-/archive/pixman-0.38.0/pixman-pixman-0.38.0.tar.bz2
tar xjf pixman-pixman-0.38.0.tar.bz2
cd pixman-pixman-0.38.0
bash autogen.sh
./configure --prefix="${VIRTUAL_ENV}"
make
make install
cd "${VIRTUAL_ENV}" && rm -rf pixman-pixman-0.38.0
cd "${CWD}"
grep '=py' "${VIRTUAL_ENV}/Epiclomal/conda_packages.txt"| cut -d= -f 1-2| sed 's/=/==/g' | while IFS= read -r package
do
  if [[ $(avail_wheels "${package}" | wc -l) -gt 2 ]]
  then
    pip install --no-index "${package}"
  else
    pip install "${package}"
  fi
done
pip install --no-index joblib==0.14.0 setuptools setuptools_scm --upgrade
pip install snakemake==5.5.4
export R_LIBS_USER=${VIRTUAL_ENV}/R
sed -i "4 i export R_LIBS_USER=${VIRTUAL_ENV}/R" "${VIRTUAL_ENV}/bin/activate"
RP_URL="https://cran.r-project.org/src/contrib/Archive"
grep '^r-' "${VIRTUAL_ENV}/Epiclomal/conda_packages.txt" | cut -d- -f2|cut -d= -f1-2|sed 's/=/_/g' | while IFS= read -r Rpackage
do
  URL="${RP_URL}/${Rpackage%%_*}/${Rpackage}.tar.gz"
  wget -P "${VIRTUAL_ENV}/R" "${URL}"
  R CMD INSTALL -l "${VIRTUAL_ENV}/R" "${VIRTUAL_ENV}/R/${Rpackage}.tar.gz"
done
rm -rf "${VIRTUAL_ENV}/R/*.tar.gz"
cd "${VIRTUAL_ENV}/Epiclomal"
python setup.py install
git clone https://bitbucket.org/jerry00/densitycut_dev.git
R CMD build densitycut_dev/
R CMD INSTALL -l "${VIRTUAL_ENV}/R" densitycut_0.0.1.tar.gz
wget -P "${VIRTUAL_ENV}/R" https://cran.r-project.org/src/contrib/Archive/bigstatsr/bigstatsr_0.2.1.tar.gz
R CMD INSTALL -l "${VIRTUAL_ENV}/R" "${VIRTUAL_ENV}/R/bigstatsr_0.2.1.tar.gz"
wget -P "${VIRTUAL_ENV}/R" "https://cran.r-project.org/src/contrib/Archive/NbClust/NbClust_2.0.4.tar.gz"
R CMD INSTALL -l "${VIRTUAL_ENV}/R" "${VIRTUAL_ENV}/R/NbClust_2.0.4.tar.gz"
R CMD build REpiclomal
R CMD INSTALL -l "${VIRTUAL_ENV}/R" REpiclomal_1.0.tar.gz
cd "$CWD"