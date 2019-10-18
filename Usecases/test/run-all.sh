#!/bin/sh

test_python_script()
{
  # test_python_script crue-2vars-symbolic.py
  pythonscript=$1
  cp $pythonscript /tmp
  python /tmp/$pythonscript
}

test_ipython_notebook()
{
  # test_ipython_notebook Logistique-calage.ipynb
  ipythonnotebook=$1
  basefilename=$(basename -- "$ipythonnotebook")
  basefilename="${basefilename%.*}"
  nbfile="$basefilename.ipynb"
  pyfile="$basefilename.py"
  jupyter nbconvert --to python $nbfile
  mv $pyfile /tmp
  python /tmp/$pyfile
}

filename=$(basename -- "$fullfile")

set -xe
# Run tests
cd ..
# axial-stressed-beam
cd axial-stressed-beam
test_ipython_notebook axial_stressed_beam.ipynb
cd ..
# cantilever_beam
cd cantilever-beam
test_ipython_notebook cantilever_beam.ipynb
cd ..
# cas-perrin
cd cas-perrin
test_ipython_notebook Sensibilite-Exemple-Perrin.ipynb
cd ..
# chaboche
cd chaboche
test_ipython_notebook Calibration-Chaboche.ipynb
test_python_script chaboche-genere-data.py
test_python_script chaboche-NLLS.py
cd ..
# chute-verticale
cd chute-verticale
test_ipython_notebook Chute-verticale.ipynb
test_python_script chute-verticale.py
test_python_script chute-verticale-vs-coefficient.py
cd ..
# crue-calage
cd crue-calage
test_ipython_notebook Calage-crue-OT.ipynb
test_ipython_notebook Calage-crue-OT-lineaire.ipynb
test_python_script crue-4vars-genere-data.py
cd ..
# crue-propagation
cd crue-propagation
test_python_script crue-2vars-symbolic.py
test_python_script crue-4vars-stochastic.py
test_python_script crue-8I3O-python.py
test_python_script crue-8vars-symbolic.py
test_python_script crue-propagation.py
cd ..
# fiabilite-RS
cd fiabilite-RS
test_python_script cas-RS.py
cd ..
# fleche-tube
cd fleche-tube
test_ipython_notebook Deflection-tube.ipynb
test_python_script fleche_tube-symbolic-FR.py
cd ..
# gsobol
cd gsobol
test_python_script gsobollib.py
test_python_script sensitivity-confidence-gsobol.py
test_python_script sensitivity-convergence-gsobol.py
cd ..
# ishigami
cd ishigami
test_python_script ishigami-AS.py
test_ipython_notebook La_fonction_Ishigami.ipynb
cd ..
# logistique-calage
cd logistique-calage
test_python_script logistique-calage-generate-data.py
test_ipython_notebook Logistique-calage.ipynb
cd ..
# logistique-champs
cd logistique-champs
test_python_script logistic-OT.py
test_ipython_notebook logistic-example.ipynb
cd ..
# morris
cd morris
test_ipython_notebook Morris-function.ipynb
test_ipython_notebook Morris-function-stand-alone.ipynb
cd ..
# oscillateur-nonlineaire
cd oscillateur-nonlineaire
test_ipython_notebook nonlinear-oscillator.ipynb
cd ..
# produit
cd produit
test_ipython_notebook Fonction_produit.ipynb
cd ..

