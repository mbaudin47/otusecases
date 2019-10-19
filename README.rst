otusecases
==========

A repository of OpenTURNS test cases

What is it?
-----------

This project contains OpenTURNS use cases and datasets.

- The datasets can be used to perform data analysis with OpenTURNS. 
The goal is to have a set of .csv files which can be easily imported 
with OpenTURNS for example to fit a marginal distribution, 
to fit a copula or to create a datamodel. 

- Each use-case is made of a function and the distribution of its 
inputs. 
Each use-case has a specific methodological goal: 
central dispersion, reliability, sensitivity analysis or 
calibration. 
In general, each use-case is presented in a Notebook: 
the equations of the function, the references (if any), 
the simplest possible study which shows how to use the test-case. 
The use-cases can be used to benchmark a method, but the 
Notebook only shows the simplest (e.g. the Monte-Carlo method): 
the goal of each example is *not* to actuall show the benchmark 
results. 
Generally, each use case has a specific goal (e.g. central 
dispersion), but some use cases can be derived to achieve different  
methodological goals. 

It is based on `OpenTURNS <http://www.openturns.org>`_.

The scripts are tested with OT 1.13.

Overview
--------

* Axial stressed beam (reliability)
* Cantilever beam (reliability)
* Perrin case (sensitivity analyis)
* Chaboche model (calibration)
* Viscous vertical fall (calibration)
* Flooding case (central dispersion, calibration)
* R-S case (reliability)
* Deflection of a tube (reliability)
* G-Sobol' function (sensitivity analysis)
* Ishigami function (sensitivity analysis)
* Logistic model (calibration, stochastic process)
* Morris function (sensitivity analysis)
* Nonlinear oscillator (reliability)
* Product function (sensitivity analysis)

How to install?
---------------

Requirements
~~~~~~~~~~~~

The dependencies are: 

- Python >= 2.7 or >= 3.3
- `OpenTURNS <http://www.openturns.org>`_ >= 1.13
- `Jupyter Notebook <https://jupyter.org>`_

Installation
~~~~~~~~~~~~

Using the latest python version is prefered! Then to install::

    git clone https://github.com/mbaudin47/otusecases.git
    cd otusecases
