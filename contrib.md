# Contributing to ESVM

First of all, thank you very much for your interest in contributing to ESVM! 

Hans Bethe once said : “We need science education to produce scientists, but we need it equally to create literacy in the public. Man has a fundamental urge to comprehend the world about him, and science gives today the only world picture which we can consider as valid. It gives an understanding of the inside of the atom and of the whole universe, or the peculiar properties of the chemical substances and of the manner in which genes duplicate in biology. An educated layman can, of course, not contribute to science, but can enjoy and participate in many scientific discoveries which as constantly made. Such participation was quite common in the 19th century, but has unhappily declined. Literacy in science will enrich a person's life.”

But Richard P. Feynmann also added : “Physics is like sex: sure, it may give some practical results, but that's not why we do it.”

Being fully in agreement with Richard P. feynman while following the advice of Hans Bethe, it seemed to me natural to develop ESVM with Github in order to try to understand the complexity of plasmas, numerical analysis and High Parallel Computing (HPC) and to make it accessible to the largest public as possible at the same time.

Here can be found how to :
- Report or fix a bug
- Propose or submit a new feature

# Report a bug using Github's [issues](https://github.com/michaeltouati/ESVM/issues)

For problems related with plotting, please make sure first that LaTeX, dvipng and ghostscript are each working and on your PATH for the Matplotlib Python package to be able to render tex fonts.
If the bug persists or is related to another problem, follow these steps :
1) Go to 'Issues' on the ESVM main repository branch
2) Click on 'New issue'
4) Describe the bug the more clear and concise as possible in the title starting the title starting with "Bug :"
5) Describe with the more details as possible the bug and
6) Click on 'Submit new issue'

OS: [e.g. Ubuntu 20,04]

Fortran compiler [e.g. gfortran 11.2.0]

Python version [e.g. Python 3.7.11]

Matplotlib Python package version [e.g. 3.4.3]

Numpy Python package version [e.g., 1.21.2]

# Propose a new feature using Github's [issues](https://github.com/michaeltouati/ESVM/issues)

If you want to propose a new feature for ESVM, follow these steps :
1) Go to 'Issues' on the ESVM main repository branch
2) Click on 'New issue'
4) Describe the new feature request the more clear and concise as possible in the title starting "Feature request :"
5) Describe with the more details as possible the new feature request and
6) Click on 'Submit new issue'

# Fix a bug or submit a new feature

ESVM Uses the Fork and pull model [Github Flow](https://guides.github.com/introduction/flow/index.html).
Pull requests are the best way to propose changes to the codebase (we use [Github Flow](https://guides.github.com/introduction/flow/index.html)). We actively welcome your pull requests:

1. Fork the repo and create your branch from `master`.
2. If you've added code that should be tested, add tests.

Please, try to keep the code Fortran 90 standard compliant by : 
- respecting 2 spaces for indentation rather than tabs
- not using object oriented Fortran 2003 features
- etc ...

3. If you've added parameters in the input-deck, please add their description in the input-deck following the same style.
4. Ensure the code compiles in debug mode 
5. Ensure the code pass the provided tests and add a test if a new solver, a new boundary condition, etc... is added.
6. Issue that pull request!

# Any contributions you make will be under the GNU General Public License v3.0

When you submit code changes, your submissions are understood to be under the same [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) that covers ESVM. 

# License
By contributing, you agree that your contributions will be licensed under its GNU General Public License v3.0.
