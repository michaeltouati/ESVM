# Contribute to ESVM 

1) download ESVM from https://github.com/michaeltouati/ESVM or type on your terminal : 

git clone https://github.com/michaeltouati/ESVM.git

2) If it's been a long time you did the step 1., make sure you have the latest ESVM updates by typing : 

git pull

3) Once you've have modified the ESVM source file(s) to fix a bug or propose a new feature, type from your local ESVM directory : 

make clean_all 

4) Type then : 

git add . 

if you've added new files

5) Type : 

git commit -m "MY COMMIT MESSAGE" -m -a 

describing the more clear and concise as possible the goal of the contribution in "MY COMMIT MESSAGE"

6) Ask for the merge of your local ESVM branch to the main branch by typing : 

git push

# Report issues or problems with ESVM 

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

# Seek support to add a new feature in ESVM

1) Go to 'Issues' on the ESVM main repository branch
2) Click on 'New issue'
4) Describe the new feature request the more clear and concise as possible in the title starting "Feature request :"
5) Describe with the more details as possible the new feature request and
6) Click on 'Submit new issue'
