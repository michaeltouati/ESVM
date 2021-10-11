#!/usr/bin/env bash
#######################################################################
##                                                                   ##
##             ElectroStatic Vlasov-Maxwell (ESVM) code              ##
##                                                                   ##
## Copyright © 2015 Michaël J TOUATI                                 ##
##                                                                   ##
## This file is part of ESVM.                                        ##
##                                                                   ##
## ESVM is free software: you can redistribute it and/or modify      ##
## it under the terms of the GNU General Public License as published ##
## by the Free Software Foundation, either version 3 of the License, ##
## or (at your option) any later version.                            ##
##                                                                   ##
## ESVM is distributed in the hope that it will be useful,           ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of    ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
## GNU General Public License for more details.                      ##
##                                                                   ##
## You should have received a copy of the GNU General Public License ##
## along with ESVM. If not, see <https://www.gnu.org/licenses/>.     ##
##                                                                   ##
#######################################################################
## Initial commit written by Michaël J TOUATI - Dec. 2015

TEST_DIR=' Maxwell/Ampere'
TEST_DIR+=' Maxwell/Poisson'
TEST_DIR+=' Parallelization/OpenMP'
TEST_DIR+=' Boundary-conditions/Periodic'
TEST_DIR+=' Boundary-conditions/Absorbing'
TEST_DIR+=' Vlasov-linear/Donor-cell'
TEST_DIR+=' Vlasov-linear/Lax-Wendroff'
TEST_DIR+=' Vlasov-linear/Beam-Warming'
TEST_DIR+=' Vlasov-linear/Fromm'
TEST_DIR+=' Vlasov-nonlinear/Minmod'
TEST_DIR+=' Vlasov-nonlinear/Superbee'
TEST_DIR+=' Vlasov-nonlinear/Van-Leer'
TEST_DIR+=' Vlasov-nonlinear/MUSCL1'
TEST_DIR+=' Vlasov-nonlinear/MUSCL2'
TEST_DIR+=' Academic-cases/Landau'
TEST_DIR+=' Academic-cases/Wakefield'
TEST_DIR+=' Academic-cases/Two-stream-insta'
# TEST_DIR+=' New-features/'

for tst in ${TEST_DIR}; do \
	cp ${tst}/input-deck . && ../.././esvm > output && cp output ${tst}/ ; \
	echo ${tst}'/output file generated'
done
rm -rf results/ output input-deck
