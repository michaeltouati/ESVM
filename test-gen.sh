#!/usr/bin/env bash

cp test-cases/Test/Ampere/input-deck . && ./esvm > output && cp output test-cases/Test/Ampere/ 
cp test-cases/Test/Poisson/input-deck . && ./esvm > output && cp output test-cases/Test/Poisson/ 
cp test-cases/Test/OpenMP/input-deck . && ./esvm > output && cp output test-cases/Test/OpenMP/ 
cp test-cases/Test/Linear-advection-schemes/Donor-cell/input-deck . && ./esvm > output && cp output test-cases/Test/Linear-advection-schemes/Donor-cell/ 
cp test-cases/Test/Linear-advection-schemes/Lax-Wendroff/input-deck . && ./esvm > output && cp output test-cases/Test/Linear-advection-schemes/Lax-Wendroff/ 
cp test-cases/Test/Linear-advection-schemes/Beam-Warming/input-deck . && ./esvm > output && cp output test-cases/Test/Linear-advection-schemes/Beam-Warming/ 
cp test-cases/Test/Linear-advection-schemes/Fromm/input-deck . && ./esvm > output && cp output test-cases/Test/Linear-advection-schemes/Fromm/ 
cp test-cases/Test/Non-linear-advection-schemes/Minmod/input-deck . && ./esvm > output && cp output test-cases/Test/Non-linear-advection-schemes/Minmod/ 
cp test-cases/Test/Non-linear-advection-schemes/Superbee/input-deck . && ./esvm > output && cp output test-cases/Test/Non-linear-advection-schemes/Superbee/ 
cp test-cases/Test/Non-linear-advection-schemes/Van-Leer/input-deck . && ./esvm > output && cp output test-cases/Test/Non-linear-advection-schemes/Van-Leer/ 
cp test-cases/Test/Non-linear-advection-schemes/MUSCL1/input-deck . && ./esvm > output && cp output test-cases/Test/Non-linear-advection-schemes/MUSCL1/ 
cp test-cases/Test/Non-linear-advection-schemes/MUSCL2/input-deck . && ./esvm > output && cp output test-cases/Test/Non-linear-advection-schemes/MUSCL2/ 