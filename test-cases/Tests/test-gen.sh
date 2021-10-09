#!/usr/bin/env bash
# Written by Dr Michaël J TOUATI
mv input-deck input-deck-old
cp test-cases/Tests/Ampere/input-deck . && ./esvm > output && cp output test-cases/Tests/Ampere/ 
cp test-cases/Tests/Poisson/input-deck . && ./esvm > output && cp output test-cases/Tests/Poisson/ 
cp test-cases/Tests/OpenMP/input-deck . && ./esvm > output && cp output test-cases/Tests/OpenMP/ 
cp test-cases/Tests/Boundary-conditions/Periodic/input-deck . && ./esvm > output && cp output test-cases/Tests/Boundary-conditions/Periodic/ 
cp test-cases/Tests/Boundary-conditions/Absorbing/input-deck . && ./esvm > output && cp output test-cases/Tests/Boundary-conditions/Absorbing/
cp test-cases/Tests/Linear-advection-schemes/Donor-cell/input-deck . && ./esvm > output && cp output test-cases/Tests/Linear-advection-schemes/Donor-cell/ 
cp test-cases/Tests/Linear-advection-schemes/Lax-Wendroff/input-deck . && ./esvm > output && cp output test-cases/Tests/Linear-advection-schemes/Lax-Wendroff/ 
cp test-cases/Tests/Linear-advection-schemes/Beam-Warming/input-deck . && ./esvm > output && cp output test-cases/Tests/Linear-advection-schemes/Beam-Warming/ 
cp test-cases/Tests/Linear-advection-schemes/Fromm/input-deck . && ./esvm > output && cp output test-cases/Tests/Linear-advection-schemes/Fromm/ 
cp test-cases/Tests/Non-linear-advection-schemes/Minmod/input-deck . && ./esvm > output && cp output test-cases/Tests/Non-linear-advection-schemes/Minmod/ 
cp test-cases/Tests/Non-linear-advection-schemes/Superbee/input-deck . && ./esvm > output && cp output test-cases/Tests/Non-linear-advection-schemes/Superbee/ 
cp test-cases/Tests/Non-linear-advection-schemes/Van-Leer/input-deck . && ./esvm > output && cp output test-cases/Tests/Non-linear-advection-schemes/Van-Leer/ 
cp test-cases/Tests/Non-linear-advection-schemes/MUSCL1/input-deck . && ./esvm > output && cp output test-cases/Tests/Non-linear-advection-schemes/MUSCL1/ 
cp test-cases/Tests/Non-linear-advection-schemes/MUSCL2/input-deck . && ./esvm > output && cp output test-cases/Tests/Non-linear-advection-schemes/MUSCL2/ 
cp test-cases/Tests/Landau/input-deck . && ./esvm > output && cp output test-cases/Tests/Landau/ 
cp test-cases/Tests/Wakefield/input-deck . && ./esvm > output && cp output test-cases/Tests/Wakefield/ 
cp test-cases/Tests/Two-stream-instability/input-deck . && ./esvm > output && cp output test-cases/Tests/Two-stream-instability/ 
rm -f output
rm -f input-deck
mv input-deck-old input-deck