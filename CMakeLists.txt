cmake_minimum_required(VERSION 2.8)
PROJECT(r-stochbb)

SET(doc_SOURCES man/Container-class.Rd man/independent.Rd man/Density-class.Rd man/invgamma.Rd
  man/ExactSampler-class.Rd man/kolmogorov.Rd man/Var-class.Rd man/logLikelihood.Rd man/affine.Rd
  man/maximum.Rd man/atomic.Rd man/minimum.Rd man/chain.Rd man/mixture.Rd man/compound.Rd
  man/normal.Rd man/condchain.Rd man/stochbb-package.Rd man/conditional.Rd man/studt.Rd
  man/delta.Rd man/uniform.Rd man/derived.Rd man/weibull.Rd man/gamma.Rd)
SET(mod_SOURCES src/Makevars src/module.cc)
SET(R_SOURCES R/stochbb.R)
SET(meta_SOURCES DESCRIPTION NAMESPACE)

add_custom_target(doc SOURCES ${doc_SOURCES})
add_custom_target(mod SOURCES ${mod_SOURCES})
add_custom_target(R SOURCES ${R_SOURCES})
add_custom_target(meta SOURCES ${meta_SOURCES})

