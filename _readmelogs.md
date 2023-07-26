# Installing the packages
These issues may be solved while installing the packege using the toml file in poetry.
## Scipy 
An older version of scipy was required hence the following lines where added to the toml file

scipy = { version = "^1.11.1", python = ">=3.11,<3.13" }

## Abacusutils
This package is used to acess and process the abacus summit data. While installing corffunc installation may througout an error.
Can be rectified by installing the required dependencies like gsl manually and then pip installing the same.

Also try installing from the source using github repo.
https://github.com/manodeep/Corrfunc/tree/master