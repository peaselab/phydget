---

## phydget
***PhyDGET: Phylogenetic Differential Gene Expression Tool<br>\linebreak Author: James B. Pease***

**Parameters**


``--data`` (required) = input expression data filepath (csv format) (type=file path, default=None)

``--out`` (required) = output file path (csv format) (type=file path, default=None)

``--tree`` (required) = input tree file path (Nexus format) (type=file path, default=None)

``--bt-burnin/--btburnin`` = BayesTrait number of burn-in steps. (type=integer, default=1000000)

``--bt-exec/--btexec`` = BayesTrait executable path (type=None, default=BayesTraitV3)

``--bt-iter/--btiter`` = BayesTrait number of iterations usedper stone in the stepping stone sampling. (type=integer, default=10000000)

``--bt-priors-alpha/--btpriorsalpha`` = BayesTrait distribution type and prior range for alpha. (type=None, default=('uniform', -10, 30))

``--bt-priors-sigma/--btpriorssigma`` = BayesTrait distribution type and and prior range for sigma^2. (type=None, default=('uniform', 0, 60))

``--bt-priors-vrbl/--btpriorsvrbl`` = BayesTrait distribution and prior range for variable rates branch length differential. This is advanced and requires BayesTraitV4. (type=None, default=None)

``--bt-stoneiter/--btstoneiter`` = BayesTrait number of iterations usedper stone in the stepping stone sampling. (type=integer, default=20000)

``--bt-stones/--bt-stones`` = BayesTrait number of stones usedin the stepping stone sampling. (type=integer, default=200)


``--keep-files/--keepfiles`` = Keep all temporary files (flag, default=False)


``--temp-dir/--tempdir`` = temporary folder for files (type=None, default=PhyDGETtmp)

``--temp-prefix/--temp-prefix`` = Temporary Directory Prefix (type=None, default=PhyDGETRun)

``--test-gene/--testgene`` = Enter exacty gene name from first column ofinput csv file to do a test run on a single gene. (type=None, default=None)

``--threads`` = Number of threads for parallelization (type=integer, default=2)

``--tip-values/--tipvalues`` = Values to place at the tips (see manual for details). (type=None, default=all)
Choices: ('all', 'amean', 'hmean', 'gmean', 'median', 'middle')

``--transform`` = Data transformation type (see manual for details). (type=None, default=log2cpm)
Choices: ('none', 'log2', 'cpm', 'log2cpm')


``--verbose`` = extra screen output (flag, default=False)


