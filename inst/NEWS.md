## Version 1.11.1

---

### CHANGES

- Parallel reading is added to readBismark-method.

### BUG FIX

- Some bugs are fixed in methHMMCMC-method to avoid creating infinity.

## Version 1.9.2

---

### CHANGES SINCE LAST RELEASE

- Made compatible with R 4 following changes in some codes in R.
- Citation is updated.


## Version 1.3.1

---

### CHANGES AND BUG FIX

- A new citation is added.
- Several bugs are fixed.


## Version 0.99.17

---

### CHANGES SINCE LAST RELEASE

- A new function, compareDMCs() is added for plotting DMCs and average predicted
    methylation for for group comparison.
- The manhattanDMCs() function is modified to plot regions of interests.


## Version 0.99.16

---

### CHANGES SINCE LAST RELEASE

- The function findDMCs() is fixed for memory management.
- The function methEM() is fixed for creating NaN in computing binomial
    probabilities.


## Version 0.99.15

---

### CHANGES SINCE LAST RELEASE

- An option is added to findDMC that lets the user enters a function to
    calculate weights for GLM using the variances obtianed from MCMC.


## Version 0.99.10 - 0.99.14

---

### CHANGES SINCE LAST RELEASE

- Several issues are fixed.
- Some modification are done to expedite the speed of processes.
- Instead of the parallel package, the BioParallel package is used.


## Version 0.99.9

---

### CHANGES SINCE LAST RELEASE

- Some parameters in examples are changed to make the running time faster.
- The BSData-method is changed to cBSData-method.
- The BSDMCs-method is changed to cBSDMCs-method.


## Version 0.99.0

---

### CHANGES SINCE LAST RELEASE

- We have refined some of the codes to speed up the running the package.
