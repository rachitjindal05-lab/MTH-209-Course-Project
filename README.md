# MTH-209-Course-Project
**gradR: Theoretical Gradient Descent and Optimization Diagnostic Dashboard**

**Description:**

I built a package called **gradR** along with a 1 minute video reel explaining the concept and runme.R file to run and test the package.

The **gradR** package is an interactive, educational tool designed for visualizing mathematical optimization. Built on a Shiny dashboard, it allows users to dynamically explore first-order gradient descent paths and evaluate second-order Hessian conditions for saddle point diagnostics. Furthermore, gradR allows users to simulate Stochastic Gradient Descent (SGD) with dynamic learning rate decay, providing a robust platform to analyze how mathematical optimization algorithms navigate complex loss landscapes.

## Prerequisites

Before installing `gradR`, ensure you have the required dependency packages installed from CRAN. You can install all of them at once by running the following command in your R console:
To install the package directly from GitHub, you will first need to have the `devtools` package installed. You can install it by running the following command in your R console:

```R
# Install required dependencies
install.packages(c("dplyr", "ggplot2", "plotly", "shiny"))
# Install devtools if you haven't already
install.packages("devtools")
```

## Getting Started

You can install it directly from the .tar.gz file once the prerequisites are met.

### Installation
Since the package is currently distributed as a source bundle, you can install it directly from the `.tar.gz` file. 

```R
# Install the package directly from the bundled file
install.packages("/path/to/gradR_1.0.0.tar.gz", repos = NULL, type = "source")

# Load the library
library(gradR)

# Launch the interactive diagnostic dashboard
visualize_gradientdescent()
```
