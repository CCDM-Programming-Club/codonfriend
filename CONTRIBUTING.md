# Contributing to codonfriend

This package is developed as part of a regular "coding club" within our research centre.
Anyone is welcome to get involved, but maybe create an issue on GitHub to introduce yourself first!


## R packaging

This package is setup following the guidelines in Hadley Wickham's "R packages" book <http://r-pkgs.had.co.nz/>.
If you are unsure about why things are where they are or where you should put things, please consult that book first.


## Further info on R

I highly recommend that anyone with some familiarity with programming have a look at Hadley Wickham's "Advanced R" book <http://adv-r.had.co.nz/>.
It explains a lot about R's internal datastructures, what the difference between `[` and `[[` is, and an introduction to how classes work in R.


## Initial setup

1) Make sure you have [git](https://git-scm.com/) and [R (version >= 3.4.0)](https://www.r-project.org/) installed.
I would also suggest installing or updating [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop) if you don't already have another favourite IDE for R.

Windows users will need to install Rtools (<https://cran.r-project.org/bin/windows/Rtools/>) as well (not an R package, follow the link).


2) [Setup ssh access for your github account](https://help.github.com/articles/connecting-to-github-with-ssh/).
This means that you can clone, pull, and push to github without having to enter your password all of the time.
Windows users might have to [configure Rstudio to find your ssh key](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN).

3) Open up R or Rstudio and install the required packages by entering the following into the console.
Note that when our end users install the package, the dependencies will be automatically installed.
But during development we need to install them manually.

```r
# required for development
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))

# Install other packages used in the library
devtools::install_dev_deps(dependencies=TRUE)
```

NOTE: pay close attention to the console output as these are installing as some of these packages require external software to be installed.
You may need to install packages like `xml2`, `libcurl` etc for the packages to install correctly.
Often R gives helpful output when things are missing. If you are having trouble, create an issue on GitHub.


4) If you are using RStudio click the "new project" button, select "Version Control" then "Git", and then enter the details for the project repo.
If you've setup your SSH credentials make sure you enter the SSH URL to the repository i.e. `git@github.com:CCDM-Programming-Club/codonfriend.git`.
If you've already cloned the git repository, pull the latest version of it from github (`git pull origin master`) and open the project in Rstudio by clicking `File>Open Project`, navigate to the Git repo, and select `codonfriend.Rproj`.

5) Now you are ready to work!


NB. If you have lots of old packages you'll probably need to update them.
You might also like to see if you can get set up in a Packrat virtual environment.


## Checking the package

"Checking" the package just makes sure that everything is configured properly, and that any dependencies you have are installed etc.

In Rstudio you can click `Check` in the Build tab in the top-right hand window pane.
Otherwise run the checks using devtools.

```r
devtools::check()
```

You'll get a couple of "notes" about non-standard files etc.
But otherwise everything should be fine.

Try to rerun checks fairly regularly to catch problems quickly.


## Build the package

R packages are basically just glorified `RData` objects.
Building the package during development prepares that RData object, so that you can import it like you would load any other package.

In RStudio you can click `Install and Restart` in the Build tab in the top-right hand window pane.
Otherwise you can build with devtools.

```r
devtools::install()
library("codonfriend")
```

You'll now have access to any exported functions from your library in the R console.


## Run the tests

We will run some tests using [testthat](http://testthat.r-lib.org/) to help make sure that our code works and is correct.
This will help us check that changes we make don't break each others' or earlier work.

There is a software development principal called "Test Driven Development" that tends to force you to write good, maintainable, tidy code (which is something we should all strive for).
Here is a blog post about TDD in R if you're interested <http://rstudio-pubs-static.s3.amazonaws.com/278724_4d8935a2955c49d9934e2113c737e70e.html>.

The basic idea is that you write a test first, thinking about how you want to call the function, what output it should have given some input, and what "edge" cases your function might encounter that cause it to behave unexpectedly.
You then write the actual function __after__ the test.

So:

1) Write a test that fails (because you haven't written the code yet, or you haven't adjusted your code to handle it yet)
2) Write or update the function to make that test pass.
3) Run the tests.
4) Restart from 1.

It's also a good way to leave reminders to do something, since the failing test won't pass until you finish the task.


To run the unit tests in RStudio you can click `More>Test Package` in the Build tab in the top-right hand window pane.
Or using devtools you can run:

```r
devtools::test()
```

In the console you should then see whether tests pass or fail, and for failures you'll get details of how the test failed.


## Style guide

Having a consistent coding style is important for readability and comprehensibility.
We each develop our own little shortcuts for conveying information, but are rarely consistent with each other or even ourselves over time.
Please adhere to the style-guide layed out in the "Advanced R" book <http://adv-r.had.co.nz/Style.html> and defer to googles style guide where Hadley doesn't specify it <https://google.github.io/styleguide/Rguide.xml>.
If you see that something doesn't adhere to that standard, please change it and submit a pull request.

To help keep style consistent you can run a [linter](https://en.wikipedia.org/wiki/Lint_(software)) on the code using `devtools` and `lintr`.

```r
devtools::lint()
```


## Exporting functions from the package

By default all functions in R packages are private, meaning that we can't access them outside of the package's source code.
To make a function available to users (and to our vignettes) you need to "export" it.
Technically this is specified in the `NAMESPACE` file, but in practise we use a package to populate this automatically for us.

To make a function public put the following special comment above it in the source code.

```r
#' @export
foo <- function(x, y, z) {
  ...
}
```

Learn more about the `NAMESPACE` file see <http://r-pkgs.had.co.nz/namespace.html>


## Using external packages within our package

It is considered bad practice to put `library()` or `require()` function calls in our R package's source code.
Instead this should be managed by the `NAMESPACE` file, but like the exports we also do this with special comments.

To make an external function available to a function in the package add this special comment above the function definition.

```r
#' @importFrom stringr str_match
foo <- function(x, y) {
  z <- str_match(x, y)
  return(z)
}
```

This will make the function `str_match` from the `stringr` package available to functions within our package.
The next time you build the project you'll see the function name appear in the `NAMESPACE` file.
Alternatively you can call the function by explicitly specifying the package each time using the `::` operator.
The following is equivalent to the previous code block.

```r
foo <- function(x, y) {
  z <- stringr::str_match(x, y)
  return(z)
}
```

The `::` alternative won't add an entry to the `NAMESPACE` file, and incurs a slight performance penalty.

Learn more about the `NAMESPACE` file see <http://r-pkgs.had.co.nz/namespace.html>


## Adding a dependency for the package

If we know that our package must use functions from another package, we need to tell R to install them when users install our package.
Dependencies are specified in the `DESCRIPTION` file.
Packages listed in the `Imports:` section are used in the package itself and will be installed with the package.
Packages listed in the `Suggests:` section are only needed to build the package, so includes things like `devtools` and `knitr` etc.

To add a dependency on another package, either edit the `DESCRIPTION` file or use devtools.

```r
devtools::use_package("ggplot2")
```

This will add that package to the `DESCRIPTION` file.
You can also specify whether it should be in the `Imports` or `Suggests` section by specifying the `type` argument.


## Documentation

We can automatically document how to use functions by adding special comments before the functions.
Check out this <http://r-pkgs.had.co.nz/man.html> page for details of how to write them, or follow the example of existing functions.

Generating the documentation will populate the `man` folder automatically based on these special comments.


## Vignettes

Vignettes are like tutorials.
In our project they are written as [Rmarkdown](https://rmarkdown.rstudio.com/) files.
You can run all of the code and output the vignette by opening the `.Rmd` file and clicking `Knit` at the top of the editor window.


## Example data

We currently have two example data files.
Both are from an early _Leptosphaeria maculans_ genome and are output from EMBOSS tools.
One is the codon usage table (`lmac.cut`) and the other has the codon adaptation indices for each gene (`lmac.cai`).

You can find the path to these files on any system with codonfriend installed using the commands.

```r
cai_path <- system.file("extdata", "lmac.cai", package = "codonfriend")
cut_path <- system.file("extdata", "lmac.cut", package = "codonfriend")
```

Any files that you place in the folder `inst/extdata` will be available like this.


## Using Git effectively

When you are working with other people on a Git project, you should avoid committing directly to the `master` git branch.
Instead, when you want to create a new feature or function you create a new branch from `master`, then you can work in that `feature` branch for as long as you want, and when your changes are complete you merge the `feature` and `master` branches together using a "pull request".

There are a few reasons why we do this.

1) It avoids having your work conflict with others' work.
2) The pull request should always be reviewed by someone else to make sure it looks ok. 
   E.g. is correct, passes all unit tests, has proper style etc.
3) If you really mess something up or decide you should have taken a different approach, you can just delete the branch and start again without it affecting anyone elses work.
4) It makes it easier to see when bugs were introduced (and therefor how we might fix them).


Here is a basic workflow that you might try:

1) Pull the most recent master branch, `git pull origin master`.
2) Create and checkout a new branch off of the master branch, `git checkout -b myfeature`.
   NB. you can go back to the master branch at anytime without losing `myfeature` with the command `git checkout master`.
3) Make whatever changes that you need to, you can add commits to the `myfeature` branch as you go.
4) When you've finished and are ready for other people to look at merging the feature, you can push it to GitHub using
   `git push -u origin myfeature` (make sure you've checked out the actual branch that you want to push i.e `myfeature`).
   The `-u` flag means that it will create a new branch in the GitHub repo,
   this isn't normally needed since the `master` branch already exists.
5) Within GitHub, open up your branch and click "Send Pull Request".
   Someone else will review your work, and add comments or suggestions.
   You can still make edits to your feature branch and push them to the new branch on GitHub.
   Those changes will be included in the pull request.
6) You're done! Well done!


In RStudio you can do all of the things i described for the shell from the `Git` tab in the top-right hand window pane.
You can switch branches from the drop down menu at the far right.
The button to create new branches is just next to that dropdown box (three little purple boxes in an "L" shape).

NOTE: If you don't have some of the dropdown boxes, please try updating R-studio to the latest version.
I'm using `Version 1.1.456`.


Happy Git-ing!
