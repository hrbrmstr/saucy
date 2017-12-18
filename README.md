
# saucy

Fast Symmetry Discovery Tool for Sparse Graphs

## Description

Many computational tools have recently begun to benefit from the use of
the symmetry inherent in the tasks they solve, and use general-purpose
graph symmetry tools to uncover this symmetry. Common applications
include organic chemistry, constraint solvers, logistics, optimization,
bio-informatics, and finite group theory. However, older
symmetry-finding tools often suffer quadratic runtime (or worse\!) in
the number of symmetries explicitly returned and are therefore of
limited use on very large, sparse, symmetric graphs. Over the last 10+
years, we developed a symmetry-discovery algorithm which exploits the
sparsity present not only in the input but also the output, i.e., the
symmetry generators themselves. By avoiding quadratic runtime on large
graphs, it improved state-of-the-art runtimes from several days to less
than a second. Recent improvements to our algorithm include additional
pruning that quickly solves the hard Miyazaki graphs.

## NOTE

I just wrapped code. <http://vlsicad.eecs.umich.edu/BK/SAUCY> for orig
ref. The following are the authors of the core C “library”:

  - Paul T. Darga <pdarga@umich.edu>
  - Mark Liffiton <liffiton@umich.edu>
  - Hadi Katebi <hadik@eecs.umich.edu>

## What’s Inside The Tin

The following functions are implemented:

  - `saucy`: Searching for Automorphisms in Underlying CNF, yes?
  - `shatter`: Identifies Symmetries in CNF Instances

## Installation

``` r
devtools::install_github("hrbrmstr/saucy")
```

## Usage

``` r
library(saucy)

# current verison
packageVersion("saucy")
```

    ## [1] '0.1.0'

### Saucy

Graph 1 (from example
files)

``` r
graph1 <- saucy::saucy(system.file("extdata", "graphfile", package="saucy"))

graph1
```

    ## (0 2)(3 4)
    ## (0 1)(2 4)
    ## 
    ##           input_file: graphfile
    ##             vertices: 5
    ##                edges: 5
    ##      group_size_base: 1
    ##       group_size_exp: 1
    ##               levels: 3
    ##                nodes: 7
    ##           generators: 2
    ##        total_support: 8
    ##      average_support: 4
    ##  nodes_per_generator: 3.5
    ##            bad_nodes: 0

``` r
str(graph1)
```

    ## List of 13
    ##  $ input_file         : chr "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/saucy/extdata/graphfile"
    ##  $ vertices           : int 5
    ##  $ edges              : int 5
    ##  $ group_size_base    : num 1
    ##  $ group_size_exp     : int 1
    ##  $ levels             : int 3
    ##  $ nodes              : int 7
    ##  $ generators         : int 2
    ##  $ total_support      : int 8
    ##  $ average_support    : num 4
    ##  $ nodes_per_generator: num 3.5
    ##  $ bad_nodes          : int 0
    ##  $ printed_output     : chr [1:2] "(0 2)(3 4)" "(0 1)(2 4)"
    ##  - attr(*, "class")= chr [1:2] "saucy" "list"

Graph 2 (from example
files)

``` r
graph2 <- saucy::saucy(system.file("extdata", "graphfile2", package="saucy"))

graph2
```

    ## (3 5)
    ## (3 6)(4 5)
    ## (0 1)
    ## (0 2)
    ## 
    ##           input_file: graphfile2
    ##             vertices: 7
    ##                edges: 7
    ##      group_size_base: 4.8
    ##       group_size_exp: 1
    ##               levels: 5
    ##                nodes: 13
    ##           generators: 4
    ##        total_support: 10
    ##      average_support: 2.5
    ##  nodes_per_generator: 3.25
    ##            bad_nodes: 0

``` r
str(graph2)
```

    ## List of 13
    ##  $ input_file         : chr "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/saucy/extdata/graphfile2"
    ##  $ vertices           : int 7
    ##  $ edges              : int 7
    ##  $ group_size_base    : num 4.8
    ##  $ group_size_exp     : int 1
    ##  $ levels             : int 5
    ##  $ nodes              : int 13
    ##  $ generators         : int 4
    ##  $ total_support      : int 10
    ##  $ average_support    : num 2.5
    ##  $ nodes_per_generator: num 3.25
    ##  $ bad_nodes          : int 0
    ##  $ printed_output     : chr [1:4] "(3 5)" "(3 6)(4 5)" "(0 1)" "(0 2)"
    ##  - attr(*, "class")= chr [1:2] "saucy" "list"

### Shatter

``` r
s1 <- saucy::shatter(system.file("extdata", "battleship.cnf", package="saucy"))

s1
```

    ##                    input_file: battleship.cnf
    ##                      vertices: 105
    ##                         edges: 320
    ##               group_size_base: 1.6
    ##                group_size_exp: 3
    ##                         nodes: 236
    ##                    generators: 6
    ##                     bad_nodes: 151
    ##                discovery_time: 0.000816
    ##     symmetry_breaking_clauses: 396
    ##          additional_variables: 102
    ##           additional_literals: 1362
    ##        consistency_violations: 0
    ##           sbp_generation_time: 0.000252
    ##                    total_time: 0.001068
    ##                     variables: 40
    ##                       clauses: 105
    ##            non_binary_clauses: 25
    ##                      literals: 360
    ##  symmetry_breaking_predicates: 293 symmetry breaking predicates
    ##                printed_output: character(0)

``` r
str(s1)
```

    ## List of 21
    ##  $ input_file                  : chr "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/saucy/extdata/battleship.cnf"
    ##  $ vertices                    : int 105
    ##  $ edges                       : int 320
    ##  $ group_size_base             : num 1.6
    ##  $ group_size_exp              : int 3
    ##  $ nodes                       : int 236
    ##  $ generators                  : int 6
    ##  $ bad_nodes                   : int 151
    ##  $ discovery_time              : num 0.000816
    ##  $ symmetry_breaking_clauses   : int 396
    ##  $ additional_variables        : int 102
    ##  $ additional_literals         : int 1362
    ##  $ consistency_violations      : int 0
    ##  $ sbp_generation_time         : num 0.000252
    ##  $ total_time                  : num 0.00107
    ##  $ variables                   : num 40
    ##  $ clauses                     : num 105
    ##  $ non_binary_clauses          : num 25
    ##  $ literals                    : num 360
    ##  $ symmetry_breaking_predicates: chr [1:293] "-6 31 0" "41 0" "-41 -6 -7 32 0" "-41 -6 42 0" ...
    ##  $ printed_output              : chr(0) 
    ##  - attr(*, "class")= chr [1:2] "shatter" "list"
