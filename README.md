
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GLabR

<!-- badges: start -->
<!-- badges: end -->

The goal of GLabR is to provide a centralized locations to hold
functions that are routinely useful in the Gonzalez Lab at UCSD

## Installation

You can install the development version of GLabR like so:

``` r
devtools::install_github("baynec2/GlabR")
```

## Examples

### Data normalization (batch correction)

The first use case of GlabR is to normalize data that we export from
proteome discoverer. Here, a text file containing PSMs is exported and
subsequently needs to be processed. To accomplish this, there are a few
different steps that need to be followed . These are described below.

1.  We need to combine PSMs from each of the different fractions (there
    are multiple fractions corresponding to a single sample). This is
    accomplished using the combine_psm_fractions() function.
    Essentially, this function filters out PSMs based on certain
    criteria and then sums the intensities for each protein that they
    map to.

2.  Then we need to normalize our data to account for batch corrections.
    This is accomplished using the normalize_to_bridge() function.

<!-- -->

1.  In some cases, there may not be a bridge channel included, and in
    these cases we will need to use the normalize_1plex() function
    instead

<!-- -->

3.  After this, we still might want to normalize the data further. To
    get data that is more normally distributed, we can use Leigh-ana’s
    method of box cox normalization through the la_box_cox_norm()
    function.

Here we can see an example of the overall workflow using a single plex
example set

``` r
library(GLabR)

data = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_1plex()
#> Rows: 13131 Columns: 54
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (13): Confidence, Identifying Node, PSM Ambiguity, Annotated Sequence, M...
#> dbl (40): PSMs Workflow ID, PSMs Peptide ID, # Proteins, # Missed Cleavages,...
#> lgl  (1): Checked
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
 
head(data)                   
#> # A tibble: 6 × 8
#>   Sample TMT   ProteinID  value protein_avg intermediate_norm median_of_sample_…
#>   <chr>  <chr> <chr>      <dbl>       <dbl>             <dbl>              <dbl>
#> 1 PCB002 126   A0A024R0K5 3419.       1803.             327.               183. 
#> 2 PCB002 127C  A0A024R0K5 1279.       1803.             122.                86.3
#> 3 PCB002 127N  A0A024R0K5 2368.       1803.             226.                98.9
#> 4 PCB002 128C  A0A024R0K5  591.       1803.              56.5              178. 
#> 5 PCB002 128N  A0A024R0K5 1656.       1803.             158.               129. 
#> 6 PCB002 129C  A0A024R0K5 1248.       1803.             119.               136. 
#> # … with 1 more variable: final_norm <dbl>
```

The final norm column has the completely normalized data.

Here we can see an example of the overall workflow using an example set
with multiple plexes

Note that when I was making this package, I discovered a problem in the
script that was originally used to normalize the data. Essentially, the
original script was only assigning the NAs from some columns as 1,
ultimately causing results to deviate from what was intended (NA from
all columns having 1)

Here I have kept the function that produces the same results as the
output from the original script (nonnormalizedall.txt) as
combine_psm_fractions_replica() for posterity.

Here we can see the slight differences produced from these different
functions.

``` r
#using function to replicate new process
new = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions() %>% 
  mutate(method = "new")
#> Rows: 13131 Columns: 54
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (13): Confidence, Identifying Node, PSM Ambiguity, Annotated Sequence, M...
#> dbl (40): PSMs Workflow ID, PSMs Peptide ID, # Proteins, # Missed Cleavages,...
#> lgl  (1): Checked
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

#using function to replicate old process
old = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions_replica() %>% 
  mutate(method = "old")
#> Rows: 13131 Columns: 54
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (13): Confidence, Identifying Node, PSM Ambiguity, Annotated Sequence, M...
#> dbl (40): PSMs Workflow ID, PSMs Peptide ID, # Proteins, # Missed Cleavages,...
#> lgl  (1): Checked
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

#Combining data
comb = bind_rows(new,old) %>% 
  pivot_wider(names_from = method, values_from = value)


#Plotting
p1 = comb %>% 
  ggplot(aes(new,old))+
  geom_point()+
  geom_smooth(method = "lm")

p1
#> `geom_smooth()` using formula 'y ~ x'
#> Warning: Removed 656 rows containing non-finite values (stat_smooth).
#> Warning: Removed 656 rows containing missing values (geom_point).
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

As we can see above, while this data is correlated, it is not exactly
the same. Future work should **only** use the combine_psm_fractions()
function.
