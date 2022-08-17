
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

-   In some cases, there may not be a bridge channel included, and in
    these cases we will need to use the normalize_1plex() function
    instead

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
 
head(data)                   
#> # A tibble: 6 × 8
#>   Sample TMT   ProteinID  value protein_avg intermediate_norm median_o…¹ final…²
#>   <chr>  <chr> <chr>      <dbl>       <dbl>             <dbl>      <dbl>   <dbl>
#> 1 PCB002 126   A0A024R0K5 3419.       1803.             328.       186.    232. 
#> 2 PCB002 127C  A0A024R0K5 1279.       1803.             123.        88.1   183. 
#> 3 PCB002 127N  A0A024R0K5 2368.       1803.             227.       101.    296. 
#> 4 PCB002 128C  A0A024R0K5  591.       1803.              56.7      181.     41.3
#> 5 PCB002 128N  A0A024R0K5 1656.       1803.             159.       131.    159. 
#> 6 PCB002 129C  A0A024R0K5 1248.       1803.             120.       139.    113. 
#> # … with abbreviated variable names ¹​median_of_sample_plex, ²​final_norm
```

The final norm column has the completely normalized data.

Here we can see an example of the overall workflow using an example set
with multiple plexes

``` r
data = read_delim("tests/testdata/normalize_to_bridge/PSM_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_to_bridge(bridge_channel_plex = 126)


head(data)
#> # A tibble: 6 × 8
#>   Sample   TMT   ProteinID  value bridge_values intermediate_n…¹ sampl…² final…³
#>   <chr>    <chr> <chr>      <dbl>         <dbl>            <dbl>   <dbl>   <dbl>
#> 1 DG014843 127C  A0A024R6I7  351.          350.             81.1    80.7    79.3
#> 2 DG014843 127N  A0A024R6I7  338.          350.             78.1    66.9    92.2
#> 3 DG014843 128C  A0A024R6I7  320.          350.             73.8    76.2    76.4
#> 4 DG014843 128N  A0A024R6I7  241.          350.             55.5    70.9    61.8
#> 5 DG014843 129C  A0A024R6I7  436           350.            101.     84.0    94.6
#> 6 DG014843 129N  A0A024R6I7  408.          350.             94.3    81.8    91.0
#> # … with abbreviated variable names ¹​intermediate_norm, ²​sample_plex_medians,
#> #   ³​final_norm
```

If you prefer data in a wide format ( as reported using the previous
script), you can specify data_format = “wide” within the normalize to
bridge function. This will give a column containing the final_norm
values for each Sample/TMT combination. Note that this is also an option
for the normalize_1plex function. An example is shown below

``` r
data = read_delim("tests/testdata/normalize_to_bridge/PSM_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_to_bridge(bridge_channel_plex = 126,data_format = "wide")


head(data)
#> # A tibble: 6 × 1,621
#>   ProteinID  DG014843_…¹ DG014…² DG014…³ DG014…⁴ DG014…⁵ DG014…⁶ DG014…⁷ DG014…⁸
#>   <chr>            <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 A0A024R6I7        79.3    92.2    76.4    61.8    94.6    91.0    97.4    66.4
#> 2 A0A075B6H7        73.0    66.6    62.6    71.7    61.0    53.8    69.1    85.2
#> 3 A0A075B7D0        94.4   127.     80.7   139.     73.0   129.    116.    156. 
#> 4 A0A087WVC6        76.2    82.2    60.6    78.3    76.8    79.7    73.8    72.9
#> 5 A0A096LPE2        80.7   142.     79.1    57.8    90.1    72.6    71.1    97.7
#> 6 A0A0A0MSV6        72.9    65.7    84.8    66.8    65.9    67.3    67.9    87.6
#> # … with 1,612 more variables: DG014843_131C <dbl>, DG014843_131N <dbl>,
#> #   DG014843_132C <dbl>, DG014843_132N <dbl>, DG014843_133C <dbl>,
#> #   DG014843_133N <dbl>, DG014843_134N <dbl>, DG014844_127C <dbl>,
#> #   DG014844_127N <dbl>, DG014844_128C <dbl>, DG014844_128N <dbl>,
#> #   DG014844_129C <dbl>, DG014844_129N <dbl>, DG014844_130C <dbl>,
#> #   DG014844_130N <dbl>, DG014844_131C <dbl>, DG014844_131N <dbl>,
#> #   DG014844_132C <dbl>, DG014844_132N <dbl>, DG014844_133C <dbl>, …
#> # ℹ Use `colnames()` to see all variable names
```

If we wanted to then use the box cox norm method, we could use
la_box_cox_norm() function.

``` r
data = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_to_bridge(bridge_channel_plex = 126) %>% 
  la_box_cox_norm()


head(data)
#> # A tibble: 6 × 5
#>   Sample TMT   ProteinID  final_norm box_cox_scaled_values
#>   <chr>  <chr> <chr>           <dbl>                 <dbl>
#> 1 PCB002 127C  A0A024R0K5      102.                  0.940
#> 2 PCB002 127N  A0A024R0K5      158.                  1.05 
#> 3 PCB002 128C  A0A024R0K5       25.9                 0.704
#> 4 PCB002 128N  A0A024R0K5       88.7                 0.899
#> 5 PCB002 129C  A0A024R0K5       63.3                 0.789
#> 6 PCB002 129N  A0A024R0K5       71.3                 0.840
```

Note that when using this function, it is expecting your data to be in
the long format. If not, it will not work. If you wish for the results
reported from la_box_cox_nrom to be in the wide format, you can specify
that like below.

``` r
data = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_to_bridge(bridge_channel_plex = 126) %>% 
  la_box_cox_norm(data_format = "wide")


head(data)
#> # A tibble: 6 × 5
#>   Sample TMT   ProteinID  final_norm box_cox_scaled_values
#>   <chr>  <chr> <chr>           <dbl>                 <dbl>
#> 1 PCB002 127C  A0A024R0K5      102.                  0.940
#> 2 PCB002 127N  A0A024R0K5      158.                  1.05 
#> 3 PCB002 128C  A0A024R0K5       25.9                 0.704
#> 4 PCB002 128N  A0A024R0K5       88.7                 0.899
#> 5 PCB002 129C  A0A024R0K5       63.3                 0.789
#> 6 PCB002 129N  A0A024R0K5       71.3                 0.840
```

Out of curiosity, what does the normalized data look like when comparing
the two methods.

``` r

p1 = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_to_bridge(bridge_channel_plex = 126) %>% 
  la_box_cox_norm() %>% 
  ggplot(aes(box_cox_scaled_values))+
  geom_histogram()+
  ggtitle("box cox scaled values")



p2 = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_to_bridge(bridge_channel_plex = 126) %>% 
  ggplot(aes(final_norm))+
  geom_histogram()+
  ggtitle("batch corrected")

p3 = ggpubr::ggarrange(p1,p2)

p3
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" /> Here
we can see the overall distribution of the data. Looks like the box cox
does a pretty good job of normalizing the data!

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

#using function to replicate old process
old = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions_replica() %>% 
  mutate(method = "old")

#Combining data
comb = bind_rows(new,old) %>% 
  pivot_wider(names_from = method, values_from = value)


#Plotting
p1 = comb %>% 
  ggplot(aes(new,old))+
  geom_point()+
  geom_smooth(method = "lm")

p1
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

As we can see above, while this data is correlated, it is not exactly
the same. **Future work should only use the combine_psm_fractions()
function.**

### Determining Differential Protein Abundance

Usually when analyzing proteomics data, we will be interested in
assessing the differences between two or more groups.

There are a number of ways to do this, but the most simple and most
commonly used approach to visualize the differentally abundant proteins
between tow groups would be by using a volcano plot. We can use GLabR to
easily plot a volcano plot as displayed in the code below.

First we have to load our data, and normalize it appropriately. To do
this we use the functions in GLabR that were previously described.

``` r
data = readr::read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>% 
  combine_psm_fractions() %>% 
  normalize_1plex() %>% 
  la_box_cox_norm() %>% 
  #Have to add the sample ID column manually by concating Sample and TMT.
  #Done outside of GLabR to allow the user freedom in determining naming convention and contents (. vs _ seperator etc) of what I would call Sample_ID
  dplyr::mutate(Sample_ID = paste0(Sample,".",TMT))
```

Next, we need to assign metadata that gives us information about what
each of these samples are.

``` r
# Loading metadata
md = readr::read_csv("tests/testdata/metadata.csv") %>% 
  #Removed redundant columns to prevent .x and .y columns in data_md
  dplyr::select(-Sample,-TMT)

# Appending md to data
data_md = dplyr::inner_join(data,md,by = "Sample_ID")
```

Now we can use the volcano_plot() function to plot the differences
between two types of samples using conditions contained in the metadata.
Here let’s compare the greatest disease severity to healthy controls
using the mayo score.

Note that the function has trouble dealing with weird column names, so
to avoid this use column names without spaces, :s, or other characters
that are dealt with differently in R

Also, note that you can change the p_theshold and log2fc threshold using
arguments. Any proteins above these thresholds will be plotted.

``` r
# Filtering our data to contain only the two conditions we want to test
# NEed to fix this later
f_data_md = data_md %>% 
  dplyr::filter(`Mayo_Endoscopic_Sub_Score` %in% c("Healthy_control","3: Severe disease (spontaneous bleeding, ulceration)")) 

#

data = f_data_md
column_split_by = "Mayo_Endoscopic_Sub_Score"

 stat = data %>%
      dplyr::group_by(ProteinID)
      tryCatch({
      stat = rstatix::t_test(stat,as.formula(paste("box_cox_scaled_values",'~',column_split_by)))
      stat = rstatix::adjust_pvalue(stat,p.col = "p",output.col = "p.adj_fdr", method = "fdr")
      }, error = function(cond){
        return(NA)
        })
volcano_plot(f_data_md,"Mayo_Endoscopic_Sub_Score",p_threshold = 0.05,fc_threshold = 1)
```

Here we can see our volcano plot! We can note that there are 5 proteins
that meet our criteria, and are called out by ProteinID on the plot.
Note that this function uses the t_test function from the excellent
rstatix package to determine statistical significance. Importantly, the
pvalue reported has been adjusted for multiple comparisons using fdr.

### Protein Idenfification.

The next step in the process is figuring out what these proteins
actually are/ what they do. In order to do this, we can use a
combination of extract_sig_proteins and our annotate proteins function.

extract_sig_proteins is a function that is intended to extract the
proteinIDs that are significant given the desired parameters. If the
same parameters are passed to the function as volcano plot, the proteins
will match between the two. We can see this below.

``` r
sig_proteins = extract_sig_proteins(f_data_md,column_split_by = "Mayo_Endoscopic_Sub_Score",p_threshold = 0.05,fc_threshold = 1)

sig_proteins
```

Now that we have a list of our proteins we can figure out what they are
using the annotate_proteins function.

This function uses Uniprot’s API to return results pertaining to the
proteins. There is a package called UniprotR that handles this, but it
was way too slow for long lists of protein IDs, so I made this solution.
We can annotate these proteins as follows.

``` r

annotated_proteins = annotate_proteins(sig_proteins)

annotated_proteins
```

Let’s say you wanted to get different information about these proteins.
To do that, you can specify what columns to return through the columns
argument. This argument takes a string of the column names that you want
to add separated by columns. The complete list of field names that are
accepted can be found here:
<https://www.uniprot.org/help/return_fields>.

In this example, lets say we want to get the GO biological process terms
for our proteins (column name accessed through api by “go_p”. We can do
this as follows.

``` r
# Result table will have accesssion and GO biological process info
annotated_proteins_GO_p = annotate_proteins(sig_proteins,columns ="accession,go_p")

annotated_proteins_GO_p
```

Here we can see that information!

### Phosphoproteomics

The above code is great for looking at proteomics, but sometimes we will
want to look at phospo post translational modifications.

The calc_phospho_ratio function allows us to do this! We can see an
example below:

``` r
proteomics_data = readr::read_delim("tests/testdata/psm_phospho_mod/PCB002_PSMs.txt")
phospho_data = readr::read_delim("tests/testdata/psm_phospho_mod/PCB001_PSM.txt")
metadata = readxl::read_excel("tests/testdata/psm_phospho_mod/metadata.xlsx")
#In this dataset, patient ID is the variable that we want to use to  combine the phospho and proteomic data. 
col_identifying_match = "PatientID"

phospho_ratio = calc_phospho_ratio(proteomics_data,phospho_data,metadata,col_identifying_match,1)

head(phospho_ratio)
#> # A tibble: 6 × 7
#>   PatientID ProteinID  Annotated_Sequence ptmRS phospho_box_co…¹ prote…² Phosp…³
#>   <chr>     <chr>      <chr>              <chr>            <dbl>   <dbl>   <dbl>
#> 1 0572      A0A024R0K5 [K].CETQNPVSAR.[R] NA               1.25    1.28    0.979
#> 2 C2        A0A024R0K5 [K].CETQNPVSAR.[R] NA              NA       1.11   NA    
#> 3 C4        A0A024R0K5 [K].CETQNPVSAR.[R] NA              NA       1.42   NA    
#> 4 0530      A0A024R0K5 [K].CETQNPVSAR.[R] NA               0.456   0.744   0.613
#> 5 0672      A0A024R0K5 [K].CETQNPVSAR.[R] NA               0.833   1.09    0.765
#> 6 0185      A0A024R0K5 [K].CETQNPVSAR.[R] NA              NA       0.961  NA    
#> # … with abbreviated variable names ¹​phospho_box_cox_scaled_values,
#> #   ²​proteomics_box_cox_scaled_values, ³​Phospho_Prot_ratio
```

In short, this code pairs up corresponding proteomics and
phosphoproteomics data sets and then calculates the ratio of
phosphorylation to using the box cox normalized data.

Sometimes, it will be useful to compare how many unique
peptide-sequences-PTM are present in a phosphoenriched experiment,
compared to a normal proteomics experiment. Not that in order to do
this, you must have your proteome discoverer analysis set up to report
infomration about phosphorylated peptides. You can do this by setting
dynamic settings to accomodate phospho mods in the sequest node, and
enabling the ptmRS node.

Once you have these results, you can compare them using the following
code. Note that for this to make any sense to do, the experiments should
contain the same samples.

``` r
phospho_enriched = readr::read_delim("tests/testdata/psm_phospho_mod/PCB001_PSM.txt")

proteomics = readr::read_delim("tests/testdata/phospho_venn_diagram/PCB002_proteomics_ptmRS_data.txt")

phospho_venn_diagram(proteomics,phospho_enriched)
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />
Here we can see that the phospho enriched experiment enriched for
phospho peptides, as we would expect.

### Designing Experiments

GLabR also has functions that are intended to help design experiments.
All together, these functions are designed to take a data frame of
metadata, assign a plate number, randomly assign a 96 well location,
assign a TMT plex number with an incorporated bridge channel (allowing
for the option to group samples within the same plex by a factor if
desired), and then randomly assign a TMT label for each sample in the
plex.

Note that these functions are designed to be modular. If you already
have metadata with the plate numbers that were generated by some other
method that you want to use (ie manually in excel), than no need to
randomize the plates and so on.

Right now, the 10 and 16 plex TMT options are supported.

We can see an example as follows for a case where we have 108 samples
across and would like to pair by mouse.

``` r
metadata = data.frame(sample_id = as.character(1:108),
                      mouse_id = rep(paste0("Mouse",1:27),each = 4)
                    )

#define metadata when using 16 plexes
final_metadata = metadata %>% 
  randomize_plates() %>% 
  randomize_wells() %>% 
  assign_plex_num("mouse_id") %>% 
  randomize_tmt()

head(final_metadata)
#> # A tibble: 6 × 6
#> # Groups:   plex_num [1]
#>   sample_id mouse_id plate well  plex_num tmt_label
#>   <chr>     <chr>    <dbl> <chr>    <dbl> <chr>    
#> 1 1         Mouse1       1 D3           1 132C     
#> 2 2         Mouse1       1 F3           1 128N     
#> 3 3         Mouse1       1 B2           1 130N     
#> 4 4         Mouse1       1 C1           1 134N     
#> 5 5         Mouse2       1 G3           1 131N     
#> 6 6         Mouse2       1 F9           1 133N
```

What if for some reason we wanted to use 10 plexes instead. We can do
that as follows.

-   need to fix this before stable release

### Pubchem

We can also use the annotate pubchem function to search Smiles for
synonyms.

Note that this can only find exact SMILE matches, it will not report
results for an isomer like the web portal will

Also note that this will return warnings when a search string cannot be
found in pubchems database.

lastly please not that the way I am determining what name to return
isn’t perfect. There are tons of synonyms for each compound, so
returning the proper name doesnt seem to be trivial. Here if the first
name at the top of the returned list has only numbers, I am then going
for the next entry in list. In the future I can patch this with a better
fix but this should be sufficient for now.

``` r
test = read.delim("tests/testdata/annotate_pubchem/Random_metabolomics_calls.txt") %>%
  pull(Smiles)


t = annotate_pubchem(test)

head(t)
#> # A tibble: 6 × 2
#>   input                                                                   output
#>   <chr>                                                                   <chr> 
#> 1 "[H]C([H])([H])Oc(n2)nc(c([H])c(OC([H])([H])[H])2)N([H])S(=O)(=O)c(c([… Sulph…
#> 2 "CN(C)CCOC(C1=CC=CC=C1)C1=CC=CC=C1"                                     Benad…
#> 3 "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)O"                                      Indol…
#> 4 "O=C(\\C=C\\C=C\\C1=CC2=C(OCO2)C=C1)N1CCCCC1"                           1-Pip…
#> 5 "CCCCC(=O)N(CC1=CC=C(C=C1)C1=CC=CC=C1C1=NNN=N1)[C@@H](C(C)C)C(O)=O"     Diovan
#> 6 "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC… Chola…
```

Note right now this doesn’t do a great job as sometimes the synonym
returned isn’t the one you would want ie “58-73-1” isnt helpful. Also is
rather slow. Will figure out how to fix and update in future.
