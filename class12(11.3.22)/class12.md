Class 12: RNASeq analysis
================
Vivian

In today’s class we will work with published RNA-seq experiment where
airway smooth muscle cells were treated with dexamethasone, a synthetic
glucocorticoid steroid with anti-inflammatory effects (Himes et
al. 2014).

## Data import

We will use good old ‘read.csv()’ to read the two things we need for
this analysis:

-count data -col data (metadata)

``` r
counts <- read.csv("airway_scaledcounts.csv", row.names=1)
metadata <-  read.csv("airway_metadata.csv")
```

``` r
#counts
#metadata
```

How many transcripts do I have? Q1: 38694 genes are in this dataset.

``` r
nrow(counts)
```

    [1] 38694

Lets have a look at the metadata…

``` r
head(metadata)
```

              id     dex celltype     geo_id
    1 SRR1039508 control   N61311 GSM1275862
    2 SRR1039509 treated   N61311 GSM1275863
    3 SRR1039512 control  N052611 GSM1275866
    4 SRR1039513 treated  N052611 GSM1275867
    5 SRR1039516 control  N080611 GSM1275870
    6 SRR1039517 treated  N080611 GSM1275871

First we should check the corespondence of the metadata and count data

``` r
metadata$id
```

    [1] "SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513" "SRR1039516"
    [6] "SRR1039517" "SRR1039520" "SRR1039521"

``` r
colnames(counts)
```

    [1] "SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513" "SRR1039516"
    [6] "SRR1039517" "SRR1039520" "SRR1039521"

To check that these are all in the same order, we can use ‘==’ test of
equality.

``` r
all( metadata$id == colnames(counts) )
```

    [1] TRUE

## Analysis via comparasion of CONTROL vs TREATED

The “treated” have the dex drug and the “control” do not. First I need
to be able to extract just the “control” columns in the ’counts” data
set.

``` r
sum(metadata$dex == "control")
```

    [1] 4

``` r
control.inds <- metadata$dex == "control"
control <- metadata[control.inds, ]
control$id
```

    [1] "SRR1039508" "SRR1039512" "SRR1039516" "SRR1039520"

Q2: There are 4 control cell lines.

Now I can use this to access just the “control” column in my ‘counts’
data…

``` r
control.counts <- counts[, control$id]
head(control.counts)
```

                    SRR1039508 SRR1039512 SRR1039516 SRR1039520
    ENSG00000000003        723        904       1170        806
    ENSG00000000005          0          0          0          0
    ENSG00000000419        467        616        582        417
    ENSG00000000457        347        364        318        330
    ENSG00000000460         96         73        118        102
    ENSG00000000938          0          1          2          0

Find the mean count value for each transcript/gene by finding the
‘rowMeans()’

``` r
control.mean <- rowMeans(control.counts)
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

And now find a mean value for all the “treated” columns in the same way
Q3: Approach on how to make the above cost more robust and condensed.
Q4: vector created called treated mean which follows the same procedure
for the treated samples

``` r
treated.id <- metadata[metadata$dex == "treated", "id"]
treated.mean <- rowMeans(counts[,treated.id])
```

Now I have ‘control.mean()’ and treated.mean.’ Let’s put them togetehr
for safe keeping and ease of use later.

``` r
meancounts <- data.frame(control.mean, treated.mean)
head(meancounts)
```

                    control.mean treated.mean
    ENSG00000000003       900.75       658.00
    ENSG00000000005         0.00         0.00
    ENSG00000000419       520.50       546.00
    ENSG00000000457       339.75       316.50
    ENSG00000000460        97.25        78.75
    ENSG00000000938         0.75         0.00

Q5: Let’s do a quick plot to see how our data looks

``` r
plot(meancounts)
```

![](class12_files/figure-gfm/unnamed-chunk-12-1.png)

This is very heavily skewed and over a wide range - calls out for a log
transform!

``` r
plot(meancounts, log = "xy")
```

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15032 x values <= 0 omitted
    from logarithmic plot

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15281 y values <= 0 omitted
    from logarithmic plot

![](class12_files/figure-gfm/unnamed-chunk-13-1.png)

Q5b: To make a log plot, you would need to use geom_point() Q6: To plot
the data on a log scale, the arguement is log.

We like working with log transfomred data as it could make things more
straight forward to interpret.

If we have no change:

``` r
log2(20/20)    #no change because the value is 0 
```

    [1] 0

What about if we had a doubling

``` r
log2(40/20)   
```

    [1] 1

``` r
log2(10/20)  
```

    [1] -1

``` r
log2(80/20)
```

    [1] 2

We like working wiht log2 fold-change values. Let’s calculate them for
our data.

``` r
meancounts$log2fc <- log2(meancounts$treated.mean / meancounts$control.mean)
head(meancounts)
```

                    control.mean treated.mean      log2fc
    ENSG00000000003       900.75       658.00 -0.45303916
    ENSG00000000005         0.00         0.00         NaN
    ENSG00000000419       520.50       546.00  0.06900279
    ENSG00000000457       339.75       316.50 -0.10226805
    ENSG00000000460        97.25        78.75 -0.30441833
    ENSG00000000938         0.75         0.00        -Inf

We want to filter out any genes (that is the rows) where we have ZERO
count data.

``` r
to.keep.inds <- rowSums(meancounts[,1:2] == 0) == 0 
head(to.keep.inds)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
               TRUE           FALSE            TRUE            TRUE            TRUE 
    ENSG00000000938 
              FALSE 

``` r
mycounts <- meancounts[to.keep.inds,]
nrow(mycounts)
```

    [1] 21817

A common threshold for calling genes as differentiating expressed is a
log2 fold-change of +2 or -2.

``` r
sum(mycounts$log2fc >= +2) 
```

    [1] 314

What percent is this?

``` r
round((sum(mycounts$log2fc >= +2)  / nrow(mycounts)) * 100,2)
```

    [1] 1.44

``` r
#down regulated
round((sum(mycounts$log2fc <= -2) / nrow(mycounts))* 100,2)
```

    [1] 2.22

Q7: The arr.ind arugment in the which function tells us which
genes(rows) and samples(columns) have zero counts. We are going to
ignore any genes with zero row counts and using unique() will make sure
we do not count any row twice if there are zero entries in both samples.
Q8. 1.44% of up regulated genes are greater than a 2 fc level . Q9:
2.22% of the down regulated genes are greater than a 2 fc level. Q10: I
do not trust these results yet because I need to do a statistical test
to make sure if this is significant or not regarding the changes.

We need some stats to check if the drug induced difference is
significance. \# Turn to DESeq2

Let’s turn to doing this the correct way with the DeSeq2 package

``` r
library(DESeq2)
```

The main function in the DESeq2 package is called ‘deseq()’. It wants
our count data and our colData (metadata) as input in a specific way.

``` r
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~dex)
```

    converting counts to integer mode

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

Now I can run the DESeq analysis

``` r
dds <- DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

``` r
results(dds)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 38694 rows and 6 columns
                     baseMean log2FoldChange     lfcSE      stat    pvalue
                    <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003  747.1942     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005    0.0000             NA        NA        NA        NA
    ENSG00000000419  520.1342      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457  322.6648      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460   87.6826     -0.1471420  0.257007 -0.572521 0.5669691
    ...                   ...            ...       ...       ...       ...
    ENSG00000283115  0.000000             NA        NA        NA        NA
    ENSG00000283116  0.000000             NA        NA        NA        NA
    ENSG00000283119  0.000000             NA        NA        NA        NA
    ENSG00000283120  0.974916      -0.668258   1.69456 -0.394354  0.693319
    ENSG00000283123  0.000000             NA        NA        NA        NA
                         padj
                    <numeric>
    ENSG00000000003  0.163035
    ENSG00000000005        NA
    ENSG00000000419  0.176032
    ENSG00000000457  0.961694
    ENSG00000000460  0.815849
    ...                   ...
    ENSG00000283115        NA
    ENSG00000283116        NA
    ENSG00000283119        NA
    ENSG00000283120        NA
    ENSG00000283123        NA

Now what we have got so far is the log2 fold-change and the adjusted
p-value for the significance.

``` r
res <- results(dds)
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 6 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj
                    <numeric>
    ENSG00000000003  0.163035
    ENSG00000000005        NA
    ENSG00000000419  0.176032
    ENSG00000000457  0.961694
    ENSG00000000460  0.815849
    ENSG00000000938        NA

``` r
plot(res$log2FoldChange, res$padj)
```

![](class12_files/figure-gfm/unnamed-chunk-26-1.png)

Well that plot sucked all the interesting P-values are down belwo zero.
I am going to take the log of the p-value

``` r
plot(res$log2FoldChange, log(res$padj))
```

![](class12_files/figure-gfm/unnamed-chunk-27-1.png)

We can flip the y-axis so the plot does not look upside down

``` r
plot(res$log2FoldChange, -log(res$padj))
abline(v=c(-2,+2), col = "blue")
abline(h=log(0.5), col = "red")
```

![](class12_files/figure-gfm/unnamed-chunk-28-1.png)

Making a volcano plot

``` r
# Setup our custom point color vector 
mycols <- rep("gray", nrow(res))
mycols[ abs(res$log2FoldChange) > 2 ]  <- "red" 

inds <- (res$padj < 0.01) & (abs(res$log2FoldChange) > 2 )
mycols[ inds ] <- "blue"

plot( res$log2FoldChange,  -log(res$padj), 
 col=mycols, ylab="-Log(P-value)", xlab="Log2(FoldChange)" )

# Cut-off lines
abline(v=c(-2,2), col="gray", lty=2)
abline(h=-log(0.1), col="gray", lty=2)
```

![](class12_files/figure-gfm/unnamed-chunk-29-1.png)

## Annotation of our gene set results

I will start by loading two Annotation packages from bioconductor:

``` r
library("AnnotationDbi")
library("org.Hs.eg.db")
```

The ‘mapIDs()’ function “maps” database identifies between different
databases. In other words it translates the identifiers used by one
database onto another database.

Let’s see what databases are available for human data

``` r
columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

My results are in the object ‘res’

``` r
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 6 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj
                    <numeric>
    ENSG00000000003  0.163035
    ENSG00000000005        NA
    ENSG00000000419  0.176032
    ENSG00000000457  0.961694
    ENSG00000000460  0.815849
    ENSG00000000938        NA

``` r
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),      # Our gene names
                     keytype="ENSEMBL",        # The format of our gene names
                     column="SYMBOL",          # The new format we want to add
                     multiVals="first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$uniprot <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="UNIPROT",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$genename <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 9 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj      symbol     uniprot               genename
                    <numeric> <character> <character>            <character>
    ENSG00000000003  0.163035      TSPAN6  A0A024RCI0          tetraspanin 6
    ENSG00000000005        NA        TNMD      Q9H2S6            tenomodulin
    ENSG00000000419  0.176032        DPM1      O60762 dolichyl-phosphate m..
    ENSG00000000457  0.961694       SCYL3      Q8IZE3 SCY1 like pseudokina..
    ENSG00000000460  0.815849    C1orf112  A0A024R922 chromosome 1 open re..
    ENSG00000000938        NA         FGR      P09769 FGR proto-oncogene, ..

# Pathway analysis

Pathway analysis (also known as gene set analysis or over-represented
analysis) aims to reduce the complexity of interpreting gene lists via
mapping the lists to known biological pathways, processes and functions

Some major gene sets include KEGG\< Go, etc We will use the **gage**
package for our first pathway analysis

``` r
library(pathview)
```

    ##############################################################################
    Pathview is an open source software package distributed under GNU General
    Public License version 3 (GPLv3). Details of GPLv3 is available at
    http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    formally cite the original Pathview paper (not just mention it) in publications
    or products. For details, do citation("pathview") within R.

    The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ##############################################################################

``` r
library(gage)
```

``` r
library(gageData)

data(kegg.sets.hs)
```

We can have a look at the first few pathways in the kegg human set

``` r
# Examine the first 2 pathways in this kegg set for humans
head(kegg.sets.hs, 2) 
```

    $`hsa00232 Caffeine metabolism`
    [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   

    $`hsa00983 Drug metabolism - other enzymes`
     [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"   "1551"  
     [9] "1553"   "1576"   "1577"   "1806"   "1807"   "1890"   "221223" "2990"  
    [17] "3251"   "3614"   "3615"   "3704"   "51733"  "54490"  "54575"  "54576" 
    [25] "54577"  "54578"  "54579"  "54600"  "54657"  "54658"  "54659"  "54963" 
    [33] "574537" "64816"  "7083"   "7084"   "7172"   "7363"   "7364"   "7365"  
    [41] "7366"   "7367"   "7371"   "7372"   "7378"   "7498"   "79799"  "83549" 
    [49] "8824"   "8833"   "9"      "978"   

``` r
# these are the blue cirlces as shown in the diagram
```

THe main ‘gage()’ function wants a vector as input that contains our
measure of importance - in our case the fold-change. The vector needs to
have ENTREZ ids assigned as the names of the vector.

Recall that vectors can have names - this is useful in book-jeeping so I
know what value corresponds to a certain gene for example.

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

    [1] -0.35070302          NA  0.20610777  0.02452695 -0.14714205 -1.73228897

Now we can run the analysis

``` r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

What is in this results object

``` r
attributes(keggres)
```

    $names
    [1] "greater" "less"    "stats"  

By default gage splits it’s results into “greater”, “less’ and”stats”
objects that you can examine. First we will look at teh “less” i.e. down
regulated results

``` r
# Look at the first three down (less) pathways
head(keggres$less, 3)
```

                                             p.geomean stat.mean p.val q.val
    hsa00232 Caffeine metabolism                    NA       NaN    NA    NA
    hsa00983 Drug metabolism - other enzymes        NA       NaN    NA    NA
    hsa01100 Metabolic pathways                     NA       NaN    NA    NA
                                             set.size exp1
    hsa00232 Caffeine metabolism                    0   NA
    hsa00983 Drug metabolism - other enzymes        0   NA
    hsa01100 Metabolic pathways                     0   NA

We can look in more detials at these pathways using the ‘pathview()’
function will take the KEGG pathway ID (printed above) and our vector of
importantance and annotate the pathway with our genes First I will look
at the hsa05310 pathway.

``` r
pathview(gene.data=foldchanges, pathway.id="hsa05310")
```

    Warning: None of the genes or compounds mapped to the pathway!
    Argument gene.idtype or cpd.idtype may be wrong.

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/vivianphan/Desktop/BIMM143/bimm143_github/class12(11.3.22)

    Info: Writing image file hsa05310.pathview.png

![Asthma pathway with our genes colored](hsa05310.png)
