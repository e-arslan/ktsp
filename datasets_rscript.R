#BINARY CLASS DATASETS
{
  
  
  ##"Broad patterns of gene expression revealed"
  ## Affymetrix oligonucleotide Hum6000 array
  ## Alon et all
  ## p=2000 n_0=22(N) n_1=40(T)
  {
    install.packages("plsgenomics")
    library("plsgenomics")
    data(Colon)
    data<-t(Colon$X)
    rownames(data)<-Colon$gene.names
    grp<-Colon$Y
    grp = as.factor(grp)
    
    nms <- matrix(rownames(data),ncol=1)
    rownames(data) <- apply(nms,1,function(x) strsplit(x,"a.")[[1]][2])
  }
  
  ## Molecular classification of cancer: class discovery 
  ## Golub et al 
  ## p=7129 n_0=25(AML) n_1=47(ALL)
  {
    source("http://bioconductor.org/biocLite.R")
    biocLite("golubEsets")
    library('golubEsets')
    
    # Because the ALL samples are either T-cell or B-cell, we may wish to consider the data
    # in 3 classes.
    # NOTE: The 'T.B.cell' is 'NA' for AML samples.
    # By default, we only consider the original two classes (i.e. ALL or AML)
    two_classes <- TRUE
    
    # The training data set.
    data('Golub_Train')
    x <- t(exprs(Golub_Train))
    y <- as.vector(pData(Golub_Train)$ALL.AML)
    
    if(!two_classes) {
      y[pData(Golub_Merge)$T.B.cell == "B-cell"] <- "ALL-B"
      y[pData(Golub_Merge)$T.B.cell == "T-cell"] <- "ALL-T"
    }
    
    golub_train <- list(x = x, y = factor(y))
    
    # The test data set.
    data('Golub_Test')
    x <- t(exprs(Golub_Test))
    y <- as.vector(pData(Golub_Test)$ALL.AML)
    
    if(!two_classes) {
      y[pData(Golub_Test)$T.B.cell == "B-cell"] <- "ALL-B"
      y[pData(Golub_Test)$T.B.cell == "T-cell"] <- "ALL-T"
    }
    
    golub_test <- list(x = x, y = factor(y))
    
    
    # The combined training and test data sets.
    data('Golub_Merge')
    x <- t(exprs(Golub_Merge))
    y <- as.vector(pData(Golub_Merge)$ALL.AML)
    
    if(!two_classes) {
      y[pData(Golub_Merge)$T.B.cell == "B-cell"] <- "ALL-B"
      y[pData(Golub_Merge)$T.B.cell == "T-cell"] <- "ALL-T"
    }
    
    golub <- list(x = x, y = factor(y))
    data<-t(golub$x)
    grp<-golub$y
    
  }
  
  # breast cancer
  # ] M. J. van de Vijver, Y. D. He, L. J. van't Veer, et al., "A
  # gene-expression signature as a predictor of survival in breast
  # cancer," 
  # 295, Good-praognosis : 115, poor-prognosis : 180
  # 25760 transcripts
  {
  biocLite("cancerdata")  
  library('cancerdata')
  data(VIJVER)
  data = exprs(VIJVER)
  grp = as.vector(pData(VIJVER)$X70_genes)
  }
  
  
  ## Prediction of central nervous system embryonal tumour outcome based 
  ## data set C
  ## Pomeroy et al
  ## p=7128 n_0=21(D) n_1=39(S)
  {
    source("http://bioconductor.org/biocLite.R")
    biocLite('stepwiseCM')
    library('stepwiseCM')
    data('CNS', package = 'stepwiseCM')
    
    x <- unname(t(CNS$mrna))
    colnames(x) <- as.character(CNS$gene.name$Name)
    
    # The category is the patient survival within 24 months after treatment.
    # 21 patients died (labeled 0)
    # 39 patients survived (labeled 1)
    y <- factor(CNS$class, levels = c(0, 1), labels = c("died", "survived"))
    
    pomeroy <- list(x = x, y = y)
    data=t(pomeroy$x)
    grp=pomeroy$y
  }
  

  
  
  ## Gene expression correlates of clinical prostate cancer behavior
  ## Singh et al
  ## p= 12600 n_0 = 50 (N) n_1 = 52 (T)
  {
    download.file(url = "http://datam.i2r.a-star.edu.sg/datasets/krbd/ProstateCancer/ProstateCancer.zip", destfile = "singh.zip")
    unzip("singh.zip")
    temp <- read.csv("prostate/prostate_TumorVSNormal_train.data", header = FALSE)
    prostate.df <- data.frame(labels = temp[,ncol(temp)], temp[,-ncol(temp)])
    
    singh <- list(
      x = temp[,-ncol(temp)],
      y = temp[,ncol(temp)]
    )
    
    # Removes downloaded files
    unlink("singh.zip")
    unlink("prostate/", recursive = TRUE)
    unlink("prostate*")
    
    data=t(singh$x)
    grp=singh$y
    
    
  }
  
  
  ## Diffuse large B-cell lymphoma outcome prediction by gene-expression
  ## Shipp et al
  ## p= 7129 n_0 = 58 (DLBCL) n_1 = 19 (FL)
  {
    
    download.file(url = "http://datam.i2r.a-star.edu.sg/datasets/krbd/DLBCL/DLBCL-Harvard.zip", destfile = "shipp.zip")
    unzip("shipp.zip")
    temp <- read.csv("DLBCLTumor.data", header = FALSE)
    
    shipp <- list(
      x = temp[,-ncol(temp)],
      y = temp[,ncol(temp)]
    )
    
    # Removes downloaded files
    unlink("shipp.zip")
    unlink("DLBCL*")
    
    data=t(shipp$x)
    grp=shipp$y
    
    
  }
  
  
  ## Genome-wide expression profiling of human blood reveals biomarkers 
  ## Borovecki et al.
  ## p= 22283 n_0 = 14 (control) n_1 = 17 (symptomatic)
  {
    download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/SeriesMatrix/GSE1751/GSE1751_series_matrix.txt.gz", destfile = "borovecki.txt.gz")
    borovecki <- read.table(gzfile("borovecki.txt.gz"), comment.char = "!",
                            stringsAsFactors = FALSE)
    
    x <- t(unname(data.matrix(borovecki[-1, -1])))
    colnames(x) <- borovecki[-1, 1]
    
    # Affymetrix U133A expression levels for 12 symptomatic and 5 presymptomatic
    # Huntington's disease patients versus 14 healthy controls.
    # We pool the symptomatic and presymptomatic samples into a single class.
    y <- factor(c(rep("symptomatic", 17), rep("control", 14)))
    
    borovecki <- list(x = x, y = y)
    
    file.remove("borovecki.txt.gz")
    
    data=t(borovecki$x)
    grp=borovecki$y
  }
  
  
  
  ## Genomic and transcriptional aberrations linked to breast cancer pathophysiologies.
  ## Chin et al.
  ## p= 22215 n_0 = 43 (negative) n_1 = 75 (positive) 
  {
    source('http://bioconductor.org/biocLite.R')
    biocLite('ArrayExpress')
    library('ArrayExpress')
    
    temp_dir <- tempdir()
    ae_obj <- getAE('E-TABM-158', type='processed', path=temp_dir)
    
    exprs_file <- file.path(temp_dir, ae_obj$processedFiles)
    sdrf_file <- file.path(temp_dir, ae_obj$sdrf)
    
    chin_x <- read.table(exprs_file, header=FALSE, sep="\t",
                         stringsAsFactors=FALSE)
    subject_id <- as.character(chin_x[1, -1])
    x <- unname(t(chin_x[-(1:2), -1]))
    
    sdrf <- read.table(sdrf_file, header=TRUE, sep="\t",
                       stringsAsFactors=TRUE, comment.char="")
    sdrf <- sdrf[sdrf[, 78] == "breastTumorExpression.txt", ]
    sdrf <- sdrf[, c(35, 74)]
    colnames(sdrf) <- c("class_label", "subject_id")
    
    sdrf$subject_id <- as.character(sdrf$subject_id)
    
    # We reorder the expressions and the class labels so that the subject ID's match
    sdrf <- sdrf[order(sdrf$subject_id), ]
    x <- x[order(subject_id), ]
    y <- sdrf$class_label
    
    class(x) <- 'numeric'
    
    chin <- list(x=x, y=y)
    
    data=t(chin$x)
    grp=chin$y
    
  }
  
  ## Prognostic Gene Expression Signatures Can Be Measured in Tissues
  ## Chowdary et al. (2006)
  ## p= 22283 n_0 = 62 (breast) n_1 = 42 (colon) 
  {
    source('http://bioconductor.org/biocLite.R')
    biocLite('GEOquery')
    library('GEOquery')
    geo_obj <- getGEO('GSE3726')
    
    chowdary_x <- exprs(geo_obj[[1]])
    x <- unname(t(chowdary_x))
    colnames(x) <- rownames(chowdary_x)
    y <- factor(as.vector(pData(geo_obj[[1]])$source_name_ch1))
    chowdary <- list(x = x, y = y)
    
    data<-t(chowdary$x)
    grp<-chowdary$y
    
    
  }
  
  
  ## The role of the Wnt-signaling antagonist DKK1 in the development 
  ## Tian et al. (2003)
  ## p= 12625 n_0 = 137 (MRI-lytic-lesion sample) n_1 = 36 (MRI-no-lytic-lesion sample) 
  {
    source('http://bioconductor.org/biocLite.R')
    biocLite('GEOquery')
    library('GEOquery')
    geo_obj <- getGEO('GSE755')
    
    tian_x <- exprs(geo_obj[[1]])
    
    x <- unname(t(tian_x))
    colnames(x) <- rownames(tian_x)
    y <- factor(as.vector(pData(geo_obj[[1]])$`description.2`))
    tian <- list(x = x, y = y)
    
    data <- t(tian$x)
    grp <- tian$y
    
  }
  
  
  ## A Marfan syndrome gene expression phenotype in cultured
  ## Yao et al. (2007)
  ## p= 4132 n_0 = 41 (cultured skin fibroblasts from control subjects) 
  ## n_1 = 60 (cultured skin fibroblasts from Marfan subjects) 
  {
    source('http://bioconductor.org/biocLite.R')
    biocLite('GEOquery')
    library('GEOquery')
    geo_obj <- getGEO('GSE8759')
    yao_x <- exprs(geo_obj[[1]])
    x <- unname(t(yao_x))
    colnames(x) <- rownames(yao_x)
    y <- factor(as.vector(pData(geo_obj[[1]])$source_name_ch1))
    yao <- list(x = x, y = y)
    
    data <- t(yao$x)
    grp <- yao$y
    
  }
  
  

  
  
  ## there are missing values
  ## DLBCL
  ## distinct types of diffuse large B-cell
  ## Alizadeh et al.
  ## p= 4026 n_0 = 24 (germinal centre B-like) n_1 = 23 (activated B-like) 
  {
    
    
    download.file(url = "http://datam.i2r.a-star.edu.sg/datasets/krbd/DLBCL/DLBCL-Stanford.zip", destfile = "Alizadeh.zip")
    unzip("Alizadeh.zip")
    temp <- read.csv("DLBCL-Stanford.data", header = FALSE)
    prostate.df <- data.frame(labels = temp[,ncol(temp)], temp[,-ncol(temp)])
    
    Alizadeh <- list(
      x = temp[,-ncol(temp)],
      y = temp[,ncol(temp)]
    )
    
    # Removes downloaded files
    unlink("Alizadeh.zip")
    
    
    data=t(Alizadeh$x)
    grp=Alizadeh$y
    
  }
  
  

  
  ## Performance comparison of two microarray platforms to assess differential 
  ## gene expression in human monocyte and macrophage cells.
  ## Maouche S et al.
  ## p= 26496 n_0 = 49(Human macrophage labeled with Cy5) n_1 = 47(Human monocyte labeled with Cy5)   
  {
    source('http://bioconductor.org/biocLite.R')
    biocLite('GEOquery')
    library('GEOquery')
    geo_obj <- getGEO('GSE10220')
    maouch_x <- exprs(geo_obj[[1]])
    x <- unname(t(maouch_x))
    colnames(x) <- rownames(maouch_x)
    y <- factor(as.vector(pData(geo_obj[[1]])$source_name_ch1))
    maouch <- list(x = x, y = y)
    
    data = t(maouch$x)
    grp = maouch$y
    
  }
  
  
  

  
  
  
  
  source("https://bioconductor.org/biocLite.R")
  biocLite("GSBenchMark")
  library(GSBenchMark)
  
  data(diracpathways)
  
  data(GSBenchMarkDatasets)
  print(GSBenchMark.Dataset.names)
  
  # 1st "leukemia_GSEA" = Armstrongdata
  {
    data(list=GSBenchMark.Dataset.names[[1]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  # 2nd "marfan_GDS2960" = Yao et al.(2007)
  {
    data(list=GSBenchMark.Dataset.names[[2]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  
  # 3rd 
  # Down-regulation of the interferon signaling pathway in T
  # lymphocytes from patients with metastatic melanoma
  # Critchley-Thorne (2007)
  # p = 20844 n_0 =23(Normal) n_1 = 23(Metastasis)
  {
    data(list=GSBenchMark.Dataset.names[[3]])
    data = exprsdata
    grp = phenotypes
  }
  
  # 4th
  # Molecular markers of early Parkinson's disease 
  # based on gene expression in blood
  # Scherzer (2007)
  # p = 22283 n_0 =22(Normal) n_1 = 50(Pakinsons)
  {
    data(list=GSBenchMark.Dataset.names[[4]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  # 5th
  # Gene expression profiles of prostate cancer reveal involvement 
  # of multiple molecular pathways in the metastatic process.
  # Chandran (2007)
  # p = 12558 n_0 =18(Normal) n_1 = 25(metastasis)
  {
    data(list=GSBenchMark.Dataset.names[[5]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  # 6th
  # Gene expression profiles of prostate cancer reveal involvement 
  # of multiple molecular pathways in the metastatic process.
  # Chandran (2007)
  # p = 12558 n_0 =65(primary) n_1 = 25(metastasis)
  {
    data(list=GSBenchMark.Dataset.names[[6]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  # 7th
  # Gene expression profiles of prostate cancer reveal involvement 
  # of multiple molecular pathways in the metastatic process.
  # Chandran (2007)
  # p = 12619 n_0 =18(normal) n_1 = 65(primary)
  {
    data(list=GSBenchMark.Dataset.names[[7]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  # 8th
  # highly accurate 2-gene classifier for
  # differentiating gastrointestinal
  # price et al. (2007)
  # p= 43931 n_0 = 31 (LMS) n_1 = 37 (GIST)
  {
    data(list=GSBenchMark.Dataset.names[[8]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  
  # 9th
  # Selection and validation of differentially expressed
  # genes in head and neck cancer
  # Kuriakose et al (2004)
  # p= 12625 n_0 = 22 (Normal) n_1 = 22 (HNSCC)
  {
    data(list=GSBenchMark.Dataset.names[[9]])
    data = exprsdata
    grp = phenotypes
  }
  
  
  # 10th
  # A two-gene expression ratio predicts clinical outcome 
  # in breast cancer patients treated with tamoxifen.
  # Ma et al.
  # p = 22575 n_0 = 32 (Responsive) n_1 = 28 (Non-Responsive)
  {
    data(breast_GDS807)
    data = exprsdata
    grp = phenotypes
  }
  
  
  # 11th
  # Gene expression analysis of bipolar disorder reveals downregulation 
  # of the ubiquitin cycle and alterations in synaptic genes.
  # Ryan et al. 
  # p = 22283 n_0 = 31 (Normal) n_1 = 30 (Bipolar)
  {
    data(list=GSBenchMark.Dataset.names[[11]])
    data = exprsdata
    grp = phenotypes
  }
  
}




#MULTICLASS DATASETS
{
  ## Molecular Classification of Crohn's Disease and Ulcerative Colitis
  ## Burczynski et al.
  ## p= 22283 n_0 = 59 (Crohn's disease) n_1 = 42 (normal) n_2 = 26 (ulcerative colitis)
  {
    source('http://bioconductor.org/biocLite.R')
    biocLite('GEOquery')
    library('GEOquery')
    geo_obj <- getGEO('GDS1615')
    burczynski_x <- Table(geo_obj)
    x <- unname(t(data.matrix(burczynski_x[, -c(1:2)])))
    colnames(x) <- burczynski_x[, 1]
    
    y <- Columns(geo_obj)$disease.state
    
    burczynski <- list(x = x, y = y)
    
    data=t(burczynski$x)
    grp=burczynski$y
    
    
    ## to get bi class data
    ind = which(grp == "Crohn's disease" | grp == "normal")
    grp = grp[ind]
    data = data[,ind]
    
  }
  
  # Sun et al. (2006) Glioma Data Set
  # p = 54613
  #astrocytomas      glioblastomas          non-tumor oligodendrogliomas 
  #         26                 81                 23                 50 
  {
  source('http://bioconductor.org/biocLite.R')
  biocLite('GEOquery')
  library('GEOquery')
  geo_obj <- getGEO('GDS1962')
  sun_x <- Table(geo_obj)
  x <- unname(t(data.matrix(sun_x[, -c(1:2)])))
  
  # The last columns of 'x' are entirely NA. We omit these columns.
  x <- x[, seq_len(54613)]
  colnames(x) <- sun_x[seq_len(54613), 1]
  y <- Columns(geo_obj)$disease.state
  sun <- list(x = x, y = y)
  data = t(sun$x)
  grp = sun$y
  }
  
  
}




# Lung cancer
#  A. Bhattacharjee, W. G. Richards, J. Staunton, et al., "Classification
# of human lung carcinomas by mRNA expression
# profiling reveals distinct adenocarcinoma subclasses,"
# 186, Adenocarcinomas : 139, non-adenocarcinomas : 47
# 12600 transcripts
http://www.pnas.org/content/98/24/13790/suppl/DC1


# prostate cancer
#  B.-L. Adam, Y. Qu, J. W. Davis, et al., "Serum protein
# fingerprinting coupled with a pattern-matching algorithm
# distinguishes prostate cancer from benign prostate hyperplasia
# and healthy men,"
# 326, cancer : 167, non-cancer : 159
# 45000 m/z(mass over charge) values

