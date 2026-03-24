# Cancer infiltration rates

This dataset features cancer infiltration rates and microsatellites
data.

## Usage

``` r
TxTum
```

## Format

A data frame with 106 rows and 60 variables.

- `CELTUMCO`:

  a numeric vector

- `age`:

  a numeric vector

- `sexe`:

  a numeric vector

- `HISTOADK`:

  a numeric vector

- `H2`:

  a numeric vector

- `P3`:

  a numeric vector

- `P4`:

  a numeric vector

- `E1`:

  a numeric vector

- `P5`:

  a numeric vector

- `R10`:

  a numeric vector

- `C3M`:

  a numeric vector

- `P6`:

  a numeric vector

- `RB`:

  a numeric vector

- `FL7A`:

  a numeric vector

- `P53`:

  a numeric vector

- `W2`:

  a numeric vector

- `P2`:

  a numeric vector

- `P1`:

  a numeric vector

- `W4`:

  a numeric vector

- `MT1`:

  a numeric vector

- `MT2`:

  a numeric vector

- `MT4`:

  a numeric vector

- `MT3`:

  a numeric vector

- `HLA`:

  a numeric vector

- `HLD`:

  a numeric vector

- `HLC`:

  a numeric vector

- `HLB`:

  a numeric vector

- `EA1`:

  a numeric vector

- `EA3`:

  a numeric vector

- `EA2`:

  a numeric vector

- `EA4`:

  a numeric vector

- `EB1`:

  a numeric vector

- `EB2`:

  a numeric vector

- `EB3`:

  a numeric vector

- `EB4`:

  a numeric vector

- `EGF1`:

  a numeric vector

- `EGF2`:

  a numeric vector

- `EGF3`:

  a numeric vector

- `EGF4`:

  a numeric vector

- `EGF5`:

  a numeric vector

- `EGF6`:

  a numeric vector

- `FL7B`:

  a numeric vector

- `VSFGF7`:

  a numeric vector

- `F3A`:

  a numeric vector

- `F3B`:

  a numeric vector

- `VSFGFR3`:

  a numeric vector

- `F4`:

  a numeric vector

- `Q5`:

  a numeric vector

- `VSTOP1`:

  a numeric vector

- `VSTOP2A`:

  a numeric vector

- `VSEGFR`:

  a numeric vector

- `AFRAEGFR`:

  a numeric vector

- `SRXRA`:

  a numeric vector

- `SMT`:

  a numeric vector

- `QMTAMPN`:

  a numeric vector

- `QMTDELN`:

  a numeric vector

- `SHL`:

  a numeric vector

- `SEA`:

  a numeric vector

- `SEB`:

  a numeric vector

- `QPCRFGF7`:

  a numeric vector

## References

Régression Bêta PLS. (French) \[PLS Beta regression.\], F. Bertrand, N.
Meyer, M. Beau-Faller, K. El Bayed, N. Izzie-J., M. Maumy-Bertrand,
(2013), J. SFdS, 154(3):143-159  

Partial Least Squares Regression for Beta Regression Models. F.
Bertrand, M. Maumy (2021). useR! 2021, Zurich.  

## Examples

``` r
data
#> function (..., list = character(), package = NULL, lib.loc = NULL, 
#>     verbose = getOption("verbose"), envir = .GlobalEnv, overwrite = TRUE) 
#> {
#>     fileExt <- function(x) {
#>         db <- grepl("\\.[^.]+\\.(gz|bz2|xz)$", x)
#>         ans <- sub(".*\\.", "", x)
#>         ans[db] <- sub(".*\\.([^.]+\\.)(gz|bz2|xz)$", "\\1\\2", 
#>             x[db])
#>         ans
#>     }
#>     my_read_table <- function(...) {
#>         lcc <- Sys.getlocale("LC_COLLATE")
#>         on.exit(Sys.setlocale("LC_COLLATE", lcc))
#>         Sys.setlocale("LC_COLLATE", "C")
#>         read.table(...)
#>     }
#>     stopifnot(is.character(list))
#>     names <- c(as.character(substitute(list(...))[-1L]), list)
#>     if (!is.null(package)) {
#>         if (!is.character(package)) 
#>             stop("'package' must be a character vector or NULL")
#>     }
#>     paths <- find.package(package, lib.loc, verbose = verbose)
#>     if (is.null(lib.loc)) 
#>         paths <- c(path.package(package, TRUE), if (!length(package)) getwd(), 
#>             paths)
#>     paths <- unique(normalizePath(paths[file.exists(paths)]))
#>     paths <- paths[dir.exists(file.path(paths, "data"))]
#>     dataExts <- tools:::.make_file_exts("data")
#>     if (length(names) == 0L) {
#>         db <- matrix(character(), nrow = 0L, ncol = 4L)
#>         for (path in paths) {
#>             entries <- NULL
#>             packageName <- if (file_test("-f", file.path(path, 
#>                 "DESCRIPTION"))) 
#>                 basename(path)
#>             else "."
#>             if (file_test("-f", INDEX <- file.path(path, "Meta", 
#>                 "data.rds"))) {
#>                 entries <- readRDS(INDEX)
#>             }
#>             else {
#>                 dataDir <- file.path(path, "data")
#>                 entries <- tools::list_files_with_type(dataDir, 
#>                   "data")
#>                 if (length(entries)) {
#>                   entries <- unique(tools::file_path_sans_ext(basename(entries)))
#>                   entries <- cbind(entries, "")
#>                 }
#>             }
#>             if (NROW(entries)) {
#>                 if (is.matrix(entries) && ncol(entries) == 2L) 
#>                   db <- rbind(db, cbind(packageName, dirname(path), 
#>                     entries))
#>                 else warning(gettextf("data index for package %s is invalid and will be ignored", 
#>                   sQuote(packageName)), domain = NA, call. = FALSE)
#>             }
#>         }
#>         colnames(db) <- c("Package", "LibPath", "Item", "Title")
#>         footer <- if (missing(package)) 
#>             paste0("Use ", sQuote(paste("data(package =", ".packages(all.available = TRUE))")), 
#>                 "\n", "to list the data sets in all *available* packages.")
#>         else NULL
#>         y <- list(title = "Data sets", header = NULL, results = db, 
#>             footer = footer)
#>         class(y) <- "packageIQR"
#>         return(y)
#>     }
#>     paths <- file.path(paths, "data")
#>     for (name in names) {
#>         found <- FALSE
#>         for (p in paths) {
#>             tmp_env <- if (overwrite) 
#>                 envir
#>             else new.env()
#>             if (file_test("-f", file.path(p, "Rdata.rds"))) {
#>                 rds <- readRDS(file.path(p, "Rdata.rds"))
#>                 if (name %in% names(rds)) {
#>                   found <- TRUE
#>                   if (verbose) 
#>                     message(sprintf("name=%s:\t found in Rdata.rds", 
#>                       name), domain = NA)
#>                   thispkg <- sub(".*/([^/]*)/data$", "\\1", p)
#>                   thispkg <- sub("_.*$", "", thispkg)
#>                   thispkg <- paste0("package:", thispkg)
#>                   objs <- rds[[name]]
#>                   lazyLoad(file.path(p, "Rdata"), envir = tmp_env, 
#>                     filter = function(x) x %in% objs)
#>                   break
#>                 }
#>                 else if (verbose) 
#>                   message(sprintf("name=%s:\t NOT found in names() of Rdata.rds, i.e.,\n\t%s\n", 
#>                     name, paste(names(rds), collapse = ",")), 
#>                     domain = NA)
#>             }
#>             files <- list.files(p, full.names = TRUE)
#>             files <- files[grep(name, files, fixed = TRUE)]
#>             if (length(files) > 1L) {
#>                 o <- match(fileExt(files), dataExts, nomatch = 100L)
#>                 paths0 <- dirname(files)
#>                 paths0 <- factor(paths0, levels = unique(paths0))
#>                 files <- files[order(paths0, o)]
#>             }
#>             if (length(files)) {
#>                 for (file in files) {
#>                   if (verbose) 
#>                     message("name=", name, ":\t file= ...", .Platform$file.sep, 
#>                       basename(file), "::\t", appendLF = FALSE, 
#>                       domain = NA)
#>                   ext <- fileExt(file)
#>                   if (basename(file) != paste0(name, ".", ext)) 
#>                     found <- FALSE
#>                   else {
#>                     found <- TRUE
#>                     switch(ext, R = , r = {
#>                       library("utils")
#>                       sys.source(file, chdir = TRUE, envir = tmp_env)
#>                     }, RData = , rdata = , rda = load(file, envir = tmp_env), 
#>                       TXT = , txt = , tab = , tab.gz = , tab.bz2 = , 
#>                       tab.xz = , txt.gz = , txt.bz2 = , txt.xz = assign(name, 
#>                         my_read_table(file, header = TRUE, as.is = FALSE), 
#>                         envir = tmp_env), CSV = , csv = , csv.gz = , 
#>                       csv.bz2 = , csv.xz = assign(name, my_read_table(file, 
#>                         header = TRUE, sep = ";", as.is = FALSE), 
#>                         envir = tmp_env), found <- FALSE)
#>                   }
#>                   if (found) 
#>                     break
#>                 }
#>                 if (verbose) 
#>                   message(if (!found) 
#>                     "*NOT* ", "found", domain = NA)
#>             }
#>             if (found) 
#>                 break
#>         }
#>         if (!found) {
#>             warning(gettextf("data set %s not found", sQuote(name)), 
#>                 domain = NA)
#>         }
#>         else if (!overwrite) {
#>             for (o in ls(envir = tmp_env, all.names = TRUE)) {
#>                 if (exists(o, envir = envir, inherits = FALSE)) 
#>                   warning(gettextf("an object named %s already exists and will not be overwritten", 
#>                     sQuote(o)))
#>                 else assign(o, get(o, envir = tmp_env, inherits = FALSE), 
#>                   envir = envir)
#>             }
#>             rm(tmp_env)
#>         }
#>     }
#>     invisible(names)
#> }
#> <bytecode: 0x12bd7bd60>
#> <environment: namespace:utils>
print(TxTum)
#>     CELTUMCO age sexe HISTOADK H2 P3 P4 E1 P5 R10 C3M P6 RB FL7A P53 W2 P2 P1
#> 1         80  59    1        1 NA  1 NA  1 NA   1   1  1 NA    1  NA  1  1  0
#> 2         99  67    1        1  1  1  1  1  1   1  NA  1  1    0   1  0  0  0
#> 3         20  60    1        0  0  1  0  0 NA   1  NA  0  0    0   1  0  0  0
#> 4         90  65    1        1  1  1  1  1  0  NA   0  0  0    1   1  1  0  0
#> 5         50  43    1        1  1  1  1 NA NA  NA   1  0  0   NA   1  1  1  0
#> 6         80  67    0        1  1 NA NA NA  0   1   1 NA  1    1   1  0  0  0
#> 7         40  53    1        1  1  1  0  1  1   1   1  0  1    1   0  0  0 NA
#> 8         70  66    1        0  0  0  0  1  1  NA   1 NA  1    1   1  1  1  0
#> 9         60  71    1        1  1  1  1 NA NA   1   1  1  1   NA   1 NA  1  1
#> 10        90  65    1        1  1  1  1  1  0   1   1  1  1    1   1  1  1  1
#> 11        40  77    0        0  0  0 NA NA  0   1   1 NA NA    1   0 NA NA  0
#> 12        70  62    1       NA  1 NA  1  1  1  NA   1  1  1    1   1 NA  0  1
#> 13        40  76    0        1 NA  1 NA  1  0  NA   0  0  0    0   1 NA  1  0
#> 14        80  51    1        0  1  0  1  1 NA   1   1  1  1   NA   0 NA  0  0
#> 15        30  73    1       NA  1  1  0 NA  0   0   1 NA NA   NA   1  0  0  1
#> 16        90  69    0        0  0  0  0  1  1   1  NA  1  1    1   1  0  0  0
#> 17        30  65    1        1  1  1  1  1  1   0   1  1  1    0  NA  1  1  0
#> 18        30  62    1        1  1 NA  1  1 NA   1   1  1  1    1   1 NA  1 NA
#> 19        30  41    1        1  1  1 NA  0 NA   1   1  0  1    0   1  0  0  0
#> 20        65  80    1        1  1  1  0  1  0   1  NA  0  1   NA  NA NA NA  0
#> 21        45  61    1       NA NA  1  1  0 NA   1  NA  1  1    1   1 NA  1 NA
#> 22        40  74    1        1 NA NA  1  1  0   1   1  0 NA    0   1 NA NA  1
#> 23        20  56    1        0  0  1 NA  1  0  NA   0 NA  1    0   1  1  0  0
#> 24        40  45    1       NA  1  1  1 NA  0   1   0 NA  0    0   1 NA  0 NA
#> 25        70  59    0        0  1  1  1 NA  0   1   1  1 NA    1  NA  1  1  0
#> 26        70  62    0        0  0  0  0  1  1   1  NA  1  1    1   1  1  1  1
#> 27        40  52    1        1  1  1  1 NA NA  NA   1  1  0   NA   1  1  1  0
#> 28        50  75    1        1  1  1  1  0  1   0   0  0  1    0   1  1 NA  1
#> 29        80  55    1        1 NA NA NA  1  0   1  NA  1  1    1   1 NA  0  0
#> 30        50  51    1        1  1  1  1 NA  0   1   1  1  1    0   1  1  1  0
#> 31        45  68    1        1  1  1  1  1  0   1  NA NA  0    1   1  0 NA  0
#> 32        90  64    1        1  1  1 NA  1  1   1   1 NA NA    1   1  1  1  0
#> 33        90  69    1        1 NA  1 NA  1 NA   1  NA  1  1   NA  NA  1  1  0
#> 34        50  68    0        1 NA  1  1  1  1   1   1  1  1    0   1  1  1  0
#> 35        60  70    1        0  0  0 NA  0  0   1   1  1  1    1   1 NA  1  1
#> 36        70  49    1        1 NA NA  1  0  1  NA   0  0  1    1   1  0  0  0
#> 37        50  54    1        1 NA  1  1 NA  1   1   1  0  1    1   1  0  0  1
#> 38        20  61    0        0  0 NA  0  1  0   1  NA NA  1    0   0  1  1  0
#> 39        70  44    1        1 NA  1  1  0  1   0   1  1  1    0   1 NA  1  0
#> 40        60  73    1        1 NA NA  1  1  0   1   1  1 NA    1   1  1  1  1
#> 41        20  51    1        1  1  0  0  0  0   0   1  1  0   NA  NA  0  0  1
#> 42        50  60    1        1  1 NA  1  1 NA  NA  NA  1  1   NA   1 NA  1 NA
#> 43        50  65    1        1 NA  1  1  0  0   0  NA  1  1    1   1 NA  1  0
#> 44        30  60    1        1  1  1  0  1  0   1  NA  1  0    1   1  0 NA  0
#> 45        40  49    1        0 NA  0 NA  0  1   1   1  1  1    1   1  1  1  0
#> 46        50  44    0       NA  0  0  0  1  0   0   0  0 NA    0   1  0  0  1
#> 47        60  44    1        1  1  1  1  1  0   1   1  1 NA    1   1  1  1  1
#> 48        90  64    1        1 NA NA  1  1  0   1   1  1  1    1   1  0  0  1
#> 49        70  52    1        1  1  1  1  1  0   1  NA  1 NA    1   1  1 NA NA
#> 50        30  54    0        0  0 NA  1 NA  1  NA   0  0  0    1   0  1  1 NA
#> 51        80  64    1        1 NA  1  1  1  1   1  NA  1  1    0   1 NA  0 NA
#> 52        30  76    1        0  0  0  0  1  0   1  NA  1  1    0   1  0  1 NA
#> 53        50  70    1        1  0  1 NA  1  0   1   1  1  1    0  NA NA  0  0
#> 54        90  51    1        1 NA  1  1  1  0   1  NA  1  1    1   1  1  1  1
#> 55        50  60    1        1 NA  1 NA  1  1   1   1  1  1    0   1  1  1  0
#> 56        40  76    1        1  1  1  1  0  1   0  NA  1  1    1   1  0  0 NA
#> 57        40  64    0       NA NA NA  0  1  1   1   1  1  0    1   0  1  1  1
#> 58        30  72    1        1  1  0 NA  0  1   1   1 NA  1   NA   1  1 NA  0
#> 59        50  52    1        1  1 NA  1  1  0  NA   1  1  0    0   1  1  1  0
#> 60        50  75    1        1  1  1 NA  1  1   1   0 NA  1    1   0 NA  0  0
#> 61        80  74    0       NA NA NA  1 NA  0   1  NA  1  1    1   1  0  0  0
#> 62        20  52    1        1  1  1 NA  1 NA  NA  NA  1 NA    1   1  0  1  0
#> 63        80  71    1        1  1  1  1  1 NA   1   1  1  1    1   1 NA  1 NA
#> 64        40  48    1        1  0  0 NA  1  0   0   1  1  1   NA   1  0  0  0
#> 65        40  70    1        1 NA  1  1  1  1   1  NA  1  1   NA  NA  1  1  0
#> 66        30  55    1        1 NA  1  0  1 NA   1  NA  0  1   NA   1 NA  0  1
#> 67        60  53    1       NA  1 NA  1  1 NA   0   0  0  1    0   1  1 NA  0
#> 68        20  43    1       NA  1  0  1  0  0  NA   0  0  0   NA  NA NA  0  0
#> 69        20  61    1        1  1 NA  0  0 NA  NA   0  0  1    1  NA  0  0  0
#> 70        90  68    1       NA  1  1  1 NA  0   1  NA  1 NA    1   1 NA  0 NA
#> 71        99  62    1        1 NA  1  1  1  0   0   0  0  1    1   1  0  0  0
#> 72        60  79    1        0 NA  0  0 NA  0   0   1  0  0   NA   0  1  1  0
#> 73        70  52    0        0  0  0 NA NA  1   0   0  0  1   NA   1  0  0  1
#> 74        70  67    1        1  1 NA NA  1  1   1  NA  1  1   NA   1  0  0  1
#> 75        30  61    1        1  1  1 NA NA  1   0  NA  0  0    1   1  1  0  1
#> 76        20  63    0       NA  1  1  1 NA  0   1   1  1  1    0   1  0 NA NA
#> 77        30  50    1        1  1 NA NA  1 NA  NA   1  1  1    1   1  1 NA  1
#> 78        30  68    1        0  0  0  1  0  1   0   0 NA  0    0   1 NA  0 NA
#> 79        30  66    1        1  1  1  1  1 NA  NA  NA  0  1   NA  NA  0  0  0
#> 80        35  41    0        1 NA  1  1  1  1   1   0  0  1   NA   1 NA  0 NA
#> 81        50  58    0       NA  1  1  1  1 NA  NA   1  1  1   NA   1 NA NA  1
#> 82        20  50    1        0  1  0  0  0  1  NA  NA  1  1   NA   0  0  0  0
#> 83        40  70    1        0  0  0  0 NA NA   0   0  0  0    0   1 NA  0  1
#> 84        50  63    1       NA  1 NA NA  1 NA   1  NA NA NA   NA   1  0  0  0
#> 85        30  70    1        1  1  1  1  0  0   0  NA  0 NA    1   0  0  0  1
#> 86        50  79    1        1  1  1 NA  1 NA  NA   1  1  0    1   1  0  0  0
#> 87        50  52    1        1  1  1 NA  1  0   1  NA  1 NA    0   1 NA  1  1
#> 88        50  60    1        1  1  1  1  1  0   1   1  0  0    1  NA  0  1 NA
#> 89        50  71    1        1 NA NA  1  1  1   1   1  1 NA    1   1 NA  0  1
#> 90        50  66    1        1 NA  1  1  1  1  NA   0  0  1   NA   0  0 NA  1
#> 91        20  62    1       NA  1  1  0 NA NA  NA   0  0  1    0   1  0  0 NA
#> 92        50  77    1        0 NA  0 NA  0  1   1   1  1  1    0   1  1  1  1
#> 93        60  73    1        1  1  1  1 NA  0   1  NA  0  1   NA   1  1 NA  1
#> 94        50  55    1        1 NA NA  1 NA NA   0   1  1  1   NA   1  1 NA  1
#> 95        30  53    0        1  1  1  1  1  1   1  NA  1 NA    1   0  0  1  0
#> 96        30  73    1        1 NA NA  1 NA  1  NA  NA  1  0   NA   1  1  1 NA
#> 97        50  49    0        1 NA  1  1  1  1  NA  NA  1  1    1   1 NA  1  1
#> 98        50  44    0        1 NA  1  1  0 NA   1   1  1  1    0   1  1  1  0
#> 99        30  65    0        1  1  1  1  0  0   1   1  1  1   NA   0  1  1  1
#> 100       80  52    1        1  1  1  1  1  1  NA   1  1  1    1   1  0  0  1
#> 101       40  67    1        1  1  1 NA  1  0   1   1  1  1    1  NA  1  1 NA
#> 102       30  67    1        0  0 NA  0  1  0   1  NA  0  1   NA   0  1  1  1
#> 103       30  55    0        0 NA  0  0  1  1  NA   1 NA  1   NA  NA NA  0  1
#> 104       30  50    1        1  1  1  1  1 NA  NA   1 NA  1    0   1  0  0 NA
#> 105       20  50    1        1  1  1  0  0  0   1  NA  1  0    1   1  1  1  0
#> 106       30  74    0        0  0  0  1  1 NA   0  NA  0 NA    1   0 NA  0  0
#>     W4 MT1 MT2 MT4 MT3 HLA HLD HLC HLB EA1 EA3 EA2 EA4 EB1 EB2 EB3 EB4 EGF1
#> 1    0  NA   0   0   1  NA   0   0   0   1   1   1   1  NA   1   1   1    1
#> 2   NA   1  NA   1   1   1   1   1   1   0   0  NA   0  NA  NA  NA   1    1
#> 3    0   0   0   0  NA   0   0   0  NA  NA  NA   0   0   1   1  NA   1   NA
#> 4   NA  NA  NA   0   0   0   0   0  NA   0  NA  NA   0   0   0   0   0   NA
#> 5    0   0  NA   0  NA   0   0  NA  NA   0   0   0  NA   0   0   0  NA   NA
#> 6   NA  NA   1   1  NA   1   1   1   1  NA   0  NA  NA  NA  NA   0   1   NA
#> 7    0  NA   1   1   0   1   1   1   1   0   0   0  NA   1  NA   1  NA    1
#> 8    0   1   1   1  NA  NA   1  NA   1   1   1   1   1  NA   1   1   1   NA
#> 9   NA  NA   1   1   1   0  NA  NA  NA   1   1   1   1  NA   1   1  NA   NA
#> 10   1   0   0   0  NA   0  NA  NA   0   1  NA   1   1  NA   0   0   0    1
#> 11  NA  NA   0   0   0   0   0  NA   0   0   0  NA  NA  NA   1   1  NA   NA
#> 12   0  NA  NA   0   0  NA   0   1   0  NA  NA   0   0   1  NA   1  NA    0
#> 13   0   1  NA   1   1   1   1   1  NA   1  NA   1  NA  NA  NA   1  NA   NA
#> 14   0   1  NA  NA   1   1   1   1   1  NA  NA   1   1   1   1   1   1   NA
#> 15  NA   0  NA  NA  NA   0  NA   0   0   0   0   0   0  NA  NA  NA   1    0
#> 16   0   1  NA   1  NA   1   1   0   1  NA   1  NA   1  NA   1   1   1    1
#> 17   0   0  NA   0  NA   0  NA  NA   0   0  NA   0  NA   0   0   0  NA   NA
#> 18   0   1  NA  NA   1  NA  NA  NA   0   1   1   1   1   0   0   0   0   NA
#> 19   0   0   0  NA  NA   0   0  NA   0  NA  NA  NA  NA   0   0  NA   0    0
#> 20   1   1  NA  NA   1   1   0  NA   0   1   1   1   1   1   1   1   1    1
#> 21   1   0  NA   0  NA   0   0  NA  NA  NA   0   0  NA   0   0  NA   0    0
#> 22  NA   0   0   0   0   0   0   0  NA  NA   0  NA   0   0  NA  NA   0   NA
#> 23   0   0   0   0  NA   0  NA   0   0   0   0   0   0  NA   0   0   0    1
#> 24   1   1   1   1   1  NA   1  NA   1   0   0   0  NA   0   0   0   0   NA
#> 25  NA   1   1   1  NA  NA   1   1   1   1   1   1  NA   1   1  NA  NA   NA
#> 26   1   0   0   0  NA   0  NA   0  NA   1   1   1   1  NA  NA   1   1   NA
#> 27   0   0  NA   0  NA   0   0  NA  NA   0   0   0  NA   0   0   1  NA   NA
#> 28   0  NA   1   1  NA   1   1   1   1  NA   0   0   0  NA   1   1   1   NA
#> 29   0   1   1   1   1   1   1   1  NA   1   1  NA  NA   1   1   1   1   NA
#> 30  NA   1   1   1  NA   1  NA  NA   1  NA   1   1   1   1   1  NA   1    1
#> 31  NA   1   1  NA   1  NA  NA  NA   1  NA  NA   1   1   1   1   1   1    1
#> 32  NA  NA   0   0  NA  NA   0   0   0  NA   1   1   1   1   1   1   1   NA
#> 33   0   1   1   1  NA  NA  NA  NA   1  NA   1   1   1  NA   1  NA   1    1
#> 34   0   1   1   1  NA   1   1   1  NA   1  NA   1   1   1  NA   1   1    1
#> 35  NA  NA  NA  NA   1   1   1  NA   1   0   0   0   0   0   0  NA  NA   NA
#> 36   0  NA   0   0   0   0  NA  NA   0   1   1  NA   1   1  NA  NA   1    1
#> 37   0  NA  NA   0  NA   0   1  NA  NA  NA   0  NA   0   0   0  NA   0    1
#> 38  NA   1   1   1   1  NA   1   1   1  NA   0   0   0   0  NA  NA   0    0
#> 39   0   1   1  NA  NA  NA  NA  NA   0   1  NA   1   1   0   0   0  NA    0
#> 40   0   1   1  NA   0  NA  NA  NA   0  NA   0   0  NA   0  NA   0   0    0
#> 41  NA  NA   0  NA   0   0   0  NA   0  NA   0   0   0   0   0  NA   0    0
#> 42   1  NA  NA   0  NA   1   1  NA   1   1   1   1  NA  NA  NA   0   0   NA
#> 43  NA  NA   1   1  NA   1  NA  NA   1   1   1   0  NA   1  NA  NA   1    1
#> 44   0   0  NA  NA   0   0   0  NA  NA   0   0   0   0  NA   0  NA  NA   NA
#> 45   0   0   0   0   0   0   0  NA   0   0  NA   0   1   1  NA   1  NA    1
#> 46  NA   1  NA   1  NA   1  NA  NA   1  NA  NA   0  NA   0   0   0   0   NA
#> 47   1   0   0   0  NA  NA   0  NA  NA  NA  NA   1  NA   0  NA   0   0   NA
#> 48  NA   0   0  NA   0   0   0  NA   0   1   0   1   1  NA   1   1   1    0
#> 49   1  NA   1   1   1  NA  NA  NA  NA   1  NA   1  NA   1  NA  NA   1    1
#> 50  NA   1   1   1  NA  NA   0   0   1   1   1   1   1  NA   1   1   1    0
#> 51   0  NA  NA  NA   0   0   0  NA   0  NA  NA  NA   1  NA  NA   1   1    0
#> 52   0   0  NA  NA   0   0   0  NA  NA   1   1   1   1   0  NA   1   1   NA
#> 53   0  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA    1
#> 54   1  NA   1   1  NA   1   1  NA   1   1   1   1   1  NA  NA  NA  NA    0
#> 55   0   0   1   0   0  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA    0
#> 56  NA   1   1   1  NA   1   1  NA  NA  NA   0   0   0   0   0   0   0    1
#> 57  NA   1  NA   1  NA   0  NA   0   0  NA   1   1   1   1   1  NA  NA   NA
#> 58  NA   0   0   0   0  NA   0  NA   0  NA  NA  NA   1   0   0   0   0    1
#> 59   0  NA   0   0  NA   0   0  NA   0  NA   0   0   0   1  NA   1  NA    1
#> 60  NA  NA   0   0   0   0   0  NA   0   0   0   0   0   1   1   1  NA    0
#> 61  NA   0  NA   0   0   0   0   0   0  NA  NA   1   1   0   0   0  NA    0
#> 62   0   1   1   1  NA   1   1  NA   1  NA   1   1   1   1   1   1  NA   NA
#> 63  NA  NA   1   1  NA   1   1  NA   1  NA  NA  NA  NA  NA   1   1   1   NA
#> 64  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA   1   1   1   1   NA
#> 65   0   0   0   0   0   0   0  NA   0  NA  NA   1  NA   0   0   0  NA    0
#> 66   0   1   1  NA   1  NA  NA  NA   1   0   0   0   0   0  NA   0   0    0
#> 67  NA  NA   1   1  NA   1   1  NA   1   0   0   0  NA   0   1   1   1   NA
#> 68   0   0   0   0   0   0   0  NA   0   1  NA  NA  NA  NA  NA   0   1    0
#> 69   0   0   0   0  NA   0  NA  NA   0  NA   0  NA   0   0   0  NA  NA    0
#> 70   1   0  NA  NA   1   1   1  NA  NA   1   1   1  NA   0   0   0   0   NA
#> 71  NA   1   1   1  NA   1   1  NA   1   1   1   1   1   1   1  NA  NA   NA
#> 72   0  NA  NA  NA  NA   0   0  NA   0  NA   0   0  NA   0   0  NA   0    0
#> 73   1   1  NA  NA  NA   1   1   0   1  NA   0  NA   0  NA   0   0  NA   NA
#> 74   1   1  NA   1   1  NA  NA   1   1  NA   1  NA   0   1   1   1   1   NA
#> 75   1   1  NA  NA  NA  NA   1  NA   0   0   0   0   0   0  NA   0   0    1
#> 76   0  NA  NA  NA  NA  NA  NA  NA   0  NA   0   0   0   0   0   0   0    0
#> 77  NA   1   1   1  NA   0   0   0   1   1  NA   1   1   1   1   1  NA    0
#> 78   1  NA  NA   0  NA   0   0   0   0   1   1   1   1   0   0  NA   0   NA
#> 79   0  NA   0  NA   1   1   0  NA   1   0   0   0  NA  NA  NA   1  NA    1
#> 80  NA  NA   0   0   0  NA   0   0   0   1  NA  NA   1  NA   0   0   0   NA
#> 81   1  NA  NA  NA   0  NA  NA   1   1   1  NA   1   1  NA  NA   1  NA    1
#> 82   0   1  NA   1  NA   1   1  NA   1  NA   0  NA  NA   0   0   0   0   NA
#> 83   1   0   0   0  NA  NA   0  NA  NA   1  NA  NA   1   0   0   0   0    0
#> 84  NA   0   0  NA  NA  NA  NA  NA  NA   1  NA  NA  NA  NA   1   1   1   NA
#> 85   0   1  NA   1  NA  NA   1  NA   1   1   1   1   1   0   0   0  NA    1
#> 86   0  NA  NA   1   1  NA   1   1   1  NA  NA  NA  NA   0   0  NA  NA    1
#> 87   1  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA    1
#> 88   0  NA   0   0   0  NA   0  NA   0  NA   1   1   1  NA   0  NA   0   NA
#> 89   1  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA   NA
#> 90   1  NA  NA  NA  NA   1   1  NA   1  NA   1   1   1   0  NA   0   0    1
#> 91  NA  NA   0   0  NA   1  NA  NA   1   0   0   0  NA   0   0   0   1    0
#> 92  NA   1  NA   0   0   0  NA  NA   0  NA   1   1  NA   1   1   1  NA   NA
#> 93   1   1   1   1   1   1   1  NA  NA  NA  NA   1  NA  NA   1   1   1   NA
#> 94  NA   0   0   0   0   0  NA  NA   0  NA  NA  NA   1   0  NA   0   0    1
#> 95   0  NA   1   1   1   1   1  NA  NA   1  NA  NA  NA  NA  NA   1   1   NA
#> 96   1   1  NA  NA  NA   0   0  NA   0   0   0   0   0   1   1   1   1    1
#> 97  NA   0   0   1   1   0   1   0   0   0   1   1  NA   0   0   0   0    0
#> 98   0   0   0   0   0  NA   1  NA   0   0   0   0   0   0  NA   0   0   NA
#> 99   1  NA   1  NA   1  NA  NA  NA   0  NA   0  NA   0   0  NA   0  NA    1
#> 100  0   0   0   0  NA  NA   1  NA   1  NA   0  NA   0  NA   1   1   0    0
#> 101  0   0   0   0  NA  NA   0   0   0   0  NA  NA   0  NA   1  NA   1    0
#> 102  1   1  NA  NA   1  NA  NA  NA   0   1   1   1  NA   1   1   1   1   NA
#> 103  1  NA  NA   0  NA  NA  NA   0   0   1   1   1   1   1   1   0   1    0
#> 104  0  NA   1   1  NA   1   1  NA   1   1   1   1   1   1   1   1   1   NA
#> 105  0   0   0   0   0   0  NA  NA  NA   1   1   1   1   1   1  NA  NA    1
#> 106  0  NA   1   1   1   1   1  NA   1   1   1  NA   1   1   1  NA  NA    1
#>     EGF2 EGF3 EGF4 EGF5 EGF6 FL7B VSFGF7 F3A F3B VSFGFR3 F4 Q5 VSTOP1 VSTOP2A
#> 1      1   NA    1    1    1   NA      1   1   1       1  0  0      0       1
#> 2      1    1    1   NA   NA    0      0  NA   1       1  0  0      0       0
#> 3      0    0    0    0    0   NA      0   0   0       0 NA  0      0       0
#> 4     NA   NA    1   NA   NA    1      1   1  NA       1 NA  0      0       0
#> 5      0    0   NA   NA   NA    0      0   0   0       0 NA  0      0       1
#> 6      1    1   NA   NA   NA    0      1   1  NA       1  0  0      0       0
#> 7      1   NA    1    1    1    1      1  NA  NA      NA  0  0      0       0
#> 8     NA    1    1    1   NA   NA      1  NA   1       1  0 NA      0       1
#> 9      0    0    0    0    0    1      1  NA   1       1 NA  1      1       1
#> 10     1   NA    1    1    1    1      1   1   1       1  1 NA      1       1
#> 11     0    0    0    0    0    1      1   0   0       0  0  0      0      NA
#> 12     0   NA    0   NA   NA    1      1   1   1       1  0  0      0       0
#> 13     1   NA   NA   NA   NA    0      0   0   0       0  1  1      1       1
#> 14     1   NA    1    1    1   NA     NA   1   1       1 NA NA      0       0
#> 15     0   NA    0    0    0    0      0  NA   1       1  1 NA      1       0
#> 16     1    1    1    0    0    1      1   1   1       1  0  0      0       0
#> 17    NA   NA    0    0    0    0      0   0   0       0  0 NA      0       1
#> 18    NA    0    0    0    0   NA      1   0   0       0  0  0      0       1
#> 19     0    0    0    0    0   NA      1  NA  NA      NA  0 NA      0       0
#> 20     1   NA    1    1    1    0      0   0   0       0 NA  0     NA      NA
#> 21     0    0   NA   NA   NA    1      1   1   1       1 NA  1      1       1
#> 22    NA    0    0    0    0   NA      0  NA   1       1  1 NA      1      NA
#> 23     0    0   NA    0    0    1      1   0   0       0  0  0      0       1
#> 24     1    1    1    1    1    0      0   1   1       1  1 NA      1       0
#> 25    NA    1   NA    1    1    1      1  NA  NA      NA NA  1      1       1
#> 26    NA   NA   NA   NA   NA   NA      1   1   0       1  1  1      1       1
#> 27     0    0   NA   NA   NA   NA     NA   0   1       0 NA  0      0       1
#> 28     1   NA    1   NA   NA    0      0   1   1       1  0 NA      0       1
#> 29     1    1   NA    1    1   NA      1   1   1       1  0  0      0       0
#> 30     1    1   NA    1    1    0      0  NA  NA      NA  0 NA      0       1
#> 31     1    1    0   NA   NA   NA      1  NA   0       0  0  0      0       0
#> 32    NA    0    0    0    0   NA      1   1  NA       1  1  1      1       1
#> 33     0    1    0    0    0    0      0   1  NA       1 NA NA      0       0
#> 34     1    1    1    1    1   NA      0   1  NA       1  0  0      0       1
#> 35     1    1    1    1    1    1      1   0   0       0  1  1      1       1
#> 36     0   NA   NA    0    0    1      1   0   0       0  0  0      0       0
#> 37     1   NA    0    0    0   NA      1   1  NA       1 NA  0      0       0
#> 38     1   NA    0    1    0    0      0   1   1       1  0  0      0       1
#> 39     0   NA    0   NA   NA   NA      0   1   0       1  0  0      0       1
#> 40     0    0    0    0    0    1      1   1   1       1  0  0      0       1
#> 41     0    0    0    0    0   NA     NA  NA   0       0  0 NA      0       0
#> 42     0   NA    0    0    0   NA     NA   1   1       1 NA  1      1       1
#> 43     1    1    1   NA   NA   NA      1   1   1       1  1 NA      1       1
#> 44     0   NA    0    0    0    1      1   1   1       1  0  0      0       0
#> 45     0   NA    0    0    0   NA      1  NA   1       1  0  0      0       1
#> 46     0   NA    0    0    0    0      0   1   1       1  0  0      0       0
#> 47     0   NA   NA    0    0   NA      1   1   1       1  0  1      1       1
#> 48     0    0    0   NA   NA    1      1   0   0       0  1  1      1       0
#> 49    NA    1    1   NA   NA    1      1   1   1       1  1  1      1       1
#> 50     0    0    0    0    0    0      1  NA   1       1 NA  0      0       1
#> 51    NA    0   NA    0   NA    0      0   0  NA       0  0 NA      0       0
#> 52     0   NA    0    0    0    0      0  NA   0       0  0 NA      0       1
#> 53     1    1    1    1    1   NA     NA   0   0       0  0  0      0       0
#> 54     0    0    0   NA   NA    1      1   0   0       0  1  1      1       1
#> 55     0    0    0   NA   NA    0      1   1   1       1  0  0      0       1
#> 56    NA    1    1    1    1    0      1   0   0       0  1  0      0       0
#> 57     1   NA    1    1    1    1      1   1   1       1  1  1      1       1
#> 58     1   NA    1   NA   NA    0      1   1   1       1  0  0      0       1
#> 59    NA    1    1   NA   NA    1      1   1   1       1  0 NA      0       1
#> 60     0    0    0    0    0    1      1   1   1       1 NA  0      0       0
#> 61     0    0    0    0    0    1      1  NA   1       1  0 NA      0       0
#> 62     1    1    1   NA   NA    1      1   0  NA       0  0  0      0       1
#> 63     0    0    0    0    0    1      1   1   1       1 NA NA      0       1
#> 64    NA    1    1   NA   NA    1      1   1   1       1  1  1      1       0
#> 65     0    0    0    0    0   NA     NA   1   1       1  1  0      0       1
#> 66     1    0    0   NA   NA    0      0   1  NA       1 NA  0      0       0
#> 67    NA    1    1    1    1    0      0   0   0       0  0  0      0       1
#> 68     1   NA    0    0    0   NA     NA   0   0       0  0  0      0       0
#> 69    NA    0    0    0    0    0      1  NA   0       0  0  0      0       0
#> 70     0   NA    0    0    0    0      1   1   1       1  1  1      1       0
#> 71     1   NA    1   NA   NA    1      1   0  NA       0  0 NA      0       0
#> 72     0    0    0    0    0    1      1   0   0       0 NA  0      0       1
#> 73     0   NA    1    1    1    1      1   0  NA       0  1 NA      1       0
#> 74     1    1    1    1    1    1      1   1   1       1  1  1      1       0
#> 75     1    1    1    1    1    1      1  NA  NA      NA  1  1      1       1
#> 76    NA    1    1    1    1    0      0  NA   0       0  0  0      0       0
#> 77     0    0   NA    0    0    1      1   0  NA       0  0  0      0       1
#> 78     0    0    0    0    0    0      0  NA   0       0 NA NA      1       0
#> 79    NA    1    1    1    1    0      0   1   1       1  0  0      0       0
#> 80     0    0    0    0    0   NA     NA   0  NA       0  1 NA      1       0
#> 81    NA    1    1    1    1   NA     NA  NA   1       1  1 NA      1      NA
#> 82    NA    0    0   NA   NA    1      1  NA   1       1  0  0      0       0
#> 83    NA    0   NA   NA   NA   NA      0   0  NA       0  1  1      1       0
#> 84    NA    1    1   NA   NA   NA     NA   1   1       1  1  1      1       0
#> 85     1    1    1    1    1    1      1   0   0       0 NA  1      1       0
#> 86     1    1   NA    0    1    1      1   1   1       1 NA  0      0       0
#> 87     0    0    0    0    0    1      1   1   1       1  1  1      1       1
#> 88     0    0   NA    0    0    1      1   1   1       1 NA  0      0       1
#> 89     1    1    1    1    1    1      1   0  NA       0  1  1      1       0
#> 90    NA    1    1    1    1    1      1   1   1       1  1 NA      1       0
#> 91     0    1    0    0    0   NA     NA   0   0       0  1  0      0       0
#> 92    NA    1    1    1    1    1      1   1   1       1  1 NA      1       1
#> 93     1    1    1    1    1    1      1   1   0       1  1  1      1       1
#> 94    NA    1    1    1    1    1      1   1   1       1  1  1      1       1
#> 95     1   NA    1    1    1   NA      1   0  NA       0 NA NA      0       1
#> 96     1   NA    1   NA   NA    1      1   1  NA       1  1  1      1       1
#> 97     0    0    0    0    0    1      1   1  NA       1  1  1      1       1
#> 98     0    1    1    1    1    0      0   0   0       0  0  0      0       1
#> 99     1    1    1    1    1    0      0   0   0       0  1  1      1       1
#> 100    0    0    0   NA   NA    1      1   1   1       1  0  0      0       0
#> 101    0   NA    0    0    0    1      1   0   0       0  0  0      0       1
#> 102    0    0    0    0    0    1      1  NA   1       1  1  1      1       1
#> 103    0    0    0   NA   NA    1      1  NA  NA      NA  1  1      1       0
#> 104    1   NA    1    1    1    0      0  NA   1       1  0  0      0       0
#> 105    0   NA    0   NA   NA   NA      1  NA   1       1 NA NA     NA       1
#> 106   NA    1    1   NA   NA   NA      1   0   0       0 NA NA     NA       0
#>     VSEGFR AFRAEGFR SRXRA SMT QMTAMPN QMTDELN SHL SEA SEB QPCRFGF7
#> 1        1        1     1   0       0       0   0   1   1        0
#> 2        1        1     1   1      NA       1   1   0   1        1
#> 3        0        1     0   0       0       0   0   0   1        1
#> 4        1        1     0   0       0       0   0   0   0        0
#> 5        0        1     1   0       0       0   0   0   0        1
#> 6        1        1     1   1      NA       1   1   0   0        1
#> 7        1        1     1   1       1      NA   1   0   1        1
#> 8        1        0     1   1       1      NA   1   1   1        0
#> 9        0        1     1   1      NA       1   0   1   1        0
#> 10       1        1     1   0       0       0   0   1   0        1
#> 11       0        0     1   0       0       0   0   0   1        0
#> 12       0        1     1   0       0       0   0   0   1        0
#> 13       1        1     0   1       0       0   1   1   1        0
#> 14       1        1     1   1       0       0   1   1   1        0
#> 15       0        1     1   0       0       0   0   0   1        0
#> 16       1        0     1   1       1      NA   1   1   1        1
#> 17       0        1     1   0       0       0   0   0   0        1
#> 18       0        1     1   1       0       0   0   1   0        0
#> 19       0        1     0   0       0       0   0  NA   0        0
#> 20       1        1     0   1       1      NA   1   1   1        1
#> 21       0        1     1   0       0       0   0   0   0        0
#> 22       0        1     1   0       1      NA   0   0   0        1
#> 23       0        1     0   0       0       0   0   0   0        1
#> 24       1        1     0   1       1      NA   1   0   0        1
#> 25       1        1     1   1       0       0   1   1   1        0
#> 26      NA        0     1   0       0       0   0   1   1        0
#> 27       0        1     1   0       0       0   0   0   1        0
#> 28       1        1     0   1       0       0   1   0   1        1
#> 29       1        1     1   1      NA       1   1   1   1        1
#> 30       1        1     1   1       1      NA   1   1   1        1
#> 31       1        1    NA   1       0       0   0   1   1        1
#> 32       0        1     1   0       0       0   0   1   1        1
#> 33       1        1     1   1       0       0   1   1   1        0
#> 34       1        1     1   1       1      NA   1   1   1        0
#> 35       1        0     1   1       1      NA   1   0   0        1
#> 36       0        1     0   0       0       0   0   1   1        1
#> 37       1        1     1   0       0       0   1   0   0        0
#> 38       0        0    NA   1       0       0   1   0   0        1
#> 39       0        1     1   1       0       0   0   1   0        1
#> 40       0        1     1   1       0       0   0   0   0        0
#> 41       0        1     1   0       0       0   0   0   0        0
#> 42       0        1     1   0       0       0   1   1   0        1
#> 43       1        1     1   1      NA       1   1   1   1        0
#> 44       0        1     1   0       0       0   0   0   0        1
#> 45       0        0     1   0       0       0   0   1   1        0
#> 46       0        0     0   1       0       0   1   0   0        0
#> 47       0        1     1   0       0       0   0   1   0        0
#> 48       0        1     1   0       1      NA   0   1   1        1
#> 49       1        1     1   1       1      NA   0   1   1        1
#> 50       0        1     0   1       0       0   0   1   1        1
#> 51       0        1     1   0       0       0   0   1   1        0
#> 52       0        0     1   0       0       0   0   1   1        1
#> 53       1       NA    NA   1       1      NA   1   0   1        1
#> 54       0        1     1   1       0       0   1   1  NA        1
#> 55       0       NA    NA   1       1      NA   0   1   1        0
#> 56       1        1     1   1      NA       1   1   0   0        1
#> 57       1        0     1   1       1      NA   0   1   1        2
#> 58       1        1     1   0       0       0   0   1   0        1
#> 59       1        1     1   0       0       0   0   0   1        1
#> 60       0        1     0   0       0       0   0   0   1        1
#> 61       1        1     1   0       0       0   0   1   0        1
#> 62       1        1     1   1       1      NA   1   1   1        0
#> 63       0        1     1   1       0       0   1  NA   1        1
#> 64       1        1     1   0       0       0   0   0   1        0
#> 65       0        1     1   0       0       0   0   1   0        1
#> 66       0        1     0   1      NA       1   1   0   0        0
#> 67       1        1     0   1      NA       1   1   0   1        0
#> 68       0        1     0   0       0       0   0   1   0        1
#> 69       0        1     0   0       0       0   0   0   0        1
#> 70       0        1     1   1       1      NA   1   1   0        1
#> 71       1        1     0   1       1      NA   1   1   1        1
#> 72       0        0     1  NA       0       0   0   0   0        0
#> 73       1        0     0   1       0       0   0   0   0        1
#> 74       1        1     1   1       0       0  NA   1   1        1
#> 75       1        1     0   1       0       0   1   0   0        1
#> 76       1        1     1  NA       0       0   0   0   0        1
#> 77       0        1     1   1       1      NA   0   1   1       NA
#> 78       0        1     0   0       0       0   0   1   0       NA
#> 79       1        1     0   0       0       0   1   0   1        1
#> 80       0        1     0   0       0       0   0   1   0        1
#> 81       1        1     1   0       0       0   1   1   1        0
#> 82       0        1     1   1       0       0   1   0   0        0
#> 83       0        0     0   0       0       0   0   1   0        1
#> 84       1        1    NA   0       0       0  NA   1   1        0
#> 85       1        1     0   1       0       0   1   1   0        1
#> 86       1        1     1   1       0       0   1  NA   0        1
#> 87       0       NA    NA   1       1      NA   0   0   1        0
#> 88       0        1     1   0       0       0   0   1   0        1
#> 89       1       NA    NA   1       0       0   1   1   1        0
#> 90       1        1     0  NA       0       0   1   1   0        0
#> 91       1        1     0   0       0       0   1   0   0        0
#> 92       1        0     1   0       0       0   0   1   1        1
#> 93       1        1     0   1       1      NA   1   1   1        0
#> 94       1        1     1   0       0       0   0   1   0        0
#> 95       1        1     1   1       0       0   1   1   1        0
#> 96       1        1     1   1      NA       1   0   0   1        1
#> 97       0        1     1   1       1      NA   1   1   0        1
#> 98       1        1     1   0       0       0   1   0   0       NA
#> 99       1        1     1   1       0       0   0   0   0        0
#> 100      0        1     1   0       0       0   1   0   1        1
#> 101      0        1     1   0       0       0   0   0   1        1
#> 102      0        0     0   1       0       0   0   1   1        1
#> 103      0        0     1   0       0       0   0   1   1        1
#> 104      1        1     1   1       0       0   1   1   1        1
#> 105      0        1     1   0       0       0   0   1   1        1
#> 106      1        1     0   1       1      NA   1   1   1       NA
summary(TxTum)
#>     CELTUMCO          age             sexe           HISTOADK     
#>  Min.   :20.00   Min.   :41.00   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:30.00   1st Qu.:52.00   1st Qu.:1.0000   1st Qu.:1.0000  
#>  Median :50.00   Median :62.00   Median :1.0000   Median :1.0000  
#>  Mean   :50.17   Mean   :61.26   Mean   :0.7925   Mean   :0.7609  
#>  3rd Qu.:68.75   3rd Qu.:69.00   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :99.00   Max.   :80.00   Max.   :1.0000   Max.   :1.0000  
#>                                                   NA's   :14      
#>        H2               P3              P4              E1      
#>  Min.   :0.0000   Min.   :0.000   Min.   :0.000   Min.   :0.00  
#>  1st Qu.:1.0000   1st Qu.:0.500   1st Qu.:0.000   1st Qu.:0.75  
#>  Median :1.0000   Median :1.000   Median :1.000   Median :1.00  
#>  Mean   :0.7568   Mean   :0.747   Mean   :0.725   Mean   :0.75  
#>  3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:1.000   3rd Qu.:1.00  
#>  Max.   :1.0000   Max.   :1.000   Max.   :1.000   Max.   :1.00  
#>  NA's   :32       NA's   :23      NA's   :26      NA's   :22    
#>        P5              R10              C3M               P6        
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:1.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :0.0000   Median :1.0000   Median :1.0000   Median :1.0000  
#>  Mean   :0.4557   Mean   :0.7625   Mean   :0.7246   Mean   :0.6593  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :27       NA's   :26       NA's   :37       NA's   :15      
#>        RB              FL7A             P53              W2        
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.000   Min.   :0.0000  
#>  1st Qu.:1.0000   1st Qu.:0.0000   1st Qu.:1.000   1st Qu.:0.0000  
#>  Median :1.0000   Median :1.0000   Median :1.000   Median :1.0000  
#>  Mean   :0.7614   Mean   :0.6364   Mean   :0.837   Mean   :0.5333  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.000   Max.   :1.0000  
#>  NA's   :18       NA's   :29       NA's   :14      NA's   :31      
#>        P2               P1              W4              MT1        
#>  Min.   :0.0000   Min.   :0.000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :0.0000   Median :0.000   Median :0.0000   Median :1.0000  
#>  Mean   :0.4945   Mean   :0.407   Mean   :0.3333   Mean   :0.5231  
#>  3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :15       NA's   :20      NA's   :34       NA's   :41      
#>       MT2           MT4            MT3            HLA              HLD       
#>  Min.   :0.0   Min.   :0.00   Min.   :0.00   Min.   :0.0000   Min.   :0.000  
#>  1st Qu.:0.0   1st Qu.:0.00   1st Qu.:0.00   1st Qu.:0.0000   1st Qu.:0.000  
#>  Median :0.5   Median :0.00   Median :0.00   Median :0.0000   Median :1.000  
#>  Mean   :0.5   Mean   :0.48   Mean   :0.48   Mean   :0.4559   Mean   :0.507  
#>  3rd Qu.:1.0   3rd Qu.:1.00   3rd Qu.:1.00   3rd Qu.:1.0000   3rd Qu.:1.000  
#>  Max.   :1.0   Max.   :1.00   Max.   :1.00   Max.   :1.0000   Max.   :1.000  
#>  NA's   :42    NA's   :31     NA's   :56     NA's   :38       NA's   :35     
#>       HLC              HLB              EA1              EA3        
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :0.0000   Median :0.0000   Median :1.0000   Median :0.0000  
#>  Mean   :0.4242   Mean   :0.4744   Mean   :0.6102   Mean   :0.4857  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :73       NA's   :28       NA's   :47       NA's   :36      
#>       EA2              EA4             EB1              EB2        
#>  Min.   :0.0000   Min.   :0.000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :1.0000   Median :1.000   Median :0.0000   Median :1.0000  
#>  Mean   :0.5753   Mean   :0.597   Mean   :0.4286   Mean   :0.5211  
#>  3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :33       NA's   :39      NA's   :36       NA's   :35      
#>       EB3              EB4              EGF1             EGF2       
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :1.0000   Median :1.0000   Median :1.0000   Median :0.0000  
#>  Mean   :0.5405   Mean   :0.5429   Mean   :0.5333   Mean   :0.4444  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :32       NA's   :36       NA's   :46       NA's   :25      
#>       EGF3             EGF4             EGF5           EGF6       
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.00   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.00   1st Qu.:0.0000  
#>  Median :1.0000   Median :0.0000   Median :0.00   Median :0.0000  
#>  Mean   :0.5068   Mean   :0.4831   Mean   :0.44   Mean   :0.4384  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.00   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.00   Max.   :1.0000  
#>  NA's   :33       NA's   :17       NA's   :31     NA's   :33      
#>       FL7B            VSFGF7            F3A              F3B        
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :1.0000   Median :1.0000   Median :1.0000   Median :1.0000  
#>  Mean   :0.6133   Mean   :0.7158   Mean   :0.5926   Mean   :0.6098  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :31       NA's   :11       NA's   :25       NA's   :24      
#>     VSFGFR3           F4              Q5             VSTOP1      
#>  Min.   :0.00   Min.   :0.000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.00   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :1.00   Median :0.000   Median :0.0000   Median :0.0000  
#>  Mean   :0.61   Mean   :0.439   Mean   :0.3718   Mean   :0.3883  
#>  3rd Qu.:1.00   3rd Qu.:1.000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.00   Max.   :1.000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :6      NA's   :24      NA's   :28       NA's   :3       
#>     VSTOP2A           VSEGFR          AFRAEGFR          SRXRA       
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:1.0000   1st Qu.:0.0000  
#>  Median :1.0000   Median :1.0000   Median :1.0000   Median :1.0000  
#>  Mean   :0.5196   Mean   :0.5143   Mean   :0.8431   Mean   :0.7071  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :4        NA's   :1        NA's   :4        NA's   :7       
#>       SMT            QMTAMPN          QMTDELN            SHL        
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :1.0000   Median :0.0000   Median :0.0000   Median :0.0000  
#>  Mean   :0.5243   Mean   :0.2268   Mean   :0.1071   Mean   :0.4423  
#>  3rd Qu.:1.0000   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
#>  NA's   :3        NA's   :9        NA's   :22       NA's   :2       
#>       SEA              SEB            QPCRFGF7    
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.000  
#>  Median :1.0000   Median :1.0000   Median :1.000  
#>  Mean   :0.5728   Mean   :0.5619   Mean   :0.598  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.000  
#>  Max.   :1.0000   Max.   :1.0000   Max.   :2.000  
#>  NA's   :3        NA's   :1        NA's   :4      
```
