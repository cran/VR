datanames <-
c("abbey", "accdeaths", "Aids2", "Animals", "anorexia", "austres",
"bacteria", "beav1", "beav2", "biopsy", "birthwt", "Boston",
"cabbages", "caith", "Cars93", "cats", "cement", "chem", "coop",
"cpus", "crabs", "Cushings", "DDT", "deaths", "drivers", "eagles",
"epil", "farms", "fdeaths", "fgl", "forbes", "GAGurine", "galaxies",
"gehan", "genotype", "geyser", "gilgais", "hills", "housing",
"immer", "Insurance", "leuk", "lh", "mammals", "mcycle", "mdeaths",
"Melanoma", "menarche", "michelson", "minn38", "motors", "muscle",
"newcomb", "nlschools", "nottem", "npk", "npr1", "oats", "OME",
"painters", "petrol", "phones", "Pima.te", "Pima.tr", "Pima.tr2",
"quine", "Rabbit", "road", "rock", "rotifer", "Rubber", "ships",
"shoes", "shrimp", "shuttle", "Sitka", "Sitka89", "Skye", "snails",
"SP500", "steam", "stormer", "survey", "synth.te", "synth.tr",
"topo", "Traffic", "UScereal", "UScrime", "VA", "waders", "whiteside",
"wtloss")

for(i in datanames)
    eval(substitute(obj <- delay(MASS.data.load(obj)), list(obj=i)))
rm(datanames,i)

MASS.data.load <- function(i)
{
    file <- file.path(system.file("data", package="MASS"),
                      paste(i, ".rda", sep=""))
    zfile <- zip.file.extract(file, "Rdata.zip")
    ## this *has* to be a temp environment, not the package.
    env <- new.env()
    load(zfile, envir=env)
    get(i, env)
}
