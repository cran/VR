rts <- function(units = NULL, ...) ts(...)
# file MASS/data.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
cement <- read.table("cement.dat", header=T)
cpus <- read.table("cpus.dat", sep=",")
forbes <- read.table("forbes.dat", header=T)
hills <- read.table("hills.dat", sep="\t", header=T)
road <- read.table("road.dat", sep="\t", header=T)
Rubber <- read.table("rubber.dat", header=T)
trees <- read.table("trees.dat", header=T)
mammals <- read.table("mammals.dat",sep="\t", header=T)
ships <- read.table("ships.dat", header=T)
mcycle <- read.table("mcycle.dat", header=T)
rock <- read.table("rock.dat")
rotifer <- read.table("rotifer.dat", header=T)
leuk <- read.table("leuk.dat")
gehan <- read.table("gehan.dat")
motors <- read.table("motors.dat", header=T)
painters <- read.table("painters.dat", header=T)
shuttle <- read.table("shuttle.dat", sep=",", header=T)
wtloss <- read.table("wtloss.dat")
cats <- read.table("cats.dat")
quine <- read.table("quine.dat")
oats <- read.table("oats.dat")
immer <- read.table("immer.dat")
steam <- read.table("steam.dat")
stormer <- read.table("stormer.dat")
michelson1 <- read.table("michels.dat")
attach(michelson1)
michelson <- data.frame(Speed, Run=factor(Run), Expt=factor(Expt))
invisible(detach())
rm(michelson1)
faithful <- read.table("faithful.dat")
galaxies <- scan("galaxies.dat")
birthwt <- read.table("birthwt.dat")
minn38 <- read.table("minn38.dat")
phones <- scan("phones.dat",list(year=0, calls=0))
   phones$calls <- 10*phones$calls
Animals <- read.table("animals.dat",sep=",", header=T)
crabs <- read.table("crabs.dat", header=T)
survey <- read.table("survey.dat", sep=",", header=T, row.names=1)
coop <- read.table("coop.dat")
Melanoma <- read.table("melanoma.dat")
Aids2 <- read.table("aids2.dat", header=T)
Tlev <- c("hs", "hsid", "id", "het", "haem",  "blood", "mother", "other")
Aids2 <- data.frame(state=c("NSW", "Other", "QLD", "VIC")[Aids2$state],
   sex=Aids2$sex, diag=Aids2$diag, death=Aids2$death, status=Aids2$status,
   T.categ=factor(Tlev[Aids2$T.categ], levels=Tlev), age=Aids2$age)
rm(Tlev)


abbey <- scan("misc.dat", skip=32, nmax=31)
chem <- scan("misc.dat", skip=36, nmax=24)
newcomb <- scan("misc.dat", skip=39, nmax=66)
shrimp <- scan("misc.dat", skip=51, nmax=18)
DDT <- scan("misc.dat", skip=54, nmax=15)
shoes <- scan("misc.dat",list(A=0, B=0), skip=57, nmax=20)


topo <- read.table("topo.dat", header=T)
shkap <- read.table("shkap.dat", header=T)
npr1 <- read.table("npr.dat")
Rabbit <- read.table("Rabbit.dat")
petrol <- read.table("petrol.dat")
sitka <- scan("Sitka.dat")
Sitka <- data.frame(size = sitka,
           Time = rep(c(152,174,201,227,258), 79),
           tree = rep(1:79, rep(5, 79)),
           treat = rep(c("ozone", "control"), c(54, 25)*5))
rm(sitka)
sitka <- scan("Sitka89.dat")
Sitka89 <- data.frame(size = sitka,
           Time = rep(c(469, 496, 528, 556, 579, 613, 639, 674), 79),
           tree = rep(1:79, rep(8, 79)),
           treat = rep(c("ozone", "control"), c(54, 25)*8))
rm(sitka)
Cushings <- read.table("Cushings.dat", header=T, row.names=1)
fgl <- read.table("fglass.dat")
fgl$type <- factor(fgl$type, labels=c("WinF", "WinNF", "Veh",  "Con", "Tabl", "Head"))
fgl$RI <- 1000*fgl$RI - 1518

beav1 <- read.table("beaver1.dat", header=T)
beav2 <- read.table("beaver2.dat", header=T)


menarche <- read.table("menarche.dat", header=T)
UScereal <- read.table("UScereal.dat", quote="", header=T, sep=",")
biopsy <- read.table("biopsy.dat", sep=",", na.strings="?",header=T)
biopsy$class <- factor(biopsy$class, labels=c("benign", "malignant"))
biopsy$ID <- as.character(biopsy$ID)
gilgais <- read.table("gilgais.dat", header=T, row.names=1)
UScrime <- read.table("UScrime.dat")
Skye <- read.table("Skye.dat", header=T)
GAGurine <- read.table("GAGurine.dat")
snails <- read.table("snails.dat", header=T)
genotype <- read.table("genotype.dat", header=T)
Cars93 <- read.table("Cars93.dat", sep=",", header=T)
Traffic <- read.table("traffic.dat", header=T)
Insurance <- read.table("insur.dat", header=T)
Insurance$Group <- ordered(Insurance$Group,
                           labels=c("<1l", "1-1.5l", "1.5-2l", ">2l"))
Insurance$Age <- ordered(Insurance$Age,
                         labels=c("<25", "25-29", "30-35", ">35"))
Insurance$District <- factor(Insurance$District)
OME <- read.table("OME.dat", header=T)
Boston <- read.table("Boston.dat", header=T)
synth.te <- read.table("synth.te", header=T)
synth.tr <- read.table("synth.tr", header=T)
Pima.te <- read.table("pima.te", header=T)
Pima.tr <- read.table("pima.tr", header=T)
Pima.tr2 <- read.table("pima.tr2", header=T)
farms <- read.table("farms.dat", header=T)
waders <- read.table("waders.dat", header=T)
caith <- read.table("Fisher.dat")
npk1 <- read.table("npk.dat", header=T)
npk <- data.frame(block=factor(rep(1:6, rep(4,6))),
      N=factor(npk1$N), P=factor(npk1$P), K=factor(npk1$K),
      yield=npk1$yield)
rm(npk1)
cabbages <- read.table("cabbage.dat")
whiteside <- read.table("Whiteside.dat")
whiteside$Insul <- factor(whiteside$Insul, levels=c("Before", "After"))
anorexia <- read.table("anorexia.dat")
housing <- read.table("housing.dat")
housing$Sat <- ordered(housing$Sat, levels=c("Low", "Medium", "High"))
housing$Infl <- factor(housing$Infl, levels=c("Low", "Medium", "High"))
housing$Cont <- factor(housing$Cont, levels=c("Low", "High"))
housing$Type <- factor(housing$Type, levels=c("Tower", "Apartment",
                                       "Atrium", "Terrace"))
muscle <- read.table("muscle.dat")
eagles <- read.table("eagles.dat", header=T)

tmp <- scan("austres.dat", what=list(num=0,year=0, mon=""))
nottem <- rts(scan("nottem.dat"), start=1920, frequency=12, units="months")
lh <- rts(scan("misc.dat", nmax=48), units="10mins")
deaths <- rts(scan("misc.dat", skip=4, nmax=72), start=1974,
  frequency=12, units="months")
accdeaths <- rts(scan("misc.dat", skip=11, nmax=72), start=1973,
  frequency=12, units="months")
mdeaths <- rts(scan("misc.dat", skip=18, nmax=72), start=1974,
  frequency=12, units="months")
fdeaths <- rts(scan("misc.dat", skip=25, nmax=72), start=1974,
  frequency=12, units="months")
drivers <- rts(scan("drivers.dat"), start=1969, frequency=12, units="months")
austres <- rts(tmp$num, start=c(1971,2), freq=4, units="quarters")
rm(tmp)

rm(rts)
for (f in ls()) 
      save(list = f, file = paste(f, ".rda", sep = ""), ascii = TRUE)
#      dump(list = f, file = paste(f, ".R", sep = ""))
