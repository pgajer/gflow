## Generating a blue-to-yellow-to-red color palette
## Defining the colors for the lowest, mid, and upper range values
colors <- c("blue", "yellow", "red")
## Creating a colorRampPalette with the specified colors
ppalette <- colorRampPalette(colors)
## Specifying the number of colors you want in the palette
n.colors <- 100
blue.yellow.red.color.palette <- ppalette(n.colors)

rainbow.palette <- rainbow(12,start=1/6,end=0)

Linf.cst.col.tbl <- c(
  Lactobacillus_crispatus = "red1",
  Lactobacillus_gasseri = "chartreuse",
  Lactobacillus_iners = "darkorange2",
  BVAB1 = "aquamarine4",
  Atopobium_vaginae = "orange",
  Gardnerella_vaginalis = "royalblue",
  Sneathia_sanguinegens = "limegreen",
  g_Anaerococcus = "blue",
  g_Corynebacterium_1 = "gold",
  g_Streptococcus = "brown",
  g_Enterococcus = "deeppink",
  g_Bifidobacterium = "darkorchid",
  Lactobacillus_jensenii = "yellow"
)

cst13.col.tbl <- c()
cst13.col.tbl["I-A"] <- "red1"
cst13.col.tbl["I-B"] <- "palevioletred2"
cst13.col.tbl["II"] <- "chartreuse"
cst13.col.tbl["III-A"] <- "darkorange2" # aquamarine4"
cst13.col.tbl["III-B"] <- "orange1"
cst13.col.tbl["IV-A"] <- "aquamarine4"
cst13.col.tbl["IV-B"] <- "royalblue"
cst13.col.tbl["IV-C0"] <- "blue"
cst13.col.tbl["IV-C1"] <- "brown"
cst13.col.tbl["IV-C2"] <- "deeppink"
cst13.col.tbl["IV-C3"] <- "darkorchid"
cst13.col.tbl["IV-C4"] <- "cyan"
cst13.col.tbl["IV-C"] <- "palevioletred4" # "deeppink"
cst13.col.tbl["V"] <- "yellow"

CST13s <- sort(unique(names(cst13.col.tbl)))

##
cst7AB.col.tbl <- c()
cst7AB.col.tbl["I-A"] <- "red1"
cst7AB.col.tbl["I-B"] <- "palevioletred2"
cst7AB.col.tbl["II"] <- "chartreuse"
cst7AB.col.tbl["III-A"] <- "darkorange2" # aquamarine4"
cst7AB.col.tbl["III-B"] <- "orange1"
cst7AB.col.tbl["IV-A"] <- "aquamarine4"
cst7AB.col.tbl["IV-B"] <- "royalblue"
cst7AB.col.tbl["IV-C"] <- "palevioletred4"
cst7AB.col.tbl["V"] <- "yellow"

## different colors
cst7AB.d.col.tbl <- c()
cst7AB.d.col.tbl["I-A"] <- "red1"
cst7AB.d.col.tbl["I-B"] <- "palevioletred2"
cst7AB.d.col.tbl["II"] <- "chartreuse"
cst7AB.d.col.tbl["III-A"] <- "darkorange2" # aquamarine4"
cst7AB.d.col.tbl["III-B"] <- "orange1"
cst7AB.d.col.tbl["IV-A"] <- "blue"
cst7AB.d.col.tbl["IV-B"] <- "#00CCFF"
cst7AB.d.col.tbl["V"] <- "yellow"

##
cst7.col.tbl <- c()
cst7.col.tbl["I"] <- "red1"
cst7.col.tbl["II"] <- "chartreuse"
cst7.col.tbl["III"] <- "darkorange2"
cst7.col.tbl["IV-A"] <- "aquamarine4"
cst7.col.tbl["IV-B"] <- "royalblue"
cst7.col.tbl["IV-C"] <- "palevioletred4" # "deeppink"
cst7.col.tbl["V"] <- "yellow"

CST7s <- sort(unique(names(cst7.col.tbl)))

##
cst5.col.tbl <- c()
cst5.col.tbl["I"] <- "red1"
cst5.col.tbl["II"] <- "chartreuse"
cst5.col.tbl["III"] <- "darkorange2"
cst5.col.tbl["IV"] <- "aquamarine4"
cst5.col.tbl["V"] <- "yellow"

CST5s <- sort(unique(names(cst5.col.tbl)))

bin.cst.col.tbl <- c("red1", "royalblue") ##c("chartreuse", "royalblue")
names(bin.cst.col.tbl) <- c("Lb", "IV")

cst2.col.tbl <- c("red", "aquamarine4")
names(cst2.col.tbl) <- c("Lb", "IV")

cst3.col.tbl <- c("red", "darkorange2", "aquamarine4")
names(cst3.col.tbl) <- c("I|II|V", "III", "IV")

bin1.col.tbl <- c()
bin1.col.tbl['1'] <- "green"
bin1.col.tbl['2'] <- "red"

bin.col.tbl <- c()
bin.col.tbl['0'] <- "gray80"
bin.col.tbl['1'] <- "red"

tri.col.tbl <- c()
tri.col.tbl['0'] <- "gray80"
tri.col.tbl['1'] <- "red"
tri.col.tbl['2'] <- "cyan"
