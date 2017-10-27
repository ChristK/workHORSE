# Preample ---------------------------------------------------
if (Sys.info()[1] == "Linux") {
  if (system("whoami", T )== "mdxasck2") {
    setwd("~/IMPACTncd Liverpool/")
    clusternumber <- ifelse(clusternumber > 20, 20, clusternumber)  
  } else {
    setwd(paste("/home/", 
                system("whoami", T), 
                "/Dropbox/PhD/Models/IMPACTncd Liverpool/", 
                sep = "", 
                collapse = ""))
  }
} else if (Sys.info()[1] == "Darwin") {
  setwd("~/dropbox/PhD/Models/IMPACTncd Liverpool/")
  sav.dir <- "~/dropbox/PhD/Models/IMPACTncd Liverpool/Output/Graphs"
} else {
  get.dropbox.folder <- function() {
    if (!require(RCurl)) 
      stop("You need to install RCurl package.")
    if (Sys.info()["sysname"] != "Windows") 
      stop("Currently, 'get.dropbox.folder' works for Windows and Linux only. Sorry.")
    db.file <- paste(Sys.getenv("LOCALAPPDATA"), "\\Dropbox\\host.db", sep = "")
    base64coded <- readLines(db.file, warn = F)[2]
    base64(base64coded, encode = F)
  }
  setwd(paste0(get.dropbox.folder(), "/PhD/Models/IMPACTncd Liverpool/"))
  sav.dir <- paste0(get.dropbox.folder(), "/PhD/Models/IMPACTncd Liverpool/Output/Graphs")
}

require(Cairo)
if (Sys.info()[1] == "Windows") Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.20/bin/gswin64c.exe")
require(extrafont)
# font_import() # to be run only once in each system
#loadfonts(device = "win", quiet = T)
loadfonts()
#fonttable()
require(data.table)
require(ggplot2)
require(ggrepel)
require(scales)
require(ggthemes)
require(Rcpp)
require(RColorBrewer)

# output quantiles
output_quants <- 
  function(x, nam = get_object_name(x), pop_extrapol = pop.fraction, prob = c(0.025, .2, 0.25, 0.5, 0.75, 0.8, 0.975)) {
    out <- as.list(round(quantile(x, prob,
                                  na.rm = T) / pop_extrapol))
    names(out) <- paste0(nam, "_", names(out))
    return(out)
  }

output_quants_lbl <- 
  function(x, nam = get_object_name(x), pop_extrapol = pop.fraction, prob = c(0.5, 0.025, 0.975)) {
    out <- as.list(round(quantile(x, prob,
                                  na.rm = T) / pop_extrapol))
    out <- list(paste0(prettyNum(out[1], ","), " (", prettyNum(out[2], ","), " to ", prettyNum(out[3], ","), ")"))
    names(out) <- paste0(nam)
    return(out)
  }


qimd_labeller <- labeller(qimd = c("1" = "QIMD 1\nleast deprived",
                                   "2" = "QIMD 2",
                                   "3" = "QIMD 3",
                                   "4" = "QIMD 4",
                                   "5" = "QIMD 5\nmost deprived"))
gbp <- dollar_format(prefix = "£", suffix = "")

# normalise a vector to 0,1 range
sourceCpp("functions.cpp", cacheDir = "./CppCache/")
identical_elements <- 
  function(x, tol = .Machine$double.eps ^ 0.5) {
    stopifnot(is.numeric(x))
    fequal(x, tol)
  }

normalise <-
  function(x, ...){
    stopifnot(is.numeric(x))
    if (identical_elements(x)) return(1)
    else return(fnormalise(x))
  }

my_theme <-   theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = rel(1), vjust = 0),
        axis.title = element_text(size = rel(1)),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(size = rel(0.8), vjust = 0.5, angle = 90),
        axis.text.y = element_text(size = rel(0.8), vjust = 0),
        strip.text.x = element_text(rel(.8)),
        strip.text.y = element_text(rel(.8)),
        strip.background = element_rect(colour = "lavenderblush4", fill = "lavenderblush2"),
        legend.title = element_blank(), # no legent title
        legend.position = "right",
        #legend.justification = c(0,1), legend.position = c(0,1),
        legend.key = element_rect(size = 1, colour = "transparent", fill = "transparent"),
        legend.key.size = unit(1, "lines")) 

# Load data
results <- fread("./Output/results.csv")
results <- results[scenario == "sc0"]
if (!exists("pop.fraction"))
{
  tt <- readLines("./Output/simulation parameters.txt")
  pop.fraction <-
    as.numeric(substring(tt[[grep(glob2rx("Population fraction = *"), tt)]], 23))
}

results[scenario == "sc0", scenario := "User Scenario"]
results[, scenario := factor(scenario)]
results[, cumul.net.utility.disc.median := median(cumul.net.utility.disc, na.rm = T), by = .(year, scenario)]
results[, cumul.net.cost.disc.median := median(cumul.net.cost.disc, na.rm = T), by = .(year, scenario)]

study.year <- 2031 # 15 years post intervention
study.sc <- "User Scenario"
colourCount <- length(unique(levels(results$scenario)))
getPalette <- colorRampPalette(brewer.pal(9, "Set1")) # create longer colour palette by interpolation

ggplot(results[year == study.year & scenario %in% study.sc],
       aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2, show.legend = F) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 3, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% study.sc,
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme  + theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))
