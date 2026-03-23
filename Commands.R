unlink("NAMESPACE")  # delete NAMESPACE as roxygen2 will generate new one


usethis::use_r("IPWM")
devtools::document()




devtools::load_all(reset = TRUE)






devtools::build()

devtools::check()
devtools::install(".", force = TRUE)

library(IPWMethod)
?IPWM_build





build()      # Builds .tar.gz
install()    # Installs locally
load_all()   # Loads package into current R session


# build example dataset
usethis::use_data_raw() # detailed commands in the data_raw file


data("sc")
data("ref_survey_1")
data("ref_survey_2")

?sc
?ref_survey_1
?ref_survey_2

str(sc)
str(ref_survey_1)
str(ref_survey_2)


# for others to use
install.packages("devtools")
devtools::install_github("Jiakun0611/IPWMethod")

# check

remove.packages("IPWMethod") # better to restart R session

browseURL(system.file("Doc", "formula_cheatsheet.pdf", package = "IPWMethod"))




path <- "."
hits <- list.files(path, recursive = TRUE, pattern = "\\.R$", full.names = TRUE)
out <- lapply(hits, function(f) {
  x <- readLines(f, warn = FALSE)
  w <- grep("design_matrix\\s*\\(", x)
  if (length(w)) data.frame(file=f, line=w, code=x[w], stringsAsFactors=FALSE)
})
do.call(rbind, out)

