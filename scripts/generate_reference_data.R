devtools::load_all()

test_cases = c(
  "Botswana2017",
  "Botswana2018",
  "DominicanRepublic2017",
  "Mozambique_Maputo_Cidade2018",
  "Netherlands2017"
)
for (test_case in test_cases) {
  print(test_case)
  pjnz <- system.file("extdata/testpjnz",
                      paste(test_case, ".PJNZ", sep = ""),
                      package = "eppasm")
  bw <- prepare_spec_fit(pjnz, proj.end = 2021.5)
  dir.create("tests/testthat/reference", FALSE, TRUE)
  saveRDS(bw, paste("tests/testthat/reference/", test_case, ".rds", sep = ""))
}
