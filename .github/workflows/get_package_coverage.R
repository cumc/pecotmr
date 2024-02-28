#!/usr/bin/env Rscript

library(covr)
library(jsonlite)
library(magrittr)
library(readr)

coverage_output <- package_coverage()
print(coverage_output)
coverage_pct <- percent_coverage(coverage_output)
coverage_list <- list(
    lines = list(pct = coverage_pct),
    statements = list(pct = coverage_pct),
    functions = list(pct = coverage_pct),
    branches = list(pct = coverage_pct)
)
coverage_list_final <- list(total = coverage_list)
dir.create("coverage")
toJSON(coverage_list_final) %>% write_lines("coverage/coverage-summary.json")
