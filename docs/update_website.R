# Update website

rmarkdown::render(here::here("docs", "index.Rmd"))
rmarkdown::render(here::here("docs", "project_background.Rmd"))
rmarkdown::render(here::here("docs", "covariates_overview.Rmd"))
rmarkdown::render(here::here("docs", "recruitment_overview.Rmd"))
rmarkdown::render(here::here("docs", "ml_predictions.Rmd"))
rmarkdown::render(here::here("docs", "multivariate_analyses.Rmd"))