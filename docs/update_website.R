# Update website

rmarkdown::render(here::here("docs", "index.Rmd"))
rmarkdown::render(here::here("docs", "project_background.Rmd"))
rmarkdown::render(here::here("docs", "covariates_overview.Rmd"))
rmarkdown::render(here::here("docs", "recruitment_overview.Rmd"))
rmarkdown::render(here::here("docs", "ml_predictions.Rmd"))
rmarkdown::render(here::here("docs", "multivariate_analyses.Rmd"))
rmarkdown::render(here::here("docs", "time_series_analyses.Rmd"))
rmarkdown::render(here::here("docs", "ml_predictions.Rmd")) # note that if significant updates are made to this script, it will take a long time to run (otherwise, progress is cached)
rmarkdown::render(here::here("docs", "groundfish_distributions.Rmd"))
rmarkdown::render(here::here("docs", "trawl_survey_indices.Rmd"))