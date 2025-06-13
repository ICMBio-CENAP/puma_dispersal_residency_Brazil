# to render documents

# 00 organize raw GPS data
# only raw script

# 01 load and prepare GPS data for analysis 
rmarkdown::render("code/01_dispersal_load_data.R", output_format = "all", output_dir = "notebooks/")

# 02 create data package
rmarkdown::render("code/02_create_data_package.R", output_format = "all", output_dir = "notebooks/")

# 03 explore gps data to analysis
rmarkdown::render("code/03_plot_explore.R", output_format = "all", output_dir = "notebooks/")

# 04 dispersal fit for 1 individual
rmarkdown::render("code/04_arima_fit_1ind.R", output_format = "all", output_dir = "notebooks/")

# 05 dispersal fit for all individuals
rmarkdown::render("code/05_arima_fit_multiple_animals_pkg.R", output_format = "all",
                  output_dir = "notebooks/")

# 06 understand dispersal
rmarkdown::render("code/06_arima_understand_dispersal.R", output_format = "all",
                  output_dir = "notebooks/")

# 07 understand movement
rmarkdown::render("code/07_understand_movement.R", output_format = "all",
                  output_dir = "notebooks/")

# 08 understand residency
rmarkdown::render("code/08_understand_residency.R", output_format = "all",
                  output_dir = "notebooks/")

# 09 literature
rmarkdown::render("code/09_compare_literature.R", output_format = "all",
                  output_dir = "notebooks/")
