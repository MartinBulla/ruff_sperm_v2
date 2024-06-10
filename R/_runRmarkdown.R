 # =============================================================
 # ❗ Runs relative to the project's root directory and
 # takes R and Rmd scripts and renders them into html exported
 # within ./Outputs/
 # =============================================================

require(rmarkdown)

rmarkdown::render("R/Out_velocity-videos_anonym.Rmd", output_dir = "Outputs", output_file = "Velocity_video_examples_anonym.html") # rmarkdown::render('R/Out_velocity_videos_examples.R', output_dir = 'Outputs',  output_file = 'Velocity_video_examples_test.html')
rmarkdown::render("R/Out_velocity-videos.Rmd", output_dir = "Outputs", output_file = "Velocity_video_examples.html") # rmarkdown::render('R/Out_velocity_videos_examples.R', output_dir = 'Outputs',  output_file = 'Velocity_video_examples_test.html')

rmarkdown::render('Protocols/Preregistration_v3.Rmd', output_dir = 'Protocols')

rmarkdown::render('Protocols/protocol_sperm.Rmd', output_dir = 'Protocols')
#rmarkdown::render('R/protocol_sperm.Rmd', output_dir = 'Protocols', word_document(toc = TRUE))

# END