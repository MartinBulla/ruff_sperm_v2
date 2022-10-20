require(rmarkdown)

rmarkdown::render('R/Out_velocity-videos.Rmd', output_dir = 'Outputs',  output_file = 'Velocity_video_examples.html')
rmarkdown::render('R/EXP_morpho_mot.R', output_dir = 'Output',  output_file = 'morpho_mot.html')
rmarkdown::render('R/EXP_morpho.R', output_dir = 'Output',  output_file = 'morpho.html')

rmarkdown::render('R/motility_html.R', output_dir = 'Output',  output_file = 'motility.html')
#rmarkdown::render('R/motility.R', output_dir = 'Output')
rmarkdown::render('R/EXP_test-40.Rmd', output_dir = 'Output')

rmarkdown::render('R/Methods.Rmd', output_dir = 'Output')

rmarkdown::render('R/Preregistration_v3.Rmd', output_dir = 'Protocols')

rmarkdown::render('R/Preregistration_v2.Rmd', output_dir = 'Protocols')

rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols')
#rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols', word_document(toc = TRUE))
#rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols', word_document(highlight = "zenburn"))

rmarkdown::render('R/protocol_sperm.Rmd', output_dir = 'Protocols')
#rmarkdown::render('R/protocol_sperm.Rmd', output_dir = 'Protocols', word_document(toc = TRUE))