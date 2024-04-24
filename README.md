# Acute-HBV-models
Shiny app for acute HBV models

To see the mean profiles for the different models: put both the ui.R and server.R file into a folder named "Acute_HBV_meanprofile". 
In that folder, create a folder named "Dataset" and download the combined datasets for ALT and HBV DNA in that folder.
In the server.R file, don't forget to set the directory where you have your observation files saved.
Set the Working Direction to the location of the folder (but not the folder itself) and type in: shiny::runApp('Acute_HBV_meanprofile')

To see the 7 patients profiles for the different models: put both the ui.R and server.R file into a folder named "Acute_HBV_patients". 
In that folder, create a folder named "Dataset" and download the 7 patients datasets for HBV DNA in that folder.
In the server.R file, don't forget to set the directory where you have your observation files saved.
Set the Working Direction to the location of the folder (but not the folder itself) and type in: shiny::runApp('Acute_HBV_patients')
