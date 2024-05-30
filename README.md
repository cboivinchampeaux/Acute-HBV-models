# Acute-HBV-models
Shiny app for acute HBV models

Make sure that you R and Rstudio versions are up to date.
Make sure to install the different packages that are in the Shinyapp - Acute HBV models.Rmd file if not already installed: <code>install.packages("Packagename")</code>

To see the mean profiles for the different models: 
1. Download the ZIP.
2. In the server.R file, don't forget to set the correct directory for the observations file in:
   <code>obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_meanprofiles/Dataset/HBVDNA_observations.csv")</code>
3. In the Shinyapp - Acute HBV models.Rmd file run shiny::runApp('Acute_HBV_meanprofiles')

To see the 6 patients profiles for the different models:  
1. Download the ZIP.
2. In the server.R file, don't forget to set the correct directory for all the patients' observations files, such as in:
   <code>patient1_obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_individualprofiles/Dataset/patient1-observationV.csv")</code>
3. In the Shinyapp - Acute HBV models.Rmd file run shiny::runApp('Acute_HBV_individualprofiles')

