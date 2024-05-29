# Acute-HBV-models
Shiny app for acute HBV models

To see the mean profiles for the different models: 
1. Download the ZIP.
2. In the server.R file, don't forget to set the correct directory of the observation file:
   obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_meanprofiles/Dataset/HBVDNA_observations.csv")
3. In the Shinyapp - Acute HBV models.Rmd file run shiny::runApp('Acute_HBV_meanprofiles')

To see the 6 patients profiles for the different models:  
1. Download the ZIP.
2. In the server.R file, don't forget to set the correct directory of the observation file:
   patient1_obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_individualprofiles/Dataset/patient1-observationV.csv")
   patient2_obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_individualprofiles/Dataset/patient2-observationV.csv")
   patient3_obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_individualprofiles/Dataset/patient3-observationV.csv")
   patient5_obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_individualprofiles/Dataset/patient5-observationV.csv")
   patient6_obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_individualprofiles/Dataset/patient6-observationV.csv")
   patient7_obsV <- read.csv("C:/Users/YourUser/.../Acute-HBV-models-main/Shiny app/Acute_HBV_individualprofiles/Dataset/patient7-observationV.csv")
3. In the Shinyapp - Acute HBV models.Rmd file run shiny::runApp('Acute_HBV_individualprofiles')

