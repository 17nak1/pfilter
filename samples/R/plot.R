#########################################            data
cureentfolder <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))

resr <- read.csv(file.path(cureentfolder, "resr.txt"))
resjs <- read.csv(file.path(cureentfolder, "resjs.txt"))
#########################################       make pomp

boxplot(resr$x, resjs$x)

