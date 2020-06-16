
cureentfolder <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))

setwd(cureentfolder)
res <- c()
for (val in c(1:50000))
{
  res <- c(res,rpois(1,0.5))
}

write.csv(res, file.path(cureentfolder, "rpoisr.txt"))


resr <- read.csv(file.path(cureentfolder, "rpoisr.txt"))
resjs <- read.csv(file.path(cureentfolder, "rpoisjs.txt"))
#########################################       make pomp

boxplot(resr$x, resjs$x)



