readCode <- function(pop_name){
  fileName <- paste0(pop_name,".txt")
  code <- readChar(fileName, file.info(fileName)$size)
  code
}


