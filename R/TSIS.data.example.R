#' Save TSIS.data data into csv files
#'
#' @param dir directory to save data. If the directory does not exist, a new folder will be created with the provided name.
#'
#' @return Three csv tables of corresponding data in TSIS.data.
#' @export
#' @examples
#' TSIS.data.example(dir='data')


TSIS.data.example<-function(dir='.'){
  # if(!file.exists(dir))
  #   dir.create(dir)
  
  download.file('https://github.com/wyguo/TSIS/raw/master/data/example_data.zip',
                destfile = 'example data.zip',quiet = T)
  unzip('example data.zip',exdir = '.')
  invisible(file.remove('example data.zip'))
  
  return('Done!')
}
