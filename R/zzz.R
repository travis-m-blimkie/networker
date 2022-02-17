.onAttach <- function(...) {
  packageStartupMessage(
    stringr::str_wrap(
      paste0(
        "Thanks for using networker v", utils::packageVersion("networker"),
        "! If you encounter any bugs or problems, please submit an issue at ",
        "the Github page: ",
        "https://github.com/travis-m-blimkie/networker/issues"
      ),
      width = getOption("width")
    )
  )
}
