.onAttach <- function(...) {
  packageStartupMessage(
    appendLF = TRUE,
    stringr::str_wrap(
      paste0(
        "Thanks for using networker v", utils::packageVersion("networker"),
        "! If you encounter any bugs or problems, please submit an issue at ",
        "the Github page:"
      ),
      width = getOption("width")
    ),
    "\nhttps://github.com/travis-m-blimkie/networker/issues\n"
  )
}
