resolve_na_action <- function(na.action) {

  # NULL -> follow options like lm()
  if (is.null(na.action)) {
    na.action <- getOption("na.action")
  }

  # 1) function input (recommended form)
  if (is.function(na.action)) {
    if (identical(na.action, stats::na.omit))    return("omit")
    if (identical(na.action, stats::na.exclude)) return("exclude")
    if (identical(na.action, stats::na.fail))    return("fail")
    if (identical(na.action, stats::na.pass))    return("pass")
    stop("Unsupported na.action function. Use na.omit/na.exclude/na.fail/na.pass.", call. = FALSE)
  }

  # 2) character input
  if (is.character(na.action) && length(na.action) == 1L) {
    x <- tolower(na.action)

    if (x %in% c("na.omit", "omit"))       return("omit")
    if (x %in% c("na.exclude", "exclude")) return("exclude")
    if (x %in% c("na.fail", "fail"))       return("fail")
    if (x %in% c("na.pass", "pass"))       return("pass")

    stop("Invalid na.action string. Use 'na.omit','na.exclude','na.fail','na.pass' (or omit/exclude/fail/pass).",
         call. = FALSE)
  }

  stop("Invalid na.action type. Use a function (na.omit/na.exclude/na.fail/na.pass), a string, or NULL.",
       call. = FALSE)
}
