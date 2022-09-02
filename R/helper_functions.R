sumna = function(x) {
  if(all(is.na(x))) NA else sum(x, na.rm = TRUE)
}
`%!in%` = Negate(`%in%`)
