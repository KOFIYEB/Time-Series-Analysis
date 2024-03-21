minmax <- function(x) {
  xmin = min(x)
  xmax = max(x)
  dminmax = xmax - xmin
  xstd = (x - xmin) / dminmax
  return(list(xstd = xstd,
              dminmax = dminmax,
              xmin = xmin)
  )
}