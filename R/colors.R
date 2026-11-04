#' Colour palettes used by getBio plotting functions
#'
#' @format Character vectors of colours.
#' \describe{
#'   \item{color_used}{Combined ggsci palettes (npg, igv, uchicago, futurama, aaas), used for Seurat plots.}
#'   \item{color_crc}{41 hand-picked named R colours, used for group/cluster colouring.}
#' }
#' @name getBio-colors
NULL

#' @rdname getBio-colors
#' @export
color_used <- c(ggsci::pal_npg()(10), ggsci::pal_igv()(9), ggsci::pal_uchicago("light")(9),
                ggsci::pal_futurama()(12), ggsci::pal_aaas()(10))[-8]

#' @rdname getBio-colors
#' @export
color_crc = c('darkcyan','mediumpurple1','darkblue','aquamarine','orchid4',
              'brown4','mediumorchid1','maroon3','antiquewhite3',
              'lightcoral','deeppink','royalblue','darkmagenta',
              'cornflowerblue','plum1','moccasin','chartreuse4',
              'chocolate2','indianred3','darkslateblue','mediumturquoise',
              'lightpink','deeppink4','deepskyblue','mediumaquamarine',
              'yellow4','deepskyblue4','steelblue','cyan2','brown1',
              'lightpink4','khaki','hotpink3','burlywood4','darkslategray',
              'lightblue','peru','lemonchiffon2','blueviolet','sienna','magenta')
