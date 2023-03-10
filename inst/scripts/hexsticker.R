library(hexSticker)
library(ggplot2)

p <- ggplot(aes(x = mpg, y = wt), data = mtcars) + geom_point()
p <- p + theme_void() + theme_transparent()
outfile <- tempfile(fileext=".png")
s <- sticker(p, package="hexSticker", filename=outfile)

plot(s)