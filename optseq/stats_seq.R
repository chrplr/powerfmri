pdf("durations_of_null_events.pdf")
par(mfrow=c(2,2))
for (f in Sys.glob('*.par'))
{
    a <- read.table(f, sep='', dec='.', quote='')
    dotchart(a$V3[a$V5=="NULL"], main=f)
}
graphics.off()
