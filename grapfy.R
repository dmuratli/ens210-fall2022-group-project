
bambu <- read.table('most_c_3.csv', sep=',', h=F)
names(bambu) <- c('position', 'characteristic', 'percentage')

> bambu_ic <- ggplot(bambu, aes(position, percentage))
> bambu_ic + geom_label( aes( label = characteristic ), nudge_x=1, nudge_y =1 )