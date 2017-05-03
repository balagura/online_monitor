## Author V. Balagura, balagura@cern.ch (19.11.2012)

for (pack in options("defaultPackages")[[1]]) suppressPackageStartupMessages(require(pack, character.only = TRUE)) # load default
## libs if running from analog of .Rprofile (no harm otherwise)

library(ggplot2, quietly = TRUE)
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
library(RGtk2, quietly = TRUE) # if put before pipe, pipe hangs
library(e1071, quietly = TRUE) # for kurtosis
## setwidth is no more available in CRAN
## library(setwidth, quietly = TRUE) # to automatically adjust R output to terminal width
library(cairoDevice, quietly = TRUE)

if ('colorout' %in% rownames(installed.packages())) library(colorout)

online.monitor.dir <- Sys.getenv("ONLINE_MONITOR_DIR")
pedestal.suppression <- Sys.getenv("ONLINE_MONITOR_PEDESTAL_SUPPRESSION") # can be '' (empty)

file.with.user.defined.plots <- '.online_monitor_user_defined_plots.RData'

online <- ifelse(system('hostname',intern=TRUE) == 'llrcaldaq.in2p3.fr', TRUE, FALSE)
FSM <- FALSE

file.suffix <- '_by_dif'
dif.glob <- 1
file.name <- Sys.getenv("ONLINE_MONITOR_FILE") # can be '' (empty)

ev <- NULL
ev.chip <- NULL
hits <- NULL
sca.glob <- 'All' # 1
chip.glob <- 'All' # 0
cut.glob <- ''


cut <- function(all.chips=FALSE, all.scas=FALSE) {
    res <- if (cut.glob == '') 'TRUE' else paste0('(',cut.glob,')') # brackets for OR: (cut1 | cut2) & chip==1
    if (all.chips == FALSE & chip.glob != 'All') {
	res <- paste0(res,' & chip==',chip.glob)
    }
    if (all.scas == FALSE & sca.glob != 'All') {
	res <- paste0(res,' & sca==',sca.glob)
    }
    res
}
cut.expr <- function(all.chips=FALSE, all.scas=FALSE) { parse(text=cut(all.chips, all.scas)) }

pdf.file.dir <- paste0(online.monitor.dir,'/plots')

map <- fread(paste0(online.monitor.dir,'/fev10_chip_channel_x_y_mapping.txt')) # data.table with chip:channel <-> X-Y mapping
## load(file=paste0(online.monitor.dir,'/map.RData'))
map[,i:=channel]

load.raw <- function(file) {
    hits <<- NULL
    ev <<- NULL
    ev.chip <<- NULL

    cat('Loading file ',file,'\n')
    df.names <- c('acq','bx','sca','chip','i','adc','trig')
    df.classes <- c(rep('integer',6), 'logical')
    hits.local <- data.table(read.table(pipe(paste0(online.monitor.dir,'/raw/raw ',file,' 0 high_gain_triggers ', pedestal.suppression)),
				     col.names=df.names, colClasses=df.classes, comment.char=''),
			  key='acq,chip,sca,bx')
    if (nrow(hits.local) > 0) {
	hits.local[,bx.cor:=cumsum( {
	    dbx = c(0,diff(bx))
	    dbx < 0 | ( dbx == 0 & c(0,diff(sca))>0 )
	       })*4096+bx, by=list(acq,chip)] # For each spill,chip: sort in SCA: bx.cor = bx + N*4096, where N - number of detected oveflows,
					# ie. number of times when BX goes down, eg. maximally: from 4095 to 0.
					# Note, if bx jumps by >4096 crossings, this is not correct (but this is the best one can do with only 12 bits for BX).
					# dbx < 0 | ( dbx == 0 & c(0,diff(sca))>0 ): treat a jump by 4096 BXs as a special case:
					# in this case BX stays the same (dbx==0) but SCA changes

	## bx.cor correction is done per chip; BX in different chips may still be compared only modulo(4096): eg. in the same event the
	## jump may be detected in one chip, but not in the other. Ie. across chips one uses bx (== bx.cor modulo(4096)), within the chip - bx.cor.

	ev.local <- hits.local[,list(n.trig=sum(trig)), keyby=list(acq,bx)]  # Find retriggers across all chips in dif
	ev.local[,bx.group:=cumsum( c(0,diff(bx)) > 4 ), by=list(acq)]  # c(0,diff(bx)) = distances to previous BX;
	## nice trick to group bxid's differing by at most 3 ("successive")
	ev.local[,`:=`(nbx=.N, ibx=1:.N), by=list(acq,bx.group)]   # finally, number of "successive" BXs in each group, ibx>1 means retriggerings

	## In FEV10, sometimes there are events with negative signals (eg. ADC==4). Pedestals in such events
	## are shifted. To avoid errors, pedestals are calculated only if there is no any channel with "negative" signal in the same event.
	##
	## Algorithm: first, find approximate pedestal values (biased by "negative" signals) from untriggered and not retriggered data per chip, channel, SCA.
        ## Then, take a median pedestal inside every chip and, finally, find the minimum across the chips.
	## In the following "unbiased" pedestal calculation, only entries above (0.75 * this value) will be considered.
	## Note, for the chip with the minimal pedestals, this assumes that all its pedestals > 0.75 * (chip average).
	negative.adc.threshold <- 0.75 * min(merge(hits.local,
						   ev.local[,list(acq,bx,nbx)], # add nbx to hits.local to require nbx==1
						   by=c('acq','bx'))[trig==FALSE & nbx==1
						       ][, list(ped=as.double(median(adc))), by=list(chip,i,sca)
							 ][,list(ped=median(ped)),by=list(chip)]$ped)
	## cat('Negitive threshold:',negative.adc.threshold, '\n')

	## same per chip (name with .chip), events for one chip are sorted first in SCA
	ev.local.chip <- hits.local[,list(n.trig.chip=sum(trig),
					  n.neg.trig.chip=sum(adc<negative.adc.threshold)), keyby=list(acq,chip,sca,bx.cor,bx)]
	ev.local.chip[,bx.group.chip:=cumsum( c(0,diff(bx.cor)) > 4 ), by=list(acq,chip)]
	ev.local.chip[,`:=`(nbx.chip=.N, ibx.chip=1:.N), by=list(acq,chip,bx.group.chip)]

	## do the same, but taking into account only events with at least one trigger
	## (in SKIROC, if trigger arrives at the clock edge, both BX and BX+1 are written, but the latter without any triggers, if there is no retriggering)
	any.event.wo.trig <- any(ev.local.chip$n.trig.chip == 0)
	if (any.event.wo.trig) {
	    ev.local.trig <- ev.local.chip[n.trig.chip>0, list(acq,chip,sca,bx.cor,bx)] # remove events wo triggers, name with suffix .trig;
	    ## note: this automatically means per chip
	    ev.local.trig[,bx.group.trig:=cumsum( c(0,diff(bx.cor)) > 4 ), by=list(acq,chip)]  # find retriggers per chip AND only for events with triggers
	    ev.local.trig[,`:=`(nbx.trig=.N, ibx.trig=1:.N), by=list(acq,chip,bx.group.trig)]
	    ev.local.chip <- ev.local.trig[ev.local.chip]
	} else ev.local.chip[,`:=`(bx.group.trig=bx.group.chip, nbx.trig=nbx.chip, ibx.trig=ibx.chip)]

	ev.local.chip <- merge(ev.local, ev.local.chip, by=c('acq','bx'), all=TRUE)
	hits.local <- merge(hits.local, ev.local.chip, by=c('acq','chip','sca','bx.cor','bx'))

#        hits.local <- hits.local[ibx==1] # remove retriggers

	hits.local <- merge(hits.local, map[,list(chip,i,x,y)], by=c('chip','i'), all.x=TRUE)

	## subtract pedestals, find them from trig==FALSE & ibx==1 & (no channel in the same chip with ADC < negative.adc.threshold), otherwise set to NA
        ## change in Feb.2016: require nbx==1 as well to remove second pedestal peak due to retrigger, then ibx==1 automatically
	hits.local[,a := {
		       unbiased.pedestal <- trig==FALSE & n.neg.trig.chip == 0 & nbx==1
		       as.double(adc) - if(any(unbiased.pedestal)) median(adc[unbiased.pedestal]) else NA
		     }, by=list(chip,i,sca)]
	setkey(hits.local, acq, chip, bx.cor, i)
	setkey(ev.local.chip, acq, chip, bx.cor)
	setkey(ev.local, acq, bx)

	hits <<- hits.local
	ev <<- ev.local
	ev.chip <<- ev.local.chip
    } else display.no.data()
}

win <- gtkWindow(show = FALSE)
gtkWindowSetTitle(win, 'Online_monitor  (questions+comments -> balagura@cern.ch)')
graphics <- gtkDrawingArea()
invisible(asCairoDevice(graphics))

daq.state.selection <- gtkComboBoxNewText()
plot.selection <- gtkComboBoxNewText()
chip.selection <- gtkComboBoxNewText()
dif.selection <- gtkComboBoxNewText()
sca.selection <- gtkComboBoxNewText()
cut.selection <- gtkEntry()
pdf <- gtkButton(label='PDF')
help.button <- gtkButton(label='Help')
open.file <- gtkButton(label='Open file')
reload <- gtkButton(label='Reload')

plots <-list()
plots[['64 ADC']] <- function() {
    d <- hits[eval(cut.expr())]
    d <- rbind(d[,all.trig:='all'],d[trig==TRUE][,all.trig:='trig'])
    bin <- max(as.integer(range(d$adc) / 200), 1)
    ggplot(d)+
	geom_histogram(aes(x=adc,color=all.trig,fill=all.trig),
		       position='identity',alpha=0.3,binwidth=bin)+facet_wrap(~i)+
			   labs(x='ADC, channels',y='Counts')+
			       scale_color_hue(name='Events')+scale_fill_hue(name='Events')
}
plots[['64 (ADC-pedestal)']] <- function() {
    d <- hits[eval(cut.expr())]
    d <- rbind(d[,all.trig:='all'],d[trig==TRUE][,all.trig:='trig'])
    bin <- max(as.integer(range(d$adc) / 200), 1)
    ggplot(d)+
	geom_histogram(aes(x=a,color=all.trig,fill=all.trig),
		       position='identity',alpha=0.3,binwidth=bin)+facet_wrap(~i)+
			   labs(x='Pedestal subtracted ADC, channels',y='Counts')+
			       scale_color_hue(name='Events')+scale_fill_hue(name='Events')
}
plots[['64 pedestals']] <- function() {
    trig.channels <- hits[trig==TRUE, list(n=.N), keyby=list(chip,i)]
    ## no cut() for trig.channels selection
    d <- hits[eval(cut.expr()) & trig==FALSE & nbx==1 & n.neg.trig.chip == 0,
	      list(adc,trig,sca), keyby=list(chip,i)]
    d <- trig.channels[d]
    pedestals <- d[,list(ped=median(as.double(adc))), keyby=list(chip,i,sca)]
    limits <- d[,list(adc.mean=median(as.double(adc)),
		      adc.mad =mad   (as.double(adc))), keyby=list(chip,i)
		][,{
		    rng <- quantile(adc.mean, c(0.01, 0.99))
		    avr.width <- median(adc.mad)
		    list(left =rng[1] - 20 * avr.width,
			 right=rng[2] + 20 * avr.width)
		    }]
    d <- d[limits$left < adc & adc < limits$right]
    ## remove channels without triggers (masked) as trig=F contains hits there
    ## note, however, sometimes triggers appear even in masked channels
    bin <- max(round(diff(range(d$adc)) / 100), 1)
    pedestals <- pedestals[,list(rng=c(min(ped),max(ped))), keyby=i]
    qplot(data=d, x=adc, binwidth=bin, facets=~i,
	  xlab='Pedestal, ADC counts', ylab='Counts', color=I('black')) +
	geom_vline(data=pedestals[limits$left < rng & rng < limits$right],
		   aes(xintercept=rng), color='red')
}
plots[['Channels w/trig']] <- function() {
    trig.channels <- hits[trig==TRUE, list(n=.N), keyby=list(chip,i)] # no cut() for trig.channels selection
    d <- hits[eval(cut.expr()), list(adc, trig), keyby=list(chip,i)]
    d <- d[trig.channels]
    d <- rbind(d[,all.trig:='all'],d[trig==TRUE][,all.trig:='trig'])
    ggplot(d)+
	geom_histogram(aes(x=adc,color=all.trig,fill=all.trig),
		       position='identity',alpha=0.3,binwidth=1)+facet_wrap(~i)+
			   labs(x='ADC, channels',y='Counts')+
			       scale_color_hue(name='Events')+scale_fill_hue(name='Events')
}
plots[['N "successive" SCA']] <- function() {
    qplot(data=ev.chip[eval(cut.expr(all.scas=TRUE)) & ibx.trig==1],nbx.trig, facets=~chip, xlab='N "successive" SCA',ylab='Counts',
	main='',binwidth=1)
}
plots[['N trigs(ch), ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==TRUE & ibx==1][,list(n.trig=.N),by=list(i,chip)][n.trig>0]
    qplot(data=d, i, n.trig, facets=~chip, xlab='Channel',ylab='N')
}
plots[['N trigs(ch,chip), ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==TRUE & ibx==1][,list(n.trig=.N),by=list(i,chip)][n.trig>0]
    qplot(data=d, i,chip,fill=n.trig,
	  geom='tile',color=I('darkgreen'),xlab='Channel',ylab='Chip') + scale_fill_gradient(low="green", high="red",name='N trigs')
}
plots[['N trigs(X,Y), ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==TRUE & ibx==1][,list(n.trig=.N),by=list(x,y,i,chip)][n.trig>0]
    qplot(data=d, x,y,fill=n.trig,geom='tile',color=I('darkgreen'),xlab='X',ylab='Y') +
	scale_fill_gradient(low="green", high="red",name='N trigs') +
	    geom_text(aes(label=paste0(chip,':',i),alpha=chip),color=I('black'),fontface=I(2),size=I(3))+scale_alpha(range=c(0.3,0.8))
}
plots[['<ADC>(ch,chip)']] <- function() {
    qplot(data=hits[eval(cut.expr()), list(n.trig=sum(trig),ADC=mean(adc)), by=list(chip,i)],i,chip,fill=ADC,
	geom='tile',color=I('darkgreen'),xlab='Channel (* means there are triggers)',ylab='Chip') + scale_fill_gradient(low="green", high="red") +
    geom_text(aes(label=ifelse(n.trig>0,'*','')),color=I('black'),fontface=I(2),size=I(3))
}
plots[['<ADC>(X,Y)']] <- function() {
    qplot(data=hits[eval(cut.expr()), list(n.trig=sum(trig),ADC=mean(adc),x=x[1],y=y[1]), by=list(chip,i)],x,y,fill=ADC,
	geom='tile',color=I('darkgreen'),xlab='X (* means there are triggers)',ylab='Y') + scale_fill_gradient(low="green", high="red") +
    geom_text(aes(label=paste0(chip,':',i,ifelse(n.trig>0,',*','')),alpha=chip),color=I('black'),fontface=I(2),size=I(3))+scale_alpha(range=c(0.3,0.8))
}
plots[['<Trigger-ped.>(ch,chip), ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==TRUE &  ibx==1 & !is.na(a)][,list(mean.trig=mean(a)),by=list(i,chip)]
    if (nrow(d)>0) {
	qplot(data=d,
	      i,chip,fill=mean.trig,geom='tile',color=I('darkgreen'),xlab='Channel',ylab='Chip') + scale_fill_gradient(low="green", high="red",name='<ADC-ped.>')
    } else display.no.data()
}
plots[['<Trigger-ped.>(X,Y), ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==TRUE &  ibx==1 & !is.na(a)][,list(mean.trig=mean(a)),by=list(x,y,chip)]
    if (nrow(d)>0) {
	qplot(data=d,
	      x,y,fill=mean.trig,geom='tile',color=I('darkgreen'),xlab='X',ylab='Y') + scale_fill_gradient(low="green", high="red",name='<ADC-ped.>') +
    geom_text(aes(label=paste0(chip,':',i),alpha=chip),color=I('black'),fontface=I(2),size=I(3))+scale_alpha(range=c(0.3,0.8))
    } else display.no.data()
}
plots[['Trig-ped. sum(X,Y), ibx=1']] <- function() {
  d <- hits[eval(cut.expr()) & trig==TRUE & ibx==1 & !is.na(a)][,list(val=sum(a)),by=list(x,y,chip,i)][val>0]
  qplot(data=d, x,y,fill=val,geom='tile',color=I('darkgreen'),xlab='X',ylab='Y') +
    scale_fill_gradient(low="green", high="red",name='Trig-ped sum') +
    geom_text(aes(label=paste0(chip,':',i),alpha=chip),color=I('black'),fontface=I(2),size=I(3))+scale_alpha(range=c(0.3,0.8))
}
plots[['RMS(ch,chip)']] <- function() {
    qplot(data=hits[eval(cut.expr()), list(n.trig=sum(trig),RMS=sd(adc)), by=list(chip,i)],i,chip,fill=RMS,
	geom='tile',color=I('darkgreen'),xlab='Channel (* means there are triggers)',ylab='Chip') + scale_fill_gradient(low="green", high="red") +
    geom_text(aes(label=ifelse(n.trig>0,'*','')),color=I('black'),fontface=I(2),size=I(3))
}
plots[['RMS(X,Y)']] <- function() {
    qplot(data=hits[eval(cut.expr()), list(n.trig=sum(trig),RMS=sd(adc),x=x[1],y=y[1]), by=list(chip,i)],x,y,fill=RMS,
	geom='tile',color=I('darkgreen'),xlab='X (* means there are triggers)',ylab='Y') + scale_fill_gradient(low="green", high="red") +
    geom_text(aes(label=paste0(chip,':',i,ifelse(n.trig>0,',*','')),alpha=chip),color=I('black'),fontface=I(2),size=I(3))+scale_alpha(range=c(0.3,0.8))
}

plots[['<Pedestals>(ch)']] <- function() {
  qplot(data=hits[eval(cut.expr()) & !is.na(a)][,list(ped=adc[1]-a[1]),keyby=list(chip,i)],i,ped,xlab='Channel',ylab='Pedestal mean, ADC counts',facets=~chip)
}
plots[['<Pedestal>(ch) VS SCA']] <- function() {
    qplot(data=hits[eval(cut.expr()) & trig==FALSE & nbx==1 & n.neg.trig.chip == 0
                    ][,list(ped=mean(adc)),keyby=list(sca,chip,i)],i,ped,color=factor(sca),facets=~chip,
	xlab='Channel',ylab='Pedestal mean, ADC counts') + scale_color_hue(name='SCA')
}
plots[['Pedestal RMS(ch) per chip']] <- function() {
    d <- hits[eval(cut.expr()) & trig==FALSE & nbx==1 & n.neg.trig.chip == 0
              ][,a:=adc-mean(adc),by=list(chip,i,sca)][,N:=.N,by=list(chip,i)][N>100][,{
    rms <- sd(a)
    kur <- kurtosis(a, type=2) # type=2 is unbiased for normal distributions (see ?kurtosis from e1071 package)
    list(rms=rms,
	 rms.e=sqrt(2/(.N-1) + kur/.N)*rms/2) # var(rms^2) = rms^4 (2/(N-1) + kur/N), sd(rms) = var(rms^2) / 2 / rms, from
  }, by=list(chip,i)]                         # http://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
  qplot(data=d, i, rms, facets=~chip, xlab='Channel',ylab='RMS, ADC channels')+ geom_pointrange(data=d,aes(ymin=rms-rms.e,ymax=rms+rms.e))
}
plots[['Pedestal RMS(ch,chip)']] <- function() {
    qplot(data=hits[eval(cut.expr()) & trig==F & nbx==1 & n.neg.trig.chip == 0
                    ][, list(a=adc - mean(adc)), by=list(chip,i,sca)][,list(RMS=sd(a)),by=list(chip,i)],
	  i,chip,fill=RMS,
	  geom='tile',color=I('darkgreen'),xlab='Channel',ylab='Chip') + scale_fill_gradient(low="green", high="red")
}
plots[['Pedestal RMS(X,Y)']] <- function() {
    qplot(data=hits[eval(cut.expr()) & trig==F & nbx==1 & n.neg.trig.chip == 0
                    ][, list(a=adc - mean(adc),x,y), by=list(chip,i,sca)][,list(RMS=sd(a),x=x[1],y=y[1]),by=list(chip,i)],
	  x,y,fill=RMS,
	  geom='tile',color=I('darkgreen'),xlab='X',ylab='Y') + scale_fill_gradient(low="green", high="red") +
    geom_text(aes(label=paste0(chip,':',i),alpha=chip),color=I('black'),fontface=I(2),size=I(3))+scale_alpha(range=c(0.3,0.8))
}
plots[['MIP, ibx=1']] <- function() {
  d <- hits[eval(cut.expr()) & trig==TRUE & ibx==1 & !is.na(a)]
  rng <- quantile(d$a, probs=c(0.005,0.995))
  a.min <- as.integer(rng[1]/10)*10
  a.max <- as.integer(rng[2]/10)*10
  qplot(data=d[a.min<a & a<a.max],a,binwidth=ceiling(a.max/200), xlab='Trigger-pedestal, ADC channels',ylab='Counts')
}
plots[['MIP per chip, ibx=1']] <- function() {
  d <- hits[eval(cut.expr()) & trig==TRUE & ibx==1 & !is.na(a)]
  rng <- quantile(d$a, probs=c(0.005,0.995))
  a.min <- as.integer(rng[1]/10)*10
  a.max <- as.integer(rng[2]/10)*10
  qplot(data=d[a.min<a & a<a.max],a,binwidth=ceiling(a.max/200), facets=~chip, xlab='Trigger-pedestal, ADC channels',ylab='Counts')
}
plots[['MIP per channel, ibx=1']] <- function() {
  d <- hits[eval(cut.expr()) & trig==TRUE & ibx==1 & !is.na(a)]
  rng <- quantile(d$a, probs=c(0.005,0.995))
  a.min <- as.integer(rng[1]/10)*10
  a.max <- as.integer(rng[2]/10)*10
  qplot(data=d[a.min<a & a<a.max],a,binwidth=ceiling(a.max/200), xlab=paste0('Trigger-pedestal, ADC channels'),ylab='Counts') + facet_wrap(facets=~i, scale='free_y')
}
plots[['BX']] <- function() {
    e <- ev[eval(cut.expr())]
    bx.rng <- range(e$bx)
    bin <- max(as.integer(diff(bx.rng)/200), 1)
    qplot(data=e,bx,binwidth=bin, xlab='BX',ylab='Events')
}
plots[['BX per chip']] <- function() {
    e <- ev.chip[eval(cut.expr())]
    bx.rng <- range(e$bx)
    bin <- max(as.integer(diff(bx.rng)/200), 1)
    qplot(data=e,bx,binwidth=bin,facets=~chip, xlab='BX',ylab='Events')
}
plots[['BX (w recycles) per chip']] <- function() {
    e <- ev.chip[eval(cut.expr())]
    bx.rng <- range(e$bx.cor)
    bin <- max(as.integer(diff(bx.rng)/200), 1)
    qplot(data=e,bx.cor,binwidth=bin,facets=~chip, xlab='BX, 4096 limit partially corrected',ylab='Events')
}
plots[['BX (w recycles) per (chip,SCA)']] <- function() {
    e <- ev.chip[eval(cut.expr())]
    bx.rng <- range(e$bx.cor)
    bin <- max(as.integer(diff(bx.rng)/200), 1)
    qplot(data=e,bx.cor,binwidth=bin,facets=sca~chip, xlab='BX, 4096 limit partially corrected',ylab='Events')
}
plots[['Spill number']] <- function() {
    e <- ev[eval(cut.expr())]
    rng <- range(e$acq)
    bin <- max(as.integer(diff(rng)/200), 1)
    qplot(data=e, acq, binwidth=bin, xlab='Spill', ylab='Events')
}
plots[['SCA per chip']] <- function() {
    qplot(data=ev.chip[eval(cut.expr(all.scas = TRUE))],
	  sca, facets=~chip, binwidth=1, xlab = 'SCA', ylab='Counts')
}
plots[['SCA per chip, ibx=1']] <- function() {
    qplot(data=ev.chip[eval(cut.expr(all.scas = TRUE)) & ibx.trig==1],
	  sca, facets=~chip, binwidth=1, xlab = 'SCA', ylab='Counts')
}
plots[['Ped. per sum(a), 1 chip, ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & ibx==1 & n.neg.trig.chip == 0]
    ## if chip.glob is given, it is selected in cut.expr() above, otherwise:
    if (chip.glob == 'All') { # select the most busy chip
      d.sum <- d[,list(a=sum(a[trig==T & !is.na(a)])), by=chip]
      chip.shown <- d.sum[which.max(a)]$chip
      d <- d[chip == chip.shown]
    } else {
	chip.shown = chip.glob
    }
    d[,a.tot:=sum(a[trig==T & !is.na(a)]), by=list(acq,sca)]
    a.tot.bin <- round(max(d$a.tot) / 10, -1)
    d <- d[trig==FALSE & nbx==1]
    d[,a.tot.bin:=round(a.tot / a.tot.bin) * a.tot.bin, by=list(acq,sca)]
    d <- d[,a.tot.bin.stat:=.N, by=list(a.tot.bin)]
    res <- d[a.tot.bin.stat > 0.01*nrow(d) & a.tot.bin.stat > 50]
    bin <- max(1, round(diff(range(res$adc)) / 200))
    qplot(data=res, adc,
	  binwidth = bin,
	  xlab='Pedestal, ADC channels',
	  ylab='Counts',main=paste0('Chip ',chip.shown),color=I('black'))+
	  facet_wrap(~a.tot.bin, scale='free_y') +
	geom_vline(aes(xintercept=median(adc)),color='red')

}
## plots[['Ped. per sum(a) per chip, ibx=1']] <- function() {
##     d <- hits[eval(cut.expr()) & ibx==1 & n.neg.trig.chip == 0]
##     ## if chip.glob is given, it is selected in cut.expr() above, otherwise:
##     if (chip.glob == 'All') { # select the most busy chip
##       d.sum <- d[,list(a=sum(a[trig==T & !is.na(a)])), by=chip]
##       chip.shown <- d.sum[which.max(a)]$chip
##       d <- d[chip == chip.shown]
##     } else {
##         chip.shown = chip.glob
##     }
##     d[,a.tot:=sum(a[trig==T & !is.na(a)]), by=list(acq,sca)]
##     a.tot.bin <- round(max(d$a.tot) / 10, -1)
##     d <- d[trig==FALSE & nbx==1]
##     d[,a.tot.bin:=round(a.tot / a.tot.bin) * a.tot.bin, by=list(acq,sca)]
##     d <- d[,a.tot.bin.stat:=.N, by=list(a.tot.bin)]
##     res <- d[a.tot.bin.stat > 0.01*nrow(d) & a.tot.bin.stat > 50,
## 	     list(pedestal = median(as.double(adc))), by=list(chip,i,sca)] # by=list(chip,i,sca, a.tot.bin)]
##     bin <- max(1, round(diff(range(res$pedestal)) / 200))
##     qplot(data=res, pedestal,
## 	  binwidth = bin,
## 	  xlab='Pedestal, ADC channels',
## 	  ylab='Counts',main=paste0('Chip ',chip.shown))+
## 	  facet_wrap(~a.tot.bin, scale='free_y')

## }
plots[['MIP VS sum(a) per chip, ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==T & ibx==1 & !is.na(a)
	      & n.neg.trig.chip == 0]
    d[,a.tot:=sum(a), by=list(acq,chip,sca)]
    rng <- quantile(d$a, probs = c(0.0, 0.70))
    a.min <- 0 #as.integer(rng[1]/10) * 10
    a.max <- as.integer(rng[2]/10) * 10
    a.tot.bin <- max(d$a.tot) / 10
    d[,a.tot.bin:=round(a.tot / a.tot.bin) * a.tot.bin, by=list(acq,chip,sca)]
    d[,a.tot.n:=sum(a.min < a & a < a.max), by=chip]
    d <- d[,a.tot.bin.stat:=.N, by=list(chip,a.tot.bin)][a.tot.bin.stat > 0.01*a.tot.n]
    qplot(data=d[a.min < a & a < a.max], a.tot.bin, a,
	  geom='violin',
	  facets=~chip, group=a.tot.bin, xlab = 'sum(ADC-pedestal) over chip',
	  ylab='Trigger-pedestal, ADC channels')
}
plots[['MIP per sum(a), 1 chip, ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==T & ibx==1 & !is.na(a)
	      & n.neg.trig.chip == 0]
    ## if chip.glob is given, it is selected in cut.expr() above, otherwise:
    if (chip.glob == 'All') { # select the most busy chip
      d.sum <- d[,list(a=sum(a)), by=chip]
      chip.shown <- d.sum[which.max(a)]$chip
      d <- d[chip == chip.shown]
    } else {
      chip.shown <- chip.glob
    }
    d[,a.tot:=sum(a), by=list(acq,sca)]
    rng <- quantile(d$a, probs = c(0.0, 0.70))
    a.min <- 0 #as.integer(rng[1]/10) * 10
    a.max <- as.integer(rng[2]/10) * 10
    a.tot.bin <- round(max(d$a.tot) / 10, -1)
    d[,a.tot.bin:=round(a.tot / a.tot.bin) * a.tot.bin, by=list(acq,sca)]
    d <- d[,a.tot.bin.stat:=.N, by=list(a.tot.bin)]
    dd <- d[a.tot.bin.stat > 0.01*nrow(d)]
    qplot(data=d[a.min < a & a < a.max & a.tot.bin.stat>50], a,
	  binwidth = ceiling(a.max/200),
	  xlab='Trigger-pedestal, ADC channels',
	  ylab='Counts',main=paste0('Chip ',chip.shown))+
	  facet_wrap(~a.tot.bin, scale='free_y')

}
plots[['sum(a) per chip, ibx=1']] <- function() {
    d <- hits[eval(cut.expr()) & trig==T & ibx==1 & !is.na(a)
	      & n.neg.trig.chip == 0,
	 list(a.tot=sum(a)), by=list(acq,chip,sca)]
    qplot(data=d, a.tot,
	  binwidth=diff(range(d$a.tot))/100,
	  xlab = 'sum(ADC-pedestal) over chip', ylab='Events') +
    facet_wrap(~chip, scale='free_y')
}
plots[['N trig in chip, ibx=1']] <- function() {
     qplot(data = hits[eval(cut.expr()) & ibx == 1 & trig==T
		       & n.neg.trig.chip == 0,
		       list(n=.N), by=list(acq,chip,sca)],
	   n, facets = ~chip, binwidth = 1,
	   xlab = "N triggers in chip",ylab = "Counts")
}
standard.plot.names <- names(plots)

if (file.exists(file.with.user.defined.plots)) {
    load(file.with.user.defined.plots)
    plots <- c(plots, user.defined.plots)
}
update.plots <- function() {
    old.sensitivity <- plot.selection$sensitive
    plot.selection$sensitive <- FALSE
    gSignalHandlerDisconnect(plot.selection,plot.selection.connect)

    gtkListStoreClear(gtkComboBoxGetModel(plot.selection))
    for (pos in 1:length(plots)) gtkComboBoxInsertText(plot.selection, pos-1, names(plots)[pos])
    n.active <- which('N "successive" SCA' == names(plots)) - 1
    if (length(n.active)==0) n.active <- 0
    plot.selection$active <- n.active

    plot.selection.connect <<- gSignalConnect(plot.selection, 'changed', change.plot.selection)
    plot.selection$sensitive <- old.sensitivity

    new.plots <- names(plots)
    user.defined.plots <- plots[new.plots[!(new.plots %in% standard.plot.names)]]
    if (length(user.defined.plots) > 0) save(user.defined.plots, file=file.with.user.defined.plots)
}
## To change histograms on the fly when the monitor is running, modify plots and call update.plots()

## Finite State Machine
socket.states <- c('UNDEFINED','OFF','READY','CONFIGURED','ACQUIRING')
states <- c('all OFF','init PC, Fake Power OFF','Fake Power ON','Configured','Run') # same, more self-explanatory
for (i in 1:length(states)) gtkComboBoxInsertText(daq.state.selection, position=i, text=states[i])

start.command <- function() {
  file.filter <- gtkFileFilter()
  gtkFileFilterAddPattern(file.filter, '*.raw')
  dialog <- gtkFileChooserDialog("Save File", win, "save",
				 "gtk-cancel", GtkResponseType["cancel"],
				 "gtk-save", GtkResponseType["accept"])
  gtkFileChooserSetDoOverwriteConfirmation(dialog, TRUE) # confirm overwriting
  gtkFileChooserSetCurrentFolder(dialog, dirname(file.name))
  gtkFileChooserAddFilter(dialog, file.filter)
  run <- FALSE
  if (dialog$run() == GtkResponseType["accept"]) {
    file.name <<- file.path(data.dir, basename(paste(dialog$getFilename(),file.suffix,dif.glob,'.raw',sep='')))
    run <- TRUE
  }
  dialog$destroy()
  if (run) {
    socket.command.one.arg('start_acq_det',basename(file.name))
    Sys.sleep(5)
    from <- list.files(data.dir,pattern=paste0('running_data',file.suffix,'[0-9]+.raw'),full.names = TRUE)
    to <- sub('running_data',sub(paste0(file.suffix,dif.glob,'.raw'),'',file.name),basename(from))
    file.copy(from, to, overwrite = TRUE)
    reload.clicked()
  }
}
xml.file <- '/opt/calicoes/ecal.xml'
state.commands.up   <- c(function() socket.command.one.arg('initialize_det',xml.file),
			 function() socket.command        ('power_on_det'),
			 function() socket.command.one.arg('configure_det',xml.file),
			 start.command)
state.commands.down <- c(function() {}, #socket.command        ('deinitialize_det'),
			 function() {}, #socket.command        ('power_off_det'),
			 function() socket.command        ('invalidate_det'),
			 function() socket.command        ('stop_acq_det'))

change.daq.state.selection <- function(ptr) {
  strt <- which(socket.command('get_state_det') == socket.states)
  end  <- which(gtkComboBoxGetActiveText(daq.state.selection) == states)
  if (strt == end) return() # strt == 'Configured' with many possible parameters is not considered
  else if (strt < end) {print(strt:(end-1)); for(i in strt:(end-1)) print((state.commands.up  [[i]])()) }
  else if (strt > end) {print((strt-1):end); for(i in (strt-1):end) print((state.commands.down[[i]])()) }
  new.state <- which(socket.command('get_state_det') == socket.states)
  daq.state.selection$active <- new.state-1
}

if (online & FSM) {
  ## Communicating via sockets
  ## ports <- data.table(read.table('/opt/calicoes/ports.txt',sep='=',col.names=c('name','port')),key='name')
  ## socket <- make.socket('localhost', ports[J('DET_PORT')]$port)
  ## on.exit(close.socket(socket))
  socket.command <- function(command) {
    write.socket(socket, paste('<cmd name="',command,'">','</cmd>\n', sep=''))
    return(get.socket.result())
  }
  socket.command.one.arg <- function(command, arg) {
    write.socket(socket, paste('<cmd name="',command,'">',
			       '<param>',arg,'</param>',
			       '</cmd>\n', sep=''))
    return(get.socket.result())
  }
  get.socket.result <- function() {
    res <- read.socket(socket)
    num <- (sub('<res ret_code="([0-9])">.*</res>\n','\\1',res) != '0')
    txt <- sub('<res ret_code="[0-9]">(.*)</res>\n','\\1',res)
    if (num != TRUE) stop(txt)
    return(txt)
  }
} else {
  ## For debugging without real sockets:
  state.glob <- 1
  state.list <- c(initialize_det=2,power_on_det=3,configure_det=4,start_acq_det=5,deinitialize_det=1,power_off_det=2,invalidate_det=3,stop_acq_det=4)
  socket.command <- function(command) {
    if (command=='get_state_det') return(socket.states[state.glob])
    state.glob <<- as.numeric(state.list[command])
    cat('Command',command,'\n')
  }
  socket.command.one.arg <- function(command, arg) {
    state.glob <<- as.numeric(state.list[command])
    cat('Command',command,arg,'\n')
  }
  get.socket.result <- function() return(names(state.list)[state.glob])
}
## ------------------------------
gtkComboBoxInsertText(chip.selection, position=0, text=paste('CHIP','All'))
gtkComboBoxInsertText(dif.selection, position=0, text='DIF')

gtkComboBoxInsertText(sca.selection, position=0, text=paste('SCA','All'))
for (i in 1:15) gtkComboBoxInsertText(sca.selection, position=i, text=paste('SCA',i))

gtkEntrySetText(cut.selection, '<Cut, eg. adc>500 & trig==TRUE>')

change.plot.selection <- function(ptr) {
    cmd <- gtkComboBoxGetActiveText(plot.selection)
    print(cmd)
    print(plots[[cmd]])   # print code
    if (is.null(hits)) display.no.data() else print(plots[[cmd]]()) # call appropriate function
}
change.chip.selection <- function(ptr) {
  chip.glob <<- sub('CHIP ','',gtkComboBoxGetActiveText(chip.selection))
  cat('CHIP ',chip.glob,'\n')
  change.plot.selection(0)
}
change.dif.selection <- function(ptr) {
    if (dif.selection$sensitive == TRUE) {
	dif.glob <<- as.numeric(sub('DIF ','',gtkComboBoxGetActiveText(dif.selection)))
	new.file <- sub('[0-9]+.raw$',paste0(dif.glob,'.raw'),file.name)
	print(new.file)
	if (file.exists(new.file)) {
	    file.name <<- new.file
	    reload.clicked()
	} else {
	    display.no.data()
	    update.dif.selection() # if file for this dif was deleted: update available files=difs,
					# if none is available: dif.selection will be disabled
	}
    }
}
update.chip.selection <- function() {
    chip.selection$sensitive <- FALSE
    gSignalHandlerDisconnect(chip.selection,chip.selection.connect)

    gtkListStoreClear(gtkComboBoxGetModel(chip.selection))
    gtkComboBoxInsertText(chip.selection, position=0, text=paste('CHIP','All'))
    chips <- sort(unique(ev.chip$chip))
    if (length(chips)>0) {
	for (pos in 1:length(chips)) gtkComboBoxInsertText(chip.selection, pos, paste('CHIP',chips[pos]))
	chip.selection$active <- 0
	chip.selection$sensitive <- TRUE
    }
    chip.selection.connect <<- gSignalConnect(chip.selection, 'changed', change.chip.selection)
}
update.dif.selection <- function() { # fill dif.selection with available dif numbers (from X in <this dir>/<file>_by_difX.raw).
    pat <- sub(paste0('(',file.suffix,')[0-9]+(\\.raw)$'), '\\1[0-9]+\\2', basename(file.name))
					## eg. xxx_by_dif3.raw -> xxx_by_dif[0-9]+.raw
    files <- list.files(dirname(file.name), pattern=pat)
    dif.numbers <- sub(paste0('.*',file.suffix,'([0-9]+)\\.raw'),'\\1',files)
    difs <- sort(as.integer(dif.numbers))

    dif.selection$sensitive <- FALSE
    gSignalHandlerDisconnect(dif.selection,dif.selection.connect)

    gtkListStoreClear(gtkComboBoxGetModel(dif.selection))
    if (length(difs)>0) {
	for (pos in 1:length(difs)) gtkComboBoxInsertText(dif.selection, pos-1, paste0('DIF ',difs[pos]))
	dif.selection$active <- which(difs == sub(paste0('.*',file.suffix,'([0-9]+)\\.raw$'), '\\1', basename(file.name)))-1

	dif.selection$sensitive <- TRUE
    }
    dif.selection.connect <<- gSignalConnect(dif.selection, 'changed', change.dif.selection)
}
change.sca.selection <- function(ptr) {
  sca.glob <<- sub('SCA ','',gtkComboBoxGetActiveText(sca.selection))
  cat('SCA ',sca.glob,'\n')
  change.plot.selection(0)
}
change.cut.selection <- function(ptr) {
  cut.glob <<- gtkEntryGetText(cut.selection)
  cat('Cut: ',cut.glob,'\n')
  change.plot.selection(0)
}
pdf.clicked <- function(ptr) {
  file.filter <- gtkFileFilter()
  gtkFileFilterAddPattern(file.filter, '*.pdf')
  dialog <- gtkFileChooserDialog("Open File", win, "open",
				 "gtk-cancel", GtkResponseType["cancel"],
				 "gtk-open", GtkResponseType["accept"])
  gtkFileChooserSetCurrentFolder(dialog, pdf.file.dir)
  gtkFileChooserAddFilter(dialog, file.filter)
  gtkFileChooserSetAction(dialog, 'save')
  if (dialog$run() == GtkResponseType["accept"]) {
    pdf.file.dir <<- dirname(dialog$getFilename())
    ggsave(dialog$getFilename(), width=12, height=8)
  }
  dialog$destroy()
}
help.button.clicked <- function(ptr) {
    cat(help.text)
    dialog <- gtkDialog(title='Online_monitor help', parent=win, flags="destroy-with-parent", buttons="gtk-close", GtkResponseType["close"], show=FALSE)
    gSignalConnect(dialog, "response", gtkWidgetDestroy)
    scroll <- gtkScrolledWindow()
    gtkScrolledWindowSetPolicy(scroll, vscrollbar.policy=GtkPolicyType['automatic'], hscrollbar.policy=GtkPolicyType['automatic'])
    dialog$vbox$packStart(scroll, expand = TRUE, fill = TRUE, padding = 0)

    label <- gtkLabel(help.text)
    label$wrap <- TRUE
    label$selectable <- TRUE
    gtkLabelSetMaxWidthChars(label, 90)
    ##    txt <- gtkTextView()
    ##    gtkTextBufferSetText(txt$buffer, help.text)
    ##    txt$editable <- FALSE
    ##  ##  txt$sensitive <- FALSE
    gtkScrolledWindowAddWithViewport(scroll, label)
    dialog$setSizeRequest(700,800)
    dialog$show()
}
open.new.file <- function() {
    update.dif.selection()
    invisible(lapply(gtkContainerGetChildren(hbox), function(x) x$sensitive = TRUE)) # enable all controls
    reload.clicked()
}
open.file.clicked <- function(ptr) {
  file.filter <- gtkFileFilter()
  gtkFileFilterAddPattern(file.filter, '*.raw')
  dialog <- gtkFileChooserDialog("Open File", win, "open",
				 "gtk-cancel", GtkResponseType["cancel"],
				 "gtk-open", GtkResponseType["accept"])
  gtkFileChooserSetCurrentFolder(dialog, dirname(file.name))
  gtkFileChooserAddFilter(dialog, file.filter)
  update <- FALSE
  if (dialog$run() == GtkResponseType["accept"]) {
    file.name <<- dialog$getFilename()
    update <- TRUE
  }
  dialog$destroy()
  if (update==TRUE) open.new.file()
}
reload.clicked <- function(ptr) {
    gtkButtonSetLabel(open.file, sub('.raw$','',basename(file.name)))
    print(qplot(x=1,y=1,geom='text',label='Processing file...',color=I('red'),xlab='',ylab='',size=I(20)))
    load.raw(file.name)
    update.chip.selection()
    change.plot.selection(0)
}
display.no.data <- function() {
    print(qplot(x=1,y=1,geom='text',label='No data',color=I('red'),xlab='',ylab='',size=I(20)))
}
vbox <- gtkVBox()
hbox <- gtkHBox()

vbox$packStart(hbox, expand = FALSE, fill = FALSE, padding = 0)
vbox$packStart(graphics, expand = TRUE, fill = TRUE, padding = 0)

## hbox$packStart(daq.state.selection,  expand = TRUE, fill = TRUE, padding = 0)
hbox$packStart(plot.selection,       expand = FALSE, fill = FALSE, padding = 0)
hbox$packStart(chip.selection,       expand = FALSE, fill = FALSE, padding = 0)
hbox$packStart(sca.selection,        expand = FALSE, fill = FALSE, padding = 0)
hbox$packStart(cut.selection,        expand = TRUE, fill = TRUE, padding = 0)
hbox$packStart(pdf,                  expand = FALSE, fill = FALSE, padding = 0)
hbox$packStart(help.button,          expand = FALSE, fill = FALSE, padding = 0)
hbox$packStart(dif.selection,        expand = FALSE, fill = FALSE, padding = 0)
hbox$packStart(open.file,            expand = TRUE, fill = TRUE, padding = 0)
hbox$packStart(reload,               expand = FALSE, fill = FALSE, padding = 0)

win$add(vbox)
win$setDefaultSize(1000,1000)

daq.state.selection$active <- which(socket.command('get_state_det') == socket.states) - 1
chip.selection$active <- 0 # 'All', otherwise chip.glob+1
sca.selection$active <- 0 # 'All', otherwise sca.glob
dif.selection$active <- 0 # 'All', otherwise sca.glob

invisible(lapply(gtkContainerGetChildren(hbox), function(x) if (! (x==open.file | x==help.button)) x$sensitive = FALSE)) # grey out everything except open.file

daq.state.selection.connect <- gSignalConnect(daq.state.selection, 'changed', change.daq.state.selection)
plot.selection.connect      <- gSignalConnect(plot.selection, 'changed', change.plot.selection)
chip.selection.connect      <- gSignalConnect(chip.selection, 'changed', change.chip.selection)
dif.selection.connect       <- gSignalConnect(dif.selection, 'changed', change.dif.selection)
sca.selection.connect       <- gSignalConnect(sca.selection, 'changed', change.sca.selection)
cut.selection.connect       <- gSignalConnect(cut.selection, 'activate', change.cut.selection)
pdf.connect                 <- gSignalConnect(pdf, 'clicked', pdf.clicked)
help.button.connect         <- gSignalConnect(help.button, 'clicked', help.button.clicked)
open.file.connect           <- gSignalConnect(open.file, 'clicked', open.file.clicked)
reload.connect              <- gSignalConnect(reload, 'clicked', reload.clicked)
# To disconnect:
# gSignalHandlerDisconnect(plot.selection,plot.selection.connect)
win$showAll()
display.no.data()
update.plots()
if (file.name != '') open.new.file()

help.text <- paste0(
'
This is R+GTK online_monitor program for analysis of SiW ECAL technological
prototype data. Start it with "<path>/online_monitor.sh [file.raw]", which
invokes online_monitor.R from the same "<path>" directory.

The plot is chosen from the list in the top left corner. Logical "AND" between
next 3 fields (CHIP, SCA and Cut) defines the selection.  The latter ("Cut"
entry) can contain an arbitrary R expression from the names of "hits" or
"ev(.per.chip)" data.tables (see below).  After typing a new "Cut" one MUST
terminate it with "ENTER".  Examples of cuts for "hits":
  chip==0 & i<10
  nbx==2 & ibx==2 & n.trig.chip==FALSE
  (!is.na(a) & a>30) | trig==TRUE
  acq %% 256 == 0  # modulo 256

"PDF" button allows to save the plot as a PDF file. "Open file" selects a new
file in a raw format. "DIF" button shows all available DIFs with similar file
names.  "Reload" rereads the same file.

In addition to GUI widgets, in the console one has a normal R session, which
may be used to produce more plots and to perform detailed ineteractive
analysis. One can therefore use this program both for simple online checks and
for complicated offline studies.

The list of predefined plots can be easily extended by populating the "plots"
list in the console R session, eg.
  plots[["My new plot"]] = function() qplot(data=hits, a, binwidth=5)
  update.plots()
The last function, update.plots(), synchronizes the "plots" structure with the
top-left list of selectable plots. Then, "My new plot" will appear in the list
of predefined plots and it will show the overall ADC pedestal subtracted
spectrum.

If the user defines his own plots, they are automatically stored in a file
', file.with.user.defined.plots,
'
in the current directory (in a structure "user.defined.plots") and then
reloaded in the next call, if online_monitor.sh is started from the same
directory.

When the plot is shown, the corresponding code which has been executed, is
printed in the console. It shows, what exactly is plotted and may be used as
an example to populate the "plots" list with your own functions.

To apply the cut selected in CHIP+SCA+Cut GUI fields, use eval(cut.expr())
expression, eg. instead of "hits" write hits[eval(cut.expr())]. To see the cut
in a text form, type "cut()". "cut.expr()" makes an R "expression" from this
text, which may be evaluated with "eval()". To ignore the CHIP or SCA cut, but
leave the rest (eg. when you plot SCA along an x-axis) use the flags
all.chips/all.scas=TRUE, eg. hits[eval(cut.expr(all.chips=TRUE,
all.scas=TRUE))].

When the raw file is opened, it is processed and the following data.tables
(see ?data.table) are created (note, they are fully kept in memory):
  hits - one raw per hit,
  ev - one raw per event,
  ev.chip - one raw per event and per chip.
To see what is inside, just type "hits", "ev" or "ev.chip" in the
console. This will print 5 first and last raws. To see more, eg. the first 50
hits, type hits[1:50].  To get help on any R command type ?command.

Only high gain ADC values are stored. Not triggered data is prescaled by
"pedestal.suppression" factor which is equal by default to 1/20. One may
change it on the command line: "online_monitor.sh --pedestal-suppression 0.5"
or just "online_monitor.sh -p 0.5", or change in the R session before opening
a new file: "pedestal.suppression=0.5".

"chip" counting starts from 0. "i"=0..63 means channel in the chip. "acq","bx"
are the spill and the bunch crossing numbers. "sca"=1..15.  Note, that
normally "bx" should increase with SCA SKIROC memory. The opposite indicates
that "bx" counter exceeded 12 bits limit (>4095) and was recycled to zero.
This is automatically detected by the program which then increments "bx" value
by 4096 (note, in reality it can be i*4096, i=1,2,...). Partially corrected in
this way value is stored in "bx.cor". The correction is done per chip, so
"bx.cor" is supposed to be used within one chip only. Across the chips, one
needs to use "bx" (without corrections), ie. "bx.cor" modulo 4096. This is
because in an event, depending on previous triggers, the correction may be
done in one chip, but not in the other.

"adc" is raw ADC value, "a" is pedestal subtracted (if pedestal position could
not be determined, a=NA meaning "not available"). Pedestal is determined per
SCA. "trig" is the trigger flag, "n.trig" - number of triggers in the event
(across all chips).  To find retriggers with almost consecutive bunch crossing
(BX), the latter are clustered in groups with the gaps inside <= 2 (ie. BX,
BX+3, BX+6, BX+9).  Here, all "bx" from 4 chips are merged together to define
BX.  The cluster size is denoted as "nbx", retrigger means nbx>1. Inside the
cluster, the bunch crossings are numbered by ibx=1..nbx.  The cluster
consecutive number is denoted by bx.group. For 4 chips, bx.group <= 15*4, it
can be 15*4 only if all nbx=1 and all 15 SCA in 4 chips are filled without
overlaps.

Similar variables with ".chip" suffix (n.trig.chip, nbx.chip, ibx.chip,
bx.group.chip) have the same meanings, but all chips are considered
independently and their "bx" are not merged. So, if BX in chip 0 is followed
by BX+1 in chip 1, this is not considered as a retrigger.  By definition,
n.trig.chip<=64 (while n.trig<=4*64), bx.group.chip<=15, nbx<=15.

Finally, variables with ".trig" suffix (nbx.trig, ibx.trig, bx.group.trig) are
analogous to ".chip" variables, but they are determined when all events
without triggers are removed. It is known, that both BX+1 and BX are readout
when signal arrives at BX clock edge.  In quiet conditions there should be no
triggers in BX+1. ".trig" variables are determined when this known effect is
removed and only problematic retriggers due to noises (with real triggers in
BX+1) are left.

"x","y" - position of the connected pixel.
Orientation of x,y axes, if looked at FEV8 with chips on the top:
      y      DIF
      | Chip3 Chip0
      | Chip2 Chip1
      |______________x
Channels with 2 or 4 connected pixels are on the top and on the right in this
picture.

Cntrl-D or quit() in the R console terminates the program.

The R interactive console responds to user input and has its own "event loop".
GTK+ GUI also has the "event loop" which is executed when R event loop is
idle.  If, however, R session is blocked (eg. when the user is reading a help
page, for example, has typed ?qplot to get a help on a "qplot" command), GTK+
GUI is hanged, this is completely normal.

If you have any questions or comments, please, do not hesitate to contact the
author, Vladislav Balagura (balagura@llr.in2p3.fr). I am especially interested
to hear which plots you create yourself or you may need. If they may be of
common interest to many people, please, inform me and I shall include them in
the next versions.
And, of course, Use R!

')
