`%notin%` <- Negate(`%in%`)

#' epidecodeR object - a S4 class object
#'
#' @slot t data.frame.
#' @slot e data.frame.
#' @slot eventcounts numeric.
#' @slot grptables list.
#' @slot grpcounts integer.
#' @slot sign.test data.frame.
#'
#' @return
#' @export
#'
#' @examples
setClass("epidecodeR", slots=list(t="data.frame", e="data.frame", eventcounts="numeric", grptables="list", grpcounts="integer",sign.test="data.frame"))

#' Analysis function for generating epidecodeR object. This function distributes dysregulated genes into user defined groups and calculates cumulative probabilities and ANOVA test statistics for significance testing in difference of log2FC means between groups of dysregulated genes
#'
#' @param events char - Name of events file. This can be a txt file with two columns: 1) id & 2) counts of events in the gene. Optionally, users can provide a 3+ column .bed file. The count of events per gene in fourth column are calculated to determine degree of events per gene; Default NULL
#' @param deg char - Name of dysregulated genes file. This file is a three column file consisting of column 1: id (Make sure ID type matches between events and deg); column 2) log2foldchange; 3) P value of signficance of fold change; Default NULL
#' @param gtf_file char - Name of compressed gtf file. Use gtf file if .bed file used as events input and users wish to count events per gene from bed file by comparing coordinates in bed to gene coordinates in gtf to assign events to genes. Note: For coordinates overlapping to multiple features in gtf, only one feature is assigned to the coordinate, which is choosen arbitrarily; Default NULL
#' @param id_type char - Name of id type used to count events per gene. ID type must match between events and DEG file. For example, if 'gene_name' is used as ID type in DEG file, same ID type must be used to assign coordinates to genes. In case the DEG list contains two ID types merged e.g. 'ENSMUSG00000035299.16|Mid1' users can give merge as parameter for id_type; Default gene_name
#' @param boundaries numeric - Number of base pairs to include within boundries of genes for event counting. This option adds # of bases to start and end of coordinates of the genes to include promotor regions within gene for overlap with bed files and event counting; Default 0
#' @param pval numeric - P value cut-off of dysregulated genes in DEG file to be considered for distribution into groups. Default: 0.05
#' @param param numeric - Defines the number and size of groups of dysregulated genes. Allowed values are param = 1 [0 events: 1+ events]; param = 2 [0 events: 1-N event: (N+1)+ event]; param = 3 [0 events; 1 event; 2-N events; (N+1)+ events]; N is user defined limit of the group provided using ints parameter
#' @param ints vector - A vector of intervals defining limits of the degree of group for param = 2 and param = 3. e.g. c(1, 4) or c(2, 5): For param = 2, Default :c(1,4) and for param = 3, Default: c(2,5)
#' @import dplyr
#' @importFrom methods new
#' @importFrom methods setClass
#' @importFrom stats aov
#' @importFrom stats na.omit
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom utils read.table
#' @return An epidecodeR object containing tables of theoretical and empirical cumulative probabilities of the log2FC (quantiles), tables of genes distributed into user defined groups, counts of genes per user defined groups, table of one-way ANOVA significance testing for difference in mean log2FC of groups of genes
#' @export
#' @exportClass epidecodeR
#' @examples
#' events<-system.file("extdata", "con_peak.bed", package="epidecodeR")
#' deg<-system.file("extdata", "deg.txt", package="epidecodeR")
#' epiobj <- epidecodeR(events = events, deg = deg, pval=0.05, param = 3, ints=c(2,5))
epidecodeR<-function (events, deg, gtf_file, id_type, boundaries, pval, param, ints) {

  ## Check inputs and throw errors
  if (missing(events) || is.null(events)) {
    message ("Error: events file is mandatory. Please check users manual for accepted formats!")
    ## Throw error and exit
  }
  if (missing(deg) || is.null(deg)) {
    message ("Error: DEG file is mandatory. Please check users manual for accepted formats!")
    ## Throw error and exit
  }
  if (missing(gtf_file)) {
    gtf_file<-NULL
  }
  if (missing(id_type) || is.null(id_type)) {
    id_type="gene_name"
  }
  if (missing(boundaries) || is.null(boundaries)) {
    boundaries = 0
  }
  if (missing(pval) || is.null(pval)) {
    pval=0.05
  }
  if (missing(param) || is.null(param)) {
    param = 2
    ints = c(1, 1)
  } else if (param == 1 | param == 2 | param == 3) {
    if (missing(ints) || is.null(ints) || length(ints) != 2) {
      if (param==2) {
        ints = c(1, 4)
      } else if (param == 3) {
        ints = c(2, 4)
      }
    }
  } else {
    message ("Invalid group selected! Allowed values are
             param = 1;
             param = 2;
             param = 3")
    ## Throw error
  }

  ## Import events file
  file<-read.table(events, header = FALSE, row.names = NULL, stringsAsFactors = F, sep = "\t", fill = TRUE)
  if (length(colnames(file))==2) { ## The file type is list of genes and event counts
    if (! is.numeric(file[1,2]) | is.null(file[1,2]) | is.na(file[1,2])) { ## Check if header is present and if present remove first row
      file<-file[-1,]
    }
    eventcounts<-as.numeric(file[,2]) ## store counts in the eventcounts variable
    names(eventcounts)<-file[,1] ## store gene names for counts in the eventcounts variable
  } else if ((length(colnames(file))>=3)==T) {
    if (is.null(gtf_file)) { ## Input is bed and GTF not required; continue with event count per gene on column 4
      eventcounts<-unlist(lapply (unique (file[,4]), function (x) length(which (file[,4]==x))))
      names(eventcounts)<-unique(file[,4])
    } else {
      bed<-with(file, GenomicRanges::GRanges(seqnames = file$V1, IRanges(start = file$V2, end = file$V3)))
      gtf<-rtracklayer::import(gtf_file)
      #gtf<-gtf[gtf$type=="gene",]
      #rownames(gtf)<-c(1:length(gtf[,1]))
      gtf<-gtf[,c("type", "gene_id", "gene_name")]
      gtf<-gtf[GenomicRanges::elementMetadata(gtf)[,"type"]=="gene"]
      #gtf<-makeGRangesFromDataFrame(gtf, keep.extra.columns=T)
      gtf<-gtf+boundaries
      overlap<-suppressWarnings(GenomicRanges::findOverlaps(bed, gtf, select = "arbitrary"))
      names(overlap)<-c(1:length(overlap))
      overlap<-na.omit(overlap)
      q<-as.data.frame(bed[as.numeric(names(overlap))])
      s<-as.data.frame(gtf[overlap])
      if (id_type == "merge") {
        df<-unique(data.frame(chr=q$seqnames, start=q$start, end=q$end, id=paste(s[,"gene_id"],s[,"gene_name"], sep="|"), width=q$end-q$start, strand=s$strand))
      } else {
        df<-unique(data.frame(chr=q$seqnames, start=q$start, end=q$end, id=s[,id_type], width=q$end-q$start, strand=s$strand))
      }

      eventcounts<-unlist(lapply (unique (df$id), function (x) length(which (df$id==x))))
      names(eventcounts)<-unique(df$id)
    }
  } else {
    message ("Error in input file! Input should be tab separated gene list and event counts or BED file of events. Exiting...")
  }

  ## Import DEG list
  del<-read.table(deg, header = T, row.names = NULL, stringsAsFactors = F, sep = "\t", fill = T)
  #del<-del[,1:3]
  #colnames(del)<-c("ID", "Log2FC", "Pval")
  if (! is.numeric(del[,2]) && ! is.numeric(del[,3])) {
    message ("DEG list not in proper format!
             Input should be three columns separated by tab with column 1 as id, column 2 as log2FC and column 3 as P value")
  }
  del<-na.omit(del)
  del<-del[del[,3]<=pval,]

  grp1<-del[del[,1] %notin% names(eventcounts),]
  cdfdf<-data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=mean(grp1[,2]), sd=sd(grp1[,2])), plot.it = F), grp="0")
  ecdfdf<-data.frame(EnvStats::ecdfPlot(grp1[,2], plot.it = F), grp="0")

  grptables<-list("0"=grp1)

  if (param == 1) {
    grp2<-del[del[,1] %in% names(eventcounts[eventcounts>0]),]
    if (dim(grp2)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=mean(grp2[,2]), sd=sd(grp2[,2])), plot.it = F), grp="1+"))
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp2[,2], plot.it = F), grp="1+"))
    } else if (dim(grp2)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp2[,2], plot.it = F), grp="1+"))
    } else {
      message("No genes map to \"1+\" group!")
    }
    grptables[["1+"]]<-grp2
    grpcounts=c("0"=length(grp1[,1]), "1+"=length(grp2[,1]))
  }


  if (param == 2) {
    grp2name=paste0("1to", ints[2], sep = "")
    grp2<-del[del[,1] %in% names(eventcounts[eventcounts>=1 & eventcounts<=ints[2]]),]
    if (dim (grp2)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=mean(grp2[,2]), sd=sd(grp2[,2])), plot.it = F), grp=grp2name))
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else if (dim(grp2)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else {
      message(paste ("No genes map to \"",grp2name,"\" group", sep = ""))
    }

    grp3name=paste0(ints[2]+1, "+", sep = "")
    grp3<-del[del[,1] %in% names(eventcounts[eventcounts>ints[2]]),]
    if (dim (grp3)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=mean(grp3[,2]), sd=sd(grp3[,2])), plot.it = F), grp=grp3name))
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else if (dim(grp3)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else {
      message(paste ("No genes map to \"",grp3name,"\" group", sep=""))
    }
    #grptables[["1"]]<-grp2
    grptables[[grp2name]]<-grp2
    #grptables[["2+"]]<-grp3
    grptables[[grp3name]]<-grp3
    #grpcounts=c("0"=length(grp1[,1]), "1"=length(grp2[,1]), "2+"=length(grp3[,1]))
    grpcounts=c(length(grp1[,1]), length(grp2[,1]), length(grp3[,1]))
    names(grpcounts)<-names(grptables)
  }


  if (param == 3) {
    grp2name="1"
    grp2<-del[del[,1] %in% names(eventcounts[eventcounts==1]),]
    if (dim (grp2)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=mean(grp2[,2]), sd=sd(grp2[,2])), plot.it = F), grp=grp2name))
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else if (dim (grp2)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else {
      message("No genes map to \"1\" group")
    }

    grp3name=paste0("2to", ints[2], sep="")
    grp3<-del[del[,1] %in% names(eventcounts[eventcounts>1 & eventcounts<=ints[2]]),]
    if (dim (grp3)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=mean(grp3[,2]), sd=sd(grp3[,2])), plot.it = F), grp=grp3name))
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else if (dim (grp3)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else {
      message (paste("No genes map to \"", grp3name, "\" group"))
    }

    grp4name=paste0(ints[2]+1, "+", sep="")
    grp4<-del[del[,1] %in% names(eventcounts[eventcounts>ints[2]]),]
    if (dim (grp4)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=mean(grp4[,2]), sd=sd(grp4[,2])), plot.it = F), grp=grp4name))
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp4[,2], plot.it = F), grp=grp4name))
    } else if (dim(grp4)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(EnvStats::ecdfPlot(grp4[,2], plot.it = F), grp=grp4name))
    } else {
      message(paste("No genes map to \"", grp4name, "\" group"))
    }
    grptables[["1"]]<-grp2
    grptables[[grp3name]]<-grp3
    grptables[[grp4name]]<-grp4
    grpcounts=c(length(grp1[,1]), length(grp2[,1]), length(grp3[,1]), length(grp4[,1]))
    names(grpcounts)<-names(grptables)
  }

  test<-data.frame(aov(Order.Statistics~grp, data = ecdfdf) %>% rstatix::tukey_hsd())

  #epidecodeR<-setClass("epidecodeR", slots=list(t="data.frame", e="data.frame", eventcounts="numeric", grptables="list", grpcounts="integer",sign.test="data.frame"))
  ro<-new("epidecodeR", t=cdfdf, e=ecdfdf, eventcounts=eventcounts, grptables=grptables, grpcounts=grpcounts, sign.test=test)
  return(ro)
}


#' Title
#'
#' @param objdf Data frame of the object from epidecodeR object
#' @param title Title of the plot
#' @param lim Limits of x-axis
#' @param xlab X-axis name
#' @param ylab Y-axis name
#' @import ggplot2
#' @return returns a generic ggplot to be filled by makeplot function
#'
#' @examples
#' \dontrun{
#' plottingfunc(objdf=obj@t, title=title, lim=lim, xlab=xlab, ylab=ylab)
#' }
plottingfunc<-function (objdf, title, lim, xlab, ylab) {
  p<-ggplot2::ggplot(objdf) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_vline(xintercept = 0, color="grey", alpha=0.5, linetype="dashed") +
    ggplot2::geom_hline(yintercept = 0.5, color="grey", alpha=0.5, linetype="dashed") +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::xlim(lim) +
    ggplot2::ylim(0,1) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(p)
}

#' Generate CDF plot using epidecodeR object generated using epidecodeR function
#'
#' @param obj epidecodeR object - epidecodeR object generated using epidecodeR function
#' @param type char - Type of CDF plot to generate; Accepted values 't': theoretical CDF plot; 'e': empirical CDF plot; 'both': Creates both theoretical and empirical plots.  Default: both
#' @param lim vector - Upper and lower limits of log2FC for X-axis
#' @param title char - Title of the plot
#' @param xlab char - X-axis label
#' @param ylab char - Y-axis label
#' @import ggplot2
#' @return A CDF plot
#' @export
#'
#' @examples events<-system.file("extdata", "con_peak.bed", package="epidecodeR")
#' deg<-system.file("extdata", "deg.txt", package="epidecodeR")
#' epiobj <- epidecodeR(events = events, deg = deg, pval=0.05, param = 3, ints=c(2,5))
#' p<-makeplot(epiobj, lim = c(-2,2), title = "H3K27ac", xlab = "siSDS3/WT")
makeplot<-function (obj, type, lim, title, xlab, ylab) {
  if (missing(obj)) {
    message ("No epidecodeR object provided!")
  }
  if (missing(type)) {
    type="both"
  }
  if (missing(lim)) {
    lim=c(-1,1)
  }
  if (missing(title)) {
    title=""
  }
  if (missing(xlab)) {
    xlab="Log2FC"
  }
  if (missing(ylab)) {
    ylab="Cumulative Probabilities"
  }
  pt<-plottingfunc(objdf=obj@t, title=title, lim=lim, xlab=xlab, ylab=ylab)
  pe<-plottingfunc(objdf=obj@e, title=title, lim=lim, xlab=xlab, ylab=ylab)
  pb<-plottingfunc(objdf = NULL, title=title, lim=lim, xlab=xlab, ylab=ylab)

  for (i in 1:length(obj@grpcounts)) {
    if (obj@grpcounts[i]==0) {
      pt<-pt+ggplot2::geom_blank(data=data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=0, sd=1), plot.it = F), grp=names(obj@grpcounts[i])), aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
      pe<-pe+ggplot2::geom_blank(data=data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=0, sd=1), plot.it = F), grp=names(obj@grpcounts[i])), aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
      pb<-pb+ggplot2::geom_blank(data=data.frame(EnvStats::cdfPlot(distribution = "norm", param.list = list(mean=0, sd=1), plot.it = F), grp=names(obj@grpcounts[i])), aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
    } else if (obj@grpcounts[i]==1) {
      message(paste("Theoretical plot for ", names(obj@grpcounts[i]), " group not included since only one gene found to belong to the group!", sep = ""))
    }
  }
  pt<-pt+ggplot2::geom_line(aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
  pt<-pt+ggplot2::scale_color_discrete(name="Groups (genes)", breaks=names(obj@grpcounts), labels=paste(names(obj@grpcounts)," (",obj@grpcounts,")", sep = "" ), drop=FALSE)
  pe<-pe+ggplot2::geom_point(aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp), size=1) +
         ggplot2::geom_line(aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp))
  pe<-pe+ggplot2::scale_color_discrete(name="Groups (genes)", breaks=names(obj@grpcounts), labels=paste(names(obj@grpcounts)," (",obj@grpcounts,")", sep = "" ), drop=FALSE)
  pb<-pb+ggplot2::geom_line(data=obj@t, aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
  pb<-pb+ggplot2::geom_point(data=obj@e, aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp),size=1) +
         ggplot2::geom_line(data=obj@e, aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp), linetype="dashed")
  pb<-pb+ggplot2::scale_color_discrete(name="Groups (genes)", breaks=names(obj@grpcounts), labels=paste(names(obj@grpcounts)," (",obj@grpcounts,")", sep = "" ), drop=FALSE)
  if (type == "t") {
    p<-pt
  } else if (type == "e") {
    p<-pe
  } else if (type == "both") {
    p<-pb
  }
  return(p)
}


#' Generates boxplot of distribution of log2FC of dysregulated genes and adjusted P value of signficance test of difference in mean log2FC between groups computed using one-way ANOVA test
#'
#' @param obj epidecodeR object - epidecodeR object generated using epidecodeR function
#' @param title char - Title of the plot
#' @param ylab char - Y-axis label
#' @import ggplot2
#'
#' @return Boxplot of distribution of log2FC of dysregulated genes between groups
#' @export
#'
#' @examples events<-system.file("extdata", "con_peak.bed", package="epidecodeR")
#' deg<-system.file("extdata", "deg.txt", package="epidecodeR")
#' epiobj <- epidecodeR(events = events, deg = deg, pval=0.05, param = 3, ints=c(2,5))
#' pt<-plot_test(epiobj, title = "H3K27ac", ylab = "siSDS3/WT")
plot_test<-function (obj, title, ylab) {

    if (missing(title) || is.null(title)) {
      title = ""
    }
    if (missing(ylab) || is.null(ylab)) {
      ylab="Log2FoldChange"
    }
    bplot<-ggplot2::ggplot()+
    ggplot2::ggtitle(title) +
    ggplot2::xlab("Groups") +
    ggplot2::ylab(ylab) +
    ggplot2::theme_classic()

  over<-NULL

  for (i in 1:length(obj@grpcounts)) {
    if (obj@grpcounts[i]==0) {
      bplot<-bplot+ggplot2::geom_blank(data=data.frame(x=rnorm(100), grp=names(obj@grpcounts[i])), aes(grp, x, group=grp, color=grp))
    } else {
      bplot<-bplot+ggplot2::geom_boxplot(data=(subset(obj@e, obj@e$grp %in% names(obj@grpcounts[i]))), aes(grp, Order.Statistics, group=grp, color=grp), width=0.4)
    }
  }

  #for (i in 1:dim(obj@sign.test)[1]) {
  #  over[i]=max(obj@e$Order.Statistics)+2*i
  #}

  if (length (unique(obj@sign.test$p.adj.signif))==1 && unique(obj@sign.test$p.adj.signif)=="ns") {
    bplot<-bplot +
      ggplot2::scale_color_discrete(name="Groups (genes)", breaks=names(obj@grpcounts), labels=paste(names(obj@grpcounts)," (",obj@grpcounts,")", sep = "" ), drop=FALSE) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text=ggplot2::element_text(size=12), axis.title=ggplot2::element_text(size=14,face="bold"))
  } else {
    bplot<-bplot +
      ggplot2::scale_color_discrete(name="Groups (genes)", breaks=names(obj@grpcounts), labels=paste(names(obj@grpcounts)," (",obj@grpcounts,")", sep = "" ), drop=FALSE) +
      ggpubr::stat_pvalue_manual(obj@sign.test, hide.ns = T, label = "p.adj", y.position = max(obj@e$Order.Statistics)*1.2, step.increase = 0.2) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text=ggplot2::element_text(size=12), axis.title=ggplot2::element_text(size=14,face="bold"))
  }


  return(bplot)
}

