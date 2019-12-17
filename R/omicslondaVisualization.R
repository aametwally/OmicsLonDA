#' Visualize Longitudinal Feature
#'
#' Visualize Longitudinal Feature
#'
#' @param se_object SummarizedExperiment object contains omics count/level matrix 
#' and  metadata 
#' @param text feature name
#' @param unit time interval unit
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures
#' (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @return null
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @importFrom methods is
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' library(SummarizedExperiment)
#' data("omicslonda_data_example")
#' omicslonda_se_object_adjusted = adjustBaseline(
#'                  se_object = omicslonda_data_example$omicslonda_se_object)
#' omicslonda_test_object = omicslonda_se_object_adjusted[1,]
#' visualizeFeature(se_object = omicslonda_test_object, text = "Feature_1",
#'                  unit = "days", ylabel = "Normalized Count", 
#'                  col = c("blue", "firebrick"), prefix = tempfile())
#' @export
visualizeFeature = function (se_object = NULL, text = "featureName",
                                unit = "days", ylabel = "Normalized Count", 
                                col = c("blue", "firebrick"), prefix = "Test")
{
    message("Visualizing Feature = ", text)
    
    ### validate se_object
    stopifnot(is(se_object, "SummarizedExperiment"))
    stopifnot(all(c("Subject", "Time", "Group") %in% colnames(colData(se_object))))
    ## validate col
    stopifnot(length(col) == 2)
    
    if (!dir.exists(prefix)){
        dir.create(file.path(prefix))
    }
    
    df = data.frame(colData(se_object))
    df$Count = as.vector(assay(se_object))
    group.levels = sort(unique(as.character(df$Group)))
    
    p = ggplot(df, aes(.data$Time, .data$Count, colour = .data$Group, group = interaction(.data$Group,
                                                         .data$Subject)))
    p = p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1,
        alpha = 0.7) + theme_bw() + ggtitle(paste("Feature = ", text,
        sep = "")) + labs(y = ylabel, x = sprintf("Time (%s)", unit)) +
        scale_colour_manual(values = col,
        labels = c(group.levels[1], group.levels[2])) +
        theme(axis.text.x = element_text(colour="black", size=12, angle=0,
                                        hjust=0.5, vjust=0.5, face="bold"),
            axis.text.y = element_text(colour="black", size=12, angle=0,
                                        hjust=0.5, vjust=0.5, face="bold"),
            axis.title.x = element_text(colour="black", size=15, angle=0,
                                        hjust=.5, vjust=0.5, face="bold"),
            axis.title.y = element_text(colour="black", size=15, angle=90,
                                        hjust=.5, vjust=.5, face="bold"),
            legend.text=element_text(size=15, face="plain"),
            legend.title = element_blank(), 
            plot.title = element_text(hjust = 0.5)) +
        theme(legend.position="top") + scale_x_continuous(breaks = waiver())


    ggsave(filename=paste(prefix, "/", "Feature_", text, ".jpg", sep=""),
            dpi = 1200, height = 10, width = 15, units = 'cm')
}






#' Visualize the feature trajectory with the fitted Splines
#'
#' Plot the longitudinal features along with the fitted splines
#'
#' @param se_object SummarizedExperiment object contains omics count/level matrix 
#' and  metadata 
#' @param omicslonda_object The returned object from omicslonda analysis 
#' @param fit.method The fitting method (ssgaussian)
#' @param text feature name
#' @param unit time unit used in Time vector (hours, days, weeks, months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures
#' (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @return null
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @importFrom methods is
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' library(SummarizedExperiment)
#' data(omicslonda_data_example)
#' omicslonda_se_object_adjusted = adjustBaseline(
#'                  se_object = omicslonda_data_example$omicslonda_se_object)
#' omicslonda_test_object = omicslonda_se_object_adjusted[1,]
#' points = seq(1, 500, length.out = 500)
#' res = omicslonda(se_object = omicslonda_test_object, n.perm = 10,
#'                  fit.method = "ssgaussian", points = points, text = "Feature_1",
#'                  parall = FALSE, pvalue.threshold = 0.05, 
#'                  adjust.method = "BH", time.unit = "days",
#'                  ylabel = "Normalized Count",
#'                  col = c("blue", "firebrick"), prefix = "OmicsLonDA_example")
#' visualizeFeatureSpline(se_object = omicslonda_test_object, omicslonda_object = 
#'                  omicslonda_data_example$omicslonda_results, fit.method = "ssgaussian",
#'                  text = "Feature_1", unit = "days",
#'                  ylabel = "Normalized Count", 
#'                  col = c("blue", "firebrick"),
#'                  prefix = tempfile())
#' @export
visualizeFeatureSpline = function (se_object = NULL, omicslonda_object = NULL, fit.method = "ssgaussian",
                                    text = "FeatureName", unit = "days",
                                    ylabel = "Normalized Count", 
                                    col = c("blue", "firebrick"),
                                    prefix = "Test")
{ 
    message("Visualizing Fitted Smoothings Splines for each Group of Feature = ", text)
    
    ### validate se_object
    stopifnot(is(se_object, "SummarizedExperiment"))
    stopifnot(all(c("Subject", "Time", "Group") %in% colnames(colData(se_object))))
    ## validate col
    stopifnot(length(col) == 2)
    ## validate fit.method
    stopifnot(fit.method %in% c("ssgaussian"))
    
    
    if (!dir.exists(prefix)){
        dir.create(file.path(prefix))
    }
    
    model = omicslonda_object$model
    df = data.frame(colData(se_object))
    df$Count = as.vector(assay(se_object))
    group.levels = sort(unique(as.character(df$Group)))
    
    dd.null = model$dd.null
    dd.0 = model$dd.0
    dd.1 = model$dd.1
    
    ln = factor(c(rep("longdash", nrow(dd.0)), rep("longdash", nrow(dd.1))))
    dm = rbind(dd.0[,c("Time", "Count", "Group", "Subject")], dd.1[,c("Time",
                                                "Count", "Group", "Subject")])
    dm$lnn=ln
    df_subset = df[,c("Time", "Count", "Group", "Subject")]
    
    
    ## Hack for matching colors to each group
    x = data.frame(group = c(group.levels, "fit.0", "fit.1"), 
               label = c(group.levels[1], group.levels[2],
                          paste(group.levels[1], ".fit", sep=""),
                         paste(group.levels[2], ".fit", sep="")),
               color = c(col, col))
    x.ordered <- x[order(x$color),]
    x.ordered$group = as.character(x.ordered$group)
    x.ordered$label = as.character(x.ordered$label)
    x.ordered$color = as.character(x.ordered$color)
    
    
    p = ggplot()
    p = p + theme_bw()  + 
        
        geom_point(data = df_subset, aes(.data$Time,
                                        .data$Count, 
                                        colour = .data$Group, group =
                                        interaction(.data$Group,
                                                    .data$Subject)),
                size=1, alpha=0.1) +
        geom_line(data = df_subset, aes(.data$Time,
                                        .data$Count, 
                                        colour = .data$Group, group =
                                        interaction(.data$Group,
                                                    .data$Subject)),
                size=1, alpha=0.1) +
        geom_line(data= dm, aes(.data$Time, .data$Count, 
                                colour = .data$Group,
                                group = interaction(.data$Group,
                                                    .data$Subject),
                                linetype=.data$lnn), size=2, alpha=0.8) +
        ggtitle(paste("Feature = ", text, sep = "")) + labs(y = ylabel,
                                                x = sprintf("Time (%s)", unit)) +
        scale_colour_manual(breaks = x.ordered$group,
                            labels = x.ordered$label,
                            values = x.ordered$color) +  
        theme(axis.text.x=element_text(colour="black", size=12, angle=0,
                                        hjust=0.5, vjust=0.5, face="bold"),
            axis.text.y = element_text(colour="black", size=12, angle=0,
                                        hjust=0.5, vjust=0.5, face="bold"),
            axis.title.x=element_text(colour="black", size=15, angle=0,
                                    hjust=.5, vjust=0.5, face="bold"),
            axis.title.y=element_text(colour="black", size=15, angle=90,
                                    hjust=.5,
                                    vjust=.5, face="bold"), 
            legend.text=element_text(size=15, face="plain"),
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
        theme(legend.position="top") + scale_x_continuous(breaks = waiver()) +
        guides(linetype=FALSE, size =FALSE)
    
    ggsave(filename=paste(prefix, "/", "Feature_", text, "_CurveFitting_",
            fit.method, ".jpg", sep=""), dpi = 1200, height = 10, width = 15,
    units = 'cm')

}





#' Visualize test statistics empirical distribution
#'
#' Visualize sest statistics empirical distribution
#'
#' @param omicslonda_object The returned object from omicslonda analysis 
#' @param text Feature name
#' @param fit.method fitting method
#' @param prefix prefix to be used to create directory for the analysis results
#' @return null
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' library(SummarizedExperiment)
#' data(omicslonda_data_example)
#' visualizeTestStatHistogram(omicslonda_object = omicslonda_data_example$omicslonda_results, 
#'                  text = "Feature_1", 
#'                  fit.method = "ssgaussian", prefix = tempfile())
#' @export
visualizeTestStatHistogram = function(omicslonda_object = NULL, text = "FeatureName", 
                                      fit.method = "ssgaussian", prefix = "Test"){
    message("Visualizing test statstic null distribution of Feature = ", text)
    
    ## validate fit.method
    stopifnot(fit.method %in% c("ssgaussian"))
    
    if (!dir.exists(prefix)){
        dir.create(file.path(prefix))
    }
    
    null_dist = omicslonda_object$distribution
    
    jpeg(filename = paste(prefix, "/", "Feature_", text, "_testStat_distribution_ALL_only", fit.method, ".jpg", sep = ""), 
         res = 1200, height = 7, width = 10, units = 'cm')
    
    hist(null_dist, xlab = "testStat", 
         breaks = 100, col = "gray", border = "gray", 
         main = paste("test statistic null dist of ", text, sep=""), 
         freq = TRUE)
    dev.off()
}





#' Visualize significant time interval
#'
#' Visualize significant time interval
#'
#' @param omicslonda_object The returned object from omicslonda analysis
#' @param fit.method Fitting method (ssgaussian)
#' @param text Feature name
#' @param unit time unit used in the Time vector
#' (hours, days, weeks, months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures
#' (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @return null
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples
#' library(SummarizedExperiment)
#' data(omicslonda_data_example)
#' visualizeArea(omicslonda_object = omicslonda_data_example$omicslonda_results, 
#'                  fit.method = "ssgaussian",
#'                  text = "Feature_1", unit = "days", 
#'                  ylabel = "Normalized Count", col =
#'                  c("blue", "firebrick"), prefix = tempfile())
#' @export
visualizeArea = function(omicslonda_object = NULL, fit.method = "ssgaussian",
                            text = "FeatureName", unit = "days",
                            ylabel = "Normalized Count", col =
                            c("blue", "firebrick"), prefix = "Test")
{
    message("Visualizing Significant Intervals of Feature = ", text)
    
    ## validate col
    stopifnot(length(col) == 2)
    ## validate fit.method
    stopifnot(fit.method %in% c("ssgaussian"))

    if (!dir.exists(prefix)){
        dir.create(file.path(prefix))
    }

    model.ss = omicslonda_object$model
    start = omicslonda_object$start
    end = omicslonda_object$end
    df = omicslonda_object$df
    group.levels = sort(unique(as.character(df$Group)))

    Time = 0 ## This line is to pass the CRAN checks for the aes in ggplot2
    sub.11 = list()
    sub.10 = list()
    xx = NULL
    for(i in seq_len(length(start)))
    {
        sub.11[[i]] = subset(model.ss$dd.1, Time >= start[i] & Time <= end[i])
        sub.10[[i]] = subset(model.ss$dd.0, Time >= start[i] & Time <= end[i])
        cmd=sprintf('geom_ribbon(data=sub.10[[%d]], aes(
                ymin = sub.11[[%d]]$Count, ymax = Count), colour= "grey3",
                fill="grey69", alpha = 0.6)', i, i)
        if (i != 1)
        {
        xx = paste(xx, cmd, sep = "+")
        } else
        {
        xx = cmd
        }
    }

    dd.0 = model.ss$dd.0
    dd.1 = model.ss$dd.1

    dm = rbind(dd.0, dd.1)
    p1 = 'ggplot(dm, aes(Time, Count, colour = Group, group =
                interaction(Group, Subject))) +
    theme_bw() + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1,
    alpha = 0.5) + ggtitle(paste("Feature = ", text, sep = "")) +
                                labs(y = ylabel,
    x = sprintf("Time (%s)", unit)) + scale_colour_manual(values = col,
    breaks = c("fit.0", "fit.1"),
    labels = c(paste(group.levels[1], ".fit", sep = ""),
    paste(group.levels[2], ".fit", sep = ""))) +
    theme(axis.text.x = element_text(colour = "black", size = 12, angle = 0,
    hjust = 0.5, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, angle = 0,
    hjust = 0.5, vjust = 0.5, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 15, angle = 0,
    hjust = 0.5, vjust = 0.5, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 15, angle = 90,
    hjust = 0.5, vjust = 0.5, face = "bold"),
    legend.text = element_text(size = 15, face="plain"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "top") + scale_x_continuous(breaks = waiver())'
    p2 = xx
    p3 = paste(p1, p2, sep="+")
    p = eval(parse(text = p3))


    ggsave(filename=paste(prefix, "/", "Feature_", text,
                            "_SignificantInterval_", fit.method, ".jpg", sep=""),
                            dpi = 1200, height = 10, width = 15,
                            units = 'cm')
}
