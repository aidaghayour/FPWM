#' A function for generating the forked Position Weight Matrix
#'
#' This function takes the generated class object and plots a forked position weight matrix.
#' @param GraphDataObj is an object of S4 class with modified and converted data ready to be plotted.
#' @param Methylation is a logical value. If it set on TRUE, Methylation level chart will also be plotted. If Flase, only sequence logos will be shown.
#' @export
FPWMPlotter <- function(TheObject, Methylation = TRUE)
{
  
  
  df <- data.frame(TheObject@forked)
  colnames(df) = c('Positions', 'A', "C", "G", "T")
  
  positions <- df[['Positions']]
  ###################################################### Grids
  W= 0.25
  H = 0.15
  N = length(TheObject@id)
  
  
  #################################
  H = (1-H)/N
  W=H*1.2
  
  grid::grid.rect(
    gp = grid::gpar(col="grey"),
    name = "Parent",
    x = grid::unit(.2, "npc"),
    y = grid::unit(.5, "npc"),
    width = grid::unit(W, "npc"),
    height = grid::unit(H, "npc"))
  vpparent <-
    grid::viewport(
      x = grid::unit(.2, "npc"),
      y = grid::unit(.5, "npc"),
      width = grid::unit(W, "npc"),
      grid::unit(H, "npc"))
  
  
  if (N %% 2 == 0){
    grid::grid.rect(
      gp = grid::gpar(col="grey"),
      name = paste0("c",floor(N/2)),
      x = grid::grobX("Parent", 0)+ grid::unit(W/1.3, "npc") ,
      y = grid::grobY("Parent", 0)+ grid::unit(H/1.5, "npc"),
      width = grid::unit(W, "npc"),
      height = grid::unit(H, "npc"))
    
    assign(paste0("vp",floor(N/2)),
           grid::viewport(
             x = grid::grobX("Parent", 0)+ grid::unit(W/1.3, "npc") ,
             y = grid::grobY("Parent", 0)+ grid::unit(H/1.5, "npc"),
             width = grid::unit(W, "npc"),
             height = grid::unit(H, "npc")))###
    
    grid::grid.rect(
      gp = grid::gpar(col="grey"),
      name = paste0("c",(floor(N/2)+1)),
      x = grid::grobX("Parent", 0)+ grid::unit(W/1.3, "npc") ,
      y = grid::grobY("Parent", 0)- grid::unit(H/1.5, "npc"),
      width = grid::unit(W, "npc"),
      height = grid::unit(H, "npc"))
    
    assign(paste0("vp",(floor(N/2)+1)),
           grid::viewport(
             x = grid::grobX("Parent", 0)+ grid::unit(W/1.3, "npc") ,
             y = grid::grobY("Parent", 0)- grid::unit(H/1.5, "npc"),
             width = grid::unit(W, "npc"),
             height = grid::unit(H, "npc")))##
    if ( N != 2 ){
      for (i in c((floor(N/2)+1):(N-1))){
        grid::grid.rect(
          gp = grid::gpar(col="grey"),
          name = paste0("c",(i+1)),
          x = grid::grobX(paste0("c",i), 90) ,
          y = grid::grobY(paste0("c",i), 270)- grid::unit(H/1.5, "npc"),
          width = grid::unit(W, "npc"),
          height = grid::unit(H, "npc"))
        
        assign(paste0("vp",(i+1)),
               grid::viewport(
                 x = grid::grobX(paste0("c",i), 90) ,
                 y = grid::grobY(paste0("c",i), 270)- grid::unit(H/1.5, "npc"),
                 width = grid::unit(W, "npc"),
                 height = grid::unit(H, "npc")))
      }
      
      for (i in c(floor(N/2):2)){
        grid::grid.rect(
          gp = grid::gpar(col="grey"),
          name = paste0("c",(i-1)),
          x = grid::grobX(paste0("c",i), 90) ,
          y = grid::grobY(paste0("c",i), 90)+ grid::unit(H/1.5, "npc"),
          width = grid::unit(W, "npc"),
          height = grid::unit(H, "npc"))
        
        assign(paste0("vp",(i-1)),
               grid::viewport(
                 x = grid::grobX(paste0("c",i), 90) ,
                 y = grid::grobY(paste0("c",i), 90)+ grid::unit(H/1.5, "npc"),
                 width = grid::unit(W, "npc"),
                 height = grid::unit(H, "npc")))
      }
    }
  }
  if (N %% 2 == 1)
  {
    grid::grid.rect(
      gp = grid::gpar(col="grey"),
      name = paste0("c",(floor(N/2)+1)),
      x = grid::grobX("Parent", 0)+ grid::unit(W/1.3, "npc") ,
      y = grid::grobY("Parent", 0),
      width = grid::unit(W, "npc"),
      height = grid::unit(H, "npc"))
    
    assign(paste0("vp",(floor(N/2)+1)),
           grid::viewport(
             x = grid::grobX("Parent", 0)+ grid::unit(W/1.3, "npc") ,
             y = grid::grobY("Parent", 0),
             width = grid::unit(W, "npc"),
             height = grid::unit(H, "npc")))
    ###
    
    for (i in c((floor(N/2)+1):(N-1))){
      grid::grid.rect(
        gp = grid::gpar(col="grey"),
        name = paste0("c",(i+1)),
        x = grid::grobX(paste0("c",i), 90) ,
        y = grid::grobY(paste0("c",i), 270)- grid::unit(H/1.5, "npc"),
        width = grid::unit(W, "npc"),
        height = grid::unit(H, "npc"))
      
      assign(paste0("vp",(i+1)),
             grid::viewport(
               x = grid::grobX(paste0("c",i), 90) ,
               y = grid::grobY(paste0("c",i), 270)- grid::unit(H/1.5, "npc"),
               width = grid::unit(W, "npc"),
               height = grid::unit(H, "npc")))
    }
    
    for (i in c((floor(N/2)+1):2)){
      grid::grid.rect(
        gp = grid::gpar(col="grey"),
        name = paste0("c",(i-1)),
        x = grid::grobX(paste0("c",i), 90) ,
        y = grid::grobY(paste0("c",i), 90)+ grid::unit(H/1.5, "npc"),
        width = grid::unit(W, "npc"),
        height = grid::unit(H, "npc"))
      
      assign(paste0("vp",(i-1)),
             grid::viewport(
               x = grid::grobX(paste0("c",i), 90) ,
               y = grid::grobY(paste0("c",i), 90)+ grid::unit(H/1.5, "npc"),
               width = grid::unit(W, "npc"),
               height = grid::unit(H, "npc")))
    }
    
  }
  ####################
  for (i in c(1:N)){
    grid::grid.segments(
      grid::grobX("Parent", 0),
      grid::grobY("Parent", 0),
      grid::grobX(paste0("c",i), 180),
      grid::grobY(paste0("c",i), 180),
      gp = grid::gpar(5,col="grey"))}
  ###################
  for (i in c(1:length(TheObject@score))) {
    
    
    grid::grid.text(gp=grid::gpar(fontsize=8, col="grey39"),
                    base::paste(base::toString(
                      base::formatC(TheObject@score[[i]], digits = 1, format = "f")
                    ), "%"),
                    x = grid::grobX(paste0("c",i), 180) - grid::unit(.2, "inches"),
                    y = grid::grobY(paste0("c",i), 0),
                    name = paste0("c",i)
    )}
  
  ##################
  if (Methylation == TRUE)
  {
    Betas = list()
    for (i in c(1:length(TheObject@betalevel))){
      Betas[[i]] <- data.frame(TheObject@betalevel[[i]])
    }
    M1 <- data.frame(TheObject@parentbeta)
    M1 <- M1[order(as.numeric(as.character(M1$X1))),]
    totals1 <-
      as.vector(by(as.numeric(as.character(M1$X2)),
                   as.numeric(as.character(M1$X1)),
                   sum))
    labels1 <-
      unlist(lapply(as.character(totals1), function(x)
        c(rep("", 2), x)))
    posi1 <- rep(totals1 + 1, each = 3)
    M1 <- data.frame(M1, pos = posi1, labels = labels1)
    M1$X1 <- as.character(M1$X1)
    M1$X1 <- factor(M1$X1, levels = unique(M1$X1))
    
    for (i in c(1:length(TheObject@id))) {
      M2 <- data.frame(TheObject@betalevel[[i]])
      M2 <- M2[order(as.numeric(as.character(M2$X1))), ]
      totals2 <-
        as.vector(by(as.numeric(as.character(M2$X2)),
                     as.numeric(as.character(M2$X1)),
                     sum))
      assign(paste0("labels",(i+1)),
        unlist(lapply(as.character(totals2), function(x)
          c(rep("", 2), x))))
      
      assign(paste0("posi",(i+1)),rep(totals2 + 1, each = 3))
      
      M2 <- data.frame(M2, pos = get(paste0("posi",(i+1))), labels = get(paste0("labels",(i+1))))
      
      M2$X1 <- as.character(M2$X1)
      M2$X1 <- factor(M2$X1, levels = unique(M2$X1))
      assign(paste0("N",i), M2)
    }
    
    
    barplot_color <- c("dodgerblue1", "darkorange1", "darkgreen")
    ######################
    
    lim = c()
    for (i in c(2:(length(TheObject@id)+1))){
      lim[i] <- max(as.numeric(as.character(get(paste0("N", (i-1)))$labels)), na.rm = TRUE)
    }
    lim[1] <- max(as.numeric(as.character(M1$labels)), na.rm = TRUE)
    
    MaxLim <- round(max(lim) * 1.6)
    
    
    BSP1 <-
      ggplotify::as.grob(
        ggplot2::ggplot(M1, ggplot2::aes(
          x = X1, y = X2, fill = X3
        ))
        + ggplot2::geom_bar(stat = "identity", position =
                              "stack")
        + ggplot2::geom_text(ggplot2::aes(y = posi1, label =
                                            labels1), vjust = 0, size=2)
        + ggplot2::theme(
          legend.text = ggplot2::element_text(size = 7),
          axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          legend.title = ggplot2::element_blank(),
          legend.background = ggplot2::element_blank(),
          legend.box.background = ggplot2::element_rect(colour = "black"),
          legend.key.size = grid::unit(0.3, "line"),
          legend.position = c(0.2, 1),
          plot.margin = ggplot2::margin(
            t = 10,
            r = 20,
            b = 0,
            l = 19,
            unit = "pt"
          ),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank()
        ) + ggplot2::scale_x_discrete(breaks = NULL) + ggplot2::scale_y_discrete(breaks =
                                                                                   NULL, limits = as.character(0:MaxLim)) + ggplot2::scale_fill_manual(values = barplot_color)
      )
    for (i in c(1:length(TheObject@id))) {
    assign(paste0("BSp",i),
      ggplotify::as.grob(
        ggplot2::ggplot(get(paste0("N", i)), ggplot2::aes(
          x = X1, y = X2, fill = X3
        ))
        + ggplot2::geom_bar(stat = "identity", position =
                              "stack")
        + ggplot2::geom_text(ggplot2::aes(y = get(paste0("posi",(i+1))), label =
                                            get(paste0("labels",(i+1)))), vjust = 0, size=2)
        + ggplot2::theme(
          legend.position = "none",
          axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(
            t = 10,
            r = 20,
            b = 0,
            l = 19,
            unit = "pt"
          ),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank()
        ) + ggplot2::scale_x_discrete(breaks = NULL) + ggplot2::scale_y_discrete(breaks =
                                                                                   NULL, limits = as.character(0:MaxLim)) + ggplot2::scale_fill_manual(values = barplot_color)
      ))}
    
  }
  
  ##################
  sp = TheObject@sp
  df['Positions'] <- NULL
  p1 <-
    ggplotify::as.grob(
      ggplot2::ggplot() + ggseqlogo::geom_logo(t(df[1:sp,])) + ggseqlogo::theme_logo() +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),axis.text.x =ggplot2::element_text(size=7)) + ggplot2::ylab('') + ggplot2::scale_y_continuous(limits = c(0, 2))
    )
  if (Methylation == TRUE)
  {
    pp1 <- cowplot::plot_grid(BSP1, p1, ncol = 1, align = "v")
    pp1 <- ggplotify::as.grob(pp1)
    grid::pushViewport(vpparent)
    grid::grid.draw(pp1)
  }else{
    p1 <- ggplotify::as.grob(p1)
    grid::pushViewport(vpparent)
    grid::grid.draw(p1)
  }
  
  grid::popViewport()
  
  
  st <- as.numeric(tail(TheObject@forked['PO'],n=1)) - sp
  c <- 1
  for (i in seq((sp+1),dim(df)[1],st)) {
    p2 <-
      ggplot2::ggplot() + ggseqlogo::geom_logo(t(df[i:((st+i)-1),])) + ggseqlogo::theme_logo() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), legend.position = "none",axis.text.x =ggplot2::element_text(size=7)) + ggplot2::ylab('') + ggplot2::scale_y_continuous(limits = c(0, 2))
    p2$scales$scales[[1]] <-
      ggplot2::scale_x_continuous(breaks = seq(1, st , by = 1),
                                  labels = base::as.character(c((sp+1):(sp+st))))
    
    if (Methylation == TRUE)
    {
      p2 <- ggplotify::as.grob(p2)
      
      pp2 <- cowplot::plot_grid(get(paste0("BSp",c)), p2, ncol = 1, align = "v")
      pp2 <- ggplotify::as.grob(pp2)
      grid::pushViewport(get(paste0("vp",c)))
      grid::grid.draw(pp2)
      grid::popViewport()
    }else{
      p2 <- ggplotify::as.grob(p2)
      grid::pushViewport(get(paste0("vp",c)))
      grid::grid.draw(p2)
      grid::popViewport()
    }
    c = c+1
  }
  
 
}
