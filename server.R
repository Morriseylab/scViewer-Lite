library(shiny)
library(shinyBS)
library(RColorBrewer)
library(Biobase)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(Seurat)
library(cowplot)
library(data.table)
library(tibble)
library(scExtras)
library(shinythemes)
library(patchwork)

cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7",
            "#8B4484", "#D3D93E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
            "#8A7C64", "#599861","#000099","#FFCC66","#99CC33","#CC99CC","#666666", "#695F74","#0447F9",
            "#89134F","#2CF7F0","#F72C35","#A5B617","#B05927","#B78ED8")

shinyServer(function(input, output, session) {
  ###################################################
  ###################################################
  ####### Display project list and load data  #######
  ###################################################
  ###################################################

  #Load Rdata
  fileload <- reactive({
    #profvis({
    withProgress(session = session, message = 'Loading Data',detail = 'Please Wait...',{
      projects = input$projects
      inFile = 'data/Reference_core.RDS'
      
      scrna=readRDS(inFile)
      return(scrna)
    })#})
  }) %>% bindCache(input$projects)


  ###################################################
  ###################################################
  ####### Compare Tsne plot with controls  ##########
  ###################################################
  ###################################################
  output$title <- renderText({
    title = "Human Reference"
    return(title)
    })
  
  
  #Genelist for tsne plot A
  output$gene1aui = renderUI({
    scrna=fileload()
    #options=sort(rownames(GetAssayData(object = scrna)))
    features=scrna@assays$RNA@meta.features
    options= sort(features$symbol)
     withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('gene1a', label='Gene Name for plot A',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  })
  
  #Genelist for tsne plot B
  output$gene2aui = renderUI({
    scrna=fileload()
    #options=sort(rownames(GetAssayData(scrna)))
    features=scrna@assays$RNA@meta.features
    options= sort(features$symbol)
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
      selectInput('gene2a', label='Gene Name for plot B',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  }) 
  
  #Conditional panel. When highlight cells is chosen, if category is 'cluster', generate drop-down with cluster names and if
  # category is 'geneexp', display error that says cannot highlight
  output$subsaui = renderUI({
    scrna=fileload()
    clusts=levels(Idents(object=scrna))
    if(input$categorya2=="clust"){
      selectInput("selclust","Select a Cluster",clusts)
    }else if(input$categorya2=="geneexp"){
      validate(need(input$categorya2!="geneexp","Cannot subselect gene expression values"))
    }
  })
  
  output$subsbui = renderUI({
    scrna=fileload()
    clusts=levels(Idents(object=scrna))
    if(input$categoryb2=="clust"){
      selectInput("selclustb","Select a Cluster",clusts)
    }else if(input$categoryb2=="geneexp"){
      validate(need(input$categoryb2!="geneexp","Cannot subselect gene expression values"))
    }
  })
  
  #Dimensionality reduction options for left plot
  output$umapa = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapa","Dimensionality Reduction",dimr,selected = "umap")})
  })
  
  #Dimensionality reduction options for right plot
  output$umapb = renderUI({
    scrna=fileload()
    dimr=names(scrna@reductions)
    withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
      selectInput("umapb","Dimensionality Reduction",dimr,selected = "umap")})
  })
  #Based on all use input, generate plots using the right category and dimensionality reduction methods
  comptsne2 = reactive({
    scrna=fileload()
    features=scrna@assays$RNA@meta.features
    gene1a = input$gene1a
    gene1a = features$id[features$symbol==gene1a]
    gene2a = input$gene2a
    gene2a = features$id[features$symbol==gene2a]
    if(input$categorya2 =="clust" & input$subsa==F){
      plot1=DimPlot(object = scrna,reduction=input$umapa,group.by = "ident",label = input$checklabel1,  pt.size = input$pointa2,label.size = 7, cols=cpallette, repel=T)+ theme(legend.position = 'bottom')+ coord_equal() & NoAxes()
    }else if(input$categorya2 =="clust" & input$subsa==TRUE){
      cells=names(Idents(object=scrna)[Idents(object=scrna)==input$selclust])
      plot1=DimPlot(object = scrna,reduction=input$umapa,cells.highlight=cells,cols.highlight="red", cols="grey",group.by = "ident",label = F,  pt.size = input$pointa2, repel=T)+ theme(legend.position = 'bottom')+ 
        scale_color_manual(labels = c("Unselected", input$selclust), values = c("grey", "red")) + coord_equal() & NoAxes()
    }else if(input$categorya2=="geneexp"){
      plot1=FeaturePlot2(object = scrna,reduction=input$umapa, features = gene1a, cols = c("grey", "blue"),pt.size = input$pointa2)+ coord_equal() & NoAxes()
    }
    
    if(input$categoryb2 =="clust" & input$subsb==F){
      plot2=DimPlot(object = scrna,reduction=input$umapb,group.by = "ident",label = input$checklabel2, pt.size = input$pointa2,label.size = 7, cols=cpallette, repel=T)+ theme(legend.position = 'bottom') + coord_equal() & NoAxes()
    }else if(input$categoryb2 =="clust" & input$subsb==TRUE){
      cells=names(Idents(object=scrna)[Idents(object=scrna)==input$selclustb])
      plot2=DimPlot(object = scrna,reduction=input$umapb,cells.highlight=cells, cols.highlight="red",group.by = "ident",label = F,  pt.size = input$pointa2, cols="grey", repel=T)+ theme(legend.position = 'bottom') +
        scale_color_manual(labels = c("Unselected", input$selclustb), values = c("grey", "red"))+ coord_equal() & NoAxes()
    }else if(input$categoryb2=="geneexp"){
      plot2=FeaturePlot2(object = scrna,reduction=input$umapb, features = gene2a, cols = c("grey", "blue"),pt.size = input$pointa2) + coord_equal() & NoAxes()
    }
    plot1 = plot1+ ggtitle("Plot A") 
    plot2 = plot2+ ggtitle("Plot B") 
    if(input$categorya2=="geneexp" & input$categoryb2!="geneexp"){
      plot3=VlnPlot(object = scrna, features = gene1a,pt.size=0,cols=cpallette)+xlab("")+ theme(legend.position = 'none')
    }else if(input$categorya2!="geneexp" & input$categoryb2=="geneexp"){
      plot3=VlnPlot(object = scrna, features = gene2a,pt.size=0,cols=cpallette)+xlab("")+ theme(legend.position = 'none')
    }else if(input$categorya2=="geneexp" & input$categoryb2=="geneexp"){
      plot3=VlnPlot(object = scrna, features = c(gene1a,gene2a),pt.size=0,cols=cpallette)+xlab("")+ theme(legend.position = 'none')
    }else{
      plot3 =plot_spacer()
    }
    #plot3 = plot3+ labs(tag = "C") 
    p=(plot1|plot2)/plot3
    p = p+ plot_layout(heights = c(2, 1))
    ggdraw(p)
  }) #%>% bindCache(input$gene1a, input$gene2a)
  
  #render final plot
  output$comptsne2 = renderPlot({
     withProgress(session = session, message = 'Generating Plot...',detail = 'Please Wait...',{
      comptsne2()
    })
  })
  
  #Handler to download plot
  output$downloadtsneplot <- downloadHandler(
    filename = function() {
      paste0("LungMapref","_CompareTsne.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=14,height = 8,useDingbats=FALSE)
      plot(comptsne2())
      dev.off()
    })
})
