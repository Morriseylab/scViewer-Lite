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

shinyServer(function(input, output, session) {
  ###################################################
  ###################################################
  ####### Display project list and load data  #######
  ###################################################
  ###################################################
 
  #Get Project list and populate drop-down
  output$projectlist = renderUI({
    excel=read.csv("data/param.csv")
    prj=as.character(excel$projects)
    prj=sort(as.character(prj))
    selectInput("projects","Select a project",as.list(prj),selected = prj[2])
  })

    #Load Rdata
  fileload <- reactive({
    #profvis({
     withProgress(session = session, message = 'Loading Data',detail = 'Please Wait...',{
    #   for (i in 1:10) {
    #     incProgress(1/10)
    #     Sys.sleep(0.1)
    #   }
        projects = input$projects
        inFile = paste('data/',as.character(projects),'.RDS',sep = '')
        scrna=readRDS(inFile)
    return(scrna)
    })#})
  }) %>% bindCache(input$projects)

  #GEt color palette
  getcpalette = reactive({
    scrna = fileload()
    cpallette = scrna@misc$celltypeColorPal
    return(cpallette)
  })
  ###################################################
  ###################################################
  ####### Compare Tsne plot with controls  ##########
  ###################################################
  ###################################################
  output$title <- renderText({
    title = input$projects
    return(title)
    })
  
  
  #Genelist for tsne plot A
  output$gene1aui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(object = scrna)))
     withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
    #     for (i in 1:10) {
    #       incProgress(1/10)
    #       Sys.sleep(0.5)
    #     }
      selectInput('gene1a', label='Gene Name for plot A',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  })
  
  #Genelist for tsne plot B
  output$gene2aui = renderUI({
    scrna=fileload()
    options=sort(rownames(GetAssayData(scrna)))
    withProgress(session = session, message = 'Generating gene list...',detail = 'Please Wait...',{
   #      for (i in 1:10) {
   #        incProgress(1/10)
   #        Sys.sleep(0.5)
   #      }
      selectInput('gene2a', label='Gene Name for plot B',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  }) 
  
  #Based on all use input, generate plots using the right category and dimensionality reduction methods
  comptsne2 = reactive({
    scrna=fileload()
    metadata=as.data.frame(scrna@meta.data)
    met= sapply(metadata,is.numeric)
    tsnea=input$tsnea2
    tsneb=input$tsneb2
    feature=names(met[met==TRUE])
    tsne=names(met[met==FALSE])
    cpallette = getcpalette()
    
    if(input$categorya2 =="clust" ){
      plot1=DimPlot(object = scrna,reduction="umap",group.by = "ident",label = input$checklabel1, pt.size = input$pointa2,label.size = 7, cols=cpallette) + coord_equal() & NoAxes()
    }else if(input$categorya2=="geneexp"){
       plot1=FeaturePlot2(object = scrna,reduction="umap", features = input$gene1a, cols = c("grey", "blue"),pt.size = input$pointa2)  + coord_equal() & NoAxes()
    } 

    if(input$categoryb2 =="clust" ){
      plot2=DimPlot(object = scrna,reduction="umap",group.by = "ident",label = input$checklabel2,pt.size = input$pointa2,label.size = 7, cols=cpallette)  + coord_equal() & NoAxes()
     }else if(input$categoryb2=="geneexp"){
        plot2=FeaturePlot2(object = scrna,reduction="umap", features = input$gene2a, cols = c("grey", "blue"),pt.size = input$pointa2)  + coord_equal() & NoAxes()
     }
    plot1 = plot1+ ggtitle("Plot A") 
    plot2 = plot2+ ggtitle("Plot B") 
    if(input$categorya2=="geneexp" & input$categoryb2!="geneexp"){
      plot3=VlnPlot(object = scrna, features = input$gene1a,pt.size=0,cols=cpallette)+xlab("")
    }else if(input$categorya2!="geneexp" & input$categoryb2=="geneexp"){
      plot3=VlnPlot(object = scrna, features = input$gene2a,pt.size=0,cols=cpallette)+xlab("")
    }else if(input$categorya2=="geneexp" & input$categoryb2=="geneexp"){
      plot3=VlnPlot(object = scrna, features = c(input$gene1a,input$gene2a),pt.size=0,cols=cpallette)+xlab("")
    }else{
      plot3 =plot_spacer()
    }
    #plot3 = plot3+ labs(tag = "C") 
    p=(plot1|plot2)/plot3
    p = p+ plot_layout(heights = c(2, 1))
      p2 <- add_sub(p, paste(input$projects,"_CompareTsne",sep=""), x = 0.87,vpadding = grid::unit(1, "lines"),size=11)
    ggdraw(p2)
  }) #%>% bindCache(input$gene1a, input$gene2a)
  
  #render final plot
  output$comptsne2 = renderPlot({
     withProgress(session = session, message = 'Generating Plot...',detail = 'Please Wait...',{
    #   for (i in 1:10) {
    #     incProgress(1/10)
    #     Sys.sleep(0.5)
    #   }
      comptsne2()
    })
  })
  
  #Handler to download plot
  output$downloadtsneplot <- downloadHandler(
    filename = function() {
      paste0(input$projects,"_CompareTsne.pdf",sep="")
    },
    content = function(file){
      pdf(file,width=14,height = 8,useDingbats=FALSE)
      plot(comptsne2())
      dev.off()
    })
})
