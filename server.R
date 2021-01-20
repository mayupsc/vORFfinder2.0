library(Gviz)
library(ggplot2)
library(ggtree)
library(tidytree)
library(treeio)
library(Biostrings)
library(DECIPHER)

taskID=system("echo $(date \"+%Y%m%d%H%M%S\")",intern = T)

server <- function(input, output, session) {
  shinyjs::runjs("$('#dataset_user').parent().removeClass('btn-default').addClass('btn-danger');")
  
  ##-- upload data
  data_user <- eventReactive(input$dataset_user, {
    print(input$dataset_user)
    input_file <- input$dataset_user$name
    if(is.null(input_file)) {
      return(NULL)
    } else {
      ext <- tools::file_ext(input_file)
      
      if(ext != "zip") {
        sendSweetAlert(
          width = "1000px",
          session = session,
          title = "Error...",
          text = "We can only accept .zip files at the moment.",
          type = "error"
        )
        input_file <- NULL
      } 
    }
    
    return(input_file)
  })
  
  ##-- orf arguments
  
  output$select_orf_render <- renderUI({
    
    orf_argv <- tagList(
      
      material_row(
        
        material_row(
          
          ##-- ORF length
          HTML("<h6 style='color:#009688;'> I.Mininal ORF length</p>"),
          br(),
          material_text_box(
            input_id = "orf_len",
            label = " ", 
            value = 50,
            color = "#9fa8da"
            )
          ),
        br(),
          
          ##-- ORF start codon
        material_row(
          HTML("<h6 style='color:#009688;'> II.start codon</p>"),
          br(),
          radioButtons(
            inputId = "start_codon",
            label = " ",
            choices = c(
              "ATG only" = "0",
              "ATG and alternative initiation codons" = "1",
              "any sense codon" = "2"
            ),
            #color = "#9fa8da",
            selected = "0"

          )
        ),
          
          ##-- Genetic Code
      material_row(
        HTML("<h6 style='color:#009688;'> III.Genetic Code</p>"),
          br(),
          material_dropdown(
            input_id = "genetic_code",
            label = " ",
            choices = c(
              "1.Standard" = "1",
              "2.Vertebrate Mitochondrial" = "2",
              "3.Yeast Mitochondiral" = "3",
              "4.Mold, Protozoan, and Coelenterate Mitochondrial and the Mycoplasma/Spiroplasma" ="4",
              "5.Invertebrate Mitochondrial" = "5",
              "6.Ciliate, Dasycladacean and Hexamita Nuclear" = "6",
              "9.Echinoderm and Flatworm Mitochondrial" ="9",
              "10.Euplotid Nuclear" = "10",
              "11.Bacterial, Archaeal and Plant Plastid" ="11",
              "12.Alternative Yeast Nuclear"="12",
              "13.Ascidian Mitochondrial"="13",
              "14.Alternative Flatworm Mitochondrial"="14",
              "16.Chlorophycean Mitochondrial"="16",
              "21.Trematode Mitochondrial"="21",
              "22.Scenedesmus obliquus Mitochondrial"="22",
              "23.Thraustochytrium Mitochondrial"="23",
              "24.Pterobranchchia Mitochondrial"="24",
              "25.Candidate Division SR1 and Gracilibacteria"="25",
              "26.Pachysolen tannophilus Nuclear"="26",
              "27.Karyorelict Nuclear"="27",
              "28.Condylostoma Nuclear"="28",
              "29.Mesodinium Nuclear"="29",
              "30.Peritrich Nuclear"="30",
              "31.Blastocrithidia Nuclear"="31"
            ),
            selected = "1",
            color = "#9fa8da"
          )
        ),
          br(),
          ##-- Nested ORF
          material_row(
            material_checkbox(input_id = "ignore_nest", 
                            label = shiny::tags$a("IV.Ignore nested ORFs",strong('bold')), 
                            initial_value = FALSE, 
                            color = "#9fa8da")
        )
      )
    )
    return(orf_argv)

  })
  

  observe({
    argv <- data_user()
    if(!is.null(argv)) {
      shinyjs::enable(id = "upload")
    }
  })

  
  runCalc <- eventReactive(input$upload, {
  
    input_file <- data_user()
    cat(paste0('input file is:',input_file,'\n'))
    cat(paste0('orf length is:',input$orf_len,'\n'))
    cat(paste0('start codon is:',input$start_codon,'\n'))
    cat(paste0('genetic code is:',input$genetic_code,'\n'))
    cat(paste0('ignore nest is:',input$ignore_nest,'\n'))
    cmd <- paste('bash script/calc.sh',taskID,input$dataset_user$datapath,input$orf_len,input$start_codon,input$genetic_code,input$ignore_nest,sep = " ")

    
    run_log <- gsub('.*/','',system(cmd,intern=T))
    #run_log2 <- paste0(system(cmd,intern=T),collapse = "\n")
    return(run_log)
    
  })
  
  output$runCalc_log <- renderPrint({
    ##-- Run calc.sh
    run_log <- runCalc()
    writeLines(run_log)
    
  })
  
  
  observe({
    run_log <- runCalc()
    print(run_log)
    if(!is.null(run_log)) {
      shinyjs::enable(id = "runlog")
    }
  })
  
  
  observeEvent(input$runlog, ignoreInit = T, {
    session$sendCustomMessage(type = "shinymaterialJS", shinyjs::js$select_material_sidenav_tab("tree_and_viewer"))
  })
  
  
  observeEvent(input$blastp, ignoreInit = T, {
    session$sendCustomMessage(type = "shinymaterialJS", shinyjs::js$select_material_sidenav_tab("blast"))
  })
  

  output$select_node_render <- renderUI({

    phylip_tree_file <- paste0(taskID,'/database/constree')
    tree <- read.newick(phylip_tree_file)
    #req(input$upload_tree, tree())

    tree_argv <- tagList(
      material_column(
        material_row(
          HTML("<h6 style='color:#009688;'> Select virus </p>"),
          br(),
          selectizeInput(
            inputId = "select_node",
            label = "",
            choices = tree$tip.label,
            width = "100%"
          )
        ),
        material_row(
          HTML("<h6 style='color:#009688;'> Select text size </p>"),
          br(),
          numericInput(
            inputId = "subtree_text_size",
            label = "",
            min = 2,
            value = 4
          )
        ),
        material_row(
          HTML("<h6 style='color:#009688;'>Select plot heigth</p>"),
          br(),
          numericInput(
            inputId = "subtree_plot_height",
            label = "",
            value = 450
          )
        ),
        material_row(
          HTML("<h6 style='color:#009688;'>Select plot width</p>"),
          br(),
          numericInput(
            inputId = "subtree_width_multiply",
            label = "",
            value = 1.2,
            min = 1,
            step = 0.1
          )
        )
      )
    )
    return(tree_argv)
  })
  
  
  
  
  output$select_virus_orf_render <- renderUI({

    req(input$select_node)
    orf_seq <- read.table(paste0(taskID,'/ORF/',input$select_node,'.orfviewer.txt'),sep = "\t",stringsAsFactors = F,header = T)[,1:8]
    #orf_seq <- read.table('ORF/AF206674.1.orfviewer.txt',sep = "\t",stringsAsFactors = F,header = T)
    orf_ID <- unique(gsub(':.[0-9]*:.[0-9]*','',orf_seq$seq.name))
    
    virus_orf_argv <- tagList(
      material_column(
        material_row(
          HTML("<h6 style='color:#009688;'>Select virus orf</p>"),
          br(),
          selectizeInput(
            inputId = "select_virus_orf",
            label = "",
            choices = orf_ID,
            width = "100%"
          )
        ),
        material_row(
          HTML("<h6 style='color:#009688;'>Select text size</p>"),
          br(),
          numericInput(
            inputId = "orfTree_text_size",
            label = "",
            min = 2,
            value = 4
          )
        ),
        material_row(
          HTML("<h6 style='color:#009688;'>Select plot height</p>"),
          br(),
          numericInput(
            inputId = "orfTree_plot_height",
            label = "",
            value = 450
          )
        ),
        material_row(
          HTML("<h6 style='color:#009688;'>Select plot width</p>"),
          br(),
          numericInput(
            inputId = "orfTree_width_multiply",
            label = "",
            value = 1.2,
            min = 1,
            step = 0.1
          )
        )
      )
    )
    return(virus_orf_argv)
  })
  


  output$subtree <- renderPlot({

    req(input$select_node,
        input$subtree_width_multiply,
        input$subtree_text_size,
        input$subtree_plot_height)
    
    phylip_tree_file <- paste0(taskID,'/database/constree')
    tree <- read.newick(phylip_tree_file)

    # getting the subtree phylo or treedata object
    sub_tree <- tree_subset(tree, node = input$select_node, levels_back = length(tree$tip.label))

    # creating the plot
     p <- ggtree(sub_tree ,aes(color = group))+
      geom_tiplab(size = input$subtree_text_size,align = T) +
      #theme_tree2() +
      scale_color_manual(values = c(`1` = "#009688", `0` = "black"))+
      theme(
        legend.position = 'none'
      )

    p + lims(x = c(0, max(p$data$x) * input$subtree_width_multiply))
  })

  # creating the ui element for the subtree
  output$subtree_render <- renderUI({

    plotOutput("subtree", height = input$subtree_plot_height)
  })
  
  
  
  output$orf_viewer <- renderPlot({
    
    req(input$select_node)
    orf_seq <- read.table(paste0(taskID,'/ORF/',input$select_node,'.orfviewer.txt'),sep = "\t",stringsAsFactors = F,header = T)
    #orf_seq <- read.table(paste0('ORF/AF206674.1.orfviewer.txt'),sep = "\t",stringsAsFactors = F,header = T)
    virus_length <- read.table(paste0(taskID,"/database/virus_length.txt"),sep = "\t",stringsAsFactors = F,header = T)
    usr_species_length <- virus_length$length[virus_length$virus==input$select_node]
    axisTrack <- GenomeAxisTrack(range = IRanges(start=1, end = usr_species_length, names="virus genome"))
    aTrack <- AnnotationTrack(start = orf_seq$start2,width = orf_seq$width, chromosome = "chrX", strand = orf_seq$strand,
                              id = orf_seq$symbol, name = "ORF Viewer", transcriptAnnotation = "symbol")

    p <- plotTracks(list(axisTrack, aTrack), from = 1, to = usr_species_length, showId = FALSE,
               fontcolor = "black", add53 = TRUE, add35 = TRUE, cex = 1,
               labelPos = "above", cex.id = 1.2, col.id = "black",
               #panel.only=TRUE,
               fill.range="#b0bec5",
               col="#b0bec5",
               fill="#80cbc4",
               main = input$select_node, featureAnnotation = "id", fontcolor.feature = "black")
    ### ORF positions

  })


  output$orf_table <- DT::renderDataTable({

    orf_seq <- read.table(paste0(taskID,'/ORF/',input$select_node,'.orfviewer.txt'),sep = "\t",stringsAsFactors = F,header = T)
    DT::datatable(
      cbind(' ' = '&oplus;', orf_seq), 
      escape = -0, 
      #extensions = "Buttons",
      options = list(
        pageLength = 10,
        dom = 'Bfrtip',
        columnDefs = list(
          list(visible = FALSE, targets = c(3)),
          list(orderable = FALSE, className = 'details-control', targets = 1)
        )
      ),
      callback = DT::JS("
                    table.column(1).nodes().to$().css({cursor: 'pointer'});
                    var format = function(d) {
                    return '<div style=\"background-color:#eee; padding: .5em;\"> orf_sequence ' +
                    d[3]   + '</div>';
                    };
                    table.on('click', 'td.details-control', function() {
                    var td = $(this), row = table.row(td.closest('tr'));
                    if (row.child.isShown()) {
                    row.child.hide();
                    td.html('&oplus;');
                    } else {
                    row.child(format(row.data())).show();
                    td.html('&CircleMinus;');
                    }
                    });")
        )
    })
  
  ##-- display orf-hit phyloTree
  
  output$orf_phylo <- renderPlot({
    
    #req(input$select_virus_orf,input$orfTree_width_multiply,input$orfTree_text_size,input$orfTree_plot_height)
    
    blast <- read.table(paste0(taskID,"/blastp/",input$select_virus_orf,".blastp"),sep = "\t",stringsAsFactors = F,header = T)
    #blast <- read.table("blastp/",sep = "\t",stringsAsFactors = F,header = T)
    colnames(blast) <- c('queryID','subjectID','identity','Alignment_length','mismatch','gap_opens','q.start','q.end','s.start','s.end','evalue','bit_score')
    
    query_virus <- unique(gsub('ORF.[0-9]*_','',gsub(':.[0-9]*:.[0-9]*','',blast$queryID)))
    hit_virus <- unique(gsub('ORF.[0-9]*_','',gsub(':.[0-9]*:.[0-9]*','',blast$subjectID)))
    
    phylip_tree_file <- paste0(taskID,'/database/constree')
    tree <- read.newick(phylip_tree_file)
    df <- data.frame('virus'=tree$tip.label,'group'='None',stringsAsFactors = F)
    df$group[df$virus %in% hit_virus] <- 'hit.virus'
    df$group[df$virus %in% query_virus] <- 'self'
    groupInfo <- split(df$virus, df$group)
    
    p_tree <- ggtree(tree)+geom_tiplab(size = input$orfTree_text_size,align = T)
    
    groupOTU(p_tree, groupInfo) +
      aes(color=group) +
      scale_color_manual(values = c('self' = "#009688", 'hit.virus' = '#FF6F00','None' = "#757575"))+
      theme(legend.position="right",
            legend.title = element_blank(),
            legend.text = element_text(size = (input$orfTree_text_size+10),face = 'bold')
            #text = element_text(size = input$orfTree_text_size)
            )+ 
      lims(x = c(0, max(p_tree$data$x) * input$orfTree_width_multiply))
    
  })
  
  
  output$orf_phylo_render <- renderUI({
    plotOutput("orf_phylo", height = input$orfTree_plot_height)
  })
  
  
  ##-- display orf blastp results  
  
  output$orf_msa <- renderUI({
    
    protein_seq_msa_file <- paste0(taskID,'/blastp/',input$select_virus_orf,'.msa.fa')
    html_file <- gsub('.msa.fa','.html',protein_seq_msa_file)
    protein_sequences <- readAAStringSet(protein_seq_msa_file, format="fasta")
    
    patterns = c("-", alphabet(protein_sequences))
    BrowseSeqs(protein_sequences, colorPatterns=T,colors = rainbow(length(patterns)),htmlFile = html_file,openURL = F)
    
    return(includeHTML(paste0(taskID,"/blastp/",input$select_virus_orf,".html")))
  })

 ##-- add download button since DT failed to download full results
   output$downloadTSV <- downloadHandler(
    filename = function(){
      paste0(input$select_node,".orf.viewer.txt")
      },
    content = function(con){
      file.copy(paste0(taskID,'/ORF/',input$select_node,".orfviewer.txt"),con)
    }
  )
  
}
