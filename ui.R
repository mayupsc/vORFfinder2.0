library(Gviz)
library(ggplot2)
library(ggtree)
library(tidytree)
library(treeio)
library(shiny)
library(shinymaterial)
library(shinyjs)
library(Biostrings)
library(DECIPHER)


jsCode <- 
  "shinyjs.select_material_sidenav_tab = function(tab_id){
$('.shiny-material-side-nav-tab-content').hide();
$('.shiny-material-side-nav-tab-content').trigger('hide');
$('.shiny-material-side-nav-tab-content').trigger('hidden');
$('.shiny-material-side-nav-tab').removeClass('active');
$('#' + tab_id).show();
$('#' + tab_id).trigger('show');
$('#' + tab_id).trigger('shown');
$('#' + tab_id + '_tab_id').addClass('active');
$('#side_nav_tabs_click_info').trigger('click');
}"


material_page(
  
  title = "Virus ORF Alignment",
  include_nav_bar = FALSE,
  #nav_bar_color = "red lighten-2",
  background_color = "#ffffff",
  
  ##-- Enabling shinyjs
  shinyjs::useShinyjs(),
  shinyjs::extendShinyjs(
    text = jsCode,
    functions = c("select_material_sidenav_tab")
  ),
  ##-- Enabling rintrojs
  #introjsUI(),
  ##-- Global CSS
  shiny::includeCSS("www/styles_global.css"),
  ##-- Github button
  #HTML("<script async defer src='https://buttons.github.io/buttons.js'>"),
  ##-- Sidebar
  material_side_nav(
    fixed = TRUE, 
    image_source = "img/huangtan.png",
    
    material_side_nav_tabs(
      side_nav_tabs = c(
        "HOME" = "home",
        "Upload your data" = "data_upload",
        "PhyloTree & ORF Viewer" = "tree_and_viewer",
        "BLAST" = 'blast',
        "Documentation" = 'Documentation'
      ),
      #icons = c("home", "cloud_upload", "arrow_downward", "arrow_downward","import_contacts"),
      color='teal'
      
    ),
    br(),
    shiny::tags$a(
      href = "https://github.com/mayupsc/virus_orf_finder",
      shiny::tags$i(
        class = "fa fa-github", style = 'font-size:30px; color: black; display: list-item; padding-left: 10px; position: fixed; bottom: 70px;'
      )
    )
  ),
  
  ##-- HOME ----
  material_side_nav_tab_content(
    includeCSS('www/css/style.css'),
    side_nav_tab_id = "home",
    br(),
    shiny::tags$h1("Introduction"),
    
    material_row(
      material_column(
        offset = 2,
        width = 9,
        shiny::tags$img(src = "img/workflow.jpg", height="100%", width="100%", align="top")
        
   )
      )
  ),

  ##-- Upload ----
  material_side_nav_tab_content(
    includeCSS('www/css/style.css'),
    side_nav_tab_id = "data_upload",
    br(),
    shiny::tags$h1("data upload & ORFfinder parameters settings"),
    br(),

      material_row(
        material_column(
          width = 5,
          offset = 3,
          material_row(
            HTML("<h6 style='color:#009688;'>Upload your data</p>"),
            ##-- Upload data
            fileInput(inputId = "dataset_user",
                      label = "The dataset must be compressed in .zip format and no space in file name"
                      )
            ),
          br(),
          material_row(
            uiOutput("select_orf_render"),
            color = '#009688'
            ),
          shinyjs::disabled(material_button(input_id = "upload", label = "Next", icon = "play_arrow", color = "#26a69a"))
          ),

        material_column(
          width = 6,
          offset = 3,


          conditionalPanel(condition = "input.upload",
                           br(),
                           material_card(
                             title = "DATA PROCESS",
                             material_row(
                               material_column(
                                 offset = 1,
                               shinycssloaders::withSpinner(
                                 verbatimTextOutput("runCalc_log"),
                                 type = 4,
                                 color = '#009688'
                                 )
                               )
                             )
                             ),
                             shinyjs::disabled(material_button(input_id = "runlog", label = "Go To ORF Viewer", icon = "play_arrow", color = "#26a69a"))
                             )
                           )
        )
    ),

  ##-- orf_database ----

  material_side_nav_tab_content(
    includeCSS('www/css/style.css'),
    side_nav_tab_id = "tree_and_viewer",
    br(),
    shiny::tags$h1("PhyloTree & ORF Viewer"),

    material_card(
      material_row(
        material_column(
          offset = 1,
          width = 3,
          title = "",
          uiOutput("select_node_render")
        ),

        material_column(
          width = 6,
          shinycssloaders::withSpinner(
            uiOutput("subtree_render"),
            type = 1,
            color = '#009688'
          )
        )
      )
    ),
    material_row(
      material_column(
        offset = 1,
        width = 9,
        material_card(
          title = shiny::tags$h3("ORF viewer"),
          shinycssloaders::withSpinner(
            plotOutput("orf_viewer"),
            type = 1,
            color = '#009688'
          )
        )
      )
    ),
    material_row(
      material_column(
        offset = 1,
        width = 9,
        withSpinner(
          DT::dataTableOutput("orf_table"),
          type = 5,
          color = '#009688'
        )
      )
    ),
    material_row(
      material_column(
        offset = 5,
        material_button(input_id = "blastp", label = "Go To BLASTP", icon = "arrow_downward", color = "#26a69a")
      )
    )
  ),

  material_side_nav_tab_content(
    includeCSS('www/css/style.css'),
    side_nav_tab_id = "blast",
    br(),
    shiny::tags$h1("blastp"),

    material_card(
      title = shiny::tags$h3("Phylo Tree"),
      material_row(
        material_column(
          offset = 1,
          width = 4,
          title = "",
          uiOutput("select_virus_orf_render")
        ),

        material_column(
          width = 6,
          shinycssloaders::withSpinner(
            uiOutput('orf_phylo_render'),
            type = 1,
            color = '#009688'
          )
        )
      )
    ),
    br(),
    br(),

    material_row(
      material_card(
        title = shiny::tags$h3("ORF Alignment"),
        material_row(
          material_column(
            offset = 1,
            width = 10,
            shinycssloaders::withSpinner(
              uiOutput('orf_msa'),
              type = 1,
              color = '#009688'
          )
          )
        )
      )
    )
  ),
  
  
  ##-- HOME ----
  material_side_nav_tab_content(
    includeCSS('www/css/style.css'),
    side_nav_tab_id = "Documentation",
    br(),
    shiny::tags$h1("Documentation"),
    br(),
    material_row(
      material_column(
        offset = 2,
        width = 9,
        shiny::tags$body(
          h4("Method"),
          p("The platform was constructed by shiny(Chang, Cheng, Allaire, Xie, & McPherson, 2020),a R package that used to building interactive web apps. 
             Phylogenetic trees embedded in our platform were calculated by phylip(Baum, 1989) and visualized by ggtree (Yu, Smith, Zhu, Guan, & Lam, 2017). 
             NCBI ORFfinder(https://www.ncbi.nlm.nih.gov/orffinder/) was used to identify orfs for each virus while NCBI blast(Altschul, Gish, Miller, Myers, & Lipman, 1990) was used to search similar orfs.",strong("bold"),'text.'),
          shiny::tags$h5("Reference"),
          p("Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2"),
          p("Baum, B. R. (1989). PHYLIP: Phylogeny Inference Package. Version 3.2. Joel Felsenstein . The Quarterly Review of Biology, 64(4), 539–541. https://doi.org/10.1086/416571"),
          p("Chang, W., Cheng, J., Allaire, J., Xie, Y., & McPherson, J. (2020). Package ‘ shiny ’: Web Application Framework for R, 238."),
          p("Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). Ggtree: an R Package for Visualization and Annotation of Phylogenetic Trees With Their Covariates and Other Associated Data. Methods in Ecology and Evolution, 8(1), 28–36. https://doi.org/10.1111/2041-210X.12628")
                      
        ),
        br()
    
                       
                       
        
      )
    )
  ),


  
  ##-- Footer ----
  div(class = "footer",
      div(includeHTML("www/html/footer.html"))
  )
)
