############# Marine Melanoma Gene Signature shiny app (aggregated) ############
## This app creates dynamic visualizations of features of the Karras, et al 2022
## Melanoma Datasets. This app takes as input a list of genes of interest in
## the columns of CSV format file

## load libraries
library("dplyr")
library("tibble")
library("ggplot2")
library("tidyr")
library("GSVA")
library("patchwork")
library("ggpubr")
library("plotly")
library("heatmaply")
library("shiny")
library("shinydashboard")

############################ load datasets #####################################
fp <- paste0("path/to/shiny/counts/")
readRDS(file = paste0(fp, "mal_agg_counts_fixed.rds")) -> df1

##  define a vector of colors to match the aesthetics of the Pozniak 2022 paper
c("#1B9E77", "#A07125", "#B16548", "#8068AE", "#D03792", "#A66753", "#7FA718",
  "#D9AA04", "#BF8B12", "#927132", "#666666") -> manual_color_scale 

############################# User Interface ###################################
ui <- dashboardPage(skin = "green",
        dashboardHeader(title = "Human Data"),
        dashboardSidebar(## start dashboardSidebar
          fileInput(inputId = "file1",
                    label = "Upload CSV file containing genes of interest",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          h5("CSV file parameters", align = "center"),
          radioButtons(inputId = "sep",
                       label = "Separator",
                       choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                       selected = ","),
          radioButtons(inputId = "quote",
                       label = "Quote",
                       choices = c(None = "", "Double Quote" = '"',
                                   "Single Quote" = "'"),
                       selected = '"'),
          checkboxInput(inputId = "header",
                        label = "Header",
                        value = TRUE)),
        ## end dashboardSidebar
        dashboardBody(
          tabsetPanel(type = "tabs",
              tabPanel(title = "Genes of Interest",
                     fluidRow(column(width = 12,
                     h5(paste0("This displays the list of uploaded genes of ",
                               "interest used to calculate the signature score")),
                     DT::dataTableOutput("table1"), align = "left"))),
              tabPanel(title = "Gene Score Ranking",
                     fluidRow(column(width = 12,
                     h5(paste0("This displays cell types from samples ranked by",
                               " highest signature score"), align = "left"),
                     DT::dataTableOutput("table2"), align = "left"))),
              tabPanel(title = "Violin plots",
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays signature scores for each ",
                                       "major cell type across timepoints"), align = "left"),
                                       plotlyOutput("violin1"),
                                       align = "left")),
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays signature scores for each ",
                                       "timepoint"), align = "left"),
                                       plotlyOutput("violin2"),
                                       align = "left")),
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays signature scores combined ",
                                       "across timepoints for each major cell type"), align = "left"),
                                       plotlyOutput("violin3"),
                                       align = "left"))),
              tabPanel(title = "Heatmaps",
                       fluidRow(column(width = 12,
                                       h5(paste0("This heatmap displays the Spearman correlation of genes from ",
                                       "gene signatures across cell types"),
                                       align = "left"),
                                       plotlyOutput("heat1"),
                                       align = "left")),
                       fluidRow(column(width = 12,
                                       h5(paste0("This heatmap displays the Spearman correlation across ",
                                       "combined cell type and timepoint"),
                                       align = "left"),
                                       plotlyOutput("heat2"),
                                       align = "left")),
                        fluidRow(column(width = 12,
                                       h5(paste0("This heatmap displays the Spearman correlation across ",
                                       "genes in the uploaded gene set calculated across all samples"),
                                       align = "left"),
                                       plotlyOutput("heat3"),
                                       align = "left"))),
              tabPanel(title = "Correlation plots",
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays the Spearman correlation between ",
                                       "the two gene signature scores for each sample"),
                                       align = "left"),
                                       plotlyOutput("corr1"),
                                       align = "left")),
                        # fluidRow(column(width = 12,
                        #                h5(paste0("This plot displays the Spearman correlation between ",
                        #                "the two gene signature scores for each time point (cell types combined)"),
                        #                 align = "left"),
                        #                plotlyOutput("corr2"),
                        #                align = "left")),
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays the Spearman correlation between ",
                                       "the two gene signature scores by cell type and time point"),
                                       align = "left"),
                                       plotlyOutput("corr2"),
                                       align = "left")))
          ) ## end tabsetPanel
        ) ## end dashboardBody 
      ) ## end dashboardPage
## end UI

################################## Server ######################################
## Define server logic 
server <- function(input, output) {
 
######################## READ Gene set of Interest table #######################
   output$table1 <- DT::renderDataTable({
       req(input$file1)
       tryCatch({
          read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_signature_1 = names(.)[1],
                          gene_signature_2 = names(.)[2]) -> g0
                },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(g0)
  })

  #################### table of ranked samples ######################
   output$table2 <- DT::renderDataTable({
     req(input$file1)
       tryCatch({
      # Read Gene of Interest table
      read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

      ## coerce genes signature columns to a list with no NA values
      gl <- list()
      for (i in 1:ncol(g1)) {
        g1 %>% 
          dplyr::select(i) %>%
          dplyr::pull() -> gl[[i]]
        gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
      }
      names(gl) <- names(g1)
         
      ## manipulate pan cancer data to apply gene sig score calculation
      df1 %>%
        as_tibble() %>% ## gets rid of row names
        tibble::column_to_rownames(var = "gene") %>%
        as.data.frame() -> df1_rows

      ## run GSVA and filter results
      GSVA::gsva(expr = as.matrix(df1_rows),
                 gset.idx.list = gl) %>%
        t() %>% as.data.frame() %>%
        tibble::rownames_to_column() %>%
        dplyr::rename(cell_name = names(.)[1]) %>%
        dplyr::arrange(desc(gene_sig_score_1)) %>%
        dplyr::mutate(gene_sig_score_1 = round(gene_sig_score_1, 3),
                      gene_sig_score_2 = round(gene_sig_score_2, 3)) %>%
        dplyr::mutate(sample = gsub("^[^_]*.", "",
                           gsub("^[^_]*.", "",
                                gsub("_[^_]+$", "",
                                     gsub("_[^_]+$", "", cell_name))))) %>%
       dplyr::mutate(timepoint = dplyr::case_when(
             grepl("BT", cell_name) ~ "BT",
             grepl("OT", cell_name) ~ "OT")) %>%
       dplyr::mutate(cell_type = dplyr::case_when(
             grepl("AntigenPresentation", cell_name) ~ "Antigen Presentation",
             grepl("Interferon", cell_name) ~ "Interferon Alpha Beta Response",
             grepl("Melanocytic", cell_name) ~ "Melanocytic",
             grepl("Mesenchymallike", cell_name) ~ "Mesenchymal like",
             grepl("Mitochondrial", cell_name) ~ "Mitochondrial low quality",
             grepl("Mitotic", cell_name) ~ "Mitotic",
             grepl("NeuralCrestlike", cell_name) ~ "Neural Crest like",
             grepl("PatientspecificA", cell_name) ~ "Patient specific A",
             grepl("PatientspecificB", cell_name) ~ "Patient specific B",
             grepl("Hypoxia", cell_name) ~ "Stress Hypoxia Response",
             grepl("p53", cell_name) ~ "Stress p53 Response")) %>%
        dplyr::mutate(aggregated_sample_cell_type = gsub("_[^_]+$", "",
                                     gsub("_[^_]+$", "", cell_name))) %>%
        dplyr::select(aggregated_sample_cell_type, sample, cell_type, timepoint,
                      gene_sig_score_1, gene_sig_score_2) %>%
        as.data.frame() -> df2
                },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(df2)
  })

################################# VIOLIN PLOT 1 ################################

  output$violin1 <- renderPlotly({
      req(input$file1)
      tryCatch({
      ## Read Gene of Interest table
      read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

      ## coerce genes signature columns to a list with no NA values
      gl <- list()
      for (i in 1:ncol(g1)) {
        g1 %>% 
          dplyr::select(i) %>%
          dplyr::pull() -> gl[[i]]
        gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
      }
      names(gl) <- names(g1)

      ## manipulate pan cancer data to apply gene sig score calculation
      df1 %>%
        dplyr::as_tibble() %>% ## gets rid of row names
        tibble::column_to_rownames(var = "gene") %>%
        as.data.frame() -> df1_rows

      ## run GSVA and filter results
      GSVA::gsva(expr = as.matrix(df1_rows),
                 gset.idx.list = gl) %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
          dplyr::rename(gene_sig = names(.)[1]) %>%
          tidyr::gather(key = "cell_name", value = "expression", -gene_sig) %>%
          dplyr::mutate(sample = gsub("^[^_]*.", "",
                                      gsub("^[^_]*.", "",
                                           gsub("_[^_]+$", "",
                                                gsub("_[^_]+$", "",
                                                     cell_name))))) %>%
          dplyr::mutate(timepoint = dplyr::case_when(
            grepl("BT", cell_name) ~ "BT",
            grepl("OT", cell_name) ~ "OT")) %>%
          dplyr::mutate(cell_type = dplyr::case_when(
            grepl("Antigen", cell_name) ~ "Antigen Presentation",
            grepl("Interferon", cell_name) ~ "Interferon Alpha Beta Response",
            grepl("Melanocytic", cell_name) ~ "Melanocytic",
            grepl("Mesenchymallike", cell_name) ~ "Mesenchymal like",
            grepl("Mitochondrial", cell_name) ~ "Mitochondrial low quality",
            grepl("Mitotic", cell_name) ~ "Mitotic",
            grepl("NeuralCrestlike", cell_name) ~ "Neural Crest like",
            grepl("PatientspecificA", cell_name) ~ "Patient specific A",
            grepl("PatientspecificB", cell_name) ~ "Patient specific B",
            grepl("Hypoxia", cell_name) ~ "Stress Hypoxia Response",
            grepl("p53", cell_name) ~ "Stress p53 Response")) %>%
          dplyr::mutate(cell_type = factor(cell_type, levels = unique(cell_type))) %>%
          ggplot(aes(x = cell_type, y = expression, fill = cell_type)) +
          geom_violin(alpha = 0.8) +
          geom_point(position = position_jitter(seed = 1, width = 0.2),
                     size = 0.4, alpha = 0.8) +
          theme_classic() +
          xlab("Cell type (timepoints combined)") +
          ylab("Gene set signature score") +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Violin plot of gene signature scores by cell type") +
          scale_fill_manual(values = manual_color_scale) +
          facet_grid(~gene_sig) -> p1
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p1)
   })

################################# VIOLIN PLOT 2 ################################

  output$violin2 <- renderPlotly({
      req(input$file1)
      tryCatch({
      ## Read Gene of Interest table
      read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

      ## coerce genes signature columns to a list with no NA values
      gl <- list()
      for (i in 1:ncol(g1)) {
        g1 %>% 
          dplyr::select(i) %>%
          dplyr::pull() -> gl[[i]]
        gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
      }
      names(gl) <- names(g1)

      ## manipulate pan cancer data to apply gene sig score calculation
      df1 %>%
        as_tibble() %>% ## gets rid of row names
        tibble::column_to_rownames(var = "gene") %>%
        as.data.frame() -> df1_rows

      ## run GSVA and filter results
      GSVA::gsva(expr = as.matrix(df1_rows),
                 gset.idx.list = gl) %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
          dplyr::rename(gene_sig = names(.)[1]) %>%
          tidyr::gather(key = "cell_name", value = "expression", -gene_sig) %>%
          dplyr::mutate(sample = gsub("^[^_]*.", "",
                                      gsub("^[^_]*.", "",
                                           gsub("_[^_]+$", "",
                                                gsub("_[^_]+$", "", cell_name))))) %>%
          dplyr::mutate(timepoint = dplyr::case_when(
            grepl("BT", cell_name) ~ "BT",
            grepl("OT", cell_name) ~ "OT")) %>%
          dplyr::mutate(cell_type = dplyr::case_when(
            grepl("Antigen", cell_name) ~ "Antigen Presentation",
            grepl("Interferon", cell_name) ~ "Interferon Alpha Beta Response",
            grepl("Melanocytic", cell_name) ~ "Melanocytic",
            grepl("Mesenchymallike", cell_name) ~ "Mesenchymal like",
            grepl("Mitochondrial", cell_name) ~ "Mitochondrial low quality",
            grepl("Mitotic", cell_name) ~ "Mitotic",
            grepl("NeuralCrestlike", cell_name) ~ "Neural Crest like",
            grepl("PatientspecificA", cell_name) ~ "Patient specific A",
            grepl("PatientspecificB", cell_name) ~ "Patient specific B",
            grepl("Hypoxia", cell_name) ~ "Stress Hypoxia Response",
            grepl("p53", cell_name) ~ "Stress p53 Response")) %>%
          dplyr::mutate(cell_type = factor(cell_type, levels = unique(cell_type))) %>%
          ggplot(aes(x = timepoint, y = expression, fill = timepoint)) +
          geom_violin(alpha = 0.8) +
          geom_point(position = position_jitter(seed = 1, width = 0.2),
                     size = 0.4, alpha = 0.8) +
          theme_classic() +
          xlab("Timepoint") +
          ylab("Gene set signature score") +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle("Violin plot of gene signature scores by timepoint") +
          facet_grid(rows = vars(gene_sig)) -> p2
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p2)
   })

################################# VIOLIN PLOT 3 ################################

  output$violin3 <- renderPlotly({
      req(input$file1)
      tryCatch({
      ## Read Gene of Interest table
      read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

      ## coerce genes signature columns to a list with no NA values
      gl <- list()
      for (i in 1:ncol(g1)) {
        g1 %>% 
          dplyr::select(i) %>%
          dplyr::pull() -> gl[[i]]
        gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
      }
      names(gl) <- names(g1)

      ## manipulate pan cancer data to apply gene sig score calculation
      df1 %>%
        as_tibble() %>% ## gets rid of row names
        tibble::column_to_rownames(var = "gene") %>%
        as.data.frame() -> df1_rows

      ## run GSVA and filter results
      GSVA::gsva(expr = as.matrix(df1_rows),
                 gset.idx.list = gl) %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
          dplyr::rename(gene_sig = names(.)[1]) %>%
          tidyr::gather(key = "cell_name", value = "expression", -gene_sig) %>%
          dplyr::mutate(sample = gsub("^[^_]*.", "",
                                      gsub("^[^_]*.", "",
                                           gsub("_[^_]+$", "",
                                                gsub("_[^_]+$", "",
                                                     cell_name))))) %>%
          dplyr::mutate(timepoint = dplyr::case_when(
            grepl("BT", cell_name) ~ "BT",
            grepl("OT", cell_name) ~ "OT")) %>%
          dplyr::mutate(cell_type = dplyr::case_when(
            grepl("Antigen", cell_name) ~ "Antigen Presentation",
            grepl("Interferon", cell_name) ~ "Interferon Alpha Beta Response",
            grepl("Melanocytic", cell_name) ~ "Melanocytic",
            grepl("Mesenchymallike", cell_name) ~ "Mesenchymal like",
            grepl("Mitochondrial", cell_name) ~ "Mitochondrial low quality",
            grepl("Mitotic", cell_name) ~ "Mitotic",
            grepl("NeuralCrestlike", cell_name) ~ "Neural Crest like",
            grepl("PatientspecificA", cell_name) ~ "Patient specific A",
            grepl("PatientspecificB", cell_name) ~ "Patient specific B",
            grepl("Hypoxia", cell_name) ~ "Stress Hypoxia Response",
            grepl("p53", cell_name) ~ "Stress p53 Response")) %>%
          dplyr::mutate(cell_type = factor(cell_type, levels = unique(cell_type))) %>%
          ggplot(aes(x = timepoint, y = expression, fill = cell_type)) +
          geom_violin(alpha = 0.8) +
          geom_point(position = position_jitter(seed = 1, width = 0.2),
                     size = 0.4, alpha = 0.8) +
          theme_classic() +
          xlab("Cell type across timepoint") +
          ylab("Gene set signature score") +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle(paste0("Violin plot of gene signature scores by cell type ",
                         "and timepoint")) +
          scale_fill_manual(values = manual_color_scale) +
          facet_grid(rows = vars(gene_sig)) -> p3
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p3)
   })

################################# HEATMAP 1 ################################
  output$heat1 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## Read Gene of Interest table
      read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

      ## coerce genes signature columns to a list with no NA values
      gl <- list()
      for (i in 1:ncol(g1)) {
        g1 %>% 
          dplyr::select(i) %>%
          dplyr::pull() -> gl[[i]]
        gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
      }
      names(gl) <- names(g1)

      ## manipulate pan cancer data to apply gene sig score calculation
      df1 %>%
        dplyr::filter(gene %in% unlist(g1)[!is.na(unlist(g1))]) %>%
        dplyr::select(gene, contains("OT"), contains("BT")) %>%
        tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
        dplyr::mutate(sample = gsub("^[^_]*.", "",
                                    gsub("^[^_]*.", "",
                                         gsub("_[^_]+$", "",
                                              gsub("_[^_]+$", "",
                                                   cell_name))))) %>%
       dplyr::mutate(timepoint = dplyr::case_when(
              grepl("BT", cell_name) ~ "BT",
              grepl("OT", cell_name) ~ "OT")) %>%
       dplyr::mutate(cell_type = dplyr::case_when(
              grepl("AntigenPresentation", cell_name) ~ "Antigen Presentation",
              grepl("Interferon", cell_name) ~ "Interferon Alpha Beta Response",
              grepl("Melanocytic", cell_name) ~ "Melanocytic",
              grepl("Mesenchymallike", cell_name) ~ "Mesenchymal like",
              grepl("Mitochondrial", cell_name) ~ "Mitochondrial low quality",
              grepl("Mitotic", cell_name) ~ "Mitotic",
              grepl("NeuralCrestlike", cell_name) ~ "Neural Crest like",
              grepl("PatientspecificA", cell_name) ~ "Patient specific A",
              grepl("PatientspecificB", cell_name) ~ "Patient specific B",
              grepl("Hypoxia", cell_name) ~ "Stress Hypoxia Response",
              grepl("p53", cell_name) ~ "Stress p53 Response")) %>%
        dplyr::group_by(gene, cell_type) %>%
        dplyr::summarise(mean_exp = mean(expression)) %>%
        tidyr::spread(key = "cell_type", value = "mean_exp") %>%
        tibble::column_to_rownames(var = "gene") %>%
        dplyr::mutate_all(log) %>%
        as.data.frame() -> r1
      is.na(r1) <- sapply(r1, is.infinite)
      r1[is.na(r1)] <- 0

        ## heatmaply heatmap
        heatmaply::heatmaply(
              x = r1,
              main = "Heatmap of log gene counts",
              xlab = "Cell types",
              ylab = "Genes of interest",
              scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                low = "blue",
                high = "firebrick")) -> p4
              },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p4)
   })
  
################################# HEATMAP 1 ################################
  output$heat2 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## Read Gene of Interest table
      read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

      ## coerce genes signature columns to a list with no NA values
      gl <- list()
      for (i in 1:ncol(g1)) {
        g1 %>% 
          dplyr::select(i) %>%
          dplyr::pull() -> gl[[i]]
        gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
      }
      names(gl) <- names(g1)

      ## manipulate pan cancer data to apply gene sig score calculation
      df1 %>%
        dplyr::filter(gene %in% unlist(g1)[!is.na(unlist(g1))]) %>%
        dplyr::select(gene, contains("OT"), contains("BT")) %>%
        tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
        dplyr::mutate(sample = gsub("^[^_]*.", "",
                                    gsub("^[^_]*.", "",
                                         gsub("_[^_]+$", "",
                                              gsub("_[^_]+$", "",
                                                   cell_name))))) %>%
       dplyr::mutate(timepoint = dplyr::case_when(
              grepl("BT", cell_name) ~ "BT",
              grepl("OT", cell_name) ~ "OT")) %>%
       dplyr::mutate(cell_type = dplyr::case_when(
              grepl("AntigenPresentation", cell_name) ~ "Antigen Presentation",
              grepl("Interferon", cell_name) ~ "Interferon Alpha Beta Response",
              grepl("Melanocytic", cell_name) ~ "Melanocytic",
              grepl("Mesenchymallike", cell_name) ~ "Mesenchymal like",
              grepl("Mitochondrial", cell_name) ~ "Mitochondrial low quality",
              grepl("Mitotic", cell_name) ~ "Mitotic",
              grepl("NeuralCrestlike", cell_name) ~ "Neural Crest like",
              grepl("PatientspecificA", cell_name) ~ "Patient specific A",
              grepl("PatientspecificB", cell_name) ~ "Patient specific B",
              grepl("Hypoxia", cell_name) ~ "Stress Hypoxia Response",
              grepl("p53", cell_name) ~ "Stress p53 Response")) %>%
        dplyr::mutate(cell_sample = paste0(cell_type, " ", timepoint)) %>%
        dplyr::group_by(gene, cell_sample) %>%
        dplyr::summarise(mean_exp = mean(expression)) %>%
        tidyr::spread(key = "cell_sample", value = "mean_exp") %>%
        tibble::column_to_rownames(var = "gene") %>%
        dplyr::mutate_all(log) %>%
        as.data.frame() -> r2
        is.na(r2) <- sapply(r2, is.infinite)
        r2[is.na(r2)] <- 0

        ## heatmaply heatmap
        heatmaply::heatmaply(
              x = r2,
              main = "Heatmap of log gene counts",
              xlab = "Cell types by time point",
              ylab = "Genes of interest",
              scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                low = "blue",
                high = "firebrick")) -> p5
              },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p5)
   })
  
################################# HEATMAP 2 ####################################
  output$heat3 <- renderPlotly({
      req(input$file1)
      tryCatch({
      ## Read Gene of Interest table
      read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

      ## coerce genes signature columns to a list with no NA values
      gl <- list()
      for (i in 1:ncol(g1)) {
        g1 %>% 
          dplyr::select(i) %>%
          dplyr::pull() -> gl[[i]]
        gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
      }
      names(gl) <- names(g1)

      ## manipulate pan cancer data to apply gene sig score calculation
      df1 %>%
        dplyr::filter(gene %in% unlist(g1)[!is.na(unlist(g1))]) %>%
        dplyr::select(gene, contains("OT"), contains("BT")) %>%
        tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
        dplyr::mutate(sample = gsub("^[^_]*.", "",
                                    gsub("^[^_]*.", "",
                                         gsub("_[^_]+$", "",
                                              gsub("_[^_]+$", "",
                                                   cell_name))))) %>%
       dplyr::mutate(timepoint = dplyr::case_when(
              grepl("BT", cell_name) ~ "BT",
              grepl("OT", cell_name) ~ "OT")) %>%
       dplyr::mutate(cell_type = dplyr::case_when(
              grepl("AntigenPresentation", cell_name) ~ "Antigen Presentation",
              grepl("Interferon", cell_name) ~ "Interferon Alpha Beta Response",
              grepl("Melanocytic", cell_name) ~ "Melanocytic",
              grepl("Mesenchymallike", cell_name) ~ "Mesenchymal like",
              grepl("Mitochondrial", cell_name) ~ "Mitochondrial low quality",
              grepl("Mitotic", cell_name) ~ "Mitotic",
              grepl("NeuralCrestlike", cell_name) ~ "Neural Crest like",
              grepl("PatientspecificA", cell_name) ~ "Patient specific A",
              grepl("PatientspecificB", cell_name) ~ "Patient specific B",
              grepl("Hypoxia", cell_name) ~ "Stress Hypoxia Response",
              grepl("p53", cell_name) ~ "Stress p53 Response")) %>%
        dplyr::mutate(cell_sample = paste0(cell_type, " ", timepoint)) %>%
        dplyr::group_by(gene, cell_sample) %>%
        dplyr::summarise(mean_exp = mean(expression)) %>%
        tidyr::spread(key = "cell_sample", value = "mean_exp") %>%
        tibble::column_to_rownames(var = "gene") %>%
        dplyr::mutate_all(log) %>%
        as.data.frame() -> r2
        is.na(r2) <- sapply(r2, is.infinite)
        r2[is.na(r2)] <- 0

        ## heatmaply heatmap
        r2 %>%
          t() %>%
          cor() %>%
        heatmaply::heatmaply(
            main = "Correlation Heatmap of log gene counts across samples",
            xlab = "Genes of interest",
            ylab = "Genes of interest",
            scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
              low = "blue",
              high = "firebrick")) -> p6
              },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p6)
   })
  
################################ SCATTER PLOT 1 ################################

  output$corr1 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## Read Gene of Interest table
         read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

        ## coerce genes signature columns to a list with no NA values
        gl <- list()
        for (i in 1:ncol(g1)) {
          g1 %>%
            dplyr::select(i) %>%
            dplyr::pull() -> gl[[i]]
          gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
        }
        names(gl) <- names(g1)

        ## manipulate pan cancer data to apply gene sig score calculation
        df1 %>%
           as_tibble() %>% ## gets rid of row names
           tibble::column_to_rownames(var = "gene") %>%
           as.data.frame() -> df1_rows

        ## run GSVA and filter results
        GSVA::gsva(expr = as.matrix(df1_rows),
                   gset.idx.list = gl) %>%
        t() %>% as.data.frame() %>%
        tibble::rownames_to_column() %>%
        dplyr::rename(sample_name = names(.)[1]) %>%
        dplyr::mutate(timepoint = dplyr::case_when(
          grepl("BT", sample_name) ~ "BT",
          grepl("OT", sample_name) ~ "OT")) %>%
        dplyr::mutate(cell_type = dplyr::case_when(
          grepl("AntigenPresentation", sample_name) ~ "Antigen Presentation",
          grepl("Interferon", sample_name) ~ "Interferon Alpha Beta Response",
          grepl("Melanocytic", sample_name) ~ "Melanocytic",
          grepl("Mesenchymallike", sample_name) ~ "Mesenchymal like",
          grepl("Mitochondrial", sample_name) ~ "Mitochondrial low quality",
          grepl("Mitotic", sample_name) ~ "Mitotic",
          grepl("NeuralCrestlike", sample_name) ~ "Neural Crest like",
          grepl("PatientspecificA", sample_name) ~ "Patient specific A",
          grepl("PatientspecificB", sample_name) ~ "Patient specific B",
          grepl("Hypoxia", sample_name) ~ "Stress Hypoxia Response",
          grepl("p53", sample_name) ~ "Stress p53 Response")) %>%
        dplyr::mutate(stem_name = paste0(cell_type, " ", timepoint)) %>%
        dplyr::mutate(sample = gsub("^[^_]*.", "",
                                    gsub("^[^_]*.", "",
                                         gsub("_[^_]+$", "",
                                              gsub("_[^_]+$", "", sample_name))))) %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(mean_gene1 = mean(gene_sig_score_1),
                         mean_gene2 = mean(gene_sig_score_2)) %>%
        ggplot(aes(x = mean_gene1, y = mean_gene2)) +
        geom_point(size = 1, color = "black") +
        stat_smooth(method = "lm", se = TRUE, fill = "gray", color = "darkgray",
                    formula = y ~ poly(x, 1, raw = TRUE)) +
        ggpubr::stat_cor(method = "spearman") +
        theme_classic() +
        xlab("Signature Score 1") +
        ylab("Signature Score 2") +
        theme(axis.title = element_text(size = 14),
              legend.position = "none") +
        ggtitle(paste0("Correlation between signature scores across all samples")) -> p7
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p7)
   })
  
################################ SCATTER PLOT 1 (original) #####################

  # output$corr1 <- renderPlotly({
  #     req(input$file1)
  #     tryCatch({
  #        ## Read Gene of Interest table
  #        read.csv(input$file1$datapath,
  #              header = input$header,
  #              sep = input$sep,
  #              quote = input$quote) %>%
  #           dplyr::select(1:2) %>%
  #           dplyr::rename(gene_sig_score_1 = names(.)[1],
  #                         gene_sig_score_2 = names(.)[2]) -> g1
  # 
  #       ## coerce genes signature columns to a list with no NA values
  #       gl <- list()
  #       for (i in 1:ncol(g1)) {
  #         g1 %>% 
  #           dplyr::select(i) %>%
  #           dplyr::pull() -> gl[[i]]
  #         gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
  #       }
  #       names(gl) <- names(g1)
  # 
  #       ## manipulate pan cancer data to apply gene sig score calculation
  #       df1 %>%
  #          as_tibble() %>% ## gets rid of row names
  #          tibble::column_to_rownames(var = "gene") %>%
  #          as.data.frame() -> df1_rows
  # 
  #       ## run GSVA and filter results
  #       GSVA::gsva(expr = as.matrix(df1_rows),
  #                  gset.idx.list = gl) %>%
  #       t() %>% as.data.frame() %>%
  #       tibble::rownames_to_column() %>%
  #       dplyr::rename(sample_name = names(.)[1]) %>%
  #       dplyr::mutate(timepoint = dplyr::case_when(
  #         grepl("BT", sample_name) ~ "BT",
  #         grepl("OT", sample_name) ~ "OT")) %>%
  #       dplyr::mutate(cell_type = dplyr::case_when(
  #         grepl("AntigenPresentation", sample_name) ~ "Antigen Presentation",
  #         grepl("Interferon", sample_name) ~ "Interferon Alpha Beta Response",
  #         grepl("Melanocytic", sample_name) ~ "Melanocytic",
  #         grepl("Mesenchymallike", sample_name) ~ "Mesenchymal like",
  #         grepl("Mitochondrial", sample_name) ~ "Mitochondrial low quality",
  #         grepl("Mitotic", sample_name) ~ "Mitotic",
  #         grepl("NeuralCrestlike", sample_name) ~ "Neural Crest like",
  #         grepl("PatientspecificA", sample_name) ~ "Patient specific A",
  #         grepl("PatientspecificB", sample_name) ~ "Patient specific B",
  #         grepl("Hypoxia", sample_name) ~ "Stress Hypoxia Response",
  #         grepl("p53", sample_name) ~ "Stress p53 Response")) %>%
  #       dplyr::mutate(stem_name = paste0(cell_type, " ", timepoint)) %>%
  #       dplyr::filter(!grepl("Patient", cell_type)) %>%
  #         ggplot(aes(x = gene_sig_score_1, y = gene_sig_score_2)) +
  #         geom_point(size = 1, aes(color = stem_name)) +
  #         stat_smooth(method = "lm", se = TRUE, fill = "gray",
  #                     color = "darkgray", formula = y ~ poly(x, 1, raw = TRUE)) +
  #         ggpubr::stat_cor(method = "spearman") +
  #         theme_classic() +
  #         xlab("Signature Score 1") +
  #         ylab("Signature Score 2") +
  #         theme(axis.title = element_text(size = 14),
  #               legend.position = "none") +
  #         ggtitle(paste0("Correlation between signature scores across cell types")) +
  #         facet_grid(cols = vars(cell_type)) -> p6
  #              },
  #     error = function(e) {## return a safeError if a parsing error occurs
  #       stop(safeError(e))
  #     })
  #     return(p6)
  #  })
################################ SCATTER PLOT 2 ################################
  output$corr2 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## Read Gene of Interest table
         read.csv(input$file1$datapath,
               header = input$header,
               sep = input$sep,
               quote = input$quote) %>%
            dplyr::select(1:2) %>%
            dplyr::rename(gene_sig_score_1 = names(.)[1],
                          gene_sig_score_2 = names(.)[2]) -> g1

        ## coerce genes signature columns to a list with no NA values
        gl <- list()
        for (i in 1:ncol(g1)) {
          g1 %>% 
            dplyr::select(i) %>%
            dplyr::pull() -> gl[[i]]
          gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
        }
        names(gl) <- names(g1)

        ## manipulate pan cancer data to apply gene sig score calculation
        df1 %>%
           as_tibble() %>% ## gets rid of row names
           tibble::column_to_rownames(var = "gene") %>%
           as.data.frame() -> df1_rows

        ## run GSVA and filter results
        GSVA::gsva(expr = as.matrix(df1_rows),
                   gset.idx.list = gl) %>%
        t() %>% as.data.frame() %>%
        tibble::rownames_to_column() %>%
        dplyr::rename(sample_name = names(.)[1]) %>%
        dplyr::mutate(timepoint = dplyr::case_when(
          grepl("BT", sample_name) ~ "BT",
          grepl("OT", sample_name) ~ "OT")) %>%
        dplyr::mutate(cell_type = dplyr::case_when(
          grepl("AntigenPresentation", sample_name) ~ "Antigen Presentation",
          grepl("Interferon", sample_name) ~ "Interferon Alpha Beta Response",
          grepl("Melanocytic", sample_name) ~ "Melanocytic",
          grepl("Mesenchymallike", sample_name) ~ "Mesenchymal like",
          grepl("Mitochondrial", sample_name) ~ "Mitochondrial low quality",
          grepl("Mitotic", sample_name) ~ "Mitotic",
          grepl("NeuralCrestlike", sample_name) ~ "Neural Crest like",
          grepl("PatientspecificA", sample_name) ~ "Patient specific A",
          grepl("PatientspecificB", sample_name) ~ "Patient specific B",
          grepl("Hypoxia", sample_name) ~ "Stress Hypoxia Response",
          grepl("p53", sample_name) ~ "Stress p53 Response")) %>%
        dplyr::mutate(stem_name = paste0(cell_type, " ", timepoint)) %>%
        dplyr::filter(!grepl("Patient", cell_type)) %>%
        ggplot(aes(x = gene_sig_score_1, y = gene_sig_score_2)) +
        geom_point(size = 1, color = "black") +
        stat_smooth(method = "lm", se = TRUE, fill = "gray", color = "darkgray",
                    formula = y ~ poly(x, 1, raw = TRUE)) +
        ggpubr::stat_cor(method = "spearman") +
        theme_classic() +
        xlab("Signature Score 1") +
        ylab("Signature Score 2") +
        theme(axis.title = element_text(size = 14),
              legend.position = "none") +
        ggtitle(paste0("Correlation between signature scores across timepoints")) +
        facet_grid(rows = vars(timepoint)) -> p8
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p8)
   })
  
################################ SCATTER PLOT 3 ################################
  # output$corr3 <- renderPlotly({
  #     req(input$file1)
  #     tryCatch({
  #        ## Read Gene of Interest table
  #        read.csv(input$file1$datapath,
  #              header = input$header,
  #              sep = input$sep,
  #              quote = input$quote) %>%
  #           dplyr::select(1:2) %>%
  #           dplyr::rename(gene_sig_score_1 = names(.)[1],
  #                         gene_sig_score_2 = names(.)[2]) -> g1
  # 
  #       ## coerce genes signature columns to a list with no NA values
  #       gl <- list()
  #       for (i in 1:ncol(g1)) {
  #         g1 %>% 
  #           dplyr::select(i) %>%
  #           dplyr::pull() -> gl[[i]]
  #         gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
  #       }
  #       names(gl) <- names(g1)
  # 
  #       ## manipulate pan cancer data to apply gene sig score calculation
  #       df1 %>%
  #          as_tibble() %>% ## gets rid of row names
  #          tibble::column_to_rownames(var = "gene") %>%
  #          as.data.frame() -> df1_rows
  # 
  #       ## run GSVA and filter results
  #       GSVA::gsva(expr = as.matrix(df1_rows),
  #                  gset.idx.list = gl) %>%
  #       t() %>% as.data.frame() %>%
  #       tibble::rownames_to_column() %>%
  #       dplyr::rename(sample_name = names(.)[1]) %>%
  #       dplyr::mutate(timepoint = dplyr::case_when(
  #         grepl("BT", sample_name) ~ "BT",
  #         grepl("OT", sample_name) ~ "OT")) %>%
  #       dplyr::mutate(cell_type = dplyr::case_when(
  #         grepl("AntigenPresentation", sample_name) ~ "Antigen Presentation",
  #         grepl("Interferon", sample_name) ~ "Interferon Alpha Beta Response",
  #         grepl("Melanocytic", sample_name) ~ "Melanocytic",
  #         grepl("Mesenchymallike", sample_name) ~ "Mesenchymal like",
  #         grepl("Mitochondrial", sample_name) ~ "Mitochondrial low quality",
  #         grepl("Mitotic", sample_name) ~ "Mitotic",
  #         grepl("NeuralCrestlike", sample_name) ~ "Neural Crest like",
  #         grepl("PatientspecificA", sample_name) ~ "Patient specific A",
  #         grepl("PatientspecificB", sample_name) ~ "Patient specific B",
  #         grepl("Hypoxia", sample_name) ~ "Stress Hypoxia Response",
  #         grepl("p53", sample_name) ~ "Stress p53 Response")) %>%
  #       dplyr::mutate(stem_name = paste0(cell_type, " ", timepoint)) %>%
  #         dplyr::filter(!grepl("Patient", cell_type)) %>%
  #       ggplot(aes(x = gene_sig_score_1, y =gene_sig_score_2)) +
  #       geom_point(size = 1, aes(color = stem_name)) +
  #       stat_smooth(method = "lm", se = TRUE, fill = "gray", color = "darkgray",
  #                   formula = y ~ poly(x, 1, raw = TRUE)) +
  #       ggpubr::stat_cor(method = "spearman") +
  #       theme_classic() +
  #       xlab("Signature Score 1") +
  #       ylab("Signature Score 2") +
  #       theme(axis.title = element_text(size = 14),
  #             legend.position = "none") +
  #       ggtitle(paste0("Correlation between signature scores across cell types",
  #                      "and timepoints")) +
  #       facet_grid(rows = vars(timepoint),
  #                  cols = vars(cell_type)) -> p8
  #              },
  #     error = function(e) {## return a safeError if a parsing error occurs
  #       stop(safeError(e))
  #     })
  #     return(p8)
  #  })
}

## Run the application 
shinyApp(ui = ui, server = server)
