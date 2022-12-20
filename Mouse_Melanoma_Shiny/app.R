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
fp <- paste0("~/Documents/tmp/2022/202210_joanna/pozniak_melanoma/counts/")
readRDS(file = paste0(fp, "NRAS13_all_43k_agg_counts.rds")) -> df1

############################# User Interface ###################################
ui <- dashboardPage(skin = "green",
        dashboardHeader(title = "Mouse Data"),
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
                                       align = "left"))),
              tabPanel(title = "Heatmaps",
                       fluidRow(column(width = 12,
                                       h5(paste0("This heatmap displays the Spearman correlation of genes from ",
                                       "gene signatures across cell types"),
                                       align = "left"),
                                       plotlyOutput("heat1"),
                                       align = "left")) #,
                        # fluidRow(column(width = 12,
                        #                h5(paste0("This heatmap displays the Spearman correlation across ",
                        #                "genes in the uploaded gene set calculated across all samples"),
                        #                align = "left"),
                        #                plotlyOutput("heat3"),
                        #                align = "left"))
                       ),
              tabPanel(title = "Correlation plots",
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays the Spearman correlation between ",
                                       "the two gene signature scores for each sample"),
                                       align = "left"),
                                       plotlyOutput("corr1"),
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
                              gsub("_[^_]+$", "", cell_name))) %>%
        dplyr::mutate(cell_type = dplyr::case_when(
          grepl("B-cell", cell_name) ~ "B-cell",
          grepl("CAF", cell_name) ~ "CAF",
          grepl("DC", cell_name) ~ "DC",
          grepl("EC", cell_name) ~ "EC",
          grepl("Malignant", cell_name) ~ "Malignant",
          grepl("Monocyte", cell_name) ~ "Monocyte/macrophage",
          grepl("Pericyte", cell_name) ~ "Pericyte",
          grepl("NKcell", cell_name) ~ "T/NKcell")) %>%
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
                              gsub("_[^_]+$", "", cell_name))) %>%
          dplyr::mutate(cell_type = dplyr::case_when(
            grepl("B-cell", cell_name) ~ "B-cell",
            grepl("CAF", cell_name) ~ "CAF",
            grepl("DC", cell_name) ~ "DC",
            grepl("EC", cell_name) ~ "EC",
            grepl("Malignant", cell_name) ~ "Malignant",
            grepl("Monocyte", cell_name) ~ "Monocyte/macrophage",
            grepl("Pericyte", cell_name) ~ "Pericyte",
            grepl("NKcell", cell_name) ~ "T/NKcell")) %>%
        ggplot(aes(x = cell_type, y = expression, fill = cell_type)) +
              geom_violin(alpha = 0.8) +
              geom_point(position = position_jitter(seed = 1, width = 0.2),
                         size = 0.4, alpha = 0.8) +
              theme_classic() +
              theme(legend.position = "none",
                    axis.text.x = element_text(angle = 45, hjust = 1)) +
              ggtitle("Violin plot of gene signature scores by cell type") +
          facet_grid(~gene_sig) -> p1
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p1)
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
        tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
        dplyr::mutate(sample = gsub("^[^_]*.", "",
                              gsub("_[^_]+$", "", cell_name))) %>%
        dplyr::mutate(cell_type = dplyr::case_when(
          grepl("B-cell", cell_name) ~ "B-cell",
          grepl("CAF", cell_name) ~ "CAF",
          grepl("DC", cell_name) ~ "DC",
          grepl("EC", cell_name) ~ "EC",
          grepl("Malignant", cell_name) ~ "Malignant",
          grepl("Monocyte", cell_name) ~ "Monocyte/macrophage",
          grepl("Pericyte", cell_name) ~ "Pericyte",
          grepl("NKcell", cell_name) ~ "T/NKcell")) %>%
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
################################# HEATMAP 2 ####################################
  # output$heat3 <- renderPlotly({
  #     req(input$file1)
  #     tryCatch({
  #     ## Read Gene of Interest table
  #     read.csv(input$file1$datapath,
  #              header = input$header,
  #              sep = input$sep,
  #              quote = input$quote) %>%
  #           dplyr::select(1:2) %>%
  #           dplyr::rename(gene_sig_score_1 = names(.)[1],
  #                         gene_sig_score_2 = names(.)[2]) -> g1
  # 
  #     ## coerce genes signature columns to a list with no NA values
  #     gl <- list()
  #     for (i in 1:ncol(g1)) {
  #       g1 %>% 
  #         dplyr::select(i) %>%
  #         dplyr::pull() -> gl[[i]]
  #       gl[[i]] <- gl[[i]][!is.na(gl[[i]])]
  #     }
  #     names(gl) <- names(g1)
  # 
  #     ## manipulate pan cancer data to apply gene sig score calculation
  #     df1 %>%
  #       dplyr::filter(gene %in% unlist(g1)[!is.na(unlist(g1))]) %>%
  #       tidyr::gather(key = "cell_name", value = "expression", -gene) %>%
  #       dplyr::mutate(sample = gsub("^[^_]*.", "",
  #                                   gsub("_[^_]+$", "", cell_name))) %>%
  #       dplyr::mutate(cell_type = dplyr::case_when(
  #         grepl("B-cell", cell_name) ~ "B-cell",
  #         grepl("CAF", cell_name) ~ "CAF",
  #         grepl("DC", cell_name) ~ "DC",
  #         grepl("EC", cell_name) ~ "EC",
  #         grepl("Malignant", cell_name) ~ "Malignant",
  #         grepl("Monocyte", cell_name) ~ "Monocyte/macrophage",
  #         grepl("Pericyte", cell_name) ~ "Pericyte",
  #         grepl("NKcell", cell_name) ~ "T/NKcell")) %>%
  #       dplyr::group_by(gene, cell_type) %>%
  #       dplyr::summarise(mean_exp = mean(expression)) %>%
  #       tidyr::spread(key = "cell_type", value = "mean_exp") %>%
  #       tibble::column_to_rownames(var = "gene") %>%
  #       dplyr::mutate_all(log) %>%
  #       as.data.frame() -> r2
  #     is.na(r2) <- sapply(r2, is.infinite)
  #     r2[is.na(r2)] <- 0
  # 
  #       ## heatmaply heatmap
  #       r2 %>%
  #         t() %>%
  #         cor() %>%
  #       heatmaply::heatmaply(
  #           main = "Correlation Heatmap of log gene counts across samples",
  #           xlab = "Genes of interest",
  #           ylab = "Genes of interest",
  #           scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
  #             low = "blue",
  #             high = "firebrick")) -> p6
  #             },
  #     error = function(e) {## return a safeError if a parsing error occurs
  #       stop(safeError(e))
  #     })
  #     return(p6)
  #  })
  
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
          dplyr::mutate(sample = gsub("^[^_]*.", "",
                                      gsub("_[^_]+$", "", sample_name))) %>%
          dplyr::mutate(cell_type = dplyr::case_when(
            grepl("B-cell", sample_name) ~ "B-cell",
            grepl("CAF", sample_name) ~ "CAF",
            grepl("DC", sample_name) ~ "DC",
            grepl("EC", sample_name) ~ "EC",
            grepl("Malignant", sample_name) ~ "Malignant",
            grepl("Monocyte", sample_name) ~ "Monocyte/macrophage",
            grepl("Pericyte", sample_name) ~ "Pericyte",
            grepl("NKcell", sample_name) ~ "T/NKcell")) %>%
          dplyr::group_by(sample) %>%
          dplyr::summarise(mean_gene1 = mean(gene_sig_score_1),
                           mean_gene2 = mean(gene_sig_score_2)) %>%
          ggplot(aes(x = mean_gene1, y = mean_gene2)) +
          geom_point(size = 1, color = "black") +
          stat_smooth(method = "lm", se = TRUE, fill = "gray", color = "darkgray",
                      formula = y ~ poly(x, 1, raw = TRUE)) +
          ggpubr::stat_cor(method = "spearman", label.x = 0, label.y = 0.5) +
          theme_classic() +
          xlab("Signature Score 1") +
          ylab("Signature Score 2") +
          theme(axis.title = element_text(size = 14),
                legend.position = "none") +
        ggtitle(paste0("Correlation between signature scores across all samples")
                ) -> p7
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p7)
   })
}

## Run the application 
shinyApp(ui = ui, server = server)
