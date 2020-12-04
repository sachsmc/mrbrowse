library(data.table)
library(ggplot2)
library(openxlsx)

load("ctdf.rda")
load("qtldf.rda")
load("metadata.rda")
load("qtlivres.rda")

server <- function(input, output, session) {


  observe({

    expselct <- input$exposurect
    if(length(expselct) > 0) {
      xtmp <- ctdf[exposure %in% expselct]

      updateSelectizeInput(session, "outcomect",
                         choices = sort(unique(xtmp$outcome)),
                         selected = input$outcomect)

    } else {

      updateSelectizeInput(session, "outcomect",
                           choices = sort(unique(ctdf$outcome)))

    }
  })

  observe({

      outselct <- input$outcomect
      xtmp2 <- ctdf[outcome == outselct]
      if(length(input$exposurect) == 0) {
        updateSelectizeInput(session, "exposurect",
                           choices = sort(unique(xtmp2$exposure)))
      } else {
        updateSelectizeInput(session, "exposurect",
                             choices = sort(unique(xtmp2$exposure)),
                             selected = input$exposurect)
      }

  })

  output$metaTab <- renderDataTable({
    x1 <- data.frame(Info = "Select an exposure")
    if(length(input$exposurect) > 0) {
      #shorte <- sapply(strsplit(input$exposurect, "_", fixed = TRUE), "[[", 1)
      x1 <- metadata[Phenotype %in% input$exposurect]
      x1$Publication <- sprintf("<a href='https://pubmed.ncbi.nlm.nih.gov/%s' target='_blank'>%s</a>", x1$pmid, x1$pmid)
      x1$main_page <- sprintf("<a href='%s' target='_blank'>%s</a>", x1$main_page, x1$main_page)
      x1$download_page <- sprintf("<a href='%s' target='_blank'>%s</a>", x1$download_page, x1$download_page)
    }
    as.data.frame(x1)

  }, escape = FALSE, options = list(searching  = FALSE, paging = FALSE))


  PlotCT <- reactive({

    if(length(input$exposurect) > 0) {
      x <- ctdf[exposure %in% input$exposurect &
                  outcome == input$outcomect]

      ggplot(x[!is.na(method) & nsnp > 1], aes(x = b, xmin = b - 1.96 * se, xmax = b + 1.96 * se,
                                               y = method)) + geom_pointrange() +
        geom_vline(xintercept = 0, linetype = 3) +
        theme_bw() + xlab("Coefficient (95% confidence interval)") +
        facet_wrap(~ exposure) + ggtitle("Combined SNPs")
    }

  })

  output$mrPlotCT <- renderPlot({

    PlotCT()

  })

  output$downloadPlotCT <- downloadHandler(
    filename = function() {
      paste(paste(input$exposurect, collapse = "-"),
            "-effect-on-", input$outcomect, "-combined-snps.pdf", sep = "")
    },
    content = function(file) {
      ggsave(filename = file, plot = PlotCT(),
             width = as.numeric(input$wp1),
             height = as.numeric(input$hp1), units = "cm")
    }
  )

  PlotCT2 <- reactive({


    if(length(input$exposurect) > 0) {
      x <- ctdf[exposure %in% input$exposurect &
                  outcome == input$outcomect]
      x <- x[!is.na(snp) &
          snp != "All - Inverse variance weighted"]
      x[, snp := reorder(snp, b)]
      ggplot(x, aes(x = b, xmin = b - 1.96 * se, xmax = b + 1.96 * se,
                                                                y = snp)) + geom_pointrange() +
        geom_vline(xintercept = 0, linetype = 3) +
        theme_bw() + xlab("Coefficient (95% confidence interval)") +
        facet_wrap(~ exposure) + ggtitle("Individual SNPs")
    }
  })

  output$mrPlotCT2 <- renderPlot({

    PlotCT2()

  })

  output$downloadPlotCT2 <- downloadHandler(
    filename = function() {
      paste(paste(input$exposurect, collapse = "-"),
            "-effect-on-", input$outcomect, "-individual-snps.pdf", sep = "")
    },
    content = function(file) {
      ggsave(filename = file, plot = PlotCT2(),
             width = as.numeric(input$wp2),
             height = as.numeric(input$hp2), units = "cm")
    }
  )

  snpTabCT <- reactive({

    x <- ctdf[outcome == input$outcomect]
    if(!is.null(input$exposurect)) {
      x <- x[exposure %in% input$exposurect]
    }
    x[, snplink := ifelse(!is.na(snp) & snp != "All - Inverse variance weighted",
                          sprintf("<a href='https://www.ncbi.nlm.nih.gov/snp/%s' target='_blank'>%s</a>",
                                  snp, snp), "")]
    #x[, c("snp")] <- NULL
    x

  })

  output$snpTabCT <- renderDataTable({

    as.data.frame(snpTabCT())

  }, escape = FALSE)

  output$downloadTabCT <- downloadHandler(
    filename = function() {
      paste(paste(input$exposurect, collapse = "-"),
            "-effect-on-", input$outcomect, "-table.xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(snpTabCT(), file = file)
    }
  )


  ## qtl tab

  observe({

    expselqtl <- input$exposureqtl
    if(length(expselqtl) > 0) {
      xtmp <- qtldf[protein %in% expselqtl]
      updateSelectizeInput(session, "outcomeqtl",
                           choices = sort(unique(xtmp$outcome)),
                           selected = input$outcomeqtl)

    } else {

      updateSelectizeInput(session, "outcomeqtl",
                           choices = sort(unique(qtldf$outcome)))

    }
  })

  observe({

    outselqtl <- input$outcomeqtl
    xtmp <- qtldf[outcome == outselqtl & mark == input$cispan]
    updateSelectizeInput(session, "exposureqtl",
                         choices = sort(unique(xtmp$protein)),
                         selected = input$exposureqtl)


  })

  observe({

    cispanme <- input$cispan
    xtmp <- qtldf[outcome == input$outcomeqtl & mark == input$cispan]
    updateSelectizeInput(session, "exposureqtl",
                         choices = sort(unique(xtmp$protein)),
                         selected = input$exposureqtl)

  })

  PlotQTL <- reactive({

    if(length(input$exposureqtl) > 0) {
      x <- qtldf[protein %in% input$exposureqtl &
                   outcome == input$outcomeqtl]

      ggplot(x[test_set == input$cispan], aes(x = b, xmin = b - 1.96 * se, xmax = b + 1.96 * se,
                                              y = method)) + geom_pointrange() +
        geom_vline(xintercept = 0, linetype = 3) +
        theme_bw() + xlab("Coefficient (95% confidence interval)") +
        facet_wrap(~ exposure) + ggtitle("Combined SNPs")
    }

  })

  output$mrPlotQTL <- renderPlot({

    PlotQTL()

  })

  output$downloadPlotQTL <- downloadHandler(
    filename = function() {
      paste(paste(input$exposureqtl, collapse = "-"),
            "-effect-on-", input$outcomeqtl, "-combined-snps.pdf", sep = "")
    },
    content = function(file) {
      ggsave(filename = file, plot = PlotQTL(),
             width = as.numeric(input$wp3),
             height = as.numeric(input$hp3), units = "cm")
    }
  )


  PlotQTL2 <- reactive({

    if(length(input$exposureqtl) > 0) {
      x <- qtldf[protein %in% input$exposureqtl &
                   outcome == input$outcomeqtl]
      x <- x[mark == input$cispan]
      x[, snp := reorder(snp, b)]
      ggplot(x, aes(x = b, xmin = b - 1.96 * se, xmax = b + 1.96 * se,
                                          y = snp)) + geom_pointrange() +
        geom_vline(xintercept = 0, linetype = 3) +
        theme_bw() + xlab("Coefficient (95% confidence interval)") +
        facet_wrap(~ exposure) + ggtitle("Individual SNPs")
    }

  })
  output$mrPlotQTL2 <- renderPlot({

    PlotQTL2()

  })

  output$downloadPlotQTL2 <- downloadHandler(
    filename = function() {
      paste(paste(input$exposureqtl, collapse = "-"),
            "-effect-on-", input$outcomeqtl, "-individual-snps.pdf", sep = "")
    },
    content = function(file) {
      ggsave(filename = file, plot = PlotQTL2(),
             width = as.numeric(input$wp4),
             height = as.numeric(input$hp4), units = "cm")
    }
  )

  qtlTabCT <- reactive({


    x <- qtldf[outcome == input$outcomeqtl]
    if(!is.null(input$exposureqtl)) {
      x <- x[protein %in% input$exposureqtl]
    }
    x <- x[is.na(mark) | mark == input$cispan]
    x <- x[is.na(test_set) | test_set == input$cispan]
    x

  })

  output$qtlTabCT <- renderDataTable({

    as.data.frame(qtlTabCT())

  }, escape = FALSE)


  output$downloadqtlTabCT <- downloadHandler(
    filename = function() {
      paste(paste(input$exposureqtl, collapse = "-"), "-effect-on-", input$outcomeqtl, "-table.xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(qtlTabCT(), file = file)
    }
  )


  ivTabQT <- reactive({

    x2 <- qtlivres
    if(!is.null(input$exposureqtl)) {
      x2 <- qtlivres[ProteinCode %in% input$exposureqtl]
    }
    x2[, snplink := sprintf("<a href='https://www.ncbi.nlm.nih.gov/snp/%s' target='_blank'>%s</a>", SNP, SNP)]
    x2

  })


  output$ivTabQT <- renderDataTable({

    as.data.frame(ivTabQT())


  }, escape = FALSE)

  output$downloadivTabQT <- downloadHandler(
    filename = function() {
      paste(paste(input$exposureqtl, collapse = "-"), "-snp-ivs-table.xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(ivTabQT(), file = file)
    }
  )

}
