## memeplotr Shiny application -------------------------------------------------
## Launched by memeplotr::run_memeplotr(). Every panel is a thin wrapper around
## the package's plotting functions, and the "R code" tab prints the exact call
## that reproduces whatever is on screen, so the app is a way into the script
## rather than a dead end.

library(shiny)
library(memeplotr)
library(ggplot2)

.has <- function(p) requireNamespace(p, quietly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
ex_files <- list(
  protein = c(xml = "example_meme.xml",     nwk = "example_tree.nwk"),
  dna     = c(xml = "example_meme_dna.xml", nwk = "example_tree_dna.nwk"))
ex_path <- function(which, what) {
  f <- ex_files[[which]][[what]]
  p <- system.file("extdata", f, package = "memeplotr")
  # An empty path means this file is not in the installed package -- almost
  # always a stale or partial install, where app.R is newer than inst/extdata.
  if (!nzchar(p)) {
    stop(sprintf(paste0("Example file '%s' is missing from the installed ",
                        "memeplotr (version %s, library %s). Reinstall the ",
                        "package: remove.packages(\"memeplotr\") then ",
                        "install.packages(<tarball>, repos = NULL, ",
                        "type = \"source\"), and restart R."),
                 f, as.character(utils::packageVersion("memeplotr")),
                 dirname(system.file(package = "memeplotr"))), call. = FALSE)
  }
  p
}

pal_choices    <- c("okabe_ito", "bright", "muted", "pastel", "set3", "earth", "grey")
shape_choices  <- c("arrow", "rect", "roundrect", "hexagon", "capsule", "chevron", "notch")
scheme_choices <- c("auto", "chemistry", "hydrophobicity", "taylor", "nucleotide", "mono")

card_or_div <- function(...) if (.has("bslib")) bslib::card(...) else div(...)

ui <- function() {
  sidebar_body <- list(
    radioButtons("src", "MEME output", inline = TRUE,
                 c("Example" = "example", "Upload" = "upload")),
    conditionalPanel(
      "input.src == 'example'",
      radioButtons("exset", "Example data", inline = TRUE,
                   c("Protein (MEME)" = "protein", "DNA (MEME)" = "dna"))),
    conditionalPanel(
      "input.src == 'upload'",
      fileInput("xml", "MEME / STREME / DREME XML", accept = c(".xml", ".txt")),
      helpText("meme.xml or streme.xml from your run directory.")),
    radioButtons("tsrc", "Phylogeny", inline = TRUE,
                 c("Example" = "example", "Upload" = "upload", "None" = "none")),
    conditionalPanel(
      "input.tsrc == 'upload'",
      fileInput("nwk", "Newick file", accept = c(".nwk", ".newick", ".tre", ".tree", ".txt"))),
    hr(),
    selectInput("order", "Sequence order",
                c("input file" = "input", "tree" = "tree", "motif count" = "count",
                  "sequence length" = "length", "name" = "name")),
    selectInput("source", "Site set",
                c("scanned sites" = "sites", "contributing sites" = "contributing")),
    selectizeInput("motifs", "Motifs (blank = all)", choices = NULL, multiple = TRUE),
    selectizeInput("seqs", "Sequences (blank = all)", choices = NULL, multiple = TRUE),
    sliderInput("pmax", "max -log10(p) filter off at 0", min = 0, max = 20, value = 0, step = 0.5),
    hr(),
    selectInput("shape", "Block shape", shape_choices),
    selectInput("palette", "Palette", pal_choices),
    checkboxInput("custom_cols", "Pick motif colours", FALSE),
    conditionalPanel("input.custom_cols == true", uiOutput("colour_ui")),
    sliderInput("height", "Block height", 0.2, 1, 0.6, 0.05),
    sliderInput("alpha", "Fill opacity", 0.1, 1, 1, 0.05),
    checkboxInput("strand_shape", "Point blocks by strand", FALSE),
    checkboxInput("label", "Label blocks", FALSE),
    checkboxInput("length_label", "Show sequence lengths", FALSE),
    checkboxInput("backbone", "Draw backbone", TRUE),
    hr(),
    checkboxInput("cladogram", "Tree as cladogram", FALSE),
    checkboxInput("align_tips", "Tip guide lines", TRUE),
    sliderInput("tree_w", "Tree width share", 0.2, 2, 1, 0.1),
    hr(),
    textInput("title", "Title", ""),
    textInput("subtitle", "Subtitle", ""),
    numericInput("w", "Export width (in)", 11, 3, 40, 0.5),
    numericInput("h", "Export height (in)", 5.5, 2, 40, 0.5),
    numericInput("dpi", "Export dpi", 300, 72, 1200, 10),
    downloadButton("dl_png", "PNG"),
    downloadButton("dl_pdf", "PDF"),
    downloadButton("dl_svg", "SVG"),
    hr(),
    downloadButton("dl_sites", "Sites (CSV)"),
    downloadButton("dl_motifs", "Motifs (CSV)")
  )

  main <- tabsetPanel(
    id = "tab",
    tabPanel("Motif map", card_or_div(plotOutput("p_map", height = "560px"))),
    tabPanel("Tree + map", card_or_div(plotOutput("p_tree", height = "560px"))),
    tabPanel(
      "Logos",
      fluidRow(
        column(4, selectizeInput("logo_motifs", "Motifs", choices = NULL, multiple = TRUE)),
        column(3, selectInput("logo_method", "Height", c("bits", "probability"))),
        column(3, selectInput("logo_scheme", "Colours", scheme_choices)),
        column(2, checkboxInput("logo_cons", "Consensus axis", TRUE))),
      fluidRow(
        column(3, checkboxInput("logo_rc", "Reverse complement", FALSE)),
        column(3, checkboxInput("logo_ssc", "Small-sample correction", FALSE)),
        column(3, numericInput("logo_ncol", "Columns", 1, 1, 4, 1)),
        column(3, numericInput("logo_outline", "Outline width", 0, 0, 2, 0.1))),
      card_or_div(plotOutput("p_logo", height = "520px"))),
    tabPanel(
      "Summary",
      fluidRow(
        column(3, selectInput("sum_kind", "Plot",
                              c("heatmap", "counts", "co-occurrence", "positions"))),
        column(3, uiOutput("sum_opt1")),
        column(3, uiOutput("sum_opt2")),
        column(3, checkboxInput("sum_label", "Cell labels", FALSE))),
      card_or_div(plotOutput("p_sum", height = "520px"))),
    tabPanel("Tables",
             radioButtons("tbl", NULL, inline = TRUE,
                          c("Sites" = "sites", "Motifs" = "motifs", "Sequences" = "sequences")),
             if (.has("DT")) DT::DTOutput("tbl_out") else tableOutput("tbl_out")),
    tabPanel("R code",
             card_or_div(
               p(class = "text-muted",
                 "Paste this into a script to reproduce the current figure."),
               verbatimTextOutput("code")))
  )

  if (.has("bslib")) {
    bslib::page_sidebar(
      title = "memeplotr - MEME motif figures with ggplot2",
      theme = bslib::bs_theme(version = 5, preset = "shiny"),
      sidebar = do.call(bslib::sidebar, c(list(width = 330, open = "open"), sidebar_body)),
      main)
  } else {
    fluidPage(titlePanel("memeplotr"),
              sidebarLayout(do.call(sidebarPanel, c(list(width = 4), sidebar_body)),
                            mainPanel(width = 8, main)))
  }
}

server <- function(input, output, session) {

  nz <- function(x) if (is.null(x) || !nzchar(x)) NULL else x

  # Uploading your own MEME output almost never pairs with the bundled demo
  # tree, so switch the phylogeny off rather than silently matching 0 tips.
  observeEvent(input$src, {
    updateRadioButtons(session, "tsrc",
                       selected = if (identical(input$src, "upload")) "none" else "example")
  }, ignoreInit = TRUE)

  res <- reactive({
    path <- if (identical(input$src, "upload")) {
      req(input$xml); input$xml$datapath
    } else {
      ex_path(input$exset %||% "protein", "xml")
    }
    out <- tryCatch(read_meme(path), error = function(e) conditionMessage(e))
    if (is.character(out)) {
      # A toast that vanishes after 10s reads as "the app just does not load",
      # so keep it dismissible-only AND stamp the message into every panel.
      msg <- paste("Could not read XML:", out)
      showNotification(msg, type = "error", duration = NULL)
      validate(need(FALSE, msg))
    }
    out
  })

  tree_path <- reactive({
    # The bundled example tree describes the bundled example data only. Pairing
    # it with an uploaded XML matches zero tips, so treat that combination as
    # "no tree" server-side; the radio above is also switched to None.
    if (identical(input$tsrc, "example") && identical(input$src, "upload")) return(NULL)
    switch(input$tsrc,
           example = ex_path(input$exset %||% "protein", "nwk"),
           upload = { req(input$nwk); input$nwk$datapath },
           none = NULL)
  })

  phy <- reactive({
    p <- tree_path(); if (is.null(p)) return(NULL)
    tryCatch(read_tree(p), error = function(e) {
      showNotification(paste("Could not read tree:", conditionMessage(e)),
                       type = "error", duration = 10)
      NULL
    })
  })

  observeEvent(res(), {
    r <- res(); req(r)
    m <- r$motifs$motif
    s <- unique(r$sequences$sequence)
    updateSelectizeInput(session, "motifs", choices = m, selected = character())
    updateSelectizeInput(session, "seqs", choices = s, selected = character())
    updateSelectizeInput(session, "logo_motifs", choices = m,
                         selected = utils::head(m, 3))
  })

  sel_motifs <- reactive(if (length(input$motifs)) input$motifs else NULL)
  sel_seqs   <- reactive(if (length(input$seqs)) input$seqs else NULL)
  pmax_val   <- reactive(if (isTRUE(input$pmax > 0)) 10^(-input$pmax) else NULL)

  output$colour_ui <- renderUI({
    r <- res(); req(r)
    m <- sel_motifs() %||% r$motifs$motif
    base <- motif_palette(length(m), input$palette)
    lapply(seq_along(m), function(i) {
      id <- paste0("col_", i)
      if (.has("colourpicker")) {
        colourpicker::colourInput(id, m[i], value = base[i], showColour = "both")
      } else {
        textInput(id, m[i], value = base[i])
      }
    })
  })

  user_cols <- reactive({
    if (!isTRUE(input$custom_cols)) return(NULL)
    r <- res(); req(r)
    m <- sel_motifs() %||% r$motifs$motif
    v <- vapply(seq_along(m), function(i) input[[paste0("col_", i)]] %||% NA_character_,
                character(1))
    if (anyNA(v)) return(NULL)
    stats::setNames(v, m)
  })

  map_args <- reactive({
    list(source = input$source, order = input$order, shape = input$shape,
         palette = input$palette, height = input$height, alpha = input$alpha,
         strand_shape = input$strand_shape, label = input$label,
         length_label = input$length_label, backbone = input$backbone,
         motifs = sel_motifs(), sequences = sel_seqs(), pvalue_max = pmax_val(),
         title = nz(input$title), subtitle = nz(input$subtitle),
         colours = user_cols())
  })

  the_map <- reactive({
    r <- res(); req(r)
    a <- map_args()
    if (identical(a$order, "tree")) a$tree <- phy()
    do.call(gg_motif_map, c(list(r), a))
  })

  the_tree_map <- reactive({
    r <- res(); req(r); p <- phy()
    if (is.null(p)) return(NULL)
    a <- map_args()
    a$order <- NULL; a$tree <- NULL
    do.call(gg_tree_motifs,
            c(list(r, p), a,
              list(cladogram = input$cladogram, align_tips = input$align_tips,
                   widths = c(input$tree_w, 2.4))))
  })

  output$p_map  <- renderPlot(the_map(), res = 110)
  output$p_tree <- renderPlot({
    p <- the_tree_map()
    validate(need(!is.null(p), "Load a Newick tree to use this panel."))
    p
  }, res = 110)

  logo_args <- reactive({
    list(method = input$logo_method, scheme = input$logo_scheme,
         show_consensus = input$logo_cons, reverse_complement = input$logo_rc,
         small_sample_correction = input$logo_ssc,
         outline_width = input$logo_outline,
         outline = if (isTRUE(input$logo_outline > 0)) "grey20" else NA)
  })

  output$p_logo <- renderPlot({
    r <- res(); req(r)
    m <- if (length(input$logo_motifs)) input$logo_motifs else utils::head(r$motifs$motif, 3)
    do.call(gg_motif_logos,
            c(list(r), list(motifs = m, ncol = input$logo_ncol), logo_args()))
  }, res = 110)

  output$sum_opt1 <- renderUI({
    switch(input$sum_kind,
           heatmap = selectInput("sum_value", "Fill", c("presence", "count")),
           counts = selectInput("sum_position", "Bars", c("stack", "dodge")),
           `co-occurrence` = selectInput("sum_measure", "Measure", c("count", "jaccard")),
           positions = selectInput("sum_geom", "Geom", c("density", "histogram", "jitter")))
  })
  output$sum_opt2 <- renderUI({
    switch(input$sum_kind,
           heatmap = checkboxInput("sum_mcols", "Colour by motif", FALSE),
           counts = NULL,
           `co-occurrence` = checkboxInput("sum_upper", "Full matrix", FALSE),
           positions = selectInput("sum_scale", "Positions", c("absolute", "relative")))
  })

  the_summary <- reactive({
    r <- res(); req(r)
    common <- list(motifs = sel_motifs(), sequences = sel_seqs(),
                   pvalue_max = pmax_val(), title = nz(input$title),
                   subtitle = nz(input$subtitle))
    ord <- list(order = input$order,
                tree = if (identical(input$order, "tree")) phy() else NULL,
                keep_empty = identical(input$order, "tree"))
    switch(
      input$sum_kind,
      heatmap = do.call(gg_motif_heatmap, c(list(r), common, ord, list(
        value = input$sum_value %||% "presence", label = input$sum_label,
        motif_colours = isTRUE(input$sum_mcols), palette = input$palette))),
      counts = do.call(gg_motif_counts, c(list(r), common, ord, list(
        position = input$sum_position %||% "stack", palette = input$palette,
        colours = user_cols()))),
      `co-occurrence` = do.call(gg_motif_cooccurrence, c(list(r), common, list(
        measure = input$sum_measure %||% "count", label = input$sum_label,
        upper = isTRUE(input$sum_upper)))),
      positions = do.call(gg_motif_positions, c(list(r), common, list(
        geom = input$sum_geom %||% "density", scale = input$sum_scale %||% "absolute",
        palette = input$palette, colours = user_cols())))
    )
  })
  output$p_sum <- renderPlot(the_summary(), res = 110)

  current <- reactive({
    switch(input$tab,
           "Motif map" = the_map(),
           "Tree + map" = the_tree_map(),
           "Logos" = {
             r <- res()
             m <- if (length(input$logo_motifs)) input$logo_motifs else utils::head(r$motifs$motif, 3)
             do.call(gg_motif_logos, c(list(r), list(motifs = m, ncol = input$logo_ncol),
                                       logo_args()))
           },
           "Summary" = the_summary(),
           the_map())
  })

  tbl <- reactive({
    r <- res(); req(r)
    d <- switch(input$tbl, sites = r$sites, motifs = r$motifs, sequences = r$sequences)
    as.data.frame(d)
  })
  if (.has("DT")) {
    output$tbl_out <- DT::renderDT(
      tbl(), options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  } else {
    output$tbl_out <- renderTable(utils::head(tbl(), 50))
  }

  save_plot <- function(file, device) {
    ggplot2::ggsave(file, current(), device = device, width = input$w,
                    height = input$h, dpi = input$dpi, limitsize = FALSE)
  }
  output$dl_png <- downloadHandler(
    function() "memeplotr-figure.png", function(f) save_plot(f, "png"))
  output$dl_pdf <- downloadHandler(
    function() "memeplotr-figure.pdf", function(f) save_plot(f, "pdf"))
  output$dl_svg <- downloadHandler(
    function() "memeplotr-figure.svg",
    function(f) {
      if (!.has("svglite")) {
        showNotification("Install the svglite package for SVG export.", type = "error")
        return(NULL)
      }
      save_plot(f, "svg")
    })
  output$dl_sites <- downloadHandler(
    function() "memeplotr-sites.csv",
    function(f) utils::write.csv(as.data.frame(res()$sites), f, row.names = FALSE))
  output$dl_motifs <- downloadHandler(
    function() "memeplotr-motifs.csv",
    function(f) utils::write.csv(as.data.frame(res()$motifs), f, row.names = FALSE))

  fmt <- function(x) {
    if (is.null(x)) return("NULL")
    if (is.character(x)) {
      if (!is.null(names(x))) {
        return(paste0("c(", paste(sprintf('"%s" = "%s"', names(x), x), collapse = ", "), ")"))
      }
      if (length(x) > 1) return(paste0('c("', paste(x, collapse = '", "'), '")'))
      return(sprintf('"%s"', x))
    }
    if (is.logical(x) || is.numeric(x)) {
      if (length(x) > 1) return(paste0("c(", paste(x, collapse = ", "), ")"))
      return(format(x))
    }
    deparse(x)
  }
  call_text <- function(fn, args) {
    args <- args[!vapply(args, is.null, logical(1))]
    paste0(fn, "(\n  x",
           if (length(args)) paste0(",\n  ", paste(sprintf("%s = %s", names(args),
                                                           vapply(args, fmt, character(1))),
                                                   collapse = ",\n  ")),
           "\n)")
  }

  output$code <- renderText({
    r <- res(); req(r)
    xml_src <- if (identical(input$src, "upload") && !is.null(input$xml)) {
      sprintf('"%s"', input$xml$name)
    } else {
      sprintf('system.file("extdata", "%s", package = "memeplotr")',
              ex_files[[input$exset %||% "protein"]][["xml"]])
    }
    nwk_src <- if (identical(input$tsrc, "upload") && !is.null(input$nwk)) {
      sprintf('"%s"', input$nwk$name)
    } else if (identical(input$tsrc, "example")) {
      sprintf('system.file("extdata", "%s", package = "memeplotr")',
              ex_files[[input$exset %||% "protein"]][["nwk"]])
    } else {
      NULL
    }
    head_txt <- paste0(
      "library(memeplotr)\n\nx <- read_meme(", xml_src, ")\n",
      if (!is.null(nwk_src)) paste0("tree <- read_tree(", nwk_src, ")\n") else "")

    a <- map_args()
    a <- a[!vapply(a, is.null, logical(1))]
    body_txt <- switch(
      input$tab,
      "Motif map" = {
        if (identical(input$order, "tree")) a$tree <- quote(tree)
        call_text("gg_motif_map", a)
      },
      "Tree + map" = {
        a$order <- NULL
        paste0("gg_tree_motifs(\n  x, tree,\n  ",
               paste(sprintf("%s = %s", names(a), vapply(a, fmt, character(1))),
                     collapse = ",\n  "),
               sprintf(",\n  cladogram = %s, align_tips = %s, widths = c(%s, 2.4)\n)",
                       input$cladogram, input$align_tips, input$tree_w))
      },
      "Logos" = {
        m <- if (length(input$logo_motifs)) input$logo_motifs else utils::head(r$motifs$motif, 3)
        call_text("gg_motif_logos", c(list(motifs = m, ncol = input$logo_ncol), logo_args()))
      },
      "Summary" = {
        fn <- switch(input$sum_kind, heatmap = "gg_motif_heatmap", counts = "gg_motif_counts",
                     `co-occurrence` = "gg_motif_cooccurrence", positions = "gg_motif_positions")
        sa <- list(motifs = sel_motifs(), sequences = sel_seqs(), pvalue_max = pmax_val(),
                   title = nz(input$title), subtitle = nz(input$subtitle))
        if (input$sum_kind %in% c("heatmap", "counts")) {
          sa$order <- input$order
          if (identical(input$order, "tree")) { sa$tree <- quote(tree); sa$keep_empty <- TRUE }
        }
        if (identical(input$sum_kind, "heatmap")) {
          sa$value <- input$sum_value; sa$label <- input$sum_label
          sa$motif_colours <- isTRUE(input$sum_mcols)
        }
        if (identical(input$sum_kind, "counts")) sa$position <- input$sum_position
        if (identical(input$sum_kind, "co-occurrence")) {
          sa$measure <- input$sum_measure; sa$label <- input$sum_label
          sa$upper <- isTRUE(input$sum_upper)
        }
        if (identical(input$sum_kind, "positions")) {
          sa$geom <- input$sum_geom; sa$scale <- input$sum_scale
        }
        call_text(fn, sa)
      },
      call_text("gg_motif_map", a))
    body_txt <- gsub("quote\\(tree\\)|tree\\(\\)", "tree", body_txt)
    paste0(head_txt, "\n", "p <- ", body_txt, "\n\nggsave(\"figure.png\", p, width = ",
           input$w, ", height = ", input$h, ", dpi = ", input$dpi, ")\n")
  })
}

shinyApp(ui(), server)
