# Purpose: Script that contains Shiny code
# Author: Anzar Alvi
# Date: November 27th, 2025
# Version: 0.1.0
# Bugs and Issues: NA

# References
# Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y,
# Allen J, McPherson J, Dipert A, Borges B (2025). shiny: Web Application
# Framework for R. R package version 1.12.0, https://shiny.posit.co/.
#
# R Core Team (2025). _R: A Language and Environment for Statistical
# Computing_. R Foundation for Statistical Computing, Vienna, Austria.
# https://www.R-project.org/.
#
# Shiny - File Upload. (2014, July 29). Shiny.
# https://shiny.posit.co/r/gallery/widgets/file-upload/
#
# Shiny - Including HTML, text, and Markdown files. (2014, July 29). Shiny.
# https://shiny.posit.co/r/gallery/application-layout/including-html-text-and-markdown-files/
#
# Shiny - Tabsets. (2014, July 30). Shiny.
# https://shiny.posit.co/r/gallery/application-layout/tabsets/

library(shiny)

# Define UI for QCseq app
ui <- fluidPage(

  # App title
  titlePanel("Single-cell RNA-seq QC Analysis and Visualization (QCseq)"),

  # Sidebar layout with input and output definitions
  sidebarLayout(

    # Sidebar panel for inputs
    sidebarPanel(

      h4("Input Data"),

      tags$p("Upload a Seurat object (.rds) containing your scRNA-seq dataset
             with an 'RNA' assay and a counts layer, or use the example dataset
             included with the QCseq package."),

      # Input: Seurat object (.rds)
      fileInput("seurat_file", "Seurat object (.rds)",
                multiple = FALSE,
                accept = c(".rds")),

      # Checkbox: use example data
      checkboxInput("use_example",
                    "Use example dataset 'sample_scRNAseq' from QCseq",
                    value = TRUE),

      tags$p(HTML(
        "The example dataset <code>sample_scRNAseq</code> is stored in
         <code>data/sample_scRNAseq.rda</code> and can be loaded in R with:<br>
         <code>data(\"sample_scRNAseq\", package = \"QCseq\")</code>."
      )),

      br(),

      h4("QC Summary Settings"),

      tags$p("These thresholds control how the QC summary classifies the
             overall dataset quality (poor, moderate, or good)."),

      numericInput("lower_threshold",
                   "Lower QC score threshold (poor ↔ moderate):",
                   value = 0.4, min = 0, max = 1, step = 0.05),

      numericInput("upper_threshold",
                   "Upper QC score threshold (moderate ↔ good):",
                   value = 0.75, min = 0, max = 1, step = 0.05),

      br(),

      h4("Visualization Settings"),

      tags$p("Histogram bar colour can be customized using a hex colour code.
             PCA plot axes correspond to principal components of QC metrics."),

      textInput("bar_colour",
                "Histogram bar colour (hex):",
                value = "#4C72B0"),

      numericInput("pc_x",
                   "Principal component on x-axis (PCx):",
                   value = 1, min = 1, max = 4, step = 1),

      numericInput("pc_y",
                   "Principal component on y-axis (PCy):",
                   value = 2, min = 1, max = 4, step = 1)

    ),

    # Main panel for displaying outputs
    mainPanel(
      # Output: Tabset w/ welcome, references, summary, table, and plots
      tabsetPanel(
        type = "tabs",

        # ---- Welcome tab ----
        tabPanel(
          "Welcome",
          HTML('
    <div style="max-width: 800px; line-height: 1.6;">

      <h3>Welcome to the QCseq Single-cell RNA-seq QC App</h3>

      <h4>About QCseq</h4>

      <p>
        QCseq is an R package for computing quality control (QC) metrics for
        single-cell RNA sequencing (scRNA-seq) data. The QC metrics focus on
        per-cell properties of the gene expression matrix and are computed from
        a Seurat object with an RNA assay and counts layer.
      </p>

      <p>
        The main QC metrics include:
      </p>

      <ul>
        <li><strong>n_genes</strong> – number of detected genes per cell.</li>
        <li><strong>n_transcripts</strong> – total transcripts per cell.</li>
        <li><strong>percent_mit</strong> – percentage of mitochondrial genes.</li>
        <li><strong>percent_ribo</strong> – percentage of ribosomal genes.</li>
      </ul>

      <p>
        These metrics are combined into a unified QC rank score using robust
        median absolute deviation (MAD) scaling and rank-based aggregation.
        The unified score <code>qc_score_rank</code> provides a single index
        of data quality for each cell.
      </p>

      <h4>Example dataset: <code>sample_scRNAseq</code></h4>

      <p>
        QCseq ships with a small example dataset named
        <code>sample_scRNAseq</code>, stored in
        <code>data/sample_scRNAseq.rda</code>. It is a Seurat object that
        contains a subset of scRNA-seq data and is used to illustrate the
        QC workflow.
      </p>

      <p>
        You can load it directly in R with:
      </p>

      <pre><code>data("sample_scRNAseq", package = "QCseq")
qc_results &lt;- QCseq::QCscore(sample_scRNAseq)
head(qc_results)</code></pre>

      <p>
        In this Shiny app, you can use this dataset by checking the option
        <strong>Use example dataset &apos;sample_scRNAseq&apos; from QCseq</strong>
        in the left sidebar.
      </p>

      <h4>How to Use This App</h4>

      <ol>
        <li>
          <strong>Choose input data:</strong>
          Upload a Seurat <code>.rds</code> file or use the example dataset
          by keeping the checkbox selected.
        </li>

        <li>
          <strong>Adjust QC thresholds (optional):</strong>
          The lower and upper thresholds control how the overall dataset
          quality is labeled (poor, moderate, or good).
        </li>

        <li>
          <strong>Customize visualizations:</strong>
          Set the histogram bar colour using a hex code, and choose which
          principal components to display on the PCA plot.
        </li>

        <li>
          Explore the tabs:
          <ul>
            <li><strong>QC summary</strong> – Plain language summary
                generated by <code>QCsummary()</code>.</li>
            <li><strong>QC metrics table</strong> – Per cell metrics and
                QC score from <code>QCscore()</code>.</li>
            <li><strong>QC score histogram</strong> – Distribution of
                <code>qc_score_rank</code> using <code>QChistogram()</code>.</li>
            <li><strong>PCA of QC metrics</strong> – PCA view of QC metrics
                using <code>QCpca()</code>, colored by QC score.</li>
          </ul>
        </li>
      </ol>

      <h4>Notes</h4>

      <ul>
        <li>The input must be a Seurat object with an RNA assay and a counts
            layer, mirroring the requirements of <code>QCcompute()</code>.</li>
        <li>Cells with extreme QC values will be reflected in
            <code>qc_score_rank</code>, which can help identify low-quality
            cells for downstream filtering.</li>
      </ul>

    </div>
          ')
        ),
        # ---- References tab ----
        tabPanel(
          "References",
          HTML('
    <div style="max-width: 800px; line-height: 1.6;">

      <h3>References</h3>

      <ul>
        <li>
          Andrews, T. S., &amp; Hemberg, M. (2018).
          <em>Identifying cell populations with scRNASeq.</em>
          Molecular Aspects of Medicine, 59, 114–122.
          <a href="https://doi.org/10.1016/j.mam.2017.07.002" target="_blank">
            https://doi.org/10.1016/j.mam.2017.07.002
          </a>.
        </li>

        <li>
          Bacher, R., Chu, L.-F., Argus, C., Bolin, J. M., Knight, P., Thomson, J.,
          Knight, P., Thomson, J., Stewart, R., &amp; Kendziorski, C. (2021).
          <em>Enhancing biological signals and detection rates in single-cell RNA-seq
          experiments with cDNA library equalization.</em>
          Nucleic Acids Research, 50(2), e12–e12.
          <a href="https://doi.org/10.1093/nar/gkab1071" target="_blank">
            https://doi.org/10.1093/nar/gkab1071
          </a>.
        </li>

        <li>
          Chang, W., Cheng, J., Allaire, J., Sievert, C., Schloerke, B., Aden-Buie, G.,
          Xie, Y., Allen, J., McPherson, J., Dipert, A., &amp; Borges, B. (2025).
          <em>shiny: Web Application Framework for R.</em>
          R package version 1.12.0.
          <a href="https://shiny.posit.co/" target="_blank">
            https://shiny.posit.co/
          </a>.
        </li>

        <li>
          Hao, Y., Stuart, T. A., Kowalski, M. H., Choudhary, S., Hoffman, P.,
          Hartman, A., Srivastava, A., Molla, G., Shaista Madad, Fernandez-Granda, C.,
          &amp; Rahul atija. (2023).
          <em>Dictionary learning for integrative, multimodal and scalable
          single-cell analysis.</em>
          Nature Biotechnology, 42, 293–304.
          <a href="https://doi.org/10.1038/s41587-023-01767-y" target="_blank">
            https://doi.org/10.1038/s41587-023-01767-y
          </a>.
        </li>

        <li>
          Lu, J., Sheng, Y., Qian, W., Pan, M., Zhao, X., &amp; Ge, Q. (2023).
          <em>scRNA-seq data analysis method to improve analysis performance.</em>
          IET Nanobiotechnology, 17(3), 246–256.
          <a href="https://doi.org/10.1049/nbt2.12115" target="_blank">
            https://doi.org/10.1049/nbt2.12115
          </a>.
        </li>

        <li>
          McCarthy, D. J., Campbell, K. R., Lun, A. T. L., &amp; Wills, Q. F. (2017).
          <em>Scater: pre-processing, quality control, normalisation and visualisation
          of single-cell RNA-seq data in R.</em>
          Bioinformatics, 33, 1179–1186.
          <a href="https://doi.org/10.1093/bioinformatics/btw777" target="_blank">
            https://doi.org/10.1093/bioinformatics/btw777
          </a>.
        </li>

        <li>
          Wickham, H. (2016).
          <em>ggplot2: Elegant Graphics for Data Analysis.</em>
          Springer-Verlag New York.
          <a href="https://ggplot2.tidyverse.org" target="_blank">
            https://ggplot2.tidyverse.org
          </a>.
        </li>

        <li>
          Ziebell, F. (2024).
          <em>RNAseqQC: Quality Control for RNA-Seq Data.</em>
          R package version 0.2.1.
          <a href="https://doi.org/10.32614/CRAN.package.RNAseqQC" target="_blank">
            https://doi.org/10.32614/CRAN.package.RNAseqQC
          </a>.
        </li>

        <li>
          R Core Team. (2025).
          <em>R: A Language and Environment for Statistical Computing.</em>
          Vienna, Austria: R Foundation for Statistical Computing.
          <a href="https://www.R-project.org/" target="_blank">
            https://www.R-project.org/
          </a>.
        </li>

        <li>
          Shiny - File Upload. (2014, July 29).
          <em>Shiny Gallery: File upload example.</em>
          <a href="https://shiny.posit.co/r/gallery/widgets/file-upload/"
             target="_blank">
            https://shiny.posit.co/r/gallery/widgets/file-upload/
          </a>.
        </li>

        <li>
          Shiny - Including HTML, text, and Markdown files. (2014, July 29).
          <em>Shiny Gallery: Including HTML, text, and Markdown files.</em>
          <a href="https://shiny.posit.co/r/gallery/application-layout/including-html-text-and-markdown-files/"
             target="_blank">
            https://shiny.posit.co/r/gallery/application-layout/including-html-text-and-markdown-files/
          </a>.
        </li>

        <li>
          Shiny - Tabsets. (2014, July 30).
          <em>Shiny Gallery: Tabsets.</em>
          <a href="https://shiny.posit.co/r/gallery/application-layout/tabsets/"
             target="_blank">
            https://shiny.posit.co/r/gallery/application-layout/tabsets/
          </a>.
        </li>

        <li>
          OpenAI. (2025).
          <em>ChatGPT. ChatGPT 5; OpenAI.</em>
          <a href="https://chatgpt.com/" target="_blank">
            https://chatgpt.com/
          </a>.
        </li>
      </ul>

    </div>
          ')
        ),

        # ---- QC summary tab ----
        tabPanel(
          "QC summary",
          h3("Overall QC summary"),
          p("Summary of dataset-level QC based on the unified QC score and
             thresholds specified in the sidebar:"),
          br(),
          verbatimTextOutput("qc_summary")
        ),

        # ---- QCcompute tab ----
        tabPanel(
          "QCcompute",
          h3("Independent QC metrics (QCcompute)"),
          p("This table shows the per-cell QC metrics computed by QCseq::QCcompute(), including
             n_genes, n_transcripts, percent_mit, and percent_ribo for each cell."),
          br(),
          tableOutput("qc_compute_table")
        ),

        # ---- QCscore tab ----
        tabPanel(
          "QCscore",
          h3("QC metrics + unified QC score (QCscore)"),
          p("This table shows the output of QCseq::QCscore(), which augments the QCcompute metrics
             with the unified QC rank score qc_score_rank for each cell."),
          br(),
          tableOutput("qc_score_table")
        ),

        # ---- QC histogram tab ----
        tabPanel(
          "QC score histogram",
          p("Distribution of QC rank scores across all cells."),
          plotOutput("qc_histogram")
        ),

        # ---- PCA plot tab ----
        tabPanel(
          "PCA of QC metrics",
          p("PCA scatter plot of QC metrics, colored by QC rank score.
             Axes correspond to the selected principal components (PCx, PCy)."),
          plotOutput("qc_pca_plot")
        )

      )
    )
  )
)

# Define server logic
server <- function(input, output) {

  # ---- Load Seurat object (uploaded or example) ----
  seurat_data <- reactive({
    # Case 1: user uploaded a file
    if (!is.null(input$seurat_file)) {
      obj <- readRDS(input$seurat_file$datapath)

      # Case 2: user wants to use example dataset
    } else if (isTRUE(input$use_example)) {
      data("sample_scRNAseq", package = "QCseq", envir = environment())
      obj <- get("sample_scRNAseq", envir = environment())

      # Case 3: neither example nor upload -> no data
    } else {
      return(NULL)
    }

    # Basic checks: must be Seurat with RNA assay
    req(inherits(obj, "Seurat"))
    req("RNA" %in% Seurat::Assays(obj))

    obj
  })

  # ---- QCcompute: independent QC metrics ----
  qc_compute <- reactive({
    req(seurat_data())
    QCseq::QCcompute(seurat_data())
  })

  # ---- QCscore: QC metrics + unified score ----
  qc_results <- reactive({
    req(seurat_data())
    QCseq::QCscore(seurat_data())
  })

  # ---- QC summary text ----
  output$qc_summary <- renderText({
    req(seurat_data())

    QCseq::QCsummary(
      expr_matrix     = seurat_data(),
      lower_threshold = input$lower_threshold,
      upper_threshold = input$upper_threshold
    )
  })

  # ---- QCcompute table ----
  output$qc_compute_table <- renderTable({
    req(qc_compute())
    qc_compute()
  }, rownames = FALSE)

  # ---- QCscore table ----
  output$qc_score_table <- renderTable({
    req(qc_results())
    qc_results()
  }, rownames = FALSE)

  # ---- Histogram of QC scores ----
  output$qc_histogram <- renderPlot({
    req(seurat_data())
    QCseq::QChistogram(
      expr_matrix = seurat_data(),
      bar_colour  = input$bar_colour
    )
  })

  # ---- PCA plot of QC metrics ----
  output$qc_pca_plot <- renderPlot({
    req(seurat_data())
    QCseq::QCpca(
      expr_matrix = seurat_data(),
      pc_x        = input$pc_x,
      pc_y        = input$pc_y
    )
  })

}

# Create Shiny app
shinyApp(ui, server)

#[END]
