library("shiny")
library("shinyWidgets")
library("plotly")
library("sansSouci")

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      numericInput("m", "Number of items", 2^3*20, min = 10, max = 1000),
      numericInput("K", "Number of regions", 8, min = 1 , max = 25),
      #      numericInput("nbRegsignal", "Number of regions with signal", 3, min = 1 , max = 10),
      numericRangeInput(
        inputId = "sig_pos", label = "Signal position",
        value = c(1, 40)
      ),
      numericInput("mu", "Signal strengh", 2, min = 1 , max = 10),
      numericInput("alpha", "Target confidence level", 0.1, min = 0, max = 1),
      # numericInput("sss", "Random seed", 1, min = 1, max = 100000),
    ),
    
    mainPanel(
      h3("Upper bounds on the number of true nulls in selection"),
      tableOutput("bound"),
      plotlyOutput("plot", inline = TRUE),
    )
  ),
  hr(),
  
  h3("Credits"),
  br(),
  
  p(strong("Reference:"), 
    span("G. Durand, G. Blanchard, P. Neuvial and E. Roquain (2018)."),
    a(href = "https://hal.archives-ouvertes.fr/hal-01829037v2", "Post hoc false positive control for structured hypotheses."),
    span("In revision for"),
    em("Scandinavian Journal of Statistics.")),

  p(strong("Code:"),
    span("available from the"),
  a(href = "https://github.com/pneuvial/sanssouci", "R package sansSouci")),

  p(strong("Project page:"), 
  a("https://www.math.univ-toulouse.fr/~pneuvial/sanssouci"))
)
server <- function(input, output, session) {
  
  simu <- reactive({
    #set.seed(311)
    m <- input$m
    sig_pos <- input$sig_pos
    K <- input$K
    s <- round(m/K)
    mu <- input$mu
    # set.seed(input$sss)
    signal <- seq(from = sig_pos[1], to = sig_pos[2])
    theta <- rep(0, m)
    theta[signal] <- 1
    mu = theta*mu
    #pvalues <- 1-pnorm(rnorm(m, mean = mu))
    pvalues <- gen.p.values(m = m, mu = mu, rho = 0)
    stat <- qnorm(1 - pvalues)
    ## tree
    dd <- dyadic.from.window.size(m, s, method = 2)
    leaf_list <- dd$leaf_list
    C <- dd$C
    
    list(p.values = pvalues, stat = stat, signal = paste0("H", theta), tree = dd)
  })    
  
  calc_bounds <- reactive({
    alpha <- input$alpha
    m <- input$m
    K <- input$K
    s <- m/K
    
    sim <- simu()
    pvalues <- sim$p.values
    signal <- sim$signal
    
    C <- sim$tree$C
    leaf_list <- sim$tree$leaf_list
    
    d <- event_data("plotly_selected")
    if (is.null(d)) {
      S <-  which(signal == "H1")
    } else {
      S <- as.numeric(d$x)
    }
    
    # Simes
    V_Simes <- length(S) - posthocBySimes(pvalues, S, alpha)
    
    # # Holm-Bonferroni (+tree)    
    # ZL <- zetas.tree(C, leaf_list, zeta.HB, pvalues, alpha = alpha)
    # V_HB <- V.star(S, C, ZL, leaf_list)
    # 
    # tree    
    ZL <- zetas.tree(C, leaf_list, zeta.DKWM, pvalues, alpha = alpha)
    V_tree <- V.star(S, C, ZL, leaf_list)
    
    # partition
    C0 <- C[length(C)]
    ZL <- zetas.tree(C0, leaf_list, zeta.DKWM, pvalues, alpha = alpha)
    V_part <- V.star(S, C0, ZL, leaf_list)
    tab <- data.frame("Oracle" = sum(signal[S] == "H0"),
                      "Simes" = V_Simes,
                      # "Holm (tree)" = V_HB,
                      "Partition" = V_part,
                      "Tree" = V_tree,
                      check.names = FALSE)  
    list(S = S, tab = tab)
  })
  
  output$sel <- renderText({
    bounds <- calc_bounds()
    S <- bounds$S
    paste("Selection:", paste(range(S), collapse=":"))
  })
  
  output$plot <- renderPlotly({
    m <- input$m
    K <- input$K
    s <- m/K
    sim <- simu()
    pvalues <- sim$p.values
    stat <- sim$stat
    signal <- sim$signal
    sig_pos <- input$sig_pos
    
    res <- calc_bounds()
    tab <- res$tab
    S <- res$S
    # ttl <- expression(sprintf("S = %s:%s, V_{Simes}(S) = %s, V[part](S) = %s V[tree](S) = %s",
    #                sig_pos[1],
    #                sig_pos[2],
    #                tab["Simes"],
    #                tab["Partition"],
    #                tab["Tree"]))
    
    reg_pos <- seq(from = s, to = m, by = s) + 0.5
    reg_dat <- data.frame(id = seq(along = reg_pos), pos = reg_pos)
    
    data <- data.frame(
      x = 1:m, 
      stat = stat,
      mlogp = -log(pvalues),
      signal = signal)

    ylim <- range(data$stat[S])
    
    p <- ggplot(data, aes(x = x, y = stat, group = signal)) +
      geom_vline(data = reg_dat, aes(xintercept = pos), 
                 color = "lightgray", linetype = "dashed", size = 1) +
      geom_rect(xmin = min(S), xmax = max(S), 
                ymin = ylim[1], ymax = ylim[2], 
                # ymin = -10, ymax = 10, 
                fill = "blue", alpha = 0.2, col = "blue") +
      geom_point(aes(shape = signal)) +
      theme_classic() + 
      scale_shape_manual(values=c(19, 1)) +
      xlab("Position") + 
      ylab("Statistic") +
      xlim(c(1, m))
    
    ggplotly(p) %>% layout(dragmode = "select")
  })
  
  
  
  output$bound <- renderTable({
    res <- calc_bounds()
    tab <- res$tab
  })
}
shinyApp(ui, server)